/* Based on https://github.com/wesbarnett/lennardjones" by James W. Barnett
 * Copyright (C) 2015 James W. Barnett <jbarnet4@tulane.edu>
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * The full license is located in a text file titled "LICENSE" in the root
 * directory of the source.
 *
 */

#include "LJ.hh"

System::System(configuration c, int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile) :
conf(c)
{
#ifdef USE_ONEAPI
    q = { sycl::cpu_selector{} };
#endif
}
    
void System::Initialize(configuration c, int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile)
{
    cout << setprecision(6) << fixed << right;
    this->prev_time_point = high_resolution_clock::now();
    this->x.resize(natoms);
    this->v.resize(natoms);
    this->f.resize(natoms);

    this->dt = dt;
    this->halfdt = 0.5*dt;
    this->halfdt2 = 0.5*dt*dt;

    this->natoms = natoms;
    this->inatomsm1 = 1.0/(natoms - 1);
    this->i2natoms = 1.0/(2.0*(double)natoms);
    this->i3natoms = 1.0/(3.0*(double)natoms);

    this->rcut2 = rcut*rcut;
    this->ecut = 4.0 * (1.0/pow(rcut,12) - 1.0/pow(rcut,6));
    this->halfecut = ecut/2.0;
    this->etail = 8.0/3.0 * M_PI * rho * (1.0/(3.0*pow(rcut,9)) - 1.0/pow(rcut,3));
    this->ptail = 16.0/3.0 * M_PI * rho*rho * (2.0/3.0*1.0/pow(rcut,9) - 1.0/pow(rcut,3));

    this->nsample = 0;
    this->rho = rho;
    this->nsteps = nsteps;
    this->xd = xdrfile_open(xtcfile.c_str(), "w");
    this->rhokB = rho*kB;

    // Calculate box dimensions based on density and number of atoms.
    double box_side = pow(natoms/rho,1.0/3.0);
    this->box[X] = box_side;
    this->box[Y] = box_side;
    this->box[Z] = box_side;
    info("Box is {:01.1f} in each dimension.", box_side);

    this->vol = volume(box);
    this->nlist = NeighborList(natoms, rlist);
    this->rdf = Rdf(rdf_nbins, box, rdf_outfile);
    this->tstat = Thermostat(reft, coll_freq, dt);
    this->vel = Velocity(v_nbins, v_max, v_min, v_outfile);

    info("Computing random positions and velocities for {} atoms...", natoms);
    // Draw from a uniform distribution centered at the origin
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> disx(-box[X]/2.0,box[X]/2.0);
    uniform_real_distribution<double> disy(-box[Y]/2.0,box[Y]/2.0);
    uniform_real_distribution<double> disz(-box[Z]/2.0,box[Z]/2.0);
    normal_distribution<double> dis_vel(0.0, sqrt(temp));
    Vector sumv(0.0, 0.0, 0.0);
    double sumv2 = 0.0;
    const double mindist2 = mindist*mindist;
    int i = 0;
    while (i < natoms)
    {

retrypoint:

        this->x[i][X] = disx(gen);
        this->x[i][Y] = disy(gen);
        this->x[i][Z] = disz(gen);

        for (int j = 0; j < i; j++)
        {

            // Too close to to other points?
            if (distance2(x[i], x[j], box) < mindist2)
            {
                if (i > maxtries)
                {
                    cout << "ERROR: Exceeded maximum number of tries in generating initial configuration." << endl;
                }
                goto retrypoint;
            }

        }

        // Point accepted if we're here
        
        this->v[i][X] = dis_vel(gen);
        this->v[i][Y] = dis_vel(gen);
        this->v[i][Z] = dis_vel(gen);

        sumv += this->v[i];
        sumv2 += dot(this->v[i], this->v[i]);

        i++;

    } 
    info("Computed random positions and velocities for {} atoms in {}ms.", natoms, GetTime());
    System::ResetTimer();

    sumv /= this->natoms;
    sumv2 /= this->natoms;
    double fs = sqrt(3.0*temp/sumv2);

    sumv2 = 0.0;
    for (int i = 0; i < this->natoms; i++)
    {
        this->v[i] = (this->v[i] - sumv) * fs;
        sumv2 += dot(this->v[i], this->v[i]);
    }
    this->temp = sumv2 / (3.0 * this->natoms);
    this->ke = 0.5 * sumv2 / this->natoms;

    PdbFile pdb(pdbfile.c_str());
    pdb.write_header(pdbfile, "oneMD LJ Simulator", "First frame");
    for (int i = 0; i < natoms; i++)
    {
        pdb.write_line(i+1, "Ar", "LIG", 1, x[i], 1.00, 0.00);
    }
    pdb.close();
    info("Created .pdb file at {}.", pdbfile);
}

void System::Print(int step)
{
    cout << setw(14) << step;
    cout << setw(14) << (step*this->dt);
    cout << setw(14) << this->temp;
    cout << setw(14) << this->press;
    cout << setw(14) << this->ke;
    cout << setw(14) << this->pe;
    cout << setw(14) << this->pe+this->ke;
    cout << setw(14) << this->GetTime() << "ms" << endl;
    return;
}

void System::PrintHeader()
{
    cout << setw(14) << "Step";
    cout << setw(14) << "Time";
    cout << setw(14) << "Temp";
    cout << setw(14) << "Pressure";
    cout << setw(14) << "KE";
    cout << setw(14) << "PE";
    cout << setw(14) << "Tot. Energy";
    cout << setw(14) << "Clock Time" << endl;
    return;
}

void System::PrintAverages()
{
    cout << "AVERAGES & CONSTANTS (" << this->nsample << " steps sampled out of " << this->nsteps << " total steps)" << endl;
    cout << setw(20) << "Number: " << setw(14) << this->natoms << endl;
    cout << setw(20) << "Density: " << setw(14) << this->rho << endl;
    cout << setw(20) << "Volume: " << setw(14) << this->vol << endl;
    cout << setw(20) << "Temperature: " << setw(14) << this->Temperature.GetAvg() << " +/- " << setw(14) << this->Temperature.GetError() << endl;
    cout << setw(20) << "Pressure: " << setw(14) << this->Pressure.GetAvg() << " +/- " << setw(14) << this->Pressure.GetError() << endl;
    cout << setw(20) << "Kinetic Energy: " << setw(14) << this->KineticEnergy.GetAvg() << " +/- " << setw(14) << this->KineticEnergy.GetError() << endl;
    cout << setw(20) << "Potential Energy: " << setw(14) << this->PotentialEnergy.GetAvg() << " +/- " << setw(14) << this->PotentialEnergy.GetError() << endl;
    cout << setw(20) << "Total Energy: " << setw(14) << this->TotalEnergy.GetAvg() << " +/- " << setw(14) << this->TotalEnergy.GetError() << endl;
    return;
}

void System::SampleRdf()
{
    this->rdf.sample(this->x, this->box);
    return;
}

void System::NormalizeRdf()
{
    this->rdf.normalize(this->natoms, this->box);
    return;
}

void System::OutputRdf()
{
    this->rdf.output();
    return;
}

void System::SampleVel()
{
    this->vel.sample(this->v);
    return;
}

void System::NormalizeVel()
{
    this->vel.normalize(natoms);
    return;
}

void System::OutputVel()
{
    this->vel.output();
    return;
}

void System::Sample()
{
    this->nsample++;
    this->KineticEnergy.Sample(this->ke);
    this->PotentialEnergy.Sample(this->pe);
    this->Pressure.Sample(this->press);
    this->Temperature.Sample(this->temp);
    this->TotalEnergy.Sample(this->ke+this->pe);
    return;
}

void System::WriteXTC(int step)
{
    rvec *x_xtc;
    matrix box_xtc;

    // Convert to "nanometer" (even though we are in reduced units)
    box_xtc[X][X] = this->box[X]*0.1;
    box_xtc[X][Y] = 0.0;
    box_xtc[X][Z] = 0.0;
    box_xtc[Y][X] = 0.0;
    box_xtc[Y][Y] = this->box[Y]*0.1;
    box_xtc[Y][Z] = 0.0;
    box_xtc[Z][X] = 0.0;
    box_xtc[Z][Y] = 0.0;
    box_xtc[Z][Z] = this->box[Z]*0.1;

    x_xtc = new rvec[this->x.size()];
#pragma omp for
#ifdef MSVC
    for (int i = 0; i < x.size(); i++)
#else
    for (unsigned int i = 0; i < x.size(); i++)
#endif
    {
        // Shift all the points to the center of the box
        this->x[i] = pbc(this->x[i], this->box);

        // Convert to "nanometers"
        x_xtc[i][X] = this->x[i][X]*0.1;
        x_xtc[i][Y] = this->x[i][Y]*0.1;
        x_xtc[i][Z] = this->x[i][Z]*0.1;
    }

    write_xtc(this->xd, this->x.size(), step, this->dt*step, box_xtc, x_xtc, 1000);

    return;
}

void System::CloseXTC()
{
    xdrfile_close(this->xd);
}

void System::ErrorAnalysis(int nblocks)
{
    this->KineticEnergy.ErrorAnalysis(nblocks);
    this->PotentialEnergy.ErrorAnalysis(nblocks);
    this->Pressure.ErrorAnalysis(nblocks);
    this->Temperature.ErrorAnalysis(nblocks);
    this->TotalEnergy.ErrorAnalysis(nblocks);
    return;
}

void System::NormalizeAverages()
{
    this->KineticEnergy.Normalize();
    this->PotentialEnergy.Normalize();
    this->Pressure.Normalize();
    this->Temperature.Normalize();
    this->TotalEnergy.Normalize();
    return;
}

int System::GetTime()
{
    auto t = duration_cast<milliseconds>(high_resolution_clock::now() - prev_time_point).count();
    return t;
}

void System::ResetTimer()
{
    info("Resetting timer.");
    System::prev_time_point = high_resolution_clock::now();
}

void System::UpdateNeighborListHostCPU()
{
    this->nlist.UpdateHostCPU(this->x, this->box);
    return;
}

void System::CalcForceHostCPU()
{

    int ncut = 0;
    double pe = 0.0;
    for (int i = 0; i < this->natoms; i++)
    {
        this->f[i] = 0.0;
    }

    #pragma omp parallel
    {

        vector <Vector> f_thread(natoms);
        for (int i = 0; i < this->natoms; i++)
        {
            f_thread[i] = 0.0;
        }

        // Uses neighbor lists to calculate forces and energies. We didn't
        // double count the atoms on the neighbor list, so we have to look at
        // each atom's list. The last atom never has it's own list since it will
        // always be on at least one other atom's list (or it is too far away to
        // interact with any other atom)
        #pragma omp for schedule(guided, CHUNKSIZE) reduction(+:ncut,pe)
        for (int i = 0; i < this->natoms-1; i++)
        {

            for (int neighb = 0; neighb < nlist.GetSize(i); neighb++)
            {

                int j = nlist.GetNeighbor(i, neighb);
                Vector dr = pbc(x[i] - x[j], box);
                double r2 = dot(dr,dr);

                if (r2 <= this->rcut2)
                {

                    double r2i = 1.0/r2;
                    double r6i = pow(r2i,3);
                    Vector fr = 48.0 * r2i * r6i * (r6i - 0.5) * dr;
                    
                    // We have to count the force both on atom i from j and on j
                    // from i, since we didn't double count on the neighbor
                    // lists
                    f_thread[i] += fr;
                    f_thread[j] -= fr;

                    pe += 4.0*r6i*(r6i-1.0) - this->ecut;
                    ncut++;

                }

            }

        }

        #pragma omp critical
        {

            for (int i = 0; i < natoms; i++)
            {
                this->f[i] += f_thread[i];
            }

        }

    }

    double vir = 0.0;
    #pragma omp parallel for schedule(guided, CHUNKSIZE) reduction(+:vir)
    for (int i = 0; i < natoms; i++)
    {
        vir += dot(f[i], x[i]);
    }
    vir *= oneSixth;
    this->press = this->rhokB * this->temp + vir/this->vol + this->ptail;

    ncut *= inatomsm1;
    this->pe = pe/this->natoms + this->etail + halfecut*(double)ncut; 

    return;

}

// Velocity Verlet integrator in two parts
void System::IntegrateHostCPU(int a, bool tcoupl)
{

    if (a == 0) 
    {
        #pragma omp parallel for schedule(guided, CHUNKSIZE)
        for (int i = 0; i < this->natoms; i++)
        {
            this->x[i] += this->v[i]*this->dt + this->f[i]*this->halfdt2;
            this->v[i] += this->f[i]*this->halfdt;
        }
    }
    else if (a == 1)
    {

        double sumv2 = 0.0;

        #pragma omp parallel for schedule(guided, CHUNKSIZE) reduction(+:sumv2)
        for (int i = 0; i < natoms; i++)
        {
            this->v[i] += this->f[i]*this->halfdt;
            sumv2 += dot(this->v[i], this->v[i]);
        }

        if (tcoupl == true)
        {
            tstat.DoCollisions(v);
        }

        this->temp = sumv2 * this->i3natoms;
        this->ke = sumv2 * this->i2natoms;

    }

    return;
}

#ifdef USE_ONEAPI
void System::UpdateNeighborListCPU()
{
    throw new NotImplementedException();
}

void System::CalcForceCPU()
{
    throw new NotImplementedException();
}

void System::IntegrateCPU(int a, bool tcoupl)
{
    throw new NotImplementedException();
}

void System::UpdateNeighborListGPU()
{
    throw new NotImplementedException();
}

void System::CalcForceGPU()
{
    throw new NotImplementedException();
}

void System::IntegrateGPU(int a, bool tcoupl)
{
    throw new NotImplementedException();
}

void System::UpdateNeighborListFPGA()
{
    throw new NotImplementedException();
}

void System::CalcForceFPGA()
{
    throw new NotImplementedException();
}

void System::IntegrateFPGA(int a, bool tcoupl)
{
    throw new NotImplementedException();
}
#endif

LJ::LJ(configuration config, Device device) : 
    Simulator("LJ", config, device),
    conf(config),
    sys(config, config.natoms, config.nsteps, config.rho, config.rcut, config.rlist, config.temp, config.dt, config.mindist, config.maxtries, config.pdbfile, config.reft, config.coll_freq, config.xtcfile, config.rdf_nbins, config.rdf_outfile, config.v_nbins, config.v_max, config.v_min, config.v_outfile)
{
#ifdef USE_ONEAPI
    q = { sycl::cpu_selector{} };
#endif
}

bool LJ::Initialize() 
{
    info("Lennard-Jones simulation of {} atoms in a cubic box of dimension {:01.1f} for {} steps.\nOriginal code by James W. Barnett https://github.com/wesbarnett/lennardjones", 
        conf.natoms, sys.box[0], conf.nsteps);
    if (!conf.debug)
    {
        info("Use --debug to print this simulator's complete configuration.");
    }
    else
    {
        cout << endl;
        cout << setw(30) << left << "[ setup ]" << endl;
        cout << setw(30) << left << "maxtries = " << setw(30) << left << conf.maxtries << endl;
        cout << setw(30) << left << "mindist = "  << setw(30) << left << conf.mindist << endl;
        cout << setw(30) << left << "dt = " << conf.dt << endl;
        cout << endl;
        cout << setw(30) << left << "[ runcontrol ]" << endl;
        cout << setw(30) << left << "nsteps = " << setw(30) << left << conf.nsteps << endl;
        cout << setw(30) << left << "eql_steps = " << setw(30) << left << conf.eql_steps << endl;
        cout << setw(30) << left << "nsample = " << setw(30) << left << conf.step_sample << endl;
        cout << setw(30) << left << "nblocks = " << setw(30) << left << conf.nblocks << endl;
        cout << endl;
        cout << setw(30) << left << "[ system ]" << endl;
        cout << setw(30) << left << "natoms = " << setw(30) << left << conf.natoms << endl;
        cout << setw(30) << left << "rho = " << setw(30) << left << conf.rho << endl;
        cout << setw(30) << left << "inittemp = " << setw(30) << left << conf.temp << endl;
        cout << setw(30) << left << "rcut = " << setw(30) << left << conf.rcut << endl;
        cout << setw(30) << left << "rlist = " << setw(30) << left << conf.rlist << endl;
        cout << setw(30) << left << "nlist = " << setw(30) << left << conf.nlist << endl;
        cout << endl;
        cout << setw(30) << left << "[ output ]" << endl;
        cout << setw(30) << left << "pdbfile = " << setw(30) << left << conf.pdbfile << endl;
        cout << setw(30) << left << "xtcfile = " << setw(30) << left << conf.xtcfile << endl;
        cout << setw(30) << left << "nxtc = " << setw(30) << left << conf.nxtc << endl;
        cout << setw(30) << left << "nlog = " << setw(30) << left << conf.nlog << endl;
        cout << endl;
        cout << setw(30) << left << "[ temperature ]" << endl;
        cout << setw(30) << left << "reft = " << setw(30) << left << conf.reft << endl;
        cout << setw(30) << left << "coupl = " << setw(30) << left << conf.tcouplstr << endl;
        cout << setw(30) << left << "coll_freq = " << setw(30) << left << conf.coll_freq << endl;
        cout << endl;
        cout << setw(30) << left << "[ rdf ]" << endl;
        cout << setw(30) << left << "sample = " << setw(30) << left << conf.dordfstr << endl;
        cout << setw(30) << left << "nbins = " << setw(30) << left << conf.rdf_nbins << endl;
        cout << setw(30) << left << "outfile = " << setw(30) << left << conf.rdf_outfile << endl;
        cout << setw(30) << left << "freq = " << setw(30) << left << conf.rdf_freq << endl;
        cout << endl;
        cout << setw(30) << left << "[ velocity ]" << endl;
        cout << setw(30) << left << "sample = " << setw(30) << left << conf.dovelstr << endl;
        cout << setw(30) << left << "min = " << setw(30) << left << conf.v_min << endl;
        cout << setw(30) << left << "max = " << setw(30) << left << conf.v_max << endl;
        cout << setw(30) << left << "nbins = " << setw(30) << left << conf.v_nbins << endl;
        cout << setw(30) << left << "outfile = " << setw(30) << left << conf.v_outfile << endl;
        cout << setw(30) << left << "freq = " << setw(30) << left << conf.v_freq << endl;
        cout << endl;
    }
    auto config = conf;
    sys.Initialize(config, config.natoms, config.nsteps, config.rho, config.rcut, config.rlist, config.temp, config.dt, config.mindist, config.maxtries, config.pdbfile, config.reft, config.coll_freq, config.xtcfile, config.rdf_nbins, config.rdf_outfile, config.v_nbins, config.v_max, config.v_min, config.v_outfile);
    return true;
}

void LJ::Update (int nd, int np, double pos[], double vel[], double f[], double acc[], double mass, double dt)
{}

void LJ::Compute (int nd, int np, double pos[], double vel[], double mass, double f[], double &pot, double &kin)
{}

void LJ::HostCPURun() 
{
    info("Using OpenMP to parallelize calculating forces on each atom from neighboring atoms and for velocity Verlet integration.");
    info("Press Ctrl-C to stop simulation.");
    auto sim_start = std::chrono::high_resolution_clock::now();
    sys.UpdateNeighborListHostCPU();
    sys.CalcForceHostCPU();
    sys.PrintHeader();
    sys.Print(0);
    for (int step = 1; step < conf.nsteps; step++)
    {
        // Main part of algorithm
        sys.IntegrateHostCPU(0, conf.tcoupl);
        sys.CalcForceHostCPU();
        sys.IntegrateHostCPU(1, conf.tcoupl);
        // Update the neighbor list this step?
        if (step % conf.nlist == 0)
        {
            sys.UpdateNeighborListHostCPU();
        }
        // Sample the RDF this step?
        if (( conf.dordf == true) && (step % conf.rdf_freq == 0) && (step > conf.eql_steps))
        {
            sys.SampleRdf();
        }

        // Sample the velocity distribution this step?
        if (( conf.dovel == true) && (step % conf.v_freq == 0) && (step > conf.eql_steps))
        {
            sys.SampleVel();
        }

        // Do other sampling this step?
        if ( (step % conf.step_sample) == 0 && (step > conf.eql_steps) )
        {
            sys.Sample();
        }

        // Print to the log this step?
        if (step % conf.nlog == 0)
        {
            sys.Print(step);
        }

        // Write to the xtc file this step?
        if (step % conf.nxtc == 0)
        {
            sys.WriteXTC(step);
        }

    }
    auto sim_time = duration_cast<microseconds>(high_resolution_clock::now() - sim_start).count() / 1000.0;
    sys.CloseXTC();
    info("Simulation total time is {:03.0f}ms.", sim_time);
    if (conf.dordf == true)
    {
        sys.NormalizeRdf();
        sys.OutputRdf();
    }

    if (conf.dovel == true)
    {
        sys.NormalizeVel();
        sys.OutputVel();
    }
    sys.ErrorAnalysis(conf.nblocks);
    sys.NormalizeAverages();
    sys.PrintAverages();
} 

#ifdef USE_ONEAPI
void LJ::CPURun() 
{
    auto device_name = q.get_device().get_info<sycl::info::device::name>();
    info("Using CPU device: {}.", device_name);
    info("Press Ctrl-C to abort simulation.");
    auto sim_start = std::chrono::high_resolution_clock::now();
    sys.UpdateNeighborListHostCPU();
    sys.CalcForceHostCPU();
    sys.PrintHeader();
    sys.Print(0);
    for (int step = 1; step < conf.nsteps; step++)
    {
        // Main part of algorithm
        sys.IntegrateHostCPU(0, conf.tcoupl);
        sys.CalcForceHostCPU();
        sys.IntegrateHostCPU(1, conf.tcoupl);
        // Update the neighbor list this step?
        if (step % conf.nlist == 0)
        {
            sys.UpdateNeighborListHostCPU();
        }
        // Sample the RDF this step?
        if (( conf.dordf == true) && (step % conf.rdf_freq == 0) && (step > conf.eql_steps))
        {
            sys.SampleRdf();
        }

        // Sample the velocity distribution this step?
        if (( conf.dovel == true) && (step % conf.v_freq == 0) && (step > conf.eql_steps))
        {
            sys.SampleVel();
        }

        // Do other sampling this step?
        if ( (step % conf.step_sample) == 0 && (step > conf.eql_steps) )
        {
            sys.Sample();
        }

        // Print to the log this step?
        if (step % conf.nlog == 0)
        {
            sys.Print(step);
        }

        // Write to the xtc file this step?
        if (step % conf.nxtc == 0)
        {
            sys.WriteXTC(step);
        }

    }
    auto sim_time = duration_cast<microseconds>(high_resolution_clock::now() - sim_start).count() / 1000.0;
    sys.CloseXTC();
    info("Simulation total time is {:03.0f}ms.", sim_time);
    if (conf.dordf == true)
    {
        sys.NormalizeRdf();
        sys.OutputRdf();
    }

    if (conf.dovel == true)
    {
        sys.NormalizeVel();
        sys.OutputVel();
    }
    sys.ErrorAnalysis(conf.nblocks);
    sys.NormalizeAverages();
    sys.PrintAverages();
} 
#endif
