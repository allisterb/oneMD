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

#ifndef __LJ_H__
#define __LJ_H__

#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <vector>

#include "../Vec3.hh"
#include "../NeighborList.hh"
#include "../PdbFile.hh"
#include "../Rdf.hh"
#include "../Thermostat.hh"
#include "../ThermodynamicVariable.hh"
#include "../CubicBox.hh"
#include "../Velocity.hh"
#include "../Simulator.hh"
#include "xdrfile_xtc.h"

#include "spdlog/spdlog.h"

using namespace std;
using namespace spdlog;

const double kB = 1.3806485279; // Boltzmann's Constant (J / K)
const double oneSixth = 1.0/6.0;

class System 
{
    private:
        simulator_config conf;
        double dt;              // time step
        double ecut;            // potential energy at cutoff
        double entot;           // instantaneous total energy (ke + pe)
        double etail;           // energy tail correction
        double halfdt;          // 0.5 * dt
        double halfdt2;         // 0.5 * dt*dt
        double halfecut;        // 0.5 * ecut
        double inatomsm1;       // 1.0/natoms - 1.0
        double i2natoms;        // 1.0(2.0*natoms)
        double i3natoms;        // 1.0(3.0*natoms)
        double ke;              // instantaneous kinetic energy
        double pe;              // instantaneous potential energy
        double press;           // instantaneous pressure
        double ptail;           // pressure tail correction
        double rcut2;           // rcut*rcut
        double rho;             // density (constant)
        double rhokB;           // density * boltzmann's constant
        double temp;            // instantaneous temperature
        double vol;             // volume (constant)
        int natoms;             // number of atoms in system
        int nsample;            // counter of number of samples
        int nsteps;             // number of steps for simulation to perform
        NeighborList nlist;
        Rdf rdf;
        ThermodynamicVariable KineticEnergy;
        ThermodynamicVariable PotentialEnergy;
        ThermodynamicVariable Pressure;
        ThermodynamicVariable Temperature;
        ThermodynamicVariable TotalEnergy;
        Thermostat tstat;
        vector <Vec3> f; // forces
        vector <Vec3> v; // velocities
        vector <Vec3> x; // positions
        Velocity vel;
        XDRFILE *xd;
        std::chrono::time_point<std::chrono::high_resolution_clock> prev_time_point;
        #ifdef USE_ONEAPI
        sycl::queue q;
        #endif
    public:
        System(simulator_config c, int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile);
        void Initialize(simulator_config c, int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile);
        void CalcForceHostCPU();
        void UpdateNeighborListHostCPU();
        void IntegrateHostCPU(int a, bool tcoupl);
        void CloseXTC();
        void ErrorAnalysis(int nblocks);
        void NormalizeAverages();
        void NormalizeRdf();
        void NormalizeVel();
        void OutputVel();
        void OutputRdf();
        void Print(int step);
        void PrintAverages();
        void PrintHeader();
        void Sample();
        void SampleRdf();
        void SampleVel();
        void WriteXTC(int step);
        int GetTime();
        void ResetTimer();
        CubicBox box;
#ifdef USE_ONEAPI
    void CalcForceCPU();
    void UpdateNeighborListCPU();
    void IntegrateCPU(int a, bool tcoupl);
    void CalcForceGPU();
    void UpdateNeighborListGPU();
    void IntegrateGPU(int a, bool tcoupl);
    void CalcForceFPGA();
    void UpdateNeighborListFPGA();
    void IntegrateFPGA(int a, bool tcoupl);
#endif
};

class LJ : public Simulator {
  private:
#ifdef USE_ONEAPI
    sycl::queue q;
#endif
  public:
    LJ(simulator_config c, const Device);
    bool Initialize();
    void Compute (int, int, double[], double[], double, double[], double&, double&);
    void Update (int, int, double[], double[], double[], double[], double, double);
    double Distance (int, double[], double[], double[]);
    void HostCPURun();
#ifdef USE_ONEAPI
    void CPURun();

#endif
    System sys;
    simulator_config conf;
};
#endif

#ifdef USE_ONEAPI
template<typename T>
SYCL_EXTERNAL void PrintDebug(sycl::stream s, const char* name, T m);
#endif