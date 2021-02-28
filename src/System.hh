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

#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#ifdef USE_ONEAPI
#include "dpc_common.hpp"
#include <CL/sycl.hpp>
using namespace oneapi;

// oneDPL headers should be included before standard headers
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/iterator>
#endif

#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <vector>

#include "Vec3.hh"
#include "NeighborList.hh"
#include "PdbFile.hh"
#include "Rdf.hh"
#include "Thermostat.hh"
#include "ThermodynamicVariable.hh"
#include "CubicBox.hh"
#include "Velocity.hh"
#include "xdrfile_xtc.h"
#include "system_configuration.h"
#include "spdlog/spdlog.h"

using namespace std;
using namespace spdlog;

constexpr double kB = 1.3806485279; // Boltzmann's Constant (J / K)
constexpr double oneSixth = 1.0/6.0;

class System 
{
    private:
        system_configuration config;  // system configuration
        const double dt;              // time step
        const double ecut;            // potential energy at cutoff
        const double entot;           // instantaneous total energy (ke + pe)
        const double etail;           // energy tail correction
        const double halfdt;          // 0.5 * dt
        const double halfdt2;         // 0.5 * dt*dt
        const double halfecut;        // 0.5 * ecut
        const double inatomsm1;       // 1.0/natoms - 1.0
        const double i2natoms;        // 1.0(2.0*natoms)
        const double i3natoms;        // 1.0(3.0*natoms)
        const double ke;              // instantaneous kinetic energy
        double pe;              // instantaneous potential energy
        double press;           // instantaneous pressure
        const double ptail;           // pressure tail correction
        const double rcut2;           // rcut*rcut
        const double rho;             // density (constant)
        const double rhokB;           // density * boltzmann's constant
        double temp;            // instantaneous temperature
        double vol;             // volume (constant)
        const int natoms;             // number of atoms in system
        int nsample;            // counter of number of samples
        const int nsteps;             // number of steps for simulation to perform
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
        System(system_configuration conf, int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile);
        void Initialize(system_configuration conf, int natoms, int nsteps, double rho, double rcut, double rlist, double temp, double dt, double mindist, double maxtries, string pdbfile, double reft, double coll_freq, string xtcfile, int rdf_nbins, string rdf_outfile, int v_nbins, double v_max, double v_min, string v_outfile);
        void UpdateNeighborListHostCPU();
        void CalcForceHostCPU();
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
};
#endif //__SYSTEM_H__