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

#pragma once

#include "common.hh"
#ifdef USE_ONEAPI
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

#include "spdlog/spdlog.h"

#include "Vec3.hh"
#include "CubicBox.hh"
#include "NeighborList.hh"
#include "PdbFile.hh"
#include "Rdf.hh"
#include "Thermostat.hh"
#include "ThermodynamicVariable.hh"
#include "Velocity.hh"
#include "xdrfile_xtc.h"
#include "simulator_config.h"

using namespace std;
using namespace spdlog;

constexpr double kB = 1.3806485279; // Boltzmann's Constant (J / K)
constexpr double oneSixth = 1.0/6.0;

class System2 
{
    protected:
        const simulator_config config;// system configuration
        const CubicBox box;
        const double rho;             // density (constant)
        const double rhokB;           // density * boltzmann's constant
        const double ecut;            // potential energy at cutoff
        const double etail;           // energy tail correction
        const double ptail;           // pressure tail correction
        const double rcut2;           // rcut*rcut
        const double rlist;           // n
        const int natoms;             // number of atoms in system
        const double dt;              // time step
        const double halfdt;          // 0.5 * dt
        const double halfdt2;         // 0.5 * dt*dt
        const double halfecut;        // 0.5 * ecut
        const double inatomsm1;       // 1.0/natoms - 1.0
        const double i2natoms;        // 1.0(2.0*natoms)
        const double i3natoms;        // 1.0(3.0*natoms)
        const int nsteps;             // number of steps for simulation to perform
        const int nsample;            // counter of number of samples
        double ke;                    // instantaneous kinetic energy
        double pe;                    // instantaneous potential energy
        double entot;                 // instantaneous total energy (ke + pe)
        double press;                 // instantaneous pressure
        double temp;                  // instantaneous temperature
        double vol;                   // volume (constant)
        std::array<int, 2> nlist;     // neighbor list
        vector <Vec3> forces;         // forces
        vector <Vec3> velocities;     // velocities
        vector <Vec3> positions;      // positions
        Rdf rdf;
        ThermodynamicVariable KineticEnergy;
        ThermodynamicVariable PotentialEnergy;
        ThermodynamicVariable Pressure;
        ThermodynamicVariable Temperature;
        ThermodynamicVariable TotalEnergy;
        Thermostat tstat;
        Velocity vel;
        XDRFILE *xd;
        std::chrono::time_point<std::chrono::high_resolution_clock> prev_time_point;
        #ifdef USE_ONEAPI
        sycl::queue q;
        #endif
    public:
        System(simulator_config conf);
        void Initialize();
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
        void SampleRdf()s;
        void SampleVel();
        void WriteXTC(int step);
        int GetTime();
        void ResetTimer();
        
    class friend Simulator;
};

#ifdef USE_ONEAPI
SYCL_EXTERN void PrintDebug(sycl:stream s, const* char m);
#endif