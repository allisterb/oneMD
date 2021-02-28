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

#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#ifdef USE_ONEAPI
// oneDPL headers should be included before standard headers
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/iterator>
#endif

#include <string>
#include <chrono>
#include "Device.hh"
#include "simulator_config.h"
#include "spdlog/spdlog.h"

#ifdef USE_ONEAPI
#include "dpc_common.hpp"
#include <CL/sycl.hpp>
using namespace oneapi;
#endif

using namespace std;
using namespace std::chrono; 
using namespace spdlog;

class NotImplementedException : public std::logic_error
{
  public: 
    NotImplementedException();
};

class Simulator {
  private:
  protected:
  public:
    Simulator(const string _name, const int _nd, const int _np, const int _ts, const double _ts_delta, const Device _device);
    Simulator(const string _name, simulator_config config, const Device _device);
    virtual bool Initialize() = 0;
    virtual void Compute (int nd, int np, double pos[], double vel[], double mass, double f[], double &pot, double &kin) = 0;
    virtual void Update (int nd, int np, double pos[], double vel[], double f[], double acc[], double mass, double dt) = 0;
    virtual void HostCPURun() = 0;
#ifdef USE_ONEAPI
    virtual void CPURun();
    virtual void GPURun();
    virtual void FPGARun();
#endif
    const string name;
    const int np;
    const int nd;
    const int ts;
    const double ts_delta;
    const Device device;
    virtual ~Simulator();
    static simulator_config default_config();
    static int config_ini_handler(void* config, const char* section, const char* name, const char* value);
};

#endif // __SIMULATOR_H__