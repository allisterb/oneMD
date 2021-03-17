# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <sstream>

#include "Figlet.hh"
#include "spdlog/spdlog.h"
#include "tclap/CmdLine.h"
#include "tclap/UnlabeledValueArg.h"
#include "Enum.h"

#include "Util.hpp"
#include "simulators/JB.hh"
#include "simulators/LJ.hh"

using namespace std;
using namespace spdlog;
using namespace TCLAP;

int main (int argc, char *argv[])
{ 
  Figlet::small.print("oneMD");
  unique_ptr<Simulator> sim(nullptr);
  auto config = Simulator::default_config();
  try 
  {  
    CmdLine cmd("oneMD is a data-parallel molecular dynamics simulator for Intel oneAPI.", ' ', "0.1", true);
    UnlabeledValueArg<string> simArg("simulator","Name of simulator to run. Defaults to LJ (Lennard-Jones) atoms in cubic box.", false, "lj", "string", cmd);
    ValueArg<string> devArg("", "device","Name of hardware device or parallel computing framework to run simulation on. Can be openmp, cpu, gpu, or fpga. Default is openmp. Alternatively use the bool flag selectors e.g. -0 or -1.", false, "OPENMP", "string");
    ValueArg<string> configArg("","config","Name of configuration file for simulation,",false,"","string");
    ValueArg<int> ndArg("", "nd","Number of dimensions for simulation from 1-3. Default is 3.",false, 3,"integer");
    ValueArg<int> npArg("", "np","Number of particles for simulation. Default is 100",false, 100,"integer");
    ValueArg<int> tsArg("", "ts","Number of time steps for simulation. Default is 10000.",false, 10000,"integer");
    ValueArg<double> tsdeltaArg("", "dt","Timestep delta in seconds for simulation. Default is 0.005.",false, 0.005, "double");
    SwitchArg debugArg("d","debug","Enable debug logging.", cmd, false);
    ValueArg<int> debugLogLevelArg("l", "debug-level","Debug logging level. Default is 1.",false, 1,"integer");
    SwitchArg openmpDeviceArg("0","openmp","Select the OpenMP parallel processing framework.", cmd, false);
    SwitchArg cpuDeviceArg("1","cpu","Select the SYCL CPU device.", cmd, false);
    cmd.add(devArg);
    cmd.add(ndArg);
    cmd.add(npArg);
    cmd.add(tsArg);
    cmd.add(tsdeltaArg);
    cmd.add(debugLogLevelArg);
    cmd.parse(argc, argv);
    Simulator::debug_log = debugArg.getValue();
    Simulator::debug_log_level = debugLogLevelArg.getValue();
    if (Simulator::debug_log) {
      set_level(level::debug);
      config.debug = true;
      info("Debug-level logging enabled. Debug log level is {}.", Simulator::debug_log_level);
    }
    auto sim_name = Util::upper(simArg.getValue());
    auto nd = ndArg.getValue();
    auto np = npArg.getValue();
    auto ts = tsArg.getValue();
    auto ts_delta = tsdeltaArg.getValue();
    auto device = Device::_from_string(Util::upper(devArg.getValue()).c_str());
    if (openmpDeviceArg.getValue())
    {
      device = Device::OPENMP;
    }
    else if (cpuDeviceArg.getValue())
    {
      device = Device::CPU;
    }
    auto config_name = configArg.getValue();
    if (config_name != "")
    { 
      info("Using configuration file {} for options not specified in CLI.", config_name);
    }
    else
    {
      info("No configuration file specified. Using default configuration values for options not specified in CLI.");
    }
    config.natoms = np;
    config.nsteps = ts;
    config.dt = ts_delta;
    config.device = device;
    if (sim_name == "JB")
    {
      sim = make_unique<JB> (nd, np, ts, ts_delta, device);
    }
    else if (sim_name == "LJ")
    {
      sim =  make_unique<LJ>(config, device);
    }
	}
  catch (ArgException &e) 
  { 
    error("Error parsing option {0}: {1}.", e.argId(), e.error());
    return 1;
  }
  catch (std::exception &e) 
  { 
    error("Runtime error parsing options: {0}", e.what());
    return 1; 
  }

  try
  {
    switch (config.device)
    {
      case Device::OPENMP:
        sim->Initialize();
        sim->HostCPURun();
        return 0;
#ifdef USE_ONEAPI
      case Device::CPU:
        sim->Initialize();
        sim->CPURun();
        return 0;
      default:
        error("The {} device is not yet implemented.", config.device._to_string());
        return 2;
#else
      default:
        error("The {} device is not enabled in this build of oneMD. Only host_cpu acceleration via OpenMP is available. \n"
        "Build oneMD using the DPC++ compiler inside a oneAPI environment to enable this device.", config.device._to_string());
        return 2;
#endif
    }
  }
  catch (std::exception &e) 
  { 
    error("Runtime error in simulator {}: {}", sim->name, e.what());
    return 1; 
  }
}