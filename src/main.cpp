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

#include "Util.hh"
#include "simulators/JB.hh"
#include "simulators/LJ.hh"

using namespace std;
using namespace spdlog;
using namespace TCLAP;

int main (int argc, char *argv[] )
{ 
  Figlet::small.print("oneMD");
	unique_ptr<Simulator> sim(nullptr);
  auto config = Simulator::default_config();
  try 
  {  
    CmdLine cmd("oneMD data-parallel molecular dynamics simulator.", ' ', "0.1", true);
    UnlabeledValueArg<string> simArg("simulator","Name of simulator to run.", true, "", "string", cmd);
    ValueArg<int> ndArg("", "nd","Number of dimensions for simulation from 1 - 3.",false, 2,"integer");
    ValueArg<int> npArg("n", "np","Number of particles for simulation.",false, 100,"integer");
    ValueArg<int> tsArg("t", "ts","Number of time steps for simulation.",false, 1000,"integer");
    ValueArg<double> tsdeltaArg("", "dt","Timestep delta in seconds.",false, 0.005,"integer");
    ValueArg<string> devArg("e", "device","Name of hardware device, accelerator or library to run simulation on.", false, "HOST_CPU", "string");
    ValueArg<string> configArg("c","config","Name of configuration file for simulation,",false,"","string");
    SwitchArg debugArg("d","debug","Enable debug-level logging.", cmd, false);
    SwitchArg hostCPUDeviceArg("0","host-cpu","Select the host CPU device.", cmd, false);
    SwitchArg cpuDeviceArg("1","cpu","Select the SYCL CPU device.", cmd, false);
    cmd.add(ndArg);
    cmd.add(npArg);
    cmd.add(tsArg);
    cmd.add(tsdeltaArg);
    cmd.add(devArg);
    cmd.parse(argc, argv);
    auto debugLog = debugArg.getValue();
    if (debugLog) {
      set_level(level::debug);
      config.debug = true;
      info("Debug-level logging enabled.");
    }
    auto sim_name = Util::upper(simArg.getValue());
    auto nd = ndArg.getValue();
    auto np = npArg.getValue();
    auto ts = tsArg.getValue();
    auto ts_delta = tsdeltaArg.getValue();
    auto device = Device::_from_string(Util::upper(devArg.getValue()).c_str());
    if (hostCPUDeviceArg.getValue())
    {
      device = Device::HOST_CPU;
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
      case Device::HOST_CPU:
        sim->Initialize();
        sim->HostCPURun();
        return 0;
#ifdef USE_ONEAPI
      default:
        error("The {} device is not implemented.", config.device._to_string());
        return 2;
#else
      default:
        error("The {} device is not enabled in this build of oneMD. Build oneMD inside a oneAPI environment to enable this device.", config.device._to_string());
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