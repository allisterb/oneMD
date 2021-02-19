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
#include "enum.h"

#include "Util.hh"
#include "simulators/jburkardt.hh"

using namespace std;
using namespace spdlog;
using namespace TCLAP;

int main ( int argc, char *argv[] )
{ 
  Figlet::small.print("oneMD");
	try 
  {  
    CmdLine cmd("oneMD data-parallel molecular dynamics simulator.", ' ', "0.1", true);
    UnlabeledValueArg<string> simArg("simulator","Name of simulator to run.", true, "cpu", "string", cmd);
    UnlabeledValueArg<int> ndArg("dimensions","Number of dimensions for simulation from 1 - 3.",true, 1,"integer", cmd);
    UnlabeledValueArg<int> npArg("particles","Number of particles for simulation.",true, 1,"integer", cmd);
    UnlabeledValueArg<int> tsArg("timesteps","Number of time steps for simulation.",true, 1,"integer", cmd);
    UnlabeledValueArg<float> tsdeltaArg("ts_delta","Timestep delta in seconds.",true, 0.2f,"integer", cmd);
    UnlabeledValueArg<string> devArg("device","Name of hardware device, accelerator or library to run simulation on.", false, "CPU", "string", cmd);
    SwitchArg debugArg("d","debug","Enable debug-level logging.", cmd, false);
    cmd.parse( argc, argv );
    auto debugLog = debugArg.getValue();
    if (debugLog) {
      set_level(level::debug);
      info("Debug-level logging enabled.");
    }
    auto simulator = simArg.getValue();
    auto nd = ndArg.getValue();
    auto np = npArg.getValue();
    auto ts = tsArg.getValue();
    auto tsDelta = tsdeltaArg.getValue();
    auto device = Device::_from_string(Util::upper(devArg.getValue()).c_str());
    auto jb = JB(nd, np, ts, tsDelta, device);
    jb.Initialize();
    jb.Run();
    return 0;
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
}