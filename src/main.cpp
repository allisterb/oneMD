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
#include "Util.hh"
#include "Simulator.hh"
#include "enum.h"

using namespace std;
using namespace spdlog;
using namespace TCLAP;

int
figlet_demo() {
  Figlet::standard.print("Fractions");
  for ( int i = 2; i <= 4; ++i ) {
    ostringstream ss;
    ss << "5/" << i << " = " << 5.0/i;
    Figlet::small.print(ss.str().c_str());
  }
  cout << "ALL DONE!\n";
  return 0;
}

int main ( int argc, char *argv[] );

int main ( int argc, char *argv[] )
{ 
  Figlet::small.print("oneMD");
	try {  
    CmdLine cmd("oneMD data-parallel molecular dynamics simulator.", ' ', "0.1", true);
    UnlabeledValueArg<string> simArg("simulator","Name of simulator to run.", true, "cpu", "string", cmd);
    UnlabeledValueArg<string> devArg("device","Name of hardware device, accelerator or library to run simulation on.", true, "cpu", "string", cmd);
    UnlabeledValueArg<int> ndArg("dimensions","Number of dimensions for simulation from 1 - 3.",true, 1,"integer", cmd);
    UnlabeledValueArg<int> npArg("particles","Number of particles for simulation.",true, 1,"integer", cmd);
    UnlabeledValueArg<int> tsArg("timesteps","Number of time steps for simulation.",true, 1,"integer", cmd);
    UnlabeledValueArg<float> tsdeltaArg("ts_delta","Timestep delta in seconds.",true, 0.2f,"integer", cmd);
    SwitchArg debugArg("d","debug","Enable debug-level logging.", cmd, false);
    
    cmd.parse( argc, argv );
    auto debugLog = debugArg.getValue();
    if (debugLog) {
      set_level(level::debug);
      info("Debug-level logging enabled.");
    }
    auto simulator = simArg.getValue();
    auto device = Device::_from_string(Util::upper(devArg.getValue()).c_str());
    auto nd = ndArg.getValue();
    auto np = npArg.getValue();
    auto ts = tsArg.getValue();
    auto tsDelta = tsdeltaArg.getValue();
	}
  catch (ArgException &e) { 
    error("Error parsing option {0}: {1}.", e.argId(), e.error());
    return 1;
  }
  catch (std::exception &e) { 
    error("Runtime error parsing options: {0}", e.what());
    return 1; 
  }
  return 0;
}