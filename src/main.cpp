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
  Figlet::standard.print("oneMD");
	try {  
    CmdLine cmd("oneMD data-parallel molecular dynamics simulator.", ' ', "0.1", true);
    UnlabeledValueArg<string> simArg("simulator","Name of simulator to run.", true, "cpu", "string", cmd);
    UnlabeledValueArg<string> devArg("device","Name of hardware device, accelerator or library to run simulation on.", true, "CPU", "string", cmd);
    UnlabeledValueArg<int> ndArg("dimensions","Number of dimensions for simulation from 1 - 3.",true, 1,"integer", cmd);
    UnlabeledValueArg<int> npArg("particles","Number of particles for simulation.",true, 1,"integer", cmd);
    UnlabeledValueArg<int> tsArg("timesteps","Number of time steps for simulation.",true, 1,"integer", cmd);
    UnlabeledValueArg<float> tsdeltaArg("ts_delta","Timestep delta in seconds.",true, 0.2f,"integer", cmd);
    SwitchArg debugArg("d","debug","Enable debug-level logging.", cmd, false);
    
    cmd.parse( argc, argv );
    auto debugLog = debugArg.getValue();
    auto simulator = simArg.getValue();
    auto device = Device::CPU;
    auto maybe_device = Device::_from_string(toupper(devArg.getValue(),locale()).c_str());
    //if (maybe_device) {
    //  device = maybe_device.value();
    //}
    //else
    //{
      //error("Invalid device selected: %s",devArg.getValue());
      //return 1;
    //}
    //switch (std::toupper(debugArg.getValue()))
    //{
      //

    //}
    auto nd = ndArg.getValue();
    auto np = npArg.getValue();
    auto ts = tsArg.getValue();
    auto tsDelta = tsdeltaArg.getValue();

    //if (debugLog) {
    //  set_level(level::debug);
    //  info("Debug-level logging enabled.");
    //}

	}
  catch (ArgException &e) { 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
  }
  return 0;
}