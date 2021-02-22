# oneMD
oneMD is a molecular dynamics simulator designed to take advantage of SIMD and GPU and FPGA hardware data parallelism and acceleration. Unlike other molecular dynamics simulators like NAMD or omegagene which rely on hardware acceleration via OpenMP or CUDA, oneMD is a new simulator written in DPC++ designed especially to utilize diverse kinds of hardware acceleration across multiple architectures .

The molecular data and trajectories are saved in standard [PDB](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction) and .xtc formats that can be viewed in programs like [VMD](https://www.ks.uiuc.edu/Research/vmd/).

## Requirements
oneMD is cross-platform and cross-arch can be built with any well-known C++ compiler including MSVC. If you are using it on your local machine you will need a oneAPI environment setup with the DPC++ compiler to build with GPU and FPGA acceleration. If you don't use a oneAPI environment then only OpenMP acceleration support will be built. 
## Installation
* Clone the repo including all sub-modules:
 ````
 git clone --recurse-submodules https://github.com/allisterb/oneMD.git
 ````
* Run `build.sh` on DevCloud or Linux, or `build` on Windows. This will place a `omd` executable in the project folder. Run
````
omd --help
````
from the project folder and you should see the oneMD help screen.
````
u61437@s001-n058:~/oneMD$ ./omdd --help                                   
                 __  __ ___                                               
 ___  _ _   ___ |  \/  |   \                                              
/ _ \| ' \ / -_)| |\/| | |) |                                             
\___/|_||_|\___||_|  |_|___/                                              
                                                                          
                                                                          
USAGE:                                                                    
                                                                          
   ./omdd  [-e <string>] [--dt <integer>] [-t <integer>] [-n <integer>]   
           [--nd <integer>] [-d] [--] [--version] [-h] <string>           
                                                                          
                                                                          
Where:                                                                    
                                                                          
   -e <string>,  --device <string>                                        
     Name of hardware device, accelerator or library to run simulation on.
                                                                          
...                               
````
Note that there seems to be a bug when oneMD is built in Release mode on Windows. Run `build-debug` on Windows if the simulation hangs on Windows and use the debug executable `omdd`.