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
./omd --help
                 __  __ ___
 ___  _ _   ___ |  \/  |   \
/ _ \| ' \ / -_)| |\/| | |) |
\___/|_||_|\___||_|  |_|___/


USAGE:

   ./omd  [-l <integer>] [--dt <double>] [--ts <integer>] [--np <integer>]
          [--nd <integer>] [--device <string>] [-1] [-0] [-d] [--]
          [--version] [-h] <string>


Where:

   -l <integer>,  --debug-level <integer>
     Debug logging level. Default is 1.

   --dt <double>
     Timestep delta in seconds for simulation. Default is 0.005.

   --ts <integer>
     Number of time steps for simulation. Default is 10000.

   --np <integer>
     Number of particles for simulation. Default is 100

   --nd <integer>
     Number of dimensions for simulation from 1-3. Default is 3.

   --device <string>
     Name of hardware device to run simulation on. Can be host_cpu, cpu,
     gpu, or fpga. Default is host_cpu. Alternatively use the bool flag
     selectors e.g. -0 or -1.

   -1,  --cpu
     Select the SYCL CPU device.

   -0,  --host-cpu
     Select the host CPU device.

   -d,  --debug
     Enable debug logging.                               
````