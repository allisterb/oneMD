#include "Simulator.hh"
#include "xdrfile.h"

NotImplementedException::NotImplementedException() : std::logic_error("Not yet implemented") { };

Simulator::Simulator(const string _name, const int _nd, const int _np, const int _ts, const float _ts_delta, const Device _device) :
  name(_name),
  nd(_nd),
  np(_np),
  ts(_ts),
  ts_delta(_ts_delta),
  device(_device)
{}

Simulator::Simulator(const string _name, configuration config, const Device _device) :
  name(_name),
  nd(3),
  np(config.natoms),
  ts(config.nsteps),
  ts_delta(config.dt),
  device(_device)
{
  info("config.natoms = {}. config.nsteps = {}. config.dt = {}. config.device={}", config.natoms, config.nsteps, config.dt, config.device._to_string());
}

Simulator::~Simulator() {}

configuration Simulator::default_config()
{ 
  return {
    .debug = false,
    .device = Device::HOST_CPU,
    .mindist = 1.0,
    .maxtries = 10e6,
    .dt = 0.005,
    .nsteps = 5000000,
    .eql_steps = 10000,
    .step_sample = 1000,
    .nblocks = 5,
    .natoms = 108,
    .rho = 0.5,
    .temp = 1.0,
    .rcut = 2.5,
    .rlist = 3.5,
    .nlist = 10,
    .pdbfile ="init.pdb",
    .xtcfile = "traj.xtc",
    .nxtc = 1000,
    .nlog = 1000,
    .tcouplstr = "no",
    .tcoupl = false,
    .coll_freq = 0.001,
    .reft = 1.0,
    .dordfstr = "no",
    .dordf = false,
    .rdf_nbins = 100,
    .rdf_outfile = "rdf.dat",
    .rdf_freq = 1000,
    .dovelstr = "no",
    .dovel = false,
    .v_max = 10.0,
    .v_min = -10.0,
    .v_nbins = 100,
    .v_outfile = "vel_dist.dat",
    .v_freq = 1000
  };
}

int Simulator::config_ini_handler(void* c, const char* section, const char* name,
                   const char* value)
{
  configuration* pconfig = (configuration*) c;
  #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0  
  return 1;
}
#ifdef USE_ONEAPI
void Simulator::CPURun()
{
  throw new NotImplementedException();
}

void Simulator::GPURun()
{
  throw new NotImplementedException();
}

void Simulator::FPGARun()
{
  throw new NotImplementedException();
}
#endif