#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include <string>
#include "Enum.h"

#include "spdlog/spdlog.h"

using namespace std;

BETTER_ENUM(Device, int, CPU, FPGA);

typedef struct
{
    const double mindist;
    const double maxtries;
    const double dt;
    const int nsteps;
    const int eql_steps;
    const int step_sample;
    const int nblocks;
    const int natoms;
    const double rho;
    double temp;
    const double rcut;
    const double rlist;
    const int nlist;
    const string pdbfile;
    const string xtcfile;
    const int nxtc;
    const int nlog;
    const std::string tcouplstr;
    bool tcoupl;
    double coll_freq;
    double reft;
    string dordfstr;
    bool dordf;
    int rdf_nbins;
    string rdf_outfile;
    int rdf_freq;
    string dovelstr;
    bool dovel; 
    double v_max;
    double v_min;
    int v_nbins;
    string v_outfile;
    int v_freq;
} configuration;

class Simulator {
  private:
  protected:
  public:
    Simulator(string name, const int _nd, const int _np, const int _ts, const float _ts_delta, const Device _device);
    Simulator(const string _name, configuration config, const Device _device);
    virtual bool Initialize() = 0;
    virtual void Compute (int nd, int np, double pos[], double vel[], double mass, double f[], double &pot, double &kin) = 0;
    virtual void Update (int nd, int np, double pos[], double vel[], double f[], double acc[], double mass, double dt) = 0;
    virtual void Run() = 0;
    const string name;
    const int np;
    const int nd;
    const int ts;
    const float ts_delta;
    const Device device;
    static configuration default_config();
    static int config_ini_handler(void* config, const char* section, const char* name, const char* value);
};

#endif // __SIMULATOR_H__