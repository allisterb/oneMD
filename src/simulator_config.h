#ifndef __SIMULATOR_CONFIG_H__
#define __SIMULATOR_CONFIG_H__

#include <string>
#include "Device.hh"

typedef struct
{
    bool debug;
    Device device;
    double mindist;
    double maxtries;
    double dt;
    int nsteps;
    int eql_steps;
    int step_sample;
    int nblocks;
    int natoms;
    double rho;
    double temp;
    double rcut;
    double rlist;
    int nlist;
    std::string pdbfile;
    std::string xtcfile;
    int nxtc;
    int nlog;
    bool tcoupl;
    double coll_freq;
    double reft;
    std::string dordfstr;
    bool dordf;
    int rdf_nbins;
    std::string rdf_outfile;
    int rdf_freq;
    bool dovel; 
    double v_max;
    double v_min;
    int v_nbins;
    std::string v_outfile;
    int v_freq;
} simulator_config;

#endif //__SIMULATOR_CONFIG_H__