#pragma once

#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include <string>
#include "enum.h"

using namespace std;

BETTER_ENUM(Device, int, CPU, FPGA);

class Simulator {
  private:
  protected:
  public:
    const string name;
    const int np;
    const int nd;
    const int ts;
    const float ts_delta;
    const Device device;
    Simulator(const string, const int np, const int nd, const int ts, const float ts_delta, const Device _device);
    virtual void Compute ( int np, int nd, double pos[], double vel[], 
        double mass, double f[], double &pot, double &kin ) = 0;
    virtual void Update ( int np, int nd, double pos[], double vel[], double f[], 
      double acc[], double mass, double dt ) = 0;
    virtual double Distance ( int nd, double r1[], double r2[], double dr[] ) = 0;
    //double cpu_time ( );
    //void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
    //void timestamp ( );
};

Simulator::Simulator(const string _name, const int _np, const int _nd, const int _ts, const float _ts_delta, const Device _device) :
name(_name),
np(_np),
nd(_nd),
ts(_ts),
ts_delta(_ts_delta),
device(_device)
{
  
}

#endif // __SIMULATOR_H__