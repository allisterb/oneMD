#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

#include <string>
#include "Enum.h"

#include "spdlog/spdlog.h"

using namespace std;

BETTER_ENUM(Device, int, CPU, FPGA);

class Simulator {
  private:
  protected:
  public:
    Simulator(string name, const int _nd, const int _np, const int _ts, const float _ts_delta, const Device _device);
    virtual bool Initialize() = 0;
    virtual void Compute (int nd, int np, double pos[], double vel[], double mass, double f[], double &pot, double &kin) = 0;
    virtual void Update (int nd, int np, double pos[], double vel[], double f[], double acc[], double mass, double dt) = 0;
    virtual double Distance (int nd, double r1[], double r2[], double dr[]) = 0;
    virtual void Run() = 0;
    const string name;
    const int np;
    const int nd;
    const int ts;
    const float ts_delta;
    const Device device;
};

#endif // __SIMULATOR_H__