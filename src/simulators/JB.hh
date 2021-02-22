#ifndef __JB_H__
#define __JB_H__

#include "../Simulator.hh"

using namespace std;
using namespace spdlog;

class JB : public Simulator {
  private:
    double *acc;
    double ctime;
    double dt;
    double e0;
    double *force;
    double kinetic;
    double mass = 1.0;
    double *pos;
    double potential;
    int step;
    int step_num;
    int step_print;
    int step_print_index;
    int step_print_num;
    double *vel;
    double cpu_time ();
    void r8mat_uniform_ab (int, int, double, double, int&, double[]);
    void timestamp ();

  public:
    JB(const int, const int, const int, const float, const Device);
    bool Initialize();
    void Compute (int, int, double[], double[], double, double[], double&, double&);
    void Update (int, int, double[], double[], double[], double[], double, double);
    double Distance (int, double[], double[], double[]);
    void CPURun();
};

#endif // __JB_H__

