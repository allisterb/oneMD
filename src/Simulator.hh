#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__

class Simulator {
  private:
  protected:
  public:
    Simulator( int np, int nd, double pos[], double vel[], double acc[] );
    virtual void Compute ( int np, int nd, double pos[], double vel[], 
        double mass, double f[], double &pot, double &kin ) = 0;
    virtual void Update ( int np, int nd, double pos[], double vel[], double f[], 
      double acc[], double mass, double dt ) = 0;
    virtual double Distance ( int nd, double r1[], double r2[], double dr[] ) = 0;
    double cpu_time ( );
    void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
    void timestamp ( );

};

#endif
