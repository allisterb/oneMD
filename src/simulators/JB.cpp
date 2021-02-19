# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
# include <sstream>

#include "JB.hh"

JB::JB(const int _nd, const int _np, const int _ts, const float _ts_delta, const Device _device) :
  Simulator("JB", _nd, _np, _ts, _ts_delta, _device),
  step_num(_ts),
  dt(_ts_delta)
{}

void JB::r8mat_uniform_ab (int m, int n, double a, double b, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A, B, the limits of the pseudorandom values.
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return;
}

void JB::timestamp ()

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

double JB::cpu_time ()

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}

bool JB::Initialize()
{
  int i;
  int j;
  int seed;
  acc = new double[nd*np];
  force = new double[nd*np];
  pos = new double[nd*np];
  vel = new double[nd*np];
  //
  //  Set the positions.
  //
  seed = 123456789;
  JB::r8mat_uniform_ab (nd, np, 0.0, 10.0, seed, pos);
  //
  //  Set the velocities.
  //
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      vel[i+j*nd] = 0.0;
    }
  }
  //
  //  Set the accelerations.
  //
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      acc[i+j*nd] = 0.0;
    }
  }
  JB::Compute(nd, np, pos, vel, mass, force, potential, kinetic);
  e0 = potential + kinetic;
  return true;
}

void JB::Run() 
{
  //
  //  This is the main time stepping loop:
  //    Compute forces and energies,
  //    Update positions, velocities, accelerations.
  //
  cout << "\n";
  cout << "  At each step, we report the potential and kinetic energies.\n";
  cout << "  The sum of these energies should be a constant.\n";
  cout << "  As an accuracy check, we also print the relative error\n";
  cout << "  in the total energy.\n";
  cout << "\n";
  cout << "      Step      Potential       Kinetic        (P+K-E0)/E0\n";
  cout << "                Energy P        Energy K       Relative Energy Error\n";
  cout << "\n";

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10; 
  ctime = this->cpu_time ();
  for (step = 0; step <= step_num; step++)
  {
    JB::Update(nd, np, pos, vel, force, acc, mass, dt);
    JB::Compute(nd, np, pos, vel, mass, force, potential, kinetic);
    if (step == step_print)
    {
      cout << "  " << setw(8) << step
          << "  " << setw(14) << potential
          << "  " << setw(14) << kinetic
          << "  " << setw(14) << ( potential + kinetic - e0 ) / e0 << "\n";
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }

  }
  //
  //  Report timing.
  //
  ctime = this->cpu_time () - ctime;
  cout << "\n";
  cout << "  Elapsed cpu time " << ctime << " seconds.\n";
  //
  //  Free memory.
  //
  delete [] acc;
  delete [] force;
  delete [] pos;
  delete [] vel;
  //
  //  Terminate.
  //
  cout << "\n";
  cout << "MD\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  this->timestamp ();
}

void JB::Update (int nd, int np, double pos[], double vel[], double f[], double acc[], double mass, double dt)

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE updates positions, velocities and accelerations.
//
//  Discussion:
//
//    The time integration is fully parallel.
//
//    A velocity Verlet algorithm is used for the updating.
//
//    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
//    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
//    a(t+dt) = f(t) / m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input/output, double POS[ND*NP], the positions.
//
//    Input/output, double VEL[ND*NP], the velocities.
//
//    Input, double F[ND*NP], the forces.
//
//    Input/output, double ACC[ND*NP], the accelerations.
//
//    Input, double MASS, the mass of each particle.
//
//    Input, double DT, the time step.
//
{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;

  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }

  return;
}

void JB::Compute (int nd, int np, double pos[], double vel[], double mass, double f[], double &pot, double &kin)

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE computes the forces and energies.
//
//  Discussion:
//
//    The computation of forces and energies is fully parallel.
//
//    The potential function V(X) is a harmonic well which smoothly
//    saturates to a maximum value at PI/2:
//
//      v(x) = ( sin ( min ( x, PI2 ) ) )**2
//
//    The derivative of the potential is:
//
//      dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
//            = sin ( 2.0 * min ( x, PI2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2008
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, double POS[ND*NP], the position of each particle.
//
//    Input, double VEL[ND*NP], the velocity of each particle.
//
//    Input, double MASS, the mass of each particle.
//
//    Output, double F[ND*NP], the forces.
//
//    Output, double &POT, the total potential energy.
//
//    Output, double &KIN, the total kinetic energy.
//
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pot = 0.0;
  kin = 0.0;

  for ( k = 0; k < np; k++ )
  {
    //
    //  Compute the potential energy and forces.
    //
    for ( i = 0; i < nd; i++ )
    {
      f[i+k*nd] = 0.0;
    }

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = this->Distance ( nd, pos+k*nd, pos+j*nd, rij );
        //
        //  Attribute half of the potential energy to particle J.
        //
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pot = pot + 0.5 * pow ( sin ( d2 ), 2 );

        for ( i = 0; i < nd; i++ )
        {
          f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
        }
      }
    }
//
//  Compute the kinetic energy.
//
    for ( i = 0; i < nd; i++ )
    {
      kin = kin + vel[i+k*nd] * vel[i+k*nd];
    }
  }

  kin = kin * 0.5 * mass;
  
  return;
}

double JB::Distance(int nd, double r1[], double r2[], double dr[])

//****************************************************************************80
//
//  Purpose:
//
//    DIST computes the displacement (and its norm) between two particles.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, double R1[ND], R2[ND], the positions.
//
//    Output, double DR[ND], the displacement vector.
//
//    Output, double D, the Euclidean norm of the displacement.
//
{
  double d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}