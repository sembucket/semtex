///////////////////////////////////////////////////////////////////////////////
// odeint.C: a code to integrate a set of ODES to produce tabulated results.
//
// Based on dumb RK4 (Numerical Recipes), and straightforward:
// everything is hard coded.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>

#include <utility.h>
#include <lapack.h>

#include <Stack.h>
#include <Array.h>

#include "nr77.h"

using namespace std;


static void ODEsys (const double& x   ,
		    const double* y   ,
		          double* dydx)
// ---------------------------------------------------------------------------
// This is the set of N ODEs; here, the autonomous Lorenz system.
// ---------------------------------------------------------------------------
{
  static const double sigma = 10.0;
  static const double b     = 8.0/3.0;
  static const double r     = 28.0;

  dydx[0] =  sigma * (y[1] - y[0]);
  dydx[1] =  r * y[0] - y[1] - y[0]*y[2];
  dydx[2] = -b * y[2] + y[0]*y[1];
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Carry manage the integration and I/O.
// ---------------------------------------------------------------------------
{
  const int    Ns = 50000;
  const int    Nv = 3;
  const double dt = 0.003;
  vector<real> y    (Nv);
  vector<real> dydx (Nv);
  int          i;
  double       t;

  // -- set ICs.
  
  t    = 0.0;
  y[0] = 0.0;
  y[1] = 0.1;
  y[2] = 0.1;

  for (i = 0; i < Ns; i++) {
    (*ODEsys)    (t, y(), dydx());
    Recipes::rk4 (y(), dydx(), Nv, t, dt, y(), ODEsys);
    t += dt;
    cout << t << '\t' << y[0] << '\t' << y[1] << '\t' << y[2] << endl;
  }
  
  return EXIT_SUCCESS;
}


