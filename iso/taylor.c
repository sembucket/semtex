/*****************************************************************************
 * taylor.c: routines for Taylor & Taylor--Green vortices.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


void Taylor2D (CVF  IC, const int*  Dim, const int  code)
/* ------------------------------------------------------------------------- *
 * Generate initial conditions corresponding to the 2D Taylor Flow
 *
 * u = -cos(x) sin (y) exp(-2 \nu t)
 * v =  sin(x) cos (y) exp(-2 \nu t)
 * w =  0
 * p = -0.25 [cos(2x) + cos(2y)] exp(-4 \nu t).
 * 
 * Input code = 0, 1, 2 generates a cyclic permutation of the velocity.
 *
 * Generate initial conditions for t = 0, components in PHYSICAL space.
 * ------------------------------------------------------------------------- */
{
  const    int    N = Dim[1];
  register real  *u = &IC[1][0][0][0].Re;
  register real  *v = &IC[2][0][0][0].Re;
  register real  *w = &IC[3][0][0][0].Re;
  register int    i, j, k;

  /* -- Fill up the cube. */

  switch (code) {
  case 0:
    for (i = 0; i < N; i++) {
      const double x = 2.0 * M_PI * i / (double) N;
      for (j = 0; j < N; j++) {
	const double y = 2.0 * M_PI * j / (double) N;
	for (k = 0; k < N; k++) {
	  const double z = 2.0 * M_PI * k / (double) N;
	  
	  u[k + N * (j + i * N)] = -cos(x) * sin(y);
	  v[k + N * (j + i * N)] =  sin(x) * cos(y);
	  w[k + N * (j + i * N)] =  0.0;
	}
      }
    }
    break;

  case 1:
    for (i = 0; i < N; i++) {
      const double x = 2.0 * M_PI * i / (double) N;
      for (j = 0; j < N; j++) {
	const double y = 2.0 * M_PI * j / (double) N;
	for (k = 0; k < N; k++) {
	  const double z = 2.0 * M_PI * k / (double) N;
	  
	  u[k + N * (j + i * N)] =  0.0;
	  v[k + N * (j + i * N)] = -cos(y) * sin(z);
	  w[k + N * (j + i * N)] =  sin(y) * cos(z);;
	}
      }
    }
    break;

  case 2:
    for (i = 0; i < N; i++) {
      const double x = 2.0 * M_PI * i / (double) N;
      for (j = 0; j < N; j++) {
	const double y = 2.0 * M_PI * j / (double) N;
	for (k = 0; k < N; k++) {
	  const double z = 2.0 * M_PI * k / (double) N;
	  
	  u[k + N * (j + i * N)] =  sin(z) * cos(x);
	  v[k + N * (j + i * N)] =  0.0;
	  w[k + N * (j + i * N)] = -cos(z) * sin(x);
	}
      }
    }
    break;

  default :
    message ("Taylor2D", "permutation key must be 0, 1 or 2", ERROR);
    break;
  }

  return;
}


void Taylor2D_error (CVF IC, const int* Dim, const Param* I, const int code)
/* ------------------------------------------------------------------------- *
 * Replace the velocity field by its error at the time indicated in I.
 * Code indicates which velocity component is zero.
 *
 * Velocity components are supplied in physical space.
 * ------------------------------------------------------------------------- */
{
  const double    decay = exp (-2.0 * I -> time / I -> Re);
  const int       N     = Dim[1];
  double          x, y, z;
  real            uvw;
  register real  *u = & IC[1][0][0][0].Re;
  register real  *v = & IC[2][0][0][0].Re;
  register real  *w = & IC[3][0][0][0].Re;
  register int    i, j, k;

  /* -- Fill up the cube. */

  switch (code) {
  case 0:
    for (i = 0; i < N; i++) {
      x = 2.0 * M_PI * i / (double) N;
      for (j = 0; j < N; j++) {
	y = 2.0 * M_PI * j / (double) N;
	for (k = 0; k < N; k++) {
	  z = 2.0 * M_PI * k / (double) N;

	  uvw = -cos(x) * sin(y) * decay;
	  u[k + N * (j + i * N)] -= uvw;

	  uvw =  sin(x) * cos(y) * decay;
	  v[k + N * (j + i * N)] -= uvw;
	}
      }
    }
    break;

  case 1:
    for (i = 0; i < N; i++) {
      x = 2.0 * M_PI * i / (double) N;
      for (j = 0; j < N; j++) {
	y = 2.0 * M_PI * j / (double) N;
	for (k = 0; k < N; k++) {
	  z = 2.0 * M_PI * k / (double) N;
	  
	  uvw = -cos(y) * sin(z) * decay;
	  v[k + N * (j + i * N)] -= uvw;

	  uvw =  sin(y) * cos(z) * decay;
	  w[k + N * (j + i * N)] -= uvw;
	}
      }
    }
    break;

  case 2:
    for (i = 0; i < N; i++) {
      x = 2.0 * M_PI * i / (double) N;
      for (j = 0; j < N; j++) {
	y = 2.0 * M_PI * j / (double) N;
	for (k = 0; k < N; k++) {
	  z = 2.0 * M_PI * k / (double) N;

	  uvw =  sin(z) * cos(x) * decay;
	  u[k + N * (j + i * N)] -= uvw;

	  uvw = -cos(z) * sin(x) * decay;
	  w[k + N * (j + i * N)] -= uvw;
	}
      }
    }
    break;

  default :
    message ("Taylor2D_error", "permutation must be 0, 1 or 2", ERROR);
    break;
  }
}


void  TaylorGreen (CVF IC, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Generate initial conditions of the 3D Taylor--Green vortex, in
 * PHYSICAL space.
 *
 *   u =  sin x cos y cos z
 *   v = -cos x sin y cos z
 *   w =  0
 * ------------------------------------------------------------------------- */
{
  const int N    = Dim[1];
  const int Npts = Dim[1] * Dim[2] * Dim[3];

  register real *u = & IC[1][0][0][0].Re;
  register real *v = & IC[2][0][0][0].Re;
  register real *w = & IC[3][0][0][0].Re;

  register int i, j, k;

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	u[k + N * (j + i * N)] =  sin(x) * cos(y) * cos(z);
	v[k + N * (j + i * N)] = -cos(x) * sin(y) * cos(z);
	w[k + N * (j + i * N)] =  0.0;
      }
    }
  }
}


void  TaylorGreenNL_error (CVF IC, const int* Dim)
/* ------------------------------------------------------------------------- *
 * Generate the nonlinear terms in the Navier--Stokes equations for 
 * the initial conditions of the 3D Taylor--Green vortex, in
 * PHYSICAL space.  Add them to input, which should be the negative of the
 * nonlinear terms.
 *
 *   u =  sin x cos y cos z
 *   v = -cos x sin y cos z
 *   w =  0
 *
 *   d(uu)/dx = 2 cos^2(y) sin(x) cos(x)            cos^2(z)
 *   d(uv)/dy = (sin^2(y) - cos^2(y)) cos(x) sin(x) cos^2(z)
 *   d(vu)/dx = (sin^2(x) - cos^2(x)) cos(y) sin(y) cos^2(z)
 *   d(vv)/dy = 2 cos^2(x) cos(y) sin(y)            cos^2(z)
 * ------------------------------------------------------------------------- */
{
  const    int   N    = Dim[1];
  const    int   Npts = Dim[1] * Dim[2] * Dim[3];
  register real *u    = &IC[1][0][0][0].Re;
  register real *v    = &IC[2][0][0][0].Re;
  register real *w    = &IC[3][0][0][0].Re;
  register int   i, j, k;
  real           UUx, UVy, VUx, VVy, cz2;

  for (i = 0; i < N; i++) {
    const double x = 2.0 * M_PI * i / (double) N;
    for (j = 0; j < N; j++) {
      const double y = 2.0 * M_PI * j / (double) N;
      for (k = 0; k < N; k++) {
	const double z = 2.0 * M_PI * k / (double) N;
	
	UUx = 2.0 * cos(x) * sin(x) * SQR(cos(y))           * SQR(cos(z));
	UVy = (SQR(sin(y)) - SQR(cos(y))) * cos(x) * sin(x) * SQR(cos(z));
	VUx = (SQR(sin(x)) - SQR(cos(x))) * cos(y) * sin(y) * SQR(cos(z));
	VVy = 2.0 * cos(y) * sin(y) * SQR(cos(x))           * SQR(cos(z));

	u[k + N * (j + i * N)] += (UUx + UVy);
	v[k + N * (j + i * N)] += (VUx + VVy);
      }
    }
  }
}


real  Brachet (const real t)
/* ------------------------------------------------------------------------- *
 * Return the Generalized Enstrophy of order 1 for the inviscid Taylor-
 * Green vortex, as estimated in Ref [5] (Table 5).
 *
 * There seems to be a normalizing factor adrift here.  The normalizing
 * factor has been adjusted from 0.5 (as per the paper) to 2.0
 * so that the series result agrees with the computation at time zero.
 * ------------------------------------------------------------------------- */
{
  register int   r;
  const    int   Ntab = 41;
  double         t2r, omega, t2 = SQR (t);
  static double  A[]  = {
     0.75000000000000000000000000000000E+00,
     0.78124999999999999999999999999999E-01,
     0.59185606060606060606060606060601E-02,
    -0.27843606425347977381470204436711E-03,
     0.61048213608414040976589547497500E-04,
    -0.96361343406379096996080171679830E-05,
     0.12376451432939409337723444617600E-05,
    -0.13001673399485842050532230751800E-06,
     0.11507664940214605511186628852000E-07,
    -0.40689927930502918781244653498000E-09,
    -0.15722647048717447228204712700000E-09,
     0.56457081462142300009592373250000E-10,
    -0.13696928782856063280154415300000E-10,
     0.29920198954228640442791535000000E-11,
    -0.61991963183363252307722688000000E-12,
     0.12450945726913139691007882000000E-12,
    -0.24774479828999593771452870000000E-13,
     0.49252841340330882023774810000000E-14,
    -0.97719041829306565745462900000000E-15,
     0.19375444767106104429043700000000E-15,
    -0.38505990845347852748213000000000E-16,
     0.76701287226861612998270000000000E-17,
    -0.15301069631879348281350000000000E-17,
     0.30585009158637265652730000000000E-18,
    -0.61279386047170433163800000000000E-19,
     0.12302069798992531554700000000000E-19,
    -0.24743842723469146812000000000000E-20,
     0.49871740662854581144000000000000E-21,
    -0.10071233331935838450000000000000E-21,
     0.20373933972348344470000000000000E-22,
    -0.41289304714858419800000000000000E-23,
     0.83820874296160258400000000000000E-24,
    -0.17043406330511116000000000000000E-24,
     0.34707362620876885000000000000000E-25,
    -0.70784176151137040000000000000000E-26,
     0.14456490670569290000000000000000E-26,
    -0.29564188944707500000000000000000E-27,
     0.60538014508785700000000000000000E-28,
    -0.12411584017048000000000000000000E-28,
     0.25476148544746000000000000000000E-29,
    -0.52351228744960000000000000000000E-30 };

  omega = A[0];
  t2r   = t2;
  for (r = 1; r < Ntab; r++) {
    omega += t2r * A[r];
    t2r   *= t2;
  }
  omega *= 0.5;

  return  omega;
}
