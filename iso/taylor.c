/*****************************************************************************
 * Generate initial conditions corresponding to the 2D Taylor Flow
 *
 * u = -cos(x) sin (y) exp(-2 PI^2 \nu t)
 * v =  sin(x) cos (y) exp(-2 PI^2 \nu t)
 * w =  0
 * p = -0.25 [cos(2x) + cos(2y)] exp(-4 PI^2 \nu t)
 * 
 * Input code = 0, 1, 2 generates a cyclic permutation of the velocity.
 *
 * Also: calculate the error in a field at any time.
 * 
 * $Id$
 *****************************************************************************/

#include "globals.h"


void Taylor2D (const int*         Dim ,
	       CVF  IC  ,
	       const int             code)
/* -------------------------------------------------------------------------
 * Generate initial conditions for t = 0, components in physical space.
 * Return prescaled for FFT.
 * ------------------------------------------------------------------------- */
{
  const int N    = Dim[1];
  const int Npts = Dim[1] * Dim[2] * Dim[3];

  /* -- Scaling factor for the FFT. */
  
  const real DFT_Fact = 1.0 / Npts;

  /* -- Fast pointers to the data. */

  register real *u = & IC[1][0][0][0].Re;
  register real *v = & IC[2][0][0][0].Re;
  register real *w = & IC[3][0][0][0].Re;

  register int i, j, k;

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

    for (i = 0; i < Npts*2; i++) {
      u[i] *= DFT_Fact;
      v[i] *= DFT_Fact;
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

    for (i = 0; i < Npts*2; i++) {
      v[i] *= DFT_Fact;
      w[i] *= DFT_Fact;
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

    for (i = 0; i < Npts*2; i++) {
      u[i] *= DFT_Fact;
      w[i] *= DFT_Fact;
    }

    break;
  default :
    message ("Taylor2D", "permutation key must be 0, 1 or 2", ERROR);
    break;
  }

  return;
}


void Taylor2D_error (const int*         Dim ,
		     CVF  IC  ,
		     const header*         I   ,
		     const int             code)
/* ------------------------------------------------------------------------
 * Replace the velocity field by its error at the time indicated by
 * header information.  Code indicates which velocity component is zero.
 *
 * Velocity components are supplied in physical space.
 * ------------------------------------------------------------------------ */
{
  const double  t     = I -> N_Step * I -> Delta_T;
  const double  decay = exp (-2.0 * M_PI * M_PI * I -> K_Visc * t);
  const int     N     = Dim[1];
  double        x, y, z;
  real          uvw;

  /* -- Fast pointers to data. */

  register real *u = & IC[1][0][0][0].Re;
  register real *v = & IC[2][0][0][0].Re;
  register real *w = & IC[3][0][0][0].Re;

  register int i, j, k;

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
