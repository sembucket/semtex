/*****************************************************************************
 * filter.c
 *
 * Routines for computing spectral filters.
 *****************************************************************************/

static char
RCSid_fourier[] = "$Id$";

#include <math.h>
#include <femdef.h>
#include <stdio.h>
#include <alplib.h>


void bvdFilter (const integer N     ,
		const integer order ,
		const integer lag   ,
		const real    attn  ,
		real*         filter)
// ---------------------------------------------------------------------------
// Load filter with the Boyd--Vandeven (i.e. erfc) filter [0, N] of given
// order (p) and lag (s).
//
// NB: N should be one less than the number of coefficients to
// which the filter will be applied.
// 
// Input parameter attn gives the attenuation at high
// wavenumbers. 0<=attn<=1, with attn = 1 giving complete attenuation.
//
// Reference:
// J.G. Levin & M. Iskandarani & D.B. Haidvogel (1997), A spectral
// filtering procedure for eddy-resolving simulations with a spectral
// element ocean model, JCP V137, 130--154.
// ---------------------------------------------------------------------------
{
  integer    i;
  real       arg, theta, chi, omega;
  const real EPS = EPSSP;

  for (i = 0; i < lag; i++)
    filter[i] = 1.0;
  
  for (i = lag; i <= N; i++) {
    theta = (i - lag) / (real) (N - lag);
    omega = fabs(theta) - 0.5;
    if ((fabs (theta - 0.5)) < EPS) 
      chi = 1.0;
    else {
      arg = 1.0 - 4.0 * SQR (omega);
      chi = sqrt (-log (arg) / (4.0 * SQR (omega)));
    }
    filter[i] = (1.0 - attn) + attn * 0.5 * erfc (2.0*sqrt(order)*chi*omega);
  }
}
