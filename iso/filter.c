/*****************************************************************************
 * filter.c: Fourier-space filter functions.
 *
 * Copyright (C) 1999 Hugh Blackburn
 *****************************************************************************/

#include "iso.h"


void bvdFilter (const int  N     ,
		const int  order ,
		const int  lag   ,
		const real attn  ,
		real*      filter)
/* ------------------------------------------------------------------------- *
 * Load filter with the Boyd--Vandeven (i.e. erfc) filter [0, N] of given
 * order (p) and lag (s).
 *
 * NB: N should be one less than the number of coefficients to
 * which the filter will be applied.
 * 
 * Input parameter attn gives the attenuation at high
 * wavenumbers. 0<=attn<=1, with attn = 1 giving complete attenuation.
 *
 * Reference:
 * J.G. Levin & M. Iskandarani & D.B. Haidvogel (1997), A spectral
 * filtering procedure for eddy-resolving simulations with a spectral
 * element ocean model, JCP V137, 130--154.
 * ------------------------------------------------------------------------- */
{
  int        i;
  real       arg, theta, chi, omega;
  const real EPS = 6.0e-7;

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


void ispectrum (CVF        U ,
		const int  k0,
		const real u0)
/* ------------------------------------------------------------------------- *
 * Supposing the spectrum of U as input is Gaussian white noise, filter
 * so that the spherically-averaged output energy spectrum is
 *
 *   E(k) = 16 sqrt (2/PI) u0^2/k0 (k/k0)^4 exp [-2 (k/k0)^2]
 * 
 * Taking account of the normalisation of 2PI k^2 that is needed when
 * looking after the full spectrum, and the fact that we also need to
 * take the square root when multiplying our complex spectral
 * coefficients, this means we need to factor coefficients by
 *
 *  2 (2/PI)^(1/4) PI^(-1) u0 (k0)^(2/5) exp [-(k/k0)^2]
 * ------------------------------------------------------------------------- */
{
  const real    A   = pow (2, 1.25) * pow (M_PI, -0.75) * pow (k0, -1.5) * u0;
  const real    k02 = k0 * k0;
  register int  c, k1, k2, k3, b1, b2;
  register real kk, tp;

  for (c = 1; c <= 3; c++) {
    U[c][0][0][0].Re = U[c][0][0][0].Im = 0.0;

    for (k1 = 1; k1 < K; k1++) {
      b1 = N - k1;
      kk = k1 * k1;
      tp = A * k1 * exp (-kk/k02);
      U[c][k1][ 0][ 0].Re *= tp;
      U[c][k1][ 0][ 0].Im *= tp;
      U[c][ 0][k1][ 0].Re *= tp;
      U[c][ 0][k1][ 0].Im *= tp;
      U[c][ 0][ 0][k1].Re *= tp;
      U[c][ 0][ 0][k1].Im *= tp;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	kk = k1 * k1 + k2 * k2;
	tp = A * sqrt (kk) * exp (-kk/k02);
	U[c][ 0][k1][k2].Re *= tp;
	U[c][ 0][k1][k2].Im *= tp;
	U[c][ 0][b1][k2].Re *= tp;
	U[c][ 0][b1][k2].Im *= tp;
	U[c][k1][ 0][k2].Re *= tp;
	U[c][k1][ 0][k2].Im *= tp;
	U[c][b1][ 0][k2].Re *= tp;
	U[c][b1][ 0][k2].Im *= tp;
	U[c][k1][k2][ 0].Re *= tp;
	U[c][k1][k2][ 0].Im *= tp;
	U[c][b1][k2][ 0].Re *= tp;
	U[c][b1][k2][ 0].Im *= tp;
	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  kk = k1 * k1 + k2 * k2 + k3 * k3;
	  tp = A * sqrt (kk) * exp (-kk/k02);
	  U[c][k1][k2][k3].Re *= tp;
	  U[c][k1][k2][k3].Im *= tp;
	  U[c][b1][k2][k3].Re *= tp;
	  U[c][b1][k2][k3].Im *= tp;
	  U[c][k1][b2][k3].Re *= tp;
	  U[c][k1][b2][k3].Im *= tp;
	  U[c][b1][b2][k3].Re *= tp;
	  U[c][b1][b2][k3].Im *= tp;
	}
      }
    }
  }
}

