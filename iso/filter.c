/*****************************************************************************
 * filter.c: Fourier-space filter functions.
 *
 * Copyright (C) 1999 Hugh Blackburn
 *
 * $Id$
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


void ispectrum (CVF        U  ,
		const real TKE,
		const real k0 ,
		real       (*Ek)(const real, const real, const real))
/* ------------------------------------------------------------------------- *
 * Supposing the spectrum of U as input is Gaussian white noise, filter
 * so that the spherically-averaged output energy spectrum is
 *
 *   E(k) = given function of TKE, k0, and k
 * ------------------------------------------------------------------------- */
{
  register int  c, k1, k2, k3, b1, b2;
  register real k, tp;

  for (c = 1; c <= 3; c++) {
    U[c][0][0][0].Re = U[c][0][0][0].Im = 0.0;

    for (k1 = 1; k1 < K; k1++) {
      b1 = N - k1;
      k  = sqrt(k1 * k1);
      tp = sqrt (Ek (TKE, k0, k)) / k;
      U[c][k1][ 0][ 0].Re *= tp;
      U[c][k1][ 0][ 0].Im *= tp;
      U[c][ 0][k1][ 0].Re *= tp;
      U[c][ 0][k1][ 0].Im *= tp;
      U[c][ 0][ 0][k1].Re *= tp;
      U[c][ 0][ 0][k1].Im *= tp;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	k  = sqrt (k1 * k1 + k2 * k2);
	tp = sqrt (Ek (TKE, k0, k)) / k;
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
	  k  = sqrt (k1 * k1 + k2 * k2 + k3 * k3);
	  tp = sqrt (Ek (TKE, k0, k)) / k;
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


real unknown01 (const real TKE,
		const real k0 ,
		const real k  )
/* ------------------------------------------------------------------------- *
 * Return spectral magnitude at wavenumber k for spectrum
 *
 *   E(k) = 16 sqrt (2/PI) u0^2 / k0 (k/k0)^4 exp [ -2 (k/k0)^2 ]
 *
 * which (reportedly) has TKE = 3/2 u0^2, epsilon = 15/4 mu u0^2 k0^2,
 * and Taylor microscale lambda = 2 / k0.
 * ------------------------------------------------------------------------- */
{
  const real u02  = 2.0 / 3.0 * TKE;
  const real A    = 16.0 * sqrt (2.0 / M_PI) * u02 / k0;
  const real kok0 = k / k0;

  return A * SQR(SQR(kok0)) * exp (-2.0 * SQR(kok0));
}


real vonKarman (const real TKE,
		const real k0 ,
		const real k  )
/* ------------------------------------------------------------------------- *
 * Return spectral magnitude at wavenumber k for von Karman spectrum
 *
 *   E(k) = 2  TKE  L Cvk (k L)^4 / [1 + (k L)^2)]^p
 *
 * where L = 2PI / k0, p = 17/6 gives -5/3 law at high k, and
 * Cvk = Gamma(p)/[Gamma(5/2)Gamma(p-5/2)] = 1.64646372716 for p = 17/6.
 * ------------------------------------------------------------------------- */
{
  const real L   = 1.0 / (2.0 * M_PI * k0);
  const real kL  = k * L;
  const real p   = 17.0 / 6.0;
  const real Cvk = 1.64646372716;

  return 2.0 * TKE * Cvk * L * pow (kL, 4.0) / pow ((1.0 + SQR(kL)), p);
}

