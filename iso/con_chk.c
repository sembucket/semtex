/*****************************************************************************
 * con_chk.c: exercise 3D convolution routines, which should do
 * de-aliased convolution by working in truncated box, and with
 * shifted grids, against "plain vanilla" de-aliasing using zero
 * padding.
 *
 * Copyright (C) 1999 Hugh Blackburn
 * 
 * $Id$
 *****************************************************************************/

#include "iso.h"

int N, K, FourKon3;

#define SIZE 8

static char prog[] = "con_chk";


int main()
/* ------------------------------------------------------------------------- *
 * Initial (random) data sets are .
 * ------------------------------------------------------------------------- */
{
  CF      U, V, U_, V_, F, F_, UU, VV, WW;
  real    *u, *v, *w, emax, test;
  complex *Wtab, *Stab, *WTAB;
  char    s[STR_MAX];
  int     seed = 1, k1, k2, k3, b1, b2, B1, B2, L, i, j, k, ntot;

  /* -- Allocation. */

  N = SIZE; K = SIZE / 2; L = 2 * SIZE; FourKon3 = (4 * K) / 3;

  cbox (0, N-1, 0, N-1, 0, K-1, &U ); /* -- Standard sized fields. */
  cbox (0, N-1, 0, N-1, 0, K-1, &V );
  cbox (0, N-1, 0, N-1, 0, K-1, &U_);
  cbox (0, N-1, 0, N-1, 0, K-1, &V_);
  cbox (0, N-1, 0, N-1, 0, K-1, &F );
  cbox (0, N-1, 0, N-1, 0, K-1, &F_);

  cbox (0, L-1, 0, L-1, 0, N-1, &UU); /* -- Double-sized for zero-padding. */
  cbox (0, L-1, 0, L-1, 0, N-1, &VV);
  cbox (0, L-1, 0, L-1, 0, N-1, &WW);

  /* -- Set up tables. */

  Wtab = cvector (0, K-1);
  WTAB = cvector (0, N-1);
  Stab = cvector (-(N-1), N-1);
  preFFT   (Wtab, K);
  preFFT   (WTAB, N);
  preShift (Stab, N);

  /* -- Initialize U, V with random numbers. */

  zeroF (U);
  zeroF (V);

  u = &U[0][0][0].Re;
  v = &V[0][0][0].Re;

  ntot = N * N * N;
#if 0
  for (i = 0; i < ntot; i++) u[i] = ran2PI (&seed);
  for (i = 0; i < ntot; i++) v[i] = ran2PI (&seed);
#else
  for (i = 0; i < ntot; i++) u[i] = 0.0;
  for (i = 0; i < ntot; i++) v[i] = 0.0;
#endif

  U[1][0][0].Re = V[1][0][0].Re = 1.0;

#if 0
  truncateF (U);
  truncateF (V);
#endif

  /* -- Zero the padded data areas UU & VV. */

  u = &UU[0][0][0].Re;
  v = &VV[0][0][0].Re;

  ntot = L * L * L;

  for (i = 0; i < ntot; i++) u[i] = 0.0;
  for (i = 0; i < ntot; i++) v[i] = 0.0;

  /* -- Load the zero-padded data areas from U & V. */

  UU[ 0][ 0][ 0] = U[ 0][ 0][ 0];
  VV[ 0][ 0][ 0] = V[ 0][ 0][ 0];

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    B1 = L - k1;

    UU[k1][ 0][ 0] = U[k1][ 0][ 0];
    UU[ 0][k1][ 0] = U[ 0][k1][ 0];
    UU[ 0][ 0][k1] = U[ 0][ 0][k1];

    VV[k1][ 0][ 0] = V[k1][ 0][ 0];
    VV[ 0][k1][ 0] = V[ 0][k1][ 0];
    VV[ 0][ 0][k1] = V[ 0][ 0][k1];

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      B2 = L - k2;

      UU[ 0][k1][k2] = U[ 0][k1][k2];
      UU[ 0][B1][k2] = U[ 0][b1][k2];
      UU[k1][ 0][k2] = U[k1][ 0][k2];
      UU[B1][ 0][k2] = U[b1][ 0][k2];
      UU[k1][k2][ 0] = U[k1][k2][ 0];
      UU[B1][k2][ 0] = U[b1][k2][ 0];
      
      VV[ 0][k1][k2] = V[ 0][k1][k2];
      VV[ 0][B1][k2] = V[ 0][b1][k2];
      VV[k1][ 0][k2] = V[k1][ 0][k2];
      VV[B1][ 0][k2] = V[b1][ 0][k2];
      VV[k1][k2][ 0] = V[k1][k2][ 0];
      VV[B1][k2][ 0] = V[b1][k2][ 0];

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	UU[k1][k2][k3] = U[k1][k2][k3];
	UU[B1][k2][k3] = U[b1][k2][k3];
	UU[k1][B2][k3] = U[k1][b2][k3];
	UU[B1][B2][k3] = U[b1][b2][k3];

	VV[k1][k2][k3] = V[k1][k2][k3];
	VV[B1][k2][k3] = V[b1][k2][k3];
	VV[k1][B2][k3] = V[k1][b2][k3];
	VV[B1][B2][k3] = V[b1][b2][k3];
      }
    }
  }

  /* -- F <-- "shifted-grid + truncated" convolution of U & V. */
  
  copyF  (U_, U);
  copyF  (V_, V);
  shift  (U_, Stab, FORWARD);
  shift  (V_, Stab, FORWARD);
  rc3DFT (U,  Wtab, INVERSE);
  rc3DFT (V,  Wtab, INVERSE);
  rc3DFT (U_, Wtab, INVERSE);
  rc3DFT (V_, Wtab, INVERSE);

  convolve (U, V, U_, V_, F, F_, Wtab, Stab);

  rc3DFT  (F, Wtab, FORWARD);
  scaleFT (F);

  /* -- WW <-- "zero padded" convolution of UU & VV. */

  N = 2 * SIZE; K = SIZE; FourKon3 = (4 * K) / 3;

  rc3DFT (UU, WTAB, INVERSE);
  rc3DFT (VV, WTAB, INVERSE);

  u = &UU[0][0][0].Re;
  v = &VV[0][0][0].Re;
  w = &WW[0][0][0].Re;

  ntot = L * L * L;

  for (i = 0; i < ntot; i++) w[i] = u[i] * v[i];

  rc3DFT  (WW, WTAB, FORWARD);
  scaleFT (WW);

  N = SIZE;
  K = SIZE / 2;

  /* -- Subtract WW from F. */

  F[ 0][ 0][ 0].Re -= WW[ 0][ 0][ 0].Re;
  F[ 0][ 0][ 0].Im -= WW[ 0][ 0][ 0].Im;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    B1 = L - k1;

    F[k1][ 0][ 0].Re -= WW[k1][ 0][ 0].Re;
    F[ 0][k1][ 0].Re -= WW[ 0][k1][ 0].Re;
    F[ 0][ 0][k1].Re -= WW[ 0][ 0][k1].Re;

    F[k1][ 0][ 0].Im -= WW[k1][ 0][ 0].Im;
    F[ 0][k1][ 0].Im -= WW[ 0][k1][ 0].Im;
    F[ 0][ 0][k1].Im -= WW[ 0][ 0][k1].Im;

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      B2 = L - k2;

      F[ 0][k1][k2].Re -= WW[ 0][k1][k2].Re;
      F[ 0][b1][k2].Re -= WW[ 0][B1][k2].Re;
      F[k1][ 0][k2].Re -= WW[k1][ 0][k2].Re;
      F[b1][ 0][k2].Re -= WW[B1][ 0][k2].Re;
      F[k1][k2][ 0].Re -= WW[k1][k2][ 0].Re;
      F[b1][k2][ 0].Re -= WW[B1][k2][ 0].Re;

      F[ 0][k1][k2].Im -= WW[ 0][k1][k2].Im;
      F[ 0][b1][k2].Im -= WW[ 0][B1][k2].Im;
      F[k1][ 0][k2].Im -= WW[k1][ 0][k2].Im;
      F[b1][ 0][k2].Im -= WW[B1][ 0][k2].Im;
      F[k1][k2][ 0].Im -= WW[k1][k2][ 0].Im;
      F[b1][k2][ 0].Im -= WW[B1][k2][ 0].Im;

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	F[k1][k2][k3].Re -= WW[k1][k2][k3].Re;
	F[b1][k2][k3].Re -= WW[B1][k2][k3].Re;
	F[k1][b2][k3].Re -= WW[k1][B2][k3].Re;
	F[b1][b2][k3].Re -= WW[B1][B2][k3].Re;

	F[k1][k2][k3].Im -= WW[k1][k2][k3].Im;
	F[b1][k2][k3].Im -= WW[B1][k2][k3].Im;
	F[k1][b2][k3].Im -= WW[k1][B2][k3].Im;
	F[b1][b2][k3].Im -= WW[B1][B2][k3].Im;
      }
    }
  }

  /* -- Find biggest error. */

  emax = 0.0;

  emax = MAG(F[ 0][ 0][ 0]);
  i = j = k = 0;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    if ((test=MAG(F[k1][ 0][ 0]))>emax) {emax=test; i=k1; j= 0; k= 0;}
    if ((test=MAG(F[ 0][k1][ 0]))>emax) {emax=test; i= 0; j=k1; k= 0;}
    if ((test=MAG(F[ 0][ 0][k1]))>emax) {emax=test; i= 0; j= 0; k=k1;}
    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      if ((test=MAG(F[ 0][k1][k2]))>emax) {emax=test; i= 0; j=k1; k=k2;}
      if ((test=MAG(F[ 0][b1][k2]))>emax) {emax=test; i= 0; j=b1; k=k2;}
      if ((test=MAG(F[k1][ 0][k2]))>emax) {emax=test; i=k1; j= 0; k=k2;}
      if ((test=MAG(F[b1][ 0][k2]))>emax) {emax=test; i=b1; j= 0; k=k2;}
      if ((test=MAG(F[k1][k2][ 0]))>emax) {emax=test; i=k1; j=k2; k= 0;}
      if ((test=MAG(F[b1][k2][ 0]))>emax) {emax=test; i=b1; j=b2; k= 0;}
      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	if ((test=MAG(F[k1][k2][k3]))>emax) {emax=test; i=k1; j=k2; k=k3;}
	if ((test=MAG(F[b1][k2][k3]))>emax) {emax=test; i=b1; j=k2; k=k3;}
	if ((test=MAG(F[k1][b2][k3]))>emax) {emax=test; i=k1; j=b2; k=k3;}
	if ((test=MAG(F[b1][b2][k3]))>emax) {emax=test; i=b1; j=b2; k=k3;}
      }
    }
  }

  printf ("Maximum error %g at k1 = %3d  k2 = %3d  k3 = %3d\n", emax, i, j, k);

  return (EXIT_SUCCESS);
}


  
