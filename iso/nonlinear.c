/*****************************************************************************
 * nonlinear.c: all computations to produce nonlinear terms in Navier--Stokes.
 *
 * Copyright (C) 1992, 1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define PEQ(Z1,c,Z2)  (Z1)->Re += (c) * (Z2)->Im; (Z1)->Im -= (c) * (Z2)->Re


void nonlinear (/* input     */  CVF             U   ,
		/* output    */  CVF             G   ,
		/* workspace */  CF              F   ,
		                 CF              F_  ,
		                 CVF             U_  ,
		/* using     */  const complex*  Wtab,
		                 const complex*  Stab)
/* ------------------------------------------------------------------------- *
 * Compute the Fourier transform of the nonlinear product -d(UjUm)/dxm,
 * Gj = -ikmFjm, using the pseudospectral method described by Orszag [1].
 *
 * In the old version, the terms were projected onto a divergence-free
 * space on the fly but now this is pulled out to be a separate step.
 *
 * The contributions are calculated only in the octodecahedral (isotropi-
 * cally truncated) space described by Orszag.
 *
 * Below is a underside view of the Fourier 1/2-box, showing Nyquist planes
 * and faces. also the labeling scheme for faces (b<==>negative wavenumber).
 * Note the box contains only positive k3 wavenumbers.
 *
 *                      | k1
 *                      |
 *                      |
 *	                +-------------+-+----------+
 *                     /|             | |          |
 *                    / |             | |          |
 *                   /  |     bk0     |N|    N     |
 *                  /   |             | |          |
 *                 /    |             | |          |
 *                +     |             | |          |
 *                | b0k +-------------+-+----------+
 *                |    /|      N      |N|    N     |
 *                |   / +-------------+-+----------+
 *                |  / /|             | |          |
 *                | /N/ |             | |          |
 *                |/ /  |             | |          |
 *                + /   |     kk0     |N|    N     |
 *                |/    |             | |          |
 *                +     |             | |          |
 *                |     |             | |          |
 *                | k0k |             | |          |
 *                |     +-------------+-+----------+---  k2
 *                |    /             / /          /
 *                |   /             / /          /
 *                |  /     0kk     /N/    0bk   /
 *                | /             / /          /
 *                |/             / /          /
 *                +-------------+-+--------- +
 *               /
 *              /
 *             / k3
 *
 * ------------------------------------------------------------------------- */
{
  register int     c, k1, b1, k2, b2, k3;
  register complex *f, *g;

  /* -- Make velocity field in PHYSICAL space on unshifted & shifted grid
   *    for alias control.
   */

  for (c = 1; c <= 3; c++) {
    copyF  (U_[c], U[c]);
    shift  (U_[c], Stab, FORWARD);
    rc3DFT (U [c], Wtab, INVERSE);
    rc3DFT (U_[c], Wtab, INVERSE);
  }

  /* -- Convolve u_hat[1] with itself to make F11 and distribute.
   *    G1 += -i k1 F11.
   */

  convolve (U[1], U[1], U_[1], U_[1], F, F_, Wtab, Stab);
  
  for (k1 = 1; k1 < K; k1++) {

    /* -- Axes.  Contributions on k1 axis only. */

    b1 = N - k1;

    f = F[k1][ 0]; g = G[1][k1][ 0]; PEQ (g,  k1, f);
    f = F[b1][ 0]; g = G[1][b1][ 0]; PEQ (g, -k1, f);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {

      /* -- Faces.  Contributions on k1--k2 & k1--k3 faces. */

      b2 = N - k2;

      f = F[k1][k2];    g = G[1][k1][k2];    PEQ (g,  k1, f);
      f = F[b1][k2];    g = G[1][b1][k2];    PEQ (g, -k1, f);
      f = F[k1][ 0]+k2; g = G[1][k1][ 0]+k2; PEQ (g,  k1, f);
      f = F[b1][ 0]+k2; g = G[1][b1][ 0]+k2; PEQ (g, -k1, f);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	/* -- Interior contributions. */

	f = F[k1][k2]+k3; g = G[1][k1][k2]+k3; PEQ (g,  k1, f);
	f = F[k1][b2]+k3; g = G[1][k1][b2]+k3; PEQ (g,  k1, f);
	f = F[b1][k2]+k3; g = G[1][b1][k2]+k3; PEQ (g, -k1, f);
	f = F[b1][b2]+k3; g = G[1][b1][b2]+k3; PEQ (g, -k1, f);

      }
    }
  }

  /* -- Convolve u_hat[2] with itself to make F22.
   *    G2 += -i k2 F22.
   */

  convolve (U[2], U[2], U_[2], U_[2], F,  F_, Wtab, Stab);
  
  for (k1 = 1; k1 < K; k1++) {

    /* -- Axes.  Contributions on k2 axis only. */

    b1 = N - k1;

    f = F[ 0][k1]; g = G[2][ 0][k1]; PEQ (g,  k1, f);
    f = F[ 0][b1]; g = G[2][ 0][b1]; PEQ (g, -k1, f);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {

      /* -- Faces.  Contributions on k1--k2 & k2--k3 faces only. */

      b2 = N - k2;

      f = F[k1][k2];    g = G[2][k1][k2];    PEQ (g,  k2, f);
      f = F[b1][k2];    g = G[2][b1][k2];    PEQ (g,  k2, f);
      f = F[ 0][k2]+k1; g = G[2][ 0][k2]+k1; PEQ (g,  k2, f);
      f = F[ 0][b2]+k1; g = G[2][ 0][b2]+k1; PEQ (g, -k2, f);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	/* -- Interior contributions. */

	f = F[k1][k2]+k3; g = G[2][k1][k2]+k3; PEQ (g,  k2, f);
	f = F[k1][b2]+k3; g = G[2][k1][b2]+k3; PEQ (g, -k2, f);
	f = F[b1][k2]+k3; g = G[2][b1][k2]+k3; PEQ (g,  k2, f);
	f = F[b1][b2]+k3; g = G[2][b1][b2]+k3; PEQ (g, -k2, f);

      }
    }
  }
  
  /* -- Convolve u_hat[3] with itself to make F33 and distribute.
   *    G3 += -i k3 F33.
   */

  convolve (U[3], U[3], U_[3], U_[3], F,  F_, Wtab, Stab);
  
  for (k1 = 1; k1 < K; k1++) {

    /* -- Axes.  Contributions on k3 axis only. */

    b1 = N - k1;

    f = F[ 0][ 0]+k1; g = G[3][ 0][ 0]+k1; PEQ (g,  k1, f);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {

      /* -- Faces.  Contributions on k1--k3 & k2--k3 faces only. */

      b2 = N - k2;

      f = F[k1][ 0]+k2; g = G[3][k1][ 0]+k2; PEQ (g,  k2, f);
      f = F[b1][ 0]+k2; g = G[3][b1][ 0]+k2; PEQ (g,  k2, f);
      f = F[ 0][k1]+k2; g = G[3][ 0][k1]+k2; PEQ (g,  k2, f);
      f = F[ 0][b1]+k2; g = G[3][ 0][b1]+k2; PEQ (g,  k2, f);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	f = F[k1][k2]+k3; g = G[3][k1][k2]+k3; PEQ (g,  k3, f);
	f = F[k1][b2]+k3; g = G[3][k1][b2]+k3; PEQ (g,  k3, f);
	f = F[b1][k2]+k3; g = G[3][b1][k2]+k3; PEQ (g,  k3, f);
	f = F[b1][b2]+k3; g = G[3][b1][b2]+k3; PEQ (g,  k3, f);

      }
    }
  }
  
  /* -- Convolve u_hat[1] with u_hat[2] to make F12 (and, by symmetry, F21).
   *    G2 += -i k1 F21;
   *    G1 += -i k2 F12;
   */

  convolve (U[1], U[2], U_[1], U_[2], F,  F_, Wtab, Stab);
  
  for (k1 = 1; k1 < K; k1++) {

    /* -- Axes. */

    b1 = N - k1;

    f = F[k1][ 0]; g = G[2][k1][ 0]; PEQ (g,  k1, f);
    f = F[b1][ 0]; g = G[2][b1][ 0]; PEQ (g, -k1, f);
    f = F[ 0][k1]; g = G[1][ 0][k1]; PEQ (g,  k1, f);
    f = F[ 0][b1]; g = G[1][ 0][b1]; PEQ (g, -k1, f);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {

      /* -- Faces.  G1 on k1--k2 & k1--k3; G2 on k1--k2 && k2--k3 */

      b2 = N - k2;

      f = F[k1][k2];    g = G[2][k1][k2];    PEQ (g,  k1, f);
      f = F[b1][k2];    g = G[2][b1][k2];    PEQ (g, -k1, f);
      f = F[k1][ 0]+k2; g = G[2][k1][ 0]+k2; PEQ (g,  k1, f);
      f = F[b1][ 0]+k2; g = G[2][b1][ 0]+k2; PEQ (g, -k1, f);
      
      f = F[k1][k2];    g = G[1][k1][k2];    PEQ (g,  k2, f);
      f = F[b1][k2];    g = G[1][b1][k2];    PEQ (g,  k2, f);
      f = F[ 0][k2]+k1; g = G[1][ 0][k2]+k1; PEQ (g,  k2, f);
      f = F[ 0][b2]+k1; g = G[1][ 0][b2]+k1; PEQ (g, -k2, f);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	/* -- Interior. */

	f = F[k1][k2]+k3; g = G[2][k1][k2]+k3; PEQ (g,  k1, f);
	f = F[k1][b2]+k3; g = G[2][k1][b2]+k3; PEQ (g,  k1, f);
	f = F[b1][k2]+k3; g = G[2][b1][k2]+k3; PEQ (g, -k1, f);
	f = F[b1][b2]+k3; g = G[2][b1][b2]+k3; PEQ (g, -k1, f);
	
	f = F[k1][k2]+k3; g = G[1][k1][k2]+k3; PEQ (g,  k2, f);
	f = F[b1][k2]+k3; g = G[1][b1][k2]+k3; PEQ (g,  k2, f);
	f = F[k1][b2]+k3; g = G[1][k1][b2]+k3; PEQ (g, -k2, f);
	f = F[b1][b2]+k3; g = G[1][b1][b2]+k3; PEQ (g, -k2, f);

      }
    }
  }
  
  /* -- Convolve u_hat[2] with u_hat[3] to make F23 (and, by symmetry, F32).
   *    G3 += -i k2 F32;
   *    G2 += -i k3 F23;
   */

  convolve (U[2], U[3], U_[2], U_[3], F,  F_, Wtab, Stab);
  
  for (k1 = 1; k1 < K; k1++) {

    /* -- Axes: contributions on k2 & k3 axes. */

    b1 = N - k1;

    f = F[ 0][k1]; g = G[3][ 0][k1]; PEQ (g,  k1, f);
    f = F[ 0][b1]; g = G[3][ 0][b1]; PEQ (g, -k1, f);

    f = F[ 0][ 0]+k1; g = G[2][ 0][ 0]+k1; PEQ (g,  k1, f);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {

      /* -- Faces.  G2 on k1--k2 & k2--k3; G3 on k1--k3 && k2--k3. */

      b2 = N - k2;

      f = F[k1][k2];    g = G[3][k1][k2];    PEQ (g,  k2, f);
      f = F[b1][k2];    g = G[3][b1][k2];    PEQ (g,  k2, f);
      f = F[ 0][k2]+k1; g = G[3][ 0][k2]+k1; PEQ (g,  k2, f);
      f = F[ 0][b2]+k1; g = G[3][ 0][b2]+k1; PEQ (g, -k2, f);

      f = F[ 0][k1]+k2; g = G[2][ 0][k1]+k2; PEQ (g,  k2, f);
      f = F[ 0][b1]+k2; g = G[2][ 0][b1]+k2; PEQ (g,  k2, f);
      f = F[k1][ 0]+k2; g = G[2][k1][ 0]+k2; PEQ (g,  k2, f);
      f = F[b1][ 0]+k2; g = G[2][b1][ 0]+k2; PEQ (g,  k2, f);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	f = F[k1][k2]+k3; g = G[3][k1][k2]+k3; PEQ (g,  k2, f);
	f = F[b1][k2]+k3; g = G[3][b1][k2]+k3; PEQ (g,  k2, f);
	f = F[k1][b2]+k3; g = G[3][k1][b2]+k3; PEQ (g, -k2, f);
	f = F[b1][b2]+k3; g = G[3][b1][b2]+k3; PEQ (g, -k2, f);
	
	f = F[k1][k2]+k3; g = G[2][k1][k2]+k3; PEQ (g,  k3, f);
	f = F[b1][k2]+k3; g = G[2][b1][k2]+k3; PEQ (g,  k3, f);
	f = F[k1][b2]+k3; g = G[2][k1][b2]+k3; PEQ (g,  k3, f);
	f = F[b1][b2]+k3; g = G[2][b1][b2]+k3; PEQ (g,  k3, f);

      }
    }
  }
  
  /* -- Convolve u_hat[1] with u_hat[3] to make F13 (and, by symmetry, F31).
   *    G3 += -i k1 F31;
   *    G1 += -i k3 F13;
   */

  convolve (U[1], U[3], U_[1], U_[3], F,  F_, Wtab, Stab);
  
  for (k1 = 1; k1 < K; k1++) {

    b1 = N - k1;

    f = F[k1][ 0];    g = G[3][k1][ 0];    PEQ (g,  k1, f);
    f = F[b1][ 0];    g = G[3][b1][ 0];    PEQ (g, -k1, f); 
    f = F[ 0][ 0]+k1; g = G[1][ 0][ 0]+k1; PEQ (g,  k1, f);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {

      b2 = N - k2;

      f = F[k1][k2];    g = G[3][k1][k2];    PEQ (g,  k1, f);
      f = F[b1][k2];    g = G[3][b1][k2];    PEQ (g, -k1, f);
      f = F[k1][ 0]+k2; g = G[3][k1][ 0]+k2; PEQ (g,  k1, f);
      f = F[b1][ 0]+k2; g = G[3][b1][ 0]+k2; PEQ (g, -k1, f);
      
      f = F[ 0][k1]+k2; g = G[1][ 0][k1]+k2; PEQ (g,  k2, f);
      f = F[ 0][b1]+k2; g = G[1][ 0][b1]+k2; PEQ (g,  k2, f);
      f = F[k1][ 0]+k2; g = G[1][k1][ 0]+k2; PEQ (g,  k2, f);
      f = F[b1][ 0]+k2; g = G[1][b1][ 0]+k2; PEQ (g,  k2, f);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {

	f = F[k1][k2]+k3; g = G[3][k1][k2]+k3; PEQ (g,  k1, f);
	f = F[k1][b2]+k3; g = G[3][k1][b2]+k3; PEQ (g,  k1, f);
	f = F[b1][k2]+k3; g = G[3][b1][k2]+k3; PEQ (g, -k1, f);
	f = F[b1][b2]+k3; g = G[3][b1][b2]+k3; PEQ (g, -k1, f);
	
	f = F[k1][k2]+k3; g = G[1][k1][k2]+k3; PEQ (g,  k3, f);
	f = F[b1][k2]+k3; g = G[1][b1][k2]+k3; PEQ (g,  k3, f);
	f = F[k1][b2]+k3; g = G[1][k1][b2]+k3; PEQ (g,  k3, f);
	f = F[b1][b2]+k3; g = G[1][b1][b2]+k3; PEQ (g,  k3, f);

      }
    }
  }

  /* -- Transform U back to FOURIER space. */

  for (c = 1; c <= 3; c++) {
    rc3DFT  (U[c], Wtab, FORWARD);
    scaleFT (U[c]);
    scaleFT (G[c]);
  }
}


void convolve (const CF       U   ,
	       const CF       V   ,
	       const CF       U_  ,
	       const CF       V_  ,
	       CF             W   ,
	       CF             W_  ,
	       const complex* Wtab,
	       const complex* Stab)
/* ------------------------------------------------------------------------- *
 * Compute the isotropically-truncated convolution sum of u & v, return in w.
 *                                                                          
 * Method as decribed in Section 6 of article by Steven Orszag, Stud. Appl. 
 * Math., VLN4, Dec. 1971, pp. 293--327.  We are only storing half the data 
 * of the full-complex storage version, since we have the transforms of real
 * data. Since we will never have to do anything with the Nyquist-plane data
 * (which will remain zero for all time), the storage scheme for the conjug-
 * ate-symmetric data is straightforward.  The only complication arises in  
 * describing the octodechedron of Fig. 1 of the article, taking account of 
 * the fact that we store negative wavenumbers in the upper ends of our box.
 *
 * The input data are the truncation will already have   
 * been carried out on u & v, and we try to restrict as many operations as  
 * possible to the truncated region to save computation.                    
 * ------------------------------------------------------------------------- */
{
  register int     _i, _j, _k, _l, _m, _n;
  register int     k1, b1, k2, b2, k3;
  register complex *A = &W[0][0][0], *B = &W_[0][0][0];
  const int        Npts = N * N * N;

  {
    register real
      *u = &U [0][0][0].Re, *u_ = &U_[0][0][0].Re,
      *v = &V [0][0][0].Re, *v_ = &V_[0][0][0].Re,
      *w = &W [0][0][0].Re, *w_ = &W_[0][0][0].Re;
    for (k1 = 0; k1 < Npts; k1++) w [k1] = u [k1] * v [k1];
    for (k1 = 0; k1 < Npts; k1++) w_[k1] = u_[k1] * v_[k1];
  }

  rc3DFT (W,  Wtab, FORWARD);
  rc3DFT (W_, Wtab, FORWARD);
  shift  (W_, Stab, INVERSE);

  A[0].Re = 0.5 * (A[0].Re + B[0].Re);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    _i = rm(k1,0,0);
    _j = rm(0,k1,0);
    _k = rm(0,0,k1);
    A[_i].Re = 0.5 * (A[_i].Re + B[_i].Re);
    A[_i].Im = 0.5 * (A[_i].Im + B[_i].Im);
    A[_j].Re = 0.5 * (A[_j].Re + B[_j].Re);
    A[_j].Im = 0.5 * (A[_j].Im + B[_j].Im);
    A[_k].Re = 0.5 * (A[_k].Re + B[_k].Re);
    A[_k].Im = 0.5 * (A[_k].Im + B[_k].Im);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      _i = rm(0,k1,k2); _j = rm(0,b1,k2);
      _k = rm(k1,0,k2); _l = rm(b1,0,k2);
      _m = rm(k1,k2,0); _n = rm(b1,k2,0);
      A[_i].Re = 0.5 * (A[_i].Re + B[_i].Re);
      A[_i].Im = 0.5 * (A[_i].Im + B[_i].Im);
      A[_j].Re = 0.5 * (A[_j].Re + B[_j].Re);
      A[_j].Im = 0.5 * (A[_j].Im + B[_j].Im);
      A[_k].Re = 0.5 * (A[_k].Re + B[_k].Re);
      A[_k].Im = 0.5 * (A[_k].Im + B[_k].Im);
      A[_l].Re = 0.5 * (A[_l].Re + B[_l].Re);
      A[_l].Im = 0.5 * (A[_l].Im + B[_l].Im);
      A[_m].Re = 0.5 * (A[_m].Re + B[_m].Re);
      A[_m].Im = 0.5 * (A[_m].Im + B[_m].Im);
      A[_n].Re = 0.5 * (A[_n].Re + B[_n].Re);
      A[_n].Im = 0.5 * (A[_n].Im + B[_n].Im);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	_i = rm(k1,k2,k3); _j = rm(b1,k2,k3);
	_k = rm(k1,b2,k3); _l = rm(b1,b2,k3);
	A[_i].Re = 0.5 * (A[_i].Re + B[_i].Re);
	A[_i].Im = 0.5 * (A[_i].Im + B[_i].Im);
	A[_j].Re = 0.5 * (A[_j].Re + B[_j].Re);
	A[_j].Im = 0.5 * (A[_j].Im + B[_j].Im);
	A[_k].Re = 0.5 * (A[_k].Re + B[_k].Re);
	A[_k].Im = 0.5 * (A[_k].Im + B[_k].Im);
	A[_l].Re = 0.5 * (A[_l].Re + B[_l].Re);
	A[_l].Im = 0.5 * (A[_l].Im + B[_l].Im);
      }
    }
  }
}


void shift (CF             U,
	    const complex* Stab,
	    const int      Drn)
/* ------------------------------------------------------------------------- *
 * Phase shift in FOURIER space <==> interpolate to shifted grid in PHYSICAL.
 * ------------------------------------------------------------------------- */
{
  register int     k1, b1, k2, b2, k3;
  register real    tempRe;
  register complex W, *u = &U[0][0][0];
  const int        SGN = (Drn == FORWARD) ? 1 : -1;

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    
    W = Stab[SGN*k1];
    SHIFT (u[rm(k1,0,0)], W);
    SHIFT (u[rm(0,k1,0)], W);
    SHIFT (u[rm(0,0,k1)], W);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;

      W = Stab[SGN*(k1+k2)];
      SHIFT (u[rm(0,k1,k2)], W);
      SHIFT (u[rm(k1,0,k2)], W);
      SHIFT (u[rm(k1,k2,0)], W);

      W = Stab[SGN*(k2-k1)];
      SHIFT (u[rm(0,b1,k2)], W);
      SHIFT (u[rm(b1,0,k2)], W);
      SHIFT (u[rm(b1,k2,0)], W);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	W = Stab[SGN*(k1+k2+k3)]; SHIFT (u[rm(k1,k2,k3)], W);
	W = Stab[SGN*(k3+k1-k2)]; SHIFT (u[rm(k1,b2,k3)], W);
	W = Stab[SGN*(k3+k2-k1)]; SHIFT (u[rm(b1,k2,k3)], W);
	W = Stab[SGN*(k3-k1-k2)]; SHIFT (u[rm(b1,b2,k3)], W);  
      }
    }
  }
}
