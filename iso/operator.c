/*****************************************************************************
 * operator.c: generate one field from another.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define ROTATE(Z) tempRe = (Z).Re; (Z).Re = -(Z).Im; (Z).Im = tempRe


void deriv (const CVF V        ,
	    const int Component,
	    CF        D        ,
	    const int Direction)
/* ------------------------------------------------------------------------- *
 * Given a field of complex velocity Fourier coefficients (V), return in D
 * the Fourier coefficients of the derivatives of a specific Component with
 * respect to a specific Direction.
 * ------------------------------------------------------------------------- */
{
  register int  k1, k2, k3, b1, b2;
  register real tempRe;
  register CF   Ui = V[Component];

  D[ 0][ 0][ 0].Re = D[ 0][ 0][ 0].Im = 0.0;

  switch (Direction) {
  case 1:

    for (k1 = 1; k1 < K; k1++) {	/* Axes. */
      D[k1][ 0][ 0].Re = k1 * Ui[k1][ 0][ 0].Re;
      D[k1][ 0][ 0].Im = k1 * Ui[k1][ 0][ 0].Im;
      ROTATE (D[k1][ 0][ 0]);

      D[ 0][k1][ 0].Re = D[ 0][k1][ 0].Im = 0.0;
      D[ 0][ 0][k1].Re = D[ 0][ 0][k1].Im = 0.0;
    }

    for (k1 = 1; k1 < K; k1++) {	/* Faces. */
      b1 = N - k1;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	D[ 0][k1][k2].Re = D[ 0][k1][k2].Im = 0.0;
	D[ 0][b1][k2].Re = D[ 0][b1][k2].Im = 0.0;

	D[k1][ 0][k2].Re =  k1 * Ui[k1][ 0][k2].Re;
	D[k1][ 0][k2].Im =  k1 * Ui[k1][ 0][k2].Im;
	D[b1][ 0][k2].Re = -k1 * Ui[b1][ 0][k2].Re;
	D[b1][ 0][k2].Im = -k1 * Ui[b1][ 0][k2].Im;
	ROTATE (D[k1][ 0][k2]);
	ROTATE (D[b1][ 0][k2]);

	D[k1][k2][ 0].Re =  k1 * Ui[k1][k2][ 0].Re;
	D[k1][k2][ 0].Im =  k1 * Ui[k1][k2][ 0].Im;
	D[b1][k2][ 0].Re = -k1 * Ui[b1][k2][ 0].Re;
	D[b1][k2][ 0].Im = -k1 * Ui[b1][k2][ 0].Im;
	ROTATE (D[k1][k2][ 0]);
	ROTATE (D[b1][k2][ 0]);
      }
    }
      
    for (k1 = 1; k1 < K; k1++) {	/* Internal. */
      b1 = N - k1;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  D[k1][k2][k3].Re =  k1 * Ui[k1][k2][k3].Re;
	  D[k1][k2][k3].Im =  k1 * Ui[k1][k2][k3].Im;
	  D[b1][k2][k3].Re = -k1 * Ui[b1][k2][k3].Re;
	  D[b1][k2][k3].Im = -k1 * Ui[b1][k2][k3].Im;
	  D[k1][b2][k3].Re =  k1 * Ui[k1][b2][k3].Re;
	  D[k1][b2][k3].Im =  k1 * Ui[k1][b2][k3].Im;
	  D[b1][b2][k3].Re = -k1 * Ui[b1][b2][k3].Re;
	  D[b1][b2][k3].Im = -k1 * Ui[b1][b2][k3].Im;

	  ROTATE (D[k1][k2][k3]);
	  ROTATE (D[b1][k2][k3]);
	  ROTATE (D[k1][b2][k3]);
	  ROTATE (D[b1][b2][k3]);
	}
      }
    }
    break;

  case 2:

    for (k1 = 1; k1 < K; k1++) {	/* Axes. */
      D[k1][ 0][ 0].Re = D[k1][ 0][ 0].Im = 0.0;

      D[ 0][k1][ 0].Re = k1 * Ui[ 0][k1][ 0].Re;
      D[ 0][k1][ 0].Im = k1 * Ui[ 0][k1][ 0].Im;
      ROTATE (D[ 0][k1][ 0]);

      D[ 0][ 0][k1].Re = D[ 0][ 0][k1].Im = 0.0;
    }

    for (k1 = 1; k1 < K; k1++) {	/* Faces. */
      b1 = N - k1;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	D[ 0][k1][k2].Re =  k1 * Ui[ 0][k1][k2].Re;
	D[ 0][k1][k2].Im =  k1 * Ui[ 0][k1][k2].Im;
	D[ 0][b1][k2].Re = -k1 * Ui[ 0][b1][k2].Re;
	D[ 0][b1][k2].Im = -k1 * Ui[ 0][b1][k2].Im;
	ROTATE (D[ 0][k1][k2]);
	ROTATE (D[ 0][b1][k2]);
	
	D[k1][ 0][k2].Re = D[k1][ 0][k2].Im = 0.0;
	D[b1][ 0][k2].Re = D[b1][ 0][k2].Im = 0.0;

	D[k1][k2][ 0].Re =  k2 * Ui[k1][k2][ 0].Re;
	D[k1][k2][ 0].Im =  k2 * Ui[k1][k2][ 0].Im;
	D[b1][k2][ 0].Re =  k2 * Ui[b1][k2][ 0].Re;
	D[b1][k2][ 0].Im =  k2 * Ui[b1][k2][ 0].Im;
	ROTATE (D[k1][k2][ 0]);
	ROTATE (D[b1][k2][ 0]);
      }
    }
      
    for (k1 = 1; k1 < K; k1++) {	/* Internal. */
      b1 = N - k1;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  D[k1][k2][k3].Re =  k2 * Ui[k1][k2][k3].Re;
	  D[k1][k2][k3].Im =  k2 * Ui[k1][k2][k3].Im;
	  D[b1][k2][k3].Re =  k2 * Ui[b1][k2][k3].Re;
	  D[b1][k2][k3].Im =  k2 * Ui[b1][k2][k3].Im;
	  D[k1][b2][k3].Re = -k2 * Ui[k1][b2][k3].Re;
	  D[k1][b2][k3].Im = -k2 * Ui[k1][b2][k3].Im;
	  D[b1][b2][k3].Re = -k2 * Ui[b1][b2][k3].Re;
	  D[b1][b2][k3].Im = -k2 * Ui[b1][b2][k3].Im;

	  ROTATE (D[k1][k2][k3]);
	  ROTATE (D[b1][k2][k3]);
	  ROTATE (D[k1][b2][k3]);
	  ROTATE (D[b1][b2][k3]);
	}
      }
    }
    break;
    
  case 3:

    for (k1 = 1; k1 < K; k1++) {	/* Axes. */
      D[k1][ 0][ 0].Re = D[k1][ 0][ 0].Im = 0.0;
      D[ 0][k1][ 0].Re = D[ 0][k1][ 0].Im = 0.0;

      D[ 0][ 0][k1].Re = k1 * Ui[ 0][ 0][k1].Re;
      D[ 0][ 0][k1].Im = k1 * Ui[ 0][ 0][k1].Im;
      ROTATE (D[ 0][ 0][k1]);
    }

    for (k1 = 1; k1 < K; k1++) {	/* Faces. */
      b1 = N - k1;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	D[ 0][k1][k2].Re = k2 * Ui[ 0][k1][k2].Re;
	D[ 0][k1][k2].Im = k2 * Ui[ 0][k1][k2].Im;
	D[ 0][b1][k2].Re = k2 * Ui[ 0][b1][k2].Re;
	D[ 0][b1][k2].Im = k2 * Ui[ 0][b1][k2].Im;
	ROTATE (D[ 0][k1][k2]);
	ROTATE (D[ 0][b1][k2]);
	
	D[k1][ 0][k2].Re = k2 * Ui[k1][ 0][k2].Re;
	D[k1][ 0][k2].Im = k2 * Ui[k1][ 0][k2].Im; 
	D[b1][ 0][k2].Re = k2 * Ui[b1][ 0][k2].Re;
	D[b1][ 0][k2].Im = k2 * Ui[b1][ 0][k2].Im;
	ROTATE (D[k1][ 0][k2]);
	ROTATE (D[b1][ 0][k2]);

	D[k1][k2][ 0].Re = D[k1][k2][ 0].Im = 0.0;
	D[b1][k2][ 0].Re = D[b1][k2][ 0].Im = 0.0;
      }
    }
      
    for (k1 = 1; k1 < K; k1++) {	/* Internal. */
      b1 = N - k1;
      for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	b2 = N - k2;
	for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	  D[k1][k2][k3].Re = k3 * Ui[k1][k2][k3].Re;
	  D[k1][k2][k3].Im = k3 * Ui[k1][k2][k3].Im;
	  D[b1][k2][k3].Re = k3 * Ui[b1][k2][k3].Re;
	  D[b1][k2][k3].Im = k3 * Ui[b1][k2][k3].Im;
	  D[k1][b2][k3].Re = k3 * Ui[k1][b2][k3].Re;
	  D[k1][b2][k3].Im = k3 * Ui[k1][b2][k3].Im;
	  D[b1][b2][k3].Re = k3 * Ui[b1][b2][k3].Re;
	  D[b1][b2][k3].Im = k3 * Ui[b1][b2][k3].Im;

	  ROTATE (D[k1][k2][k3]);
	  ROTATE (D[b1][k2][k3]);
	  ROTATE (D[k1][b2][k3]);
	  ROTATE (D[b1][b2][k3]);
	}
      }
    }
    break;

  default:
    message ("deriv", "Direction for derivative not in [1..3]", ERROR);
    break;
  }
}


void curl (const CVF V ,
	   CVF       W ,
	   CF        WK)
/* ------------------------------------------------------------------------- *
 * Given the Fourier coefficients of a vector field V, compute the Fourier
 * coefficients of its curl, W.
 * ------------------------------------------------------------------------- */
{
  zeroVF (W);

  deriv (V, 3, WK, 2); setF (W[1], WK);
  deriv (V, 2, WK, 3); subF (W[1], WK);

  deriv (V, 1, WK, 3); setF (W[2], WK);
  deriv (V, 3, WK, 1); subF (W[2], WK);

  deriv (V, 2, WK, 1); setF (W[3], WK);
  deriv (V, 1, WK, 2); subF (W[3], WK);
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
 * Compute the isotropically-truncated convolution sum of U & V,
 * return in W.
 *
 * Inputs U & V are the physical space values of the two scalar fields
 * to be convolved, and U_ & V_ are their values represented on a grid
 * that is shifted in physical space.  The return value W, is the
 * Fourier transform of the convolution of U & V.
 *                                                                          
 * Method as decribed in Section 6 of article by Steven Orszag,
 * Stud. Appl.  Math., VLN4, Dec. 1971, pp. 293--327.  We are only
 * storing half the data of the full-complex storage version, since we
 * have the transforms of real data. Since we will never have to do
 * anything with the Nyquist-plane data (which will remain zero for
 * all time), the storage scheme for the conjugate-symmetric data is
 * straightforward.  The only complication arises in describing the
 * octodechedron of Fig. 1 of the article, taking account of the fact
 * that we store negative wavenumbers in the upper ends of our box.
 * 
 * It is assumed that the truncation to the octodecahedral space will
 * already have been carried out on inputs U & V, and we restrict as
 * many operations as possible to the truncated region to save
 * computation.
 * ------------------------------------------------------------------------- */
{
  const int        Npts = N * N * N;
  register int     _i, _j, _k, _l, _m, _n;
  register int     k1, b1, k2, b2, k3;
  register complex *A = &W [0][0][0],    *B  = &W_[0][0][0];
  register real    *u = &U [0][0][0].Re, *u_ = &U_[0][0][0].Re,
                   *v = &V [0][0][0].Re, *v_ = &V_[0][0][0].Re,
                   *w = &W [0][0][0].Re, *w_ = &W_[0][0][0].Re;
  
  for (k1 = 0; k1 < Npts; k1++) w_[k1] = u [k1] * v [k1];

  rc3DFT (W_, Wtab, FORWARD);

  A[0].Re = B[0].Re;
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    _i = RM(k1,0,0);
    _j = RM(0,k1,0);
    _k = RM(0,0,k1);
    A[_i] = B[_i];
    A[_j] = B[_j];
    A[_k] = B[_k];

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      _i = RM(0,k1,k2); _j = RM(0,b1,k2);
      _k = RM(k1,0,k2); _l = RM(b1,0,k2);
      _m = RM(k1,k2,0); _n = RM(b1,k2,0);
      A[_i] = B[_i];
      A[_j] = B[_j];
      A[_k] = B[_k];
      A[_l] = B[_l];
      A[_m] = B[_m];
      A[_n] = B[_n];

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	_i = RM(k1,k2,k3); _j = RM(b1,k2,k3);
	_k = RM(k1,b2,k3); _l = RM(b1,b2,k3);
	A[_i] = B[_i];
	A[_j] = B[_j];
	A[_k] = B[_k];
	A[_l] = B[_l];
      }
    }
  }

  for (k1 = 0; k1 < Npts; k1++) w_[k1] = u_[k1] * v_[k1];

  rc3DFT (W_, Wtab, FORWARD);
  shift  (W_, Stab, INVERSE);

  A[0].Re = 0.5 * (A[0].Re + B[0].Re);

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    _i = RM(k1,0,0);
    _j = RM(0,k1,0);
    _k = RM(0,0,k1);
    A[_i].Re = 0.5 * (A[_i].Re + B[_i].Re);
    A[_i].Im = 0.5 * (A[_i].Im + B[_i].Im);
    A[_j].Re = 0.5 * (A[_j].Re + B[_j].Re);
    A[_j].Im = 0.5 * (A[_j].Im + B[_j].Im);
    A[_k].Re = 0.5 * (A[_k].Re + B[_k].Re);
    A[_k].Im = 0.5 * (A[_k].Im + B[_k].Im);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;
      _i = RM(0,k1,k2); _j = RM(0,b1,k2);
      _k = RM(k1,0,k2); _l = RM(b1,0,k2);
      _m = RM(k1,k2,0); _n = RM(b1,k2,0);
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
	_i = RM(k1,k2,k3); _j = RM(b1,k2,k3);
	_k = RM(k1,b2,k3); _l = RM(b1,b2,k3);
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
  const int        SGN = (Drn == FORWARD) ? 1 : -1;
  register int     k1, b1, k2, b2, k3;
  register real    tempRe;
  register complex W, *u = &U[0][0][0];

  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    
    W = Stab[SGN*k1];
    SHIFT (u[RM(k1,0,0)], W);
    SHIFT (u[RM(0,k1,0)], W);
    SHIFT (u[RM(0,0,k1)], W);

    for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
      b2 = N - k2;

      W = Stab[SGN*(k1+k2)];
      SHIFT (u[RM(0,k1,k2)], W);
      SHIFT (u[RM(k1,0,k2)], W);
      SHIFT (u[RM(k1,k2,0)], W);

      W = Stab[SGN*(k2-k1)];
      SHIFT (u[RM(0,b1,k2)], W);
      SHIFT (u[RM(b1,0,k2)], W);
      SHIFT (u[RM(b1,k2,0)], W);

      for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	W = Stab[SGN*(k1+k2+k3)]; SHIFT (u[RM(k1,k2,k3)], W);
	W = Stab[SGN*(k3+k1-k2)]; SHIFT (u[RM(k1,b2,k3)], W);
	W = Stab[SGN*(k3+k2-k1)]; SHIFT (u[RM(b1,k2,k3)], W);
	W = Stab[SGN*(k3-k1-k2)]; SHIFT (u[RM(b1,b2,k3)], W);  
      }
    }
  }
}
