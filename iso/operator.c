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
