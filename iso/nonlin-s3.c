/*****************************************************************************
 * nonlin-s3.c: all computations to produce nonlinear terms in
 * Navier--Stokes.  This version incorporates Shah-Ferziger "S^3" SGSS
 * model.
 *
 * Copyright (C) 1992, 1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

#define PEQ(Z1,c,Z2) (Z1)->Re += (c) * (Z2)->Im; (Z1)->Im -= (c) * (Z2)->Re


void nonlinear (CVF            U   ,
		CVF            G   ,
		CVF            U_  ,
		CF             F   ,
		CF             F_  ,
		const complex* Wtab,
		const complex* Stab)
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
  static real      *expand = 0, *shrink = 0;

  if (!expand) {		/* -- Set up DFT filter coefficients.  Once. */
    k1 = K; k2 = k1 - 1;
    expand = (real*) malloc (2 * (k1 + 1) * sizeof (real));
    shrink = expand + k1 + 1;
    
    bvdFilter (k2, F_ORDER, (int) (F_ROLL * k2), F_ATTEN, shrink);
    for (c = 0; c <= k2; c++) expand[c] = 1.0 / shrink[c];
  }

  /* -- "Stimulate" the small scales uf U. */

  for (c = 1; c <= 3; c++) filterF (U[c], expand, NULL);

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

  /* -- "De-stimulate" the small scales of U, and G. */

  for (c = 1; c <= 3; c++) {
    filterF (U[c], shrink, 0);
    filterF (G[c], shrink, 0);
  }
}
