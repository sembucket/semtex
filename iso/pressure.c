/*****************************************************************************
 * pressure.c: compute pressure field.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 *
 * $Id$
 ****************************************************************************/

#include "iso.h"


void pressure(/* input     */ CVF      V  ,
                              CVF      G  ,
              /* output    */ CF       P  ,
              /* workspace */ CVF      WK ,
                              CF       VV ,
              /* using     */ complex* Wtab,
	                      complex* Stab)
/* ------------------------------------------------------------------------- *
 * From the Fourier coefficients of the velocity field (V) and the nonlinear
 * terms in the evolution equations (G), compute the Fourier coefficients of
 * the pressure field (P) (strictly, the pressure divided by the density).
 *
 * We use the identity:
 *    G_i = -i*k_i*P - i*k_j*FT[invFT(V_j)*invFT(V_i)]
 *
 * Since we can apparently arbitrarily use either:
 *        P = {i*G_1 - k_j*FT[invFT(V_j)*invFT(V_1)]} / k_1
 * or     P = {i*G_2 - k_j*FT[invFT(V_j)*invFT(V_2)]} / k_2 
 * or     P = {i*G_3 - k_j*FT[invFT(V_j)*invFT(V_3)]} / k_3
 * to compute P, we choose to use the first of these three, except on the
 * k1=0 face, where we use the second equation.  On the k3-axis, we use the
 * third equation.
 * Note that along the k1 axis, G_1 is zero, along the k2 axis G_2 is zero
 * and along the k3 axis G_3 is zero since the velocity Fourier components
 * are constrained to evolve in planes normal to the wavenumber vector at
 * each point in wavespace.
 *
 * The pressure derivative Fourier coefficients should come out parallel to
 * the wavenumber vector at each point in wavespace (no check made).
 * ------------------------------------------------------------------------- */
{
  int k1, k2, k3, b1, b2;

  truncate (P);			/* Set all high-wavenumber stuff zero. */
  P[0][0][0].Re = 0.0;		/* Zero the mean pressure.             */

  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      b2 = N-k2;
      P[0][k1][k2].Re = -G[2][0][k1][k2].Im;
      P[0][k1][k2].Im =  G[2][0][k1][k2].Re;
      P[0][b1][k2].Re = -G[2][0][b1][k2].Im;
      P[0][b1][k2].Im =  G[2][0][b1][k2].Re;

      P[k1][0][k2].Re = -G[1][k1][0][k2].Im;
      P[k1][0][k2].Im =  G[1][k1][0][k2].Re;
      P[b1][0][k2].Re = -G[1][b1][0][k2].Im;
      P[b1][0][k2].Im =  G[1][b1][0][k2].Re;

      P[k1][k2][0].Re = -G[1][k1][k2][0].Im;
      P[k1][k2][0].Im =  G[1][k1][k2][0].Re;
      P[b1][k2][0].Re = -G[1][b1][k2][0].Im;
      P[b1][k2][0].Im =  G[1][b1][k2][0].Re;

      for (k3=1; k3<K && k2+k3INSIDE && k1+k3INSIDE; k3++) {
	P[k1][k2][k3].Re = -G[1][k1][k2][k2].Im;
	P[k1][k2][k3].Im =  G[1][k1][k2][k3].Re;
	P[b1][k2][k3].Re = -G[1][b1][k2][k3].Im;
	P[b1][k2][k3].Im =  G[1][b1][k2][k3].Re;
	P[k1][b2][k3].Re = -G[1][k1][b2][k3].Im;
	P[k1][b2][k3].Im =  G[1][k1][b2][k3].Re;
	P[b1][b2][k3].Re = -G[1][b1][b2][k3].Im;
	P[b1][b2][k3].Im = -G[1][b1][b2][k3].Re;
      }
    }
  }

  convolve(V[1], V[1], VV, WK, Wtab, Stab);
  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    P[k1][0][0].Re = -VV[k1][0][0].Re;
    P[k1][0][0].Im = -VV[k1][0][0].Im;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      b2 = N-k2;
      P[k1][0][k2].Re -= k1 * VV[k1][0][k2].Re;
      P[k1][0][k2].Im -= k1 * VV[k1][0][k2].Im;
      P[b1][0][k2].Re += k1 * VV[k1][0][k2].Re;
      P[b1][0][k2].Im += k1 * VV[k1][0][k2].Im;
      P[k1][k2][0].Re -= k1 * VV[k1][k2][0].Re;
      P[k1][k2][0].Im -= k1 * VV[k1][k2][0].Im;
      P[b1][k2][0].Re += k1 * VV[b1][k2][0].Re;
      P[b1][k2][0].Im += k1 * VV[b1][k2][0].Im;     
      for (k3=1; k3<K && k2+k3INSIDE && k1+k3INSIDE; k3++) {
	P[k1][k2][k3].Re -= k1 * VV[k1][k2][k2].Im;
	P[k1][k2][k3].Im -= k1 * VV[k1][k2][k3].Re;
	P[b1][k2][k3].Re += k1 * VV[b1][k2][k3].Im;
	P[b1][k2][k3].Im += k1 * VV[b1][k2][k3].Re;
	P[k1][b2][k3].Re -= k1 * VV[k1][b2][k3].Im;
	P[k1][b2][k3].Im -= k1 * VV[k1][b2][k3].Re;
	P[b1][b2][k3].Re += k1 * VV[b1][b2][k3].Im;
	P[b1][b2][k3].Im += k1 * VV[b1][b2][k3].Re;
      }
    }
  }
  
  convolve(V[1], V[2], VV, WK, Wtab, Stab);
  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      b2 = N-k2;
      P[k1][k2][0].Re -= k2 * VV[k1][k2][0].Re;
      P[k1][k2][0].Im -= k2 * VV[k1][k2][0].Im;
      P[b1][k2][0].Re -= k2 * VV[b1][k2][0].Re;
      P[b1][k2][0].Im -= k2 * VV[b1][k2][0].Im;
      for (k3=1; k3<K && k2+k3INSIDE && k1+k3INSIDE; k3++) {
	P[k1][k2][k3].Re -= k2 * VV[k1][k2][k2].Im;
	P[k1][k2][k3].Im -= k2 * VV[k1][k2][k3].Re;
	P[b1][k2][k3].Re -= k2 * VV[b1][k2][k3].Im;
	P[b1][k2][k3].Im -= k2 * VV[b1][k2][k3].Re;
	P[k1][b2][k3].Re += k2 * VV[k1][b2][k3].Im;
	P[k1][b2][k3].Im += k2 * VV[k1][b2][k3].Re;
	P[b1][b2][k3].Re += k2 * VV[b1][b2][k3].Im;
	P[b1][b2][k3].Im += k2 * VV[b1][b2][k3].Re;
      }
    }
  }
  
  convolve(V[1], V[3], VV, WK, Wtab, Stab);
  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      b2 = N-k2;
      P[k1][0][k2].Re -= k2 * VV[k1][0][k2].Re;
      P[k1][0][k2].Im -= k2 * VV[k1][0][k2].Im;
      P[b1][0][k2].Re -= k2 * VV[b1][0][k2].Re;
      P[b1][0][k2].Im -= k2 * VV[b1][0][k2].Im;
      for (k3=1; k3<K && k2+k3INSIDE && k1+k3INSIDE; k3++) {
	P[k1][k2][k3].Re -= k3 * VV[k1][k2][k2].Im;
	P[k1][k2][k3].Im -= k3 * VV[k1][k2][k3].Re;
	P[b1][k2][k3].Re -= k3 * VV[b1][k2][k3].Im;
	P[b1][k2][k3].Im -= k3 * VV[b1][k2][k3].Re;
	P[k1][b2][k3].Re -= k3 * VV[k1][b2][k3].Im;
	P[k1][b2][k3].Im -= k3 * VV[k1][b2][k3].Re;
	P[b1][b2][k3].Re -= k3 * VV[b1][b2][k3].Im;
	P[b1][b2][k3].Im -= k3 * VV[b1][b2][k3].Re;
      }
    }
  }
  
  convolve(V[2], V[2], VV, WK, Wtab, Stab);
  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    P[0][k1][0].Re = -VV[0][k1][0].Re;
    P[0][k1][0].Im = -VV[0][k1][0].Im;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      P[0][k1][k2].Re -= k1 * VV[0][k1][k2].Re;
      P[0][k1][k2].Im -= k1 * VV[0][k1][k2].Im;
      P[0][b1][k2].Re += k1 * VV[0][b1][k2].Re;
      P[0][b1][k2].Im += k1 * VV[0][b1][k2].Im;
    }
  }
  
  convolve(V[2], V[3], VV, WK, Wtab, Stab);
  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      P[0][k1][k2].Re -= k2 * VV[0][k1][k2].Re;
      P[0][k1][k2].Im -= k2 * VV[0][k1][k2].Im;
      P[0][b1][k2].Re -= k2 * VV[0][b1][k2].Re;
      P[0][b1][k2].Im -= k2 * VV[0][b1][k2].Im;
    }
  }  
  
  convolve(V[3], V[3], VV, WK, Wtab, Stab);
  for (k1=1; k1<K; k1++) {
    P[0][0][k1].Re = -VV[0][0][k1].Re;
    P[0][0][k1].Im = -VV[0][0][k1].Im;
  }

  for (k1=1; k1<K; k1++) {
    b1 = N-k1;
    for (k2=1; k2<K && k1+k2INSIDE; k2++) {
      b2 = N-k2;
      P[0][k1][k2].Re /=  k1;
      P[0][k1][k2].Im /=  k1;
      P[0][b1][k2].Re /= -k1;
      P[0][b1][k2].Im /= -k1;
      P[k1][0][k2].Re /=  k1;
      P[k1][0][k2].Im /=  k1;
      P[b1][0][k2].Re /= -k1;
      P[b1][0][k2].Im /= -k1;
      P[k1][k2][0].Re /=  k1;
      P[k1][k2][0].Im /=  k1;
      P[b1][k2][0].Re /= -k1;
      P[b1][k2][0].Im /= -k1;
      for (k3=1; k3<K && k2+k3INSIDE && k1+k3INSIDE; k3++) {
	P[k1][k2][k3].Re /=  k1; 
	P[k1][k2][k3].Im /=  k1; 
	P[b1][k2][k3].Re /= -k1; 
	P[b1][k2][k3].Im /= -k1; 
	P[k1][b2][k3].Re /=  k1; 
	P[k1][b2][k3].Im /=  k1;
	P[b1][b2][k3].Re /= -k1; 
	P[b1][b2][k3].Im /= -k1; 
      }
    }
  }
}
  
