/*****************************************************************************
 * integrate.c: time advancement of Fourier coefficients.
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 * 
 * $Id$
 *****************************************************************************/

#include "iso.h"


void integrate (CVF          U,
                const CVF*   G,
		const Param* I)
/* ------------------------------------------------------------------------- *
 * Update the velocity field Fourier coefficients.  Time advancement
 * is explicit, using the Adams-Bashforth family of schemes, up to
 * AB3.  An integrating-factor method is used for the viscous terms
 * (strictly, I_factor below is the inverse of the factor).
 *
 * See Canuto  et al. ([4]) \S\,4.4.2.
 * ------------------------------------------------------------------------- */
{
  const int     order = CLAMP (I -> step + 1, 1, I -> norder);
  const real    Neg_nu_dt = -I -> dt * I -> kinvis;
  register int  c, k1, b1, k2, b2, k3;
  register real kSqrd, I_factor, C1, C2, C3;
  register CF   G0, G1, G2, V;

  if (order == 1) {			/* -- Take an Euler step. */

    for (c = 1; c <= 3; c++) {
      G0 = G[0][c];
      V  = U[c];

      for (k1 = 1; k1 < K; k1++) {
	b1       = N - k1;
	kSqrd    = k1 * k1;	/* Axes. */
	I_factor = exp (kSqrd * Neg_nu_dt);
	C1       = I -> dt * I_factor;
	
	V[k1][ 0][ 0].Re *= I_factor;
	V[k1][ 0][ 0].Im *= I_factor;
	V[ 0][k1][ 0].Re *= I_factor;
	V[ 0][k1][ 0].Im *= I_factor;
	V[ 0][ 0][k1].Re *= I_factor;
	V[ 0][ 0][k1].Im *= I_factor;

	V[k1][ 0][ 0].Re += C1 * G0[k1][ 0][ 0].Re;
	V[k1][ 0][ 0].Im += C1 * G0[k1][ 0][ 0].Im;
	V[ 0][k1][ 0].Re += C1 * G0[ 0][k1][ 0].Re;
	V[ 0][k1][ 0].Im += C1 * G0[ 0][k1][ 0].Im;
	V[ 0][ 0][k1].Re += C1 * G0[ 0][ 0][k1].Re;
	V[ 0][ 0][k1].Im += C1 * G0[ 0][ 0][k1].Im;

	for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	  b2       = N - k2;
	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = exp (kSqrd * Neg_nu_dt);
	  C1       = I -> dt * I_factor;
	  
	  V[ 0][k1][k2].Re *= I_factor;   /* i = 0 */
	  V[ 0][k1][k2].Im *= I_factor;
	  V[ 0][b1][k2].Re *= I_factor;
	  V[ 0][b1][k2].Im *= I_factor;

	  V[ 0][k1][k2].Re += C1 * G0[ 0][k1][k2].Re;
	  V[ 0][k1][k2].Im += C1 * G0[ 0][k1][k2].Im;
	  V[ 0][b1][k2].Re += C1 * G0[ 0][b1][k2].Re; 
	  V[ 0][b1][k2].Im += C1 * G0[ 0][b1][k2].Im; 

	  V[k1][ 0][k2].Re *= I_factor;  /* j = 0 */
	  V[k1][ 0][k2].Im *= I_factor;
	  V[b1][ 0][k2].Re *= I_factor;
	  V[b1][ 0][k2].Im *= I_factor;

	  V[k1][ 0][k2].Re += C1 * G0[k1][ 0][k2].Re;
	  V[k1][ 0][k2].Im += C1 * G0[k1][ 0][k2].Im;
	  V[b1][ 0][k2].Re += C1 * G0[b1][ 0][k2].Re; 
	  V[b1][ 0][k2].Im += C1 * G0[b1][ 0][k2].Im; 

	  V[k1][k2][ 0].Re *= I_factor;  /* k = 0 */
	  V[k1][k2][ 0].Im *= I_factor;
	  V[b1][k2][ 0].Re *= I_factor;
	  V[b1][k2][ 0].Im *= I_factor;

	  V[k1][k2][ 0].Re += C1 * G0[k1][k2][ 0].Re;
	  V[k1][k2][ 0].Im += C1 * G0[k1][k2][ 0].Im;
	  V[b1][k2][ 0].Re += C1 * G0[b1][k2][ 0].Re; 
	  V[b1][k2][ 0].Im += C1 * G0[b1][k2][ 0].Im; 

	  for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++) {
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = exp (kSqrd * Neg_nu_dt);
	    C1       = I -> dt * I_factor;
	    
	    V[k1][k2][k3].Re *= I_factor;
	    V[k1][k2][k3].Im *= I_factor;
	    V[b1][k2][k3].Re *= I_factor;
	    V[b1][k2][k3].Im *= I_factor;
	    V[k1][b2][k3].Re *= I_factor;
	    V[k1][b2][k3].Im *= I_factor;
	    V[b1][b2][k3].Re *= I_factor;
	    V[b1][b2][k3].Im *= I_factor;
	    
	    V[k1][k2][k3].Re += C1 * G0[k1][k2][k3].Re; 
	    V[k1][k2][k3].Im += C1 * G0[k1][k2][k3].Im;
	    V[b1][k2][k3].Re += C1 * G0[b1][k2][k3].Re; 
	    V[b1][k2][k3].Im += C1 * G0[b1][k2][k3].Im;
	    V[k1][b2][k3].Re += C1 * G0[k1][b2][k3].Re; 
	    V[k1][b2][k3].Im += C1 * G0[k1][b2][k3].Im;
	    V[b1][b2][k3].Re += C1 * G0[b1][b2][k3].Re; 
	    V[b1][b2][k3].Im += C1 * G0[b1][b2][k3].Im;
	  }
	}
      }
    }

  } else if (order == 2) {	/* -- Adams--Bashforth 2 timestep. */

    for (c = 1; c <= 3; c++) {
      G0 = G[0][c];
      G1 = G[1][c];
      V  = U[c];

      for (k1 = 1; k1 < K; k1++) {
	b1       = N - k1;
	kSqrd    = k1 * k1;	/* Axes. */
	I_factor = exp (kSqrd * Neg_nu_dt);
	C1       = 0.5 * I -> dt * I_factor;
	C2       = -I_factor * C1;
	C1      *= 3.0;
	
	V[k1][ 0][ 0].Re *= I_factor;
	V[k1][ 0][ 0].Im *= I_factor;
	V[ 0][k1][ 0].Re *= I_factor;
	V[ 0][k1][ 0].Im *= I_factor;
	V[ 0][ 0][k1].Re *= I_factor;
	V[ 0][ 0][k1].Im *= I_factor;
	
	V[k1][ 0][ 0].Re += C1 * G0[k1][ 0][ 0].Re +
	                    C2 * G1[k1][ 0][ 0].Re ;
	V[k1][ 0][ 0].Im += C1 * G0[k1][ 0][ 0].Im +
	                    C2 * G1[k1][ 0][ 0].Im ;
	V[ 0][k1][ 0].Re += C1 * G0[ 0][k1][ 0].Re +
	                    C2 * G1[ 0][k1][ 0].Re ;
	V[ 0][k1][ 0].Im += C1 * G0[ 0][k1][ 0].Im +
	                    C2 * G1[ 0][k1][ 0].Im ;
	V[ 0][ 0][k1].Re += C1 * G0[ 0][ 0][k1].Re +
	                    C2 * G1[ 0][ 0][k1].Re ;
	V[ 0][ 0][k1].Im += C1 * G0[ 0][ 0][k1].Im +
	                    C2 * G1[ 0][ 0][k1].Im ;

	for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	  b2       = N-k2;
	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = exp (kSqrd * Neg_nu_dt);
	  C1       = 0.5 * I -> dt * I_factor;
	  C2       = -I_factor * C1;
	  C1      *= 3.0;
	  
	  V[ 0][k1][k2].Re *= I_factor;   /* i = 0 */
	  V[ 0][k1][k2].Im *= I_factor;
	  V[ 0][b1][k2].Re *= I_factor;
	  V[ 0][b1][k2].Im *= I_factor;
	  
	  V[ 0][k1][k2].Re += C1 * G0[ 0][k1][k2].Re +
	                      C2 * G1[ 0][k1][k2].Re ;
	  V[ 0][k1][k2].Im += C1 * G0[ 0][k1][k2].Im +
	                      C2 * G1[ 0][k1][k2].Im ;
	  V[ 0][b1][k2].Re += C1 * G0[ 0][b1][k2].Re +
	                      C2 * G1[ 0][b1][k2].Re ;
	  V[ 0][b1][k2].Im += C1 * G0[ 0][b1][k2].Im +
	                      C2 * G1[ 0][b1][k2].Im ;
	  
	  V[k1][ 0][k2].Re *= I_factor;  /* j = 0 */
	  V[k1][ 0][k2].Im *= I_factor;
	  V[b1][ 0][k2].Re *= I_factor;
	  V[b1][ 0][k2].Im *= I_factor;
	  
	  V[k1][ 0][k2].Re += C1 * G0[k1][ 0][k2].Re +
	                      C2 * G1[k1][ 0][k2].Re ;
	  V[k1][ 0][k2].Im += C1 * G0[k1][ 0][k2].Im +
	                      C2 * G1[k1][ 0][k2].Im ;
	  V[b1][ 0][k2].Re += C1 * G0[b1][ 0][k2].Re +
	                      C2 * G1[b1][ 0][k2].Re ;
	  V[b1][ 0][k2].Im += C1 * G0[b1][ 0][k2].Im +
	                      C2 * G1[b1][ 0][k2].Im ;
	  
	  V[k1][k2][ 0].Re *= I_factor;  /* k = 0 */
	  V[k1][k2][ 0].Im *= I_factor;
	  V[b1][k2][ 0].Re *= I_factor;
	  V[b1][k2][ 0].Im *= I_factor;
	  
	  V[k1][k2][ 0].Re += C1 * G0[k1][k2][ 0].Re +
	                      C2 * G1[k1][k2][ 0].Re ;
	  V[k1][k2][ 0].Im += C1 * G0[k1][k2][ 0].Im +
	                      C2 * G1[k1][k2][ 0].Im ;
	  V[b1][k2][ 0].Re += C1 * G0[b1][k2][ 0].Re +
	                      C2 * G1[b1][k2][ 0].Re ;
	  V[b1][k2][ 0].Im += C1 * G0[b1][k2][ 0].Im +
	                      C2 * G1[b1][k2][ 0].Re ;

	  for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++){
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = exp (kSqrd * Neg_nu_dt);
	    C1       = 0.5 * I -> dt * I_factor;
	    C2       = -I_factor * C1;
	    C1      *= 3.0;
	    
	    V[k1][k2][k3].Re *= I_factor;
	    V[k1][k2][k3].Im *= I_factor;
	    V[b1][k2][k3].Re *= I_factor;
	    V[b1][k2][k3].Im *= I_factor;
	    V[k1][b2][k3].Re *= I_factor;
	    V[k1][b2][k3].Im *= I_factor;
	    V[b1][b2][k3].Re *= I_factor;
	    V[b1][b2][k3].Im *= I_factor;
	    
	    V[k1][k2][k3].Re += C1 * G0[k1][k2][k3].Re +
	                        C2 * G1[k1][k2][k3].Re ;
	    V[k1][k2][k3].Im += C1 * G0[k1][k2][k3].Im +
	                        C2 * G1[k1][k2][k3].Im ;
	    V[b1][k2][k3].Re += C1 * G0[b1][k2][k3].Re +
	                        C2 * G1[b1][k2][k3].Re ;
	    V[b1][k2][k3].Im += C1 * G0[b1][k2][k3].Im +
	                        C2 * G1[b1][k2][k3].Im ;
	    V[k1][b2][k3].Re += C1 * G0[k1][b2][k3].Re +
	                        C2 * G1[k1][b2][k3].Re ;
	    V[k1][b2][k3].Im += C1 * G0[k1][b2][k3].Im +
	                        C2 * G1[k1][b2][k3].Im ;
	    V[b1][b2][k3].Re += C1 * G0[b1][b2][k3].Re +
	                        C2 * G1[b1][b2][k3].Re ;
	    V[b1][b2][k3].Im += C1 * G0[b1][b2][k3].Im +
	                        C2 * G1[b1][b2][k3].Im ;
	  } 
	} 
      }
    }

  } else {			/*  -- AB3 step. */

    for (c = 1; c <= 3; c++) {
      G0 = G[0][c];
      G1 = G[1][c];
      G2 = G[2][c];
      V  = U[c];

      for (k1 = 1; k1 < K; k1++) {
	b1       = N - k1;
	kSqrd    = k1 * k1;	/* Axes. */
	I_factor = exp (kSqrd * Neg_nu_dt);
	C1       = I -> dt * I_factor;
	C2       = I_factor * C1;
	C3       = I_factor * C2;
	C1      *=  23.0 / 12.0;
	C2      *= -16.0 / 12.0;
	C3      *=   5.0 / 12.0;
	
	V[k1][ 0][ 0].Re *= I_factor;
	V[k1][ 0][ 0].Im *= I_factor;
	V[ 0][k1][ 0].Re *= I_factor;
	V[ 0][k1][ 0].Im *= I_factor;
	V[ 0][ 0][k1].Re *= I_factor;
	V[ 0][ 0][k1].Im *= I_factor;
	
	V[k1][ 0][ 0].Re += C1 * G0[k1][ 0][ 0].Re +
	                    C2 * G1[k1][ 0][ 0].Re +
	                    C3 * G2[k1][ 0][ 0].Re ;
	V[k1][ 0][ 0].Im += C1 * G0[k1][ 0][ 0].Im +
	                    C2 * G1[k1][ 0][ 0].Im +
	                    C3 * G1[k1][ 0][ 0].Im ;
	V[ 0][k1][ 0].Re += C1 * G0[ 0][k1][ 0].Re +
	                    C2 * G1[ 0][k1][ 0].Re +
	                    C3 * G1[ 0][k1][ 0].Re ;
	V[ 0][k1][ 0].Im += C1 * G0[ 0][k1][ 0].Im +
	                    C2 * G1[ 0][k1][ 0].Im +
	                    C3 * G1[ 0][k1][ 0].Im ;
	V[ 0][ 0][k1].Re += C1 * G0[ 0][ 0][k1].Re +
	                    C2 * G1[ 0][ 0][k1].Re +
	                    C3 * G1[ 0][ 0][k1].Re ;
	V[ 0][ 0][k1].Im += C1 * G0[ 0][ 0][k1].Im +
	                    C2 * G1[ 0][ 0][k1].Im +
	                    C3 * G1[ 0][ 0][k1].Im ;

	for (k2 = 1; k2 < K && k1+k2 INSIDE; k2++) {
	  b2       = N-k2;
	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = exp (kSqrd * Neg_nu_dt);
	  C1       = I -> dt * I_factor;
	  C2       = I_factor * C1;
	  C3       = I_factor * C2;
	  C1      *=  23.0 / 12.0;
	  C2      *= -16.0 / 12.0;
	  C3      *=   5.0 / 12.0;
	  
	  V[ 0][k1][k2].Re *= I_factor;   /* i = 0 */
	  V[ 0][k1][k2].Im *= I_factor;
	  V[ 0][b1][k2].Re *= I_factor;
	  V[ 0][b1][k2].Im *= I_factor;
	  
	  V[ 0][k1][k2].Re += C1 * G0[ 0][k1][k2].Re +
	                      C2 * G1[ 0][k1][k2].Re +
	                      C3 * G1[ 0][k1][k2].Re ;
	  V[ 0][k1][k2].Im += C1 * G0[ 0][k1][k2].Im +
	                      C2 * G1[ 0][k1][k2].Im +
	                      C3 * G1[ 0][k1][k2].Im ;
	  V[ 0][b1][k2].Re += C1 * G0[ 0][b1][k2].Re +
	                      C2 * G1[ 0][b1][k2].Re +
	                      C3 * G1[ 0][b1][k2].Re ;
	  V[ 0][b1][k2].Im += C1 * G0[ 0][b1][k2].Im +
	                      C2 * G1[ 0][b1][k2].Im +
	                      C3 * G1[ 0][b1][k2].Im ;
	  
	  V[k1][ 0][k2].Re *= I_factor;  /* j = 0 */
	  V[k1][ 0][k2].Im *= I_factor;
	  V[b1][ 0][k2].Re *= I_factor;
	  V[b1][ 0][k2].Im *= I_factor;
	  
	  V[k1][ 0][k2].Re += C1 * G0[k1][ 0][k2].Re +
	                      C2 * G1[k1][ 0][k2].Re +
	                      C3 * G1[k1][ 0][k2].Re ;
	  V[k1][ 0][k2].Im += C1 * G0[k1][ 0][k2].Im +
	                      C2 * G1[k1][ 0][k2].Im +
	                      C3 * G1[k1][ 0][k2].Im ;
	  V[b1][ 0][k2].Re += C1 * G0[b1][ 0][k2].Re +
	                      C2 * G1[b1][ 0][k2].Re +
	                      C3 * G1[b1][ 0][k2].Re ;
	  V[b1][ 0][k2].Im += C1 * G0[b1][ 0][k2].Im +
	                      C2 * G1[b1][ 0][k2].Im +
	                      C3 * G1[b1][ 0][k2].Im ;
	  
	  V[k1][k2][ 0].Re *= I_factor;  /* k = 0 */
	  V[k1][k2][ 0].Im *= I_factor;
	  V[b1][k2][ 0].Re *= I_factor;
	  V[b1][k2][ 0].Im *= I_factor;
	  
	  V[k1][k2][ 0].Re += C1 * G0[k1][k2][ 0].Re +
	                      C2 * G1[k1][k2][ 0].Re +
	                      C3 * G1[k1][k2][ 0].Re ;
	  V[k1][k2][ 0].Im += C1 * G0[k1][k2][ 0].Im +
	                      C2 * G1[k1][k2][ 0].Im +
	                      C3 * G1[k1][k2][ 0].Im ;
	  V[b1][k2][ 0].Re += C1 * G0[b1][k2][ 0].Re +
	                      C2 * G1[b1][k2][ 0].Re +
	                      C3 * G1[b1][k2][ 0].Re ;
	  V[b1][k2][ 0].Im += C1 * G0[b1][k2][ 0].Im +
	                      C2 * G1[b1][k2][ 0].Re +
	                      C3 * G1[b1][k2][ 0].Re ;

	  for (k3 = 1; k3 < K && k2+k3 INSIDE && k1+k3 INSIDE; k3++){
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = exp (kSqrd * Neg_nu_dt);
	    C1       = I -> dt * I_factor;
	    C2       = I_factor * C1;
	    C3       = I_factor * C2;
	    C1      *=  23.0 / 12.0;
	    C2      *= -16.0 / 12.0;
	    C3      *=   5.0 / 12.0;
	    
	    V[k1][k2][k3].Re *= I_factor;
	    V[k1][k2][k3].Im *= I_factor;
	    V[b1][k2][k3].Re *= I_factor;
	    V[b1][k2][k3].Im *= I_factor;
	    V[k1][b2][k3].Re *= I_factor;
	    V[k1][b2][k3].Im *= I_factor;
	    V[b1][b2][k3].Re *= I_factor;
	    V[b1][b2][k3].Im *= I_factor;
	    
	    V[k1][k2][k3].Re += C1 * G0[k1][k2][k3].Re +
	                        C2 * G1[k1][k2][k3].Re +
	                        C3 * G1[k1][k2][k3].Re ;
	    V[k1][k2][k3].Im += C1 * G0[k1][k2][k3].Im +
	                        C2 * G1[k1][k2][k3].Im +
	                        C3 * G1[k1][k2][k3].Im ;
	    V[b1][k2][k3].Re += C1 * G0[b1][k2][k3].Re +
	                        C2 * G1[b1][k2][k3].Re +
	                        C3 * G1[b1][k2][k3].Re ;
	    V[b1][k2][k3].Im += C1 * G0[b1][k2][k3].Im +
	                        C2 * G1[b1][k2][k3].Im +
	                        C3 * G1[b1][k2][k3].Im ;
	    V[k1][b2][k3].Re += C1 * G0[k1][b2][k3].Re +
	                        C2 * G1[k1][b2][k3].Re +
	                        C3 * G1[k1][b2][k3].Re ;
	    V[k1][b2][k3].Im += C1 * G0[k1][b2][k3].Im +
	                        C2 * G1[k1][b2][k3].Im +
	                        C3 * G1[k1][b2][k3].Im ;
	    V[b1][b2][k3].Re += C1 * G0[b1][b2][k3].Re +
	                        C2 * G1[b1][b2][k3].Re +
	                        C3 * G1[b1][b2][k3].Re ;
	    V[b1][b2][k3].Im += C1 * G0[b1][b2][k3].Im +
	                        C2 * G1[b1][b2][k3].Im +
	                        C3 * G1[b1][b2][k3].Im ;
	  } 
	} 
      }
    }
  } 
} 
