/*****************************************************************************
 * integrate.c: time advancement of Fourier coefficients.
 * 
 * $Id$
 *****************************************************************************/

#include "iso.h"


void  integrate (/* update */ CVF           U    ,
                 /* using  */ const CVF     G    ,
 		              const CVF     G_old,
                              const Param*  I    ,
                              const int*    Dim  )
/* ------------------------------------------------------------------------- *
 * Update the velocity field Fourier coefficients.  Time advancement is
 * explicit, using an Euler step at the first timestep, and Adams-Bashforth
 * second order subsequently.  An integrating-factor method is used for the
 * viscous terms (strictly, I_factor below is the inverse of the factor).
 *
 * See Canuto et al. ([4]) \S\,4.4.2.
 * ------------------------------------------------------------------------- */
{
  register int  c, k1, b1, k2, b2, k3;
  const    int  N         = Dim[1];
  const    int  K         = Dim[3];
  const    int  FOURKon3  = (4 * K) / 3;
  const    real Neg_nu_dt =  -I -> dt / I -> Re;
  real          kSqrd;

  
  if (! I -> step) {			/* -- Do an Euler step. */
    register real  I_factor, C1;

    for (c = 1; c <= 3; c++)
      for (k1 = 1; k1 < K; k1++) {
	b1 = N - k1;

	kSqrd    = k1 * k1;	/* Axes. */
	I_factor = exp (kSqrd * Neg_nu_dt);
	C1       = I -> dt * I_factor;
	
	U[c][k1][ 0][ 0].Re *= I_factor;
	U[c][k1][ 0][ 0].Im *= I_factor;
	U[c][ 0][k1][ 0].Re *= I_factor;
	U[c][ 0][k1][ 0].Im *= I_factor;
	U[c][ 0][ 0][k1].Re *= I_factor;
	U[c][ 0][ 0][k1].Im *= I_factor;

	U[c][k1][ 0][ 0].Re += C1 * G[c][k1][ 0][ 0].Re;
	U[c][k1][ 0][ 0].Im += C1 * G[c][k1][ 0][ 0].Im;
	U[c][ 0][k1][ 0].Re += C1 * G[c][ 0][k1][ 0].Re;
	U[c][ 0][k1][ 0].Im += C1 * G[c][ 0][k1][ 0].Im;
	U[c][ 0][ 0][k1].Re += C1 * G[c][ 0][ 0][k1].Re;
	U[c][ 0][ 0][k1].Im += C1 * G[c][ 0][ 0][k1].Im;

	for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
	  b2 = N - k2;
	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = exp (kSqrd * Neg_nu_dt);
	  C1       = I -> dt * I_factor;
	  
	  U[c][ 0][k1][k2].Re *= I_factor;   /* i = 0 */
	  U[c][ 0][k1][k2].Im *= I_factor;
	  U[c][ 0][b1][k2].Re *= I_factor;
	  U[c][ 0][b1][k2].Im *= I_factor;

	  U[c][ 0][k1][k2].Re += C1 * G[c][ 0][k1][k2].Re;
	  U[c][ 0][k1][k2].Im += C1 * G[c][ 0][k1][k2].Im;
	  U[c][ 0][b1][k2].Re += C1 * G[c][ 0][b1][k2].Re; 
	  U[c][ 0][b1][k2].Im += C1 * G[c][ 0][b1][k2].Im; 

	  U[c][k1][ 0][k2].Re *= I_factor;  /* j = 0 */
	  U[c][k1][ 0][k2].Im *= I_factor;
	  U[c][b1][ 0][k2].Re *= I_factor;
	  U[c][b1][ 0][k2].Im *= I_factor;

	  U[c][k1][ 0][k2].Re += C1 * G[c][k1][ 0][k2].Re;
	  U[c][k1][ 0][k2].Im += C1 * G[c][k1][ 0][k2].Im;
	  U[c][b1][ 0][k2].Re += C1 * G[c][b1][ 0][k2].Re; 
	  U[c][b1][ 0][k2].Im += C1 * G[c][b1][ 0][k2].Im; 

	  U[c][k1][k2][ 0].Re *= I_factor;  /* k = 0 */
	  U[c][k1][k2][ 0].Im *= I_factor;
	  U[c][b1][k2][ 0].Re *= I_factor;
	  U[c][b1][k2][ 0].Im *= I_factor;

	  U[c][k1][k2][ 0].Re += C1 * G[c][k1][k2][ 0].Re;
	  U[c][k1][k2][ 0].Im += C1 * G[c][k1][k2][ 0].Im;
	  U[c][b1][k2][ 0].Re += C1 * G[c][b1][k2][ 0].Re; 
	  U[c][b1][k2][ 0].Im += C1 * G[c][b1][k2][ 0].Im; 

	  for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++){
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = exp (kSqrd * Neg_nu_dt);
	    C1       = I -> dt * I_factor;
	    
	    U[c][k1][k2][k3].Re *= I_factor;
	    U[c][k1][k2][k3].Im *= I_factor;
	    U[c][b1][k2][k3].Re *= I_factor;
	    U[c][b1][k2][k3].Im *= I_factor;
	    U[c][k1][b2][k3].Re *= I_factor;
	    U[c][k1][b2][k3].Im *= I_factor;
	    U[c][b1][b2][k3].Re *= I_factor;
	    U[c][b1][b2][k3].Im *= I_factor;
	    
	    U[c][k1][k2][k3].Re += C1 * G[c][k1][k2][k3].Re; 
	    U[c][k1][k2][k3].Im += C1 * G[c][k1][k2][k3].Im;
	    U[c][b1][k2][k3].Re += C1 * G[c][b1][k2][k3].Re; 
	    U[c][b1][k2][k3].Im += C1 * G[c][b1][k2][k3].Im;
	    U[c][k1][b2][k3].Re += C1 * G[c][k1][b2][k3].Re; 
	    U[c][k1][b2][k3].Im += C1 * G[c][k1][b2][k3].Im;
	    U[c][b1][b2][k3].Re += C1 * G[c][b1][b2][k3].Re; 
	    U[c][b1][b2][k3].Im += C1 * G[c][b1][b2][k3].Im;
	  }
	}
      }
    
  } else {			/* -- Adams-Bashforth 2 timestep. */

    register real  I_factor, C1, C2;

    for (c = 1; c <= 3; c++)
      for (k1 = 1; k1 < K; k1++) {
	b1 = N - k1;

	kSqrd    = k1 * k1;	/* Axes. */
	I_factor = exp (kSqrd * Neg_nu_dt);
	C1       = 0.5 * I -> dt * I_factor;
	C2       = -I_factor * C1;
	C1      *= 3.0;
	
	U[c][k1][ 0][ 0].Re *= I_factor;
	U[c][k1][ 0][ 0].Im *= I_factor;
	U[c][ 0][k1][ 0].Re *= I_factor;
	U[c][ 0][k1][ 0].Im *= I_factor;
	U[c][ 0][ 0][k1].Re *= I_factor;
	U[c][ 0][ 0][k1].Im *= I_factor;
	
	U[c][k1][ 0][ 0].Re +=   C1 *     G[c][k1][ 0][ 0].Re
	                       + C2 * G_old[c][k1][ 0][ 0].Re;
	U[c][k1][ 0][ 0].Im +=   C1 *     G[c][k1][ 0][ 0].Im
	                       + C2 * G_old[c][k1][ 0][ 0].Im;
	U[c][ 0][k1][ 0].Re +=   C1 *     G[c][ 0][k1][ 0].Re
	                       + C2 * G_old[c][ 0][k1][ 0].Re;
	U[c][ 0][k1][ 0].Im +=   C1 *     G[c][ 0][k1][ 0].Im
	                       + C2 * G_old[c][ 0][k1][ 0].Im;
	U[c][ 0][ 0][k1].Re +=   C1 *     G[c][ 0][ 0][k1].Re
	                       + C2 * G_old[c][ 0][ 0][k1].Re;
	U[c][ 0][ 0][k1].Im +=   C1 *     G[c][ 0][ 0][k1].Im
	                       + C2 * G_old[c][ 0][ 0][k1].Im;

	for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
	  b2 = N-k2;

	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = exp (kSqrd * Neg_nu_dt);
	  C1       = 0.5 * I -> dt * I_factor;
	  C2       = -I_factor * C1;
	  C1      *= 3.0;
	  
	  U[c][ 0][k1][k2].Re *= I_factor;   /* i = 0 */
	  U[c][ 0][k1][k2].Im *= I_factor;
	  U[c][ 0][b1][k2].Re *= I_factor;
	  U[c][ 0][b1][k2].Im *= I_factor;
	  
	  U[c][ 0][k1][k2].Re +=   C1 *     G[c][ 0][k1][k2].Re
	                         + C2 * G_old[c][ 0][k1][k2].Re;
	  U[c][ 0][k1][k2].Im +=   C1 *     G[c][ 0][k1][k2].Im
	                         + C2 * G_old[c][ 0][k1][k2].Im;
	  U[c][ 0][b1][k2].Re +=   C1 *     G[c][ 0][b1][k2].Re
	                         + C2 * G_old[c][ 0][b1][k2].Re;
	  U[c][ 0][b1][k2].Im +=   C1 *     G[c][ 0][b1][k2].Im
	                         + C2 * G_old[c][ 0][b1][k2].Im;
	  
	  U[c][k1][ 0][k2].Re *= I_factor;  /* j = 0 */
	  U[c][k1][ 0][k2].Im *= I_factor;
	  U[c][b1][ 0][k2].Re *= I_factor;
	  U[c][b1][ 0][k2].Im *= I_factor;
	  
	  U[c][k1][ 0][k2].Re +=   C1 *     G[c][k1][ 0][k2].Re
	                         + C2 * G_old[c][k1][ 0][k2].Re;
	  U[c][k1][ 0][k2].Im +=   C1 *     G[c][k1][ 0][k2].Im
	                         + C2 * G_old[c][k1][ 0][k2].Im;
	  U[c][b1][ 0][k2].Re +=   C1 *     G[c][b1][ 0][k2].Re
	                         + C2 * G_old[c][b1][ 0][k2].Re;
	  U[c][b1][ 0][k2].Im +=   C1 *     G[c][b1][ 0][k2].Im
	                         + C2 * G_old[c][b1][ 0][k2].Im;
	  
	  U[c][k1][k2][ 0].Re *= I_factor;  /* k = 0 */
	  U[c][k1][k2][ 0].Im *= I_factor;
	  U[c][b1][k2][ 0].Re *= I_factor;
	  U[c][b1][k2][ 0].Im *= I_factor;
	  
	  U[c][k1][k2][ 0].Re +=   C1 *     G[c][k1][k2][ 0].Re
	                         + C2 * G_old[c][k1][k2][ 0].Re;
	  U[c][k1][k2][ 0].Im +=   C1 *     G[c][k1][k2][ 0].Im
	                         + C2 * G_old[c][k1][k2][ 0].Im;
	  U[c][b1][k2][ 0].Re +=   C1 *     G[c][b1][k2][ 0].Re 
	                         + C2 * G_old[c][b1][k2][ 0].Re;
	  U[c][b1][k2][ 0].Im +=   C1 *     G[c][b1][k2][ 0].Im
	                         + C2 * G_old[c][b1][k2][ 0].Re;

	  for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++){
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = exp (kSqrd * Neg_nu_dt);
	    C1       = 0.5 * I -> dt * I_factor;
	    C2       = -I_factor * C1;
	    C1      *= 3.0;
	    
	    U[c][k1][k2][k3].Re *= I_factor;
	    U[c][k1][k2][k3].Im *= I_factor;
	    U[c][b1][k2][k3].Re *= I_factor;
	    U[c][b1][k2][k3].Im *= I_factor;
	    U[c][k1][b2][k3].Re *= I_factor;
	    U[c][k1][b2][k3].Im *= I_factor;
	    U[c][b1][b2][k3].Re *= I_factor;
	    U[c][b1][b2][k3].Im *= I_factor;
	    
	    U[c][k1][k2][k3].Re +=   C1 *     G[c][k1][k2][k3].Re
	                           + C2 * G_old[c][k1][k2][k3].Re;
	    U[c][k1][k2][k3].Im +=   C1 *     G[c][k1][k2][k3].Im
	                           + C2 * G_old[c][k1][k2][k3].Im;
	    U[c][b1][k2][k3].Re +=   C1 *     G[c][b1][k2][k3].Re
	                           + C2 * G_old[c][b1][k2][k3].Re;
	    U[c][b1][k2][k3].Im +=   C1 *     G[c][b1][k2][k3].Im
	                           + C2 * G_old[c][b1][k2][k3].Im;
	    U[c][k1][b2][k3].Re +=   C1 *     G[c][k1][b2][k3].Re
	                           + C2 * G_old[c][k1][b2][k3].Re;
	    U[c][k1][b2][k3].Im +=   C1 *     G[c][k1][b2][k3].Im
	                           + C2 * G_old[c][k1][b2][k3].Im;
	    U[c][b1][b2][k3].Re +=   C1 *     G[c][b1][b2][k3].Re
	                           + C2 * G_old[c][b1][b2][k3].Re;
	    U[c][b1][b2][k3].Im +=   C1 *     G[c][b1][b2][k3].Im
	                           + C2 * G_old[c][b1][b2][k3].Im;
	  } 
	} 
      } 
  } 
} 
