/*===========================================================================
 * RCS Information:
 * ----------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 *===========================================================================*/

#include "globals.h"

void integrate( /* update */ complex_vector_field  Vel,
                /* using  */ complex_vector_field  G,
 		             complex_vector_field  G_old,
                             header                Run_info,
                             ivector               Dimension,
                /* switch */ int                   Start)
/*===========================================================================*/
/* Update the velocity field Fourier coefficients.  Time advancement is      */
/* explicit, using an Euler step at the first timestep, and Adams-Bashforth  */
/* second order subsequently.  An integrating-factor method is used for the  */
/* viscous terms (strictly, I_factor below is the inverse of the factor).    */
/*===========================================================================*/
{
  int   c, k1, b1, k2, b2, k3, K, N, FOURKon3;
  float kSqrd, Neg_nu_dt, I_factor, C1, C2;


  N = Dimension[1];
  K = Dimension[3];
  FOURKon3 = (4*K) / 3;
  
  Neg_nu_dt = -Run_info.K_Visc * Run_info.Delta_T;


  if (Start) {
    /*-----------------------------------------------------------------------*/
    /* Do an Euler step.                                                     */
    /*-----------------------------------------------------------------------*/
    for (c=1; c<=3; c++)
      for (k1=1; k1<K; k1++) {
	b1 = N-k1;

	kSqrd    = k1*k1;	/* Axes. */
	I_factor = expf(kSqrd * Neg_nu_dt);
	C1       = Run_info.Delta_T * I_factor;
	
	Vel[c][k1][0][0].Re *= I_factor;
	Vel[c][k1][0][0].Im *= I_factor;
	Vel[c][0][k1][0].Re *= I_factor;
	Vel[c][0][k1][0].Im *= I_factor;
	Vel[c][0][0][k1].Re *= I_factor;
	Vel[c][0][0][k1].Im *= I_factor;
	
	Vel[c][k1][0][0].Re += C1 * G[c][k1][0][0].Re;
	Vel[c][k1][0][0].Im += C1 * G[c][k1][0][0].Im;
	Vel[c][0][k1][0].Re += C1 * G[c][0][k1][0].Re;
	Vel[c][0][k1][0].Im += C1 * G[c][0][k1][0].Im;
	Vel[c][0][0][k1].Re += C1 * G[c][0][0][k1].Re;
	Vel[c][0][0][k1].Im += C1 * G[c][0][0][k1].Im;

	for (k2=1; k2<K && k1+k2<=FOURKon3; k2++) {
	  b2 = N-k2;

	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = expf(kSqrd * Neg_nu_dt);
	  C1       = Run_info.Delta_T * I_factor;
	  
	  Vel[c][0][k1][k2].Re *= I_factor;   /* i = 0 */
	  Vel[c][0][k1][k2].Im *= I_factor;
	  Vel[c][0][b1][k2].Re *= I_factor;
	  Vel[c][0][b1][k2].Im *= I_factor;
	  
	  Vel[c][0][k1][k2].Re += C1 * G[c][0][k1][k2].Re;
	  Vel[c][0][k1][k2].Im += C1 * G[c][0][k1][k2].Im;
	  Vel[c][0][b1][k2].Re += C1 * G[c][0][b1][k2].Re; 
	  Vel[c][0][b1][k2].Im += C1 * G[c][0][b1][k2].Im; 
	  
	  Vel[c][k1][0][k2].Re *= I_factor;  /* j = 0 */
	  Vel[c][k1][0][k2].Im *= I_factor;
	  Vel[c][b1][0][k2].Re *= I_factor;
	  Vel[c][b1][0][k2].Im *= I_factor;
	  
	  Vel[c][k1][0][k2].Re += C1 * G[c][k1][0][k2].Re;
	  Vel[c][k1][0][k2].Im += C1 * G[c][k1][0][k2].Im;
	  Vel[c][b1][0][k2].Re += C1 * G[c][b1][0][k2].Re; 
	  Vel[c][b1][0][k2].Im += C1 * G[c][b1][0][k2].Im; 
	  
	  Vel[c][k1][k2][0].Re *= I_factor;  /* k = 0 */
	  Vel[c][k1][k2][0].Im *= I_factor;
	  Vel[c][b1][k2][0].Re *= I_factor;
	  Vel[c][b1][k2][0].Im *= I_factor;
	  
	  Vel[c][k1][k2][0].Re += C1 * G[c][k1][k2][0].Re;
	  Vel[c][k1][k2][0].Im += C1 * G[c][k1][k2][0].Im;
	  Vel[c][b1][k2][0].Re += C1 * G[c][b1][k2][0].Re; 
	  Vel[c][b1][k2][0].Im += C1 * G[c][b1][k2][0].Im; 

	  for (k3=1; k3<K && k2+k3<=FOURKon3 && k1+k3<=FOURKon3; k3++) {
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = expf(kSqrd * Neg_nu_dt);
	    C1       = Run_info.Delta_T * I_factor;
	    
	    Vel[c][k1][k2][k3].Re *= I_factor;
	    Vel[c][k1][k2][k3].Im *= I_factor;
	    Vel[c][b1][k2][k3].Re *= I_factor;
	    Vel[c][b1][k2][k3].Im *= I_factor;
	    Vel[c][k1][b2][k3].Re *= I_factor;
	    Vel[c][k1][b2][k3].Im *= I_factor;
	    Vel[c][b1][b2][k3].Re *= I_factor;
	    Vel[c][b1][b2][k3].Im *= I_factor;
	    
	    Vel[c][k1][k2][k3].Re += C1 * G[c][k1][k2][k3].Re; 
	    Vel[c][k1][k2][k3].Im += C1 * G[c][k1][k2][k3].Im;
	    Vel[c][b1][k2][k3].Re += C1 * G[c][b1][k2][k3].Re; 
	    Vel[c][b1][k2][k3].Im += C1 * G[c][b1][k2][k3].Im;
	    Vel[c][k1][b2][k3].Re += C1 * G[c][k1][b2][k3].Re; 
	    Vel[c][k1][b2][k3].Im += C1 * G[c][k1][b2][k3].Im;
	    Vel[c][b1][b2][k3].Re += C1 * G[c][b1][b2][k3].Re; 
	    Vel[c][b1][b2][k3].Im += C1 * G[c][b1][b2][k3].Im;
	  }
	}
      }
    
  } else {
    /*-----------------------------------------------------------------------*/
    /* Adams-Bashforth 2 timestep.                                           */
    /*-----------------------------------------------------------------------*/
    for (c=1; c<=3; c++)
      for (k1=1; k1<K; k1++) {
	b1 = N-k1;

	kSqrd    = k1*k1;	/* Axes. */
	I_factor = expf(kSqrd * Neg_nu_dt);
	C1       = 0.5 * Run_info.Delta_T * I_factor;
	C2       = -I_factor * C1;
	C1      *= 3.0;
	
	Vel[c][k1][0][0].Re *= I_factor;
	Vel[c][k1][0][0].Im *= I_factor;
	Vel[c][0][k1][0].Re *= I_factor;
	Vel[c][0][k1][0].Im *= I_factor;
	Vel[c][0][0][k1].Re *= I_factor;
	Vel[c][0][0][k1].Im *= I_factor;
	
	Vel[c][k1][0][0].Re += C1 *     G[c][k1][0][0].Re
	                     + C2 * G_old[c][k1][0][0].Re;
	Vel[c][k1][0][0].Im += C1 *     G[c][k1][0][0].Im
	                     + C2 * G_old[c][k1][0][0].Im;
	Vel[c][0][k1][0].Re += C1 *     G[c][0][k1][0].Re
	                     + C2 * G_old[c][0][k1][0].Re;
	Vel[c][0][k1][0].Im += C1 *     G[c][0][k1][0].Im
	                     + C2 * G_old[c][0][k1][0].Im;
	Vel[c][0][0][k1].Re += C1 *     G[c][0][0][k1].Re
	                     + C2 * G_old[c][0][0][k1].Re;
	Vel[c][0][0][k1].Im += C1 *     G[c][0][0][k1].Im
	                     + C2 * G_old[c][0][0][k1].Im;

	for (k2=1; k2<K && k1+k2<=FOURKon3; k2++) {
	  b2 = N-k2;

	  kSqrd    = k1*k1 + k2*k2; /* Faces. */
	  I_factor = expf(kSqrd * Neg_nu_dt);
	  C1       = 0.5 * Run_info.Delta_T * I_factor;
	  C2       = -I_factor * C1;
	  C1      *= 3.0;
	  
	  Vel[c][0][k1][k2].Re *= I_factor;   /* i = 0 */
	  Vel[c][0][k1][k2].Im *= I_factor;
	  Vel[c][0][b1][k2].Re *= I_factor;
	  Vel[c][0][b1][k2].Im *= I_factor;
	  
	  Vel[c][0][k1][k2].Re += C1 *     G[c][0][k1][k2].Re
	                        + C2 * G_old[c][0][k1][k2].Re;
	  Vel[c][0][k1][k2].Im += C1 *     G[c][0][k1][k2].Im
	                        + C2 * G_old[c][0][k1][k2].Im;
	  Vel[c][0][b1][k2].Re += C1 *     G[c][0][b1][k2].Re
	                        + C2 * G_old[c][0][b1][k2].Re;
	  Vel[c][0][b1][k2].Im += C1 *     G[c][0][b1][k2].Im
	                        + C2 * G_old[c][0][b1][k2].Im;
	  
	  Vel[c][k1][0][k2].Re *= I_factor;  /* j = 0 */
	  Vel[c][k1][0][k2].Im *= I_factor;
	  Vel[c][b1][0][k2].Re *= I_factor;
	  Vel[c][b1][0][k2].Im *= I_factor;
	  
	  Vel[c][k1][0][k2].Re += C1 *     G[c][k1][0][k2].Re
	                        + C2 * G_old[c][k1][0][k2].Re;
	  Vel[c][k1][0][k2].Im += C1 *     G[c][k1][0][k2].Im
	                        + C2 * G_old[c][k1][0][k2].Im;
	  Vel[c][b1][0][k2].Re += C1 *     G[c][b1][0][k2].Re
	                        + C2 * G_old[c][b1][0][k2].Re;
	  Vel[c][b1][0][k2].Im += C1 *     G[c][b1][0][k2].Im
	                        + C2 * G_old[c][b1][0][k2].Im;
	  
	  Vel[c][k1][k2][0].Re *= I_factor;  /* k = 0 */
	  Vel[c][k1][k2][0].Im *= I_factor;
	  Vel[c][b1][k2][0].Re *= I_factor;
	  Vel[c][b1][k2][0].Im *= I_factor;
	  
	  Vel[c][k1][k2][0].Re += C1 *     G[c][k1][k2][0].Re
	                        + C2 * G_old[c][k1][k2][0].Re;
	  Vel[c][k1][k2][0].Im += C1 *     G[c][k1][k2][0].Im
	                        + C2 * G_old[c][k1][k2][0].Im;
	  Vel[c][b1][k2][0].Re += C1 *     G[c][b1][k2][0].Re 
	                        + C2 * G_old[c][b1][k2][0].Re;
	  Vel[c][b1][k2][0].Im += C1 *     G[c][b1][k2][0].Im
	                        + C2 * G_old[c][b1][k2][0].Re;

	  for (k3=1; k3<K && k2+k3<=FOURKon3 && k1+k3<=FOURKon3; k3++) {
	    kSqrd    = k1*k1 + k2*k2 + k3*k3; /* Internal. */
	    I_factor = expf(kSqrd * Neg_nu_dt);
	    C1       = 0.5 * Run_info.Delta_T * I_factor;
	    C2       = -I_factor * C1;
	    C1      *= 3.0;
	    
	    Vel[c][k1][k2][k3].Re *= I_factor;
	    Vel[c][k1][k2][k3].Im *= I_factor;
	    Vel[c][b1][k2][k3].Re *= I_factor;
	    Vel[c][b1][k2][k3].Im *= I_factor;
	    Vel[c][k1][b2][k3].Re *= I_factor;
	    Vel[c][k1][b2][k3].Im *= I_factor;
	    Vel[c][b1][b2][k3].Re *= I_factor;
	    Vel[c][b1][b2][k3].Im *= I_factor;
	    
	    Vel[c][k1][k2][k3].Re += C1 *     G[c][k1][k2][k3].Re
	                           + C2 * G_old[c][k1][k2][k3].Re;
	    Vel[c][k1][k2][k3].Im += C1 *     G[c][k1][k2][k3].Im
	                           + C2 * G_old[c][k1][k2][k3].Im;
	    Vel[c][b1][k2][k3].Re += C1 *     G[c][b1][k2][k3].Re
	                           + C2 * G_old[c][b1][k2][k3].Re;
	    Vel[c][b1][k2][k3].Im += C1 *     G[c][b1][k2][k3].Im
	                           + C2 * G_old[c][b1][k2][k3].Im;
	    Vel[c][k1][b2][k3].Re += C1 *     G[c][k1][b2][k3].Re
	                           + C2 * G_old[c][k1][b2][k3].Re;
	    Vel[c][k1][b2][k3].Im += C1 *     G[c][k1][b2][k3].Im
	                           + C2 * G_old[c][k1][b2][k3].Im;
	    Vel[c][b1][b2][k3].Re += C1 *     G[c][b1][b2][k3].Re
	                           + C2 * G_old[c][b1][b2][k3].Re;
	    Vel[c][b1][b2][k3].Im += C1 *     G[c][b1][b2][k3].Im
	                           + C2 * G_old[c][b1][b2][k3].Im;
	  } 
	} 
      } 
  } 
} 
