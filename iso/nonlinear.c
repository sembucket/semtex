/*****************************************************************************
 * nonlinear.c: all computations to produce nonlinear terms in Navier--Stokes.
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


#define ROTATE(Z)   tempRe = (Z.Re); (Z.Re) = -(Z.Im); (Z.Im) = tempRe;
#define SHIFT(Z)    tempRe = (Z.Re); \
                    (Z.Re) = tempRe*cosA - (Z.Im)*sinA; \
                    (Z.Im) = (Z.Im)*cosA + tempRe*sinA;

                         
static void  convolve (const CF,       const CF,       CF,   CVF,
		       const complex*, const complex*, const int*);



void nonlin (/* input     */  const CVF       U   ,
	     /* output    */  CVF             G   ,
	     /* workspace */  CF              F   ,
                              CVF             WK  ,
	     /* using     */  const complex*  Wtab,
                              const complex*  Stab,
                              const int*      Dim )
/* ------------------------------------------------------------------------- *
 * Compute the nonlinear convolution sums created from the combinations of
 * the three velocity-component Fourier coefficients, to produce the terms
 * which contribute to the divergence-free projections described in e.g.,
 * the RHS of Orszag eqs. (2.10) and (2.12), or of Canuto et
 * al. eq. (7.2.11), also Lesieur's eq. (IV-2-7).  The method is to take
 * velocity components two at a time, compute their convolution (F), and
 * add the appropriate terms into the RHS for each velocity component,
 * stored in G[1], G[2], G[3] respectively.  Results are returned for a
 * time-step.
 *
 * The contributions are calculated only in the octodecahedral (isotropi-
 * cally truncated) space described by Orszag.
 * ------------------------------------------------------------------------- */
{
  int               c;
  register int      k1, b1, k2, b2, k3;
  register real     C1, C2, C3, kSqrd, kkonk2, tempRe;
  register complex  F_hat;
  const    int      N        = Dim[1];
  const    int      K        = Dim[3];
  const    int      FOURKon3 = (4 * K) / 3;

  /* -- Convolve u_hat[1] with itself to make F11 and distribute.
   *    Note that nothing happens on the axes, or on the k1=0 face for F11. */

  convolve (U[1], U[1], F, WK, Wtab, Stab, Dim);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2     = N - k2;
      kSqrd  = k1*k1 + k2*k2;
      kkonk2 = k1*k1 / kSqrd;
      C1 =     k1 * (kkonk2 - 1.0);
      C2 =     k2 *  kkonk2;
      
      F_hat = F[k1][ 0][k2];
      G[1][k1][ 0][k2].Re  = C1 * F_hat.Re;
      G[1][k1][ 0][k2].Im  = C1 * F_hat.Im;
      G[3][k1][ 0][k2].Re  = C2 * F_hat.Re;
      G[3][k1][ 0][k2].Im  = C2 * F_hat.Im;
      
      F_hat = F[k1][k2][ 0];
      G[1][k1][k2][ 0].Re = C1 * F_hat.Re;
      G[1][k1][k2][ 0].Im = C1 * F_hat.Im;
      G[2][k1][k2][ 0].Re = C2 * F_hat.Re;
      G[2][k1][k2][ 0].Im = C2 * F_hat.Im;
      
      C1 = -C1;
      
      F_hat = F[b1][ 0][k2];
      G[1][b1][ 0][k2].Re = C1 * F_hat.Re;
      G[1][b1][ 0][k2].Im = C1 * F_hat.Im;
      G[3][b1][ 0][k2].Re = C2 * F_hat.Re;
      G[3][b1][ 0][k2].Im = C2 * F_hat.Im;
      
      F_hat = F[b1][k2][0];
      G[1][b1][k2][0].Re = C1 * F_hat.Re;
      G[1][b1][k2][0].Im = C1 * F_hat.Im;
      G[2][b1][k2][0].Re = C2 * F_hat.Re;
      G[2][b1][k2][0].Im = C2 * F_hat.Im;

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd  = k1*k1 + k2*k2 + k3*k3;
	kkonk2 = k1*k1 / kSqrd;
	C1     = k1 * (kkonk2 - 1.0);
	C2     = k2 *  kkonk2;
	C3     = k3 *  kkonk2;
	
	F_hat = F[k1][k2][k3];
	G[1][k1][k2][k3].Re = C1 * F_hat.Re;
	G[1][k1][k2][k3].Im = C1 * F_hat.Im;
	G[2][k1][k2][k3].Re = C2 * F_hat.Re;
	G[2][k1][k2][k3].Im = C2 * F_hat.Im;
	G[3][k1][k2][k3].Re = C3 * F_hat.Re;
	G[3][k1][k2][k3].Im = C3 * F_hat.Im;
	
	C1 = -C1;
	
	F_hat = F[b1][k2][k3];
	G[1][b1][k2][k3].Re = C1 * F_hat.Re;
	G[1][b1][k2][k3].Im = C1 * F_hat.Im;
	G[2][b1][k2][k3].Re = C2 * F_hat.Re;
	G[2][b1][k2][k3].Im = C2 * F_hat.Im;
	G[3][b1][k2][k3].Re = C3 * F_hat.Re;
	G[3][b1][k2][k3].Im = C3 * F_hat.Im;
	
	C2 = -C2;
	
	F_hat = F[b1][b2][k3];
	G[1][b1][b2][k3].Re = C1 * F_hat.Re;
	G[1][b1][b2][k3].Im = C1 * F_hat.Im;
	G[2][b1][b2][k3].Re = C2 * F_hat.Re;
	G[2][b1][b2][k3].Im = C2 * F_hat.Im;
	G[3][b1][b2][k3].Re = C3 * F_hat.Re;
	G[3][b1][b2][k3].Im = C3 * F_hat.Im;
	
	C1 = -C1;
	
	F_hat = F[k1][b2][k3];
	G[1][k1][b2][k3].Re  = C1 * F_hat.Re;
	G[1][k1][b2][k3].Im  = C1 * F_hat.Im;
	G[2][k1][b2][k3].Re  = C2 * F_hat.Re;
	G[2][k1][b2][k3].Im  = C2 * F_hat.Im;
	G[3][k1][b2][k3].Re  = C3 * F_hat.Re;
	G[3][k1][b2][k3].Im  = C3 * F_hat.Im;
      }
    }
  }

  /* -- Convolve u_hat[2] with itself to make F22.
   *    Nothing happens on the axes or on k2=0 face. */

  convolve (U[2], U[2], F, WK, Wtab, Stab, Dim);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2     = N - k2;
      kSqrd  = k1*k1 + k2*k2;
      kkonk2 = k1*k1 / kSqrd;
      C1     = k1 * (kkonk2 - 1.0);
      C2     = k2 *  kkonk2;
      
      F_hat = F[0][k1][k2];
      G[2][0][k1][k2].Re = C1 * F_hat.Re;
      G[2][0][k1][k2].Im = C1 * F_hat.Im;
      G[3][0][k1][k2].Re = C2 * F_hat.Re;
      G[3][0][k1][k2].Im = C2 * F_hat.Im;
      
      C1 = -C1;
      
      F_hat = F[0][b1][k2];
      G[2][0][b1][k2].Re = C1 * F_hat.Re;
      G[2][0][b1][k2].Im = C1 * F_hat.Im;
      G[3][0][b1][k2].Re = C2 * F_hat.Re;
      G[3][0][b1][k2].Im = C2 * F_hat.Im;
      
      kkonk2 = k2*k2 / kSqrd;
      C1     = k1 *  kkonk2;
      C2     = k2 * (kkonk2 - 1.0);
      
      F_hat = F[k1][k2][0];
      G[1][k1][k2][0].Re += C1 * F_hat.Re;
      G[1][k1][k2][0].Im += C1 * F_hat.Im;
      G[2][k1][k2][0].Re += C2 * F_hat.Re;
      G[2][k1][k2][0].Im += C2 * F_hat.Im;
      
      F_hat = F[b1][k2][0];
      G[1][b1][k2][0].Re -= C1 * F_hat.Re;
      G[1][b1][k2][0].Im -= C1 * F_hat.Im;
      G[2][b1][k2][0].Re += C2 * F_hat.Re;
      G[2][b1][k2][0].Im += C2 * F_hat.Im;

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd  = k1*k1 + k2*k2 + k3*k3;
	kkonk2 = k2*k2 / kSqrd;
	C1     = k1 *  kkonk2;
	C2     = k2 * (kkonk2 - 1.0);
	C3     = k3 *  kkonk2;
	
	F_hat = F[k1][k2][k3];
	G[1][k1][k2][k3].Re += C1 * F_hat.Re;
	G[1][k1][k2][k3].Im += C1 * F_hat.Im;
	G[2][k1][k2][k3].Re += C2 * F_hat.Re;
	G[2][k1][k2][k3].Im += C2 * F_hat.Im;
	G[3][k1][k2][k3].Re += C3 * F_hat.Re;
	G[3][k1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][k2][k3];
	G[1][b1][k2][k3].Re -= C1 * F_hat.Re;
	G[1][b1][k2][k3].Im -= C1 * F_hat.Im;
	G[2][b1][k2][k3].Re += C2 * F_hat.Re;
	G[2][b1][k2][k3].Im += C2 * F_hat.Im;
	G[3][b1][k2][k3].Re += C3 * F_hat.Re;
	G[3][b1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[k1][b2][k3];
	G[1][k1][b2][k3].Re += C1 * F_hat.Re;
	G[1][k1][b2][k3].Im += C1 * F_hat.Im;
	G[2][k1][b2][k3].Re -= C2 * F_hat.Re;
	G[2][k1][b2][k3].Im -= C2 * F_hat.Im;
	G[3][k1][b2][k3].Re += C3 * F_hat.Re;
	G[3][k1][b2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][b2][k3];
	G[1][b1][b2][k3].Re -= C1 * F_hat.Re;
	G[1][b1][b2][k3].Im -= C1 * F_hat.Im;
	G[2][b1][b2][k3].Re -= C2 * F_hat.Re;
	G[2][b1][b2][k3].Im -= C2 * F_hat.Im;
	G[3][b1][b2][k3].Re += C3 * F_hat.Re;
	G[3][b1][b2][k3].Im += C3 * F_hat.Im;
      }
    }
  }
  
  /* -- Convolve u_hat[3] with itself to make F33 and distribute.
   *    Nothing happens on axes or k3=0 face.                      */

  convolve (U[3], U[3], F, WK, Wtab, Stab, Dim);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2     = N - k2;
      kSqrd  = k1*k1 + k2*k2;
      kkonk2 = k2*k2 / kSqrd;
      C1     = k1 *  kkonk2;
      C2     = k2 * (kkonk2 - 1.0);
      
      F_hat = F[0][k1][k2];
      G[2][0][k1][k2].Re += C1 * F_hat.Re;
      G[2][0][k1][k2].Im += C1 * F_hat.Im;
      G[3][0][k1][k2].Re += C2 * F_hat.Re;
      G[3][0][k1][k2].Im += C2 * F_hat.Im;
      
      F_hat = F[0][b1][k2];
      G[2][0][b1][k2].Re -= C1 * F_hat.Re;
      G[2][0][b1][k2].Im -= C1 * F_hat.Im;
      G[3][0][b1][k2].Re += C2 * F_hat.Re;
      G[3][0][b1][k2].Im += C2 * F_hat.Im;
      
      F_hat = F[k1][0][k2];
      G[1][k1][0][k2].Re += C1 * F_hat.Re;
      G[1][k1][0][k2].Im += C1 * F_hat.Im;
      G[3][k1][0][k2].Re += C2 * F_hat.Re;
      G[3][k1][0][k2].Im += C2 * F_hat.Im;
      
      F_hat = F[b1][0][k2];
      G[1][b1][0][k2].Re -= C1 * F_hat.Re;
      G[1][b1][0][k2].Im -= C1 * F_hat.Im;
      G[3][b1][0][k2].Re += C2 * F_hat.Re;
      G[3][b1][0][k2].Im += C2 * F_hat.Im;

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd  = k1*k1 + k2*k2 + k3*k3;
	kkonk2 = k3*k3 / kSqrd;
	C1     = k1 *  kkonk2;
	C2     = k2 *  kkonk2;
	C3     = k3 * (kkonk2 - 1.0);
	
	F_hat = F[k1][k2][k3];
	G[1][k1][k2][k3].Re += C1 * F_hat.Re;
	G[1][k1][k2][k3].Im += C1 * F_hat.Im;
	G[2][k1][k2][k3].Re += C2 * F_hat.Re;
	G[2][k1][k2][k3].Im += C2 * F_hat.Im;
	G[3][k1][k2][k3].Re += C3 * F_hat.Re;
	G[3][k1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][k2][k3];
	G[1][b1][k2][k3].Re -= C1 * F_hat.Re;
	G[1][b1][k2][k3].Im -= C1 * F_hat.Im;
	G[2][b1][k2][k3].Re += C2 * F_hat.Re;
	G[2][b1][k2][k3].Im += C2 * F_hat.Im;
	G[3][b1][k2][k3].Re += C3 * F_hat.Re;
	G[3][b1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[k1][b2][k3];
	G[1][k1][b2][k3].Re += C1 * F_hat.Re;
	G[1][k1][b2][k3].Im += C1 * F_hat.Im;
	G[2][k1][b2][k3].Re -= C2 * F_hat.Re;
	G[2][k1][b2][k3].Im -= C2 * F_hat.Im;
	G[3][k1][b2][k3].Re += C3 * F_hat.Re;
	G[3][k1][b2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][b2][k3];
	G[1][b1][b2][k3].Re -= C1 * F_hat.Re;
	G[1][b1][b2][k3].Im -= C1 * F_hat.Im;
	G[2][b1][b2][k3].Re -= C2 * F_hat.Re;
	G[2][b1][b2][k3].Im -= C2 * F_hat.Im;
	G[3][b1][b2][k3].Re += C3 * F_hat.Re;
	G[3][b1][b2][k3].Im += C3 * F_hat.Im;
      }
    }
  }
  
  /* --Convolve u_hat[1] with u_hat[2] to make F12 (and, by symmetry, F21). */

  convolve (U[1], U[2], F, WK, Wtab, Stab, Dim);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    C1 = -k1;
    
    G[1][0][k1][0].Re = C1 * F[0][k1][0].Re;
    G[1][0][k1][0].Im = C1 * F[0][k1][0].Im;
    G[2][k1][0][0].Re = C1 * F[k1][0][0].Re;
    G[2][k1][0][0].Im = C1 * F[k1][0][0].Im;

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2    =  N - k2;
      kSqrd =  k1*k1 + k2*k2;
      C1    = -k1;
      C2    =  k2 * (2.0*k1*k1/kSqrd - 1.0);
      C3    =  k1 * (2.0*k2*k2/kSqrd - 1.0);
      
      F_hat = F[0][k1][k2];
      G[1][0][k1][k2].Re = C1 * F_hat.Re;
      G[1][0][k1][k2].Im = C1 * F_hat.Im;
      
      F_hat = F[k1][0][k2];
      G[2][k1][0][k2].Re = C1 * F_hat.Re;
      G[2][k1][0][k2].Im = C1 * F_hat.Im;
      
      C1 = k1;
      
      F_hat = F[0][b1][k2];
      G[1][0][b1][k2].Re = C1 * F_hat.Re;
      G[1][0][b1][k2].Im = C1 * F_hat.Im;
      
      F_hat = F[b1][0][k2];
      G[2][b1][0][k2].Re = C1 * F_hat.Re;
      G[2][b1][0][k2].Im = C1 * F_hat.Im;
      
      F_hat = F[k1][k2][0];
      G[1][k1][k2][0].Re += C2 * F_hat.Re;
      G[1][k1][k2][0].Im += C2 * F_hat.Im;
      G[2][k1][k2][0].Re += C3 * F_hat.Re;
      G[2][k1][k2][0].Im += C3 * F_hat.Im;
      
      F_hat = F[b1][k2][0];
      G[1][b1][k2][0].Re += C2 * F_hat.Re;
      G[1][b1][k2][0].Im += C2 * F_hat.Im;
      G[2][b1][k2][0].Re -= C3 * F_hat.Re;
      G[2][b1][k2][0].Im -= C3 * F_hat.Im;

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd = k1*k1 + k2*k2 + k3*k3;
	C1    = k2 * (2.0*k1*k1/kSqrd - 1.0);
	C2    = k1 * (2.0*k2*k2/kSqrd - 1.0);
	C3    = 2.0*k1*k2*k3/kSqrd;
	
	F_hat = F[k1][k2][k3];
	G[1][k1][k2][k3].Re += C1 * F_hat.Re;
	G[1][k1][k2][k3].Im += C1 * F_hat.Im;
	G[2][k1][k2][k3].Re += C2 * F_hat.Re;
	G[2][k1][k2][k3].Im += C2 * F_hat.Im;
	G[3][k1][k2][k3].Re += C3 * F_hat.Re;
	G[3][k1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][k2][k3];
	G[1][b1][k2][k3].Re += C1 * F_hat.Re;
	G[1][b1][k2][k3].Im += C1 * F_hat.Im;
	G[2][b1][k2][k3].Re -= C2 * F_hat.Re;
	G[2][b1][k2][k3].Im -= C2 * F_hat.Im;
	G[3][b1][k2][k3].Re -= C3 * F_hat.Re;
	G[3][b1][k2][k3].Im -= C3 * F_hat.Im;
	
	F_hat = F[k1][b2][k3];
	G[1][k1][b2][k3].Re -= C1 * F_hat.Re;
	G[1][k1][b2][k3].Im -= C1 * F_hat.Im;
	G[2][k1][b2][k3].Re += C2 * F_hat.Re;
	G[2][k1][b2][k3].Im += C2 * F_hat.Im;
	G[3][k1][b2][k3].Re -= C3 * F_hat.Re;
	G[3][k1][b2][k3].Im -= C3 * F_hat.Im;
	
	F_hat = F[b1][b2][k3];
	G[1][b1][b2][k3].Re -= C1 * F_hat.Re;
	G[1][b1][b2][k3].Im -= C1 * F_hat.Im;
	G[2][b1][b2][k3].Re -= C2 * F_hat.Re;
	G[2][b1][b2][k3].Im -= C2 * F_hat.Im;
	G[3][b1][b2][k3].Re += C3 * F_hat.Re;
	G[3][b1][b2][k3].Im += C3 * F_hat.Im;
      }
    }
  }
  
  /* -- Convolve u_hat[2] with u_hat[3] to make F23 (and, by symmetry, F32). */

  convolve (U[2], U[3], F, WK, Wtab, Stab, Dim);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    C1 = -k1;
    
    G[2][0][0][k1].Re = C1 * F[0][0][k1].Re;
    G[2][0][0][k1].Im = C1 * F[0][0][k1].Im;
    G[3][0][k1][0].Re = C1 * F[0][k1][0].Re;
    G[3][0][k1][0].Im = C1 * F[0][k1][0].Im;

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2    =  N - k2;
      kSqrd =  k1*k1 + k2*k2;
      C1    =  k2 * (2.0*k1*k1/kSqrd - 1.0);
      C2    =  k1 * (2.0*k2*k2/kSqrd - 1.0);
      C3    = -k2;
      
      F_hat = F[0][k1][k2];
      G[2][0][k1][k2].Re += C1 * F_hat.Re;
      G[2][0][k1][k2].Im += C1 * F_hat.Im;
      G[3][0][k1][k2].Re += C2 * F_hat.Re;
      G[3][0][k1][k2].Im += C2 * F_hat.Im;
      
      F_hat = F[0][b1][k2];
      G[2][0][b1][k2].Re += C1 * F_hat.Re;
      G[2][0][b1][k2].Im += C1 * F_hat.Im;
      G[3][0][b1][k2].Re -= C2 * F_hat.Re;
      G[3][0][b1][k2].Im -= C2 * F_hat.Im;
      
      F_hat = F[k1][0][k2];
      G[2][k1][0][k2].Re += C3 * F_hat.Re;
      G[2][k1][0][k2].Im += C3 * F_hat.Im;
      
      F_hat = F[b1][0][k2];
      G[2][b1][0][k2].Re += C3 * F_hat.Re;
      G[2][b1][0][k2].Im += C3 * F_hat.Im;
      
      F_hat = F[k1][k2][0];
      G[3][k1][k2][0].Re = C3 * F_hat.Re;
      G[3][k1][k2][0].Im = C3 * F_hat.Im;
      
      F_hat = F[b1][k2][0];
      G[3][b1][k2][0].Re = C3 * F_hat.Re;
      G[3][b1][k2][0].Im = C3 * F_hat.Im;

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd = k1*k1 + k2*k2 + k3*k3;
	C1    = 2.0*k1*k2*k3/kSqrd;
	C2    = k3 * (2.0*k2*k2/kSqrd - 1.0);
	C3    = k2 * (2.0*k3*k3/kSqrd - 1.0);
	
	F_hat = F[k1][k2][k3];
	G[1][k1][k2][k3].Re += C1 * F_hat.Re;
	G[1][k1][k2][k3].Im += C1 * F_hat.Im;
	G[2][k1][k2][k3].Re += C2 * F_hat.Re;
	G[2][k1][k2][k3].Im += C2 * F_hat.Im;
	G[3][k1][k2][k3].Re += C3 * F_hat.Re;
	G[3][k1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][k2][k3];
	G[1][b1][k2][k3].Re -= C1 * F_hat.Re;
	G[1][b1][k2][k3].Im -= C1 * F_hat.Im;
	G[2][b1][k2][k3].Re += C2 * F_hat.Re;
	G[2][b1][k2][k3].Im += C2 * F_hat.Im;
	G[3][b1][k2][k3].Re += C3 * F_hat.Re;
	G[3][b1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[k1][b2][k3];
	G[1][k1][b2][k3].Re -= C1 * F_hat.Re;
	G[1][k1][b2][k3].Im -= C1 * F_hat.Im;
	G[2][k1][b2][k3].Re += C2 * F_hat.Re;
	G[2][k1][b2][k3].Im += C2 * F_hat.Im;
	G[3][k1][b2][k3].Re -= C3 * F_hat.Re;
	G[3][k1][b2][k3].Im -= C3 * F_hat.Im;
	
	F_hat = F[b1][b2][k3];
	G[1][b1][b2][k3].Re += C1 * F_hat.Re;
	G[1][b1][b2][k3].Im += C1 * F_hat.Im;
	G[2][b1][b2][k3].Re += C2 * F_hat.Re;
	G[2][b1][b2][k3].Im += C2 * F_hat.Im;
	G[3][b1][b2][k3].Re -= C3 * F_hat.Re;
	G[3][b1][b2][k3].Im -= C3 * F_hat.Im;
      }
    }
  }
  
  /* -- Convolve u_hat[1] with u_hat[3] to make F13 (and, by symmetry, F31). */

  convolve (U[1], U[3], F, WK, Wtab, Stab, Dim);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    C1 = -k1;
    
    G[1][0][0][k1].Re = C1 * F[0][0][k1].Re;
    G[1][0][0][k1].Im = C1 * F[0][0][k1].Im;
    G[3][k1][0][0].Re = C1 * F[k1][0][0].Re;
    G[3][k1][0][0].Im = C1 * F[k1][0][0].Im;

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2    = N - k2;
      kSqrd = k1*k1 + k2*k2;
      C1    = k2;
      C2    = k2 * (2.0*k1*k1/kSqrd - 1.0);
      C3    = k1 * (2.0*k2*k2/kSqrd - 1.0);
      
      F_hat = F[0][k1][k2];
      G[1][0][k1][k2].Re -= C1 * F_hat.Re;
      G[1][0][k1][k2].Im -= C1 * F_hat.Im;
      
      F_hat = F[0][b1][k2];
      G[1][0][b1][k2].Re -= C1 * F_hat.Re;
      G[1][0][b1][k2].Im -= C1 * F_hat.Im;
      
      F_hat = F[k1][0][k2];
      G[1][k1][0][k2].Re += C2 * F_hat.Re;
      G[1][k1][0][k2].Im += C2 * F_hat.Im;
      G[3][k1][0][k2].Re += C3 * F_hat.Re;
      G[3][k1][0][k2].Im += C3 * F_hat.Im;
      
      F_hat = F[b1][0][k2];
      G[1][b1][0][k2].Re += C2 * F_hat.Re;
      G[1][b1][0][k2].Im += C2 * F_hat.Im;
      G[3][b1][0][k2].Re -= C3 * F_hat.Re;
      G[3][b1][0][k2].Im -= C3 * F_hat.Im;
      
      C1 = k1;
      
      F_hat = F[k1][k2][0];
      G[3][k1][k2][0].Re -= C1 * F_hat.Re;
      G[3][k1][k2][0].Im -= C1 * F_hat.Im;
      
      F_hat = F[b1][k2][0];
      G[3][b1][k2][0].Re += C1 * F_hat.Re;
      G[3][b1][k2][0].Im += C1 * F_hat.Im;

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	kSqrd = k1*k1 + k2*k2 + k3*k3;
	C1    = k3 * (2.0*k1*k1/kSqrd - 1.0);
	C2    = 2.0*k1*k2*k3/kSqrd;
	C3    = k1 * (2.0*k3*k3/kSqrd - 1.0);
	
	F_hat = F[k1][k2][k3];
	G[1][k1][k2][k3].Re += C1 * F_hat.Re;
	G[1][k1][k2][k3].Im += C1 * F_hat.Im;
	G[2][k1][k2][k3].Re += C2 * F_hat.Re;
	G[2][k1][k2][k3].Im += C2 * F_hat.Im;
	G[3][k1][k2][k3].Re += C3 * F_hat.Re;
	G[3][k1][k2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][k2][k3];
	G[1][b1][k2][k3].Re += C1 * F_hat.Re;
	G[1][b1][k2][k3].Im += C1 * F_hat.Im;
	G[2][b1][k2][k3].Re -= C2 * F_hat.Re;
	G[2][b1][k2][k3].Im -= C2 * F_hat.Im;
	G[3][b1][k2][k3].Re -= C3 * F_hat.Re;
	G[3][b1][k2][k3].Im -= C3 * F_hat.Im;
	
	F_hat = F[k1][b2][k3];
	G[1][k1][b2][k3].Re += C1 * F_hat.Re;
	G[1][k1][b2][k3].Im += C1 * F_hat.Im;
	G[2][k1][b2][k3].Re -= C2 * F_hat.Re;
	G[2][k1][b2][k3].Im -= C2 * F_hat.Im;
	G[3][k1][b2][k3].Re += C3 * F_hat.Re;
	G[3][k1][b2][k3].Im += C3 * F_hat.Im;
	
	F_hat = F[b1][b2][k3];
	G[1][b1][b2][k3].Re += C1 * F_hat.Re;
	G[1][b1][b2][k3].Im += C1 * F_hat.Im;
	G[2][b1][b2][k3].Re += C2 * F_hat.Re;
	G[2][b1][b2][k3].Im += C2 * F_hat.Im;
	G[3][b1][b2][k3].Re -= C3 * F_hat.Re;
	G[3][b1][b2][k3].Im -= C3 * F_hat.Im;
      }
    }
  }
  
  /* -- Finally, we multiply the whole lot by (0, i), which is a PI/2
   *    rotation in the complex plane.
   *    The traverse amounts to a concise account of non-Nyquist data. */

  for (c = 1; c <= 3; c++)
    for (k1 = 1; k1 < K; k1++) {
      b1 = N - k1;
      ROTATE (G[c][k1][ 0][ 0]);
      ROTATE (G[c][ 0][k1][ 0]);
      ROTATE (G[c][ 0][ 0][k1]);

      for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
	b2 = N - k2;
	ROTATE (G[c][ 0][k1][k2]);  ROTATE (G[c][ 0][b1][k2]);
	ROTATE (G[c][k1][ 0][k2]);  ROTATE (G[c][b1][ 0][k2]);
	ROTATE (G[c][k1][k2][ 0]);  ROTATE (G[c][b1][k2][ 0]);

	for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	  ROTATE (G[c][k1][k2][k3]);
	  ROTATE (G[c][b1][k2][k3]);
	  ROTATE (G[c][k1][b2][k3]);
	  ROTATE (G[c][b1][b2][k3]);
	}
      }
    }
}


static void  convolve (/* input     */ const CF        u   ,
	 	                       const CF        v   ,
		       /* return    */ CF              w   ,
		       /* workspace */ CVF             WK  ,
		       /* using     */ const complex*  Wtab,
		                       const complex*  Stab,
                                       const int*      Dim )
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
 * By the time the data are passed here, the truncation will already have   
 * been carried out on u & v, and we try to restrict as many operations as  
 * possible to the truncated region to save computation.                    
 * ------------------------------------------------------------------------- */
{
  register int    k1, b1, k2, b2, k3;
  register real   cosA, sinA, tempRe, NormFact;
  register real  *UHead, *VHead, *WHead, *uHead, *vHead;
  CF              U,      V,      W;

  const int N        = Dim[1];
  const int K        = Dim[3];
  const int FOURKon3 = (4 * K) / 3;
  const int Npts     = 2 * Dim[1] * Dim[2] * Dim[3];

  U = WK[1];
  V = WK[2];
  W = WK[3];

  UHead = &U[0][0][0].Re; 
  VHead = &V[0][0][0].Re; 
  WHead = &W[0][0][0].Re; 

  uHead = &u[0][0][0].Re;
  vHead = &v[0][0][0].Re;

  /* -- Parts without phase shifts: make U_hat & V_hat [eq. (6.2)]. */

  for (k1 = 0; k1 < Npts; k1++) {
    UHead[k1] = uHead[k1];
    VHead[k1] = vHead[k1];
  }
  rc3DFT (U, Dim, Wtab, INVERSE);
  rc3DFT (V, Dim, Wtab, INVERSE);
  
  /* w_hat, from eq. (6.4), normalization neglected.  Only shift into w the
   * retained Fourier modes of w_hat.  This controls roundoff in high modes.
   * Note that FOURKon3 = (4*K)/3 is always rounded down so test <=, not <.  */

  for (k1 = 0; k1 < Npts; k1++)
    WHead[k1] = UHead[k1] * VHead[k1];
  rc3DFT (W, Dim, Wtab, FORWARD);
  
  w[0][0][0].Re = W[0][0][0].Re; /* Mean value. */
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    w[k1][0][0] = W[k1][0][0];            /* Axes. */
    w[0][k1][0] = W[0][k1][0];
    w[0][0][k1] = W[0][0][k1];
    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2 = N - k2;
      w[0][k1][k2] = W[0][k1][k2];        /* Faces. */
      w[0][b1][k2] = W[0][b1][k2];
      w[k1][0][k2] = W[k1][0][k2];
      w[b1][0][k2] = W[b1][0][k2];
      w[k1][k2][0] = W[k1][k2][0];
      w[b1][k2][0] = W[b1][k2][0];
      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	w[k1][k2][k3] = W[k1][k2][k3];    /* Internal. */
	w[b1][k2][k3] = W[b1][k2][k3];
	w[k1][b2][k3] = W[k1][b2][k3];
	w[b1][b2][k3] = W[b1][b2][k3];
      }
    }
  }
  
  /* Done with the unshifted computations.  Now we do the data on shifted
   * grid starting with eq (6.3).  U=U_tilde, V=V_tilde.                  */

  for (k1 = 0; k1 < Npts; k1++) {
    UHead[k1] = uHead[k1];
    VHead[k1] = vHead[k1];
  }
  
  for (k1 = 1; k1 < K; k1++) {
    b1   = N - k1;
    cosA = Stab[k1].Re;             /* Axes. */
    sinA = Stab[k1].Im;
    SHIFT (U[k1][0][0]);   SHIFT (V[k1][0][0]);
    SHIFT (U[0][k1][0]);   SHIFT (V[0][k1][0]);
    SHIFT (U[0][0][k1]);   SHIFT (V[0][0][k1]);

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2   = N - k2;
      cosA = Stab[k1+k2].Re;        /* Faces. */
      sinA = Stab[k1+k2].Im;
      SHIFT (U[0][k1][k2]);
      SHIFT (V[0][k1][k2]);
      SHIFT (U[k1][0][k2]);
      SHIFT (V[k1][0][k2]);
      SHIFT (U[k1][k2][0]);
      SHIFT (V[k1][k2][0]);
      
      cosA = Stab[k2-k1].Re;
      sinA = Stab[k2-k1].Im;
      SHIFT (U[0][b1][k2]);
      SHIFT (V[0][b1][k2]);
      SHIFT (U[b1][0][k2]);
      SHIFT (V[b1][0][k2]);
      SHIFT (U[b1][k2][0]);
      SHIFT (V[b1][k2][0]);

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	cosA = Stab[k1+k2+k3].Re;     /* Internal. */
	sinA = Stab[k1+k2+k3].Im;
	SHIFT (U[k1][k2][k3]);  
	SHIFT (V[k1][k2][k3]);
	cosA = Stab[k3+k1-k2].Re;
	sinA = Stab[k3+k1-k2].Im;
	SHIFT (U[k1][b2][k3]);
	SHIFT (V[k1][b2][k3]);
	cosA = Stab[k3+k2-k1].Re;
	sinA = Stab[k3+k2-k1].Im;
	SHIFT (U[b1][k2][k3]);  
	SHIFT (V[b1][k2][k3]);
	cosA = Stab[k3-k1-k2].Re;
	sinA = Stab[k3-k1-k2].Im;
	SHIFT (U[b1][b2][k3]);  
	SHIFT (V[b1][b2][k3]);
      }
    }
  }
  
  /* -- Inverse FFT in eq. (6.3). */

  rc3DFT (U, Dim, Wtab, INVERSE);
  rc3DFT (V, Dim, Wtab, INVERSE);
  
  /* -- Product in space of U_tilde*V_tilde. */

  for (k1 = 0; k1 < Npts; k1++)
    WHead[k1] = UHead[k1] * VHead[k1];
  
  /* -- Eq. (6.5), excepting phase shifts and normalization. */

  rc3DFT (W, Dim, Wtab, FORWARD);
  
  /* Now we phase-shift w_tilde.  Note that we can do this after forward
   * transformation since exp[-iPI/K*k.(1/2, 1/2, 1/2)] is independent of j
   * in eq. (6.5).  Restrict phase shifts to untruncated region.            */

  for (k1 = 1; k1 < K; k1++) {
    b1   = N - k1;
    cosA = Stab[-k1].Re;            /* Axes. */
    sinA = Stab[-k1].Im;
    SHIFT (W[k1][0][0]);
    SHIFT (W[0][k1][0]);
    SHIFT (W[0][0][k1]);

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2   = N - k2;
      cosA = Stab[-k1-k2].Re;       /* Faces. */
      sinA = Stab[-k1-k2].Im;
      SHIFT (W[0][k1][k2]);
      SHIFT (W[k1][0][k2]);
      SHIFT (W[k1][k2][0]);
      
      cosA = Stab[k1-k2].Re;
      sinA = Stab[k1-k2].Im;
      SHIFT (W[0][b1][k2]);
      SHIFT (W[b1][0][k2]);
      SHIFT (W[b1][k2][0]);

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	cosA = Stab[-k1-k2-k3].Re;    /* Internal */
	sinA = Stab[-k1-k2-k3].Im;
	SHIFT (W[k1][k2][k3]);  
	cosA = Stab[k2-k1-k3].Re;
	sinA = Stab[k2-k1-k3].Im;
	SHIFT (W[k1][b2][k3]);
	cosA = Stab[k1-k2-k3].Re;
	sinA = Stab[k1-k2-k3].Im;
	SHIFT (W[b1][k2][k3]);  
	cosA = Stab[k1+k2-k3].Re;
	sinA = Stab[k1+k2-k3].Im;
	SHIFT (W[b1][b2][k3]);  
      }
    }
  }
  
  /* Eq (6.6), on unfiltered region, including normalization.  The normaliz-
   * ing factor includes the factor of 1/2 from eq (6.6), but also two
   * factors of 2 from the real-data inverse FFTs in the k-direction.        */

  NormFact = 1.0 / (N*N*K);
  
  w[0][0][0].Re = NormFact * (w[0][0][0].Re + W[0][0][0].Re);
  
  for (k1 = 1; k1 < K; k1++) {
    b1 = N - k1;
    w[k1][0][0].Re = NormFact * (w[k1][0][0].Re + W[k1][0][0].Re);
    w[k1][0][0].Im = NormFact * (w[k1][0][0].Im + W[k1][0][0].Im);
    w[0][k1][0].Re = NormFact * (w[0][k1][0].Re + W[0][k1][0].Re);
    w[0][k1][0].Im = NormFact * (w[0][k1][0].Im + W[0][k1][0].Im);
    w[0][0][k1].Re = NormFact * (w[0][0][k1].Re + W[0][0][k1].Re);
    w[0][0][k1].Im = NormFact * (w[0][0][k1].Im + W[0][0][k1].Im);

    for (k2 = 1; k2 < K && k1+k2 <= FOURKon3; k2++) {
      b2 = N - k2;
      w[0][k1][k2].Re = NormFact * (w[0][k1][k2].Re + W[0][k1][k2].Re);
      w[0][k1][k2].Im = NormFact * (w[0][k1][k2].Im + W[0][k1][k2].Im);
      w[0][b1][k2].Re = NormFact * (w[0][b1][k2].Re + W[0][b1][k2].Re);
      w[0][b1][k2].Im = NormFact * (w[0][b1][k2].Im + W[0][b1][k2].Im);
      w[k1][0][k2].Re = NormFact * (w[k1][0][k2].Re + W[k1][0][k2].Re);
      w[k1][0][k2].Im = NormFact * (w[k1][0][k2].Im + W[k1][0][k2].Im);
      w[b1][0][k2].Re = NormFact * (w[b1][0][k2].Re + W[b1][0][k2].Re);
      w[b1][0][k2].Im = NormFact * (w[b1][0][k2].Im + W[b1][0][k2].Im);
      w[k1][k2][0].Re = NormFact * (w[k1][k2][0].Re + W[k1][k2][0].Re);
      w[k1][k2][0].Im = NormFact * (w[k1][k2][0].Im + W[k1][k2][0].Im);
      w[b1][k2][0].Re = NormFact * (w[b1][k2][0].Re + W[b1][k2][0].Re);
      w[b1][k2][0].Im = NormFact * (w[b1][k2][0].Im + W[b1][k2][0].Im);

      for (k3 = 1; k3 < K && k2+k3 <= FOURKon3 && k1+k3 <= FOURKon3; k3++) {
	w[k1][k2][k3].Re = NormFact * (w[k1][k2][k3].Re + W[k1][k2][k3].Re);
	w[k1][k2][k3].Im = NormFact * (w[k1][k2][k3].Im + W[k1][k2][k3].Im);
	w[b1][k2][k3].Re = NormFact * (w[b1][k2][k3].Re + W[b1][k2][k3].Re);
	w[b1][k2][k3].Im = NormFact * (w[b1][k2][k3].Im + W[b1][k2][k3].Im);
	w[k1][b2][k3].Re = NormFact * (w[k1][b2][k3].Re + W[k1][b2][k3].Re);
	w[k1][b2][k3].Im = NormFact * (w[k1][b2][k3].Im + W[k1][b2][k3].Im);
	w[b1][b2][k3].Re = NormFact * (w[b1][b2][k3].Re + W[b1][b2][k3].Re);
	w[b1][b2][k3].Im = NormFact * (w[b1][b2][k3].Im + W[b1][b2][k3].Im);
      }
    }
  }
}
