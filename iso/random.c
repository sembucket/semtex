/*****************************************************************************
 * random.c: random number generation, randomize a CVF with a uniform 
 * spectral density (uniform PSD at each Fourier component).
 *
 * Copyright (C) 1992-1999 Hugh Blackburn
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"

                  
void randomise (CVF  U   ,
		int* seed,
		CF   work)
/* ------------------------------------------------------------------------- *
 * Create a random spectrum with uniform spectral density in U, but also
 * ensure that it is divergence-free, and only set on octodecahedral space.
 * ------------------------------------------------------------------------- */
{
  const int ntot = N * N * N;
  int       c, i, j, k, bi, bj;
  real*     u = &U[1][0][0][0].Re;
  real*     v = &U[2][0][0][0].Re;
  real*     w = &U[3][0][0][0].Re;

  for (i = 0; i < ntot; i++) u[i] = ran2PI (seed);
  for (i = 0; i < ntot; i++) v[i] = ran2PI (seed);
  for (i = 0; i < ntot; i++) w[i] = ran2PI (seed);

  truncateVF (U);
  projectVF  (U, work);
}


#if 1
    
real ran2PI (int* idum)
/* ------------------------------------------------------------------------- *
 * Generate IUD random variates on (0, 2PI).  This is a doctoring of ran1()
 * from Numerical Recipes in C, 1st ed.
 * ------------------------------------------------------------------------- */
{
#define M1    259200
#define IA1   7141
#define IC1   54773
#define RM1   (1.0/M1)
#define M2    134456
#define IA2   8121
#define IC2   28411
#define RM2   (1.0/M2)
#define M3    243000
#define IA3   4561
#define IC3   51349

  static   long  ix1, ix2, ix3;
  static   float r[98];
  static   int   iff = 0;
  float          temp;
  register int   j;
    
  if (*idum < 0 || iff == 0) {
    iff=1;
    ix1=(IC1-(*idum)) % M1;
    ix1=(IA1*ix1+IC1) % M1;
    ix2=ix1 % M2;
    ix1=(IA1*ix1+IC1) % M1;
    ix3=ix1 % M3;
    for (j=1;j<=97;j++) {
      ix1=(IA1*ix1+IC1) % M1;
      ix2=(IA2*ix2+IC2) % M2;
      r[j]=(ix1+ix2*RM2)*RM1;
    }
    *idum=1;
  }
  ix1=(IA1*ix1+IC1) % M1;
  ix2=(IA2*ix2+IC2) % M2;
  ix3=(IA3*ix3+IC3) % M3;
  j=1 + ((97*ix3)/M3);
  if (j > 97 || j < 1) message ("ran2PI", "This cannot happen.", ERROR);
  temp=r[j];
  r[j]=(ix1+ix2*RM2)*RM1;

  return 2.0 * M_PI * temp;
    
#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3
}

#else

real ran2PI (long* idum)
/* ------------------------------------------------------------------------- *
 * Long period (>2x10^18) random number generator of L'Ecuyer with Bays-
 * Durham shuffle and added safeguards.  Returns a uniform random deviate
 * between 0.0 & 1.0 (exclusive of endpoints).  Call with idum a negative
 * integer to initialize; thereafter, do not alter idum between successive
 * deviates in a sequence.  RNMX should approximate the largest floating
 * value that is less than 1.
 * ------------------------------------------------------------------------- */
{

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-13
#define RNMX  (1.0-EPS)

  int          j;
  long         k;
  static long  idum2=123456789;
  static long  iy=0;
  static long  iv[NTAB];
  double       temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum = 1;
    else              *idum = -(*idum);
    idum2 = (*idum);
    for (j=NTAB+7; j>=0; j--) {
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) -k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if (*idum < 0) *idum += IM1;

  k = idum2 / IQ2;
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;

  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;

  if ((temp=AM*iy) > RNMX) return 2.0 * M_PI * RNMX;
  else                     return 2.0 * M_PI * temp;

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
}

#endif


#if 0

/* Rogallo-style ICs. */

static void doalphabeta   (complex*, complex*, real, int*);
static void setcomponents (CVF, int, int, int, int, int, int,
			   complex, complex, real, real);


void randomise (int  seed,
		CVF  U   )
/* ------------------------------------------------------------------------- *
 * Create a random spectrum with uniform spectral density in U, but also
 * ensure that it is divergence-free, and only set on octodecahedral space.
 *
 * k12 = sqrt(k1*k1 + k2*k2), kk = sqrt(k1*k1 + k2*k2 + k3*k3).  See
 * Rogallo (1981), but note we use a different form on the axes, since
 * we know that the velocity Fourier components must be perpendicular
 * to the axes, and the form derived by Rogallo blows up on the k3
 * axis.
 * ------------------------------------------------------------------------- */
{
  int     c, i, bi, j, bj, k;
  real    k12, kk;
  complex alpha, beta;

  for (i = 1; i < K; i++) {		/* -- Axes. */
    kk = i;
    doalphabeta(&alpha, &beta, kk, &seed);
    U[1][i][0][0].Re =          (U[1][i][0][0].Im = 0.0);
    U[2][i][0][0].Re = alpha.Re; U[2][i][0][0].Im = alpha.Im;
    U[3][i][0][0].Re = beta.Re;  U[3][i][0][0].Im = beta.Im;
    doalphabeta(&alpha, &beta, kk, &seed);
    U[1][0][i][0].Re = alpha.Re; U[2][0][i][0].Im = alpha.Im;
    U[2][0][i][0].Re =          (U[2][0][i][0].Im = 0.0);
    U[3][0][i][0].Re = beta.Re;  U[3][0][i][0].Im = beta.Im;
    doalphabeta(&alpha, &beta, kk, &seed);
    U[1][0][0][i].Re = alpha.Re; U[1][0][0][i].Im = alpha.Im;
    U[2][0][0][i].Re = beta.Re;  U[2][0][0][i].Im = beta.Im;
    U[3][0][0][i].Re =          (U[3][0][0][i].Im = 0.0);
  }

  for (i = 1; i < K; i++) {		/* -- Faces. */
    bi = N-i;
    for (j = 1; j < K; j++) {
      kk = sqrt (i*i + j*j);
      k12 = i;                         /* i = 0 */
      doalphabeta(&alpha, &beta, kk, &seed);
      setcomponents(U,   0,  i, j,   0,  i, j,   alpha, beta, k12, kk);
      doalphabeta(&alpha, &beta, kk, &seed);
      setcomponents(U,   0, bi, j,   0, -i, j,   alpha, beta, k12, kk);
                                         /* j = 0 */
      doalphabeta(&alpha, &beta, kk, &seed);
      setcomponents(U,    i, 0, j,    i, 0, j,   alpha, beta, k12, kk);
      doalphabeta(&alpha, &beta, kk, &seed);
      setcomponents(U  , bi, 0, j,   -i, 0, j,   alpha, beta, k12, kk); 
      k12 = sqrt (i*i + j*j);          /* k = 0 */
      doalphabeta(&alpha, &beta, kk, &seed);
      setcomponents(U,    i, j, 0,    i, j, 0,   alpha, beta, k12, kk);
      doalphabeta(&alpha, &beta, kk, &seed);
      setcomponents(U,   bi, j, 0,   -i, j, 0,   alpha, beta, k12, kk);
    }
  }

  for (i = 1; i < K; i++) {		/* -- Internal. */
    bi = N - i;
    for (j = 1; j < K; j++) {
      bj = N - j;
      k12 = sqrt (i*i + j*j);
      for (k = 1; k < K; k++) {
	kk = sqrt (i*i + j*j + k*k);
	doalphabeta(&alpha, &beta, kk, &seed);
	setcomponents(U,    i,  j, k,     i,  j, k, alpha, beta, k12, kk);
	doalphabeta(&alpha, &beta, kk, &seed);
	setcomponents(U,   bi,  j, k,    -i,  j, k, alpha, beta, k12, kk);
	doalphabeta(&alpha, &beta, kk, &seed);
	setcomponents(U,    i, bj, k,     i, -j, k, alpha, beta, k12, kk);
	doalphabeta(&alpha, &beta, kk, &seed);
	setcomponents(U,   bi, bj, k,    -i, -j, k, alpha, beta, k12, kk);
      }
    }
  }

  truncateVF (U);		/* -- Filter off high frequencies. */
}


static void doalphabeta (complex* alpha,
			 complex* beta ,
			 real     kk   ,
			 int*     seed )
/* ------------------------------------------------------------------------- *
 *
 * ------------------------------------------------------------------------- */
{
  static real theta1, theta2, phi;

  theta1      = ran2PI (seed);
  theta2      = ran2PI (seed);
  phi         = ran2PI (seed);
  alpha -> Re = cos (theta1) * cos (phi) / kk;
  alpha -> Im = sin (theta1) * cos (phi) / kk;
  beta  -> Re = cos (theta2) * sin (phi) / kk;
  beta  -> Im = sin (theta2) * sin (phi) / kk;
}

                
static void setcomponents(/* update     */ CVF     U  ,
			  /* location   */ int     i  ,
		                           int     j  ,
		                           int     k  ,
			  /* wavenumber */ int     k1 ,
		                           int     k2 ,
		                           int     k3 ,
			  /* using      */ complex A  ,
		                           complex B  ,
		                           real    k12,
		                           real    kk )
/* ------------------------------------------------------------------------- *
 * This should (check) make U divergence-free.
 * ------------------------------------------------------------------------- */
{
  static real denom;
  
  denom = k12*kk;

  U[1][i][j][k].Re = (A.Re*kk*k2 + B.Re*k1*k3) / denom;
  U[1][i][j][k].Im = (A.Im*kk*k2 + B.Im*k1*k3) / denom; 

  U[2][i][j][k].Re = (B.Re*k2*k3 - A.Re*kk*k1) / denom;
  U[2][i][j][k].Im = (B.Im*k2*k3 - A.Im*kk*k1) / denom;

  U[3][i][j][k].Re = -B.Re * k12 / kk;
  U[3][i][j][k].Im = -B.Im * k12 / kk;
}

#endif
