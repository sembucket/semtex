/*****************************************************************************
 *                       RANDOM NUMBER GENERATION                            *
 *                                                                           *
 * The following set of routines provides several types of random number     * 
 * generation.  The only routines provided as part of VECLIB generate a      *
 * set of numbers distributed uniformly on (0,1), but an extension to        *
 * normally-distributed RVs has been implemented.                            *
 *                                                                           *
 * Routines come from Numerical Recipes.                                     *
 *****************************************************************************/

#include <math.h>
#include <time.h>


static double  UD     (double, double);
static double  GD     (double, double);
static double  ran2   (long *);
static double  gasdev (long *);
static long    iseed  = 0;




 
double  dranu(void)
/* ========================================================================= *
 * Provide a single random number UD on (0, 1).                              * 
 * ========================================================================= */
{
  return UD(0.0, 1.0);
}





float  sranu(void)
/* ========================================================================= *
 * Provide a single random number UD on (0, 1).                              * 
 * ========================================================================= */
{
  return (float) UD(0.0, 1.0);
}





double  dnormal(double mean, double sdev)
/* ========================================================================= *
 * Provide a single random number, Normal(mean, sdev).                       * 
 * ========================================================================= */
{
  return GD(mean, sdev);
}





float  snormal(float mean, float sdev)
/* ========================================================================= *
 * Provide a single random number, Normal(mean, sdev).                       * 
 * ========================================================================= */
{
  return (float) GD((double) mean, (double) sdev);
}





void  dvrandom(int n, double *x, int incx)
/* ========================================================================= *
 * Randomize vector x, UD on (0, 1).                                         *
 * ========================================================================= */
{
  register int i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = UD(0.0, 1.0);
    x += incx;
  }
}





void  svrandom(int n, float *x, int incx)
/* ========================================================================= *
 * Randomize vector x, UD on (0, 1).                                         *
 * ========================================================================= */
{
  register int i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = (float) UD(0.0, 1.0);
    x += incx;
  }
}





void  dvnormal(int n, double mean, double sdev, double *x, int incx)
/* ========================================================================= *
 * Randomize vector x, Normal(mean, sdev).                                   *
 * ========================================================================= */
{
  register int i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = GD(mean, sdev);
    x += incx;
  }
}





void  svnormal(int n, float mean, float sdev, float *x, int incx)
/* ========================================================================= *
 * Randomize vector x, Normal(mean, sdev).                                   *
 * ========================================================================= */
{
  register int i;

  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i=0; i<n; i++) {
    *x = (float) GD((double) mean, (double) sdev);
    x += incx;
  }
}





void  raninit(int flag)
/* ========================================================================= *
 * Initialize random number generator.                                       *
 * ========================================================================= */
{
  iseed = (flag < 0) ? time(NULL) : flag;

  (void) ran2(&iseed);
}





/*****************************************************************************
 * Remaining routines are only accessible internally.                        *
 *****************************************************************************/


static double UD(double low, double high)
/* ========================================================================= *
 * Return a random number UD on (low, high).                                 *
 * ========================================================================= */
{
  return ran2(&iseed) * (high-low) + low;
}





static double GD(double mean, double sdev)
/* ========================================================================= *
 * Return normally distributed random deviate with specified mean and        *
 * standard deviation.                                                       *
 * ========================================================================= */
{
  return gasdev(&iseed) * sdev + mean;
}





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

static double ran2(long *idum)
/* ========================================================================= *
 * Long period (>2x10^18) random number generator of L'Ecuyer with Bays-     *
 * Durham shuffle and added safeguards.  Returns a uniform random deviate    *
 * between 0.0 & 1.0 (exclusive of endpoints).  Call with idum a negative    *
 * integer to initialize; thereafter, do not alter idum between successive   *
 * deviates in a sequence.  RNMX should approximate the largest floating     *
 * value that is less than 1.                                                *
 * ========================================================================= */
{
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

  if ((temp=AM*iy) > RNMX) return RNMX;
  else                     return temp;
}





double gasdev(long *idum)
/* ========================================================================= *
 * Returns a normally distributed deviate with zero mean and unit variance,  *
 * using ran2(idum) as the source of uniform deviates.                       *
 * ========================================================================= */
{
  static int     iset = 0;
  static double  gset;
  double         fac, rsq, v1, v2;


  if (iset == 0) {
    do {
      v1 = 2.0*ran2(idum) - 1.0;
      v2 = 2.0*ran2(idum) - 1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0*log(rsq) / rsq);
    gset = v1 * fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}
