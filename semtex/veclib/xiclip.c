/*****************************************************************************
 * xiclip: inverted clip to interval [alpha,beta]:
 *   if x[i] < (alpha+beta)/2 y[i] = MIN(x[i],alpha) else y[i]=MAX(x[i],beta)
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>

#if defined(__uxp__)
#pragma global novrec
#pragma global noalias
#endif

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

void diclip (integer n, const double alpha, const double beta,
	     const double* x, integer incx,
	           double* y, integer incy)
{
  register integer i;
  register double  xtmp;
  const double     mval = 0.5*(alpha + beta);

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    xtmp = x[i*incx];
    y[i*incy] = (xtmp < mval) ? MIN(xtmp, alpha) : MAX(xtmp, beta);
  }
}


void iiclip (integer n, const integer alpha, const integer beta,
	     const integer* x, integer incx,
	           integer* y, integer incy)
{
  register integer i, xtmp;
  const integer    mval = (alpha + beta)/2;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    xtmp = x[i*incx];
    y[i*incy] = (xtmp <= mval) ? MIN(xtmp, alpha) : MAX(xtmp, beta);
  }
}


void siclip (integer n, const float alpha, const float beta,
	     const float* x, integer incx,
	           float* y, integer incy)
{
  register integer i;
  register float   xtmp;
  const float      mval=0.5*(alpha + beta);

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) {
    xtmp = x[i*incx];
    y[i*incy] = (xtmp < mval) ? MIN(xtmp, alpha) : MAX(xtmp, beta);
  }
}

#undef MIN
#undef MAX
