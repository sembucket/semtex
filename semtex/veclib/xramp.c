/*****************************************************************************
 * xramp:  ramp function:  x[i] = alpha + i*beta.
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void dramp (integer n, double alpha, double beta, double* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha + i * beta;
}


void iramp (integer n, integer alpha, integer beta, integer* x, integer incx)
{
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha + i * beta;
}


void sramp (integer n, float alpha, float beta, float* x, integer incx)
{
  register integer i;
  
  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) x[i*incx] = alpha + i * beta;
}
