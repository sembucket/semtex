/*****************************************************************************
 * xsum:  sum = 0;  sum += x[i];
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>


double dsum (int_t n, const double* x, int_t incx)
{
  register int_t i;
  register double  sum = 0.0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) sum += x[i*incx];
  
  return sum;
}


int_t isum (int_t n, const int_t* x, int_t incx)
{
  register int_t i, sum = 0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) sum += x[i*incx];
  
  return sum;
}


float ssum (int_t n, const float* x, int_t incx)
{
  register int_t i;
  register float   sum = 0.0F;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) sum += x[i*incx];
  
  return sum;
}
