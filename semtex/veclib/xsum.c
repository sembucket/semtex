/*****************************************************************************
 * xsum:  sum = 0;  sum += x[i];
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>


double dsum (integer n, const double* x, integer incx)
{
  register integer i;
  register double  sum = 0.0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) sum += x[i*incx];
  
  return sum;
}


integer isum (integer n, const integer* x, integer incx)
{
  register integer i, sum = 0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) sum += x[i*incx];
  
  return sum;
}


float ssum (integer n, const float* x, integer incx)
{
  register integer i;
  register float   sum = 0.0F;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) sum += x[i*incx];
  
  return sum;
}
