/*****************************************************************************
 * icount:  number of non-zero values in x.
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


integer icount (integer n, const integer* x, integer incx)
{
  register integer i, sum = 0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++ ) sum += (x[i*incx]) ? 1 : 0;

  return sum;
}
