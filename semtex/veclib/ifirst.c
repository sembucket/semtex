/*****************************************************************************
 * ifirst:  index of first non-zero value in x.
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


integer ifirst (integer n, const integer* x, integer incx)
{ 
  register integer i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) if (x[i*incx]) return i;
  
  return -1;
}
