/*****************************************************************************
 * xseq:  y[i] = alpha == x[i].
 *
 * $Id$
 *****************************************************************************/

#include <femdef.h>


void iseq (integer n, integer alpha,
	   const integer* x, integer incx,
	         integer* y, integer incy)
{
  register integer  i;

  x += (incx<0) ? (-n+1)*incx : 0;
  y += (incy<0) ? (-n+1)*incy : 0;

  for (i = 0; i < n; i++) y[i*incy] = (x[i*incx] == alpha) ? 1 : 0;
}
