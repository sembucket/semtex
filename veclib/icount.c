/*****************************************************************************
 * icount:  number of non-zero values in x.
 *
 * $Id: icount.c,v 9.1 2019/05/30 06:36:13 hmb Exp $
 *****************************************************************************/

#include <cfemdef.h>


int_t icount (int_t n, const int_t* x, int_t incx)
{
  register int_t i, sum = 0;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++ ) sum += (x[i*incx]) ? 1 : 0;

  return sum;
}
