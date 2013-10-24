/*****************************************************************************
 * lany: return 1 if any x are true: iany = 0; if (x[i]) lany = 1.
 *
 * $Id$
 *****************************************************************************/

#include <cfemdef.h>


int_t lany (int_t n, const int_t* x, int_t incx)
{ 
  register int_t i;

  x += (incx<0) ? (-n+1)*incx : 0;

  for (i = 0; i < n; i++) if (x[i*incx]) return 1;
  
  return 0;
}
