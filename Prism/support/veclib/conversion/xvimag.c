/*
 * xvimag: y = IMAG(x)
 */

#include "veclib/veclib.h"

void zvimag(int n, zcomplex *x, int incx, double *y, int incy)
{
  while( n-- ) {
    *y = x->i;
    x += incx;
    y += incy;
  }
  return;
}
