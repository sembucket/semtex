/*****************************************************************************
 * xzero:  x[i] = 0.                                                         *
 *****************************************************************************/

#include <string.h>


void dzero(int n, double *x, int incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (double));

  else {
    register int   i;

    x += (incx < 0) ? (-n + 1)*incx : 0;

    for (i = 0; i < n; i++) {
      *x = 0.0;
      x += incx;
    }
  }
}





void izero(int n, int *x, int incx)
{
  if (incx == 1)
    memset(x, '\0', n * sizeof (int));

  else {
    register int   i;

    x += (incx < 0) ? (-n + 1)*incx : 0;

    for (i=0; i<n; i++) {
      *x = 0;
      x += incx;
    }
  }
}





void szero(int n, float *x, int incx)
{



  if (incx == 1)
    memset(x, '\0', n * sizeof (float));

  else {
    register int   i;

    x += (incx < 0) ? (-n + 1)*incx : 0;

    for (i = 0; i < n; i++) {
      *x = 0.0F;
      x += incx;
    }
  }
}
