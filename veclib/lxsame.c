/*****************************************************************************
 * lxsame: return 1 if all elements of x any y match, else zero.
 *
 * Floating point versions use absolute tolerances on the allowable difference.
 *****************************************************************************/

#include <math.h>

#define EPSSP   6.0e-7
#define EPSDP   6.0e-14


int lisame (int n, const int* x, int incx, const int* y, int incy)
{ 
  register int i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) {
    if (*x != *y) return 0;
    x += incx;
    y += incy;
  }
  
  return 1;
}


int ldsame (int n, const  double *x, int incx, const double* y, int incy)
{ 
  register int i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) {
    if (fabs (*x - *y) > EPSDP) return 0;
    x += incx;
    y += incy;
  }
  
  return 1;
}


int lssame (int n, const  float *x, int incx, const float* y, int incy)
{ 
  register int i;

  x += (incx < 0) ? (-n + 1) * incx : 0;
  y += (incy < 0) ? (-n + 1) * incy : 0;

  for (i = 0; i < n; i++) {
#ifdef __GNUC__			/* -- No support for single-precision math. */
    if (fabs  (*x - *y) > EPSSP) return 0;
#else
    if (fabsf (*x - *y) > EPSSP) return 0;
#endif
    x += incx;
    y += incy;
  }
  
  return 1;
}
