/*****************************************************************************
 * xspline: Cubic spline interpolation.
 *
 * Two routines to perform cubic spline interpolation.  The function
 * "spline" must be called once to compute the spline coefficients, and then
 * "splint" can be called any number of times to evaluate the spline.
 *
 * Input:
 *    n           number of control points
 *    x[0..n-1]   list of nodes
 *    y[0..n-1]   function values at the nodes, y = f(x)
 *    yp1         function derivative f'(x) at the first node
 *    yp2         function derivative f'(x) at the last node
 *
 * If yp1 or ypn is > 1.e30, "spline" uses "natural" spline end conditions
 * (d2y/dx2 = 0) at the respective end.
 *
 * Output:
 *    y2[0..n-1]  the spline coefficients
 *
 * These are routines from Numerical Recipes, modified to use base-0
 * indexed arrays.
 *
 * $Id$
 *****************************************************************************/

#include <stdio.h>
#include "alplib.h"



void dspline(int n, double yp1, double ypn, const double *x, const double *y,
	                                                           double *y2)
{
  double   h = x[1] - x[0], 
          *u = dvector(0, n-2);
  double   p, qn, sig, un, hh;
  register i, k;

  if (yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u [0] = (3.0 / h) * ((y[1]-y[0]) / h - yp1);
  }

  for (i = 1; i < n-1; i++) {
    hh     = 1.0 / (x[i+1] - x[i-1]);
    sig    = h * hh;
    p      = sig * y2[i-1] + 2.0;
    y2[i]  = (sig - 1.0) / p;
    u [i]  = (y[i] - y[i-1]) / h;
    u [i]  = (y[i+1] - y[i]) / (h = x[i+1] - x[i]) - u[i];
    u [i]  = (6.0 * u[i] * hh - sig * u[i-1]) / p;
  }

  if (ypn > 0.99e30)
    qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0 / h) * (ypn - (y[n-1]-y[n-2])/h);
  }

  y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);
  for (k = n-2; k; k--) y2[k] = y2[k] * y2[k+1] + u[k];

  freeDvector(u, 0);
  return;
}





double dsplint(int n, double x, const double *xa, const double *ya,
	                                          const double *y2a)
{
  register k;
  register double h, b, a;
  static   klo = -1, khi = -1;

  /* check the results of the previous search */

  if (klo < 0 || khi > n || xa[khi] < x || xa[klo] > x) {
    klo = 0;
    khi = n-1;
  }
    
  while (khi-klo > 1) {    /* search for the bracketing interval */
    k = (khi+klo) >> 1;
    if (xa[k] > x) 
      khi = k; 
    else 
      klo = k;
  }

  h  = xa[khi]-xa[klo];
  a  = (xa[khi]-x)/h;
  b  = (x-xa[klo])/h;

  return a * ya[klo] + b * ya[khi] + 
    ((a*a*a-a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * (h*h)/6.0;
}			  





void sspline(int n, float yp1, float ypn, const float *x, const float *y,
	                                                        float *y2)
{
  float    h = x[1] - x[0], 
          *u = svector(0, n-2);
  float    p, qn, sig, un, hh;
  register i, k;

  if (yp1 > 0.99e30)
    y2[0] = u[0] = 0.0F;
  else {
    y2[0] = -0.5F;
    u [0] = (3.0F / h) * ((y[1]-y[0]) / h - yp1);
  }

  for (i = 1; i < n-1; i++) {
    hh     = 1. / (x[i+1] - x[i-1]);
    sig    = h * hh;
    p      = sig * y2[i-1] + 2.;
    y2[i]  = (sig - 1.) / p;
    u [i]  = (y[i] - y[i-1]) / h;
    u [i]  = (y[i+1] - y[i]) / (h = x[i+1]-x[i]) - u[i];
    u [i]  = (6.0F * u[i] * hh - sig * u[i-1])/p;
  }

  if (ypn > 0.99e30)
    qn = un = 0.0F;
  else {
    qn = 0.5F;
    un = (3.0F / h)*(ypn - (y[n-1]-y[n-2])/h);
  }

  y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.);
  for (k = n-2; k; k--) y2[k] = y2[k] * y2[k+1] + u[k];

  freeSvector(u, 0);
}





float ssplint(int n, float x, const float *xa, const float *ya,
	                                       const float *y2a)
{
  register k;
  register float h, b, a;
  static   klo = -1, khi = -1;

  /* check the results of the previous search */

  if (klo < 0 || khi > n || xa[khi] < x || xa[klo] > x) {
    klo = 0;
    khi = n-1;
  }
    
  while (khi-klo > 1) {    /* search for the bracketing interval */
    k = (khi+klo) >> 1;
    if (xa[k] > x) 
      khi = k; 
    else 
      klo = k;
  }

  h  = xa[khi]-xa[klo];
  a  = (xa[khi]-x)/h;
  b  = (x-xa[klo])/h;

  return a * ya[klo] + b * ya[khi] + 
    ((a*a*a-a) * y2a[klo] + (b*b*b - b) * y2a[khi]) * (h*h)/6.0F;
}			  






