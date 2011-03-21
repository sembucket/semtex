// ***************************************************************************
// Routines to compute aerofoil coordinates.
// ***************************************************************************

// $Id$

#include <math.h>
#include <Utility.h>


void NACA00 (double t, double& s, double& x, double& y)
// ---------------------------------------------------------------------------
// Adaption of routine from Fletcher to compute NACA-00't' aerofoil.
//
// The first time the routine is called, an internal table is created
// that describes how the surface coordinate (curve length, s) varies
// with 0 <= x <= 1 along the upper surface of the airfoil.  The
// total curve length is returned in s.
//
// Subsequent calls interpolate in this table to find x, given the surface
// coordinate s.  Then y is calculated from the formula that describes the
// profile, and x & y are returned.
//
// NACA coefficients modified slightly from those given in Mises' book
// in order to get closure at trailing edge.  This is done by reducing
// quartic coefficient slightly.
//
// Input parameter t is dimensionless thickness.
// ---------------------------------------------------------------------------
{
  const  int      TAB_MAX  = 128;
  const  double   a[] = {1.4845, -0.6300, -1.7580, 1.4215, -0.5180};
  static double*  xa;
  static double*  sa;

  if (!xa) {			// -- Make lookup table.
    int     i,   ip;
    double  dum, dx, fl, flp;

    xa = dvector (0, TAB_MAX);
    sa = dvector (0, TAB_MAX);

    sa[0] = 0.0;
    xa[0] = 0.0;
    xa[1] = sqr (a[0] / (1.0 / t - a[1]));
    sa[1] = 0.5 * M_PI * xa[1];
    dum   = 3.0 * a[3] + 4.0 * a[4] * xa[1];
    dum   = 2.0 * a[2] + dum * xa[1];
    dum   = t * (0.5 * a[0] / sqrt (xa[1]) + a[1] + dum * xa[1]);
    flp   = sqrt (1.0 + sqr (dum));

    for (i = 1; i < TAB_MAX; i++) {
      ip     = i + 1;
      xa[ip] = ip / (double) TAB_MAX;
      dx     = xa[ip] - xa[i];
      fl     = flp;
      dum    = 3.0 * a[3] + 4.0 * a[4] * xa[ip];
      dum    = 2.0 * a[2] + dum * xa[ip];
      dum    = t * (0.5 * a[0] / sqrt (xa[ip]) + a[1] + dum * xa[ip]);
      flp    = sqrt (1.0 + sqr (dum));
      sa[ip] = sa[i] + 0.5 * (fl + flp) * dx;
    }
    s = sa[TAB_MAX];
    return;
  }
  
  // -- Use linear interpolation in table.

  if (s > sa[TAB_MAX])
    message ("NACA00", "input s exceeds curve length", ERROR);
    
  int      i, im;
  double   dum;
  
  for (i = 1; i <= TAB_MAX; i++) {
    if (s > sa[i]) continue;
    im  = i - 1;
    x   = xa[im] + (xa[i] - xa[im]) * (s - sa[im]) / (sa[i] - sa[im]);
    x   = (x < EPSSP) ? EPSSP : x;
    dum = a[3] + a[4] * x;
    dum = a[2] + dum  * x;
    dum = a[1] + dum  * x;
    y   = t * (a[0] * sqrt (x) + dum * x);
    return;
  }
}
