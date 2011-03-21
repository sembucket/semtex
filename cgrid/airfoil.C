/*****************************************************************************
 * airfoil.C: routines for airfoils given by analytic functions.
 *
 * $Id$
 *****************************************************************************/


#include <math.h>
#include <Utility.h>
#include <Grid.h>


void NACAdsdx (double x, double* s, double* dsdx, void* params)
// ---------------------------------------------------------------------------
// This routine calculates the rate of change of chord length dsdx at
// given (dimensionless) location along zero line for NACA-00't' series
// foils, and is designed to be called by RK evaluation routines.
//
// Input 's' is not needed, but is included so routine conforms to
// conventions used by RK routines.  dsdx is base-1 indexed.
// ---------------------------------------------------------------------------
{
  NACAparam*  N = (NACAparam*) params;
  double      D;

  D =  4.0 * N -> a[4];
  D =  3.0 * N -> a[3] + D * x;
  D =  2.0 * N -> a[2] + D * x;
  D =        N -> a[1] + D * x;
  D =  0.5 * N -> a[0] / sqrt (x) + D;
  D *= N -> t;

  dsdx[1] = sqrt (1.0 + sqr (D));
}


double NACAfoil (double x, void* params)
// ---------------------------------------------------------------------------
// Return y value for NACA00't' x coordinate.
// ---------------------------------------------------------------------------
{
  NACAparam*  N = (NACAparam*) params;
  double      y;

  y = N -> a[4];
  y = N -> a[3] + y * x;
  y = N -> a[2] + y * x;
  y = N -> a[1] + y * x;
  y = N -> t * (N -> a[0] * sqrt (x) + y * x);

  return y;
}
