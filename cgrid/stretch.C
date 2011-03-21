/****************************************************************************
 * Stretching functions.
 ****************************************************************************/

// $Id$

#include <math.h>


double stretch1 (double P, double Q, double estar)
// ---------------------------------------------------------------------------
// Stretching function, eq. 13.44 from Fletcher.
//
// Dimensionless parameter estar (between zero & one) gets mapped
// into [0, 1] by function.
// ---------------------------------------------------------------------------
{
  double s;

  s = Q * (1.0 - estar);
  s = 1.0 - tanh (s) / tanh (Q);
  s = P * estar + (1.0 - P) * s;

  return s;
}
