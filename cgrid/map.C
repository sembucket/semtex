// ***************************************************************************
// Generate dz/dc along branch cut for airfoil grid generation.
// ***************************************************************************

// $Id$

#include <iostream.h>
#include <math.h>
#include <Utility.h>


extern double stretch1 (double, double,  double);
extern double NACA00   (double, double&, double&, double&);



int main ()
// ===========================================================================
//
// ===========================================================================
{
  const int     NC   = 16;
  const int     NO   = 16;
  const int     NY   = 8;
  const double  P    = 0.1;
  const double  Q    = 2.0;
  const double  Lo   = 4.0;	// Outflow length.
  const double  Lc   = 2.0;	// Cross-flow dimension.
  const double  AcLo = 2.0634;	// ArcCosh(Lo).
  const double  T    = 0.5;	// Airfoil thickness;

  int           i;
  double        estar, ratse, r, t, chord;

  double*       cx = dvector (0, NC+NO);
  double*       cy = dvector (0, NC+NO);
  double*       x  = dvector (0, NC+NO);
  double*       y  = dvector (0, NY);

  // -- Generate control grid.

  for (i = 0; i < NC; i++) {
    estar = i * 1.0 / NC;
    ratse = 1.0 - estar;
    r     = stretch1 (P, Q, estar);
    t     = 1.0 - stretch1 (P, Q, ratse);
    cx[i]  = ratse * r + estar * t;
  }

  for (i = 0; i <= NO; i++) {
    estar = i * 1.0 / NO;
    r = 1.0 + cosh(AcLo * estar) * stretch1 (P, Q, estar);
    cx[NC + i] =  r;
  }

  for (i = 0; i <= NY; i++) {
    estar = i * Lc / NY;
    r     = exp (estar) - 1.0;
    cy[i] = r;
  }

  // -- Calculate position along upper side of J = 1 line.

  NACA00 (T, chord, x[0], y[0]);

  cout << chord << endl;
/*  
  for (i = 0; i <= NC; i++)
    NACA00 (T, cx[i] / chord, x[i], y[i]);
  for (i = 1; i <= NO; i++) {
    x[NC + i] = cx[i];
    y[NC + i] = 0.0;
  }
  
  // -- Print out input for cgrid.f

  cout << 2 * (NC + NO) - 1 << "  " << NY << "  " << 0.10 << "  " << 0 << endl;

  for (i = 0; i <  NO + NC; i++)
    cout << x[NC + NO - i] << "  " << -y[NC + NO - i] << endl;
  for (i = 0; i <= NC + NO; i++)
    cout << x[i]           << "  " <<  y[i]           << endl;
*/  
  return 0;
}

