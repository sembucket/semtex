// ***************************************************************************
// cgen: generate input for FORTRAN C-grid generation routines.
//
// This code makes input for a grid wrapped around NACA airfoil of unit chord.
//
// $Id$
// ***************************************************************************

#include <cstdlib>
#include <cmath>

#include <iostream>

using namespace std;

#include <utility.h>
#include <veclib.h>

#include <Grid.h>


const NACAparam  NACA0012 = { 0.12, 1.4845, -0.63, -1.758, 1.4215, -0.518 };
const int        NC       = 32;	 // Number of points along 1/2 body boundary.
const int        NO       = 16;	 // Number of points to outflow boundary.
const int        NY       = 16;   // Number of cross-flow control points.
const double     b4       = 0.12;

static void  controlGrids  (int, int, int, double*, double*);
static void  airfoilSpline (int, const double*, double*, double*);
static void  printHeader   (int, int);
static void  printBoundary (int, int, const double*,
			              const double*, const double*);
static void  printJacobian (int, int, const double*, const double*);


int main ()
// ---------------------------------------------------------------------------
// Generate stretched coordintes for control grid, and airfoil surface
// interpolation.  Then output data for use by cgrid.f.
//
// For now, all choices of mesh size & extent are hard-coded.  Aerofoil
// is a symmetrical NACA section, slightly modified to obtain closure
// at TE.
// ---------------------------------------------------------------------------
{
  double*  cx = dvector (0, NC+NO-2);   // Control grid streamwise points.
  double*  cy = dvector (0, NY-1);      // Control grid cross-flow points.
  double*  ax = dvector (0, NC-1);      // Body boundary x-points.
  double*  ay = dvector (0, NC-1);      // Body boundary y-points.

  controlGrids  (NC,  NO,  NY,  cx,  cy);
  airfoilSpline (NC,  cx,  ax,  ay);
  printHeader   (NC + NO,  NY);
  printBoundary (NC,  NO,  cx,  ax,  ay);
  printJacobian (NC + NO - 1,   NY,  cx,  cy);

  return EXIT_SUCCESS;
}


static void  controlGrids (int NC, int NO, int NY, double* cx, double* cy)
// ---------------------------------------------------------------------------
// Generate control grid.  Body has unit length.
// ---------------------------------------------------------------------------
{
  const double  P    = 0.1;	// Stretch parameter.
  const double  Q    = 2.0;	// Stretch parameter.
  const double  Lo   = 4.0;	// Outflow length.
  const double  Lc   = 2.0;	// Cross-flow dimension.
  const double  AcLo = 2.0634;	// ArcCosh(Lo).

  int     i;
  double  estar, ratse, r, t;

  // -- Bilinear blend of two stretching functions along body.

  for (i = 0; i < NC; i++) {
    estar = i * 1.0 / (NC - 1);
    ratse = 1.0 - estar;
    r     = stretch1 (P, Q, estar);
    t     = 1.0 - stretch1 (P, Q, ratse);
    cx[i]  = ratse * r + estar * t;
  }

  // -- Cosh-stretched stretching function along outflow.

  for (i = 0; i < NO; i++) {
    estar = i * 1.0 / (NO - 1);
    r     = 1.0 + cosh(AcLo * estar) * stretch1 (P, Q, estar);
    cx[NC + i - 1] =  r;
  }

  // -- Exponential stretching in cross-flow direction.

  for (i = 0; i < NY; i++) {
    estar = i * Lc / (NY - 1);
    r     = cosh (estar) - 1.0;
    cy[i] = r;
  }
}


static void  airfoilSpline (int NC, const double* cx, double* ax, double* ay)
// ---------------------------------------------------------------------------
// Integrate ODE for NACA0012 upper surface to generate tabular data for
// surface coordinate and chord coordinate.  Then generate cubic spline
// coefficients to interpolate from surface to chord coordinate.  Finally
// load airfoil upper-surface x & y coordinates ax & ay according to
// control grid locations supplied in cx.
// ---------------------------------------------------------------------------
{
  const int  PMAX  = 513;
  double*    s     = dvector (1, 1);
  double**   state = dmatrix (1, PMAX, 1, 2);

  int        nok, nbad, nmax;
  double     arclength;

  s[1] = 0.0;

  odeint (s, 1, EPSm12, 1.0, EPSm12, EPSm12, EPSDP, &nok, &nbad, PMAX,
	  &nmax, state, EPSm3, NACAdsdx, (void*) &NACA0012, rkqs);

  arclength = state[nmax][2];

  double* x = dvector (0, nmax-1);
  double* r = dvector (0, nmax-1);
  double* p = dvector (0, nmax-1);

  Veclib::copy   (nmax, &state[1][1],   2, x, 1);
  Veclib::copy   (nmax, &state[1][2],   2, r, 1);
  Veclib::spline (nmax, 1.0e35, 1.0e35, r, x, p);

  ax[0] = ay[0] = 0.0;

  for (int i = 1; i < NC; i++) {
    ax[i] = Veclib::splint (nmax, cx[i] * arclength, r, x, p);
    ay[i] = NACAfoil (ax[i], (void*) &NACA0012);
  }

  ax[NC-1] = 1.0;
  ay[NC-1] = 0.0;

  freeDvector (s, 1);
  freeDmatrix (state, 1, 1);
  freeDvector (x, 0);
  freeDvector (r, 0);
  freeDvector (p, 0);
}


static void  printHeader (int NX, int NY)
// ---------------------------------------------------------------------------
// Refer to cgrid.f documentation.
//
// Dissipation coefficient b4 < 0.125;
// ---------------------------------------------------------------------------
{ 
  cout.setf(ios::fixed);
  cout << NX + NX - 3 << "  " << NY << " " << b4  << "  1"  << endl;
  cout << endl;
}


static void  printBoundary (int NC, int NO, const double* cx,
			                    const double* ax, const double* ay)
// ---------------------------------------------------------------------------
// Refer to cgrid.f documentation.
//
// These are the x--y locations along the j = 1 line.
// ---------------------------------------------------------------------------
{
  int  i;

  for (i = 0; i < NO; i++)
    cout << cx[NC + NO - 2 - i] << "  " << 0.0 << endl;

  for (i = 1; i < NC; i++)
    cout << ax[NC - 1 - i] << "  " << -ay[NC - 1 - i] << endl;

  for (i = 1; i < NC; i++)
    cout << ax[i]          << "  " <<  ay[i]          << endl;

  for (i = 1; i < NO; i++)
    cout << cx[NC - 1 + i]     << "  " << 0.0 << endl;

  cout << endl;
}


static void  printJacobian (int NX, int NY, const double* cx, const double* cy)
// ---------------------------------------------------------------------------
// Refer to cgrid.f documentation.
//
// These are the cell areas at each control grid point, in column-major order.
// Along edges, use areas in the inwards-projected direction.
// At interior points, use average of surrounding cell areas.
// ---------------------------------------------------------------------------
{
  int     i, j;
  double  area;

  // -- LH edge.

  area = (cx[NX - 1] - cx[NX - 2]) * (cy[1] - cy[0]);
  cout << area << endl;

  for (j = 1; j < NY - 1; j++) {
    area = (cx[NX - 1] - cx[NX - 2]) * 0.5 * (cy[j + 1] - cy[j - 1]);
    cout << area << endl;
  }

  area = (cx[NX - 1] - cx[NX - 2]) * (cy[NY - 1] - cy[NY - 2]);
  cout << area << endl;

  // -- Internal locations up to reflection line.

  for (i = 1; i < NX - 1; i++) {
    area = 0.5 * (cx[NX - i] - cx[NX - i - 2]) * (cy[1] - cy[0]);
    cout << area << endl;

    for (j = 1; j < NY - 1; j++) {
      area = 0.25 * (cx[NX - i] - cx[NX - i - 2]) * (cy[j + 1] - cy[j - 1]);
      cout << area << endl;
    }
    
    area = 0.5 * (cx[NX - i] - cx[NX - i - 2]) * (cy[NY - 1] - cy[NY - 2]);
    cout << area << endl;
  }

  // -- Reflection line.

  area = (cx[1] - cx[0]) * (cy[1] - cy[0]);
  cout << area << endl;

  for (j = 1; j < NY - 1; j++) {
    area = (cx[1] - cx[0]) * 0.5 * (cy[j + 1] - cy[j - 1]);
    cout << area << endl;
  }

  area = (cx[1] - cx[0]) * (cy[NY - 1] - cy[NY - 2]);
  cout << area << endl;

  // -- Internal locations up to end line.

  for (i = 1; i < NX - 1; i++) {
    area = 0.5 * (cx[i + 1] - cx[i - 1]) * (cy[1] - cy[0]);
    cout << area << endl;

    for (j = 1; j < NY - 1; j++) {
      area = 0.25 * (cx[i + 1] - cx[i - 1]) * (cy[j + 1] - cy[j - 1]);
      cout << area << endl;
    }
    
    area = 0.5 * (cx[i + 1] - cx[i - 1]) * (cy[NY - 1] - cy[NY - 2]);
    cout << area << endl;
  }

  // -- RH edge.

  area = (cx[NX - 1] - cx[NX - 2]) * (cy[1] - cy[0]);
  cout << area << endl;

  for (j = 1; j < NY - 1; j++) {
    area = (cx[NX - 1] - cx[NX - 2]) * 0.5 * (cy[j + 1] - cy[j - 1]);
    cout << area << endl;
  }

  area = (cx[NX - 1] - cx[NX - 2]) * (cy[NY - 1] - cy[NY - 2]);
  cout << area << endl;

}
