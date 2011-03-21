///////////////////////////////////////////////////////////////////////////////
// plate: generate mesh for a flat plate, with stretched coordinates.
//
// Usage: plate [session]
//
// Default action is to print output suitable for use by sm macro meshplot
// on cout.  If a session file is nonminated, output suitable for fem input
// is also generated, written to session.
//
// Plate LE is at coordinate origin.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#include <utility.h>
#include <veclib.h>

const int    NC = 10;	// Number of points along semi-chord.
const double LC = 0.5;	// Semi-chord length.

const int    NY = 10;	// Number of points in cross flow direction.
const double LY = 1.0;	// Mesh cross-flow semi-length.

const int    NI = 10;    // Number of points in inlet direction.
const double LI = 1.0;	// Mesh inlet to plate LE length.

const int    NO = 20;   // Number of points in outflow direction.
const int    LO = 2.0;	// Mesh TE to outflow length.

extern double stretch1 (double, double, double);


int main (int argc, char** argv)
// ---------------------------------------------------------------------------
// Generate stretched coordinates for grid, print up.
// ---------------------------------------------------------------------------
{
  char*     session = 0;
  ofstream  file;
  if (argc > 1) {
    session = argv[1];
    file.open (session);
  }

  int      i, j, k;
  double   P, Q;
  double   r, estar;

  double*  cx = dvector (0, NC-1);
  double*  cy = dvector (0, NY-1);
  double*  ix = dvector (0, NI-1);
  double*  ox = dvector (0, NO-1);
  
  int      nx    = 2 * (NC - 1) + (NI - 1) + (NO - 1);
  int      ny    = 2 * (NY - 1);
  int      nel   = nx * ny;
  int      nvert = (nx + 1) * (ny + 1);

  double** xgrid = dmatrix (0, ny, 0, nx);
  double** ygrid = dmatrix (0, ny, 0, nx);
  int**    igrid = imatrix (0, ny, 0, nx);

  // -- Stretching function along inlet length.

  P = 0.02;
  Q = 3.5;

  for (i = 0; i < NI; i++) {
    estar = i * 1.0 / (NI - 1);
    r     = LI * (1.0 - stretch1 (P, Q, estar));
    ix[i] = -LI + r;
  }

  // -- Stretching function along plate semi-chord.

  P = 0.1;
  Q = 3.0;

  for (i = 0; i < NC; i++) {
    estar = i * 1.0 / (NC - 1);
    r     = LC * stretch1 (P, Q, estar);
    cx[i] = r;
  }

  // -- Stretching function along outlet length.

  P = 0.05;
  Q = 3.5;

  for (i = 0; i < NO; i++) {
    estar = i * 1.0 / (NO - 1);
    r     = LO * stretch1 (P, Q, estar);
    ox[i] = 2.0 * LC + r;
  }

  // -- Exponential stretching in cross-flow direction.

  P = 0.05;
  Q = 4.0;

  for (i = 0; i < NY; i++) {
    estar = i * 1.0 / (NY - 1);
    r     = LY * stretch1 (P, Q, estar);
    cy[i] = r;
  }

  // -- Print up for sm macros.

  cout << "2 2 1 " << nel << " NR NS NZ NEL" << endl;

  double  x1, x2, y1, y2;

  for (j = 0; j < NY - 1; j++) {
    y1 = cy[j];
    y2 = cy[j + 1];

    for (i = 0; i < NI - 1; i++) { // -- Inlet.
      x1 = ix[i];
      x2 = ix[i + 1];

      cout << x1 << " " <<  y1 << endl;
      cout << x2 << " " <<  y1 << endl;
      cout << x1 << " " <<  y2 << endl;
      cout << x2 << " " <<  y2 << endl;

      cout << x1 << " " << -y2 << endl;
      cout << x2 << " " << -y2 << endl;
      cout << x1 << " " << -y1 << endl;
      cout << x2 << " " << -y1 << endl;
    }

    for (i = 0; i < NC - 1; i++) { // -- First semi-chord.
      x1 = cx[i];
      x2 = cx[i + 1];

      cout << x1 << " " <<  y1 << endl;
      cout << x2 << " " <<  y1 << endl;
      cout << x1 << " " <<  y2 << endl;
      cout << x2 << " " <<  y2 << endl;

      cout << x1 << " " << -y2 << endl;
      cout << x2 << " " << -y2 << endl;
      cout << x1 << " " << -y1 << endl;
      cout << x2 << " " << -y1 << endl;
    }

    for (i = NC - 1; i > 0; i--) { // -- Second semi-chord.
      x1 = 2.0 * LC - cx[i];
      x2 = 2.0 * LC - cx[i-1];

      cout << x1 << " " <<  y1 << endl;
      cout << x2 << " " <<  y1 << endl;
      cout << x1 << " " <<  y2 << endl;
      cout << x2 << " " <<  y2 << endl;

      cout << x1 << " " << -y2 << endl;
      cout << x2 << " " << -y2 << endl;
      cout << x1 << " " << -y1 << endl;
      cout << x2 << " " << -y1 << endl;
    }

    for (i = 0; i < NO - 1; i++) { // -- Inlet.
      x1 = ox[i];
      x2 = ox[i + 1];

      cout << x1 << " " <<  y1 << endl;
      cout << x2 << " " <<  y1 << endl;
      cout << x1 << " " <<  y2 << endl;
      cout << x2 << " " <<  y2 << endl;

      cout << x1 << " " << -y2 << endl;
      cout << x2 << " " << -y2 << endl;
      cout << x1 << " " << -y1 << endl;
      cout << x2 << " " << -y1 << endl;
    }
  }

  if (!session) return EXIT_SUCCESS;

  // -- Fill grid arrays.

  for (j = 0, k = 1; j <= nx; j++)
    for (i = 0; i <= ny; i++)
      igrid[i][j] = k++;

  for (j = 0; j <= nx; j++) {
    for (i = 0; i < NY; i++) {
      ygrid[NY+i-1][j] =  cy[i];
      ygrid[NY-i-1][j] = -cy[i];
    }
  }

  for (i = 0; i <= ny; i++) {
    for (j = 0; j < NI; j++)
      xgrid[i][j] = ix[NI - j - 1];
    for (j = 0; j < NC; j++) {
      xgrid[i][NI + j - 1]           =            cx[j];
      xgrid[i][NI + NC + NC - j - 3] = 2.0 * LC - cx[j];
    }
    for (j = 0; j < NO; j++)
      xgrid[i][NI + NC + NC + j - 3] = ox[j];
  }

  // -- Print up a session file.

  file << "PROBLEM   { navierstokes geometry 2D-cartesian }" << endl << endl;

  file << "PARAMETER {"         << endl;
  file << "0 options"           << endl << endl;
  file << "6 integer"           << endl;
  file << "2    N_VAR"          << endl;
  file << "2    N_TIME"         << endl;
  file << "5    N_POLY"         << endl;
  file << "100  N_STEP"         << endl;
  file << "100  IO_FLD"         << endl;
  file << "10   IO_HIS"         << endl << endl;
  file << "3 floating point"    << endl;
  file << "0.005 DELTAT"        << endl;
  file << "100.0 Re"            << endl;
  file << "1/Re  KINVIS"        << endl;
  file << "}"                   << endl << endl;

  file << "BOUNDARY {"          << endl;
  file << "3 boundary segments" << endl;
  file << "1 essential 0.0 0.0" << endl;
  file << "2 outflow"           << endl;
  file << "3 wall"              << endl;
  file << "}"                   << endl << endl;

  file << "MESH {" << endl;

  // -- Vertices.

  file << nvert << " vertices" << endl;
  
  for (j = 0; j <= nx; j++)
    for (i = 0; i <= ny; i++)
      file
	<< igrid[i][j]
	  << setw (14)
	    << xgrid[i][j]
	      << setw (14)
		<< ygrid[i][j]
		  << setw (14)
		    << 0.0
		      << endl;
  
  // -- Elements.

  file << endl << nel << " elements" << endl;

  for (j = 0, k = 1; j < nx; j++)
    for (i = 0; i < ny; i++)
      file
	<< k++
	  << "\t4\t"
	    << igrid[i][j]
	      << ' '
		<< igrid[i][j+1]
		  << ' '
		    << igrid[i+1][j+1]
		      << ' '
			<< igrid[i+1][j]
			  << endl;

  // -- Boundaries.

  file << endl << "3 boundaries" << endl;
  
  file << "1\t1\t" << nx + ny + nx + 1;
  for (j = nx; j >= 0; j--)
    file << ' ' << igrid[ny][j];
  for (i = ny - 1; i > 0; i--)
    file << ' ' << igrid[i][0];
  for (j = 0; j <= nx; j++)
    file << ' ' << igrid[0][j];
  file << endl;

  file << "2\t2\t" << ny + 1;
  for (i = 0; i <= ny; i++)
    file << ' ' << igrid[i][nx];
  file << endl;

  file << "3\t3\t" << 4 * (NC - 1);
  for (j = NI - 1; j < NI + NC + NC - 4; j++)
    file << ' ' << igrid[NY-1][j];
  for (j = NI + NC + NC - 3; j > NI - 2; j--)
    file << ' ' << igrid[NY-1][j];
  file << endl;

  file << endl << "0 curves" << endl;
  
  file << "}" << endl;

  file.close ();

  return EXIT_SUCCESS;
}
