#include <cmath>
#include <iostream>

using namespace std;

#include <utility.h>
#include <veclib.h>

#include <Grid.h>


int main()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int  PMAX     = 513;
  double*    s        = dvector (1,    1);
  double**   state    = dmatrix (1, PMAX, 1, 2);
  NACAparam  NACA0012 = {0.12, 1.4845, -0.63, -1.758, 1.4215, -0.518};
  int        nok, nbad, nmax;

  s[1] = 0.0;

  odeint (s, 1, EPSm12, 1.0, EPSm12, EPSm12, EPSDP, &nok, &nbad, PMAX,
	  &nmax, state, EPSm3, NACAdsdx, &NACA0012, rkqs);

  cout << "nok: " << nok << " nbad: " << nbad << " nmax: " << nmax << endl;

  double arclength = state[nmax][2];
  cout << "arclength: " << arclength << endl;

  double* x = dvector (0, nmax-1);
  double* r = dvector (0, nmax-1);
  double* p = dvector (0, nmax-1);

  Veclib::copy   (nmax, &state[1][1], 2, x, 1);
  Veclib::copy   (nmax, &state[1][2], 2, r, 1);
  Veclib::spline (nmax, 1.0e35, 1.0e35, r, x, p);

  for (int i = 0; i <= 100; i++) {
    double xp = Veclib::splint (nmax, i*arclength/100.0, r, x, p);
    double yp = NACAfoil (xp, &NACA0012);

    cout << xp << "  " << yp << endl;
  }

  return 0;
}
