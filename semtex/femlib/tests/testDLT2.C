//////////////////////////////////////////////////////////////////////////////
// testDLT2.C: test run 2D Discrete Legendre Transforms.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <Array.h>
#include <Veclib.h>
#include <Femlib.h>
#include <Blas.h>
#include <Utility.h>

#define SIZE 6


int main ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  int            i, j, k, l, p, q;
  const int      np = SIZE, np2 = np * np, np4 = np2 * np2;
  double         ci, *x, *y, *z, *t, *FW, *FT, *BW, *BT, *CF, *CB;
  double         ts, tf;
  const double   *w, *tab;

  x   = new double [8*np2+np*(np+1)+2*np4];
  y   = x   + np2;
  z   = y   + np2;
  t   = z   + np2;
  FW  = t   + np2;
  FT  = FW  + np2;
  BW  = FT  + np2;
  BT  = BW  + np2;
  CF  = BT  + np2;
  CB  = CF  + np4;
  tab = CB  + np4;

  // -- Get GLL grid points & weights, prepare table of Legendre polys.

  Femlib::legCoef (np, &tab);
  Femlib::quad    (LL, np, np, 0, 0, &w, 0, 0, 0, 0);

  // -- Create forward & inverse tensor-product transform matrices.

  for (i = 0; i < np; i++) {
    ci = tab[Veclib::row_major (np, i, np)];
    for (j = 0; j < np; j++) {
      FW[Veclib::row_major (i, j, np)] =
	ci * w[j] * tab[Veclib::row_major (i, j, np)];
      BW[Veclib::row_major (i, j, np)] =
	tab[Veclib::row_major (j, i, np)];
    }
  }

  // -- And their transposes.

  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      FT[Veclib::row_major (i, j, np)] = FW[Veclib::row_major (j, i, np)];
      BT[Veclib::row_major (i, j, np)] = BW[Veclib::row_major (j, i, np)];
    }
  }

  // -- Now manufacture 2D forward, inverse DLT matrices.

  for (k = 0, i = 0; i < np; i++)
    for (j = 0; j < np; j++, k++)
      for (l = 0, p = 0; p < np; p++)
	for (q = 0; q < np; q++, l++) {
	  CF[Veclib::row_major (k, l, np2)]
	    = FW[Veclib::row_major(i,p,np)] * FT[Veclib::row_major(q,j,np)];
	  CB[Veclib::row_major (k, l, np2)]
	    = BW[Veclib::row_major(i,p,np)] * BT[Veclib::row_major(q,j,np)];
	}

  // -- Create original data.

  Veclib::vrandom (np2, x, 1);

  // -- Forward DLT.
/*
  Blas::mxm (FW, np, x,  np, t, np);
  Blas::mxm (t,  np, FT, np, y, np);
*/
  Blas::mxv (CF, np2, x, np2, y);

  // -- Inverse DLT.
/*
  Blas::mxm (BW, np, y,  np, t, np);
  Blas::mxm (t,  np, BT, np, z, np);
*/
  Blas::mxv (CB, np2, y, np2, z);

  // -- Print everything up.

  cout << "Problem size: " << np << " x " << np << endl;

  for (i = 0; i < np2; i++)
    cout << setw(14) << x[i] << setw(14) << y[i] << setw(14) << z[i] << endl;

  // -- Timing.

  ts = dclock();

  for (i = 0; i < 100000; i++) {
    Blas::mxm (FW, np, x,  np, t, np);
    Blas::mxm (t,  np, FT, np, y, np);
  }

  tf = dclock();

  cout << "Time for 100000 T-P transforms:      " << tf - ts << endl;

  ts = dclock();

  for (i = 0; i < 100000; i++)
    Blas::mxv (CF, np2, x, np2, y);

  tf = dclock();

  cout << "Time for 100000 unrolled transforms: " << tf - ts << endl;

  return (EXIT_SUCCESS);
}
