///////////////////////////////////////////////////////////////////////////////
// svv.C: provide SVV-stabilized differentiation operator matrices.
//
// The operating convention is similar to what the equivalent femlib
// routines provide: we only return pointers if the input pointers are
// non-null.
//
// The tokens SVV_EPSx and SVV_Mx (where 0 < N_P < SVV_MN and 0 _N_Z/2
// < SVV_MZ) should be pre-defined for both x-y and z (Fourier)
// directions, see routines coeffs and coeffs_z below.
//
// Copyright (c) 2004 <--> $Date$, Hugh Blackburn
//
// [1] C Xu & R Pasquetti (2004), 'Stabilized spectral element computations
//     of high Reynolds number incompressible flows', JCP 196, 680-704.
// [2] HM Blackburn & S Schmidt (2003), 'Spectral element filtering
//     techniques for large eddy simulation with dynamic estimation',
//     JCP 186, 610-629.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <map>

#include "dns.h"
    
namespace SVV {

typedef struct { real_t* dv; real_t* dt; } vector_pair;

const real_t* coeffs (const int_t np)
// ---------------------------------------------------------------------------
// Return from internal storage a pointer to an array S of SVV filter
// weights for x and y direction.  We use the tokens SVV_MN and
// SVV_EPSN (see [1]) to define this. Specifically, we create
//
//   S = 1 + eps_N/nu * Q
//
// If the relevant vector of weights doesn't exist, create it first.
// ---------------------------------------------------------------------------
{
  static map<int_t, real_t*>    cmap;
  map<int_t, real_t*>::iterator c = cmap.find (np);

  if (c == cmap.end()) {
    const real_t eps = Femlib:: value ("SVV_EPSN");
    const int_t  mN  = Femlib:: value ("SVV_MN");
    const real_t nu  = Femlib:: value ("KINVIS");

    const int_t  N = np - 1;
    real_t*      svvcoeff = new real [np];
    int_t        i;

    cmap[np] = svvcoeff;
    Veclib::zero (np, svvcoeff, 1);
    for (i = mN+1; i < np; i++) svvcoeff[i] = exp (-sqr ((N-i)/(1.0*mN-i)));
    for (i = 0; i < np; i++)    svvcoeff[i] = 1.0 + eps/nu * svvcoeff[i];
  }

  return cmap[np];
}


const real_t* coeffs_z (const int_t numModes)
// ---------------------------------------------------------------------------
// Return from internal storage a pointer to an array S of SVV filter
// weights for the z (Fourier) direction. We use the tokens SVV_MZ and
// SVV_EPSZ (see [1]) to define this. Specifically, we create
//
//   S = 1 + eps_N/nu * Q.
//
// If the relevant vector of weights doesn't exist, create it first.
// ---------------------------------------------------------------------------
{
  static map<int_t, real_t*>    czmap;
  map<int_t, real_t*>::iterator c = czmap.find (numModes);

  if (c == czmap.end()) {
    const real_t eps  = Femlib:: value ("SVV_EPSZ");
    const int_t  mN   = Femlib:: value ("SVV_MZ");
    const real_t nu   = Femlib:: value ("KINVIS");

    const int_t  N    = Geometry::nMode() - 1;
    const int_t  base = Geometry::baseMode();

    real_t*      svvcoeff = new real [numModes];
    int_t        i;

    czmap[numModes] = svvcoeff;
    Veclib::zero (numModes, svvcoeff, 1);
    for (i = base; i < base + numModes; i++) 
      if (i > mN) svvcoeff[i-base] = exp (-sqr ((N-i)/(1.0*mN -i))) ;
    for (i = 0; i < numModes; i++)  svvcoeff[i] = 1.0 + eps/nu * svvcoeff[i];
  }

  return czmap[numModes];
}


void operators (const int_t    np ,
		const real_t** SDV,
		const real_t** SDT)
// ---------------------------------------------------------------------------
// Return from internal storage pointers to (flat) arrays that provide
// the SVV-stabilized derivative operator SDV and its transpose SDT
// which will be used in computing Helmholtz matrices (or their
// equivalent operators).  Again, if these derivative operators do not
// pre-exist, we create them and leave on internal static storage.
//
//   SDV = [M^{-1}][diag(1+SVV_EPSN/KINVIS*Q)][M][DV],
//
// where [DV] is the standard spectral element derivative
// operator. [M] is the discrete polynomial transform matrix that
// takes a set of nodal values into a polyomial space, see
// e.g. Appendix A of [2].  Below we have compile-in support for
// transforms either to Legendre polynomial space (default) or
// modal/hp/Jacobi polynomial space.
// ---------------------------------------------------------------------------
{
  static map<int_t, vector_pair> dmap;
  map <int_t, vector_pair>::iterator d = dmap.find (np);

  if (d == dmap.end()) {
    real_t*        dv = new real [sqr (np)];
    real_t*        dt = new real [sqr (np)];
    const real_t*  S  = SVV::coeffs (np);
    vector<real_t> sqrtS (np);
    int_t          i, j;
    const real_t   *DV;  	// -- Lagrange interpolant operator matrix.
    const real_t   *MF;		// -- Forward polynomial transform matrix.
    const real_t   *MI;		// -- Inverse polynomial transform matrix.

    // -- Gather up the various matrices.

    for (i = 0; i < np; i++) sqrtS[i] = sqrt (S[i]);
    Femlib::quadrature (0, 0, &DV, 0, np, 'L', 0.0, 0.0);

#if 1  // -- Legendre polynomial transform.
    Femlib::legTran (np, &MF, 0, &MI, 0, 0, 0);
#else  // -- Modal polynomial transform.
    Femlib::modTran (np, &MF, 0, &MI, 0, 0, 0);
#endif

    // -- Multiply them.

    for (i = 0; i < np; i++)
      Veclib::smul (np, sqrtS[i], MF + i*np, 1, dv + i*np, 1);
    Blas::mxm (MI, np, dv, np, dt, np);
    Blas::mxm (dt, np, DV, np, dv, np);

    // -- Make the transpose.

    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++)
	dt[Veclib::row_major (j, i, np)] = dv[Veclib::row_major (i, j, np)];

    // -- Set the map storage.

    dmap[np].dv = dv;
    dmap[np].dt = dt;
  }

  if (SDV) *SDV = const_cast<const real_t*> (dmap[np].dv);
  if (SDT) *SDT = const_cast<const real_t*> (dmap[np].dt);
}
}

#if 0
// -- Self-contained testing.
//g++ svv.C -I. -I../semtex/include -L../semtex/lib/Darwin -L/sw/lib -lfem -lvec -framework Accelerate -lg2c
int main (int    argc,
	  char** argv)
{
  const int     np = 6;
  int_t         i, j;
  const real_t* c;
  const real_t* S;
  const real_t* I;

  Femlib::initialize(&argc, &argv);
  Femlib::ivalue ("SVV_MN",   3);
  Femlib:: value ("SVV_EPSN", 0.2);

  c = SVV::coeffs(np);

  for (i = 0; i < np; i++) cout << c[i] << "  ";
  cout << endl;

  Femlib::quadrature (0, 0, &S, 0, np, 'L', 0.0, 0.0);

  cout << "-- The original matrix D" << endl;
  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      cout << S[Veclib::row_major (i, j, np)] << "  ";
    }
    cout << endl;
  }

  SVV::operators (np, &S, &I);

  cout << "-- With SVV modification" << endl;
  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      cout << S[Veclib::row_major (i, j, np)] << "  ";
    }
    cout << endl;
  }

  cout << "-- Transpose" << endl;
  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      cout << I[Veclib::row_major (i, j, np)] << "  ";
    }
    cout << endl;
  }

  return EXIT_SUCCESS;
}
#endif
