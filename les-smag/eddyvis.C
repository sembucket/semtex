///////////////////////////////////////////////////////////////////////////////
// eddyvis.C: calculate eddy viscosity for LES.
//
// Copyright (c) Hugh Blackburn 1998--1999.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <les.h>

static integer DIM, CYL, ITR_MAX;
static real    EPS2;

static void strainRate  (const Domain*, AuxField**,AuxField**);
static void viscoModel  (const Domain*, AuxField**, AuxField**, AuxField*);
static real RNG_quartic (const real, const real, const real);


void eddyViscosity (const Domain* D ,
		    AuxField**    Us,
		    AuxField**    Uf,
		    AuxField*     EV)
// ---------------------------------------------------------------------------
// Compute the resolved strain-rate field and EV, the eddy viscosity field
// from it.
//
// On entry, D contains velocity (and pressure) fields, and the first
// levels of Us & Uf are free.  The diagonal entries of the rate of
// strain tensor are left in Us, and the off-diagonal entries are left
// in Uf (here for 3D):
//
//                      / Us[0]  Uf[0]  Uf[1] \
//                S =   |    .   Us[1]  Uf[2] |
//                ~     \    .       .  Us[2] /
// ---------------------------------------------------------------------------
{
  DIM     = Geometry::nDim();
  CYL     = Geometry::system() == Geometry::Cylindrical;
  ITR_MAX = (integer) Femlib::value ("STEP_MAX");
  EPS2    = sqr (EPSSP);

  if ((int) Femlib::value ("RNG"))
    ROOTONLY EV -> addToPlane (0, Femlib::value ("KINVIS"));

  strainRate  (D, Us, Uf);
  viscoModel  (D, Us, Uf, EV);
  D -> u[0] -> smooth (EV);

  ROOTONLY EV -> addToPlane (0, Femlib::value ("-KINVIS"));

#if defined(DEBUG)
  *EV = 0.0;
  ROOTONLY EV -> addToPlane (0, Femlib::value ("-REFVIS"));
#endif
}


static void strainRate (const Domain* D ,
			AuxField**    Us,
			AuxField**    Uf)
// ---------------------------------------------------------------------------
// On entry D contains the velocity fields Ui and the first-level
// areas of Us and Uf are free.  Construct the symmetric strain-rate
// tensor terms, leave the diagonal terms Sii (unsummed) in Us and the
// off-diagonal terms Sij (i != j) in Uf.
// 
// For Cartesian coordinates
//           
//             / du/dx  1/2(du/dy + dv/dx)  1/2(du/dz + dw/dx) \
//       S =   |    .           dv/dy       1/2(dv/dz + dw/dy) |,
//       ~     \    .             .                 dw/dz      /
//
// while for cylindrical coordinates
//           
//         / du/dx  1/2(du/dy + dv/dx)  1/2(1/y*du/dz +    dw/dx)    \
//   S =   |    .           dv/dy       1/2(1/y*dv/dz + dw/dy - w/y) |.
//   ~     \    .             .             1/y*dw/dz +     1/y*v    /
// ---------------------------------------------------------------------------
{
  register integer i, j;

  if (CYL) {			// -- Cylindrical geometry.

    AuxField* tp1 = Us[0];
    AuxField* tp2 = Us[1];
  
    for (i = 0; i < DIM; i++)
      for (j = 0; j < DIM; j++) {
	if (j == i) continue;
	if (i == 2 && j == 1) {
	  (*tp1 = *D -> u[2]) . gradient (1);
	  (*tp2 = *D -> u[2]) . divR();
	  *tp1 -= *tp2;
	} else {
	  (*tp1 = *D -> u[i]) . gradient (j);
	  if (j == 2) tp1 -> divR();
	}
	if   (j > i) *Uf[i + j - 1]  = *tp1;
	else         *Uf[i + j - 1] += *tp1;
      }
  
    for (i = 0; i < DIM; i++)
      for (j = i + 1; j < DIM; j++)
	*Uf[i + j - 1] *= 0.5;
  
    // -- Diagonal.

    for (i = 0; i < DIM; i++) {
      (*Us[i] = *D -> u[i]) . gradient (i);
      if (i == 2) (*Us[2] += *D -> u[1]) . divR();
    }

  } else {			// -- Cartesian geometry.

    // -- Off-diagonal terms.

    AuxField* tmp = Us[0];

    for (i = 0; i < DIM; i++)
      for (j = 0; j < DIM; j++) {
	if (j == i) continue;
	(*tmp = *D -> u[i]) . gradient (j);
	if   (j > i) *Uf[i + j - 1]  = *tmp;
	else         *Uf[i + j - 1] += *tmp;
      }
      
    for (i = 0; i < DIM; i++)
      for (j = i + 1; j < DIM; j++)
	*Uf[i + j - 1] *= 0.5;

    // -- Diagonal.

    for (i = 0; i < DIM; i++)
      (*Us[i] = *D -> u[i]) . gradient (i);
  }
}


static void viscoModel (const Domain* D ,
			AuxField**    Us,
			AuxField**    Uf,
			AuxField*     EV)
// ---------------------------------------------------------------------------
// On entry the first-level areas of Us & Uf contain the components of
// the strain-rate tensor S and EV contains the old values of eddy viscosity.
// On exit, by default EV contains the Smagorinsky eddy viscosity field
//
//                 \nu_S = (Cs \Delta)^2 |S|, where
//
// |S| = sqrt {2[(S11)^2 + (S22)^2 + (S33)^2 + 2(S12)^2 + 2(S13)^2 + 2(S23)^2]}
//
// NB (1): the outer factor of 2, sometimes omitted in derivations.
// NB (2): the products in |S| are only dealiased for serial operation.
//
// For RNG, the "decreed" value of C_SMAG = 0.1114, RNG_C = 75, RNG_BIG = 500.
// ---------------------------------------------------------------------------
{
  register integer i, j, k;
  const integer    nP     = Geometry::nPlane();
  const integer    NP     = Geometry::planeSize();
  const integer    nPR    = Geometry::nProc();
  const integer    nZ     = Geometry::nZProc();
  const integer    nZ32   = (nPR > 1) ? nZ : (3 * nZ) >> 1;
  const integer    nTot32 = nZ32 * NP;
  const real       Cs     = Femlib::value ("C_SMAG");
  const real       molvis = Femlib::value ("REFVIS");

  vector<real>     work (2 * nTot32 + nP);
  real*            tmp    = work();
  real*            sum    = tmp + nTot32;
  real*            delta  = sum + nTot32;

  EV -> lengthScale       (delta);
  D  -> u[0] -> smooth (1, delta);

  Veclib::zero (nTot32, sum, 1);
  
  for (i = 0; i < DIM; i++) {
    for (j = i + 1; j < DIM; j++) {
      Uf[i + j - 1] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
    }
  }
  Blas::scal (nTot32, 2.0, sum, 1);

  for (i = 0; i < DIM; i++) {
    Us[i] -> transform32 (INVERSE, tmp);
    Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
  }
  Blas::scal (nTot32, 2.0, sum, 1);

  Veclib::vsqrt (nTot32, sum, 1, sum, 1);

  // -- At this point we have S in physical space.

  if ((int) Femlib::value ("RNG")) {

    register real *S, *E, D2;
    register real cutoff;
    const real    C       = Femlib::value ("RNG_C");
    const real    BIG     = Femlib::value ("RNG_BIG");
    const real    molvis3 = molvis * molvis * molvis;
    const real    Cs2     = Cs  * Cs;
    const real    Cs4     = Cs2 * Cs2;
  
    EV -> transform32 (INVERSE, tmp);

    for (k = 0; k < nZ32; k++) {
      S = sum + k * NP;
      E = tmp + k * NP;
      
      for (i = 0; i < nP; i++) {
	D2 = delta[i] * delta[i];
	cutoff  = Cs4  * D2 * D2 / molvis3;
	cutoff *= S[i] * S[i] * E[i];

	if      (cutoff <= C) 
	  S[i] = molvis;
	else if (cutoff > BIG)
	  S[i] = molvis + Cs2 * D2 * S[i];
	else
	  S[i] = molvis * RNG_quartic (E[i]/molvis, C-1., -cutoff*E[i]/molvis);
      }
    }

  } else {			// -- Smagorinsky.

    Blas::scal   (nP, Cs, delta, 1);
    Veclib::vmul (nP, delta, 1, delta, 1, delta, 1);

    for (k = 0; k < nZ32; k++)
      Veclib::svvtp (nP, molvis, delta, 1, sum + k * NP, 1, sum + k * NP, 1);
  }

  D -> u[0] -> smooth (nZ32, sum);

  // -- Transform back to Fourier space.

  EV -> transform32 (FORWARD, sum);
}


static real RNG_quartic (const real x0,
			 const real a1,
			 const real a0)
// ---------------------------------------------------------------------------
// Solve quartic equation x^4 + a1 x + a0 = 0 by Newton-Raphson.
// Input value x0 is an initial guess for the ratio of the turbulent and
// molecular viscosities.  Minimum returned value for x is 1.0.
// ---------------------------------------------------------------------------
{
  const char routine[] = "RNG_quartic";
  register real x = x0, dx, f, df, x3;
  register int  i = 0;

  for (i = 0; i < ITR_MAX; i++) {
    x3  = x * x * x;
    df  = a1 + 4.0 * x3;
    f   = a0 + x * (a1 + x3);
    dx  = -f / df;
    x  += dx;
    if (dx * dx < EPS2) break;
  }

  if (i == ITR_MAX) message (routine, "failed to converge", ERROR);
#if defined(DEBUG)
  if (x < 1.0) message (routine, "invalid solution",   WARNING);
#endif

  return max (1.0, x);
}
