///////////////////////////////////////////////////////////////////////////////
// SGSS.C: calculate subgrid-scale stress for LES.
//
// Copyright (c) 1999 Hugh Blackburn
//
// This is where turbulence modelling takes place: the task is to compute
// the subgrid-scale stress tensor from the resolved-scale velocity field.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <les.h>

static void filters    (const integer, const integer,
			real*, real*, real*, real*);
static void strainRate (const Domain*, AuxField***, AuxField***);


void SGSS (const Domain* D ,
	   AuxField***   Us,
	   AuxField***   Uf,
	   matrix<real>& Ut)
// ---------------------------------------------------------------------------
// Compute the subgrid-scale stress estimate.
//
// On entry, D contains old velocity (and pressure) fields, and the
// first levels of Us & Uf are free.
//
// ---------------------------------------------------------------------------
{
  integer           i, j, ij;
  const integer     DIM    = Geometry::nDim();
  const integer     np     = Geometry::nP();
  const integer     nZ     = Geometry::nZ();
  const integer     nZP    = Geometry::nZProc();
  const integer     nP     = Geometry::planeSize();
  const integer     nPP    = Geometry::nBlock();
  const integer     nPR    = Geometry::nProc();
  const integer     nTot   = Geometry::nTotProc();
  const integer     nZ32   = (nPR > 1) ? nZP : (3 * nZ) >> 1;
  const integer     nTot32 = nZ32 * nP;
  vector<real>      work ((2 * DIM) * nTot32 + nZ + 2 + 2 * np);
  vector<real*>     uu32 (DIM);
  vector<real*>     uv32 (DIM);
  vector<AuxField*> uu   (DIM);
  vector<AuxField*> uv   (DIM);
  AuxField*         utmp = D -> u[DIM];
  real              *expand_ft, *shrink_ft, *expand_lt, *shrink_lt;

  strainRate (D, Us, Uf);

  // -------------------------------------------------------------------------

  for (i = 0; i < DIM; i++) {
    uu  [i] = Us[i][0];
    uv  [i] = Uf[i][0];
    uu32[i] = work() +  i        * nTot32;
    uv32[i] = work() + (i + DIM) * nTot32;
  }

  // -- Create filters.  Fourier filters are zero-phase (all real).

  expand_ft = uv32[DIM-1] + nTot32;
  shrink_ft = expand_ft + (nZ >> 1) + 1;
  expand_lt = shrink_ft + (nZ >> 1) + 1;
  shrink_lt = expand_lt + np;

  filters (nZ >> 1 + 1, np, expand_ft, shrink_ft, expand_lt, shrink_lt);

  // -- Build 1st terms on RHS; attenuated products of amplified variables.

  // -- Complete transform to basis space, and amplify the small scales.

  for (i = 0; i < DIM; i++) {
    *uu[i] = *D -> u[i];
    uu [i] -> DLT2D  (+1, expand_lt);
    uu [i] -> DFfilt (    expand_ft);
  }

  // -- Transform to physical space.

  for (i = 0; i < DIM; i++) {
    uu[i] -> DLT2D       (-1);
    uu[i] -> transform32 (-1, uu32[i]);
  }

  // -- Off-diagonal product terms.

  for (i = 0; i < DIM; i++)
    for (j = i + 1; j < DIM; j++) {
      ij = i + j - 1;
      Veclib::vmul (nTot32, uu32[i], 1, uu32[j], 1, uv32[ij], 1);
    }

  // -- Diagonal product terms.

  for (i = 0; i < DIM; i++)
    Veclib::vmul (nTot32, uu32[i], 1, uu32[i], 1, uu32[i], 1);

  // -- Transform to basis space, attenuate fine scales, invert DLT.

  for (i = 0; i < DIM; i++) {
    uu[i] -> transform32 (+1, uu32[i]);
    uu[i] -> DFfilt      (    shrink_ft);
    uu[i] -> DLT2D       (+1, shrink_lt);
    uu[i] -> DLT2D       (-1);

    uv[i] -> transform32 (+1, uv32[i]);
    uv[i] -> DFfilt      (    shrink_ft);
    uv[i] -> DLT2D       (+1, shrink_lt);
    uv[i] -> DLT2D       (-1);
  }

  // -- 2nd terms on RHS: products in physical space for unfiltered fields.

  // -- Transform to physical space.

  for (i = 0; i < DIM; i++)
    D -> u[i] -> transform32 (-1, uu32[i]);

  // -- Off-diagonal product terms.

  for (i = 0; i < DIM; i++)
    for (j = i + 1; j < DIM; j++) {
      ij = i + j - 1;
      Veclib::vmul (nTot32, uu32[i], 1, uu32[j], 1, uv32[ij], 1);
    }

  // -- Diagonal product terms.

  for (i = 0; i < DIM; i++)
    Veclib::vmul (nTot32, uu32[i], 1, uu32[i], 1, uu32[i], 1);

  // -- Transform products back to Fourier, and subtract off to get SGSS.

  for (i = 0; i < DIM; i++) {
    tmp -> transform32 (+1, uu32[i]);
    *uu[i] -= *tmp;

    tmp -> transform32 (+1, uv32[i]);
    *uv[i] -= *tmp;
  }
}


static void filters (const integer nFourier ,
		     const integer nLegendre,
		     real*         expand_ft,
		     real*         shrink_ft,
		     real*         expand_lt,
		     real*         shrink_lt)
// ---------------------------------------------------------------------------
// Create erfc-type filters and their inverses.
//
// The Fourier filters are of length NZ / 2 + 1 (so notionally they take
// the Nyquist frequency into account too), and are zero-phase (i.e. only
// real coefficients).  The Legendre filters are long enough to fit
// along the side of an element (i.e. Geometry::nP() long).
//
// At present, erfc or Boyd--Vandeven filters are computed.  These are
// defined by three parameters (supplied in the TOKENS section).
//   F_ROLL,  real number in 0--1 that defines where filter roll-off begins;
//   F_ORDER, integer that defines filter order (e.g. 2 or more);
//   F_ATTEN, real number in 0--1 that defines degree of attenuation.
//
// The parameters are assumed to be the same for both Legendre and
// Fourier filters shrink_lt and shrink_ft.  The filters expand_lt
// and expand_ft are their pointwise inverses.
// ---------------------------------------------------------------------------
{
  static integer order = 0;
  static real    lag   = 0.0;
  static real    atten = 0.0;
  integer        i, n, nm;

  if (!order) {
    order = (integer) Femlib::value ("F_ORDER");
    lag   =           Femlib::value ("F_ROLL" );
    atten =           Femlib::value ("F_ATTEN");
  }

  nm = nFourier  - 1;
  Femlib::erfcFilter (nm, order, (integer) (lag * nm), atten, shrink_ft);
  for (i = 0; i <= nm; i++) expand_ft[i] = 1.0 / shrink_ft[i];
  
  nm = nLegendre - 1;
  Femlib::erfcFilter (nm, order, (integer) (lag * nm), atten, shrink_lt);
  for (i = 0; i <= nm; i++) expand_lt[i] = 1.0 / shrink_lt[i];
}



static void strainRate (const Domain* D ,
			AuxField***   Us,
			AuxField***   Uf)
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

    AuxField* tp1 = Us[0][0];
    AuxField* tp2 = Us[1][0];
  
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
	if   (j > i) *Uf[i + j - 1][0]  = *tp1;
	else         *Uf[i + j - 1][0] += *tp1;
      }
  
    for (i = 0; i < DIM; i++)
      for (j = i + 1; j < DIM; j++)
	*Uf[i + j - 1][0] *= 0.5;
  
    // -- Diagonal.

    for (i = 0; i < DIM; i++) {
      (*Us[i][0] = *D -> u[i]) . gradient (i);
      if (i == 2) (*Us[2][0] += *D -> u[1]) . divR();
    }

  } else {			// -- Cartesian geometry.

    // -- Off-diagonal terms.

    AuxField* tmp = Us[0][0];

    for (i = 0; i < DIM; i++)
      for (j = 0; j < DIM; j++) {
	if (j == i) continue;
	(*tmp = *D -> u[i]) . gradient (j);
	if   (j > i) *Uf[i + j - 1][0]  = *tmp;
	else         *Uf[i + j - 1][0] += *tmp;
      }
      
    for (i = 0; i < DIM; i++)
      for (j = i + 1; j < DIM; j++)
	*Uf[i + j - 1][0] *= 0.5;

    // -- Diagonal.

    for (i = 0; i < DIM; i++)
      (*Us[i][0] = *D -> u[i]) . gradient (i);
  }
}
