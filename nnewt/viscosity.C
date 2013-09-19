///////////////////////////////////////////////////////////////////////////////
// viscosity.C: calculate non-Newtonian viscosity for nnewt.
//
// Copyright (c) 1998 <--> $Date$,
//   Murray Rudman, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <nnewt.h>

static int_t NCOM, NDIM;

static void strainRate  (const Domain*, AuxField**, AuxField**);
static void viscoModel  (const Domain*, AuxField**, AuxField**, AuxField*);


void viscosity (const Domain* D  ,
		AuxField**    Us ,
		AuxField**    Uf ,
		AuxField*     NNV)
// ---------------------------------------------------------------------------
// Compute the strain-rate field and NNV, the non-Newtonian viscosity
// field from it.
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
  NDIM = Geometry::nDim();
  NCOM = (D -> u) . size() - 1;

  strainRate  (D, Us, Uf);
  viscoModel  (D, Us, Uf, NNV);

  ROOTONLY NNV -> addToPlane (0, Femlib::value ("-KINVIS"));

#if defined (NOMODEL)
  *NNV = 0.0;
  ROOTONLY NNV -> addToPlane (0, Femlib::value ("-REFVIS"));
#endif
}


static void strainRate (const Domain* D ,
			AuxField**    Us,
			AuxField**    Uf)
// ---------------------------------------------------------------------------
// On entry D contains the velocity fields Ui and the first-level
// areas of Us and Uf are free.  Construct the symmetric strain-rate
// tensor terms.
// 
// For Cartesian coordinates
//           
//             / du/dx  1/2(du/dy + dv/dx)  1/2(du/dz + dw/dx) \
//       S =   |    .           dv/dy       1/2(dv/dz + dw/dy) |,
//       ~     \    .             .                 dw/dz      /
//
// Cartesian is straightforward: leave the diagonal terms in Us and
// the off-diagonal terms in Uf.  This works fine for 2D/2C, 2D/3C, 3D/3C.
// I.e. stored in (3D/3C shown, other cases are obvious restrictions):
//
//             / Us[0]  Uf[0]  Uf[1] \
//       S =   |    .   Us[1]  Uf[2] |,
//       ~     \    .    .     Us[2] /
//
// while for cylindrical coordinates
//           
//         / du/dx  1/2(du/dy + dv/dx)  1/2(1/y*du/dz +    dw/dx)    \
//   S =   |    .           dv/dy       1/2(1/y*dv/dz + dw/dy - w/y) |.
//   ~     \    .             .             1/y*dw/dz +     v/y      /
//
//
// Cylindrical is more complicated owing to the cross-coupling terms,
// especially v/y.  Always there is enough total storage in Us and Uf,
// we just have to pack it in different ways.
//
// 2D/2C stored like:
//             / Us[0]  Uf[0]   0    \
//       S =   |    .   Us[1]   0    |,
//       ~     \    .    .     Uf[1] /
//
// 2D/3C stored like:
//             / Us[0]  Uf[0]  Uf[1] \
//       S =   |    .   Us[1]  Uf[2] |,
//       ~     \    .    .     Us[2] /
//
// 3D/3C as for Cartesian.
// ---------------------------------------------------------------------------
{
  int_t i, j;

  if (Geometry::cylindrical()) {

    AuxField* tp1 = Us[0];
    AuxField* tp2 = Us[1];

    if (NDIM == 2 && NCOM == 2) {

      (*Uf[0] = *D -> u[0]) . gradient (1);
      (*tp1   = *D -> u[1]) . gradient (0);
      (*Uf[0] += *tp1) *= 0.5;

      (*Us[0] = *D -> u[0]) . gradient (0);
      (*Us[1] = *D -> u[1]) . gradient (1);

      // -- "Diagonal" hoop-stress term v/y placed into Uf[1].

      (*Uf[1] = *D -> u[1]) . divY();

      // -- Then any on-axis terms are replaced using l'Hopital's rule.

      D -> u[0] -> overwriteForGroup ("axis", Us[1], Uf[1]);

    } else if (NDIM == 2 && NCOM == 3) {

      (*Uf[0] = *D -> u[0]) . gradient (1);
      (*tp1   = *D -> u[1]) . gradient (0);
      (*Uf[0] += *tp1) *= 0.5;
      
      (*Uf[1] = *D -> u[2]) . gradient (0) *= 0.5;

      (*Uf[2] = *D -> u[2]) . gradient (1);
      (*tp1 = *D -> u[2]) . divY ();
      D -> u[0] -> overwriteForGroup ("axis", Uf[2], tp1);
      (*Uf[2] -= *tp1) *= 0.5;	

      // -- Note Us gets an "off-diagonal" term.

      (*Us[0] = *D -> u[0]) . gradient (0);
      (*Us[1] = *D -> u[1]) . gradient (1);
      (*Us[2] = *D -> u[1]) . divY ();
      D -> u[0] -> overwriteForGroup ("axis", Us[1], Us[2]);

    } else {
      
      // -- The "standard" 3D/3C case.

      for (i = 0; i < NCOM; i++)
	for (j = 0; j < NDIM; j++) {
	  if (j == i) continue;
	  if (i == 2 && j == 1) {
	    (*tp1 = *D -> u[2]) . gradient (1);
	    (*tp2 = *D -> u[2]) . divY();
	    D -> u[0] -> overwriteForGroup ("axis", tp1, tp2);
	    *tp1 -= *tp2;
	  } else {
	    (*tp1 = *D -> u[i]) . gradient (j);
	    if (j == 2) tp1 -> divY();
	  }
	  if   (j > i) *Uf[i + j - 1]  = *tp1;
	  else         *Uf[i + j - 1] += *tp1;
	}
  
      for (i = 0; i < NCOM; i++) *Uf[i] *= 0.5;
  
      // -- Diagonal.

      for (i = 0; i < NDIM; i++) {
	(*Us[i] = *D -> u[i]) . gradient (i);
	if (i == 2) {
	  Us[2] -> divY();
	  (*tp1 = *D -> u[1]) . divY();
	  D -> u[0] -> overwriteForGroup ("axis", Us[1], tp1);
	  *Us[2] += *tp1;
	}
      }
    }

  } else {			// -- Cartesian geometry.

    // -- Off-diagonal terms.

    AuxField* tmp = Us[0];

    for (i = 0; i < NCOM; i++) {
      for (j = 0; j < NDIM; j++) {
	if (i == j) continue;
	(*tmp = *D -> u[j]) . gradient (i);
	if (j > i) *Uf[i + j - 1]  = *tmp;
	else       *Uf[i + j - 1] += *tmp;
      }
    }
      
    for (i = 0; i < NCOM; i++) *Uf[i] *= 0.5;

    // -- Diagonal.

    for (i = 0; i < NDIM; i++) (*Us[i] = *D -> u[i]) . gradient (i);
  }
}


static void viscoModel (const Domain* D  ,
			AuxField**    Us ,
			AuxField**    Uf ,
			AuxField*     NNV)
// ---------------------------------------------------------------------------
// On entry the first-level areas of Us & Uf contain the components of
// the strain-rate tensor S and NNV contains the old values of
// non-Newtonian viscosity.  On exit, by default NNV contains (for
// example) the Cross model non-Newtonian viscosity field:
//
// NNV = (\nu_0 + \nu_\inf (K|S|)^N)/ (1+(K|S|)^N), where
//
// |S| = sqrt{2[(S11)^2 + (S22)^2 + (S33)^2 + 2(S12)^2 + 2(S13)^2 + 2(S23)^2]},
//
// \nu_0 = KINVIS from session file (i.e. REFVIS below), \nu_\inf is
// the infinite shear rate viscosity and K and N are the cross model
// parameters stored in CROSS_K and CROSS_N.  Other non-Newtonian
// models can be defined as desired.
//
// As noted in the header to NS.C, NNV = -KINVIS for debugging (NB: Fourier!).
//
// NB: the products in |S| are only dealiased for single-processor operation.
// ---------------------------------------------------------------------------
{
  const int_t    nP     = Geometry::nPlane();
  const int_t    NP     = Geometry::planeSize();
  const int_t    nPR    = Geometry::nProc();
  const int_t    nZ     = Geometry::nZProc();
  const int_t    nZ32   = Geometry::nZ32();
  const int_t    nTot32 = nZ32 * NP;
  const real_t   RefVis = Femlib::value ("KINVIS");

  static vector<real_t> work (2 * nTot32);

  int_t   i, j;
  real_t* tmp = &work[0];
  real_t* sum = tmp + nTot32;

  Veclib::zero (nTot32, sum, 1);

#if 0 // -- Dump 2x2 strain rate tensor components for examination.

  *D -> ua[0] = *Us[0];
  *D -> ua[1] = *Us[1];
  *D -> ua[2] = *Uf[0];
  *D -> ua[3] = *Uf[1];
  return;
#endif

  if (Geometry::cylindrical()) {

    if (NDIM == 2 && NCOM == 2) {

      Uf[0] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);

      Blas::scal (nTot32, 2.0, sum, 1);

      Uf[1] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
      Us[0] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
      Us[1] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);

      Blas::scal (nTot32, 2.0, sum, 1);

    } else {			// -- Other two cases (2D/3C and 3D/3C).
      
      for (i = 0; i < NCOM; i++) {
	Uf[i] -> transform32 (INVERSE, tmp);
	Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
      }
      Blas::scal (nTot32, 2.0, sum, 1);

      for (i = 0; i < NCOM; i++) {
	Us[i] -> transform32 (INVERSE, tmp);
	Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
      }
      Blas::scal (nTot32, 2.0, sum, 1);

    }

  } else {			// -- Cartesian.

    if (NCOM == 2) {
      Uf[0] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
    } else {
      for (i = 0; i < NCOM; i++) {
	Uf[i] -> transform32 (INVERSE, tmp);
	Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
      }
    }
    Blas::scal (nTot32, 2.0, sum, 1);

    for (i = 0; i < NDIM; i++) {
      Us[i] -> transform32 (INVERSE, tmp);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
    }
    Blas::scal (nTot32, 2.0, sum, 1);      
  }

  Veclib::vsqrt (nTot32, sum, 1, sum, 1);

#if 0 // -- Dump strain rate in pressure for debugging (and return).
  D -> u[NCOM] -> transform32 (FORWARD, sum);
  // *D -> u[NCOM] = 1.; // -- Test.
#else

  // -- Smooth shear rate.

  D -> u[0] -> smooth (nZ32, sum);
  
  // -- At this point we have |S| in physical space (in sum).

  if (Femlib::ivalue ("PowerLaw")) {

    // -- Power law model: visc = K*S**(N-1).

    const real_t K         = Femlib::value ("PL_K");
    const real_t N         = Femlib::value ("PL_N") - 1.0;
    const real_t Gamma_min = Femlib::value ("PL_ZERO");  // Limit shear.

    Veclib::clipup (nTot32, Gamma_min, sum, 1, tmp, 1);
    Veclib::spow   (nTot32, N, tmp, 1, tmp, 1);
    Veclib::smul   (nTot32, K, tmp, 1, tmp, 1);
    
  } else if (Femlib::ivalue ("HB")) {

    // -- Herschel-Bulkley model: visc = Tau_Y/S + K*S**(N-1).

    const real_t Tau_Y     = Femlib::value ("YIELD_STRESS");
    const real_t K         = Femlib::value ("HB_K");
    const real_t N         = Femlib::value ("HB_N") - 1.0;
    const real_t Gamma_min = Femlib::value ("HB_ZERO");  // Limit shear.

    Veclib::clipup (nTot32, Gamma_min, sum, 1, sum, 1);
    Veclib::spow   (nTot32, N, sum, 1, tmp, 1);
    Veclib::smul   (nTot32, K, tmp, 1, tmp, 1);
    Veclib::sdiv   (nTot32, Tau_Y, sum, 1, sum, 1);
    Veclib::vadd   (nTot32, sum, 1, tmp, 1, tmp, 1);

  } else if (Femlib::ivalue ("CAR_YAS")) {

    // -- Carreau-Yasuda model - visc = nu_inf + (nu_0-nu_inf)/(1+(K*S)^N)^A.

    const real_t K      = Femlib::value ("CY_K");
    const real_t N      = Femlib::value ("CY_N");
    const real_t A      = Femlib::value ("CY_A");
    const real_t mu_0   = Femlib::value ("VISC_ZERO");
    const real_t mu_inf = Femlib::value ("VISC_INF");
    const real_t mu_dif = mu_0 - mu_inf;
    const real_t one    = 1.0;

    Veclib::smul (nTot32, K,   sum, 1, tmp, 1); 
    Veclib::spow (nTot32, N,   tmp, 1, tmp, 1);
    Veclib::sadd (nTot32, one, tmp, 1, tmp, 1);
    Veclib::spow (nTot32, A,   tmp, 1, tmp, 1);
    Veclib::sdiv (nTot32, one, tmp, 1, tmp, 1);
    Veclib::smul (nTot32, mu_dif, tmp, 1, tmp, 1); 
    Veclib::sadd (nTot32, mu_inf, tmp, 1, tmp, 1);

  } else {

    // -- Cross model - visc = nu_inf + (nu_0-nu_inf)/(1+K*(S^N)).

    const real_t K      = Femlib::value ("CROSS_K");
    const real_t N      = Femlib::value ("CROSS_N");
    const real_t mu_0   = Femlib::value ("VISC_ZERO");
    const real_t mu_inf = Femlib::value ("VISC_INF");
    const real_t mu_dif = mu_0 - mu_inf;
    const real_t one    = 1.0;

    Veclib::spow (nTot32, N,   sum, 1, tmp, 1);
    Veclib::smul (nTot32, K,   tmp, 1, tmp, 1); 
    Veclib::sadd (nTot32, one, tmp, 1, tmp, 1);
    Veclib::sdiv (nTot32, one, tmp, 1, tmp, 1);
    Veclib::smul (nTot32, mu_dif, tmp, 1, tmp, 1); 
    Veclib::sadd (nTot32, mu_inf, tmp, 1, tmp, 1);
    
  }

  // -- Smooth non-Newtonian viscosity.

  //  D -> u[0] -> smooth (nZ32, tmp);
  
  // -- Transform back to Fourier space.
  
  NNV -> transform32 (FORWARD, tmp);
#endif
}


