///////////////////////////////////////////////////////////////////////////////
// eddyvis.C: calculate eddy viscosity for LES.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>

static integer DIM, CYL, C3D;

static void strainRate  (const Domain*,AuxField***,AuxField***);
static void Smagorinsky (AuxField***, AuxField***, AuxField*);


void eddyViscosity (const Domain* D ,
		    AuxField***   Us,
		    AuxField***   Uf,
		    AuxField*     EV)
// ---------------------------------------------------------------------------
// Compute the resolved strain-rate field and E, the eddy viscosity field
// from it.
//
// On entry, D contains velocity (and pressure) fields, and the first
// levels of Us & Uf are free.  The diagonal entries of the rate of
// strain tensor are left in Us, and the off-diagonal entries are left
// in Uf (here for 3D):
//
//                      / Us[0][0]  Uf[0][0]  Uf[1][0] \
//                S =   |    .      Us[1][0]  Uf[2][0] |
//                ~     \    .         .      Us[2][0] /
//
// As a final step, subtract off the difference between the kinematic
// and reference viscosities stored in REFVIS.
// ---------------------------------------------------------------------------
{
  DIM = Geometry::nDim();
  CYL = Geometry::system() == Geometry::Cylindrical;
  C3D = CYL && DIM == 3;

  strainRate  (D,  Us, Uf);
  Smagorinsky (Us, Uf, EV);

#if !defined(DEBUG)
  ROOTONLY EV -> addToPlane (0, -Femlib::value ("REFVIS"));
#endif
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
  integer i, j;

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


static void Smagorinsky (AuxField*** Us,
			 AuxField*** Uf,
			 AuxField*   EV)
// ---------------------------------------------------------------------------
// On entry the first-level areas of Us & Uf contain the components of
// the strain-rate tensor S.  On exit EV contains the Smagorinsky eddy
// viscosity field \nu_S = (Cs \Delta)^2 |S|, where
//
// |S| = sqrt [(S11)^2 + (S22)^2 + (S33)^2 + 2(S12)^2 + 2(S13)^2 + 2(S23)^2].
//
// As noted in the header to NS.C, EV = -KINVIS for debugging (NB: Fourier!).
//
// NB: the products in |S| are only dealiased for single-processor operation.
// ---------------------------------------------------------------------------
{
#if defined(DEBUG)
  *EV = 0.0;
  ROOTONLY EV -> addToPlane (0, -0.5 * Femlib::value ("KINVIS"));

#else

  integer       i, j;
  const integer nP     = Geometry::planeSize();
  const integer nPR    = Geometry::nProc();
  const integer nZ     = Geometry::nZProc();
  const integer nZ32   = (nPR > 1) ? nZ : (3 * nZ) >> 1;
  const integer nTot32 = nZ32 * nP;

  vector<real>  work (2 * nTot32);
  real*         tmp = work();
  real*         sum = tmp + nTot32;

  Veclib::zero (nTot32, sum, 1);
  
  for (i = 0; i < DIM; i++) {
    for (j = i + 1; j < DIM; j++) {
      Uf[i + j - 1][0] -> transform32 (tmp, -1);
      Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
    }
    Blas::scal (nTot32, 2.0, sum, 1);
    Us[i][0] -> transform32 (tmp, -1);
    Veclib::vvtvp (nTot32, tmp, 1, tmp, 1, sum, 1, sum, 1);
  }

  Veclib::vsqrt (nTot32, sum, 1, sum, 1);
  EV -> transform32 (sum, +1);
  EV -> Smagorinsky ();

#endif
}

