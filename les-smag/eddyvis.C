///////////////////////////////////////////////////////////////////////////////
// eddyvis.C: calculate eddy viscosity for LES.  Cartesian.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <les.h>

static integer DIM, NZ;

static void strainRate  (const Domain*,AuxField***,AuxField***);
static void Smagorinsky (AuxField***, AuxField***, AuxField*);


void eddyViscosity (const Domain* D ,
		    AuxField***   Us,
		    AuxField***   Uf,
		    AuxField*     EV)
// ---------------------------------------------------------------------------
// Compute eddy viscosity field E from resolved-scale velocities.
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
  NZ  = Geometry::nZ();

  strainRate  (D,  Us, Uf);
  Smagorinsky (Us, Uf, EV);

  *EV -= Femlib::value ("REFVIS");
}


static void strainRate (const Domain* D ,
			AuxField***   Us,
			AuxField***   Uf)
// ---------------------------------------------------------------------------
// On entry D contains the velocity fields Ui and the first-level areas of
// Us and Uf are free.  Construct the symmetric strain-rate tensor terms,
// leave the diagonal terms Sii (unsummed) in Us and the off-diagonal terms
// Sij (i != j) in Uf.
//           
//             / du/dx  1/2(du/dy + dv/dx)  1/2(du/dz + dw/dx) \
//       S =   |    .           dv/dy       1/2(dv/dz + dw/dy) |
//       ~     \    .             .                 dw/dz      /
// ---------------------------------------------------------------------------
{
  integer i, j;

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


static void Smagorinsky (AuxField*** Us,
			 AuxField*** Uf,
			 AuxField*   EV)
// ---------------------------------------------------------------------------
// On entry the first-level areas of Us & Uf contain the components of the
// strain-rate tensor S.  On exit EV contains the Smagorinsky eddy viscosity
// field \nu_S = (Cs \Delta)^2 |S|, where
// |S| = sqrt [(S11)^2 + (S22)^2 + (S33)^2 + 2(S12)^2 + 2(S13)^2 + 2(S23)^2].
//
// As noted in the header to NS.C, EV = -KINVIS for debugging (NB: Fourier!).
// ---------------------------------------------------------------------------
{
#if defined(DEBUG)
  (*EV = 0.0) . addToPlane (0, -0.5 * Femlib::value ("KINVIS"));

#else
  integer       i, j;
  const integer nP     = Geometry::planeSize();
  const integer nZ32   = (3 * NZ) >> 1;
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

