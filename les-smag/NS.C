///////////////////////////////////////////////////////////////////////////////
// NS.C: Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
//             This version has Smagorinsky eddy-viscosity LES.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <NS.h>

#ifdef __DECCXX
  #pragma define_template roll<AuxField*>
#endif

typedef ModalMatrixSystem ModeSys;
static  int               DIM;

static void strainRate    (const Domain*,AuxField***,AuxField***);
static void eddyViscosity (AuxField***,AuxField***,AuxField*);
static void nonLinear     (Domain*,AuxField***,AuxField***,AuxField*,Vector&);
static void waveProp      (Domain*,const AuxField***,const AuxField***);
static void setPForce     (const Domain*,const AuxField***,AuxField***);
static void project       (const Domain*,AuxField***,AuxField***);
static void setUForce     (const Domain*,AuxField***);

static ModeSys** preSolve (const Domain*);
static void      Solve    (Field*, AuxField*, ModeSys*, const int, const int);


void NavierStokes (Domain*   D,
		   Analyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  DIM = D -> nField() - 1;

  int        i, j;
  const real dt       =       Femlib::value ("D_T");
  const int  nOrder   = (int) Femlib::value ("N_TIME");
  const int  nStep    = (int) Femlib::value ("N_STEP");
  const int  nZ       = (int) Femlib::value ("N_Z");

  Field*     Pressure = D -> u[DIM];
  ModeSys**  MMS      = preSolve (D);

  // -- Set up multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [DIM];
  AuxField*** Uf = new AuxField** [DIM];

  for (i = 0; i < DIM; i++) {
    Us[i] = new AuxField* [nOrder];
    Uf[i] = new AuxField* [nOrder];
    for (j = 0; j < nOrder; j++) {
      Us[i][j] = new AuxField (D -> Esys, nZ);
      Uf[i][j] = new AuxField (D -> Esys, nZ);
    }
  }

  // -- Set up multi-level storage for pressure BCS.

  PBCmgr::build (Pressure);

  // -- Set up spatially-constant forcing terms.

  Vector ff = {0.0, 0.0, 0.0};

  ff.x = Femlib::value ("FFX");
  ff.y = Femlib::value ("FFY");
  ff.z = Femlib::value ("FFZ");

  // -- Set up eddy viscosity storage.

  AuxField* varVis = new AuxField (D -> Esys, nZ);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Compute spatially-varying kinematic viscosity.

    strainRate    (D, Us, Uf);
    eddyViscosity (Us, Uf, varVis);
    *varVis -= Femlib::value ("KINVIS");

    // -- Unconstrained forcing substep.

    nonLinear (D, Us, Uf, varVis, ff);
    waveProp  (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      (const AuxField***) Us, (const AuxField***) Uf, 1);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce (D, (const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Viscous correction substep.

    setUForce (D, Uf);
    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    }

    // -- Process results of this step.

    A -> analyse (Us);
  }
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
//           / dU1/dx1  1/2(dU1/dx2 + dU2/dx1)  1/2(dU1/dx3 + dU3/dx1) \
//   Sij =   |    .              dU2/dx2        1/2(dU2/dx3 + dU3/dx2) |  (3D)
//           \    .                 .                     dU3/dx3      /
// ---------------------------------------------------------------------------
{
  int i, j;

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


static void eddyViscosity (AuxField*** Us ,
			   AuxField*** Uf ,
			   AuxField*   nuT)
// ---------------------------------------------------------------------------
// On entry the first-level areas of Us & Uf contain the components of the
// strain-rate tensor Sij.  Construct in nuT the Smagorinsky eddy-viscosity
// field (Cs \Delta)^2 |S| where (3D, symmetry)
//
// |S| = sqrt [(S11)^2 + (S22)^2 + (S33)^2 + 2(S12)^2 + 2(S13)^2 + 2(S23)^2].
// ---------------------------------------------------------------------------
{
  int          i, j;
  const int    nZ     = Geometry::nZ();
  const int    nP     = Geometry::planeSize();
  const int    nZ32   = (3 * nZ) >> 1;
  const int    nTot32 = nZ32 * nP;

  vector<real> work (2 * nTot32);
  real*        tmp    = work();
  real*        sum    = tmp + nTot32;

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
  nuT -> transform32 (sum, +1);
  nuT -> Smagorinsky ();
}


static void nonLinear (Domain*     D ,
		       AuxField*** Us,
		       AuxField*** Uf,
		       AuxField*   vV,
		       Vector&     ff)
// ---------------------------------------------------------------------------
// Compute nonlinear + forcing terms in Navier--Stokes equations
//                       N(u) + 2*div(vV*S) + ff.
//
// Here N(u) represents the nonlinear advection terms in the N--S equations
// transposed to the RHS, div(vV*S) is the divergence of the non-constant
// component of the eddy-viscosity SGS terms and ff is a vector of body
// force per unit mass.
//
// On entry, D contains the old velocity (and pressure) fields, the
// lowest levels of Us & Uf contain the components of the old strain-rate
// tensor and vV contains the spatially-varying viscosity.  On exit, the
// velocity field storage areas of D are free, the zeroth level of Us contains
// the old velocities and the zeroth level of Uf contains the most recent
// explcit forcing terms.
// Velocity field data areas of D and first level of Us are swapped, then
// the next stage of nonlinear forcing terms N(u) - a are computed from
// velocity fields and left in the first level of Uf.
//
// Nonlinear terms N(u) are computed in skew-symmetric form (Zang 1991)
//                 ~ ~
//                   N = -0.5 ( u . grad u + div uu )
//                   ~          ~        ~       ~~
//
// i.e., in Cartesian component form
//
//              N  = -0.5 ( u  d(u ) / dx  + d(u u ) / dx ).
//               i           j    i      j      i j      j
//
// Data are transformed to physical space for most of the operations, with
// the Fourier transform extended using zero padding for dealiasing.  For
// gradients in the Fourier direction however, the data must be transferred
// back to Fourier space.
// ---------------------------------------------------------------------------
{
  int          i, j;
  vector<real> F (DIM);

  F[0] = ff.x; F[1] = ff.y; if (DIM == 3) F[2] = ff.z;

  const int         nZ     = Geometry::nZ();
  const int         nP     = Geometry::planeSize();
  const int         nTot   = nZ * nP;
  const int         nZ32   = (3 * nZ) >> 1;
  const int         nTot32 = nZ32 * nP;
  const real        EPS    = (sizeof(real) == sizeof(double)) ? EPSDP : EPSSP;
  vector<real>      work ((2 * DIM + 1) * nTot32);
  vector<real*>     u32 (DIM);
  vector<real*>     n32 (DIM);
  vector<AuxField*> U   (DIM);
  vector<AuxField*> N   (DIM);
  Field*            master = D -> u[0];
  real*             tmp    = work() + 2 * DIM * nTot32;

  // -- Divergence of SGS.

  for (i = 0; i < DIM; i++) {
    u32[i] = work() +  i        * nTot32;
    n32[i] = work() + (i + DIM) * nTot32;
  }

  vV -> transform32 (u32[0], -1);

  for (i = 0; i < DIM; i++) {

    // -- Diagonal.

    Us[i][0] -> transform32 (n32[i], -1);

    // -- Off-diagonal.

    for (j = i + 1; j < DIM; j++) {
      Uf[i + j - 1][0] -> transform32 (tmp, -1);
      Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);
      Veclib::vadd (nTot32, tmp, 1, n32[j], 1, n32[j], 1);
    }

    // -- Convert rates of strain into stresses, account for factor -0.5 below.

    Veclib::svvtt (nTot32, -4.0, u32[0], 1, n32[i], 1, n32[1], 1);

    if (i == 2) {
      Femlib::DFTr (n32[2], nZ32, nP, +1);
      Veclib::zero (nTot32 - nTot, n32[2] + nTot, 1);
      master -> gradient (nZ, n32[2], 2);
      Femlib::DFTr (n32[2], nZ32, nP, -1);
    } else {
      master -> gradient (nZ32, n32[i], i);
    }
  }

  // -- Nonlinear terms.

  for (i = 0; i < DIM; i++) {
    AuxField::swapData  (D -> u[i], Us[i][0]);
    U[i] = Us[i][0];
    N[i] = Uf[i][0];
    U[i] -> transform32 (u32[i], -1);
  }
  
  for (i = 0; i < DIM; i++) {
    for (j = 0; j < DIM; j++) {
      
      // -- Perform n_i = u_j d(u_i) / dx_j.

      Veclib::copy (nTot32, u32[i], 1, tmp,  1);
      if (j == 2) {
	Femlib::DFTr (tmp, nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, tmp, j);
	Femlib::DFTr (tmp, nZ32, nP, -1);
      } else {
	master -> gradient (nZ32, tmp, j);
      }
      Veclib::vvtvp (nTot32, u32[j], 1, tmp,  1, n32[i], 1, n32[i], 1);

      // -- Perform n_i += d(u_i u_j) / dx_j.

      Veclib::vmul  (nTot32, u32[i], 1, u32[j], 1, tmp,  1);
      if (j == 2) {
	Femlib::DFTr (tmp, nZ32, nP, +1);
	Veclib::zero (nTot32 - nTot, tmp + nTot, 1);
	master -> gradient (nZ, tmp, j);
	Femlib::DFTr (tmp, nZ32, nP, -1);
      } else {
	master -> gradient (nZ32, tmp, j);
      }
      Veclib::vadd (nTot32, tmp, 1, n32[i], 1, n32[i], 1);

    }
    N[i]   -> transform32 (n32[i], +1);
    master -> smooth (N[i]);

    if (fabs (F[i]) > EPS) N[i] -> addToPlane (0, -2.0*F[i]);
    *N[i] *= -0.5;
  }
}


static void waveProp (Domain*           D ,
		      const AuxField*** Us,
		      const AuxField*** Uf)
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// On entry, the most recent velocity fields are in Us, and the most
// recent nonlinear terms in Uf.  The intermediate velocity field u^ is
// computed and left in D's velocity areas. 
// ---------------------------------------------------------------------------
{
  int               i, q;
  vector<AuxField*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  int  Je = (int) Femlib::value ("N_TIME");
  Je = min (D -> step, Je);

  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i < DIM; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[i][q]);
      H[i] -> axpy ( beta [q]    , *Uf[i][q]);
    }
}


static void setPForce (const Domain*     D ,
		       const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
//
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf as a result of gradient operation.
// ---------------------------------------------------------------------------
{
  int         i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {
   *Uf[i][0] = *Us[i][0];
    Uf[i][0] -> gradient (i);
  }
  
  for (i = 1; i < DIM; i++) *Uf[0][0] += *Uf[i][0];

  *Uf[0][0] /= dt;
}


static void project (const Domain* D ,
		     AuxField***   Us,
		     AuxField***   Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in lowest level of Us.  Constrain velocity field:
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in lowest level of Uf.
//
// A frame-swapping operation takes place in the first time level of the
// Fourier direction of Uf.  This returns to original place the swapping done
// by setPForce.
// ---------------------------------------------------------------------------
{
  int        i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM];
    Uf[i][0] -> gradient (i);
  
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);
  }
}


static void setUForce (const Domain* D ,
		       AuxField***   Uf)
// ---------------------------------------------------------------------------
// On entry, intermediate velocity storage u^^ is in lowest levels of Uf.
// Multiply by -1.0 / (D_T * KINVIS) to create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  int        i;
  const real alpha = -1.0 / Femlib::value ("D_T*KINVIS");

  for (i = 0; i < DIM; i++) *Uf[i][0] *= alpha;
}


static ModeSys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystems for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSystem pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocity components,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  char                 name;
  int                  i;
  const int            nSys   = D -> Nsys.getSize();
  const int            nZ     = Geometry::nZ();
  const int            nModes = (nZ + 1) >> 1;
  const int            base   = 0;
  const int            itLev  = (int) Femlib::value ("ITERATIVE");
  const int            nOrder = (int) Femlib::value ("N_TIME");
  const real           beta   = Femlib::value ("BETA");
  ModeSys**            M      = new ModeSys* [DIM + 1];
  vector<Element*>&    E      = ((Domain*) D) -> Esys;
  const NumberSystem** N      = new const NumberSystem* [3];

  // -- Velocity systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real   lambda2 = alpha[0] / Femlib::value ("KINVIS*D_T");

    for (i = 0; i < DIM; i++) {
      name = D -> u[i] -> name();
      D -> setNumber (name, N);
      M[i] = new ModalMatrixSystem (lambda2, beta, name, base, nModes, E, N);
    }
  } else
    for (i = 0; i < DIM; i++) M[i] = 0;

  // -- Pressure system.

  if (itLev < 2) {
    name = D -> u[DIM] -> name();
    D -> setNumber (name, N);
    M[DIM] = new ModalMatrixSystem (0.0, beta, name, base, nModes, E, N);
  } else
    M[DIM] = 0;

  return M;
}


static void Solve (Field*    U     ,
		   AuxField* Force ,
		   ModeSys*  M     ,
		   const int step  ,
		   const int nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  char routine[] = "Solve";

  const int  iterative = M == 0;
  const char name      = U -> name();
  const int  velocity  = name == 'u' || name == 'v' || name == 'w';
  const int  pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    if   (iterative) U -> solve (Force, 0.0);
    else             U -> solve (Force, M);
    return;
  }

  if (velocity) {
    if (iterative || step < nOrder) {
      const int    Je = min (step, nOrder);
      vector<real> alpha (Je + 1);
      Integration::StifflyStable (Je, alpha());
      const real   lambda2 = alpha[0] / Femlib::value ("D_T*KINVIS");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
