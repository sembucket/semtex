///////////////////////////////////////////////////////////////////////////////
// integrate.C: integrate unsteady momentum equations forward in time.
//
// Copyright (c) Hugh Blackburn 1999--2000
//
// This version incorporates eddy-viscosity SGS stress, with the
// dynamic algorithm.  After computation of the spatially-variable
// eddy viscoisty, it is handled within the framework of the
// implicit-explicit time integration by splitting the SGS stress into
// two parts: one which can be dealt with implicitly as diffusion of
// momentum with a spatially-constant kinematic viscosity KINVIS, the
// other dealt with explicitly as the divergence of a stress.  This
// way the velocity BCs are still set correctly in the viscous
// substep.
//
// For testing, define DEBUG during compilation.  This will set the
// kinematic viscosity equal to twice the input value and formulate
// the SGS stresses as those supplied by the product of the strain
// rate tensor with a spatially-constant NEGATIVE viscosity -KINVIS:
// the results should be the very similar over a short time period to
// a DNS with total viscosity KINVIS.
//
// No attempt to incorporate dealiasing in z-direction.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"

typedef ModalMatrixSys  Msys;
static  integer         DIM, CYL, C3D;
   
static void   nonLinear (Domain*, AuxField***, AuxField***,
			 AuxField*, vector<real>&);
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField***, AuxField***);
static void   project   (const Domain*, AuxField***, AuxField***);
static Msys** preSolve  (const Domain*);
static void   Solve     (Field*, AuxField*, Msys*,const integer,const integer);


void integrate (Domain*      D,
		const real*  FourierMask,
		const real*  LegendreMask,
		LESAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' ('w') and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  DIM = 3;			// -- Guaranteed if we got this far.
  CYL = Geometry::system() == Geometry::Cylindrical;

  integer       i, j, k;
  const real    dt     =           Femlib::value ("D_T"   );
  const integer nOrder = (integer) Femlib::value ("N_TIME");
  const integer nStep  = (integer) Femlib::value ("N_STEP");
  const integer nZ     = Geometry::nZProc();
  const integer nTot   = Geometry::nTotProc();

  // -- Create global matrix systems: rename viscosities beforehand.

#if defined(DEBUG)
  Femlib::value ("REFVIS", Femlib::value ("2.0 * KINVIS"));
#endif

  if (Femlib::value ("REFVIS") > 0.0) {
    real kinVis = Femlib::value ("REFVIS");
    real refVis = Femlib::value ("KINVIS");
    Femlib::value ("KINVIS", kinVis);
    Femlib::value ("REFVIS", refVis);
  } 

  Msys** MMS = preSolve (D);

  // -- Create extra storage needed for computation of SGSS, nonlinear terms.
  //    First DIM*nOrder*2 of these are used for Us & Uf.
  
  matrix<real> Ut (25, Geometry::nTotProc());

  // -- Create & initialize multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [(size_t) DIM];
  AuxField*** Uf = new AuxField** [(size_t) DIM];

  for (k = 0, i = 0; i < DIM; i++) {
    Us[i] = new AuxField* [(size_t) nOrder];
    Uf[i] = new AuxField* [(size_t) nOrder];
    for (j = 0; j < nOrder; j++) {
      *(Us[i][j] = new AuxField (Ut(k++), nZ, D -> elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (Ut(k++), nZ, D -> elmt)) = 0.0;
    }
  }

  // -- Create multi-level storage for pressure BCs.

  Field* Pressure = D -> u[DIM];
  PBCmgr::build (Pressure);

  // -- Create spatially-constant forcing terms.

  vector<real> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], +1);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);

    // -- Compute nonlinear terms and strain-rate tensors.

    nonLinear (D, Us, Uf, Ut);

    // -- Compute modified eddy viscosity and associated SGGS.

    SGSS (D, Us, Uf, Ut);

    // -- Add divergence of SGSS to nonlinear terms.

    turbModel (D, Us, Uf, Ut);

    // -- Fourier Transform velocity fields & nonlinear terms, add forcing.

    transform (D, Uf, Ut, ff);

    // -- Take unconstrained forcing substep.

    waveProp (D, (const AuxField***) Us, (const AuxField***) Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure,
		      (const AuxField***) Us, (const AuxField***) Uf, 1);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < DIM; i++) {
      AuxField::swapData (D -> u[i], Us[i][0]);
      roll (Uf[i], nOrder);
    }
    setPForce ((const AuxField***) Us, Uf);
    Solve     (Pressure, Uf[0][0], MMS[DIM], D -> step, nOrder);
    project   (D, Us, Uf);

    // -- Update multilevel velocity storage.

    for (i = 0; i < DIM; i++) {
      *Us[i][0] = *D -> u[i];
      roll (Us[i], nOrder);
    }

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [1][0], Uf [2][0], +1);
      AuxField::couple (D -> u[1], D -> u[2], +1);
    }
    for (i = 0; i < DIM; i++)
      Solve (D -> u[i], Uf[i][0], MMS[i], D -> step, nOrder);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], -1);

    // -- Process results of this step.

    A -> analyse (Us);
  }

  // -- Dump ratio eddy/molecular viscosity to file visco.fld.

  ofstream          evfl;
  vector<AuxField*> visco (1);

  visco[0] = EV;

  ROOTONLY {
    evfl.open ("visco.fld", ios::out);
    EV -> addToPlane (0, Femlib::value ("KINVIS"));
  }

  (*EV /= Femlib::value ("REFVIS")) . transform (-1);

  writeField (evfl, D -> name, D -> step, D -> time, visco);

  ROOTONLY evfl.close();
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
  integer           i, q;
  vector<AuxField*> H (DIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < DIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, (integer) Femlib::value ("N_TIME"));
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


static void setPForce (const AuxField*** Us,
		       AuxField***       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in first time-level of Us.
// Create div u^ / D_T in the first dimension, first level storage
// of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < DIM; i++) {
   *Uf[i][0] = *Us[i][0];
    Uf[i][0] -> gradient (i);
  }

  if (C3D)
    Uf[2][0] -> divR();

  for (i = 1; i < DIM; i++) *Uf[0][0] += *Uf[i][0];

  if (CYL) {
    *Uf[1][0]  = *Us[1][0];
    Uf[1][0] ->  divR();
    *Uf[0][0] += *Uf[1][0];
  }

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
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to create
// forcing for viscous step.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < DIM; i++) {

   *Uf[i][0] = *D -> u[DIM];
    Uf[i][0] -> gradient (i);

    if (C3D) Uf[2][0] -> divR();

    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Field::swapData (Us[i][0], Uf[i][0]);

    *Uf[i][0] *= alpha;
  }
}


static Msys** preSolve (const Domain* D)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSyss for each Field of D.  If iterative solution
// is selected for any Field, the corresponding ModalMatrixSys pointer
// is set to zero.
//
// ITERATIVE == 1 selects iterative solvers for velocity components,
// ITERATIVE == 2 adds iterative solver for pressure as well.
// ---------------------------------------------------------------------------
{
  const integer           nmodes = Geometry::nModeProc();
  const integer           base   = Geometry::baseMode();
  const integer           itLev  = (integer) Femlib::value ("ITERATIVE");
  const integer           nOrder = (integer) Femlib::value ("N_TIME");
  const real              beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [(size_t) (DIM + 1)];
  integer                 i;

  // -- Velocity systems.

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

    for (i = 0; i < DIM; i++)
      M[i] = new Msys (lambda2, beta, base, nmodes, E, D -> b[i]);
  } else
    for (i = 0; i < DIM; i++)
      M[i] = 0;

  // -- Pressure system.

  if (itLev < 2)
    M[DIM] = new Msys (0.0, beta, base, nmodes, E, D -> b[DIM]);
  else
    M[DIM] = 0;

  return M;
}


static void Solve (Field*        U     ,
		   AuxField*     Force ,
		   Msys*         M     ,
		   const integer step  ,
		   const integer nOrder)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for U, using F as a forcing Field.  Iterative
// or direct solver selected on basis of field type, step, time order
// and command-line arguments.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Solve";
  const integer iterative = M == 0;
  const char    name      = U -> name();
  const integer velocity  = name == 'u' || name == 'v' || name == 'w';
  const integer pressure  = name == 'p';

  if (!(velocity || pressure))
    message (routine, "input field type not recognized", ERROR);

  if (pressure) {
    if   (iterative) U -> solve (Force, 0.0);
    else             U -> solve (Force, M);
    return;
  }

  if (velocity) {
    if (iterative || step < nOrder) {
      const integer    Je = min (step, nOrder);
      vector<real> alpha (Je + 1);
      Integration::StifflyStable (Je, alpha());
      const real   lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
      
      U -> solve (Force, lambda2);

    } else
      U -> solve (Force, M);

    return;
  }
}
