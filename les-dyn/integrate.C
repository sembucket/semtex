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
// For testing, define NOMODEL during compilation.  This will set the
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
static  integer         NORD, NDIM, CYL, C3D;
   
static void   waveProp  (Domain*, const AuxField***, const AuxField***);
static void   setPForce (const AuxField**, AuxField**);
static void   project   (const Domain*, AuxField**, AuxField**);
static Msys** preSolve  (const Domain*);
static void   Solve     (Domain*, const integer, AuxField*, Msys*);
static void   pushdown  (AuxField***, const integer, const integer);

void integrate (Domain*      D,
		LESAnalyser* A)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Fields 'u', 'v' 'w' and
// constraint Field 'p'.
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  NDIM = 3;			// -- Guaranteed if we got this far.
  NORD = (integer) Femlib::value ("N_TIME");
  CYL  = Geometry::system() == Geometry::Cylindrical;

  integer       i, j;
  const real    dt    =           Femlib::value ("D_T"   );
  const integer nStep = (integer) Femlib::value ("N_STEP");
  const integer nZ    = Geometry::nZProc();
  const integer nTot  = Geometry::nTotProc();

  // -- Create global matrix systems: rename viscosities beforehand.

#if !defined (NOMODEL)
  if (Femlib::value ("REFVIS") > 0.0) {
    real kinVis = Femlib::value ("REFVIS");
    real refVis = Femlib::value ("KINVIS");
    Femlib::value ("KINVIS", kinVis);
    Femlib::value ("REFVIS", refVis);
  } 
#endif

  Msys** MMS = preSolve (D);

  // -- Create extra storage needed for computation of SGSS, nonlinear
  //    terms.  Last 2*NDIM*NORD of these are used for Us & Uf, first 17
  //    are used for SGSS modelling work.
  
  matrix<real> Ut (17 + 2*NDIM*NORD, Geometry::nTotProc());

  // -- Create & initialise multi-level storage for velocities and forcing.

  AuxField*** Us = new AuxField** [(size_t) 2 * NORD];
  AuxField*** Uf = Us + NORD;

  for (i = 0; i < NORD; i++) {
    Us[i] = new AuxField* [(size_t) 2 * NDIM];
    Uf[i] = Us[i] + NDIM;
    for (j = 0; j < NDIM; j++) {
      *(Us[i][j] = new AuxField (Ut(17 +    2*i *NDIM+j), nZ, D->elmt)) = 0.0;
      *(Uf[i][j] = new AuxField (Ut(17 + (1+2*i)*NDIM+j), nZ, D->elmt)) = 0.0;
    }
  }

  // -- Create multi-level storage for pressure BCs.

  Field* Pressure = D -> u[NDIM];
  PBCmgr::build (Pressure);

  // -- Create spatially-constant forcing terms.

  vector<real> ff (3);

  ff[0] = Femlib::value ("FFX");
  ff[1] = Femlib::value ("FFY");
  ff[2] = Femlib::value ("FFZ");

  // -- Apply coupling to radial & azimuthal velocity BCs.

  if (C3D) Field::coupleBCs (D -> u[1], D -> u[2], FORWARD);

  // -- Timestepping loop.

  while (D -> step < nStep) {
 
    D -> step += 1; 
    D -> time += dt;
    Femlib::value ("t", D -> time);
#if 1

    for (i = 0; i < NDIM; i++) {
      lowpass (D -> udat[i]);
    }

#else

    // -- Compute nonlinear terms + divergence(SGSS) + body forces.

    nonLinear (D, Ut, ff);

    // -- Unconstrained forcing substep.

    waveProp (D, (const AuxField***)Us, (const AuxField***)Uf);

    // -- Pressure projection substep.

    PBCmgr::maintain (D -> step, Pressure, 
		      (const AuxField**)Us[0],
		      (const AuxField**)Uf[0]);
    Pressure -> evaluateBoundaries (D -> step);
    for (i = 0; i < NDIM; i++) {
      *Pressure  = *(const AuxField*) D -> u[i];
      *D -> u[i] = *(const AuxField*) Us [0][i];
      *Us [0][i] = *(const AuxField*) Pressure;
    }
    pushdown  (Uf, NORD, NDIM);
    setPForce ((const AuxField**)Us[0], Uf[0]);
    Solve     (D, NDIM, Uf[0][0], MMS[NDIM]);
    project   (D, Us[0], Uf[0]);

    // -- Update multilevel velocity storage.

    for (i = 0; i < NDIM; i++) *Us[0][i] = *D -> u[i];
    pushdown (Us, NORD, NDIM);

    // -- Viscous correction substep.

    if (C3D) {
      AuxField::couple (Uf [0][1], Uf [0][2], FORWARD);
      AuxField::couple (D -> u[1], D -> u[2], FORWARD);
    }
    for (i = 0; i < NDIM; i++) Solve (D, i, Uf[0][i], MMS[i]);
    if (C3D)
      AuxField::couple (D -> u[1], D -> u[2], INVERSE);

    // -- Process results of this step.
#endif
    A -> analyse (Us[0]);

  }
#if 0
  // -- Dump ratio eddy/molecular viscosity to file visco.fld.

  dynamic (D, Ut, 0);
  AuxField* EV = new AuxField (Ut(15), nZ, D->elmt, 'e');

  ofstream          evfl;
  vector<AuxField*> visco (1);
  visco[0] = EV;

  ROOTONLY evfl.open ("visco.fld", ios::out);

  *EV /= Femlib::value ("REFVIS");

  writeField (evfl, D -> name, D -> step, D -> time, visco);

  ROOTONLY evfl.close();
#endif
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
  vector<AuxField*> H (NDIM);	// -- Mnemonic for u^{Hat}.

  for (i = 0; i < NDIM; i++) {
     H[i] = D -> u[i];
    *H[i] = 0.0;
  }

  const integer Je = min (D -> step, NORD);
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);
  
  Integration::StifflyStable (Je, alpha());
  Integration::Extrapolation (Je, beta ());
  Blas::scal (Je, Femlib::value ("D_T"), beta(),  1);

  for (i = 0; i < NDIM; i++)
    for (q = 0; q < Je; q++) {
      H[i] -> axpy (-alpha[q + 1], *Us[q][i]);
      H[i] -> axpy ( beta [q]    , *Uf[q][i]);
    }
}


static void setPForce (const AuxField** Us,
		       AuxField**       Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in Us.  Create div u^ / D_T
// in the first dimension of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt = Femlib::value ("D_T");

  for (i = 0; i < NDIM; i++) (*Uf[i] = *Us[i]) . gradient (i);

  if (C3D) Uf[2] -> divR();

  for (i = 1; i < NDIM; i++) *Uf[0] += *Uf[i];

  if (CYL) *Uf[0] += (*Uf[1] = *Us[1]) . divR();

  *Uf[0] /= dt;
}


static void project (const Domain* D ,
		     AuxField**    Us,
		     AuxField**    Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in Us.  Constrain velocity field:
//                    u^^ = u^ - D_T * grad P;
// u^^ is left in Uf.
//
// After creation of u^^, it is scaled by -1.0 / (D_T * KINVIS) to
// create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  integer    i;
  const real dt    = Femlib::value ("D_T");
  const real alpha = -1.0 / Femlib::value ("D_T * KINVIS");

  for (i = 0; i < NDIM; i++) {

    (*Uf[i] = *D -> u[NDIM]) . gradient (i);

    if (C3D) Uf[2] -> divR();

    *Uf[i] *= -dt;
    *Uf[i] += *Us[i];
    *Uf[i] *= alpha;
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
  const real              beta   = Femlib::value ("BETA");
  const vector<Element*>& E      = D -> elmt;
  Msys**                  M      = new Msys* [(size_t) (NDIM + 1)];
  integer                 i;
  vector<real>            alpha (Integration::OrderMax + 1);
  Integration::StifflyStable (NORD, alpha());
  const real              lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");

  // -- Velocity systems.

  for (i = 0; i < NDIM; i++)
    M[i] = new Msys (lambda2, beta, base, nmodes, E, D -> b[i],
		     (itLev < 1) ? DIRECT : JACPCG);

  // -- Pressure system.

  M[NDIM] = new Msys (0.0,    beta, base, nmodes, E, D -> b[NDIM], 
		      (itLev < 2) ? DIRECT : JACPCG);

  return M;
}


static void Solve (Domain*       D,
		   const integer i,
		   AuxField*     F,
		   Msys*         M)
// ---------------------------------------------------------------------------
// Solve Helmholtz problem for D->u[i], using F as a forcing Field.
// Iterative or direct solver selected on basis of field type, step,
// time order and command-line arguments.
// ---------------------------------------------------------------------------
{
  const integer step = D -> step;

  if (i < NDIM && step < NORD) { // -- We need a temporary matrix system.
    const integer Je      = min (step, NORD);    
    const integer base    = Geometry::baseMode();
    const integer nmodes  = Geometry::nModeProc();
    vector<real>  alpha (Je + 1);
    Integration::StifflyStable (Je, alpha());
    const real    lambda2 = alpha[0] / Femlib::value ("D_T * KINVIS");
    const real    beta    = Femlib::value ("BETA");

    Msys* tmp = new Msys
      (lambda2, beta, base, nmodes, D -> elmt, D -> b[i], JACPCG);
    D -> u[i] -> solve (F, tmp);
    delete tmp;

  } else D -> u[i] -> solve (F, M);
}


static void pushdown (AuxField***   U ,
		      const integer nr,
		      const integer nc)
// ---------------------------------------------------------------------------
// Pushdown time-level stacks of velocity, forcing storage, without
// changing pointers (i.e. by copying data).
// ---------------------------------------------------------------------------
{
  integer i, j;

  for (i = nr - 1; i; i--)
    for (j = 0; j < nc; j++)
      *U[i][j] = *U[i-1][j];
}
