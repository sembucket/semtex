///////////////////////////////////////////////////////////////////////////////
// transport.C: Unsteady advection--diffusion solver, Cartesian cords.
//
// Copyright (C) 1999 Hugh Blackburn
//
// Here there is only a single velocity component, and we solve
// for the transport of the scalar field 'c'.  The single component
// of fluid velocity, 'w', is in the Fourier/z direction; this is a 
// linear transport problem:
//
//   dc     dc
//   -- + w -- = D div grad c.
//   dt     dz
//
// Time split is discussed in: Karniadakis, Israeli & Orszag,
// "High-order splitting methods for the incompressible Navier--Stokes
// equations", JCP 9(2). 1991.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "Sem.h"

typedef ModalMatrixSys Msys;

static void  advection (AuxField*, AuxField*, const real*);
static void  extrap    (AuxField*, AuxField**, AuxField**,
			const int, const real);
static Msys* preSolve  (const Domain*, const real);
static void  Solve     (Field*, Field*, const AuxField*, const AuxField*,
		        AuxField*, AuxField*, Msys*, Msys*,
			MixPatch*, const int);

static int np, nz, nel1, nel2;
static Geometry::CoordSys space;

#define SET1 Geometry::set (np, nz, nel1, space)
#define SET2 Geometry::set (np, nz, nel2, space)


void transport (Domain*     d1   ,
		Domain*     d2   ,
		Analyser*   a1   ,
		Analyser*   a2   ,
		const real* gas  ,
		MixPatch*   patch)
// ---------------------------------------------------------------------------
// On entry, d1 & d2 contain storage for scalar Field 'c', pre-transformed in
// the z direction.  Input w contains the z-component velocity field.
// The first domain is assumed to be the mobile phase, and the matching
// gas data supplies the velocity field for advection.
//
// Us is multi-level auxillary Field storage for scalar history and 
// Uf is multi-level auxillary Field storage for advection history terms.
// ---------------------------------------------------------------------------
{
  const real dt    = Femlib::value ("D_T");
  const int  Je    = (int) Femlib::value ("N_TIME");
  const int  nStep = (int) Femlib::value ("N_STEP");
  const real DC1   = Femlib::value ("DC_GAS");
  const real DC2   = Femlib::value ("DC_FIX");
  int        i, order;

  // -- Set file-scope global variables.

  np    = Geometry::nP();
  nz    = Geometry::nZ();
  nel1  = d1 -> elmt.getSize();
  nel2  = d2 -> elmt.getSize();
  space = Geometry::Cartesian;

  // -- Create global matrix systems.

  SET1; Msys* MMS1 = preSolve (d1, DC1);
  SET2; Msys* MMS2 = preSolve (d2, DC2);

  // -- Coupled diffusion substep work areas.

  SET1; AuxField* tmp1 = new AuxField (d1 -> elmt, nz);
  SET2; AuxField* tmp2 = new AuxField (d2 -> elmt, nz);

  // -- Create & initialise multi-level storage for old fields and advection.

  AuxField** Us1 = new AuxField* [(size_t) Je];
  AuxField** Uf1 = new AuxField* [(size_t) Je];
  AuxField** Us2 = new AuxField* [(size_t) Je];

  for (i = 0; i < Je; i++) {
    SET1;
    *(Us1[i] = new AuxField (d1 -> elmt, nz)) = 0.0;
    *(Uf1[i] = new AuxField (d1 -> elmt, nz)) = 0.0;
    SET2;
    *(Us2[i] = new AuxField (d2 -> elmt, nz)) = 0.0;
  }

  // -- Dump startup analysis information.

  SET1; a1 -> analyse();
  SET2; a2 -> analyse();

  // -- Timestepping loop.

  while (d1 -> step < nStep) {
 
    d1 -> step += 1; d1 -> time += dt;
    d2 -> step += 1; d2 -> time += dt;

    Femlib::value ("t", d1 -> time);
    
    order = clamp (d1 -> step, 0, Je);

    // -- Explicit update for advection and old time levels on d1.

    SET1;
    AuxField::swapData (d1 -> u[0], Us1[0]);
    advection (Us1[0],  Uf1[0], gas);
    extrap    (d1 -> u[0], Us1, Uf1, order, DC1);
    roll      (Us1, Je);
    roll      (Uf1, Je);
    AuxField::swapData (d1 -> u[0], Us1[0]);

    // -- Explicit update for old time levels on d2 (no advection).

    SET2;
    AuxField::swapData (d2 -> u[0], Us2[0]);
    extrap (d2 -> u[0], Us2, 0, order, DC2);
    roll   (Us2, Je);
    AuxField::swapData (d2 -> u[0], Us2[0]);

    // -- Diffusion substep, coupling two domains over patch.

    Solve (d1 -> u[0], d2 -> u[0],
	   Us1[0],     Us2[0],
	   tmp1,       tmp2,
	   MMS1,       MMS2,
	   patch,      order);

    // -- Print diagnostic information.

    SET1; a1 -> analyse();
    SET2; a2 -> analyse();
  }
}


static void advection (AuxField*   C,
		       AuxField*   A,
		       const real* w)
// ---------------------------------------------------------------------------
// Compute the Fourier transform of the advection terms
// 
//       ^
//      dc
//    w --  = i beta k w c
//      dz                 k
//         k
//
// and leave in A.  This is a linear operation, no dealiasing needed.
// Old time level of scalar is passed in as C.
// ---------------------------------------------------------------------------
{
  static const int nzp = Geometry::nZProc();
  int              i;

  for (i = 0; i < nzp; i++)
    A -> setPlane (i, w);

  (A -> times (*A, *C)) . gradient (2);
}


static void extrap (AuxField*  C  ,
		    AuxField** Us ,
		    AuxField** Uf ,
		    const int  step,
		    const real DiffusionCoeff)
// ---------------------------------------------------------------------------
// Compute the RHS explicit approximation of forcing terms.
//
// On entry, the most recent scalar fields are in Us, and the most
// recent advection terms in Uf.  The forcing field is
// computed and left in C's storage areas.
// ---------------------------------------------------------------------------
{
  static real  dt = Femlib::value ("D_T");
  int          q;
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> beta  (Integration::OrderMax);

  Integration::StifflyStable (step, alpha());
  Integration::Extrapolation (step, beta ());
  Blas::scal (step + 1, 1.0 / (dt * DiffusionCoeff), alpha(), 1);
  Blas::scal (step,     1.0 /       DiffusionCoeff,  beta(),  1);

  *C = 0.0;

          for (q = 0; q < step; q++) C -> axpy (alpha[q + 1], *Us[q]);

  if (Uf) for (q = 0; q < step; q++) C -> axpy (beta [q]    , *Uf[q]);

}


static Msys* preSolve (const Domain* D             ,
		       const real    DiffusionCoeff)
// ---------------------------------------------------------------------------
// Set up ModalMatrixSystem for D -> u[0].  If iterative solution
// (ITERATIVE == 1) is selected, return 0.
// ---------------------------------------------------------------------------
{
  const int               nmodes = Geometry::nModeProc();
  const int               base   = Geometry::baseMode();
  const int               itLev  = (int) Femlib::value ("ITERATIVE");
  const int               nOrder = (int) Femlib::value ("N_TIME");
  const real              beta   = Femlib::value ("BETA");
  const real              dt     = Femlib::value ("D_T");
  Msys*                   M;
  const vector<Element*>& E      = ((Domain*) D) -> elmt;

  if (itLev < 1) {
    vector<real> alpha (Integration::OrderMax + 1);
    Integration::StifflyStable (nOrder, alpha());
    const real lambda2 = alpha[0] / (dt * DiffusionCoeff);

    M = new Msys (lambda2, beta, base, nmodes, E, D -> b[0]);

  } else M = 0;

  return M;
}


static void Solve (Field*          c1   ,
		   Field*          c2   ,
		   const AuxField* f1   ,
		   const AuxField* f2   ,
		   AuxField*       t1   ,
		   AuxField*       t2   ,
		   Msys*           M1   ,
		   Msys*           M2   ,
		   MixPatch*       patch,
		   const int       je   )
// ---------------------------------------------------------------------------
// Coupled Helmholtz solution for fields c1 & c2.  Termination is less
// efficient than it could be, since we continue to iterate until all
// modes have converged, instead of stopping each mode as it is done.
// ---------------------------------------------------------------------------
{
  static const int  verbose = (int) Femlib::value ("VERBOSE");
  static const int  JE      = (int) Femlib::value ("N_TIME");
  static const int  MAXITN  = (int) Femlib::value ("STEP_MAX");
  static const real TOLREL  = Femlib::value ("TOL_REL");
  static const real DCON1   = Femlib::value ("DC_GAS");
  static const real DCON2   = Femlib::value ("DC_FIX");
  static const real DT      = Femlib::value ("D_T");
  static const int  nmodes  = Geometry::nModeProc();

  const int    DIRECT = !(je < JE || M1 == 0);		     
  vector<real> alpha (Integration::OrderMax + 1);
  vector<real> L2new1 (nmodes), L2old1 (nmodes);
  vector<real> L2new2 (nmodes), L2old2 (nmodes);
  vector<int>  test1  (nmodes), test2  (nmodes);
  integer      i, k, converged;

  Integration::StifflyStable (je, alpha());
  const real lambda2_1 = alpha[0] / (DT * DCON1);
  const real lambda2_2 = alpha[0] / (DT * DCON2);

  L2old1 = 0.0;
  L2old2 = 0.0;

  for (converged = 0, i = 0; i < MAXITN && !converged; i++) {
    c1 -> setPatch (patch);
    c2 -> setPatch (patch);

    c1 -> getPatch (patch);
    c2 -> getPatch (patch);

    SET1; *t1 = *f1;
    SET2; *t2 = *f2;

    if (DIRECT) {

      SET1; c1 -> solve (t1, M1);
      SET2; c2 -> solve (t2, M2);

    } else {

      SET1; c1 -> solve (t1, lambda2_1);
      SET2; c2 -> solve (t2, lambda2_2);

    }

    SET1;
    for (k = 0; k < nmodes; k++) {
      L2new1[k] = c1 -> mode_L2 (k);
      test1 [k] = fabs ((L2new1[k] - L2old1[k])/L2new1[k]) < TOLREL;
      L2old1[k] = L2new1[k];
      VERBOSE cout << "mode: " << k << ", energy1 = " << L2old1[k];
    }

    SET2;
    for (k = 0; k < nmodes; k++) {
      L2new2[k] = c2 -> mode_L2 (k);
      test2 [k] = fabs ((L2new2[k] - L2old2[k])/L2new2[k]) < TOLREL;
      L2old2[k] = L2new2[k];
      VERBOSE cout  << "mode: " << k << ", energy2 = " << L2old2[k];
    }
    VERBOSE cout << endl;

    for (converged = 1, k = 0; converged && k < nmodes; k++)
      converged = converged && test1[k] && test2[k];
  }
  
  cout << i << " BC iterations" << endl;
}
