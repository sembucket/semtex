///////////////////////////////////////////////////////////////////////////////
// newton.C: solve a time-invariant NS flow by a Newton iteration.
//
// Copyright (C) 2002 Hugh Blackburn.
//
// The iterative method was devised by Laurette Tuckerman, see Refs
// [1] & [2], and is based on time-stepping with both linearised and
// full Navier--Stokes integrators.  The matrix-free linear systems
// are here solved using the Bi-Conjugate-Gradients-Stabilized
// algorithm, with code from the Templates package, Ref [3], or
// the Bi-Conjugate-Gradients-Squared algorithm, with code from the
// NSPCG package, Ref [4] (both codes available through netlib).
//
// USAGE
// -----
// newton [options] basesession pertsession
//   *session: specifies name of semtex session file.
//   options:
//   -h       ... print this message
//   -v       ... set verbose
//   -m <num> ... set maximum number of Newton iterations
//   -c <num> ... set convergence tolerance on Newton iteration
//   -n <num> ... set maximum number of BiCGS iterations per Newton step
//   -t <num> ... set BiCGS convergence tolerance
// 
// FILES
// -----
// A number of semtex files are required for each solution system
//   *session:     semtex session file
//   *session.num: computed automatically if not supplied
//   *session.rst: restart file (only relevant for base flow)
//
// Note that the base and perturbation session files *should be
// identical* except for boundary conditions (and possibly the USER
// section) -- typically the boundary conditions for the perturbation
// should be all zero. The TOKENS get read and the parser is reset,
// each time (in order), so at least the geometric descriptors in each
// session file have to match: it's safe to ensure the TOKENS match,
// but those in pertsession override those in basesession.
//
// REFERENCES
// ----------
// [1] C.K. Mamun & L.S. Tuckerman (1995), "Asymmetry and Hopf Bifurcation
//     in spherical Couette flow, Phys Fluids V7N1, 80--91.
// [2] L.S. Tuckerman & D. Barkley (2000), "Bifurcation analysis for
//     timesteppers", in Numerical Methods for Bifurcation Problems,
//     ed E. Doedel & L.S. Tuckerman, Springer. 453--466.
// [3] R. Barrett, M. Berry, T. Chan, J. Demmel, J. Donato, J. Dongarra, 
//     V. Eijkhout, R. Pozo, C. Romine & H. van der Vorst (1993),
//     "Templates for the Solution of Linear Systems: Building Blocks for
//     Iterative Methods", SIAM Publications.
// [4] T. C. Oppe, W. D. Joubert & D. R. Kincaid (1988), "NSPCG User's Guide",
//     Center for Numerical Analysis, University of Texas at Austin.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "newt.h"
#include <new.h>

static char prog[] = "newton";

// -- Linear system convergence control.

const real    INITOL  = 1.0e-3;
const real    SHRINK  = 0.67;
const integer NSTABLE = 3;

// -- Duplicates for base and perturbation systems.

static Domain* BaseDomain;
static Domain* PertDomain;
static FEML*   BaseFile;
static FEML*   PertFile;
static BCmgr*  BaseBCman;
static BCmgr*  PertBCman;

// -- These are shared/utilised by both systems.

static Mesh*            mesh;
static vector<Element*> elmt;

// -- File-scope routines.

static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs    (int, char**, int&, int&, real&, real&, int&,
			char*&, char*&);
static int  preprocess (const char*, const char*);
static void initVec    (real*);
static void NS_update  (const real*, real*);

// -- Routines passed to linear system solvers.

#if defined (NSPCG)
void matvec (const real*, const integer*, const real*, const integer*,
	     const integer&, const real*, real*);
void ident  (const real*, const integer*, const real*, const integer*,
	     const integer&, const real*, real*);
#else
void matvec (const real&, const real*, const real&, real*);
void ident  (real*, const real*);
#endif


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#if !defined (__DECCXX)
  ios::sync_with_stdio();
#endif

  int  maxiLsys = 100, maxiNewt = 20, i, itn;
  real tol, pretol, tolLsys = 1.0e-6, tolNewt = 1.0e-6, rnorm;
  int  verbose  = 0, converged = 0, ier;
  char *BaseSession, *PertSession;

#if defined (NSPCG)
  int  iparm[25];
  real rparm[16], ubar[1];
#endif

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, maxiLsys, maxiNewt, tolLsys, tolNewt, verbose,
	   BaseSession, PertSession);

  // -- Echo execution parameters.

  cout << "-- Newton convergence tol  : " << tolNewt  << endl;
  cout << "          iteration limit  : " << maxiNewt << endl;
  cout << "-- BiCGS  convergence tol  : " << tolLsys  << endl;
  cout << "          iteration limit  : " << maxiLsys << endl;

  // -- Allocate storage.
  
  const int ntot = preprocess (BaseSession, PertSession);

#if defined (NSPCG)
  int       nw   = 9 * ntot;
  const int wdim = 3 * ntot + nw;
#else
  const int wdim = 10 * ntot;
#endif

  vector<real> work (wdim);
  Veclib::zero (wdim, work(), 1);

  real* u    = work();
  real* U    = u  + ntot;
  real* dU   = U  + ntot;
  real* lwrk = dU + ntot; 

  cout.setf (ios::scientific, ios::floatfield); cout.precision (2);

  initVec (U);
  pretol = INITOL;
  
  // -- Newton iteration.

  for (i = 1; !converged && i <= maxiNewt; i++) {
  
    NS_update (U, dU);
    Veclib::vsub (ntot, dU, 1, U, 1, dU, 1);
    
    itn = maxiLsys;
    tol = pretol;

    Veclib::zero (ntot, u, 1);

#if defined (NSPCG)
    F77NAME (dfault) (iparm, rparm);
    iparm[1] = itn; iparm[2] = 0; rparm[0] = pretol;

    F77NAME (bcgsw) (matvec, ident, ident, 0, 0, 0, 0,
		     ntot, u, ubar, dU, lwrk, nw, iparm, rparm, ier);

    itn = iparm[1]; tol = rparm[6];
#else
    F77NAME (bicgstab) (ntot, dU, u, lwrk, ntot, itn, tol, matvec, ident, ier);
#endif

    if (ier < 0) {
      cout << "WARNING: error return from iterative solver" << endl;
      break;
    }

    rnorm     = sqrt (Blas::nrm2 (ntot, u, 1) / Blas::nrm2 (ntot, U, 1));
    converged = rnorm < tolNewt;

    cout << "Iteration "     << setw(3) << i
	 << ", tol: "        << pretol
         << ", BiCGS itns: " << setw(3) << itn
	 << ", resid: "      << tol
	 << ", Rnorm: "      << rnorm 
	 << endl;

    if (itn < NSTABLE)		// -- Tighten tolerance.
      pretol *= (pretol < tolLsys) ? 1.0 : SHRINK;
    else if (itn == maxiLsys)	// -- Loosen tolerance (and try again).
      pretol /= SHRINK;

    if (itn < maxiLsys)		// -- Accept adjustment to solution.
      Veclib::vsub (ntot, U, 1, u, 1, U, 1);
  }

  if (converged)
    cout << "Writing converged solution."                     << endl;
  else
    cout << "Not converged: writing final solution estimate." << endl;

  BaseDomain -> dump();

  Femlib::finalize();
  return (EXIT_SUCCESS);
}


static void getargs (int    argc    ,
		     char** argv    , 
		     int&   maxiLsys,
		     int&   maxiNewt,
		     real&  tolLsys ,
		     real&  tolNewt ,
		     int&   verbose ,
		     char*& Bsession,
		     char*& Psession)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "newton [options] basesession pertsession\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-m <num> ... set maximum number of Newton iterations\n"
    "-c <num> ... set convergence tolerance on Newton iteration\n"
    "-n <num> ... set maximum number of BiCGS iterations per Newton step\n"
    "-t <num> ... set BiCGS convergence tolerance\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      verbose = 1;
      Femlib::value ("VERBOSE", 
		     static_cast<integer>(Femlib::value("VERBOSE") + 1));
      break;
    case 'c':
      if (*++argv[0]) tolNewt  = atof (  *argv);
      else { --argc;  tolNewt  = atof (*++argv); }
      break;
    case 'm':
      if (*++argv[0]) maxiNewt = atoi (  *argv);
      else { --argc;  maxiNewt = atoi (*++argv); }
      break;
    case 'n':
      if (*++argv[0]) maxiLsys = atoi (  *argv);
      else { --argc;  maxiLsys = atoi (*++argv); }
      break;
    case 't':
      if (*++argv[0]) tolLsys  = atof (  *argv);
      else { --argc;  tolLsys  = atof (*++argv); }
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc != 2) message (prog, "no session files", ERROR);
  else { Bsession = *argv; Psession = *++argv; }
}


static int preprocess (const char* Bsession,
		       const char* Psession)
// ---------------------------------------------------------------------------
// Create objects needed for semtex execution, given the session file names.
//
// Return length of an solution vector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real* z;
  integer     i, np, nel, nz, nf, cyl;

  BaseFile = new FEML (Bsession);
  PertFile = new FEML (Psession);

  mesh = new Mesh (BaseFile);

  cyl  = static_cast<integer>(Femlib::value ("CYLINDRICAL"));
  np   = static_cast<integer>(Femlib::value ("N_POLY"));
  nz   = static_cast<integer>(Femlib::value ("N_Z"));
  nel  = mesh -> nEl();

  Geometry::set (np, nz, nel, (cyl)?Geometry::Cylindrical:Geometry::Cartesian);

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

  BaseBCman  = new BCmgr  (BaseFile, elmt);
  PertBCman  = new BCmgr  (PertFile, elmt);
  BaseDomain = new Domain (BaseFile, elmt, BaseBCman);
  PertDomain = new Domain (PertFile, elmt, PertBCman);

  nf = BaseDomain -> nField();

  // -- Tie the two systems together.

  for (i = 0; i < nf; i++) {
    PertDomain -> U[i]    = BaseDomain -> u[i];
    PertDomain -> Udat[i] = BaseDomain -> udat[i];
  }

  // -- Over-ride any CHKPOINT flag in session file.

  Femlib::value ("CHKPOINT", 1);

  BaseDomain -> restart();
  BaseDomain -> report ();

  return (nf - 1) * Geometry::planeSize() * Geometry::nZ();
}


static void initVec (real* tgt)
// ---------------------------------------------------------------------------
// Initialise tgt from initial guess stored in base flow.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int NC = BaseDomain -> nField() - 1;
  const int NP = Geometry::planeSize();
  const int NZ = Geometry::nZ();

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      BaseDomain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void NS_update (const real* src,
		       real*       tgt)
// ---------------------------------------------------------------------------
// Generate tgt by integrating NS problem for input velocity src.
// Then reinstate src as the base flow velocity field.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int NC = BaseDomain -> nField() - 1;
  const int NP = Geometry::planeSize();
  const int NZ = Geometry::nZ();

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      BaseDomain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);

  BaseDomain -> step = 0;
  BaseDomain -> time = 0.0;

  integrate (BaseDomain, nonlinear);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++) {
      BaseDomain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
      BaseDomain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);
    }
}

#if defined (NSPCG)

void matvec (const real*    dum0,
	     const integer* dum1,
	     const real*    dum2,
	     const integer* dum3,
	     const integer& dum4,
	     const real*    src ,
	     real*          tgt )
// ---------------------------------------------------------------------------
// Integrate linear NS problem for input velocity x, generating y = LNS(x)-x.
// This routine is used by NPSCG BCGS solver, as operator SUBA.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int NC   = PertDomain -> nField() - 1;
  const int NP   = Geometry::planeSize();
  const int NZ   = Geometry::nZ();
  const int ntot = NC * NP * NZ;

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      PertDomain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);

  PertDomain -> step = 0;
  PertDomain -> time = 0.0;

  integrate (PertDomain, linear);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      PertDomain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);

  Veclib::vsub (ntot, tgt, 1, src, 1, tgt, 1);
}


void ident (const real*    dum0,
	    const integer* dum1,
	    const real*    dum2,
	    const integer* dum3,
	    const integer& dum4,
	    const real*    src ,
	    real*          tgt )
// ---------------------------------------------------------------------------
// This is a preconditioner routine supplied to BCGSW, as COPY.
// As we don't have a preconditioner for the problem, this is an identity. 
// ---------------------------------------------------------------------------
{
  const int NC   = PertDomain -> nField() - 1;
  const int NP   = Geometry::planeSize();
  const int NZ   = Geometry::nZ();
  const int ntot = NC * NP * NZ;

  Veclib::copy (ntot, src, 1, tgt, 1);
}

#else

void matvec (const real& alpha,
	     const real* x    ,
	     const real& beta ,
	     real*       y    )
// ---------------------------------------------------------------------------
// Integrate linear NS problem for input velocity x, generate
//   y = alpha [LNS(x) - x] + beta y.
// This routine is used by BICGSTAB solver, as operator MATVEC.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int NC   = PertDomain -> nField() - 1;
  const int NP   = Geometry::planeSize();
  const int NZ   = Geometry::nZ();
  const int ntot = NC * NP * NZ;
  vector<real> work (ntot);

  Veclib::copy (ntot, y, 1, work(), 1);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      PertDomain -> u[i] -> setPlane (k, x + (i*NZ + k)*NP);

  PertDomain -> step = 0;
  PertDomain -> time = 0.0;

  integrate (PertDomain, linear);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      PertDomain -> u[i] -> getPlane (k, y + (i*NZ + k)*NP);

  Veclib::vsub (ntot, y, 1, x, 1, y, 1);
  Blas::scal   (ntot, alpha, y, 1);
  Blas::axpy   (ntot, beta, work(), 1, y, 1);
}


void ident (real*       tgt,
	    const real* src)
// ---------------------------------------------------------------------------
// This is a preconditioner routine supplied to BICGSTAB, as PSOLVE.
// As we don't have a preconditioner for the problem, this is an identity. 
// ---------------------------------------------------------------------------
{
  const int NC   = PertDomain -> nField() - 1;
  const int NP   = Geometry::planeSize();
  const int NZ   = Geometry::nZ();
  const int ntot = NC * NP * NZ;

  Veclib::copy (ntot, src, 1, tgt, 1);
}

#endif



