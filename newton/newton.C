///////////////////////////////////////////////////////////////////////////////
// newton.C: solve a time-invariant NS flow by a Newton iteration.
//
// Copyright (c) 2002 <--> $Date$, Hugh Blackburn
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
// [1] CK Mamun & LS Tuckerman (1995), "Asymmetry and Hopf Bifurcation
//     in spherical Couette flow, Phys Fluids 7(1), 80--91.
// [2] LS Tuckerman & D Barkley (2000), "Bifurcation analysis for
//     timesteppers", in Numerical Methods for Bifurcation Problems,
//     ed E. Doedel & L.S. Tuckerman, Springer. 453--466.
// [3] R Barrett, M Berry, T Chan, J Demmel, J Donato, J Dongarra, 
//     V Eijkhout, R Pozo, C Romine & H van der Vorst (1993),
//     "Templates for the Solution of Linear Systems: Building Blocks for
//     Iterative Methods", SIAM Publications.
// [4] TC Oppe, WD Joubert & DR Kincaid (1988), "NSPCG User's Guide",
//     Center for Numerical Analysis, University of Texas at Austin.
// [5] HM Blackburn (2002) "Three-dimensional instability and state
//     selection in an oscillatory axisymmetric swirling flow",
//     Phys Fluids 14(11): 3983--3996.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "newt.h"

static char prog[] = "newton";

// -- Linear system convergence control.

const real_t INITOL  = 1.0e-3;
const real_t SHRINK  = 0.67;
const int_t  NSTABLE = 3;

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

static void  memExhaust () { message ("new", "free store exhausted", ERROR); }
static void  getargs    (int, char**, int_t&, int_t&, real_t&, real_t&,
			 int_t&, char*&, char*&);
static int_t preprocess (const char*, const char*);
static void  initVec    (real_t*);
static void  NS_update  (const real_t*, real_t*);

// -- Routines passed to linear system solvers.

#if defined (NSPCG)
void matvec (const real_t*, const int_t*, const real_t*, const int_t*,
	     const int_t&, const real_t*, real_t*);
void ident  (const real_t*, const int_t*, const real_t*, const int_t*,
	     const int_t&, const real_t*, real_t*);
#else  // -- TEMPLATES code, the default.
void matvec (const real_t&, const real_t*, const real_t&, real_t*);
void ident  (real_t*, const real_t*);
#endif


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  int_t  maxiLsys = 100, maxiNewt = 20, i, itn;
  real_t tol, pretol, tolLsys = 1.0e-6, tolNewt = 1.0e-6, rnorm;
  int_t  verbose  = 0, ier;
  bool   converged = false;
  char   *BaseSession, *PertSession;

#if defined (NSPCG)
  int_t  iparm[25];
  real_t rparm[16], ubar[1];
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
  
  const int_t ntot = preprocess (BaseSession, PertSession);

#if defined (NSPCG)
  int_t       nw   = 9 * ntot;
  const int_t wdim = 3 * ntot + nw;
#else
  const int_t wdim = 10 * ntot;
#endif

  vector<real_t> work (wdim);
  Veclib::zero (wdim, &work[0], 1);

  real_t* u    = &work[0];
  real_t* U    = u  + ntot;
  real_t* dU   = U  + ntot;
  real_t* lwrk = dU + ntot; 

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


static void getargs (int     argc    ,
		     char**  argv    , 
		     int_t&  maxiLsys,
		     int_t&  maxiNewt,
		     real_t& tolLsys ,
		     real_t& tolNewt ,
		     int_t&  verbose ,
		     char*&  Bsession,
		     char*&  Psession)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "newton [options] basesession pertsession\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-m <num> ... set maximum number of Newton iterations, "
    "default=20\n"
    "-c <num> ... set convergence tolerance on Newton iteration, "
    "default=1e-6\n"
    "-n <num> ... set maximum number of BiCGS iterations per Newton step, "
    "default=100\n"
    "-t <num> ... set BiCGS convergence tolerance, "
    "default=1e-6\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      verbose = 1;
      Femlib::ivalue ("VERBOSE", Femlib::ivalue("VERBOSE") + 1);
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


static int_t preprocess (const char* Bsession,
			 const char* Psession)
// ---------------------------------------------------------------------------
// Create objects needed for semtex execution, given the session file names.
//
// Return length of an solution vector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real_t* z;
  int_t         i, np, nel, nz, nf, cyl;

  BaseFile = new FEML (Bsession);
  PertFile = new FEML (Psession);

  mesh = new Mesh (BaseFile);

  cyl  = Femlib::ivalue ("CYLINDRICAL");
  np   = Femlib::ivalue ("N_POLY");
  nz   = Femlib::ivalue ("N_Z");
  nel  = mesh -> nEl();

  Geometry::set (np, nz, nel, (cyl)?Geometry::Cylindrical:Geometry::Cartesian);

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

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

  Femlib::ivalue ("CHKPOINT", 1);

  BaseDomain -> restart();
  BaseDomain -> report ();

  return (nf - 1) * Geometry::planeSize() * Geometry::nZ();
}


static void initVec (real_t* tgt)
// ---------------------------------------------------------------------------
// Initialise tgt from initial guess stored in base flow.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t NC = BaseDomain -> nField() - 1;
  const int_t NP = Geometry::planeSize();
  const int_t NZ = Geometry::nZ();

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      BaseDomain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void NS_update (const real_t* src,
		       real_t*       tgt)
// ---------------------------------------------------------------------------
// Generate tgt by integrating NS problem for input velocity src.
// Then reinstate src as the base flow velocity field.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t NC = BaseDomain -> nField() - 1;
  const int_t NP = Geometry::planeSize();
  const int_t NZ = Geometry::nZ();

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

void matvec (const real_t* dum0,
	     const int_t*  dum1,
	     const real_t* dum2,
	     const int_t*  dum3,
	     const int_t&  dum4,
	     const real_t* src ,
	     real_t*       tgt )
// ---------------------------------------------------------------------------
// Integrate linear NS problem for input velocity x, generating y = LNS(x)-x.
// This routine is used by NPSCG BCGS solver, as operator SUBA.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t NC   = PertDomain -> nField() - 1;
  const int_t NP   = Geometry::planeSize();
  const int_t NZ   = Geometry::nZ();
  const int_t ntot = NC * NP * NZ;

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


void ident (const real_t* dum0,
	    const int_t*  dum1,
	    const real_t* dum2,
	    const int_t*  dum3,
	    const int_t&  dum4,
	    const real_t* src ,
	    real_t*       tgt )
// ---------------------------------------------------------------------------
// This is a preconditioner routine supplied to BCGSW, as COPY.
// As we don't have a preconditioner for the problem, this is an identity. 
// ---------------------------------------------------------------------------
{
  const int_t NC   = PertDomain -> nField() - 1;
  const int_t NP   = Geometry::planeSize();
  const int_t NZ   = Geometry::nZ();
  const int_t ntot = NC * NP * NZ;

  Veclib::copy (ntot, src, 1, tgt, 1);
}

#else  // -- Default use of TEMPLATES routine.

void matvec (const real_t& alpha,
	     const real_t* x    ,
	     const real_t& beta ,
	     real_t*       y    )
// ---------------------------------------------------------------------------
// Integrate linear NS problem for input velocity x, generate
//   y = alpha [LNS(x) - x] + beta y.
// This routine is used by BICGSTAB solver, as operator MATVEC.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t NC   = PertDomain -> nField() - 1;
  const int_t NP   = Geometry::planeSize();
  const int_t NZ   = Geometry::nZ();
  const int_t ntot = NC * NP * NZ;
  vector<real_t> work (ntot);

  Veclib::copy (ntot, y, 1, &work[0], 1);

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
  Blas::axpy   (ntot, beta, &work[0], 1, y, 1);
}


void ident (real_t*       tgt,
	    const real_t* src)
// ---------------------------------------------------------------------------
// This is a preconditioner routine supplied to BICGSTAB, as PSOLVE.
// As we don't have a preconditioner for the problem, this is an identity. 
// ---------------------------------------------------------------------------
{
  const int_t NC   = PertDomain -> nField() - 1;
  const int_t NP   = Geometry::planeSize();
  const int_t NZ   = Geometry::nZ();
  const int_t ntot = NC * NP * NZ;

  Veclib::copy (ntot, src, 1, tgt, 1);
}

#endif



