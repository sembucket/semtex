///////////////////////////////////////////////////////////////////////////////
// newton.C: solve a time-invariant NS flow by a Newton iteration.
//
// Copyright (C) 2002 Hugh Blackburn.
//
// The iterative method was derived by Laurette Tuckerman, see Refs
// [1] & [2], and is based on time-stepping with both linearised and
// full Navier--Stokes integrators.  The matrix-free linear systems
// are here solved using the Bi-Conjugate-Gradients-Stabilized,
// algorithm, with code from the Templates package, Ref [3].
//
// USAGE
// -----
// newton [options] session
//   session: specifies name of semtex session file.
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
// A number of semtex files are required --
//   session:     semtex session file
//   session.num: computed automatically if not supplied
//   session.rst: restart file (initialised with white noise if not supplied) 
//
// REFERENCES
// ----------
// [1] C.K. Mamun & L.S. Tuckerman (1995), "Asymmetry and Hopf Bifurcation
//     in spherical Couette flow, Phys Fluids V7N1, 80--91.
// [2] L.S. Tuckerman & D. Barkley (2000), "Bifurcation analysis for
//     timesteppers", in Numerical Methods for Bifurcation Problems,
//     ed E. Doedel & L.S. Tuckerman, Springer. 453--466.
// [3] Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
//     Eijkhout, Pozo, Romine, and van der Vorst (1993), "Templates for the 
//     Solution of Linear Systems: Building Blocks for Iterative 
//     Methods", SIAM Publications.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "newt.h"
#include <new.h>

static char             prog[] = "newton";
static char*            session;
static Domain*          domain;
static Analyser*        analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs    (int, char**, int&, int&, real&, real&, int&, char*&);
static int  preprocess (const char*);
static void initVec    (real*);
static void NS_update  (const real*, real*);

void matvec (const real&, const real*, const real&, real*);
void ident  (real*, const real*);


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

  int  maxiLsys = 100,    maxiNewt = 20, i, itn;
  real tolLsys  = 1.0e-6, tolNewt  = 1.0e-6, rnorm, tol;
  int  verbose  = 0, converged = 0, info;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, maxiLsys, maxiNewt, tolLsys, tolNewt, verbose, session);

  // -- Allocate storage.
  
  const int ntot = preprocess (session);
  const int wdim = 10 * ntot;

  vector<real> work (wdim);
  Veclib::zero (wdim, work(), 1);

  real* u    = work();
  real* U    = u  + ntot;
  real* dU   = U  + ntot;
  real* lwrk = dU + ntot; 

  initVec (U);

  // -- Newton iteration.

  for (i = 1; !converged && i <= maxiNewt; i++) {
  
    NS_update (U, dU);
    Veclib::vsub (ntot, dU, 1, U, 1, dU, 1);

    itn = maxiLsys; tol = tolLsys; Veclib::zero (ntot, u, 1);

    F77NAME(bicgstab) (ntot, dU, u, lwrk, ntot, itn, tol, matvec, ident, info);
    
    if (info < 0) message (prog, "error return from bicgstab", ERROR);

    rnorm     = sqrt (Blas::nrm2 (ntot, u, 1) / Blas::nrm2 (ntot, U, 1));
    converged = rnorm < tolNewt;

    cout << "Iteration "           << setw(3) << i 
	 << ", BiCGS iterations: " << setw(3) << itn
	 << ", resid: "            << setw(6) << tol
	 << ", Rnorm: "            << setw(6) << rnorm << endl;

    Veclib::vsub (ntot, U, 1, u, 1, U, 1);
  }

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
		     char*& session )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "newton [options] session\n"
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
      else { --argc;  tolNewt  = atoi (*++argv); }
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

  if   (argc != 1) message (prog, "no session file",   ERROR);
  else             session = *argv;
}


static int preprocess (const char* session)
// ---------------------------------------------------------------------------
// Create objects needed for semtex execution, given the session file name.
//
// Return length of an solution vector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real* z;
  integer     i, np, nel, nz, cyl;

  file = new FEML (session);
  mesh = new Mesh (file);

  cyl  = static_cast<integer>(Femlib::value ("CYLINDRICAL"));
  np   = static_cast<integer>(Femlib::value ("N_POLY"));
  nz   = static_cast<integer>(Femlib::value ("N_Z"));
  nel  = mesh -> nEl();

  Geometry::set (np, nz, nel, (cyl)?Geometry::Cylindrical:Geometry::Cartesian);

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

  bman    = new BCmgr    (file, elmt);
  domain  = new Domain   (file, elmt, bman);
  analyst = new Analyser (domain, file);

  // -- Over-ride any CHKPOINT flag in session file.

  Femlib::value ("CHKPOINT", 1);

  domain -> restart ();
  domain -> report  ();

  return (domain->nField() - 1) * Geometry::planeSize() * Geometry::nZ();
}


static void initVec (real* tgt)
// ---------------------------------------------------------------------------
// Initialise tgt from initial guess stored in base flow.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int NC = domain -> nField() - 1;
  const int NP = Geometry::planeSize();
  const int NZ = Geometry::nZ();

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      domain -> U[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void NS_update (const real* src,
		       real*       tgt)
// ---------------------------------------------------------------------------
// Generate tgt by integrating NS problem for input velocity src.
// Also set base flow to input velocity.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int NC = domain -> nField() - 1;
  const int NP = Geometry::planeSize();
  const int NZ = Geometry::nZ();

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++) {
      domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);
      domain -> U[i] -> setPlane (k, src + (i*NZ + k)*NP);
    }

  domain -> step = 0;
  domain -> time = 0.0;

  integrate (domain, analyst, nonlinear);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


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
  const int NC   = domain -> nField() - 1;
  const int NP   = Geometry::planeSize();
  const int NZ   = Geometry::nZ();
  const int ntot = NC * NP * NZ;
  vector<real> work (ntot);

  Veclib::copy (ntot, y, 1, work(), 1);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> setPlane (k, x + (i*NZ + k)*NP);

  domain -> step = 0;
  domain -> time = 0.0;

  integrate (domain, analyst, linear);

  for (i = 0; i < NC; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, y + (i*NZ + k)*NP);

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
  const int NC   = domain -> nField() - 1;
  const int NP   = Geometry::planeSize();
  const int NZ   = Geometry::nZ();
  const int ntot = NC * NP * NZ;

  Veclib::copy (ntot, src, 1, tgt, 1);
}

