///////////////////////////////////////////////////////////////////////////////
// arnoldi.C: compute leading eigenvalues and eigenvectors for
// stability analysis based on linearised (Navier--Stokes) operators.
// The base flow can be either steady or periodic in time, two or
// three component, cylindrical or Cartesian, but must be
// two-dimensional.
// 
// Based on code floK by Dwight Barkley & Ron Henderson.
//
// Copyright (C) 1999-2002 Hugh Blackburn.
//
// The eigenpairs computed in the subspace are related to the Ritz
// estimates of those in the original space in a simple way: the
// eigenvalues are the same, and the Ritz eigenvectors are related to
// the subspace eigenvectors through a linear transformation (see
// Saad, p.175).
//
// USAGE
// -----
// arnoldi [options] session
//   session: specifies name of semtex session file.
//   options:
//   -h       ... print this message
//   -v       ... set verbose
//   -k <num> ... set dimension of subspace (maximum number of pairs) to num
//   -m <num> ... set maximum number of iterations         (m >= k)
//   -n <num> ... compute num eigenvalue/eigenvector pairs (n <= k)
//   -t <num> ... set eigenvalue tolerance to num [Default = 1e-6].
// 
// FILES
// -----
// A number of semtex files are required --
//   session:     semtex session file
//   session.num: computed automatically if not supplied
//   session.bse: base flow, containing N_SLICE field dumps
//   session.rst: restart file (initialised with white noise if not supplied) 
//
// REFERENCES
// ----------
// [1]  D. Barkley & R.D. Henderson (1996), "Three-dimensional Floquet
//      stability analysis of the wake of a circular cylinder",
//      J. Fluid Mech V322, 215--241.
// [2]  Y. Saad (1991), "Numerical methods for large eigenvalue problems"
//      Wiley.
// [3]  L.S. Tuckerman & D. Barkley (2000), "Bifurcation analysis for
//      timesteppers", in Numerical Methods for Bifurcation Problems,
//      ed E. Doedel & L.S. Tuckerman, Springer. 453--466.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "stab.h"
#include <new.h>

static char             prog[] = "arnoldi";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void memExhaust () { message ("new", "free store exhausted", ERROR); }
static void getargs    (int, char**, int&, int&, int&, int&, real&, char*&);
static int  preprocess (const char*);

static void EV_init    (real*);
static void EV_update  (const real*, real*);
static void EV_small   (real**, const int, const int, 
		        real*, real*, real*, real&, const int); 
static int  EV_test    (const int, const int, real*, real*, real*,
		        const real, const real, const int);
static void EV_sort    (real*, real*, real*, real*, const int);
static void EV_post    (real**, real**, const int, const int, 
			const int, const real*, const real*, 
			const real*, const int);
static void EV_big     (real**, real**, const int, const int, const int,
			const real*, const real*, const real*);


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

  int  kdim = 2, nvec = 2, nits = 2, verbose = 0, converged = 0;
  real norm, resnorm, evtol = 1.0e-6;
  int  i, j;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, kdim, nits, nvec, verbose, evtol, session);

  // -- Check, echo parameter values.

  if (kdim < 1)    message (prog, "param error: KDIM must be > 1",     ERROR);
  if (nvec < 1)    message (prog, "param error: NVEC must be > 1",     ERROR);
  if (nits < kdim) message (prog, "param error: NITS must be >= KDIM", ERROR);
  if (kdim < nvec) message (prog, "param error: NVEC must be <= KDIM", ERROR);

  cout << "-- Krylov dimension        : " << kdim  << endl;
  cout << "   Convergence dimension   : " << nvec  << endl;
  cout << "   Convergence tolerance   : " << evtol << endl;
  cout << "   Maximum iterations      : " << nits  << endl;

  // -- Allocate eigenproblem storage.
  
  const int ntot = preprocess (session);
  const int wdim = kdim + kdim + (kdim * kdim) + 2*ntot*(kdim + 1);

  vector<real> work (wdim);
  Veclib::zero (wdim, work(), 1);

  real*  wr   = work();
  real*  wi   = wr   + kdim;
  real*  zvec = wi   + kdim;
  real*  kvec = zvec + kdim * kdim;
  real*  tvec = kvec + ntot * (kdim + 1);
  real** Kseq = new real* [kdim + 1];
  real** Tseq = new real* [kdim + 1];

  for (i = 0; i <= kdim; i++) {
    Kseq[i] = kvec + i * ntot;
    Tseq[i] = tvec + i * ntot;
  }

  // -- Load starting vector.

  EV_init (Kseq[0]);
  norm = Blas::nrm2 (ntot, Kseq[0], 1);
  Blas::scal (ntot, 1.0/norm, Kseq[0], 1);

  // -- Fill initial Krylov sequence -- during which convergence may occur.

  for (i = 1; !converged && i <= kdim; i++) {
    EV_update (Kseq[i - 1], Kseq[i]);
    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small  (Tseq, ntot, i, zvec, wr, wi, resnorm, verbose);
    converged = EV_test (i, i, zvec, wr, wi, resnorm, evtol, min(i, nvec));
    converged = max (converged, 0); // -- Only exit on evtol.
  }

  // -- Carry out iterative solution.

  for (i = kdim + 1; !converged && i <= nits; i++) {

    // -- Normalise Krylov sequence & update.

    norm = Blas::nrm2 (ntot, Kseq[1], 1);
    for (j = 1; j <= kdim; j++) {
      Blas::scal   (ntot, 1.0/norm, Kseq[j], 1);
      Veclib::copy (ntot, Kseq[j], 1, Kseq[j - 1], 1);
    }
    EV_update (Kseq[kdim - 1], Kseq[kdim]);
    
    // -- Get subspace eigenvalues, test for convergence.

    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, kdim, zvec, wr, wi, resnorm, verbose); 

    converged = EV_test (i, kdim, zvec, wr, wi, resnorm, evtol, nvec);
  }

  EV_post (Tseq, Kseq, ntot, min(--i, kdim), nvec, zvec, wr, wi, converged);

  Femlib::finalize();
  return (EXIT_SUCCESS);
}


static void EV_init (real* tgt)
// ---------------------------------------------------------------------------
// Load initial vector from domain velocity fields.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int ND = Geometry::nPert();
  const int NP = Geometry::planeSize();
  const int NZ = Geometry::nZ();
    
  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void EV_update  (const real* src,
			real*       tgt)
// ---------------------------------------------------------------------------
// Generate tgt by applying linear operator (here, a linearised
// Navier--Stokes integrator) to src.
// ---------------------------------------------------------------------------
{
  int       i, k;
  const int ND = Geometry::nPert();
  const int NP = Geometry::planeSize();
  const int NZ = Geometry::nZ();

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);

  domain -> step = 0;
  integrate (domain, analyst);

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void EV_small (real**    Kseq   ,
		      const int ntot   ,
		      const int kdim   ,
		      real*     zvec   ,
		      real*     wr     ,
		      real*     wi     ,
		      real&     resnorm,
		      const int verbose)
// ---------------------------------------------------------------------------
// Here we take as input the Krylov sequence Kseq =
//          x,
//   A      x,
//   A^2    x,
//   A^3    x,
//    ...
//   A^kdim x;
//
// and from it find Q (an orthonormal basis), by a Gram--Schmidt
// algorithm that produces the decomposition Kseq = Q R.  Kseq is
// destroyed in the process, and replaced by an othonormal basis for
// Kseq.
//
// Then we compute Hessenberg matrix H = Q* A Q (using Q* A = R), find
// its eigenvalues (wr, wi) and (right) eigenvectors, zvec, related to
// those of A.
//
// The residual norm for each eigenvector is related to H(kdim+1, kdim),
// which is passed back for convergence testing.
// ---------------------------------------------------------------------------
{
  char         routine[] = "EV_small";
  const int    kdimp = kdim + 1;
  int          i, j, ier, lwork = 10 * kdim;
  vector<real> work (kdimp * kdimp + kdim * kdim + lwork);
  real         *R     = work(),
               *H     = R + kdimp * kdimp,
               *rwork = H + kdim * kdim;

  Veclib::zero (kdimp * kdimp, R, 1);

  // -- Modified G--S orthonormalisation.

  for (i = 0; i < kdimp; i++) {
    real gsc = Blas::nrm2 (ntot, Kseq[i], 1);
    if (gsc == 0.0)
      message (routine, "basis vectors linearly dependent", ERROR);

    R[Veclib::col_major (i, i, kdimp)] = gsc;
    Blas::scal (ntot, 1.0 / gsc, Kseq[i], 1);
    
    for (j = i + 1; j < kdimp; j++) {
      gsc = Blas::dot (ntot,  Kseq[j], 1, Kseq[i], 1);
      Blas::axpy (ntot, -gsc, Kseq[i], 1, Kseq[j], 1);
      R[Veclib::col_major (i, j, kdimp)] = gsc;
    }
  }

  // -- QR decomposition completed.  Print up R as diagnostic.

  if (verbose) {
    cout << "R =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	cout << setw (14) << R[Veclib::col_major (i, j, kdimp)];
      cout << endl;
    }
  }

  // -- H(i, j) = (q_i, A q_j) = 1 / R(j, j) * (R(i, j + 1) - H(i, l).R(l, j)),
  //    with the last inner product taken over l < j.

  for (i = 0; i < kdim; i++) {
    for (j = 0; j < kdim; j++) {
      H[Veclib::col_major (i, j, kdim)] =
	R[Veclib::col_major (i, j + 1, kdimp)]
	- Blas::dot (j, H + i, kdim, R + j * kdimp, 1);
      H[Veclib::col_major (i, j, kdim)] /= R [Veclib::col_major (j, j, kdimp)];
    }
  }

  // -- Print up H as diagnostic.

  if (verbose) {
    cout << "H =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	cout << setw (14) << H[Veclib::col_major (i, j, kdim)];
      cout << endl;
    }
  }

  // -- Find eigenpairs of H using LAPACK routine dgeev (q.v.).
  //    For complex-conjugate eigenvalues, the corresponding
  //    complex eigenvector is a real-imaginary pair of rows
  //    of zvec.

  F77NAME(dgeev) ("N","V",kdim,H,kdim,wr,wi,0,1,zvec,kdim,rwork,lwork,ier);

  if (ier) message (routine, "error return from dgeev", ERROR);

  // -- Print up eigenvectors as diagnostic.

  if (verbose) {
    cout << "zvec =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	cout << setw (14) << zvec[Veclib::col_major (i, j, kdim)];
      cout << endl;
    }
  }
  
  // -- Compute residual information.

  resnorm = fabs (R[Veclib::col_major (kdim,     kdim,     kdimp)] /
		  R[Veclib::col_major (kdim - 1, kdim - 1, kdimp)] );
}


static int EV_test (const int  itrn   ,
		    const int  kdim   ,
		    real*      zvec   ,
		    real*      wr     ,
		    real*      wi     ,
		    const real resnorm,
		    const real evtol  ,
		    const int  nvec   )
// ---------------------------------------------------------------------------
// Test convergence of eigenvalues and print up diagnostic information.
//
// Return value:
//   nvec:  all of the first nvec eigenvalue estimates have converged;
//  -1,-2:  the residuals aren't shrinking;
//   0:     neither of the above is true: not converged.
// ---------------------------------------------------------------------------
{
  int          i, idone = 0;
  vector<real> work (kdim);
  real         re_ev, im_ev, abs_ev, ang_ev, re_Aev, im_Aev;
  real*        resid = work();
  const real   period = Femlib::value ("D_T * N_STEP");
  static real  min_max1, min_max2;
 
  if (min_max1 == 0.0) min_max1 = 1000.0;
  if (min_max2 == 0.0) min_max2 = 1000.0;

  // -- Sort subspace eigenvectors by residual.

  for (i = 0; i < kdim; i++) {
    resid[i] = resnorm * fabs (zvec[kdim - 1 + i * kdim])
      / sqrt (Blas::dot (kdim, zvec + i * kdim, 1, zvec + i * kdim, 1));
    if (wi[i] < 0.0) resid[i - 1] = resid[i] = hypot (resid[i - 1], resid[i]);
  }
  EV_sort (zvec, wr, wi, resid, kdim);

  // -- Stopping test.

  if      (resid[nvec - 1] < evtol)                              idone = nvec;
  else if (min_max1 < 0.01 && resid[nvec - 1] > 10.0 * min_max1) idone = -1;
  else if (min_max2 < 0.01 && resnorm         > 10.0 * min_max2) idone = -2;

  min_max1 = min (min_max1, resid[nvec - 1]);
  min_max2 = min (min_max2, resnorm);

  // -- Print diagnostic information.

  cout << "-- Iteration = " << itrn << ", H(k+1, k) = " << resnorm << endl;

  cout.precision(4);
  cout.setf(ios::scientific, ios::floatfield);

  cout << "EV  Magnitude   Angle       Growth      Frequency   Residual"
       << endl;

  for (i = 0; (i < kdim) && (i < 10); i++) {
    re_ev  = wr[i];
    im_ev  = wi[i];
    abs_ev = hypot (re_ev, im_ev);
    ang_ev = atan2 (im_ev, re_ev);
    re_Aev = log (abs_ev) / period;
    im_Aev = ang_ev       / period;
    cout << setw(2)  << i
         << setw(12) << abs_ev
      	 << setw(12) << ang_ev
	 << setw(12) << re_Aev
	 << setw(12) << im_Aev
	 << setw(12) << resid[i]
	 << endl;
  }

  cout.precision(6);
  cout.setf(ios::fixed);
  
  return idone;
}


static void EV_sort (real*     evec,
		     real*     wr  ,
		     real*     wi  ,
		     real*     test,
		     const int dim )
// ---------------------------------------------------------------------------
// Insertion sort to rearrange eigenvalues and eigenvectors to ascending
// order according to vector test.  See equivalent Numerical Recipes routine.
// ---------------------------------------------------------------------------
{
  int          i, j;
  vector<real> work (dim);
  real         wr_tmp, wi_tmp, te_tmp, *z_tmp = work();

  for (j = 1; j < dim; j++) {
    wr_tmp = wr  [j];
    wi_tmp = wi  [j];
    te_tmp = test[j];
    Veclib::copy (dim, evec + j * dim, 1, z_tmp, 1);
    i = j - 1;
    while (i >= 0 && test[i] > te_tmp) {
      wr  [i + 1] = wr  [i];
      wi  [i + 1] = wi  [i];
      test[i + 1] = test[i];
      Veclib::copy (dim, evec + i * dim, 1, evec + (i + 1) * dim, 1);
      i--;
    }
    wr  [i + 1] = wr_tmp;
    wi  [i + 1] = wi_tmp;
    test[i + 1] = te_tmp;
    Veclib::copy (dim, z_tmp, 1, evec + (i + 1) * dim, 1);
  }
}


static void EV_post (real**      Tseq,
		     real**      Kseq, 
		     const int   ntot, 
		     const int   kdim, 
		     const int   nvec,
		     const real* zvec, 
		     const real* wr  , 
		     const real* wi  , 
		     const int   icon)
// ---------------------------------------------------------------------------
// Carry out postprocessing of estimates, depending on value of icon
// (as output by EV_test).
//
// icon == 0:
//   Solution has not converged.  The final value of the Krylov sequence 
//   is presumed to have already been output within the integration loop.
// icon == -1, -2:
//   Solution has not converged, but the residuals are getting no smaller.
//   In this case we write out the starting vector for the last kdim
//   iterates.
// icon == nvec:
//   All the eigenvalue estimates under test have converged, and they are 
//   output to separate files.
// ---------------------------------------------------------------------------
{
  const char* routine = "EV_post";
  char        msg[StrMax], nom[StrMax];
  int         i, j, k;
  const int   ND = Geometry::nPert();
  const int   NP = Geometry::planeSize();
  const int   NZ = Geometry::nZ();
  const real* src;
  ofstream    file;

  if (icon == 0) {

    cout << prog
	 << ": not converged, writing final Krylov vector."
	 << endl;

    // -- At present, this output is dealt with automatically in integrate().
    //    But we do it again here; we may want to modify Domain::dump().

    src = Kseq[kdim];
    for (i = 0; i < ND; i++)
      for (k = 0; k < NZ; k++)
	domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);
    strcat    (strcpy (nom, domain -> name), ".fld");
    file.open (nom, ios::out); file << *domain; file.close();

  } else if (icon < 0) {
    
    cout << prog
	 << ": minimum residual reached, writing initial Krylov vector."
	 << endl;
    
    src = Kseq[0];
    for (i = 0; i < ND; i++)
      for (k = 0; k < NZ; k++)
	domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);
    strcat    (strcpy (nom, domain -> name), ".fld");
    file.open (nom, ios::out); file << *domain; file.close();

  } else if (icon == nvec) {

    cout << prog
	 << ": converged, writing "
	 << icon
	 << " eigenvectors."
	 << endl;
    
    EV_big (Tseq, Kseq, ntot, kdim, icon, zvec, wr, wi);
    for (j = 0; j < icon; j++) {
      src = Kseq[j];
      for (i = 0; i < ND; i++)
	for (k = 0; k < NZ; k++)
	  domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);
      sprintf   (msg, ".eig.%1d", j);
      strcat    (strcpy (nom, domain -> name), msg);
      file.open (nom, ios::out); file << *domain; file.close();
    }

  } else {
    sprintf (msg, "input convergence value %1d: not recognised", icon);
    message (routine, msg, ERROR);
  }
}


static void EV_big (real**      bvecs,
	            real**      evecs,
		    const int   ntot ,
		    const int   kdim ,
		    const int   nvec ,
		    const real* zvecs,
		    const real* wr   ,
		    const real* wi   )
// ---------------------------------------------------------------------------
// Compute the Ritz eigenvector estimates of the linear operator using
// the eigenvalues and eigenvectors of H computed in the subspace.
// 
// Input
// -----
// bvecs: orthonormal basis of the Krylov sequence 
//        (produced in EV_small), dimensions ntot * kdim.
// zvecs: eigenvectors of H, dimensions kdim * kdim.
// wr,wi: eigenvalues  of H, each kdim long.
//
// Output
// ------
// evecs: eigenvector estimates on the original space, ntot * kdim.
//
// The Ritz estimates are computed as (see Saad, p.175):
//
//        [evecs]            = [bvecs]            [zvecs]
//               ntot x kdim          ntot x kdim        kdim x kdim
//
// Potentially, convergence has been achieved on nvec eigenpairs
// before the nominated Krylov dimension has been filled, so the
// computations are only done for the first nvec eigenvectors.
// ---------------------------------------------------------------------------
{
  real norm, wgt;
  int  i, j;

  // -- Generate big e-vectors.

  for (j = 0; j < nvec; j++) {
    Veclib::zero (ntot, evecs[j], 1);
    for (i = 0; i < nvec; i++) {
      wgt = zvecs[Veclib::col_major (i, j, kdim)];
      Blas::axpy (ntot, wgt, bvecs[i], 1, evecs[j], 1);
    }
  }

  // -- Normalize big e-vectors.

  for (i = 0; i < nvec; i++)
    if (wi[i] == 0.0) {
      norm = Blas::nrm2 (ntot, evecs[i], 1);
      Blas::scal (ntot, 1.0/norm, evecs[i], 1);
    } else if (wi[i] > 0.0) {
      norm  = sqr (Blas::nrm2 (ntot, evecs[i],   1));
      norm += sqr (Blas::nrm2 (ntot, evecs[i+1], 1));
      norm = sqrt (norm);
      Blas::scal (ntot, 1.0/norm, evecs[i],   1);
      Blas::scal (ntot, 1.0/norm, evecs[i+1], 1);
      i++;
    }

#if BIG_RESIDS  // -- This is not yet recoded: HMB 8/2/2002.

  // Compute residuals of big vectors directly (i.e. A(u)-mu(u)).
  real resid;
 
  for (i = 0; i < nwrt; i++) {
    V_copy(ntot, evecs[i], bvecs[0]);
    A_op(bvecs[0]);
    tsteps -= nsteps;
    ftime  -= nsteps*dt;
    if(wi[i]==0.) {
      V_axpy(ntot, -wr[i], evecs[i], bvecs[0]);
      resid = V_nrm2(ntot, bvecs[0]);
      printf("big resid(%d) = %g \n", i, resid);
    }
    else if(wi[i] > 0.) {
      V_axpy(ntot, -wr[i], evecs[i],   bvecs[0]);
      V_axpy(ntot,  wi[i], evecs[i+1], bvecs[0]);
      resid = V_nrm2(ntot, bvecs[0])/V_nrm2(ntot, evecs[i]);
      printf("big resid(%d) = %g \n", i, resid);
    }
    else {
      V_axpy(ntot, -wr[i], evecs[i],   bvecs[0]);
      V_axpy(ntot,  wi[i], evecs[i-1], bvecs[0]);
      resid = V_nrm2(ntot, bvecs[0])/V_nrm2(ntot, evecs[i]);
      printf("big resid(%d) = %g \n", i, resid);
    }
  }
#endif
}


static void getargs (int    argc   ,
		     char** argv   ,
		     int&   kdim   , 
		     int&   maxit  ,
		     int&   neval  ,
		     int&   verbose,
		     real&  evtol  ,
		     char*& session)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "arnoldi [options] session\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-k <num> ... set dimension of subspace (maximum number of pairs) to num\n"
    "-m <num> ... set maximum number of iterations         (m >= k)\n"
    "-n <num> ... compute num eigenvalue/eigenvector pairs (n <= k)\n"
    "-t <num> ... set eigenvalue tolerance to num [Default 1e-6]\n";

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
    case 'k':
      if (*++argv[0]) kdim  = atoi (  *argv);
      else { --argc;  kdim  = atoi (*++argv); }
      break;
    case 'm':
      if (*++argv[0]) maxit = atoi (  *argv);
      else { --argc;  maxit = atoi (*++argv); }
      break;
    case 'n':
      if (*++argv[0]) neval = atoi (  *argv);
      else { --argc;  neval = atoi (*++argv); }
      break;
    case 't':
      if (*++argv[0]) evtol = atof (  *argv);
      else { --argc;  evtol = atof (*++argv); }
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
// Return length of an eigenvector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real* z;
  integer     i, np, nel, npert;

  file  = new FEML (session);
  mesh  = new Mesh (file);

  np    =  static_cast<integer>(Femlib::value ("N_POLY"));
  nel   = mesh -> nEl();
  npert = file -> attribute ("FIELDS", "NUMBER") - 1;
  Geometry::set (nel, npert);

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

  bman    = new BCmgr        (file, elmt);
  domain  = new Domain       (file, elmt, bman);
  analyst = new StabAnalyser (domain, file);

  // -- Over-ride any CHKPOINT flag in session file.

  Femlib::value ("CHKPOINT", 1);

  domain -> restart ();
  domain -> loadBase();
  domain -> report  ();

  return Geometry::nPert() * Geometry::planeSize() * Geometry::nZ();
}
