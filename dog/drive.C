///////////////////////////////////////////////////////////////////////////////
// drive.C: compute leading eigenvalues and eigenvectors for stability
// analysis based on linearised (Navier--Stokes) operators.
// Optionally compute the optimal transient growth initial condition,
// final condition, or the adjoint modes. The base flow can be either
// steady or periodic in time, two or three component, cylindrical or
// Cartesian, but must be two-dimensional.
// 
// Originally based on code "floK" by Dwight Barkley & Ron Henderson.
//
// Copyright (c) 1999 <--> $Date$, Hugh Blackburn
//
// The eigenpairs computed in the subspace are related to the Ritz
// estimates of those in the original space in a simple way: the
// eigenvalues are the same, and the Ritz eigenvectors are related to
// the subspace eigenvectors through a linear transformation (see
// Saad, p.175).
//
// USAGE
// -----
// dog [options] session
//   session: specifies name of semtex session file.
//   options:
//   -h       ... print this message
//   -v       ... set verbose
//   -a||g||s ... solve adjoint or optimal growth or optimal shrink problem
//   -k <num> ... set dimension of subspace (maximum number of pairs) to num
//   -m <num> ... set maximum number of iterations         (m >= k)
//   -n <num> ... compute num eigenvalue/eigenvector pairs (n <= k)
//   -t <num> ... set eigenvalue tolerance to num [Default = 1e-6].
//
#ifdef FLIP
//
// Floquet analysis for reflection/time translation (RT) symmetric
// base flows.  The idea here is that the mapping of an instability
// applied by RT-symmetric flows with period T is like the square of
// two mappings of period T/2 (see [4]). Here we explicitly deal only
// with the 1/2-period map, but have to apply a symmetry operation to
// the perturbations at the end of every 1/2-period before Krylov
// analysis. The mapping vector is precomputed by flipmap.C, which
// contains the rule for transforming the perturbation velocity field.
//
// The user has to ensure that the integration time in the session
// file is T/2 (just as, for dog, it must be T).
//
// Usage is the same as for dog, except the code is called dog-H.
// (H is for "half-period".)
//
#endif
// 
// FILES
// -----
// A number of semtex files are required --
//   session:     semtex session file
//   session.num: computed automatically if not supplied
//   session.bse: base flow, containing N_SLICE field dumps
//   session.rst: restart file (initialised with white noise if not supplied)
#ifdef FLIP
//   session.map: contains the transformation to be applied (see above)
#endif
//
// REFERENCES
// ----------
// [1]  D Barkley & RD Henderson (1996), "Three-dimensional Floquet
//      stability analysis of the wake of a circular cylinder",
//      J Fluid Mech V322, 215--241.
// [2]  Y Saad (1991), "Numerical methods for large eigenvalue problems"
//      Wiley.
// [3]  LS Tuckerman & D Barkley (2000), "Bifurcation analysis for
//      timesteppers", in Numerical Methods for Bifurcation Problems,
//      ed E Doedel & LS Tuckerman, Springer. 453--466.
// [4]  HM Blackburn, F Marques & JM Lopez (2005), "Symmetry breaking of 
//      two-dimensional time-periodic wakes", J Fluid Mech V522, 395--411.
// [5]  D Barkley, HM Blackburn & SJ Sherwin (2008), "Direct optimal
//      growth analysis for timesteppers", IJNMF V57, 1435--1458.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <stab.h>

#ifdef FLIP
static char             prog[] = "dog-H";
static char             generator;
static vector<int_t>    positive, negative;
static void loadmap     (const char*);
static void mirror      (real_t*);
#else
static char             prog[] = "dog";
#endif
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void  getargs    (int, char**, problem_t&, int_t&, int_t&, int_t&,
			 int_t&, real_t&, char*&);
static int_t preprocess (const char*, bool&);

static void  EV_init    (real_t*);
static void  EV_update  (const problem_t, const real_t*, real_t*);
static void  EV_small   (real_t**, const int_t, const real_t*, const int_t,
			 real_t*, real_t*, real_t*, real_t&, const int_t, 
			 ofstream&); 
static int_t EV_test    (const int_t, const int_t, real_t*, real_t*, real_t*,
			 const real_t, const real_t, const int_t, ofstream&);
static void  EV_sort    (real_t*, real_t*, real_t*, real_t*, const int_t);
static void  EV_post    (const problem_t, real_t**, real_t**,
			 const int_t, const int_t, const int_t,
			 const real_t*, const real_t*, 
			 const real_t*, const int_t, ofstream&);
static void  EV_big     (const problem_t, real_t**, real_t**,
			 const int_t, const int_t, const int_t,
			 const real_t*, const real_t*,
			 const real_t*, ofstream&);

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  int_t     kdim = 2, nvec = 2, nits = 2, verbose = 0, converged = 0;
  real_t    resnorm, evtol = 1.0e-6;
  int_t     i, j, k;
  char      buf[StrMax];
  ofstream  runinfo;
  problem_t task = PRIMAL;
  bool      restart = false;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, task, kdim, nits, nvec, verbose, evtol, session);

  // -- Check, echo parameter values.

  if (kdim < 1)    message (prog, "param error: KDIM must be > 1",     ERROR);
  if (nvec < 1)    message (prog, "param error: NVEC must be > 1",     ERROR);
  if (nits < kdim) message (prog, "param error: NITS must be >= KDIM", ERROR);
  if (kdim < nvec) message (prog, "param error: NVEC must be <= KDIM", ERROR);

  // -- Install lookup copies for reporting purposes.

  Femlib::ivalue ("KRYLOV_KDIM", kdim);
  Femlib::ivalue ("KRYLOV_NVEC", nvec);
  Femlib::ivalue ("KRYLOV_NITS", nits);
  Femlib::value  ("KRYLOV_KTOL", evtol);

  strcat (strcpy (buf, session), ".evl");
  runinfo.open (buf, ios::out);

  // -- Set up to run with semtex.
  
  const int_t ntot = preprocess (session, restart);

#if defined (ARPACK)		// -- Eigensolution by ARPACK.

  const int_t done = 99, lworkl = 3*kdim*kdim + 6*kdim;
  int_t       ido, info, iparam[11], ipntr[14], select[kdim];

  iparam [0] = 1;		// -- Shifting will be handled by ARPACK.
  iparam [1] = 0; 		// -- Not used.
  iparam [2] = nits;		// -- Input: maximum, output: number done.
  iparam [3] = 1;		// -- Blocksize, ARPACK say = 1.
  iparam [4] = 0;		// -- Output, number of converged values.
  iparam [5] = 0;		// -- Not used.
  iparam [6] = 1;		// -- Mode: solve A x = lambda x.
  iparam [7] = 0; 		// -- For user shifts, not used here.
  iparam [8] = 0;		// -- Output, number of Op x operations.
  iparam [9] = 0;		// -- Output, not used here.
  iparam[10] = 0;		// -- Output, number of re-orthog steps.

  // -- Allocate and zero storage.

  vector<real_t> work(3*ntot + lworkl + ntot*kdim + ntot +
		      2*(nvec+1) + 3*kdim + ntot*(nvec+1));
  work.clear();

  real_t*        workd  = &work[0];
  real_t*        workl  = workd + 3*ntot;
  real_t*        v      = workl + lworkl;
  real_t*        resid  = v + ntot*kdim;
  real_t*        dr     = resid + ntot;
  real_t*        di     = dr + nvec + 1;
  real_t*        workev = di + nvec + 1;
  real_t*        z      = workev + 3*kdim;

  // -- Either read in a restart, or set random IC. 

  EV_init (resid);

  // -- Set up for reverse communication.

  F77NAME(dnaupd) (ido=0, "I", ntot, "LM", nvec, evtol, resid, kdim, 
		   v, ntot, iparam, ipntr, workd, workl, lworkl, info=1);

  // -- IRAM iteration.

  i = 0;
  while (ido != done) {
    EV_update (task, workd+ipntr[0]-1, workd+ipntr[1]-1);

    F77NAME(dnaupd) (ido, "I", ntot, "LM", nvec, evtol, resid, kdim,
		     v, ntot, iparam, ipntr, workd, workl, lworkl, info);

    runinfo << "ARPACK iteration: " << ++i << endl;
  }

  if (info < 0) message (prog, "ARPACK error", ERROR);

  runinfo << "Converged in " << iparam[8] << " iterations" << endl;

  // -- Post-process to obtain eigenvalues and Ritz eigenvectors.

  F77NAME(dneupd) (1, "A", select, dr, di, z, ntot, 0, 0, workev,
		   "I", ntot, "LM", nvec, evtol, resid, kdim,
		   v, ntot, iparam, ipntr, workd, workl, lworkl, info);

  // -- Print up eigenvalues.

  real_t       re_ev, im_ev, abs_ev, ang_ev, re_Aev, im_Aev;
  const real_t period = Femlib::value ("D_T * N_STEP");

  runinfo.precision(4);
  runinfo.setf(ios::scientific, ios::floatfield);

  runinfo << "EV  Magnitude   Angle       Growth      Frequency" << endl;

  for (j = 0; j < nvec; i++) {
    re_ev  = dr[j];
    im_ev  = di[j];
    abs_ev = hypot (re_ev, im_ev);
    ang_ev = atan2 (im_ev, re_ev);
    re_Aev = log (abs_ev) / period;
    im_Aev = ang_ev       / period;
    runinfo << setw(2)  << j
         << setw(12) << abs_ev
      	 << setw(12) << ang_ev
	 << setw(12) << re_Aev
	 << setw(12) << im_Aev
	 << endl;
  }

  // -- Print up eigenvectors.

  for (j = 0; j < nvec; j++) {
    char     msg[StrMax], nom[StrMax];
    real_t*  src = z + j * ntot;
    ofstream file;
    for (i = 0; i < Geometry::nPert(); i++)
      for (k = 0; k < Geometry::nZ(); k++)
	domain -> u[i] -> setPlane 
	  (k, src + (i*Geometry::nZ() + k)*Geometry::planeSize());
    sprintf   (msg, ".eig.%1d", j);
    strcat    (strcpy (nom, domain -> name), msg);
    file.open (nom, ios::out); file << *domain; file.close();
  }


#else                           // -- Eigensolution by DB algorithm.

  // -- Allocate and zero eigenproblem storage.

  const int_t wdim = kdim + kdim + kdim*kdim + (2*ntot + 1)*(kdim + 1);

  vector<real_t> work (wdim);
  work.clear();

  real_t*  alpha = &work[0];	             // -- Scale factors.
  real_t*  wr    = alpha + (kdim + 1);       // -- Eigenvalues (real part).
  real_t*  wi    = wr   + kdim;	             //                (imag part).
  real_t*  zvec  = wi   + kdim;	             // -- Subspace eigenvectors.
  real_t*  kvec  = zvec + kdim * kdim;       // -- Krylov sequence (flat).
  real_t*  tvec  = kvec + ntot * (kdim + 1); // -- Orthonormalised equivalent.
  real_t** Kseq  = new real_t* [kdim + 1];   // -- Handles for Krylov sequence.
  real_t** Tseq  = new real_t* [kdim + 1];   // -- Handles for ortho  sequence.

  for (i = 0; i <= kdim; i++) {
    Kseq[i] = kvec + i * ntot;
    Tseq[i] = tvec + i * ntot;
  }

  // -- Load starting vector.

  EV_init (Kseq[0]);

  if (!restart) {

    // -- If we had a random initial condition, apply operator
    //    once to enforce BCs, incompressibility, etc.

    EV_update (task, Kseq[0], Kseq[1]);
    Veclib::copy (ntot, Kseq[1], 1, Kseq[0], 1);
  }

  alpha[0] = sqrt (Blas::nrm2 (ntot, Kseq[0], 1));
  Blas::scal (ntot, 1.0/alpha[0], Kseq[0], 1);

  // -- Fill initial Krylov sequence, during which convergence may occur.

  for (i = 1; !converged && i <= kdim; i++) {

    EV_update (task, Kseq[i - 1], Kseq[i]);

    alpha[i] = sqrt (Blas::nrm2 (ntot, Kseq[i], 1));
    Blas::scal (ntot, 1.0/alpha[i], Kseq[i], 1);

    Veclib::copy (ntot * (i + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, alpha, i, zvec, wr, wi, resnorm, verbose, runinfo);
    converged = EV_test (i,i,zvec,wr,wi,resnorm,evtol,min(i, nvec), runinfo);
    converged = max (converged, 0); // -- Only exit on evtol.
  }

  // -- Carry out iterative solution with a full sequence.

  for (i = kdim + 1; !converged && i <= nits; i++) {

    // -- Update Krylov sequence. 

    for (j = 1; j <= kdim; j++) {
      alpha[j - 1] = alpha[j];
      Veclib::copy (ntot, Kseq[j], 1, Kseq[j - 1], 1);
    }

    EV_update (task, Kseq[kdim - 1], Kseq[kdim]);

    // -- Compute new scale factor.

    alpha[kdim] = sqrt (Blas::nrm2 (ntot, Kseq[kdim], 1));
    Blas::scal (ntot, 1.0/alpha[kdim], Kseq[kdim], 1);
    
    // -- Get subspace eigenvalues, test for convergence.

    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, alpha, kdim, zvec, wr, wi, resnorm, verbose,runinfo);

    converged = EV_test (i, kdim, zvec, wr, wi, resnorm, evtol, nvec, runinfo);
  }

  EV_post(task,Tseq,Kseq,ntot,min(--i,kdim),nvec,zvec,wr,wi,converged,runinfo);

#endif

  runinfo.close();
  Femlib::finalize();
  return (EXIT_SUCCESS);
}


static void EV_init (real_t* tgt)
// ---------------------------------------------------------------------------
// Load initial vector from domain velocity fields.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t ND = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  const int_t NZ = Geometry::nZ();
    
  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);
}


static void EV_update  (const problem_t task,
                        const real_t*   src,
			real_t*         tgt)
// ---------------------------------------------------------------------------
// Generate tgt by applying linear operator (here, a linearised
// Navier--Stokes integrator) to src.  Src and tgt could be the same
// vector.
//
// If RT-flipping, apply flip-map.
// ---------------------------------------------------------------------------
{
  int_t       i, k;
  const int_t ND = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  const int_t NZ = Geometry::nZ();

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);

  switch (task) {
  case PRIMAL:			// -- Forward in time.
    integrate (linAdvect , domain, analyst); break;

  case ADJOINT:			// -- Backward in time.
    integrate (linAdvectT, domain, analyst); break;

  case GROWTH:			// -- Forward, then backward.
    integrate (linAdvect , domain, analyst);
    integrate (linAdvectT, domain, analyst); break;

  case SHRINK:			// -- Backward, then forward.
    integrate (linAdvectT, domain, analyst);
    integrate (linAdvect , domain, analyst); break;

  default:
    message ("EV_update", "Impossible task", ERROR); break;
  }

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++)
      domain -> u[i] -> getPlane (k, tgt + (i*NZ + k)*NP);

#ifdef FLIP
  mirror (tgt);
#endif
}


static void EV_small (real_t**      Kseq   ,
		      const int_t   ntot   ,
		      const real_t* alpha  ,
		      const int_t   kdim   ,
		      real_t*       zvec   ,
		      real_t*       wr     ,
		      real_t*       wi     ,
		      real_t&       resnorm,
		      const int_t   verbose,
		      ofstream&     runinfo   )
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
// destroyed in the process, and replaced by an orthonormal basis for
// Kseq (Q).
//
// Then we compute Hessenberg matrix H = Q* A Q (using Q* A = R), find
// its eigenvalues (wr, wi) and (right) eigenvectors, zvec, related to
// those of A.
//
// The residual norm for each eigenvector is related to R(kdim+1, kdim),
// which says how much the (k+1) component lies outside Q.
// ---------------------------------------------------------------------------
{
  char           routine[] = "EV_small";
  const int_t    kdimp = kdim + 1;
  int_t          i, j, ier, lwork = 10 * kdim;
  vector<real_t> work (kdimp * kdimp + kdim * kdim + lwork);
  real_t         *R     = &work[0],
                 *H     = R + kdimp * kdimp,
                 *rwork = H + kdim * kdim;

  Veclib::zero (kdimp * kdimp, R, 1);

  // -- Modified G--S orthonormalisation.

  for (i = 0; i < kdimp; i++) {
    real_t gsc = Blas::nrm2 (ntot, Kseq[i], 1);
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
    runinfo << "R =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	runinfo << setw (14) << R[Veclib::col_major (i, j, kdimp)];
      runinfo << endl;
    }
  }

  // -- H(i, j) = (q_i, A q_j) 
  //            = 1 / R(j, j) * (alpha(j + 1)*R(i, j + 1) - H(i, l).R(l, j)),
  //    with the last inner product taken over l < j.

  for (i = 0; i < kdim; i++) {
    for (j = 0; j < kdim; j++) {
      H[Veclib::col_major (i, j, kdim)] =
	alpha[j + 1] * R[Veclib::col_major (i, j + 1, kdimp)]
	- Blas::dot (j, H + i, kdim, R + j * kdimp, 1);
      H[Veclib::col_major (i, j, kdim)] /= R [Veclib::col_major (j, j, kdimp)];
    }
  }

  // -- Print up H as diagnostic.

  if (verbose) {
    runinfo << "H =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	runinfo << setw (14) << H[Veclib::col_major (i, j, kdim)];
      runinfo << endl;
    }
  }

  // -- Find eigenpairs of H using LAPACK routine dgeev.
  //    For complex-conjugate eigenvalues, the corresponding
  //    complex eigenvector is a real-imaginary pair of rows
  //    of zvec.

  F77NAME(dgeev) ("N","V",kdim,H,kdim,wr,wi,0,1,zvec,kdim,rwork,lwork,ier);

  if (ier) message (routine, "error return from dgeev", ERROR);

  // -- Print up (unsorted) eigenvalues and eigenvectors as diagnostic.

  if (verbose) {
    runinfo << "eval =" << endl;
    for (i = 0; i < kdim; i++) runinfo << setw (14) << wr[i]; runinfo << endl;
    for (i = 0; i < kdim; i++) runinfo << setw (14) << wi[i]; runinfo << endl;
    runinfo << "zvec =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	runinfo << setw (14) << zvec[Veclib::col_major (i, j, kdim)];
      runinfo << endl;
    }
  }
  
  // -- Compute residual information.

  resnorm = alpha[kdim] * fabs 
    (R[Veclib::col_major (kdim,     kdim,     kdimp)] /
     R[Veclib::col_major (kdim - 1, kdim - 1, kdimp)] );
}


static int_t EV_test (const int_t  itrn   ,
		      const int_t  kdim   ,
		      real_t*      zvec   ,
		      real_t*      wr     ,
		      real_t*      wi     ,
		      const real_t resnorm,
		      const real_t evtol  ,
		      const int_t  nvec   ,
		      ofstream&    runinfo)
// ---------------------------------------------------------------------------
// Test convergence of eigenvalues and print up diagnostic information.
//
// Return value:
//   nvec:  all of the first nvec eigenvalue estimates have converged;
//  -1,-2:  the residuals aren't shrinking;
//   0:     neither of the above is true: not converged.
//
// Hold off on testing residual shrinkage until we're in steady-state
// operation.
// ---------------------------------------------------------------------------
{
  int_t          i, idone = 0;
  vector<real_t> work (kdim);
  real_t         re_ev, im_ev, abs_ev, ang_ev, re_Aev, im_Aev;
  real_t*        resid = &work[0];
  const real_t   period = Femlib::value ("D_T * N_STEP");
  static real_t  min_max1, min_max2;
 
  if (min_max1 == 0.0) min_max1 = 1000.0;
  if (min_max2 == 0.0) min_max2 = 1000.0;

  // -- Sort subspace eigenvectors by residual.

  for (i = 0; i < kdim; i++) {
    resid[i] = resnorm *
      fabs (zvec[kdim - 1 + i*kdim]) / Blas::nrm2 (kdim, zvec + i*kdim, 1);
    if (wi[i] < 0.0) resid[i - 1] = resid[i] = hypot (resid[i - 1], resid[i]);
  }
  EV_sort (zvec, wr, wi, resid, kdim);

  // -- Stopping test.

  if      (resid[nvec - 1] < evtol * hypot (wr[0], wi[0]))       idone = nvec;
  else if (min_max1 < 0.01 && resid[nvec - 1] > 10.0 * min_max1) idone = -1;
  else if (min_max2 < 0.01 && resnorm         > 10.0 * min_max2) idone = -2;

  min_max1 = (itrn > kdim) ? min (min_max1, resid[nvec - 1]) : min_max1;
  min_max2 = (itrn > kdim) ? min (min_max2, resnorm)         : min_max2;

  // -- Print diagnostic information.

  runinfo << "-- Iteration = " << itrn << ", H(k+1, k) = " << resnorm << endl;

  runinfo.precision(4);
  runinfo.setf(ios::scientific, ios::floatfield);

  runinfo << "EV  Magnitude   Angle       Growth      Frequency   Residual"
       << endl;

  for (i = 0; i < kdim; i++) {
    re_ev  = wr[i];
    im_ev  = wi[i];
    abs_ev = hypot (re_ev, im_ev);
    ang_ev = atan2 (im_ev, re_ev);
    re_Aev = log (abs_ev) / period;
    im_Aev = ang_ev       / period;
    runinfo << setw(2)  << i
         << setw(12) << abs_ev
      	 << setw(12) << ang_ev
	 << setw(12) << re_Aev
	 << setw(12) << im_Aev
	 << setw(12) << resid[i]
	 << endl;
  }

  runinfo.precision(6);
  runinfo.setf(ios::fixed);
  
  return idone;
}


static void EV_sort (real_t*     evec,
		     real_t*     wr  ,
		     real_t*     wi  ,
		     real_t*     test,
		     const int_t dim )
// ---------------------------------------------------------------------------
// Insertion sort to rearrange eigenvalues and eigenvectors to ascending
// order according to vector test.  See equivalent Numerical Recipes routine.
// ---------------------------------------------------------------------------
{
  int_t          i, j;
  vector<real_t> work (dim);
  real_t         wr_tmp, wi_tmp, te_tmp, *z_tmp = &work[0];

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


static void EV_post (const problem_t task,
		     real_t**        Tseq,
		     real_t**        Kseq, 
		     const int_t     ntot, 
		     const int_t     kdim, 
		     const int_t     nvec,
		     const real_t*   zvec, 
		     const real_t*   wr  , 
		     const real_t*   wi  , 
		     const int_t     icon,
		     ofstream&       runinfo)
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
  const char*   routine = "EV_post";
  char          msg[StrMax], nom[StrMax];
  int_t         i, j, k;
  const int_t   ND = Geometry::nPert();
  const int_t   NP = Geometry::planeSize();
  const int_t   NZ = Geometry::nZ();
  const real_t* src;
  ofstream      file;

  if (icon == 0) {

    runinfo << prog
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
    
    runinfo << prog
	 << ": minimum residual reached, writing initial Krylov vector."
	 << endl;
    
    src = Kseq[0];
    for (i = 0; i < ND; i++)
      for (k = 0; k < NZ; k++)
	domain -> u[i] -> setPlane (k, src + (i*NZ + k)*NP);
    strcat    (strcpy (nom, domain -> name), ".fld");
    file.open (nom, ios::out); file << *domain; file.close();

  } else if (icon == nvec) {

    runinfo << prog
	 << ": converged, writing "
	 << icon
	 << " eigenvectors."
	 << endl;
    
    EV_big (task, Tseq, Kseq, ntot, kdim, icon, zvec, wr, wi, runinfo);
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


static void EV_big (const problem_t task ,
		    real_t**        bvecs,
	            real_t**        evecs,
		    const int_t     ntot ,
		    const int_t     kdim ,
		    const int_t     nvec ,
		    const real_t*   zvecs,
		    const real_t*   wr   ,
		    const real_t*   wi   ,
		    ofstream&       runinfo )
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
// computations are only done for the first nvec eigenvectors, but use
// all kdim vectors in the Krylov sequence.  If this is not yet full,
// that's OK because the sequence was zeroed at allocation.
// ---------------------------------------------------------------------------
{
  real_t norm, resid, wgt;
  int_t  i, j;

  // -- Generate big e-vectors.

  for (j = 0; j < nvec; j++) {
    Veclib::zero (ntot, evecs[j], 1);
    for (i = 0; i < kdim; i++) {
      wgt = zvecs[Veclib::col_major (i, j, kdim)];
      Blas::axpy (ntot, wgt, bvecs[i], 1, evecs[j], 1);
    }
  }

  // -- Normalize big e-vectors.

  for (i = 0; i < nvec; i++)
    if (wi[i] == 0.0) {		// -- This a real mode.
      norm = Blas::nrm2 (ntot, evecs[i], 1);
      Blas::scal (ntot, 1.0/norm, evecs[i], 1);
    } else {			// -- It's complex.
      norm  = sqr (Blas::nrm2 (ntot, evecs[i],   1));
      norm += sqr (Blas::nrm2 (ntot, evecs[i+1], 1));
      norm = sqrt (norm);
      Blas::scal (ntot, 1.0/norm, evecs[i],   1);
      Blas::scal (ntot, 1.0/norm, evecs[i+1], 1);
      i++;
    }

  if (Femlib::ivalue ("BIG_RESIDS") > 0) {

    // -- Compute residuals of big eigenvectors directly, i.e. A(u)-lambda(u).

    for (i = 0; i < nvec; i++) {

      EV_update (task, evecs[i], bvecs[0]);

      if (wi[i] == 0.0)
	Blas::axpy (ntot, -wr[i], evecs[i], 1, bvecs[0], 1);
      else if(wi[i] > 0.) {
	Blas::axpy (ntot, -wr[i], evecs[i],   1, bvecs[0], 1);
	Blas::axpy (ntot,  wi[i], evecs[i+1], 1, bvecs[0], 1);
      } else {
	Blas::axpy (ntot, -wr[i], evecs[i],   1, bvecs[0], 1);
	Blas::axpy (ntot,  wi[i], evecs[i-1], 1, bvecs[0], 1);
      }

      resid = Blas::nrm2 (ntot, bvecs[0], 1) / Blas::nrm2 (ntot, evecs[i], 1);
    
      runinfo << "Big-space residual of eigenvector " <<i<< ": " << resid << endl;
    }
  }
}


static void getargs (int        argc   ,
		     char**     argv   ,
		     problem_t& task   ,
		     int_t&     kdim   , 
		     int_t&     maxit  ,
		     int_t&     neval  ,
		     int_t&     verbose,
		     real_t&    evtol  ,
		     char*&     session)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "dsa(-H) [options] session\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-a||g||s ... solve adjoint or optimal growth or optimal shrink, problem\n"
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
      Femlib::ivalue ("VERBOSE", Femlib::ivalue("VERBOSE") + 1);
      break;
    case 'a':
      task = ADJOINT;
      break;
    case 'g':
      task = GROWTH;
      break;
    case 's':
      task = SHRINK;
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

  // -- Here is a minor hack, installs TASK in parser.

  Femlib::ivalue ("TASK", task);

  // -- While Fourier temporal interpolation is the default for base
  //    flow reconstruction, we can switch to 4-point (cubic) Lagrange
  //    interpolation by setting token LAGRANGE_INT. Here we install
  //    it in the parser table but set it to be disabled (0).

  Femlib::ivalue ("LAGRANGE_INT", 0);
}


static int_t preprocess (const char* session,
			 bool&       restart)
// ---------------------------------------------------------------------------
// Create objects needed for semtex execution, given the session file name.
//
// Return length of an eigenvector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real_t* z;
  int_t         i, np, nel, npert;

  // -- Set default additional tokens.

  Femlib::value  ("BASE_PERIOD", 0.0);
  Femlib::value  ("T_OFFSET",    0.0);
  Femlib::ivalue ("BIG_RESIDS",  0);

  // -- Start up dealing with session file.

  file  = new FEML (session);
  mesh  = new Mesh (file);

  np    = Femlib::ivalue ("N_P");
  nel   = mesh -> nEl();
  npert = file -> attribute ("FIELDS", "NUMBER") - 1;
  Geometry::set (nel, npert);

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  bman   = new BCmgr  (file, elmt);
  domain = new Domain (file, elmt, bman);

  // -- Load restart and base flow data.

  restart = domain -> restart ();
  domain -> loadBase();
  domain -> report  ();

  analyst = new StabAnalyser (domain, file);

  // -- Over-ride any CHKPOINT flag in session file.

  Femlib::value ("CHKPOINT", 1);

#ifdef FLIP
  // -- Load RT-symmetry map.

  loadmap (session);
#endif

  return Geometry::nPert() * Geometry::planeSize() * Geometry::nZ();
}


#ifdef FLIP
static void loadmap (const char* session)
// ---------------------------------------------------------------------------
// Load symmetry mapping information from session.map.
// ---------------------------------------------------------------------------
{
  const int_t np  = Femlib::ivalue ("N_P");
  const int_t nel = mesh -> nEl();
  char        buf[StrMax], err[StrMax];
  ifstream    file;
  int_t       i, NR, NS, NEL, NMAP;
  
  file.open (strcat (strcpy (buf, session), ".map"));

  if (!file) {
    sprintf (err, "cannot find map file %s", buf);
    message (prog, err, ERROR);
  }

  file >> NR >> NS >> NEL >> NEL;
  file.ignore (StrMax, '\n');

  if (NR != np || NS != np || NEL != nel)
    message (prog, "map file doesn't conform with session file", ERROR);
  file >> generator;
  file >> NMAP;

  cout
    << "   RT-flip mapping         : "
    << NMAP
    << " points by "
    << generator
    << " reflection" << endl;
  
  positive.resize (NMAP);
  negative.resize (NMAP);

  for (i = 0; i < NMAP; i++) file >> positive[i] >> negative[i];

  if (!file)
    message (prog, "bad (premature end of?) map file", ERROR);

  file.close();
}


static void mirror (real_t* tgt)
// ---------------------------------------------------------------------------
// Apply RT-flip-map. Note that to avoid holes in the mapping, the
// gather and scatter vectors have to be used in the order shown.
// ---------------------------------------------------------------------------
{
  int_t          i, k;
  const int_t    ND = Geometry::nPert();
  const int_t    NP = Geometry::planeSize();
  const int_t    NZ = Geometry::nZ();
  const int_t    NM = positive.size();
  static real_t* tmp;

  if (!tmp) tmp = new real_t [NP];

  // -- First, the reflection.

  for (i = 0; i < ND; i++)
    for (k = 0; k < NZ; k++) {
      Veclib::copy (NP, tgt + (i*NZ+k)*NP, 1, tmp, 1);
      Veclib::gathr_scatr (NM,tmp,&negative[0],&positive[0],tgt + (i*NZ+k)*NP);
    }
  
  // -- Then the sign change.

  if (generator == 'x')		// -- Change sign of 'u'.
    for (k = 0; k < NZ; k++)
      Veclib::neg (NP, tgt + (0*NZ+k)*NP, 1);
  else				// -- Change sign of 'v'.
    for (k = 0; k < NZ; k++)
      Veclib::neg (NP, tgt + (1*NZ+k)*NP, 1);

  // -- Update simulation time by T/2.

  domain -> time += 0.5 * domain -> period;
}
#endif
