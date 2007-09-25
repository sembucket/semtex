///////////////////////////////////////////////////////////////////////////////
// arnoldi.C: compute leading eigenvalues and eigenvectors of (real)
// arbitrary sparse matrix using Arnoldi iteration.  Sparse matrices
// are read and stored internally using Harwell--Boeing (HB) format.
// 
// Based on code floK by Dwight Barkley
//
// USAGE
// -----
// arnoldi [options] matfile
//   matfile: specifies name of HB format matrix input file
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
// Input is obtained from HB format matrix file.  Length of eigenvectors
// is set to match the size of matrix specified in file.
// 
// Output is written in ASCII to stdout.  The output consists of eigen pairs,
// in blocks separated by a blank line.  Each block commences with the real
// and imaginary parts of the eigenvalue, followed by the eigenvector 
// components.
//
// REFERENCES
// ----------
// [1]  D. Barkley & R.D. Henderson (1996), "Three-dimensional Floquet
//      stability analysis of the wake of a circular cylinder",
//      J. Fluid Mech V322, 215--241.
// [2]  Y. Saad (1991), "Numerical methods for large eigenvalue problems"
//      Wiley.
// [3]  W.J. Stewart & A. Jennings, (1981), "A simultaneous iteration
//      algorithm for real matrices", ACM Trans Math Soft V7N2, 184--198.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <cstdarg>		/* System C headers.  */
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cctype>
#include <cstring>
#include <climits>
#include <cfloat>
#include <cassert>

#include <iostream>		/* System C++ headers. */
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

#include <cfemdef.h>		// -- Semtex headers.
#include <blas.h>
#include <lapack.h>
#include <utility.h>
#include <veclib.h>

#include <iohb.h>		// -- HB IO routines from NIST.

#define F77name(x) x ## _

extern "C" {
  void F77name(dspmvc)		// -- HB matrix-vector product.
    (const int_t&   trans ,
     const int_t&   n     ,
     const int_t&   m     ,
     const double*  a     ,
     const int_t*   colptr,
     const int_t*   rowind,
     const double*  x     ,
     const int_t&   ldx   ,
     double*        y     ,
     const int_t&   ldy   );

#if defined (ARPACK)
  void F77name(dnaupd) 		// -- ARPACK reverse-communications interface.
    (int_t&         ido   ,
     const char*    bmat  ,
     const int_t&   n     ,
     const char*    which ,
     const int_t&   nev   ,
     const real_t&  tol   ,
     real_t*        resid ,
     const int_t&   ncv   ,
     real_t*        v     ,
     const int_t&   ldv   ,
     int_t*         iparam,
     int_t*         ipntr ,
     real_t*        workd ,
     real_t*        workl ,
     const int_t&   lworkl,
     int_t&         info
     );
#endif
}

static char prog[] = "arnoldi";

static void  getargs  (int, char**,int_t&,int_t&,int_t&,int_t&,real_t&,char*&);
static void  EV_small (real_t**, const int_t, const real_t*, const int_t, 
		       real_t*, real_t*, real_t*, real_t&, const int_t); 
static int_t EV_test  (const int_t, const int_t, real_t*, real_t*, real_t*,
		       const real_t, const real_t, const int_t);
static void  EV_sort  (real_t*, real_t*, real_t*, real_t*, const int_t);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char   *HBfile, *HBtype;
  int_t  *HBrptr, *HBcptr;
  real_t *HBval;
  int_t   HBnr, HBnc, HBnz, HBnrhs;

  int_t   ntot, kdim = 2, nvec = 2, nits = 2, verbose = 0;
  real_t  evtol = 1.0e-6;

  getargs (argc, argv, kdim, nits, nvec, verbose, evtol, HBfile);

  // -- Check parameter values.

  if (kdim < 1)    message (prog, "param error: KDIM must be > 1",     ERROR);
  if (nvec < 1)    message (prog, "param error: NVEC must be > 1",     ERROR);
  if (nits < kdim) message (prog, "param error: NITS must be >= KDIM", ERROR);
  if (kdim < nvec) message (prog, "param error: NVEC must be <= KDIM", ERROR);

  // -- Open Harwell--Boeing file, extract data.

  readHB_info (HBfile, &HBnr, &HBnc, &HBnz, &HBtype, &HBnrhs);
  if (verbose) {
    cerr << "Matrix in file " << HBfile << "is " << HBnr << " x " << HBnc
	 << " with " << HBnz << " nozero entries, type " << HBtype << endl;
    cerr << HBnrhs << " right-hand-sides available" << endl;
  }
  readHB_newmat_double (HBfile, &HBnr, &HBnc, &HBnz, &HBcptr, &HBrptr, &HBval);

  ntot = HBnr;

#if defined (ARPACK)		// -- Solution using ARPACK.

  const int_t done = 99, lworkl = 3*kdim*kdim + 6*kdim;
  int_t       ido, info, iparam[11], ipntr[14];

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

  // -- Allocate storage.

  vector<real_t> work (3*ntot + lworkl + ntot*kdim + ntot);
  real_t*        workd = &work[0];
  real_t*        workl = workd + 3*ntot;
  real_t*        v     = workl + lworkl;
  real_t*        resid = v + ntot*kdim;

  // -- Set up for reverse communication.

  cout << "set up ... " << flush; 

  F77name(dnaupd) (ido=0, "I", ntot*ntot, "LM", nvec, evtol, resid, kdim, 
		   v, ntot, iparam, ipntr, workd, workl, lworkl, info);

  cout << "done" << endl;

  // -- IRAM iteration.

  while (ido != done) {

    cout << "call Aop ... " << flush ;
    F77name(dspmvc) (0, HBnr, 1, HBval, HBcptr, HBrptr,
		     workd+ipntr[0]-1, HBnr, workd+ipntr[1]-1, HBnr);

    cout << "done" << endl;
    
    F77name(dnaupd) (ido, "I", ntot*ntot, "LM", nvec, evtol, resid, kdim,
		     v, ntot, iparam, ipntr, workd, workl, lworkl, info);

    cout << "Info: " << info << "  resid: " << resid << endl;
  }

#else                           // -- Solution using Dwight's algorithm.

  // -- Allocate eigenproblem storage.
  
  int_t           i, j, converged = 0;
  real_t          resnorm;
  int_t           nwork  = kdim + kdim + kdim*kdim + (2*ntot + 1)*(kdim + 1);

  vector<real_t>  work (nwork);
  Veclib::zero (nwork, &work[0], 1);

  real_t*         alpha = &work[0];
  real_t*         wr    = alpha + (kdim + 1);
  real_t*         wi    = wr   + kdim;
  real_t*         zvec  = wi   + kdim;
  real_t*         kvec  = zvec + kdim * kdim;
  real_t*         tvec  = kvec + ntot * (kdim + 1);
  real_t**        Kseq  = new real_t* [kdim + 1];
  real_t**        Tseq  = new real_t* [kdim + 1];

  for (i = 0; i <= kdim; i++) {
    Kseq[i] = kvec + i * ntot;
    Tseq[i] = tvec + i * ntot;
  }

  // -- Generate random initial guess, zero evals.

  Veclib::vrandom (ntot, Kseq[0], 1);
  alpha[0] = sqrt(Blas::nrm2 (ntot, Kseq[0], 1));
  Blas::scal (ntot, 1.0/alpha[0], Kseq[0], 1);

  // -- Fill initial Krylov sequence.

  for (i = 1; !converged && i <= kdim; i++) {

    F77NAME(dspmvc) (0, HBnr, 1, HBval, HBcptr, HBrptr,
		     Kseq[i - 1], HBnr, Kseq[i], HBnr);

    alpha[i] = sqrt(Blas::nrm2(ntot, Kseq[i], 1));
    Blas::scal (ntot, 1.0/alpha[i], Kseq[i], 1);

    Veclib::copy (ntot * (i + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, alpha, i, zvec, wr, wi, resnorm, verbose);
    converged = EV_test (i, i, zvec, wr, wi, resnorm, evtol, min(i, nvec));
    converged = max (converged, 0);
  }

  // -- Carry out iterative solution.

  for (i = kdim + 1; !converged && i <= nits; i++) {

    for (j = 1; j <= kdim; j++) {
      alpha[j - 1] = alpha[j];
      Veclib::copy (ntot, Kseq[j], 1, Kseq[j - 1], 1);
    }

      // -- Matrix-vector product.

    F77name(dspmvc) (0, HBnr, 1, HBval, HBcptr, HBrptr,
		     Kseq[kdim - 1], HBnr, Kseq[kdim], HBnr);

    alpha[kdim] = sqrt (Blas::nrm2 (ntot, Kseq[kdim], 1));
    Blas::scal (ntot, 1.0/alpha[kdim], Kseq[kdim], 1);

    // -- Get subspace eigenvalues, test for convergence.

    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, alpha, kdim, zvec, wr, wi, resnorm, verbose); 

    converged = EV_test (i, kdim, zvec, wr, wi, resnorm, evtol, nvec);
  }

  if      (!converged)
    message (prog, "not converged", ERROR);
  else if (converged == nvec)
    message (prog, ": all estimates converged",  REMARK);
  else
    message (prog, ": minimum residual reached", REMARK);

#endif

  return (EXIT_SUCCESS);
}


static void getargs (int     argc,
		     char**  argv ,
		     int_t&  kdim , 
		     int_t&  maxit,
		     int_t&  neval,
		     int_t&  verb ,
		     real_t& evtol,
		     char*&  mfile)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "arnoldi [options] matfile\n"
    "matfile: specifies name of HB format matrix input file\n"
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
      verb = 1;
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

  if   (argc != 1) message (prog, "no matrix file",   ERROR);
  else             mfile = *argv;
}


static void EV_small (real_t**      Kseq   ,
		      const int_t   ntot   ,
		      const real_t* alpha  , 
		      const int_t   kdim   ,
		      real_t*       zvec   ,
		      real_t*       wr     ,
		      real_t*       wi     ,
		      real_t&       resnorm,
		      const int_t   verbose)
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
// destroyed in the process.
//
// Then we compute Hessenberg matrix H = Q* A Q (using Q* A = R), find
// its eigenvalues (wr, wi) and eigenvectors, zvec, related to those of A.
//
// The residual norm for each eigenvector is related to H(kdim+1, kdim),
// which is passed back for convergence testing.
// ---------------------------------------------------------------------------
{
  char           routine[] = "EV_small";
  const int_t    kdimp = kdim + 1;
  int_t          i, j, ier, lwork = 10 * kdim;
  vector<real_t> work (kdimp * kdimp + kdim * kdim + lwork);
  real_t         *R = &work[0], *H = R + kdimp*kdimp, *rwork = H + kdim*kdim;

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
    cout << "R =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	cout << setw (14) << R[Veclib::col_major (i, j, kdimp)];
      cout << endl;
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
    cout << "H =" << endl;
    for (i = 0; i < kdim; i++) {
      for (j = 0; j < kdim; j++) 
	cout << setw (14) << H[Veclib::col_major (i, j, kdim)];
      cout << endl;
    }
  }

  // -- Find eigenpairs of H.

  F77name(dgeev) ("N","V",kdim,H,kdim,wr,wi,0,1,zvec,kdim,rwork,lwork,ier);

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
		      const int_t  nvec   )
// ---------------------------------------------------------------------------
// Return value:
//   nvec:  all of the first nvec eigenvalue estimates have converged;
//  -1,-2:  the residuals aren't shrinking;
//   0:     neither of the above is true: not converged.
// ---------------------------------------------------------------------------
{
  int_t          i, idone;
  vector<real_t> work (kdim);
  real_t         re_ev, im_ev, max_resid, *resid = &work[0];
  static real_t  min_max1, min_max2;
 
  if (min_max1 == 0.0) min_max1 = 1000.0;
  if (min_max2 == 0.0) min_max2 = 1000.0;

  // -- Sort subspace eigenvectors by residual.

  for (i = 0; i < kdim; i++) {
    resid[i] = resnorm * fabs (zvec[kdim - 1 + i * kdim])
      / sqrt (Blas::dot (kdim, zvec + i * kdim, 1, zvec + i * kdim, 1));
    if (wi[i] < 0.0) resid [i - 1] = resid[i] = hypot (resid[i - 1],resid [i]);
  }
  EV_sort (zvec, wr, wi, resid, kdim);

  // -- Stopping test.

  if      (resid[nvec - 1] < evtol)
    idone = nvec;
  else if (min_max1 < 0.001 && resid[nvec - 1] > 100.0 * min_max1 ||
	   min_max2 < 0.001 && resnorm         > 100.0 * min_max2 )
    idone = -1;
  else
    idone = 0;

  min_max1 = min (min_max1, resid[nvec - 1]);
  min_max2 = min (min_max2, resnorm);

  // -- Print diagnostic information.

  cout << "-- Iteration = " << itrn << ", H(k+1, k) = " << resnorm << endl;

  for (i = 0; i < kdim; i++) {
    re_ev = wr[i];
    im_ev = wi[i];
    cout << "Eigval(" << i << ") = ("
	 << setw(14) << re_ev << ", " 
	 << setw(14) << im_ev
	 << ")       resid = " << resid[i] << endl;
  }
  
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

  
