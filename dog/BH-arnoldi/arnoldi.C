///////////////////////////////////////////////////////////////////////////////
// arnoldi.C: compute leading eigenvalues and eigenvectors of (real)
// arbitrary sparse matrix using Arnoldi iteration.  Sparse matrices
// are read and stored internally using Harwell--Boeing (HB) format.
// 
// Based on code by Dwight Barkley & Ron Henderson.
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
#include <strstream>
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
    (const integer& trans ,
     const integer& n     ,
     const integer& m     ,
     const double*  a     ,
     const integer* colptr,
     const integer* rowind,
     const double*  x     ,
     const integer& ldx   ,
     double*        y     ,
     const integer& ldy   );
}

static char prog[] = "arnoldi";

static void getargs  (int, char**, int&, int&, int&, int&, real&, char*&);
static void EV_small (real**, const int, const int, 
		      real*, real*, real*, real&, const int); 
static int  EV_test  (const int, const int, real*, real*, real*,
		      const real, const real, const int);
static void EV_sort  (real*, real*, real*, real*, const int);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char   *HBfile, *HBtype;
  int    *HBrptr,  *HBcptr;
  double *HBval;
  int     HBnr, HBnc, HBnz, HBnrhs;
  int     ntot, kdim = 2, nvec = 2, nits = 2, verbose = 0;
  real    evtol = 1.0e-6;

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

  // -- Allocate eigenproblem storage.
  
  integer       i, itrn, converged = 0;
  real          norm, resnorm;
  vector<real>  work (kdim + kdim + kdim * kdim + 2 * ntot * (kdim + 1));
  real*         wr   = &work[0];
  real*         wi   = wr   + kdim;
  real*         zvec = wi   + kdim;
  real*         kvec = zvec + kdim * kdim;
  real*         tvec = kvec + ntot * (kdim + 1);
  real**        Kseq = new real* [kdim + 1];
  real**        Tseq = new real* [kdim + 1];

  for (i = 0; i <= kdim; i++) {
    Kseq[i] = kvec + i * ntot;
    Tseq[i] = tvec + i * ntot;
  }

  // -- Generate random initial guess, zero evals.

  Veclib::vrandom (ntot, Kseq[0], 1);
  norm = Blas::nrm2 (ntot, Kseq[0], 1);
  Blas::scal (ntot, 1.0/norm, Kseq[0], 1);

  // -- Fill initial Krylov sequence.

  for (i = 1; i <= kdim; i++) {
    F77NAME(dspmvc) (0, HBnr, 1, HBval, HBcptr, HBrptr,
		     Kseq[i - 1], HBnr, Kseq[i], HBnr);
    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, i, zvec, wr, wi, resnorm, verbose);
    EV_test  (i, i, zvec, wr, wi, resnorm, evtol, i);
  }

  // -- Carry out iterative solution.

  for (itrn = kdim; !converged && itrn < nits; itrn++) {

    if (itrn != kdim) {		// -- Roll vectors, normalise.
      norm = Blas::nrm2 (ntot, Kseq[1], 1);
      for (i = 1; i <= kdim; i++) {
	Blas::scal   (ntot, 1.0/norm, Kseq[i], 1);
	Veclib::copy (ntot, Kseq[i], 1, Kseq[i - 1], 1);
      }

      // -- Matrix-vector product.

      F77name(dspmvc) (0, HBnr, 1, HBval, HBcptr, HBrptr,
		       Kseq[kdim - 1], HBnr, Kseq[kdim], HBnr);
    }

    // -- Get subspace eigenvalues, test for convergence.

    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, kdim, zvec, wr, wi, resnorm, verbose); 

    converged = EV_test (itrn, kdim, zvec, wr, wi, resnorm, evtol, nvec);
  }

  if      (!converged)
    message (prog, "not converged", ERROR);
  else if (converged == nvec)
    message (prog, ": all estimates converged",  REMARK);
  else
    message (prog, ": minimum residual reached", REMARK);

  return (EXIT_SUCCESS);
}


static void getargs (int    argc,
		     char** argv ,
		     int&   kdim , 
		     int&   maxit,
		     int&   neval,
		     int&   verb ,
		     real&  evtol,
		     char*& mfile)
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
// destroyed in the process.
//
// Then we compute Hessenberg matrix H = Q* A Q (using Q* A = R), find
// its eigenvalues (wr, wi) and eigenvectors, zvec, related to those of A.
//
// The residual norm for each eigenvector is related to H(kdim+1, kdim),
// which is passed back for convergence testing.
// ---------------------------------------------------------------------------
{
  char         routine[] = "EV_small";
  const int    kdimp = kdim + 1;
  int          i, j, ier, lwork = 10 * kdim;
  vector<real> work (kdimp * kdimp + kdim * kdim + lwork);
  real         *R = &work[0], *H = R + kdimp * kdimp, *rwork = H + kdim * kdim;

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
// Return true if converged.
// ---------------------------------------------------------------------------
{
  int          i, idone;
  vector<real> work (kdim);
  real         re_ev, im_ev, max_resid, *resid = &work[0];
  static real  min_max1, min_max2;
 
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
  real         wr_tmp, wi_tmp, te_tmp, *z_tmp = &work[0];

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

  
