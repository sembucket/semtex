///////////////////////////////////////////////////////////////////////////////
//flok.C: compute leading eigenvalues and eigenvectors for stability
//analysis based on linearised Navier--Stokes operators.
// 
// Based on code by Dwight Barkley & Ron Henderson.
//
// USAGE
// -----
// flok [options] session
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
//
// REFERENCES
// ----------
// [1]  D. Barkley & R.D. Henderson (1996), "Three-dimensional Floquet
//      stability analysis of the wake of a circular cylinder",
//      J. Fluid Mech V322, 215--241.
// [2]  Y. Saad (1991), "Numerical methods for large eigenvalue problems"
//      Wiley.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "stab.h"

static char prog[] = "flok";

static void getargs    (int, char**, int&, int&, int&, int&, real&, char*&);
static int  preprocess (const char*);

static void EV_init    (real*);
static void EV_update  (const real*, real*);
static void EV_small   (real**, const int, const int, 
		        real*, real*, real*, real&, const int); 
static int  EV_test    (const int, const int, real*, real*, real*,
		        const real, const real, const int);
static void EV_sort    (real*, real*, real*, real*, const int);
static int  EV_big     (real**, real**, int, int, real*, real*, real*);

static char*            session;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;
static BoundarySys*     bsys;
static Domain*          domain;
static STABAnalyser*    analyst;


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  int  kdim = 2, nvec = 2, nits = 2, verbose = 0, converged = 0;
  real norm, resnorm, evtol = 1.0e-6;
  int  i, itrn;

  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, kdim, nits, nvec, verbose, evtol, session);

  // -- Check parameter values.

  if (kdim < 1)    message (prog, "param error: KDIM must be > 1",     ERROR);
  if (nvec < 1)    message (prog, "param error: NVEC must be > 1",     ERROR);
  if (nits < kdim) message (prog, "param error: NITS must be >= KDIM", ERROR);
  if (kdim < nvec) message (prog, "param error: NVEC must be <= KDIM", ERROR);
  
  const int ntot = preprocess (session);

  // -- Allocate eigenproblem storage.

  vector<real> work (kdim + kdim + (kdim * kdim) + 3*ntot*(kdim + 1));
  real*        wr   = work();
  real*        wi   = wr   + kdim;
  real*        zvec = wi   + kdim;
  real*        kvec = zvec + kdim * kdim;
  real*        tvec = kvec + ntot * (kdim + 1);
  real*        evec = tvec + ntot * (kdim + 1);
  real**       Kseq = new real* [kdim + 1];
  real**       Tseq = new real* [kdim + 1];
  real**       Eseq = new real* [kdim + 1];

  for (i = 0; i <= kdim; i++) {
    Kseq[i] = kvec + i * ntot;
    Tseq[i] = tvec + i * ntot;
    Eseq[i] = evec + i * ntot;
  }

  // -- Load starting vector.

  EV_init (Kseq[0]);
  norm = Blas::nrm2 (ntot, Kseq[0], 1);
  Blas::scal (ntot, 1.0/norm, Kseq[0], 1);

  // -- Fill initial Krylov sequence.

  for (i = 1; i <= kdim; i++) {
    EV_update (Kseq[i - 1], Kseq[i]);
    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small  (Tseq, ntot, i, zvec, wr, wi, resnorm, verbose);
    EV_test   (i, i, zvec, wr, wi, resnorm, evtol, i);
  }

  // -- Carry out iterative solution.

  for (itrn = kdim; !converged && itrn <= nits; itrn++) {

    if (itrn != kdim) {
      norm = Blas::nrm2 (ntot, Kseq[1], 1);
      for (i = 1; i <= kdim; i++) {
	Blas::scal   (ntot, 1.0/norm, Kseq[i], 1);
	Veclib::copy (ntot, Kseq[i], 1, Kseq[i - 1], 1);
      }
      EV_update (Kseq[kdim - 1], Kseq[kdim]);
    }
    
    // -- Get subspace eigenvalues, test for convergence.

    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, kdim, zvec, wr, wi, resnorm, verbose); 

    converged = EV_test (itrn, kdim, zvec, wr, wi, resnorm, evtol, nvec);
  }
 
  if (!converged)
    message (prog, "not converged", ERROR);
  else if (converged == nvec) {
    message (prog, ": all estimates converged",  REMARK);
  } else
    message (prog, ": minimum residual reached", REMARK);

  return (EXIT_SUCCESS);
}


static void EV_init (real* tgt)
// ---------------------------------------------------------------------------
// Load initial vector from domain velocity fields.
// ---------------------------------------------------------------------------
{
  int       i;
  const int DIM = Geometry::nVec();
  const int NP  = Geometry::planeSize();
    
  for (i = 0; i < DIM; i++) domain -> u[i] -> getPlane (0, tgt + i * NP);
}


static void EV_update  (const real* src,
			real*       tgt)
// ---------------------------------------------------------------------------
// Generate tgt by applying linear operator to src.
// ---------------------------------------------------------------------------
{
  integer       j;
  const integer NDIM = Geometry::nVec();
  const integer NP   = Geometry::planeSize();

  for (j = 0; j < NDIM; j++) domain -> u[j] -> setPlane (0, src + j*NP);

  domain -> step = 0;
  integrate (domain, analyst);

  for (j = 0; j < NDIM; j++) domain -> u[j] -> getPlane (0, tgt + j*NP);
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
  real         re_ev, im_ev, abs_ev, *resid = work();
  real         re_Aev, im_Aev, period;
  static real  min_max1, min_max2;
 
  if (min_max1 == 0.0) min_max1 = 1000.0;
  if (min_max2 == 0.0) min_max2 = 1000.0;

  // -- Sort subspace eigenvectors by residual.

  for (i = 0; i < kdim; i++) {
    resid[i] = resnorm * fabs (zvec[kdim - 1 + i * kdim])
      / sqrt (Blas::dot (kdim, zvec + i * kdim, 1, zvec + i * kdim, 1));
    if (wi[i] < 0.0) resid [i - 1] = resid[i] = hypot (resid[i - 1], resid[i]);
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

  cout.precision(5);
  cout.setf(ios::fixed, ios::floatfield);

  period = Femlib::value ("D_T") * (integer) Femlib::value ("N_STEP");

  cout << "Eigval   Re\t    Im\t       ABS()"
       << "\t  Resid\t     a\t        b" << endl;
  for (i = 0; (i < kdim) & (i < 10); i++) {
    re_ev = wr[i];
    im_ev = wi[i];
    abs_ev = sqrt(re_ev*re_ev + im_ev*im_ev);
    re_Aev = 1/(2*period) * log(re_ev*re_ev + im_ev*im_ev);
    im_Aev = 1/period *atan2(im_ev, re_ev);
    cout <<i << " =\t"
	 << setw(8) << re_ev << ",  " 
	 << setw(8) << im_ev << ",  "
         << setw(8) << abs_ev << ",  "
	 << setw(8) << resid[i] << ",  "
	 << setw(8) << re_Aev << ",  " 
	 << setw(8) << im_Aev << endl;
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


static int EV_big (real**  K_Seq,   // Krylov subspace matrix
	           real**  evecs,   // eigenvectors
		   int     ntot,
		   int     kdim,
		   real*   z_vec,
		   real*   wr,
		   real*   wi)
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  real norm, wgt;
  int i, j;

  /* generate big e-vector  */

  for(j=0;j<kdim;j++){
    Veclib::zero (ntot, evecs[j], 1);

    for(i=0;i<kdim;i++){
      wgt = z_vec[i+j*kdim];
      Blas::axpy (ntot, wgt, K_Seq[i], 1, evecs[j], 1);
    }
  }

  /* normalize big e-vectors */
  for(i=0;i<kdim;i++){
    if(wi[i]==0.0) {
      norm = Blas::nrm2 (ntot, evecs[i], 1);
      Blas::scal (ntot, 1.0/norm, evecs[i], 1);
    }
    else if(wi[i] > 0.0) {
      norm  = pow(Blas::nrm2(ntot, evecs[i], 1),  2.);
      norm += pow(Blas::nrm2(ntot, evecs[i+1], 1),2.);
      norm = sqrt(norm);
      Blas::scal (ntot, 1./norm, evecs[i], 1 );
      Blas::scal (ntot, 1./norm, evecs[i+1], 1);
      i++;
    }

  }

  // this next bit requires use of A_op which hasn't yet been done
  // J.E 3/3/2000
  
#if BIG_RESIDS
  // Compute residuals of big vectors directly
  real resid;
 
  for(i=0;i<nwrt;i++){
    V_copy(ntot, evecs[i], tvecs[0]);
    A_op(tvecs[0]);
    tsteps -= nsteps;
    ftime  -= nsteps*dt;
    if(wi[i]==0.) {
      V_axpy(ntot, -wr[i], evecs[i], tvecs[0]);
      resid = V_nrm2(ntot, tvecs[0]);
      printf("big resid(%d) = %g \n", i, resid);
    }
    else if(wi[i] > 0.) {
      V_axpy(ntot, -wr[i], evecs[i],   tvecs[0]);
      V_axpy(ntot,  wi[i], evecs[i+1], tvecs[0]);
      resid = V_nrm2(ntot, tvecs[0])/V_nrm2(ntot, evecs[i]);
      printf("big resid(%d) = %g \n", i, resid);
    }
    else {
      V_axpy(ntot, -wr[i], evecs[i],   tvecs[0]);
      V_axpy(ntot,  wi[i], evecs[i-1], tvecs[0]);
      resid = V_nrm2(ntot, tvecs[0])/V_nrm2(ntot, evecs[i]);
      printf("big resid(%d) = %g \n", i, resid);
    }
  }
#endif

  return 0;

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
  char usage[] = "flok [options] session\n"
    "session: specifies name of session file\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-k <num> ... set dimension of subspace (maximum number of pairs) to num\n"
    "-m <num> ... set maximum number of iterations         (m >= k)\n"
    "-n <num> ... compute num eigenvalue/eigenvector pairs (n <= k)\n"
    "-t <num> ... set eigenvalue tolerance to num [Default 1e-6]\n"
    "-chk     ... checkpoint field dumps\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      verbose = 1;
      Femlib::value ("VERBOSE",   (integer) Femlib::value ("VERBOSE")   + 1);
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
// Create objects needed for semtex execution, given the session file
// name.
//
// Return length of an eigenvector.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys space;
  const real*        z;
  integer            i, np, nz, nel, nvec;

  file  = new FEML (session);
  mesh  = new Mesh (file);

  nel   = mesh -> nEl();
  np    =  (integer) Femlib::value ("N_POLY");
  nz    =  (integer) Femlib::value ("N_Z");
  space = ((integer) Femlib::value ("CYLINDRICAL")) ? 
                     Geometry::Cylindrical : Geometry::Cartesian;
  nvec  = file -> attribute ("FIELDS", "NUMBER") - 1;
  Geometry::set (np, nz, nel, space, nvec);

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

  bman    = new BCmgr        (file, elmt);
  domain  = new Domain       (file, elmt, bman);
  analyst = new STABAnalyser (domain, file);

  domain -> restart ();
  domain -> loadbase();
  domain -> report  ();

  return Geometry::nVec() * Geometry::planeSize();
}
