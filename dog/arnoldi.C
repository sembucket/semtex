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


#include <stab.h>
#include <unistd.h>

#define F77name(x) x ## _

extern "C" {
  void F77name(dgeev)		// -- Lapack eigensystem routine.
    (const char*    N    ,      // computes all eigenvalues of matrix ??
     const char*    V    ,      // and returns as ??
     const integer& dim1 ,
     double*        H    ,
     const integer& dim2 ,
     double*        wr   ,
     double*        wi   ,
     double* f1   ,
     const integer& f2   ,
     double*        Hvec ,
     const integer& dim3 ,
     double*          rwork,
     const integer& lwork,
     integer&       ier  );
}

static ofstream   eig_strm; // File for eigenvalues.

static char prog[] = "arnoldi";
static void getargs  (int, char**, int&, int&, int&, int&, real&,
		      int&, char*&);
static void EV_small (real**, const int, const int, 
		      real*, real*, real*, real&, const int); 
static int  EV_test  (const int, const int, real*, real*, real*,
		      const real, const real, const int, const int,
		      int, const real);
static void EV_sort  (real*, real*, real*, real*, const int);
static int  EV_big (real**, real**, int, int, real*, real*, real*);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, BoundarySys*&, Domain*&);
void NavierStokes (Domain*, STABAnalyser*);
void check_input(int, int, int);
void arnoldi_report(int, int, int, int);
void write_restart (char*, int, int, int, int, int, real, real, real*);
void create_eigstrm(ofstream&, char*);

// ------------------ MAGIC NUMBERS (default values) ------------------------

const int    kdim_def = 2;       
const int    nvec_def = 2;
const int    nits_def = 2;
const real  evtol_def = 1.0e-12;

#if BIG_RESIDS
    const int        nwrt = 3;       
#endif

// --------------------------------------------------------------------------
// ------------------- MAIN -------------------------------------------------
// --------------------------------------------------------------------------

int main (int    argc,
	  char** argv)
{
  int          kdim = kdim_def;
  int          nvec = nvec_def;
  int          nits = nits_def;
  real        evtol = evtol_def;

  int              verbose = 0;
  int              restart = 0;      // FALSE - no restart file
  char             res_name[StrMax];
  char             stp_cmd[StrMax];
  int              stop_cnd = 0;     // default = false
  ifstream         rst_strm;         // Krylov matrix restart file.
  char*            session;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  BoundarySys*     bsys;
  Domain*          domain;
  STABAnalyser*     adjunct;
  integer       i, it_start, itrn, j, converged = 0;  
  real          norm, resnorm;
  int           Total_step = 0;  // default value is overiden by restart


  Femlib::initialize (&argc, &argv);

  Femlib::value ("DUMP"      , 1 );  // -- toggle field dump for each NS call
  Femlib::value ("IO_STEP"   , 1 );  // -- 0 = no step output
  Femlib::value ("IO_EIG"    , 2 );  // -- # of eigenvalues to output
  Femlib::value ("MAX_BASE_F", 50);  // -- max # of base files
  Femlib::value ("KDIM"      , 20);  // -- dimension of krylov matrix
  Femlib::value ("NVEC"      , 5 );  // -- number of eigenvectors to test
  Femlib::value ("IO_KRY"    , 10);  // -- how often to write a krylov file

  getargs (argc, argv, kdim, nits, nvec, verbose, evtol, restart, session);
  
  if (restart) {
    // check for file and open if exists.
    rst_strm.open(strcat (strcpy(res_name,session), ".kry"));

    if (rst_strm) {
      cout << "Reading Krylov matrix data from " << res_name << endl;
      
      rst_strm >> kdim;  // cannot be overriden -> Kseq dimensions.

      // get nvec and nits if command line values not equal to default
      int trash;  // holds restart integers - if not needed.
      rst_strm >> (nits > nits_def ? trash : nits); 
      rst_strm >> (nvec > nvec_def ? trash : nvec);
      rst_strm >> Total_step;
      rst_strm >> evtol;
      
    }
    else message ("arnoldi", "no restart file found", ERROR);
  }
  else {
    check_input(kdim, nvec, nits);
  }

  preprocess (session, file, mesh, elmt, bman, bsys, domain);

  adjunct = new STABAnalyser (domain, file);

  if (!restart) domain -> restart();

  domain -> loadbase();

  if (restart) rst_strm >> domain->time;

  domain -> report();

  create_eigstrm(eig_strm, domain-> name);  
    		 
  int  DIM =  domain -> nField() -1;
  
  cout << "DIM ....." << DIM << endl;
  int ntot = (domain -> nField() -1) * Geometry::planeSize();

  arnoldi_report(ntot, kdim, nits, DIM);
  
  // -- Allocate eigenproblem storage.

  vector<real>  work (kdim + kdim + (kdim * kdim) +
		      3*ntot*(kdim + 1) );

  real*         wr   = work();
  real*         wi   = wr   + kdim;
  real*         zvec = wi   + kdim;
  real*         kvec = zvec + kdim * kdim;
  real*         tvec = kvec + ntot * (kdim + 1);
  real*         evec = tvec + ntot * (kdim + 1);
  real**        Kseq = new real* [kdim + 1];
  real**        Tseq = new real* [kdim + 1];
  real**        Eseq = new real* [kdim + 1];

  for (i = 0; i <= kdim; i++) {
    Kseq[i] = kvec + i * ntot;
    Tseq[i] = tvec + i * ntot;
    Eseq[i] = evec + i * ntot;
  }


  // ****************************************************************
  // ************** construct or load krylov matrix *****************
  // ****************************************************************

  if (!restart) { // ************ construction ************

    // -- Generate random initial guess, this random is 0 - 1 !?
    Veclib::vrandom (ntot, Kseq[0], 1); 
    // shift to range (-0.5 to 0.5)
    Veclib::sadd (ntot, -0.5, Kseq[0], 1, Kseq[0], 1); 

    norm = Blas::nrm2 (ntot, Kseq[0], 1);
    Blas::scal (ntot, 1.0/norm, Kseq[0], 1);

    // -- Fill initial Krylov sequence.
    for (i = 1; i <= kdim; i++) {

      // copy Kseq[i-1] across to domain structure 
      for (j = 0; j < DIM; j++) {
	domain -> u[j]->setPlane(0, j * Geometry::planeSize() + Kseq[i-1]);
      }

      // call to Linear NS routine to fill Kseq.
      domain -> step = 0;
      NavierStokes (domain, adjunct);
      Total_step += domain->step;

      // copy out fields from NS to Kseq[i]
      for(j = 0; j < DIM; j++)
	domain -> u[j]->getPlane(0, j* Geometry::planeSize() + Kseq[i]);
 
      Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);

      EV_small (Tseq, ntot, i, zvec, wr, wi, resnorm, verbose);
      EV_test  (i, i, zvec, wr, wi, resnorm, evtol, i, Total_step,
		restart, domain->time);
    }
    // write restart file now and at end of iterations.
    cout << "matrix constructed ... writing restart";

    write_restart(session, kdim, nits, nvec, ntot, Total_step,
		  evtol, domain->time, kvec);
  }
  else { // *************** restarting from krylov matrix *********

    cout << "loading krylov matrix.....";

    // binary read of date -> make sure endl is removed.
    while (rst_strm.get() != '\n') continue;
    rst_strm.read((char *) kvec, (ntot*(kdim+1)*sizeof(real)));

    // check input is sensible (enabled if debugging)
    // for (int q = 0 ;q <10 ;q++) cout << setw(20) << kvec[q] << endl;

    rst_strm.close();
    cout << "done" << endl << endl;

    // compare with last dump from previous run to validate restart.

    cout << "Restart file eigenvalue status:" << endl; 
    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, kdim, zvec, wr, wi, resnorm, verbose); 
    converged = EV_test (0, kdim, zvec, wr, wi, resnorm, evtol, nvec,
			 Total_step, restart, domain->time);
  }

  // construct command for stop file test.
  strcat(strcat( strcpy(stp_cmd, "test -e "), session),".stp");

  // ********************************************************************
  // **************** Iterative Loop ************************************
  // ********************************************************************

  // -- Carry out iterative solution.

  it_start = (restart ? 0 : kdim);

  for (itrn = ++it_start; !converged && !stop_cnd && itrn <= nits; itrn++) {

    if ((itrn > kdim) && (itrn & 10))
      write_restart(session, kdim, nits, nvec, ntot, Total_step,
		    evtol, domain->time, kvec);

    stop_cnd = !system(stp_cmd);

    if (itrn != kdim || restart) {	     // -- Roll vectors, normalise.
      norm = Blas::nrm2 (ntot, Kseq[1], 1);
      for (i = 1; i <= kdim; i++) {
	Blas::scal   (ntot, 1.0/norm, Kseq[i], 1);
	Veclib::copy (ntot, Kseq[i], 1, Kseq[i - 1], 1);
      }
      
      for(j = 0; j < DIM; j++)
	domain -> u[j]->setPlane(0, j * Geometry::planeSize() + Kseq[kdim-1]);

      // setup and call Linear NS op (pass vector address)
      // domain -> time = 0.0;
      domain -> step = 0;
      NavierStokes (domain, adjunct);
      Total_step += domain->step;

      for(j = 0; j < DIM; j++)
	domain -> u[j]->getPlane(0, j* Geometry::planeSize() + Kseq[kdim]);
    }
    
    // -- Get subspace eigenvalues, test for convergence.

    Veclib::copy (ntot * (kdim + 1), kvec, 1, tvec, 1);
    EV_small (Tseq, ntot, kdim, zvec, wr, wi, resnorm, verbose); 

    converged = EV_test (itrn, kdim, zvec, wr, wi, resnorm, evtol, nvec,
			 Total_step, restart, domain->time);
  }

  if (stop_cnd) cout << "encountered stop file... stopping program" << endl;

  // ******************************************************************
  // ************** output results to screen and file *****************
  // ******************************************************************


  if (itrn > kdim) write_restart(session, kdim, nits, nvec, ntot,
				 Total_step, evtol, domain-> time, kvec);

  // - calculate matching eigenvectors.
  cout << "calculating eigenvectors with length "<< ntot << endl;
  EV_big (Tseq, Eseq, ntot, kdim, zvec, wr, wi);

  cout << "dumping leading eigenvector to field file";
  domain -> dump ();
  cout << endl;
 
  if      (!converged)
    message (prog, "not converged", ERROR);
  else if (converged == nvec) {
    message (prog, ": all estimates converged",  REMARK);
    
  }
  else
    message (prog, ": minimum residual reached", REMARK);

  // dump leading eigenvector to field file


  return (EXIT_SUCCESS);
}

// --------------------------------------------------------------------------
// --------------- EV_small -------------------------------------------------
// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------
// --------------- EV_test --------------------------------------------------
// --------------------------------------------------------------------------

static int EV_test (const int  itrn   ,
		    const int  kdim   ,
		    real*      zvec   ,
		    real*      wr     ,
		    real*      wi     ,
		    const real resnorm,
		    const real evtol  ,
		    const int  nvec,
		    const int  step,
		    int restart,
		    const real time)
// ---------------------------------------------------------------------------
// Return true if converged.
// ---------------------------------------------------------------------------
{
  int          i, idone;
  vector<real> work (kdim);
  real         re_ev, im_ev, abs_ev, *resid = work();
  real         re_Aev, im_Aev, period;
  static real  min_max1, min_max2;
  int   nEig = (integer) Femlib::value ("IO_EIG");
 
 
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

  cout.precision(5);
  // cout.setf(ios::showpoint);
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

    if ((i < nEig) && (itrn > nEig || restart)) {
      eig_strm << setprecision(2) << setw(10) << time
	       << setw(8) << step
	       << setw(7) << i
	       << setprecision(6) << setw(15) << re_ev 
	       << setw(15) << im_ev
	       << setw(15) << abs_ev
	       << setw(15) << resid[i]
	       << setw(15) << re_Aev
	       << setw(15) << im_Aev << endl;

    }
  }
  
  // return 0;
  return idone;
}

// --------------------------------------------------------------------------
// --------------- EV_sort --------------------------------------------------
// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------
// --------------- EV_big ---------------------------------------------------
// --------------------------------------------------------------------------

static int EV_big (real**  K_Seq,   // Krylov subspace matrix
	           real**  evecs,   // eigenvectors
		   int     ntot,
		   int     kdim,
		   real*   z_vec,
		   real*   wr,
		   real*   wi)

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

// --------------------------------------------------------------------------
// --------------- GETARGS --------------------------------------------------
// --------------------------------------------------------------------------

static void getargs (int    argc,
		     char** argv ,
		     int&   kdim , 
		     int&   maxit,
		     int&   neval,
		     int&   verb ,
		     real&  evtol,
		     int&   restart,
		     char*& session)

  /* Parse command-line arguments. */

{
  char usage[] = "arnoldi [options] session\n"
    "session: specifies name of session file\n"
    "options:\n"
    "-h       ... print this message\n"
    "-v       ... set verbose\n"
    "-k <num> ... set dimension of subspace (maximum number of pairs) to num\n"
    "-m <num> ... set maximum number of iterations         (m >= k)\n"
    "-n <num> ... compute num eigenvalue/eigenvector pairs (n <= k)\n"
    "-t <num> ... set eigenvalue tolerance to num [Default 1e-6]\n"
    "-chk     ... checkpoint field dumps\n"
    "-r       ... restart from session.kry file. must supply session.\n"
    "               uses old values of n and m unless overiden.\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      // verb = 1;
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
    case 'r':
      restart = 1;    // TRUE
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc != 1) message (prog, "no session file",   ERROR);
  else             session = *argv;
}

// --------------------------------------------------------------------------
// --------------- PREPROCESS -----------------------------------------------
// --------------------------------------------------------------------------

static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			BoundarySys*&     bsys   ,
			Domain*&          domain )

  /* Create objects needed for execution, given the session file name.
     They are listed in order of creation. */

{
  const integer      verbose = (integer) Femlib::value ("VERBOSE");
  Geometry::CoordSys space;
  const real*        z;
  integer            i, np, nz, nel;

  // -- Initialise problem and set up mesh geometry.

  VERBOSE cout << "Building mesh ..." << endl;

  file = new FEML (session);
  mesh = new Mesh (file);

  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.

  VERBOSE cout << "Setting geometry ... ";

  nel   = mesh -> nEl();
  np    =  (integer) Femlib::value ("N_POLY");
  nz    =  (integer) Femlib::value ("N_Z");
  space = ((integer) Femlib::value ("CYLINDRICAL")) ? 
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);

  VERBOSE cout << "done" << endl;

  // -- Build all the elements.

  VERBOSE cout << "Building elements ... ";

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, mesh, z, np);

  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.

  VERBOSE cout << "Building boundary condition manager ..." << endl;

  bman = new BCmgr (file, elmt);

  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.

  VERBOSE cout << "Building domain ..." << endl;

  domain = new Domain (file, elmt, bman);

  VERBOSE cout << "done" << endl;
}

// --------------------------------------------------------------------------
// ---------- check_input ---------------------------------------------------
// --------------------------------------------------------------------------

void check_input(int kdim, int nvec, int nits)

  /* check arguments and exit if incorrect */

{
  if (kdim < 1)    message (prog, "param error: KDIM must be > 1",     ERROR);
  if (nvec < 1)    message (prog, "param error: NVEC must be > 1",     ERROR);
  if (nits < kdim) message (prog, "param error: NITS must be >= KDIM", ERROR);
  if (kdim < nvec) message (prog, "param error: NVEC must be <= KDIM", ERROR);
}

// --------------------------------------------------------------------------
// ---------- arnoldi_report ------------------------------------------------
// --------------------------------------------------------------------------

void arnoldi_report (int ntot, int kdim, int nits, int DIM)

  /* display info on arnoldi matrices */

{
  cout << endl;
  cout << "-------- Arnoldi Info" << endl;
  cout << "plane size          : " << Geometry::planeSize() << endl;
  cout << "physical dimensions : " << DIM << endl;
  cout << "ntot                : " << ntot << endl;
  cout << "kdim                : " << kdim << endl;
  cout << "nits                : " << nits << endl;
  cout << endl;

}

void write_restart (char*   session,
		    int         kdim,
		    int         nits,
		    int         nvec,
		    int         ntot,
		    int         NS_step,
		    real        evtol,
		    real        NS_time,
		    real*       Kvec)
{
  char  kry_name[StrMax];

  ofstream kry_strm;

  kry_strm.open (strcat(strcpy(kry_name,session), ".kry"));
  if (!kry_strm) message ("arnoldi->write_restart",
			"cannot open session.kry for output",
			ERROR);

  // ascii write of required variables
  kry_strm << kdim << endl 
	   << nits << endl 
	   << nvec << endl
	   << NS_step << endl
	   << evtol << endl
	   << NS_time << endl;

  // binary write of memory areas.
  kry_strm.write((char *) Kvec, ntot*(kdim+1)*sizeof(real));
 
  // uncomment this + statement in read area to check binary read/write
  // for (int pl = 0; pl < 15 ; pl++) cout << Kvec[pl] << endl;
  
  kry_strm.close();

  cout << "....done" << endl;
}

void create_eigstrm(ofstream& strm, char* name)
{
  char   str[StrMax], sys_cmd[StrMax];
  int new_file;
  
  strcpy(str, strcat( strcpy(str,name), ".eig")); 
  strcpy(sys_cmd, strcat(strcpy(sys_cmd,"test -e "),str) );
  new_file = system(sys_cmd);
  
  // open eigenvalue file if required (IO_EIG > 0 )
  if ((integer) Femlib::value ("IO_EIG")) {
    
    strm.open (str, (new_file ? ios::out : ios::app));
    if (!strm) cerr << "can't open eigenvalue file" << endl;

    strm.setf (ios::scientific, ios::floatfield);
    if (new_file) {
      strm << '#';
      strm << setw(9) << "time"
	   << setw(8) << "step"
	   << setw(7) << "eig #"
	   << setw(15) << "real e val" 
	   << setw(15) << "imag e val"
	   << setw(15) << "abs e val"
	   << setw(15) << "resid"
	   << setw(15) << "Re A(eval)"
	   << setw(15) << "Im A(eval)" << endl;
      strm << '#';
    
      for (int q = 0; q < (9+8+7+6*15);q++) strm << "-";
      
      strm << endl;
    }

  }
}
