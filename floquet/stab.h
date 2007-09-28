#ifndef STAB_H
#define STAB_H

//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>

typedef enum { PRIMAL, ADJOINT, GROWTH, SHRINK, EVD, SVD } problem_t;

// -- EVD => eigenvalue decomposition.
// -- SVD => singular value decomposition.

class Krylov
// ===========================================================================
// Data structure to hold storage for Krylov sequences etc.
// ===========================================================================
{
public:
  Krylov (int_t kdim, int_t ndim, problem_t job);

  int_t     K;		// -- Krylov dimension.
  int_t     N;		// -- Full space dimension.
  problem_t task;	// -- Problem under solution.
  
  real_t*   wr;		// -- K subspace eigenvalues (real part).
  real_t*   wi;		// -- K subspace eigenvalues (imag part).
  real_t*   svec;	// -- K x K subspace eigenvectors.
  real_t*   qvec;	// -- M x K+1 orthonormal basis for Krylov sequence.
  real_t*   rvec;       // -- Upper-triangular R (K+1 x K+1): Kseq = Q R
  real_t*   kvec;       // -- N x K+1 Krylov sequence.
  real_t*   jvec;       // -- N x K+1 transpose's Krylov sequence (SVD only).
  real_t**  Kseq;       // -- Handles for individual kvecs.
  real_t**  Jseq;       // -- Handles for individual jvecs.
  real_t**  Q;          // -- Handles for individual qvecs.
};


class StabAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  StabAnalyser (Domain*, FEML*);
  void analyse (AuxField**, AuxField** = 0);

private:
  vector<HistoryPoint*> base_history; // -- Locations, etc. of history points.
  ofstream              bhs_strm;     // -- File for base history points.
};

void integrate  (void (*)(Domain*, AuxField**, AuxField**),
		 Domain*, StabAnalyser*);
void linAdvect  (Domain*, AuxField**, AuxField**);
void linAdvectT (Domain*, AuxField**, AuxField**);

// -- Fortran interfaces to things not included in semtex.

extern "C" {
#if defined (ARPACK)
  void F77NAME(dnaupd) 		// -- ARPACK reverse-communications interface.
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
     int_t&         info  );
  void F77NAME(dneupd)		// -- Postprocessing.
    (const int_t&   rvec  ,
     const char*    howmny,
     const int*     select,
     real_t*        dr    ,
     real_t*        di    ,
     real_t*        z     ,
     const int_t&   ldz   ,
     const real_t&  sigmar,
     const real_t&  sigmai,
     real_t*        workev,
     const char*    bmat  ,	// -- Remainder unchanged after dnaupd.
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
     int_t&         info  );     
#endif
}


#endif
