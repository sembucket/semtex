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

#endif
