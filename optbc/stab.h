#ifndef STAB_H
#define STAB_H

//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>
#include <data2df.h>

typedef enum { PRIMAL, ADJOINT, GROWTH, SHRINK, EVD, SVD } problem_t;

// -- EVD => eigenvalue decomposition.
// -- SVD => singular value decomposition.

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
		 Domain*, StabAnalyser*, vector<real_t*>& ,vector<real_t*>&);
void integrate  (void (*)(Domain*, AuxField**, AuxField**,  vector<AuxField*>,  vector<AuxField*>, int_t),
				 Domain*, StabAnalyser*, vector<real_t*>& ,vector<real_t*>&, vector<real_t*>&, real_t*, real_t*, int_t);
void linAdvect  (Domain*, AuxField**, AuxField**);
void linAdvectT (Domain*, AuxField**, AuxField**);
void linAdvect_unsteady  (Domain*, AuxField**, AuxField**, vector<AuxField*>, vector<AuxField*>, int_t);
void linAdvectT_unsteady (Domain*, AuxField**, AuxField**, vector<AuxField*>, vector<AuxField*>, int_t);
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