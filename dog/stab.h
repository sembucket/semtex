#ifndef STAB_H
#define STAB_H

//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>

typedef enum { PRIMAL, ADJOINT, GROWTH } problem_t;

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

void integrate (const problem_t, Domain*, StabAnalyser*);

extern "C" {
  void F77NAME(dgeev)		// -- Lapack eigensystem routine.
    (const char*    N    ,
     const char*    V    ,
     const int_t&   dim1 ,
     double*        H    ,
     const int_t&   dim2 ,
     double*        wr   ,
     double*        wi   ,
     double*        f1   ,
     const int_t&   f2   ,
     double*        Hvec ,
     const int_t&   dim3 ,
     double*        rwork,
     const int_t&   lwork,
     int_t&         ier  );
}

#endif
