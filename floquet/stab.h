#ifndef STAB_H
#define STAB_H
//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "sem.h"

class StabAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  StabAnalyser (Domain*, FEML*);
  void analyse (AuxField**);

private:
  vector<HistoryPoint*> base_history; // -- Locations, etc. of history points.
  ofstream              bhs_strm;     // -- File for base history points.
};

void integrate (Domain*, StabAnalyser*);

extern "C" {
  void F77NAME(dgeev)		// -- Lapack eigensystem routine.
    (const char*    N    ,
     const char*    V    ,
     const integer& dim1 ,
     double*        H    ,
     const integer& dim2 ,
     double*        wr   ,
     double*        wi   ,
     double*        f1   ,
     const integer& f2   ,
     double*        Hvec ,
     const integer& dim3 ,
     double*        rwork,
     const integer& lwork,
     integer&       ier  );
}

#endif
