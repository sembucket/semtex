#ifndef STAB_H
#define STAB_H
//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "Sem.h"

class STABAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  STABAnalyser (Domain*, FEML*);
  void analyse (AuxField**);

private:
  vector<HistoryPoint*> base_history; // -- Locations, etc. of history points.
  ofstream              bhs_strm;     // -- File for base history points.
};

void integrate (Domain*, STABAnalyser*);

#define F77name(x) x ## _

extern "C" {
  void F77name(dgeev)		// -- Lapack eigensystem routine.
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
