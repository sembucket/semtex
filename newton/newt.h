#ifndef NEWT_H
#define NEWT_H
//////////////////////////////////////////////////////////////////////////////
// newt.h: header file for Newton's method steady state NS solver.
//
// http://www.netlib.org/templates/index.html
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "Sem.h"

class NewtAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  NewtAnalyser (Domain*, FEML*);
  void analyse (AuxField**);

private:
  vector<HistoryPoint*> base_history; // -- Locations, etc. of history points.
  ofstream              bhs_strm;     // -- File for base history points.
};


typedef void (*Advection) (Domain*, AuxField**, AuxField**);
void integrate (Domain*, Analyser*, Advection);
//Advection linear, nonlinear;
void nonlinear (Domain*, AuxField**, AuxField**);
void linear    (Domain*, AuxField**, AuxField**);


extern "C" {
  void F77NAME(bicgstab)	// -- Templates Bi-Conj-Grad-Stab solver.
    (const integer& N    ,
     const real*    B    ,
     real*          X    ,
     real*          WORK ,
     const integer& LDW  ,
     integer&       ITER ,
     real&          RESID,
     void (*MATVEC) (const real&, const real*, const real&, real*),
     void (*PSOLVE) (real*, const real*),
     integer&       INFO );
}

#endif
