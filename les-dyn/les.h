//////////////////////////////////////////////////////////////////////////////
// les.h: header file for LES solver.
//
// Copyright (c) Hugh Blackburn 1998--2000
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


class LESAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  LESAnalyser  (Domain*, FEML*);
  void analyse (AuxField**);

private:
  ofstream flx_strm;
};

// -- filter.C:

void initFilters ();
void lowpass     (real*);

// -- integrate.C:

void integrate (Domain*, LESAnalyser*);

// -- nonlinear.C:

void nonLinear (Domain*, matrix<real>&, vector<real>&);
