#ifndef LES_SMAG_H
#define LES_SMAG_H
//////////////////////////////////////////////////////////////////////////////
// les.h: header file for LES solver.
//
// Copyright (c) 1998 <--> $Date$, Hugh Blackburn
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "sem.h"


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


// -- In eddyvis.C:

void eddyViscosity (const Domain*, AuxField**, AuxField**, AuxField*);
#endif
