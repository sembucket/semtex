//////////////////////////////////////////////////////////////////////////////
// nnewt.h: header file for NNEWT solver.
//
// Copyright (c) 1998 <--> $Date$,
//   Murray Rudman, Hugh Blackburn
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "sem.h"


class nnewtAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  nnewtAnalyser (Domain*, FEML*);
  void analyse  (AuxField**, AuxField**, AuxField*);

private:
  ofstream flx_strm;
};


// -- In nnewtvis.C:

void viscosity (const Domain*, AuxField**, AuxField**, AuxField*);
