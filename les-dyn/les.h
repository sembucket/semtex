//////////////////////////////////////////////////////////////////////////////
// les.h: header file for LES solver.
//
// Copyright (c) Hugh Blackburn 1998--1999
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
  void analyse (AuxField***);

private:
  ofstream flx_strm;
};
