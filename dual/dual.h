#ifndef DUAL_H
#define DUAL_H
//////////////////////////////////////////////////////////////////////////////
// dual.h: header file for direct numerical simulation solver.
//
// Copyright (c) 1999 Hugh Blackburn
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "Sem.h"


class DualAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  DualAnalyser  (Domain*, FEML*);
  void analyse (AuxField**);

private:
  ofstream flx_strm;
};
#endif
