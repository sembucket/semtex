#ifndef STAB_H
#define STAB_H
//////////////////////////////////////////////////////////////////////////////
// stab.h: header file for linearised NS stability solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

class STABAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  STABAnalyser  (Domain*, FEML*);
  void analyse (AuxField**);

private:
  vector<HistoryPoint*> base_history  ; // Locations, etc. of history points.
  ofstream              bhs_strm ; // file for base history points.
};
#endif
