//////////////////////////////////////////////////////////////////////////////
// buoy.h: header file for N--S with heat transfer & buoyancy.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


class ScatAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  ScatAnalyser (Domain*, FEML*);
  void analyse (AuxField**);

private:
  ofstream flx_strm;
};
