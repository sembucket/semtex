//////////////////////////////////////////////////////////////////////////////
// buoy.h: header file for N--S with heat transfer & buoyancy.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


class BuoyAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  BuoyAnalyser (Domain&, FEML&);
  void analyse (AuxField***);

private:
  ofstream flx_strm;
};
