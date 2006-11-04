//////////////////////////////////////////////////////////////////////////////
// scat.h: header file for N--S with heat transfer & buoyancy.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>


class ScatAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  ScatAnalyser (Domain*, FEML*);
  void analyse (AuxField**, AuxField**);

private:
  ofstream flx_strm;
};
