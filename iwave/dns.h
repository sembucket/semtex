//////////////////////////////////////////////////////////////////////////////
// dns.h: header file for direct numerical simulation solver.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>


class DNSAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  DNSAnalyser  (Domain*, FEML*);
  void analyse (AuxField**, AuxField**);

private:
  ofstream flx_strm;
};

void nonLinear (Domain*, AuxField**, AuxField**);
