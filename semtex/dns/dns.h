#ifndef DNS_H
#define DNS_H
//////////////////////////////////////////////////////////////////////////////
// dns_h: header file for direct numerical simulation solver.
//
// Copyright (c) 1994,2003 Hugh Blackburn
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
  void analyse (AuxField**);

private:
  ofstream flx_strm;
};

void nonLinear (Domain*, AuxField**, AuxField**, vector<real>&);
#endif
