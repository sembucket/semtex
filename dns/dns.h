#ifndef DNS_H
#define DNS_H
//////////////////////////////////////////////////////////////////////////////
// dns.h: header file for direct numerical simulation solver.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


class DNSAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  DNSAnalyser  (Domain*, FEML*);
  void analyse (AuxField***);

private:
  ofstream flx_strm;
};
#endif
