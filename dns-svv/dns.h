#ifndef DNS_H
#define DNS_H
//////////////////////////////////////////////////////////////////////////////
// dns.h: header file for direct numerical simulation solver.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include "sem.h"


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

namespace SVV {
  const real_t* coeffs (const int_t, const int_t);
  void operators (const int_t, const real_t**, const real_t**);
}

#endif
