#ifndef LES_H
#define LES_H
//////////////////////////////////////////////////////////////////////////////
// les.h: header file for LES solver.
//
// Copyright (c) Hugh Blackburn 1998
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
  void analyse (AuxField**);

private:
  ofstream _flx_strm;
};


class SumIntegrator
// ===========================================================================
// Implement first-order system smoothing of dynamic Smag estimate.
// ===========================================================================
{
friend ifstream& operator >> (ifstream&, SumIntegrator&);
friend ofstream& operator << (ofstream&, SumIntegrator&);
public:
  SumIntegrator  (Domain*);
  ~SumIntegrator () { delete [] _work; }

  void update (real*);
  void dump   ();

private:
  const Domain* _domain;
  AuxField*     _Lmix2 ;
  real*         _work  ;
  real          _BB    ;
  real          _AA    ;
  integer       _ntot  ;
  integer       _nz    ;
};


// -- filter.C:

void initFilters ();
void lowpass     (real*);

// -- integrate.C:

void integrate (Domain*, LESAnalyser*, SumIntegrator*);

// -- nonlinear.C:

void nonLinear (Domain*, SumIntegrator*, vector<real*>&, vector<real>&);
void dynamic   (Domain*, vector<real*>&, const int = 1);
#endif
