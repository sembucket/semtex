#ifndef LES_H
#define LES_H
//////////////////////////////////////////////////////////////////////////////
// les.h: header file for LES solver.
//
// Copyright (c) Hugh Blackburn 1998
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>


class LESAnalyser : public Analyser
// ===========================================================================
// Implement step-by-step processing and output control for flow solver.
// ===========================================================================
{
public:
  LESAnalyser  (Domain*, FEML*);
  void analyse (AuxField**, AuxField**);

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

  void update (real_t*);
  void dump   ();

private:
  const Domain* _domain;
  AuxField*     _Lmix2 ;
  real_t*       _work  ;
  real_t        _BB    ;
  real_t        _AA    ;
  int_t         _ntot  ;
  int_t         _nz    ;
};


// -- filter.C:

void initFilters ();
void lowpass     (real_t*);

// -- integrate.C:

void integrate (Domain*, LESAnalyser*, SumIntegrator*);

// -- nonlinear.C:

void nonLinear (Domain*, SumIntegrator*, vector<real_t*>&, vector<real_t>&);
void dynamic   (Domain*, vector<real_t*>&, const bool = true);
#endif
