///////////////////////////////////////////////////////////////////////////////
// boundary.C: implement Boundary class functions.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// SYNOPSIS
// --------
// Boundaries correspond to domain edges that have boundary conditions
// applied (as opposed to periodic edges).  The ordering of internal
// storage for condition values and geometric factors corresponds to
// CCW traverse of 2D element edges.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem_h>


void Boundary::evaluate (const integer plane,
			 const integer step ,
			 real*         tgt  ) const
// ---------------------------------------------------------------------------
// Load boundary condition storage area with numeric values.
// ---------------------------------------------------------------------------
{
  _bcond -> evaluate (_np, _id, plane, _elmt, _side, step, _nx, _ny, tgt);
}


void Boundary::set (const real*    src,
		    const integer* b2g,
		    real*          tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to over-ride (set) values
// in globally-numbered tgt.  This will only take place on essential BCs.
//
// b2g is a pointer to the global node numbers for the appropriate
// element's edge nodes.
// ---------------------------------------------------------------------------
{
  _bcond -> set (_side, b2g, src, tgt);
}


void Boundary::sum (const real*    src,
		    const integer* b2g,
		    real*          wrk,
		    real*          tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to add in the boundary-integral
// terms generated in constructing the weak form of the MWR into globally-
// numbered tgt.  This will only take place on natural BCs.
//
// b2g is a pointer to the global node numbers for the appropriate
// element's edge nodes.  wrk is a work array, np long.
// ---------------------------------------------------------------------------
{
  _bcond -> sum (_side, b2g, src, _area, wrk, tgt);
}


void Boundary::augmentSC (const integer  nband ,
			  const integer  nsolve,
			  const integer* b2g   ,
			  real*          work  ,
			  real*          H     ) const
// ---------------------------------------------------------------------------
// Add in diagonal terms <K, w> to (banded LAPACK) H on mixed BCs.
// Work array must be np long.
// ---------------------------------------------------------------------------
{
  _bcond -> augmentSC (_side, nband, nsolve, b2g + bOff(), _area, work, H);
}


void Boundary::augmentOp (const integer* b2g,
			  const real*    src,
			  real*          tgt) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  _bcond -> augmentOp (_side, b2g + bOff(), _area, src, tgt);
}


void Boundary::augmentDg (const integer* b2g,
			  real*          tgt) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  _bcond -> augmentDg (_side, b2g + bOff(), _area, tgt);
}


void Boundary::print () const
// ---------------------------------------------------------------------------
// (Debugging) utility to print internal information.
// ---------------------------------------------------------------------------
{
  char info[StrMax];

  cout << "** Boundary id: " << _id + 1 << " -> ";
  cout << _elmt ->  ID() + 1 << "." << _side + 1;
  cout << " (Element id.side)" << endl;
  
  _bcond -> describe (info);

  cout << info << endl;

  cout << "  " << _np << " (number of points along edge)" << endl;
  cout << "         nx             ny             area";
  cout << endl;
  
  printVector (cout, "rrr", _np, _nx, _ny, _area);
}



