//////////////////////////////////////////////////////////////////////////////
// condition.C: functions used to evaluate & apply BCs.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// All classes inherit the (semi-abstract) base class Condition.
// Owing to the different behaviour of essential and natural BCs, a
// number of routines take no action; due also to the generality of
// the base class function calls, many routines do not use all
// parameters.
//
// The two basic kinds of conditions are essential and natural.
// Essential conditions are set/imposed directly, while natural
// conditions are applied by summation of integral terms, owing to the
// fact that they derive from integration by parts.  So function "set"
// is only used by essential-type conditions while function "sum" is
// only used by natural-type conditions.
//
// The HOPBC condition are a special case of natural conditions used
// for the pressure Poisson equation: they are derived from the
// momentum equations and computed using an extrapolative process
// described in [1].
//
// A further general class of available condition is the Robin or
// mixed type, of form dc/dn + K( c - C ) = 0. See [2]. 

// The convective BC is a variant of this, see [3]. The wave equation
//
// u_t + V u_x = 0
//
// becomes in discrete form with backwards Euler for temporal
// derivative
//
// [u^(n+1) - u^(n)]*DT^{-1} + V u_x^(n+1) = 0
//
// which is rearranged as a Robin BC:
//
// u_x^(n+1) + (V*DT)^{-1} [ u^(n+1) - u^(n) ] = 0
//
// So in the typical Robin form
//
// u_n + K ( u - C ) = 0
//
// K = V*DT and C = u^(n)
//
// Convective BCs are good at containing blowups related to inflow
// across the outflow boundary.  Guess V to be the mean convection
// velocity.
//
// NB: if you invent a new BC which is of essential/Dirichlet type,
// you also need to edit mesh.C/buildMask() so that enumerate sets the
// global numbering mask to 1 for this type. Otherwise the BC won't be
// correctly applied.
//
// References:
// 1. Karniadakis, Israeli & Orszag (1991) "High-order splitting methods
//    for the incompressible Navier--Stokes equations", JCP 9(2).
// 2. Blackburn (2001) "Dispersion and diffusion in coated tubes of 
//    arbitrary cross-section", Comput Chem Eng 25(2/3).
// 3. Orlanski (1976) "A simple boundary condition for unbounded 
//    hyperbolic flows", JCP 21.
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


void Essential::evaluate (const int_t    np  ,
			  const int_t    id  ,
			  const int_t    nz  ,
			  const Element* E   ,
			  const int_t    side,
			  const int_t    step,
			  const real_t*  nx  ,
			  const real_t*  ny  ,
			  real_t*        tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, _value, tgt, 1);
}


void Essential::set (const int_t   side,
		     const int_t*  bmap,
		     const real_t* src ,
		     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt. 
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void Essential::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "essential:\t%g", _value);
}


EssentialFunction::EssentialFunction (const char* f)
// ---------------------------------------------------------------------------
// Essential condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Essential class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
}


void EssentialFunction:: evaluate (const int_t    np  ,
				   const int_t    id  ,
				   const int_t    nz  ,
				   const Element* E   ,
				   const int_t    side,
				   const int_t    step,
				   const real_t*  nx  ,
				   const real_t*  ny  ,
				   real_t*        tgt ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// ---------------------------------------------------------------------------
{
  E -> sideEval (side, tgt, _function);
}


void EssentialFunction::set (const int_t   side,
			     const int_t*  bmap,
			     const real_t* src ,
			     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt
// (as for Essential class).
// ---------------------------------------------------------------------------
{
  const int_t  nm    = Geometry::nP() - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::scatr (nm, src, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] = src[nm];
  else             tgt[start[nm]] = src[nm];  
}


void EssentialFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "essential-function:\t%s", _function);
}


void Natural::evaluate (const int_t    np  ,
			const int_t    id  ,
			const int_t    nz  ,
			const Element* E   ,
			const int_t    side,
			const int_t    step,
			const real_t*  nx  ,
			const real_t*  ny  ,
			real_t*        tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, _value, tgt, 1);
}


void Natural::sum (const int_t   side  ,
		   const int_t*  bmap  ,
		   const real_t* src   ,
		   const real_t* weight,
		   real_t*       work  ,
		   real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms (from src) into globally-numbered tgt.
// Work vector must be np long.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void Natural::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "natural:\t%g", _value);
}


NaturalFunction::NaturalFunction (const char* f)
// ---------------------------------------------------------------------------
// Natural condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Natural class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
}


void NaturalFunction::evaluate (const int_t    np  ,
				const int_t    id  ,
				const int_t    nz  ,
				const Element* E   ,
				const int_t    side,
				const int_t    step,
				const real_t*  nx  ,
				const real_t*  ny  ,
				real_t*        tgt ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// ---------------------------------------------------------------------------
{
  E -> sideEval (side, tgt, _function);
}


void NaturalFunction::sum (const int_t   side  ,
			   const int_t*  bmap  ,
			   const real_t* src   ,
			   const real_t* weight,
			   real_t*       work  ,
			   real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "natural-function:\t%s", _function);
}


void NaturalHOPBC::evaluate (const int_t    np  ,
			     const int_t    id  ,
			     const int_t    nz  ,
			     const Element* E   ,
			     const int_t    side,
			     const int_t    step,
			     const real_t*  nx  ,
			     const real_t*  ny  ,
			     real_t*        tgt ) const
// ---------------------------------------------------------------------------
// Load external via a call to PBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  PBCmgr::evaluate (id, np, nz, step, nx, ny, tgt); 
}


void NaturalHOPBC::sum (const int_t   side  ,
			const int_t*  bmap  ,
			const real_t* src   ,
			const real_t* weight,
			real_t*       work  ,
			real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::vmul (np, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];  
}


void NaturalHOPBC::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "high-order-pressure");
}


Mixed::Mixed (const char* v)
// ---------------------------------------------------------------------------
// The format for a Mixed BC is: "field = mulvalue;refvalue", where in
// dc/dn + K( c - C ) = 0, K is mulvalue and C is refvalue.  The
// separator must be ';' or ',' (and without white space). Each of the
// two supplied values is expected to evaluate to a real_t constant.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Mixed::Mixed";
  const char sep[] = ";,";
  char       buf[StrMax], *tok;

  strcpy (buf, v);
  _K_ = Femlib::value (tok = strtok (buf,  sep));
  _C_ = Femlib::value (tok = strtok (NULL, sep));
}


void Mixed::evaluate (const int_t    np  ,
		      const int_t    id  ,
		      const int_t    nz  ,
		      const Element* E   ,
		      const int_t    side,
		      const int_t    step,
		      const real_t*  nx  ,
		      const real_t*  ny  ,
		      real_t*        tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constants.
// This takes no part in the solution 
// ---------------------------------------------------------------------------
{
#if defined (DEBUG)
  Veclib::fill (np, _K_ * _C_, tgt, 1);
#endif
}


void Mixed::sum (const int_t   side  ,
		 const int_t*  bmap  ,
		 const real_t* src   ,
		 const real_t* weight,
		 real_t*       work  ,
		 real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*C to RHS forcing.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_ * _C_, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];
}


void Mixed::augmentSC (const int_t   side  ,
		       const int_t   nband ,
		       const int_t   nsolve,
		       const int_t*  bmap  ,
		       const real_t* area  ,
		       real_t*       work  ,
		       real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void Mixed::augmentOp (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       const real_t* src ,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm] * src[i];
}


void Mixed::augmentDg (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm];
}


void Mixed::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "mixed:\t%g\t%g", _K_, _C_);
}


Convective::Convective (const char* v)
// ---------------------------------------------------------------------------
// The format for a Convective BC is: "field = velocity", e.g. "u =
// V".  But we store K = 1.0/(V*dt).  The reference value C = 0.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Convective::Convective";
  char str[StrMax];

  _K_ = (strtod (v, 0));
  _K_ = 1.0 / (_K_ * Femlib::value ("D_T"));
}


void Convective::sum (const int_t   side  ,
		      const int_t*  bmap  ,
		      const real_t* src   ,
		      const real_t* weight,
		      real_t*       work  ,
		      real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*c^{n} (where c^{n} is supplied as src) to RHS forcing.
// ---------------------------------------------------------------------------
{ 
  const int_t  np    = Geometry::nP();
  const int_t  nm    = np - 1;
  const int_t* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::svvtt (np, _K_, src, 1, weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];
}


void Convective::augmentSC (const int_t   side  ,
			    const int_t   nband ,
			    const int_t   nsolve,
			    const int_t*  bmap  ,
			    const real_t* area  ,
			    real_t*       work  ,
			    real_t*       tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i, k;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}


void Convective::augmentOp (const int_t   side, 
			    const int_t*  bmap,
			    const real_t* area,
			    const real_t* src ,
			    real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm] * src[i];
}


void Convective::augmentDg (const int_t   side, 
			    const int_t*  bmap,
			    const real_t* area,
			    real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    nm    = np - 1;
  const int_t*   start = bmap;
  register int_t i;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  for (i = 0; i < nm; i++)
    tgt[start[i]] += _K_ * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_ * area[nm];
}


void Convective::extract (const int_t   np     ,
			  const real_t* src    ,
			  const int_t   doffset,
			  const int_t   dskip  , 
			  real_t*       tgt    ) const
// ---------------------------------------------------------------------------
// Load field data from appropriate edge into external BC storage tgt.
// This is the same operation as Edge::get, but here will only happen for
// a convective BC edge.
// ---------------------------------------------------------------------------
{
  Veclib::copy (np, src + doffset, dskip, tgt, 1);
}


void Convective::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "convective:\t%g", _K_);
}
