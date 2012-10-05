//////////////////////////////////////////////////////////////////////////////
// condition.C: functions used to evaluate & apply BCs.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
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
// along with Semtex (see the file COPYING); if not, write to the Freem
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
// --
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
// described in KIO91.
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


void Essential::get (const int_t   side,
		     const int_t*  bmap,
		     const real_t* src ,
		     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Gather globally-numbered src into external value storage area tgt. 
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

  Veclib::gathr (np, src, start, tgt);
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

void EssentialFunction::get (const int_t   side,
			     const int_t*  bmap,
			     const real_t* src ,
			     real_t*       tgt ) const
// -----------------------------------------------------------------------------
// Gather globally-numbered src into external value storage area tgt. 
// -----------------------------------------------------------------------------
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

  Veclib::gathr (np, src, start, tgt);
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
// Add boundary-integral terms into globally-numbered tgt.
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
// v is not used here. _K_ and _C_ are zeroed and will be updated later
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Mixed::Mixed";
  const char  sep[] = ";,";
  const int_t np  = Geometry::nP();
  const int_t nb  = Geometry::nBnode();
   
  _K_  = new real_t [static_cast<size_t>(np)];
  _C_  = new real_t [static_cast<size_t>(np)];
  k_func = new real_t [static_cast<size_t>(nb)];
  Veclib::zero (np, _C_, 1);
  Veclib::zero (np, _K_, 1); 

  // TESTING::: Veclib::zero (np, _K_, 1);
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
// ---------------------------------------------------------------------------
{
  Veclib::vmul(np, _K_, 1,  _C_, 1, tgt, 1);
}


void Mixed::sum (const int_t   side  ,
		 const int_t*  bmap  ,
		 const real_t* src   ,
		 const real_t* weight,
		 real_t*       work  ,
		 real_t*       tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*C (supplied as src) to RHS forcing.
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

  Veclib::vmul (np, _K_, 1, area, 1, work, 1);

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
    tgt[start[i]] += _K_[i] * area[i] * src[start[i]];
  
  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_[nm] * area[nm] * src[i];
  //	printf(								"%f         \n", *( _K_+1));
}



void Mixed::augmentDg (const int_t   side, 
		       const int_t*  bmap,
		       const real_t* area,
		       real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs. Add in diagonal terms <K, w> to globally-numbered tgt.
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
    tgt[start[i]] += _K_[i] * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_[nm] * area[nm];

}


void Mixed::switchK (const int_t   np     ,
		     const real_t* K_adjoint_timei    ,
		     const real_t* K_adjoint_timej    , //useless
		     const int_t   doffset,
		     const int_t   dskip  ,
		     const real_t* nx    ,
		     const real_t* ny, 
		     const bool forwards  ) const
// ---------------------------------------------------------------------------
// Load field data from appropriate edge into _K_.
// Only happen for Toutlfow BC edge.
// ---------------------------------------------------------------------------
{   

  Veclib::copy(np, k_func+doffset*np, 1, _K_, 1);
}

void Mixed::switchK (const int_t num) const
// ---------------------------------------------------------------------------
// Change _K_ based on position of current Boundary object. Used if _K_ is a
// spatially varying function.
// ---------------------------------------------------------------------------
{
  const int_t np = Geometry::nP();

  Veclib::copy(np, k_func+num*np, 1, _K_, 1);

}
void Mixed::takeC ( const int_t np,
		    const real_t* p_adjoint ) const
// ---------------------------------------------------------------------------
// update _C_
// ---------------------------------------------------------------------------
{   
  Veclib::vdiv (np, p_adjoint, 1, _K_, 1, _C_, 1);
} 

void Mixed::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "mixed:\t%g\t%g", _K_[0], _C_[0]);
}

void  Controlbc::set (const int_t   side,
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

void Controlbc::get (const int_t   side,
		     const int_t*  bmap,
		     const real_t* src ,
		     real_t*       tgt ) const
// ---------------------------------------------------------------------------
// Gather globally-numbered src into external value storage area tgt. 
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

  Veclib::gathr (np, src, start, tgt);
}

real_t Controlbc::controlnorm( const real_t* weight,
			       real_t*       uci   ) const
//Integrate over the control bc to get the norm.
{
  const int_t np = Geometry::nP();
  int_t       i;
  real_t      norm = 0.0;
  for (i = 0; i < np; i++) 
    norm += *(uci+i) * *(uci+i) * *(weight+i);
  return norm;
}



real_t Controlbc::controllength( const real_t* weight) const
//Integrate over the control bc to get the norm.
{
const int_t   np    = Geometry::nP();
int_t i;
real_t length = 0.0;
   for (i = 0; i < np; i++) length +=  weight[i];
 return length;
}




void Controlbc::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "controlbc:\t 0");
}



void Toutflow::evaluate (const int_t    np  ,
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
// ---------------------------------------------------------------------------
{
   Veclib::zero (np, tgt, 1);
}


void Toutflow::augmentSC (const int_t   side  ,
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
  Veclib::vmul (np, _K_, 1, area, 1, work, 1);

  for (i = 0; i < nm; i++)
    if ((k = start[i]) < nsolve)
      tgt[Lapack::band_addr (k, k, nband)] += work[i];

  i = (side == 3) ? bmap[0] : start[nm];
  if (i < nsolve) tgt[Lapack::band_addr (i, i, nband)] += work[nm];
}



void Toutflow::augmentOp (const int_t   side, 
			  const int_t*  bmap,
			  const real_t* area,
			  const real_t* src ,
			  real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are Toutflow BCs.  Add in diagonal terms
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
    tgt[start[i]] += _K_[i] * area[i] * src[start[i]];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_[nm] * area[nm] * src[i];
}


void Toutflow::augmentDg (const int_t   side, 
			  const int_t*  bmap,
			  const real_t* area,
			  real_t*       tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are Toutflow
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
    tgt[start[i]] += _K_[i] * area[i];

  i = (side == 3) ? bmap[0] : start[nm];
  tgt[i] += _K_[nm]* area[nm];

}




void Toutflow::switchK (const int_t   np     ,
			const real_t* Ux    ,
			const real_t* Uy    ,
			const int_t   doffset,
			const int_t   dskip  ,
			const real_t* nx    ,
			const real_t* ny, 
			const bool forwards  ) const
// ---------------------------------------------------------------------------
// Load field data from appropriate edge into _K_.
// Only happen for Toutlfow BC edge.
// ---------------------------------------------------------------------------
{  real_t* Uxedge;
   real_t* Uyedge;
   Uxedge= new real_t [static_cast<size_t>(np)];
   Uyedge= new real_t [static_cast<size_t>(np)];
   
 if (!forwards){
   Veclib::copy (np, Ux + doffset, dskip, Uxedge, 1);
   Veclib::copy (np, Uy + doffset, dskip, Uyedge, 1);
   
   Veclib::vvtvvtp(np, Uxedge, 1, nx, 1, Uyedge, 1, ny, 1, _K_, 1);
   Veclib::smul(np, 1.0/Femlib::value ("KINVIS"), _K_, 1, _K_, 1);	 
//	 *(_K_+0)=100; 	*(_K_+1)=60; 	*(_K_+2)=160;	*(_K_+3)=80;	*(_K_+4)=100; 
 }
 else{
   Veclib::zero (np, _K_, 1);
//	 	 *(_K_+0)=1; 	*(_K_+1)=1; 	*(_K_+2)=1;	*(_K_+3)=1;	*(_K_+4)=1; 
 }
//	 printf("%e,  %e,  %e,  %e,  %e, \n", *(_K_+0),*(_K_+1),*(_K_+2),*(_K_+3),*(_K_+4) );
	
   delete [] Uxedge;
   delete [] Uyedge;
}



void Toutflow::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "toutflow:-UnRe for backwards");
}
