///////////////////////////////////////////////////////////////////////////////
// condition.C: functions used to evaluate & apply BCs.
//
// Copyright (c) 1994,2003 Hugh Blackburn
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

#include <Sem.h>


Essential::Essential (const char* v,
		      const char* g) :
// ---------------------------------------------------------------------------
// Essential condition constructor, for a condition with a constant real
// value, also sets boundary group name in base class.
// ---------------------------------------------------------------------------
  _value (strtod (v, 0))
{
  strcpy ((_grp = new char [strlen (g) + 1]), g);
}


void Essential::evaluate (const integer  np  ,
			  const integer  id  ,
			  const integer  nz  ,
			  const Element* E   ,
			  const integer  side,
			  const integer  step,
			  const real*    nx  ,
			  const real*    ny  ,
			  real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, _value, tgt, 1);
}


void Essential::set (const integer  side,
		     const integer* bmap,
		     const real*    src ,
		     real*          tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt. 
// ---------------------------------------------------------------------------
{
  const integer  nm    = Geometry::nP() - 1;
  const integer* start = bmap;
  
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


void Essential::sum (const integer  side  ,
		     const integer* bmap  ,
		     const real*    src   ,
		     const real*    weight,
		     real*          work  ,
		     real*          tgt   ) const
// ---------------------------------------------------------------------------
// To be used for natural BC integral terms, has no effect on essential BCs.
// ---------------------------------------------------------------------------
{ }


void Essential::augmentSC (const integer  side  ,
			   const integer  nband ,
			   const integer  nsolve,
			   const integer* bmap  ,
			   const real*    area  ,
			   real*          work  ,
			   real*          tgt   )  const
// ---------------------------------------------------------------------------
// Do nothing, this is for mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void Essential::augmentOp (const integer  side, 
			   const integer* bmap,
			   const real*    area,
			   const real*    src ,
			   real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void Essential::augmentDg (const integer  side, 
			   const integer* bmap,
			   const real*    area,
			   real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void Essential::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, essential:\t%g", _grp, _value);
}


EssentialFunction::EssentialFunction (const char* f,
				      const char* g)
// ---------------------------------------------------------------------------
// Essential condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Essential class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
  strcpy ((_grp      = new char [strlen (g) + 1]), g);
}


void EssentialFunction:: evaluate (const integer  np  ,
				   const integer  id  ,
				   const integer  nz  ,
				   const Element* E   ,
				   const integer  side,
				   const integer  step,
				   const real*    nx  ,
				   const real*    ny  ,
				   real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// ---------------------------------------------------------------------------
{
  E -> sideEval (side, tgt, _function);
}


void EssentialFunction::set (const integer  side,
			     const integer* bmap,
			     const real*    src ,
			     real*          tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt
// (as for Essential class).
// ---------------------------------------------------------------------------
{
  const integer  nm    = Geometry::nP() - 1;
  const integer* start = bmap;
  
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


void EssentialFunction::sum (const integer  side  ,
			     const integer* bmap  ,
			     const real*    src   ,
			     const real*    weight,
			     real*          work  ,
			     real*          tgt   ) const
// ---------------------------------------------------------------------------
// Take no action on essential BC.
// ---------------------------------------------------------------------------
{ }


void EssentialFunction::augmentSC (const integer  side  ,
				   const integer  nband ,
				   const integer  nsolve,
				   const integer* bmap  ,
				   const real*    area  ,
				   real*          work  ,
				   real*          tgt   )  const
// ---------------------------------------------------------------------------
// Do nothing for essential BCs.
// ---------------------------------------------------------------------------
{ }


void EssentialFunction::augmentOp (const integer  side, 
				   const integer* bmap,
				   const real*    area,
				   const real*    src ,
				   real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void EssentialFunction::augmentDg (const integer  side, 
				   const integer* bmap,
				   const real*    area,
				   real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void EssentialFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, essential-function:\t%s", _grp,_function);
}


Natural::Natural (const char* v,
		  const char* g) :
// ---------------------------------------------------------------------------
// Used to apply natural type BCs using a constant real value.
// ---------------------------------------------------------------------------
  _value (strtod (v, 0))
{
  strcpy ((_grp = new char [strlen (g) + 1]), g);
}


void Natural::evaluate (const integer  np  ,
			const integer  id  ,
			const integer  nz  ,
			const Element* E   ,
			const integer  side,
			const integer  step,
			const real*    nx  ,
			const real*    ny  ,
			real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, _value, tgt, 1);
}


void Natural::set (const integer  side,
		   const integer* bmap,
		   const real*    src ,
		   real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since natural BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void Natural::sum (const integer  side  ,
		   const integer* bmap  ,
		   const real*    src   ,
		   const real*    weight,
		   real*          work  ,
		   real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// Work vector must be np long.
// ---------------------------------------------------------------------------
{ 
  const integer  np    = Geometry::nP();
  const integer  nm    = np - 1;
  const integer* start = bmap;
  
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


void Natural::augmentSC (const integer  side  ,
			 const integer  nband ,
			 const integer  nsolve,
			 const integer* bmap  ,
			 const real*    area  ,
			 real*          work  ,
			 real*          tgt   )  const
// ---------------------------------------------------------------------------
// Do nothing.
// ---------------------------------------------------------------------------
{ }


void Natural::augmentOp (const integer  side, 
			 const integer* bmap,
			 const real*    area,
			 const real*    src ,
			 real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void Natural::augmentDg (const integer  side, 
			 const integer* bmap,
			 const real*    area,
			 real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void Natural::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, natural:\t%g", _grp, _value);
}


NaturalFunction::NaturalFunction (const char* f,
				  const char* g)
// ---------------------------------------------------------------------------
// Natural condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Natural class.
// ---------------------------------------------------------------------------
{
  strcpy ((_function = new char [strlen (f) + 1]), f);
  strcpy ((_grp      = new char [strlen (g) + 1]), g);
}


void NaturalFunction::evaluate (const integer  np  ,
				const integer  id  ,
				const integer  nz  ,
				const Element* E   ,
				const integer  side,
				const integer  step,
				const real*    nx  ,
				const real*    ny  ,
				real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// ---------------------------------------------------------------------------
{
  E -> sideEval (side, tgt, _function);
}


void NaturalFunction::set (const integer  side,
			   const integer* bmap,
			   const real*    src ,
			   real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since natural BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void NaturalFunction::sum (const integer  side  ,
			   const integer* bmap  ,
			   const real*    src   ,
			   const real*    weight,
			   real*          work  ,
			   real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const integer  np    = Geometry::nP();
  const integer  nm    = np - 1;
  const integer* start = bmap;
  
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


void NaturalFunction::augmentSC (const integer  side  ,
				 const integer  nband ,
				 const integer  nsolve,
				 const integer* bmap  ,
				 const real*    area  ,
				 real*          work  ,
				 real*          tgt   )  const
// ---------------------------------------------------------------------------
// Do nothing.
// ---------------------------------------------------------------------------
{ }


void NaturalFunction::augmentOp (const integer  side, 
				 const integer* bmap,
				 const real*    area,
				 const real*    src ,
				 real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void NaturalFunction::augmentDg (const integer  side, 
				 const integer* bmap,
				 const real*    area,
				 real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void NaturalFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, natural-function:\t%s", _grp, _function);
}


NaturalHOPBC::NaturalHOPBC (const char* g)
// ---------------------------------------------------------------------------
// Construct by initializing base Condition class.
// ---------------------------------------------------------------------------
{ 
  strcpy ((_grp = new char [strlen (g) + 1]), g);
}


void NaturalHOPBC::evaluate (const integer  np  ,
			     const integer  id  ,
			     const integer  nz  ,
			     const Element* E   ,
			     const integer  side,
			     const integer  step,
			     const real*    nx  ,
			     const real*    ny  ,
			     real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external via a call to PBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  PBCmgr::evaluate (id, np, nz, step, nx, ny, tgt); 
}


void NaturalHOPBC::set (const integer  side,
			const integer* bmap,
			const real*    src ,
			real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since natural BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void NaturalHOPBC::sum (const integer  side  ,
			const integer* bmap  ,
			const real*    src   ,
			const real*    weight,
			real*          work  ,
			real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  const integer  np    = Geometry::nP();
  const integer  nm    = np - 1;
  const integer* start = bmap;
  
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


void NaturalHOPBC::augmentSC (const integer  side  ,
			      const integer  nband ,
			      const integer  nsolve,
			      const integer* bmap  ,
			      const real*    area  ,
			      real*          work  ,
			      real*          tgt   )  const
// ---------------------------------------------------------------------------
// Do nothing.
// ---------------------------------------------------------------------------
{ }


void NaturalHOPBC::augmentOp (const integer  side, 
			      const integer* bmap,
			      const real*    area,
			      const real*    src ,
			      real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void NaturalHOPBC::augmentDg (const integer  side, 
			      const integer* bmap,
			      const real*    area,
			      real*          tgt ) const
// ---------------------------------------------------------------------------
// Do nothing, this is for Mixed BCs only.
// ---------------------------------------------------------------------------
{ }


void NaturalHOPBC::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, high-order-pressure", _grp);
}


Mixed::Mixed (const char* v,
	      const char* g)
// ---------------------------------------------------------------------------
// The format for a Mixed BC is: "field = mulvalue, refvalue".
// Each of these is expected to evaluate to a real constant.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Mixed::Mixed";
  const char sep[] = ";,";
  char       buf[StrMax], *tok;

  strcpy ((_grp = new char [strlen (g) + 1]), g);

  strcpy (buf, v);
  tok = strtok (buf, sep);

  _K_ = Femlib::value (tok);

  if (_K_ <= EPSSP) {
    sprintf (buf, "transfer coefficient K must be positive (%s)", tok);
    message (routine, buf, ERROR);
  }

  _C_ = Femlib::value (tok = strtok (0, sep));
}


void Mixed::evaluate (const integer  np  ,
		      const integer  id  ,
		      const integer  nz  ,
		      const Element* E   ,
		      const integer  side,
		      const integer  step,
		      const real*    nx  ,
		      const real*    ny  ,
		      real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constants.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, _K_ * _C_, tgt, 1);
}


void Mixed::set (const integer  side,
		 const integer* bmap,
		 const real*    src ,
		 real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since mixed BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void Mixed::sum (const integer  side  ,
		 const integer* bmap  ,
		 const real*    src   ,
		 const real*    weight,
		 real*          work  ,
		 real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*C (supplied as src) to RHS forcing.
// ---------------------------------------------------------------------------
{ 
  const integer  np    = Geometry::nP();
  const integer  nm    = np - 1;
  const integer* start = bmap;
  
  switch (side) {
  case 1: start += nm;           break;
  case 2: start += nm + nm;      break;
  case 3: start += nm + nm + nm; break;
  default: break;
  }

  Veclib::smul (np, _K_ * _C_,  weight, 1, work, 1);

  Veclib::scatr_sum (nm, work, start, tgt);
  if   (side == 3) tgt[bmap [ 0]] += work[nm];
  else             tgt[start[nm]] += work[nm];
}


void Mixed::augmentSC (const integer  side  ,
		       const integer  nband ,
		       const integer  nsolve,
		       const integer* bmap  ,
		       const real*    area  ,
		       real*          work  ,
		       real*          tgt   )  const
// ---------------------------------------------------------------------------
// Add <K, w> terms to (banded LAPACK) matrix tgt.
// ---------------------------------------------------------------------------
{
  const integer    np    = Geometry::nP();
  const integer    nm    = np - 1;
  const integer*   start = bmap;
  register integer i, k;
  
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


void Mixed::augmentOp (const integer  side, 
		       const integer* bmap,
		       const real*    area,
		       const real*    src ,
		       real*          tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  const integer    np    = Geometry::nP();
  const integer    nm    = np - 1;
  const integer*   start = bmap;
  register integer i;
  
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


void Mixed::augmentDg (const integer  side, 
		       const integer* bmap,
		       const real*    area,
		       real*          tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  const integer    np    = Geometry::nP();
  const integer    nm    = np - 1;
  const integer*   start = bmap;
  register integer i;
  
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
  sprintf (tgt, "boundary-group: %s, mixed:\t%g\t%g", _grp, _K_, _C_);
}
