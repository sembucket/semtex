///////////////////////////////////////////////////////////////////////////////
// condition.C: functions used to evaluate & apply BCs.
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
//
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>


Essential::Essential (const char* v,
		      const char* g) : value (strtod (v, 0))
// ---------------------------------------------------------------------------
// Essential condition constructor, for a condition with a constant real
// value, also sets boundary group name in base class.
// ---------------------------------------------------------------------------
{
  grp = strdup (g);
}


void Essential::evaluate (const int      np  ,
			  const int      id  ,
			  const int      nz  ,
			  const Element* E   ,
			  const int      side,
			  const int      step,
			  const real*    nx  ,
			  const real*    ny  ,
			  real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, value, tgt, 1);
}


void Essential::set (const Element* E   ,
		     const int      side,
		     const int*     bmap,
		     const real*    src ,
		     real*          tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt. 
// ---------------------------------------------------------------------------
{
  E -> sideSet (side, bmap, src, tgt);
}


void Essential::sum (const Element* E     ,
		     const int      side  ,
		     const int*     bmap  ,
		     const real*    src   ,
		     const real*    weight,
		     real*          tgt   ) const
// ---------------------------------------------------------------------------
// To be used for natural BC integral terms, has no effect on essential BCs.
// ---------------------------------------------------------------------------
{ }


void Essential::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, essential:\t%g", grp, value);
}


EssentialFunction::EssentialFunction (const char* f,
				      const char* g) : function (strdup (f))
// ---------------------------------------------------------------------------
// Essential condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Essential class.
// ---------------------------------------------------------------------------
{
  grp = strdup (g);
}


void EssentialFunction:: evaluate (const int      np  ,
				   const int      id  ,
				   const int      nz  ,
				   const Element* E   ,
				   const int      side,
				   const int      step,
				   const real*    nx  ,
				   const real*    ny  ,
				   real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// ---------------------------------------------------------------------------
{
  E -> sideEval (side, tgt, function);
}


void EssentialFunction::set (const Element* E   ,
			     const int      side,
			     const int*     bmap,
			     const real*    src ,
			     real*          tgt ) const
// ---------------------------------------------------------------------------
// Scatter external value storage area src into globally-numbered tgt
// (as for Essential class).
// ---------------------------------------------------------------------------
{
  E -> sideSet (side, bmap, src, tgt);
}


void EssentialFunction::sum (const Element* E     ,
			     const int      side  ,
			     const int*     bmap  ,
			     const real*    src   ,
			     const real*    weight,
			     real*          tgt   ) const
// ---------------------------------------------------------------------------
// Take no action on essential BC.
// ---------------------------------------------------------------------------
{ 

}


void EssentialFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, essential-function:\t%s", grp, function);
}


Natural::Natural (const char* v,
		  const char* g) : value (strtod (v, 0))
// ---------------------------------------------------------------------------
// Used to apply natural type BCs using a constant real value.
// ---------------------------------------------------------------------------
{
  grp = strdup (g);
}


void Natural::evaluate (const int      np  ,
			const int      id  ,
			const int      nz  ,
			const Element* E   ,
			const int      side,
			const int      step,
			const real*    nx  ,
			const real*    ny  ,
			real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constant.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, value, tgt, 1);
}


void Natural::set (const Element* E   ,
		   const int      side,
		   const int*     bmap,
		   const real*    src ,
		   real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since natural BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void Natural::sum (const Element* E     ,
		   const int      side  ,
		   const int*     bmap  ,
		   const real*    src   ,
		   const real*    weight,
		   real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  E -> sideDsSum (side, bmap, src, weight, tgt);
}


void Natural::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, natural:\t%g", grp, value);
}


NaturalFunction::NaturalFunction (const char* f,
				  const char* g) : function (strdup (f))
// ---------------------------------------------------------------------------
// Natural condition that sets value by interpreting function string,
// which can be an explicit function of x, y, z & t as well as any
// installed token values.  Otherwise the same as Natural class.
// ---------------------------------------------------------------------------
{
  grp = strdup (g);
}


void NaturalFunction::evaluate (const int      np  ,
				const int      id  ,
				const int      nz  ,
				const Element* E   ,
				const int      side,
				const int      step,
				const real*    nx  ,
				const real*    ny  ,
				real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external tgt by interpreting function.
// ---------------------------------------------------------------------------
{
  E -> sideEval (side, tgt, function);
}


void NaturalFunction::set (const Element* E   ,
			   const int      side,
			   const int*     bmap,
			   const real*    src ,
			   real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since natural BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ 

}


void NaturalFunction::sum (const Element* E     ,
			   const int      side  ,
			   const int*     bmap  ,
			   const real*    src   ,
			   const real*    weight,
			   real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  E -> sideDsSum (side, bmap, src, weight, tgt);
}


void NaturalFunction::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, natural-function:\t%g", grp, function);
}


NaturalHOPBC::NaturalHOPBC (const char* g)

// ---------------------------------------------------------------------------
// Construct by initializing base Condition class.
// ---------------------------------------------------------------------------
{ 
  grp = strdup (g);
}


void NaturalHOPBC::evaluate (const int      np  ,
			     const int      id  ,
			     const int      nz  ,
			     const Element* E   ,
			     const int      side,
			     const int      step,
			     const real*    nx  ,
			     const real*    ny  ,
			     real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external via a call to PBCmgr to compute terms.
// ---------------------------------------------------------------------------
{
  PBCmgr::evaluate (id, np, nz, step, nx, ny, tgt); 
}


void NaturalHOPBC::set (const Element* E   ,
			const int      side,
			const int*     bmap,
			const real*    src ,
			real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since natural BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void NaturalHOPBC::sum (const Element* E     ,
			const int      side  ,
			const int*     bmap  ,
			const real*    src   ,
			const real*    weight,
			real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.
// ---------------------------------------------------------------------------
{ 
  E -> sideDsSum (side, bmap, src, weight, tgt);
}


void NaturalHOPBC::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, high-order-pressure", grp);
}


Mixed::Mixed (const char* v,
	      const char* g)
// ---------------------------------------------------------------------------
// The format for a Mixed BC is: "field = mulvalue, refvalue".
// ---------------------------------------------------------------------------
{
  char buf[StrMax], *tok, *sep = "=,";

  grp = strdup (g);

  strcpy (buf, v);
  tok = strtok (buf, sep);
  K = atof (tok = strtok (0, sep));
  C = atof (tok = strtok (0, sep));
}


void Mixed::evaluate (const int      np  ,
		      const int      id  ,
		      const int      nz  ,
		      const Element* E   ,
		      const int      side,
		      const int      step,
		      const real*    nx  ,
		      const real*    ny  ,
		      real*          tgt ) const
// ---------------------------------------------------------------------------
// Load external value storage area tgt with installed constants.
// ---------------------------------------------------------------------------
{
  Veclib::fill (np, K * C, tgt, 1);
}


void Mixed::set (const Element* E   ,
		 const int      side,
		 const int*     bmap,
		 const real*    src ,
		 real*          tgt ) const
// ---------------------------------------------------------------------------
// Take no action, since mixed BCs are applied in the weak sense.
// ---------------------------------------------------------------------------
{ }


void Mixed::sum (const Element* E     ,
		 const int      side  ,
		 const int*     bmap  ,
		 const real*    src   ,
		 const real*    weight,
		 real*          tgt   ) const
// ---------------------------------------------------------------------------
// Add boundary-integral terms into globally-numbered tgt.  This is
// used to add K*C (supplied as src) to RHS forcing.
// ---------------------------------------------------------------------------
{ 
  E -> sideDsSum (side, bmap, src, weight, tgt);
}


void Mixed::describe (char* tgt) const
// ---------------------------------------------------------------------------
// Load descriptive/diagnostic material into tgt.
// ---------------------------------------------------------------------------
{
  sprintf (tgt, "boundary-group: %s, mixed:\t%g\t%g", grp, K, C);
}
