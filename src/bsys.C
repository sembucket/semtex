///////////////////////////////////////////////////////////////////////////////
// bsys.C: BoundarySys class functions.
//
// Copyright (C) 1999 Hugh Blackburn
//
// The information to be returned by class functions are the global
// numbering scheme and vector of Boundary*'s for a given Field and
// Fourier mode number.  There is one BoundarySys for each Field, but
// a possible modal dependence for the appropriate BCs (in fact, only
// for 3D cylindrical coordinate systems in which the axis appears).
//
// Use of cylindrical coordinates is flagged by the Geometry class
// variable.  In the case where the number of space dimensions is also
// 3, the number of boundary frames and numbering systems is set to 3,
// for the 0th, 1st and 2nd (and higher) modes, irrespective of the
// number of Fourier modes actually used.  See BCmgr.C.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>


BoundarySys::BoundarySys (BCmgr*                  bcmgr,
			  const vector<Element*>& elmt ,
			  const char              name ) :
// ---------------------------------------------------------------------------
// Construct internal storage for boundary systems for all modes.
// ---------------------------------------------------------------------------
  field_name (name),
  nbound     (bcmgr -> nBCedges())
{
  const char              routine[] = "BoundarySys::BoundarySys";
  const integer           np = Geometry::nP();
  ListIterator<BCtriple*> edge (bcmgr -> getBCedges());
  char                    group;
  BCtriple*               BCT;
  integer                 i, j, k, offset;

  number = new NumberSys* [3];
  for (i = 0; i < 3; i++) number[i] = bcmgr -> getNumberSys (field_name, i);

  boundary = new vector<Boundary*> [3];

  if (!nbound) { for (i = 0; i < 3; i++) boundary[i].setSize (0); return; }

  // -- Construct vectors of Boundary pointers using BCmgr.

  for (i = 0; i < 3; i++) boundary[i].setSize (nbound);
 
  // -- Mode 0 boundaries, and default settings for other modes.

  edge.reset();
  for (offset = 0, i = 0; i < nbound; i++, edge.next(), offset += np) {
    BCT   = edge.current();
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    boundary[0][i] = boundary[1][i] = boundary[2][i] =
      new Boundary
      (i, offset, bcmgr -> groupInfo (group),
       bcmgr -> getCondition (group, field_name, 0), elmt[j], k);
  }

  if (!(Geometry::system() == Geometry::Cylindrical && Geometry::nDim() == 3))
    return;

  // -- Mode 1 boundaries, adjusted on axis.

  edge.reset();
  for (offset = 0, i = 0; i < nbound; i++, edge.next(), offset += np) {
    BCT   = edge.current();
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    if (strstr (bcmgr -> groupInfo (group), "axis"))
      boundary[1][i] = new Boundary
	(i,offset,"axis",bcmgr->getCondition(group,field_name,1),elmt[j],k);
  }

  // -- Mode 2 boundaries,adjusted on axis.
  
  edge.reset();
  for (offset = 0, i = 0; i < nbound; i++, edge.next(), offset += np) {
    BCT   = edge.current();
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    if (strstr (bcmgr -> groupInfo (group), "axis"))
      boundary[2][i] = new Boundary
	(i,offset,"axis",bcmgr->getCondition(group,field_name,2),elmt[j],k);
  }
}


const vector<Boundary*>& BoundarySys::BCs (integer mode) const
// ---------------------------------------------------------------------------
// Return appropriate vector of Boundary*'s for field name, according
// to Fourier mode.  Mode number is the actual number, counting from
// zero, not modulo number of modes on this process.
// ---------------------------------------------------------------------------
{
  return boundary [clamp (mode, 0, 2)];
}


const NumberSys* BoundarySys::Nsys (integer mode) const
// ---------------------------------------------------------------------------
// Return appropriate NumberSystem* for field name, according to
// Fourier mode.
// ---------------------------------------------------------------------------
{
  return number [clamp (mode, 0, 2)];
}


const real* BoundarySys::Imass (const integer mode) const
// ---------------------------------------------------------------------------
// Return appropriate inverse mass matrix for field name, according to
// Fourier mode.
// ---------------------------------------------------------------------------
{
  return number [clamp (mode, 0, 2)] -> imass();
}
