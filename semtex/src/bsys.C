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
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


BoundarySys::BoundarySys (BCmgr*                  bcmgr,
			  const vector<Element*>& elmt ,
			  const char              name ) :
// ---------------------------------------------------------------------------
// Construct internal storage for boundary systems for all modes.
// ---------------------------------------------------------------------------
  field_name (name),
  nbound     (bcmgr -> nBCedges()),
  mixed      (0)
{
  const integer           np = Geometry::nP();
  ListIterator<BCtriple*> edge (bcmgr -> getBCedges());
  BCtriple*               BCT;
  const Condition*        C;
  const char*             S;
  char                    buf[StrMax], group;
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
    S     = bcmgr -> groupInfo    (group);
    C     = bcmgr -> getCondition (group, field_name, 0);
    
    C -> describe (buf);
    if (strstr (buf, "mixed")) mixed = 1;
    
    boundary[0][i] =
    boundary[1][i] =
    boundary[2][i] = new Boundary (i, S, C, elmt[j], k);
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
    S     = bcmgr -> groupInfo    (group);
    C     = bcmgr -> getCondition (group, field_name, 1);

    if (strstr (S, "axis"))
      boundary[1][i] = new Boundary (i, "axis", C, elmt[j], k);
  }

  // -- Mode 2 boundaries,adjusted on axis.
  
  edge.reset();
  for (offset = 0, i = 0; i < nbound; i++, edge.next(), offset += np) {
    BCT   = edge.current();
    group = BCT -> group;
    j     = BCT -> elmt;
    k     = BCT -> side;
    S     = bcmgr -> groupInfo    (group);
    C     = bcmgr -> getCondition (group, field_name, 2);

    if (strstr (bcmgr -> groupInfo (group), "axis"))
      boundary[2][i] = new Boundary (i, "axis", C, elmt[j], k);
  }
}


const vector<Boundary*>& BoundarySys::BCs (const integer mode) const
// ---------------------------------------------------------------------------
// Return appropriate vector of Boundary*'s for field name, according
// to Fourier mode.  Mode number is the actual number, counting from
// zero, not modulo number of modes on this process.
// ---------------------------------------------------------------------------
{
  return boundary [clamp (mode, (const integer)0, (const integer)2)];
}


const NumberSys* BoundarySys::Nsys (const integer mode) const
// ---------------------------------------------------------------------------
// Return appropriate NumberSystem* for field name, according to
// Fourier mode.
// ---------------------------------------------------------------------------
{
  return number [clamp (mode, (const integer)0, (const integer)2)];
}


const real* BoundarySys::Imass (const integer mode) const
// ---------------------------------------------------------------------------
// Return appropriate inverse mass matrix for field name, according to
// Fourier mode.
// ---------------------------------------------------------------------------
{
  return number [clamp (mode, (const integer)0, (const integer)2)] -> imass();
}
