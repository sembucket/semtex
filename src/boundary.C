/*****************************************************************************
 * boundary.C:  Boundary edge routines.
 *****************************************************************************/

static char RCSid[]= "$Id$";

#include "Fem.h"


Boundary::Boundary (int ident, Element* E, int sideNo,  BC* cond)
// ---------------------------------------------------------------------------
// Constructor.  Allocate new memory for value & geometric factors.
// ---------------------------------------------------------------------------
{
  id        = ident;
  elmt      = E;
  side      = sideNo;
  condition = cond;
  value     = rvector (E -> nKnot());
  nx        = rvector (E -> nKnot());
  ny        = rvector (E -> nKnot());
  area      = rvector (E -> nKnot());

  E -> sideOffset (side, offset, skip);
  E -> sideGeom   (side, nx, ny, area);
}


Boundary::Boundary (const Boundary& B, const List<Element*>& E)
// ---------------------------------------------------------------------------
// Make a copy of an existing Boundary edge, with new data storage area
// and with element pointer into the list of associated new Elements (input).
//
// The boundary condition is set to the input condition.
// ---------------------------------------------------------------------------
{
  char routine[] = "Boundary::Boundary(Boundary&, List<Element*>&)";

  memcpy (this, &B, sizeof (Boundary));

  int found = 0;
  for (ListIterator<Element*> k(E); !found && k.more(); k.next())
    if (found = k.current () -> ID () == elmt -> ID ()) {
      elmt      = k.current ();
      value     = rvector (elmt -> nKnot());
      condition = B.condition;
    }

  if (!found) message (routine, "can't find a match to new element id", ERROR);
}


void  Boundary::evaluate (int step)
// ---------------------------------------------------------------------------
// Load boundary condition storage area with numeric values.
// ---------------------------------------------------------------------------
{
  char routine[] = "Boundary::evaluate";
  int  np        = elmt -> nKnot ();

  switch (condition -> kind) {
  case ESSENTIAL: case NATURAL: case VALUE: case FLUX:
    Veclib::fill (np, condition -> value, value, 1);
    break;
  case ESSENTIAL_FN: case NATURAL_FN: 
    elmt -> sideEval (side, value, condition -> string);
    break;
  case OUTFLOW: case WALL:
    Veclib::fill (np, 0.0, value, 1);
    break;
  case HOPBC:
    PBCmanager::evaluate (id, np, step, value, nx, ny);
    break;
  default:
    message (routine, "illegal BC kind", ERROR);
    break;
  }
}


void  Boundary::print () const
// ---------------------------------------------------------------------------
// (Debugging) utility to print internal information.
// ---------------------------------------------------------------------------
{
  char routine[] = "Boundary::print";

  cout << "** Boundary id: " << id  << " -> ";
  cout <<     elmt ->  ID () << "." << side;
  cout << " (Element id.side)" << endl;
  
  switch (condition -> kind) {
  case ESSENTIAL:
    cout << "ESSENTIAL:    " << condition -> value  << endl;
    break;
  case ESSENTIAL_FN:
    cout << "ESSENTIAL_FN: " << condition -> string << endl;
    break;
  case NATURAL:
    cout << "NATURAL:      " << condition -> value  << endl;
    break;
  case NATURAL_FN:
    cout << "NATURAL_FN:   " << condition -> string << endl;
    break;
  case WALL:
    cout << "WALL"                                  << endl;
    break;
  case OUTFLOW:
    cout << "OUTFLOW"                               << endl;
    break;
  case HOPBC:
    cout << "HOPBC"                                 << endl;
    break;
  default:
    message (routine, "unknown boundary condition kind", ERROR);
    break;
  }
  cout << "  " << elmt -> nKnot() << " (number of points along edge)" << endl;
  cout << "         nx             ny             area           value";
  cout << endl;
  
  printVector (cout, "rrrr", elmt -> nKnot(), nx, ny, area, value);
}


void  Boundary::enforce (real* target) const
// ---------------------------------------------------------------------------
// This is for ESSENTIAL BCs: load BC data into globally-numbered target.
// ---------------------------------------------------------------------------
{
  elmt -> sideScatr (side, value, target);
}


void  Boundary::dsSum (real* target) const
// ---------------------------------------------------------------------------
// This is for NATURAL BCs: direct-stiffness-sum weighted BC data into target.
// ---------------------------------------------------------------------------
{
  elmt -> sideDsSum (side, value, area, target);
}


int  Boundary::isWall () const
// ---------------------------------------------------------------------------
// Return 1 for wall-type BC, 0 otherwise.
// ---------------------------------------------------------------------------
{
  return condition -> kind == WALL;
}


int  Boundary::isEssential () const
// ---------------------------------------------------------------------------
// Return 1 for essential-type BC, 0 for natural.
// ---------------------------------------------------------------------------
{
  if (   condition -> kind == ESSENTIAL
      || condition -> kind == ESSENTIAL_FN
      || condition -> kind == VALUE
      || condition -> kind == WALL        )
    return 1;
  else
    return 0;
}


void  Boundary::mask (int* gmask) const
// ---------------------------------------------------------------------------
// Mask global node-number vector for this boundary.
// ---------------------------------------------------------------------------
{
  elmt -> sideMask (side, gmask);
}


void Boundary::curlCurl (const real*  U ,  const real*  V ,
			 real*        wx,  real*        wy) const
// ---------------------------------------------------------------------------
// Evaluate dw/dx & dw/dy (where w is the z-component of vorticity) from
// element velocity fields, according to the side of the element on which
// they are required.  u & v are Boundaries on the two velocity Fields.
// U & V are pointers to the data storage areas for U & V velocity Fields.
//
// NB: sense of traverse in wx & wy is BLAS-conformant (according to sign
// of skip on relevant edge).
// ---------------------------------------------------------------------------
{
  const int  ntot   = elmt -> nTot ();
  const int  offset = elmt -> nOff ();

  real*  w   = rvector (ntot);
  real*  vx  = rvector (ntot);
  real*  uy  = rvector (ntot);

  Veclib::copy (ntot, U + offset, 1, uy, 1);
  Veclib::copy (ntot, V + offset, 1, vx, 1);
  elmt -> grad (vx, uy);

  // -- Vorticity, w = dv/dx - du/dy.

  Veclib::vsub (ntot, vx, 1, uy, 1, w, 1);

  // -- find dw/dx & dw/dy on appropriate edge.

  elmt -> sideGrad (side, w, wx, wy);

  freeVector (w);
  freeVector (vx);
  freeVector (uy);
}


void  Boundary::resetKind (const BC* new_other, const BC* new_outflow)
// ---------------------------------------------------------------------------
// Examine & reset BC kind for this Boundary.
// This routine is intended for resetting pressure BCs.
// ---------------------------------------------------------------------------
{
  if   (condition -> kind == OUTFLOW) condition = (BC*) new_outflow;
  else                                condition = (BC*) new_other;
}


Vector  Boundary::normalTraction (const real* src, const int skp) const
// ---------------------------------------------------------------------------
// Compute normal tractive force on this boundary segment, using src as
// a pressure stress.
// ---------------------------------------------------------------------------
{
  register int  i, np    = nKnot ();
  Vector        Force    = {0.0, 0.0, 0.0};
  real*         pressure = rvector (np);

  Veclib::copy (np, src, skp, pressure, 1);

  for (i = 0; i < np; i++) {
    Force.x += nx[i] * pressure[i] * area[i];
    Force.y += ny[i] * pressure[i] * area[i];
  }

  freeVector (pressure);

  return Force;
}


Vector  Boundary::tangentTraction (const real*  dx       ,
				   const real*  dy       ,
				   const int    skp      ,
				   const real   mu       ,
				   const int    component) const
// ---------------------------------------------------------------------------
// Compute 1st or 2nd component of viscous stress on this boundary segment.
// ---------------------------------------------------------------------------
{
  register int   i, np = nKnot ();
  Vector         Force = {0.0, 0.0, 0.0};
  real*          ux    = rvector (np);
  real*          uy    = rvector (np);

  Veclib::copy (np, dx, skp, ux, 1);
  Veclib::copy (np, dy, skp, uy, 1);
  
  switch (component) {
  case 1:
    for (i = 0; i < np; i++) {
      Force.x += (2.0*ux[i]*nx[i] + uy[i]*ny[i]) * area[i];
      Force.y +=                    uy[i]*nx[i]  * area[i];
    }
    break;
  case 2:
    for (i = 0; i < np; i++) {
      Force.x +=                    ux[i]*ny[i]  * area[i];
      Force.y += (2.0*uy[i]*ny[i] + ux[i]*nx[i]) * area[i];
    }
    break;
  default:
    break;
  }
  
  freeVector (ux);
  freeVector (uy);

  Force.x *= -mu;
  Force.y *= -mu;

  return Force;
}
