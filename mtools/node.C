///////////////////////////////////////////////////////////////////////////////
// node.cc
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <qmesh.h>


ostream& operator << (ostream& s,
		      Node&    n)
// ---------------------------------------------------------------------------
// Print up node information.
// ---------------------------------------------------------------------------
{
  s << setw (5) << n.id << " " << setw (10) << n.ideal;
  
  switch (n.kind) {
  case Node::INTERIOR_FIXED:         s << "  1 "; break;
  case Node::LOOP_OFFSET_FIXED:      s << "  2 "; break;
  case Node::LOOP_OFFSET_MOBILE:     s << "  3 "; break;
  case Node::LOOP_BOUNDARY_FIXED:    s << "  4 "; break;
  case Node::LOOP_BOUNDARY_MOBILE:   s << "  5 "; break;
  case Node::DOMAIN_OFFSET_FIXED:    s << "  6 "; break;
  case Node::DOMAIN_OFFSET_MOBILE:   s << "  7 "; break;
  case Node::DOMAIN_BOUNDARY_FIXED:  s << "  8 "; break;
  case Node::DOMAIN_BOUNDARY_MOBILE: s << "  9 "; break;
  default:                           s << "  0 "; break;
  }

  s << setw (10) << n.loc.x << "  " << setw (10) << n.loc.y;

  return s;
}


void  Node::setPos (const Point& p)
// ---------------------------------------------------------------------------
// Move a node to a specified location, but only if it is free to do so.
// ---------------------------------------------------------------------------
{ 
  if (kind == INTERIOR) {
    loc.x = p.x; 
    loc.y = p.y;
  } 
}


void Node::xadd (Node* N)
// ---------------------------------------------------------------------------
// Traverse list of contacting nodes.  Add in N if not already present.
// ---------------------------------------------------------------------------
{
  contact.xadd (N);
}


Point Node::centroid () const
// ---------------------------------------------------------------------------
// Return centroid of connected Nodes (weighted).
// ---------------------------------------------------------------------------
{
  char      err[StrMax], routine[] = "Node::centroid";
  const int M = contact.length ();

  Point P, C, D;
  Node* N;
  real  w, den;
  
  P   = 0.0;
  D   = 0.0;
  den = 0.0;
  w   = 0.0;

  if (!M) {
    sprintf (err, "Node %1d has empty list of contacting Nodes", id);
    error (routine, err, WARNING);
    return P;
  }

  for (ListIterator<Node*> n (contact); n.more (); n.next ()) {
    N = n.current ();

    P += N -> pos ();		                // -- Plain centroid.
    w += 1.0;
/*
    w    = N -> pos().distance (this -> pos());	// -- Length-weighted.
    P    = N -> pos() - this -> pos();
    D   += w * P;
    den += w;
*/
  }

  P *= 1.0 / (w);
/*
  P = this -> pos() + 1.0 / den * D;
*/
  return P;
}
