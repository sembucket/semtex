///////////////////////////////////////////////////////////////////////////////
// global.C: routines for the Global class.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>

real        Global::refCoeff  = 0.0;
real        Global::C1        = 0.5;
real        Global::C2        = 0.3;
real        Global::C3        = 0.2;
real        Global::gblSize   = 1.0;
int         Global::nodeIdMax = 0;
int         Global::loopIdMax = 0;
int         Global::verbose   = 0;
List<Node*> Global::nodeList;
List<Quad*> Global::quadList;


real Global::limits (Point& Pmin,
		     Point& Pmax)
// ---------------------------------------------------------------------------
// Establish Pmin & Pmax that define bounding box for nodes, return hypotenuse.
// ---------------------------------------------------------------------------
{
  ListIterator<Node*> n (nodeList);
  Node*               N;
  real                X, Y, xmin, ymin, xmax, ymax;

  xmin = ymin =  FLT_MAX;
  xmax = ymax = -FLT_MAX;

  for (n.reset(); n.more(); n.next()) {
    N = n.current();
    X = N -> pos () . x;
    Y = N -> pos () . y;
    if      (X < xmin) xmin = X;
    else if (X > xmax) xmax = X;
    if      (Y < ymin) ymin = Y;
    else if (Y > ymax) ymax = Y;
  }

  Pmin.x = xmin;  Pmin.y = ymin;
  Pmax.x = xmax;  Pmax.y = ymax;

  return gblSize = hypot (Pmax.x - Pmin.x, Pmax.y - Pmin.y);
}


Node* Global::exist (const Node* N)
// ---------------------------------------------------------------------------
// Check if a Node corresponding to N has already been created.
// Return pointer to the existing Node if it has, else zero.
// ---------------------------------------------------------------------------
{
  char           err[StrMax], routine[] = "Global::exist";
  int            found = 0;
  register Node* oldNode;
  const Point    P    = N -> pos();
  const real     size = lengthScale();
  const real     TOL  = 0.0001;
  
  ListIterator<Node*> n (nodeList);

  for (n.reset(); !found && n.more(); n.next()) {
    oldNode = n.current();
    found   = oldNode -> pos().distance (P) / size < TOL;
  }

  if (found) {
    sprintf (err, "position for Node %1d exists, deleting", N -> ID());
    message (routine, err, WARNING);
    return oldNode;
  }
  
  return 0;
}
