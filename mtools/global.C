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
list<Node*> Global::nodeList;
list<Quad*> Global::quadList;


real Global::limits (Point& Pmin,
		     Point& Pmax)
// ---------------------------------------------------------------------------
// Establish Pmin & Pmax that define bounding box for nodes, return hypotenuse.
// ---------------------------------------------------------------------------
{
  Node*                 N;
  list<Node*>::iterator n;
  real                  X, Y, xmin, ymin, xmax, ymax;

  xmin = ymin =  1.0e30;
  xmax = ymax = -1.0e30;

  for (n = nodeList.begin(); n != nodeList.end(); n++) {
    N = *n;
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
// Check if a Node corresponding in position to N has already been created.
// Return pointer to the existing Node if it has, else zero.
// ---------------------------------------------------------------------------
{
  char           err[StrMax], routine[] = "Global::exist";
  int            found = 0;
  register Node* oldNode;
  const Point    P    = N -> pos();
  const real     size = lengthScale();
  const real     TOL  = 0.0001;
  
  list<Node*>::iterator n;

  for (n = nodeList.begin(); !found && n != nodeList.end(); n++) {
    oldNode = *n;
    found   = oldNode -> pos().distance (P) / size < TOL;
  }

  if (found) {
    sprintf (err, "position for Node %1d exists, deleting", N -> ID());
    message (routine, err, WARNING);
    return oldNode;
  }
  
  return 0;
}
