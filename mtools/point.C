///////////////////////////////////////////////////////////////////////////////
// point.C: routines for the Point class.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>


real Point::distance (const Point& p) const
// ---------------------------------------------------------------------------
// Return distance between *this and p.
// ---------------------------------------------------------------------------
{
  real dx = x - p.x;
  real dy = y - p.y;

  return hypot (dx, dy);
}


real Point::angle (const Point& p) const
// ---------------------------------------------------------------------------
// Return angle of line from *this to p, in range -PI--PI.
// ---------------------------------------------------------------------------
{
  double dx = p.x - x;
  double dy = p.y - y;
  
  return atan2 (dy, dx);
}


real Point::angle (const Point& p1,
		   const Point& p2) const
// ---------------------------------------------------------------------------
// Return angle in radians which *this subtends with p1 & p2 in CCW rotation.
// ---------------------------------------------------------------------------
{
  double dx, dy, t1, t2, dir, c, s, a;

  dx = p1.x - x;
  dy = p1.y - y;
  t1 = atan2 (dy, dx);

  dx = p2.x - x;
  dy = p2.y - y;
  t2 = atan2 (dy, dx);

  dir = this -> turn (p1, p2);

  if (dir < 0.0 && t2 < t1) {	      // -- Angular difference lies in 0--PI.
    c = cos (t2 - t1);
    s = sin (t2 - t1);
    a = atan2 (s, c);
  } else if (dir >= 0.0 && t2 < t1) { // -- Angular difference lies in PI--2PI.
    c = cos (t1 - t2);
    s = sin (t1 - t2);
    a = atan2 (s, c);
    a = (a <= 0.0) ? TWOPI + a : TWOPI - a;
  } else			      // -- Other cases don't need fixing.
    a = t2 - t1;

  return a;
}


real Point::turn (const Point& p1,
		  const Point& p2) const
// ---------------------------------------------------------------------------
// The value returned is positive if the traverse from p1 through *this and
// on to p2 makes a CCW turn, otherwise negative.  See Ref. [2], \S\,7.8.3.
// ---------------------------------------------------------------------------
{
  real dx1 = x - p1.x;
  real dy1 = y - p1.y;
  real dx2 = p2.x - x;
  real dy2 = p2.y - y;

  return dx1 * dy2 - dy1 * dx2;
}


int Point::ccw (const Point& p1,
		const Point& p2) const
// ---------------------------------------------------------------------------
// Similar to Point::turn except that +1 returned for CCW turn, -1, 0 etc.
// ---------------------------------------------------------------------------
{
  real dx1, dx2, dy1, dy2, xy1, xy2;
  
  dx1 = p1.x - x;   dy1 = p1.y - y;
  dx2 = p2.x - x;   dy2 = p2.y - y;
  xy1 = dx1  * dy2; xy2 = dy1  * dx2;

  if (               xy1  >   xy2               ) return +1;
  if (               xy1  <   xy2               ) return -1;
  if (     dx1*dx2 < 0.0  ||  dy1*dy2 < 0.0     ) return -1;
  if ( dx1*dx1 + dy1*dy1  <   dx2*dx2 + dy2*dy2 ) return +1;
  
  return 0; 
}


Point& Point::operator += (const Point& rhs)
// ---------------------------------------------------------------------------
// Add x & y locations of rhs to *this.
// ---------------------------------------------------------------------------
{
  x += rhs.x;
  y += rhs.y;

  return *this;
}

Point& Point::operator -= (const Point& rhs)
// ---------------------------------------------------------------------------
// Subtract x & y locations of rhs from *this.
// ---------------------------------------------------------------------------
{
  x -= rhs.x;
  y -= rhs.y;

  return *this;
}


Point& Point::operator = (const real& val)
// ---------------------------------------------------------------------------
// Set x & y to val;
// ---------------------------------------------------------------------------
{
  x = y = val;

  return *this;
}


Point& Point::operator *= (const real& alpha)
// ---------------------------------------------------------------------------
// Multiply x * y values by alpha;
// ---------------------------------------------------------------------------
{
  x *= alpha;
  y *= alpha;

  return *this;
}


ostream& operator << (ostream& s, const Point& p) { s << p.x <<" "<< p.y; return s; }
istream& operator >> (istream& s, Point& p) { s >> p.x      >> p.y; return s; }


Point operator + (const Point& lhs,
		  const Point& rhs)
// ---------------------------------------------------------------------------
// Add rhs to lhs.
// ---------------------------------------------------------------------------
{
  Point result (lhs);

  return result += rhs;
}


Point operator - (const Point& lhs,
		  const Point& rhs)
// ---------------------------------------------------------------------------
// Subtract rhs from lhs.
// ---------------------------------------------------------------------------
{
  Point result (lhs);

  return result -= rhs;
}


Point operator * (const real&  alpha,
		  const Point& rhs  )
// ---------------------------------------------------------------------------
// Return alpha * rhs;
// ---------------------------------------------------------------------------
{
  Point result (rhs);

  return result *= alpha;
}


Point intersect (const Point& start1, const Point& end1,
		 const Point& start2, const Point& end2)
// ---------------------------------------------------------------------------
// Return point of intersection for two lines defined by start & end points.
// ---------------------------------------------------------------------------
{
  real a1, b1 = 1.0e30, c1;
  real a2, b2 = 1.0e30, c2;
  real det, sinsq;

  if (start1.x == end1.x) {
    a1 = 1.0; b1 = 0.0; c1 = -start1.x;
  }
  if (start2.x == end2.x) {
    a2 = 1.0; b2 = 0.0; c2 = -start2.x;
  }
  
  if (b1 > 0.0) {
    a1 = end1.y   - start1.y;
    b1 = start1.x - end1.x;
    c1 = start1.y * end1.x - end1.y * start1.x;
  }
  
  if (b2 > 0.0) {
    a2 = end2.y   - start2.y;
    b2 = start2.x - end2.x;
    c2 = start2.y * end2.x - end2.y * start2.x;
  }
  
  det   = a1 * b2 - a2 * b1;
  sinsq = sqr (det) / ((sqr(a1) + sqr(b1)) * (sqr(a2) + sqr(b2)));
  
  if (sinsq < EPSSP * EPSSP)	// -- Point at infinity.
    return Point (1.0e30, 1.0e30);
  else
    return Point ((b1 * c2 - b2 * c1) / det, (c1 * a2 - c2 * a1) / det);
}


Point Point::relative (const Point& P     ,
		       const real&  length,
		       const real&  angle ) const
// ---------------------------------------------------------------------------
// Return the Point which lies at length from current Point, at angle
// relative to the line which joins *this to P.  Input angle is assumed
// to lie between 0 & TWOPI (i.e. positive).
// ---------------------------------------------------------------------------
{
  real  dx, dy, theta;
  Point Z;

  dx    = P.x - x;
  dy    = P.y - y;
  theta = atan2 (dy, dx) + angle;
  
  Z.x = x + length * cos (theta);
  Z.y = y + length * sin (theta);

  return Z;
}


int cross (const Point& start1, const Point& end1,
	   const Point& start2, const Point& end2)
// ---------------------------------------------------------------------------
// Return 1 if line *segments* defined by start & end Point cross over.
// See Ref. [3], Ch 24.
// ---------------------------------------------------------------------------
{
  return ((start1.ccw (end1, start2) * start1.ccw (end1, end2)) <= 0)
      && ((start2.ccw (end2, start1) * start2.ccw (end2, end1)) <= 0);
}


int cull (const Point& Pi, const Point& Pk,
	  const Point& Pj, const Point& Pl)
// ---------------------------------------------------------------------------
// Check potential visibility of Pk from Pi.
// Pj & Pl are Points which lie CW & CCW respectively of Pk in Loop.
// Return 1 if Pk is obscured from Pi, otherwise 0.
//
// This code is closely related to back-face culling techniques in hidden
// surface removal.  We employ cross products of vectors, using
//    T1: the turn from Pk to Pi through Pj,
//    T2: the turn from Pk to Pi through Pl,
//    T3: the turn from Pj to Pl through Pk.
// Given the properties of vectors, there cannot be two zero values amongst
// T1, T2, T3, although there can be either 3 zeros or one.  The three-zero
// case has to be tested on distances between points, but the rest can be
// determined according to a table of rules.
//
// Visibility truth table taken on sign of turns: (return the compliment)
//                                                
//    +====+====+====+===+    +====+====+====+===+    +====+====+====+===+
//    | T1 | T2 | T3 |   |    | T1 | T2 | T3 |   |    | T1 | T2 | T3 |   |
//    +----+----+----+---+    +----+----+----+---+    +----+----+----+---+
//    |  + |  0 |  + | . |    |  0 |  0 |  0 | ? |    |  + |  0 |  - | . |
//    |  0 |  + |  + | . |    |  + |  - |  0 | . |    |  0 |  + |  - | 1 |
//    |  - |  0 |  + | . |    |  - |  + |  0 | 1 |    |  - |  0 |  - | 1 |
//    |  0 |  - |  + | . |    |  + |  + |  0 | 1 |    |  0 |  - |  - | . |
//    |  + |  - |  + | . |    |  - |  - |  0 | 1 |    |  + |  - |  - | . |
//    |  - |  + |  + | 1 |    +====+====+====+===+    |  - |  + |  - | 1 |
//    |  + |  + |  + | . |                            |  + |  + |  - | 1 |
//    |  - |  - |  + | . |                            |  - |  - |  - | 1 |
//    +====+====+====+===+                            +====+====+====+===+
// ---------------------------------------------------------------------------
{
  char routine[] = "cull";

  const real EPS = EPSSP;
  const real TOL = 0.01;

  const real dj = Pj.distance (Pk);
  const real dl = Pl.distance (Pk);
  const real di = Pi.distance (Pk);

  if (dj < EPS || dl < EPS || di < EPS) {
//    if (Global::verbose) error (routine, "two points are coincident", WARNING);
//    if (verbose) error (routine, "two points are coincident", WARNING);
    return 0;
  }

  Point pk (0.0, 0.0);
  Point pj = 1.0 / dj * (Pj - Pk);
  Point pl = 1.0 / dl * (Pl - Pk);
  Point pi = 1.0 / di * (Pi - Pk);

  const real T1 = pj.turn (pk, pi);
  const real T2 = pl.turn (pk, pi);
  const real T3 = pk.turn (pj, pl);
  
  if (fabs (T3) < EPS) {
    if (fabs (T1) < TOL || fabs (T2) < TOL)
      return (di < dj && di < dl) ? 0 : 1;
    else if (T1 > 0.0 && T2 < 0.0)
      return 1;
    else 
      return 0;

  } else if (T3 > 0.0) {
    if (T1 < -TOL && T2 > TOL)
      return 0;
    else
      return 1;
    
  } else if (T3 < 0.0) {
    if (fabs (T2) < EPS && T1 > 0.0 ||
	fabs (T1) < EPS && T2 < 0.0 ||
	      T1  > EPS && T2 < EPS )
      return 1;
    else
      return 0;

  } else error (routine, "can't get here!", ERROR);

  return 0;
}


Point unitNormal (const Point& P1,
		  const Point& P2)
// ---------------------------------------------------------------------------
// Return the unit normal vector for the line which joins P1 to P2.
// Vector is oriented so that it points to the right during traverse
// from P1 to P2.  This corresponds to an outward direction during a CCW
// traverse of a convex Loop.
// ---------------------------------------------------------------------------
{
  Point p;
  real  dx, dy, len;

  len = P2.distance (P1);

  if (len < EPSSP)
    p = 1.0e30;
  else {
    dx  = P2.x - P1.x;
    dy  = P2.y - P1.y;

    p.x =  dy / len;
    p.y = -dx / len;
  }

  return p;
}


real Point::dot (const Point& p) const
// ---------------------------------------------------------------------------
// Dot product of points treated as vectors.
// ---------------------------------------------------------------------------
{
  real dp = p.x * x + p.y * y;
  
  return dp;
}
