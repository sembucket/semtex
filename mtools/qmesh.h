///////////////////////////////////////////////////////////////////////////////
// qmesh.h: header file for quadrilateral mesh generator.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#ifndef QMESH_H
#define QMESH_H

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>

#include <sm_options.h>
#include <sm_declare.h>

#include <Array.h>
#include <List.h>

typedef float real;
enum lev {WARNING, ERROR, REMARK};
void message (const char*, const char*, const lev&);

const  int    StrMax = 256;
const  double EPSSP  = 6.0e-7;
const  double EPSDP  = 6.0e-14;
const  double TWOPI  = 6.28318530717958647692;

extern int verbose;
extern int graphics;

template<class T> inline T sqr(T x)      { return x * x;            }
template<class T> inline T sgn(T x)      { return (x < 0) ? -1 : 1; }
template<class T> inline T min(T a, T b) { return (a < b) ?  a : b; }
template<class T> inline T max(T a, T b) { return (a > b) ?  a : b; }

// -- Splitting parameters: "rules of thumb" described in Reference [1].

const real C1     = 0.5;	// -- Angle  weight factor.
const real C2     = 0.3;	// -- Length weight factor.
const real C3     = 0.2;	// -- Error  weight factor.

const int  InsMax = 64;	        // -- Max Nodes to insert on splitting line.


class Point
// ===========================================================================
// 2D points.
// ===========================================================================
{
friend ostream& operator << (ostream&, Point&);
friend istream& operator >> (istream&, Point&);
public:
  Point () : x (0.0), y (0.0) { }
  Point (const real& a, const real& b) { x = a; y = b; }
  Point (const Point& p) { x = p.x; y = p.y; }

  real  magnitude () const { return hypot (x, y); }
  real  distance  (const Point&) const;
  real  angle     (const Point&) const;
  real  angle     (const Point&, const Point&) const;
  real  turn      (const Point&, const Point&) const;

  Point& operator += (const Point&);
  Point& operator -= (const Point&);
  Point& operator  = (const real&);
  Point& operator *= (const real&);

  real x;
  real y;
};
Point operator + (const Point&, const Point&);
Point operator - (const Point&, const Point&);
Point operator * (const real&,  const Point&);
Point intersect  (const Point&, const Point&, const Point&, const Point&);


class Node
// ===========================================================================
// A Node has a Point, an index number and an ideal element length scale.
// ===========================================================================
{
friend ostream& operator << (ostream&, Node&);
public:
  Node (const int& i, const Point& p, const real& s, const int& b = 0) 
    : id (i), loc (p), ideal (s), boundary (b) { }

  const int&   ID       () const { return id;    }
  const Point& pos      () const { return loc;   }
  const real&  prefSize () const { return ideal; }

  void  xadd     (Node*);
  void  setPos   (const Point& p) { if (!boundary){loc.x = p.x; loc.y = p.y;} }
  Point centroid () const;

private:
  int         id;
  Point       loc;
  real        ideal;
  int         boundary;
  List<Node*> contact;
};


class Loop
// ===========================================================================
// A Loop has a list of Nodes and possibly a splitting line and two sub-Loops.
// ===========================================================================
{
friend istream& operator >> (istream&, Loop&);
friend ostream& operator << (ostream&, Loop&);
public:
  Loop ();
  Loop (vector<Node*>&, vector<Node*>&, Node*, Node*, const int&);

  const int& ID () const { return id; }

  void  split   ();
  void  connect ();
  void  smooth  ();
  void  draw    () const;

  void  limits (Point&, Point&)               const;
  int   points (vector<real>&, vector<real>&) const;
  int   line   (vector<real>&, vector<real>&) const;
  
private:
  int                id;
  static int         node_id_max;
  static int         loop_id_max;
  static List<Node*> node_list;

  vector<Node*> nodes;
  vector<Node*> splitline;
  Loop*         left;
  Loop*         right;

  real  lengthScale  () const;
  Point centroid     () const;

  void  visibleNodes (List<Node*>&, const int&)                          const;
  void  bestSplit    (List<Node*>*, int&, int&, int&)                    const;
  void  bestLine     (List<Node*>&, const real&,
		      const int&, int&, int&, real&)                     const;
  real  spaceNodes   (const Node*, const Node*, const real&, const int&) const;
  void  insertNodes  (const Node*, const Node*, const int&);

  void  splitSix     (int&, int&);
};


// -- Routines in graphics.cc:

void initGraphics  ();
void stopGraphics  ();
void eraseGraphics ();
void drawBox  (const Loop*);
void drawLoop (const Loop*);

#endif
