///////////////////////////////////////////////////////////////////////////////
// Qmesh.h: header file for quadrilateral mesh generator.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#ifndef QMESH_H
#define QMESH_H

#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include <Array.h>
#include <List.h>

typedef  double real;
enum lev {WARNING, ERROR, REMARK};

const  int    StrMax = 256;
const  double EPSSP  = 6.0e-7;
const  double EPSDP  = 6.0e-14;
const  double TWOPI  = 6.28318530717958647692;

template<class T> inline T sqr(T x)      { return x * x;            }
template<class T> inline T sgn(T x)      { return (x < 0) ? -1 : 1; }
template<class T> inline T min(T a, T b) { return (a < b) ?  a : b; }
template<class T> inline T max(T a, T b) { return (a > b) ?  a : b; }

const int InsMax = 512;	// -- Max Nodes to insert on splitting line.


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
  real  dot       (const Point&) const;
  real  angle     (const Point&) const;
  real  angle     (const Point&, const Point&) const;
  real  turn      (const Point&, const Point&) const;
  int   ccw       (const Point&, const Point&) const;
  Point relative  (const Point&, const real&, const real&) const;

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
int   cross      (const Point&, const Point&, const Point&, const Point&);
int   cull       (const Point&, const Point&, const Point&, const Point&);
Point unitNormal (const Point&, const Point&);


class Node
// ===========================================================================
// A Node has a Point, an ID number and an ideal element length scale.
//
// In addition it knows if it lies on a domain boundary (hence, can't be
// shifted during mesh smoothing), is a boundary node flagged for offset,
// or is an interior node.
//
// A List of contacting Nodes is computed after all Nodes are created,
// which is used during Laplacian smoothing.
// ===========================================================================
{
friend ostream& operator << (ostream&, Node&);
public:
  enum  nodekind {
    INTERIOR,			/* -- Numeric code 0            */
    INTERIOR_FIXED,		/* -- Numeric code 1            */
    LOOP_OFFSET_FIXED,		/* -- Numeric code 2            */
    LOOP_OFFSET_MOBILE,		/* -- Numeric code 3            */
    LOOP_BOUNDARY_FIXED,	/* -- Numeric code 4            */
    LOOP_BOUNDARY_MOBILE,	/* -- Numeric code 5            */
    DOMAIN_OFFSET_FIXED,	/* -- Numeric code 6            */
    DOMAIN_OFFSET_MOBILE,	/* -- Numeric code 7 (not used) */
    DOMAIN_BOUNDARY_FIXED,	/* -- Numeric code 8            */
    DOMAIN_BOUNDARY_MOBILE	/* -- Numeric code 9 (not used) */
  };

  Node (const int&           i,
	const Point&         p,
	const real&          s,
	const Node::nodekind b = INTERIOR)
    : id (i), loc (p), ideal (s), kind (b) { }

  const int&   ID       () const { return id;    }
  const Point& pos      () const { return loc;   }
  const real&  prefSize () const { return ideal; }

  void  renumber (const int i) { id = i; }
  void  xadd     (Node*);
  void  sever    () { contact.clear(); }
  void  setPos   (const Point&);
  void  setKind  (const nodekind& k) { kind  = k; }
  void  setID    (const int&      i) { id    = i; }
  void  setSize  (const real&     s) { ideal = s; }
  Point centroid () const;
  int   interior () const { return
			      kind == INTERIOR               ||
			      kind == INTERIOR_FIXED;
  }
  int   offset   () const { return
			      kind == LOOP_OFFSET_FIXED      ||
			      kind == LOOP_OFFSET_MOBILE     ||
			      kind == DOMAIN_OFFSET_FIXED    ||
			      kind == DOMAIN_OFFSET_MOBILE;
  }
  int   mobile   () const { return
			      kind == INTERIOR               ||
			      kind == LOOP_BOUNDARY_MOBILE   ||
			      kind == LOOP_OFFSET_MOBILE     ||
			      kind == DOMAIN_BOUNDARY_MOBILE ||
			      kind == DOMAIN_OFFSET_MOBILE;
  }
  char  classify () const {
    if      (kind == INTERIOR_FIXED      ||
	     kind == LOOP_OFFSET_FIXED   ||
	     kind == LOOP_BOUNDARY_FIXED  ) return 'F';
    else if (kind == DOMAIN_OFFSET_FIXED ||
	     kind == DOMAIN_BOUNDARY_FIXED) return 'U';
    return 'M';
  }
  int   adjncy   () const { return contact.length(); }

private:
  int         id;
  Point       loc;
  real        ideal;
  nodekind    kind;
  List<Node*> contact;
};


class Quad
// ===========================================================================
// A Quad is a vector of four Node*s.
// ===========================================================================
{
public:
  Node* vertex[4];
};


class Global
// ===========================================================================
// Global maintains the global list of all allocated Nodes for the problem,
// and other geometric data.
// ===========================================================================
{
public:
  static Node*       exist (const Node*);
  static int         nodeIdMax;
  static int         loopIdMax;
  static List<Node*> nodeList;
  static List<Quad*> quadList;

  static real        limits (Point&, Point&);
  static real        lengthScale () { return gblSize; }

  static real        refCoeff;	// -- Refinement coefficient per Ref. [1].
  static const real  C1 = 0.5;	// -- Angle  weight factor.
  static const real  C2 = 0.3;	// -- Length weight factor.
  static const real  C3 = 0.2;	// -- Error  weight factor.

  static int         verbose;

private:
  static real gblSize;
};


class Loop
// ===========================================================================
// A Loop has a list of Nodes and possibly a splitting line and two sub-Loops.
// ===========================================================================
{
friend istream& operator >> (istream&, Loop&);
friend ostream& operator << (ostream&, Loop&);
public:
  Loop (const int&);
  Loop (vector<Node*>&, vector<Node*>&, const int&, const int&, const int&);

  const int&  ID   () const { return id; }
  real        area () const;

  void  offset  ();
  void  split   ();
  void  connect ();
  void  smooth  ();
 
  void  limits   (Point&, Point&)                 const;
  int   points   (vector<float>&, vector<float>&) const;
  int   line     (vector<float>&, vector<float>&) const;
  void  quads    (List<Quad*>&)                   const;
  void  drawQuad (const int& = 0)                 const;

  void  nodeLabel (const int&, char*)             const;
  void  printMesh (ostream&)                      const;

private:
  int                id;

  vector<Node*> nodes;
  vector<Node*> splitline;
  Loop*         left;
  Loop*         right;

  Node* exist        (const Node*, List<Node*>&)                         const;
  real  lengthScale  ()                                                  const;
  Point centroid     ()                                                  const;

  void  visibleNodes (List<Node*>&, const int&)                          const;
  void  bestSplit    (List<Node*>*, int&, int&, int&)                    const;
  void  bestLine     (List<Node*>&, const real&,
		      const int&, int&, int&, real&)                     const;
  real  spaceNodes   (const Node*, const Node*, const real&, const int&) const;
  void  insertNodes  (const Node*, const Node*, const int&);

  void  splitSix     (int&, int&);

  void  printQuad    (ostream&, int&)                                    const;
  void  countElements(int&)                                              const;
};


// -- Routines in graphics.cc:

extern int graphics;

void initGraphics  (const char*);
void stopGraphics  ();
void eraseGraphics ();
void drawBox       ();
void drawBox       (const Loop*);
void drawLoop      (const Loop*, const int& = 0);
void drawMesh      (List<Quad*>&, const int = 0);
void hardCopy      (List<Quad*>&);
void qpause        ();
void message       (const char*, const char*, const lev&);

// -- Routines in misc.cc:

extern int verbose;

void      error     (const char*, const char*, const lev&);
istream&  seekBlock (istream&, const char*);
istream&  endBlock  (istream&);
char*     upperCase (char*);

#endif
