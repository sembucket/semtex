#ifndef MESH_H
#define MESH_H
///////////////////////////////////////////////////////////////////////////////
// mesh: header file for Mesh and related classes.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>

#include <cfemdef.h>
#include <feml.h>

class Curve;


class Mesh
// ===========================================================================
// Mesh class is used as a storage for inter-element connectivities,
// BC tags, element-corner vertex locations and curved boundaries.
// It plays no direct part in the eventual solution of the discrete
// equations.  Presently, meshes are 2D only.
// 
// Mesh provides seven externally-visible facilities: 
// (1) Return number of elements in mesh,
// (2) Generation of element-boundary boundary-to-global mapping vector.
// (3) Generation of element-boundary Essential BC mask vector.
// (4) Generation of element mesh knot points.
// (5) Print up NEKTON-style .rea information.
// (6) Return x--y extent of mesh, based on traverse of element vertices.
// (7) Return true if an element touches the z axis (cylindrical coords).
//
// (2--4) can be independently computed for arbitrary polynomial orders.
// ===========================================================================
{
public:
  class Node;
  class Elmt;
  class Side;

  enum IDstatus {UNSET = -1};

  Mesh (FEML*, const bool = true);

  integer nEl         () const { return _elmtTable.size(); }
  integer buildMap    (const integer, integer*);
  void    buildMask   (const integer, const char, integer*);
  void    meshElmt    (const integer, const integer, const real*, const real*,
		       real*, real*) const;
  void printNek       () const;
  void extent         (Point&, Point&) const;

  bool isAxial (const integer id) const {
    return _elmtTable[id] -> side [0] -> axial && 
      (static_cast<integer>(Femlib::value("CYLINDRICAL")) == 1); 
  };

  static void showGlobalID (Mesh&);
  static void showAssembly (Mesh&);

  class Node {
  public:
    integer ID;
    integer gID;
    Point   loc;
    Node*   periodic;
  };

  class Elmt {
    friend class Mesh;
  private:
    vector<Node*> node;
    vector<Side*> side;
  public:
    integer ID;
    integer nNodes   () const    { return node.size(); }
    Node*   ccwNode  (integer i) { return node[(i         +1) % nNodes()]; }
    Node*   cwNode   (integer i) { return node[(i+nNodes()-1) % nNodes()]; }
    Point   centroid () const;
  };
  
  class Side {
    friend class Mesh;
  private:
    Side (integer id, Node* n1, Node* n2) :
      ID(id), startNode(n1), endNode(n2), mateElmt(0) {
      gID.resize (0); mateSide = 0;
      axial = (n1 -> loc.y == 0.0) && (n2 -> loc.y == 0.0);
    }
    Side () {}
  public:
    integer         ID;
    bool            axial;
    vector<integer> gID;
    Node*           startNode;
    Node*           endNode;
    Elmt*           thisElmt;
    Elmt*           mateElmt;	// -- Doubles as a flag for union:
    union { Side* mateSide; char group; };
    void connect (const integer, integer&);
  };

private:
  FEML&          _feml;
  vector<Node*>  _nodeTable;
  vector<Elmt*>  _elmtTable;
  vector<Curve*> _curveTable;

  void surfaces      ();
  void curves        ();
  void assemble      ();
  void checkAssembly ();
  void checkAxial    ();
  void chooseNode    (Node*, Node*);
  void fixPeriodic   ();

  void describeGrp (char, char*)          const;
  void describeBC  (char, char, char*)    const;
  bool matchBC     (const char, const char, const char);

  void meshSide (const integer, const integer, const integer,
		 const real*, Point*)                          const;
};


class Curve
// ===========================================================================
// Base class for curved edges.
// ===========================================================================
{
public:
  virtual ~Curve () { }

  virtual void compute  (const integer, const real*, Point*) const = 0;
  virtual void printNek ()                                   const = 0;

  bool ismatch (const integer e, const integer s) const {
    return (e == curveSide -> thisElmt -> ID && s == curveSide -> ID);
  }

protected:
  Mesh::Side* curveSide;
};


class CircularArc : public Curve
// ===========================================================================
// Prototype curved edge class.   To make new kinds, retain the same
// interface structure and put new instances in Mesh::curves.
// ===========================================================================
{
public:
  CircularArc (const integer, Mesh::Side*, const real);
  virtual void compute  (const integer, const real*, Point*) const;
  virtual void printNek () const;

private:
  Point   centre;
  real    radius;
  real    semiangle;
  integer convexity;
};


class spline2D
// ===========================================================================
// Definition of spline2D struct used by Spline class.
// ===========================================================================
{
public:
  const char*  name  ;		// name of file containing knot points
  vector<real> x     ;		// knot x-coordinates
  vector<real> y     ;		// knot y-coordinates
  vector<real> sx    ;		// spline x-coefficients
  vector<real> sy    ;		// spline y-coefficients
  vector<real> arclen;	        // arclength along curve at knot locations
  integer      pos   ;		// index of last confirmed position
};


class Spline : public Curve
// ===========================================================================
// A 2D splined curve defined in a named file.
// ===========================================================================
{
public:
  Spline (const integer, Mesh::Side*, const char*);
  virtual void compute  (const integer, const real*, Point*) const;
  virtual void printNek () const;

private:
  char*     _name;
  spline2D* _geom;
  real      _startarc;
  real      _endarc;

  spline2D* getGeom  (const char*);
  real      getAngl  (const real&);
  integer   closest  (const Point&);
  real      arcCoord ();
};


#endif
