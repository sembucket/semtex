/*****************************************************************************
 * Routines to deal with meshes based on vertices.
 *
 * Here is an example mesh description.  Note that all sections of mesh
 * description are required.
 *
 * mesh 
 * {
 * 9	vertices
 * #id  x       y       z
 * 1	0.0	0.0	0.0
 * 2	2.0	0.0	0.0
 * 3	4.0	0.0	0.0
 * 4	0.0	0.5	0.0
 * 5	2.0	0.5	0.0
 * 6	4.0	0.5	0.0
 * 7	0.0	1.0	0.0
 * 8	2.0	1.0	0.0
 * 9	4.0	1.0	0.0
 *
 * 4	elements
 * #id  n_vert  vertices
 * 1	4	1 2 5 4
 * 2	4	2 3 6 5
 * 3	4	4 5 8 7
 * 4	4	5 6 9 8
 *
 * 4	boundaries
 * #id  BCtag   n_vert  vertices
 * 1	2	3	1 2 3
 * 2	2	3	7 8 9
 * 3	1	3	7 4 1
 * 4	3	3	3 6 9
 *
 * 1	curve
 * #id   kind   elmt   side   param   param
 * 1	 c	1      2      1.0     +
 * }
 *****************************************************************************/

// $Id$

#include "Fem.h"

#ifdef __DECCXX
  #pragma define_template Array1d<int>
  #pragma define_template Array1d<Mesh::Node*>
  #pragma define_template Array1d<Mesh::Side*>
#endif


static const int UNSET = -1;
static const int NIL   =  0;


class Curve {
  // ========================================================================
  // Abstract base class for curved edges.
  // ========================================================================
public:
  virtual  ~Curve () { }
  virtual  void    compute (int, const real*, Point*) const = 0;
  int              ismatch (int e, int s) const {
    return  (e == curveSide -> thisElmt -> ID && s == curveSide -> ID);
  }

protected:
  Mesh::Side*  curveSide;
};


class CircularArc : public Curve {
  // ========================================================================
  // Prototype curved edge class.   To make new kinds, retain the same
  // interface structure and put new instances in Mesh::createCurve below.
  // ========================================================================
public:
  CircularArc (const char*, const Mesh&);
  virtual void  compute (int, const real*, Point*) const;

private:
  Point  centre;
  real   radius;
  real   semiangle;
};


Curve*  Mesh::createCurve (const char* s)
// ---------------------------------------------------------------------------
// Wrapper function to return a base class pointer for each kind of curve.
//
// This routine needs modification each time a new curve type is implemented.
//
// Start of curve string is:
// curveID  curveKind  elementID  sideID
// After this can follow as many parameters as needed to define curve.
// ---------------------------------------------------------------------------
{
  char    routine[] = "Mesh::createCurve";
  char    err[StrMax];

  Curve*  C = 0;
  char    kind;
  int     id;

  istrstream strm (s, StrMax);

  strm >> id >> kind;

  if (strm.bad ()) {
    sprintf (err, "couldn't parse curve ident from: %s", s);
    message (routine, err, ERROR);
  }
  
  switch (kind) {  // -- Add new cases for new kinds.
  case 'c':
    C = new CircularArc (s, *this);
    break;
  default:
    sprintf (err, "can't parse a known curve type from: %s", s);
    message (routine, err, ERROR);
    break;
  }

  return C;
}


CircularArc::CircularArc (const char* s, const Mesh& m)
// ---------------------------------------------------------------------------
// Constructor for CircularArc.  Input string is to contain parameters for
// arc as follows:
//
// id    kind   elmt   side   radius  convexity
// 1     c      1      2      1.0     +
//
// Convexity == "+" ==> arc increases area enclosed by element (cf line),
// Convexity == "-" ==> arc decreases area enclosed by element.
// ---------------------------------------------------------------------------
{
  char   routine[] = "CircularArc::CircularArc";
  char   err[StrMax];
  int    verbose  = option ("VERBOSE");

  int    cid, elmt, side, convexity;

  istrstream strm (s, StrMax);

  if (verbose == 2) message (routine, s, REMARK);

  strm >> cid >> err[0] >> elmt >> side >> radius >> err[0];

  if (strm.bad ()) {
    sprintf (err, "couldn't parse circular arc parameters from: %s", s);
    message (routine, err, ERROR);
  } else if (radius < 0.0) {
    sprintf (err, "curve %1d: negative radius: %f", cid, radius);
    message (routine, err, ERROR);
  } else if (err[0] == '+') {
    convexity = 1;
  } else if (err[0] == '-') {
    convexity = -1;
  } else {
    message (routine, "element convexity not \"+\" or \"-\"", ERROR);
  }

  // -- Traverse mesh sides to find this one.

  curveSide = m.getSide (elmt, side);

  if (!curveSide) {
    sprintf (err, "can't find element %1d side %1d in mesh", elmt, side);
    message (routine, err, ERROR);
  }

  // -- Locate centre on bisector of nodes according to convexity.

  Point P1 = curveSide -> startNode -> location ();
  Point P2 = curveSide -> endNode   -> location ();
  Point unitNormal, link, midpoint, centroid = {0.0, 0.0, 0.0};
  real  dx, dy, l, sign = 0.0;

  midpoint.x   = 0.5 * (P2.x + P1.x);
  midpoint.y   = 0.5 * (P2.y + P1.y);
  dx           =        P2.x - P1.x;
  dy           =        P2.y - P1.y;
  l            = hypot (dx, dy);
  unitNormal.x = -dy / l;
  unitNormal.y =  dx / l;

  if (2.0 * radius < l) {
    sprintf (err, "curve %1d:\ncircle, radius %f, can't span nodes %1d & %1d",
	     cid, radius, 
	     curveSide -> startNode -> ID, curveSide -> endNode -> ID);
    message (routine, err, ERROR);
  } else {
    semiangle = asin (0.5*l / radius);
  }

  centroid = curveSide -> thisElmt -> centroid ();
  
  link.x = centroid.x - midpoint.x;
  link.y = centroid.y - midpoint.y;

  // -- Sign +1 if centre lies in direction of centroid from side midpoint.

  sign = link.x * unitNormal.x + link.y * unitNormal.y;
  sign = convexity * sign / fabs (sign);

  centre.x = midpoint.x + sign * cos (semiangle) * radius * unitNormal.x;
  centre.y = midpoint.y + sign * cos (semiangle) * radius * unitNormal.y;
}


void  CircularArc::compute (int np, const real* spacing, Point* knot) const
// ---------------------------------------------------------------------------
// Distribute np knots along arc according to spacing on -1, 1.
// ---------------------------------------------------------------------------
{
  Point  P1 = curveSide -> startNode -> location ();
  Point  P2 = curveSide -> endNode   -> location ();
  real   theta1, theta2, dtheta, phi;
  int    nm = np - 1;

  theta1 = atan2 (P1.y - centre.y, P1.x - centre.x);
  theta2 = atan2 (P2.y - centre.y, P2.x - centre.x);
  dtheta = theta2 - theta1;

  if (abs (dtheta) > 2.0*semiangle + EPSSP)
    dtheta += (dtheta < 0.0) ? TWOPI : -TWOPI;

  knot[0].x  = P1.x;  knot[0].y  = P1.y;
  knot[nm].x = P2.x;  knot[nm].y = P2.y;

  for (int i(1); i < nm; i++) {
    phi = theta1 + dtheta * 0.5 * (spacing[i] + 1.0);
    knot[i].x = centre.x + radius * cos (phi);
    knot[i].y = centre.y + radius * sin (phi);
  }
}





class NodeReader {
public:
  NodeReader (Mesh& m , istream& strm) : mesh(m), in(strm) { }
  int          getElmtID ();
  int          getSize();
  Mesh::Node*  getNode();
private:
  Mesh&     mesh;
  istream&  in;
};


int  NodeReader::getElmtID ()
// ---------------------------------------------------------------------------
// Return Elmt ID tag
// ---------------------------------------------------------------------------
{
  int id;
  in >> id;

  return id;
}



int  NodeReader::getSize ()
// ---------------------------------------------------------------------------
// Return the number of nodes in an element.
// ---------------------------------------------------------------------------
{
  int size;
  in >> size;

  return size;
}


Mesh::Node*  NodeReader::getNode ()
// ---------------------------------------------------------------------------
// Read an element node number & return a pointer to that node.
// ---------------------------------------------------------------------------
{
  char  routine[] = "NodeReader::getNode";

  int   nodeNum;
  in >> nodeNum;

  int          found(0);
  Mesh::Node*  N;

  for (ListIterator<Mesh::Node*> i(mesh.nodeTable); !found && i.more();
       i.next())  {
    N = i.current ();
    found = nodeNum == N -> ID;
  }
  
  if (found) {
    return N;
  } else {
    char s[StrMax];
    ostrstream ost (s, StrMax);
    ost << "Node number " << nodeNum << " not defined";
    message (routine, s, ERROR);
  }
  
  return 0;
}


void Mesh::assemble ()
// ---------------------------------------------------------------------------
// Traverse Elmt list and fill in edge-based connectivity information.
// ---------------------------------------------------------------------------
{
  // -- First, build element sides.

  for (ElmtsOfMesh i(*this); i.more(); i.next()) {
    Elmt& E        = i.current();
    int   numNodes = E.nNodes();
    int   j        = 0;

    for (NodesOfElmt k(E); k.more(); k.next(), j++) {
      Node& n1 = k.current();
      Node& n2 = k.ccwNeighbour();
      Side* S  = new Side (&n1, &n2);
      S -> setLength (0);
      S -> ID  = j + 1;
      S -> gID = UNSET;
      S -> thisElmt = &E;
      E.insertSide (S, j);
    }
  }

  // -- Now traverse list of Elmts and build side--side connections.

  for (i.reset(); i.more(); i.next()) {
    Elmt& E = i.current();

    for (SidesOfElmt k(E); k.more(); k.next()) {
      Side& thisSide = k.current();
      thisSide.mateElmt = NIL;
      thisSide.mateSide = NIL;
      int   found       = 0;

      for (ElmtsOfMesh m(*this); !found && m.more(); m.next()) {
	Elmt& M = m.current();
	
	for (SidesOfElmt j(M); !found && j.more(); j.next()) {
	  Side& mateSide = j.current();

	  found = (thisSide.startNode == mateSide.endNode &&
		   mateSide.startNode == thisSide.endNode );

	  if (found) {
	    thisSide.mateElmt = &M;
	    thisSide.mateSide = &mateSide;
	  }
	}
      }
    }
  }
}


void   Mesh::installTag (istream& strm, int tag, int nvert)
// ---------------------------------------------------------------------------
// In strm there should be nvert vertex indices.  Install boundary condition
// tag in each element side that is indicated by vertex numbers.
// ---------------------------------------------------------------------------
{
  char routine[] = "Mesh::installTag";
  char s[StrMax];

  if (nvert < 2) {
    sprintf (s, "BC tag %1d, number of vertices < 2: %1d", tag, nvert);
    message (routine, s, ERROR);
  }

  int i1, i2;
  strm >> i1;

  for (int i = 1; i < nvert; i++) {
    strm >> i2;
    for (ListIterator<Elmt*> e(elmtList); e.more(); e.next()) {
      Elmt& E = *e.current();
      for (SidesOfElmt s(E); s.more(); s.next()) {
	Side& S = s.current();
	if ((S.startNode -> ID == i1 && S.endNode -> ID == i2) ||
	    (S.startNode -> ID == i2 && S.endNode -> ID == i1)  ) {
	  if (S.mateElmt) {
	    S.mateSide -> mateElmt = 0;
	    S.mateSide -> BConTag    = tag;
	  }
	  S.mateElmt = 0;
	  S.BConTag    = tag;
	}
      }
    }
    i1 = i2;
  }
}


istream& operator >> (istream& strm, Mesh::Node& n)
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  strm >> n.ID >> n.where.x >> n.where.y >> n.where.z;
  n.gID = UNSET;

  return strm;
}
  

void  operator >> (NodeReader& reader, Mesh::Elmt& e)
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  e.ID = reader.getElmtID();

  int  nNodesInElement = reader.getSize();

  e.nodeArray.setSize(nNodesInElement);
  e.sideArray.setSize(nNodesInElement);
  for (int i(0); i < nNodesInElement; i++) e.nodeArray[i] = reader.getNode ();
}


void  operator << (ostream& strm, Mesh::Elmt& e)
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  int    nNodesInElement = e.nNodes ();
  strm << "Element: " << " " << e.ID << ", Vertices:";

  for (Mesh::NodesOfElmt i(e); i.more(); i.next()) {
    strm << " " << i.current().ID;
  }
  
  strm << ", Mating:";

  for (Mesh::SidesOfElmt s(e); s.more(); s.next()) {
    Mesh::Side& S = s.current ();
    strm << " " << S.ID << "->";
      if (S.mateElmt == NIL)
	strm << "bt" << S.BConTag;
      else
	strm << S.mateElmt -> ID << "." << S.mateSide -> ID;
  }

  cout << endl;
}


istream&  operator >> (istream& strm, Mesh& m)
// ---------------------------------------------------------------------------
// Read in all the information that declares a mesh.  Set up edge--edge
// connectivity but don't do any global node numbering.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Mesh::operator >>";
  char  s[StrMax],  err[StrMax];
  int   verbose   = option ("VERBOSE");

  // -- Look for number of mesh vertices, then input table.

  int          i, nVert;
  Mesh::Node*  N;

  strm >> nVert;
  strm.getline (s, StrMax);
  upperCase    (s);
  if (strstr (s, "VERT")) {
    if (nVert < 3) {
      sprintf (err, "need at least 3 vertices, got: %1d", nVert);
      message (routine, err, ERROR);
    }
  } else {
    sprintf (err, "expected number of vertices, got: %s", s);
    message (routine, err, ERROR);
  }

  for (i = 0; i < nVert; i++) {
    N = new Mesh::Node;
    strm >> *N;
    m.nodeTable.add (N);
  }

  // -- Look for number of elements, then input list.

  int          nElmt;
  Mesh::Elmt*  E;
  NodeReader   reader (m, strm);

  strm >> nElmt;
  strm.getline (s, StrMax);
  upperCase    (s);
  if (strstr (s, "ELEMENT")) {
    if (nElmt < 1) {
      sprintf (err, "need at least one element, got: %1d", nElmt);
      message (routine, err, ERROR);
    }
  } else {
    sprintf (err, "expected number of elements, got: %s", s);
    message (routine, err, ERROR);
  }

  for (i = 0; i < nElmt; i++) {
    E = new Mesh::Elmt;
    reader >> *E;
    m.elmtList.add (E);
  }

  // -- Fill in element--element connectivity information.

  m.assemble ();

  // -- Set boundary tag information.

  int nBound, iTag;

  strm >> nBound;
  strm.getline (s, StrMax);
  upperCase    (s);
  if (strstr (s, "BOUNDAR")) {
    if (nBound < 1) {
      sprintf (err, "need at least one declared boundary, got: %1d", nBound);
      message (routine, err, ERROR);
    }
  } else {
    sprintf (err, "expected number of boundaries, got: %s", s);
    message (routine, err, ERROR);
  }

  for (i = 0; i < nBound; i++) {
    strm >> iTag >> iTag >> nVert;
    m.installTag (strm, iTag, nVert);
  }

  // -- Now check that every element side is accounted for.

  m.checkAssembly ();

  // -- Look for curved side information and load it if required.
  //    Curve input is loaded line-by-line from the file.

  int     nCurve;
  Curve*  C;

  strm >> nCurve;
  strm.getline (s, StrMax);
  upperCase    (s);
  if (strstr (s, "CURVE")) {
    if (nCurve < 0) {
      sprintf (err, "can't have a negative number of curves: %1d", nBound);
      message (routine, err, ERROR);
    }
  } else {
    sprintf (err, "expected number of curves, got: %s", s);
    message (routine, err, ERROR);
  }

  for (i = 0; i < nCurve; i++) {
    strm.getline (s, StrMax);
    C = m.createCurve  (s);
    m.curveList.add (C);
  }

  return strm;
}


void  Mesh::checkAssembly ()
// ---------------------------------------------------------------------------
// All element sides have to either mate an adjoining element or fall on
// a boundary.  Check it out.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Mesh::checkAssembly";
  char  err[StrMax];

  for (ElmtsOfMesh e(*this); e.more(); e.next()) {
    for (SidesOfElmt s(e.current()); s.more(); s.next()) {
      Side& S = s.current();
      if (S.mateSide == NIL) {
	sprintf (err, "Element %1d Side %1d not set", S.thisElmt -> ID, S.ID);
	message (routine, err, ERROR);
      }
    }
  }
}


void operator << (ostream& strm, Mesh::Side& s)
// ---------------------------------------------------------------------------
// Print Elmt connectivity information in SEM-compatible form.
// ---------------------------------------------------------------------------
{
  strm << setw(5) << s.thisElmt -> ID;
  strm << setw(5) << s.ID;
  strm << setw(5) << 1;
  strm << setw(5) << s.startNode -> gID;
  strm << endl;

  int ni = s.gID.numElts();

  if (ni)
    for (int i(0); i < ni; i++) {
      strm << setw(5) << s.thisElmt -> ID;
      strm << setw(5) << s.ID;
      strm << setw(5) << i + 2;
      strm << setw(5) << s.gID[i];
      strm << endl;
    }

  strm << setw(5) << s.thisElmt -> ID;
  strm << setw(5) << s.ID;
  strm << setw(5) << ni + 2;
  strm << setw(5) << s.endNode -> gID;
  strm << endl;
}


void Mesh::printAssembly (Mesh& m)
// ---------------------------------------------------------------------------
// Print edge--edge connectivity information.
// ---------------------------------------------------------------------------
{
  cout << m.elmtList.length () << " elements" << endl;

  for (ElmtsOfMesh  i(m); i.more(); i.next()) cout << i.current ();
}


void Mesh::printConnections (Mesh& m)
// ---------------------------------------------------------------------------
// Print knot connectivity information (global node numbers).
// ---------------------------------------------------------------------------
{
  cout << "# " << m.elmtList.length() << " NEL" << endl;
  cout << "# 0 BANDWIDTH" << endl;

  for (ElmtsOfMesh e(m); e.more(); e.next())
    for (SidesOfElmt s(e.current()); s.more(); s.next())
      cout << s.current();
}


void Mesh::Side::connect (int ni, int& gid)
// ---------------------------------------------------------------------------
// Fill in connectivity for this element side, updating global number gid.
// ---------------------------------------------------------------------------
{
  if (startNode -> gID == UNSET) startNode -> gID = gid++;

  if (ni) {  // -- Do side-internal gIDs.
    int i, k;

    if (mateElmt != NIL) {
      Side& otherSide = *mateSide;

      if (otherSide.gID[0] == UNSET)
	for (i = 0; i < ni; i++) gID[i] = gid++;
      else
	for (i = 0, k = ni - 1; i < ni; i++, k--) gID[i] = otherSide.gID[k];

    } else
      for (i = 0; i < ni; i++)   gID[i] = gid++;
  }

  if (endNode -> gID == UNSET)   endNode -> gID = gid++;
}


int  Mesh::connectSC (int np)
// ---------------------------------------------------------------------------
// Generate connectivity (i.e. global knot numbers) for a mesh with
// np knot points (i.e. Lagrange knots) along each element side,
// ignoring internal points (i.e. connectivity for static-condensation form).
// Return maximum global node number.  NB: global numbers generated here
// start at 0.
//
// NB: np >= 2.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Mesh::connectSC";

  if (np < 2) message (routine, "need at least 2 knots", ERROR);

  // -- Create element-side based gID storage, if required, & initialize gIDs.

  int ni   = np - 2;

  for (ElmtsOfMesh e(*this); e.more(); e.next()) {
    for (SidesOfElmt s(e.current()); s.more(); s.next()) {
      Side& S = s.current();
      S.gID.setSize(ni);      
      S.startNode -> gID = UNSET;
      S.endNode   -> gID = UNSET;
      if (ni)      S.gID = UNSET;
    }
  }

  // -- Generate connectivity information.

  int gid = 0;

  for (e.reset(); e.more(); e.next())
    for (SidesOfElmt s(e.current()); s.more(); s.next())
      s.current().connect (ni, gid);

  return gid;
}


Mesh::Side*  Mesh::getSide (int elmt, int side) const
// ---------------------------------------------------------------------------
// Return pointer to element Side indicated by input idents.
// ---------------------------------------------------------------------------
{
  Side*  CS = 0;

  for (ElmtsOfMesh i(*this); !CS && i.more(); i.next())
    for (SidesOfElmt k(i.current()); !CS && k.more(); k.next()) {
      Side& S = k.current();
      if (elmt == S.thisElmt -> ID && side == S.ID) CS = &S;
    }

  return CS;
}


void  Mesh::Elmt::getGid (int* b)
// ---------------------------------------------------------------------------
// Fill b with Elmt-boundary global node numbers.
// ---------------------------------------------------------------------------
{
  int start, nint;

  for (SidesOfElmt s(*this); s.more(); s.next()) {
    Side& S  = s.current ();
    nint     = S.gID.numElts();
    start    = (S.ID - 1) * (nint + 1);
    b[start] = S.startNode -> gID;
    for (register int i(0); i < nint; i++) b[start + i + 1] = S.gID[i];
  }
}


Point  Mesh::Elmt::centroid () const
// ---------------------------------------------------------------------------
// Return point that is centroid of element Node points.
// ---------------------------------------------------------------------------
{
  Point  C = {0.0, 0.0, 0.0};

  for (Mesh::NodesOfElmt n(*this); n.more(); n.next()) {
    Point P = n.current().location();
    C.x += P.x;
    C.y += P.y;
  }

  C.x /= nodeArray.numElts ();
  C.y /= nodeArray.numElts ();

  return C;
}


void   Mesh::meshSide (int          np     ,
		       int          elmt   ,
		       int          side   ,
		       const real*  spacing,
		       Point*       knot   ) const
// ---------------------------------------------------------------------------
// If a curved side can be identified for the nominated element and side,
// compute the points using appropriate routine.  Otherwise compute points
// along a straight side.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Mesh::meshSide";

  if (np < 2) message (routine, "must have at least two points", ERROR);

  for (ListIterator<Curve*> c(curveList); c.more(); c.next()) {
    Curve& C = *c.current();
    if (C.ismatch (elmt, side)) {
      C.compute (np, spacing, knot);
      return;
    }
  }

  // -- Default: straight line.

  Point  P1, P2;

  for (ElmtsOfMesh e(*this); e.more(); e.next()) {
    for (Mesh::SidesOfElmt s(e.current()); s.more(); s.next()) {
      Side& S = s.current();
      if (S.thisElmt -> ID == elmt && S.ID == side) {
	P1 = S.startNode -> location ();
	P2 = S.endNode   -> location ();
      }
    }
  }

  real  dx = P2.x - P1.x;
  real  dy = P2.y - P1.y;

  for (int i(0); i < np; i++) {
    knot [i].x = P1.x + dx * 0.5 * (spacing[i] + 1.0);
    knot [i].y = P1.y + dy * 0.5 * (spacing[i] + 1.0);
  }
  return;
}
