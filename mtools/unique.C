///////////////////////////////////////////////////////////////////////////////
// unique.C: fix a list of nodes and elements produced by mesh
// digitisation (e.g. by digmesh, Murray's code), where the vertices
// may lack uniqueness.
//
// Usage: unique [-h] [-t tol] [file]
//
// Where tol is a positional radius tolerance (default value: DefTol)
// and file is an ASCII file of form
//
//  47     Nodes
//  1.9120   0.4040   0.0000
//  1.9790   0.0337   0.0000
//  ...
//  ...
//
//  196     Elements
//     1     4     1     2     3     4
//     2     4     5     1     4     6
//     3     4     7     5     6     8
//     4     4     9     7     8    10
//     ...
//     ...
//
// By default we deal here only with quad meshes, so the second index
// in the element list above is always 4.
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>		// C headers.
#include <limits.h>
#include <stdio.h>
#include <math.h>

#include <iostream>		// C++ headers.
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
using namespace std;

typedef  double real;
enum lev {WARNING, ERROR, REMARK};

const char* prog   = "unique";
const int   StrMax = 256;
const real  DefTol = 1.0e-3;


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
// A Node has an ID and a Point location.
// ===========================================================================
{
friend ostream& operator << (ostream&, Node&);
public:
  Node (const int& i, const Point& p) : _id (i), _loc (p) { }
  const int&   ID  () const { return _id;  }
  const Point& pos () const { return _loc; }
private:
  int   _id;
  Point _loc;
};


ostream& operator << (ostream& s,
		      Node&    n)
// ---------------------------------------------------------------------------
// Print up node information.
// ---------------------------------------------------------------------------
{
  s << setw  (5) << n._id << setw (10) << n._loc.x << setw (10) << n.loc.y;
  return s;
}


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// We have a table of Node*'s (the length of which is given in the
// input file), and each element has an array of indices into this
// table, but Nodes pointed to by various elements of the table are
// not necessarily unique.  Each time a new node definition is read in
// from the input file, we traverse a list of Node*'s and check if we
// already have this stored (to within positional tolerance): if not,
// a new Node is created.  Finally we print up the list of unique
// Nodal positions, and the Element listing using the IDs in the
// Node*'s.  I.e. the table of Node*'s gives a level of indirection
// through which we are able to give unique re-indexation.
// ---------------------------------------------------------------------------
{
  const char* usage = "Usage: %s [-h] [-t tol] [file]";
  char        buf[StrMax];
  ifstream    infile;
  real        PosTol = DefTol;
  real        x, y;
  int         num, i, id1, id2, id3, id4, gid = 0, found;
  Point       P;
  Node*       N;

  list<Node*>                 uniq;
  list<Node*>::const_iterator u;
  vector<Node*>               nodetab;
  vector<int*>*               elmttab;

  // -- Parse command line.
  
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 't':
      if (*++argv[0]) PosTol = atof (*argv);
      else { --argc;  PosTol = atof (*++argv); }
      break;
    default:
      sprintf (buf, usage, prog);
      cerr << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if   (argc == 1) infile.open   (*argv, ios::in);
  else             infile.attach (0);

  if (!infile) {
    sprintf (buf, "unable to open file: %s", *argv);
    error   (prog, buf, ERROR);
  }
  
  // -- Read Node locations, generate list of unique Node*'s.
  
  infile >> num >> buf;
  if (strstr (buf, "Nodes"))
    message (prog, "Failed reading number of Nodes from file", ERROR);
  nodetab.resize (num);
  
  // -- First Node *is* unique!

  file >> P.x >> P.y; file.getline (buf, StrMax);
  unique.push_front (N = new Node (gid++, P));
  nodetab[0] = N;

  // -- Continue to end of Nodes.
  
  for (i = 1; i < num; i++) {
    file >> P.x >> P.y; file.getline (buf, StrMax);
    for (found = 0, u = uniq.begin(); !found && u != uniq.end(); u++)
      found = u -> distance (P) < PosTol;
    if  (found) nodetab[i] = u;
    else        unique.push_front (N = nodetab[i] = new Node (gid++, P));
  }

  // -- Now read in Elements.

  // -- Print up revised table of unique Nodes.

  // -- Print up revised table of Element Node indices.
    
  return EXIT_SUCCESS;
}
