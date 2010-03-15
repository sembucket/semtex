///////////////////////////////////////////////////////////////////////////////
// qmesh.C: Driver for qmesh: 2D quadrilateral unstructured mesh generator.
//
// Usage: qmesh [options] [file]
//   [options]:
//   -h       ... print this message
//   -g       ... disable X11 (and PostScript hardcopy) graphics
//   -m       ... output suitable for subsequent merging
//   -s <num> ... number of smoothing passes
//   -r <num> ... refinement coefficient (0--1)
//   -v[v...] ... increase verbosity level
//
// Window limits for graphics are set in file limits.sm, if it exists;
// if not, limits are generated from input data.
// The file limits.sm contains (in order): xmin xmax ymin ymax.
//
// Starting from an initial list of nodes and loops declared in the
// input file, first recursively subdivide each loop.  The loop
// subdivision algorithm used closely follows that given in Ref. [1].
// When this process is finished, each loop is subdivided into a
// binary tree of four-noded loops. 
//
// This information is transferred to a list of quads, which is
// "improved" by culling nodes and quads, after which the remainder is
// smoothed and printed up.
//
// -- Pseudocode description of loop subdivision:
// read in nodes which define a closed loop in the plane;
// 
// generate boundary offset nodes, subdividing original loop each time
// a 4-noded subloop is formed;
//
// while (any loop has more than 4 nodes) {
//   for (each loop with more than six nodes) {
//     for (each node in loop) generate list of visible nodes;
//     choose best splitting line;
//     subdivide splitting line;
//     split loop into two loops;
//   }
//   for (each six-noded loop) {
//     classify;
//     split into two subloops;
//   }
// }
// -- End.
//
// A binary tree is used to maintain the loop/subloop structure, see
// file loop.C.
//
// By convention, CCW is direction of loop traverses and positive angles.      
//
// References:
// ----------
// [1] J. A. Talbert & A. R. Parkinson,  1990.  Development of an automatic,
//     two-dimensional finite element mesh generator using quadrilateral
//     elements and Bezier curve boundary definition.  IJNME V29, 1551--1567.
// [2] F. S. Hill, Jr., 1990.  Computer Graphics.  Collier Macmillan.
// [3] R. Sedgewick, 1990.  Algorithms in C.  Addison-Wesley.
// [4] J. Z. Zhu, O. C. Zienkiewicz, E. Hinton & J. Wu, 1991.  A new approach
//     to the development of automatic quadrilateral mesh generation.
//     IJMME V32, 849--866.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>

char prog[] = "qmesh";

// -- Local routines.

static void     getArgs      (int, char**, int&, int&, istream*&);
static istream& operator >>  (istream&, List<Node*>&);
static istream& operator >>  (istream&, List<Loop*>&);
static int      loopDeclared (istream& s);
static void     connect      (List<Quad*>&);
static void     renumber     (List<Node*>&);
static void     deleteNodes  (List<Quad*>&);
static void     deleteQuads  (List<Quad*>&);
static void     smooth       (List<Node*>&);
static void     printNodes   (ostream&, List<Node*>&, const int);
static void     printMesh    (ostream&, List<Quad*>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for mesh generator.
// ---------------------------------------------------------------------------
{
  istream*    input;
  int         i, nsmooth = 0, merger = 0;
  List<Loop*> initial;
  List<Quad*> elements;

  getArgs (argc, argv, nsmooth, merger, input);

  // -- Load all predeclared nodes and loops.

  *input >> Global::nodeList;
  *input >> initial;

  // -- Start drawing of subdivision process.

  if (graphics) {
    initGraphics ("x11");
    drawBox      ();
  }

  ListIterator<Loop*> I (initial);
  Loop*               L;

  // -- Subdivide predeclared loops until all are quads.

  for (I.reset(); I.more(); I.next()) {
    L = I.current();

    if (graphics) drawLoop (L);
    L -> offset  ();
    if (graphics) drawLoop (L);
    L -> split   ();
    L -> quads   (elements);
  }

  // -- Join up all quads into an element mesh.

  connect (elements);

  // -- Try to improve mesh by Node and Quad elimination, see Ref. [4].

  do {
    i = Global::nodeList.length();

    deleteNodes (elements);
    connect     (elements);
    deleteQuads (elements);
    connect     (elements);

  } while (i != Global::nodeList.length());

  renumber (Global::nodeList);

  if (graphics) { eraseGraphics(); drawMesh (elements); }

  // -- Laplacian smoothing.

  for (i = 0; i < nsmooth; i++) {
    smooth (Global::nodeList); if (graphics) drawMesh (elements);
  }

  // -- Plot and output final mesh.

  if (graphics) {
    eraseGraphics ();
    drawBox       ();
    drawMesh      (elements);
    hardCopy      (elements);
    stopGraphics  ();
  }

  printNodes (cout, Global::nodeList, merger);
  printMesh  (cout, elements);

  return EXIT_SUCCESS;
}


static void getArgs (int       argc   ,
		     char**    argv   ,
		     int&      nsmooth,
		     int&      merger ,
		     istream*& input  )
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char buf[StrMax], c;
  char usage[]   =
    "Usage: %s [options] [file]\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -g       ... disable X11 (and PostScript hardcopy) graphics\n"
    "  -m       ... output suitable for subsequent merging\n"
    "  -s <num> ... number of smoothing passes\n"
    "  -r <num> ... refinement coefficient (0--1)\n"
    "  -v[v...] ... increase verbosity level\n"
    "\n"
    "Window limits can be set in file limits.sm: xmin xmax ymin ymax.\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'g':
      graphics = 0;
      break;
    case 'm':
      merger = 1;
      break;
    case 'r':
      if (*++argv[0])
	Global::refCoeff = atof(*argv);
      else {
	--argc;
	Global::refCoeff = atof(*++argv);
      }
      if (Global::refCoeff < 0.0 || Global::refCoeff > 1.0)
	error (prog, "refinement coefficient not in range 0--1", ERROR);
      if (Global::refCoeff > 0.7)
	error (prog, ": refinement coefficient is large", WARNING);
      break;
    case 's':
      if (*++argv[0])
	nsmooth = atoi(*argv);
      else {
	--argc;
	nsmooth = atoi(*++argv);
      }
      break;
    case 'v':
      do Global::verbose++; while (*++argv[0] == 'v');
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
      
  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) error (prog, "unable to open input file", ERROR);
  } else input = &cin;
}


static istream& operator >> (istream&     s,
			     List<Node*>& n)
// ---------------------------------------------------------------------------
// Node information for all predeclared loops exists in the following
// format:
//
// NN BOUNDARY NODES        { NN gives number of following nodes.
// 1  0.3  B  0.0 0.0       { Tag, prefsize, kind, x, y.
// 2  0.4  B  1.0 0.0       { Allowed kinds: B <==> fixed location boundary.
// ..                       {                O <==> fixed, but create offset.
// NN 0.1  B  0.0 0.1       {                I <==> interior (moveable).
//
// Read all this into n.  The first line is all on one line but following
// are free-format.  Kind information (B, O, I) upper or lower case.
//
// Node tags must be supplied in increasing order starting at 1,
// increment 1.
// ---------------------------------------------------------------------------
{
  char  routine[] = "operator >> (istream&, Loop<Node*>&)";
  char  err[StrMax], buf[StrMax];

  int            i, id, kind, npts, found;
  Point          pnt, P1, P2;
  real           size;
  register Node* N;

  if (!(s >> npts)) {
    sprintf (err, "problem reading initial number of nodes");
    message (routine, err, ERROR);
  }

  s.getline (buf, StrMax); upperCase (buf);
  if (!strstr (buf, "BOUNDARY NODES")) {
    sprintf (err, "can't locate nodes in loop: %s", buf);
    message (routine, err, ERROR);
  }

  for (i = 0; i < npts; i++) {
    s >> id >> size >> kind >> pnt;

    if (!s) {
      sprintf (err, "problem reading point %1d", i+1);
      message (routine, err, ERROR);
    } else if (++Global::nodeIdMax != id) {
      sprintf (err, "node %1d specified out of order (%1d)",
	       id, Global::nodeIdMax);
      message (routine, err, ERROR);
    }

    switch (kind) {
    case 0: N = new Node (id, pnt, size, Node::INTERIOR);              break;
    case 1: N = new Node (id, pnt, size, Node::INTERIOR_FIXED);        break;
    case 2: N = new Node (id, pnt, size, Node::LOOP_OFFSET_FIXED);     break;
    case 3: N = new Node (id, pnt, size, Node::LOOP_OFFSET_MOBILE);    break;
    case 4: N = new Node (id, pnt, size, Node::LOOP_BOUNDARY_FIXED);   break;
    case 5: N = new Node (id, pnt, size, Node::LOOP_BOUNDARY_MOBILE);  break;
    case 6: N = new Node (id, pnt, size, Node::DOMAIN_OFFSET_FIXED);   break;
    case 8: N = new Node (id, pnt, size, Node::DOMAIN_BOUNDARY_FIXED); break;
    default:
      sprintf (err, "read unused/unknown Node kind specifier: %1d", kind);
      message (routine, err, ERROR);
      break;
    }

    if (!(Global::exist (N)))
      n.add (N);
    else {
      sprintf (err, "Node %1d already allocated, check input", N -> ID());
      message (routine, err, ERROR);
    }    
  }

  Global::limits (P1, P2);

  return s;
}


static istream& operator >> (istream&     s,
			     List<Loop*>& l)
// ---------------------------------------------------------------------------
// Administer reading all the loops.
// ---------------------------------------------------------------------------
{
  char  routine[] = "operator >> (istream&, List<Loop*>&)";
  int   numnodes;
  Loop* L;
  
  while (numnodes = loopDeclared (s)) {
    L = new Loop (numnodes);
    s >> *L;
    l.add (L);
  }
  
  return s;
}


static int loopDeclared (istream& s)
// ---------------------------------------------------------------------------
// Test: is a loop declared for input?
//
// MM NODE LOOP             { MM gives number of following node tag numbers.
// 1 2 3 4 5 ... NN ...     { Loop is assumed closed by return to start tag.
//
// MM specifiers must be on lines of their own, but input is
// otherwise free-format.
//
// The first tag does not need to be re-specified at the end of the node loop
// tag list: MM is the number of nodes in the loop not including the return
// to the first node (e.g. a quad loop would have MM = 4, not 5).
//
// MM must be even.
// ---------------------------------------------------------------------------
{
  char routine[] = "loopDeclared";
  char buf[StrMax], err[StrMax];
  int  n;

  if (!(s >> n)) return 0;

  s.getline (buf, StrMax); upperCase (buf);
  if (!strstr (buf, "NODE LOOP")) {
    sprintf (err, "can't locate number of nodes in loop: %s", buf);
    message (routine, err, ERROR);
  }

  if (n & 1) {
    sprintf (err, "loop has odd number of points: %1d", n);
    message (routine, err, WARNING);
  }

  return n;
}


static void connect (List<Quad*>& elements)
// ---------------------------------------------------------------------------
// Visit all quads and for each node, add information about the nodes it
// is connected to.
// ---------------------------------------------------------------------------
{
  char routine[] = "connect";

  ListIterator<Quad*> q (elements);
  ListIterator<Node*> n (Global::nodeList);

  register int        i, found1, found2;
  register Quad*      Q;
  register Node*      N;
  register Node*      N1;
  register Node*      N2;
  
  for (; n.more(); n.next()) n.current() -> sever();

  for (; q.more(); q.next()) {
    Q = q.current();
    for (i = 0; i < 4; i++) {
      N1 = Q -> vertex[i];
      N2 = Q -> vertex[(i + 1) % 4];
      for (found1 = 0, found2 = 0, n.reset();
	   !(found1 && found2) && n.more();
	   n.next()) {
	N = n.current();
	if (!found1) if (found1 = (N == N1)) N -> xadd (N2);
	if (!found2) if (found2 = (N == N2)) N -> xadd (N1);
      }
    }
  }

  if (Global::verbose)
    for (n.reset(); n.more(); n.next()) {
      N = n.current();
      cout << routine << ": node ID: " << N -> ID()
	   << ", adjacency = " << N -> adjncy() << endl;
    }
}


static void smooth (List<Node*>& nodes)
// ---------------------------------------------------------------------------
// Laplacian smoothing.  Visit each Node, move non-boundary nodes to
// centroid of connected Nodes.
// ---------------------------------------------------------------------------
{
  ListIterator<Node*> n (nodes);
  register Node*      N;
  Point               cen;

  for (; n.more(); n.next()) { N  = n.current();
    cen = N -> centroid ();
    N -> setPos (cen);
  }
}


static void printNodes (ostream&     strm  ,
			List<Node*>& nodes ,
			const int    merger)
// ---------------------------------------------------------------------------
// Print Node information in FEML format.  BUT, if merger != 0, add
// information to indicate if output nodes are to be fixed (F), mobile
// (M) or unique (U).  This information will be used for possible
// subsequent merging and smoothing of output files.
// ---------------------------------------------------------------------------
{
  int                 i = 0;
  Node*               N;
  ListIterator<Node*> n (nodes);

  strm << "<NODES NUMBER=" << nodes.length() << ">" <<endl;
  for (n.reset(); n.more(); n.next()) {
    N = n.current ();
    strm << setw (5)  << ++i
	 << setw (16) << N -> pos().x 
	 << setw (16) << N -> pos().y
	 << setw (16) << 0.0;

    if (merger)
      switch (N -> classify()) { 
      case 'F': strm << "\tF"; break;
      case 'U': strm << "\tU"; break;
      default:  strm << "\tM"; break;
      }

    strm << endl;
  }
  strm << "</NODES>" << endl << endl;
}


static void printMesh (ostream&     strm,
		       List<Quad*>& mesh)
// ---------------------------------------------------------------------------
// Print Quad information in FEML format.
// ---------------------------------------------------------------------------
{
  int                 i, id = 0;
  ListIterator<Quad*> q (mesh);
  Quad*               Q;
  
  strm << "<ELEMENTS NUMBER=" << mesh.length() << ">" <<endl;
  for (q.reset(); q.more(); q.next()) {
    Q = q.current();
    strm << setw(5) << ++id << "  <Q>";
    for (i = 0; i < 4; i++) 
      strm << setw (5) << Q -> vertex[i] -> ID();
    strm << "  </Q>" << endl;
  }
  strm << "</ELEMENTS>" << endl;
}


static void deleteNodes (List<Quad*>& mesh)
// ---------------------------------------------------------------------------
// Improve mesh by node elimination.  See \S 3.1.1 in Ref [4].
// ---------------------------------------------------------------------------
{
  char  routine[] = "deleteNodes", err[StrMax];
  int   i, i1, i2, found;
  Node  *N;
  Quad  *Q, *Q1, *Q2;

  do {
    ListIterator<Node*> n (Global::nodeList);
    ListIterator<Quad*> q (mesh);

    for (found = 0; !found && n.more(); n.next()) {
      N = n.current();
      found = N -> adjncy() == 2 && N -> interior();
      if (found) {
	for (Q1 = 0, Q2 = 0; !(Q1 && Q2) && q.more(); q.next()) {
	  Q = q.current();
	  for (i = 0; i < 4; i++)
	    if (Q -> vertex[i] == N) {
	      if   (!Q1) { Q1 = Q; i1 = i; }
	      else       { Q2 = Q; i2 = i; }
	    }
	}
	if (!(Q1 && Q2)) {
	  sprintf (err, "node %1d marked, can't find two elements", N -> ID());
	  message (routine, err, ERROR);
	}
	Q1 -> vertex[i1] = Q2 -> vertex[(i2 + 2) % 4];
      }
    }
    if (found) {
      if (Global::verbose) 
	cout << routine << ": node " << N -> ID() << " deleted" << endl;
      Global::nodeList.remove (N);
      mesh            .remove (Q2);
    }
  } while (found);
}


static void deleteQuads (List<Quad*>& mesh)
// ---------------------------------------------------------------------------
// Improve mesh by element elimination.  See \S 3.1.2 in Ref [4].
// ---------------------------------------------------------------------------
{
  char  routine[] = "deleteQuads";
  int   i, i1, i2, found;
  Node  *N1, *N2;
  Quad  *Q, *P, *Q1, *Q2;

  do {
    ListIterator<Quad*> q (mesh);
    ListIterator<Quad*> p (mesh);

    for (found = 0; !found && q.more(); q.next()) {
      Q = q.current();
      for (i = 0; !found && i < 2; i++) {
	N1 = Q -> vertex [i];
	N2 = Q -> vertex [(i + 2) % 4];

	found = ((N1 -> adjncy() == 3) && N1 -> interior() &&
		 (N2 -> adjncy() == 3) && N2 -> interior());
      }
      if (found) {
	for (Q1 = 0, Q2 = 0; !(Q1 && Q2) && p.more(); p.next()) {
	  P = p.current();
	  if (P != Q) {
	    for (i = 0; i < 4; i++)
	      if (P -> vertex[i] == N2) P -> vertex[i] = N1;
	  }
	}
      }
    }
    if (found) {
      if (Global::verbose) 
	cout << routine << ": node " << N2 -> ID() << " deleted" << endl;
      Global::nodeList.remove (N2);
      mesh            .remove (Q);
    }
  } while (found);
}


static void renumber (List<Node*>& nodes)
// ---------------------------------------------------------------------------
// After "improving" mesh, there may be some holes in the Node ID numbers.
// This routine fixes that.
// ---------------------------------------------------------------------------
{
  register int        id;
  ListIterator<Node*> n (nodes);
  
  for (id = 0; n.more(); n.next()) n.current() -> renumber (++id);
}

