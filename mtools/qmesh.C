///////////////////////////////////////////////////////////////////////////////
// qmesh.C: Driver for qmesh: 2D quadrilateral unstructured mesh generator.
//
// Usage: qmesh [options] [file]
//   [options]:
//   -h       ... print this message
//   -g       ... disable X11 (and PostScript hardcopy) graphics
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
// binary tree of four-noded loops.  Then this information is
// transferred to a list of quads, which are smoothed and printed up.
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
// file loop.cc
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
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <qmesh.h>


char prog[] = "qmesh";

real        Global::refCoeff    = 0.0;
real        Global::gblSize     = 1.0;
int         Global::nodeIdMax   = 0;
int         Global::loopIdMax   = 0;
List<Node*> Global::nodeList;

// -- Local routines.

static void     getArgs      (int, char**, int&, ifstream&);
static istream& operator >>  (istream&, List<Node*>&);
static istream& operator >>  (istream&, List<Loop*>&);
static int      loopDeclared (istream& s);
static void     connect      (List<Quad*>&);
static void     smooth       (List<Node*>&);
static void     printNodes   (ostream&, List<Node*>&);
static void     printMesh    (ostream&, List<Quad*>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for mesh generator.
// ---------------------------------------------------------------------------
{
  int         i, nsmooth = 0;
  ifstream    infile;
  List<Loop*> initial;
  List<Quad*> elements;

  getArgs (argc, argv, nsmooth, infile);

  // -- Load all predeclared nodes and loops.

  infile >> Global::nodeList;
  infile >> initial;
  infile.close ();

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

  // --  Join up all quads into mesh, and smooth.

  connect (elements);
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

  printNodes (cout, Global::nodeList);
  printMesh  (cout, elements);

  return EXIT_SUCCESS;
}


static void getArgs (int       argc   ,
		     char**    argv   ,
		     int&      nsmooth,
		     ifstream& infile )
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
      do verbose++; while (*++argv[0] == 'v');
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if   (argc == 1) infile.open   (*argv, ios::in);
  else             infile.attach (0);

  if (!infile) {
    sprintf (buf, "unable to open file: %s", *argv);
    error   (prog, buf, ERROR);
  }
}


static istream& operator >> (istream&     s,
			     List<Node*>& n)
// ---------------------------------------------------------------------------
// Node information for all predeclared loops exists in the following
// format:
//
// NN boundary nodes        { NN gives number of following nodes.
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

  int            i, id, npts, found;
  Point          pnt;
  real           size;
  char           kind;
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

    switch (toupper (kind)) {
    case 'B': N = new Node (id, pnt, size, Node::BOUNDARY); break;
    case 'I': N = new Node (id, pnt, size, Node::INTERIOR); break;
    case 'O': N = new Node (id, pnt, size, Node::OFFSET  ); break;
    default:
      sprintf (err, "read unknown Node kind specifier: %c", kind);
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
// MM node loop             { MM gives number of following node tag numbers.
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
// Return 1 if it has, else zero.
// ---------------------------------------------------------------------------
{
  char           err[StrMax], routine[] = "Global::exist";
  int            found = 0;
  register Node* oldNode;
  const Point    P    = N -> pos();
  const real     size = lengthScale();
  const real     TOL  = 0.001;
  
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


static void connect (List<Quad*>& elements)
// ---------------------------------------------------------------------------
// Visit all quads and for each node, add information about the nodes it
// is connected to.
// ---------------------------------------------------------------------------
{
  ListIterator<Quad*> q (elements);
  ListIterator<Node*> n (Global::nodeList);

  register int        i, found1, found2;
  register Quad*      Q;
  register Node*      N;
  register Node*      N1;
  register Node*      N2;

  for (q.reset(); q.more(); q.next()) {
    Q = q.current();
    for (i = 0; i < 4; i++) {
      N1 = Q -> vertex[i];
      N2 = Q -> vertex[(i + 1) % 4];
      for (found1 = 0, found2 = 0, n.reset(); n.more(); n.next()) {
	N = n.current();
	if (!found1) if (found1 = (N == N1)) N -> xadd (N2);
	if (!found2) if (found2 = (N == N2)) N -> xadd (N1);
	if (found1 && found2) break;
      }
    }
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


static void printNodes (ostream&     strm ,
			List<Node*>& nodes)
// ---------------------------------------------------------------------------
// Print Node information in FEML format.
// ---------------------------------------------------------------------------
{
  Node*               N;
  ListIterator<Node*> n (nodes);

  strm << "<NODES NUMBER=" << nodes.length() << ">" <<endl;
  for (n.reset(); n.more(); n.next()) {
    N = n.current();
    strm << setw (5)  << N -> ID()
	 << setw (16) << N -> pos().x 
	 << setw (16) << N -> pos().y
	 << setw (16) << 0.0 << endl;
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
  strm << "</ELEMENTS>" << endl << endl;
}

