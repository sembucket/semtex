///////////////////////////////////////////////////////////////////////////////
// merge.C: merge two files of nodes and loops into one.
//
// usage: merge file1 file2
//
// Method: read in the predefined Node section of first file.  Then
// read in the same information from the second file, but xadd the new
// Nodes to the list of old Nodes, and at the same time build a lookup
// table of indices of the new Nodes and pointers to those in the
// extended list.  Print up the data for the extended list of Nodes,
// followed by the Node loop information in the first file and the
// second file.  It is assumed that within each file (taken separately)
// there are no redundant node indentifiers.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <qmesh.h>

char prog[] = "merge";

real        Global::refCoeff  = 0.0;
real        Global::gblSize   = 1.0;
int         Global::nodeIdMax = 0;
int         Global::loopIdMax = 0;
int         Global::verbose   = 0;
List<Node*> Global::nodeList;

Node*       Global::exist (const Node*);
static void buildTable    (istream&, int, int, vector<Node*>&);
static void printNodes    (ostream&, List<Node*>&);
static void processLoops  (ostream&, istream&,  int, vector<Node*>&);
static int  loopDeclared  (istream&);



int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char          usage[] = "Usage: merge file1 file2\n";
  char          buf[StrMax], err[StrMax];
  int           n1, n2;
  vector<Node*> lut;

  // -- Open files.

  if (argc != 3) message (prog, usage, ERROR);
  ifstream f1 (argv[1]), f2 (argv[2]);
  if (!f1) {
    sprintf (err, "can't open file: %s", argv[1]);
    message (prog, err, ERROR);
  }
  if (!f2) {
    sprintf (err, "can't open file: %s", argv[2]);
    message (prog, err, ERROR);
  }

  // -- How many Nodes are predeclared in total?

  if (!(f1 >> n1)) {
    sprintf (err, "problem reading initial number of nodes in %s", *argv[1]);
    message (prog, err, ERROR);
  }
  f1.getline (buf, StrMax); upperCase (buf);
  if (!strstr (buf, "BOUNDARY NODES")) {
    sprintf (err, "%s: can't locate nodes in loop: %s", *argv[1], buf);
    message (prog, err, ERROR);
  }

  if (!(f2 >> n2)) {
    sprintf (err, "problem reading initial number of nodes in %s", *argv[1]);
    message (prog, err, ERROR);
  }
  f2.getline (buf, StrMax); upperCase (buf);
  if (!strstr (buf, "BOUNDARY NODES")) {
    sprintf (err, "%s: can't locate nodes in loop: %s", *argv[1], buf);
    message (prog, err, ERROR);
  }

  lut.setSize (n1 + n2);

  // -- Input nodes to global node list, build LUT, then print list.

  buildTable (f1, 0,  n1, lut);
  buildTable (f2, n1, n2, lut);

  printNodes (cout, Global::nodeList);

  // -- Input/output loops with revised node numbers.

  processLoops (cout, f1,  0, lut);
  processLoops (cout, f2, n1, lut);

  return (EXIT_SUCCESS);
}


static void buildTable (istream&       strm  ,
			const int      offset,
			const int      nnodes,
			vector<Node*>& table )
// ---------------------------------------------------------------------------
// Xadd new (previously undeclared) nodes to global list and LUT.
// ---------------------------------------------------------------------------
{
  char  routine[] = "buildTable";
  char  err[StrMax], buf[StrMax];
  int   i, id;
  real  size;
  char  kind;
  Point pnt;
  Node  *N, *O;

  for (i = 0; i < nnodes; i++) {
    strm >> id >> size >> kind >> pnt;

    if (!strm) {
      sprintf (err, "problem reading point %1d", i + 1);
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

    if ((O = Global::exist (N)) == N) {
      delete N;
      id = ++Global::nodeIdMax;
      switch (toupper (kind)) {
      case 'B': N = new Node (id, pnt, size, Node::BOUNDARY); break;
      case 'I': N = new Node (id, pnt, size, Node::INTERIOR); break;
      case 'O': N = new Node (id, pnt, size, Node::OFFSET  ); break;
      }
      Global::nodeList.add (N);
    }
    table [i + offset] = O;
  }
}


Node* Global::exist (const Node* N)
// ---------------------------------------------------------------------------
// Check if a Node corresponding to N has already been created.
// Return pointer to old Node if it has, else pointer to new Node.
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
  
  return N;
}


static void printNodes (ostream&     ostrm,
			List<Node*>& nodes)
// ---------------------------------------------------------------------------
// Print up list of unique Nodes from both files.
// ---------------------------------------------------------------------------
{
  ListIterator<Node*> n (nodes);

  ostrm << nodes.length() << " boundary nodes" << endl;
  for (; n.more(); n.next()) ostrm << *n.current() << endl;
  ostrm << endl;
}


static void processLoops (ostream&       ostrm ,
			  istream&       istrm , 
			  int            offset,
			  vector<Node*>& table)
// ---------------------------------------------------------------------------
// Read in loop information, output revised Node IDs.
// ---------------------------------------------------------------------------
{
  int i, j, n;

  while (n = loopDeclared (istrm)) {
    ostrm << n << " node loop" << endl;
    for (i = 0; i < n; i++) {
      istrm >> j;
      ostrm << setw(5) << table[j + offset - 1] -> ID();
      if (!((i + 1) % 10)) ostrm << endl;
    }
    if (n % 10) ostrm << endl;
    ostrm << endl;
  }
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
