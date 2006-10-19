///////////////////////////////////////////////////////////////////////////////
// merge.C: merge two files of nodes and elements into one, smooth, output.
//
// usage: merge file1 [file2]
//
// Method: read in the predefined Node section of first file.  Then
// read in the same information from the second file, but xadd the new
// Nodes to the list of old Nodes, and at the same time build a lookup
// table of indices of the new Nodes and pointers to those in the
// extended list.  Print up the data for the extended list of Nodes,
// followed by the Element information in the first file and the
// second file.  It is assumed that within each file (taken
// separately) there are no redundant node indentifiers.
//
// If file2 is not specified, merge acts on the single input file,
// i.e. just acts to smooth the mesh.
//
// The file formats are (nearly) standard FEML, with NODES and
// ELEMENTS sections.  Each node has an additional character label: U,
// F or M.  U ==> a (fixed) mesh boundary node, F ==> a fixed mesh
// interior node and M ==> a moveable node.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>

char prog[] = "merge";

static const int NSMOOTH = 0;

static int  getNum      (istream*, const char*);
static void buildTable  (istream*, const int, const int, Node**);
static void getElements (istream*, Node**, Quad**, const int);
static void connect     (List<Node*>&, Quad**, const int);
static void smooth      (List<Node*>&, const int);
static void printUp     (ostream&, const List<Node*>&, Quad**, const int);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char     usage[] = "Usage: merge file1 file2\n";
  char     err[StrMax];
  int      nnod1 = 0, nnod2 = 0, nel1 = 0, nel2 = 0;
  ifstream *f1, *f2 = 0;
  Node**   nodeLUT;
  Quad**   Elmt;

  // -- Open files.

  switch (argc) {
  case 3:			// -- Standard usage: merge and smooth.
    f2 = new ifstream (argv[2]);
    if (!*f2) {
      sprintf (err, "can't open file: %s", argv[2]);
      message (prog, err, ERROR);
    }
    // -- Fall through.
  case 2:			// -- Just apply a smoothing to single file.
    f1 = new ifstream (argv[1]);
    if (!*f1) {
      sprintf (err, "can't open file: %s", argv[1]);
      message (prog, err, ERROR);
    }
    break;
  default:
    message (prog, usage, ERROR);
    break;
  }

  // -- Find number of Nodes in each file.

  nnod1 = getNum (f1, "NODES");
  nnod2 = getNum (f2, "NODES");

  nodeLUT = new Node* [nnod1 + nnod2];

  // -- Input nodes to global node list, build LUT.

  buildTable (f1, 0,     nnod1, nodeLUT);
  buildTable (f2, nnod1, nnod2, nodeLUT);

  // -- Find number of elements in each file.

  nel1 = getNum (f1, "ELEMENTS");
  nel2 = getNum (f2, "ELEMENTS");

  Elmt = new Quad* [nel1 + nel2];

  // -- Input element information, install Node*s from LUT.

  getElements (f1, nodeLUT,         Elmt,        nel1);
  getElements (f2, nodeLUT + nnod1, Elmt + nel1, nel2);

  initGraphics ("x11");
  drawBox      ();
  drawMesh     (Global::quadList);

  // -- Process combined element information.

  connect (Global::nodeList, Elmt, nel1 + nel2);
  smooth  (Global::nodeList, NSMOOTH);

  eraseGraphics ();
  drawBox       ();
  drawMesh      (Global::quadList);

  Global::verbose = 1;

  hardCopy      (Global::quadList);

  stopGraphics  ();

  // -- Print up.

  printUp (cout, Global::nodeList, Elmt, nel1 + nel2);

  return (EXIT_SUCCESS);
}


static int getNum (istream*    strm   ,
		   const char* keyword)
// ---------------------------------------------------------------------------
// Look for keyword NUMBER=num, extract & return num, leave stream after EOL.
// ---------------------------------------------------------------------------
{
  char routine[] = "getNum";
  char buf[StrMax], err[StrMax], key[StrMax], sep[] = "=>", *tok;
  int  num;

  if (!strm) return 0;

  sprintf (key, "%s NUMBER", keyword);

  while (strm -> getline (buf, StrMax)) {
    upperCase (buf);
    if (strstr (buf, key)) break;
  }

  if (!strstr (buf, key)) {
    sprintf (err, "can't locate number of %s in file", keyword);
    message (routine, err, ERROR);
  }

  tok = strtok (buf, sep);
  tok = strtok (0,   sep);
  sscanf (tok, "%d", &num);
  
  return num;
}


static void buildTable (istream*       strm  ,
			const int      offset,
			const int      nnodes,
			Node**         table )
// ---------------------------------------------------------------------------
// Xadd new (previously undeclared) nodes to global list and LUT.
// ---------------------------------------------------------------------------
{
  char  routine[] = "buildTable";
  char  err[StrMax];
  int   i, id;
  real  z;
  char  kind;
  Point pnt;
  Node  *N, *O;

  if (!strm) return;

  for (i = 0; i < nnodes; i++) {
    *strm >> id >> pnt >> z >> kind;

    if (!*strm) {
      sprintf (err, "problem reading point %1d", i + 1);
      message (routine, err, ERROR);
    }

    z = 0.0;

    switch (toupper (kind)) {
    case 'U': N = new Node (id, pnt, z, Node::DOMAIN_BOUNDARY_FIXED); break;
    case 'F': N = new Node (id, pnt, z, Node::INTERIOR_FIXED);        break;
    case 'M': N = new Node (id, pnt, z, Node::INTERIOR);              break;
    default:
      sprintf (err, "read unknown Node kind specifier: %c", kind);
      message (routine, err, ERROR);
      break;
    }

    if (O = Global::exist (N))
      table [i + offset] = O;
    else {
      delete N;
      id = ++Global::nodeIdMax;
      switch (toupper (kind)) {
      case 'U': N = new Node (id, pnt, z, Node::DOMAIN_BOUNDARY_FIXED); break;
      case 'F': N = new Node (id, pnt, z, Node::INTERIOR_FIXED);        break;
      case 'M': N = new Node (id, pnt, z, Node::INTERIOR);              break;
      }
      Global::nodeList.add (N);
      table [i + offset] = N;
    }
  }
}


static void getElements (istream*  strm   ,
			 Node**    nodeLUT,
			 Quad**    Elmt   ,
			 const int nElmt  )
// ---------------------------------------------------------------------------
// Read in declared elements, match up vertex IDs with Node*'s.
// ---------------------------------------------------------------------------
{
  char         routine[] = "getElements", buf[StrMax], err[StrMax];
  register int i, j, id;
  Quad*        Q;

  if (!strm) return;

  for (i = 0; i < nElmt; i++) {
    *strm >> id >> buf;
    if (!*strm || !strstr (buf, "<Q>")) {
      sprintf (err, "problem reading element info, expected <Q>, got %s", buf);
      message (routine, err, ERROR);
    }
    Global::quadList.add (Elmt[i] = Q = new Quad);
    for (j = 0; j < 4; j++) {
      *strm >> id;
      Q -> vertex[j] = nodeLUT[--id];
    }
    *strm >> buf;
  }
}


static void connect (List<Node*>& Nodes,
		     Quad**       Elmt ,
		     const int    nElmt)
// ---------------------------------------------------------------------------
// We know the connectivity at Element level, but not at the mesh level.
// Traverse Elements and update connectivity in Node vector.
// ---------------------------------------------------------------------------
{
  register int        i, j, found1, found2;
  register Node       *N, *N1, *N2;
  ListIterator<Node*> n (Nodes);

  for (i = 0; i < nElmt; i++)
    for (j = 0; j < 4; j++) {
      N1 = Elmt[i] -> vertex[j];
      N2 = Elmt[i] -> vertex[(j + 1) % 4];
      found1 = found2 = 0;
      for (n.reset(); n.more(); n.next()) {
	N = n.current();
	if (!found1) if (found1 = (N == N1)) N -> xadd (N2);
	if (!found2) if (found2 = (N == N2)) N -> xadd (N1);
	if ( found1  &&  found2 ) break;
      }
    }
}


static void smooth (List<Node*>& Nodes,
		    const int    nPass)
// ---------------------------------------------------------------------------
// Run through all Nodes and move them to centroid of connected Nodes.
// (But only if they're Interior Nodes.)
// ---------------------------------------------------------------------------
{
  register int        k;
  register Node*      N;
  ListIterator<Node*> n (Nodes);
  Point               cen;
  
  for (k = 0; k < nPass; k++)
    for (n.reset(); n.more(); n.next()) {
      N   = n.current();
      cen = N -> centroid();
      N -> setPos (cen);
    }
}


static void printUp (ostream&            strm ,
		     const List<Node*>&  Nodes,
		     Quad**              Elmt ,
		     const int           nElmt)
// ---------------------------------------------------------------------------
// The output format is the same as the input format to ensure that merge
// can operate on its own output.
// ---------------------------------------------------------------------------
{
  register int        i, j;
  ListIterator<Node*> n (Nodes);
  register Node*      N;

  // -- Print up Nodes.

  strm << setprecision (6);
  strm.setf (ios::scientific, ios::floatfield);

  strm << "<NODES NUMBER=" << Nodes.length() << ">" << endl;
  for (i = 0; n.more(); n.next()) {
    N = n.current();
    strm << setw  (5) << ++i
	 << setw (15) << N -> pos().x
	 << setw (15) << N -> pos().y
	 << "  0.0  " << N -> classify()
	 << endl;
  }

  strm << "</NODES>" << endl << endl;
  
  // -- Print elements.

  strm << "<ELEMENTS NUMBER=" << nElmt << ">" << endl;
  for (i = 0; i < nElmt; i++) {
    strm << setw (5) << i + 1;
    strm << "\t<Q>";
    for (j = 0; j < 4; j++) strm << setw(5) << Elmt[i] -> vertex[j] -> ID();
    strm << "\t</Q>" << endl;
  }
  strm << "</ELEMENTS>" << endl << endl;
}
