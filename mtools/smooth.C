///////////////////////////////////////////////////////////////////////////////
// smooth.cc:  read qmesh output (format) and apply Laplacian smoothing
// to all interior nodes.  This allows hand-generated meshes to be smoothed.
// Meshes don't have to be composed of quads.
// Print on cout in same format.
//
// smooth [-h] [-s N] [file]
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <qmesh.h>

#ifdef GRAPHICS
  #include <sm_options.h>
  #include <sm_declare.h>
#endif


static char  prog[] = "smooth";

static void getArgs     (int, char**, int&, ifstream&);
static int  getVertices (ifstream&, vector<Node*>&);
static int  getElements (ifstream&, vector<Node*>&, vector<Node*>*&);
static void connect     (vector<Node*>&, vector<Node*>*&, const int);
static void smooth      (vector<Node*>&, const int);
static void printUp     (ostream&, const vector<Node*>&,
			 const vector<Node*>*&, const int);

static void smStart (const char*);
static void smBox   ();
static void smDraw  (const vector<Node*>*&, const int);
static void smStop  ();


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver for other routines.
// ---------------------------------------------------------------------------
{
  ifstream       infile;
  int            nPass = 1;	// -- Number of smoothing passes.
  int            nVert, nEl;	// -- Numbers read from input.
  vector<Node*>  vertices;
  vector<Node*>* elements;

  getArgs (argc, argv, nPass, infile);

  seekBlock (infile, "Mesh");

  nVert = getVertices (infile, vertices);
  nEl   = getElements (infile, vertices, elements);

  endBlock (infile);

  smStart ("x11");
  smBox   ();

  connect  (vertices, elements, nEl);
  smooth   (vertices, nPass);

  smDraw  (elements, nEl);
  smStop  ();

  printUp  (cout, vertices, elements, nEl);
  
  return EXIT_SUCCESS;
}


static void getArgs (int       argc  ,
		     char**    argv  ,
		     int&      npass ,
		     ifstream& infile)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last (optional) argument is the name of an input file.
// ---------------------------------------------------------------------------
{
  char buf[StrMax], c;
  char usage[] = "Usage: %s [options] [file]\n"
    "  [options]:\n"
    "  -h   ... print this message\n"
    "  -s N ... use N smoothing passes [Default: 1]\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 's':
      if (*++argv[0])
        npass = atoi(*argv);
      else {
        --argc;
        npass = atoi(*++argv);
      }
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
}


static int getVertices (ifstream&      S,
			vector<Node*>& V)
// ---------------------------------------------------------------------------
// Read in declared vertices and create corresponding Nodes.
//
// Note that we check that input vertex IDs are given in order, since
// later the element node IDs can be used to index into V without clashes.
// ---------------------------------------------------------------------------
{
  char  routine[] = "getVertices", err[StrMax], buf[StrMax], *c;
  int   i, id, num, refId;
  Node* N;
  Point P, R;
  real  size, x, y, z;

  S >> num >> buf;
  refId = num; 
  
  c = buf;
  while (*c = toupper (*c)) c++;
  if (!strstr (buf, "VERTICES"))
    error (routine, "problem reading number of vertices", ERROR);

  V.setSize (num);
  for (i = 0; i < num; i++) {

    // -- Input unreflected Nodes.

    S >> id >> size >> buf >> x >> y >> z;
    if (id != i + 1)
      error (routine, "input list of Vertices out of order: check", ERROR);

    P.x = x; P.y = y;
    if (buf[0] == 'B')
      N = new Node (id, P, 1.0, Node::BOUNDARY);
    else if (buf[0] == 'I')
      N = new Node (id, P, 1.0);
    else {
      sprintf (err, "unknown Node type: %s", buf);
      error   (routine, err, ERROR);
    }
    V[i] = N;
  }

  return num;
}


static int getElements (ifstream&       S,
			vector<Node*>&  V,
			vector<Node*>*& E)
// ---------------------------------------------------------------------------
// Read in declared elements, match up vertex IDs with Node*'s.
// ---------------------------------------------------------------------------
{
  char         routine[] = "getElements", buf[StrMax], *c;
  register int i, j;
  int          id, M, J;
  const int    N = V.getSize();

  S >> M >> buf;
  
  c = buf;
  while (*c = toupper (*c)) c++;
  if (!strstr (buf, "ELEMENT"))
    error (routine, "problem reading number of elements", ERROR);

  E = new vector<Node*> [M];
  
  for (i = 0; i < M; i++) {
    S >> id >> J;
    E[i].setSize (J);
    for (j = 0; j < J; j++) {
      S >> id;
      E[i][j] = V[id - 1];
    }
  }

  return M;
}


static void connect (vector<Node*>&  V  ,
		     vector<Node*>*& E  ,
		     const int       nEl)
// ---------------------------------------------------------------------------
// We know the connectivity at Element level, but not at the mesh level.
// Traverse Elements and update connectivity in Node vector.
// ---------------------------------------------------------------------------
{
  register int   i, j, k, found1, found2;
  int            J;
  const int      K = V.getSize();
  register Node* N;
  register Node* N1;
  register Node* N2;

  for (i = 0; i < nEl; i++) {
    J = E[i].getSize();
    for (j = 0; j < J; j++) {
      N1 = E[i][j];
      N2 = E[i][(j + 1) % J];
      found1 = found2 = 0;
      for (k = 0; k < K; k++) {
	N = V[k];
	if (!found1) if (found1 = (N == N1)) N -> xadd (N2);
	if (!found2) if (found2 = (N == N2)) N -> xadd (N1);
	if ( found1  &&  found2 ) break;
      }
    }
  }
}


static void smooth (vector<Node*>& V    ,
		    const int      nPass)
// ---------------------------------------------------------------------------
// Run through all Nodes and move them to centroid of connected Nodes.
// (But only if they're Interior Nodes.)
// ---------------------------------------------------------------------------
{
  register int   i, k;
  const int      I = V.getSize();
  register Node* N;
  Point          cen;
  
  for (k = 0; k < nPass; k++)
    for (i = 0; i < I; i++) {
      N   = V[i];
      cen = N -> centroid();
      N -> setPos (cen);
    }
}


static void printUp (ostream&              S  ,
		     const vector<Node*>&  V  ,
		     const vector<Node*>*& E  ,
		     const int             nEl)
// ---------------------------------------------------------------------------
// The output format is the same as the input format to ensure that reflect
// can operate on its own output.
// ---------------------------------------------------------------------------
{
  register int i, j;
  const int    N = V.getSize();
  int          J;

  // -- Print up Nodes.

  S << "Mesh {" << endl;

  S << N << "  Vertices" << endl;
  for (i = 0; i < N; i++)
    S << *V[i] << setw (10) << 0.0 << endl;

  S << endl;
  
  // -- Print elements.

  S << nEl << "  Elements" << endl;
  for (i = 0; i < nEl; i++) {
    S << i + 1;
    J = E[i].getSize();
    S << setw (5) << J;
    for (j = 0; j < J; j++) S << setw(5) << E[i][j] -> ID();
    S << endl;
  }
  
  // -- Print up boundary Nodes (this will need hand editing).

  S << endl << "1  Boundary";

  S << endl;
  S << "1  1  ";
  j = 0;
  for (i = 0; i < N; i++) j += (V[i] -> interior()) ? 0 : 1;
  S << setw (5) << j;
  for (i = 0; i < N; i++)
    if (!V[i] -> interior()) S << setw (5) << V[i] -> ID();

  S << endl;
  S << endl << "0  Curves" << endl;
  
  S << "}" << endl;
}


static void smStart (const char* device)
// ---------------------------------------------------------------------------
// Do whatever is needed to start up.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  vector<float> pp (1);
  pp[0] = 41.0;

  if (sm_device( (char*) device)) {
    cerr << "unable to initialize plotting device" << endl;
    exit (EXIT_FAILURE);
  }
 
  sm_graphics ();
  sm_location (3000, 31000, 3000, 31000);
  sm_defvar   ("TeX_strings", "1");
  sm_expand   (1.0);
  sm_ptype    (pp(), 1);
  sm_lweight  (1);
  sm_erase    ();
  sm_window   (1, 1, 1, 1);
#endif
}


static void smStop ()
// ---------------------------------------------------------------------------
// Shut down graphics.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  char a;

  sm_alpha ();

  cerr << "Press <return> to exit" << endl;
  sm_redraw (0);
  cin.get(a);
#endif
}


static void smBox ()
// ---------------------------------------------------------------------------
// Underloaded version of drawBox.  Look for limits in file "limits.sm",
// otherwise use [-10, -10], [10, 10].
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  Point Pmax, Pmin, Centre, Range;
  float xmin, xmax, ymin, ymax;

  ifstream file ("limits.sm");

  if (!file) {
    xmin = -10.0;
    xmax =  10.0;
    ymin = -10.0;
    ymax =  10.0;
  } else
    file >> xmin >> xmax >> ymin >> ymax;

  sm_limits (xmin, xmax, ymin, ymax);
  sm_box    (1, 2, 0, 0);
  sm_gflush ();
  sm_box    (1, 2, 0, 0);
  sm_expand (1.6);
  sm_gflush ();
  sm_gflush ();
#endif
}


static void smDraw (const vector<Node*>*& E  ,
		    const int             nEl)
// ---------------------------------------------------------------------------
// Draw elements of mesh.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  int           i, j, J;
  vector<float> x(32);
  vector<float> y(32);
  
  for (i = 0; i < nEl; i++) {
    J = E[i].getSize();
    for (j = 0; j < J; j++) {
      x[j] = E[i][j] -> pos().x;
      y[j] = E[i][j] -> pos().y;
    }
    sm_conn   (x(), y(), J);
    sm_draw   (x[0], y[0]);
    sm_points (x(), y(), J);
  }
      
  sm_gflush ();
  sm_gflush ();
#endif
}

