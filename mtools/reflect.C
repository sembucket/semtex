///////////////////////////////////////////////////////////////////////////////
// reflect.C:  read qmesh output and generate planar reflection about
// either x or y axis, or generate rotation about named point.  Print on cout.
//
// reflect [-x || -y] [-r x0 y0 ang] [-h] [file]
//
// Either X or Y reflections change the sense of rotation around elements.
// Angular rotations are taken CCW, measured in degrees.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <qmesh.h>


enum   Action { UNDEFINED, X, Y, R };
static char       prog[] = "reflect";
static const real D2R    = 1.0 / 57.29577951308232087721;

static void getArgs     (int, char**, Action&, Point&, real&, ifstream&);
static int  getVertices (ifstream&, vector<Node*>&, const Action,
			 const Point&, const real);
static int  getElements (ifstream&, vector<Node*>&, vector<Node*>*&);
static void printUp     (ostream&, const vector<Node*>&,
			 const vector<Node*>*&, const int, const Action);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver for other routines.
// ---------------------------------------------------------------------------
{
  ifstream       infile;
  Action         mirror = UNDEFINED;
  Point          axis;
  real           angle;
  int            nVert, nEl;	// -- Numbers read from input.
  vector<Node*>  vertices;
  vector<Node*>* elements;

  getArgs (argc, argv, mirror, axis, angle, infile);

  seekBlock (infile, "Mesh");

  nVert = getVertices (infile, vertices, mirror, axis, angle);
  nEl   = getElements (infile, vertices, elements);

  endBlock (infile);

  printUp  (cout, const_cast<const vector<Node*>&>(vertices),
	    const_cast<const vector<Node*>*&>(elements), nEl, mirror);
  
  return EXIT_SUCCESS;
}


static void getArgs (int       argc  ,
		     char**    argv  ,
		     Action&   mirror,
		     Point&    axis  ,
		     real&     angle ,
		     ifstream& infile)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last (optional) argument is the name of an input file.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[] = "Usage: %s [-x || -y] [-r x0 y0 ang] [-h] [file]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'x':
      mirror = X;
      break;
    case 'y':
      mirror = Y;
      break;
    case 'r':
      mirror = R;
      --argc;
      axis.x = atof (*++argv);
      --argc;
      axis.y = atof (*++argv);
      --argc;
      angle  = atof (*++argv) * D2R;
      break;
    default:
      sprintf (buf, usage, prog);
      cerr << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (mirror == UNDEFINED) {
    sprintf (buf, usage, prog);
    cout << buf;
    exit (EXIT_FAILURE);
  }
  
  if   (argc == 1) infile.open   (*argv, ios::in);
  else             infile.attach (0);

  if (!infile) {
    sprintf (buf, "unable to open file: %s", *argv);
    error   (prog, buf, ERROR);
  }
}


static int getVertices (ifstream&      S,
			vector<Node*>& V,
			const Action   A,
			const Point&   O,
			const real     t)
// ---------------------------------------------------------------------------
// Read in declared vertices and create corresponding Nodes.
//
// Vertices for reflected Nodes are also generated in the upper locations
// of V, so there is a rule for getting from the original Nodes to the
// reflected Nodes and vice versa.  If the reflected Nodes land in the
// same locations (within rounding) as the originals, the Node pointers are
// the same as the unreflected ones.  Later on we have to construct
// a list of the unique Nodes to print out in the header though.
//
// Note that we check that input vertex IDs are given in order, since
// later the element node IDs can be used to index into V without clashes.
// ---------------------------------------------------------------------------
{
  char  routine[] = "getVertices", err[StrMax], buf[StrMax], *c;
  int   i, id, num, refId;
  Node* N;
  Point P, Pd;
  real  size, x, y, z, cc = cos (t), ss = sin (t);

  S >> num >> buf;
  refId = num; 
  
  c = buf;
  while (*c = toupper (*c)) c++;
  if (!strstr (buf, "VERTICES"))
    error (routine, "problem reading number of vertices", ERROR);

  V.setSize (num+num);
  for (i = 0; i < num; i++) {

    // -- Input unreflected Nodes.

    S >> id >> size >> buf >> x >> y >> z;
    if (id != i + 1)
      error (routine, "input list of Vertices out of order: check", ERROR);

    P.x = x; P.y = y;
    if (buf[0] == 'B')
      N = new Node (id, P, 1.0, Node::DOMAIN_BOUNDARY_FIXED);
    else if (buf[0] == 'I')
      N = new Node (id, P, 1.0);
    else {
      sprintf (err, "unknown Node type: %s", buf);
      error   (routine, err, ERROR);
    }
    V[i] = N;
   
    // -- Generate reflected/rotated Nodes.

    switch (A) {
    case X:
      Pd.x = -P.x; Pd.y =  P.y;
      break;
    case Y:
      Pd.x =  P.x; Pd.y = -P.y;
      break;
    case R:
      Pd.x = P.x*cc - P.y*ss + O.x*(1.0-cc) + O.y*ss;
      Pd.y = P.x*ss + P.y*cc + O.y*(1.0-cc) - O.x*ss;
      break;
    }

    if (Pd.distance (P) > EPSSP)
      N = new Node (++refId, Pd, 1.0,
		    (N -> interior()) ?
		    Node::INTERIOR : Node::DOMAIN_BOUNDARY_FIXED);

    V[i + num] = N;
  }

  return num;
}


static int getElements (ifstream&       S,
			vector<Node*>&  V,
			vector<Node*>*& E)
// ---------------------------------------------------------------------------
// Read in declared elements, match up vertex IDs with Node*'s.
// Generate reflected elements.
//
// Note that we know how many elements there'll be (twice the number declared
// in input) but not the number of unique vertices, because some may retain
// their identity under reflection.
// ---------------------------------------------------------------------------
{
  char         routine[] = "getElements", buf[StrMax], *c;
  register int i, j;
  int          id, M, J;
  const int    N = V.getSize() >> 1;

  S >> M >> buf;
  
  c = buf;
  while (*c = toupper (*c)) c++;
  if (!strstr (buf, "ELEMENT"))
    error (routine, "problem reading number of elements", ERROR);

  E = new vector<Node*> [M + M];
  
  for (i = 0; i < M; i++) {
    S >> id >> J;
    E[i    ].setSize (J);
    E[i + M].setSize (J);
    for (j = 0; j < J; j++) {
      S >> id;
      E[i    ][j] = V[    id - 1];
      E[i + M][j] = V[N + id - 1];
    }
  }

  return M;
}


static void printUp (ostream&              S  ,
		     const vector<Node*>&  V  ,
		     const vector<Node*>*& E  ,
		     const int             nEl,
		     const Action          A  )
// ---------------------------------------------------------------------------
// The output format is the same as the input format to ensure that reflect
// can operate on its own output.
// ---------------------------------------------------------------------------
{
  List<Node*>  U;
  register int i, j;
  const int    N = V.getSize() >> 1;
  Node*        n;
  int          J;

  // -- Create & print up list of unique Nodes.

  for (i = 0; i < N + N; i++) U.xadd (V[i]);

  S << "Mesh {" << endl;

  S << U.length() << "  Vertices" << endl;
  for (ListIterator<Node*> u(U); u.more(); u.next()) {
    n = u.current();
    S << *n << setw (10) << 0.0 << endl;
  }

  S << endl;
  
  // -- Print the original elements.

  S << nEl + nEl << "  Elements" << endl;
  for (i = 0; i < nEl; i++) {
    S << i + 1;
    J = E[i].getSize();
    S << setw (5) << J;
    for (j = 0; j < J; j++) S << setw(5) << E[i][j] -> ID();
    S << endl;
  }
  
  // -- Print up reflected elements.

  for (i = 0; i < nEl; i++) {
    S << nEl + i + 1;
    J = E[i + nEl].getSize();
    S << setw(5) << J;
    for (j = 0; j < J; j++)
      S << setw(5)
	<< ((A == R) ? E[i + nEl][j] -> ID() : E[i + nEl][J - j - 1] -> ID());
    S << endl;
  }

  // -- Print up boundary Nodes (this will need hand editing).

  S << endl << "2  Boundaries";

  S << endl;
  S << "1  1  ";
  j = 0;
  for (i = 0; i < N; i++) j += (V[i] -> interior()) ? 0 : 1;
  S << setw (5) << j;
  for (i = 0; i < N; i++)
    if (!V[i] -> interior()) S << setw (5) << V[i] -> ID();

  S << endl;
  S << "2  1  ";
  j = 0;
  for (i = N; i < N + N; i++) j += (V[i] -> interior()) ? 0 : 1;
  S << setw (5) << j;
  for (i = N + N - 1; i >= N; i--)
    if (!V[i] -> interior()) S << setw (5) << V[i] -> ID();  

  S << endl;
  S << endl << "0  Curves" << endl;
  
  S << "}" << endl;
}

