///////////////////////////////////////////////////////////////////////////////
// rectmesh.C: create a session file for a rectangular mesh of elements.
//
// Usage:
// -----
// rectmesh [file]
//
// Files:
// -----
// Input consists of a list of x, followed by y, locations of element
// boundaries, one per line.  A blank line separates x from y
// locations.  Output consists of a (2D) session file with an element
// order of 7, and "wall" group boundaries around the domain border.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <sem_h>

static char prog[] = "rectmesh";
static void getargs (int, char**, istream*&);
static void header  ();


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char                    line[STR_MAX];
  istream*                input;
  real                    x, y;
  stack <real>            X, Y;
  vector<vector<Point*> > vertex;
  int                     Nx = 0, Ny = 0;
  int                     i, j, k;

  getargs (argc, argv, input);

  // -- Read x, then y locations onto two stacks.

  while (input -> getline(line, STR_MAX).gcount() > 1) {
    istrstream (line, strlen(line)) >> x;
    X.push (x);
    Nx++;
  }

  while (input -> getline(line, STR_MAX)) {
    istrstream (line, strlen(line)) >> y;
    Y.push (y);
    Ny++;
  }
  
  // -- Insert into vertex matrix.

  vertex.resize (Ny);
  for (i = 0; i < Ny; i++) vertex[i].resize(Nx);
  for (i = 0; i < Ny; i++)
    for (j = 0; j < Nx; j++) {
      vertex[i][j] = new Point;
      vertex[i][j] -> z = 0.0;
    }

  j = Nx;
  while (j--) { vertex[0][j] -> x = X.top(); X.pop(); }
  i = Ny;
  while (i--) { vertex[i][0] -> y = Y.top(); Y.pop(); }

  for (i = 0; i < Ny; i++)
    for (j = 0; j < Nx; j++) {
      vertex[i][j] -> x = vertex[0][j] -> x;
      vertex[i][j] -> y = vertex[i][0] -> y;
    }

  header();

  // -- Print up vertex list.

  cout << "<NODES NUMBER=" << Nx*Ny << ">" << endl;

  for (k = 0, i = 0; i < Ny; i++)
    for (j = 0; j < Nx; j++)
      cout << setw(5)  << ++k << "\t"
	   << setw(15) << vertex[i][j] -> x
	   << setw(15) << vertex[i][j] -> y
	   << setw(15) << vertex[i][j] -> z
	   << endl;
  
  cout << "</NODES>" << endl;

  // -- Print up elements.

  cout << endl << "<ELEMENTS NUMBER=" << (Nx-1)*(Ny-1) << ">" << endl;

  for (k = 0, i = 0; i < (Ny - 1); i++)
    for (j = 0; j < (Nx - 1); j++)
      cout << setw(5) << ++k << "\t" << "<Q>"
	   << setw(5) << j +  i      * Nx + 1
	   << setw(5) << j +  i      * Nx + 2
	   << setw(5) << j + (i + 1) * Nx + 2
	   << setw(5) << j + (i + 1) * Nx + 1
	   << "    </Q>" << endl;
    
  cout << "</ELEMENTS>" << endl;

  // -- Print up surfaces.

  cout << endl << "<SURFACES NUMBER=" << 2*((Nx-1)+(Ny-1)) << ">" << endl;

  for (k = 0, j = 0; j < (Nx - 1); j++)
    cout << setw(5) << ++k << setw(5) 
	 <<  j + 1
	 << "    1"
	 << "    <B> w </B>" << endl;

  for (i = 0; i < (Ny - 1); i++)
    cout << setw(5) << ++k 
	 << setw(5) << (i + 1) * (Nx - 1)
	 << "    2"
	 << "    <B> w </B>" << endl;

  for (j = Nx - 1; j > 0; j--)
    cout << setw(5) << ++k 
	 << setw(5) <<  j + (Nx - 1) * (Ny - 2)
	 << "    3"
	 << "    <B> w </B>" << endl;

  for (i = Ny - 1; i > 0; i--)
    cout << setw(5) << ++k 
	 << setw(5) << (i - 1) * (Nx - 1) + 1   
	 << "    4"
	 << "    <B> w </B>" << endl;
    
  cout << "</SURFACES>" << endl;

  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: rectmesh [options] [file]\n"
    "  options:\n"
    "  -h ... print this message\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open geometry file", ERROR);
  } else input = &cin;
}


static void header ()
// ---------------------------------------------------------------------------
// Output header information.  Mesh is valid (eg for meshpr) without this.
// ---------------------------------------------------------------------------
{
  cout << "<FIELDS>\n\tu\tv\tp\n</FIELDS>" << endl << endl;

  cout << "<GROUPS NUMBER=1>\n\t1\tw\twall\n</GROUPS>" << endl << endl;

  cout << "<BCS NUMBER=1>\n\t1\tw\t3" << endl;
  cout << "\t\t\t<D>\tu = 0.0\t</D>"  << endl;
  cout << "\t\t\t<D>\tv = 0.0\t</D>"  << endl;
  cout << "\t\t\t<H>\tp = 0.0\t</H>"  << endl;
  cout << "</BCS>" << endl << endl;

  cout << "<TOKENS>" << endl;
  cout << "\tFFX       = 0.0"       << endl;
  cout << "\tRNG       = 0"         << endl;
  cout << "\tC_SMAG    = 0.1"       << endl;
  cout << "\tKINVIS    = 2e-6"      << endl;
  cout << "\tREFVIS    = 10*KINVIS" << endl;        
  cout << "\tD_T       = 0.005"     << endl;
  cout << "\tN_STEP    = 100"       << endl;
  cout << "\tN_TIME    = 2"         << endl;
  cout << "\tN_POLY    = 7"         << endl;
  cout << "\tN_Z       = 1"         << endl;
  cout << "\tBETA      = 1.0"       << endl;
  cout << "\tIO_CFL    = 50"        << endl;
  cout << "\tIO_FLD    = 1000"      << endl;
  cout << "\tAVERAGE   = 0"         << endl;
  cout << "\tCHKPOINT  = 1"         << endl;
  cout << "\tITERATIVE = 1"         << endl;
  cout << "</TOKENS>" << endl << endl;
}
