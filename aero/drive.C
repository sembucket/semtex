//////////////////////////////////////////////////////////////////////////////
// drive.C: control finite-element flow solver.
//
// Usage:
// =====
// aero [options] session
//   options:
//   -h       ... print usage prompt
//   -i[i]    ... use iterative solver for viscous [and pressure] steps
//   -v[v...] ... increase verbosity level
//   -chk     ... checkpoint field dumps
//
// Files:
// =====
// Set-up and execution of the problem is controlled by a session file,
// which is block-structured.  Each block has a (pre-defined) name and
// is delimited by '{' & '}'.  Contents of session files are case-insensitive,
////except* for function strings and named floating-point parameters.
// Required blocks are: PROBLEM, PARAMETER, BOUNDARY and MESH.
// Characters between block names and the opening '{' are ignored, and
// can serve as comments.
//
// (Body motion parameters are set up in a file named session.bdy.
// See comments in body.C for file structure.)
// 
// PROBLEM block
// -------------
// This determines the type of problem to be solved and the kind of
// geometry used (Cartesian, cylindrical, ...).
//
// (1)  Unsteady (Navier)--Stokes problems.
// (2)  The default (and only current) geometry is 2D Cartesian.
//
// Example for a Navier--Stokes problem:
//
// problem {
// navierstokes
// [geometry   2D-Cartesian]
// [forcing    function_1 .. function_DIM]
// }
//
// (1)  Unsteady Navier--Stokes problem.  An unsteady Stokes problem can also
//      be solved: at present this is done by conditional compilation (define
//      "STOKES" when compiling NS.C).
// (2)  The kinematic viscosity "KINVIS" should be set in the floating point
//      parameter section.  Default value, pre-installed, is unity.
// (3)  Forcing functions describe spatially-distributed force per unit mass
//      in the various components of the Navier--Stokes equations (in order).
//      DIM == 2 for 2D Cartesian.  Use 'x', 'y', and 't' as variables.
// (4)  Solution algorithm is the 'stiffly-stable' method.  This may be
//      inappropriate for low element orders.  
//
// PARAMETER block:
// ---------------
// In order, there should be option, integer and floating-point parameters.
// The number of each type is given at the start of each sub-block.
//
// BOUNDARY block:
// --------------
// Boundary conditions for integer-tagged, non-overlapping boundary segments.
// The number of segments is given at the start of the block.
//
// MESH block:
// ----------
// Commences with a list of vertices, followed by a list of elements with
// corner vertices given as tags in vertex list.  The non-overlapping
// boundary sectors are supplied as lists of vertices, together with an
// integer tag which ties back to the BOUNDARY block.  Last comes
// a list of curved edge specifiers, which are defined on a element-edge
// basis.
//
// Program development by:
// ======================
// Hugh Blackburn
// CSIRO
// Division of Building, Construction and Engineering
// P.O. Box 56
// Highett, Vic 3190
// Australia
// hmb@dbce.csiro.au
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <aero.h>
#include <new.h>


static char prog[]  = "aero";
static void memExhaust () { message ("new", "free store exhausted", ERROR); }

static void getArgs (int, char**, char*&);
static void setUp   (ifstream&);


int main (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#ifndef __DECCXX
  ios::sync_with_stdio ();
#endif

  ifstream*  input = new ifstream;
  char*      session;
  char       s[StrMax];

  // -- Initialization section.
  
  cout << prog << ": aeroelastic Navier--Stokes solver"  << endl;
  cout << "      (c) Hugh Blackburn 1995, 1996." << endl << endl;

  initialize    ();
  getArgs       (argc, argv, session);
  input -> open (session);
  setUp         (*input);

  // -- Get mesh information, save BC and parameter information.

  Mesh*  M = preProcess (*input);
  input -> close ();
  
  // -- Set up domain with single field, 'u'.

  Domain*  D = new Domain (*M, session, iparam ("N_POLY"));
  D -> u[0] -> setName ('u');

  // -- Add remaining velocity fields.

  const int DIM = iparam ("N_VAR");
  SystemField*  newField;
  for (int i = 1; i < DIM; i++) {
    newField = new SystemField (*D -> u[0]);
    D -> addField (newField);
    D -> u[i] -> setName ('u' + i);
  }

  // -- And constraint field 'p'.

  SystemField* Pressure = new SystemField (*D -> u[0], 1);
  D -> addField (Pressure);
  Pressure -> setName ('p');  
  PBCmanager::build (*Pressure);
  Pressure -> connect (*M, iparam ("N_POLY"));

  // -- Startup.

  D -> restart ();

  // -- Seek body information, construct body.

  input -> open (strcat (strcpy (s, session), ".bdy"));
  Body*  B = new Body (*input);
  input -> close ();

  B -> force   (*D);

  Analyser*  A = new Analyser (*D, *B);

  // -- Solve.

  NavierStokes (D, B, A);

  return EXIT_SUCCESS;
}


static void getArgs (int argc, char** argv, char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char routine[] = "getArgs";
  char buf[StrMax], c;
  char usage[]   =
    "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h        ... print this message\n"
    "  -i[i]     ... use iterative solver for viscous [& pressure] steps\n"
    "  -v[v...]  ... increase verbosity level\n"
    "  -chk      ... checkpoint field dumps\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      do
	setOption ("ITERATIVE", option ("ITERATIVE") + 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	setOption ("VERBOSE",   option ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv)) {
	setOption ("CHKPOINT", 1);
      } else {
	fprintf (stdout, usage, prog);
	exit (EXIT_FAILURE);	  
      }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if (argc != 1) {
    sprintf (buf, usage, prog);
    cout << buf;
    exit (EXIT_FAILURE);
  } else
    session = *argv;
}


static void setUp (ifstream& file)
// ---------------------------------------------------------------------------
// Return a function to solve, fill in optional function-string arrays.
// ---------------------------------------------------------------------------
{
  if (file.fail ()) message (prog, "couldn't open session file", ERROR);

  char  routine[] = "setProblem";
  char  s[StrMax], err[StrMax];

  seekBlock (file, "problem");

  file >> s;
  upperCase (s);

  if (strcmp (s, "NAVIERSTOKES") == 0) {

    setOption ("PROBLEM",  NAVIERSTOKES );
    setOption ("GEOMETRY", CART2D       );
    
    while (file >> s) {
      if (s[0] == '}') break;
      upperCase (s);

      if (strcmp (s, "GEOMETRY") == 0) {
	file >> s;
	upperCase (s);
	if (strstr (s, "2D") && strstr (s, "CART")) {
	  setIparam ("N_VAR",    TWO_COMPONENT);
	  setOption ("GEOMETRY", CART2D       );
	} else {
	  sprintf (err, "can't set geometry to: %s", s);
	  message (routine, err, ERROR);
	}
	
      } else if (strcmp (s, "FORCING") == 0) { // -- Read, but take no action.
	int  i, nstrng = iparam ("N_VAR");
	for (i = 0; i < nstrng; i++) file >> err;
	
      } else {
	sprintf (err, "Navier--Stokes problem option? : %s", s);
	message (routine, err, ERROR);
      }
    }
    
    if (s[0] != '}') endBlock (file);
    
  } else {
    sprintf (err, "couldn't recognize a problem type in string: %s", s);
    message (routine, err, ERROR);
  }
}
