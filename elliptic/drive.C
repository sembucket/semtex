/*****************************************************************************
 * drive.C: compute solution to elliptic problem, compare to exact solution.
 *
 * Usage:
 * =====
 * elliptic [options] session
 *   options:
 *   -h       ... print usage prompt
 *   -v[v...] ... increase verbosity level
 *
 * Files:
 * =====
 * Set-up and execution of the problem is controlled by a session file,
 * which is block-structured.  Each block has a (pre-defined) name and
 * is delimited by '{' & '}'.  Contents of session files are case-insensitive,
 * *except* for function strings and named floating-point parameters.
 * Required blocks are: PROBLEM, PARAMETER, BOUNDARY and MESH.
 * Characters between block names and the opening '{' are ignored, and
 * can serve as comments.
 * 
 * PROBLEM block
 * -------------
 * This determines the type of problem to be solved and the kind of
 * geometry used (Cartesian, cylindrical, ...).
 *
 * (1)  Currently can solve Helmholtz (Possion, Laplace) and unsteady
 *      (Navier)--Stokes problems.
 * (2)  The default (and only current) geometry is 2D Cartesian.
 *
 * Example for a Helmholtz problem:
 *
 * problem {
 * helmholtz
 * [geometry   2D-Cartesian]
 * [forcing    function]
 * [exact      function]
 * }
 *
 * (1)  For a Helmholtz problem, the parameter "LAMBDA2" should be set in the
 *      floating point parameters section.  The default value (pre-installed)
 *      is zero, which converts the Helmholtz problem to a Poisson problem.
 * (2)  The forcing function describes (temporo-spatially-varying) forcing.
 *      If omitted, and LAMBDA2 == 0.0, the problem degenerates to Laplace.
 * (3)  Forcing function can use variables 'x' and 'y' for spatial variation
 *      and 't' for temporal variation.  Function should contain no spaces,
 *      which terminate scanning.
 *
 * PARAMETER block:
 * ---------------
 * In order, there should be option, integer and floating-point parameters.
 * The number of each type is given at the start of each sub-block.
 *
 * BOUNDARY block:
 * --------------
 * Boundary conditions for integer-tagged, non-overlapping boundary segments.
 * The number of segments is given at the start of the block.
 *
 * MESH block:
 * ----------
 * Commences with a list of vertices, followed by a list of elements with
 * corner vertices given as tags in vertex list.  The non-overlapping
 * boundary sectors are supplied as lists of vertices, together with an
 * integer tag which ties back to the BOUNDARY block.  Last comes
 * a list of curved edge specifiers, which are defined on a element-edge
 * basis.
 *
 * Program development by:
 * ======================
 * Hugh Blackburn
 * CSIRO
 * Division of Building, Construction and Engineering
 * P.O. Box 56
 * Highett, Vic 3190
 * Australia
 * hmb@dbce.csiro.au
 *
 *****************************************************************************/

static char  RCSid[] = "$Id$";

#include <Fem.h>
#include <new.h>

static char  prog[]  = "elliptic";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }
static const int MATCH = 0;

void         Helmholtz (Domain*, Mesh*, char*);
static void  getArgs   (int, char**, char*&);
static void  setUp     (ifstream&, char*&, char*&);


int  main (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);
#ifndef __DECCXX
  ios::sync_with_stdio ();
#endif

  ifstream  file;
  char*     session;
  char*     forcing = 0;
  char*     exact   = 0;

  initialize ();

  getArgs (argc, argv, session);

  file.open (session);
  setUp   (file, forcing, exact);

  Mesh*    M = preProcess (file);
  Domain*  D = new Domain (*M, session, iparam ("N_POLY"));

  file.close ();

  D -> openFiles ();

  Helmholtz (D, M, forcing);

  if (exact) D -> u[0] -> errors (exact);

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
    "  -v[v...]  ... increase verbosity level\n";
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      do
	setOption ("VERBOSE", option ("VERBOSE") + 1);
      while (*++argv[0] == 'v');
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
  
  if   (argc != 1) message (routine, "no session definition file", ERROR);
  else             session = *argv;
}


static void setUp (ifstream& file, char*& force, char*& exact)
// ---------------------------------------------------------------------------
// Return a function to solve, fill in optional strings for distributed 
// forcing and exact solution, if known.
// ---------------------------------------------------------------------------
{
  if (file.fail ()) message (prog, "couldn't open session file", ERROR);

  char   routine[] = "setProblem";
  char   s[StrMax], err[StrMax];

  seekBlock (file, "problem");

  file >> s;
  upperCase (s);

  if (strcmp (s, "HELMHOLTZ") == 0) {
    setIparam ("N_VAR",    SCALAR   );
    setOption ("PROBLEM",  HELMHOLTZ);
    setOption ("GEOMETRY", CART2D   );

    while (file >> s) {
      if (s[0] == '}') break;
      upperCase (s);

      if (strcmp (s, "GEOMETRY") == MATCH) {
	file >> s;
	upperCase (s);
	if (strstr (s, "2D") && strstr (s, "CART")) {
	  setOption ("GEOMETRY", CART2D);
	} else {
	  sprintf (err, "can't set geometry to: %s", s);
	  message (routine, err, ERROR);
	}

      } else if (strcmp (s, "FORCING") == MATCH) {
	force = new char [StrMax];
	file >> force;

      } else if (strcmp (s, "EXACT") == MATCH) {
	exact = new char [StrMax];
	file >> exact;
	
      } else {
	sprintf (err, "Helmholtz problem option? : %s", s);
	message (routine, s, ERROR);
      }
    }
    
    if (s[0] != '}') endBlock (file);
    
  } else {
    sprintf (err, "couldn't recognize a problem type in string: %s", s);
    message (routine, err, ERROR);
  }
}



