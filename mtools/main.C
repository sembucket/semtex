///////////////////////////////////////////////////////////////////////////////
// main.C: Driver for qmesh: 2D quadrilateral unstructured mesh generator.
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
// The algorithm used closely follows that given in Ref. [1]. 
//
// -- Pseudocode description:
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
// 
// smooth mesh;
//
// print up;
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


// -- Global variables.

char prog[] = "qmesh";

real refcoeff = 0.0;

// -- Local routines.

static void getArgs (int, char**, int&, ifstream&);


int main (int argc, char** argv)
// ---------------------------------------------------------------------------
// Driver routine for mesh generator.
// ---------------------------------------------------------------------------
{
  ifstream infile;
  int      i, nsmooth = 0;

  getArgs (argc, argv, nsmooth, infile);

  Loop* L = new Loop;

  infile >> *L;
  infile.close ();

  if (graphics) { initGraphics ("x11"); drawBox (L); drawLoop (L);}

  L -> offset  ();
  if (graphics) drawLoop (L);

  L -> split   ();
  L -> connect ();

  for (i = 0; i < nsmooth; i++) {
    L -> smooth (); 
    if (graphics) L -> drawQuad ();
  }

  if (graphics) {
    eraseGraphics ();
    drawBox       (L);
    L -> drawQuad ();
    hardCopy      (L);
    stopGraphics  ();
  }

  L -> printMesh (cout);

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
	refcoeff = atof(*argv);
      else {
	--argc;
	refcoeff = atof(*++argv);
      }
      if (refcoeff < 0.0 || refcoeff > 1.0)
	error (prog, "refinement coefficient not in range 0--1", ERROR);
      if (refcoeff > 0.7)
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
