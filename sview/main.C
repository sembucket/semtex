///////////////////////////////////////////////////////////////////////////////
//
// SYNOPSIS
// --------
// sview: Interactive display of isosurfaces.  This version works with
// NEKTON/PRISM/SEMTEX multi-element 3D data, and uses OpenGL and GLUT.
//
// USAGE
// -----
// sview [options] meshfile fieldfile
// meshfile  ... specifies name of NEKTON-format mesh file.
// fieldfile ... specifies name of NEKTON-format (binary) field file.
// options:
//   -h      ... print this message
// 
//
// REFERENCES
// ----------
// [1] Lorensen & Cline (1987), paper in SIGGRAPH 87.
//
// [2] Neider, J., Davis, T. & Woo, M. (1993), OpenGL Programming Guide,
//     Addison--Wesley.
//
// [3] Kilgard, M.J. (1996), OpenGL Programming for the X Window System,
//     Addison--Wesley.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>

static char  prog[] = "sview";

// -- Global variables needed for graphics routines.

Flag  State;
Iso** Surface;
Sem*  Mesh;

static void getargs  (int, char**, ifstream&, ifstream&);


void main (int    argc,
	   char** argv)
// ---------------------------------------------------------------------------
// Driver.
//
// There are two levels of user interaction.  When GLUT has control of
// the drawing window, interaction occurs through mouse motion events
// and keystrokes in the window.  When ESC is pressed, the drawing
// window is iconised and the user types commands for isosurface
// generation into the window from which sview was initialised.
// Control is returned to GLUT by a request to draw isosurfaces.
//
// ---------------------------------------------------------------------------
{
  ifstream mfile, ffile;

  // -- Initialise.

  glutInit (&argc, argv);
  getargs  (argc,  argv, mfile, ffile);

  State.drawbox = GL_TRUE;
  State.rotate  = GL_TRUE;
  State.xrot    = 0.0;
  State.yrot    = 0.0;
  State.zrot    = 0.0;
  State.xtrans  = 0.0;
  State.ytrans  = 0.0;
  State.ztrans  = 0.0;

  // -- Load external file information.

  Mesh = loadMesh (mfile);
  
  State.radius = 4.0 * State.length;

  // -- Set up windowing.

  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow    ("sview");

  // -- Define callbacks for events inside drawing window.

  glutReshapeFunc  (reshape);
  glutDisplayFunc  (display);
  glutKeyboardFunc (keyboard);
  glutSpecialFunc  (speckeys);
  //  glutMotionFunc   (motion);
  glutIdleFunc     (NULL);

  // -- Transfer control to GLUT (no return).

  initGraphics ();

  glutMainLoop ();
}


static void getargs (int       argc ,
		     char**    argv ,
		     ifstream& mfile,
		     ifstream& ffile)
// ---------------------------------------------------------------------------
// Process command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = 
    "sview [options] meshfile fieldfile\n"
    "meshfile  ... specifies name of NEKTON-format mesh file.\n"
    "fieldfile ... specifies name of NEKTON-format (binary) field file.\n"
    "options:\n"
    "-h        ... print this message\n";
  char err[StrMax];

  while (--argc && **++argv == '-')
    switch (**++argv) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (err, "illegal option: %c\n", **argv);
      message (prog, err, ERROR);
      break;
    }

  if (argc == 2) {
    mfile.open (argv[0], ios::in);
    ffile.open (argv[1], ios::in);
    if (!mfile) message (prog, "couldn't open mesh file",  ERROR);
    if (!ffile) message (prog, "couldn't open field file", ERROR);
  } else {
    cerr << usage;
    exit (EXIT_FAILURE);
  }    
}


void message (const char* routine,
	      const char* text   ,
	      const lev&  level  )
// ---------------------------------------------------------------------------
// Error message handler.
// ---------------------------------------------------------------------------
{
  switch (level) {
  case WARNING:
    cerr << "WARNING: " << routine << ": " << text << endl;
    return;
  case REMARK:
    cerr << text << endl;
    return;
  case ERROR:
    cerr << "ERROR: " << routine << ": " << text << endl;
    exit (EXIT_FAILURE);
    break;
  }
}
