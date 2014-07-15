///////////////////////////////////////////////////////////////////////////////
// sview: Interactive display of isosurfaces.  This version works with
// NEKTON/PRISM/SEMTEX multi-element 3D data, and uses OpenGL and GLUT.
//
// Copyright (c) 1999 <--> $Date$, Hugh Blackburn
//
// USAGE
// -----
// sview [options] meshfile fieldfile
// meshfile  ... specifies name of NEKTON-format mesh file.
// fieldfile ... specifies name of NEKTON-format (binary) field file.
// options:
//   -b        ... start without mesh box
//   -c        ... set cylindrcal coordinates
//   -d        ... dump a TIFF image to file "sview.tif" on quitting
//   -h        ... print this message
//   -p <file> ... read particle locations from file
//   -s <file> ... read commands from file
//   -w        ... white background [Default: black]
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

#include "sview.h"

static char prog[] = "sview";
using namespace std;

// -- Global variables needed for graphics routines.

Flag         State;
Iso**        Surface;
Iso**        Display;
Sem*         Mesh;
Data*        Fields;
vector<Pnt*> Point;

static void getargs (int, char**, char*&, char*&, char*&, char*&);


int main (int    argc,
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
  int  i;
  char *mfile, *ffile, *script = 0, *pfile = 0;
  char start[] = 
    "-- sview : isosurface viewer for spectral element meshes --\n";

  // -- Set graphics state defaults.

  State.drawbox = GL_TRUE;
  State.drawiso = GL_FALSE;
  State.drawpar = GL_FALSE;
  State.rotate  = GL_TRUE;
  State.blackbk = GL_TRUE;
  State.alpha   = GL_FALSE;
  State.alias   = GL_FALSE;
  State.fog     = GL_FALSE;
  State.cylind  = GL_FALSE;
  State.dump    = GL_FALSE;

  State.xrot    = 0.0;
  State.yrot    = 0.0;
  State.zrot    = 0.0;
  State.xtrans  = 0.0;
  State.ytrans  = 0.0;
  State.ztrans  = 0.0;
  State.radius  = 1.0;
  State.wangle  = 45.0;

  // -- Initialise.

  glutInit (&argc, argv);
  glutInitWindowPosition(0, 0);

  getargs  (argc,  argv, mfile, ffile, pfile, script);

  cout << start;

  Surface = new Iso* [IsoMax];
  Display = new Iso* [IsoMax];
  for (i = 0; i < IsoMax; i++) {
    Surface[i] = 0;
    Display[i] = 0;
  }

  // -- Load external file information.

  Mesh   = loadMesh  (mfile);
  Fields = setFields (ffile);
  loadPnts (pfile);

  State.radius = 1.0 * State.length;

  // -- Set up windowing.

  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow    ("sview");

  // -- Define callbacks for events inside drawing window.

  glutReshapeFunc  (reshape);
  glutDisplayFunc  (display);
  glutKeyboardFunc (keyboard);
  glutSpecialFunc  (speckeys);
  glutIdleFunc     (NULL);

  initGraphics();

  if (script) processScript (script);

  // -- Transfer control to GLUT (no return).

  glutMainLoop();

  return EXIT_SUCCESS;
}


static void getargs (int    argc ,
		     char** argv ,
		     char*& mfile,
		     char*& ffile,
		     char*& pfile,
		     char*& sfile)
// ---------------------------------------------------------------------------
// Process command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = 
    "sview [options] meshfile fieldfile\n"
    "meshfile  ... specifies name of NEKTON-format mesh file.\n"
    "fieldfile ... specifies name of NEKTON-format (binary) field file.\n"
    "options:\n"
    "-b        ... start without mesh box\n"
    "-c        ... set cylindrical coordinates\n"
    "-d        ... dump a TIFF image to file \"sview.tif\" on quitting\n"
    "-h        ... print this message\n"
    "-p <file> ... specifies name of SEMTEX format point data file.\n"
    "-s <file> ... read commands from file\n"
    "-w        ... white background\n"
    "-g x y    ... set window size *** BUG FIX FOR PROBLEM IN GLUTINIT ***\n";

  char err[StrMax];
  int  x, y;

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'b':
      State.drawbox = GL_FALSE;
      break;
    case 'c':
      State.cylind = GL_TRUE;
      break;
    case 'd':
      State.dump = GL_TRUE;
      break;
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 's':
      --argc;
      sfile = *++argv;
      break;
    case 'p':
      --argc;
      pfile = *++argv;
      State.drawpar = GL_TRUE;
      break;
    case 'w':
      State.blackbk = !State.blackbk;
      break;

    case 'g':
      --argc; --argc;
      x = atoi (*++argv);
      y = atoi (*++argv);
      glutInitWindowSize(x, y);
      break;

    default:
      sprintf (err, "illegal option: %c\n", **argv);
      message (prog, err, ERROR);
      break;
    }

  if (argc == 2) {
    mfile = argv[0];
    ffile = argv[1];
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


void processScript (const char *name)
// ---------------------------------------------------------------------------
// Process sview commands from file name.
// ---------------------------------------------------------------------------
{
  const char routine[] = "processScript";
  char       buf[StrMax], command;
  ifstream   fp(name);

  if (!fp) {
    sprintf (buf, "couldn't open script file -- %s", name);
    message (routine, buf, WARNING);
  } else {
    sprintf (buf, "-- Reading commands from %s", name);
    message (routine, buf, REMARK);

    while (fp >> command) {
      fp.getline     (buf, StrMax);
      processCommand (command, buf);
    }

    fp.close();
  }
}


void quit ()
// ---------------------------------------------------------------------------
// This gets called on every legal exit path.
// ---------------------------------------------------------------------------
{
  if (State.dump) {
    writetiff ("sview.tif", "Isosurface", COMPRESSION_PACKBITS);
    cout << "Wrote file sview.tif" << endl;
  }
  cerr << "-- sview : normal termination" << endl;
  exit (EXIT_SUCCESS);
}