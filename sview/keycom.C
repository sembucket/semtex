///////////////////////////////////////////////////////////////////////////////
// keycom.C: process commands typed by user into command window to generate
// isosurfaces.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>

#include <Sview.h>


void commandLine ()
// ---------------------------------------------------------------------------
// Read commands from command line in window from which programme
// was initiated.
// ---------------------------------------------------------------------------
{
  char routine[] ="commandLine";
  char help[] =
    "-- sview help menu --\n"
    "a                     : add default isosurface to end of storage list\n"
    "d                     : display default [0] isosurface\n"
    "d <n1> [<n2> .. <nx>] : display numbered isosurfaces\n"
    "f <c>                 : set field to c [Default: first field]\n"
    "l                     : list available isosurfaces\n"
    "m <val>               : make isosurface level val with current field\n"
    "n                     : invert normals of default isosurface\n";
  char  buf[StrMax], command;
  float value;

  cout << "-- Enter isosurface commands, h for usage prompt." << endl;
  cout << "> ";
  cin >> command;
  cin.getline (buf, StrMax);

  switch (command) {
  case 'h':
    cout << help;
    break;
  case 'd':
    if (!Surface[0]) {
      message (routine, "no default isosurface defined yet", WARNING);
      break;
    }
    State.drawiso = GL_TRUE;
    State.drawbox = GL_FALSE;

    glutPostRedisplay ();
    glutIdleFunc      (0);
    glutShowWindow    ();
    break;
  case 'm':
    istrstream (buf, strlen (buf)) >> value;
    Surface[0] = makeSurf (Mesh   -> nel  ,
			   Mesh   -> xgrid,
			   Mesh   -> ygrid,
			   Mesh   -> zgrid,
			   Fields -> elmt ,
			   Mesh   -> idim ,
			   Mesh   -> jdim ,
			   Mesh   -> kdim ,
			   value,    0    );
    break;
  case 'n':
    if (!Surface[0]) {
      message (routine, "no default isosurface defined yet", WARNING);
      break;
    }
    flipNorms (Surface[0]);
    break;
  default:
    break;
  }

}
