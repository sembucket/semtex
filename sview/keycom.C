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
  char help[] =
    "-- sview help menu --\n"
    "a                     : add default isosurface to end of storage list\n"
    "d                     : display default [0] isosurface\n"
    "d <n1> [<n2> .. <nx>] : display numbered isosurfaces\n"
    "f <c>                 : set field to c [Default: first field]\n"
    "l                     : list available isosurfaces\n"
    "m <val>               : make isosurface level val with current field\n";
  char c;

  cout << "-- Enter isosurface commands, h for usage prompt." << endl;
  cout << "> ";
  cin  >> c;

  switch (c) {
  case 'h':
    cout << help;
    break;
  default:
    break;
  }

  glutIdleFunc   (NULL);
  glutShowWindow ();
}
