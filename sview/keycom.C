///////////////////////////////////////////////////////////////////////////////
// keycom.C: process commands typed by user into command window to generate
// isosurfaces.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>

static void catalogue ();
static int  addSurf   ();


void commandLine ()
// ---------------------------------------------------------------------------
// Read commands from command line in window from which programme
// was initiated.
// ---------------------------------------------------------------------------
{
  char buf[StrMax], command;

  cout << "(h = help) > ";
  cin >> command;
  cin.getline (buf, StrMax);

  processCommand (command, buf);
}


void processCommand (const char  command,
		     const char* buf    )
// ---------------------------------------------------------------------------
// Do command, with arguments from buf if appropriate.
// ---------------------------------------------------------------------------
{
  char routine[] ="processCommand";
  char help[] =
    "-- sview help menu --\n"
    "a                     : add default isosurface to end of storage\n"
    "d <n1> [<n2> .. <nx>] : display numbered isosurfaces [Default: 0]\n"
    "f <c>                 : set field to <c> [Default: first field]\n"
    "h                     : show this menu\n"
    "l                     : list available isosurfaces\n"
    "m <val>               : make isosurface level val with current field\n"
    "n <n>                 : invert normal of <n>th isosurface [Default: 0]\n"
    "o <r>                 : set radius of viewpoint to <r>\n"
    "p <x> <y> <z>         : specify coordinate rotation\n"
    "q                     : quit\n"
    "r <n>                 : remove/delete surface[<n>] from storage\n"
    "t <x> <y> <z>         : specify coordinate translation\n"
    "v <angle>             : angle of view, degrees (0 < v < 90)\n"
    "z <zoomfactor>        : specify magnification\n"
    "?                     : report current t/p/z state values\n";
  char   name;
  float  value;
  int    i, j, N;

  switch (command) {

  case 'a':
    addSurf ();
    break;

  case 'd':
    if (!Surface[0]) {
      message (routine, "no default isosurface defined yet", WARNING);
      break;
    }

    for (i = 0; i < IsoMax; i++) Display[i] = 0;

    if (!strlen (buf))
      Display[0] = Surface[0];	// -- Default.
    else {			// -- Parse list of indices.
      istringstream strm (buf);
      int           count = 0;

      while (strm >> i) {
	if (++count >= IsoMax) {
	  message (routine, "too many surfaces requested", WARNING);
	  break;
	} else if ((i+1) > countSurf (Surface)) {
	  message (routine, "requested surface not in store", WARNING);
	} else {
	  Display[count-1] = Surface[i];
	}
      }
    }

    State.drawiso = GL_TRUE;

    glutPostRedisplay ();
    glutIdleFunc      (0);
    glutShowWindow    ();
    break;

  case 'f':
    { istringstream ibuf (buf); ibuf >> name; loadData (Fields, name); }
    break;

  case 'h':
    cout << help;
    break;

  case 'l':
    catalogue ();
    break;

  case 'm':
    istringstream (buf) >> value;

    if (Surface[0]) {
      delete [] Surface[0] -> info;
      delete [] Surface[0] -> pxyz;
      delete [] Surface[0] -> nxyz;
      delete [] Surface[0] -> plist;
      delete    Surface[0];
      Surface[0] = 0;
    }

    Surface[0] = makeSurf (Mesh   -> nel  ,
			   Mesh   -> xgrid, Mesh -> ygrid, Mesh -> zgrid,
			   Fields -> elmt ,
			   Mesh   -> idim , Mesh -> jdim , Mesh -> kdim ,
			   Fields -> current,
			   value,    0      );
    break;

  case 'n':
    if (!Surface[0]) {
      message (routine, "no default isosurface defined yet", WARNING);
      break;
    }
    if (!strlen (buf)) flipNorms (Surface[0]);
    else {
      istringstream (buf) >> i;

      if (i >= 0 && i < (N = countSurf (Surface))) flipNorms (Surface[i]);
    }
    break;

  case 'p': {
    istringstream strm(buf);
    double val[3];

    for (i = 0; i < 3; i++) strm >> val[i];

    if (!strm)
      message (routine, "insufficient arguments to set angles", WARNING);
    else {
      State.xrot = val[0];
      State.yrot = val[1];
      State.zrot = val[2];
    }

    break;
  }

  case 'q':
    quit();
    break;

  case 'r':
    if (!strlen (buf)) {
      message (routine, "no index flagged for deletion", WARNING);
    } else {
      istringstream (buf) >> i;
      
      if (i == 0)
	message (routine, "can't delete default surface [0]", WARNING);

      else if (i < (N = countSurf (Surface))) {
	Iso* kill = Surface[i];
	for (j = i + 1; j < N; j++)
	  Surface[j - 1] = Surface[j];
	Surface[N - 1] = 0;

	free (kill -> info);
	free (kill -> pxyz);
	free (kill -> nxyz);
	free (kill -> plist);
	free (kill);

      } else
      	message (routine, "index flagged for deletion unavailable", WARNING);

      for (i = 1; i < IsoMax; i++) Display[i] = 0;
      Display[0] = Surface[0];

    }
    break;

  case 't': {
    istringstream strm(buf);
    float         val[3];

    for (i = 0; i < 3; i++) strm >> val[i];

    if (!strm)
      message (routine, "insufficient arguments to set positions", WARNING);
    else {
      State.xtrans = val[0];
      State.ytrans = val[1];
      State.ztrans = val[2];
    }

    break;
  }

  case 'v': {
    istringstream strm(buf);
    float         val;

    strm >> val;

    if (!strm)
      message (routine, "insufficient arguments to set view angle", WARNING);
    else if (val < 0.0 || val > 90.0)
      message (routine, "0 < angle < 90", WARNING);
    else
      State.wangle = val;
    break;
  }

  case 'z': {
    istringstream strm(buf);
    float         val;

    strm >> val;

    if (!strm)
      message (routine, "insufficient arguments to set radius", WARNING);
    else
      State.radius = (1./val) * State.length;
    break;
  }

  case '?': {
    cout << "-- Camera" << endl;
    cout << "   t " 
	 << State.xtrans << ' '
	 << State.ytrans << ' '
	 << State.ztrans << endl;
    cout << "   p " 
	 << State.xrot   << ' '
	 << State.yrot   << ' '
	 << State.zrot   << endl;
    cout << "   v "
	 << State.wangle << endl;
    cout << "   z "
	 << State.length / State.radius << endl;
      
    cout << "-- Surfaces" << endl;

    const int N = countSurf(Display);
    if (!N)
      cout << "   none" << endl;
    else {
      int i;
      for (i = 0; i < N; i++) {
	cout << "   " << Display[i]->info << endl;
      }
    }
    break;
  }

  default:
    break;
  }
}


int countSurf (Iso** list)
// ---------------------------------------------------------------------------
// How many surfaces (including the default, 0) are stored in global array?
// ---------------------------------------------------------------------------
{
  int i = 0;

  while (i < IsoMax && list[i]) i++;

  return i;
}


static void catalogue ()
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  int       i;
  const int N = countSurf (Surface);

  cout << "Available fields are: " << Fields -> name << endl;
  if (!N) 
    cout << "No surfaces defined" << endl;
  for (i = 0; i < N; i++)
    cout << "Surface[" << i << "]: " << Surface[i] -> info << endl;
}


static int addSurf ()
// ---------------------------------------------------------------------------
// Add default surface to end of Surface array.
// ---------------------------------------------------------------------------
{
  char routine[] = "addSurf";
  int       i = 0;
  const int N = countSurf (Surface);

  if (N == IsoMax)
    message (routine, "storage area is full", WARNING);

  if (N) Surface[N] = copySurf (Surface[0]);

  return i;
}


