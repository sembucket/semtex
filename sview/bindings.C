/*
 * Python bindings for sview
 */

#include <Python.h>
#include "Sview.h"

/* Provide external C linkage for our initialization routine */

extern "C" {
  void initSview();
}

// -- Global variables needed for graphics routines.

Flag  State;
Iso** Surface;
Iso** Display;
Sem*  Mesh;
Data* Fields;

/* ------------------------------------------------------------------------- */

static PyObject *sv_mainloop (PyObject *self, PyObject *args)
{
  // -- Transfer control to GLUT (no return).

  initGraphics ();
  glutMainLoop ();

  Py_INCREF(Py_None);
  return Py_None;
}

/* ------------------------------------------------------------------------- */

static PyObject *sv_meshin (PyObject *self, PyObject *args)
{
  char *file;

  if (!PyArg_ParseTuple(args,"s", &file))
    return 0;

  Mesh = loadMesh(file);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *sv_flowin (PyObject *self, PyObject *args)
{
  char *file;

  if (!PyArg_ParseTuple(args,"s", &file))
    return 0;

  Fields = setFields (file);

  Py_INCREF(Py_None);
  return Py_None;
}


/* ------------------------------------------------------------------------- */

static PyObject *sv_help (PyObject *self, PyObject *args)
{
  static char *help = 
    "-- sview help menu --\n"
    "a                     : add default isosurface to end of storage\n"
    "d <n1> [<n2> .. <nx>] : display numbered isosurfaces [Default: 0]\n"
    "f <c>                 : set field to <c> [Default: first field]\n"
    "h                     : show this menu\n"
    "l                     : list available isosurfaces\n"
    "m <val>               : make isosurface level val with current field\n"
    "n <n>                 : invert normal of <n>th isosurface [Default: 0]\n"
    "p <x> <y> <z>         : specify coordinate rotation\n"
    "q                     : quit\n"
    "r <n>                 : remove/delete surface[<n>] from storage\n"
    "t <x> <y> <z>         : specify coordinate translation\n"
    "z <zoomfactor>        : specify magnification\n"
    "?                     : report current t/p/z state values\n";
  fputs (help, stderr);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *noop (PyObject *self, PyObject *args)
{
  cerr << "----- not implemented" << endl;

  Py_INCREF(Py_None);
  return Py_None;
}

/* ------------------------------------------------------------------------- */


static PyMethodDef SviewMethods[] = {
  {"mainLoop",  sv_mainloop, METH_VARARGS},
  {"help",      sv_help,     METH_VARARGS},
  {"meshin",    sv_meshin,   METH_VARARGS},
  {"flowin",    sv_flowin,   METH_VARARGS},

  {"add",       noop,        METH_VARARGS},
  {"display",   noop,        METH_VARARGS},
  {"field",     noop,        METH_VARARGS},
  {"list",      noop,        METH_VARARGS},
  {"make",      noop,        METH_VARARGS},
  {"invert",    noop,        METH_VARARGS},
  {"remove",    noop,        METH_VARARGS},
  {0, 0, 0}
};


void initSview()
{
  PyObject *m = Py_InitModule("Sview", SviewMethods);

  int  i;
  char *mfile, *ffile;
  char start[] = 
    "-- sview : isosurface viewer for spectral element meshes --\n"
    "           OpenGL version CSIRO 1999\n";

  // -- Initialise.

  int argc = 1;
  char *argv[] = { "sview", NULL };

  glutInit (&argc, argv);

  cout << start;


  Surface = new Iso* [IsoMax];
  Display = new Iso* [IsoMax];
  for (i = 0; i < IsoMax; i++) {
    Surface[i] = 0;
    Display[i] = 0;
  }

  // -- Load external file information.

  mfile = "tg.mesh";
  ffile = "tg.fld";
#if 0
  Mesh   = loadMesh  (mfile);
  Fields = setFields (ffile);

  // -- Set graphics state defaults.

  State.drawbox = GL_TRUE;
  State.drawiso = GL_FALSE;
  State.rotate  = GL_TRUE;
  State.blackbk = GL_TRUE;
  State.noalias = GL_FALSE;

  State.xrot    = 0.0;
  State.yrot    = 0.0;
  State.zrot    = 0.0;
  State.xtrans  = 0.0;
  State.ytrans  = 0.0;
  State.ztrans  = 0.0;
  State.radius  = 1.0 * State.length;

  // -- Set up windowing.

  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow    ("sview");

  // -- Define callbacks for events inside drawing window.

  glutReshapeFunc  (reshape);
  glutDisplayFunc  (display);
  glutKeyboardFunc (keyboard);
  glutSpecialFunc  (speckeys);
  glutIdleFunc     (NULL);

  // -- Transfer control to GLUT (no return).

  initGraphics ();

  glutMainLoop ();
#endif
}

// ---------------------------------------------------------------------------
// Error message handler.
// ---------------------------------------------------------------------------

void message (const char* routine,
	      const char* text   ,
	      const lev&  level  )
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


