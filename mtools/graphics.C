///////////////////////////////////////////////////////////////////////////////
// graphics.cc
//
// All drawing is done using Super Mongo graphics package.
// Contact: Robert Lupton, rhl@astro.princeton.edu
//
// GRAPHICS needs to be defined during compilation otherwise all calls
// are to empty stubs.  This allows compilation on systems that don't
// have SM graphics installed.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <qmesh.h>

#ifdef GRAPHICS
  #include <sm_options.h>
  #include <sm_declare.h>
#endif


int graphics = 1;


void initGraphics (const char* device)
// ---------------------------------------------------------------------------
// Do whatever is needed to start up.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  char routine[] = "initGraphics";

  vector<float> pp (1);
  pp[0] = 41.0;

  if (sm_device ((char*) device))
    message (routine, "unable to initialize plotting device", ERROR);

  sm_graphics ();
  sm_location (3000, 31000, 3000, 31000);
  sm_defvar   ("TeX_strings", "1");
  sm_expand   (1.0);
  sm_ptype    (pp(), 1);
  sm_lweight  (1);
  sm_erase    ();
  sm_window   (1, 1, 1, 1, 1, 1);
#endif
}


void stopGraphics ()
// ---------------------------------------------------------------------------
// Shut down graphics.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  char a;

  sm_alpha ();
  cerr << "Press <return> to exit" << endl;
  sm_redraw (0);
  cin.get(a);
#endif
}


void eraseGraphics ()
// ---------------------------------------------------------------------------
// Clear window.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  sm_erase ();
#endif
}


void drawBox ()
// ---------------------------------------------------------------------------
// If the file "limits.sm" can be found, open it and attempt to read
// xmin xmax ymin ymax from it.  Otherwise, determine them from Node list.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  Point Pmax, Pmin, Centre, Range;
  float xmin, xmax, ymin, ymax;
  real  AR;

  ifstream file ("limits.sm");

  if (!file) {

    Global::limits (Pmin, Pmax);

    Centre.x = 0.5 * (Pmin.x + Pmax.x);
    Centre.y = 0.5 * (Pmin.y + Pmax.y);
    Range.x  = 0.5 * (Pmax.x - Pmin.x);
    Range.y  = 0.5 * (Pmax.y - Pmin.y);
    AR       = Range.y / Range.x;

    const real EXPAND = 1.1;

    if (AR >= 1.0) {
      xmin = Centre.x - EXPAND * AR * Range.x;
      xmax = Centre.x + EXPAND * AR * Range.x;
      ymin = Centre.y - EXPAND      * Range.y;
      ymax = Centre.y + EXPAND      * Range.y;
    } else {
      xmin = Centre.x - EXPAND      * Range.x;
      xmax = Centre.x + EXPAND      * Range.x;
      ymin = Centre.y - EXPAND / AR * Range.y;
      ymax = Centre.y + EXPAND / AR * Range.y;
    } 

  } else {
    file >> xmin >> xmax >> ymin >> ymax;
  }

  sm_limits (xmin, xmax, ymin, ymax);
  sm_box    (1, 2, 0, 0);
  sm_expand (1.6);
  sm_gflush ();
#endif
}


void drawLoop (const Loop* L      ,
	       const int&  numbers)
// ---------------------------------------------------------------------------
// Draw L.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  int npts;
  vector<float> x;
  vector<float> y;
  
  npts = L -> points (x, y);

  sm_conn (x(), y(), npts);
  sm_draw (x[0], y[0]);

  if (numbers) {
    register int i;
    char         label[StrMax];

    sm_expand (0.3);
    for (i = 0; i < npts; i++) {
      sm_relocate (x[i], y[i]);
      L -> nodeLabel (i, label);
      sm_label (label);
    }
    sm_expand (1.6);
  } else
    sm_points (x(), y(), npts);

  sm_gflush ();
#endif
}


void drawMesh (List<Quad*>& mesh,
	       const int    nums)
// ---------------------------------------------------------------------------
// Draw the mesh (with replicated mating edges --- the easy option!).
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  register int        i;
  float               x[4], y[4];
  ListIterator<Quad*> q (mesh);
  Quad*               Q;
  char                label[StrMax];

  for (; q.more(); q.next()) {
    Q = q.current();
    for (i = 0; i < 4; i++) {
      x[i] = Q -> vertex[i] -> pos().x;
      y[i] = Q -> vertex[i] -> pos().y;
    }

    sm_conn (x, y, 4);
    sm_draw (x[0], y[0]);

    if (nums) {
      sm_expand (0.3);
      for (i = 0; i < 4; i++) {
	sprintf     (label, "%1d", Q -> vertex[i] -> ID());
	sm_relocate (x[i], y[i]);
	sm_label    (label);
      }
      sm_expand (1.6);
    }
  }

  sm_gflush ();
#endif
}


void hardCopy (List<Quad*>& mesh)
// ---------------------------------------------------------------------------
// Draw m in PostScript file.
// ---------------------------------------------------------------------------
{
#ifdef GRAPHICS
  sm_gflush    ();
  sm_hardcopy  ();
  sm_alpha     ();

  initGraphics ("postfile mesh.eps");
  drawBox      ();
  drawMesh     (mesh, Global::verbose);

  sm_hardcopy  ();
  sm_alpha     ();
#endif
}


void qpause ()
// ---------------------------------------------------------------------------
// Wait for input.
// ---------------------------------------------------------------------------
{
   char a;
#ifdef GRAPHICS
  sm_alpha ();
#endif
  cerr << "Press <return> to continue" << endl;
#ifdef GRAPHICS
  sm_redraw (0);
#endif
  cin.get(a);
}

