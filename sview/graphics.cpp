///////////////////////////////////////////////////////////////////////////////
// graphics.C: All drawing is done through OpenGL commands and GLUT callbacks.
//
// Copyright (c) 1999 <--> $Date$, Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "sview.h"

using namespace std;

// -- Materials definitions, in display order ---
//    blue, red, yellow, green, purple, white, turquoise, grey-blue

static GLfloat diffuse[IsoMax][4] = {
  {0.2,  0.2,  1.0,  1.0},
  {0.95, 0.2,  0.0,  1.0},
  {1.0,  1.0,  0.0,  1.0},
  {0.0,  0.9,  0.3,  1.0},
  {1.0,  1.0,  0.0,  1.0},
  {1.0,  1.0,  1.0,  1.0},
  {0.0,  0.9,  0.9,  1.0},
  {0.5,  0.7,  1.0,  1.0}
};

static GLfloat fogColor[4] = {0.5, 0.5, 0.5, 1.0};

static void drawMesh   ();
static void drawSurf   ();
static void drawPoints ();
static void skeleton   ();
static void polarView  (GLdouble, GLdouble, GLdouble, GLdouble);

static void drawSpecial();


void keyboard (unsigned char key,
	       int           x  ,
	       int           y  )
// ---------------------------------------------------------------------------
// GLUT callback for keyboard events within graphics window.
// ---------------------------------------------------------------------------
{
  switch (key) {
  case 27:			// -- ESC.
    glutIdleFunc     (commandLine);
    glutIconifyWindow();
    break;
  case 'q':
    quit();
    break;
  case '+':
    State.radius *= 0.95;
    skeleton();
    break;
  case '-':
  case '_':
    State.radius /= 0.95;
    skeleton();
    break;
  case 'a':
    if (State.alpha = !State.alpha) {
      glEnable    (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    } else {
      glDisable   (GL_BLEND);
    }
    break;
  case 'A':
    if (State.alias = !State.alias) {
      glDisable (GL_LINE_SMOOTH);
      glDisable (GL_POINT_SMOOTH);
      glDisable (GL_POLYGON_SMOOTH);
    } else {
      glEnable  (GL_LINE_SMOOTH);
      glEnable  (GL_POINT_SMOOTH);
      glEnable  (GL_POLYGON_SMOOTH);
    }
    break;
  case 'b':
    State.drawbox = !State.drawbox;
    break;
  case 'c':
    {
      ifstream table ("colours.sv");
      if (table) {
	cout << "-- Loading colour definitions from colours.sv [";
	int k = 0;
	double r, g, b, a;
	while (k < IsoMax && table >> r >> g >> b >> a) {
	  diffuse[k][0] = r;
	  diffuse[k][1] = g;
	  diffuse[k][2] = b;
	  diffuse[k][3] = a;
	  k++;
	}
	cout << k << "]" << endl;
      }
      table.close();
    }
    break;
  case 'd':
    //    writetiff ("sview.tif", "Isosurface", COMPRESSION_PACKBITS);
    writetiff ("sview.tif", "Isosurface", COMPRESSION_LZW);
    cout << "Wrote file sview.tif" << endl;
    break;
  case 'f': 
    if (State.fog = !State.fog) {
      glEnable (GL_FOG);
      glHint   (GL_FOG_HINT, GL_NICEST);
      glFogi   (GL_FOG_MODE, GL_EXP2);
      glFogfv  (GL_FOG_COLOR, fogColor);
      glFogf   (GL_FOG_DENSITY, 0.01);
    } else {
      glDisable (GL_FOG);
    }
    break;
  case 'i':
    if   (State.blackbk = !State.blackbk) glClearColor (0.0, 0.0, 0.0, 0.0);
    else                                  glClearColor (1.0, 1.0, 1.0, 0.0);
    break;
  case 'k':
    {
      ofstream knobs ("knobs.sv");
      knobs << "v"
	    << " " << State.wangle << endl;
      knobs << "z"
	    << " " << State.length / State.radius << endl;
      knobs << "p"
	    << " " << State.xrot
	    << " " << State.yrot
	    << " " << State.zrot   << endl;
      knobs << "t" 
	    << " " << State.xtrans
	    << " " << State.ytrans
	    << " " << State.ztrans << endl;
      knobs.close();
      cout << "Wrote file knobs.sv" << endl;
    }
    break;
  case 'l':
    processScript ("knobs.sv");
    break;
  case 'n':
    if (State.wangle > 0.0)  State.wangle = MAX (State.wangle - 5.0, 0.0);
    skeleton();
    break;
  case 'P':
    State.drawpar = !State.drawpar;
    break;
  case 'r':
    State.xrot    = 0.0;
    State.yrot    = 0.0;
    State.zrot    = 0.0;
    State.xtrans  = 0.0;
    State.ytrans  = 0.0;
    State.ztrans  = 0.0;
    State.radius  = 1.0 * State.length;
    break;
  case 'w':
    if (State.wangle < 90.0) State.wangle = MIN (State.wangle + 5.0, 90.0);
    skeleton();
    break;
  }
  State.drawiso = GL_TRUE;
  glutPostRedisplay();
}


void speckeys (int key,
	       int x  ,
	       int y  )
// ---------------------------------------------------------------------------
// GLUT callback for special key events within graphics window.
// ---------------------------------------------------------------------------
{
  switch (key) {
  case GLUT_KEY_HOME:
    State.rotate = !State.rotate;
    break;
  case GLUT_KEY_LEFT:
    if (State.rotate)
      State.yrot += 5.0;
    else
      State.xtrans -= 0.01 * State.length;
    skeleton();
    break;
  case GLUT_KEY_RIGHT:
    if (State.rotate)
      State.yrot -= 5.0;
    else
      State.xtrans += 0.01 * State.length;
    skeleton();
    break;
  case GLUT_KEY_UP:
    if (State.rotate)
      State.xrot -= 5.0;
    else
      State.ytrans += 0.01 * State.length;
    skeleton();
    break;
  case GLUT_KEY_DOWN:
    if (State.rotate)
      State.xrot += 5.0;
    else
      State.ytrans -= 0.01 * State.length;
    skeleton();
    break;
  case GLUT_KEY_PAGE_UP:
    if (State.rotate)
      State.zrot -= 5.0;
    else
      State.ztrans += 0.01 * State.length;
    skeleton();
    break;
  case GLUT_KEY_PAGE_DOWN:
    if (State.rotate)
      State.zrot += 5.0;
    else
      State.ztrans -= 0.01 * State.length;
    skeleton();
    break;
  default:
    return;
  }
  
  State.drawiso = GL_TRUE;
  glutPostRedisplay();
}


void display ()
// ---------------------------------------------------------------------------
// GLUT callback for display of graphics window.
// ---------------------------------------------------------------------------
{
  GLdouble AR =
    glutGet (GLUT_WINDOW_WIDTH)/ (GLdouble) glutGet (GLUT_WINDOW_HEIGHT);

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix ();
  polarView    (State.radius, State.zrot, State.xrot, State.yrot);
  glTranslated (State.xtrans, State.ytrans, State.ztrans);

  glMatrixMode   (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective (State.wangle, AR, 0.01 * State.length, 100 * State.length);
  glMatrixMode   (GL_MODELVIEW);

  if (State.drawbox) drawMesh ();
  if (State.drawpar) drawPoints();
  if (State.drawiso) drawSurf ();

  drawSpecial();

  glPopMatrix     ();
  glutSwapBuffers ();

  if (State.dump) quit();
}


void reshape (int w,
	      int h)
// ---------------------------------------------------------------------------
// GLUT callback for reshaping of graphics window.
// ---------------------------------------------------------------------------
{
  GLdouble AR = w / (GLdouble) h;
  glViewport (0, 0, w, h);

  glMatrixMode   (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective (State.wangle, AR, 0.01 * State.length, 100.0 * State.length);
  glMatrixMode   (GL_MODELVIEW);
}


void initGraphics ()
// ---------------------------------------------------------------------------
// Set up drawing defaults.
// ---------------------------------------------------------------------------
{
  // -- Attempt to load surface colour definition table from file.

  ifstream table ("colours.sv");

  if (table) {
    cout << "-- Loading colour definitions from colours.sv [";
    int    k = 0;
    double r, g, b, a;
    while (k < IsoMax && table >> r >> g >> b >> a) {
      diffuse[k][0] = r;
      diffuse[k][1] = g;
      diffuse[k][2] = b;
      diffuse[k][3] = a;
      k++;
    }
    cout << k << "]" << endl;
  }

  table.close();

  // -- Set window background colour and linewidth for skeleton.

  if   (State.blackbk) glClearColor (0.0, 0.0, 0.0, 0.0);
  else                 glClearColor (1.0, 1.0, 1.0, 0.0);

  glLineWidth (1.5);
  glPointSize (2.5);

  // -- Isosurface shading.

  glShadeModel   (GL_SMOOTH);
  glEnable       (GL_DEPTH_TEST);

  // -- Lighting.

  glEnable       (GL_LIGHTING);
  glEnable       (GL_LIGHT0);
  glLightModelf  (GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  
  // -- Blending (enabled by default).

  glEnable       (GL_BLEND);
  glBlendFunc    (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // -- Anti-aliasing (disabled by default);

  glDisable       (GL_LINE_SMOOTH);
  glDisable       (GL_POINT_SMOOTH);
  glDisable       (GL_POLYGON_SMOOTH);

  // -- Shared material properties:
  //    both surfaces have white highlights,
  //    but only "front" surface is shiny.

  GLfloat specular[] = {1.0, 1.0, 1.0, 1.0};

  //  glMaterialf    (GL_FRONT,          GL_SHININESS, 80.0); 
  glMaterialf    (GL_FRONT_AND_BACK,          GL_SHININESS, 80.0); 
  glMaterialfv   (GL_FRONT_AND_BACK, GL_SPECULAR,  specular);

  // -- Perspective.

  glMatrixMode   (GL_PROJECTION);
  gluPerspective (State.wangle, 1.0, 0.01 * State.length, 100 * State.length);
  glMatrixMode   (GL_MODELVIEW);

  // -- Issue prompt, also indicating initialisation is OK.

  cout <<
    "!  press ESC key in display window to enter surface manipulation mode\n"
    "!  press 'q' key in display window to exit sview" << endl;

}


// ===========================================================================
// Remaining routines have file scope only.
// ===========================================================================


static void drawMesh ()
// ---------------------------------------------------------------------------
// Draw outlines of spectral elements: render 12 sides of each element.
// ---------------------------------------------------------------------------
{
  register int i, j, k;
  const int    nel = Mesh -> nel;
  int          nr, ns, nt, nskip, ntot;
  float        *x, *y, *z;

  glDisable (GL_LIGHTING);
  glColor3f (0.7, 0.7, 0.7);

  for (k = 0; k < nel; k++) {

    x  = Mesh -> xgrid[k];
    y  = Mesh -> ygrid[k];
    z  = Mesh -> zgrid[k];

    nr = Mesh -> idim[k];
    ns = Mesh -> jdim[k];
    nt = Mesh -> kdim[k];
    
    nskip = nr * ns;
    ntot  = nskip * nt;

    // -- 1.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nr; i++)
      { j = i;                               glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 2.

    glBegin (GL_LINE_STRIP); for (i = 0; i < ns; i++)
      { j = i * nr;                          glVertex3f (x[j], y[j], z[j]); }
    glEnd ();
    
    // -- 3.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nr; i++)
      { j = nskip - i - 1;                   glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 4.

    glBegin (GL_LINE_STRIP); for (i = 0; i < ns; i++)
      { j = nr * (i + 1) - 1;                glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 5.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nt; i++)
      { j = i * nskip;                       glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 6.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nt; i++)
      { j = i * nskip + nr - 1;              glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 7.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nt; i++)
      { j = (i + 1) * nskip - 1;             glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 8.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nt; i++)
      { j = (i + 1) * nskip - nr;            glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 9.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nr; i++)
      { j = ntot - nskip + i;                glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 10.

    glBegin (GL_LINE_STRIP); for (i = 0; i < ns; i++)
      { j = ntot - nskip + i * nr;           glVertex3f (x[j], y[j], z[j]); }
    glEnd ();
    
    // -- 11.

    glBegin (GL_LINE_STRIP); for (i = 0; i < ns; i++)
      { j = ntot - nskip + (i + 1) * nr - 1; glVertex3f (x[j], y[j], z[j]); }
    glEnd ();

    // -- 12.

    glBegin (GL_LINE_STRIP); for (i = 0; i < nr; i++)
      { j = ntot - i - 1;                    glVertex3f (x[j], y[j], z[j]); }
    glEnd ();
  }
}


static void drawSurf ()
// ---------------------------------------------------------------------------
// Draw the isosurfaces selected for display.
// ---------------------------------------------------------------------------
{
  const int    N = countSurf (Display);
  register int i, v0, v1, v2;
  int          j, npoly;
  float        *norm, *vert;
  int          *pind;

  glEnable (GL_LIGHTING);  

  for (j = 0; j < N; j++) {

    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse [j]); 
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, diffuse [j]);

    npoly = Display[j] -> npoly;
    norm  = Display[j] -> nxyz;
    vert  = Display[j] -> pxyz;
    pind  = Display[j] -> plist;

    for (i = 0; i < npoly; i++) {
      v0 = 3 * pind[3 * i];
      v1 = 3 * pind[3 * i + 1];
      v2 = 3 * pind[3 * i + 2];
      glBegin (GL_POLYGON);
      glNormal3fv (norm + v0); glVertex3fv (vert + v0);
      glNormal3fv (norm + v1); glVertex3fv (vert + v1);
      glNormal3fv (norm + v2); glVertex3fv (vert + v2);
      glEnd ();
    }
  }
}


static void drawPoints ()
// ---------------------------------------------------------------------------
// Draw set of isolated points.
// ---------------------------------------------------------------------------
{
  const int    N = Point.size();
  register int i;
  Pnt*         datum;

  glDisable (GL_LIGHTING);  

  if   (State.blackbk) glColor3f (1.0, 1.0, 1.0);
  else                 glColor3f (0.0, 0.0, 0.0);
  
  glBegin (GL_POINTS);
  for (i = 0; i < N; i++) {
    datum = Point[i];
    glVertex3f (datum -> x, datum -> y, datum -> z);
  }
  glEnd();
}


static void polarView (GLdouble distance ,
		       GLdouble twist    ,
		       GLdouble elevation,
		       GLdouble azimuth  )
// ---------------------------------------------------------------------------
// Viewing transformation for spherical polar coordinates.
// ---------------------------------------------------------------------------
{
  glTranslated (0.0, 0.0, -distance);
  glRotated    (-twist,     0.0, 0.0, 1.0);
  glRotated    (-elevation, 1.0, 0.0, 0.0);
  glRotated    (azimuth,    0.0, 1.0, 0.0);
}


static void skeleton ()
// ---------------------------------------------------------------------------
// Draw the outline of spectral elements.
// ---------------------------------------------------------------------------
{
  glClear         (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix    ();
  polarView       (State.radius, State.zrot, State.xrot, State.yrot);
  glTranslated    (State.xtrans, State.ytrans, State.ztrans);
  drawMesh        ();
  glPopMatrix     ();
  glutSwapBuffers ();
}


static void drawSpecial ()
// ---------------------------------------------------------------------------
// Any special drawing code (to be hand edited).
// ---------------------------------------------------------------------------
{
#if 0
  int           i;
  const int     NC = 256;
  static float* coord = new float [2*NC];
  float         x, *yc = coord, *zc = yc + NC;

  glDisable (GL_LIGHTING);

  if   (State.blackbk) glColor3f (1.0, 1.0, 1.0);
  else                 glColor3f (0.0, 0.0, 0.0);

  for (i = 0; i < NC; i++) {
    yc[i] = cos (2.0*i*M_PI/NC);
    zc[i] = sin (2.0*i*M_PI/NC);
  }

  x = 1.25;
  
  glBegin (GL_LINE_STRIP);
  for (i = 0; i < NC; i++)
   glVertex3f (x, yc[i], zc[i]); 
  glVertex3f (x, yc[0], zc[0]);
  glEnd();

  x = -1.25;
  
  glBegin (GL_LINE_STRIP);
  for (i = 0; i < NC; i++)
   glVertex3f (x, yc[i], zc[i]); 
  glVertex3f (x, yc[0], zc[0]);
  glEnd();
#endif
}