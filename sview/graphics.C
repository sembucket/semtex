///////////////////////////////////////////////////////////////////////////////
// graphics.C: All drawing is done through OpenGL commands and GLUT callbacks.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>

static void drawMesh  ();
static void drawSurf  ();
static void skeleton  ();
static void polarView (GLdouble, GLdouble, GLdouble, GLdouble);
static void initLight ();


void keyboard (unsigned char key,
	       int           x  ,
	       int           y  )
// ---------------------------------------------------------------------------
// GLUT callback for keyboard events within graphics window.
// ---------------------------------------------------------------------------
{
  switch (key) {
  case 27:
    glutIdleFunc     (commandLine);
    glutIconifyWindow();
    break;
  case 'q':
    cerr << "-- sview : normal termination" << endl;
    exit (EXIT_SUCCESS);
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
  case 'b':
    State.drawbox = !State.drawbox;
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
  case GLUT_KEY_INSERT:
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
  glClear         (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix    ();
  polarView       (State.radius, State.zrot, State.xrot, State.yrot);
  glTranslated    (State.xtrans, State.ytrans, State.ztrans);

  if (State.drawbox)               drawMesh ();
  if (State.drawiso && Surface[0]) drawSurf ();

  glPopMatrix     ();

  glutSwapBuffers ();
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
  gluPerspective (45.0, AR, 0.1 * State.length, 100.0 * State.length);
  glMatrixMode   (GL_MODELVIEW);
}


void initGraphics ()
// ---------------------------------------------------------------------------
// Set up drawing defaults.
// ---------------------------------------------------------------------------
{
  char help[] =
    "   press ESC key in display window to enter surface manipulation mode\n"
    "   press 'q' key in display window to exit sview\n";
  
  cout << help;

  glClearColor   (0.0, 0.0, 0.0, 0.0);
  glShadeModel   (GL_SMOOTH);
  glEnable       (GL_DEPTH_TEST);

  initLight();

  glMatrixMode   (GL_PROJECTION);
  gluPerspective (45.0, 1.0, 0.0 * State.length, 100 * State.length);
  glMatrixMode   (GL_MODELVIEW);
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
  register int i, v0, v1, v2;
  const int    npoly = Surface[0] -> npoly;
  float*       norm  = Surface[0] -> nxyz;
  float*       vert  = Surface[0] -> pxyz;
  int*         pind  = Surface[0] -> plist;
  
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


static void initLight ()
// ---------------------------------------------------------------------------
// Initialize lighting and materials.
// ---------------------------------------------------------------------------
{
  GLfloat ambient[]             = {0.1, 0.1,   0.1,  1.0};
  GLfloat diffuse[]             = {0.5, 1.0,   1.0,  1.0};
  GLfloat position0[]           = {0.0, 0.0,  20.0,  0.0};
  GLfloat position1[]           = {0.0, 0.0, -20.0,  0.0};
  GLfloat front_mat_shininess[] = {60.0};
  GLfloat front_mat_specular[]  = {0.2, 0.2,   0.2,  1.0};
  GLfloat front_mat_diffuse[]   = {0.5, 0.28,  0.38, 1.0};
  GLfloat lmodel_ambient[]      = {1.0, 1.0,   1.0,  1.0};
  GLfloat lmodel_twoside[]      = {GL_FALSE};

  glLightfv      (GL_LIGHT0, GL_AMBIENT,  ambient);
  glLightfv      (GL_LIGHT0, GL_DIFFUSE,  diffuse);
  glLightfv      (GL_LIGHT0, GL_POSITION, position0);
  glEnable       (GL_LIGHT0);

  glLightfv      (GL_LIGHT1, GL_AMBIENT,  ambient);
  glLightfv      (GL_LIGHT1, GL_DIFFUSE,  diffuse);
  glLightfv      (GL_LIGHT1, GL_POSITION, position1);
  glEnable       (GL_LIGHT1);

  glLightModelfv (GL_LIGHT_MODEL_AMBIENT,  lmodel_ambient);
  glLightModelfv (GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
  glEnable       (GL_LIGHTING);

  glMaterialfv   (GL_FRONT_AND_BACK, GL_SHININESS, front_mat_shininess);
  glMaterialfv   (GL_FRONT_AND_BACK, GL_SPECULAR,  front_mat_specular);
  glMaterialfv   (GL_FRONT_AND_BACK, GL_DIFFUSE,   front_mat_diffuse);
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
