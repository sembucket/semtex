///////////////////////////////////////////////////////////////////////////////
// graphics.C: All drawing is done through OpenGL commands and GLUT callbacks.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sview.h>

static void drawMesh  ();
static void polarView (GLdouble, GLdouble, GLdouble, GLdouble);


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
  case '+':
    State.radius *= 0.95;
    glutPostRedisplay();
    break;
  case '-':
    State.radius /= 0.95;
    glutPostRedisplay();
    break;
  case 'q':
    cerr << "sview: normal termination" << endl;
    exit (EXIT_SUCCESS);
  }
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
      State.xtrans -= 0.05 * State.xmax;
    break;
  case GLUT_KEY_RIGHT:
    if (State.rotate)
      State.yrot -= 5.0;
    else
      State.xtrans += 0.05 * State.xmax;
    break;
  case GLUT_KEY_UP:
    if (State.rotate)
      State.xrot -= 5.0;
    else
      State.ytrans += 0.05 * State.ymax;
    break;
  case GLUT_KEY_DOWN:
    if (State.rotate)
      State.xrot += 5.0;
    else
      State.ytrans -= 0.05 * State.ymax;
    break;
  case GLUT_KEY_PAGE_UP:
    if (State.rotate)
      State.zrot += 5.0;
    else
      State.ztrans += 0.05 * State.zmax;
    break;
  case GLUT_KEY_PAGE_DOWN:
    if (State.rotate)
      State.zrot -= 5.0;
    else
      State.ztrans -= 0.05 * State.zmax;
    break;
  default:
    return;
  }
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

  drawMesh        ();

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


void motion (int x,
	     int y)
// ---------------------------------------------------------------------------
// GLUT callback for mouse motion events within graphics window.
// ---------------------------------------------------------------------------
{
  glutPostRedisplay ();
}


void initGraphics ()
// ---------------------------------------------------------------------------
// Set up drawing defaults.
// ---------------------------------------------------------------------------
{
  char help[] = 
    "-- sview: isosurface viewer for spectral element meshes --\n"
    "          OpenGL version CSIRO 1999\n"
    "press ESC key within display window to enter surface manipulation mode\n"
    "press 'q' key within display window to exit sview\n";

  cout << help;

  glClearColor   (0.0, 0.0, 0.0, 0.0);
  glShadeModel   (GL_FLAT);
  glEnable       (GL_DEPTH_TEST);

  glMatrixMode   (GL_PROJECTION);
  gluPerspective (45.0, 1.0, 0.1 * State.length, 100 * State.length);

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
  register int i, k;
  const int    np    = Mesh -> np;
  const int    nz    = Mesh -> nz;
  const int    nel   = Mesh -> nel;
  const int    npm   = np - 1;
  const int    np2   = sqr (np);
  const int    np2m  = np2 - 1;
  const int    npnpm = np * (np - 1);
  const int    nzm   = nz - 1;
  float*       x;
  float*       y;
  float*       z = Mesh -> zmesh;

  glColor3f (0.7, 0.7, 0.7);

  for (k = 0; k < nel; k++) {

    x = Mesh -> xmesh + k * np2;
    y = Mesh -> ymesh + k * np2;

    // -- 1.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[i], y[i], z[0]);
    glEnd ();

    // -- 2.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[npm + i * np], y[npm + i * np], z[0]);
    glEnd ();
    
    // -- 3.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[np2m - i], y[np2m - i], z[0]);
    glEnd ();

    // -- 4.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[npnpm - i * np], y[npnpm - i * np], z[0]);
    glEnd ();

    // -- 5.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < nz; i++)
      glVertex3f (x[0], y[0], z[i]);
    glEnd ();

    // -- 6.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < nz; i++)
      glVertex3f (x[npm], y[npm], z[i]);
    glEnd ();

    // -- 7.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < nz; i++)
      glVertex3f (x[np2m], y[np2m], z[i]);
    glEnd ();

    // -- 8.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < nz; i++)
      glVertex3f (x[npnpm], y[npnpm], z[i]);
    glEnd ();

    // -- 9.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[i], y[i], z[nzm]);
    glEnd ();

    // -- 10.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[npm + i * np], y[npm + i * np], z[nzm]);
    glEnd ();
    
    // -- 11.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[np2m - i], y[np2m - i], z[nzm]);
    glEnd ();

    // -- 12.

    glBegin (GL_LINE_STRIP); 
    for (i = 0; i < np; i++)
      glVertex3f (x[npnpm - i * np], y[npnpm - i * np], z[nzm]);
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

