#ifndef SVIEW_H
#define SVIEW_H
///////////////////////////////////////////////////////////////////////////////
// Prototypes, classes and constants for sview.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

#include <GL/glut.h>

enum  lev    {WARNING, ERROR, REMARK};
const int    StrMax = 256;
const int    IsoMax = 8;
const double EPSSP  = 6.0e-7;
const double EPSDP  = 6.0e-14;
const double TWOPI  = 6.28318530717958647692;

template<class T> inline T sqr(T x)      { return x * x;            }
template<class T> inline T sgn(T x)      { return (x < 0) ? -1 : 1; }
template<class T> inline T min(T a, T b) { return (a < b) ?  a : b; }
template<class T> inline T max(T a, T b) { return (a > b) ?  a : b; }

typedef struct flag {		/* Global state variables.              */
  GLboolean drawbox;		/* Toggle drawing spectral elements.    */
  GLboolean rotate ;
  GLdouble  radius ;
  GLdouble  xrot   ;
  GLdouble  yrot   ;
  GLdouble  zrot   ;
  GLdouble  xtrans ;
  GLdouble  ytrans ;
  GLdouble  ztrans ;
  GLfloat   xmin   ;
  GLfloat   xmax   ;
  GLfloat   ymin   ;
  GLfloat   ymax   ;
  GLfloat   zmin   ;
  GLfloat   zmax   ;
  GLdouble  length ;
} Flag;

typedef struct sem {		/* Spectral element mesh information.   */
  int    np   ;			/* Number of points on element edge.    */
  int    nz   ;			/* Number of z planes.                  */
  int    nel  ;			/* Number of elements.                  */
  int    ntot ;			/* np*np*nel.                           */
  float* xmesh;			/* Row-major array of x-mesh locations. */
  float* ymesh;			/* Row-major array of y-mesh locations. */
  float* zmesh;			/* Z-mesh locations.                    */
} Sem;

typedef struct iso {		/* Isosurface wireframe information.    */
  int    id          ;		/* Identifier.                          */
  char   name[StrMax];		/* Description.                         */
  int    disp_indx   ;		/* Order for display rendering.         */
  int    nvert       ;		/* Number of triangle vertices.         */
  int    npoly       ;		/* Number of triangles.                 */
  int    npntspl     ;		/* Number of points per polygon (?3).   */
  int    normsgn     ;		/* Toggle for surface normals.          */
  float  xmin, xmax  ;		/* Limits in physical space.  X.        */
  float  ymin, ymax  ;		/* Y.                                   */
  float  zmin, zmax  ;		/* Z.                                   */
  float* pxyz        ;		/* Triangles.                           */
  int*   plist       ;		/* Some related descriptor?             */
} Iso;

/* -- Global variables needed for graphics routines. */

extern Flag  State;
extern Iso** Surface;
extern Sem*  Mesh;

/* -- Routines in main.C: */

void message (const char*, const char*, const lev&);

/* -- Routines in semIO.C: */

Sem* loadMesh (ifstream&);

/* -- Routines in keycom.C: */

void commandLine ();

/* -- Routines in graphics.C: */

void reshape      (int, int);
void display      ();
void keyboard     (unsigned char, int, int);
void speckeys     (int, int, int);
void motion       (int, int);
void initGraphics ();

#endif
