#ifndef SVIEW_H
#define SVIEW_H
///////////////////////////////////////////////////////////////////////////////
// Prototypes, classes and constants for sview.
//
// Copyright (C) 1999 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <iostream.h>
#include <strstream.h>
#include <fstream.h>

#include <GL/glut.h>

#include <tiffio.h>		/* Sam Leffler's libtiff library. */

enum  lev    {WARNING, ERROR, REMARK};
const int    StrMax = 256;
const int    IsoMax = 8;
const int    FldMax = 32;

template<class T> inline T sqr(T x)      { return x * x;            }
template<class T> inline T sgn(T x)      { return (x < 0) ? -1 : 1; }
template<class T> inline T min(T a, T b) { return (a < b) ?  a : b; }
template<class T> inline T max(T a, T b) { return (a > b) ?  a : b; }

typedef struct flag {		/* Global state variables.              */
  GLboolean drawbox;		/* Toggle drawing spectral elements.    */
  GLboolean drawiso;		/* Toggle drawing enabled isosurfaces   */
  GLboolean rotate ;		/* Toggle rotation/translation of view. */
  GLboolean blackbk;		/* Toggle black/white background.       */
  GLboolean noalias;		/* Toggle antaliasing of polygons.      */
  GLboolean cylind ;		/* True for cylindrical coordinates.    */
  GLboolean dump   ;		/* Flag dump of TIFF image.             */
  GLdouble  radius ;		/* Positioning of viewing point.        */
  GLdouble  xrot   ;		/*                                      */
  GLdouble  yrot   ;		/*                                      */
  GLdouble  zrot   ;		/*                                      */
  GLdouble  xtrans ;		/*                                      */
  GLdouble  ytrans ;		/*                                      */
  GLdouble  ztrans ;		/*                                      */
  GLfloat   xmin   ;		/* Mesh extents.                        */
  GLfloat   xmax   ;		/*                                      */
  GLfloat   ymin   ;		/*                                      */
  GLfloat   ymax   ;		/*                                      */
  GLfloat   zmin   ;		/*                                      */
  GLfloat   zmax   ;		/*                                      */
  GLdouble  length ;		/* A mesh length scale from extents.    */
} Flag;

typedef struct sem {		/* Spectral element mesh information.   */
  int     nel  ;		/* Number of elements.                  */
  int     nrep ;		/* Number of periodic z-repetitions.    */
  int*    idim ;		/* "r" dim for each element (equal!).   */
  int*    jdim ;		/* "s" dim for each element (equal!).   */
  int*    kdim ;		/* "t" dim for each element (equal!).   */
  float** xgrid;		/* x locations for each element.        */
  float** ygrid;		/* y    ""     ""   ""     ""           */
  float** zgrid;		/* z    ""     ""   ""     ""           */
} Sem;

typedef struct data {		/* Spectral element data information.   */
  int       nfields        ;	/* Number of fields in file.            */
  char      name[StrMax]   ;	/* Their one-character names.           */
  ifstream  file           ;	/* File handle.                         */
  streampos fldPosn[FldMax];	/* Locations of the various fields.     */
  char      current	   ;	/* Name of current field.               */
  float**   elmt           ;	/* Equivalent storage, by element.      */
} Data;

typedef struct iso {		/* Isosurface wireframe information.    */
  char*  info;			/* String that describes surface.       */
  int    nvert;			/* Number of triangle vertices.         */
  int    npoly;			/* Number of triangles.                 */
  float* pxyz ;			/* Vertices.                            */
  float* nxyz ;			/* Normals at vertices.                 */
  int*   plist;			/* Vertex indices for each triangle.    */
} Iso;

/* -- Global variables needed for graphics routines. */

extern Flag  State;		/* Local (non-OpenGL/GLUT) state variables. */
extern Iso** Surface;		/* Array of stored isosurfaces.             */
extern Iso** Display;		/* Array of isosurfaces chosen for display. */
extern Sem*  Mesh;		/* Element nodal location data.             */
extern Data* Fields;		/* Scalar field data structure/retrieval.   */

/* -- External routines in main.C: */

void message       (const char*, const char*, const lev&);
void processScript (const char*);
void quit          ();

/* -- Routines in semIO.C: */

Sem*  loadMesh  (const char*);
Data* setFields (const char*);
int   loadData  (Data*, char);

/* -- Routines in keycom.C: */

void commandLine    ();
void processCommand (const char, const char*);
int  countSurf      (Iso**);

/* -- Routines in graphics.C: */

void reshape      (int, int);
void display      ();
void keyboard     (unsigned char, int, int);
void speckeys     (int, int, int);
void initGraphics ();

/* -- Routines in isowire.C: */

Iso* makeSurf  (int, float**, float**, float**, float**,
		int*, int*, int*, char, float, int);
Iso* copySurf  (Iso*);
void flipNorms (Iso*);

/* -- Routines in image.C: */

int writetiff (char*, char*, int);

#endif
