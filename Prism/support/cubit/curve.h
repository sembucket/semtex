#ifndef CURVE_H
#define CURVE_H

#include <stdarg.h>
#include "element.h"
#include "edge.h"


/* ------------------------------------------------------------------------- */

/* Basic curve types */

struct line {                 /* straight line segment */
  struct {
    double x, y;              /* endpoints */
  } pos[2];
};

struct arc {                  /* circular arc */
  double x, y;                /* center */
  double radius;              /* radius */
  struct {
    double start;             /* limits */
    double range;
  } theta;
};

typedef struct spline_t {    /* Bezier cubic spline */
  struct {
    double a, b, c, d;
  } x, y;
} spline_t;

typedef struct geofile_t {   /* geometry file */
  int     npts;              /* number of points */
  int     pos;               /* memory for searches along the curve */
  char   *name;              /* file name */
  double *x, *sx;            /* coordinates and spline coefficients */
  double *y, *sy;
  double *arclen;            /* arclength along the curve */
} geofile_t;

/* ------------------------------------------------------------------------- */

typedef enum CurveType {      /* ----------- Curve Types ------------- */
  Line       =  'L',          /* straight line                         */
  Arc        =  'A',          /* circular arc                          */
  Spline     =  'S',          /* Bezier cubic spline section           */
  File       =  'F'           /* coordinates from a file               */
} CurveType;
    
typedef struct curve {        /* ------------ Curved Side ------------ */
  CurveType    type;          /* Type of curve                         */
  union {
    struct line      line;    /* Straight line                         */
    struct arc       arc;     /* Circular arc                          */
    struct spline_t  spline;  /* Bezier cubic spline                   */
    struct geofile_t file;    /* Geometry file                         */
  } info;     
} Curve;

/* Prototypes */

Curve *Curve_alloc  (CurveType type, ...);
void   Curve_free   (Curve *curve);
void   Curve_attach (Curve *curve, Element *elmt, Edge *edge);

#endif
