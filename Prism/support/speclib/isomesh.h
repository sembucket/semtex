#ifndef ISOMESH_H
#define ISOMESH_H

/* Isoparametric Mesh Generation
 *
 * $Id$
 *
 * Copyright (c) 1994 Ronald Dean Henderson
 *
 * Isoparametric Mesh Generation, prototypes and type declarations
 * ------------------------------------------------------------------------- */


/* Macro's and #define's */

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif
#
#define distance(p1,p2)   (sqrt(pow(p2.x-p1.x,2.) + pow(p2.y-p1.y,2.)))
#define genSetup(p1,p2,p,np,start,skip) \
  np = p->edge->np; start = p->edge->start; skip = p->edge->skip; \
  p1 = edgPoint (p->elmt, p->id); \
  p2 = edgPoint (p->elmt,(p->id + 1)%4);
#

/* --------     C U R V E D    S I D E    T Y P E S     --------- */

typedef struct lintype {     /* ........... Straight ............ */
  double   x0, y0        ;   /* Starting coordinates              */
  double   x1, y1        ;   /* Ending coordinates                */
} C_Strait;                  /* --------------------------------- */

typedef struct arctype {     /* .............. Arc .............. */
  double   xc, yc        ;   /* Coordinates of the center         */
  double   radius        ;   /* Radius of the curve               */
} C_Arc;                     /* --------------------------------- */

typedef struct partype {     /* .......... Parametric ........... */
  char *   function      ;   /* Parametric form c(t)              */
} C_Parametric;              /* --------------------------------- */

typedef struct wavtype {     /* ............. Wave .............. */
  double   wavenum       ;   /* Wave-number for the curve         */
  double   magnitude     ;   /* Magnitude of the wave             */
} C_Wave;                    /* --------------------------------- */

typedef struct filtype {     /* ............. File .............. */
  char *   name          ;   /* The name of the file to read      */
  double   xoffset       ;   /* Possible offsets of the           */
  double   yoffset       ;   /*    geometry for the current edge  */
} C_File;                    /* --------------------------------- */

typedef union curvinfo {     /* ....... Curved Side Infos ....... */
    C_Strait        line ;   /* straight line                     */
    C_Arc           arc  ;   /* an arc with specified curvature   */
    C_Parametric    map  ;   /* a parametric form                 */
    C_Wave          wave ;   /* a sinusoidal wave                 */
    C_File          file ;   /* coordinates from a file           */
} CurveInfo;                 /* --------------------------------- */

typedef enum {               /* ....... Curved Side Types ....... */
    T_Strait             ,   /* line                              */
    T_Arc                ,   /* arc                               */
    T_Parametric         ,   /* parametric function               */
    T_SinWave            ,   /* sin-wave                          */
    T_File                   /* points (spline fitted)            */
} CurveType;

typedef struct curve {       /* .... CURVED SIDE Definition ..... */
  CurveType    type [4]  ;   /* Array containing curve types      */
  CurveInfo    info [4]  ;   /* Array containing curve defs       */
} Curve;                     /* --------------------------------- */

typedef struct point {       /* A 2-D point */
  double     x, y     ;      /* coordinates */
} Point;

typedef struct vector {      /* A 2-D vector */
  double     x, y     ;      /* components   */
  double     length   ;      /* length       */
} Vector;

typedef struct geomf  {      /* Curve defined in a file */
  int           npts  ;      /* number of points        */
  int           pos   ;      /* last confirmed position */
  char         *name  ;      /* File/curve name         */
  double       *x, *y ;      /* coordinates             */
  double       *sx,*sy;      /* spline coefficients     */
  double       *arclen;      /* arclen along the curve  */
  struct geomf *next  ;      /* link to the next        */
} Geometry;

/* Prototypes */

void genxy  (Element *U, Curve *curv);
void geomap (Element *U);

#endif
