/*
 * Curve implementation
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cubit.h"
#include "curve.h"
#include "element.h"
#include "edge.h"
#include "isomesh.h"

Curve *Curve_alloc (CurveType type, ...)
{
  Curve *curve = (Curve*) calloc(1,sizeof(Curve));
  va_list ap;

  va_start(ap, type);

  switch (type) {
  case Arc:
    curve->type = Arc;
    curve->info.arc.radius = va_arg(ap,double);
    break;

  case Spline:
    curve->type = Spline;
    curve->info.spline.x.b = va_arg(ap,double);
    curve->info.spline.y.b = va_arg(ap,double);
    curve->info.spline.x.c = va_arg(ap,double);
    curve->info.spline.y.c = va_arg(ap,double);
    break;

  case File:
    curve->type = File;
    curve->info.file = geofile_alloc(va_arg(ap,char*));
    break;

    /* Anthing else defaults to a straight edge */

  default:
    fprintf (stderr, "curve: type [%d] not implemented\n", (int) type);
  case Line:
    curve->type = Line;
    break;
  }

  va_end(ap);
  return curve;
}

void Curve_free (Curve *curve) 
{
  switch (curve->type) {
  case Arc:
    break;
  case Line:
    break;
  case Spline:
    break;
  case File:
    free (curve->info.file.name);
    free (curve->info.file.x);
    free (curve->info.file.sx);
    free (curve->info.file.y);
    free (curve->info.file.sy);
    free (curve->info.file.arclen);
    break;
  default:
    break;
  }

  free (curve);
}

/* Attach the curve to the edge of an element */

void Curve_attach (Curve *curve, Element *elmt, Edge *edge)
{
  if (edge->curve) {
    fprintf (stderr, "curve: edge is already attached!\n");
    return;
  } else
    edge->curve = curve;
  
  switch (curve->type) {
  case Arc:
    curve->info.arc = make_arc(elmt,edge,curve->info.arc.radius);
    break;

  case Spline: {
    double pos[4];
    pos[0] = curve->info.spline.x.b;
    pos[1] = curve->info.spline.y.b;
    pos[2] = curve->info.spline.x.c;
    pos[3] = curve->info.spline.y.c; 
    curve->info.spline = spline_alloc(elmt,edge,pos);
    break;
  }

  case File: {
    double xb[_MAX_NB];
    double yb[_MAX_NB];
    const int np = edge->np;
    _geofile_fit (&curve->info.file, elmt, edge, np, xb, yb);
    ecopy (np, xb, 1, *elmt->xmesh + edge->start, edge->skip);
    ecopy (np, yb, 1, *elmt->ymesh + edge->start, edge->skip);
    break;
  }

  default:
    break;
  }
}
