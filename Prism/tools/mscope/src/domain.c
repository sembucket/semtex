/*
 * Some functions for manipulating the computational domain
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <string.h>
#include <math.h>

#include "cubit/cubit.h"
#include "mscope.h"

extern Domain Geometry;

/* ------------------------------------------------------------------------- */

int Domain_init (const char *name) 
{
  Domain_reset();

  Geometry.name     = strdup(name);
  Geometry.param    = param_alloc(32);
  Geometry.force    = NULL;
  Geometry.solution = NULL;
  Geometry.ic       = NULL;
  Geometry.user     = NULL;
  Geometry.fields   = NULL;
  Geometry.history  = NULL;
  Geometry.mesh     = NULL;
  Geometry.BC       = NULL;
  Geometry.A        = NULL;

  return 0;
}

int Domain_reset()
{
  int i;

  if (Geometry.name)
    free (Geometry.name);
  if (Geometry.param)
    param_free (Geometry.param);
  if (Geometry.force)
    keyword_free (Geometry.force);
  if (Geometry.solution)
    keyword_free (Geometry.solution);
  if (Geometry.ic)
    keyword_free (Geometry.ic);
  if (Geometry.history)
    keyword_free (Geometry.history);
  if (Geometry.fields)
    keyword_free (Geometry.fields);
  if (Geometry.user)
    keyword_free (Geometry.user);

  if (Geometry.mesh)
    Mesh_free (Geometry.mesh);
  if (Geometry.BC)
    BC_free (Geometry.BC);
  if (Geometry.A)
    Matrix_free (Geometry.A);

  for (i = 0; i < Geometry.nfields; i++)
    Field_free (Geometry.solVector[i]);

  memset(&Geometry, '\0', sizeof(Domain));

  /* Reset certain "special" parameters */

  dparam_set ("XSCALE", 1.);
  dparam_set ("XSHIFT", 0.);
  dparam_set ("YSCALE", 1.);
  dparam_set ("YSHIFT", 0.);
  dparam_set ("LZ",     1.);
  dparam_set ("BETA",   2.*M_PI);

  return 0;
}

int Domain_check() {
  return Geometry.name != NULL;
}

int Domain_require()
{
  int status = Domain_check();
  if (status==0)
    fprintf (stderr, "There is no computational domain\n");
  return status;
}

/* Compute the bounding box for the domain */

int Domain_bbox()
{
  Element *elmt;

  Point min;
  Point max;

  min.x = 1000.;
  min.y = 1000.;
  max.x = -min.x;
  max.y = -min.y;

  for (elmt = Geometry.mesh->head; elmt; elmt = elmt->next) {
    int    np =  elmt->nr * elmt->ns;
    double *x = *elmt->xmesh;
    double *y = *elmt->ymesh;
      
    do {
      min.x = MIN(min.x, *x);
      min.y = MIN(min.y, *y);
      max.x = MAX(max.x, *x);
      max.y = MAX(max.y, *y);
    } while (x++, y++, --np);
  }

  dparam_set("XMIN", min.x);
  dparam_set("XMAX", max.x);
  dparam_set("YMIN", min.y);
  dparam_set("YMAX", max.y);

  return 0;
}

Field* Domain_addField (char type) {
  int i = Geometry.nfields++;
  Field *u = Geometry.solVector[i] = Field_alloc(Geometry.mesh);
  FIELD_TYPE(u) = type;
  return u;
}

Field* Domain_getField (char type) {
  Field *u = Domain_chkField(type);
  if (!u) fprintf(stderr, "no such field type -- %c\n", type);
  return u;
}

Field* Domain_chkField (char type) 
{
  Field *u = NULL;
  int i;

  for (i = 0; i < Geometry.nfields; i++) {
    if (FIELD_TYPE(Geometry.solVector[i])==type) {
      u = Geometry.solVector[i];
      break;
    }
  }

  return u;
}






