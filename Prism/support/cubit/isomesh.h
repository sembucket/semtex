#ifndef ISOMESH_H
#define ISOMESH_H

/*
 * Copyright (c) 1995 Ronald Dean Henderson
 *
 * Isoparametric Mesh Generation, prototypes and type declarations
 * 
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include "element.h"
#include "edge.h"

/* Normal sequence is to call the functions in this order: */

void shape   (Element *U);
void curve   (Element *U, Edge *edge);   
void blend   (Element *U);
void map     (Element *U);
void normals (Element *U);

/* Shape-generating functions */

struct line      make_line   (Element *, Edge *);
struct arc       make_arc    (Element *, Edge *, double radius);
struct spline_t  spline_alloc (Element *, Edge *, double pos[]);

struct geofile_t geofile_alloc (const char *name);

void spline_node (const spline_t *s, int node, double *x, double *y);
void spline_cut  (const spline_t *s, spline_t *s1, spline_t *s2);

#endif
