/* Implementation of measure_t
 *
 * Copyright (c) 1997 Ronald D. Henderson and Caltech
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "prism/prism.h"
#include "prism/domain.h"
#include "prism/measure.h"

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "speclib/field.h"
#include "speclib/edge.h"

#define SQR(a) ((a)*(a))

/* ------------------------------------------------------------------------- */

static void separator (FILE *fp, int width) {
  int n = width-2;
  fputs ("# ", fp);
  while (n--) fputc('-',fp);
  fputc ('\n', fp);
}

/* Integrate over a single (xy)-plane */

static double integrate (const Field *u)
{
  const int npts = u->nr * u->ns;
  double sum = 0.;

  Element *elmt;
  int i;

  for (elmt = Field_head(u); elmt; elmt = elmt->next) {
    const double *mass = *elmt->mass;
    const double *valu = *elmt->field;
    for (i = 0; i < npts; i++)
      sum += mass[i] * valu[i] * valu[i];
  }

  return sqrt(sum);
}

/* ------------------------------------------------------------------------- */

measure_t *measure_alloc (struct domain *d)
{
  static char *routine = "measure";

  measure_t *m = (measure_t*) calloc(1,sizeof(measure_t));

  const int nz   = FIELD_NZ(d->U);
  const int npts = FIELD_NR(d->U)*FIELD_NS(d->U);
  const int nel  = Field_count(d->U);

  int i, k, col;

  if (!m) Prism_error("%s: couldn't allocate a measure_t", routine);

  m->body      = BC_get(d->Ubc,'W');
  m->energy    = (double*) calloc(MAX(nz,2),sizeof(double));
  m->viscosity = dparam("KINVIS");
  m->d         = d;

  for (k = 0; k < nel; k++)
    m->area += dsum(npts, *d->U[k].mass, 1);

  /* open the file and write the header */

  ROOTONLY {
    char fname[FILENAME_MAX];

    sprintf (fname, "%s.mea", d->name);
    m->fp = fopen(fname,"w");

    if (m->fp == NULL)
      Prism_error("%s: unable to open the output file -- %s\n", 
		  routine, fname);

    fprintf (m->fp, 
	     "# Measurements file\n"
	     "#\n"
	     "# Contents by column...\n"
	     "#   1  time\n");

    fprintf (m->fp, "#   %d  pressure -- FX\n", col=2);
    for (i = 1; i < DIM; i++)
      fprintf (m->fp, "#   %d              F%c\n", ++col, 'X'+i);

    fprintf (m->fp, "#   %d  viscous  -- FX\n", ++col);
    for (i= 1; i < DIM; i++)
      fprintf (m->fp, "#   %d              F%c\n", ++col, 'X'+i);

    if (dparam("FLOWRATE")>0.)
      fprintf (m->fp, "#   %d  applied pressure drop\n", ++col);
    
    fprintf (m->fp, 
	     "#   %d+ normalized energy of the kth Fourier mode\n", ++col);

    fprintf (m->fp,
	     "#\n"
	     "# Multiply by 2 if you want to convert forces reported in this\n"
	     "# file to force coefficients\n"
	     "#\n"
	     "# Normalization for energy (computational domain size) = %g\n", 
	     m->area);

    separator (m->fp, 60);
  }

  return m;
}

void measure_free (measure_t *m)
{
  if (m->fp) { fclose(m->fp); }

  free (m->body);
  free (m->energy);
  free (m);
}

/* ------------------------------------------------------------------------- */

static void pressure (measure_t *m)
{
  const Field *p = m->d->P;
  int i;

  Bedge *bc;
  double tmp[_MAX_NORDER];

  for (i = 0; i < DIM; i++)
    MEASURE_FP(m,i) = 0.;

  /* Pressure can only exert a force in (x,y) */

  for (bc = m->body; bc; bc = bc->next) {
    edge_gathr (bc->edge, p[bc->elmt->id].base[0], tmp);
    MEASURE_FP(m,0) -= edge_integrate_nx (bc->edge, tmp);
    MEASURE_FP(m,1) -= edge_integrate_ny (bc->edge, tmp);
  }
}

static void viscous (measure_t *m)
{
  Domain *d = m->d;

  const Field *u = d->U;
  const Field *v = d->V;
  const Field *w = d->W;

  Field *dudx = d->Uf[0];
  Field *dudy = d->Vf[0];

  Bedge *bc;
  double dx[_MAX_NORDER];
  double dy[_MAX_NORDER];

  int i;

  for (i = 0; i < DIM; i++)
    MEASURE_FV(m,i) = 0.;

#if DIM==2
  Field_setFrameMulti (0, 4, u, v,    dudx, dudy);
#else
  Field_setFrameMulti (0, 5, u, v, w, dudx, dudy);
#endif

  for (bc = m->body; bc; bc = bc->next) {
    const int id = bc->elmt->id;

    /* .......... X component .......... */

    Element_grad (&u[id], &dudx[id], &dudy[id]);
    edge_gathr (bc->edge, dudx[id].base[0], dx);
    edge_gathr (bc->edge, dudy[id].base[0], dy);

    MEASURE_FV(m,0) += 
      edge_integrate_ny(bc->edge, dy) + 2.*edge_integrate_nx(bc->edge, dx);
    MEASURE_FV(m,1) += 
      edge_integrate_nx(bc->edge, dy);

    /* .......... Y component .......... */

    Element_grad (&v[id], &dudx[id], &dudy[id]);
    edge_gathr (bc->edge, dudx[id].base[0], dx);
    edge_gathr (bc->edge, dudy[id].base[0], dy);

    MEASURE_FV(m,0) +=
      edge_integrate_ny(bc->edge, dx);
    MEASURE_FV(m,1) +=
      edge_integrate_nx(bc->edge, dx) + 2.*edge_integrate_ny(bc->edge, dy);

#if DIM==3

    /* .......... Z component .......... */

    Element_grad (&w[id], &dudx[id], &dudy[id]);
    edge_gathr (bc->edge, dudx[id].base[0], dx);
    edge_gathr (bc->edge, dudy[id].base[0], dy);

    MEASURE_FV(m,2) +=
      edge_integrate_nx(bc->edge, dx) + edge_integrate_ny(bc->edge, dy);
#endif

  }

  for (i = 0; i < DIM; i++)
    MEASURE_FV(m,i) *= m->viscosity;
}

static void forces (measure_t *m) {
  pressure(m);
  viscous (m);
}

static void energy (measure_t *m)
{
  Domain *d      = m->d;
  double *energy = m->energy;

  const double scal = 0.5 / m->area;

  const Field *u = d->U;
  const Field *v = d->V;
  const Field *w = d->W;

  const int nz = FIELD_NZ(u);
  int k;

  for (k = 0; k < nz; k++) {
    double eu, ev, ew;
#if DIM==2
    eu = integrate(u);
    ev = integrate(v);
    ew = 0.;
#else
    Field_setFrameMulti (k, 3, u, v, w);
    eu = integrate(u);
    ev = integrate(v);
    ew = integrate(w);
#endif
    energy[k] = scal * (eu*eu + ev*ev + ew*ew);
  }

  /* Make sure the energy in the imaginary part of mode 0 is set to zero. *
   * This corresponds to the Nyquist frequency and may or may not be      *
   * cleared by the FFT's.                                                */

  ROOTONLY { energy[1] = 0.; }
}

void measure_analyze (measure_t *m)
{
  Domain *d = m->d;

  const int    nz     = FIELD_NZ(d->U);
  const int    nmodes = MAX(nz >> 1,1);
  const double time   = dparam("TIME");
  int i, k;

  double tmp[_MAX_NZ];

  forces(m);
  energy(m);

  ROOTONLY {
    fprintf (m->fp, "%g ", time);

    for (i = 0; i < DIM; i++)
      fprintf (m->fp, "%g ", MEASURE_FP(m,i));
    for (i = 0; i < DIM; i++)
      fprintf (m->fp, "%g ", MEASURE_FV(m,i));

    if (dparam("FLOWRATE")>0.)
      fprintf (m->fp, "%g ", dparam("PDROP"));
  }

  for (i = 0, k = 0; i < nz; i += 2, k++)
    tmp[k] = MEASURE_ENERGY(m,i)+MEASURE_ENERGY(m,i+1);

#ifdef PARALLEL
#define MSGTYPE(i) (1000+i)

  if (!(i = comm_rank())) {
    for (k = 0; k < nmodes; k++)
      fprintf (m->fp, "%g ", tmp[k]);
    for (i = 1; i < comm_size(); i++) {
      comm_recv (MSGTYPE(i), tmp, nmodes*sizeof(double));
      for (k = 0; k < nmodes; k++)
	fprintf (m->fp, "%g ", tmp[k]);
    }
    fputc ('\n', m->fp);
  } else
    comm_send (MSGTYPE(i), tmp, nmodes*sizeof(double), 0);

#undef MSGTYPE
#else
  for (k = 0; k < nmodes; k++)
    fprintf (m->fp, "%g ", tmp[k]);
  fputc ('\n', m->fp);
#endif
}

#undef SQR

