/* Implementation of Edge
 *
 * Copyright (c) 1997 Ronald D. Henderson and Caltech
 * $Id$
 * ------------------------------------------------------------------------- */

#include "speclib.h"
#include "edge.h"

void edge_gathr (const Edge *edge, const double *u, double *ue) {
  const int np   = edge->np;
  const int skip = edge->skip;
  int i, j;
  for (i = 0, j = edge->start; i < np; i++, j += skip)
    ue[i] = u[j];
}

void edge_scatr (const Edge *edge, const double *ue, double *u) {
  const int np   = edge->np;
  const int skip = edge->skip;
  int i, j;
  for (i = 0, j = edge->start; i < np; i++, j += skip)
    u[j] = ue[i];
}

double edge_integrate_s (const Edge *edge, const double *f)
{
  const int np = edge->np;
  const double *area = edge->area;

  double *weight, sum = 0.;
  int i;

  getops(np, 0, &weight, 0, 0);
  for (i = 0; i < np; i++)
    sum += weight[i] * area[i] * f[i];
  return sum;
}

double edge_integrate_n (const Edge *edge, const double *u, const double *v)
{
  const int np = edge->np;
  const double *area = edge->area;
  const double *nx   = edge->unx;
  const double *ny   = edge->uny;

  double *weight, sum = 0.;
  int i;

  getops(np, 0, &weight, 0, 0);
  for (i = 0; i < np; i++)
    sum += weight[i] * area[i] * (nx[i]*u[i] + ny[i]*v[i]);
  return sum;
}

double edge_integrate_nx (const Edge *edge, const double *u)
{
  const int np = edge->np;
  const double *area = edge->area;
  const double *nx   = edge->unx;

  double *weight, sum = 0.;
  int i;

  getops(np, 0, &weight, 0, 0);
  for (i = 0; i < np; i++)
    sum += weight[i] * area[i] * nx[i] * u[i];
  return sum;
}

double edge_integrate_ny (const Edge *edge, const double *u)
{
  const int np = edge->np;
  const double *area = edge->area;
  const double *ny   = edge->uny;

  double *weight, sum = 0.;
  int i;

  getops(np, 0, &weight, 0, 0);
  for (i = 0; i < np; i++)
    sum += weight[i] * area[i] * ny[i] * u[i];
  return sum;
}

