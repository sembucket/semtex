/*
 * Aliases
 *
 * $Revision$
 *
 * Author:       R. D. Henderson
 * ------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "veclib/veclib.h"

#include "cubit.h"
#include "vertex.h"
#include "edge.h"
#include "alias.h"
#include "element.h"

/* eta() -- The mortar test function of order N-2, eqn. (2.29) */

#define SIGN(n)   (1.-2.*((n) % 2))

static double eta (int q, double z, double *zgll, int n)
{
  const double alpha = n * (n + 1.);
  const double dz    = zgll[q] - z;

  if (fabs(dz*dz) < FLT_EPSILON)
    return SIGN(n-q-1) * pnd2leg(z,n);
  if (fabs(z*z - 1.) < FLT_EPSILON)
    return SIGN(n-q) * pow(z, n-1.) * alpha * 0.5 / dz;
  else
    return SIGN(n-q) * pndleg(z,n) / dz;
}

#undef SIGN

/* ------------------------------------------------------------------------ *
 * form_q() -- Form the projection operator                                 *
 *                                                                          *
 * Input:   n   ...   order of the polynomial along the edges (nr|ns - 1)   *
 *          m   ...   order of the polynomial along the mortar              *
 *          s0  ...   mortar offset                                         *
 *                                                                          *
 * The "gamma_x" arguments are the length of the mortar (p), the integra-   *
 * tion strip (s), and the element edge (k).  The mortar offset is the      *
 * signed distance from the edge of the ELEMENT to the edge of the MORTAR   *
 * (see below).                                                             *
 *                                                                          *
 *                        | o    o     o     o    o |  element edge         *
 *                                                                          *
 *   |<------  s0  -------|=========================|  integration strip    *
 *       (i.e, s0 < 0)                                                      *
 *   |..........................................................| mortar    *
 *                                                                          *
 * ------------------------------------------------------------------------ */

static double **form_q 
(int n, int m, double s0, double gamma_p, double gamma_s, double gamma_k)
{
  const double tol = FLT_EPSILON;
  const int np  = n + 1;
  const int mp  = m + 1;

  double **q = dmatrix(0, n, 0, m);
  double *z1 = dvector(0, n);
  double *z2 = dvector(0, n);

  double *zs, *zm, *w, **d, fac;
  int i, j, k;
  
  /* Get collocation points & weights for the integration */

  getops (mp, &zm, &w, &d, &d);
  getops (np, &zs, &w, &d, &d);
  dzero  (np*mp, *q, 1);

  for (i = 0; i < np; i++) {
    z1[i] = (1. + zs[i]) * gamma_s / gamma_k - 1.;
    z2[i] = (1. + zs[i]) * gamma_s / gamma_p - 1.;
  }

  if (s0*s0 > tol*gamma_p) {
    if (s0 >= 0.) {
      fac = 2. * fabs(s0) / gamma_k;
      dsadd(np, fac, z1, 1, z1, 1);
    } else {
      fac = 2. * fabs(s0) / gamma_p;
      dsadd(np, fac, z2, 1, z2, 1);
    }
  }

  /* Compute the P-operator from eqns. (2.32-2.33) */

  fac = gamma_s / 2.;
  for (i = 1; i < n; i++)
    for (j = 0; j <= m; j++)
      for (k = 0; k <= m; k++)
	q[i][j] += fac * w[k] * eta(i, z1[k], zs, n) * hgll(j, z2[k], zm, mp);


  /* Enforce the vertex-pinning conditions, eqn. (2.18-2.19) */

  fac = gamma_k / 2.;
  if (s0 < tol) {
    for (i = 1; i < n; i++)
      for (j = 0; j <= m; j++)
	q[i][j] -= fac * w[0] * eta(i, z1[0], zs, n) * hgll(j, z2[0], zm, mp);
    
    for (i = 0, j = 0; j <= m; j++)        /* Vertex condition */
      q[i][j] = hgll(j, z2[i], zm, mp);
  }


  /* Convert s0 to an offset from the other end and check that one */

  s0 = gamma_k - (s0 + gamma_p);
  if (s0 < tol) {
    for (i = 1; i < n; i++)
      for (j = 0; j <= m; j++)
	q[i][j] -= fac * w[n] * eta(i, z1[n], zs, n) * hgll(j, z2[n], zm, mp);

    for (i = n, j = 0; j <= m; j++)        /* Vertex condition */
      q[i][j] = hgll (j, z2[i], zm, mp);
  }

  /* Now divide by the B-operator, eqn. (2.31) */

  for (i = 1; i < n; i++) {
    fac = 1. / (.5 * gamma_k * w[i] * eta(i, zs[i], zs, n));
    dscal (mp, fac, q[i], 1);
  }

  free (z1);
  free (z2);

  return q;
}

/* ------------------------------------------------------------------------- */

Alias *Alias_alloc (void) {
  Alias *alias = (Alias*) calloc (1, sizeof(Alias));
  assert(alias);
  return alias;
}

void Alias_free (Alias *alias) {
  free (alias);
}

void Alias_build (Alias *alias, double *pos, int type, int np)
{
  static int    n;     /* only maintain one set of matrices */
  static double **Z1;
  static double **Z2;

  if (np == 0)
    return;

  if (n == 0) {
    n  = np - 1;
    Z1 = form_q (n, n,  0., 2., 1., 1.);
    Z2 = form_q (n, n, -1., 2., 1., 1.);
  }

  memcpy (alias->pos, pos, 2 * sizeof(double));

  switch (alias->type = type) {
  case 0:
    break;
  case 1:
    alias->Z = Z1;
    break;
  case 2:
    alias->Z = Z2;
    break;
  default:
    break;
  }
}



