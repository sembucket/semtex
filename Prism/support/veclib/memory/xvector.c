/*
 * vector allocation with arbitrary indexing
 */

#include <stdio.h>
#include <stdlib.h>

double *dvector(int nl, int nh)
{
  double *v;
  
  v = (double*) malloc((nh-nl+1)*sizeof(double));
  if (!v) fprintf (stderr, "veclib: out of memory\n");
  return v-nl;
}

float *svector(int nl, int nh)
{
  float *v;

  v = (float*) malloc((nh-nl+1)*sizeof(float));
  if (!v) fprintf (stderr, "veclib: out of memory\n");

  return v-nl;
}

int *ivector(int nl, int nh)
{
  int *v;
  
  v = (int*) malloc((nh-nl+1)*sizeof(int));
  if (!v) fprintf (stderr, "veclib: out of memory\n");
  return v-nl;
}

void free_dvector(double *v, int nl)
{
  free(v+nl);
}

void free_svector(float *v, int nl)
{
  free(v+nl);
  return;
}

void free_ivector(int *v, int nl)
{
  free(v+nl);
  return;
}
