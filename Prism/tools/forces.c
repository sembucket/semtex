/*
 * Compute shear stress and pressure along the walls
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include "veclib/veclib.h"
#include "speclib/speclib.h"

char *prog   = "forces";
char *author = "R. D. Henderson";
char *rcsid  = "$Id$";
char *usage  = "usage: forces [options] -r session input > output\n";
char *help   =
"options:\n"
"-h          print this help message\n"
"\n"
"This program computes the forces due to shear stress along the walls in a \n"
"flow.  The output is similar that of Prism's built-in IO_MEA calculation, \n"
"except that the itegrated values are written to the header and the wall   \n"
"values are written to the file.  After (fx,fy,fz), each of the primitive  \n"
"variables in the file is written.\n";

/* FILE *out = stdout; Changed hmb Jan 2002 */
FILE *out = NULL;

#define SQR(a) ((a)*(a))

/* ------------------------------------------------------------------------- */

static void separator (FILE *fp, int width) {
  int n = width-2;
  fputs ("# ", fp);
  while (n--) fputc('-',fp);
  fputc ('\n', fp);
}

static void load_or_clear (FieldFile *f, Field *u) 
{
  if (FieldFile_checkType(f,FIELD_TYPE(u)))
    FieldFile_get(f,u);
  else
    Field_scal (0.0,u);
}

/* ------------------------------------------------------------------------- */

main (int argc, char *argv[])
{
  FILE      *fp;
  FieldFile *f;

  Field *u, *v, *w, *p;
  Field *dudx, *dudy;
  Field *other[FIELDFILE_MAX];
  Bedge *bclist, *walls, *bc;
  int n, nflds;

  if (argc != 3) {
    fputs (usage, stderr);
    exit (-1);
  }

  out = stdout;

  /* ---------- Read the input file ---------- */

  fp = fopen(argv[1], "r");
  
  speclib_init();
  ReadParams(fp);
  u = ReadMesh(fp);
  v = Field_dup(u); FIELD_TYPE(v) = 'v';
  w = Field_dup(u); FIELD_TYPE(w) = 'w';

  dudx = Field_dup(u);
  dudy = Field_dup(u);

  bclist = ReadBCs(fp, 0, u);
  walls  = BC_get(bclist,'W');

  fclose (fp);

  /* ---------- Read the field file ---------- */

  fp = fopen(argv[2], "r");
  f  = (FieldFile*) calloc(1, sizeof(FieldFile));
  FieldFile_read (f, fp);

  load_or_clear(f,u);
  load_or_clear(f,v);
  load_or_clear(f,w);

  nflds = FIELDFILE_COUNT(f);
  for (n = 0; n < nflds; n++) {
    other[n] = Field_dup(u);
    FIELD_TYPE(other[n]) = FIELDFILE_TYPE(f,n);
    FieldFile_get(f,other[n]);
  }

  fclose (fp);

  /* Only work with the mean flow */

  if (FIELD_NZ(u)>1) {
    Field_FFT(u, -1);
    Field_FFT(v, -1);
    Field_FFT(w, -1);

    for (n = 0; n < nflds; n++)
      Field_FFT(other[n], -1);
  }

  fprintf(out, "# Pressure and viscous forces along the walls\n");
  fprintf(out, "# x   y   fx   fy   fz");
  for (n = 0; n < nflds; n++)
    fprintf (out, "   %c", FIELD_TYPE(other[n]));
  fprintf(out, "\n");

  for (bc = walls; bc; bc = bc->next) {
    const int id = bc->elmt->id;

    const double *nx = bc->edge->unx;
    const double *ny = bc->edge->uny;
    const int   npts = EDGE_NPTS(bc->edge);
    
    int i;

    double fx[_MAX_NORDER];
    double fy[_MAX_NORDER];
    double fz[_MAX_NORDER];
    double dx[_MAX_NORDER];
    double dy[_MAX_NORDER];

    /* initialize */

    double **uo = dmatrix(0, nflds-1, 0, npts);
    for (n = 0; n < nflds; n++)
      edge_gathr(bc->edge, other[n][id].base[0], uo[n]);

    /* .......... X component .......... */

    Element_grad (&u[id], &dudx[id], &dudy[id]);
    edge_gathr (bc->edge, dudx[id].base[0], dx);
    edge_gathr (bc->edge, dudy[id].base[0], dy);

    for (i = 0; i < npts; i++) {
      fx[i] = dy[i]*ny[i] + 2.*dx[i]*nx[i];
      fy[i] = dy[i]*nx[i];
      fz[i] = 0.;
    }

    /* .......... Y component .......... */

    Element_grad (&v[id], &dudx[id], &dudy[id]);
    edge_gathr (bc->edge, dudx[id].base[0], dx);
    edge_gathr (bc->edge, dudy[id].base[0], dy);

    for (i = 0; i < npts; i++) {
      fx[i] += dx[i]*ny[i];
      fy[i] += dx[i]*nx[i] + 2.*dy[i]*ny[i];
    }

    /* .......... Z component .......... */

    Element_grad (&w[id], &dudx[id], &dudy[id]);
    edge_gathr (bc->edge, dudx[id].base[0], dx);
    edge_gathr (bc->edge, dudy[id].base[0], dy);

    for (i = 0; i < npts; i++)
      fz[i] += dx[i]*nx[i] + dy[i]*ny[i];


    /* Output */

    edge_gathr (bc->edge, *u[id].xmesh, dx);
    edge_gathr (bc->edge, *u[id].ymesh, dy);

    for (i = 0; i < npts; i++) {
      fprintf (out, "%g %g %g %g %g",
	       dx[i], dy[i], fx[i], fy[i], fz[i]);
      for (n = 0; n < nflds; n++)
	fprintf (out, " %g", uo[n][i]);
      fprintf (out, "\n");
    }

    free_dmatrix(uo, 0, 0);
  }

  return 0;
}
