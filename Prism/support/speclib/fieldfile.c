/*
 * FieldFile implementation
 *
 * $Revision$
 *
 * Author:  R. D. Henderson
 *
 * -------------------------------------------------------------------- */

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "comm/comm.h"
#include "veclib/veclib.h"
#include "speclib/speclib.h"

#include "speclib/fieldfile.h"
#include "speclib/config.h"

/* Internal prototypes */

int FieldFile_write_sp    (FieldFile *f, FILE *fp);
int FieldFile_write_mp    (FieldFile *f, FILE *fp);
int FieldFile_read_sp     (FieldFile *f, FILE *fp);
int FieldFile_read_mp     (FieldFile *f, FILE *fp);
int FieldFile_projectz_sp (FieldFile *f, int nz);
int FieldFile_projectz_mp (FieldFile *f, int nz);

/* Deallocate memory for any stored Field data */

static FieldFile *empty (FieldFile *f) 
{
  if (f->count != 0) {
    int n = 0;
    do
      free (f->data[n]);
    while (++n < f->count);
    f->count = 0;
  }

  memset (f->type, 0, FIELDFILE_STRSIZE);
  memset (f->data, 0, FIELDFILE_MAX*sizeof(double*));

  return f;
}

/* ------------------------------------------------------------------------- */

FieldFile *FieldFile_alloc () 
{
  FieldFile *f = (FieldFile*) calloc(1,sizeof(FieldFile));

  f->count     = 0;
  f->step      = iparam("STEP");
  f->time      = dparam("TIME");
  f->time_step = dparam("DT");
  f->Re        = dparam("Re");
  f->beta      = dparam("BETA");

  sprintf  (f->name,    "unknown");
  sprintf  (f->created, "unknown");
  sprintf  (f->format,  "binary-%s", ARCH);

  return empty(f);
}

void FieldFile_free (FieldFile *f) 
{
  empty(f);
  free (f);
}

int FieldFile_setSize (FieldFile *f, 
		       int count, int nr, int ns, int nz, int nel)
{
  const int npts = nr*ns*nz*nel;
  int n;

  f->count = count;
  f->nr    = nr;
  f->ns    = ns;
  f->nz    = nz;
  f->nel   = nel;

  for (n = 0; n < count; n++) {
    f->type[n] = '?';
    f->data[n] = (double*) malloc(npts*sizeof(double));
  }

  return 0;
}

int FieldFile_put (FieldFile *f, const Field *u)
{
  const int n = f->count;
  int npts;

  if (f->count == 0) {
    f->nr        = u->nr;
    f->ns        = u->ns;
    f->nz        = u->nz;
    f->nel       = Field_count(u);
  } else if (f->nr != u->nr || f->ns != u->ns || f->nz != u->nz)
    return speclib_error("fieldfile: input field is the wrong size!\n");

  npts       = f->nr * f->ns * f->nz * f->nel;
  f->type[n] = u->type;
  f->data[n] = dvector (0, npts-1);

  dcopy (npts, *u->base, 1, f->data[n], 1);
  
  return f->count++;
}

int FieldFile_get (const FieldFile *f, Field *u)
{
  const int npts  = f->nr * f->ns * f->nz * f->nel;
  const int count = f->count;
  int n;

  for (n = 0; n < count; n++) {
    if (f->type[n] == u->type) {
      dcopy (npts, f->data[n], 1, *u->base, 1);
      return n;
    }
  }
  
  speclib_error("fieldfile: field type %c not found\n", u->type);
  return 0;
}

/* Headers */

static int header_write (FieldFile *f, FILE *fp)
{
  char resolution[FIELDFILE_STRSIZE];

  /* Create the time-stamp */

  time_t tp = time((time_t*)NULL);
  strftime (f->created, FIELDFILE_STRSIZE-1, "%a %b %d %H:%M:%S %Y",
	    localtime(&tp));

  /* Create the resolution string */

  sprintf (resolution, "%d %d %d %d", f->nr, f->ns, f->nz, f->nel);

  /* Write the header info */

  fprintf (fp, "%-25s Session\n", f->name);
  fprintf (fp, "%-25s Created\n", f->created);
  fprintf (fp, "%-25s Nr Ns Nz Nelements\n", resolution);
  fprintf (fp, "%-25d Step\n", f->step);
  fprintf (fp, "%-25g Time\n", f->time);
  fprintf (fp, "%-25g Time Step\n", f->time_step);
  fprintf (fp, "%-25g Re\n", f->Re);
  fprintf (fp, "%-25g Beta\n", f->beta);
  fprintf (fp, "%-25s Fields Written\n", f->type);
  fprintf (fp, "%-25s Format\n", f->format);

  return 0;
}

/* NB: First test to see if there is another FieldFile on the input        *
 * stream.  If there is then load it, but if there isn't just return the   *
 * EOF flag and leave any data currently in the FieldFile alone.           *
 *                                                                         *
 * DO NOT CHANGE THIS BEHAVIOR.  Some applications rely on _read() giving  *
 * back an EOF with the data left intact in order to operate on the last   *
 * dump in a file.                                                         */

static int gettypes (char *s, char *t)
{
  char *p = strstr(s, "Field");
  int   n = 0;
  if (!p) speclib_error("invalid \"Field Types\" string");
  memset(t, 0, FIELDFILE_MAX);
  while (s != p && n < FIELDFILE_MAX) {
    if (isalpha(*s)) t[n++] = *s;
    s++;
  }
  return n;
}

#define READLINE(fp,fmt,args) \
  if (fscanf(fp,fmt,args)==EOF) perror("fieldfile"); fgets(buf,BUFSIZ,fp);

static int header_read (FieldFile *f, FILE *fp) 
{
  char buf[BUFSIZ];

  /* Read the session name */
  
  if (fscanf(fp, "%s", buf) == EOF)
    return FIELDFILE_EOF;          
  else
    empty (f);

  strncpy (f->name, buf, FIELDFILE_STRSIZE-1);
  fgets   (buf, BUFSIZ, fp);

  /* Read the "created" flag */

  fgets (buf, BUFSIZ, fp);
  strncpy (f->created, buf, 25);

  /* Read the the resolution parameters */

  fgets (buf, BUFSIZ, fp);
  sscanf(buf, "%d%d%d%d", &f->nr, &f->ns, &f->nz, &f->nel);

  /* Read the simulation parameters */

  READLINE (fp, "%d" , &f->step);  
  READLINE (fp, "%lf", &f->time);
  READLINE (fp, "%lf", &f->time_step);
  READLINE (fp, "%lf", &f->Re);
  READLINE (fp, "%lf", &f->beta);
  
  fgets (buf, BUFSIZ, fp);
  f->count = gettypes(buf, f->type);

  fscanf (fp, "%s", f->format);
  fgets  (buf, BUFSIZ, fp);
  
  return 0;
}

#undef READLINE

int FieldFile_header (FieldFile *f, FILE *fp) {
  return header_write(f,fp);
}

int FieldFile_write (FieldFile *f, FILE *fp)
{
#ifdef PARALLEL
  return FieldFile_write_mp (f, fp);
#else
  return FieldFile_write_sp (f, fp);
#endif
}

int FieldFile_read (FieldFile *f, FILE *fp)
{
#ifdef PARALLEL
  return FieldFile_read_mp (f,fp);
#else
  return FieldFile_read_sp (f,fp);
#endif
}

/* ------------------------------------------------------------------------- *
 * Single-processor versions of _read() and _write().                        *
 * ------------------------------------------------------------------------- */

int FieldFile_write_sp (FieldFile *f, FILE *fp)
{
  const int count = f->count;
  const int npts  = f->nr * f->ns * f->nz * f->nel;
  int i, n;

  /* Write the header */

  FieldFile_header (f, fp);     

  switch (tolower(*f->format)) {
  case 'a':
    for (i = 0; i < npts; i++) {
      for (n = 0; n < count; n++)
	if (fprintf (fp, "%#16.10g ", f->data[n][i]) < 0)
	  perror ("fieldfile: write ascii:");
      fputc ('\n', fp);
    }
    break;
    
  case 'p': {
    for (n = 0; n < count; n++)
      for (i = 0; i < npts; i++) {
	float tmp = (float) f->data[n][i];
	if (fwrite (&tmp, sizeof(float), 1, fp) != 1)
	  perror ("fieldfile: write packed:");
      }
    break;
  }
    
  case 'b':
    for (n = 0; n < count; n++)
      if (fwrite (f->data[n], sizeof(double), npts, fp) != npts)
	perror ("fieldfile: write binary");
    break;

  default:
    speclib_error("fieldfile: unknown format -- %s\n", f->format);
    break;
  }

  fflush (fp);
  return count;
}

int FieldFile_read_sp (FieldFile *f, FILE *fp)
{
  /* Read the header */

  const int status = header_read(f,fp);

  /* Read the next dump... */

  if (status != FIELDFILE_EOF) {
    const int npts  = f->nr * f->ns * f->nz * f->nel;
    const int count = f->count;
    int i, n;

    for (n = 0; n < count; n++)
      f->data[n] = dvector(0, npts-1);

    /* Load the data */

    switch (tolower(*f->format)) {
    case 'a':
      for (i = 0; i < npts; i++) {
	for (n = 0; n < count; n++)
	  if (fscanf (fp, "%lf", & f->data[n][i]) != 1)
	    perror ("fieldfile: read ascii");
      }
      while (fgetc(fp) != '\n') 
	{ /* read to end-of-line */ }
      break;
      
    case 'p':
      for (n = 0; n < count; n++) {
	for (i = 0; i < npts; i++) {
	  float tmp;
	  if (fread (&tmp, sizeof(float), 1, fp) != 1)
	    perror ("fieldfile: read packed");
	  f->data[n][i] = (double) tmp;
	}
      }
      break;
      
    case 'b':
      for (n = 0; n < count; n++)
	if (fread (f->data[n], sizeof(double), npts, fp) != npts)
	  perror ("fieldfile: read binary");
      break;
      
    default:
      speclib_error("fieldfile: unknown format -- %s\n", f->format);
      break;
    }
  }

  return status;
}

/* ---------------------------------------------------------------------- *
 * FieldFile_interp   -- Interpolate a FieldFile to a different basis     *
 * FieldFile_project  -- Interpolate a FieldFile to a different order     *
 *                                                                        *
 * This function interpolates a FieldFile to a different elemental basis. *
 * For a multiframe field, it interpolates one frame at a time.  It will  *
 * not interpolate to a different number of frames.                       *
 *                                                                        *
 * FieldFile_project is a wrapper that projects the data to another GLL   *
 * mesh of a higher or lower polynomial order.                            *
 *                                                                        *
 * Input:   f       FieldFile w/data assumed to be at the GLL points      *
 *          nx      Number of points in the (xx)-direction                *
 *          xx      Array of points over the range  -1 <= xx <= 1.        *
 *          ny      Number of points in the (yy)-direction                *
 *          yy      Array of points over the range  -1 <= yy <= 1.        *
 *                                                                        *
 * Return value: -1 for error, 1 otherwise.                               *
 * ---------------------------------------------------------------------- */

/* Interpolation function */

static double 
 *data_interp (int nr, int ns, int nz, int nel, double *data,
	       int nx, int ny, double **imr, double **ims)
{
  const int nrns = nr * ns;
  const int nxny = nx * ny;
  double *nvec = dvector (0, nxny * nz * nel - 1), *p = nvec;
  int k;

  double tmp [_MAX_NORDER*_MAX_NORDER];

  for (k = 0; k < nel*nz; k++, data += nrns, p += nxny) {
    dgemm ('N', 'N', nr, ny, ns, 1., data, nr, *ims, ns, 0., tmp, nr);
    dgemm ('T', 'N', nx, ny, nr, 1., *imr, nr,  tmp, nr, 0.,   p, nx);
  }

  return nvec;
}

void FieldFile_project (FieldFile *f, int nr, int ns)
{
  if (f->nr != nr || f->ns != ns) {
    double *zr, *zs;

    coef (nr); getops (nr, &zr, 0, 0, 0);
    coef (ns); getops (ns, &zs, 0, 0, 0);

    FieldFile_interp (f, nr, zr, ns, zs);
  }
}

void FieldFile_interp (FieldFile *f, int nx, double *xx, int ny, double *yy)
{
  int nr        = f->nr;     
  int ns        = f->ns;     
  int nz        = f->nz;
  int nel       = f->nel;
  int nfields   = strlen(f->type);
  double **data = f->data;

  double   **imr, **itmr, **ims, **itms, *zr, *zs;

  if (nx > _MAX_NORDER || ny > _MAX_NORDER)
    speclib_error("FieldFile_interp: maximum order is %d", _MAX_NORDER);

  imr     = dmatrix (0, nx-1, 0, nr-1);
  itmr    = dmatrix (0, nr-1, 0, nx-1);
  ims     = dmatrix (0, ny-1, 0, ns-1);
  itms    = dmatrix (0, ns-1, 0, ny-1);

  /* Compute the GLL points */

  coef(nr);  getops(nr, &zr, 0, 0, 0);
  coef(ns);  getops(ns, &zs, 0, 0, 0);

  /* Compute the interpolation matrices */

  igllm (imr, itmr, zr, xx, nr, nx);
  igllm (ims, itms, zs, yy, ns, ny);

  /* Loop through the fields and interpolate each one */

  while (nfields--) {
    double *p = data_interp (nr, ns, nz, nel, *data, nx, ny, imr, ims);
    free (*data);
    *data++ = p;
  }

  f->nr = nx;
  f->ns = ny;

  free_dmatrix (imr , 0, 0);
  free_dmatrix (itmr, 0, 0);
  free_dmatrix (ims , 0, 0);
  free_dmatrix (itms, 0, 0);
}

/* ----------------------------------------------------------------------- *
 * The following function provide access to the FieldFile data structure.  *
 * ----------------------------------------------------------------------- */

int FieldFile_getSize (const FieldFile *f, 
		       int* nr, int* ns, int* nz, int* nel)
{
  *nr  = f->nr;
  *ns  = f->ns;
  *nz  = f->nz;
  *nel = f->nel;

  return f->nr * f->ns * f->nz * f->nel;
}

int FieldFile_getTypeList (const FieldFile *f, char *list) {
  strncpy (list, f->type, FIELDFILE_MAX);
  return strlen(list);
}

int FieldFile_setFormat (FieldFile *f, format_t format) 
{
  switch (format) {
  case BINARY:
    sprintf (f->format, "binary-%s", ARCH);
    break;
  case PACKED:
    sprintf (f->format, "packed-%s", ARCH);
    break;
  case ASCII:
    sprintf (f->format, "ascii");
    break;
  default:
    break;
  }

  return (int) format;
}

char* FieldFile_setName (FieldFile *f, const char *name) {
  return strncpy(f->name, name, FIELDFILE_STRSIZE-1);
}

char* FieldFile_getName (const FieldFile *f, char *name) {
  return strncpy(name, f->name, FIELDFILE_STRSIZE-1);
}

int FieldFile_checkType (const FieldFile *f, char type)
{
  int i = 0;
  int n = f->count;
  int status = 0;

  for (i = 0; i < n; i++) {
    if (f->type[i] == type) {
      status = 1;
      break;
    }
  }

  return status;
}

/* ------------------------------------------------------------------------- *
 * FieldFile_projectz() -- Fourier interpolation                             *
 *                                                                           *
 * This function interpolates a field file to a different number of frames   * 
 * using Fourier interpolation.                                              * 
 * ------------------------------------------------------------------------- */

void FieldFile_projectz (FieldFile *f, int nz)
{
#ifdef PARALLEL
  FieldFile_projectz_mp (f,nz);
#else
  FieldFile_projectz_sp (f,nz);
#endif
}

static void transpose (int dir, int nxy, int nz, double *da, double *db)
{
  int i;

  if (dir > 0) {                /* transpose da[nxy][nz] -> db[nz][nxy] */
    for (i = 0; i < nxy; i++)
      dcopy (nz, da+i, nxy, db + i*nz, 1);

  } else {                      /* transpose db[nz][nzy] -> da[nxy][nz] */
    for (i = 0; i < nz; i++)
      dcopy (nxy, db + i, nz, da + i*nxy, 1);
  }
}

int FieldFile_projectz_sp (FieldFile *f, int nz)
{
  if (FIELDFILE_NZ(f) != nz) {
    const int nxy = FIELDFILE_NR(f)*FIELDFILE_NS(f)*FIELDFILE_NEL(f);
    const int nz1 = FIELDFILE_NZ(f);
    const int nz2 = nz;

    double *wrk1 = dvector(0,nxy*nz1);
    double *wrk2 = dvector(0,nxy*nz2);

    const int nflds = FIELDFILE_COUNT(f);
    int i, j, n;

    for (n = 0; n < nflds; n++) {

      /* Assemble the transpose */
      transpose (+1, nxy, nz1, FIELDFILE_DATA(f,n), wrk1);

      /* Fourier transform (unless nz1=1, i.e. 2D -> 3D projection) */

      if (nz1 > 1) {
	for (i = 0; i < nxy; i++)
	  realft(nz1>>1, wrk1+i*nz1, -1);
      }

      /* Truncate or 0-extend the data */
      if (nz2 > nz1) {
	const int nzd = nz2-nz1;
	for (j = 0; j < nxy; j++) {
	  dcopy (nz1, wrk1+j*nz1, 1, wrk2+j*nz2, 1);
	  dzero (nzd, wrk2+j*nz2+nz1, 1);
	}
      } else {
	for (j = 0; j < nxy; j++) {
	  dcopy (nz2, wrk1+j*nz1, 1, wrk2+j*nz2, 1);
	}
      }

      /* Inverse transform (unless nz2=1, i.e. 3D -> 2D projection) */
      
      if (nz2 > 1) {
	for (i = 0; i < nxy; i++)
	  realft(nz2>>1, wrk2+i*nz2, +1);
      }

      /* Reallocate and assemble the new data */
      free(FIELDFILE_DATA(f,n));
      FIELDFILE_DATA(f,n) = dvector(0,nxy*nz2);
      transpose (-1, nxy, nz2, FIELDFILE_DATA(f,n), wrk2);
    }

    FIELDFILE_NZ(f) = nz;

    free(wrk1);
    free(wrk2);
  }
  return 0;
}

#ifdef PARALLEL

/* ------------------------------------------------------------------------- *
 * _mp(): parallel I/O routines                                              *
 * ------------------------------------------------------------------------- */

#define FLD_HANDLE 1024    /* communications handle */
#define FLD_NBLOCK 1025    /* block size message */
#define FLD_DATA   1026    /* block data message */
#define FLD_CTS    1027    /* clear-to-send message */

int FieldFile_write_mp (FieldFile *f, FILE *fp)
{
  const int count = f->count;
  const int npts  = f->nr * f->ns * f->nz * f->nel;
  const int rank  = comm_rank();
  const int size  = comm_size();

  char buf[BUFSIZ]; /* communications buffer */

  int n, p;

  if (rank==0) {
    f->nz *= size;
    FieldFile_header (f, fp);
    f->nz /= size;
  }

  switch (tolower(*f->format)) {
  case 'b':
    for (n = 0; n < count; n++) {

      /* Node zero collects data and transfers it to the file */

      if (rank==0) {

	fwrite (f->data[n], sizeof(double), npts, fp);
	for (p = 1; p < size; p++) {
	  int nblock;

	  comm_send (FLD_HANDLE, buf, 0, p);

	  do {
	    comm_recv (FLD_NBLOCK, &nblock, sizeof(int));
	    if (nblock > 0) {
	      comm_recv (FLD_DATA, buf, nblock*sizeof(double));
	      fwrite (buf, sizeof(double), nblock, fp);
	      comm_send (FLD_CTS, buf, 0, p);
	    }
	  } while (nblock > 0);
	}
      }

      /* All other nodes just send data to node zero */

      else {
	double *ptr   = f->data[n];
	int    nwords = npts;
	int    nblock = MIN(BUFSIZ/sizeof(double), npts);

	comm_recv (FLD_HANDLE, buf, 0);

	while (nwords > 0) {
	  comm_send (FLD_NBLOCK, &nblock, sizeof(int), 0);
	  comm_send (FLD_DATA, ptr, nblock*sizeof(double), 0);

	  nwords -= nblock;
	  ptr    += nblock;
	  nblock  = MIN(nblock, nwords);
	  comm_recv (FLD_CTS, buf, 0);
	}

	nblock = 0;
	comm_send (FLD_NBLOCK, &nblock, sizeof(int), 0);
      }
    }
    break;

  default:
    speclib_error ("FieldFile_write_mp: only binary output in parallel");
    break;
  }

  if (rank==0) fflush(fp);

  return 0;
}

int FieldFile_read_mp (FieldFile *f, FILE *fp)
{
  const int rank = comm_rank();
  const int size = comm_size();

  char buf[BUFSIZ]; /* communications buffer */

  int status;

  if (rank==0) status = header_read(f,fp);
  comm_bcast (&status, sizeof(int), 0);

  if (status != FIELDFILE_EOF) {
    comm_bcast (f, sizeof(FieldFile), 0);
    f->nz /= size;
  }

  /* Read the next dump... */

  if (status != FIELDFILE_EOF) {
    const int npts  = f->nr * f->ns * f->nz * f->nel;
    const int count = f->count;
    int n, p;

    for (n = 0; n < count; n++)
      f->data[n] = dvector(0, npts-1);

    /* Load the data */

    switch (tolower(*f->format)) {
    case 'b':
      for (n = 0; n < count; n++) {

	/* Node zero collects data and transfers it to the file */

	if (rank==0) {

	  if (fread (f->data[n], sizeof(double), npts, fp) != npts)
	    perror ("FieldFile_read_mp: read binary");

	  for (p = 1; p < size; p++) {
	    int nwords = npts;
	    int nblock = MIN(BUFSIZ/sizeof(double), npts);

	    while (nwords > 0) {

	      if (fread (buf, sizeof(double), nblock, fp) != nblock)
		perror ("FieldFile_read_mp: read binary");

	      comm_send (FLD_NBLOCK, &nblock, sizeof(int), p);
	      comm_send (FLD_DATA, buf, nblock*sizeof(double), p);

	      nwords -= nblock;
	      nblock  = MIN(nblock,nwords);
	    }

	    nblock = 0;
	    comm_send (FLD_NBLOCK, &nblock, sizeof(int), p);
	  }
	}

      /* All other nodes just recv data from node zero */

	else {
	  double *ptr = f->data[n];
	  int nwords = 0;
	  int nblock;

	  do {
	    comm_recv (FLD_NBLOCK, &nblock, sizeof(int));
	    if (nblock > 0)
	      comm_recv (FLD_DATA, ptr, nblock*sizeof(double));
	    ptr    += nblock;
	    nwords += nblock;
	  } while (nblock > 0);

	  if (nwords != npts) {
	    speclib_error
	      ("FieldFile_read_mp: node %d expected %d words, got %d\n",
	       rank, npts, nwords);
	  }
	}
      }
      break;

    case 'a':
    case 'p':
      speclib_error("FieldFile_read_mp: only binary input in parallel");
      break;

    default:
      speclib_error("fieldfile: unknown format -- %s\n", f->format);
      break;
    }
  }

  return status;
}

int FieldFile_projectz_mp (FieldFile *f, int nz)
{
  speclib_error ("FieldFile_projectz_mp: not implemented yet\n");
  return 0;
}

#endif
