/*
 * FieldFiles
 *
 * RCS Information
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * -------------------------------------------------------------------- */

#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "cubit.h"

#include "veclib/veclib.h"

#define   DESCRIP   25     /* Point at which header descriptions begin */

static char *hdr_fmt[] = { 
  "%-25s "    "Session\n",
  "%-25s "    "Created\n",
  "%-25s "    "Nr, Ns, Nz, Elements\n",
  "%-25d "    "Step\n",
  "%-25.6g "  "Time\n",
  "%-25.6g "  "Time step\n",
  "%-25.6g "  "Reynolds number\n",
  "%-25.6g "  "Beta\n",
  "%-25s "    "Fields Written\n",
  "%-25s "    "Format\n"
};
                                   /* .......  String constants ...... */ 
#define BINARY   strings[0]	   /* 0. Binary format file            */ 
#define ASCII    strings[1]	   /* 1. ASCII format file             */ 
#define PACKED   strings[2]        /* 2. Packed binary file            */
#define REYNOLDS strings[3]	   /* 3. Reynolds number               */ 
#define BETA     strings[4]	   /* 4. Spanwise wavenumber           */ 
#define DT       strings[5]	   /* 5. Time step                     */ 
#define BINFMT   strings[6]        /* 6. Binary format string          */
#define CUBIT    strings[7]        /* 7. Cubit indentifier           */
#define CHECKPT  strings[8]        /* 8. Checkpoint option             */
#define NOGLL    strings[9]        /* 9. No GLL interpolation option   */

static char *strings [] = { 
  "binary", "ascii", "packed", "Re", "BETA", "DT",
#ifdef ARCH
  "binary-" ARCH,
#else
  "binary",
#endif
  "cubit", "checkpt", "nogll"
};


/* Private functions */

static int gettypes (char *t, char *s);
static int checkfmt (char *format);

/* ----------------------------------------------------------------------- *
 * FieldFile_alloc() - Allocate a FieldFile                                *
 * ----------------------------------------------------------------------- */

FieldFile *FieldFile_alloc (void)
{
  time_t    tp;
  char      buf[DESCRIP+1];
  FieldFile *f;

  /* Get the current date and time from the system */

  tp = time((time_t*) NULL);
  strftime (buf, DESCRIP, "%a %b %d %H:%M:%S %Y", localtime(&tp));

  /* Allocate space and initialize the structure */

  f = (FieldFile*) calloc (1, sizeof(FieldFile));
  assert (f);

  /* Initialize the FieldFile structure */

  f->created = strdup(buf);
  f->format  = BINFMT;

  return f;
}

/* ---------------------------------------------------------------------- *
 * FieldFile_free() -- Free the memory associated with a FieldFile        *
 * ---------------------------------------------------------------------- */

void FieldFile_free (FieldFile *f)
{
  int n = strlen(f->type);
  double **d = f->data;

  if (f->name)    free (f->name);
  if (f->created) free (f->created);
  if (f->format)  free (f->format);

  while (*d && n--) free (*d++);

  free (f);
  return;
}

/* ---------------------------------------------------------------------- *
 * FieldFile_write() -- Write a field file                                *
 *                                                                        *
 * This function writes a field file from a FieldFile structure.  The     *
 * format is a simple array, stored in one of the following formats:      *
 *                                                                        *
 *        ascii              ASCII text                                   *
 *        binary             64 bit real                                  *
 *        packed             32 bit real                                  *
 *                                                                        *
 * NB: Internal format is always binary (64 bit).                         *
 *                                                                        *
 * Return value: number of variables written, -1 for error.               *
 * ---------------------------------------------------------------------- */
  
int FieldFile_write (FieldFile *f, FILE *fp)
{
  int nfields, ntotz;
  int i, n;

  FieldFile_header (f, fp);     /* Write the header */


  /* Write the field files (no error checking yet) */

  ntotz   = f->nr * f->ns * f->nz * f->nel;
  nfields = strlen (f->type);

  switch (tolower(*f->format)) {
  case 'a':
    for (i = 0; i < ntotz; i++) {
      for (n = 0; n < nfields; n++)
	if (fprintf (fp, "%#16.10g ", f-> data[n][i]) != 1) {
	  perror (CUBIT);
	  return -1;
	}
      fputc ('\n', fp);
    }
    break;
    
  case 'b':
    for (n = 0; n < nfields; n++)
      if (fwrite (f->data[n], sizeof(double), ntotz, fp) != ntotz)
	perror (CUBIT);
    break;
    
  case 'p': {
    float tmp;
    for (n = 0; n < nfields; n++)
      for (i = 0; i < ntotz; i++) {
	tmp  = (float) f-> data[n][i];
	if (fwrite (&tmp, sizeof(float), 1, fp) != 1)
	  perror (CUBIT);
      }
  }

  default:
    sprintf(cubit_err_msg, "unknown format -- %s", f->format);
    cubit_err(NULL);
    break;
  }

  fflush (fp);
  return nfields;
}

/* Write a FieldFile header */

void FieldFile_header (FieldFile *f, FILE *fp)
{
  char sizbuf [DESCRIP+1];

  sprintf (sizbuf, "%d %d %d %d", f->nr, f->ns, f->nz, f->nel);

  fprintf (fp, hdr_fmt[0], f->name);
  fprintf (fp, hdr_fmt[1], f->created);
  fprintf (fp, hdr_fmt[2], sizbuf);
  fprintf (fp, hdr_fmt[3], f->step);
  fprintf (fp, hdr_fmt[4], f->time);
  fprintf (fp, hdr_fmt[5], f->time_step);
  fprintf (fp, hdr_fmt[6], f->Re);
  fprintf (fp, hdr_fmt[7], f->beta);
  fprintf (fp, hdr_fmt[8], f->type);
  fprintf (fp, hdr_fmt[9], f->format);

  return;
}

/* ---------------------------------------------------------------------- *
 * FieldFile_read() -- Read a field file                                  *
 *                                                                        *
 * This function loads a field file into a FieldFile structure.  Formats  *
 * are the same as in FieldFile_write().                                  *
 *                                                                        *
 * Return value: number of fields loaded if successful, EOF for end-of-   *
 *               file, 0 on error.                                        *
 * ---------------------------------------------------------------------- */

#define READLINE(fmt,arg) \
  if (fscanf(fp,fmt,arg)==EOF) return EOF; fgets(buf,BUFSIZ,fp)

int FieldFile_read (FieldFile *f, FILE *fp)
{
  char     buf[BUFSIZ];
  int i, n, ntotz, nfields;

  if (fscanf (fp, "%s", buf) == EOF) return EOF;          /* session name */
  f->name    = strdup (buf); fgets (buf, BUFSIZ, fp);
  fgets (buf, DESCRIP, fp);                               /* creation date */
  f->created = strdup (buf); fgets (buf, BUFSIZ, fp);
  
  /* Read the the file size parameters */

  if (fscanf (fp, "%d%d%d%d", &f->nr, &f->ns, &f->nz, &f->nel) != 4)
    return 0;
  fgets (buf, BUFSIZ, fp);

  READLINE ("%d" , &f->step);                  /* simulation parameters */
  READLINE ("%lf", &f->time);
  READLINE ("%lf", &f->time_step);
  READLINE ("%lf", &f->Re);
  READLINE ("%lf", &f->beta);
  
  nfields = gettypes (f->type, fgets(buf,BUFSIZ,fp));
  fscanf(fp, "%s", buf);
  if (strchr(f->format = strdup(buf),'-'))
    if (checkfmt (strchr(f->format,'-')+1))
      fputs ("cubit: warning: field file may be in the wrong "
	     "format for this architecture\n", stderr);
  fgets (buf, BUFSIZ, fp);

  /* Allocate memory for the data if it's not allocated already */

  ntotz = f->nr * f->ns * f->nz * f->nel;
  for (n = 0; n < nfields; n++)
    if (!(f->data [n]))
      f->data [n] = dvector (0, ntotz-1);

  /* Load the data */

  switch (tolower(*f->format)) {
  case 'a':
    for (i = 0; i < ntotz; i++) {
      for (n = 0; n < nfields; n++)
	if (fscanf (fp, "%lf", f->data[n] + i) != 1)
	  perror (CUBIT);
    }
    fgets (buf, BUFSIZ, fp);
    break;

  case 'b':
    for (n = 0; n < nfields; n++)
      if (fread (f->data[n], sizeof(double), ntotz, fp) != ntotz) {
	perror (CUBIT);
	return 0;
      }
    break;

  case 'p': {
    float tmp;
    for (n = 0; n < nfields; n++)
      for (i = 0; i < ntotz; i++) {
	if (fread (&tmp, sizeof(float), 1, fp) != 1) {
	  perror (CUBIT);
	  return 0;
	}
	f-> data[n][i] = (double) tmp;
      }
  }
    
  default:
    sprintf (cubit_err_msg, "unknown field format -- %s\n", f->format);
    cubit_err(NULL);
    break;
  }

  return nfields;
}

#undef READLINE

/* Transfer a field type array from s -> t.            *
 *                                                     *
 * NOTE: Assume that t is at least _MAX_FIELDS long!! */

static int gettypes (char *t, char *s)
{
  char    *p = strstr (s, "Fields");
  int n = 0;
  if (!p) cubit_err ("invalid \"Fields\" string");
  memset(t, '\0', _MAX_FIELDS);
  while (s != p && n < _MAX_FIELDS) {
    if (isalpha(*s)) t[n++] = *s;
    s++;
  }
  return n;
}


/* Check binary format compatibility */

#if (defined(sgi) || defined(cm5)) && !defined(ieee)
#define ieee
#endif

static int checkfmt (char *arch)
{
#if 1
  return 0;
#else
  char        **p;
  static char *fmtlist[] = {
#if defined(ieee)                    /* ..... Generic IEEE machine ..... */
    "ieee", 
#endif
#
#if defined(sgi)                    /* ........ Silicon Graphics ....... */
    "sgi", "iris4d", "SGI", "IRIX",
#endif
#
#if defined(cm5)                    /* ........ Connection Machine ..... */
    "cm5",
#endif
#
#if defined(i860)                   /* ........... Intel i860 .......... */
    "i860",
#endif
#
#if defined(_CRAY)                  /* ........... Cray Y-MP ........... */
    "cray", "CRAY",
#endif
#
     0 };   /* a NULL string pointer to signal the end of the list */

  for (p = fmtlist; *p; p++)
    if (strncmp (arch, *p, strlen(*p)) == 0)
      return 0;

  return 1;
#endif
}

#undef ieee

/* ---------------------------------------------------------------------- *
 * File_backup() -- Create a backup copy of a file                        *
 * ---------------------------------------------------------------------- */

extern int unlink (const char *path);
extern int link   (const char *path1, const char *path2);

int File_backup (char *path1)
{
  int  stat;
  char path2[FILENAME_MAX];

  sprintf (path2, "%s.bak", path1);
  unlink  (path2);                    /* unlink path2 regardless    */
  if (!(stat = link(path1, path2)))   /* try to link path1 -> path2 */
    unlink (path1);                   /* unlink path1 only if the   */
  return stat;                        /* link was sucessful         */
}
    
/* ---------------------------------------------------------------------- *
 * File_open() -- Open a file                                             *
 *                                                                        *
 * This is just a wrapper around fopen() that performs the proper synch-  *
 * ronization for a concurrent file (if necessary).                       *
 * ---------------------------------------------------------------------- */

FILE *File_open (char *path, char *mode)
{
  FILE *fp;

#if defined (i860) || defined (sim860)
  if (strstr (path, "/cfs"))
    if (mynode() == 0) {
      fp = fopen (path, mode);
      gsync();
    } else {
      gsync();
      fp = fopen (path, mode);
    }
  else
#endif
    fp = fopen (path, mode);

  return fp;
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

static double *data_interp (int nr, int ns, int nz, int nel, double *data,
			    int nx, int ny, double **imr, double **ims);

void FieldFile_project (FieldFile *f, int nr, int ns)
{
  double *zr, *zs;

  if (f->nr != nr || f->ns != ns) {
    coef (nr); getops (nr, & zr, 0, 0, 0);
    coef (ns); getops (ns, & zs, 0, 0, 0);

    FieldFile_interp (f, nr, zr, ns, zs);
  }

  return;
}

void FieldFile_interp (FieldFile *f, int nx, double *xx, int ny, double *yy)
{
  int      nr, ns, nz, nel, nfields;
  double   **imr, **itmr, **ims, **itms, *zr, *zs;
  double   **data;

  if (nx > _MAX_NORDER || ny > _MAX_NORDER) {
    sprintf (cubit_err_msg, 
	     "FieldFile_interp: Maximum order is %d", _MAX_NORDER);
    cubit_err (NULL);
    return;
  }

  nr      = f->nr;     
  ns      = f->ns;     
  nz      = f->nz;
  nel     = f->nel;
  data    = f->data;
  nfields = strlen(f->type);

  imr     = dmatrix (0, nx-1, 0, nr-1);
  itmr    = dmatrix (0, nr-1, 0, nx-1);
  ims     = dmatrix (0, ny-1, 0, ns-1);
  itms    = dmatrix (0, ns-1, 0, ny-1);

  /* Compute the GLL points */

  coef  (nr);  getops(nr, &zr, 0, 0, 0);
  coef  (ns);  getops(ns, &zs, 0, 0, 0);

  /* Compute the interpolation matrices */

  igllm (imr, itmr, zr, xx, nr, nx);
  igllm (ims, itms, zs, yy, ns, ny);

  /* Loop through the fields and interpolate each one */

  while (nfields--) 
    {
      double *p = 
	data_interp (nr, ns, nz, nel, *data, nx, ny, imr, ims);
      free (*data);
      *data++ = p;
    }

  f->nr = nx;
  f->ns = ny;

  free_dmatrix (imr , 0, 0);
  free_dmatrix (itmr, 0, 0);
  free_dmatrix (ims , 0, 0);
  free_dmatrix (itms, 0, 0);
  return;
}

/* Interpolation function */

static double 
 *data_interp (int nr, int ns, int nz, int nel, double *data,
	       int nx, int ny, double **imr, double **ims)
{
  const int   nrns = nr * ns,
          nxny = nx * ny;
  double *new  = dvector (0, nxny * nz * nel - 1), *p = new;
  int k;

  double tmp [_MAX_NORDER*_MAX_NORDER];

  for (k = 0; k < nel*nz; k++, data += nrns, p += nxny) {
    dgemm ('N', 'N', nr, ny, ns, 1., data, nr, *ims, ns, 0., tmp, nr);
    dgemm ('T', 'N', nx, ny, nr, 1., *imr, nr,  tmp, nr, 0.,   p, nx);
  }

  return new;
}

/* ---------------------------------------------------------------------- *
 * FieldFile_get() -- Methods to read information from a FieldFile        *
 * ---------------------------------------------------------------------- */

int FieldFile_getFieldCount (const FieldFile *f)
{ return strlen(f->type); }

char *FieldFile_getFieldList(const FieldFile *f, char *t)
{ return strcpy(t, f->type); }

char *FieldFile_getName (const FieldFile *f, char *s)
{ return strcpy(s, f->name); }

int FieldFile_getSize (const FieldFile *f, int *nr, int *ns, int *nz, int *nel)
{ 
  *nr  = f->nr;
  *ns  = f->ns;
  *nz  = f->nz;
  *nel = f->nel;
  return f->nr * f->ns * f->nz * f->nel;
}

/* ------------------------------------------------------------------------- *
 * FieldFile_store() -- Store one Field in a FieldFile                       *
 *                                                                           *
 * Note that FieldFile data is static and is layed out in x--y planes rather *
 * than blocks of elements.   This is to maintain compatibility with the     *
 * older fieldfile format.                                                   *
 * ------------------------------------------------------------------------- */

int FieldFile_store (FieldFile *f, const Field *u)
{
  const int nel  = Field_count(u);
  const int nrns = FIELD_NR(u)*FIELD_NS(u);
  const int nz   = FIELD_NZ(u);

  int i, k;

  char *p = strchr(f->type, FIELD_TYPE(u));
  int   n = strlen(f->type);

  if (f->nr == 0) {
    f->nr  = FIELD_NR(u);
    f->ns  = FIELD_NS(u);
    f->nz  = FIELD_NZ(u);
    f->nel = nel;
  }

  if (p) {
    fprintf (stderr, "warning: overwriting field %c in %s\n",
	     FIELD_TYPE(u), f->type);
    n = (int) (p - f->type);
  } else {
    f->type[n] = FIELD_TYPE(u);
    f->data[n] = dvector(0, nrns*nz*nel-1);
  }

  for (k = 0; k < nz; k++) {
    Element *elmt = FIELD_HEAD(u);
    int j = 0;
    while (elmt) {
      const int id = elmt->id;
      double *src  = FIELD_DATA(u,id) + k*nrns;
      double *dest = f->data[n] + j*nrns + k*nel*nrns;

      memcpy(dest, src, nrns*sizeof(double));

      elmt = elmt->next;
      j++;
    }
  }
  
  return 0;
}

int FieldFile_load (const FieldFile *f, Field *u)
{
  const int nel  = Field_count(u);
  const int nrns = FIELD_NR(u)*FIELD_NS(u);
  const int nz   = FIELD_NZ(u);

  char buf[BUFSIZ], *p;
  int  n;
  int  k;

  if (!(p = strchr(f->type, FIELD_TYPE(u)))) {
    sprintf (buf, "variable type %c not found in %c", FIELD_TYPE(u), f->type);
    cubit_err (buf);
  }

  if (f->nr  != FIELD_NR(u) ||
      f->ns  != FIELD_NS(u) ||
      f->nz  != FIELD_NZ(u) ||
      f->nel != FIELD_NELMT(u)) {
    sprintf (buf, "session and field file don't match: "
	     "[%d %d %d %d] vs. [%d %d %d %d]\n",
	     FIELD_NR(u), FIELD_NS(u), FIELD_NZ(u), nel,
	     f->nr,       f->ns,       f->nz,       f->nel);
    cubit_err (buf);
  }
      
  n = (int) (p - f->type);

  for (k = 0; k < nz; k++) {
    Element *elmt = FIELD_HEAD(u);
    int j = 0;
    while (elmt) {
      const int id = elmt->id;
      double *src  = f->data[n] + nrns*(k*nel+j);
      double *dest = FIELD_DATA(u,id) + k*nrns;

      memcpy(dest, src, nrns*sizeof(double));
      elmt = elmt->next;
      j++;
    }
  }

  return 0;
}

/* --------------------------------------------------------------------- *
 * FieldFile_setName() -- Set the session name for a FieldFile           *
 * --------------------------------------------------------------------- */

char *FieldFile_setName (FieldFile *f, const char *s)
{
  if (f->name)
    free (f->name);

  f->name = strdup(s);
  return f->name;
}

  
