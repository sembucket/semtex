#ifndef FIELDFILE_H
#define FIELDFILE_H

/* FieldFile
 *
 * $Id$
 * 
 * Copyright (c) 1998 R. D. Henderson and Caltech
 *
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>

#include "speclib/field.h"

#define FIELDFILE_MAX      31
#define FIELDFILE_STRSIZE  32
#define FIELDFILE_EOF      -1

#define FIELDFILE_NR(f)     (f->nr)
#define FIELDFILE_NS(f)     (f->ns)
#define FIELDFILE_NZ(f)     (f->nz)
#define FIELDFILE_NEL(f)    (f->nel)
#define FIELDFILE_STEP(f)   (f->step)
#define FIELDFILE_TIME(f)   (f->time)
#define FIELDFILE_DT(f)     (f->time_step)
#define FIELDFILE_RE(f)     (f->Re)
#define FIELDFILE_BETA(f)   (f->beta)
#define FIELDFILE_NAME(f)   (f->name)

#define FIELDFILE_COUNT(f)  (f->count)
#define FIELDFILE_TYPE(f,n) (f->type[n])
#define FIELDFILE_DATA(f,n) (f->data[n])

typedef enum format {
  BINARY,
  PACKED,
  ASCII
} format_t;

typedef struct fieldfile {   /* ..... Field FILE structure ...... */
  int      count           ; /* number of stored fields           */
  int      nr, ns, nz, nel ; /* mesh size parameters              */
  int      step            ; /* time step number                  */
  double   time            ; /* time                              */
  double   time_step       ; /* time step                         */
  double   Re              ; /* Reynolds number                   */
  double   beta            ; /* fundamental wavenumber for 3D     */

  /* All arrays must be fixed-length so that a FieldFile can be   *
   * broadcast to other processors w/out packing and unpacking.   */

  char     name    [FIELDFILE_STRSIZE];    /* session name        */
  char     created [FIELDFILE_STRSIZE];    /* time stamp          */
  char     format  [FIELDFILE_STRSIZE];    /* data format string  */
  char     type    [FIELDFILE_STRSIZE];    /* stored field types  */
  double  *data    [FIELDFILE_MAX];        /* stored field data   */
} FieldFile;                 /* --------------------------------- */

/* Prototypes */

FieldFile* FieldFile_alloc  ();
void       FieldFile_free   (FieldFile *f);

int FieldFile_setFormat (FieldFile *f, format_t format);
int FieldFile_setSize   (FieldFile *f, int count, int nr, int ns, int nz, int nel);
char* FieldFile_setName (FieldFile *f, const char *name);
char* FieldFile_getName (const FieldFile *f, char *name);

int FieldFile_get    (const FieldFile *f, Field *u);
int FieldFile_put    (FieldFile *f, const Field *u);    
int FieldFile_write  (FieldFile *f, FILE *fp);
int FieldFile_header (FieldFile *f, FILE *fp);
int FieldFile_read   (FieldFile *f, FILE *fp);

void FieldFile_project  (FieldFile *f, int nr, int ns);
void FieldFile_projectz (FieldFile *f, int nz);
void FieldFile_interp   (FieldFile *f, int nx, double *x, int ny, double *y);

/* Info */

int FieldFile_getTypeList (const FieldFile *f, char *list);
int FieldFile_checkType   (const FieldFile *f, char type);

#endif

