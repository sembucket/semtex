#ifndef FIELDFILE_H
#define FIELDFILE_H

#include <stdio.h>
#include "field.h"

#define FIELDFILE_NAME(f)    (f->name)
#define FIELDFILE_NR(f)      (f->nr)
#define FIELDFILE_NS(f)      (f->ns)
#define FIELDFILE_NZ(f)      (f->nz)
#define FIELDFILE_NELMT(f)   (f->nel)
#define FIELDFILE_STEP(f)    (f->step)
#define FIELDFILE_TIME(f)    (f->time)
#define FIELDFILE_DT(f)      (f->time_step)
#define FIELDFILE_RE(f)      (f->Re)
#define FIELDFILE_BETA(f)    (f->beta)
#define FIELDFILE_NFLDS(f)   (FieldFile_getFieldCount(f))

typedef struct fieldfile {   /* ..... Field FILE structure ...... */
  char*   name           ;   /* session name                      */
  char*   created        ;   /* date and time created             */
  int     nr, ns, nz, nel;   /* mesh size parameters              */

  void*   misc           ;   /* miscellaneous                     */

  int     step           ;   /* --- Got to go!                    */
  double  time           ;
  double  time_step      ;
  double  Re, beta       ;

  char*   format         ;   /* file format string                */
  char    type[_MAX_FIELDS]; /* list of variable types            */
  double* data[_MAX_FIELDS]; /* array of data values              */
} FieldFile;                 /* --------------------------------- */

/* Prototypes */

FieldFile* FieldFile_alloc  (void);
void       FieldFile_free   (FieldFile *f);
int        FieldFile_write  (FieldFile *f, FILE *fp);
void       FieldFile_header (FieldFile *f, FILE *fp);
int        FieldFile_read   (FieldFile *f, FILE *fp);
void       FieldFile_project(FieldFile *f, int nr, int ns);
void       FieldFile_interp (FieldFile *f, int nx, double *x, 
			                   int ny, double *y);
int        File_backup      (char *path);
FILE*      File_open        (char *path, char *mode);

int        FieldFile_getFieldCount (const FieldFile *f);
char*      FieldFile_getFieldList  (const FieldFile *f, char *t);
char*      FieldFile_getName       (const FieldFile *f, char *s);
int        FieldFile_getSize       (const FieldFile *f, int*,int*,int*,int*);

char*      FieldFile_setName (FieldFile *f, const char *s);

int        FieldFile_load  (const FieldFile *f, Field *U);
int        FieldFile_store (FieldFile *f, const Field *U);

#endif
