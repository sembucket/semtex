/*
 * Utilities for PRISM
 *
 * RCS Information
 * ------------------------------------
 * $Author$
 * $Date$
 * $RCSfile$
 * $Revision$
 * ---------------------------------------------------------------------- */

#ifndef   SEM_UTILS_H
#define   SEM_UTILS_H

#include <stdio.h>
#include <math.h>

#include "veclib/veclib.h"
#include "speclib/speclib.h"
#include "mesh.h"

#define MAXARGS  16            /* maximum # of other arguments */
#define UNSET    0xffffff      /* flag for an unset parameter  */

extern char *prog;             
extern char *usage, *help,     /* defined by the application */
            *rcsid, *author;          

typedef struct {               /* structure for storing fp/names */
  FILE *fp;
  char *name;
} aFile;

typedef struct {             /* .......  File structure ........ */
  aFile  in      ;           /* primary input file               */
  aFile  out     ;           /* primary output file              */
  aFile  mesh    ;           /* mesh coordinates file [optional] */
  aFile  rea     ;           /* session file          [optional] */
  aFile  mor     ;           /* connectivity file     [optional] */
} FileList;                  /* ................................ */


/* Macros for reading the argument list */

#define   nextargi(n,argc,argv)\
        { if (*++argv[0]) n = atoi(*argv); \
          else { n = atoi(*++argv); argc--; } (*argv)[1] = '\0'; }
#define   nextargf(f,argc,argv)\
        { if (*++argv[0]) f = atof(*argv);\
          else { f = atof(*++argv); argc--; } (*argv)[1] = '\0'; }
#define   nextargc(c,argc,argv)\
        { if (*++argv[0]) c = *(*argv); else { c = *(*++argv); argc--; } }
#define   nextargs(s,argc,argv)\
        { if (*++argv[0]) s = *argv; else { s = *++argv; argc--; }\
	  *argv += strlen(s)-1; } 

/* Prototypes, listed by source file */


/* Functions from sem_utils.c */

int      generic_args   (int argc, char *argv[], FileList *f);
void     error_msg      (char *msg);
Element* load_aux_field (Element *U, char type, FieldFile *f);
int      GetRS          (Element *U, double x, double y, double *r, double *s);
int      lookupF        (Element *U, double x, double y, double *r, double *s);
double   interpF        (int nr, double *hr, int ns, double *hs, double *data);

#endif
