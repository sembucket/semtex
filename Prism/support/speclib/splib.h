#ifndef SPLIB_H
#define SPLIB_H

/* SPLIB
 * 
 * $Id$
 *
 * Copyright (c) 1998 R. D. Henderson and Caltech
 * 
 * ------------------------------------------------------------------------- */

/* Points and weights */

void   zwgl    (double*,double*,int);
void   zwgll   (double*,double*,int);
void   zwgrl   (double*,double*,int);
void   jacgl   (int,double,double,double*);

/* Polynomial evaluation */

void   jacobf  (int,double*,double*,double*,double*,double*,double*,double);
void   jacgr   (int,double,double,double*);
void   jacg    (int,double,double,double*);
double pnleg   (double zp, int nz);
double pndleg  (double zp, int nz);
double pnd2leg (double zp, int nz);
double pnddleg (double zp, int nz);

/* Lagrangian interpolants */

double hgl     (int i, double zp, double *z, int nz);
double hgll    (int i, double zp, double *z, int nz);

/* Derivative operators */

void   dgll    (double**,double**,double*,int);

/* Interpolation operators */

void   genim12 (double **,double**,double*,double*,int,int);
void   genim21 (double **,double**,double*,double*,int,int);
void   genim13 (double **,double**,double*,double*,int,int);
void   genim31 (double **,double**,double*,double*,int,int);
void   igllm   (double **im12, double **im12t, 
		     double *zgll, double *zm, int nz, int mz);
		
#endif
