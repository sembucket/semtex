/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY LIBFEM.A              *
 *****************************************************************************/

/*------------------*
 * RCS Information: *
 *------------------*/

/* $Id$ */



#ifndef FEMLIB_H		/* BEGIN femlib.h DECLARATIONS */
#define FEMLIB_H

#include <stdio.h>


#define STANDARD   'f' 	/* "Standard" finite element basis functions. */
#define GLL        's'	/* Gauss-Lobatto-Legendre basis functions.    */

/* ------------------------------------------------------------------------- *
 * Routines from initial.c:                                                  *
 * ------------------------------------------------------------------------- */

void   initialize  (void);
double interpret   (char *);

void   setOption  (char *, int);
int    getOption  (char *);
void   showOption (FILE *);

void   setIparam  (char *, int);
int    getIparam  (char *);
void   showIparam (FILE *);

void   setDparam  (char *, double);
double getDparam  (char *);
void   showDparam (FILE *);

/* ------------------------------------------------------------------------- *
 * Routines from polyops.c:                                                  *
 * ------------------------------------------------------------------------- */

void   dermat_g (int K, double *zero, int I, double  *x,
		                      double **D,  double **DT);
void   intmat_g (int K, double *zero, int I, double  *x, 
		                      double **IN, double **IT);
void   dermat_k (int K, double *zero, double **D,  double **DT);

void   jacg     (int n, double alpha, double beta, double *xjac);
void   jacgr    (int n, double alpha, double beta, double *xjac);
void   jacgl    (int n, double alpha, double beta, double *xjac);

void   zwgl     (double *z, double *w, int np);
void   zwgrl    (double *z, double *w, int np);
void   zwgll    (double *z, double *w, int np);

double pnleg    (double  z, int n);
double pndleg   (double  z, int n);
double pnd2leg  (double  z, int n);

void   jacobf   (int     n     ,  double  x     ,
		 double  alpha ,  double  beta  ,
		 double *poly  ,  double *pder  ,
		 double *polym1,  double *pderm1,
		 double *polym2,  double *pderm2);

void   dgll     (int nz, double *z, double **D, double **DT);

void   uniknot  (int nk, double *k);

/* ------------------------------------------------------------------------- *
 * Routines from operators.c:                                                *
 * ------------------------------------------------------------------------- */

void  quadOps(int basis ,	/* input: element basis: STANDARD or GLL     */
	      int np    ,	/* input: number of knot points              */
	      int nq    ,	/* input: number of quadrature points        */
	      double  **kp ,	/* pointer to knot point storage             */
	      double  **qp ,	/* pointer to quadrature point storage       */
	      double  **qw ,	/* pointer to quadrature weight storage      */
	      double ***in ,	/* pointer to interpolation matrix           */
	      double ***it ,	/* pointer to transposed interpolant matrix  */
	      double ***dr ,	/* pointer to derivative matrix              */
	      double ***dt );	/* pointer to transposed derivative matrix   */

void  meshOps(int oldbasis   ,	/* input: element basis: STANDARD or GLL     */
	      int newbasis   ,	/* input: desired basis: STANDARD or GLL     */
	      int np         ,	/* input: number of knot points              */
	      int ni         ,	/* input: number of interpolant points       */
	      double  **mesh ,	/* pointer to interpolated mesh storage      */
	      double ***in   ,	/* pointer to interpolation matrix           */
	      double ***it   ,	/* pointer to transposed interpolant matrix  */
	      double ***dr   ,	/* pointer to derivative matrix              */
	      double ***dt   ); /* pointer to transposed derivative matrix   */

#endif
