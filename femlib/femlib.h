/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY LIBFEM.A
 *****************************************************************************/

/* $Id$ */



#ifndef femlibH
#define femlibH

/* ------------------------------------------------------------------------- *
 * Routines from initial.c:                                                  *
 * ------------------------------------------------------------------------- */

void   initialize ();
double interpret  (const char *);

void   vecInit    (const char*, const char *);
void   vecInterp  (int, ...);

void   setOption  (const char *, int);
int    option     (const char *);
void   showOption ();

void   setIparam  (const char *, int);
int    iparam     (const char *);
void   showIparam ();

void   setDparam  (const char *, double);
double dparam     (const char *);
void   showDparam ();

/* ------------------------------------------------------------------------- *
 * Routines from polyops.c:                                                  *
 * ------------------------------------------------------------------------- */

void   dermat_g (int K, double *zero, int I, double  *x,
		                      double **DV, double **DT);
void   intmat_g (int K, double *zero, int I, double  *x, 
		                      double **IN, double **IT);
void   dermat_k (int K, double *zero, double **DV, double **DT);

void   jacg     (int n, double alpha, double beta, double *xjac);
void   jacgr    (int n, double alpha, double beta, double *xjac);
void   jacgl    (int n, double alpha, double beta, double *xjac);

void   zwgl     (double *z, double *w, int np);
void   zwgrl    (double *z, double *w, int np);
void   zwgll    (double *z, double *w, int np);

double pnleg    (double  z, int n);
double pndleg   (double  z, int n);
double pnd2leg  (double  z, int n);

void   dgll     (int nz, double *z, double **D, double **DT);

void   uniknot  (int nk, double *k);

/* ------------------------------------------------------------------------- *
 * Routines from operators.c:                                                *
 * ------------------------------------------------------------------------- */

void  quadOps(int rule     ,	/* input: quadrature rule: GL or LL          */
	      int np       ,	/* input: number of knot points              */
	      int nq       ,	/* input: number of quadrature points        */
	      double  **kp ,	/* pointer to knot point storage             */
	      double  **qp ,	/* pointer to quadrature point storage       */
	      double  **qw ,	/* pointer to quadrature weight storage      */
	      double ***in ,	/* pointer to interpolation matrix           */
	      double ***it ,	/* pointer to transposed interpolant matrix  */
	      double ***dr ,	/* pointer to derivative matrix              */
	      double ***dt );	/* pointer to transposed derivative matrix   */

void  meshOps(int oldbasis   ,	/* input: element basis: STD or GLL          */
	      int newbasis   ,	/* input: desired basis: STD or GLL          */
	      int np         ,	/* input: number of knot points              */
	      int ni         ,	/* input: number of interpolant points       */
	      double  **mesh ,	/* pointer to interpolated mesh storage      */
	      double ***in   ,	/* pointer to interpolation matrix           */
	      double ***it   ,	/* pointer to transposed interpolant matrix  */
	      double ***dr   ,	/* pointer to derivative matrix              */
	      double ***dt   ); /* pointer to transposed derivative matrix   */

/* ------------------------------------------------------------------------- *
 * Routines from RCM.f:                                                      *
 * ------------------------------------------------------------------------- */

extern void genrcm_ (int*, int*, int*, int*, int*, int*);

#define genrcm(neqns, xadj, adjncy, perm, mask, xls) \
(_alpIreg[0] = neqns, genrcm_(_alpIreg, xadj, adjncy, perm, mask, xls))



#endif
