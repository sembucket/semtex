#ifndef femlibH
#define femlibH
/*****************************************************************************
 *         FUNCTION PROTOTYPES FOR ROUTINES IN LIBRARY LIBFEM.A
 *****************************************************************************/

/* $Id$ */

/* -- Routines from initial.c: */

void   yy_initialize (void);
double yy_interpret  (const char*);

void   yy_vec_init   (const char*, const char*);
void   yy_vec_interp (const int, ...);

void   yy_help       (void);
void   yy_show       (void);
int    yy_dump       (char*, const int);

/* -- Routines from polyops.c: */

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

int    quadComplete (int dim, int np);

/* -- Routines from operators.c: */

void  dQuadOps(const int rule,	/* input: quadrature rule: GL or LL          */
	       const int np  ,	/* input: number of knot points              */
	       const int nq  ,	/* input: number of quadrature points        */
	       double  **kp  ,	/* pointer to knot point storage             */
	       double  **qp  ,	/* pointer to quadrature point storage       */
	       double  **qw  ,	/* pointer to quadrature weight storage      */
	       double ***in  ,	/* pointer to interpolation matrix           */
	       double ***it  ,	/* pointer to transposed interpolant matrix  */
	       double ***dr  ,	/* pointer to derivative matrix              */
	       double ***dt  );	/* pointer to transposed derivative matrix   */

void  dMeshOps(const int old ,	/* input: element basis: STD or GLL          */
	       const int new ,	/* input: desired basis: STD or GLL          */
	       const int np  ,	/* input: number of knot points              */
	       const int ni  ,	/* input: number of interpolant points       */
	       double  **mesh,	/* pointer to interpolated mesh storage      */
	       double ***in  ,	/* pointer to interpolation matrix           */
	       double ***it  ,	/* pointer to transposed interpolant matrix  */
	       double ***dr  ,	/* pointer to derivative matrix              */
	       double ***dt  ); /* pointer to transposed derivative matrix   */

void  sQuadOps(const int rule,	/* input: quadrature rule: GL or LL          */
	       const int np  ,	/* input: number of knot points              */
	       const int nq  ,	/* input: number of quadrature points        */
	       float   **kp  ,	/* pointer to knot point storage             */
	       float   **qp  ,	/* pointer to quadrature point storage       */
	       float   **qw  ,	/* pointer to quadrature weight storage      */
	       float  ***in  ,	/* pointer to interpolation matrix           */
	       float  ***it  ,	/* pointer to transposed interpolant matrix  */
	       float  ***dr  ,	/* pointer to derivative matrix              */
	       float  ***dt  );	/* pointer to transposed derivative matrix   */

void  sMeshOps(const int old ,	/* input: element basis: STD or GLL          */
	       const int new ,	/* input: desired basis: STD or GLL          */
	       const int np  ,	/* input: number of knot points              */
	       const int ni  ,	/* input: number of interpolant points       */
	       float   **mesh,	/* pointer to interpolated mesh storage      */
	       float  ***in  ,	/* pointer to interpolation matrix           */
	       float  ***it  ,	/* pointer to transposed interpolant matrix  */
	       float  ***dr  ,	/* pointer to derivative matrix              */
	       float  ***dt  ); /* pointer to transposed derivative matrix   */

/* -- Routines from mapping.c: */

void edgemaps (const int np, int** emap, int** pmap);

/* -- Routines from family.c: */

void iadopt   (const int, int**   );
void dadopt   (const int, double**);
void sadopt   (const int, float** );

void iabandon (int**   );
void dabandon (double**);
void sabandon (float** );

int  FamilySize (int*, int*, int*);

/* -- Routines from RCM.f: */

void genrcm_ (int*, int*, int*, int*, int*, int*);
#define genrcm(neqns, xadj, adjncy, perm, mask, xls) \
(_alpIreg[0] = neqns, genrcm_(_alpIreg, xadj, adjncy, perm, mask, xls))

void fnroot_ (int*, int*, int*, int*, int*, int*, int*);
#define fnroot(root, xadj, adncy, mask, nlvl, xls, ls) \
(_alpIreg[0] = root, _alpIreg[1] = nlvl,                   \
 fnroot_(_alpIreg, xadj, adjncy, mask, _alpIreg + 1, xls, ls))

void rcm_ (int*, int*, int*, int*, int*, int*, int*);
#define rcm(root, xadj, adjncy, mask, perm, ccsize, deg)  \
(_alpIreg[0] = root, rcm_(_alpIreg, xadj, adjncy, mask, perm, ccsize, deg))

/* -- Routines from fftpack.f (NETLIB/FFTPACK): */

void srffti_ (int*, float*);
void srfftf_ (int*, float*, float*);
void srfftb_ (int*, float*, float*);

#define srffti(n,wsave)   (_alpIreg[0]=n, srffti_ (_alpIreg, wsave))
#define srfftf(n,r,wsave) (_alpIreg[0]=n, srfftf_ (_alpIreg, r, wsave))
#define srfftb(n,r,wsave) (_alpIreg[0]=n, srfftb_ (_alpIreg, r, wsave))

void drffti_ (int*, double*);
void drfftf_ (int*, double*, double*);
void drfftb_ (int*, double*, double*);

#define drffti(n,wsave)   (_alpIreg[0]=n, drffti_ (_alpIreg, wsave))
#define drfftf(n,r,wsave) (_alpIreg[0]=n, drfftf_ (_alpIreg, r, wsave))
#define drfftb(n,r,wsave) (_alpIreg[0]=n, drfftb_ (_alpIreg, r, wsave))

/* -- Routines from fourier.c */

void sDFTr (float*,  const int, const int, const int, const int, const int);
void dDFTr (double*, const int, const int, const int, const int, const int);

#endif
