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

void   dermat_g (const int, const double*, const int,
		 const double*, double**, double**);
void   intmat_g (const int, const double*, const int,
		 const double*, double**, double**);
void   dermat_k (const int, const double*, double**, double**);

void   jacg     (const int, const double, const double, double*);
void   jacgr    (const int, const double, const double, double*);
void   jacgl    (const int, const double, const double, double*);

void   zwgl     (double*, double*, const int);
void   zwgrl    (double*, double*, const int);
void   zwgll    (double*, double*, const int);

double pnleg    (const double, const int);
double pndleg   (const double, const int);
double pnd2leg  (const double, const int);

void   dgll     (const int, const double*, double**, double**);

void   uniknot  (const int, double*);

int    quadComplete (const int, const int);

/* -- Routines from operators.c: */

void  dQuadOps(const int       rule, /* input: quadrature rule: GL or LL     */
	       const int       np  , /* input: number of knot points         */
	       const int       nq  , /* input: number of quadrature points   */
	       const double**  kp  , /* pointer to knot point storage        */
	       const double**  qp  , /* pointer to quadrature point storage  */
	       const double**  qw  , /* pointer to quadrature weight storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  );/* pointer to transposed deriv matrix   */

void  dMeshOps(const int       old , /* input: element basis: STD or GLL     */
	       const int       new , /* input: desired basis: STD or GLL     */
	       const int       np  , /* input: number of knot points         */
	       const int       ni  , /* input: number of interpolant points  */
	       const double**  mesh, /* pointer to interpolated mesh storage */
	       const double*** in  , /* pointer to interpolation matrix      */
	       const double*** it  , /* pointer to transposed interp matrix  */
	       const double*** dr  , /* pointer to derivative matrix         */
	       const double*** dt  );/* pointer to transposed deriv matrix   */

void  sQuadOps(const int       rule, /* input: quadrature rule: GL or LL     */
	       const int       np  , /* input: number of knot points         */
	       const int       nq  , /* input: number of quadrature points   */
	       const float**   kp  , /* pointer to knot point storage        */
	       const float**   qp  , /* pointer to quadrature point storage  */
	       const float**   qw  , /* pointer to quadrature weight storage */
	       const float***  in  , /* pointer to interpolation matrix      */
	       const float***  it  , /* pointer to transposed interp matrix  */
	       const float***  dr  , /* pointer to derivative matrix         */
	       const float***  dt  );/* pointer to transposed deriv matrix   */

void  sMeshOps(const int       old , /* input: element basis: STD or GLL     */
	       const int       new , /* input: desired basis: STD or GLL     */
	       const int       np  , /* input: number of knot points         */
	       const int       ni  , /* input: number of interpolant points  */
	       const float**   mesh, /* pointer to interpolated mesh storage */
	       const float***  in  , /* pointer to interpolation matrix      */
	       const float***  it  , /* pointer to transposed interp matrix  */
	       const float***  dr  , /* pointer to derivative matrix         */
	       const float***  dt  );/* pointer to transposed deriv matrix   */

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
