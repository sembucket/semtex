/*****************************************************************************
 * PRESSURE.C: routines to deal with pressure field and boundary conditions. *
 *****************************************************************************/

static char
RCSid[] = "$Id$";


#include "Fem.h"
#include "integration.h"


typedef struct  hobc {		/* -------- High order PBC storage --------- */
  int           id   ;		/* A structure matching each Bedge (same id) */
  double      **Px   ;		/* dP / dx at domain boundary                */
  double      **Py   ;		/* dP / dy                                   */
  double      **Ux   ;		/* Velocity storage             (x-compt)    */
  double      **Uy   ;		/*                              (y-compt)    */
  struct hobc  *next ;		/* Link to next one                          */
} HOBC;				/* ----------------------------------------- */

static HOBC *BCList = NULL;


static void curlCurl (int            np    ,  int            side,
		      int            estart,  int            skip,
		      const Element *u     ,  const Element *v   ,
		      double        *wx    ,  double        *wy  );





Bedge *buildPBC (const Bedge *B, const Element *E)
/* ========================================================================= *
 * Build extra file-scope storage structures required for high order PBCs,   *
 * return a Bedge structure for pressure.                                    *
 * Input is pointer to a list of velocity Bedges.                            *
 *                                                                           *
 * Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting     *
 * methods for the incompressible Navier--Stokes equations", JCP V9N2.       *
 * ========================================================================= */
{
  char      routine[] = "buildPBC";

  Bedge    *P, *p;
  BcRecord *DefaultType;

  int       nOrder    = iparam ("TORDER");
  int       nEdge     = countBedge (B);
  int       np        = B -> np;
  int       ntot      = nEdge * np;
  int       i, q;



  P = copyBedge (B, (Element *)E);

  DefaultType = (BcRecord *) calloc (1, sizeof (BcRecord));
  DefaultType -> kind = HOPBC;
    
  if (!(BCList = (HOBC *) calloc (nEdge, sizeof (HOBC))))
    message (routine, "allocation of HOBC list failed", ERROR);

  BCList -> Px = dmatrix (0, nOrder-1, 0, ntot-1);
  BCList -> Py = dmatrix (0, nOrder-1, 0, ntot-1);
  BCList -> Ux = dmatrix (0, nOrder-1, 0, ntot-1);
  BCList -> Uy = dmatrix (0, nOrder-1, 0, ntot-1);

  if (!(BCList -> Px && BCList -> Py && BCList -> Ux && BCList -> Uy))
    message (routine, "allocation of HOBC storage failed", ERROR);

  for (p = P, i = 0; i < nEdge; i++, p = p -> next, B = B -> next) {

    if (i != nEdge - 1) (BCList + i) -> next = BCList + i + 1;

#if 1
    /* -- use "conventional" pressure outflow BC (0.0). */
    if (B -> bc -> kind == OUTFLOW) {
      p -> bc = (BcRecord *) calloc (1, sizeof (BcRecord));
      p -> bc -> kind = ESSENTIAL;
    } else
#endif
      p -> bc = DefaultType;

    (BCList + i) -> id = p -> id;

    if (i > 0) {
      (BCList + i) -> Px = (double **) calloc (nOrder, sizeof (double *));
      (BCList + i) -> Py = (double **) calloc (nOrder, sizeof (double *));
      (BCList + i) -> Ux = (double **) calloc (nOrder, sizeof (double *));
      (BCList + i) -> Uy = (double **) calloc (nOrder, sizeof (double *));

      for (q = 0; q < nOrder; q++) {
	(BCList + i) -> Px[q]  = BCList -> Px[q] + i * np;
	(BCList + i) -> Py[q]  = BCList -> Py[q] + i * np;
	(BCList + i) -> Ux[q]  = BCList -> Ux[q] + i * np;
	(BCList + i) -> Uy[q]  = BCList -> Uy[q] + i * np;
      }
    }
  }

  return P;
}





void  maintainPBC (const Bedge     *B       , 
		         int        step    ,
		         Element  **Us[DIM] ,
		         Element  **Uf[DIM] )
/* ========================================================================= *
 * Update storage for evaluation of high-order pressure boundary condition.  *
 *                                                                           *
 * No smoothing is done to high-order spatial derivatives computed here.     *
 *                                                                           *
 * We add estimates of velocity local acceleration to the linear terms;      *
 * since these must be constructed from velocity fields, there must be the   *
 * same amount of storage as for the time order of the scheme, (since, e.g.  *
 * for a first order scheme we need two levels of velocities to estimate a   *
 * time derivative: these come from the new one passed in, and the boundary  *
 * store).  Note also that this term cannot be estimated on the first step.  *
 * ========================================================================= */
{
  int      np = B -> np;
  HOBC    *p;
  double  *tmp;
  double   nu = dparam ("KINVIS");

  int      q, nOrder, Je = iparam ("TORDER");
  double   alpha[TIME_ORDER_MAX];
  double   gamma;
  double   invDt = 1.0 / dparam ("DELTAT");

  const Element *Ux = Us[0][0];
  const Element *Uy = Us[1][0];
  const Element *Nx = Uf[0][0];
  const Element *Ny = Uf[1][0];

  double *wx = dvector(0, np-1);
  double *wy = dvector(0, np-1);


  nOrder = Je;
  
  for (p = BCList; B; B = B -> next, p = p -> next) {

    /* -- Roll grad P storage area up, load new level of nonlinear terms. */

    rollUp (p -> Px, nOrder, tmp);
    Veclib::copy (np,*(Nx + B->elmt->id)->value+B->estart,B->eskip,*p -> Px,1);

    rollUp (p -> Py, nOrder, tmp);
    Veclib::copy (np,*(Ny + B->elmt->id)->value+B->estart,B->eskip,*p -> Py,1);

    /* -- Add in -nu * curl curl u. */

    curlCurl (np, B->side, B->estart, B->eskip,
	      Ux + B->elmt->id, Uy + B->elmt->id, wx, wy);

    Blas::axpy (np, -nu, wy, 1, *p -> Px, 1);
    Blas::axpy (np,  nu, wx, 1, *p -> Py, 1);

    /* -- Estimate -du/dt by backwards differentiation, add in. */

    if (step > 1 && B -> bc -> kind == HOPBC) {

      Je    =    MIN (step - 1, Je);
      tmp   =    wx;
      gamma =    Icoef[Je - 1].gamma;
      Veclib::copy (Je, Icoef[Je - 1].alpha, 1, alpha, 1);

      Veclib::copy (np,*(Ux + B->elmt->id)->value + B->estart, B->eskip,tmp,1);
      Blas::scal (np, gamma, tmp, 1);
      for (q = 0; q < Je; q++)
	Blas::axpy (np, -alpha[q], *(p -> Ux + q), 1, tmp, 1);
      Blas::axpy (np, -invDt, tmp, 1, *p -> Px, 1);

      Veclib::copy (np,*(Uy + B->elmt->id)->value + B->estart, B->eskip,tmp,1);
      Blas::scal (np, gamma, tmp, 1);
      for (q = 0; q < Je; q++)
	Blas::axpy (np, -alpha[q], *(p -> Uy + q), 1, tmp, 1);
      Blas::axpy (np, -invDt, tmp, 1, *p -> Py, 1);
    }

    /* -- Roll velocity storage area up, load new level. */

    rollUp (p -> Ux, nOrder, tmp);
    Veclib::copy  (np, *(Ux+B->elmt->id)->value+B->estart, B->eskip, *p -> Ux, 1);

    rollUp (p -> Uy, nOrder, tmp);
    Veclib::copy  (np, *(Uy+B->elmt->id)->value+B->estart, B->eskip, *p -> Uy, 1);
  }

  freeDvector (wx, 0);
  freeDvector (wy, 0);    
}





void  evaluatePBC (Bedge *B, int step)
/* ========================================================================= *
 * Load PBC storage with values obtained from HOBC multi-level storage.      *
 *                                                                           *
 * The boundary condition for evaluation is                                  *
 *                                                                           *
 *   dP       /                           du  \                              *
 *   -- = n . | N(u) - a + f + \nu*L(u) - --  |  =  n . grad P.              *
 *   dn   ~   \ ~ ~    ~   ~       ~ ~    dt  /     ~                        *
 *                                                                           *
 * Grad P is estimated at the end of the current timestep using explicit     *
 * extrapolation.                                                            *
 * ========================================================================= */
{
  int      np = B -> np;
  HOBC    *p;

  int      q, Je = iparam ("TORDER");
  double   beta[TIME_ORDER_MAX];
  
  double  *tmpX = dvector(0, np-1);
  double  *tmpY = dvector(0, np-1);


  Je = MIN (step, Je);
  Veclib::copy (Je, Icoef[Je - 1].beta, 1, beta, 1);

  for (p = BCList; B; B = B -> next, p = p -> next) {

    if (B -> bc -> kind == ESSENTIAL) {
      Veclib::zero (B -> np, B -> value, 1);
      continue;
    }

    np = B -> np;
    Veclib::zero (np, tmpX, 1);
    Veclib::zero (np, tmpY, 1);

    for (q = 0; q < Je; q++) {
      Blas::axpy (np, beta[q], *(p -> Px + q), 1, tmpX, 1);
      Blas::axpy (np, beta[q], *(p -> Py + q), 1, tmpY, 1);
    }
    
    Veclib::vmul  (np, B -> nx, 1, tmpX, 1, B -> value, 1);
    Veclib::vvtvp (np, B -> ny, 1, tmpY, 1, B -> value, 1, B -> value, 1);
  }

  freeDvector (tmpX, 0);
  freeDvector (tmpY, 0);
}





void  rollStacks (int stacksize, Element **Us[DIM],  Element **Uf[DIM])
/* ========================================================================= *
 * Service routine to maintain storage stacks.                               *
 * ========================================================================= */
{
  int       i;
  Element  *tmp;


  for (i = 0; i < DIM; i++) {
    rollUp (Us[i], stacksize, tmp);
    rollUp (Uf[i], stacksize, tmp);
  }
}





void setPForce (Domain *D, Element **Us[DIM], Element **Uf[DIM])
/* ========================================================================= *
 * On input, intermediate velocity storage u^ is in D.  Transfer this to     *
 * the first levels of Us and create div u^ / DELTAT in the first dimension, *
 * first level storage of Uf as a forcing field for discrete PPE.            *
 * ========================================================================= */
{
  int       i;
  Element  *tmp;
  double    invDt = 1.0 / dparam ("DELTAT");
  int       ntot  = D -> nEl * SQR (D -> u[0] -> np);


  for (i = 0; i < DIM; i++) {
    tmp       = D -> u[i];
    D -> u[i] = Us[i][0];
    Us[i][0]  = tmp;

    Veclib::copy (ntot, *Us[i][0] -> value, 1, *Uf[i][0] -> value, 1);
  }

  grad  (Uf[0][0], Uf[1][0]);
  Veclib::svvpt (ntot, invDt,
	 *Uf[0][0] -> value, 1, *Uf[1][0] -> value, 1, *Uf[0][0] -> value, 1);

  dsSmooth (Uf[0][0]);
}




void  constrain (Domain *D, Element **Us[DIM], Element **Uf[DIM])
/* ========================================================================= *
 * On input, new pressure field is stored in D and intermediate velocity     *
 * level u^ is stored in lowest level of Us.  Constrain velocity field:      *
 *                    u^^ = u^ - DELTAT * grad P;                            *
 * u^^ is left in lowest level of Us.                                        *
 * ========================================================================= */
{
  int     i;
  int     ntot = D -> nEl * SQR (D -> u[0] -> np);
  double  dt   = dparam ("DELTAT");


  for (i = 0; i < DIM; i++)
    Veclib::copy (ntot, *D -> u[DIM] -> value, 1, *Uf[i][0] -> value, 1);

  grad (Uf[0][0], Uf[1][0]);

  for (i = 0; i< DIM; i++) {
    Blas::axpy (ntot, -dt, *Uf[i][0] -> value, 1, *Us[i][0] -> value, 1);
    dsSmooth (Us[i][0]);
  }
}





static void curlCurl (int            np    ,  int            side,
		      int            estart,  int            skip,
		      const Element *u     ,  const Element *v   ,
		      double        *wx    ,  double        *wy  )
/* ========================================================================= *
 * Evaluate dw/dx & dw/dy (where w is the z-component of vorticity) from     *
 * element velocity fields, according to the side of the element on which    *
 * they are required.                                                        *
 *                                                                           *
 * NB: sense of traverse in wx & wy is BLAS-conformant (according to sign    *
 * of skip on relevant edge).                                                *
 * ========================================================================= */
{
  int       ntot = SQR (np);

  double  **vx = dmatrix (0, np-1, 0, np-1);
  double  **uy = dmatrix (0, np-1, 0, np-1);
  double  **w  = dmatrix (0, np-1, 0, np-1);

  double   *tmpA = dvector (0, ntot-1);
  double   *tmpB = dvector (0, ntot-1);

  double  **DV, **DT;


  quadOps (LL, np, np, NULL, NULL, NULL, NULL, NULL, &DV, &DT);

  /* -- dvdx = dvdr*drdx + dvds*dsdx */

  Blas::mxm     (*v -> value, np, *DT,         np, tmpA, np);
  Blas::mxm     (*DV,         np, *v -> value, np, tmpB, np);
  Veclib::vmul  (ntot, tmpA, 1, *v -> drdx, 1, tmpA, 1);
  Veclib::vvtvp (ntot, tmpB, 1, *v -> dsdx, 1, tmpA, 1, *vx, 1);

  /* -- dudy = dudr*drdy + duds*dsdy */

  Blas::mxm     (*u -> value, np, *DT,         np, tmpA, np);
  Blas::mxm     (*DV,         np, *u -> value, np, tmpB, np);
  Veclib::vmul  (ntot, tmpA, 1, *u -> drdy, 1, tmpA, 1);
  Veclib::vvtvp (ntot, tmpB, 1, *u -> dsdy, 1, tmpA, 1, *uy, 1);

  /* -- w = dv/dx - du/dy. */

  Veclib::vsub  (ntot, *vx, 1, *uy, 1, *w, 1);

  /* -- find dw/dx & dw/dy on appropriate edge. */

  switch (side) {
  case 0:
    Blas::gemv  ("T", np, np,  1.0, *DV, np, *w  + estart, skip, 0.0, tmpA, 1);
    Blas::gemv  ("N", np, np,  1.0, *w , np, *DV + estart, skip, 0.0, tmpB, 1);
    break;
  case 2:
    Blas::gemv  ("T", np, np, -1.0, *DV, np, *w  + estart, skip, 0.0, tmpA, 1);
    Blas::gemv  ("N", np, np, -1.0, *w , np, *DV + estart, skip, 0.0, tmpB, 1);
    break;
  case 1:
    Blas::gemv  ("T", np, np,  1.0, *w , np, *DT + estart, skip, 0.0, tmpA, 1);
    Blas::gemv  ("T", np, np,  1.0, *DV, np, *w  + estart, skip, 0.0, tmpB, 1);
    break;
  case 3:
    Blas::gemv  ("T", np, np, -1.0, *w , np, *DT + estart, skip, 0.0, tmpA, 1);
    Blas::gemv  ("T", np, np, -1.0, *DV, np, *w  + estart, skip, 0.0, tmpB, 1);
    break;
  }
   
  Veclib::vmul  (np, tmpA, 1, *u -> drdx + estart, skip, wx, 1);
  Veclib::vvtvp (np, tmpB, 1, *u -> dsdx + estart, skip, wx, 1, wx, 1);

  Veclib::vmul  (np, tmpA, 1, *u -> drdy + estart, skip, wy, 1);
  Veclib::vvtvp (np, tmpB, 1, *u -> dsdy + estart, skip, wy, 1, wy, 1);

  freeDmatrix (vx,   0, 0);
  freeDmatrix (uy,   0, 0);
  freeDmatrix (w ,   0, 0);
  freeDvector (tmpA, 0);
  freeDvector (tmpB, 0);
}
