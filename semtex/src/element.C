/*****************************************************************************
 * ELEMENT.C: Element utility routines.                                      *
 *****************************************************************************/

static char
RCSid[] = "$Id$";


#include "Fem.h"


int countElmts (const Element *E)
/* ========================================================================= *
 * Traverse & count.                                                         *
 * ========================================================================= */
{
  int n = 0;

  for (; E; E = E -> next) n++;

  return n;
}





void meshElement (Element *E     ,
		  int     *vnode ,
		  Point   *vertex)
/* ========================================================================= *
 * Generate element internal node locations.                                 *
 *                                                                           *
 * E      : pointer to current element                                       *
 * vnode  : pointer to an array of indices for locations in vertex list      *
 * vertex : vertex locations                                                 *
 *                                                                           *
 * For now the mesh has only linear edges.                                   *
 * ========================================================================= */
{
  double  *z;
  double **x  = E->xmesh,
         **y  = E->ymesh;
  int      np = E->np, nm = E->np-1, i, j;


  quadOps (LL, np, np, &z, 0, 0, 0, 0, 0, 0);
  
  /* Load element vertices. */
  
  x[ 0][ 0] = vertex[ vnode[0] ].x;
  x[ 0][nm] = vertex[ vnode[1] ].x;
  x[nm][nm] = vertex[ vnode[2] ].x;
  x[nm][ 0] = vertex[ vnode[3] ].x;

  y[ 0][ 0] = vertex[ vnode[0] ].y;
  y[ 0][nm] = vertex[ vnode[1] ].y;
  y[nm][nm] = vertex[ vnode[2] ].y;
  y[nm][ 0] = vertex[ vnode[3] ].y;

  /* Linear interpolation along element edges to get edge-internal points.  */
  
  for (i = 1; i < nm; i++) {

    x[ 0][ i] = 0.5*( (1.0 - z[i])*x[ 0][ 0] + (1.0 + z[i])*x[ 0][nm] );
    x[ i][nm] = 0.5*( (1.0 - z[i])*x[ 0][nm] + (1.0 + z[i])*x[nm][nm] );
    x[nm][ i] = 0.5*( (1.0 - z[i])*x[nm][ 0] + (1.0 + z[i])*x[nm][nm] );
    x[ i][ 0] = 0.5*( (1.0 - z[i])*x[ 0][ 0] + (1.0 + z[i])*x[nm][ 0] );

    y[ 0][ i] = 0.5*( (1.0 - z[i])*y[ 0][ 0] + (1.0 + z[i])*y[ 0][nm] );
    y[ i][nm] = 0.5*( (1.0 - z[i])*y[ 0][nm] + (1.0 + z[i])*y[nm][nm] );
    y[nm][ i] = 0.5*( (1.0 - z[i])*y[nm][ 0] + (1.0 + z[i])*y[nm][nm] );
    y[ i][ 0] = 0.5*( (1.0 - z[i])*y[ 0][ 0] + (1.0 + z[i])*y[nm][ 0] );

  }

  /* Make interior points by averaging linear & bilinear interpolations. */

  for (i = 1; i < nm; i++)
    for (j = 1; j < nm; j++) {

      x[i][j] = 0.50 * ( (1.0 - z[j])*x[i][0] + (1.0 + z[j])*x[i][nm] 
		     +   (1.0 - z[i])*x[0][j] + (1.0 + z[i])*x[nm][j] )
	     
	      - 0.25 * ( (1.0 - z[j])*(1.0 - z[i])*x[ 0][ 0]
                     +   (1.0 - z[j])*(1.0 + z[i])*x[nm][ 0]
		     +   (1.0 - z[i])*(1.0 + z[j])*x[ 0][nm]
		     +   (1.0 + z[i])*(1.0 + z[j])*x[nm][nm] );

      y[i][j] = 0.50 * ( (1.0 - z[j])*y[i][0] + (1.0 + z[j])*y[i][nm] 
		     +   (1.0 - z[i])*y[0][j] + (1.0 + z[i])*y[nm][j] )
	     
	      - 0.25 * ( (1.0 - z[j])*(1.0 - z[i])*y[ 0][ 0]
                     +   (1.0 - z[j])*(1.0 + z[i])*y[nm][ 0]
		     +   (1.0 - z[i])*(1.0 + z[j])*y[ 0][nm]
		     +   (1.0 + z[i])*(1.0 + z[j])*y[nm][nm] );
    }
}





void  mapElements (Element  *E)
/* ========================================================================= *
 * Generate geometric factors associated with mapping from 2D Cartesian to   *
 * isoparametrically-mapped space:                                           *
 *                                                                           *
 *   dxdr, drdx,   = dx/dr,  dr/dx,  "Forward Partials"                      *
 *   dxds, dsdx,   = dx/ds,  ds/dx,                                          *
 *   dydr, drdy,   = dy/dr,  dr/dy,  "Inverse Partials"                      *
 *   dyds, dsdy,   = dy/ds,  ds/dy,                                          *
 *   jac           = dx/dr * dy/ds - dx/ds * dy/dr.                          *
 *                                                                           *
 * The following relationships are used, where                               *
 *                                                                           *
 *   IN[j][k] = h_k (x_j) is a Lagrangian interpolant matrix operator, and   * 
 *   DV[j][k] = h'_k(x_j) is a Lagrangian derivative  matrix operator,       *
 *                                                                           * 
 * [ IT = transpose(IN), DT = transpose(DV) ]:                               *
 *                                                                           *
 *   [dxdr] = [IN][X][DT];    [dydr] = [IN][Y][DT],                          * 
 *   [dxds] = [DV][X][IT];    [dyds] = [DV][Y][IT].                          *
 *                                                                           *
 * For a GL quadrature rule, the inverse partials and mass matrix are        *
 * returned for spatial locations at the mesh nodes, while the forward       *
 * partials and other geometric factors are for spatial locations at the     *
 * quadrature points.  For LL rule, everything is at the nodes.              *
 * ========================================================================= */
{
  char      routine[] = "mapElements";
  char      buf[STR_MAX];
  double  **x,    **y;
  double  **dxdr, **drdx, **dxds, **dsdx;
  double  **dydr, **drdy, **dyds, **dsdy;
  double  **G1,   **G2,   **G3,   **G4;
  double  **mass;
  double  **jac,  **DV, **DT, **IN, **IT;
  double  **tM,    *tV;		            /* Temporaries.                  */
  double    *w,   **WW;		            /* Weights & w.w' outer product. */
  int       np, nq, ntot;


  if ( option ("RULE") == LL ) {

    for (; E; E = E -> next) {
      np   = E -> np;
      ntot = SQR (np);
      x    = E -> xmesh;
      y    = E -> ymesh;

      E -> dxdr = dxdr = dmatrix (0, np-1, 0, np-1);
      E -> drdx = drdx = dmatrix (0, np-1, 0, np-1);
      E -> dxds = dxds = dmatrix (0, np-1, 0, np-1);
      E -> dsdx = dsdx = dmatrix (0, np-1, 0, np-1);
      E -> dydr = dydr = dmatrix (0, np-1, 0, np-1);
      E -> drdy = drdy = dmatrix (0, np-1, 0, np-1);
      E -> dyds = dyds = dmatrix (0, np-1, 0, np-1);
      E -> dsdy = dsdy = dmatrix (0, np-1, 0, np-1);
      E -> G1   = G1   = dmatrix (0, np-1, 0, np-1);
      E -> G2   = G2   = dmatrix (0, np-1, 0, np-1);
      E -> G3   = G3   = dmatrix (0, np-1, 0, np-1);
      E -> G4   = G4   = dmatrix (0, np-1, 0, np-1);
      E -> mass = G4;

      jac = dmatrix (0, np-1, 0, np-1);
      WW  = dmatrix (0, np-1, 0, np-1);
      tV  = dvector (0, ntot-1);

      if (!(dxdr && drdx && dxds && dsdx && dydr && drdy && dyds && dsdy
	    && G1 && G2 && G3 && G4 && jac && WW && tV))
	message (routine, "allocation failure", ERROR);

      quadOps (LL, np, np, 0, 0, &w, 0, 0, &DV, &DT);
      Veclib::zero   (ntot, *WW, 1);
      Blas::  ger    (np, np, 1.0, w, 1, w, 1, *WW, np);
      
      Blas::mxm ( *x, np, *DT, np, *dxdr, np);
      Blas::mxm (*DV, np,  *x, np, *dxds, np);
      Blas::mxm ( *y, np, *DT, np, *dydr, np);
      Blas::mxm (*DV, np,  *y, np, *dyds, np);

      Veclib::vmul  (ntot,        *dxdr, 1, *dyds, 1,  tV,  1);
      Veclib::vvvtm (ntot, tV, 1, *dxds, 1, *dydr, 1, *jac, 1);

      if ((*jac)[Veclib::imin (ntot, *jac, 1)] <= EPSSP) {
	sprintf (buf, "jacobian of element %1d nonpositive", E -> id + 1);
	message (routine, buf, ERROR);
      }
	
      Veclib::vmul  (ntot, *dyds, 1, *dyds, 1,  tV, 1);
      Veclib::vvtvp (ntot, *dxds, 1, *dxds, 1,  tV, 1, *G1, 1);
      Veclib::vdiv  (ntot, *G1,   1, *jac,  1,  tV, 1);
      Veclib::vmul  (ntot,  tV,   1, *WW,   1, *G1, 1);

      Veclib::vmul  (ntot, *dydr, 1, *dydr, 1,  tV, 1);
      Veclib::vvtvp (ntot, *dxdr, 1, *dxdr, 1,  tV, 1, *G2, 1);
      Veclib::vdiv  (ntot, *G2,   1, *jac,  1,  tV, 1);
      Veclib::vmul  (ntot,  tV,   1, *WW,   1, *G2, 1);

      Veclib::vmul  (ntot, *dydr, 1, *dyds, 1,  tV, 1);
      Veclib::neg   (ntot, tV,    1);
      Veclib::vvvtm (ntot, tV,    1, *dxdr, 1, *dxds, 1, *G3, 1);
      Veclib::vdiv  (ntot, *G3,   1, *jac,  1,  tV, 1);
      Veclib::vmul  (ntot,  tV,   1, *WW,   1, *G3, 1);

      Veclib::vmul  (ntot, *jac,  1, *WW,   1, *G4, 1);

      Veclib::copy (ntot, *dyds, 1, *drdx, 1);
      Veclib::vneg (ntot, *dxds, 1, *drdy, 1);
      Veclib::vneg (ntot, *dydr, 1, *dsdx, 1);
      Veclib::copy (ntot, *dxdr, 1, *dsdy, 1);

      Veclib::vdiv (ntot, *drdx, 1, *jac, 1, *drdx, 1);
      Veclib::vdiv (ntot, *drdy, 1, *jac, 1, *drdy, 1);
      Veclib::vdiv (ntot, *dsdx, 1, *jac, 1, *dsdx, 1);
      Veclib::vdiv (ntot, *dsdy, 1, *jac, 1, *dsdy, 1);

      freeDmatrix (jac, 0, 0);
      freeDmatrix (WW,  0, 0);
      freeDvector (tV,  0);
    }
      
  } else {	                    /* RULE == GL */

    for (; E; E = E -> next) {
      np   = E -> np;
      nq   = E -> nq;
      x    = E -> xmesh;
      y    = E -> ymesh;
    
      /* Quadrature point computations. */

      ntot = SQR(nq);

      E -> dxdr = dxdr = dmatrix (0, nq-1, 0, nq-1);
      E -> dxds = dxds = dmatrix (0, nq-1, 0, nq-1);
      E -> dydr = dydr = dmatrix (0, nq-1, 0, nq-1);
      E -> dyds = dyds = dmatrix (0, nq-1, 0, nq-1);
      E -> G1   = G1   = dmatrix (0, nq-1, 0, nq-1);
      E -> G2   = G2   = dmatrix (0, nq-1, 0, nq-1);
      E -> G3   = G3   = dmatrix (0, nq-1, 0, nq-1);
      E -> G4   = G4   = dmatrix (0, nq-1, 0, nq-1);

      jac = dmatrix (0, nq-1, 0, nq-1);
      WW  = dmatrix (0, nq-1, 0, nq-1);
      tM  = dmatrix (0, np-1, 0, nq-1);
      tV  = dvector (0, ntot-1);
      
      if (!(dxdr && dxds && dydr && dyds && G1 && G2 && G3 && G4
	    && jac && WW && tM && tV))
	message (routine, "allocation failure", ERROR);

      quadOps (GL, np, nq, 0, 0, &w, &IN, &IT, &DV, &DT);
      Veclib::zero   (ntot, *WW, 1);
      Blas::  ger    (nq, nq, 1.0, w, 1, w, 1, *WW, nq);

      Blas::mxm ( *x, np, *DT, np, *tM,   nq);
      Blas::mxm (*IN, nq, *tM, np, *dxdr, nq);
      Blas::mxm ( *y, np, *DT, np, *tM,   nq);
      Blas::mxm (*IN, nq, *tM, np, *dydr, nq);
      Blas::mxm ( *x, np, *IT, np, *tM,   nq);
      Blas::mxm (*DV, nq, *tM, np, *dxds, nq);
      Blas::mxm ( *y, np, *IT, np, *tM,   nq);
      Blas::mxm (*DV, nq, *tM, np, *dyds, nq);

      Veclib::vmul  (ntot,        *dxdr, 1, *dyds, 1,  tV,  1);
      Veclib::vvvtm (ntot, tV, 1, *dxds, 1, *dydr, 1, *jac, 1);

      Veclib::vmul  (ntot, *dyds, 1, *dyds, 1,  tV, 1);
      Veclib::vvtvp (ntot, *dxds, 1, *dxds, 1,  tV, 1, *G1, 1);
      Veclib::vdiv  (ntot, *G1,   1, *jac,  1,  tV, 1);
      Veclib::vmul  (ntot,  tV,   1, *WW,   1, *G1, 1);

      Veclib::vmul  (ntot, *dydr, 1, *dydr, 1,  tV, 1);
      Veclib::vvtvp (ntot, *dxdr, 1, *dxdr, 1,  tV, 1, *G2, 1);
      Veclib::vdiv  (ntot, *G2,   1, *jac,  1,  tV, 1);
      Veclib::vmul  (ntot,  tV,   1, *WW,   1, *G2, 1);

      Veclib::vmul  (ntot, *dydr, 1, *dyds, 1,  tV, 1);
      Veclib::neg   (ntot, tV,    1);
      Veclib::vvvtm (ntot, tV,    1, *dxdr, 1, *dxds, 1, *G3, 1);
      Veclib::vdiv  (ntot, *G3,   1, *jac,  1,  tV, 1);
      Veclib::vmul  (ntot,  tV,   1, *WW,   1, *G3, 1);

      Veclib::vmul  (ntot, *jac,  1, *WW,   1, *G4, 1);

      freeDmatrix (jac, 0, 0);
      freeDmatrix (WW,  0, 0);
      freeDmatrix (tM,  0, 0);
      freeDvector (tV,  0);

      /* Node point computations. */
    
      ntot = SQR (np);

      E -> drdx = drdx = dmatrix (0, np-1, 0, np-1);
      E -> drdy = drdy = dmatrix (0, np-1, 0, np-1);
      E -> dsdx = dsdx = dmatrix (0, np-1, 0, np-1);
      E -> dsdy = dsdy = dmatrix (0, np-1, 0, np-1);
      E -> mass = mass = dmatrix (0, np-1, 0, np-1);

      dxdr = dmatrix (0, np-1, 0, np-1);
      dxds = dmatrix (0, np-1, 0, np-1);
      dydr = dmatrix (0, np-1, 0, np-1);
      dyds = dmatrix (0, np-1, 0, np-1);
      jac  = dmatrix (0, np-1, 0, np-1);
      WW   = dmatrix (0, np-1, 0, np-1);
      tV   = dvector (0, ntot-1);
      
      if (!(drdx && drdy && dsdx && dsdy && mass && dxdr && dxds
	    && dydr && dyds && jac && WW && tV))
	message (routine, "allocation failure", ERROR);

      quadOps (LL, np, np, 0, 0, &w, 0, 0, &DV, &DT);
      Veclib::zero   (ntot, *WW, 1);
      Blas::  ger    (np, np, 1.0, w, 1, w, 1, *WW, np);

      Blas::mxm ( *x, np, *DT, np, *dxdr, np);
      Blas::mxm (*DV, np,  *x, np, *dxds, np);
      Blas::mxm ( *y, np, *DT, np, *dydr, np);
      Blas::mxm (*DV, np,  *y, np, *dyds, np);

      Veclib::vmul  (ntot,        *dxdr, 1, *dyds, 1,  tV,  1);
      Veclib::vvvtm (ntot, tV, 1, *dxds, 1, *dydr, 1, *jac, 1);

      Veclib::vmul (ntot, *jac, 1, *WW, 1, *mass, 1);

      Veclib::copy (ntot, *dyds, 1, *drdx, 1);
      Veclib::vneg (ntot, *dxds, 1, *drdy, 1);
      Veclib::vneg (ntot, *dydr, 1, *dsdx, 1);
      Veclib::copy (ntot, *dxdr, 1, *dsdy, 1);

      Veclib::vdiv (ntot, *drdx, 1, *jac, 1, *drdx, 1);
      Veclib::vdiv (ntot, *drdy, 1, *jac, 1, *drdy, 1);
      Veclib::vdiv (ntot, *dsdx, 1, *jac, 1, *dsdx, 1);
      Veclib::vdiv (ntot, *dsdy, 1, *jac, 1, *dsdy, 1);

      freeDmatrix (dxdr, 0, 0);
      freeDmatrix (dxds, 0, 0);
      freeDmatrix (dydr, 0, 0);
      freeDmatrix (dyds, 0, 0);
      freeDmatrix (jac,  0, 0);
      freeDmatrix (WW,   0, 0);
      freeDvector (tV,   0);
    }
  }
}




