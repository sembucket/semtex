//////////////////////////////////////////////////////////////////////////////
// element.C: Element utility routines.
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include "Fem.h"


Element::Element (const int& ident, const int& nknot, const int& nside)
// ---------------------------------------------------------------------------
// Create a new element.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Element::Element";

  if (nside != 4)
    message (routine, "need 4 sides for quadrilateral element",  ERROR);
  else if (nknot < 2)
    message (routine, "need at least 2 knots for element edges", ERROR);
  else
    setState (ident, nknot, nside);
}


Element::Element (const Element& parent, const int& np)
// ---------------------------------------------------------------------------
// Use parent as a model for a new element, possibly of different order.
// ---------------------------------------------------------------------------
{
  if (!np) {
    memcpy (this, &parent, sizeof (Element));

  } else if (np == parent.np) {
    memcpy (this, &parent, sizeof (Element));
    bmap  = solve = 0;
    drdx = dsdx = drdy = dsdy = G1 = G2 = G3 = G4 = mass = 0;
    hbi  = hii  = 0;

  } else {
    setState (parent.id, np, parent.ns);
    xmesh = ymesh = 0;
    bmap  = solve = 0;
    drdx = dsdx = drdy = dsdy = G1 = G2 = G3 = G4 = mass = 0;
    hbi  = hii  = 0;

  }
}


void  Element::project (const Element& parent, const real* src, real* target)
// ---------------------------------------------------------------------------
// Isoparametrically project mesh from parent Element onto this Element.
// Optionally project Field data (src) from parent onto projected element.
// ---------------------------------------------------------------------------
{
  char       routine[] = "Element::project";
  const int  nold      = parent.np;

  if (!nold || !np)
    message (routine, "size of either this or parent Element not set", ERROR);
  
  if (nold == np) {
    Veclib::copy (nTot(), *parent.xmesh, 1, *xmesh, 1);
    Veclib::copy (nTot(), *parent.ymesh, 1, *ymesh, 1);
    if (src) Veclib::copy (nTot(), src, 1, target, 1);

  } else {
    real  **IN, **IT, *tmp = rvector (np*nold);

    meshOps (GLL, GLL, nold, np, 0, &IN, &IT, 0, 0);

    Blas::mxm (*IN, np, *parent.xmesh, nold, tmp,    nold);
    Blas::mxm (tmp, np, *IT,           nold, *xmesh, np  );
      
    Blas::mxm (*IN, np, *parent.ymesh, nold, tmp,    nold);
    Blas::mxm (tmp, np, *IT,           nold, *ymesh, np  );

    if (src) {
      Blas::mxm (*IN, np, src, nold, tmp,    nold);
      Blas::mxm (tmp, np, *IT, nold, target, np  );
    }
    
    freeVector (tmp);
  }
}


void  Element::setState (int ident, int nknot, int nside)
// ---------------------------------------------------------------------------
// Set element state variables.
// ---------------------------------------------------------------------------
{
  char routine[] = "Element::setState";

  id   = ident;
  np   = nknot;
  ns   = nside;
  rule = option ("RULE");

  switch (rule) {
  case LL:
    nq = np;
    break;
  case GL:
    nq = quadComplete (2, np);
    break;
  default:
    message (routine, "quadrature rule set incorrectly", ERROR);
    break;
  }

  setEdgeMaps (*this);
}


void Element::setEdgeMaps (Element& E)
// ---------------------------------------------------------------------------
// An emap is an edge-based list of indices, going CCW around the element
// edges and then traversing the interior.  This allows us to access element
// storage by traversing edges, tying in with edge-based global numbering.
//
// Pmap (the inversion of emap) is built afterwards.
//
// To convert from row-major storage to boundary-major, then back:
//   gathr (ntot, value, emap, tmp  );
//   gathr (ntot, tmp,   pmap, value)  OR  scatr (ntot, tmp, emap, value);
// ---------------------------------------------------------------------------
{
  char  routine[] = "Element::setEdgeMaps";

  if (!E.np) message (routine, "element order not set", ERROR);

  static List<mapping*> maplist;

  for (ListIterator<mapping*> m(maplist); m.more(); m .next())
    if (m.current () -> ne == E.np) {
      E.emap = m.current () -> e;
      E.pmap = m.current () -> p;
      return;
    }

  // -- Not found: build a new pair of emap & pmap.

  register int  i, j, k, n;
  
  E.emap = ivector (sqr (E.np));

  // -- Traverse exterior CCW.

  k = 0;
  n = 0;
  E.emap[0] = 0;
  for (i = 1; i < E.np; i++) E.emap[++k] = n += 1;
  for (i = 1; i < E.np; i++) E.emap[++k] = n += E.np;
  for (i = 1; i < E.np; i++) E.emap[++k] = n -= 1;
  for (i = 1; i < E.np; i++) E.emap[++k] = n -= E.np;
  
  // -- Traverse interior in row-major order.

  int nm = E.np - 1;
  for (i = 1; i < nm; i++)
    for (j = 1; j < nm; j++)
      E.emap[k++] = i * E.np + j;
  
  // -- Build pmap.
  
  E.pmap = ivector (sqr (E.np));
  for (i = 0; i < sqr (E.np); i++) E.pmap[E.emap[i]] = i;
    
  // -- Insert in the list of maps.

  mapping *M;
  maplist.add (M = new mapping);
  M -> e = E.emap;
  M -> p = E.pmap;
}


inline
void  Element::terminal (int side, int& estart, int& eskip, int& bstart) const
// ---------------------------------------------------------------------------
// Evaluate the element-edge terminal values of estart, skip, bstart.
// NB: BLAS-conformant terminal start values are delivered for negative skips.
// ---------------------------------------------------------------------------
{
  switch (side) {
  case 1:
    estart = 0;
    eskip  = 1;
    bstart = 0;
    break;
  case 2:
    estart = np - 1;
    eskip  = np;
    bstart = np - 1;
    break;
  case 3:
    estart = np * (np - 1);
    eskip  = -1;
    bstart = 2  * (np - 1);
    break;
  case 4:
    estart = 0;
    eskip  = -np;
    bstart = 3  * (np - 1);
    break;
  }
}


void  Element::install (int jump, real* mesh, int* map, int* mask)
// ---------------------------------------------------------------------------
// Make Element xmesh, ymesh, bmap & solve point to new memory, set offset
// in Field data area.  Mesh only installed if pointer is non-zero.
// ---------------------------------------------------------------------------
{
  if (jump) offset = jump;

  if (mesh) {
    xmesh = new real* [np];
    ymesh = new real* [np];
    
    real *xm = mesh;
    real *ym = mesh + sqr(np);

    for (register int j = 0; j < np; j++) {
      xmesh[j] = xm + j * np;
      ymesh[j] = ym + j * np;
    }
  }
  
  if (mask) solve = mask;
  if (map)  bmap  = map;
}


void  Element::gidInsert (const int* b)
// ---------------------------------------------------------------------------
// Insert an array with Element-boundary global node numbers into bmap.
// ---------------------------------------------------------------------------
{
  Veclib::copy (nExt(), b, 1, bmap, 1);
}


void  Element::mesh (const Mesh& M)
// ---------------------------------------------------------------------------
// Generate Element mesh internal node locations.
// ---------------------------------------------------------------------------
{
  real*         z;
  real**        x  = xmesh;
  real**        y  = ymesh;
  register int  i, j, nm = np - 1;
  Point*        P = new Point [np];

  // -- Build Element-boundary knots using Mesh functions.

  quadOps (LL, np, np, &z, 0, 0, 0, 0, 0, 0);

  for (j = 1; j <= ns; j++) {
    M.meshSide (np, id, j, z, P);
    for (i = 0; i < nm; i++) {
      switch (j) {
      case 1:
	x[0][i]       = P[i].x;
	y[0][i]       = P[i].y;
	break;
      case 2:
	x[i][nm]      = P[i].x;
	y[i][nm]      = P[i].y;
	break;
      case 3:
	x[nm][nm - i] = P[i].x;
	y[nm][nm - i] = P[i].y;
	break;
      case 4:
	x[nm - i][0]  = P[i].x;
	y[nm - i][0]  = P[i].y;
	break;
      }
    }
  }

  // --  Make interior knots by averaging linear & bilinear interpolations.

  for (i = 1; i < nm; i++)
    for (j = 1; j < nm; j++) {

      x[i][j] = 0.50 * ( (1.0 - z[j]) * x[i][0] + (1.0 + z[j]) * x[i][nm] 
		     +   (1.0 - z[i]) * x[0][j] + (1.0 + z[i]) * x[nm][j] )
	     
	      - 0.25 * ( (1.0 - z[j]) * (1.0 - z[i]) * x[ 0][ 0]
                     +   (1.0 - z[j]) * (1.0 + z[i]) * x[nm][ 0]
		     +   (1.0 - z[i]) * (1.0 + z[j]) * x[ 0][nm]
		     +   (1.0 + z[i]) * (1.0 + z[j]) * x[nm][nm] );

      y[i][j] = 0.50 * ( (1.0 - z[j]) * y[i][0] + (1.0 + z[j]) * y[i][nm] 
		     +   (1.0 - z[i]) * y[0][j] + (1.0 + z[i]) * y[nm][j] )
	     
	      - 0.25 * ( (1.0 - z[j]) * (1.0 - z[i]) * y[ 0][ 0]
                     +   (1.0 - z[j]) * (1.0 + z[i]) * y[nm][ 0]
		     +   (1.0 - z[i]) * (1.0 + z[j]) * y[ 0][nm]
		     +   (1.0 + z[i]) * (1.0 + z[j]) * y[nm][nm] );
    }

  delete [] P;
}


void  Element::map ()
// ---------------------------------------------------------------------------
// Generate geometric factors associated with mapping from 2D Cartesian to
// isoparametrically-mapped space:
//
//   dxdr, drdx,   = dx/dr,  dr/dx,  "Forward Partials"
//   dxds, dsdx,   = dx/ds,  ds/dx,
//   dydr, drdy,   = dy/dr,  dr/dy,  "Inverse Partials"
//   dyds, dsdy,   = dy/ds,  ds/dy,
//   jac           = dx/dr * dy/ds - dx/ds * dy/dr.
//
// The following relationships are used, where
//
//   IN[j][k] = h_k (x_j) is a Lagrangian interpolant matrix operator, and
//   DV[j][k] = h'_k(x_j) is a Lagrangian derivative  matrix operator,
//
// [ IT = transpose(IN), DT = transpose(DV) ]:
//
//   [dxdr] = [IN][X][DT];    [dydr] = [IN][Y][DT],
//   [dxds] = [DV][X][IT];    [dyds] = [DV][Y][IT].
//
// For a GL quadrature rule, the inverse partials and mass matrix are
// returned for spatial locations at the mesh nodes, while the forward
// partials and other geometric factors are for spatial locations at the
// quadrature points.  For LL rule, everything is at the nodes.
// ---------------------------------------------------------------------------
{
  char    routine[] = "Element::map";
  char    buf[StrMax];
  int     ntot;
  real   *x, *y;
  real  **DV,   **DT,   **IN,   **IT;
  real   *jac,   *dxdr,  *dxds,  *dydr,  *dyds,  *tM,   *tV,   *w,  *WW;

  const real EPS = (sizeof(real) == sizeof (double)) ? EPSDP : EPSSP;

  if (rule == LL) {
    ntot = nTot();
    x    = *xmesh;
    y    = *ymesh;

    drdx = rvector (ntot);
    dsdx = rvector (ntot);
    drdy = rvector (ntot);
    dsdy = rvector (ntot);
    G1   = rvector (ntot);
    G2   = rvector (ntot);
    G3   = rvector (ntot);
    G4   = rvector (ntot);
    mass = G4;
    
    dxdr = rvector (ntot);
    dxds = rvector (ntot);
    dydr = rvector (ntot);
    dyds = rvector (ntot);

    jac  = rvector (ntot);
    WW   = rvector (ntot);
    tV   = rvector (ntot);
    
    quadOps (LL, np, np, 0, 0, &w, 0, 0, &DV, &DT);
    Veclib::zero (ntot, WW, 1);
    Blas::ger    (np, np, 1.0, w, 1, w, 1, WW, np);
    
    Blas::mxm (  x, np, *DT, np, dxdr, np);
    Blas::mxm (*DV, np,   x, np, dxds, np);
    Blas::mxm (  y, np, *DT, np, dydr, np);
    Blas::mxm (*DV, np,   y, np, dyds, np);
    
    Veclib::vmul  (ntot,        dxdr, 1, dyds, 1, tV,  1);
    Veclib::vvvtm (ntot, tV, 1, dxds, 1, dydr, 1, jac, 1);
    
    if (jac[Veclib::imin (ntot, jac, 1)] <= EPSSP) {
      sprintf (buf, "jacobian of element %1d nonpositive", id);
      message (routine, buf, ERROR);
    }
    
    Veclib::vmul  (ntot, dyds, 1, dyds, 1, tV, 1);
    Veclib::vvtvp (ntot, dxds, 1, dxds, 1, tV, 1, G1, 1);
    Veclib::vdiv  (ntot, G1,   1, jac,  1, tV, 1);
    Veclib::vmul  (ntot, tV,   1, WW,   1, G1, 1);
    
    Veclib::vmul  (ntot, dydr, 1, dydr, 1, tV, 1);
    Veclib::vvtvp (ntot, dxdr, 1, dxdr, 1, tV, 1, G2, 1);
    Veclib::vdiv  (ntot, G2,   1, jac,  1, tV, 1);
    Veclib::vmul  (ntot, tV,   1, WW,   1, G2, 1);
    
    Veclib::vmul  (ntot, dydr, 1, dyds, 1, tV,   1);
    Veclib::neg   (ntot, tV,   1);
    Veclib::vvvtm (ntot, tV,   1, dxdr, 1, dxds, 1, G3, 1);
    Veclib::vdiv  (ntot, G3,   1, jac,  1, tV,   1);
    Veclib::vmul  (ntot, tV,   1, WW,   1, G3,   1);

    if (Blas::nrm2 (ntot, G3, 1) < EPS) {
      freeVector (G3);
      G3 = 0;
    }
    
    Veclib::vmul  (ntot, jac, 1, WW, 1, G4, 1);
    
    Veclib::copy (ntot, dyds, 1, drdx, 1);
    Veclib::vneg (ntot, dxds, 1, drdy, 1);
    Veclib::vneg (ntot, dydr, 1, dsdx, 1);
    Veclib::copy (ntot, dxdr, 1, dsdy, 1);
    
    Veclib::vdiv (ntot, drdx, 1, jac, 1, drdx, 1);
    Veclib::vdiv (ntot, drdy, 1, jac, 1, drdy, 1);
    Veclib::vdiv (ntot, dsdx, 1, jac, 1, dsdx, 1);
    Veclib::vdiv (ntot, dsdy, 1, jac, 1, dsdy, 1);

    freeVector (dxdr);
    freeVector (dxds);
    freeVector (dydr);
    freeVector (dyds);
    freeVector (jac);
    freeVector (WW);
    freeVector (tV);
      
  } else {  // -- rule == GL
    x = *xmesh;
    y = *ymesh;
    
    // -- Quadrature point computations.

    ntot = sqr (nq);

    G1   = rvector (ntot);
    G2   = rvector (ntot);
    G3   = rvector (ntot);
    G4   = rvector (ntot);
     
    dxdr = rvector (ntot);
    dxds = rvector (ntot);
    dydr = rvector (ntot);
    dyds = rvector (ntot);
   
    jac  = rvector (ntot);
    WW   = rvector (ntot);
    tM   = rvector (np * nq);
    tV   = rvector (ntot);

    quadOps (GL, np, nq, 0, 0, &w, &IN, &IT, &DV, &DT);
    Veclib::zero (ntot, WW, 1);
    Blas::ger    (nq, nq, 1.0, w, 1, w, 1, WW, nq);
    
    Blas::mxm (  x, np, *DT, np, tM,   nq);
    Blas::mxm (*IN, nq,  tM, np, dxdr, nq);
    Blas::mxm (  y, np, *DT, np, tM,   nq);
    Blas::mxm (*IN, nq,  tM, np, dydr, nq);
    Blas::mxm (  x, np, *IT, np, tM,   nq);
    Blas::mxm (*DV, nq,  tM, np, dxds, nq);
    Blas::mxm (  y, np, *IT, np, tM,   nq);
    Blas::mxm (*DV, nq,  tM, np, dyds, nq);
    
    Veclib::vmul  (ntot,        dxdr, 1, dyds, 1, tV,  1);
    Veclib::vvvtm (ntot, tV, 1, dxds, 1, dydr, 1, jac, 1);
    
    Veclib::vmul  (ntot, dyds, 1, dyds, 1, tV, 1);
    Veclib::vvtvp (ntot, dxds, 1, dxds, 1, tV, 1, G1, 1);
    Veclib::vdiv  (ntot, G1,   1, jac,  1, tV, 1);
    Veclib::vmul  (ntot, tV,   1, WW,   1, G1, 1);
    
    Veclib::vmul  (ntot, dydr, 1, dydr, 1, tV, 1);
    Veclib::vvtvp (ntot, dxdr, 1, dxdr, 1, tV, 1, G2, 1);
    Veclib::vdiv  (ntot, G2,   1, jac,  1, tV, 1);
    Veclib::vmul  (ntot, tV,   1, WW,   1, G2, 1);
    
    Veclib::vmul  (ntot, dydr, 1, dyds, 1, tV,   1);
    Veclib::neg   (ntot, tV,   1);
    Veclib::vvvtm (ntot, tV,   1, dxdr, 1, dxds, 1, G3, 1);
    Veclib::vdiv  (ntot, G3,   1, jac,  1, tV,   1);
    Veclib::vmul  (ntot, tV,   1, WW,   1, G3,   1);

    if (Blas::nrm2 (ntot, G3, 1) < EPS) {
      freeVector (G3);
      G3 = 0;
    }
    
    Veclib::vmul  (ntot, jac, 1, WW, 1, G4, 1);

    freeVector (dxdr);
    freeVector (dxds);
    freeVector (dydr);
    freeVector (dyds);
    freeVector (jac);
    freeVector (WW);
    freeVector (tM);
    freeVector (tV);
    
    // -- Node point computations.
    
    ntot = nTot();

    drdx = rvector (ntot);
    drdy = rvector (ntot);
    dsdx = rvector (ntot);
    dsdy = rvector (ntot);
    mass = rvector (ntot);
    
    dxdr = rvector (ntot);
    dxds = rvector (ntot);
    dydr = rvector (ntot);
    dyds = rvector (ntot);
    jac  = rvector (ntot);
    WW   = rvector (ntot);
    tV   = rvector (ntot);
    
    quadOps      (LL, np, np, 0, 0, &w, 0, 0, &DV, &DT);
    Veclib::zero (ntot, WW, 1);
    Blas::ger    (np, np, 1.0, w, 1, w, 1, WW, np);
    
    Blas::mxm (  x, np, *DT, np, dxdr, np);
    Blas::mxm (*DV, np,   x, np, dxds, np);
    Blas::mxm (  y, np, *DT, np, dydr, np);
    Blas::mxm (*DV, np,   y, np, dyds, np);
    
    Veclib::vmul  (ntot,        dxdr, 1, dyds, 1, tV,  1);
    Veclib::vvvtm (ntot, tV, 1, dxds, 1, dydr, 1, jac, 1);
    
    Veclib::vmul (ntot, jac, 1, WW, 1, mass, 1);
    
    Veclib::copy (ntot, dyds, 1, drdx, 1);
    Veclib::vneg (ntot, dxds, 1, drdy, 1);
    Veclib::vneg (ntot, dydr, 1, dsdx, 1);
    Veclib::copy (ntot, dxdr, 1, dsdy, 1);
    
    Veclib::vdiv (ntot, drdx, 1, jac, 1, drdx, 1);
    Veclib::vdiv (ntot, drdy, 1, jac, 1, drdy, 1);
    Veclib::vdiv (ntot, dsdx, 1, jac, 1, dsdx, 1);
    Veclib::vdiv (ntot, dsdy, 1, jac, 1, dsdy, 1);
    
    freeVector (dxdr);
    freeVector (dxds);
    freeVector (dydr);
    freeVector (dyds);
    freeVector (jac);
    freeVector (WW);
    freeVector (tV);
  }
}


void  Element::printBndry (const real* src) const
// ---------------------------------------------------------------------------
// (Debugging) utility to print edge connectivity & value information.
// ---------------------------------------------------------------------------
{
  int   ne = nExt ();
  char  s[StrMax];

  message ("", "#   id  emap  bmap solve value", REMARK);

  for (int i = 0; i < ne; i++) {
    sprintf
      (s, "%6d%6d%6d%6d%6d", id, emap[i], bmap[i], solve[i], src[emap[i]]);
    message ("", s, REMARK);
    }
}


void  Element::bndryInsert (const real* src, real* target) const
// ---------------------------------------------------------------------------
// Project from globally-numbered storage to boundary nodes of element
// in Field storage, i.e. target[emap[i]] = src[bmap[i]].
// To invert: target[bmap[i]] = src[emap[i]].
// ---------------------------------------------------------------------------
{
  Veclib::gathr_scatr (nExt (), src, bmap, emap, target);
}


void  Element::bndryDsSum (const real* src, real* target) const
// ---------------------------------------------------------------------------
// Direct-stiffness-sum from element boundary to globally-numbered storage,
// i.e. target[bmap[i]] += mass[emap[i]] * src[emap[i]].
// This is using in smoothing Fields along element boundaries.
// ---------------------------------------------------------------------------
{
  register int    i, b, e;
  register int    nxt = nExt();
  register real*  wt = mass;

  for (i = 0; i < nxt; i++) {
    b = bmap[i];
    e = emap[i];
    target[b] += wt[e] * src[e];
  }
}


void  Element::dsForcingSC (real* F, real* target) const
// ---------------------------------------------------------------------------
// Create statically-condensed boundary Helmholtz forcing for this element
// and insert it into target by direct stiffness summation.
//  
// NB: forcing, F, is modified.
//
// As a first step, F is negated and weighted by elemental mass matrix so
// that it contains the weak form of the forcing:
//   F = - mass * value.
//
// Elemental storage is then sorted in F so that it is ordered with boundary
// nodes foremost, i.e. it contains the partition { F | F }.
//                                                   b   i
//
// Statically-condensed boundary forcing is created in the first partition:
//                       -1                         -1
//   F   <--   F  -  h  h   F           (matrix h  h   supplied by *this)
//    b         b     bi ii  i                   bi ii
//
// and summed into the target vector.  In the summation, there is no need
// to check if the global node is to be solved for or is fixed, since
// the fixed (ESSENTIAL-BC) partition of target is overwritten later.
// ---------------------------------------------------------------------------
{
  int  ntot = nTot ();
  int  nint = nInt ();
  int  next = nExt ();

  Veclib::neg  (ntot, F, 1);
  Veclib::vmul (ntot, F, 1, mass, 1, F, 1);

  real*  tmp = rvector (ntot);
  Veclib::gathr (ntot, F, emap, tmp);
  Veclib::copy  (ntot, tmp, 1, F, 1);
  freeVector (tmp);

  if (nint) Blas::gemv ("T", nint,next, -1.0, hbi,nint, F + next,1, 1.0, F,1);

  for (register int i = 0; i < next; i++) target[bmap[i]] += F[i];
}


void  Element::dsForcing (real* F) const
// ---------------------------------------------------------------------------
// Input F is negated and weighted by elemental mass matrix so
// that it contains the weak form of the forcing:
//   F = - mass * value.
// ---------------------------------------------------------------------------
{
  int  ntot = nTot ();

  Veclib::neg  (ntot, F, 1);
  Veclib::vmul (ntot, F, 1, mass, 1, F, 1);
}


void  Element::resolveSC (const real* RHS, real* F, real* target) const
// ---------------------------------------------------------------------------
// Complete static condensation solution for internal values of Element.
//
// On entry, global-node solution values are in RHS and F contains the
// weak form of internal forcing in its top end (as installed by dsForcingSC).
//
// If u is current Element, compute internal solution according to:
//            -1      -1
//   u  <--  h  F  - h  h   u
//    i       ii i    ii ib  b
// ----------------------------------------------------------------------------
{
  int    next = nExt ();
  real*  tmp  = rvector (next);

  // -- Load globally-numbered RHS into element-boundary storage.

  Veclib::gathr (next, RHS, bmap, tmp   );
  Veclib::scatr (next, tmp, emap, target);

  // -- Complete static-condensation solution for element-internal storage.

  int nint = nInt ();
  if (nint) {
    int    info = 0;
    real*  Fi   = F + next;
    Lapack::pptrs ("U", nint, 1, hii, Fi, nint, info);
    Blas::gemv    ("N", nint, next, -1.0, hbi, nint, tmp, 1, 1.0, Fi, 1);
    Veclib::scatr (nint, Fi, emap + next, target);
  }
  
  freeVector (tmp);
}


int Element::bandwidthSC () const
// ---------------------------------------------------------------------------
// Find the global equation bandwidth of this element, excluding diagonal.
// ---------------------------------------------------------------------------
{
  register int i;
  register int Min  = INT_MAX;
  register int Max  = INT_MIN;
  register int next = nExt();
  
  for (i = 0; i < next; i++)
    if (solve[i]) {
      Min = min (bmap[i], Min);
      Max = max (bmap[i], Max);
    }

  return Max - Min;
}


void  Element::HelmholtzSC (const real&  lambda2,
			    real*        hbb    ,
			    real*        rmat   ,
			    real*        rwrk   ) 
// ---------------------------------------------------------------------------
// Compute the discrete elemental Helmholtz matrix and return the
// statically condensed form in hbb, retain the interior-exterior coupling
// matrix in hbi, and the interior resolution matrix factor in hii.
//
// Uncondensed System   -->   Statically condensed form returned in hbb.
//
//                                                           -1
//  +---------+------+       +---------+    +------+  +------+  +---------+
//  |         |      |       |         |    |      |  |      |  |         |
//  |         |      |       |         |    |      |  | hii  |  |   hib   |
//  |   hbb   | hbi  |  -->  |   hbb   | -  | hbi  |  |      |  |         |
//  |         |      |       |         |    |      |  +------+  +---------+
//  |         |      |       |         |    |      |
//  +---------+------+       +---------+    +------+
//  |         |      |
//  |   hib   | hii  |
//  |         |      |
//  +---------+------+
//
//
// Element matrices are built row-by-row, sorted to place entries for
// internal nodes first, posted into local partitions.  Then the internal
// resolution matrix hii is factorized and the static condensation completed.
// In addition, the interior-exterior partition hbi is postmultiplied
// by hii(inverse) for convenience in the resolution stage.
//
// hbb:    nExt  by nExt    matrix;  (row-major 1D storage).
// hbi:    nExt  by nInt    matrix;  (row-major 1D storage).
// hii:    nInt  by nInt    matrix;  (packed-symmetric 1D storage).
// rmat:   nKnot by nKnot   matrix;  (row-major 1D storage).
// rwrk:   nExt*(nExt+nInt) vector;
// ---------------------------------------------------------------------------
{
  char      routine[] = "Element::HelmholtzSC";
  register  int i, j, k;
  int       eq, ij = 0;
  int       ntot = nTot();
  int       next = nExt(),  nint = nInt(),  ipack;
  real    **IT, **DV, **DT;

  // -- Allocate element internal/external resolution matrices.

  if (nint) {
    ipack = ((nint + 1) * nint) >> 1;
    hii  = rvector (ipack);            // LAPACK packed-symmetric store.
    hbi  = rvector (next*nint);        // Full store.
  }

  // -- Construct hbb, hbi, hii partitions of elemental Helmholtz matrix.

  quadOps (rule, np, nq, 0, 0, 0, 0, &IT, &DV, &DT);

  for (i = 0; i < np; i++)
    for (j = 0; j < np; j++, ij++) {

      helmRow (IT, DV, DT, lambda2, i, j, rmat, rwrk, rwrk+nq, rwrk+nq+nq);

      Veclib::gathr (ntot, rmat, emap, rwrk);

      if ( (eq = pmap[ij]) < next ) {
	Veclib::copy (next, rwrk, 1, hbb + eq * next, 1);
	if (nint) Veclib::copy (nint, rwrk+next, 1, hbi + eq*nint, 1);
      } else
	for (k = eq; k < ntot; k++)
	  hii[Lapack::pack_addr (eq - next, k - next)] = rwrk[k];
    }

#ifdef DEBUG
  if (option ("VERBOSE") > 3) printMatSC (hbb);
#endif

  // -- Carry out static condensation step.

  if (nint) {
    int  info = 0;

    // -- Factor hii.

    Lapack::pptrf ("U", nint, hii, info);
    if (info) message (routine, "dpptrf failed to factor hii", ERROR);

    // -- Statically condense hbb.

    Veclib::copy  (nint*next, hbi, 1, rwrk, 1);
    Lapack::pptrs ("U", nint, next, hii, rwrk, nint, info);
    Blas::gemm    ("T", "N", next, next, nint, -1.0, hbi,
		   nint, rwrk, nint, 1.0, hbb, next); 

    // -- Create hib*hii(inverse), leave in hbi.

    Lapack::pptrs ("U", nint, next, hii, hbi, nint, info);
  }
}


void Element::HelmholtzSC (const Element& master)
// ---------------------------------------------------------------------------
// Set stored Helmholtz matrix partitions from master.
// ---------------------------------------------------------------------------
{
  hii = master.hii;
  hbi = master.hbi;
}


void  Element::printMatSC (const real* hbb) const
// ---------------------------------------------------------------------------
// (Debugging) utility to print up element matrices.
// ---------------------------------------------------------------------------
{
  char  s[8*StrMax];
  int   i, j, next = nExt(), nint = nInt(), ntot = nTot();

  sprintf (s, "-- Uncondensed Helmholtz matrices, element %1d", id);
  message ("", s, REMARK);

  sprintf (s, "-- hbb:");
  message ("", s, REMARK);

  for (i = 0; i < next; i++) {
    for (j = 0; j < next; j++)
      cout << setw (10) << hbb[Veclib::row_major (i, j, next)];
    cout << endl;
  }

  sprintf (s, "-- hii:");
  message ("", s, REMARK);

  for (i = 0; i < nint; i++) {
    for (j = 0; j < nint; j++)
      cout << setw (10) << hii[Lapack::pack_addr (i, j)];
    cout << endl;
  }

  sprintf (s, "-- hbi:");
  message ("", s, REMARK);

  for (i = 0; i < next; i++) {
    for (j = 0; j < nint; j++)
      cout << setw (10) << hbi[Veclib::row_major (i, j, nint)];
    cout << endl;
  }
}


void Element::Helmholtz (const real&  lambda2,
			 real*        h      ,
			 real*        rmat   ,
			 real*        rwrk   ) const
// ---------------------------------------------------------------------------
// Compute the discrete elemental Helmholtz matrix, return in h.
//
// This routine can be used when static condensation is not employed, and is
// included mainly to ease checking of entire element matrices.
// Node ordering produced is row-major.
//
// h:    vector, length np*np*np*np;
// rmat: vector, length np*np;
// rwrk: vector, length nq*nq + 2*nq.
// ---------------------------------------------------------------------------
{
  register  int ij = 0;
  int       ntot   = nTot ();
  real    **IT, **DV, **DT;

  quadOps (rule, np, nq, 0, 0, 0, 0, &IT, &DV, &DT);

  for (register int i = 0; i < np; i++)
    for (register int j = 0; j < np; j++, ij++) {
      helmRow (IT, DV, DT, lambda2, i, j, rmat, rwrk, rwrk+nq, rwrk+nq+nq);
      Veclib::copy (ntot, rmat, 1, h + ij * np, 1);
    }
}


void Element::helmRow (real**  IT     ,
		       real**  DV     ,
		       real**  DT     ,
		       const real&    lambda2,
		       const int&     i      ,
		       const int&     j      ,
		       real*          hij    ,
		       real*          W1     ,
		       real*          W2     ,
		       real*          W0     ) const
// ---------------------------------------------------------------------------
// Build row [i,j] of the elemental Helmholtz matrix in array hij (np x np).
//
// Three work arrays W1, W2 and W0 are nominated.  Only W1 is used if the
// element quadrature rule is LL (Lobatto), but all three are used for GL
// (Gauss) quadrature.  W1 and W2 should be at least nq long; W0 should be
// nq x nq long.
//
// For a 2D tensor product form, the elemental Helmholtz matrix is produced
/// as (sums on p & q indices assumed):
//
// h      =         G1  IN  DT  IN  DT     \
//  ij mn             pq  pi  jq  pm  nq   |
//        +         G2  DV  IT  DV  IT     |
//                    pq  pi  jq  pm  nq   |
//        +         G3  DV  IT  IN  DT      >  "STIFFNESS"
//                    pq  pi  jq  pm  nq   |
//        +         G3  IN  DT  DV  IT     |
//                    pq  pi  jq  pm  nq   /
//                2
//        + lambda  G4  IN  IT  IN  IT        "MASS"
//                    pq  pi  jq  pm  nq
//
// where the terms G1, G2, G3, G4 contain geometric mapping factors and
// quadrature weights, and the matrices IN, IT are the Lagrangian
// interpolation matrix (from the nodes to the quadrature points) and its
// transpose, while DV, DT are the Lagrangian derivative matrix & transpose.
// If the Lobatto quadrature rule is used, the interpolant matrices are
// identities (and may be input as NULL pointers).
// ---------------------------------------------------------------------------
{
  register  int m, n;
  int       ntot = sqr (nq);

  if (rule == LL) {		// -- Lobatto quadrature:

    Veclib::zero (ntot, hij, 1);

    for (n = 0; n < nq; n++) {
      Veclib::vmul (nq, DT[j], 1, DT[n], 1, W1, 1);
      hij[Veclib::row_major (i, n, nq)]  = Blas::dot (nq, G1 + i*nq, 1, W1, 1);
    }

    for (m = 0; m < nq; m++) {
      Veclib::vmul (nq, DT[i], 1, DT[m], 1, W1, 1);
      hij[Veclib::row_major (m, j, nq)] += Blas::dot (nq, G2 + j,   nq, W1, 1);
    }

    if (G3)
      for (m = 0; m < nq; m++)
	for (n = 0; n < nq; n++) {
	  hij [Veclib::row_major (m, n, nq)] +=
	    G3[Veclib::row_major (i, n, nq)] * DV[n][j] * DV[i][m];
	  hij [Veclib::row_major (m, n, nq)] +=
	    G3[Veclib::row_major (m, j, nq)] * DV[j][n] * DV[m][i];
      }

    hij[Veclib::row_major (i, j, nq)] +=
      lambda2 * G4[Veclib::row_major (i, j, nq)];
 
  } else {			// -- Gauss quadrature:

    for (m = 0; m < np; m++)
      for (n = 0; n < np; n++) {

	Veclib::zero (ntot, W0, 1);
	Veclib::vmul (nq, IT[i], 1, IT[m], 1, W1, 1);
	Veclib::vmul (nq, DT[j], 1, DT[n], 1, W2, 1);
	Blas::ger    (nq, nq, 1.0, W2, 1, W1, 1, W0, nq);
	hij[Veclib::row_major (m, n, np)] = Blas::dot (ntot, G1, 1, W0, 1);

	Veclib::zero (ntot, W0, 1);
	Veclib::vmul (nq, DT[i], 1, DT[m], 1, W1, 1);
	Veclib::vmul (nq, IT[j], 1, IT[n], 1, W2, 1);
	Blas::ger    (nq, nq, 1.0, W2, 1, W1, 1, W0, nq);
	hij[Veclib::row_major (m, n, np)] += Blas::dot (ntot, G2, 1, W0, 1);

	if (G3) {
	  Veclib::zero (ntot, W0, 1);
	  Veclib::vmul (nq, DT[i], 1, IT[m], 1, W1, 1);
	  Veclib::vmul (nq, IT[j], 1, DT[n], 1, W2, 1);
	  Blas::ger    (nq, nq, 1.0, W2, 1, W1, 1, W0, nq);
	  hij[Veclib::row_major (m, n, np)] += Blas::dot (ntot, G3, 1, W0, 1);

	  Veclib::zero (ntot, W0, 1);
	  Veclib::vmul (nq, IT[i], 1, DT[m], 1, W1, 1);
	  Veclib::vmul (nq, DT[j], 1, IT[n], 1, W2, 1);
	  Blas::ger    (nq, nq, 1.0, W2, 1, W1, 1, W0, nq);
	  hij[Veclib::row_major (m, n, np)] += Blas::dot (ntot, G3, 1, W0, 1);
	}

	Veclib::zero (ntot, W0, 1);
	Veclib::vmul (nq, IT[i], 1, IT[m], 1, W1, 1);
	Veclib::vmul (nq, IT[j], 1, IT[n], 1, W2, 1);
	Blas::ger    (nq, nq, 1.0, W2, 1, W1, 1, W0, nq);
	hij[Veclib::row_major (m, n, np)] +=
	  lambda2 * Blas::dot (ntot, G4, 1, W0, 1);
      }
  }
}


void  Element::assemble (const real*  hbb    ,
			 real*        Hp     ,
			 real*        Hc     ,
			 const int    nsolve ,
			 const int    ncons  ,
			 const int    nband  ) const
// ---------------------------------------------------------------------------
// Post elemental external node matrix hbb into global Helmholtz matrix Hp
// and constraint partition Hc using direct stiffness summation.
// ---------------------------------------------------------------------------
{
  if (!nsolve) return;

  register int  i, j, m, n;
  int           next = nExt ();

  if (nband) {		          // -- LAPACK upper triangular band storage.

    for (i = 0; i < next; i++)
      if ((m = bmap[i]) < nsolve)
	for (j = 0; j < next; j++)
	  if ((n = bmap[j]) < nsolve && n >= m)
	    Hp[Lapack::band_addr (m, n, nband)] +=
	      hbb[Veclib::row_major (i, j, next)];

  } else {			 // -- LAPACK upper triangular packed storage.

    for (i = 0; i < next; i++)
      if ((m = bmap[i]) < nsolve)
	for (j = 0; j < next; j++)
 	  if ((n = bmap[j]) < nsolve && n >= m)
	    Hp[Lapack::pack_addr (m, n)] += 
	      hbb[Veclib::row_major (i, j, next)];
  }

  if (ncons > 1)       		 // -- A constraint partition must exist...
    for (i = 0; i < next; i++)
      if ((m = bmap[i]) < nsolve)
	for (j = 0; j < next; j++)
	  if ((n = bmap[j]) >= nsolve)
	    Hc[Veclib::row_major (m, n - nsolve, ncons)] +=
	      hbb[Veclib::row_major (i, j, next)];
}


void  Element::grad (real* targA, real* targB) const
// ---------------------------------------------------------------------------
// Operate partial derivative d(target)/dxi = d_dr*drdxi + d_ds*dsdxi,
// where the appropriate component of gradient is selected by input pointers.
// Values are computed at node points.
// ---------------------------------------------------------------------------
{
  real  **DV, **DT;
  quadOps (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  int   ntot = nTot();
  real* tmpA = rvector (ntot);
  real* tmpB = rvector (ntot);
  real* target;

  if ((target = targA)) {
    Blas::mxm     (target, np, *DT,    np, tmpA, np);
    Blas::mxm     (*DV,    np, target, np, tmpB, np);
    Veclib::vmul  (ntot, tmpA, 1, drdx, 1, tmpA, 1);
    Veclib::vvtvp (ntot, tmpB, 1, dsdx, 1, tmpA, 1, target, 1);
  }

  if ((target = targB)) {
    Blas::mxm     (target, np, *DT,    np, tmpA, np);
    Blas::mxm     (*DV,    np, target, np, tmpB, np);
    Veclib::vmul  (ntot, tmpA, 1, drdy, 1, tmpA, 1);
    Veclib::vvtvp (ntot, tmpB, 1, dsdy, 1, tmpA, 1, target, 1);
  }
  
  freeVector (tmpA);
  freeVector (tmpB);
}


void  Element::connectivSC (List<int>* adjList) const
// ---------------------------------------------------------------------------
// AdjList is an array of linked lists, each of which describes the global
// nodes that have connectivity with the the current node, i.e. which make
// a contribution to the weighted-residual integral for this node.
// This routine fills in the contribution from the current element.
//
// For general finite elements, all nodes of an element are interconnected,
// while for statically-condensed elements, only the boundary nodes are
// considered (since internal nodes are not global).
//
// Essential-BC nodes are ignored, since we're only interested in mimimizing
// bandwidths of global matrices.
// ---------------------------------------------------------------------------
{
  int           next = nExt();
  int           found, gidCurr, gidMate;
  register int  i, j;
    
  for (i = 0; i < next; i++)
    if (solve[i]) {
      gidCurr = bmap[i];

      for (j = 0; j < next; j++)
	if (i != j && solve[j]) {
	  gidMate = bmap[j];

	  found = 0;
	  for (ListIterator<int>
	       a (adjList[gidCurr]); !found && a.more (); a.next ())
	    found = a.current () == gidMate;
	  if (!found) adjList[gidCurr].add (gidMate);
	}
    }
}


void  Element::sideOffset (int side, int& start, int& skip) const
// ---------------------------------------------------------------------------
// Return starting offset & skip in Field storage for this side.
// ---------------------------------------------------------------------------
{
  int estart, bstart;
  terminal (side, estart, skip, bstart);
  start = offset + estart;
}


void  Element::sideGeom (int side, real* nx, real* ny, real* area) const
// ---------------------------------------------------------------------------
// Generate unit outward normal components and change-of-variable
// Jacobian, area, for use in computation of edge integrals.
//
// We will always use Lobatto-Legendre quadrature for these integrals; 
// however, we need to do some recomputation of local forward partial
// derivatives along edges.
//
// Computed vectors have edge-traverse ordering, i.e. are made to operate
// on vectors obtained from base storage using BLAS-conformant copy.
// ---------------------------------------------------------------------------
{
  if (side < 1 || side > ns)
    message ("Element::sideGeom", "illegal side", ERROR);

  int     low, skip;
  real  **D, *w, *xr, *xs, *yr, *ys, *len;

  quadOps (LL, np, np, 0, 0, &w, 0, 0, &D, 0);

  switch (side) {
  case 1: 
    skip = 1;
    xr   = rvector (np);
    yr   = rvector (np);
    
    Blas::gemv    ("T", np, np, 1.0, *D, np, *xmesh, 1, 0.0, xr, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, *ymesh, 1, 0.0, yr, 1);
    Veclib::vmul  (np, xr, 1, xr, 1, area, 1);
    Veclib::vvtvp (np, yr, 1, yr, 1, area, 1, area, 1);
    Veclib::smul  (np, -1.0, dsdx, skip, nx, 1);
    Veclib::smul  (np, -1.0, dsdy, skip, ny, 1);
    
    freeVector (xr);
    freeVector (yr);
    break;

  case 2: 
    low  = np - 1;
    skip = np;
    xs   = rvector (np);
    ys   = rvector (np);
      
    Blas::gemv    ("T", np, np, 1.0, *D, np, *xmesh+low, np, 0.0, xs, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, *ymesh+low, np, 0.0, ys, 1);
    Veclib::vmul  (np, xs, 1, xs, 1, area, 1);
    Veclib::vvtvp (np, ys, 1, ys, 1, area, 1, area, 1);
    Veclib::copy  (np, drdx+low, skip, nx, 1);
    Veclib::copy  (np, drdy+low, skip, ny, 1);
    
    freeVector (xs);
    freeVector (ys);
    break;

  case 3:
    low  = np * (np - 1);
    skip = -1;
    xr   = rvector (np);
    yr   = rvector (np);
	
    Blas::gemv    ("T", np, np, 1.0, *D, np, *xmesh+low, 1, 0.0, xr, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, *ymesh+low, 1, 0.0, yr, 1);
    Veclib::vmul  (np, xr, 1, xr, 1, area, 1);
    Veclib::vvtvp (np, yr, 1, yr, 1, area, 1, area, 1);
    Veclib::copy  (np, dsdx+low, skip, nx, 1);
    Veclib::copy  (np, dsdy+low, skip, ny, 1);
    
    freeVector (xr);
    freeVector (yr);
    break;

  case 4:
    skip = -np;
    xs   = rvector (np);
    ys   = rvector (np);
      
    Blas::gemv    ("T", np, np, 1.0, *D, np, *xmesh, np, 0.0, xs, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, *ymesh, np, 0.0, ys, 1);
    Veclib::vmul  (np, xs, 1, xs, 1, area, 1);
    Veclib::vvtvp (np, ys, 1, ys, 1, area, 1, area, 1);
    Veclib::smul  (np, -1.0, drdx, skip, nx, 1);
    Veclib::smul  (np, -1.0, drdy, skip, ny, 1);
    
    freeVector (xs);
    freeVector (ys);
    break;
  }
  
  Veclib::vsqrt (np, area, 1, area, 1);
  Veclib::vmul  (np, area, 1, w,    1, area, 1);

  len = rvector (np);

  Veclib::vhypot (np, nx, 1, ny,  1, len, 1);
  Veclib::vdiv   (np, nx, 1, len, 1, nx,  1);
  Veclib::vdiv   (np, ny, 1, len, 1, ny,  1);

  freeVector (len);
}


void  Element::sideEval (int side, real* target, const char* func) const
// ---------------------------------------------------------------------------
// Evaluate function func along side of element, returning in target.
//
// The function can use variables "x", "y" & "t" (and any floating-point
// parameters previously set).
// ---------------------------------------------------------------------------
{
  int estart, skip, bstart;
  terminal (side, estart, skip, bstart);

  real *x = rvector (np);
  real *y = rvector (np);

  Veclib::copy (np, *xmesh + estart, skip, x, 1);
  Veclib::copy (np, *ymesh + estart, skip, y, 1);

  vecInit   ("x y", func);
  vecInterp (np, x, y, target);

  freeVector (x);
  freeVector (y);
}


void  Element::sideScatr (int side, const real* src, real* target) const
// ---------------------------------------------------------------------------
// Scatter vector src into globally-numbered target vector along this side.
// ---------------------------------------------------------------------------
{
  int estart, skip, bstart;
  terminal (side, estart, skip, bstart);

  if (side == ns) {
    Veclib::scatr (np - 1, src, bmap + bstart, target);
    target[bmap[0]] = src[np - 1];
  }
  else
    Veclib::scatr (np,     src, bmap + bstart, target);
}


void  Element::sideDsSum (int         side  ,
			  const real* src   ,
			  const real* area  ,
			  real*       target) const
// ---------------------------------------------------------------------------
// Direct-stiffness-sum (area-weighted) vector src into globally-numbered
// target on side.  This is for evaluation of ESSENTIAL BCs with
// Lobatto quadrature.
//
// Complication at side ends, since NATURAL BCs defer to ESSENTIAL BCs.
// ---------------------------------------------------------------------------
{
  register int  i, nm = np - 1;
  real*         tmp   = rvector (np);
  
  int estart, skip, bstart;
  terminal (side, estart, skip, bstart);

  Veclib::vmul (np, src, 1, area, 1, tmp, 1);

  if (solve[bstart])          target[bmap[bstart  ]] += tmp[0];

  for (i = 1; i < nm; i++)    target[bmap[bstart+i]] += tmp[i];

  if (side == ns && solve[0]) target[bmap[0       ]] += tmp[i];
  else if (solve[bstart+i])   target[bmap[bstart+i]] += tmp[i];

  freeVector (tmp);
}


void  Element::sideMask (int side, int* gmask) const
// ---------------------------------------------------------------------------
// Switch off globally-numbered solve mask (for ESSENTIAL-BC side).
// ---------------------------------------------------------------------------
{
  register int i, nm = np - 1;

  --side;

  for (i = 0; i < nm; i++) gmask[bmap[side*nm + i]] &= 0;
  gmask[bmap[nm*(side + 1) % nExt()]]               &= 0;
}


void  Element::sideMask (int side, real* gnode) const
// ---------------------------------------------------------------------------
// Switch off globally-numbered gnode at ESSENTIAL BC locations.
// ---------------------------------------------------------------------------
{
  register int i, nm = np - 1;

  --side;

  for (i = 0; i < nm; i++) gnode[bmap[side * nm + i]] = 0.0;
  gnode[bmap[nm*(side + 1) % nExt()]]             = 0.0;
}


void Element::bndryMask (real* essl, real* othr) const
// ---------------------------------------------------------------------------
// Switch off globally-numbered gnode at ESSENTIAL BC locations, or at other.
// ---------------------------------------------------------------------------
{
  register int  i, next = nExt ();

  if (essl) for (i = 0; i < next; i++) essl[bmap[i]] *= solve[i] ? 1.0 : 0.0;
  if (othr) for (i = 0; i < next; i++) othr[bmap[i]] *= solve[i] ? 0.0 : 1.0;
}


void  Element::sideGrad (int side, const real* src, real* c1, real* c2 ) const
// ---------------------------------------------------------------------------
// Using geometric factors for this Element, return the first and second
// component, c1 and/or c2, of grad src (length nTot()) along side.
//
// We have to take some special care on  sides 3 & 4, where the usual skips
// are negative: we instead use positive skips for formation of dc/dr, dc/ds,
// then a -1 skip when multiplying by dr/dx, ds/dx, etc.
// ---------------------------------------------------------------------------
{
  int estart, skip, bstart, d;
  terminal (side, estart, skip, bstart);

  real*  ddr = rvector (np);
  real*  dds = rvector (np);

  real  **DV, **DT;

  quadOps (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  // -- Make dc/dr, dc/ds along element edge.

  switch (side) {
  case 1:
    d = 1;
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("N", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, dds, 1);
    break;
  case 2:
    d = 1;
    Blas::gemv ("T", np, np, 1.0, src, np, *DT + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, dds, 1);
    break;
  case 3:
    d = -1;
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("N", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, dds, 1);
    break;
  case 4:
    d = -1;
    Blas::gemv ("T", np, np, 1.0, src, np, *DT + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, dds, 1);
    break;
  }

  // -- dc/dx = dc/dr * dr/dx + dc/ds * ds/dx.

  if (c1) {
    Veclib::vmul  (np, ddr, d, drdx + estart, skip, c1, 1);
    Veclib::vvtvp (np, dds, d, dsdx + estart, skip, c1, 1, c1, 1);
  }
  
  // -- dc/dy = dc/dr * dr/dy + dc/ds * ds/dy.

  if (c2) {
    Veclib::vmul  (np, ddr, d, drdy + estart, skip, c2, 1);
    Veclib::vvtvp (np, dds, d, dsdy + estart, skip, c2, 1, c2, 1);
  }

  freeVector (ddr);
  freeVector (dds);

/*
  int estart, skip, bstart, d;
  terminal (side, estart, skip, bstart);

  real*  ddr = rvector (np);
  real*  dds = rvector (np);

  real  **DV, **DT;

  quadOps (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  // -- Make dc/dr, dc/ds along element edge.

  switch (side) {
  case 1:
    d = 1;
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("N", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, dds, 1);
    break;
  case 2:
    d = 1;
    Blas::gemv ("T", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, dds, 1);
    break;
  case 3:
    d = -1;
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("N", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, dds, 1);
    break;
  case 4:
    d = -1;
    Blas::gemv ("T", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, dds, 1);
    break;
  }

  // -- dc/dx = dc/dr * dr/dx + dc/ds * ds/dx.

  if (c1) {
    Veclib::vmul  (np, ddr, d, drdx + estart, skip, c1, 1);
    Veclib::vvtvp (np, dds, d, dsdx + estart, skip, c1, 1, c1, 1);
  }
  
  // -- dc/dy = dc/dr * dr/dy + dc/ds * ds/dy.

  if (c2) {
    Veclib::vmul  (np, ddr, d, drdy + estart, skip, c2, 1);
    Veclib::vvtvp (np, dds, d, dsdy + estart, skip, c2, 1, c2, 1);
  }

  freeVector (ddr);
  freeVector (dds);
*/
}


void  Element::evaluate (const char* function, real* target) const
// ---------------------------------------------------------------------------
// Evaluate function over mesh points, store in target.
// Function can explicitly use "x" and "y", for which mesh values are used.
// ---------------------------------------------------------------------------
{
  vecInit   ("x y", function);
  vecInterp (nTot(), *xmesh, *ymesh, target);
}


real  Element::area () const
// ---------------------------------------------------------------------------
// Return area of element, using element quadrature rule.
// ---------------------------------------------------------------------------
{
  return  Veclib::sum (nq*nq, G4, 1);
}


real  Element::integral (const char* function) const
// ---------------------------------------------------------------------------
// Return integral of function over element, using element quadrature rule.
// ---------------------------------------------------------------------------
{
  real   intgrl;
  int    ntot = nq*nq;
  real*  tmp  = rvector (ntot);

  vecInit ("x y", function);

  switch (rule) {

  case LL:
    vecInterp (ntot, *xmesh, *ymesh, tmp);
    break;

  case GL:
    real   *x   = rvector (ntot),
           *y   = rvector (ntot),
           *wrk = rvector (nq*np);
    real  **IN, **IT;

    quadOps (GL, np, nq, 0, 0, 0, &IN, &IT, 0, 0);

    Blas::mxm (*IN, nq, *xmesh, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,    np, x,   nq);

    Blas::mxm (*IN, nq, *ymesh, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,    np, y,   nq);

    vecInterp (ntot, x, y, tmp);

    freeVector (x);
    freeVector (y);
    freeVector (wrk);
    break;
  }

  Veclib::vmul (ntot, tmp, 1, G4, 1, tmp, 1);

  intgrl = Veclib::sum (ntot, tmp, 1);
  
  freeVector (tmp);

  return intgrl;
}


real  Element::norm_inf (const real* src) const
// ---------------------------------------------------------------------------
// Return infinity-norm of element value.
// ---------------------------------------------------------------------------
{
  return  fabs (src[Blas::iamax (nTot (), src, 1)]);
}


real  Element::norm_L2 (const real* src) const
// ---------------------------------------------------------------------------
// Return L2-norm of Element value, using Element quadrature rule.
// ---------------------------------------------------------------------------
{
  register real   L2 = 0.0;
  register int    i, ntot = nq*nq;
  register real*  dA;

  switch (rule) {

  case LL:
    dA = G4;

    for (i = 0; i < ntot; i++) L2 += src[i] * src[i] * dA[i];

    break;

  case GL:
    real*            wrk = rvector (nq*np);
    register real*   u   = rvector (ntot);
    real**           IN;
    real**           IT;
    
    dA = G4;

    quadOps (GL, np, nq, 0, 0, 0, &IN, &IT, 0, 0);

    Blas::mxm (*IN, nq,  src, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,  np, u,   nq);

    for (i = 0; i < ntot; i++) L2 += u[i] * u[i] * dA[i];
    
    freeVector (wrk);
    freeVector (u);

    break;
  }

  return sqrt (L2);
}


real  Element::norm_H1 (const real* src) const
// ---------------------------------------------------------------------------
// Return Sobolev-1 norm of Element value, using Element quadrature rule.
// ---------------------------------------------------------------------------
{
  register real   H1 = 0;
  register int    i, ntot = nq*nq;
  register real*  u;
  register real*  gw;

  switch (rule) {

  case LL:
    gw = G4;

    // -- Add in L2 norm of u.

    for (i = 0; i < ntot; i++) H1 += src[i] * src[i] * gw[i];

    // -- Add in L2 norm of grad u.

    u = rvector (ntot);
    
    Veclib::copy (ntot, src, 1, u, 1);
    grad (u, 0);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * gw[i];

    Veclib::copy (ntot, src, 1, u, 1);
    grad (0, u);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * gw[i];

    freeVector (u);
    break;

  case GL:
    real*   wrk = rvector (nq*np);
    real**  IN;
    real**  IT;
    real**  DV;
    real**  DT;

    quadOps (GL, np, nq, 0, 0, 0, &IN, &IT, &DV, &DT);

    Blas::mxm (*IN, nq,  src, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,  np, u,   nq);

    // -- Add in L2 norm of u.

    gw = G4;
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * gw[i];

    // -- Add in L2 norm of grad u.

    gw = G1;
    Blas::mxm (*IN, nq,  src, np, wrk, np);
    Blas::mxm (wrk, nq, *DV,  np, u,   nq);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * gw[i];

    gw = G2;
    Blas::mxm (*DV, nq,  src, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,  np, u,   nq);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * gw[i];

    if (G3) {
      gw = G3;
      Blas::mxm (*DV, nq,  src, np, wrk, np);
      Blas::mxm (wrk, nq, *DT,  np, u,   nq);
      for (i = 0; i < ntot; i++) H1 += 2.0 * u[i] * u[i] * gw[i];
    }
    
    freeVector (wrk);
    freeVector (u);

    break;
  }

  return sqrt (H1);
}


void  Element::printMesh () const
// ---------------------------------------------------------------------------
// Print mesh information in prism-compatible form.
// ---------------------------------------------------------------------------
{
  printVector (cout, "rr", nTot(), *xmesh, *ymesh);
}


void  Element::economizeGeo ()
// ---------------------------------------------------------------------------
// Insert element geometric factors in economized storage.
// ---------------------------------------------------------------------------
{
  int  size;

  size = nTot ();

  drdx = FamilyMgr::insert (size, drdx);
  dsdx = FamilyMgr::insert (size, dsdx);
  drdy = FamilyMgr::insert (size, drdy);
  dsdy = FamilyMgr::insert (size, dsdy);

  size = sqr (nQuad ());
  
  G1   = FamilyMgr::insert (size, G1);
  G2   = FamilyMgr::insert (size, G2);
  G3   = FamilyMgr::insert (size, G3);

  if (G4 == mass) {
    mass = FamilyMgr::insert (size, mass);
    G4   = mass;
  } else {
    mass = FamilyMgr::insert (nTot (), mass);
    G4   = FamilyMgr::insert (size,    G4  );
  }
}


void  Element::economizeMat ()
// ---------------------------------------------------------------------------
// Insert element static condensation matrices in economized storage.
// ---------------------------------------------------------------------------
{
  int  nint = nInt (), next = nExt ();
  int  size;

  if (!nint) return;

  size = ((nint + 1) * nint) >> 1;
  
  hii  = FamilyMgr::insert (size, hii);

  size = next * nint;

  hbi  = FamilyMgr::insert (size, hbi);
}


void  Element::condense (const real* src, real* external, real* internal) const
// ---------------------------------------------------------------------------
// Load this element's boundary data into globally-numbered external, and
// internal data into un-numbered internal.
// ---------------------------------------------------------------------------
{
  register int  next = nExt ();

  Veclib::gathr_scatr (next,    src, emap,  bmap, external);
  Veclib::gathr       (nInt (), src, emap + next, internal);
}


void  Element::expand (real* targ, const real* external, const real* internal)
// ---------------------------------------------------------------------------
// Load this element's boundary data from globally-numbered external, and
// internal data from un-numbered internal.
// ---------------------------------------------------------------------------
{
  register int  next = nExt ();

  Veclib::gathr_scatr (next,    external, bmap,  emap, targ);
  Veclib::scatr       (nInt (), internal, emap + next, targ);
}


void Element::dsSum (const real* src, real* external, real* internal) const
// ---------------------------------------------------------------------------
// Sum this src's boundary data into globally-numbered external, and
// internal data into un-numbered internal.
// ---------------------------------------------------------------------------
{
  register int  next = nExt ();

  Veclib::gathr_scatr_sum (next,    src, emap,  bmap, external);
  Veclib::gathr_sum       (nInt (), src, emap + next, internal);
}


void Element::HelmholtzOperator (const real* src, real* tgt, 
				 const real& L2,  real* wrk) const
// ---------------------------------------------------------------------------
// Apply elemental discrete Helmholtz operator on src to make tgt.
// This routine is specific to 2D Cartesian space with GLL integration.
//
// Workspace vector wrk must hold 2 * nTot () elements.
// ---------------------------------------------------------------------------
{
  register int   ij, nt = nTot ();
  register real  tmp;
  real*          R = wrk;
  real*          S = R + nTot ();
  real**         DV;
  real**         DT;

  quadOps (rule, np, nq, 0, 0, 0, 0, 0, &DV, &DT);

  Blas::gemm ("N", "N", np, np, np, 1.0, *DT, np, src, np, 0.0, R, np);
  Blas::gemm ("N", "N", np, np, np, 1.0, src, np, *DV, np, 0.0, S, np);

  if (G3) {
    for (ij = 0; ij < nt; ij++) {
      tmp      = R [ij];
      R  [ij]  = G1[ij] * R  [ij] + G3[ij] * S  [ij];
      S  [ij]  = G2[ij] * S  [ij] + G3[ij] * tmp;
      tgt[ij]  = G4[ij] * src[ij];
    }
  } else {
    for (ij = 0; ij < nt; ij++) {
      R  [ij] *= G1[ij];
      S  [ij] *= G2[ij];
      tgt[ij]  = G4[ij] * src[ij];
    }
  }

  Blas::gemm ("N", "N", np, np, np, 1.0,  S,  np, *DT, np, L2,  tgt, np);
  Blas::gemm ("N", "N", np, np, np, 1.0, *DV, np,  R,  np, 1.0, tgt, np);
}
