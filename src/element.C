//////////////////////////////////////////////////////////////////////////////
// element.C
//////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <Sem.h>


Element::Element (const int   i   ,
		  const Mesh& M   ,
		  const real* z   ,
		  const int   nk  ,
		  const int   doff,
		  const int   boff) :

		  id          (i   ),
		  ns          (4   ),
		  np          (nk  ),
		  dOffset     (doff),
		  bOffset     (boff)
// ---------------------------------------------------------------------------
// Create a new quad element, nk X nk.  Spacing along any size generated
// by mapping z (defined on domain [-1, 1], np points) onto side.
//
// Compute information for internal storage, and economize.
// ---------------------------------------------------------------------------
{
  char routine[] = "Element::Element";
  const int nk2  = sqr (nk);

  if (nk < 2) message (routine, "need > 2 knots for element edges", ERROR);

  rule = (int) Femlib::value ("RULE");
  nq = (rule == LL) ? np : Femlib::nquad (2, np);

  Femlib::buildMaps (np, &emap, &pmap);
  
  xmesh = new real [nk2];
  ymesh = new real [nk2];
  
  M.meshElmt (id, np, z, xmesh, ymesh);
  
  Femlib::adopt (nk2, &xmesh);
  Femlib::adopt (nk2, &ymesh);
 
  map ();
}


Element::~Element ()
// ---------------------------------------------------------------------------
// Clean up internal storage using Femlib family routines.
// ---------------------------------------------------------------------------
{
  Femlib::abandon (&xmesh);
  Femlib::abandon (&ymesh);

  Femlib::abandon (&drdx);
  Femlib::abandon (&dsdx);
  Femlib::abandon (&drdy);
  Femlib::abandon (&dsdy);

  Femlib::abandon (&G1  );
  Femlib::abandon (&G2  );
  Femlib::abandon (&G3  );
  Femlib::abandon (&G4  );

  Femlib::abandon (&mass);
}


void Element::map ()
// ---------------------------------------------------------------------------
// Generate geometric factors associated with mapping from 2D Cartesian to
// isoparametrically-mapped space:
//
//   dxdr, dydr,   = dx/dr,  dy/dr,  "Forward Partials"
//   dxds, dyds,   = dx/ds,  dy/ds.
//   drdx, drdy,   = dr/dx,  dr/dy,  "Inverse Partials"
//   dsdx, dsdy,   = ds/dx,  ds/dy.
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
// For a Gauss--Legendre quadrature rule, the inverse partials and mass
// matrix are returned for spatial locations at the mesh nodes, while
// the forward partials and other geometric factors are for spatial
// locations at the quadrature points.  In general, the amount of storage
// allocated for forward and inverse partials differ.
//
// For Lobatto--Legendre rule, everything is at the nodes.
//
// The inverse partials are retained in Element storage, to be used in
// element gradient operations (e.g. dP/dx = dP/dr * dr/dx + dP/ds * ds/dx)
// while the forward partials are retained in scrambled form (in combination
// with quadrature weights) as "geometric factors" G1--G4, to be used in
// element quadrature operations.
//
// Null-mapping optimizations mentioned below occur when the element geometry
// ensures that the entries of a vector are zero to within roundoff, due
// either to the edges of elements being aligned with coordinate axes
// (as can happen for the inverse partials) or if the element is an
// undistorted (but possibly rotated) rectangle (this applies to G3).
// In these cases the associated memory is deleted and the pointers
// are set to zero, so they can serve as flags in subsequent computations.
// ---------------------------------------------------------------------------
{
  char  routine[] = "Element::map";
  char  err[StrMax];
  int   ntot;
  real  *x = xmesh, *y = ymesh;
  real  **DV, **DT,  **IN,  **IT;
  real  *jac, *dxdr, *dxds, *dydr, *dyds, *tM, *tV, *w, *WW;

  vector<real> work;
  const real   EPS = 4 * nTot()*((sizeof(real)==sizeof(double)) ? EPSDP:EPSSP);

  if (rule == LL) {		// -- Lobatto--Legendre quadrature.
    ntot = nTot();

    // -- Permanent/family allocations.

    drdx = new real [ntot];
    dsdx = new real [ntot];
    drdy = new real [ntot];
    dsdy = new real [ntot];
    G1   = new real [ntot];
    G2   = new real [ntot];
    G3   = new real [ntot];
    G4   = new real [ntot];
    
    // -- Temporaries.

    work.setSize (7 * ntot);

    dxdr = work();
    dxds = dxdr + ntot;
    dydr = dxds + ntot;
    dyds = dydr + ntot;

    jac  = dyds + ntot;
    WW   = jac  + ntot;
    tV   = WW   + ntot;
    
    Femlib::quad (LL, np, np, 0, 0, &w, 0, 0, &DV, &DT);
    Veclib::zero (ntot, WW, 1);
    Blas::ger    (np, np, 1.0, w, 1, w, 1, WW, np);
    
    Blas::mxm (  x, np, *DT, np, dxdr, np);
    Blas::mxm (*DV, np,   x, np, dxds, np);
    Blas::mxm (  y, np, *DT, np, dydr, np);
    Blas::mxm (*DV, np,   y, np, dyds, np);
    
    Veclib::vmul  (ntot,        dxdr, 1, dyds, 1, tV,  1);
    Veclib::vvvtm (ntot, tV, 1, dxds, 1, dydr, 1, jac, 1);
    
    if (jac[Veclib::imin (ntot, jac, 1)] < EPSSP) {
      sprintf (err, "jacobian of element %1d nonpositive", id);
      message (routine, err, ERROR);
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
    
    Veclib::vmul  (ntot, jac, 1, WW, 1, G4, 1);
    
    Veclib::copy (ntot, dyds, 1, drdx, 1);
    Veclib::vneg (ntot, dxds, 1, drdy, 1);
    Veclib::vneg (ntot, dydr, 1, dsdx, 1);
    Veclib::copy (ntot, dxdr, 1, dsdy, 1);
    
    Veclib::vdiv (ntot, drdx, 1, jac, 1, drdx, 1);
    Veclib::vdiv (ntot, drdy, 1, jac, 1, drdy, 1);
    Veclib::vdiv (ntot, dsdx, 1, jac, 1, dsdx, 1);
    Veclib::vdiv (ntot, dsdy, 1, jac, 1, dsdy, 1);

    // -- Calculations are done.  Do null-mapping optimizations.

    if (Blas::nrm2 (ntot, drdx, 1) < EPS) { delete [] drdx; drdx = 0; }
    if (Blas::nrm2 (ntot, drdy, 1) < EPS) { delete [] drdy; drdy = 0; }
    if (Blas::nrm2 (ntot, dsdx, 1) < EPS) { delete [] dsdx; dsdx = 0; }
    if (Blas::nrm2 (ntot, dsdy, 1) < EPS) { delete [] dsdy; dsdy = 0; }
    if (Blas::nrm2 (ntot, G3,   1) < EPS) { delete [] G3;   G3   = 0; }

    // -- Check for family redundancies.

    Femlib::adopt (ntot, &drdx);
    Femlib::adopt (ntot, &drdy);
    Femlib::adopt (ntot, &dsdx);
    Femlib::adopt (ntot, &dsdy);
    Femlib::adopt (ntot, &G1  );
    Femlib::adopt (ntot, &G2  );
    Femlib::adopt (ntot, &G3  );
    Femlib::adopt (ntot, &G4  );

    // -- Install alias for mass matrix.

    mass = G4;

  } else {  // -- rule == GL: Gauss--Legendre quadrature.
    
    // -- Quadrature point computations.

    ntot = sqr (nq);

    G1 = new real [ntot];
    G2 = new real [ntot];
    G3 = new real [ntot];
    G4 = new real [ntot];
    
    work.setSize (7 * ntot + np * nq);

    dxdr = work();
    dxds = dxdr + ntot;
    dydr = dxds + ntot;
    dyds = dydr + ntot;
   
    jac  = dyds + ntot;
    WW   = jac  + ntot;
    tM   = WW   + ntot;
    tV   = tM + np * nq;

    Femlib::quad (GL, np, nq, 0, 0, &w, &IN, &IT, &DV, &DT);
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

    Veclib::vmul (ntot, jac, 1, WW, 1, G4, 1);

    // -- Do null-mapping optimization.

    if (Blas::nrm2 (ntot, G3, 1) < EPS) { delete [] G3; G3 = 0; }
    
    // -- Check family membership.

    Femlib::adopt (ntot, &G1);
    Femlib::adopt (ntot, &G2);
    Femlib::adopt (ntot, &G3);
    Femlib::adopt (ntot, &G4);

    // -- Node point computations.
    
    ntot = nTot();

    drdx = new real [ntot];
    drdy = new real [ntot];
    dsdx = new real [ntot];
    dsdy = new real [ntot];
    mass = new real [ntot];

    work.setSize (7 * ntot);

    dxdr = work();
    dxds = dxdr + ntot;
    dydr = dxds + ntot;
    dyds = dydr + ntot;
    jac  = dyds + ntot;
    WW   = jac  + ntot;
    tV   = WW   + ntot;
    
    Femlib::quad (LL, np, np, 0, 0, &w, 0, 0, &DV, &DT);
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

    // -- Do null-mapping optimizations.

    if (Blas::nrm2 (ntot, drdx, 1) < EPS) { delete [] drdx; drdx = 0; }
    if (Blas::nrm2 (ntot, drdy, 1) < EPS) { delete [] drdy; drdy = 0; }
    if (Blas::nrm2 (ntot, dsdx, 1) < EPS) { delete [] dsdx; dsdx = 0; }
    if (Blas::nrm2 (ntot, dsdy, 1) < EPS) { delete [] dsdy; dsdy = 0; }

    // -- Family membership.

    Femlib::adopt (ntot, &drdx);
    Femlib::adopt (ntot, &drdy);
    Femlib::adopt (ntot, &dsdx);
    Femlib::adopt (ntot, &dsdy);
    Femlib::adopt (ntot, &mass);
  }
}


void Element::bndryDsSum (const int*  btog,
			  const real* src ,
			  real*       tgt ) const
// ---------------------------------------------------------------------------
// Direct-stiffness-sum from element boundary to globally-numbered storage,
// i.e. tgt[btog[i]] += mass[emap[i]] * src[emap[i]].
// This is using in smoothing Fields along element boundaries.
// ---------------------------------------------------------------------------
{
  register int   i, e;
  const int      next = nExt();
  register real* wt = mass;

  for (i = 0; i < next; i++) {
    e = emap[i];
    tgt[btog[i]] += wt[e] * src[e];
  }
}


void Element::bndryMask (const int*  bmsk,
			 real*       tgt ,
			 const real* src ,
			 const int*  btog) const
// ---------------------------------------------------------------------------
// Mask the values in (row-major) tgt according to the mask vector bmsk
// and optionally globally-numbered vector src.
//
// If src is non-zero, it is used with boundary-to-global mapping vector
// btog to impose values within tgt on locations where bmsk is non-zero;
// other locations are unaffected.
//
// If src is zero, then the values within tgt where bmsk is zero are
// set to zero (i.e. tgt itself is taken as the source).  Btog is then
// not used and may be zero also.
//
// INCORPORATES/REPLACES old routines Element::mask and Element::setEssential.
// ---------------------------------------------------------------------------
{
  register int i, e;
  const int    ntot = nTot();
  const int    next = nExt();
  const int    nint = nInt();

  if (src) {
    for (i = 0; i < next; i++) {
      e = emap[i];
      tgt[e] = (bmsk[i]) ? src[btog[i]] : tgt[e];
    }

  } else {
    vector<real>   work (ntot);
    register real* tmp = work();

    Veclib::gathr (ntot, tgt, emap, tmp);
    for (i = 0; i < next; i++)
      tmp[i] = (bmsk[i]) ? tmp[i] : 0.0;
    Veclib::zero  (nint, tmp + next, 1);
    Veclib::gathr (ntot, tmp, pmap, tgt);
  }
}


void Element::e2gSumSC (real*       F   ,
			const int*  btog,
			real*       tgt ,
			const real* hbi ) const
// ---------------------------------------------------------------------------
// Create statically-condensed boundary Helmholtz forcing for this element
// from row-major F and insert it into globally-numbered tgt by direct
// stiffness summation.
//
// NB: forcing, F, is modified.
//
// On entry, F contains the elemental weak boundary-constrained forcing
//   - M f - H g.
//
// Elemental storage is then sorted in F so that it is ordered with boundary
// nodes foremost, i.e. it contains the partition { F | F }.
//                                                   b   i
//
// Statically-condensed boundary forcing is created in the first partition:
//                       -1                         -1
//   F   <--   F  -  h  h   F           (matrix h  h   supplied as hbi)
//    b         b     bi ii  i                   bi ii
//
// and summed into the tgt vector.  In the summation, there is no need
// to check if the global node is to be solved for or is fixed, since
// the fixed (essential-BC) partition of tgt is overwritten later.
//
// REPLACES old routine Element::dsForcingSC
// ---------------------------------------------------------------------------
{
  const int    ntot = nTot();
  const int    nint = nInt();
  const int    next = nExt();
  vector<real> work (ntot);

  Veclib::gathr (ntot, F, emap, work());
  Veclib::copy  (ntot, work(), 1, F, 1);

  if (nint) Blas::gemv ("T", nint,next, -1.0, hbi,nint, F + next,1, 1.0, F,1);

  Veclib::scatr_sum (next, F, btog, tgt);
}


void Element::g2eSC (const real* RHS ,
		     const int*  btog,
		     real*       F   , 
		     real*       tgt ,
		     const real* hbi ,
		     const real* hii ) const
// ---------------------------------------------------------------------------
// Complete static condensation solution for internal values of Element.
//
// On entry, global-node solution values are in RHS and F contains the
// weak form of internal forcing in its top end (as installed by e2gSumSC).
//
// If u is current Element, compute internal solution according to:
//            -1      -1
//   u  <--  h  F  - h  h   u
//    i       ii i    ii ib  b
//
// REPLACES old routine Element::resolveSC
// ----------------------------------------------------------------------------
{
  const int    next = nExt();
  const int    nint = nInt();
  vector<real> work (next);

  // -- Load globally-numbered RHS into element-boundary storage.

  Veclib::gathr (next, RHS,    btog, work());
  Veclib::scatr (next, work(), emap, tgt   );

  // -- Complete static-condensation solution for element-internal storage.

  if (nint) {
    int   info;
    real* Fi = F + next;

    Lapack::pptrs ("U", nint, 1, hii, Fi, nint, info);
    Blas::gemv    ("N", nint, next, -1.0, hbi, nint, work(), 1, 1.0, Fi, 1);
    Veclib::scatr (nint, Fi, emap + next, tgt);
  }
}


void Element::HelmholtzSC (const real lambda2,
			   real*      hbb    ,
			   real*      hbi    ,
			   real*      hii    ,
			   real*      rmat   ,
			   real*      rwrk   ) const
// ---------------------------------------------------------------------------
// Compute the discrete elemental Helmholtz matrix and return the
// statically condensed form in hbb, the interior-exterior coupling
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
  char         routine[] = "Element::HelmholtzSC";
  register int i, j, k;
  int          eq, ij = 0;
  const int    ntot = nTot();
  const int    next = nExt();
  const int    nint = nInt();
  real         **IT, **DV, **DT;

  // -- Construct hbb, hbi, hii partitions of elemental Helmholtz matrix.

  Femlib::quad (rule, np, nq, 0, 0, 0, 0, &IT, &DV, &DT);

  for (i = 0; i < np; i++)
    for (j = 0; j < np; j++, ij++) {

      helmRow ((const real**) IT, (const real**) DV, (const real**) DT,
	       lambda2, i, j, rmat, rwrk, rwrk+nq, rwrk+nq+nq);

      Veclib::gathr (ntot, rmat, emap, rwrk);

      if ( (eq = pmap[ij]) < next ) {
	Veclib::copy (next, rwrk, 1, hbb + eq * next, 1);
	if (nint) Veclib::copy (nint, rwrk+next, 1, hbi + eq*nint, 1);
      } else
	for (k = eq; k < ntot; k++)
	  hii[Lapack::pack_addr (eq - next, k - next)] = rwrk[k];
    }

#ifdef DEBUG
  if ((int) Femlib::value ("VERBOSE") > 3) printMatSC (hbb, hbi, hii);
#endif

  // -- Carry out static condensation step.

  if (nint) {
    int info = 0;

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


void Element::printMatSC (const real* hbb,
			  const real* hbi,
			  const real* hii) const
// ---------------------------------------------------------------------------
// (Debugging) utility to print up element matrices.
// ---------------------------------------------------------------------------
{
  char s[8*StrMax];
  int  i, j, next = nExt(), nint = nInt(), ntot = nTot();

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


void Element::Helmholtz (const real lambda2,
			 real*      h      ,
			 real*      rmat   ,
			 real*      rwrk   ) const
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
  register int ij   = 0;
  const int    ntot = nTot ();
  real         **IT, **DV, **DT;
  real         *wk1 = rwrk, *wk2 = wk1 + nq, *wk3 = wk2 + nq;

  Femlib::quad (rule, np, nq, 0, 0, 0, 0, &IT, &DV, &DT);

  for (register int i = 0; i < np; i++)
    for (register int j = 0; j < np; j++, ij++) {
      helmRow ((const real**) IT, (const real**) DV, (const real**)DT,
	       lambda2, i, j, rmat, wk1, wk2, wk3);
      Veclib::copy (ntot, rmat, 1, h + ij * np, 1);
    }
}


void Element::helmRow (const real** IT     ,
		       const real** DV     ,
		       const real** DT     ,
		       const real   lambda2,
		       const int    i      ,
		       const int    j      ,
		       real*        hij    ,
		       real*        W1     ,
		       real*        W2     ,
		       real*        W0     ) const
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
  register int m, n;
  const int    ntot = sqr (nq);

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


void Element::grad (real* tgtA,
		    real* tgtB) const
// ---------------------------------------------------------------------------
// Operate partial derivative d(tgt)/dxi = d_dr*drdxi + d_ds*dsdxi,
// where the appropriate component of gradient is selected by input pointers.
// Values are computed at node points.
// ---------------------------------------------------------------------------
{
  real  **DV, **DT;
  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  const int    ntot = nTot();
  vector<real> work (ntot + ntot);

  real* tmpA = work();
  real* tmpB = tmpA + ntot;
  real* tgt;

  if ((tgt = tgtA)) {
    if (drdx && dsdx) {
      Blas::mxm     (tgt, np, *DT, np, tmpA, np);
      Blas::mxm     (*DV, np, tgt, np, tmpB, np);
      Veclib::vmul  (ntot, tmpA, 1, drdx, 1, tmpA, 1);
      Veclib::vvtvp (ntot, tmpB, 1, dsdx, 1, tmpA, 1, tgt, 1);
    } else if (drdx) {
      Blas::mxm     (tgt, np, *DT, np, tmpA, np);
      Veclib::vmul  (ntot, tmpA, 1, drdx, 1, tgt, 1);
    } else {
      Blas::mxm     (*DV, np, tgt, np, tmpB, np);
      Veclib::vmul  (ntot, tmpB, 1, dsdx, 1, tgt, 1);
    }
  }

  if ((tgt = tgtB)) {
    if (drdy && dsdy) {
      Blas::mxm     (tgt, np, *DT, np, tmpA, np);
      Blas::mxm     (*DV, np, tgt, np, tmpB, np);
      Veclib::vmul  (ntot, tmpA, 1, drdy, 1, tmpA, 1);
      Veclib::vvtvp (ntot, tmpB, 1, dsdy, 1, tmpA, 1, tgt, 1);
    } else if (drdy) {
      Blas::mxm     (tgt, np, *DT, np, tmpA, np);
      Veclib::vmul  (ntot, tmpA, 1, drdy, 1, tgt, 1);
    } else {
      Blas::mxm     (*DV, np, tgt, np, tmpB, np);
      Veclib::vmul  (ntot, tmpB, 1, dsdy, 1, tgt, 1);
    }
  }
}


inline
void Element::terminal (const int side  ,
			int&      estart,
			int&      eskip ,
			int&      bstart) const
// ---------------------------------------------------------------------------
// Evaluate the element-edge terminal values of estart, skip, bstart.
// NB: BLAS-conformant terminal start values are delivered for negative skips.
//
// Side numbering starts at zero.
// ---------------------------------------------------------------------------
{
  switch (side) {
  case 0:
    estart = 0;
    eskip  = 1;
    bstart = 0;
    break;
  case 1:
    estart = np - 1;
    eskip  = np;
    bstart = np - 1;
    break;
  case 2:
    estart = np * (np - 1);
    eskip  = -1;
    bstart = 2  * (np - 1);
    break;
  case 3:
    estart = 0;
    eskip  = -np;
    bstart = 3  * (np - 1);
    break;
  }
}


void Element::sideOffset (const int side ,
			  int&      start,
			  int&      skip ) const
// ---------------------------------------------------------------------------
// Return starting offset & skip in Field storage for this side.
// ---------------------------------------------------------------------------
{
  int estart, bstart;

  terminal (side, estart, skip, bstart);
  start = dOffset + estart;
}


void Element::sideGeom (const int side,
			real*     nx  ,
			real*     ny  ,
			real*     area) const
// ---------------------------------------------------------------------------
// Generate unit outward normal components and change-of-variable
// Jacobian, area, for use in computation of edge integrals.
//
// We will always use Lobatto-Legendre quadrature for these integrals; 
// however, we need to do some recomputation of local forward partial
// derivatives along edges.
//
// Computed vectors have CCW edge-traverse ordering, i.e. are made to operate
// on vectors obtained from base storage using BLAS-conformant copy.
// ---------------------------------------------------------------------------
{
  if (side < 0 || side >= ns)
    message ("Element::sideGeom", "illegal side", ERROR);

  register int low, skip;
  real         **D, *w, *xr, *xs, *yr, *ys, *len;
  vector<real> work (np + np);

  Femlib::quad (LL, np, np, 0, 0, &w, 0, 0, &D, 0);

  switch (side) {
  case 0: 
    skip = 1;
    xr   = work();
    yr   = xr + np;
    
    Blas::gemv    ("T", np, np, 1.0, *D, np, xmesh, 1, 0.0, xr, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, ymesh, 1, 0.0, yr, 1);
    Veclib::vmul  (np, xr, 1, xr, 1, area, 1);
    Veclib::vvtvp (np, yr, 1, yr, 1, area, 1, area, 1);
    if   (dsdx) Veclib::smul (np, -1.0, dsdx, skip, nx, 1);
    else        Veclib::zero (np,                   nx, 1);
    if   (dsdy) Veclib::smul (np, -1.0, dsdy, skip, ny, 1);
    else        Veclib::zero (np,                   ny, 1);
    
    break;

  case 1: 
    low  = np - 1;
    skip = np;
    xs   = work();
    ys   = xs + np;
      
    Blas::gemv    ("T", np, np, 1.0, *D, np, xmesh+low, np, 0.0, xs, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, ymesh+low, np, 0.0, ys, 1);
    Veclib::vmul  (np, xs, 1, xs, 1, area, 1);
    Veclib::vvtvp (np, ys, 1, ys, 1, area, 1, area, 1);
    if   (drdx) Veclib::copy (np, drdx+low, skip, nx, 1);
    else        Veclib::zero (np,                 nx, 1);
    if   (drdy) Veclib::copy (np, drdy+low, skip, ny, 1);
    else        Veclib::zero (np,                 ny, 1);
    
    break;

  case 2:
    low  = np * (np - 1);
    skip = -1;
    xr   = work();
    yr   = xr + np;
	
    Blas::gemv    ("T", np, np, 1.0, *D, np, xmesh+low, 1, 0.0, xr, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, ymesh+low, 1, 0.0, yr, 1);
    Veclib::vmul  (np, xr, 1, xr, 1, area, 1);
    Veclib::vvtvp (np, yr, 1, yr, 1, area, 1, area, 1);
    if   (dsdx) Veclib::copy (np, dsdx+low, skip, nx, 1);
    else        Veclib::zero (np,                 nx, 1);
    if   (dsdy) Veclib::copy (np, dsdy+low, skip, ny, 1);
    else        Veclib::zero (np,                 ny, 1);
    
    break;

  case 3:
    skip = -np;
    xs   = work();
    ys   = xs + np;
      
    Blas::gemv    ("T", np, np, 1.0, *D, np, xmesh, np, 0.0, xs, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, ymesh, np, 0.0, ys, 1);
    Veclib::vmul  (np, xs, 1, xs, 1, area, 1);
    Veclib::vvtvp (np, ys, 1, ys, 1, area, 1, area, 1);
    if   (drdx) Veclib::smul (np, -1.0, drdx, skip, nx, 1);
    else        Veclib::zero (np,                   nx, 1);
    if   (drdy) Veclib::smul (np, -1.0, drdy, skip, ny, 1);
    else        Veclib::zero (np,                   ny, 1);
    
    break;
  }
  
  Veclib::vsqrt (np, area, 1, area, 1);
  Veclib::vmul  (np, area, 1, w,    1, area, 1);

  len = work();

  Veclib::vhypot (np, nx, 1, ny,  1, len, 1);
  Veclib::vdiv   (np, nx, 1, len, 1, nx,  1);
  Veclib::vdiv   (np, ny, 1, len, 1, ny,  1);
}


void Element::sideEval (const int   side,
			real*       tgt ,
			const char* func) const
// ---------------------------------------------------------------------------
// Evaluate function func along side of element, returning in tgt.
//
// The function can use variables "x", "y" & "t" (and any floating-point
// parameters previously set).
// ---------------------------------------------------------------------------
{
  register  int estart, skip, bstart;
  terminal (side, estart, skip, bstart);
  
  vector<real> work(np + np);
  real         *x, *y;

  x = work();
  y = x + np;

  Veclib::copy (np, xmesh + estart, skip, x, 1);
  Veclib::copy (np, ymesh + estart, skip, y, 1);

  Femlib::prepVec  ("x y", func);
  Femlib__parseVec (np, x, y, tgt);
}


void Element::sideGrad (const int   side,
			const real* src ,
			real*       c1  ,
			real*       c2  ) const
// ---------------------------------------------------------------------------
// Using geometric factors for this Element, return the first and second
// component, c1 and/or c2, of grad src (length nTot()) along side.
//
// We have to take some special care on  sides 2 & 3, where the usual skips
// are negative: we instead use positive skips for formation of dc/dr, dc/ds,
// then a -1 skip when multiplying by dr/dx, ds/dx, etc.
// ---------------------------------------------------------------------------
{
  register  int estart, skip, bstart, d;
  terminal (side, estart, skip, bstart);
  
  vector<real> work (np + np);
  real         *ddr, *dds, **DV, **DT;

  ddr = work();
  dds = ddr + np;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  // -- Make dc/dr, dc/ds along element edge.

  switch (side) {
  case 0:
    d = 1;
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("N", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, dds, 1);
    break;
  case 1:
    d = 1;
    Blas::gemv ("T", np, np, 1.0, src, np, *DT + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, dds, 1);
    break;
  case 2:
    d = -1;
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("N", np, np, 1.0, src, np, *DV + estart, d*skip, 0.0, dds, 1);
    break;
  case 3:
    d = -1;
    Blas::gemv ("T", np, np, 1.0, src, np, *DT + estart, d*skip, 0.0, ddr, 1);
    Blas::gemv ("T", np, np, 1.0, *DV, np, src + estart, d*skip, 0.0, dds, 1);
    break;
  }

  // -- dc/dx = dc/dr * dr/dx + dc/ds * ds/dx.

  if (c1) {
    if   (drdx) Veclib::vmul  (np, ddr, d, drdx + estart, skip, c1, 1);
    else        Veclib::zero  (np, c1, 1);
    if   (dsdx) Veclib::vvtvp (np, dds, d, dsdx + estart, skip, c1, 1, c1, 1);
  }
  
  // -- dc/dy = dc/dr * dr/dy + dc/ds * ds/dy.

  if (c2) {
    if   (drdy) Veclib::vmul  (np, ddr, d, drdy + estart, skip, c2, 1);
    else        Veclib::zero  (np, c2, 1);
    if   (dsdy) Veclib::vvtvp (np, dds, d, dsdy + estart, skip, c2, 1, c2, 1);
  }
}


void Element::sideSet (const int   side,
		       const int*  bmap,
                       const real* src ,
                       real*       tgt ) const
// ---------------------------------------------------------------------------
// Load edge vector src into globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  register int estart, skip, bstart;
  const int    nm = np - 1;

  terminal (side, estart, skip, bstart);

  Veclib::scatr (nm, src, bmap + bstart, tgt);
  if   (side == ns - 1) tgt[bmap[          0]] = src[nm];
  else                  tgt[bmap[bstart + nm]] = src[nm];

}


void Element::sideDsSum (const int   side,
			 const int*  bmap,
                         const real* src ,
                         const real* area,
                         real*       tgt ) const
// ---------------------------------------------------------------------------
// Direct-stiffness-sum (area-weighted) vector src into globally-numbered
// tgt on side.  This is for evaluation of natural BCs using Gauss--
// Lobatto quadrature.
// ---------------------------------------------------------------------------
{
  const int nm  = np - 1;
  vector<real> tmp (np);
  
  int estart, skip, bstart;
  terminal (side, estart, skip, bstart);

  Veclib::vmul (np, src, 1, area, 1, tmp(), 1);

  Veclib::scatr_sum (nm, tmp(), bmap + bstart, tgt);
  if   (side == ns - 1) tgt[bmap[          0]] += tmp[nm];
  else                  tgt[bmap[bstart + nm]] += tmp[nm];
}


void Element::evaluate (const char* func,
			real*       tgt ) const
// ---------------------------------------------------------------------------
// Evaluate function over mesh points, store in tgt.
// Function can explicitly use "x" and "y", for which mesh values are used.
// ---------------------------------------------------------------------------
{
  Femlib::prepVec  ("x y", func);
  Femlib__parseVec (nTot (), xmesh, ymesh, tgt);
}


real Element::integral (const char* func) const
// ---------------------------------------------------------------------------
// Return integral of func over element, using element quadrature rule.
// ---------------------------------------------------------------------------
{
  real         intgrl;
  const int    ntot = nq * nq;
  vector<real> tmp (ntot);

  Femlib::prepVec ("x y", func);

  switch (rule) {

  case LL:
    Femlib__parseVec (ntot, xmesh, ymesh, tmp());
    break;

  case GL:
    vector<real> work (ntot + ntot + nq * np);
    real         *x, *y, *wrk, **IN, **IT;

    x   = work();
    y   = x + ntot;
    wrk = y + ntot;

    Femlib::quad (GL, np, nq, 0, 0, 0, &IN, &IT, 0, 0);

    Blas::mxm (*IN, nq, xmesh, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,   np, x,   nq);

    Blas::mxm (*IN, nq, ymesh, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,   np, y,   nq);

    Femlib__parseVec (ntot, x, y, tmp());

    break;
  }

  Veclib::vmul (ntot, tmp(), 1, G4, 1, tmp(), 1);

  intgrl = Veclib::sum (ntot, tmp(), 1);
  
  return intgrl;
}


real Element::norm_inf (const real* src) const
// ---------------------------------------------------------------------------
// Return infinity-norm of element value.
// ---------------------------------------------------------------------------
{
  return fabs (src[Blas::iamax (nTot(), src, 1)]);
}


real Element::norm_L2 (const real* src) const
// ---------------------------------------------------------------------------
// Return L2-norm of Element value, using Element quadrature rule.
// ---------------------------------------------------------------------------
{
  register int   i;
  register real  L2   = 0.0;
  const int      ntot = nq * nq;
  register real* dA   = G4;

  switch (rule) {

  case LL:
    for (i = 0; i < ntot; i++) L2 += src[i] * src[i] * dA[i];
    break;

  case GL:
    vector<real> work (ntot + nq * np);
    real         *wrk, *u, **IN, **IT;
    
    u   = work();
    wrk = u + ntot;

    Femlib::quad (GL, np, nq, 0, 0, 0, &IN, &IT, 0, 0);

    Blas::mxm (*IN, nq,  src, np, wrk, np);
    Blas::mxm (wrk, nq, *IT,  np, u,   nq);

    for (i = 0; i < ntot; i++) L2 += u[i] * u[i] * dA[i];
    
    break;
  }

  return sqrt (L2);
}


real Element::norm_H1 (const real* src) const
// ---------------------------------------------------------------------------
// Return Sobolev-1 norm of Element value, using Element quadrature rule.
// ---------------------------------------------------------------------------
{
  register real H1 = 0;
  register int  i;
  const int     ntot = nq * nq;
  
  vector<real>  work (ntot);
  register real *u, *dA;
  
  u = work();

  switch (rule) {

  case LL:
    dA = G4;

    // -- Add in L2 norm of u.

    for (i = 0; i < ntot; i++) H1 += src[i] * src[i] * dA[i];

    // -- Add in L2 norm of grad u.

    Veclib::copy (ntot, src, 1, u, 1);
    grad (u, 0);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * dA[i];

    Veclib::copy (ntot, src, 1, u, 1);
    grad (0, u);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * dA[i];

    break;

  case GL:
    vector<real> work2 (nq * np);
    real         *w, **IN, **IT, **DV, **DT;

    w = work2();

    Femlib::quad (GL, np, nq, 0, 0, 0, &IN, &IT, &DV, &DT);

    Blas::mxm (*IN, nq,  src, np, w, np);
    Blas::mxm (w,   nq, *IT,  np, u, nq);

    // -- Add in L2 norm of u.

    dA = G4;
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * dA[i];

    // -- Add in L2 norm of grad u.

    dA = G1;
    Blas::mxm (*IN, nq,  src, np, w, np);
    Blas::mxm (w,   nq, *DV,  np, u, nq);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * dA[i];

    dA = G2;
    Blas::mxm (*DV, nq,  src, np, w, np);
    Blas::mxm (w,   nq, *IT,  np, u, nq);
    for (i = 0; i < ntot; i++) H1 += u[i] * u[i] * dA[i];

    if (G3) {
      dA = G3;
      Blas::mxm (*DV, nq,  src, np, w, np);
      Blas::mxm (w,   nq, *DT,  np, u, nq);
      for (i = 0; i < ntot; i++) H1 += 2.0 * u[i] * u[i] * dA[i];
    }

    break;
  }

  return sqrt (H1);
}


void Element::e2g (const real* src     ,
		   const int*  btog    ,
		   real*       external,
		   real*       internal) const
// ---------------------------------------------------------------------------
// Src is a row-major vector, representing this Element's data.
// Load boundary data into globally-numbered vector external, and
// internal data into un-numbered vector internal.
//
// REPLACES old routine Element::condense
// ---------------------------------------------------------------------------
{
  const int next = nExt();
  const int nint = nInt();

  Veclib::gathr_scatr (next, src, emap,  btog, external);
  if (internal)
    Veclib::gathr     (nint, src, emap + next, internal);
}


void Element::e2gSum (const real* src     ,
		      const int*  btog    ,
		      real*       external,
		      real*       internal) const
// ---------------------------------------------------------------------------
// Src is a row-major vector, representing this element's data.
// Sum boundary data into globally-numbered vector external, and
// internal data into un-numbered vector internal.
//
// REPLACES old routine Element::dsForcing
// ---------------------------------------------------------------------------
{
  const int next = nExt();
  const int nint = nInt();

  Veclib::gathr_scatr_sum (next, src, emap,  btog, external);
  if (internal)
    Veclib::gathr_sum     (nint, src, emap + next, internal);
}


void Element::g2e (real*       tgt     ,
		   const int*  btog    ,
		   const real* external,
		   const real* internal) const
// ---------------------------------------------------------------------------
// Tgt is a row-major vector, representing this Element's data.
// Load boundary data from globally-numbered vector external, and
// internal data from un-numbered vector internal.
//
// REPLACES old routine Element::expand
// ---------------------------------------------------------------------------
{
  const int next = nExt();
  const int nint = nInt();

  Veclib::gathr_scatr (next, external, btog,  emap, tgt);
  if (internal)
    Veclib::scatr     (nint, internal, emap + next, tgt);
}


void Element::HelmholtzOp (const real* src,
			   real*       tgt, 
			   const real  L2 ,
			   real*       wrk) const
// ---------------------------------------------------------------------------
// Apply elemental discrete Helmholtz operator on src to make tgt.
// This routine is specific to 2D Cartesian space with GLL integration.
//
// Input workspace vector wrk must hold 2 * nTot() elements.
// ---------------------------------------------------------------------------
{
  register int  ij;
  const int     ntot = nTot();
  register real tmp, *R, *S;
  real          **DV, **DT;

  R = wrk;
  S = R + ntot;

  Femlib::quad (rule, np, nq, 0, 0, 0, 0, 0, &DV, &DT);

  Blas::gemm ("N", "N", np, np, np, 1.0, *DT, np, src, np, 0.0, R, np);
  Blas::gemm ("N", "N", np, np, np, 1.0, src, np, *DV, np, 0.0, S, np);

  if (G3) {
    for (ij = 0; ij < ntot; ij++) {
      tmp      = R [ij];
      R  [ij]  = G1[ij] * R  [ij] + G3[ij] * S  [ij];
      S  [ij]  = G2[ij] * S  [ij] + G3[ij] * tmp;
      tgt[ij]  = G4[ij] * src[ij];
    }
  } else {
    for (ij = 0; ij < ntot; ij++) {
      R  [ij] *= G1[ij];
      S  [ij] *= G2[ij];
      tgt[ij]  = G4[ij] * src[ij];
    }
  }

  Blas::gemm ("N", "N", np, np, np, 1.0,  S,  np, *DT, np, L2,  tgt, np);
  Blas::gemm ("N", "N", np, np, np, 1.0, *DV, np,  R,  np, 1.0, tgt, np);
}
