///////////////////////////////////////////////////////////////////////////////
// element.C
///////////////////////////////////////////////////////////////////////////////

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
// Create a new quad element, nk X nk.  Node spacing along any side generated
// by mapping z (defined on domain [-1, 1], np points) onto side.
//
// Compute information for internal storage, and economize.
// ---------------------------------------------------------------------------
{
  char routine[] = "Element::Element";
  const int nk2  = sqr (nk);

  if (nk < 2) message (routine, "need > 2 knots for element edges", ERROR);

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
// For Lobatto--Legendre rule, everything is at the nodes, hence the
// interpolant matrices are identities.
// 
// The inverse partials are retained in Element storage, to be used in
// element gradient operations (e.g. dP/dx = dP/dr * dr/dx + dP/ds *
// ds/dx) while the forward partials are retained in scrambled form (in
// combination with quadrature weights) as "geometric factors" G1--G4, to
// be used in element quadrature operations.  For cylindrical geometries,
// G1--G4 are multiplied by y (i.e. r) as a consequence of the fact that
// Helmholtz equations are symmetrized by premultiplication by this
// factor for cylindrical coords.
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
  char       routine[] = "Element::map", err[StrMax];
  const int  ntot = nTot();
  const real EPS  = 4 * nTot()*((sizeof(real)==sizeof(double)) ? EPSDP:EPSSP);
  const real *x   = xmesh, *y = ymesh;
  const real **DV, **DT, *w;
  real       *jac, *dxdr, *dxds, *dydr, *dyds, *tV, *WW;

  vector<real> work;

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
    
  if (jac[Veclib::imin (ntot, jac, 1)] < EPS) {
    sprintf (err, "jacobian of element %1d nonpositive", id + 1);
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
  
  Veclib::vmul  (ntot, jac,  1, WW,   1, G4, 1);

  Veclib::copy (ntot, dyds, 1, drdx, 1);
  Veclib::vneg (ntot, dxds, 1, drdy, 1);
  Veclib::vneg (ntot, dydr, 1, dsdx, 1);
  Veclib::copy (ntot, dxdr, 1, dsdy, 1);
    
  Veclib::vdiv (ntot, drdx, 1, jac, 1, drdx, 1);
  Veclib::vdiv (ntot, drdy, 1, jac, 1, drdy, 1);
  Veclib::vdiv (ntot, dsdx, 1, jac, 1, dsdx, 1);
  Veclib::vdiv (ntot, dsdy, 1, jac, 1, dsdy, 1);

  if (Geometry::system() == Geometry::Cylindrical) {
    Veclib::vmul (ntot, G1, 1, y, 1, G1, 1);
    Veclib::vmul (ntot, G2, 1, y, 1, G2, 1);
    Veclib::vmul (ntot, G3, 1, y, 1, G3, 1);
    Veclib::vmul (ntot, G4, 1, y, 1, G4, 1);
  } 

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
  register real* wt = G4;

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
			   const real betak2 ,
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
// lambda2 is the Helmholtz constant.
//
// k2 is the square of the wavenumber for the Fourier decomposition that is
// used in the azimuthal direction in cylindrical coordinates, and effectively
// serves as a flag for use of cylindrical coordinates: it should always be
// zero for Cartesian coordinates.

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
  const real   **DV, **DT;

  // -- Construct hbb, hbi, hii partitions of elemental Helmholtz matrix.

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  for (i = 0; i < np; i++)
    for (j = 0; j < np; j++, ij++) {

      helmRow ((const real**) DV, (const real**) DT,
	       lambda2, betak2, i, j, rmat, rwrk);

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
			 const real betak2 ,
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
// ---------------------------------------------------------------------------
{
  register int ij   = 0;
  const int    ntot = nTot ();
  const real   **DV, **DT;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  for (register int i = 0; i < np; i++)
    for (register int j = 0; j < np; j++, ij++) {
      helmRow ((const real**) DV, (const real**) DT,
	       lambda2, betak2, i, j, rmat, rwrk);
      Veclib::copy (ntot, rmat, 1, h + ij * np, 1);
    }
}


void Element::HelmholtzDg (const real lambda2,
			   const real betak2 ,
			   real*      diag   ,
			   real*      work   ) const
// ---------------------------------------------------------------------------
// Create the diagonal of the elemental Helmholtz matrix in diag.  The
// diagonal is sorted in emap order: i.e., boundary nodes are first.
//
// Input vector diag must be nTot() long, work must be Ntot() + nKnot() long.
// Construction is very similar to that in helmRow except that m, n = i, j.
// ---------------------------------------------------------------------------
{
  register int  i, j;
  const int     ntot = sqr (np);
  const real    EPS  = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  const real**  DT;
  register real *dg = work, *tmp = work + ntot;
  real          r2, HCon;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, 0, &DT);

  if (Geometry::system() == Geometry::Cylindrical) {
    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++, dg++) {
	r2   = sqr (ymesh[Veclib::row_major (i, j, np)]);
	HCon = (r2 > EPS) ? (betak2 / r2 + lambda2) : 0.0;
	Veclib::vmul (np, DT[j], 1, DT[j], 1, tmp, 1);
	*dg  = Blas::dot   (np, G1 + i*np, 1, tmp, 1);
	Veclib::vmul (np, DT[i], 1, DT[i], 1, tmp, 1);
	*dg += Blas::dot   (np, G2 + j,   np, tmp, 1);
	if (G3)
	  *dg += 2.0 * G3[Veclib::row_major (i, j, np)] * DT[j][j] * DT[i][i];
	*dg += HCon  * G4[Veclib::row_major (i, j, np)];
      }
  } else {
    HCon = lambda2 + betak2;
    for (i = 0; i < np; i++)
      for (j = 0; j < np; j++, dg++) {
	Veclib::vmul (np, DT[j], 1, DT[j], 1, tmp, 1);
	*dg  = Blas::dot   (np, G1 + i*np, 1, tmp, 1);
	Veclib::vmul (np, DT[i], 1, DT[i], 1, tmp, 1);
	*dg += Blas::dot   (np, G2 + j,   np, tmp, 1);
	if (G3)
	  *dg += 2.0 * G3[Veclib::row_major (i, j, np)] * DT[j][j] * DT[i][i];
	*dg += HCon  * G4[Veclib::row_major (i, j, np)];
      }
  }

  Veclib::gathr (ntot, work, emap, diag);
}


void Element::helmRow (const real** DV     ,
		       const real** DT     ,
		       const real   lambda2,
		       const real   betak2 ,
		       const int    i      ,
		       const int    j      ,
		       real*        hij    ,
		       real*        work   ) const
// ---------------------------------------------------------------------------
// Build row [i,j] of the elemental Helmholtz matrix in array hij (np x np).
//
// Lambda2 is the Helmholtz constant.
//
// k2 is the square of the wavenumber for the Fourier decomposition that is
// used in the azimuthal direction in cylindrical coordinates, and effectively
// serves as a flag for use of cylindrical coordinates: it should always be
// zero for Cartesian coordinates.
//
// Input array work should be at least np long.
//
// For a 2D tensor product form, the elemental Helmholtz matrix is produced
// as (sums on p & q indices assumed):
//
// h      = G1  IN  DT  IN  DT     \
//  ij mn     pq  pi  jq  pm  nq   |
//        + G2  DV  IT  DV  IT     |
//            pq  pi  jq  pm  nq   |
//        + G3  DV  IT  IN  DT      >                              "STIFFNESS"
//            pq  pi  jq  pm  nq   |
//        + G3  IN  DT  DV  IT     |
//            pq  pi  jq  pm  nq   /
//                
//        + G4  IN  IT  IN  IT        (k2 / sqr (r  ) + lambda2)        "MASS"
//            pq  pi  jq  pm  nq                  pq
//
// where the terms G1, G2, G3, G4 contain geometric mapping factors and
// quadrature weights, and the matrices IN, IT are the Lagrangian
// interpolation matrix (from the nodes to the quadrature points) and its
// transpose, while DV, DT are the Lagrangian derivative matrix & transpose.
//
// For Gauss--Lobatto--Legendre integration, the interpolant matrices are
// indentities, and are not required.
// ---------------------------------------------------------------------------
{
  register int m, n;
  const int    ntot = sqr (np);
  const real   r2   = sqr (ymesh[Veclib::row_major (i, j, np)]);
  const real   EPS  = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  const real   hCon = (Geometry::system() == Geometry::Cylindrical &&
		       r2 > EPS) ? (betak2 / r2 + lambda2) : betak2 + lambda2;

  Veclib::zero (ntot, hij, 1);

  for (n = 0; n < np; n++) {
    Veclib::vmul (np, DT[j], 1, DT[n], 1, work, 1);
    hij[Veclib::row_major (i, n, np)]  = Blas::dot (np, G1 + i*np, 1, work, 1);
  }

  for (m = 0; m < np; m++) {
    Veclib::vmul (np, DT[i], 1, DT[m], 1, work, 1);
    hij[Veclib::row_major (m, j, np)] += Blas::dot (np, G2 + j,   np, work, 1);
  }

  if (G3)
    for (m = 0; m < np; m++)
      for (n = 0; n < np; n++) {
	hij [Veclib::row_major (m, n, np)] +=
	  G3[Veclib::row_major (i, n, np)] * DV[n][j] * DV[i][m];
	hij [Veclib::row_major (m, n, np)] +=
	  G3[Veclib::row_major (m, j, np)] * DV[j][n] * DV[m][i];
      }

  hij[Veclib::row_major (i, j, np)] += G4[Veclib::row_major (i, j, np)] * hCon;
}


void Element::HelmholtzOp (const real  lambda2,
			   const real  betak2 ,
			   const real* src    ,
			   real*       tgt    , 
			   real*       wrk    ) const
// ---------------------------------------------------------------------------
// Apply elemental discrete Helmholtz operator on src to make tgt.
//
// Lambda2 is the Helmholtz constant.
//
// k2 is the square of the wavenumber for the Fourier decomposition that is
// used in the azimuthal direction in cylindrical coordinates, and effectively
// serves as a flag for use of cylindrical coordinates: it should always be
// zero for Cartesian coordinates.
//
// Input workspace vector wrk must hold 2 * nTot() elements.
// ---------------------------------------------------------------------------
{
  register int  ij;
  const int     ntot = nTot();
  register real tmp, r2, hCon;
  register real *R, *S, *g1, *g2, *g3, *g4;
  const real    **DV, **DT;
  const real    EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  R  = wrk;
  S  = R + ntot;
  g1 = G1;
  g2 = G2;
  g3 = G3;
  g4 = G4;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  Blas::gemm ("N", "N", np, np, np, 1.0, *DT, np, src, np, 0.0, R, np);
  Blas::gemm ("N", "N", np, np, np, 1.0, src, np, *DV, np, 0.0, S, np);

  if (Geometry::system() == Geometry::Cylindrical) {
    if (g3) {
      for (ij = 0; ij < ntot; ij++) {
	r2       = sqr (ymesh[ij]);
	hCon     = (r2 > EPS) ? (betak2 / r2 + lambda2) : 0.0;
	tmp      = R [ij];
	R  [ij]  = g1[ij] * R  [ij] + g3[ij] * S  [ij];
	S  [ij]  = g2[ij] * S  [ij] + g3[ij] * tmp;
	tgt[ij]  = g4[ij] * src[ij] * hCon;
      }
    } else {
      for (ij = 0; ij < ntot; ij++) {
	r2       = sqr (ymesh[ij]);
	hCon     = (r2 > EPS) ? (betak2 / r2 + lambda2) : 0.0;
	R  [ij] *= g1[ij];
	S  [ij] *= g2[ij];
	tgt[ij]  = g4[ij] * src[ij] * hCon;
      }
    }

  } else {			// -- Cartesian.
    hCon = betak2 + lambda2;
#ifdef __uxp__
    if (g3) {
      Veclib::copy    (ntot,       R,  1, tgt, 1);
      Veclib::vvtvvtp (ntot,       g1, 1, R,   1, g3,  1, S,   1, R, 1);
      Veclib::vvtvvtp (ntot,       g2, 1, S,   1, g3,  1, tgt, 1, S, 1);
      Veclib::svvtt   (ntot, hCon, g4, 1, src, 1, tgt, 1);
    } else {
      Veclib::vmul    (ntot,       R,  1, g1,  1, R,   1);
      Veclib::vmul    (ntot,       S,  1, g2,  1, S,   1);
      Veclib::svvtt   (ntot, hCon, g4, 1, src, 1, tgt, 1);
    }
#else
    if (g3) {
      for (ij = 0; ij < ntot; ij++) {
	tmp      = R [ij];
	R  [ij]  = g1[ij] * R  [ij] + g3[ij] * S  [ij];
	S  [ij]  = g2[ij] * S  [ij] + g3[ij] * tmp;
	tgt[ij]  = g4[ij] * src[ij] * hCon;
      }
    } else {
      for (ij = 0; ij < ntot; ij++) {
	R  [ij] *= g1[ij];
	S  [ij] *= g2[ij];
	tgt[ij]  = g4[ij] * src[ij] * hCon;
      }
    }
#endif
  }

  Blas::gemm ("N", "N", np, np, np, 1.0,  S,  np, *DT, np, 1.0, tgt, np);
  Blas::gemm ("N", "N", np, np, np, 1.0, *DV, np,  R,  np, 1.0, tgt, np);
}


void Element::grad (real* tgtA,
		    real* tgtB) const
// ---------------------------------------------------------------------------
// Operate partial derivative d(tgt)/dxi = d_dr*drdxi + d_ds*dsdxi,
// where the appropriate component of gradient is selected by input pointers.
// Values are computed at node points.
// ---------------------------------------------------------------------------
{
  const real **DV, **DT;
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
//
// For cylindrical coordinates the area variable is weighted by y (i.e. r).
// ---------------------------------------------------------------------------
{
  if (side < 0 || side >= ns)
    message ("Element::sideGeom", "illegal side", ERROR);

  register int low, skip;
  const real   **D, *w;
  real         *xr, *xs, *yr, *ys, *len;
  vector<real> work (np + np);

  Femlib::quad (LL, np, np, 0, 0, &w, 0, 0, &D, 0);

  switch (side) {
  case 0: 
    low  = 0;
    skip = 1;
    xr   = work();
    yr   = xr + np;
    
    Blas::gemv    ("T", np, np, 1.0, *D, np, xmesh+low, 1, 0.0, xr, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, ymesh+low, 1, 0.0, yr, 1);
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
    low  = 0;
    skip = -np;
    xs   = work();
    ys   = xs + np;
      
    Blas::gemv    ("T", np, np, 1.0, *D, np, xmesh+low, np, 0.0, xs, 1);
    Blas::gemv    ("T", np, np, 1.0, *D, np, ymesh+low, np, 0.0, ys, 1);
    Veclib::vmul  (np, xs, 1, xs, 1, area, 1);
    Veclib::vvtvp (np, ys, 1, ys, 1, area, 1, area, 1);

    if   (drdx) Veclib::smul (np, -1.0, drdx, skip, nx, 1);
    else        Veclib::zero (np,                   nx, 1);
    if   (drdy) Veclib::smul (np, -1.0, drdy, skip, ny, 1);
    else        Veclib::zero (np,                   ny, 1);

    break;
  }
  
  Veclib::vsqrt  (np, area, 1, area, 1);
  Veclib::vmul   (np, area, 1, w,    1, area, 1);
  if (Geometry::system() == Geometry::Cylindrical)
    Veclib::vmul (np, area, 1, ymesh+low, skip, area, 1);

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
  register int d, estart, skip, bstart;
  terminal (side, estart, skip, bstart);
  
  vector<real> work (np + np);
  const real   **DV, **DT;
  real         *ddr, *dds;

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
  const int    ntot = np * np;
  vector<real> tmp (ntot);

  Femlib::prepVec  ("x y", func);
  Femlib__parseVec (ntot, xmesh, ymesh, tmp());

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
  const int      ntot = sqr (np);
  register real* dA   = G4;

  for (i = 0; i < ntot; i++) L2 += src[i] * src[i] * dA[i];

  return sqrt (L2);
}


real Element::norm_H1 (const real* src) const
// ---------------------------------------------------------------------------
// Return Sobolev-1 norm of Element value, using Element quadrature rule.
// ---------------------------------------------------------------------------
{
  register real H1 = 0;
  register int  i;
  const int     ntot = sqr (np);
  
  vector<real>  work (ntot);
  register real *u, *dA;
  
  u = work();

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
// ---------------------------------------------------------------------------
{
  const int next = nExt();
  const int nint = nInt();

  Veclib::gathr_scatr (next, external, btog,  emap, tgt);
  if (internal)
    Veclib::scatr     (nint, internal, emap + next, tgt);
}


void Element::divR (real* src) const
// ---------------------------------------------------------------------------
// Divide src by y (i.e. r in cylindrical coordinates), take special action
// where r = 0.  This is used in taking theta component of gradient.
// ---------------------------------------------------------------------------
{
  register int   i;
  register real  rad, rinv;
  register real* y   = ymesh;
  const int      N   = nTot();
  const real     EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  for (i = 0; i < N; i++) {
    rad     = y[i];
    rinv    = (rad > EPS) ? 1.0 / rad : 0.0;
    src[i] *= rinv;
  }
}


void Element::mulR (real* src) const
// ---------------------------------------------------------------------------
// Multiply src by y (i.e. r in cylindrical coordinates).
// ---------------------------------------------------------------------------
{
  Veclib::vmul (nTot(), src, 1, ymesh, 1, src, 1);
}


void Element::sideDivR (const int   side,
			const real* src ,
			real*       tgt ) const
// ---------------------------------------------------------------------------
// Deliver in tgt the side traverse of src (elemental storage) divided by
// y (i.e. r), take special action where r = 0.
// ---------------------------------------------------------------------------
{
         int  i, base, skip;
         real r, rinv, *y;
   const real *s;
  const          real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  switch (side) {
  case 0: 
    base = 0;
    skip = 1;
    y    = ymesh;
    s    = src;
    break;
  case 1:
    base = np - 1;
    skip = np;
    y    = ymesh + base;
    s    = src   + base;
    break;
  case 2:
    base = np * np - 1;
    skip = -1;
    y    = ymesh + base;
    s    = src   + base;
    break;
  case 3:
    base = np * (np - 1);
    skip = -np;
    y    = ymesh + base;
    s    = src   + base;
    break;
  }

  for (i = 0; i < np; i++) {
    r      = y[i*skip];
    rinv   = (r > EPS) ? 1.0 / r : 0.0;
    tgt[i] = rinv * s[i*skip];
  }
}


void Element::sideDivR2 (const int   side,
			 const real* src ,
			 real*       tgt ) const
// ---------------------------------------------------------------------------
// Deliver in tgt the side traverse of src (elemental storage) divided by
// y^2 (i.e. r^2), take special action where r = 0.
// ---------------------------------------------------------------------------
{
  register       int  i, base, skip;
  register       real r, rinv2, *y;
  register const real *s;
  const          real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  switch (side) {
  case 0: 
    base = 0;
    skip = 1;
    y    = ymesh;
    s    = src;
    break;
  case 1:
    base = np - 1;
    skip = np;
    y    = ymesh + base;
    s    = src   + base;
    break;
  case 2:
    base = np * np - 1;
    skip = -1;
    y    = ymesh + base;
    s    = src   + base;
    break;
  case 3:
    base = np * (np - 1);
    skip = -np;
    y    = ymesh + base;
    s    = src   + base;
    break;
  }

  for (i = 0; i < np; i++) {
    r      = y[i*skip];
    rinv2  = (r > EPS) ? 1.0 / sqr(r) : 0.0;
    tgt[i] = rinv2 * s[i*skip];
  }
}


int Element::locate (const real x    ,
		     const real y    ,
		     real&      r    ,
		     real&      s     ,
		     const int  guess) const
// ---------------------------------------------------------------------------
// If x & y fall in this element, compute the corresponding r & s values, 
// and return 1.  Otherwise return 0.
//
// If guess = 0 (the default argument), the input value of (r, s) is used
// as an initial guess for N--R iteration.  Otherwise the  (r, s) value that
// corresponds to the closest point in the Element mesh to (x, y) is used.
// ---------------------------------------------------------------------------
{
  const int    MaxItn = 16;
  const real   EPS    = 50*((sizeof (double)==sizeof (real)) ? EPSDP : EPSSP);
  const real   DIVERG = 1.5;
  real         *J, *F, *ir, *is, *dr, *ds, *tp;
  vector<real> work (5 * np + 6);
  int          ipiv[2], info, i, j;
  
  tp = work();
  ir = tp + np;
  is = ir + np;
  dr = is + np;
  ds = dr + np;
  J  = ds + np;
  F  = J  + 4;
  
  if (guess) {
    const int    ntot = nTot();
    vector<real> tmp (2 * ntot);
    const real*  knot;
    real         *tx = tmp(), *ty = tmp() + ntot;

    Femlib::quad     (LL, np, np, &knot, 0, 0, 0, 0, 0, 0);
    Veclib::ssub     (ntot, x, xmesh, 1, tx, 1);
    Veclib::ssub     (ntot, y, ymesh, 1, ty, 1);
    Veclib::vvtvvtp  (ntot, tx, 1, tx, 1, ty, 1, ty, 1, tx, 1);
    
    i = Veclib::imin (ntot, tx, 1);
    j = i % np;
    i = (i - j) / np;

    r = knot[i];
    s = knot[j];
  }

  i = 0;
  do {
    Femlib::interp (LL, np, r, s, ir, is, dr, ds);

               Blas::gemv ("T", np, np, 1.0, xmesh, np, ir, 1, 0.0, tp, 1);
    F[0] = x - Blas::dot  (np, is, 1, tp, 1);
    J[2] =     Blas::dot  (np, ds, 1, tp, 1);
               Blas::gemv ("T", np, np, 1.0, ymesh, np, ir, 1, 0.0, tp, 1);
    F[1] = y - Blas::dot  (np, is, 1, tp, 1);
    J[3] =     Blas::dot  (np, ds, 1, tp, 1);
               Blas::gemv ("T", np, np, 1.0, xmesh, np, dr, 1, 0.0, tp, 1);
    J[0] =     Blas::dot  (np, is, 1, tp, 1);
               Blas::gemv ("T", np, np, 1.0, ymesh, np, dr, 1, 0.0, tp, 1);
    J[1] =     Blas::dot  (np, is, 1, tp, 1);
    
    Lapack::gesv (2, 1, J, 2, ipiv, F, 2, info);
    
    r += F[0];
    s += F[1];

    if (fabs (r) > DIVERG || fabs (s) > DIVERG) return 0;

  } while (++i < MaxItn && (fabs (F[0]) > EPS || fabs (F[1]) > EPS));

  return (i < MaxItn) ? 1 : 0;
}


real Element::probe (const real  r  ,
		     const real  s  ,
		     const real* src) const
// ---------------------------------------------------------------------------
// Return the value of field storage located at r, s, in this element.
// ---------------------------------------------------------------------------
{
  real         *ir, *is, *tp;
  vector<real> work (3 * np);

  ir = work();
  is = ir + np;
  tp = is + np;

  Femlib::interp   (LL, np, r, s, ir, is, 0, 0);
  Blas::gemv       ("T", np, np, 1.0, src, np, ir, 1, 0.0, tp, 1);

  return Blas::dot (np, is, 1, tp, 1);
}
