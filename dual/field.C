///////////////////////////////////////////////////////////////////////////////
// field.C: derived from AuxField, Field adds boundary conditions,
// global numbering, and the ability to solve Helmholtz problems.
//
// Copyright (C) 1994 <--> $Date$, Hugh Blackburn
//
// NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBN
//
// This version allows the user to freeze either the zeroth or first 
// Fourier mode, according to token
//
// FREEZE = -1 (default) update both modes as per standard dual.
// FREEZE = 0 Freeze velocity in mode 0 (and ignore its pressure field).
// FREEZE = 1 Freeze velocity in mode 1 (and ignore its pressure field).
// 
// NBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBNBN
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


Field::Field (BoundarySys*      B,
	      real_t*           M,
	      const int_t       N,
	      vector<Element*>& E,
	      const char        C) :
// ---------------------------------------------------------------------------
// Create storage for a new Field from scratch.
// ---------------------------------------------------------------------------
  AuxField (M, N, E, C),
  _bsys    (B)
{
  const int_t              np  = Geometry::nP();
  const int_t              nzb = Geometry::basePlane();
  const vector<Boundary*>& BC  = _bsys -> BCs (0);
  register real_t*         p;
  register int_t           i, k;

  // -- Allocate storage for boundary data.
  
  _nbound = _bsys -> nSurf();
  _nline  = _nbound * np;

  _line  = new real_t* [static_cast<size_t> (_nz)];
  _sheet = new real_t  [static_cast<size_t> (_nz * _nline)];

  for (k = 0; k < _nz; k++) _line[k] = _sheet + k*_nline;

  Veclib::zero (_nz * _nline, _sheet, 1);

  // -- Set values for boundary data, 0th Fourier mode only, enforce z = 0.

  Femlib::value ("z", 0);
  for (p = _line[0], i = 0; i < _nbound; i++, p += np)
    BC[i] -> evaluate (0, 0, p);
}


void Field::printBoundaries (const Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// (Debugging) Utility to print information contained in a Boundary list.
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = F -> _bsys -> BCs (0);
    int_t                    i;
  
    cout << "# -- Field '" << F -> name() << "' Boundary Information:" << endl;
    if (!F -> _nbound) cout << "No BCs for this Field" << endl;
    for (i = 0; i < F -> _nbound; i++) BC[i] -> print();
  }
}


void Field::evaluateBoundaries (const int_t step   ,
				const bool  notused)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
//
// This routine is only meant to be used for boundary conditions that must
// be re-evaluated at every step, such as high-order pressure BCs.
//
// NB: modified for dual: argument "unused" retained for compatibility
// with standard semtex class prototype.
// ---------------------------------------------------------------------------
{
  const int_t    np    = Geometry::nP();
  const int_t    kfund = Femlib::ivalue ("BETA");
  register int_t i, k;
  real_t*        p;

  const vector<Boundary*>& BC0 = _bsys -> BCs (0);
  const vector<Boundary*>& BCk = _bsys -> BCs (kfund);

  for (p = _line[0], i = 0; i < _nbound; i++, p += np)
    BC0[i] -> evaluate (0, step, p);

  for (k = 1; k < 3; k++) {
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BCk[i] -> evaluate (k, step, p);
  }
}


void Field::evaluateM0Boundaries (const int_t step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind, but only for Mode 0.
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> BCs (0);
    const int_t              np = Geometry::nP();
    real_t*                  p;
    register int_t           i;

    for (p = _line[0], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (0, step, p);
  }
}


void Field::addToM0Boundaries (const real_t  val,
			       const char* grp)
// ---------------------------------------------------------------------------
// Add val to zeroth Fourier mode's bc storage area on BC group "grp".
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> BCs (0);
    const int_t              np = Geometry::nP();
    real_t*                  p;
    register int_t           i;

    for (p = _line[0], i = 0; i < _nbound; i++, p += np)
      BC[i] -> addForGroup (grp, val, p);
  }
}


Field& Field::smooth (AuxField* slave)
// ---------------------------------------------------------------------------
// Smooth slave field along element boundaries using *this, with
// mass-average smoothing.  The operation is equivalent to finding
//
//            -1
//   {u} = [M]   Sum [M] {u}  ,
//      g     g     e   e   e
//
// where g ==> global, e ==> elemental, [M] ==> mass matrix, and the
// summation is a "direct stiffness summation", or matrix assembly.
//
// If slave == 0, smooth this -> data.
// ---------------------------------------------------------------------------
{
  const int_t      nel     = Geometry::nElmt();
  const int_t      npnp    = Geometry::nTotElmt();
  const int_t      next    = Geometry::nExtElmt();
  const NumberSys* N       = _bsys -> Nsys  (0);
  const real_t*    imass   = _bsys -> Imass (0);
  const int_t      nglobal = N    -> nGlobal();
  const int_t*     btog    = N    -> btog();
  const int_t*     gid;
  register int_t   i, k;
  vector<real_t>   work (nglobal);
  real_t           *src, *dssum = &work[0];

  for (k = 0; k < _nz; k++) {

    Veclib::zero (nglobal, dssum, 1);
    src = (slave) ? slave -> _plane[k] : _plane[k];
    gid = btog;

    for (i = 0; i < nel; i++, src += npnp, gid += next)
      _elmt[i] -> bndryDsSum (gid, src, dssum);

    Veclib::vmul (nglobal, dssum, 1, imass, 1, dssum, 1);
    src = (slave) ? slave -> _plane[k] : _plane[k];
    gid = btog;

    for (i = 0; i < nel; i++, src += npnp, gid += next)
      _elmt[i] -> bndryInsert (gid, dssum, src);
  }

  return *this;
}


void Field::smooth (const int_t nZ ,
		    real_t*     tgt) const
// ---------------------------------------------------------------------------
// Smooth tgt field along element boundaries using *this, with
// mass-average smoothing.  Tgt is assumed to be arranged by planes, with
// planeSize() offset between each plane of data.
// ---------------------------------------------------------------------------
{
  const int_t      nel     = Geometry::nElmt();
  const int_t      npnp    = Geometry::nTotElmt();
  const int_t      next    = Geometry::nExtElmt();
  const int_t      nP      = Geometry::planeSize();
  const NumberSys* N       = _bsys -> Nsys  (0);
  const real_t*    imass   = _bsys -> Imass (0);
  const int_t      nglobal = N    -> nGlobal();
  const int_t*     btog    = N    -> btog();
  const int_t*     gid;
  register int_t   i, k;
  vector<real_t>   work (nglobal);
  real_t           *src, *dssum = &work[0];

  for (k = 0; k < nZ; k++) {

    Veclib::zero (nglobal, dssum, 1);
    src = tgt + k * nP;
    gid = btog;

    for (i = 0; i < nel; i++, src += npnp, gid += next)
      _elmt[i] -> bndryDsSum (gid, src, dssum);

    Veclib::vmul (nglobal, dssum, 1, imass, 1, dssum, 1);
    src = tgt + k * nP;
    gid = btog;

    for (i = 0; i < nel; i++, src += npnp, gid += next)
      _elmt[i] -> bndryInsert (gid, dssum, src);
  }
}


real_t Field::scalarFlux (const Field* C)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute normal flux of field C on all "wall" group boundaries.
//
// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = C -> _bsys -> BCs (0);
  vector<real_t>           work(4 * Geometry::nP());
  real_t                   F = 0.0;
  register int_t           i;
  
  for (i = 0; i < C -> _nbound; i++)
    F += BC[i] -> scalarFlux ("wall", C -> _data, &work[0]);

  return F;
}


Vector Field::normTraction (const Field* P)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute normal tractive forces on all WALL boundaries, taking P to be
// the pressure field.
//
// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = P -> _bsys -> BCs (0);
  const int_t              nsurf = P -> _nbound;
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real_t>           work(Geometry::nP());
  register int_t           i;
  
  for (i = 0; i < nsurf; i++) {
    secF = BC[i] -> normTraction ("wall", P -> _data, &work[0]);
    F.x += secF.x;
    F.y += secF.y;
  }

  return F;
}


Vector Field::tangTraction (const Field* U,
			    const Field* V,
			    const Field* W)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute (2D) tangential viscous tractive forces on all WALL boundaries,
// treating U & V as first and second velocity components, respectively.
//
// Compute viscous tractive forces on wall from
//
//  t_i  = - T_ij * n_j       (minus sign for force exerted BY fluid ON wall),
//
// where
//
//  T_ij = viscous stress tensor (here in Cartesian coords)
//                          dU_i    dU_j
//       = RHO * KINVIS * ( ----  + ---- ) .
//                          dx_j    dx_i
//
// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& UBC =       U->_bsys->BCs(0);
  const vector<Boundary*>& WBC = (W) ? W->_bsys->BCs(0) : (vector<Boundary*>)0;
  const int_t              np      = Geometry::nP();
  const int_t              _nbound = U -> _nbound;
  const real_t             mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real_t>           work(4 * np);
  register int_t           i;

  for (i = 0; i < _nbound; i++) {
    secF = UBC[i] -> tangTraction  ("wall", U->_data, V->_data, &work[0]);
    F.x        -= mu * secF.x;
    F.y        -= mu * secF.y;
    if (W) F.z -= mu * WBC[i] -> scalarFlux ("wall", W->_data, &work[0]);
  }

  return F;
}


Field& Field::solve (AuxField*             f  ,
		     const ModalMatrixSys* MMS)
// ---------------------------------------------------------------------------
// Problem for solution is
//                                          
//                      div grad u - lambda^2 u = f,
//
// which is set up in discrete form as
//
//                       H v = - M f - H g + <h, w>.
//
// This routine creates the RHS vector from the input forcing field f
// and the Field's boundary conditions g (essential) & h (natural).
// Forcing field f's data area is overwritten/destroyed during
// processing.
//
// For DIRECT (Cholesky) solution:
//
//   The RHS vector is constructed with length of the number of
//   element-edge nodes in the problem (n_gid).  The first n_solve
//   values contain forcing terms for the free (non essential-BC)
//   nodes in the problem, derived from the forcing field and the
//   natural BCs "h", while the remaining values get loaded from
//   essential BC values, "g".
//
// For JACPCG (iterative) solution:
//
//   All vectors are ordered with globally-numbered (element-
//   boundary) nodes first, followed by all element-internal nodes.
//   The zeroing operation which occurs after each application of the
//   Helmholtz operator serves to apply the essential BCs, which are
//   zero during the iteration (see file header).
//
//   The notation follows that used in Fig 2.5 of Barrett et al.,
//   "Templates for the Solution of Linear Systems", netlib.
//
// NB: modified for dual.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Field::solve";
  const int_t np    = Geometry::nP();
  const int_t nel   = Geometry::nElmt();
  const int_t next  = Geometry::nExtElmt();
  const int_t npnp  = Geometry::nTotElmt();
  const int_t ntot  = Geometry::nPlane();
  int_t       i, k, mode;

  for (k = 0; k < _nz; k++) {	// -- Loop over planes of data.
    
    // -- Select Fourier mode, set local pointers and variables.

    if   (k == 0) mode = 0;
    else mode = (Geometry::cylindrical()) ? Femlib::ivalue ("BETA") : 1;

    if (Femlib::ivalue ("FREEZE") == 0 && mode == 0) return *this;
    if (Femlib::ivalue ("FREEZE") == 1 && mode >  0) return *this;

    const MatrixSys*         M       = (*MMS)[(k != 0)];
    const vector<Boundary*>& B       = M -> _BC;
    const NumberSys*         N       = M -> _NS;
    real_t                   lambda2 = M -> _HelmholtzConstant;
    real_t                   betak2  = M -> _FourierConstant;
    int_t                    nsolve  = M -> _nsolve;
    int_t                    nglobal = M -> _nglobal;
    int_t                    nzero   = nglobal - nsolve;
    real_t*                  forcing = f -> _plane[k];
    real_t*                  unknown = _plane     [k];
    real_t*                  bc      = _line      [k];

    switch (M -> _method) {

    case DIRECT: {
      const real_t*  H       = const_cast<const real_t*>  (M -> _H);
      const real_t** hii     = const_cast<const real_t**> (M -> _hii);
      const real_t** hbi     = const_cast<const real_t**> (M -> _hbi);
      const int_t*   b2g     = const_cast<const int_t*>   (N -> btog());
      int_t          nband   = M -> _nband;

      vector<real_t> work (nglobal + 4*npnp);
      real_t         *RHS = &work[0], *tmp = RHS + nglobal;
      int_t          info;
      
      // -- Build RHS = - M f - H g + <h, w>.

      Veclib::zero (nglobal, RHS, 1);
      
      this -> getEssential (bc, RHS, B, N);
      this -> constrain    (forcing, lambda2, betak2, RHS, N, tmp);
      this -> buildRHS     (forcing, bc, RHS, 0, hbi, nsolve, nzero, B, N,tmp);
      
      // -- Solve for unknown global-node values (if any).
      
      if (nsolve) Lapack::pbtrs("U",nsolve,nband-1,1,H,nband,RHS,nglobal,info);
      
      // -- Carry out Schur-complement solution for element-internal nodes.
      
      for (i = 0; i < nel; i++, b2g += next, forcing += npnp, unknown += npnp)
	_elmt[i] -> global2localSC (RHS,b2g,forcing,unknown,hbi[i],hii[i],tmp);

      unknown -= ntot;
    
      // -- Scatter-gather essential BC values into plane.

      Veclib::zero (nglobal, RHS, 1);
    
      this -> getEssential (bc, RHS, B,   N);
      this -> setEssential (RHS, unknown, N);
    }
    break;

    case JACPCG: {
      const int_t    StepMax = Femlib::ivalue ("STEP_MAX");  
      const int_t    npts    = M -> _npts;
      real_t         alpha, beta, dotp, epsb2, r2, rho1, rho2;
      vector<real_t> work (5 * npts + 3 * Geometry::nPlane());

      real_t* r   = &work[0];
      real_t* p   = r + npts;
      real_t* q   = p + npts;
      real_t* x   = q + npts;
      real_t* z   = x + npts;
      real_t* wrk = z + npts;

      Veclib::zero (nglobal, x, 1);

      this -> getEssential (bc,x,B,N);  
      this -> constrain    (forcing,lambda2,betak2,x,N,wrk);
      this -> buildRHS     (forcing,bc,r,r+nglobal,0,nsolve,nzero,B,N,wrk);

      epsb2  = Femlib::value ("TOL_REL") * sqrt (Blas::dot (npts, r, 1, r, 1));
      epsb2 *= epsb2;

      // -- Build globally-numbered x from element store.

      this -> local2global (unknown, x, N);
  
      // -- Compute first residual using initial guess: r = b - Ax.
      //    And mask to get residual for the zero-BC problem.

      Veclib::zero (nzero, x + nsolve, 1);   
      Veclib::copy (npts,  x, 1, q, 1);

      this -> HelmholtzOperator (q, p, lambda2, betak2, mode, wrk);

      Veclib::zero (nzero, p + nsolve, 1);
      Veclib::zero (nzero, r + nsolve, 1);
      Veclib::vsub (npts, r, 1, p, 1, r, 1);

      r2 = Blas::dot (npts, r, 1, r, 1);

      // -- PCG iteration.

      i = 0;
      while (r2 > epsb2 && ++i < StepMax) {

	// -- Preconditioner.

	Veclib::vmul (npts, M -> _PC, 1, r, 1, z, 1);

	rho1 = Blas::dot (npts, r, 1, z, 1);

	// -- Update search direction.

	if (i == 1)
	  Veclib::copy  (npts,             z, 1, p, 1); // -- p = z.
	else {
	  beta = rho1 / rho2;	
	  Veclib::svtvp (npts, beta, p, 1, z, 1, p, 1); // -- p = z + beta p.
	}

	// -- Matrix-vector product.

	this -> HelmholtzOperator (p, q, lambda2, betak2, mode, wrk);
	Veclib::zero (nzero, q + nsolve, 1);

	// -- Move in conjugate direction.

	dotp  = Blas::dot (npts, p, 1, q, 1);
	alpha = rho1 / dotp;
	Blas::axpy (npts,  alpha, p, 1, x, 1); // -- x += alpha p.
	Blas::axpy (npts, -alpha, q, 1, r, 1); // -- r -= alpha q.

	rho2 = rho1;
	r2   = Blas::dot (npts, r, 1, r, 1);
      }
  
      if (i == StepMax) message (routine, "step limit exceeded", WARNING);
  
      // -- Unpack converged vector x, impose current essential BCs.

      this -> global2local (x, unknown, N);
      
      this -> getEssential (bc, x, B,   N);
      this -> setEssential (x, unknown, N);
  
      if (Femlib::ivalue ("VERBOSE") > 1) {
	char s[StrMax];
	sprintf (s, ":%3d iterations, field '%c'", i, _name);
	message (routine, s, REMARK);
      }
    }
    break;
    }
  }
  return *this;
}


void Field::constrain (real_t*          force  ,
		       const real_t     lambda2,
 		       const real_t     betak2 ,
		       const real_t*    esstlbc,
		       const NumberSys* N      ,
		       real_t*          work   ) const
// ---------------------------------------------------------------------------
// Replace f's data with constrained weak form of forcing: - M f - H g.
// On input, essential BC values (g) have been loaded into globally-numbered
// esstlbc, other values are zero.
// ---------------------------------------------------------------------------
{
  const int_t       np    = Geometry::nP();
  const int_t       nel   = Geometry::nElmt();
  const int_t       next  = Geometry::nExtElmt();
  const int_t       npnp  = Geometry::nTotElmt();
  const int_t       ntot  = Geometry::nPlane();
  const int_t*      emask = N -> emask();
  const int_t*      btog  = N -> btog();
  register Element* E;
  register int_t    i;
  real_t           *u = work, *tmp = work + npnp;

  // -- Manufacture -(M f + H g).

  for (i = 0; i < nel; i++, emask++, btog += next, force += npnp) {
    E = _elmt[i];
    E -> weight (force);	// -- f <-- M f.
    if (*emask) {		// -- f <-- M f + H g.
      Veclib::zero      (npnp, u, 1);
      E -> global2local (u, btog, esstlbc, 0);
      E -> HelmholtzOp  (lambda2, betak2, u, u, tmp);
      Veclib::vadd      (npnp, force, 1, u, 1, force, 1);
    }
  }

  force -= ntot;
  Veclib::neg (ntot, force, 1);
}


void Field::HelmholtzOperator (const real_t* x      ,
			       real_t*       y      ,
			       const real_t  lambda2,
			       const real_t  betak2 ,
			       const int_t   mode   ,
			       real_t*       work   ) const
// ---------------------------------------------------------------------------
// Discrete 2D global Helmholtz operator which takes the vector x into
// vector y, including direct stiffness summation.  Vectors x & y have 
// global ordering: that is, with nglobal (element edge nodes, with
// redundancy removed) coming first, followed by nel blocks of element-
// internal nodes.
//
// Vector work must have length 3 * Geometry::nPlane().
// ---------------------------------------------------------------------------
{
  const int_t      np      = Geometry::nP();
  const int_t      nel     = Geometry::nElmt();
  const int_t      npnp    = Geometry::nTotElmt();
  const int_t      ntot    = Geometry::nPlane();
  const NumberSys* NS      = _bsys -> Nsys (mode);
  const int_t      nglobal = NS -> nGlobal() + Geometry::nInode();
  register int_t   i;
  const real_t     *DV, *DT;
  real_t           *P = work, *R = P + ntot, *S = R + ntot;

  Veclib::zero (nglobal,     y, 1);
  Veclib::zero (ntot + ntot, R, 1);

  // -- Add in contributions from mixed BCs while x & y are global vectors.

  if (_bsys -> mixBC()) {
    const vector<Boundary*>& BC   = _bsys -> BCs (0);
    const int_t*             bmap = NS    -> btog();

    for (i = 0; i < _nbound; i++) BC[i] -> augmentOp (bmap, x, y);
  }

  // -- Add in contributions from elemental Helmholtz operations.

  Femlib::quadrature (0, 0, &DV, 0  , np, GLJ, 0.0, 0.0);
  Femlib::quadrature (0, 0, 0  , &DT, np, GLJ, 0.0, 0.0);

  this -> global2local (x, P, NS);

  Femlib::grad2 (P, P, R, S, DV, DT, np, np, nel);

  for (i = 0; i < nel; i++, R += npnp, S += npnp, P += npnp)
    _elmt[i] -> HelmholtzKern (lambda2, betak2, R, S, P, P);
 
  P -= ntot;
  R -= ntot;
  S -= ntot;
  
  Femlib::grad2 (R, S, P, P, DT, DV, np, np, nel);

  this -> local2globalSum (P, y, NS);
}


void Field::buildRHS (real_t*                  force ,
		      const real_t*            bc    ,
		      real_t*                  RHS   ,
		      real_t*                  RHSint,
		      const real_t**           hbi   ,
		      const int_t              nsolve,
		      const int_t              nzero ,
		      const vector<Boundary*>& bnd   ,
		      const NumberSys*         N     ,
		      real_t*                  work  ) const
// ---------------------------------------------------------------------------
// Build RHS for direct or iterative solution.
//
// Iterative solution is flagged by presence of RHSint, a pointer to
// element-internal node storage.  If RHSint is zero, then hbi, a vector
// of pointers to element interior/exterior coupling matrices, must be 
// non-zero.
//
// Compute RHS vector for direct solution of Helmholtz problem as
//
//                      - M f - H g + <h, w>.
// 
// On input, force contains a plane of
//
//                          - M f - H g
//
// in element (row-major) form, and bc contains the line of BC data values
// for this plane of data: only natural BCs are used in formation of <h, w>.
// ---------------------------------------------------------------------------
{
  const int_t              np      = Geometry::nP();
  const int_t              nel     = Geometry::nElmt();
  const int_t              next    = Geometry::nExtElmt();
  const int_t              nint    = Geometry::nIntElmt();
  const int_t              npnp    = Geometry::nTotElmt();
  const int_t              nglobal = N -> nGlobal();
  const int_t*             gid;
  register const Boundary* B;
  register int_t           i, boff;

  if   (RHSint) Veclib::zero (nglobal + Geometry::nInode(), RHS, 1);
  else          Veclib::zero (nglobal,                      RHS, 1);

  // -- Add in contribution from forcing f = - M f - H g.

  for (gid = N -> btog(), i = 0; i < nel; i++, force += npnp, gid += next) {
    if (RHSint) {
      _elmt[i] -> local2globalSum   (force, gid, RHS, RHSint); RHSint += nint;
    } else
      _elmt[i] -> local2globalSumSC (force, gid, RHS, hbi[i], work);
  }

  // -- Add in <h, w>.

  for (gid = N -> btog(), i = 0; i < _nbound; i++, bc += np) {
    B    = bnd[i];
    boff = B -> bOff();

    B -> sum (bc, gid + boff, work, RHS);
  }

  // -- Zero any contribution that <h, w> made to essential BC nodes.

  Veclib::zero (nzero, RHS + nsolve, 1);
}


void Field::local2global (const real_t*    src,
			  real_t*          tgt,
			  const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (src) into globally-numbered tgt, with element-
// boundary values in the first N -> nGlobal() places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      next = Geometry::nExtElmt();
  const int_t      nint = Geometry::nIntElmt();
  const int_t      npnp = Geometry::nTotElmt();
  const int_t*     gid  = N -> btog();
  register int_t   i;
  register real_t* internal = tgt + N -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> local2global (src, gid, tgt, internal);
}


void Field::global2local (const real_t*    src,
			  real_t*          tgt,
			  const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (tgt) from src, which has globally-numbered
// (element- boundary) values in the first nglobal places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int_t          nel  = Geometry::nElmt();
  const int_t          next = Geometry::nExtElmt();
  const int_t          nint = Geometry::nIntElmt();
  const int_t          npnp = Geometry::nTotElmt();
  const int_t*         gid  = N -> btog();
  register int_t       i;
  register const real_t* internal = src + N -> nGlobal();

  for (i = 0; i < nel; i++, tgt += npnp, gid += next, internal += nint)
    _elmt[i] -> global2local (tgt, gid, src, internal);
}


void Field::local2globalSum (const real_t*    src,
			     real_t*          tgt,
			     const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Direct stiffness sum a plane of data (src) into globally-numbered
// tgt, with element- boundary values in the first N -> nGlobal()
// places, followed by element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      next = Geometry::nExtElmt();
  const int_t      nint = Geometry::nIntElmt();
  const int_t      npnp = Geometry::nTotElmt();
  const int_t*     gid  = N -> btog();
  register int_t   i;
  register real_t* internal = tgt + N -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> local2globalSum (src, gid, tgt, internal);
}


void Field::getEssential (const real_t*            src,
			  real_t*                  tgt,
			  const vector<Boundary*>& bnd,
			  const NumberSys*         N  ) const
// ---------------------------------------------------------------------------
// On input, src contains a line of BC values for the current data plane.
// Scatter current values of essential BCs into globally-numbered tgt.
//
// The construction of the essential BCs has to account for cases
// where element corners may touch the domain boundary but do not have
// an edge along a boundary.  This is done by working with a
// globally-numbered vector.
// ---------------------------------------------------------------------------
{
  const int_t              np = Geometry::nP();
  const int_t*             btog = N -> btog();
  register const Boundary* B;
  register int_t           i, boff;
  
  for (i = 0; i < _nbound; i++, src += np) {
    B    = bnd[i];
    boff = B -> bOff();
  
    B -> set (src, btog + boff, tgt);
  }
}


void Field::setEssential (const real_t*    src,
			  real_t*          tgt,
			  const NumberSys* N  )
// ---------------------------------------------------------------------------
// Gather globally-numbered src into essential BC nodes of current
// data plane.
// ---------------------------------------------------------------------------
{
  const int_t    nel  = Geometry::nElmt();
  const int_t    next = Geometry::nExtElmt();
  const int_t    npnp = Geometry::nTotElmt();
  const int_t*   emask = N -> emask();
  const int_t*   bmask = N -> bmask();
  const int_t*   btog  = N -> btog();
  register int_t i;

  for (i = 0; i < nel; i++, bmask += next, btog += next, tgt += npnp)
    if (emask[i]) _elmt[i] -> bndryMask (bmask, tgt, src, btog);
}


void Field::coupleBCs (Field*      v  ,
		       Field*      w  ,
		       const int_t dir)
// ---------------------------------------------------------------------------
// Couples/uncouple boundary condition values for the radial and
// azimuthal velocity fields in cylindrical coordinates, depending on
// indicated direction.  This action is required due to the coupling
// in the viscous terms of the N--S equations in cylindrical coords.
//
// dir == +1
// ---------
//           v~ <-- v + i w
//           w~ <-- v - i w
// dir == -1
// ---------
//           v  <-- 0.5   * (v~ + w~)
//           w  <-- 0.5 i * (w~ - v~)
//
// Since there is no coupling for the viscous terms in the 2D
// equation, do nothing for the zeroth Fourier mode.
//
// NB: modified for dual.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "Field::couple";
  const int_t    nL = v -> _nline;
  vector<real_t> work (nL);
  real_t         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
  if (dir == FORWARD) {

    Vr = v -> _line[1];
    Vi = v -> _line[2];
    Wr = w -> _line[1];
    Wi = w -> _line[2];

    Veclib::copy (nL, Vr, 1, tp, 1);
    Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
    Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
    Veclib::copy (nL, Wr, 1, tp, 1);
    Veclib::copy (nL, Wi, 1, Wr, 1);
    Veclib::vsub (nL, Vi, 1, tp, 1, Wi, 1);
    Veclib::vadd (nL, Vi, 1, tp, 1, Vi, 1);

  } else if (dir == INVERSE) {

    Vr = v -> _line[1];
    Vi = v -> _line[2];
    Wr = w -> _line[1];
    Wi = w -> _line[2];

    Veclib::copy  (nL,      Vr, 1, tp, 1);
    Veclib::svvpt (nL, 0.5, Vr, 1, Wr, 1, Vr, 1);
    Veclib::svvmt (nL, 0.5, Wr, 1, tp, 1, Wr, 1);
    Veclib::copy  (nL,      Wi, 1, tp, 1);
    Veclib::copy  (nL,      Wr, 1, Wi, 1);
    Veclib::svvpt (nL, 0.5, Vi, 1, tp, 1, Wr, 1);
    Veclib::svvmt (nL, 0.5, Vi, 1, tp, 1, Vi, 1);

  } else
    message (routine, "unknown direction given", ERROR);
}


real_t Field::modeConstant (const char   name,
			    const int_t  mode,
			    const real_t beta)
// ---------------------------------------------------------------------------
// For cylindrical coordinates & 3D, the radial and azimuthal fields
// are coupled before solution of the viscous step.  This means that
// the Fourier constant used for solution may vary from that which
// applies to the axial component.
//
// For Field v~, betak -> betak + 1 while for w~, betak -> betak - 1.
//
// For the uncoupled Fields v, w solved for the zeroth Fourier mode,
// the "Fourier" constant in the Helmholtz equations is 1.
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim()    <          3          ||
      Geometry::system() == Geometry::Cartesian || 
      name               ==         'c'         ||
      name               ==         'p'         ||
      name               ==         'u'          ) return beta * mode;

  if      (name == 'v') return (mode == 0) ? 1.0 : beta * mode + 1.0;
  else if (name == 'w') return (mode == 0) ? 1.0 : beta * mode - 1.0;
  else message ("Field::modeConstant", "unrecognized Field name given", ERROR);

  return -1.0;			// -- Never happen.
}
