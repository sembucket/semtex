///////////////////////////////////////////////////////////////////////////////
// field.C: derived from AuxField, Field adds boundary conditions,
// global numbering, and the ability to solve Helmholtz problems.
//
// Copyright (C) 1994, 2001 Hugh Blackburn
//
// HELMHOLTZ PROBLEMS
// ------------------
// Solve routines provide solution to the discrete form of the Helmholtz eqn
//                      
//                       div grad u - \lambda^2 u = f,
//
// on domain \Omega, subject to essential BCs u = g on \Gamma_g and
// natural BCs \partial u / \partial n = h on \Gamma_h, where the
// boundary \Gamma of \Omega is the union of (non-overlapping)
// \Gamma_g and \Gamma_h and n is the unit outward normal vector on
// \Gamma.  \lambda^2 is called the Helmholtz constant below.
//
// The Galerkin form, using integration by parts with weighting functions w
// which are zero on \Gamma_g, is
//
//            (grad u, grad w) + \lambda^2 (u, w) = - (f, w) + <h, w>
//
// where (a, b) = \int a.b d\Omega is an integration over the domain and
//       <a, b> = \int a.b d\Gamma is an integration over the domain boundary.
//
// The discrete (finite element) equivalent is to solve
//
//                   K.u + \lambda^2 M.u = - M.f + <h, w>
//
// or
//
//                         H.u = - M.f + <h, w>
//
// where K, M and H are respectively (assembled) "stiffness", "mass"
// and Helmholtz matrices.
//
// Some complications arise from dealing with essential boundary
// conditions, since typically the elemental matrices K^e, M^e which
// are assembled to form K and M do not account for the boundary
// requirements on w.  There are a number of ways of dealing with this
// issue: one approach is to partition H as it is formed (here F =
// -M.f + <h, w>):
//
//   +--------+-------------+ /  \     /  \
//   |        |             | |  |     |  |
//   |   Hp   |     Hc      | |u |     |F |   u : nodal values for solution.
//   |        |(constraint) | | s|     | s|    s
//   |        |             | |  |     |  |       (n_solve values)
//   +--------+-------------+ +--+     +--+
//   |        | H_ess: this | |  |  =  |  |
//   |        | partition   | |  |     |  |
//   |    T   | relates to  | |u |     |F |   u : are given essential BCs.
//   |  Hc    | essential   | | g|     | g|    g
//   |        | BC nodes    | |  |     |  |       (n_global - n_solve values)
//   |        | and is not  | |  |     |  |
//   |        | assembled.  | |  |     |  |
//   +--------+-------------+ \  /     \  /
//
// Partition out the sections of the matrix corresponding to the known
// nodal values (essential BCs), and solve instead the constrained
// problem
//
//   +--------+               /  \     /  \     +-------------+ /  \
//   |        |               |  |     |  |     |             | |  |
//   |   Hp   |               |u |     |F |     |     Hc      | |  |
//   |        |               | s|  =  | s|  -  |             | |u |
//   |        |               |  |     |  |     |             | | g|
//   +--------+               \  /     \  /     +-------------+ |  |.
//                                                              |  |
//                                                              |  |
//                                                              \  /
//
// Here n_global is the number of nodes that receive global node
// numbers, typically those on the mesh edges.  N_solve is the number
// of these nodes that have values that must be solved for,
// i.e. n_global minus the number of global nodes situated on
// essential-type boundaries.
//
// An alternative approach (USED HERE) to the constrained problem is to let
//
//                      u = v + g     (v = 0 on \Gamma_g)
//
// and solve instead
//
//                      H.v = - M.f - H.g + <h, w>
//
// (where only the partition Hp is needed for the matrix H in the
// LHS), then afterwards compose u = v + g to get the full solution.
// The advantage to this method is that the constraint partition Hc
// does not need to be assembled or stored.  The operations M.f and
// H.g can be performed on an element-by-element basis.
//
// FIELD NAMES
// -----------
// The (one character) names of field variables are significant, and have
// the following reserved meanings:
// 
// u:  First velocity component.            (Cylindrical: axial     velocity.)
// v:  Second velocity component.           (Cylindrical: radial    velocity.)
// w:  Third velocity component.            (Cylindrical: azimuthal velocity.)
// p:  Pressure divided by density.
// c:  Scalar for transport or elliptic problems.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


Field::Field (BoundarySys*      B,
	      real*             M,
	      const integer     N,
	      vector<Element*>& E,
	      const char        C) :
// ---------------------------------------------------------------------------
// Create storage for a new Field from scratch.
// ---------------------------------------------------------------------------
  AuxField (M, N, E, C),
  _bsys    (B)
{
  const integer            np  = Geometry::nP();
  const integer            npr = Geometry::nProc();
  const integer            nzb = Geometry::basePlane();
  const vector<Boundary*>& BC  = _bsys -> BCs (0);
  const real               dz  = Femlib::value ("TWOPI / BETA / N_Z");
  register real*           p;
  register integer         i, k;

  // -- Allocate storage for boundary data, round up for Fourier transform.
  
  _nbound = _bsys -> nSurf();
  _nline  = _nbound * np;
  if   (npr > 1) _nline += 2 * npr - _nline % (2 * npr);
  else           _nline += _nline % 2;

  _line  = new real* [static_cast<size_t>(_nz)];
  _sheet = new real  [static_cast<size_t>(_nz * _nline)];

  for (k = 0; k < _nz; k++) _line[k] = _sheet + k*_nline;

  Veclib::zero (_nz * _nline, _sheet, 1);

  // -- Set values for boundary data, but enforce z = 0.

  for (k = 0; k < _nz; k++) {
    Femlib::value ("z", 0.0);
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (k, 0, p);
  }

  // -- Do NOT Fourier transform boundary data; already in Fourier space.
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
    integer                  i;
  
    cout << "# -- Field '" << F -> name() << "' Boundary Information:" << endl;
    if (!F -> _nbound) cout << "No BCs for this Field" << endl;
    for (i = 0; i < F -> _nbound; i++) BC[i] -> print();
  }
}


void Field::evaluateBoundaries (const integer step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
// Note that for 3D this evaluation is done in Fourier-transformed space.
//
// This routine is only meant to be used for boundary conditions that must
// be re-evaluated at every step, such as high-order pressure BCs.
// ---------------------------------------------------------------------------
{
  const integer     np = Geometry::nP();
  vector<Boundary*> BC;
  real*             p;
  register integer  i, k;

  if (Geometry::nPert() == 2)
    BC = _bsys -> BCs (0);
  else
    BC = _bsys -> BCs (static_cast<integer>(Femlib::value ("K_FUND")));

  for (k = 0; k < _nz; k++)
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (k, step, p);
}


void Field::evaluateM0Boundaries (const integer step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind, but only for Mode 0.
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> BCs (0);
    const integer            np = Geometry::nP();
    real*                    p;
    register integer         i;

    for (p = _line[0], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (0, step, p);
  }
}


void Field::addToM0Boundaries (const real  val,
			       const char* grp)
// ---------------------------------------------------------------------------
// Add val to zeroth Fourier mode's bc storage area on BC group "grp".
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> BCs (0);
    const integer            np = Geometry::nP();
    real*                    p;
    register integer         i;

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
  const integer    nel     = Geometry::nElmt();
  const integer    npnp    = Geometry::nTotElmt();
  const integer    next    = Geometry::nExtElmt();
  const NumberSys* N       = _bsys -> Nsys  (0);
  const real*      imass   = _bsys -> Imass (0);
  const integer    nglobal = N    -> nGlobal();
  const integer*   btog    = N    -> btog();
  const integer*   gid;
  register integer i, k;
  vector<real>     work (nglobal);
  real             *src, *dssum = work();

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


void Field::smooth (const integer nZ ,
		    real*         tgt) const
// ---------------------------------------------------------------------------
// Smooth tgt field along element boundaries using *this, with
// mass-average smoothing.  Tgt is assumed to be arranged by planes, with
// planeSize() offset between each plane of data.
// ---------------------------------------------------------------------------
{
  const integer    nel     = Geometry::nElmt();
  const integer    npnp    = Geometry::nTotElmt();
  const integer    next    = Geometry::nExtElmt();
  const integer    nP      = Geometry::planeSize();
  const NumberSys* N       = _bsys -> Nsys  (0);
  const real*      imass   = _bsys -> Imass (0);
  const integer    nglobal = N    -> nGlobal();
  const integer*   btog    = N    -> btog();
  const integer*   gid;
  register integer i, k;
  vector<real>     work (nglobal);
  real             *src, *dssum = work();

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


real Field::flux (const Field* C)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute normal flux of field C on all "wall" group boundaries.
//
// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = C -> _bsys -> BCs (0);
  vector<real>             work(3 * Geometry::nP());
  real                     F = 0.0, *tmp = work();
  register integer         i;
  
  for (i = 0; i < C -> _nbound; i++)
    F += BC[i] -> flux ("wall", C -> _data, tmp);

  return F;
}


Vector Field::normalTraction (const Field* P)
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
  const integer            nsurf = P -> _nbound;
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real>             work(Geometry::nP());
  real*                    tmp = work();
  register integer         i;
  
  for (i = 0; i < nsurf; i++) {
    secF = BC[i] -> normalTraction ("wall", P -> _data, tmp);
    F.x += secF.x;
    F.y += secF.y;
  }

  return F;
}


Vector Field::tangentTraction (const Field* U,
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
  const vector<Boundary*>& UBC = U -> _bsys -> BCs(0);
  const vector<Boundary*>& WBC = (W) ? W -> _bsys -> BCs(0) : 0;
  const integer            np      = Geometry::nP();
  const integer            nbound  = U -> _nbound;
  const real               mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real>             work(3 * np);
  real                     *ddx = work(), *ddy = ddx + np;
  register integer         i;

  for (i = 0; i < nbound; i++) {
    secF = UBC[i] -> tangentTraction  ("wall", U->_data, V->_data, ddx, ddy);
    F.x        -= mu * secF.x;
    F.y        -= mu * secF.y;
    if (W) F.z -= mu * WBC[i] -> flux ("wall", W->_data, ddx);
  }

  return F;
}


void Field::normTractionV (real*        fx,
			   real*        fy,
			   const Field* P )
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute normal tractive forces at each z location on wall
// boundaries.  Fx & Fy are assumed to contain sufficient (nZProc)
// storage and to be zero on entry.  P could be in physical space or
// Fourier transformed.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC    = P -> _bsys -> BCs (0);
  const integer            np    = Geometry::nP();
  const integer            nz    = Geometry::nZProc();
  const integer            nsurf = P -> _nbound;
  Vector                   secF;
  vector<real>             work(np);
  real                     *p, *tmp = work();
  register integer         i, j;
  
  for (j = 0; j < nz; j++) {
    p = P -> _plane[j];
    for (i = 0; i < nsurf; i++) {
      secF = BC[i] -> normalTraction ("wall", p, tmp);
      fx[j] += secF.x;
      fy[j] += secF.y;
    }
  }
}


void Field::tangTractionV (real*        fx,
			   real*        fy,
			   real*        fz,
			   const Field* U ,
			   const Field* V ,
			   const Field* W )
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute tangential tractive forces at each z location on wall
// boundaries.  Fx, fy, fz assumed to contain sufficient storage and
// be zero on entry.  U, V, W could be in physical space or Fourier
// transformed.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& UBC = U -> _bsys -> BCs(0);
  const vector<Boundary*>& WBC = (W) ? W -> _bsys -> BCs(0) : 0;
  const integer            np      = Geometry::nP();
  const integer            nz      = Geometry::nZProc();
  const integer            nbound = U -> _nbound;
  const real               mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF;
  vector<real>             work(3 * np);
  real                     *ddx, *ddy, *tmp = work(), *u, *v, *w;
  register integer         i, j;

  ddx = work();
  ddy = ddx + np;

  for (j = 0; j < nz; j++) {
    u = U -> _plane[j];
    v = V -> _plane[j];
    w = (W) ? W -> _plane[j] : 0;
    for (i = 0; i < nbound; i++) {
      secF = UBC[i] -> tangentTraction ("wall", u, v, ddx, ddy);
             fx[j] -= mu * secF.x;
             fy[j] -= mu * secF.y;
      if (W) fz[j] -= mu * WBC[i] -> flux ("wall", w, tmp);
    }
  }
}


Field& Field::solve (AuxField*        f,
		     const MatrixSys* M)
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
// ---------------------------------------------------------------------------
{
  const char routine[] = "Field::solve";
  const integer nel    = Geometry::nElmt();
  const integer next   = Geometry::nExtElmt();
  const integer npnp   = Geometry::nTotElmt();
  const integer ntot   = Geometry::nPlane();
  integer       i, k;

  const vector<Boundary*>& B       = M -> _BC;
  const NumberSys*         N       = M -> _NS;
  real                     lambda2 = M -> _HelmholtzConstant;
  real                     betak2  = M -> _FourierConstant;
  integer                  nsolve  = M -> _nsolve;
  integer                  nglobal = M -> _nglobal;
  integer                  nzero   = nglobal - nsolve;

  for (k = 0; k < _nz; k++) {
    real* forcing = f -> _plane[k];
    real* unknown = _plane     [k];
    real* bc      = _line      [k];

    switch (M -> _method) {
    
    case DIRECT: {
      const real*    H       = const_cast<const real*>   (M -> _H);
      const real**   hii     = const_cast<const real**>  (M -> _hii);
      const real**   hbi     = const_cast<const real**>  (M -> _hbi);
      const integer* b2g     = const_cast<const integer*>(N -> btog());
      integer        nband   = M -> _nband;

      vector<real>   work (nglobal + npnp);
      real           *RHS = work(), *tmp = RHS + nglobal;
      integer        info;
      
      // -- Build RHS = - M f - H g + <h, w>.

      Veclib::zero (nglobal, RHS, 1);
      
      this -> getEssential (bc, RHS, B, N);
      this -> constrain    (forcing, lambda2, betak2, RHS, N);
      this -> buildRHS     (forcing, bc, RHS, 0, hbi, nsolve, nzero, B, N);
      
      // -- Solve for unknown global-node values (if any).
    
      if (nsolve) Lapack::pbtrs("U",nsolve,nband-1,1,H,nband,RHS,nglobal,info);
      
      // -- Carry out Schur-complement solution for element-internal nodes.
      
      for (i = 0; i < nel; i++, b2g += next, forcing += npnp, unknown += npnp)
	_elmt[i] -> g2eSC (RHS, b2g, forcing, unknown, hbi[i], hii[i], tmp);

      unknown -= ntot;
    
      // -- Scatter-gather essential BC values into plane.

      Veclib::zero (nglobal, RHS, 1);
    
      this -> getEssential (bc, RHS, B,   N);
      this -> setEssential (RHS, unknown, N);
    }
    break;

    case JACPCG: {
      const integer StepMax = static_cast<integer>(Femlib::value ("STEP_MAX"));
      const integer npts    = M -> _npts;
      real          alpha, beta, dotp, epsb2, r2, rho1, rho2 =0.0;
      vector<real>  work (5 * npts + 3 * Geometry::nPlane());

      real* r   = work();
      real* p   = r + npts;
      real* q   = p + npts;
      real* x   = q + npts;
      real* z   = x + npts;
      real* wrk = z + npts;

      Veclib::zero (nglobal, x, 1);

      this -> getEssential (bc, x, B, N);  
      this -> constrain    (forcing, lambda2, betak2, x, N);
      this -> buildRHS     (forcing, bc, r, r+nglobal, 0, nsolve, nzero, B, N);

      epsb2  = Femlib::value ("TOL_REL") * sqrt (Blas::dot (npts, r, 1, r, 1));
      epsb2 *= epsb2;

      // -- Build globally-numbered x from element store.

      this -> local2global (unknown, x, N);
  
      // -- Compute first residual using initial guess: r = b - Ax.
      //    And mask to get residual for the zero-BC problem.

      Veclib::zero (nzero, x + nsolve, 1);   
      Veclib::copy (npts,  x, 1, q, 1);
    
      this -> HelmholtzOperator (q, p, lambda2, betak2, wrk, 0);

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

	this -> HelmholtzOperator (p, q, lambda2, betak2, wrk, 0);
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
  
      if ((integer) Femlib::value ("VERBOSE") > 1) {
	char s[StrMax];
	sprintf (s, ":%3d iterations, field '%c'", i, _name);
	message (routine, s, REMARK);
      }
    }
    break;
  
    default:
      message (routine, "called with a method that isn't implemented", ERROR);
      break;
    }
  }
  return *this; 
}

void Field::constrain (real*            force  ,
		       const real       lambda2,
 		       const real       betak2 ,
		       const real*      esstlbc,
		       const NumberSys* N      ) const
// ---------------------------------------------------------------------------
// Replace f's data with constrained weak form of forcing: - M f - H g.
// On input, essential BC values (g) have been loaded into globally-numbered
// esstlbc, other values are zero.
// ---------------------------------------------------------------------------
{
  const integer     np    = Geometry::nP();
  const integer     nel   = Geometry::nElmt();
  const integer     next  = Geometry::nExtElmt();
  const integer     npnp  = Geometry::nTotElmt();
  const integer     ntot  = Geometry::nPlane();
  const integer*    emask = N -> emask();
  const integer*    btog  = N -> btog();
  const real        **DV, **DT;
  register Element* E;
  register integer  i;
  vector<real>      work (3 * npnp);
  real              *u, *r, *s;

  u = work();
  r = u + npnp;
  s = r + npnp;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  // -- Manufacture -(M f + H g).

  for (i = 0; i < nel; i++, emask++, btog += next, force += npnp) {
    E = _elmt[i];
    
    E -> weight (force);	// -- f <-- M f.

    if (*emask) {		// -- f <-- M f + H g.
      Veclib::zero       (3 * npnp, u, 1);
      E -> g2e           (u, btog, esstlbc, 0);
      Femlib::grad2      (u, u, r, s, *DV, *DT, np, 1);
      E -> HelmholtzKern (lambda2, betak2, r, s, u, u);
      Femlib::grad2      (r, s, u, u, *DT, *DV, np, 1);
      Veclib::vadd       (npnp, force, 1, u, 1, force, 1);
    }
  }

  force -= ntot;
  Veclib::neg (ntot, force, 1);
}


void Field::HelmholtzOperator (const real*   x      ,
			       real*         y      ,
			       const real    lambda2,
			       const real    betak2 ,
			       real*         work   ,
			       const integer mode   ) const
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
  const integer    np      = Geometry::nP();
  const integer    nel     = Geometry::nElmt();
  const integer    npnp    = Geometry::nTotElmt();
  const integer    ntot    = Geometry::nPlane();
  const NumberSys* NS      = _bsys -> Nsys (mode);
  const integer    nglobal = NS -> nGlobal() + Geometry::nInode();
  const real       **DV, **DT;
  register integer i;
  real             *P = work, *R = P + ntot, *S = R + ntot;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  Veclib::zero (nglobal,     y, 1);
  Veclib::zero (ntot + ntot, R, 1);

  // -- Add in contributions from mixed BCs while x & y are global vectors.

  if (_bsys -> mixBC()) {
    const vector<Boundary*>& BC   = _bsys -> BCs (0);
    const integer*           bmap = NS    -> btog();

    for (i = 0; i < _nbound; i++)
      BC[i] -> augmentOp (bmap, x, y);
  }

  // -- Add in contributions from elemental Helmholtz operations.

  this -> global2local (x, P, NS);

  Femlib::grad2 (P, P, R, S, *DV, *DT, np, nel);

  for (i = 0; i < nel; i++, R += npnp, S += npnp, P += npnp)
    _elmt[i] -> HelmholtzKern (lambda2, betak2, R, S, P, P);
 
  P -= ntot;
  R -= ntot;
  S -= ntot;
  
  Femlib::grad2 (R, S, P, P, *DT, *DV, np, nel);

  this -> local2globalSum (P, y, NS);
}


void Field::buildRHS (real*                    force ,
		      const real*              bc    ,
		      real*                    RHS   ,
		      real*                    RHSint,
		      const real**             hbi   ,
		      const integer            nsolve,
		      const integer            nzero ,
		      const vector<Boundary*>& bnd   ,
		      const NumberSys*         N     ) const
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
  const integer            np      = Geometry::nP();
  const integer            nel     = Geometry::nElmt();
  const integer            next    = Geometry::nExtElmt();
  const integer            nint    = Geometry::nIntElmt();
  const integer            npnp    = Geometry::nTotElmt();
  const integer            nglobal = N -> nGlobal();
  const integer*           gid;
  register const Boundary* B;
  vector<real>             work (npnp);
  register integer         i, boff;

  if   (RHSint) Veclib::zero (nglobal + Geometry::nInode(), RHS, 1);
  else          Veclib::zero (nglobal,                      RHS, 1);

  // -- Add in contribution from forcing f = - M f - H g.

  for (gid = N -> btog(), i = 0; i < nel; i++, force += npnp, gid += next) {
    if (RHSint) {
      _elmt[i] -> e2gSum   (force, gid, RHS, RHSint); RHSint += nint;
    } else
      _elmt[i] -> e2gSumSC (force, gid, RHS, hbi[i], work());
  }

  // -- Add in <h, w>.

  for (gid = N -> btog(), i = 0; i < _nbound; i++, bc += np) {
    B    = bnd[i];
    boff = B -> bOff();

    B -> sum (bc, gid + boff, work(), RHS);
  }

  // -- Zero any contribution that <h, w> made to essential BC nodes.

  Veclib::zero (nzero, RHS + nsolve, 1);
}


void Field::local2global (const real*      src,
			  real*            tgt,
			  const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (src) into globally-numbered tgt, with element-
// boundary values in the first N -> nGlobal() places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const integer   nel  = Geometry::nElmt();
  const integer   next = Geometry::nExtElmt();
  const integer   nint = Geometry::nIntElmt();
  const integer   npnp = Geometry::nTotElmt();
  const integer*  gid  = N -> btog();
  register integer i;
  register real* internal = tgt + N -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> e2g (src, gid, tgt, internal);
}


void Field::local2globalSum (const real*      src,
			     real*            tgt,
			     const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Direct stiffness sum a plane of data (src) into globally-numbered
// tgt, with element- boundary values in the first N -> nGlobal()
// places, followed by element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const integer    nel  = Geometry::nElmt();
  const integer    next = Geometry::nExtElmt();
  const integer    nint = Geometry::nIntElmt();
  const integer    npnp = Geometry::nTotElmt();
  const integer*   gid  = N -> btog();
  register integer i;
  register real*   internal = tgt + N -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> e2gSum (src, gid, tgt, internal);
}


void Field::global2local (const real*      src,
			  real*            tgt,
			  const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (tgt) from src, which has globally-numbered
// (element- boundary) values in the first nglobal places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const integer        nel  = Geometry::nElmt();
  const integer        next = Geometry::nExtElmt();
  const integer        nint = Geometry::nIntElmt();
  const integer        npnp = Geometry::nTotElmt();
  const integer*       gid  = N -> btog();
  register integer     i;
  register const real* internal = src + N -> nGlobal();

  for (i = 0; i < nel; i++, tgt += npnp, gid += next, internal += nint)
    _elmt[i] -> g2e (tgt, gid, src, internal);
}


void Field::getEssential (const real*              src,
			  real*                    tgt,
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
  const integer            np = Geometry::nP();
  const integer*           btog = N -> btog();
  register const Boundary* B;
  register integer         i, boff;
  
  for (i = 0; i < _nbound; i++, src += np) {
    B    = bnd[i];
    boff = B -> bOff();
  
    B -> set (src, btog + boff, tgt);
  }
}


void Field::setEssential (const real*      src,
			  real*            tgt,
			  const NumberSys* N  )
// ---------------------------------------------------------------------------
// Gather globally-numbered src into essential BC nodes of current
// data plane.
// ---------------------------------------------------------------------------
{
  const integer    nel   = Geometry::nElmt();
  const integer    next  = Geometry::nExtElmt();
  const integer    npnp  = Geometry::nTotElmt();
  const integer*   emask = N -> emask();
  const integer*   bmask = N -> bmask();
  const integer*   btog  = N -> btog();
  register integer i;

  for (i = 0; i < nel; i++, bmask += next, btog += next, tgt += npnp)
    if (emask[i]) _elmt[i] -> bndryMask (bmask, tgt, src, btog);
}


void Field::coupleBCs (Field*        v  ,
		       Field*        w  ,
		       const integer dir)
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
// ---------------------------------------------------------------------------
{
  if (Geometry::nPert() < 3) return;

  const char    routine[] = "Field::couple";
  const integer nL        =  v -> _nline;
  vector<real>  work (nL);
  real          *Vr, *Vi, *Wr, *Wi, *tp = work();
  
  if (dir == FORWARD) {

    if (v->_nz == 1) {		// -- Half-complex.
      Vr = v -> _line[0];
      Wi = w -> _line[0];

      Veclib::copy (nL, Vr, 1, tp, 1);
      Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
    } else {			// -- Full complex.
      Vr = v -> _line[0];
      Vi = v -> _line[1];
      Wr = w -> _line[0];
      Wi = w -> _line[1];

      Veclib::copy (nL, Vr, 1, tp, 1);
      Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nL, Wr, 1, tp, 1);
      Veclib::copy (nL, Wi, 1, Wr, 1);
      Veclib::vsub (nL, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nL, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == INVERSE) {

    if (v->_nz == 1) {		// -- Half-complex.
      Vr = v -> _line[0];
      Wr = w -> _line[0];

      Veclib::copy  (nL,      Vr, 1, tp, 1);
      Veclib::svvpt (nL, 0.5, Vr, 1, Wr, 1, Vr, 1);
      Veclib::svvmt (nL, 0.5, Wr, 1, tp, 1, Wr, 1);
    } else {			// -- Full complex.
      Vr = v -> _line[0];
      Vi = v -> _line[1];
      Wr = w -> _line[0];
      Wi = w -> _line[1];

      Veclib::copy  (nL,      Vr, 1, tp, 1);
      Veclib::svvpt (nL, 0.5, Vr, 1, Wr, 1, Vr, 1);
      Veclib::svvmt (nL, 0.5, Wr, 1, tp, 1, Wr, 1);
      Veclib::copy  (nL,      Wi, 1, tp, 1);
      Veclib::copy  (nL,      Wr, 1, Wi, 1);
      Veclib::svvpt (nL, 0.5, Vi, 1, tp, 1, Wr, 1);
      Veclib::svvmt (nL, 0.5, Vi, 1, tp, 1, Vi, 1);
    }

  } else
    message (routine, "unknown direction given", ERROR);
}


real Field::modeConstant (const char    name,
			  const integer mode,
			  const real    beta)
// ---------------------------------------------------------------------------
// For cylindrical coordinates & 3D, the radial and azimuthal fields
// are coupled before solution of the viscous step.  This means that
// the Fourier constant used for solution may vary from that which
// applies to the axial component.
//
// For Field v~, betak -> betak + 1 while for w~, betak -> betak - 1.
//
// For the uncoupled Fields v, w solved for the zeroth Fourier mode,
// the "Fourier" constant in the Helmholtz equations is +/-1.
// ---------------------------------------------------------------------------
{
  if (Geometry::nPert()   <          3          ||
      Geometry::system() == Geometry::Cartesian || 
      name               ==         'c'         ||
      name               ==         'p'         ||
      name               ==         'u'          ) return beta * mode;

  if      (name == 'v') return beta * mode + 1.0;
  else if (name == 'w') return beta * mode - 1.0;
  else message ("Field::modeConstant", "unrecognized Field name given", ERROR);

  return -1.0;			// -- Never happen.
}
