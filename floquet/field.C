///////////////////////////////////////////////////////////////////////////////
// field.C: derived from AuxField, Field adds boundary conditions,
// global numbering, and the ability to solve Helmholtz problems.
//
// Copyright (c) 1994,2003 Hugh Blackburn
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
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "Sem.h"


Field::Field (BoundarySys*      B,
	      real*             M,
	      const int         N,
	      vector<Element*>& E,
	      const char        C) :
// ---------------------------------------------------------------------------
// Create storage for a new Field from scratch.
// ---------------------------------------------------------------------------
  AuxField (M, N, E, C),
  _bsys    (B)
{
  const int                np  = Geometry::nP();
  const int                npr = Geometry::nProc();
  const int                nzb = Geometry::basePlane();
  const vector<Boundary*>& BC  = _bsys -> BCs (0);
  register real*           p;
  register int             i, k;

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
    int                      i;
  
    cout << "# -- Field '" << F -> name() << "' Boundary Information:" << endl;
    if (!F -> _nbound) cout << "No BCs for this Field" << endl;
    for (i = 0; i < F -> _nbound; i++) BC[i] -> print();
  }
}


void Field::evaluateBoundaries (const int step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
// Note that for 3D this evaluation is done in Fourier-transformed space.
//
// This routine is only meant to be used for boundary conditions that must
// be re-evaluated at every step, such as high-order pressure BCs.
// ---------------------------------------------------------------------------
{
  const int         np = Geometry::nP();
  vector<Boundary*> BC;
  real*             p;
  register int      i, k;

  if (Geometry::nPert() == 2)
    BC = _bsys -> BCs (0);
  else
    BC = _bsys -> BCs (static_cast<int>(Femlib::value ("K_FUND")));

  for (k = 0; k < _nz; k++)
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BC[i] -> evaluate (k, step, p);
}


void Field::evaluateM0Boundaries (const int step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind, but only for Mode 0.
// ---------------------------------------------------------------------------
{
  ROOTONLY {
    const vector<Boundary*>& BC = _bsys -> BCs (0);
    const int                np = Geometry::nP();
    real*                    p;
    register int             i;

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
    const int                np = Geometry::nP();
    real*                    p;
    register int             i;

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
  const int        nel     = Geometry::nElmt();
  const int        npnp    = Geometry::nTotElmt();
  const int        next    = Geometry::nExtElmt();
  const NumberSys* N       = _bsys -> Nsys  (0);
  const real*      imass   = _bsys -> Imass (0);
  const int        nglobal = N     -> nGlobal();
  const integer*   btog    = N     -> btog();
  const integer*   gid;
  register int     i, k;
  vector<real>     work (nglobal);
  real             *src, *dssum = &work[0];

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


void Field::smooth (const int nZ ,
		    real*     tgt) const
// ---------------------------------------------------------------------------
// Smooth tgt field along element boundaries using *this, with
// mass-average smoothing.  Tgt is assumed to be arranged by planes, with
// planeSize() offset between each plane of data.
// ---------------------------------------------------------------------------
{
  const int        nel     = Geometry::nElmt();
  const int        npnp    = Geometry::nTotElmt();
  const int        next    = Geometry::nExtElmt();
  const int        nP      = Geometry::planeSize();
  const NumberSys* N       = _bsys -> Nsys  (0);
  const real*      imass   = _bsys -> Imass (0);
  const int        nglobal = N    -> nGlobal();
  const integer*   btog    = N    -> btog();
  const integer*   gid;
  register int     i, k;
  vector<real>     work (nglobal);
  real             *src, *dssum = &work[0];

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
  vector<real>             work(4 * Geometry::nP());
  real                     F = 0.0;
  register int             i;
  
  for (i = 0; i < C -> _nbound; i++)
    F += BC[i] -> flux ("wall", C -> _data, &work[0]);

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
  const int                nsurf = P -> _nbound;
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real>             work(Geometry::nP());
  register int             i;
  
  for (i = 0; i < nsurf; i++) {
    secF = BC[i] -> normalTraction ("wall", P -> _data, &work[0]);
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
  const vector<Boundary*>& UBC =       U->_bsys->BCs(0);
  const vector<Boundary*>& WBC = (W) ? W->_bsys->BCs(0) : (vector<Boundary*>)0;
  const int                np      = Geometry::nP();
  const int                nbound  = U -> _nbound;
  const real               mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real>             work(4 * np);
  register int             i;

  for (i = 0; i < nbound; i++) {
    secF = UBC[i] -> tangentTraction  ("wall", U->_data, V->_data, &work[0]);
    F.x        -= mu * secF.x;
    F.y        -= mu * secF.y;
    if (W) F.z -= mu * WBC[i] -> flux ("wall", W->_data, &work[0]);
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
  const int                np    = Geometry::nP();
  const int                nz    = Geometry::nZProc();
  const int                nsurf = P -> _nbound;
  Vector                   secF;
  vector<real>             work(np);
  real                     *p;
  register int             i, j;
  
  for (j = 0; j < nz; j++) {
    p = P -> _plane[j];
    for (i = 0; i < nsurf; i++) {
      secF = BC[i] -> normalTraction ("wall", p, &work[0]);
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
  const vector<Boundary*>& UBC =       U->_bsys->BCs(0);
  const vector<Boundary*>& WBC = (W) ? W->_bsys->BCs(0) : (vector<Boundary*>)0;
  const int                np      = Geometry::nP();
  const int                nz      = Geometry::nZProc();
  const int                nbound = U -> _nbound;
  const real               mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF;
  vector<real>             work(4 * np);
  real                     *u, *v, *w;
  register int             i, j;

  for (j = 0; j < nz; j++) {
    u = U -> _plane[j];
    v = V -> _plane[j];
    w = (W) ? W -> _plane[j] : 0;
    for (i = 0; i < nbound; i++) {
      secF = UBC[i] -> tangentTraction ("wall", u, v, &work[0]);
             fx[j] -= mu * secF.x;
             fy[j] -= mu * secF.y;
      if (W) fz[j] -= mu * WBC[i] -> flux ("wall", w, &work[0]);
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
  const int  nel    = Geometry::nElmt();
  const int  next   = Geometry::nExtElmt();
  const int  npnp   = Geometry::nTotElmt();
  const int  ntot   = Geometry::nPlane();

  int        i, k;

  const vector<Boundary*>& B       = M -> _BC;
  const NumberSys*         N       = M -> _NS;
  real                     lambda2 = M -> _HelmholtzConstant;
  real                     betak2  = M -> _FourierConstant;
  int                      nsolve  = M -> _nsolve;
  int                      nglobal = M -> _nglobal;
  int                      nzero   = nglobal - nsolve;

  for (k = 0; k < _nz; k++) {
    real* forcing = f -> _plane[k];
    real* unknown = _plane     [k];
    real* bc      = _line      [k];

    switch (M -> _method) {
    
    case DIRECT: {
      const real*    H     = const_cast<const real*>   (M -> _H);
      const real**   hii   = const_cast<const real**>  (M -> _hii);
      const real**   hbi   = const_cast<const real**>  (M -> _hbi);
      const integer* b2g   = const_cast<const integer*>(N -> btog());
      int            nband = M -> _nband;

      vector<real>   work (nglobal + 4*npnp);
      real           *RHS = &work[0], *tmp = RHS + nglobal;
      int            info;
      
      // -- Build RHS = - M f - H g + <h, w>.

      Veclib::zero (nglobal, RHS, 1);
      
      this -> getEssential (bc, RHS, B, N);
      this -> constrain    (forcing, lambda2,betak2,RHS,N,tmp);
      this -> buildRHS     (forcing, bc,RHS,0,hbi,nsolve,nzero,B,N,tmp);
      
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
      const int     StepMax = static_cast<int>(Femlib::value ("STEP_MAX"));
      const int     npts    = M -> _npts;
      real          alpha, beta, dotp, epsb2, r2, rho1, rho2 = 0.0;
      vector<real>  work (5*npts + 4*Geometry::nTotElmt());

      const int mode =
	(((Geometry::problem()) == Geometry::O2_2D)?0:1) * Geometry::kFund();

      real* r   = &work[0];
      real* p   = r + npts;
      real* q   = p + npts;
      real* x   = q + npts;
      real* z   = x + npts;
      real* wrk = z + npts;

      Veclib::zero (nglobal, x, 1);

      this -> getEssential (bc, x, B, N);  
      this -> constrain    (forcing, lambda2,betak2,x,N,wrk);
      this -> buildRHS     (forcing, bc,r,r+nglobal,0,nsolve,nzero,B,N,wrk);

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
  
      if (static_cast<int>(Femlib::value ("VERBOSE")) > 1) {
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
		       const NumberSys* N      ,
		       real*            work   ) const
// ---------------------------------------------------------------------------
// Replace f's data with constrained weak form of forcing: - M f - H g.
// On input, essential BC values (g) have been loaded into globally-numbered
// esstlbc, other values are zero.
//
// Input vector work should be 4*Geometry::nTotElmt() long.
// ---------------------------------------------------------------------------
{
  const int         np    = Geometry::nP();
  const int         nel   = Geometry::nElmt();
  const int         next  = Geometry::nExtElmt();
  const int         npnp  = Geometry::nTotElmt();
  const int         ntot  = Geometry::nPlane();
  const integer*    emask = N -> emask();
  const integer*    btog  = N -> btog();
  const real        **DV, **DT;
  register Element* E;
  register int      i;
  real              *u = work, *tmp = work + npnp;

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


void Field::HelmholtzOperator (const real*x      ,
			       real*      y      ,
			       const real lambda2,
			       const real betak2 ,
			       const int  mode   ,
			       real*      work   ) const
// ---------------------------------------------------------------------------
// Discrete 2D global Helmholtz operator which takes the vector x into
// vector y, including direct stiffness summation.  Vectors x & y have 
// global ordering: that is, with nglobal (element edge nodes, with
// redundancy removed) coming first, followed by nel blocks of element-
// internal nodes.
//
// Vector work must have length 4*Geometry::nTotElmt().
// ---------------------------------------------------------------------------
{
  const int        np      = Geometry::nP();
  const int        nel     = Geometry::nElmt();
  const int        npnp    = Geometry::nTotElmt();
  const int        next    = Geometry::nExtElmt();
  const int        nint    = Geometry::nIntElmt();
  const int        ntot    = Geometry::nPlane();
  const NumberSys* NS      = _bsys -> Nsys (mode);
  const integer*   gid     = NS -> btog();
  const int        nglobal = NS -> nGlobal() + Geometry::nInode();
  const real*      xint    = x + NS -> nGlobal();
  real*            yint    = y + NS -> nGlobal();
  real             *P = work, *tmp = work + npnp;
  register int     i;
  Element*         E;

  Veclib::zero (nglobal, y, 1);

  // -- Add in contributions from mixed BCs while x & y are global vectors.

  if (_bsys -> mixBC()) {
    const vector<Boundary*>& BC   = _bsys -> BCs (0);
    const integer*           bmap = NS    -> btog();

    for (i = 0; i < _nbound; i++)
      BC[i] -> augmentOp (bmap, x, y);
  }

  // -- Add in contributions from elemental Helmholtz operations.

  for (i = 0; i < nel; i++, gid += next, xint += nint, yint += nint) {
    E = _elmt[i];
    E -> global2local    (P, gid, x, xint);
    E -> HelmholtzOp     (lambda2, betak2, P, P, tmp);
    E -> local2globalSum (P, gid, y, yint);
  }
}


void Field::buildRHS (real*                    force ,
		      const real*              bc    ,
		      real*                    RHS   ,
		      real*                    RHSint,
		      const real**             hbi   ,
		      const int                nsolve,
		      const int                nzero ,
		      const vector<Boundary*>& bnd   ,
		      const NumberSys*         N     ,
		      real*                    work  ) const
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
//
// Input vector work should be Geometry::nTotElmt() long.
// ---------------------------------------------------------------------------
{
  const int                np      = Geometry::nP();
  const int                nel     = Geometry::nElmt();
  const int                next    = Geometry::nExtElmt();
  const int                nint    = Geometry::nIntElmt();
  const int                npnp    = Geometry::nTotElmt();
  const int                nglobal = N -> nGlobal();
  const integer*           gid;
  register const Boundary* B;
  register int             i, boff;

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


void Field::local2global (const real*      src,
			  real*            tgt,
			  const NumberSys* N  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (src) into globally-numbered tgt, with element-
// boundary values in the first N -> nGlobal() places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  const int      nel  = Geometry::nElmt();
  const int      next = Geometry::nExtElmt();
  const int      nint = Geometry::nIntElmt();
  const int      npnp = Geometry::nTotElmt();
  const integer* gid  = N -> btog();
  register int   i;
  register real* internal = tgt + N -> nGlobal();

  for (i = 0; i < nel; i++, src += npnp, gid += next, internal += nint)
    _elmt[i] -> local2global (src, gid, tgt, internal);
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
  const int            nel  = Geometry::nElmt();
  const int            next = Geometry::nExtElmt();
  const int            nint = Geometry::nIntElmt();
  const int            npnp = Geometry::nTotElmt();
  const integer*       gid  = N -> btog();
  register int         i;
  register const real* internal = src + N -> nGlobal();

  for (i = 0; i < nel; i++, tgt += npnp, gid += next, internal += nint)
    _elmt[i] -> global2local (tgt, gid, src, internal);
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
  const int                np = Geometry::nP();
  const integer*           btog = N -> btog();
  register const Boundary* B;
  register int             i, boff;
  
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
  const int      nel  = Geometry::nElmt();
  const int      next = Geometry::nExtElmt();
  const int      npnp = Geometry::nTotElmt();
  const integer* emask = N -> emask();
  const integer* bmask = N -> bmask();
  const integer* btog  = N -> btog();
  register int   i;

  for (i = 0; i < nel; i++, bmask += next, btog += next, tgt += npnp)
    if (emask[i]) _elmt[i] -> bndryMask (bmask, tgt, src, btog);
}


void Field::coupleBCs (Field*    v  ,
		       Field*    w  ,
		       const int dir)
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
  if (Geometry::problem() == Geometry::O2_2D) return;

  const char   routine[] = "Field::couple";
  const int    nL        =  v -> _nline;
  vector<real> work (nL);
  real         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
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


real Field::modeConstant (const char name,
			  const int  mode,
			  const real beta)
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
  const char routine[] = "Field::modeConstant";

  if (Geometry::problem() == Geometry::O2_2D     ||
      Geometry::system()  == Geometry::Cartesian || 
      name                ==         'c'         ||
      name                ==         'p'         ||
      name                ==         'u'          ) return beta * mode;

  if      (name == 'v') return beta * mode + 1.0;
  else if (name == 'w') return beta * mode - 1.0;
  else message (routine, "unrecognized Field name given", ERROR);

  return -1.0;			// -- Never happen.
}
