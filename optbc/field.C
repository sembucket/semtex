///////////////////////////////////////////////////////////////////////////////
// field.C: derived from AuxField, Field adds boundary conditions,
// global numbering, and the ability to solve Helmholtz problems.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
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
  const int_t              npr = Geometry::nProc();
  const int_t              nzb = Geometry::basePlane();
  vector<Boundary*>        BC ;
  register real_t*         p;
  register int_t           i, k;
  
// Ask Hugh why the fundamental mode is used to initialize.
  if (!Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                         BC = _bsys -> BCs (Femlib::ivalue ("BETA"));
  
  // -- Allocate storage for boundary data, round up for Fourier transform.
  
  _nbound = _bsys -> nSurf();
  _nline  = _nbound * np;
  if   (npr > 1) _nline += 2 * npr - _nline % (2 * npr);
  else           _nline += _nline % 2;

  _line  = new real_t* [static_cast<size_t>(_nz)];
  _sheet = new real_t  [static_cast<size_t>(_nz * _nline)];

  for (k = 0; k < _nz; k++) _line[k] = _sheet + k*_nline;

  Veclib::zero (_nz * _nline, _sheet, 1);
  
  // The evaluation of bc here is removed.

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
    int_t                    i;
  
    cout << "# -- Field '" << F -> name() << "' Boundary Information:" << endl;
    if (!F -> _nbound) cout << "No BCs for this Field" << endl;
    for (i = 0; i < F -> _nbound; i++) BC[i] -> print();
  }
}


void Field::evaluateBoundaries (const int_t step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
// Note that for 3D this evaluation is done in Fourier-transformed space.
//
// This routine is only meant to be used for boundary conditions that must
// be re-evaluated at every step, such as high-order pressure BCs.
// evaluate boundaries other than control bcs.
// ---------------------------------------------------------------------------
{
  const int_t       np = Geometry::nP();
  vector<Boundary*> BC;
  real_t*           p;
  int_t             i, k;

  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));

  for (k = 0; k < _nz; k++)
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      if(!(*BC[i]->group()=='c' || *BC[i]->group()=='d') || _bsys -> field() == 'p'){
	BC[i] -> evaluate (k, step, p);
      }
}



void Field::evaluateControl (const int_t step, real_t* uci)
// ---------------------------------------------------------------------------
// evaluate control bcs.
// ---------------------------------------------------------------------------
{
  const int_t       np = Geometry::nP();
  vector<Boundary*> BC;
  real_t*           p;
  int_t             i, k;
  int_t             lengthcontrol= _bsys -> sizecontrolbc ();
  real_t            cowt, timenow, dt, adjointtime;
  real_t            timedecay = (1-exp(-timenow*timenow))*(1-exp(-adjointtime*adjointtime));
  real_t	    shift=1.83/2*1.75; //testing:::
  
  dt          = Femlib::value ("D_T");
  adjointtime = dt*(Femlib::ivalue ("N_STEP")-step);
  timenow     = dt*step;
    
  // -- Compute f(t)

  //To initialize uc set temporal dependence to 1.
  if (step == 0){
    timedecay = 1;
    cowt      = 1;
    if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
    else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA")); 
  }
  else{
    if (Femlib::ivalue  ("TIMEDECAY")==1) timedecay = exp(-3*(timenow-shift)*(timenow-shift));
    
    if (Femlib::ivalue  ("TIMEDECAY")==2) timedecay = 1; 
    
    if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
    else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));
    
    if (Femlib::ivalue ("controlbc_t_dependency") == 0){
      cowt=1;
    }
    else if (Femlib::ivalue ("controlbc_t_dependency") == 1) {
      cowt=cos( Femlib::value ("frequency_controlbc") *timenow);
    }      
  }
  
  cowt = cowt*timedecay;
  
  for (k = 0; k < _nz; k++){
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      {
	if(*BC[i]->group() == 'c'){
	  if (_nz==1){
	    Veclib::smul (np, cowt, uci, 1, p, 1);   
	    //Re{u_c+}=\hat{u}coswtf(t);Re{v_c}=\hat{v}coswtf(t);Im{w_c}=\hat{w}coswtf(t);
	  }
	  else if (k==0){
	    Veclib::smul  (np, cowt, uci, 1, p, 1);
	    //Veclib::svtvp (np, -siwt, uci+lengthcontrol, 1, p, 1, p, 1); 
	  }
	
	  else{
	    Veclib::smul  (np, cowt, uci, 1, p, 1);
	    //Veclib::svtvp (np, siwt, uci-lengthcontrol, 1, p, 1, p, 1); 
	  }
	  uci+=np;
	}	
      }
  }
}

void Field::fixControl ()
// ---------------------------------------------------------------------------
// Update boundary storage to accommodate intersection of control boundaries
// with other boundary types.
// ---------------------------------------------------------------------------
{
  const NumberSys*  N;
  const int_t       npert  = Geometry::nPert();
  int_t             sizebc = _bsys -> sizecontrolbc();
  int_t             np     = Geometry::nP();

  vector<Boundary*> BC;
  vector<real_t>    work;
  real_t            *p, *gvector;
  int_t             i, nglobal;

  if (npert == 2){
    BC = _bsys -> BCs  (0);
    N  = _bsys -> Nsys (0);
  }
  else{
    BC = _bsys -> BCs  (Femlib::ivalue("BETA")); 
    N  = _bsys -> Nsys (Femlib::ivalue("BETA"));
  }

  nglobal = N -> nGlobal();
  work.resize(nglobal);
  gvector = &work[0];
  Veclib::zero(nglobal, gvector, 1);

  // -- Write from _line storage to globally numbered vector then back
  //    to _line storage.                                           
  for(i = 0; i < _nz; i++){
    p = _line[i];
    this -> getEssential(p, gvector, BC, N);
    this -> giveEssential(gvector, p, BC, N);   
  }
}
 
void Field::getControl(real_t* uci)
// -----------------------------------------------------------------------
// Update external control boundary storage from current _line data.
// -----------------------------------------------------------------------
{  
  int_t sizebc = _bsys -> sizecontrolbc();

  int_t  j;
  real_t *p;

  // -- Write from _line storage to external control boundary storage. 
  for(j = 0; j < _nz; j++){
    p = _line[j];
    Veclib::copy(sizebc, p, 1, uci+j*sizebc, 1);
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
    int_t                    i;

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
    int_t                    i;

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
  const int_t      nglobal = N     -> nGlobal();
  const int_t*     btog    = N     -> btog();
  const int_t*     gid;
  int_t            i, k;
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
  int_t            i, k;
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


real_t Field::flux (const Field* C)
// ---------------------------------------------------------------------------
// Static member function.
//
// Compute gradient flux of field C on all "wall" group boundaries.
//
// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  const vector<Boundary*>& BC = C -> _bsys -> BCs (0);
  vector<real_t>           work(4 * Geometry::nP());
  real_t                   F = 0.0;
  int_t                    i;
  
  for (i = 0; i < C -> _nbound; i++)
    F += BC[i] -> scalarFlux ("wall", C -> _data, &work[0]);

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
  const int_t              nsurf = P -> _nbound;
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real_t>           work(Geometry::nP());
  int_t                    i;
  
  for (i = 0; i < nsurf; i++) {
    secF = BC[i] -> normTraction ("wall", P -> _data, &work[0]);
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
  const int_t              np      = Geometry::nP();
  const int_t              nbound  = U -> _nbound;
  const real_t             mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF, F = {0.0, 0.0, 0.0};
  vector<real_t>           work(4 * np);
  int_t                    i;

  for (i = 0; i < nbound; i++) {
    secF = UBC[i] -> tangTraction  ("wall", U->_data, V->_data, &work[0]);
    F.x        -= mu * secF.x;
    F.y        -= mu * secF.y;
    if (W) F.z -= mu * WBC[i] -> scalarFlux ("wall", W->_data, &work[0]);
  }

  return F;
}


void Field::normTractionV (real_t*      fx,
			   real_t*      fy,
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
  const int_t              np    = Geometry::nP();
  const int_t              nz    = Geometry::nZProc();
  const int_t              nsurf = P -> _nbound;
  Vector                   secF;
  vector<real_t>           work(np);
  real_t                   *p;
  int_t                    i, j;
  
  for (j = 0; j < nz; j++) {
    p = P -> _plane[j];
    for (i = 0; i < nsurf; i++) {
      secF = BC[i] -> normTraction ("wall", p, &work[0]);
      fx[j] += secF.x;
      fy[j] += secF.y;
    }
  }
}


void Field::tangTractionV (real_t*      fx,
			   real_t*      fy,
			   real_t*      fz,
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
  const int_t              np      = Geometry::nP();
  const int_t              nz      = Geometry::nZProc();
  const int_t              nbound = U -> _nbound;
  const real_t             mu      = Femlib::value ("RHO * KINVIS");
  Vector                   secF;
  vector<real_t>           work(4 * np);
  real_t                   *u, *v, *w;
  int_t                    i, j;

  for (j = 0; j < nz; j++) {
    u = U -> _plane[j];
    v = V -> _plane[j];
    w = (W) ? W -> _plane[j] : 0;
    for (i = 0; i < nbound; i++) {
      secF = UBC[i] -> tangTraction ("wall", u, v, &work[0]);
             fx[j] -= mu * secF.x;
             fy[j] -= mu * secF.y;
      if (W) fz[j] -= mu * WBC[i] -> scalarFlux ("wall", w, &work[0]);
    }
  }
}


Field& Field::solve (Domain* D,
		     AuxField*        f,
		     const MatrixSys* M,
		     const vector<AuxField*>& Baseflow,
		     bool forwards)
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
  const char  routine[] = "Field::solve";
  const int_t nel    = Geometry::nElmt();
  const int_t next   = Geometry::nExtElmt();
  const int_t npnp   = Geometry::nTotElmt();
  const int_t ntot   = Geometry::nPlane();
  int_t       i, k;

  const vector<Boundary*>& B       = M -> _BC;
  const NumberSys*         N       = M -> _NS;
  real_t                   lambda2 = M -> _HelmholtzConstant;
  real_t                   betak2  = M -> _FourierConstant;
  int_t                    nsolve  = M -> _nsolve;
  int_t                    nglobal = M -> _nglobal;
  int_t                    nzero   = nglobal - nsolve;

  for (k = 0; k < _nz; k++) {
    real_t* forcing = f -> _plane[k];
    real_t* unknown = _plane     [k];
    real_t* bc      = _line      [k];

    switch (M -> _method) {
    
    case DIRECT: {
      const real_t*  H     = const_cast<const real_t*> (M -> _H);
      const real_t** hii   = const_cast<const real_t**>(M -> _hii);
      const real_t** hbi   = const_cast<const real_t**>(M -> _hbi);
      const int_t*   b2g   = const_cast<const int_t*>  (N -> btog());
      int_t          nband = M -> _nband;
      vector<real_t> work (nglobal + 4*npnp);
      real_t         *RHS = &work[0], *tmp = RHS + nglobal;
      int_t          info;

      // -- Build RHS = - M f - H g + <h, w>.

      Veclib::zero (nglobal, RHS, 1);
    
      this -> getEssential (bc, RHS, B, N);
      this -> constrain    (forcing, lambda2,betak2,RHS,N,tmp);
      this -> buildRHS     (forcing, bc,RHS,0,hbi,nsolve,nzero,B,N,tmp);

      // -- Solve for unknown global-node values (if any).
      //      cout << "nsolve: " << nsolve << "nglobal: " << nglobal << endl;
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
      real_t         alpha, beta, dotp, epsb2, r2, rho1, rho2 = 0.0;
      vector<real_t> work (5*npts + 4*Geometry::nTotElmt());

      const int_t mode =
	(Geometry::problem() == Geometry::O2_2D ||
	 Geometry::problem() == Geometry::SO2_2D ) ? 0 : 1 ;

      real_t* r   = &work[0];
      real_t* p   = r + npts;
      real_t* q   = p + npts;
      real_t* x   = q + npts;
      real_t* z   = x + npts;
      real_t* wrk = z + npts;

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
      Veclib::copy (npts, x, 1, q, 1);
    
      this -> HelmholtzOperator (q, p, lambda2, betak2, mode, wrk, Baseflow, forwards);

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
	  Veclib::copy  (npts, z, 1, p, 1); // -- p = z.
	else {
	  beta = rho1 / rho2;	
	  Veclib::svtvp (npts, beta, p, 1, z, 1, p, 1); // -- p = z + beta p.
	}

	// -- Matrix-vector product.

	this -> HelmholtzOperator (p, q, lambda2, betak2, mode, wrk, Baseflow, forwards);
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
  
    default:
      message (routine, "called with a method that isn't implemented", ERROR);
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
//
// Input vector work should be 4*Geometry::nTotElmt() long.
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
  real_t            *u = work, *tmp = work + npnp;
  int_t             i;

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
			       real_t*       work   ,
			       const vector<AuxField*>&  Baseflow, 
			       bool          forwards) const
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
  const char       routine[] = "Field::HelmholtzOperator";
  const int_t      np        = Geometry::nP();
  const int_t      nel       = Geometry::nElmt();
  const int_t      npnp      = Geometry::nTotElmt();
  const int_t      next      = Geometry::nExtElmt();
  const int_t      nint      = Geometry::nIntElmt();
  const int_t      ntot      = Geometry::nPlane();
  const NumberSys* NS        = _bsys -> Nsys (mode * Femlib::ivalue ("BETA"));
  const int_t*     gid       = NS -> btog();
  const int_t      nglobal   = NS -> nGlobal() + Geometry::nInode();
  const real_t*    xint      = x + NS -> nGlobal();
  real_t*          yint      = y + NS -> nGlobal();
  real_t           *P = work, *tmp = work + npnp;
  int_t            i;
  Element*         E;

  Veclib::zero (nglobal, y, 1);

  // -- Add in contributions from mixed BCs while x & y are global vectors.

  if (_bsys -> mixBC()) {
    const vector<Boundary*>& BC   = _bsys -> BCs (0);
    const int_t*             bmap = NS    -> btog();
    int_t                    num  = 0;

    for (i = 0; i < _nbound; i++){
      if(*BC[i] -> group() =='c'){
	BC[i] -> switchK(num);
	num += 1;
      }
      BC[i] -> augmentOp (bmap, x, y);
    }
  }

  // -- Add in contributions from Toutflow BCs while x & y are global vectors.
  if (_bsys -> ToutflowBC()) {
    if (Geometry::nSlice()>1) message (routine, "Sorry, cannot support time-dependent base flow now", ERROR);
    const vector<Boundary*>& BC   = _bsys -> BCs (0);
    const int_t*           bmap = NS    -> btog();
    for (i = 0; i < _nbound; i++){
      BC[i] -> switchK (Baseflow[0]->plane(0), Baseflow[1]->plane(0), forwards);  
      BC[i] -> augmentOp (bmap, x, y);
    }
  }
  
  // -- Add in contributions from elemental Helmholtz operations.
  
  for (i = 0; i < nel; i++, gid += next, xint += nint, yint += nint) {
    E = _elmt[i];
    E -> global2local    (P, gid, x, xint);
    E -> HelmholtzOp     (lambda2, betak2, P, P, tmp);
    E -> local2globalSum (P, gid, y, yint);
  }
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
//
// Input vector work should be Geometry::nTotElmt() long.
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
  int_t                    i, boff;

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
  const int_t  nel  = Geometry::nElmt();
  const int_t  next = Geometry::nExtElmt();
  const int_t  nint = Geometry::nIntElmt();
  const int_t  npnp = Geometry::nTotElmt();
  const int_t* gid  = N -> btog();
  int_t        i;
  real_t*      internal = tgt + N -> nGlobal();

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
  const int_t   nel  = Geometry::nElmt();
  const int_t   next = Geometry::nExtElmt();
  const int_t   nint = Geometry::nIntElmt();
  const int_t   npnp = Geometry::nTotElmt();
  const int_t*  gid  = N -> btog();
  int_t         i;
  const real_t* internal = src + N -> nGlobal();

  for (i = 0; i < nel; i++, tgt += npnp, gid += next, internal += nint)
    _elmt[i] -> global2local (tgt, gid, src, internal);
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
  const int_t     np   = Geometry::nP();
  const int_t*    btog = N -> btog();
  const Boundary* B;
  int_t           i, boff;
  
  for (i = 0; i < _nbound; i++, src += np) {
    B    = bnd[i];
    boff = B -> bOff();
  
    B -> set (src, btog + boff, tgt);
  }
}


void Field::giveEssential (const real_t*            src,
			   real_t*                  tgt,
			   const vector<Boundary*>& bnd,
			   const NumberSys*         N  ) const
// ---------------------------------------------------------------------------
// On input, src is a globally numbered vector of boundary data for the
// current data plane. Gather current values of essential BCs into non-
// globally numbered vector tgt.
// ---------------------------------------------------------------------------
{
  const int_t     np   = Geometry::nP();
  const int_t*    btog = N -> btog();
  const Boundary* B;
  int_t           i, boff;
  
  for (i = 0; i < _nbound; i++, tgt += np) {
    B    = bnd[i];
    boff = B -> bOff();
  
    B -> get (src, btog + boff, tgt);
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
  const int_t  nel  = Geometry::nElmt();
  const int_t  next = Geometry::nExtElmt();
  const int_t  npnp = Geometry::nTotElmt();
  const int_t* emask = N -> emask();
  const int_t* bmask = N -> bmask();
  const int_t* btog  = N -> btog();
  int_t        i;

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
// --------------
{
  if (Geometry::problem() == Geometry::O2_2D ||
      Geometry::problem() == Geometry::SO2_2D ) return;

  const char     routine[] = "Field::couple";
  const int_t    nL        =  v -> _nline;
  vector<real_t> work (nL);
  real_t         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
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
// the "Fourier" constant in the Helmholtz equations is +/-1.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Field::modeConstant";

  if (Geometry::problem() == Geometry::O2_2D     ||
      Geometry::problem() == Geometry::SO2_2D    ||
      Geometry::system()  == Geometry::Cartesian || 
      name                ==         'c'         ||
      name                ==         'p'         ||
      name                ==         'u'          ) return beta * mode;

  if      (name == 'v') return beta * mode + 1.0;
  else if (name == 'w') return beta * mode - 1.0;
  else message (routine, "unrecognized Field name given", ERROR);

  return -1.0;			// -- Never happen.
}


real_t Field::normc(real_t* uci)
//get the norm of the control bc.
{  const int_t      np  = Geometry::nP();
  vector<Boundary*> BC;
  int_t             i, k;
  real_t            norm = 0.0;
  real_t            length= 0.0;
  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));

   for (k = 0; k < _nz; k++)    
    for ( i = 0; i < _nbound; i++)
	if(*BC[i]->group()=='c'||*BC[i]->group() =='d'){
	  norm   += BC[i] -> controlnorm (uci);
	  length += BC[i] -> controllength ();
	  uci +=np;
	  }
//scale the norm with the length of Controlbc
//norm=norm/length;
return norm;
}

real_t Field::normc_mixed(real_t* uci,real_t* lagi)
// Get the ith conponent of [uc,lag].
{ 
  const int_t       np  = Geometry::nP();
  vector<Boundary*> BC;
  int_t             i, k;
  real_t            norm = 0.0;
  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));
  
  for (k = 0; k < _nz; k++)    
    for ( i = 0; i < _nbound; i++)
      if(*BC[i]->group()=='c'){
	norm += BC[i] -> controlnorm_mixed (uci,lagi);
	uci  +=np;
	lagi +=np;
      }
  return norm;
}

void Field::add_adjoint(real_t* adjoint_gradient,
			const int_t step, 
			const int_t lengthcontrol)
// Add the gradient of the adjoint in the integration. direction gradient term.
{  
  int_t             i, j, k, id, side;
  const int_t       npnp = Geometry::nTotElmt();
  const int_t       np   = Geometry::nP();
  vector<Boundary*> BC;
  vector<real_t>    work;
  work.resize (2 * npnp);
  real_t*  ux    = new  real_t [static_cast<size_t>(npnp)];
  real_t*  uy    = new  real_t [static_cast<size_t>(npnp)];
  real_t*  gradx = new  real_t [static_cast<size_t>(np)];
  real_t*  grady = new  real_t [static_cast<size_t>(np)];
  real_t*  gradn = new  real_t [static_cast<size_t>(np)];
  real_t   timenow     = Femlib::value  ("D_T") * (Femlib::ivalue ("N_STEP")-step);
  real_t   adjointtime = Femlib::value  ("D_T")*step; 
  real_t   timedecay   = (1-exp(-timenow*timenow))*(1-exp(-adjointtime*adjointtime));
  real_t   shift       = 1.83/2*1.75; //tesing:::
 
  if (Femlib::ivalue  ("TIMEDECAY")==1) timedecay = exp(-3*(timenow-shift)*(timenow-shift)); //testing:::
  if (Femlib::ivalue  ("TIMEDECAY")==2) timedecay = 1; 
 
  real_t cowt, siwt;

 
  if ( Femlib::ivalue ("controlbc_t_dependency")==0)   {
    cowt=1; siwt=0;
  }
  
  else if ( Femlib::ivalue ("controlbc_t_dependency")==1) {
    cowt=cos( -Femlib::value ("frequency_controlbc") * timenow);
    siwt=sin( -Femlib::value ("frequency_controlbc") * timenow);
  }
 	  
  cowt=cowt*timedecay;
  siwt=siwt*timedecay; 	  
 
  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));

  for (k = 0; k < _nz; k++)      
    for (i = 0; i < _nbound; i++ )
      if(*BC[i]->group()=='c'||*BC[i]->group() =='d') {
	id = BC[i]->element()->ID();
	side = BC[i]->side();

	// Write field data for element into ux and uy.
	Veclib::copy (npnp, _plane[k]+id*npnp, 1, ux, 1);
	Veclib::copy (npnp, _plane[k]+id*npnp, 1, uy, 1);

	// Compute partial derivative of current field u WRT x and y and store in ux and uy respectively.
	BC[i] -> element() -> grad(ux, 0, &work[0]);
	BC[i] -> element() -> grad(0, uy, &work[0]);
	
	// Write partial derivatives on control boundaries to vectors gradx and grady
	switch (side) {
	case 0: for (j=0;j<np;j++)
	    {*(gradx+j)=*(ux+j);  *(grady+j)=*(uy+j);}  
	  break;
	case 1: for (j=0;j<np;j++)
	    {*(gradx+j)=*(ux+np*(j+1)-1);  *(grady+j)=*(uy+np*(j+1)-1);}
	  break;
	case 2: for (j=0;j<np;j++)
	    {*(gradx+j)=*(ux+npnp-1-j);   *(grady+j)=*(uy+npnp-1-j);}
	  break;
	case 3: for (j=0;j<np;j++)
	    {*(gradx+j)=*(ux+np*(np-1)-j*np); *(grady+j)=*(uy+np*(np-1)-j*np); }
	  break;
	default:     break;
	}
       
	//Compute n.gradu.
	BC[i] -> normal_gradient(gradx,grady,gradn);
	
	// Multiply by f(t) and add to what was previously written to adjoint_gradient at earlier time steps. Note -ive signs included as the velocity gradient term is subtracted. See eqn. 3.9 Mao, Blackburn and Sherwin 2012.
	if (_nz==1)
	  Veclib::svtvp (np, -Femlib::value("D_T * KINVIS")*cowt, gradn, 1, adjoint_gradient, 1, adjoint_gradient, 1);
	else if (k==0){
	  Veclib::svtvp (np, -Femlib::value("D_T * KINVIS")*cowt, gradn, 1, adjoint_gradient, 1, adjoint_gradient, 1);
	  Veclib::svtvp (np, -Femlib::value("D_T * KINVIS")*siwt, gradn, 1, adjoint_gradient+lengthcontrol, 1, adjoint_gradient+lengthcontrol, 1);
	}
	else{
	  Veclib::svtvp (np, -Femlib::value("D_T * KINVIS")*cowt, gradn, 1, adjoint_gradient, 1, adjoint_gradient, 1);
	  Veclib::svtvp (np,  Femlib::value("D_T * KINVIS")*siwt, gradn, 1, adjoint_gradient-lengthcontrol, 1, adjoint_gradient-lengthcontrol, 1); 
	}
       
	adjoint_gradient +=np;
      }

  delete [] ux;
  delete [] uy;
  delete [] gradx;
  delete [] grady;
  delete [] gradn;
}

void Field::add_adjoint_mixed(real_t*                  adjoint_gradient,
			      const int_t              step, 
			      const int_t              lengthcontrol,
			      const vector<AuxField*>& U)
// Compute integral term, g(u*, omega), in gradient of Lagrangian for
// mixed BC applied on adjoint "outflow" boundary.
{  
  const int_t npnp    = Geometry::nTotElmt();
  const int_t np      = Geometry::nP();

  real_t  timenow     = Femlib::value  ("D_T") * (Femlib::ivalue ("N_STEP")-step);
  real_t  adjointtime = Femlib::value  ("D_T")*step;
  real_t  timedecay   = (1-exp(-timenow*timenow))*(1-exp(-adjointtime*adjointtime));
  real_t  shift       = 1.83/2*1.75;

  vector<Boundary*> BC;
  real_t            *u, *Ux, *Uy, *p, *src, *nx, *ny, *tmp1, *tmp2;
  real_t            cowt;
  int_t             i, j, k, id, side, dOffset, dSkip;

  u  = new real_t[np];
  nx = new real_t[np];
  ny = new real_t[np];
  tmp1 = new real_t[np];
  tmp2 = new real_t[np];

  Ux = U[0] -> plane(0);
  Uy = U[1] -> plane(0);

  if (Femlib::ivalue  ("TIMEDECAY")==1) timedecay = exp(-3*(timenow-shift)*(timenow-shift));
  if (Femlib::ivalue  ("TIMEDECAY")==2) timedecay = 1;

  if ( Femlib::ivalue ("controlbc_t_dependency")==0){
    cowt=1;
  }
  else if ( Femlib::ivalue ("controlbc_t_dependency")==1){
    cowt=cos( Femlib::value ("frequency_controlbc") * timenow);
  }

  cowt=cowt*timedecay;

  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));

  for(k = 0; k < _nz; k++){
    p = _plane[k];
    for(i = 0; i < _nbound; i++){
      if(*BC[i] -> group() == 'c'){
	id      = BC[i] -> element() -> ID();
	side    = BC[i] -> side();
	dOffset = BC[i] -> dOff();
	dSkip   = BC[i] -> dSkip();
	src     = p + id*npnp;

	BC[i] -> element() -> sideGet(side, src, u);
	BC[i] -> controlnxny(nx, ny);

	Veclib::smul(np, -1, nx, 1, nx, 1);
	Veclib::smul(np, -1, ny, 1, ny, 1);
	Veclib::copy(np, Ux + dOffset, dSkip, tmp1, 1);
	Veclib::copy(np, Uy + dOffset, dSkip, tmp2, 1);
	Veclib::vvtvvtp(np, tmp1, 1, nx, 1, tmp2, 1, ny, 1, tmp1, 1);

	Veclib::vmul(np, u, 1, tmp1, 1, u, 1);

	Veclib::svtvp(np, cowt, u, 1, adjoint_gradient, 1, adjoint_gradient, 1);

	adjoint_gradient += np;
      }
    }
  }

  delete [] u;
  delete [] nx;
  delete [] ny;
  delete [] tmp1;
  delete [] tmp2;
}

void Field::add_adjoint_pressure(real_t* adjoint_px,
				 real_t* adjoint_py,
				 const int_t step,
				 const int_t lengthcontrol)
//add the gradient of the adjoint in the integration, pressure term.
{  
  const int_t npnp        = Geometry::nTotElmt();
  const int_t np          = Geometry::nP();
  real_t*     p_element   = new  real_t [static_cast<size_t>(npnp)];
  real_t*     p_edge      = new  real_t [static_cast<size_t>(np)]; 
  real_t*     nxp         = new  real_t [static_cast<size_t>(np)]; 
  real_t*     nyp         = new  real_t [static_cast<size_t>(np)]; 
  real_t      timenow     = Femlib::value  ("D_T") * (Femlib::ivalue ("N_STEP")-step);
  real_t      adjointtime = Femlib::value  ("D_T")*step; 
  real_t      timedecay   = (1-exp(-timenow*timenow))*(1-exp(-adjointtime*adjointtime));
  real_t      shift       = 1.83/2*1.75; //testing:::

  vector<Boundary*> BC;
  real_t            cowt, siwt;
  int_t             i, j, k, id, side;


  if (Femlib::ivalue  ("TIMEDECAY")==1) timedecay = exp(-3*(timenow-shift)*(timenow-shift)); //testing:::
  if (Femlib::ivalue  ("TIMEDECAY")==2) timedecay = 1; 
  
  if ( Femlib::ivalue ("controlbc_t_dependency")==0)   {
    cowt=1; siwt=0;
  }
  else if ( Femlib::ivalue ("controlbc_t_dependency")==1) {
    cowt=cos( -Femlib::value ("frequency_controlbc") * timenow);
    siwt=sin( -Femlib::value ("frequency_controlbc") * timenow);
  }
  
  cowt=cowt*timedecay;
  siwt=siwt*timedecay;
  
  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));
  
  for (k = 0; k < _nz; k++)      
    for (i = 0; i < _nbound; i++ )
      if(*BC[i]->group()=='c') {
	id =BC[i]->element()->ID();
	side=BC[i]->side();

	Veclib::copy (npnp, _plane[k]+id*npnp, 1, p_element , 1);
	
	switch (side) {
	case 0: for (j=0;j<np;j++)
	    *( p_edge+j)=*(p_element+j); 
	  break;
	case 1: for (j=0;j<np;j++)
	    *( p_edge+j)=*(p_element+np*(j+1)-1);  
	  break;
	case 2: for (j=0;j<np;j++)
	    *( p_edge+j)=*(p_element+npnp-1-j);   
	  break;
	case 3: for (j=0;j<np;j++)
	    *( p_edge+j)=*(p_element+np*(np-1)-j*np);
	  break;
	default:     break;
	}
	
	// Compute unit normal x p* on control boundary.
	BC[i] -> direction_pressure(p_edge,nxp,nyp);
	
	// Multiply by f(t) and add to what has previously been written to adjoint_gradient. Note that the pressure term only applies to the u and v components of the perturbation.
	if (_nz==1){
	  Veclib::svtvp (np, Femlib::value("D_T")*cowt, nxp, 1, adjoint_px, 1, adjoint_px, 1); 
	  Veclib::svtvp (np, Femlib::value("D_T")*cowt, nyp, 1, adjoint_py, 1, adjoint_py, 1); 
	}
	else if (k==0){
	  Veclib::svtvp (np, Femlib::value("D_T")*cowt, nxp, 1, adjoint_px, 1, adjoint_px, 1); 
	  Veclib::svtvp (np, Femlib::value("D_T")*cowt, nyp, 1, adjoint_py, 1, adjoint_py, 1); 
	  Veclib::svtvp (np, Femlib::value("D_T")*siwt, nxp, 1, adjoint_px+lengthcontrol, 1, adjoint_px+lengthcontrol, 1); 
	  Veclib::svtvp (np, Femlib::value("D_T")*siwt, nyp, 1, adjoint_py+lengthcontrol, 1, adjoint_py+lengthcontrol, 1);
	}
	else{
	  Veclib::svtvp (np, Femlib::value("D_T")*cowt, nxp, 1, adjoint_px, 1, adjoint_px, 1); 
	  Veclib::svtvp (np, Femlib::value("D_T")*cowt, nyp, 1, adjoint_py, 1, adjoint_py, 1); 
	  Veclib::svtvp (np,-Femlib::value("D_T")*siwt, nxp, 1, adjoint_px-lengthcontrol, 1, adjoint_px-lengthcontrol, 1); 
	  Veclib::svtvp (np,-Femlib::value("D_T")*siwt, nyp, 1, adjoint_py-lengthcontrol, 1, adjoint_py-lengthcontrol, 1); 
	}
	
	adjoint_px +=np;
	adjoint_py +=np;
      }
  
  delete []  p_element;
  delete []  p_edge ;
  delete []  nxp;
  delete []  nyp;
}

void Field::getInfK(real_t* ki, int_t step)
// -----------------------------------------------------------------
// Retrieve and store inflow velocity information for later use by 
// mixed BCs in the adjoint system.
// -----------------------------------------------------------------
{
  const int_t np   = Geometry::nP();
  const int_t nz   = Geometry::nZ();
  const int_t npnp = Geometry::nTotElmt();
  const int_t nc   = _bsys -> sizecontrolbc();

  vector<Boundary*> BC;
  vector<real_t>    work(2*npnp);
  real_t            *p, *src, *ux, *uy, *uxc, *uyc, *tmp;
  int_t             side, id, dOffset, dSkip, i, j, k, m;

  ux  = new real_t[npnp];
  uy  = new real_t[npnp];
  uxc = new real_t[np];
  uyc = new real_t[np];
  tmp = new real_t[np];

  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));

  for (k = 0; k < nz; k++){
    p = _plane[k];
    for (j = 0; j < _nbound; j++){
      if (*BC[j] -> group() == 'c'){
	side = BC[j] -> side();
	id   = BC[j] -> element() -> ID();
	src  = p + id*npnp;

	BC[j] -> element() -> sideGet(side, src, ki);
	/*	cout << "Inflow Velocity Data: " << endl;
	for (m = 0; m < np; m++){
	  cout << "Surface No: " << j << " "
	       << "Element ID: " << id << " "
	       << "Side: " << side << " "
	       << "Node: " << m << " "
	       << "ui: " << ki[m] << endl;
	       }*/
	

	Veclib::copy(npnp, src, 1, ux, 1);
	Veclib::copy(npnp, src, 1, uy, 1);

	BC[j] -> element() -> grad(ux, 0, &work[0]);
	BC[j] -> element() -> grad(0, uy, &work[0]);

	switch(side){
	case 0: for (i = 0; i < np; i++){
	    *(uxc+i) = *(ux+i);
	    *(uyc+i) = *(uy+i);
	  }
	  break;
	case 1: for (i = 0; i < np; i++){
	    *(uxc+i) = *(ux+np*(i+1)-1);
	    *(uyc+i) = *(uy+np*(i+1)-1);
	  }
	  break;
	case 2: for (i = 0; i < np; i++){
	    *(uxc+i) = *(ux+npnp-1-i);
	    *(uyc+i) = *(uy+npnp-1-i);
	  }
	  break;
	case 3: for (i = 0; i < np; i++){
	    *(uxc+i) = *(ux+np*(np-1)-i*np);
	    *(uyc+i) = *(uy+np*(np-1)-i*np);
	  }
	  break;
	default: break;
	}

	BC[j] -> normal_gradient(uxc, uyc, tmp);
	/*	cout << "Normal Gradient Data: " << endl;
	for (m = 0; m < np; m++){
	  cout << "Surface No: " << j << " "
	       << "Element ID: " << id << " "
	       << "Side: " << side << " "
	       << "Node: " << m << " "
	       << "ui: " << tmp[m] << endl;
	       }*/
      
	// Effectively Veclib::vdiv with check for divisor = 0.
	for(i = 0; i < np; i++){
	  if(ki[i] != 0){
	    ki[i] = tmp[i]/ki[i];
	  }
	  else {
	    ki[i] = 0;
	  }
	}

	//Veclib::smul(np, -1, ki, 1, ki, 1);

	//Veclib::vdiv(np, tmp, 1, ki, 1, ki, 1);

	// Print routine for testing.
	/*	cout << "K Data: " << endl;
	for (m = 0; m < np; m++){
	  cout << "Surface No: " << j << " "
	       << "Element ID: " << id << " "
	       << "Side: " << side << " "
	       << "Node: " << m << " "
	       << "ui: " << ki[m] << endl;
	       }*/

	//	Veclib::fill(np, 1.0, ki, 1); // Remove after testing.
	
	ki += np;
      }
    }
  }

  delete [] ux;
  delete [] uy;
  delete [] uxc;
  delete [] uyc;
  delete [] tmp;
}

void Field::updateK(real_t* ki)
// -------------------------------------------------------------------
// Write K data in domain storage to condition object storage.
// -------------------------------------------------------------------
{
  const int_t nc   = _bsys -> sizecontrolbc();

  vector<Boundary*> BC;
  real_t*           K;
  int_t             i;

  if (Geometry::nPert() == 2) BC = _bsys -> BCs (0);
  else                        BC = _bsys -> BCs (Femlib::ivalue ("BETA"));
  
  for (i = 0; i < _nbound; i++){
    if (*BC[i] -> group() == 'c'){
      K = BC[i] -> bcond() -> kfunc();
      Veclib::copy(nc, ki, 1, K, 1);
      break;
    }
  }
}

void Field::infBoundaryData(bool forwards, int_t step)
// --------------------------------------------------------------------
// Function to output fluid velocity component and normal gradient on 
// the control boundary.
// --------------------------------------------------------------------
{
  const int_t np   = Geometry::nP();
  const int_t nz   = Geometry::nZ();
  const int_t npnp = Geometry::nTotElmt();
  const int_t nc   = _bsys -> sizecontrolbc();

  vector<Boundary*> BC;
  vector<real_t>    work(2*npnp);
  real_t            *uinf, *guinf, *tmp, *p, *src, *x, *y,
                    *ux, *uy, *uxc, *uyc, *K, *kay, *Ku, *res;
  int_t             side, id, i, j, k, m, num;
  char              forward[StrMax], adjoint[StrMax];
  char              *n;

  n = new char[2];
  n[0] = _name;
  n[1] = '\0';

  sprintf(forward, "forward.dat.%d.", step);
  strcat(forward, n);
  sprintf(adjoint, "adjoint.dat.%d.", step);
  strcat(adjoint, n);

  m   = 0;
  num = 0;

  x     = new real_t[nc];
  y     = new real_t[nc];
  K     = new real_t[nc];
  Ku    = new real_t[nc];
  uinf  = new real_t[nc];
  guinf = new real_t[nc];
  res   = new real_t[nc];
  ux    = new real_t[npnp];
  uy    = new real_t[npnp];
  uxc   = new real_t[np];
  uyc   = new real_t[np];

  if (Geometry::nPert() == 2) BC = _bsys -> BCs(0);
  else                        BC = _bsys -> BCs(Femlib::ivalue("BETA"));

  for(k = 0; k < nz ; k++){
    p = _plane[k];
    for(i = 0; i < _nbound; i++){
      if(*BC[i] -> group() == 'c'){
	BC[i] -> switchK(num);
	side = BC[i] -> side();
	id   = BC[i] -> element() -> ID();
	src  = p + id*npnp;
	kay  = BC[i] -> bcond() -> getK();
	
	Veclib::copy(npnp, src, 1, ux, 1);
	Veclib::copy(npnp, src, 1, uy, 1);

        BC[i] -> element() -> grad(ux, 0, &work[0]);
        BC[i] -> element() -> grad(0, uy, &work[0]);

        switch(side){
        case 0: for (j = 0; j < np; j++){
            *(uxc+j) = *(ux+j);
            *(uyc+j) = *(uy+j);
          }
          break;
        case 1: for (j = 0; j < np; j++){
            *(uxc+j) = *(ux+np*(j+1)-1);
            *(uyc+j) = *(uy+np*(j+1)-1);
          }
          break;
        case 2: for (j = 0; j < np; j++){
            *(uxc+j) = *(ux+npnp-1-j);
            *(uyc+j) = *(uy+npnp-1-j);
          }
          break;
        case 3: for (j = 0; j < np; j++){
            *(uxc+j) = *(ux+np*(np-1)-j*np);
            *(uyc+j) = *(uy+np*(np-1)-j*np);
          }
          break;
        default: break;
        }

	BC[i] -> controlbcmesh(x+m*np, y+m*np);
	BC[i] -> element() -> sideGet(side, src, uinf+m*np);
        BC[i] -> normal_gradient(uxc, uyc, guinf+m*np);

	Veclib::copy(np, kay, 1, K+m*np, 1);
	Veclib::vmul(np, K+m*np, 1, uinf+m*np, 1, Ku+m*np, 1);

	Veclib::vadd(np, guinf+m*np, 1, Ku+m*np, 1, res+m*np, 1);
	Veclib::vdiv(np, res+m*np, 1, Ku+m*np, 1, res+m*np, 1);
	Veclib::smul(np, 100.0, res+m*np, 1, res+m*np, 1);

	num += 1;
	m += 1;
      }
    }
  }
  if(forwards){
    ofstream file(forward); 
    file << "x" << setw(15)
	 << "y" << setw(15)
	 << "u" << setw(15)
	 << "grad u" << endl;
    for(i = 0; i < nc; i++){
      file << x[i]     << setw(15)
	   << y[i]     << setw(15)
	   << uinf[i]  << setw(15)
	   << guinf[i] << endl;
    }
    file.close();
  }
  else{
    ofstream file1(adjoint);
    file1 << "Mixed BC Data: Grad u* + K x u* = 0" << endl;
    file1 << "x"       << setw(15)
	  << "y"       << setw(15)
	  << "u*"      << setw(15)
	  << "Grad u*" << setw(15)
	  << "K"       << setw(15)
	  << "K x u*"  << setw(15)
	  << "% Error"   << endl;
    for(i = 0; i < nc; i++){
      file1 << x[i]        << setw(15)
	    << y[i]        << setw(15)
	    << uinf[i]     << setw(15)
	    << guinf[i]    << setw(15)
	    << K[i]        << setw(15)
	    << Ku[i]       << setw(15)
	    << abs(res[i]) <<endl;
    }
    file1.close();
  }
  
  delete [] x;
  delete [] y;
  delete [] K;
  delete [] Ku;
  delete [] uinf;
  delete [] guinf;
  delete [] res;
  delete [] ux;
  delete [] uy;
  delete [] uxc;
  delete [] uyc;
  delete [] n;
}
