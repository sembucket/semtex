///////////////////////////////////////////////////////////////////////////////
// field.C: derived from AuxField, Field adds boundary conditions, global
// numbering, and the ability to solve Helmholtz problems.
//
// HELMHOLTZ PROBLEMS
// ------------------
// Solve routines provide solution to the discrete form of the Helmholtz eqn
//                      
//                       div grad u - \lambda^2 u = f,
//
// on domain \Omega, subject to essential BCs u = g on \Gamma_g and natural
// BCs \partial u / \partial n = h on \Gamma_h, where the boundary \Gamma
// of \Omega is the union of (non-overlapping) \Gamma_g and \Gamma_h
// and n is the unit outward normal vector on \Gamma.  \lambda^2 is
// called the Helmholtz constant below.
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
//                   K u + \lambda^2 M u = - M f + <h, w>
//
// or
//
//                         H u = - M f + <h, w>
//
// where K, M and H are respectively (assembled) "stiffness", "mass" and
// Helmholtz matrices.
//
// Some complications arise from dealing with essential boundary conditions,
// since typically the elemental matrices K^e, M^e which are assembled to
// form K and M do not account for the boundary requirements on w.  There
// are a number of ways of dealing with this issue: one approach is to
// partition H as it is formed (here F = -M f + <h, w>):
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
// nodal values (essential BCs), and solve instead the constrained problem
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
// Here n_global is the number of nodes that receive global node numbers,
// typically those on the mesh edges.  N_solve is the number of these
// nodes that have values that must be solved for, i.e. n_global minus the
// number of global nodes situated on essential-type boundaries.
//
// An alternative approach (USED HERE) to the constrained problem is to let
//
//                      u = v + g     (v = 0 on \Gamma_g)
//
// and solve instead
//
//                      H v = - M f - H g + <h, w>
//
// (where only the partition Hp is needed for the matrix H in the LHS), then
// afterwards compose u = v + g to get the full solution.  The advantage
// to this method is that the constraint partition Hc does not need
// to be assembled or stored.  The operations M f and H g can be performed
// on an element-by-element basis.
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
// BOUNDARY CONDITIONS
// -------------------
// (See also BCmgr.C.)  Boundary conditions may be arbitrarily applied to
// the various fields, EXCEPT in the case of cylindrical coordinates, see
// below.
//
// CYLINDRICAL COORDINATES
// -----------------------
// Special conditions apply in cylindrical geometries in the case where
// the axis is included in the domain.  In cylindical geometries it is
// necessary that the kind of boundary condition on the v & w components
// of velocity are identical on boundaries which don't lie on the axis.
// This is required so that those two fields, when coupled, have BCs of
// the same type away from the axis (value/function sub-types are equivalent).
//
// Summary of axial BCs for cylindrical coordinates:
//         +------------+------------------------------------------+
//         | (Coupled~) |       Axis BC for Fourier mode index     |
//         |  Variable  |         0             1           2...   |
//         +------------+------------------------------------------+
//         |   u, p, c  |   du/dr = 0       u   = 0      u  = 0    |
//         |      v~    |     v~  = 0       v~  = 0      v~ = 0    |
//         |      w~    |     w~  = 0    dw~/dr = 0      w~ = 0    |
//         +------------+------------------------------------------+
//
// In order to deal with this modal dependence on BCs, three levels of
// Boundary pointers are maintained for cylindrical Fields, corresponding
// to Fourier modes 0, 1, 2 (and higher), even if there are no axial BCs
// to be applied.
//
// The names of the numbering schemes that will be used for cylindrical
// 3D problems with axial BCs are (note case sensitivity):
//         +------------+------------------------------------------+
//         | (Coupled~) | Numbering scheme for Fourier mode index  |
//         |  Variable  |         0             1           2...   |
//         +------------+------------------------------------------+
//         |      u     |         u             U           U      |
//         |      v~    |         v             v           v      |
//         |      w~    |         w             W           w      |
//         |      p     |         p             P           P      |
//         |      c     |         c             C           C      |
//         +------------+------------------------------------------+
//
///////////////////////////////////////////////////////////////////////////////

static char 
RCSid[] = "$Id$";

#include <Sem.h>


Field::Field (FEML&                feml ,
	      BCmgr&               bcmgr,
	      vector<Element*>&    Elts ,
	      const char           name ,
	      const NumberSystem** Nmbr ) :

	      AuxField            (Elts ,
				   name ),
	      n_line              (0    )
// ---------------------------------------------------------------------------
// Create storage for a new Field from scratch.
//
// After constructing AuxField, use SURFACES information in FEML file to
// count number of surfaces on which BCs are applied, get Boundary Condition
// pointers from BCmgr, and allocate value storage areas for BCs.
//
// Example SURFACE description (NB: <P> -- periodic -- surfaces don't get BCs):
//
// <SURFACES NUMBER=6>
// #       tag     elmt    face    type
//         1       1       1       <P>     3       3       </P>
//         2       2       1       <P>     4       3       </P>
//         3       2       2       <B>     o       </B>
//         4       4       2       <B>     o       </B>
//         5       3       4       <B>     v       </B>
//         6       1       4       <B>     v       </B>
// </SURFACES>
//
// Use of cylindrical coordinates is flagged by the Geometry class variable.
// In the case where the number of space dimensions is also 3, the number
// of boundary frames and numbering systems is set to 3, for the 0th, 1st
// and 2nd (and higher) modes, irrespective of the number of Fourier modes
// actually used.  Internal variable n_bmodes flags if 1 or 3 sets are kept.
// The higher-mode boundaries only differ from the zero-mode ones on axis.
// Routines above this one are responsible for setting up the ordering of
// the required numbering systems, but the boundary conditions are set up here.
// ---------------------------------------------------------------------------
{
  char      routine[] = "Field::Field", err[StrMax], tag[StrMax];
  char      group, nextc;
  int       i, k, t, elmt, side;
  const int Nsurf = feml.attribute ("SURFACES", "NUMBER");
  const int nZ    =  Geometry::nZ();

  n_bmodes = (Geometry::system() == Geometry::Cylindrical &&
	      Geometry::nDim()   == 3)  ?  3 : 1;

  // -- Install NumberSystems.

  Nsys = new const NumberSystem* [n_bmodes];
  for (i = 0; i < n_bmodes; i++) Nsys[i] = Nmbr[i];

  // -- Construct temporary lists of Conditions by scanning SURFACE info,
  //    which has already been checked for consistency by Mesh constructor.

  List<char> grupID;
  List<int>  elmtID;
  List<int>  sideID;

  for (i = 0; i < Nsurf; i++) {

    while ((feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    // -- Get element and side number information, tag.

    feml.stream() >> t >> elmt >> side >> tag;
    elmt--;
    side--;

    if (strcmp (tag, "<B>") == 0) {
      
      feml.stream() >> group;
      
      grupID.add (group);
      elmtID.add (elmt );
      sideID.add (side );

      feml.stream() >> tag;
      if (strcmp (tag, "</B>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <B> with %s", t, tag);
	message (routine, err, ERROR);
      }
    } else
      feml.stream().ignore (StrMax, '\n');
  }

  // -- Construct vector of Boundary pointers using information just read in.
 
  Element*           E;
  Boundary*          B;
  ListIterator<char> gID (grupID);
  ListIterator<int>  eID (elmtID);
  ListIterator<int>  sID (sideID);
  
  n_bound     = grupID.length();
  boundary    = new Boundary** [n_bmodes];
  boundary[0] = new Boundary*  [n_bmodes * n_bound];

  // -- Zero-mode boundaries are used by both Cartesian & cylindrical forms.

  for (i = 0; i < n_bound; i++, gID.next(), eID.next(), sID.next()) {
    group = gID.current();
    side  = sID.current();
    elmt  = eID.current();
    E     = Elmt[elmt];
    boundary[0][i] = new Boundary (i, n_line,
				   bcmgr.descriptor (group),
				   bcmgr.retrieve   (group, field_name),
				   E, side);
    n_line += E -> nKnot();
  }

  // -- Higher mode boundaries get constructed only if 3D cylindrical.

  for (k = 1; k < n_bmodes; k++) {
    boundary[k] = boundary[0] + k * n_bound;
    n_line      = 0;
    for (gID.reset(), eID.reset(), sID.reset(), i = 0;
	 i < n_bound;
	 i++, gID.next(), eID.next(), sID.next()) {
      group = gID.current();      
      side  = sID.current();
      elmt  = eID.current();
      E     = Elmt[elmt];
      boundary[k][i] = (strstr (bcmgr.descriptor (group), "axis")) ? 
	new Boundary (i, n_line, // -- Create a new Boundary structure.
	  bcmgr.descriptor (group),
	  bcmgr.retrieve   (group, field_name, k),
	  E, side) :
	boundary[0][i];		// -- Alias an old one.
      n_line += E -> nKnot();
    }
  }

  // -- Allocate storage for boundary data.

  if (n_line & 1) n_line++;	// -- Round up to aid Fourier transform.
  
  line  = new real* [nZ];
  sheet = new real  [nZ * n_line];

  for (k = 0; k < nZ; k++) line[k] = sheet + k * n_line;

  Veclib::zero (nZ * n_line, sheet, 1);

  // -- Set values for boundary data (in physical space).

  const real dz = Femlib::value ("TWOPI / BETA") / nZ;

  for (k = 0; k < nZ; k++) {
    Femlib::value ("z", k * dz);
    for (i = 0; i < n_bound; i++) {
      B = boundary[0][i];
      t = B -> vOff();
      B -> evaluate (k, 0, line[k] + t);
    }
  }

  // -- Fourier transform boundary data.

  bTransform (+1);
}


void Field::bTransform (const int sign)
// ---------------------------------------------------------------------------
// Compute forward or backward 1D-DFT of boundary value storage areas.
//
// Normalization is carried out on forward transform, so that the zeroth
// mode's real data are the average over the homogeneous direction of the
// physical space values.  See also comments for AuxField::transform.
// ---------------------------------------------------------------------------
{
  Femlib::DFTr (sheet, Geometry::nZ(), n_line, sign);
}


void Field::printBoundaries (const Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// (Debugging) Utility to print information contained in a Boundary list.
// ---------------------------------------------------------------------------
{
  int i;
  
  cout
    << "# -- Field '"          
      << F -> name()
	<< "' Boundary Information:"
	  << endl;

  if (!F -> n_bound) cout << "No BCs for this Field" << endl;

  for (i = 0; i < F -> n_bound; i++) F -> boundary[0][i] -> print();
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
  register int       i, j, k, m, voff;
  const int          nZ = Geometry::nZ();
  register Boundary* B;

  for (k = 0; k < nZ; k++) {
    m = k >> 1;    
    j = min (m, n_bmodes - 1);
    for (i = 0; i < n_bound; i++) {
      B    = boundary[j][i];
      voff = B -> vOff();

      B -> evaluate (k, step, line[k] + voff);
    }
  }
}


void Field::evaluateM0Boundaries (const int step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind, but only for Mode 0.
// ---------------------------------------------------------------------------
{
  register int       i, voff;
  const int          nZ = Geometry::nZ();
  register Boundary* B;

  for (i = 0; i < n_bound; i++) {
    B    = boundary[0][i];
    voff = B -> vOff();

    B -> evaluate (0, step, line[0] + voff);
  }
}


void Field::addToM0Boundaries (const real  val,
			       const char* grp)
// ---------------------------------------------------------------------------
// Add val to zeroth Fourier mode's bc storage area on BC group "grp".
// ---------------------------------------------------------------------------
{
  register int       i, voff;
  register Boundary* B;

  for (i = 0; i < n_bound; i++) {
    B    = boundary[0][i];
    voff = B -> vOff();

    B -> addForGroup (grp, val, sheet + voff);
  }
}


Field& Field::smooth (AuxField* slave)
// ---------------------------------------------------------------------------
// Smooth slave field along element boundaries using *this, with
// mass-average smoothing.
//
// If slave == 0, smooth this -> data.
// ---------------------------------------------------------------------------
{
  register int      j, k, boff, doff;
  register Element* E;
  const int         nglobal = Nsys[0] -> nGlobal();
  const real*       imass   = Nsys[0] -> imass();
  const int*        btog    = Nsys[0] -> btog();
  const int         nE      = Geometry::nElmt();
  const int         nZ      = Geometry::nZ();
  vector<real>      work (nglobal);
  real              *src, *dssum = work();

  for (k = 0; k < nZ; k++) {

    Veclib::zero (nglobal, dssum, 1);
    src = (slave) ? slave -> plane[k] : plane[k];

    for (j = 0; j < nE; j++) {
      E    = Elmt[j];
      boff = E -> bOff();
      doff = E -> dOff();

      E -> bndryDsSum (btog + boff, src + doff, dssum);
    }

    Veclib::vmul (nglobal, dssum, 1, imass, 1, dssum, 1);

    for (j = 0; j < nE; j++) {
      E    = Elmt[j];
      boff = E -> bOff();
      doff = E -> dOff();

      E -> bndryInsert (btog + boff, dssum, src + doff);
    }
  }

  return *this;
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
  register int i;
  Vector       secF, F = {0.0, 0.0, 0.0};
  vector<real> work(Geometry::nP());
  
  for (i = 0; i < P -> n_bound; i++) {
    secF = P -> boundary[0][i] -> normalTraction ("wall", P -> data, work());
    F.x += secF.x;
    F.y += secF.y;
  }

  return F;
}


Vector Field::tangentTraction (const Field* U,
			       const Field* V)
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
//  T_ij = viscous stress tensor
//                          dU_i    dU_j
//       = RHO * KINVIS * ( ----  + ---- ) .
//                          dx_j    dx_i
//
// This only has to be done on the zero (mean) Fourier mode.
// ---------------------------------------------------------------------------
{
  register int i;
  const real   mu = Femlib::value ("RHO * KINVIS");
  Vector       secF, F = {0.0, 0.0, 0.0};
  vector<real> work(2 * Geometry::nP());
  real         *ddx, *ddy;

  ddx = work();
  ddy = ddx + Geometry::nP();

  for (i = 0; i < U -> n_bound; i++) {
    secF = U -> boundary[0][i] -> tangentTraction 
      ("wall", U -> data, V -> data, mu, ddx, ddy);
    F.x += secF.x;
    F.y += secF.y;
  }

  return F;
}


Field& Field::solve (AuxField*                f  ,
		     const ModalMatrixSystem* MMS)
// ---------------------------------------------------------------------------
// Carry out DIRECT solution of this Field using f as forcing; f must have
// the same structure as the Field to be solved.
//
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
//
// The RHS vector is constructed with length of the number of element-edge
// nodes in the problem (n_gid).  The first n_solve values contain forcing
// terms for the free (not essential BC) nodes in the problem, derived from
// the forcing field and the natural BCs "h", while the remaining values
// get loaded from essential BC values, "g".
//
// Forcing field f's data area is overwritten/destroyed during processing.
// ---------------------------------------------------------------------------
{
  register int        j, k, m, doff, boff;
  int                 info, nzero, nsolve, nband;
  const int*          b2g;
  const int           nglobal = Nsys[0] -> nGlobal();
  const int           nZ      = Geometry::nZ();
  const int           nE      = Geometry::nElmt();
  const MatrixSystem* M;
  const NumberSystem* N;
  const Boundary**    B;
  Element*            E;
  const real          *H, **hbi, **hii;
  vector<real>        work (nglobal);
  real                *RHS = work(), *forcing, *unknown, *bc;
  real                betak2, lambda2;

  for (k = 0; k < nZ; k++) {	// -- Loop over planes of data.
    
    if (k == 1) continue;	// -- Nyquist plane will always be set to zero.

    // -- Select Fourier mode, set local pointers and variables.

    m       = k >> 1;
    j       = min (m, n_bmodes - 1);

    B       = (const Boundary**) boundary[j];
    N       = Nsys[j];
    nband   = N -> nBand();
    b2g     = N -> btog ();
    
    M       = (*MMS) [m];
    H       = (const real*)  M -> H;
    hii     = (const real**) M -> hii;
    hbi     = (const real**) M -> hbi;
    nsolve  = M -> nsolve;
    lambda2 = M -> HelmholtzConstant;
    betak2  = M -> FourierConstant;

    nzero   = nglobal - nsolve;
    forcing = f -> plane[k];
    unknown = plane[k];
    bc      = line[k];

    // -- Build RHS = - M f - H g + <h, w>.

    Veclib::zero (nglobal, RHS, 1);
    getEssential (bc, RHS, B, N);
    constrain    (forcing, lambda2, betak2, RHS, N);
    buildRHS     (forcing, bc, RHS, 0, hbi, nsolve, nzero,B, N);

    // -- Solve for unknown global-node values (if any).

    if (nsolve)
      Lapack::pbtrs ("U", nsolve, nband - 1, 1, H, nband, RHS, nglobal, info);

    // -- Carry out Schur-complement solution for element-internal nodes.

    for (j = 0; j < nE; j++) {
      E    = Elmt[j];
      doff = E -> dOff();
      boff = E -> bOff();

      E -> g2eSC (RHS, b2g+boff, forcing+doff, unknown+doff, hbi[j], hii[j]);
    }

    // -- Scatter-gather essential BC values into plane.

    Veclib::zero (nglobal, RHS, 1);
    getEssential (bc, RHS, B,   N);
    setEssential (RHS, unknown, N);
  }

  return *this;
}


Field& Field::solve (AuxField*  f      ,
		     const real lambda2)
// ---------------------------------------------------------------------------
// Carry out ITERATIVE Conjugate Gradient solution of this Field with
// forcing from f.  Input value of L2 is the zero-mode Helmholtz constant.
//
// Problem for solution is the discrete version of the
// weak form of
//                                          2
//                       div grad u - lambda  u = f.
//
// Forcing field F's data area is overwritten/destroyed during processing.
//
// In this version, all vectors are ordered with globally-numbered (element-
// boundary) nodes first, followed by all element-internal nodes.
//
// The notation follows that used in Fig 2.5 of
//   Barrett et al., "Templates for the Solution of Linear Sytems", netlib.
// ---------------------------------------------------------------------------
{
  char         routine[] = "Field::solve";
  register int i, j, k, m;
  int          singular, nsolve, nzero;
  real         rho1, rho2, alpha, beta, r2, epsb2, betak2, dotp;
  real         *forcing, *unknown, *bc;
  const int    nZ      = Geometry::nZ();
  const int    nglobal = Nsys[0] -> nGlobal();
  const int    StepMax = (int) Femlib::value ("STEP_MAX");
  const int    npts    = nglobal + Geometry::nInode();
  const real   betaZ   = Femlib::value ("BETA");
  const real   FTINY   = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  
  const NumberSystem* N;
  const Boundary**    B;

  // -- Allocate storage.
  
  vector<real> workspace (6 * npts + 4 * Geometry::nTotElmt());

  real* r   = workspace();
  real* p   = r  + npts;
  real* q   = p  + npts;
  real* x   = q  + npts;
  real* z   = x  + npts;
  real* PC  = z  + npts;
  real* wrk = PC + npts;

  for (k = 0; k < nZ; k++) {	// -- Loop over planes of data.

    if (k == 1) continue;	// -- Nyquist planes remain zero always.

    // -- Select Fourier mode, set local pointers and variables.

    m        = k >> 1;
    j        = min (m, n_bmodes - 1);
    betak2   = sqr (Field::modeConstant (field_name, m) * betaZ);

    B        = (const Boundary**) boundary[j];
    N        = Nsys[j];
    nsolve   = N -> nSolve();
    singular = nglobal == nsolve && fabs (lambda2 + betak2) < FTINY;
    nsolve   = (singular) ? nsolve - 1 : nsolve;
    nzero    = nglobal - nsolve;

    forcing = f -> plane[k];
    unknown = plane[k];
    bc      = line[k];
    
    // -- Build diagonal preconditioner.

    jacobi (lambda2, betak2, PC, N);

    // -- f <-- - M f - H g, then (r = ) b = - M f - H g + <h, w>.

    Veclib::zero (nglobal, x, 1);
    getEssential (bc, x, B, N);  
    constrain    (forcing, lambda2, betak2, x, N);
    buildRHS     (forcing, bc, r, r + nglobal, 0, nsolve, nzero, B, N);

    epsb2  = Femlib::value ("TOL_REL") * sqrt (Blas::dot (npts, r, 1, r, 1))
           + Femlib::value ("TOL_ABS");
    epsb2 *= epsb2;

    // -- Build globally-numbered x from element store.

    local2global (unknown, x, N);
  
    // -- Compute first residual using initial guess: r = b - Ax.

    Veclib::zero (nzero, x + nsolve, 1);   
    Veclib::copy (npts,  x, 1, q, 1);

    HelmholtzOperator (q, p, lambda2, betak2, wrk, N);

    Veclib::zero (nzero, p + nsolve, 1);
    Veclib::vsub (npts, r, 1, p, 1, r, 1);

    r2 = Blas::dot (npts, r, 1, r, 1);

    // -- PCG iteration.

    i = 0;
    while (r2 > epsb2 && ++i < StepMax) {

      // -- Preconditioner.

      Veclib::vmul (npts, PC, 1, r, 1, z, 1);

      rho1 = Blas::dot (npts, r, 1, z, 1);

      // -- Update search direction.

      if (i == 1)
	Veclib::copy  (npts,             z, 1, p, 1); // -- p = z.
      else {
	beta = rho1 / rho2;	
	Veclib::svtvp (npts, beta, p, 1, z, 1, p, 1); // -- p = z + beta p.
      }

      // -- Matrix-vector product.

      HelmholtzOperator (p, q, lambda2, betak2, wrk, N);
      Veclib::zero (nzero, q + nsolve, 1);

      // -- Move in conjugate direction.

      dotp = Blas::dot (npts, p, 1, q, 1);
      dotp = (dotp > FTINY) ? dotp : FTINY;
			
      alpha = rho1 / dotp;
      Blas::axpy (npts,  alpha, p, 1, x, 1); // -- x += alpha p.
      Blas::axpy (npts, -alpha, q, 1, r, 1); // -- r -= alpha q.

      rho2 = rho1;
      r2   = Blas::dot (npts, r, 1, r, 1);
    }
  
    if (i == StepMax) message (routine, "step limit exceeded", WARNING);
  
    // -- Unpack converged vector x, impose current essential BCs.

    global2local (x, unknown, N);

    getEssential (bc, x, B,   N);
    setEssential (x, unknown, N);
  
    if ((int) Femlib::value ("VERBOSE")) {
      char s[StrMax];
      sprintf (s, ":%3d iterations, field '%c'", i, field_name);
      message (routine, s, REMARK);
    }
  }

  return *this;
}


void Field::jacobi (const real          lambda2,
		    const real          betak2 ,
		    real*               PC     ,
		    const NumberSystem* N      ) const
// ---------------------------------------------------------------------------
// Build diagonal (point Jacobi) preconditioner, PC, length npts.
//
// PC is arranged with global nodes first, followed by element-internal values.
// ---------------------------------------------------------------------------
{
  register int      i, next, nint;
  register Element* E;
  const int         nE   = Geometry::nElmt();
  const int         npts = N -> nGlobal() + Geometry::nInode();
  const int*        btog = N -> btog();
  real*             PCi  = PC + N -> nGlobal();
  vector<real>      work (2 * Geometry::nTotElmt() + Geometry::nP());
  real              *ed = work(), *ewrk = work() + Geometry::nTotElmt();
  
  Veclib::zero (npts, PC, 1);

  for (i = 0; i < nE; i++) {
    E    = Elmt[i];
    next = E -> nExt();
    nint = E -> nInt();

    E -> HelmholtzDg (lambda2, betak2, ed, ewrk);
    
    Veclib::scatr_sum (next, ed, btog, PC);
    Veclib::copy      (nint, ed + next, 1, PCi, 1);

    btog += next;
    PCi  += nint;
  }

#if 1
  Veclib::vrecp (npts, PC, 1, PC, 1);
#else  // -- Turn off preconditioner for testing.
  Veclib::fill (npts, 1.0, PC, 1);
#endif

}


void Field::constrain (real*               force  ,
		       const real          lambda2,
		       const real          betak2 ,
		       const real*         esstlbc,
		       const NumberSystem* N      ) const
// ---------------------------------------------------------------------------
// Replace f's data with constrained weak form of forcing: - M f - H g.
// On input, essential BC values (g) have been loaded into globally-numbered
// esstlbc, other values are zero.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      j, ntot, doff, boff;
  const int*        emask = N -> emask();
  const int*        bmask = N -> bmask();
  const int*        btog  = N -> btog();
  const int         nE    = Geometry::nElmt();
  real              *fdp, *in, *out, *ewrk;
  vector<real> work (4 * Geometry::nTotElmt());

  in   = work();
  out  = in  + Geometry::nTotElmt();
  ewrk = out + Geometry::nTotElmt();

  // -- Manufacture -(M f + H g).

  for (j = 0; j < nE; j++) {
    E    = Elmt[j];
    ntot = E -> nTot();
    doff = E -> dOff();
    boff = E -> bOff();
    fdp  = force + doff;

    E -> weight (fdp);		// -- f <-- M f.

    if (emask[j]) {		// -- f <-- M f + H g.

      Veclib::zero     (ntot, in, 1);
      E -> g2e         (in, btog + boff, esstlbc, 0);
      E -> HelmholtzOp (lambda2, betak2, in, out, ewrk);
      Veclib::vadd     (ntot, fdp, 1, out, 1, fdp, 1);
    }

    Veclib::neg (ntot, fdp, 1);	// -- f <-- -(M f + H g).
  }
}


void Field::HelmholtzOperator (const real*         x      ,
			       real*               y      ,
			       const real          lambda2,
			       const real          betak2 ,
			       real*               work   ,
			       const NumberSystem* N      ) const
// ---------------------------------------------------------------------------
// Discrete 2D global Helmholtz operator which takes the vector x into
// vector y, including direct stiffness summation.
//
// Vector work must have length 4 * elmt_nt_max.
// ---------------------------------------------------------------------------
{
  register Element*    E;
  register int         j, nint, boff;
  const int            nglobal = N -> nGlobal();
  const int*           btog    = N -> btog();
  const int            nE      = Geometry::nElmt();
  real*                in      = work;
  real*                out     = in  + Geometry::nTotElmt();
  real*                ewk     = out + Geometry::nTotElmt();
  register const real* x_int   = x + nglobal;
  register       real* y_int   = y + nglobal;

  Veclib::zero (nglobal + Geometry::nInode(), y, 1);

  for (j = 0; j < nE; j++) {
    E    = Elmt[j];
    nint = E -> nInt();
    boff = E -> bOff();

    E -> g2e         (in,  btog + boff, x, x_int);
    E -> HelmholtzOp (lambda2, betak2, in, out, ewk);
    E -> e2gSum      (out, btog + boff, y, y_int);

    x_int += nint;
    y_int += nint;
  }
}


void Field::buildRHS (real*               force ,
		      const real*         bc    ,
		      real*               RHS   ,
		      real*               RHSint,
		      const real**        hbi   ,
		      const int           nsolve,
		      const int           nzero ,
		      const Boundary**    bnd   ,
		      const NumberSystem* N     ) const
// ---------------------------------------------------------------------------
// Build RHS for direct or iterative solution.
//
// Iterative solution is flagged by presence of RHSint, a pointer to
// element-internal node storage.  If RHS is zero, then hbi, a vector of
// pointers to element interior/exterior coupling matrices, must be 
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
  register Element*        E;
  register const Boundary* B;
  register int             j, boff, doff;
  const int                nglobal = N -> nGlobal();
  const int*               btog    = N -> btog();
  const int                nE      = Geometry::nElmt();

  if   (RHSint) Veclib::zero (nglobal + Geometry::nInode(), RHS, 1);
  else          Veclib::zero (nglobal,                      RHS, 1);

  // -- Add in contribution from forcing f = - M f - H g.

  for (j = 0; j < nE; j++) {
    E    = Elmt[j];
    boff = E -> bOff();
    doff = E -> dOff();

    if (RHSint) {
      E -> e2gSum   (force + doff, btog + boff, RHS, RHSint);
      RHSint += E -> nInt();

    } else
      E -> e2gSumSC (force + doff, btog + boff, RHS, hbi[j]);
  }
  
  // -- Add in <h, w>.

  for (j = 0; j < n_bound; j++) {
    B    = bnd[j];
    doff = B -> vOff();
    boff = B -> bOff();

    B -> sum (bc + doff, btog + boff, RHS);
  }

  // -- Zero any contribution that <h, w> made to essential BC nodes.

  Veclib::zero (nzero, RHS + nsolve, 1);
}


void Field::local2global (const real*         src,
			  real*               tgt,
			  const NumberSystem* N  ) const
// ---------------------------------------------------------------------------
// Load a plane of data (src) into globally-numbered tgt, with element-
// boundary values in the first N -> nGlobal() places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      j, boff, doff;
  const int*        btog = N -> btog();
  const int         nE   = Geometry::nElmt();
  register real*    internal = tgt + N -> nGlobal();

  for (j = 0; j < nE; j++) {
    E    = Elmt[j];
    boff = E -> bOff();
    doff = E -> dOff();

    E -> e2g (src + doff, btog + boff, tgt, internal);
    internal += E -> nInt();
  }
}


void Field::global2local (const real*         src,
			  real*               tgt,
			  const NumberSystem* N  )
// ---------------------------------------------------------------------------
// Load a plane of data (tgt) from src, which has globally-numbered (element-
// boundary) values in the first nglobal places, followed by element-internal
// locations in emap ordering.
// ---------------------------------------------------------------------------
{
  register Element*    E;
  register int         j, boff, doff;
  const int*           btog     = N -> btog();
  const int            nE       = Geometry::nElmt();
  register const real* internal = src + N -> nGlobal();

  for (j = 0; j < nE; j++) {
    E    = Elmt[j];
    boff = E -> bOff();
    doff = E -> dOff();

    E -> g2e (tgt + doff, btog + boff, src, internal);
    internal += E -> nInt();
  }
}


void Field::getEssential (const real*         src,
			  real*               tgt,
			  const Boundary**    bnd,
			  const NumberSystem* N  ) const
// ---------------------------------------------------------------------------
// On input, src contains a line of BC values for the current data plane.
// Scatter current values of essential BCs into globally-numbered tgt.
//
// The construction of the essential BCs has to account for cases where
// element corners may touch the domain boundary but do not have an edge
// along a boundary.  This is done by working with a globally-numbered vector.
// ---------------------------------------------------------------------------
{
  register int             i, boff, voff;
  const int*               btog = N -> btog();
  register const Boundary* B;
  
  for (i = 0; i < n_bound; i++) {
    B = bnd[i];
    boff = B -> bOff();
    voff = B -> vOff();
  
    B -> set (src + voff, btog + boff, tgt);
  }
}


void Field::setEssential (const real*         src,
			  real*               tgt,
			  const NumberSystem* N  )
// ---------------------------------------------------------------------------
// Gather globally-numbered src into essential BC nodes of current data plane.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      k, boff, doff;
  const int*        emask = N -> emask();
  const int*        bmask = N -> bmask();
  const int*        btog  = N -> btog();
  const int         nE    = Geometry::nElmt();

  for (k = 0; k < nE; k++) {
    if (emask[k]) {
      E    = Elmt[k];
      boff = Elmt[k] -> bOff();
      doff = Elmt[k] -> dOff();
    
      E -> bndryMask (bmask + boff, tgt + doff, src, btog + boff);
    }
  }
}


void Field::coupleBCs (Field*    v  ,
		       Field*    w  ,
		       const int dir)
// ---------------------------------------------------------------------------
// Couples/uncouple boundary condition values for the radial and azimuthal
// velocity fields in cylindrical coordinates, depending on indicated
// direction.  This action is required due to the coupling in the viscous
// terms of the N--S equations in cylindrical coords.
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
// Since there is no coupling for the viscous terms in the 2D equation,
// do nothing for the zeroth Fourier mode.
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim() < 3) return;

  char         routine[] = "Field::couple";
  register int Re, Im, k;
  const int    nZ    = Geometry::nZ(),
               nL    = v -> n_line,
               nMode = nZ >> 1;
  vector<real> work (nL);
  real         *Vr, *Vi, *Wr, *Wi, *tp = work();
  
  if (dir == 1) {

    for (k = 1; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> line[Re];
      Vi = v -> line[Im];
      Wr = w -> line[Re];
      Wi = w -> line[Im];

      Veclib::copy (nL, Vr, 1, tp, 1);
      Veclib::vsub (nL, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nL, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nL, Wr, 1, tp, 1);
      Veclib::copy (nL, Wi, 1, Wr, 1);
      Veclib::vsub (nL, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nL, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == -1) {

    for (k = 1; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> line[Re];
      Vi = v -> line[Im];
      Wr = w -> line[Re];
      Wi = w -> line[Im];

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


int Field::modeConstant (const char name,
			 const int  mode)
// ---------------------------------------------------------------------------
// For cylindrical coordinates & 3D, the radial and azimuthal fields are
// coupled before solution of the viscous step.  This means that the Fourier
// constant used for solution may vary from that which applies to the axial
// component.  For Field v~, k -> k + 1 while for w~, k -> k - 1.
// For the uncoupled Fields v, w solved for the zeroth Fourier mode, the
// "Fourier" constant in the Helmholtz equations is 1.
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim()    <          3          ||
      Geometry::system() == Geometry::Cartesian || 
      name               ==         'c'         ||
      name               ==         'p'         ||
      name               ==         'u'          ) return mode;

  if      (name == 'v') return  mode + 1;
  else if (name == 'w') return (mode == 0) ? 1 : mode - 1;
  else message ("Field::modeConstant", "unrecognized Field name given", ERROR);

  return INT_MAX;
}
