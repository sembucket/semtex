///////////////////////////////////////////////////////////////////////////////
// field.C: derived from AuxField, Field adds boundary conditions, global
// numbering, and the ability to solve Helmholtz problems.
//
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
// since typically the elemental matrices A^e, B^e which are assembled to
// form A and B do not account for the boundary requirements on w.  There
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
// Partition out the parts of the matrix corresponding to the known nodal
// values (essential BCs), and solve instead the constrained problem
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
//                      u = v + g (v zero on \Gamma_g)
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
///////////////////////////////////////////////////////////////////////////////

static char 
RCSid[] = "$Id$";

#include <Sem.h>


Field::Field (FEML&               feml ,
	      BCmgr&              bcmgr,
	      vector<Element*>&   Elts ,
	      const NumberSystem* Nmbr,
	      const int           nz   ,
	      const char          name ) :

	      AuxField           (Elts ,
				  nz   ,
				  name ),
	      Nsys               (Nmbr ),
	      n_line             (0    )
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
// ---------------------------------------------------------------------------
{
  char              routine[] = "Field::Field", err[StrMax], tag[StrMax];
  char              group, nextc;
  int               i, k, t, elt, sid;
  const int         Ns = feml.attribute ("SURFACES", "NUMBER");
  List<Condition*>  tmpBC;
  List<int>         elmtID;
  List<int>         sideID;
  List<const char*> descript;

  // -- Construct a temporary list of Conditions by scanning SURFACE info,
  //    which has already been checked for consistency by Mesh constructor.

  for (i = 0; i < Ns; i++) {

    while ((nextc = feml.stream().peek()) == '#') // -- Skip comments.
      feml.stream().ignore (StrMax, '\n');

    // -- Get element and side number information, tag.

    feml.stream() >> t >> elt >> sid >> tag;
    elt--;
    sid--;

    if (strcmp (tag, "<B>") == 0) {
      
      feml.stream() >> group;
      
      tmpBC   .add (bcmgr.retrieve   (group, field_name));
      descript.add (bcmgr.descriptor (group));
      elmtID  .add (elt);
      sideID  .add (sid);

      feml.stream() >> tag;
      if (strcmp (tag, "</B>") != 0) {
	sprintf (err, "Surface %1d: couldn't close tag <B> with %s", t, tag);
	message (routine, err, ERROR);
      }
    } else
      feml.stream().ignore (StrMax, '\n');
  }

  // -- Construct vector of Boundary pointers using information just read in.
 
  Element*                  E;
  Boundary*                 B;
  ListIterator<Condition*>  cPt (tmpBC);
  ListIterator<int>         eID (elmtID);
  ListIterator<int>         sID (sideID);
  ListIterator<const char*> des (descript);
  
  n_bound  = tmpBC.length();
  boundary = new Boundary* [n_bound];

  for (i = 0; i < n_bound; i++,cPt.next(),eID.next(),sID.next(),des.next()) {
    sid         = sID.current();
    elt         = eID.current();
    E           = Elmt[elt];
    boundary[i] = new Boundary (i, n_line, des.current(), cPt.current(),E,sid);
    n_line     += E -> nKnot();
  }

  // -- Allocate storage for boundary data.

  line  = new real* [n_z];
  sheet = new real  [n_z * n_line];

  for (k = 0; k < n_z; k++) line[k] = sheet + k * n_line;

  Veclib::zero (n_z * n_line, sheet, 1);

  // -- Set values for boundary data (in physical space).

  const real dz = Femlib::value ("TWOPI/BETA") / nz;

  for (k = 0; k < n_z; k++) {
    Femlib::value ("z", k * dz);
    for (i = 0; i < n_bound; i++) {
      B = boundary[i];
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
  if (n_z < 2) return;

  register int i;
  const int    ntot = n_z * n_line;
  vector<real> work (3 * n_z + 15);
  real*        tmp  = work();
  real*        Wtab = tmp + n_z;
  real*        ptr  = sheet;

  Femlib::rffti (n_z, Wtab);

  switch (sign) {

  case 1:
    for (i = 0; i < n_line; i++, ptr++) {
      Veclib::copy  (n_z, ptr, n_line, tmp, 1);
      Femlib::rfftf (n_z, tmp, Wtab);
      Veclib::copy  (n_z - 2, tmp + 1, 1, ptr + 2 * n_line, n_line);
      ptr[0]      = tmp[0];
      ptr[n_line] = 0.0;
    }
    Blas::scal (ntot, 1.0 / n_z, sheet, 1);
    break;

  case -1:
    for (i = 0; i < n_line; i++, ptr++) {
      tmp[n_z - 1] = 0.0;
      tmp[0]       = ptr[0];
      Veclib::copy  (n_z - 2, ptr + 2 * n_line, n_line, tmp + 1, 1);
      Femlib::rfftb (n_z, tmp, Wtab);
      Veclib::copy  (n_z, tmp, 1, ptr, n_line);
    }
    break;

  default:
    message ("Field::bTransform", "illegal direction flag", ERROR);
    break;
  }
}


void Field::printBoundaries (const Field* F)
// ---------------------------------------------------------------------------
// Static member function.
//
// (Debugging) Utility to print information contained in a Boundary list.
// ---------------------------------------------------------------------------
{
  int  i;
  
  cout
    << "# -- Field '"          
      << F -> name()
	<< "' Boundary Information:"
	  << endl;

  if (!F -> n_bound) cout << "No BCs for this Field" << endl;

  for (i = 0; i < F -> n_bound; i++) F -> boundary[i] -> print();
}


void Field::evaluateBoundaries (const int step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
// Note that for 3D this evaluation is done in Fourier-transformed space.
// ---------------------------------------------------------------------------
{
  register int       i, k, voff;
  register Boundary* B;

  for (k = 0; k < n_z; k++)
    for (i = 0; i < n_bound; i++) {
      B    = boundary[i];
      voff = B -> vOff();

      B -> evaluate (k, step, line[k] + voff);
    }
}


void Field::addToBoundaries (const real  val,
			     const char* grp)
// ---------------------------------------------------------------------------
// Add val to zeroth Fourier mode's bc storage area on BC group "grp".
// ---------------------------------------------------------------------------
{
  register int       i, voff;
  register Boundary* B;

  for (i = 0; i < n_bound; i++) {
    B    = boundary[i];
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
  const int         nglobal = Nsys -> nGlobal();
  const real*       imass   = Nsys -> imass();
  const int*        btog    = Nsys -> btog();
  const int         N       = Nsys -> nEl();
  vector<real>      work (nglobal);
  real              *src, *dssum = work();


  for (k = 0; k < n_z; k++) {

    Veclib::zero (nglobal, dssum, 1);
    src = (slave) ? slave -> plane[k] : plane[k];

    for (j = 0; j < N; j++) {
      E    = Elmt[j];
      boff = E -> bOff();
      doff = E -> dOff();

      E -> bndryDsSum (btog + boff, src + doff, dssum);
    }

    Veclib::vmul (nglobal, dssum, 1, imass, 1, dssum, 1);

    for (j = 0; j < N; j++) {
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
  vector<real> work(P -> elmt_np_max);
  
  for (i = 0; i < P -> n_bound; i++) {
    secF = P -> boundary[i] -> normalTraction ("wall", P -> data, work());
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
  const real   mu = Femlib::value ("RHO*KINVIS");
  Vector       secF, F = {0.0, 0.0, 0.0};
  vector<real> work(2 * U -> elmt_np_max);
  real         *ddx, *ddy;

  ddx = work();
  ddy = ddx + U -> elmt_np_max;

  for (i = 0; i < U -> n_bound; i++) {
    secF = U -> boundary[i] -> tangentTraction ("wall", U -> data, V -> data,
						mu,     ddx,       ddy      );
    F.x += secF.x;
    F.y += secF.y;
  }

  return F;
}


Field& Field::solve (AuxField*                f  ,
		     const ModalMatrixSystem* MMS,
		     MatrixOperator           Opr)
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
  int                 j, k, info, nzero, nsolve, doff, boff;
  const int           nglobal = Nsys -> nGlobal();
  const int           nband   = Nsys -> nBand();
  const int           N       = Nsys -> nEl();
  const int*          b2g     = Nsys -> btog();
  const MatrixSystem* M;
  Element*            E;
  const real          *H, **hbi, **hii;
  vector<real>        work (nglobal);
  real                *RHS = work(), *forcing, *unknown, *bc;
  real                lambda2;

  for (k = 0; k < n_z; k++) {	// -- Loop over planes of data.
    
    if (k == 1) continue;	// -- Nyquist plane will always be set to zero.

    // -- Select Fourier mode, set local pointers and variables.

    j       = k >> 1;

    M       = (*MMS) [j];
    H       = (const real*)  M -> H;
    hii     = (const real**) M -> hii;
    hbi     = (const real**) M -> hbi;
    nsolve  = M -> nsolve;
    lambda2 = M -> HelmholtzConstant;

    nzero   = nglobal - nsolve;
    forcing = f -> plane[k];
    unknown = plane[k];
    bc      = line[k];

    // -- Build RHS = - M f - H g + <h, w>.

    Veclib::zero (nglobal, RHS, 1);
    getEssential (bc, RHS);
    constrain    (forcing, lambda2, RHS, Opr);
    buildRHS     (forcing, bc, RHS, 0, hbi, nsolve, nzero);

    // -- Solve for unknown global-node values (if any).

    if (nsolve)
      Lapack::pbtrs ("U", nsolve, nband - 1, 1, H, nband, RHS, nglobal, info);

    // -- Carry out Schur-complement solution for element-internal nodes.

    for (j = 0; j < N; j++) {
      E    = Elmt[j];
      doff = E -> dOff();
      boff = E -> bOff();

      E -> g2eSC (RHS, b2g+boff, forcing+doff, unknown+doff, hbi[j], hii[j]);
    }

    // -- Scatter-gather essential BC values into plane.

    Veclib::zero (nglobal, RHS, 1);
    getEssential (bc, RHS);
    setEssential (RHS, unknown);
  }

  return *this;
}


Field& Field::solve (AuxField*      f  ,
		     const real     L2 ,
		     MatrixOperator Opr)
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
// Operator points to an Element routine that applies the discrete Helmholtz
// to element-level storage to assemble contributions to storage ordered
// as above.  Default operator is for 2D Cartesian geometry.
//
// The notation follows that used in Fig 2.5 of
//   Barrett et al., "Templates for the Solution of Linear Sytems", netlib.
// ---------------------------------------------------------------------------
{
  char       routine[] = "Field::solve";
  int        i, k, singular, nsolve, nzero;
  real       rho1, rho2, alpha, beta, epsb;
  const real betaZ = Femlib::value ("BETA");
  const real FTINY = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  real       lambda2;
  real       *forcing, *unknown, *bc;

  const int StepMax = (int) Femlib::value ("STEP_MAX");
  const int nglobal = Nsys -> nGlobal();
  const int npts    = nglobal + n_elmt_inodes;

  // -- Allocate storage.
  
  vector<real> workspace (4 * npts + 4 * elmt_nt_max);

  real* r   = workspace();
  real* p   = r + npts;
  real* w   = p + npts;
  real* x   = w + npts;
  real* wrk = x + npts;

  for (k = 0; k < n_z; k++) {	// -- Loop over planes of data.

    if (k == 1) continue;	// -- Nyquist planes remain zero always.

    // -- Select Fourier mode, set local pointers and variables.

    lambda2  = L2 + sqr ((k >> 1) * betaZ);

    nsolve   = Nsys -> nSolve();
    singular = nglobal == Nsys -> nSolve() && fabs (lambda2) < FTINY;
    nsolve   = (singular) ? nsolve - 1 : nsolve;
    nzero    = nglobal - nsolve;

    forcing = f -> plane[k];
    unknown = plane[k];
    bc      = line[k];

    // -- f <-- - M f - H g, then (r = ) b = - M f - H g + <h, w>.

    Veclib::zero (nglobal, x, 1);
    getEssential (bc, x);  
    constrain    (forcing, lambda2, x, Opr);
    buildRHS     (forcing, bc, r, r + nglobal, 0, nsolve, nzero);

    epsb = Femlib::value ("TOL_REL") * sqrt (Blas::dot (npts, r, 1, r, 1));

    // -- Build globally-numbered x from element store.

    local2global (unknown, x);
  
    // -- Compute first residual using initial guess: r = b - Ax.
    
    Veclib::copy (npts, x, 1, w, 1);
    Veclib::zero (nzero, w + nsolve, 1);

    HelmholtzOperator (w, p, lambda2, wrk, Opr);

    Veclib::zero (nzero, p + nsolve, 1);
    Veclib::vsub (npts, r, 1, p, 1, r, 1);

    rho1 = Blas::dot (npts, r, 1, r, 1);

    // -- CG iteration.  Vector x will not evolve at essential BC nodes.

    for (i = 0; sqrt (rho1) > epsb && i < StepMax; i++) {

      if   (i) Veclib::svtvp (npts, beta, p, 1, r, 1, p, 1); //  p = r + beta p
      else     Veclib::copy  (npts,             r, 1, p, 1); //  p = r

      HelmholtzOperator (p, w, lambda2, wrk, Opr);
      Veclib::zero (nzero, w + nsolve, 1);

      alpha = rho1 / Blas::dot (npts, p, 1, w, 1);
      Blas::axpy (npts,  alpha, p, 1, x, 1); // -- x += alpha p
      Blas::axpy (npts, -alpha, w, 1, r, 1); // -- r -= alpha w

      rho2 = rho1;
      rho1 = Blas::dot (npts, r, 1, r, 1);
      beta = rho1 / rho2;
    }
  
    if (i == StepMax) message (routine, "step limit exceeded", WARNING);
  
    // -- Unpack converged vector x, impose current essential BCs.

    global2local (x, unknown);

    getEssential (bc, x);
    setEssential (x, unknown);
  
    if ((int) Femlib::value ("VERBOSE")) {
      char s[StrMax];
      sprintf (s, ":%3d CG iterations, field '%c'", i, field_name);
      message (routine, s, REMARK);
    }
  }

  return *this;
}


void Field::constrain (real*          force   ,
		       const real     lambda2 ,
		       const real*    esstlbc ,
		       MatrixOperator Operator) const
// ---------------------------------------------------------------------------
// Replace f's data with constrained weak form of forcing: - M f - H g.
// On input, essential BC values (g) have been loaded into globally-numbered
// esstlbc, other values are zero.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      j, ntot, doff, boff;
  const int*        emask = Nsys -> emask();
  const int*        bmask = Nsys -> bmask();
  const int*        btog  = Nsys -> btog();
  const int         N     = Nsys -> nEl();
  real              *fdp, *in, *out, *ewrk;
  vector<real> work (4 * elmt_nt_max);

  in   = work();
  out  = in  + elmt_nt_max;
  ewrk = out + elmt_nt_max;

  // -- Manufacture -(M f + H g).

  for (j = 0; j < N; j++) {
    E    = Elmt[j];
    ntot = E -> nTot();
    doff = E -> dOff();
    boff = E -> bOff();
    fdp  = force + doff;

    E -> weight (fdp);		// -- f <-- M f.

    if (emask[j]) {		// -- f <-- M f + H g.

      Veclib::zero    (ntot, in, 1);
      E ->  g2e       (in, btog + boff, esstlbc, 0);
     (E ->* Operator) (in, out, lambda2, ewrk);
      Veclib::vadd    (ntot, fdp, 1, out, 1, fdp, 1);
    }

    Veclib::neg (ntot, fdp, 1);	// -- f <-- -(M f + H g).
  }
}


void Field::HelmholtzOperator (const real*    x       ,
			       real*          y       ,
			       const real     lambda2 ,
			       real*          work    ,
			       MatrixOperator Op      ) const
// ---------------------------------------------------------------------------
// Discrete 2D global Helmholtz operator which takes the vector x into
// vector y, including direct stiffness summation.
//
// Vector work must have length 4 * elmt_nt_max.
// ---------------------------------------------------------------------------
{
  register Element*    E;
  register int         j, nint, boff;
  const int            nglobal = Nsys -> nGlobal();
  register const real* x_int   = x + nglobal;
  register       real* y_int   = y + nglobal;
  const int*           btog    = Nsys -> btog();
  const int            N       = Nsys -> nEl();
  real*                in      = work;
  real*                out     = in  + elmt_nt_max;
  real*                ewk     = out + elmt_nt_max;

  Veclib::zero (nglobal + n_elmt_inodes, y, 1);

  for (j = 0; j < N; j++) {
    E    = Elmt[j];
    nint = E -> nInt();
    boff = E -> bOff();

    E -> g2e    (in,  btog + boff, x, x_int);
   (E ->* Op)   (in,  out, lambda2, ewk);
    E -> e2gSum (out, btog + boff, y, y_int);

    x_int += nint;
    y_int += nint;
  }
}


void Field::buildRHS (real*        force ,
		      const real*  bc    ,
		      real*        RHS   ,
		      real*        RHSint,
		      const real** hbi   ,
		      const int    nsolve,
		      const int    nzero ) const
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
  register Element*  E;
  register Boundary* B;
  register int       j, boff, doff;
  const int          nglobal = Nsys -> nGlobal();
  const int*         btog    = Nsys -> btog();
  const int          N       = Nsys -> nEl();

  if   (RHSint) Veclib::zero (nglobal + n_elmt_inodes, RHS, 1);
  else          Veclib::zero (nglobal,                 RHS, 1);

  // -- Add in contribution from forcing f = - M f - H g.

  for (j = 0; j < N; j++) {
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
    B    = boundary[j];
    doff = B -> vOff();
    boff = B -> bOff();

    B -> sum (bc + doff, btog + boff, RHS);
  }

  // -- Zero any contribution that <h, w> made to essential BC nodes.

  Veclib::zero (nzero, RHS + nsolve, 1);
}


void Field::local2global (const real* src,
			  real*       tgt) const
// ---------------------------------------------------------------------------
// Load a plane of data (src) into globally-numbered tgt, with element-
// boundary values in the first N -> nGlobal() places, followed by
// element-internal locations in emap ordering.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      j, boff, doff;
  const int*        btog = Nsys -> btog();
  const int         N    = Nsys -> nEl();
  register real*    internal = tgt + Nsys -> nGlobal();

  for (j = 0; j < N; j++) {
    E    = Elmt[j];
    boff = E -> bOff();
    doff = E -> dOff();

    E -> e2g (src + doff, btog + boff, tgt, internal);
    internal += E -> nInt();
  }
}


void Field::global2local (const real* src,
			  real*       tgt)
// ---------------------------------------------------------------------------
// Load a plane of data (tgt) from src, which has globally-numbered (element-
// boundary) values in the first nglobal places, followed by element-internal
// locations in emap ordering.
// ---------------------------------------------------------------------------
{
  register Element*    E;
  register int         j, boff, doff;
  const int*           btog     = Nsys -> btog();
  const int            N        = Nsys -> nEl();
  register const real* internal = src + Nsys -> nGlobal();

  for (j = 0; j < N; j++) {
    E    = Elmt[j];
    boff = E -> bOff();
    doff = E -> dOff();

    E -> g2e (tgt + doff, btog + boff, src, internal);
    internal += E -> nInt();
  }
}


void Field::getEssential (const real* src,
			  real*       tgt) const
// ---------------------------------------------------------------------------
// On input, src contains a line of BC values for the current data plane.
// Scatter current values of essential BCs into globally-numbered tgt.
//
// The construction of the essential BCs has to account for cases where
// element corners may touch the domain boundary but do not have an edge
// along a boundary.  This is done by working with a globally-numbered vector.
// ---------------------------------------------------------------------------
{
  register Boundary* B;
  register int       i, boff, voff;
  const int*         btog = Nsys -> btog();
  
  for (i = 0; i < n_bound; i++) {
    B    = boundary[i];
    boff = B -> bOff();
    voff = B -> vOff();
  
    B -> set (src + voff, btog + boff, tgt);
  }
}


void Field::setEssential (const real* src,
			  real*       tgt)
// ---------------------------------------------------------------------------
// Gather globally-numbered src into essential BC nodes of current data plane.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      k, boff, doff;
  const int*        emask = Nsys -> emask();
  const int*        bmask = Nsys -> bmask();
  const int*        btog  = Nsys -> btog();
  const int         N     = Nsys -> nEl();

  for (k = 0; k < N; k++) {
    if (emask[k]) {
      E    = Elmt[k];
      boff = Elmt[k] -> bOff();
      doff = Elmt[k] -> dOff();
    
      E -> bndryMask (bmask + boff, tgt + doff, src, btog + boff);
    }
  }
}
