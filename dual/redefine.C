///////////////////////////////////////////////////////////////////////////////
// redefine.C: redefined and added class methods for the standard semtex
// distribution.  This must be linked before the standard class objects.
//
// Copyright (C) 2000 Hugh Blackburn
//
// For dual, data are NOT transformed to physical space for storage in
// restart files.  There are 3 planes of data: the first corresponds
// to mode 0 (which is real only), the second and third to the real an
// imaginary parts of the selected Fourier mode.
//
// NB: define "K_FUND" in the TOKENS section, which pulls in the right
// set of BCs on the axis for cylindrical geometries.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


// ===========================================================================
// From analysis.C (these are all redefinitions, with original prototypes):


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
// ---------------------------------------------------------------------------
  src (D)
{
  char str[StrMax];

  cout << setprecision (6);

  // -- Initialize averaging.

  if ((integer) Femlib::value ("AVERAGE")) {
    vector<AuxField*> extra (0);
    stats = new Statistics (D, extra);
  } else                              
    stats = 0;

  // -- Set up for output of modal energies every IO_CFL steps if 3D.

  if (Geometry::nDim() == 3) {
    strcat (strcpy (str, src -> name), ".mdl");
    ROOTONLY {
      mdl_strm.open (str, ios::out); 
      mdl_strm <<
"#         Time          Mode0          ModeC       ModeC.Re       ModeC.Im"
<< endl
	       <<
"# ------------------------------------------------------------------------" 
<< endl;
      mdl_strm.setf (ios::scientific, ios::floatfield);
      mdl_strm.precision (8);
    }
  }
}


void Analyser::analyse (AuxField** work)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const integer verbose = (integer) Femlib::value ("VERBOSE");
  const integer cflstep = (integer) Femlib::value ("IO_CFL");

  // -- Step-by-step updates.

  cout << "Step: " << src -> step << "  Time: " << src -> time << endl;

  // -- CFL, energy, divergence information.

  if (cflstep && !(src -> step % cflstep)) {
    modalEnergy ();
    estimateCFL ();
    divergence  (work);
  }

  // -- Periodic dumps and global information.
  
  const integer periodic = !(src->step %  (integer) Femlib::value("IO_HIS")) ||
                           !(src->step %  (integer) Femlib::value("IO_FLD"));
  const integer final    =   src->step == (integer) Femlib::value("N_STEP");
  const integer state    = periodic || final;

  if (state) {
     
    // -- Statistical analysis.

    if (stats) stats -> update (work);
  }

  // -- Field and statistical dumps.

  src -> dump ();
  if (stats) stats -> dump();
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const integer    DIM   = Geometry::nDim();
  register integer i;
  real             re, im, ek[4];

  Veclib::zero (4, ek, 1);

  for (i = 0; i < DIM; i++) {
    src -> u[i] -> mode_en (0, re, im);
    ek[0] += re;
    src -> u[i] -> mode_en (1, re, im);
    ek[1] += re + im;
    ek[2] += re;
    ek[3] += im;
  }

  mdl_strm << setw(10) << src -> time 
	   << setw(15) << ek[0]
	   << setw(15) << ek[1]
	   << setw(15) << ek[2]
	   << setw(15) << ek[3]
	   << endl;
}


// ===========================================================================
// From auxfield.C (there is an added function: AuxField::convolve,
// others match prototypes):


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.  Physical space.
//
// NB: function should not have z as a variable for dual: enforce z=0.
//
// This is a redefinition.
// ---------------------------------------------------------------------------
{
  const integer    nel = Geometry::nElmt();
  const integer    np2 = Geometry::nTotElmt();
  register integer i, k;
  real*            p;

  Femlib::value ("z", 0.0);

  for (k = 0; k < _nz; k++)
    for (p = _plane[k], i = 0; i < nel; i++, p += np2)
      _elmt[i] -> evaluate (function, p);
  
  return *this;
}


AuxField& AuxField::convolve (const AuxField& x, const AuxField& y)
// ---------------------------------------------------------------------------
// Set this AuxField equal to the convolution of a & b, all terms made
// in Fourier-transformed space.  This is for use with dual, and all
// AuxFields possess only two modes, three data planes.
//
//                  ^     ^  ^    ^ ^    ^ ^
//  Mode 0 result: x y  = x  y  + x y  + x y  
//                    0    -c c    0 0    c -c
//
//                  ^     ^ ^   ^ ^
//  Mode c result: x y  = x y + x y
//                    c    0 c   c 0
//
// This is a new method.
// ---------------------------------------------------------------------------
{
  const integer nP = Geometry::planeSize();

  Veclib::vmul    (nP, x._plane[0], 1, y._plane[0], 1, _plane[0], 1);
  Veclib::svvttvp (nP, 2.0, x._plane[1], 1, y._plane[1], 1,
		   _plane[0], 1, _plane[0], 1);
  Veclib::svvttvp (nP, 2.0, x._plane[2], 1, y._plane[2], 1,
		   _plane[0], 1, _plane[0], 1);

  Veclib::vvtvvtp (nP, x._plane[0], 1, y._plane[1], 1, 
		   x._plane[1], 1, y._plane[0], 1, _plane[1], 1);
  Veclib::vvtvvtp (nP, x._plane[0], 1, y._plane[2], 1, 
		   x._plane[2], 1, y._plane[0], 1, _plane[2], 1);

  return *this;
}


AuxField& AuxField::gradient (const integer dir)
// ---------------------------------------------------------------------------
// Operate on AuxField to produce the nominated index of the gradient.
// dir == 0 ==> gradient in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
//
// NB: modfied for dual, only 1 non-zero mode, nominated as "K_FUND"
// in session.
//
// This is a redefinition.
// ---------------------------------------------------------------------------
{
  const char       routine[] = "AuxField::gradient";
  const integer    nel  = Geometry::nElmt();
  const integer    np   = Geometry::nP();
  const integer    npnp = np  * np;
  const integer    ntot = nel * npnp;
  const integer    nP   = Geometry::planeSize();
  vector<real>     work;
  register real    *xr, *xs, *tmp;
  register integer i, k;
  const real       **DV, **DT;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  switch (dir) {

  case 0:
    work.setSize (2 * nP);
    xr = work();
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      tmp = _plane[k];

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (tmp, tmp, xr, xs, *DV, *DT, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, tmp += npnp)
	_elmt[i] -> gradX (xr, xs, tmp);
      
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 1:
    work.setSize (2 * nP);
    xr = work();
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      tmp = _plane[k];

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (tmp, tmp, xr, xs, *DV, *DT, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, tmp += npnp)
	_elmt[i] -> gradY (xr, xs, tmp);
      
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 2: {
    const real betak = Femlib::value ("BETA * K_FUND");

    work.setSize (nP);
    xr = work();

    Veclib::zero (nP, _data, 1);
    Veclib::copy (nP,                       _plane[1], 1, xr, 1);
    Veclib::smul (nP, -betak, _plane[2], 1, _plane[1], 1);
    Veclib::smul (nP,  betak, xr,        1, _plane[2], 1);
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }

  return *this;
}


void AuxField::gradient (const integer nP ,
			 real*         src,
			 const integer dir) const
// ---------------------------------------------------------------------------
// Use Field structure to perform gradient operations on data area
// src, according to nominated direction.  Input value nZ is the
// number of planes to operate on.  Input value nP is the size of src
// in the orthogonal direction: for dir == 0, 1, this should be the
// size of a data plane, but for dir == 2 it can be arbitrary,
// e.g. the size of a data plane or the size of a block of data which
// has planar/row structure.
//
// NB: modified for dual: nZ is always going to be 3 (nP can still be
// arbitrary).
//
// This is a new method.
// ---------------------------------------------------------------------------
{
  const char       routine[] = "AuxField::gradient";
  const integer    nel  = Geometry::nElmt();
  const integer    np   = Geometry::nP();
  const integer    npnp = np  * np;
  const integer    ntot = nel * npnp;
  register integer i, k;
  vector<real>     work;
  register real    *plane, *xr, *xs, *Re, *Im;
  const real       **DV, **DT;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);

  switch (dir) {

  case 0:
    work.setSize (2 * nP);
    xr = work();
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      plane = src + k * nP;

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (plane, plane, xr, xs, *DV, *DT, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, plane += npnp)
	_elmt[i] -> gradX (xr, xs, plane);
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 1:
    work.setSize (2 * nP);
    xr = work();
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      plane = src + k * nP;

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (plane, plane, xr, xs, *DV, *DT, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, plane += npnp)
	_elmt[i] -> gradY (xr, xs, plane);
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 2: {
    const real betak = Femlib::value ("BETA * K_FUND");

    work.setSize (nP);
    xr = work();

    Veclib::zero (nP, src, 1);

    Re = src + nP;
    Im = Re  + nP;
    Veclib::copy (nP,             Re, 1, xr, 1);
    Veclib::smul (nP, -betak,  Im, 1, Re, 1);
    Veclib::smul (nP,  betak,  xr, 1, Im, 1);
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }
}


real AuxField::mode_L2 (const integer mode) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
// Mode numbers run 0 -- n_z/2 - 1.
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const integer     nel  = Geometry::nElmt();
  const integer     npnp = Geometry::nTotElmt();
  register real     area = 0.0, Ek = 0.0, *Re, *Im;
  register integer  i;
  register Element* E;

  if (mode == 0) {
    Re = _plane[0];
    for (i = 0; i < nel; i++, Re += npnp) {
      E       = _elmt[i];
      area   += E -> area();
      Ek     += sqr (E -> norm_L2 (Re));
    }
  } else {
    Re = _plane[1];
    Im = _plane[2];
    for (i = 0; i < nel; i++, Re += npnp, Im += npnp) {
      E     = _elmt[i];
      area += E -> area();
      Ek   += sqr (E -> norm_L2 (Re));
      Ek   += sqr (E -> norm_L2 (Im));
    }
  }

  return Ek / (2.0 * area);
}


void AuxField::mode_en (const integer mode  ,
			real&         RePart,
			real&         ImPart) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
//
// NB: Modified for dual: only 2 modes exist, 0 & 1.  We return
// mode_L2 split into components produced by real and imaginary parts.
//
// New method.
// ---------------------------------------------------------------------------
{
  const integer     nel  = Geometry::nElmt();
  const integer     npnp = Geometry::nTotElmt();
  register real     area = 0.0, *Re, *Im;
  register integer  i;
  register Element* E;

  RePart = ImPart = 0.0;
  
  if (mode == 0) {
    Re = _plane[0];
    for (i = 0; i < nel; i++, Re += npnp) {
      E       = _elmt[i];
      area   += E -> area();
      RePart += sqr (E -> norm_L2 (Re));
    }
  } else {
    Re = _plane[1];
    Im = _plane[2];
    for (i = 0; i < nel; i++, Re += npnp, Im += npnp) {
      E       = _elmt[i];
      area   += E -> area();
      RePart += sqr (E -> norm_L2 (Re));
      ImPart += sqr (E -> norm_L2 (Im));
    }
  }

  RePart /= 2.0 * area;
  ImPart /= 2.0 * area;
}


void AuxField::couple (AuxField*     v  ,
		       AuxField*     w  ,
		       const integer dir)
// ---------------------------------------------------------------------------
// Couples/uncouple field data for the radial and azimuthal velocity
// fields in cylindrical coordinates, depending on indicated
// direction.  This action is required due to the coupling in the
// viscous terms of the N--S equations in cylindrical coords.
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
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Field::couple";
  const integer nP    =  Geometry::planeSize();
  vector<real>  work (nP);
  real          *Vr, *Vi, *Wr, *Wi, *tp = work();
  
  if (dir == FORWARD) {

    Vr = v -> _plane[1];
    Vi = v -> _plane[2];
    Wr = w -> _plane[1];
    Wi = w -> _plane[2];

    Veclib::copy (nP, Vr, 1, tp, 1);
    Veclib::vsub (nP, Vr, 1, Wi, 1, Vr, 1);
    Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);
    Veclib::copy (nP, Wr, 1, tp, 1);
    Veclib::copy (nP, Wi, 1, Wr, 1);
    Veclib::vsub (nP, Vi, 1, tp, 1, Wi, 1);
    Veclib::vadd (nP, Vi, 1, tp, 1, Vi, 1);

  } else if (dir == INVERSE) {

    Vr = v -> _plane[1];
    Vi = v -> _plane[2];
    Wr = w -> _plane[1];
    Wi = w -> _plane[2];
    
    Veclib::copy  (nP,      Vr, 1, tp, 1);
    Veclib::svvpt (nP, 0.5, Vr, 1, Wr, 1, Vr, 1);
    Veclib::svvmt (nP, 0.5, Wr, 1, tp, 1, Wr, 1);
    Veclib::copy  (nP,      Wi, 1, tp, 1);
    Veclib::copy  (nP,      Wr, 1, Wi, 1);
    Veclib::svvmt (nP, 0.5, Vi, 1, tp, 1, Wr, 1);
    Veclib::svvpt (nP, 0.5, Vi, 1, tp, 1, Vi, 1);

  } else
    message (routine, "unknown direction given", ERROR);
}


// ===========================================================================
// From domain.C, just a redefine I/O so there are no Fourier transforms.


void Domain::restart ()
// ---------------------------------------------------------------------------
// If a restart file "name".rst can be found, use it for input.  If
// this fails, initialize all Fields to zero ICs.
//
// Carry out forwards Fourier transformation, zero Nyquist data.
//
// Redefintion.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer nF = nField();
  char          restartfile[StrMax];
  
  ROOTONLY cout << "-- Initial condition       : ";
  ifstream file (strcat (strcpy (restartfile, name), ".rst"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << restartfile;
      cout.flush();
    }
    file >> *this;
    file.close();
  } else {
    ROOTONLY cout << "set to zero";
    for (i = 0; i < nF; i++) *u[i] = 0.0;
  }

  ROOTONLY cout << endl;
  
  Femlib::value ("t", time);
  step = 0;
}


void Domain::dump ()
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;
  ofstream output;
  
  ROOTONLY {
    const char    routine[] = "Domain::dump";
    char          dumpfl[StrMax], backup[StrMax], command[StrMax];
    const integer verbose   = (integer) Femlib::value ("VERBOSE");
    const integer chkpoint  = (integer) Femlib::value ("CHKPOINT");

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, name), ".fld");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, name), ".chk");
	if (!initial) {
	  strcat  (strcpy (backup, name), ".chk.bak");
	  sprintf (command, "mv ./%s ./%s", dumpfl, backup);
	  system  (command);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, name), ".fld");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }
    
    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }
  
  output << *this;
  output.close();
}


///////////////////////////////////////////////////////////////////////////////
// From field.C: these are all redefinitions.


Field::Field (BoundarySys*      B,
	      real*             M,
	      const integer     N,
	      vector<Element*>& E,
	      const char        C) :
// ---------------------------------------------------------------------------
// Create storage for a new Field from scratch.
//
// Redefinition.
// ---------------------------------------------------------------------------
  AuxField (M, N, E, C),
  _bsys    (B)
{
  const integer            np  = Geometry::nP();
  const integer            nzb = Geometry::basePlane();
  const vector<Boundary*>& BC  = _bsys -> BCs (0);
  register real*           p;
  register integer         i, k;

  // -- Allocate storage for boundary data.
  
  _nbound = _bsys -> nSurf();
  _nline  = _nbound * np;

  _line  = new real* [(size_t)  _nz];
  _sheet = new real  [(size_t) (_nz * _nline)];

  for (k = 0; k < _nz; k++) _line[k] = _sheet + k*_nline;

  Veclib::zero (_nz * _nline, _sheet, 1);

  // -- Set values for boundary data, 0th Fourier mode only, enforce z = 0.

  Femlib::value ("z", 0);
  for (p = _line[0], i = 0; i < _nbound; i++, p += np)
    BC[i] -> evaluate (0, 0, p);
}


void Field::evaluateBoundaries (const integer step)
// ---------------------------------------------------------------------------
// Traverse Boundaries and evaluate according to kind.
//
// This routine is only meant to be used for boundary conditions that must
// be re-evaluated at every step, such as high-order pressure BCs.
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const integer    np    = Geometry::nP();
  const integer    kfund = (integer) Femlib::value ("K_FUND");
  register integer i, k;
  real*            p;

  const vector<Boundary*>& BC0 = _bsys -> BCs (0);
  const vector<Boundary*>& BCk = _bsys -> BCs (kfund);

  for (p = _line[0], i = 0; i < _nbound; i++, p += np)
    BC0[i] -> evaluate (0, step, p);

  for (k = 1; k < 3; k++) {
    for (p = _line[k], i = 0; i < _nbound; i++, p += np)
      BCk[i] -> evaluate (k, step, p);
  }
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
// Redefinition.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Field::solve";
  const integer np    = Geometry::nP();
  const integer nel   = Geometry::nElmt();
  const integer next  = Geometry::nExtElmt();
  const integer npnp  = Geometry::nTotElmt();
  const integer ntot  = Geometry::nPlane();
  const integer kfund = (integer) Femlib::value ("K_FUND");
  integer       i, k, mode;

  for (k = 0; k < _nz; k++) {	// -- Loop over planes of data.
    
    // -- Select Fourier mode, set local pointers and variables.

    if   (k == 0) mode = 0;
    else          mode = kfund;

    const MatrixSys*         M       = (*MMS)[(k != 0)];
    const vector<Boundary*>& B       = M -> _BC;
    const NumberSys*         N       = M -> _NS;
    real                     lambda2 = M -> _HelmholtzConstant;
    real                     betak2  = M -> _FourierConstant;
    integer                  nsolve  = M -> _nsolve;
    integer                  nglobal = M -> _nglobal;
    integer                  nzero   = nglobal - nsolve;

    real*                    forcing = f -> _plane[k];
    real*                    unknown = _plane     [k];
    real*                    bc      = _line      [k];

    switch (M -> _method) {

    case DIRECT: {
      const real*    H       = (const real*)    M -> _H;
      const real**   hii     = (const real**)   M -> _hii;
      const real**   hbi     = (const real**)   M -> _hbi;
      const integer* b2g     = (const integer*) N -> btog();
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
      const integer StepMax = (integer) Femlib::value ("STEP_MAX");  
      const integer npts    = M -> _npts;
      real           alpha, beta, dotp, epsb2, r2, rho1, rho2;
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

      this -> HelmholtzOperator (q, p, lambda2, betak2, wrk, mode);

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

	this -> HelmholtzOperator (p, q, lambda2, betak2, wrk, mode);
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
    }
  }
  return *this;
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
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Field::couple";
  const integer nL = v -> _nline;
  vector<real>  work (nL);
  real          *Vr, *Vi, *Wr, *Wi, *tp = work();
  
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


//////////////////////////////////////////////////////////////////////////////
// From geometry.C: the only routine here is redefined!


void Geometry::set (const integer  NP,
		    const integer  NZ,
		    const integer  NE,
		    const CoordSys CS)
// ---------------------------------------------------------------------------
// Load values of static internal variables.
// ---------------------------------------------------------------------------
{
  static char routine[] = "Geometry::set", err[StrMax];

  pid   = (integer) Femlib::value ("I_PROC");
  nproc = (integer) Femlib::value ("N_PROC");

  np   = NP; nz = NZ; nel = NE; csys = CS;
  nzp  = nz;
  ndim = 3;

  if (nz != 3) {
    sprintf (err, "dual needs N_Z == 3 (%1d)", nz);
    message (routine, err, ERROR);
  }

  if (nproc != 1) {
    sprintf (err, "This is a serial code (1 process, not %1d)", nproc);
    message (routine, err, ERROR);
  }

  psize = nPlane();
}


//////////////////////////////////////////////////////////////////////////////
// From matrix.C, redefine constructor.


static List<MatrixSys*> MS;

ModalMatrixSys::ModalMatrixSys (const real              lambda2 ,
				const real              beta    ,
				const integer           baseMode,
				const integer           numModes,
				const vector<Element*>& Elmt    ,
				const BoundarySys*      Bsys    ,
				const SolverKind        method  )
// ---------------------------------------------------------------------------
// Generate or retrieve from internal database MS the vector of
// MatrixSys's which will be used to solve all the Fourier-mode
// discrete Helmholtz problems for the associated scalar Fields
// (called out by names).
//
// Input variables:
//   lambda2 : Helmholtz constant for the problem,	
//   beta    : Fourier length scale = TWOPI / Lz,
//   nmodes  : number of Fourier modes which will be solved,
//   Elmt    : vector of Element*'s used to make local Helmholtz matrices,
//   Bsys    : boundary system for this field.
//   method  : specify the kind of solver we want (Cholesky, PCG ...).
//
// NB: for dual, there will be only two modes, we select the mode
// numbers here.
// ---------------------------------------------------------------------------
{
  const char               name = Bsys -> field();
  const integer            kfund = (integer) Femlib::value ("K_FUND");
  integer                  found, index;
  ListIterator<MatrixSys*> m (MS);
  MatrixSys*               M;

  _fields = new char [strlen (Bsys -> Nsys (0) -> fields()) + 1];
  strcpy (_fields, Bsys -> Nsys (0) -> fields());
  _Msys.setSize (2);

  if (method == DIRECT) {
    cout << "-- Installing matrices for field '" << name << "' [";
    cout.flush();
    Femlib::synchronize();
  }

  for (index = 0; index < 2; index++) {
    const integer    mode      = index * kfund;
    const NumberSys* N         = Bsys -> Nsys (mode);
    const real       betak2    = sqr (Field::modeConstant (name, mode, beta));

    for (found = 0, m.reset(); !found && m.more(); m.next()) {
      M     = m.current();
      found = M -> match (lambda2, betak2, N, method);
    }
    if (found) {
      _Msys[index] = M;
      if (method == DIRECT) { cout << '.'; cout.flush(); }
    } else {
      _Msys[index] =
	new MatrixSys (lambda2, betak2, mode, Elmt, Bsys, method);
      MS.add (_Msys[index]);
      if (method == DIRECT) { cout << '*'; cout.flush(); }
    }
  }

  if (method == DIRECT) {
    Femlib::synchronize();
    cout << "]" << endl;
    cout.flush();
  }
}


///////////////////////////////////////////////////////////////////////////////
// From pressure.C, all redefinitions:



void PBCmgr::maintain (const integer    step   ,
		       const Field*     P      ,
		       const AuxField** Us     ,
		       const AuxField** Uf     ,
		       const integer    timedep)
// ---------------------------------------------------------------------------
// Update storage for evaluation of high-order pressure boundary condition.
// Storage order for each edge represents a CCW traverse of element boundaries.
//
// If the velocity field varies in time on HOPB field boundaries (e.g. due
// to time-varying BCs) the local fluid acceleration will be estimated
// from input velocity fields by explicit extrapolation if timedep is true.
// This correction cannot be carried out at the first timestep, since the
// required extrapolation cannot be done.  If the acceleration is known,
// (for example, a known reference frame acceleration) it is probably better
// to leave timedep unset, and to use PBCmgr::accelerate() to add in the
// accelerative term.  Note also that since grad P is dotted with n, the
// unit outward normal, at a later stage, timedep only needs to be set if
// there are wall-normal accelerative terms.
//
// Field* master gives a list of pressure boundary conditions with which to
// traverse storage areas (note this assumes equal-order interpolations).
//
// No smoothing is done to high-order spatial derivatives computed here.
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const real      nu    =           Femlib::value ("KINVIS");
  const real      invDt = 1.0     / Femlib::value ("D_T");
  const integer   nTime = (integer) Femlib::value ("N_TIME");
  const integer   kfund = (integer) Femlib::value ("K_FUND");
  const integer   nEdge = P -> _nbound;
  const integer   nZ    = P -> _nz;
  const integer   nP    =  Geometry::nP();

  const AuxField* Ux = Us[0];
  const AuxField* Uy = Us[1];
  const AuxField* Uz = (nZ > 1) ? Us[2] : 0;
  const AuxField* Nx = Uf[0];
  const AuxField* Ny = Uf[1];

  const vector<Boundary*>& BC = P -> _bsys -> BCs (0);
  register Boundary*       B;
  register integer         i, k, q;
  integer                  offset, skip, Je;

  vector<real>       work (4 * nP + Integration::OrderMax + 1);

  // -- Roll grad P storage area up, load new level of nonlinear terms Uf.

  rollv (Pnx, nTime);
  rollv (Pny, nTime);

  for (i = 0; i < nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();
    
    for (k = 0; k < nZ; k++) {
      Veclib::copy (nP, Nx -> _plane[k] + offset, skip, Pnx[0][i][k], 1);
      Veclib::copy (nP, Ny -> _plane[k] + offset, skip, Pny[0][i][k], 1);
    }
  }

  // -- Add in -nu * curl curl u.

  real  *UxRe, *UxIm, *UyRe, *UyIm, *UzRe, *UzIm, *tmp;
  real* xr    = work();
  real* xi    = xr + nP;
  real* yr    = xi + nP;
  real* yi    = yr + nP;
  real* alpha = yi + nP;

  for (i = 0; i < nEdge; i++) {
    B      = BC[i];
    offset = B -> dOff ();
    skip   = B -> dSkip();

    UxRe = Ux -> _plane[0];
    UyRe = Uy -> _plane[0];
    B -> curlCurl (0, UxRe, 0, UyRe, 0, 0, 0, xr, 0, yr, 0);
    Blas::axpy (nP, -nu, xr, 1, Pnx[0][i][0], 1);
    Blas::axpy (nP, -nu, yr, 1, Pny[0][i][0], 1);

    UxRe = Ux -> _plane[1];
    UxIm = Ux -> _plane[2];
    UyRe = Uy -> _plane[1];
    UyIm = Uy -> _plane[2];
    UzRe = Uz -> _plane[1];
    UzIm = Uz -> _plane[2];

    B -> curlCurl (kfund, UxRe,UxIm, UyRe,UyIm, UzRe,UzIm, xr,xi, yr,yi);

    Blas::axpy (nP, -nu, xr, 1, Pnx[0][i][1], 1);
    Blas::axpy (nP, -nu, xi, 1, Pnx[0][i][2], 1);
    Blas::axpy (nP, -nu, yr, 1, Pny[0][i][1], 1);
    Blas::axpy (nP, -nu, yi, 1, Pny[0][i][2], 1);
  }

  if (timedep) {

    // -- Estimate -du / dt by backwards differentiation and add in.
    
    if (step > 1) {
      Je  = min (step - 1, nTime);
      tmp = xr;
      Integration::StifflyStable (Je, alpha);
      
      for (i = 0; i < nEdge; i++) {
	B      = BC[i];
	offset = B -> dOff ();
	skip   = B -> dSkip();

	for (k = 0; k < nZ; k++) {
	  Veclib::copy (nP, Ux -> _plane[k] + offset, skip, tmp, 1);
	  Blas::scal   (nP, alpha[0], tmp, 1);
	  for (q = 0; q < Je; q++)
	    Blas::axpy (nP, alpha[q + 1], Unx[q][i][k], 1, tmp, 1);
	  Blas::axpy (nP, -invDt, tmp, 1, Pnx[0][i][k], 1);
	  
	  Veclib::copy (nP, Uy -> _plane[k] + offset, skip, tmp, 1);
	  Blas::scal   (nP, alpha[0], tmp, 1);
	  for (q = 0; q < Je; q++)
	    Blas::axpy (nP, alpha[q + 1], Uny[q][i][k], 1, tmp, 1);
	  Blas::axpy (nP, -invDt, tmp, 1, Pny[0][i][k], 1);
	}
      }
    }

    // -- Roll velocity storage area up, load new level.

    rollv (Unx, nTime);
    rollv (Uny, nTime);
      
    for (i = 0; i < nEdge; i++) {
      B      = BC[i];
      offset = B -> dOff ();
      skip   = B -> dSkip();
    
      for (k = 0; k < nZ; k++) {
	Veclib::copy (nP, Ux -> _plane[k] + offset, skip, Unx[0][i][k], 1);
	Veclib::copy (nP, Uy -> _plane[k] + offset, skip, Uny[0][i][k], 1);
      }
    }
  }
}


void PBCmgr::evaluate (const integer id   ,
		       const integer np   ,
		       const integer plane,
		       const integer step ,
		       const real*   nx   ,
		       const real*   ny   ,
		       real*         tgt  )
// ---------------------------------------------------------------------------
// Load PBC value with values obtained from HOBC multi-level storage.
//
// The boundary condition for evaluation is
//
//   dP       /                           du  \
//   -- = n . | N(u) - a + f + \nu*L(u) - --  |  =  n . grad P.
//   dn   ~   \ ~ ~    ~   ~       ~ ~    dt  /     ~
//
// Grad P is estimated at the end of the current timestep using explicit
// extrapolation, then dotted into n.
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  if (step < 1) return;

  register integer q, Je = (integer) Femlib::value ("N_TIME");
  vector<real>     work (Integration::OrderMax + 2 * np);
  real*            beta  = work();
  real*            tmpX  = beta + Integration::OrderMax;
  real*            tmpY  = tmpX + np;

  Je = min (step, Je);
  Integration::Extrapolation (Je, beta);
  Veclib::zero (2 * np, tmpX, 1);
  
  for (q = 0; q < Je; q++) {
    Blas::axpy (np, beta[q], Pnx[q][id][plane], 1, tmpX, 1);
    Blas::axpy (np, beta[q], Pny[q][id][plane], 1, tmpY, 1);
  }
    
  Veclib::vmul  (np, nx, 1, tmpX, 1, tgt, 1);
  Veclib::vvtvp (np, ny, 1, tmpY, 1, tgt, 1, tgt, 1);
}


///////////////////////////////////////////////////////////////////////////////
// From statistics.C: all redefinitions: no computation of Reynolds stresses:


Statistics::Statistics (Domain*            D    ,
			vector<AuxField*>& extra) : 
// ---------------------------------------------------------------------------
// Store averages for all Domain Fields, and any extra AuxFields
// supplied.
//
// Try to initialize from file session.avg, failing that set all
// buffers to zero.  Number of fields in file should be same as
// Statistics::avg buffer.
//
// Redefinition.
// ---------------------------------------------------------------------------
  name (D -> name),
  base (D)
{
  integer       i, j;
  const integer ND    = Geometry::nDim();
  const integer NF    = base -> u.getSize();
  const integer NE    = extra.getSize();
  const integer NR    = 0;
  const integer NT    = NF + NE + NR;
  const integer nz    = Geometry::nZProc();
  const integer ntot  = Geometry::nTotProc();
  real*         alloc = new real [(size_t) NT * ntot];

  ROOTONLY cout << "-- Initializing averaging  : ";  

  // -- Set pointers, allocate storage.

  src.setSize (NF + NE);	// -- Straight running average of these.
  avg.setSize (NT);		// -- Additional are computed from src.
  
  for (i = 0; i < NF; i++) src[     i] = (AuxField*) base -> u[i];
  for (i = 0; i < NE; i++) src[NF + i] = extra[i];

  for (j = 0, i = 0; i < NF + NE; i++, j++)
    avg[i] = new AuxField (alloc+j*ntot, nz, base -> elmt, src[i] -> name());
  for (i = 0; i < NT - NF - NE; i++, j++)
    avg[i + NF + NE] = new AuxField (alloc+j*ntot, nz, base -> elmt, 'A' + i);

  // -- Initialise averages, either from file or zero.
  //    This is much the same as Domain input routine.

  char     s[StrMax];
  ifstream file (strcat (strcpy (s, name), ".avg"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << s;
      cout.flush();
    }
    file >> *this;
    file.close();
  } else {			// -- No file, set to zero.
    ROOTONLY cout << "set to zero";
    for (i = 0; i < NT; i++) *avg[i] = 0.0;
    navg = 0;
  }

  ROOTONLY cout << endl;
}


void Statistics::update (AuxField** work)
// ---------------------------------------------------------------------------
// Update running averages, using zeroth time level of work as
// workspace. 
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer NT = avg.getSize();
  const integer ND = Geometry::nDim();
  const integer NR = 0;
  const integer NA = NT - NR;

  // -- Running averages only.

  for (i = 0; i < NA; i++) {
    *avg[i] *= (real)  navg;
    *avg[i] += *src[i];
    *avg[i] /= (real) (navg + 1);
  }

  navg++;
}


void Statistics::dump ()
// ---------------------------------------------------------------------------
// Similar to Domain::dump.
//
// Redefinition.
// ---------------------------------------------------------------------------
{
  const integer step     = base -> step;
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;

  ofstream      output;
  const integer NT = avg.getSize();
  const integer ND = Geometry::nDim();

  ROOTONLY {
    const char    routine[] = "Statistics::dump";
    const integer verbose   = (integer) Femlib::value ("VERBOSE");
    const integer chkpoint  = (integer) Femlib::value ("CHKPOINT");
    char          dumpfl[StrMax], backup[StrMax], command[StrMax];

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, name), ".avg");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, name), ".ave");
	if (!initial) {
	  strcat  (strcpy (backup, name), ".ave.bak");
	  sprintf (command, "mv ./%s ./%s", dumpfl, backup);
	  system  (command);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, name), ".avg");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }

    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }

  output << *this;

  ROOTONLY output.close();
}


