///////////////////////////////////////////////////////////////////////////////
// boundary.C: implement Boundary class functions.
//
// Copyright (c) 1994,2003 Hugh Blackburn
//
// SYNOPSIS
// --------
// Boundaries correspond to domain edges that have boundary conditions
// applied (as opposed to periodic edges).  The ordering of internal
// storage for condition values and geometric factors corresponds to
// CCW traverse of 2D element edges.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <Sem.h>


Boundary::Boundary (const integer    Ident ,
		    const char*      Bgroup,
		    const Condition* Bcondn,
		    const Element*   Elmt  ,
		    const integer    Side  ) :
// ---------------------------------------------------------------------------
// Constructor.  Allocate new memory for value & geometric factors.
// ---------------------------------------------------------------------------
  _id      (Ident ),
  _np      (Geometry::nP()),
  _bgroup  (Bgroup),
  _bcondn  (Bcondn),
  _elmt    (Elmt  ),
  _side    (Side  )
{
  const char    routine[] = "Boundary::Boundary";
  const integer npnp      = sqr (_np);
  char          err[StrMax];

  _x    = new real [(size_t) 5 * _np];
  _y    = _x  + _np;
  _nx   = _y  + _np;
  _ny   = _nx + _np;
  _area = _ny + _np;

  _doffset = _elmt -> ID() * npnp;
  switch (_side) {
  case 0: _doffset += 0;               _dskip = 1;    break;
  case 1: _doffset += (_np - 1);       _dskip = _np;  break;
  case 2: _doffset += _np * (_np - 1); _dskip = -1;   break;
  case 3: _doffset += 0;               _dskip = -_np; break;
  default:
    sprintf (err, "cannot construct side %1d", _side + 1);
    message (routine, err, ERROR);
  }
  _elmt -> sideGeom (_side, _x, _y, _nx, _ny, _area);
}


void Boundary::geometry (real* X   ,
			 real* Y   ,
			 real* Nx  ,
			 real* Ny  ,
			 real* Area) const
// ---------------------------------------------------------------------------
// Copy internal geometric info for exterior use.
// ---------------------------------------------------------------------------
{
  Veclib::copy (_np, _x, 1, X, 1);
  Veclib::copy (_np, _y, 1, Y, 1);
  if (Nx)   Veclib::copy (_np, _nx,   1, Nx,   1);
  if (Ny)   Veclib::copy (_np, _ny,   1, Ny,   1);
  if (Area) Veclib::copy (_np, _area, 1, Area, 1);
}


void Boundary::evaluate (const integer plane,
			 const integer step ,
			 real*         tgt  ) const
// ---------------------------------------------------------------------------
// Load boundary condition storage area with numeric values.
// ---------------------------------------------------------------------------
{
  _bcondn -> evaluate (_np, _id, plane, _elmt, _side, step, _nx, _ny, tgt);
}


void Boundary::set (const real*    src,
		    const integer* b2g,
		    real*          tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to over-ride (set) values
// in globally-numbered tgt.  This will only take place on essential BCs.
//
// b2g is a pointer to the global node numbers for the appropriate
// element's edge nodes.
// ---------------------------------------------------------------------------
{
  _bcondn -> set (_side, b2g, src, tgt);
}


void Boundary::sum (const real*    src,
		    const integer* b2g,
		    real*          wrk,
		    real*          tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to add in the boundary-integral
// terms generated in constructing the weak form of the MWR into globally-
// numbered tgt.  This will only take place on natural BCs.
//
// b2g is a pointer to the global node numbers for the appropriate
// element's edge nodes.  wrk is a work array, np long.
// ---------------------------------------------------------------------------
{
  _bcondn -> sum (_side, b2g, src, _area, wrk, tgt);
}


void Boundary::augmentSC (const integer  nband ,
			  const integer  nsolve,
			  const integer* b2g   ,
			  real*          work  ,
			  real*          H     ) const
// ---------------------------------------------------------------------------
// Add in diagonal terms <K, w> to (banded LAPACK) H on mixed BCs.
// Work array must be np long.
// ---------------------------------------------------------------------------
{
  _bcondn -> augmentSC (_side, nband, nsolve, b2g + bOff(), _area, work, H);
}


void Boundary::augmentOp (const integer* b2g ,
			  const real*    src ,
			  real*          tgt ) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise Helmholtz
// operations where there are mixed BCs.  Add in diagonal terms
// <K*src, w> to tgt.  Both src and tgt are globally-numbered vectors.
// ---------------------------------------------------------------------------
{
  _bcondn -> augmentOp (_side, b2g + bOff(), _area, src, tgt);
}


void Boundary::augmentDg (const integer* b2g,
			  real*          tgt) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  _bcondn -> augmentDg (_side, b2g + bOff(), _area, tgt);
}


void Boundary::print () const
// ---------------------------------------------------------------------------
// (Debugging) utility to print internal information.
// ---------------------------------------------------------------------------
{
  char info[StrMax];

  cout << "** Boundary id: " << _id + 1 << " -> ";
  cout << _elmt ->  ID() + 1 << "." << _side + 1;
  cout << " (Element id.side)" << endl;
  
  _bcondn -> describe (info);

  cout << info << endl;

  cout << "  " << _np << " (number of points along edge)" << endl;
  cout << "         nx             ny             area";
  cout << endl;
  
  printVector (cout, "rrr", _np, _nx, _ny, _area);
}


void Boundary::curlCurl (const integer k ,
			 const real*   Ur,
			 const real*   Ui,
			 const real*   Vr,
			 const real*   Vi,
			 const real*   Wr,
			 const real*   Wi,
			 real*         xr,
			 real*         xi,
			 real*         yr,
			 real*         yi) const
// ---------------------------------------------------------------------------
// Generate (the Fourier mode equivalent of) curl curl u along this boundary.
//
// Input k is the Fourier-mode index.
//
// Input pointers Ur, Ui etc correspond to the real and imaginary planes of
// data for the three components of vector field u corresponding to the
// kth Fourier mode.  The third component is treated as the transformed
// direction.
//
// Output pointers are to the (real and imaginary parts of) the first and
// second components of curl curl u along this boundary edge.  The third
// component is not computed as it is not required by the application.
//
// When k == 0, all the imaginary components, also the third velocity vector
// component pointers are not used, and may be provided as NULL values.
// This allows the same routine to be used for 2D solutions.
// ---------------------------------------------------------------------------
{
  const int npnp     = sqr (_np);
  const int elmtOff  = _elmt -> ID() * npnp;
  const int localOff = _doffset - elmtOff;

  static vector<real> work (5 * npnp + 3 * _np);
  real* gw = &work[0];
  real* ew = gw + npnp + npnp;
  real* w  = ew + _np  + _np;
  real* vx = w  + npnp;
  real* uy = vx + npnp;
  real* t  = uy + npnp;
  
  // -- Make pointers to current element storage.

  Ur += elmtOff; Ui += elmtOff;
  Vr += elmtOff; Vi += elmtOff;
  Wr += elmtOff; Wi += elmtOff;

  if (k == 0) {			// -- Zeroth mode / 2D.

    Veclib::copy (npnp, Ur, 1, uy, 1);
    Veclib::copy (npnp, Vr, 1, vx, 1);

    _elmt -> grad (vx, uy, gw);

    // -- (Z-component of) vorticity, w = dv/dx - du/dy.

    Veclib::vsub (npnp, vx, 1, uy, 1, w, 1);

    // -- Find dw/dx & dw/dy on appropriate edge.

    _elmt -> sideGrad (_side, w, yr, xr, ew);

    // -- Add in cylindrical space modification to complete x-component.

    if (Geometry::cylindrical()) {
      _elmt -> sideGet (_side, w, t);

      Veclib::vmul (_np, xr, 1, _y, 1, xr, 1);
      Veclib::vmul (_np, yr, 1, _y, 1, yr, 1);
      Veclib::vadd (_np, xr, 1, t,  1, xr, 1);
    }

    // -- Sign change to complete y-component of curl curl u.
    
    Veclib::neg (_np, yr, 1);

  } else {			// -- 3D.

    const real betaK  = k * Femlib::value ("BETA");
    const real betaK2 = sqr (betaK);

    // -- Make the equivalents of the 2D terms above.

    Veclib::copy      (npnp, Ur, 1, uy, 1);
    Veclib::copy      (npnp, Vr, 1, vx, 1);
    _elmt -> grad     (vx, uy, gw);
    Veclib::vsub      (npnp, vx, 1, uy, 1, w, 1);
    _elmt -> sideGrad (_side, w, yr, xr, ew);
    Veclib::neg       (_np, yr, 1);

    if (Geometry::cylindrical()) {
      _elmt -> sideGet (_side, w, t);
      Veclib::vmul (_np, xr, 1, _y, 1, xr, 1);
      Veclib::vmul (_np, yr, 1, _y, 1, yr, 1);
      Veclib::vadd (_np, xr, 1, t,  1, xr, 1);     
    }

    Veclib::copy      (npnp, Ui, 1, uy, 1);
    Veclib::copy      (npnp, Vi, 1, vx, 1);
    _elmt -> grad     (vx, uy, gw);
    Veclib::vsub      (npnp, vx, 1, uy, 1, w, 1);
    _elmt -> sideGrad (_side, w, yi, xi, ew);
    Veclib::neg       (_np, yi, 1);

    if (Geometry::cylindrical()) {
      _elmt -> sideGet (_side, w, t);
      Veclib::vmul (_np, xi, 1, _y, 1, xi, 1);
      Veclib::vmul (_np, yi, 1, _y, 1, yi, 1);
      Veclib::vadd (_np, xi, 1, t,  1, xi, 1);   
    }

    // -- Semi-Fourier terms based on Wr.

    Veclib::copy  (npnp, Wr, 1, vx, 1);
    Veclib::copy  (npnp, Wr, 1, uy, 1);
    _elmt -> grad (vx, uy, gw);
    if (Geometry::cylindrical()) {
      _elmt -> sideDivY (_side, Wr, t);
      Blas::axpy (_np, betaK, vx + localOff, _dskip, xi, 1);
      Blas::axpy (_np, betaK, uy + localOff, _dskip, yi, 1);
      Blas::axpy (_np, betaK, t, 1, yi, 1);
    } else {
      Blas::axpy (_np, betaK, vx + localOff, _dskip, xi, 1);
      Blas::axpy (_np, betaK, uy + localOff, _dskip, yi, 1);
    }

    // -- Semi-Fourier terms based on Wi.

    Veclib::copy  (npnp, Wi, 1, vx, 1);
    Veclib::copy  (npnp, Wi, 1, uy, 1);
    _elmt -> grad (vx, uy, gw);
    if (Geometry::cylindrical()) {
      _elmt -> sideDivY (_side, Wi, t);
      Blas::axpy (_np, -betaK, vx + localOff, _dskip, xr, 1);
      Blas::axpy (_np, -betaK, uy + localOff, _dskip, yr, 1);
      Blas::axpy (_np, -betaK, t, 1, yr, 1);
    } else {
      Blas::axpy (_np, -betaK, vx + localOff, _dskip, xr, 1);
      Blas::axpy (_np, -betaK, uy + localOff, _dskip, yr, 1);
    }

    // -- Fourier second derivatives in the third direction.

    if (Geometry::cylindrical()) {
      _elmt -> sideDivY (_side, Ur,   t);
      Blas::axpy        (_np, betaK2, t, 1, xr, 1);
      _elmt -> sideDivY (_side, Ui,   t);
      Blas::axpy        (_np, betaK2, t, 1, xi, 1);
      _elmt -> sideDivY (_side, Vr,   t);
      Blas::axpy        (_np, betaK2, t, 1, yr, 1);
      _elmt -> sideDivY (_side, Vi,   t);
      Blas::axpy        (_np, betaK2, t, 1, yi, 1);
    } else {
      Blas::axpy (_np, betaK2, Ur + localOff, _dskip, xr, 1);
      Blas::axpy (_np, betaK2, Ui + localOff, _dskip, xi, 1);
      Blas::axpy (_np, betaK2, Vr + localOff, _dskip, yr, 1);
      Blas::axpy (_np, betaK2, Vi + localOff, _dskip, yi, 1);
    }
  }
}


Vector Boundary::normalTraction (const char* grp,
				 const real* p  ,
				 real*       wrk) const
// ---------------------------------------------------------------------------
// Compute normal tractive force on this boundary segment, if it lies
// in group called grp, using p as a pressure stress field data area.
//
// Wrk is a work vector elmt_np_max long.
// ---------------------------------------------------------------------------
{
  Vector Force = {0.0, 0.0, 0.0};

  if (strcmp (grp, _bcondn -> group()) == 0) {
    register integer i;
    const integer    np = Geometry::nP();

    Veclib::copy (_np, p + _doffset, _dskip, wrk, 1);

    for (i = 0; i < _np; i++) {
      Force.x += _nx[i] * wrk[i] * _area[i];
      Force.y += _ny[i] * wrk[i] * _area[i];
    }
  }

  return Force;
}


Vector Boundary::tangentTraction (const char* grp,
				  const real* u  ,
				  const real* v  ,
				  real*       ux ,
				  real*       uy ) const
// ---------------------------------------------------------------------------
// Compute viscous stress on this boundary segment, if it lies in group grp.
// u is data area for first velocity component field, v is for second.
// Ux and uy are work vectors, each elmt_np_max long.
// ---------------------------------------------------------------------------
{
  Vector              Force = {0.0, 0.0, 0.0};
  static vector<real> work (2 * _np);

  if (strcmp (grp, _bcondn -> group()) == 0) {
    const integer    offset = _elmt -> ID() * sqr (_np);
    register integer i;

    _elmt -> sideGrad (_side, u + offset, ux, uy, &work[0]);

    for (i = 0; i < _np; i++) {
      Force.x += (2.0*ux[i]*_nx[i] + uy[i]*_ny[i]) * _area[i];
      Force.y +=                     uy[i]*_nx[i]  * _area[i];
    }

    _elmt -> sideGrad (_side, v + offset, ux, uy, &work[0]);

    for (i = 0; i < _np; i++) {
      Force.x +=                     ux[i]*_ny[i]  * _area[i];
      Force.y += (2.0*uy[i]*_ny[i] + ux[i]*_nx[i]) * _area[i];
    }
  }

  return Force;
}


real Boundary::flux (const char* grp,
		     const real* src,
		     real*       wrk) const
// ---------------------------------------------------------------------------
// Compute wall-normal flux of field src on this boundary segment,
// if it lies in group grp.  Wrk is a work vector, 4 * elmt_np_max long.
// NB: n is a unit outward normal, with no component in Fourier direction.
// NB: For cylindrical coords, it is assumed we are dealing with a scalar!
// ---------------------------------------------------------------------------
{
  register real dcdn = 0.0;
  
  if (strcmp (grp, _bcondn -> group()) == 0) {
    const real*      data = src + _elmt -> ID() * Geometry::nTotElmt();
    register integer i;
    register real    *cx = wrk, *cy = wrk + _np, *r = wrk + _np + _np;

    _elmt -> sideGrad (_side, data, cx, cy, r);
    for (i = 0; i < _np; i++)
      dcdn += (cx[i]*_nx[i] + cy[i]*_ny[i]) * _area[i];
  }

  return dcdn;
}


void Boundary::addForGroup (const char* grp,
			    const real  val,
			    real*       tgt) const
// ---------------------------------------------------------------------------
// Add val to tgt if this Boundary falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, _bcondn -> group()) == 0)
    Veclib::sadd (_np, val, tgt, 1, tgt, 1);
}


void Boundary::setForGroup (const char* grp,
			    const real  val,
			    real*       tgt) const
// ---------------------------------------------------------------------------
// Set tgt to val if this Boundary falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, _bcondn -> group()) == 0)
    Veclib::fill (_np, val, tgt, 1);
}


void Boundary::get (const real* src,
		    real*       tgt) const
// ---------------------------------------------------------------------------
// Load np-long tgt (representing storage along edge of element) from
// element-wise data storage src.
// ---------------------------------------------------------------------------
{
  Veclib::copy (_np, src + _doffset, _dskip, tgt, 1);
}


const char* Boundary::group () const
// ---------------------------------------------------------------------------
// Return group of underlying boundary Condition.
// ---------------------------------------------------------------------------
{
  return _bcondn -> group();
}


void Boundary::mulY (real* tgt) const
// ---------------------------------------------------------------------------
// Multiply tgt by y (i.e. radius) along this edge.
// ---------------------------------------------------------------------------
{
  Veclib::vmul (_np, tgt, 1, _y, 1, tgt, 1);
}

