//////////////////////////////////////////////////////////////////////////////
// egde.C: implement element-edge operators.
//
// Copyright (c) 2003<-->$Date$, Hugh Blackburn
//
// Edges, like boundaries (to which they contribute) always belong to
// a group.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


Edge::Edge (const char*    grp ,
	    const Element* elmt,
	    const integer  side) : 
// ---------------------------------------------------------------------------
// Class constructor.
// ---------------------------------------------------------------------------
  _np   (Geometry::nP()),
  _elmt (elmt),
  _side (side)
{
  const char    routine[] = "Edge::Edge";
  const integer npnp      = sqr (_np);
  char          err[StrMax];

  strcpy ((_group = new char [static_cast<size_t>(strlen (grp) + 1)]), grp);

  _x    = new real [static_cast<size_t>(5 * _np)];
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
    sprintf (err, "cannot construct edge %1d", _side + 1);
    message (routine, err, ERROR);
  }

  _elmt -> sideGeom (_side, _x, _y, _nx, _ny, _area);

  if (Geometry::cylindrical()) Veclib::vmul (_np, _area, 1, _y, 1, _area, 1);
}


void Edge::geometry (real* X   ,
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


void Edge::curlCurl (const integer k  ,
		     const real*   Ur ,
		     const real*   Ui ,
		     const real*   Vr ,
		     const real*   Vi ,
		     const real*   Wr ,
		     const real*   Wi ,
		     real*         xr ,
		     real*         xi ,
		     real*         yr ,
		     real*         yi ,
		     real*         wrk) const
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
//
// Work vector wrk 5*np*np + 3*np long.
// ---------------------------------------------------------------------------
{
  const integer npnp     = sqr (_np);
  const integer elmtOff  = _elmt -> ID() * npnp;
  const integer localOff = _doffset - elmtOff;

  real* gw = wrk;
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
      _elmt -> sideDivY (_side, w, t);
      Veclib::vadd      (_np, xr, 1, t, 1, xr, 1);
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
      _elmt -> sideDivY (_side, w, t);
      Veclib::vadd      (_np, xr, 1, t, 1, xr, 1);
    }

    Veclib::copy      (npnp, Ui, 1, uy, 1);
    Veclib::copy      (npnp, Vi, 1, vx, 1);
    _elmt -> grad     (vx, uy, gw);
    Veclib::vsub      (npnp, vx, 1, uy, 1, w, 1);
    _elmt -> sideGrad (_side, w, yi, xi, ew);
    Veclib::neg       (_np, yi, 1);

    if (Geometry::cylindrical()) {
      _elmt -> sideDivY (_side, w, t);
      Veclib::vadd     (_np, xi, 1, t, 1, xi, 1);
    }

    // -- Semi-Fourier terms based on Wr.

    Veclib::copy  (npnp, Wr, 1, vx, 1);
    Veclib::copy  (npnp, Wr, 1, uy, 1);
    _elmt -> grad (vx, uy, gw);
    if (Geometry::cylindrical()) {
      _elmt -> sideDivY  (_side, vx,  t);
      Blas::axpy         (_np,  betaK, t, 1, xi, 1);
      _elmt -> sideDivY  (_side, uy,  t);
      Blas::axpy         (_np,  betaK, t, 1, yi, 1);
      _elmt -> sideDivY2 (_side, Wr,  t);
      Blas::axpy         (_np,  betaK, t, 1, yi, 1);
    } else {
      Blas::axpy (_np, betaK, vx + localOff, _dskip, xi, 1);
      Blas::axpy (_np, betaK, uy + localOff, _dskip, yi, 1);
    }

    // -- Semi-Fourier terms based on Wi.

    Veclib::copy  (npnp, Wi, 1, vx, 1);
    Veclib::copy  (npnp, Wi, 1, uy, 1);
    _elmt -> grad (vx, uy, gw);
    if (Geometry::cylindrical()) {
      _elmt -> sideDivY  (_side, vx,   t);
      Blas::axpy         (_np, -betaK, t, 1, xr, 1);
      _elmt -> sideDivY  (_side, uy,   t);
      Blas::axpy         (_np, -betaK, t, 1, yr, 1);
      _elmt -> sideDivY2 (_side, Wi,   t);
      Blas::axpy         (_np, -betaK, t, 1, yr, 1);
    } else {
      Blas::axpy (_np, -betaK, vx + localOff, _dskip, xr, 1);
      Blas::axpy (_np, -betaK, uy + localOff, _dskip, yr, 1);
    }

    // -- Fourier second derivatives in the third direction.

    if (Geometry::cylindrical()) {
      _elmt -> sideDivY2 (_side, Ur,   t);
      Blas::axpy         (_np, betaK2, t, 1, xr, 1);
      _elmt -> sideDivY2 (_side, Ui,   t);
      Blas::axpy         (_np, betaK2, t, 1, xi, 1);
      _elmt -> sideDivY2 (_side, Vr,   t);
      Blas::axpy         (_np, betaK2, t, 1, yr, 1);
      _elmt -> sideDivY2 (_side, Vi,   t);
      Blas::axpy         (_np, betaK2, t, 1, yi, 1);
    } else {
      Blas::axpy (_np, betaK2, Ur + localOff, _dskip, xr, 1);
      Blas::axpy (_np, betaK2, Ui + localOff, _dskip, xi, 1);
      Blas::axpy (_np, betaK2, Vr + localOff, _dskip, yr, 1);
      Blas::axpy (_np, betaK2, Vi + localOff, _dskip, yi, 1);
    }
  }
}


Vector Edge::normalTraction (const char* grp,
			     const real* p  ,
			     real*       wrk) const
// ---------------------------------------------------------------------------
// Compute normal tractive force on this boundary segment, if it lies
// in group called grp, using p as a pressure stress field data area.
//
// Wrk is a work vector elmt_np_max long.
// ---------------------------------------------------------------------------
{
  integer i;
  Vector  Force = {0.0, 0.0, 0.0};

  if (strcmp (grp, _group) == 0) {

    _elmt -> sideGet (_side, p, wrk);

    for (i = 0; i < _np; i++) {
      Force.x += _nx[i] * wrk[i] * _area[i];
      Force.y += _ny[i] * wrk[i] * _area[i];
    }
  }

  return Force;
}


Vector Edge::tangentTraction (const char* grp,
			      const real* u  ,
			      const real* v  ,
			      real*       wrk) const
// ---------------------------------------------------------------------------
// Compute viscous stress on this boundary segment, if it lies in group grp.
// u is data area for first velocity component field, v is for second.
//
// Work is a work vector, 4 * _np long.
// ---------------------------------------------------------------------------
{
  Vector Force = {0.0, 0.0, 0.0};

  if (strcmp (grp, _group) == 0) {
    register integer i;
    real         *ux = wrk + 2 * _np, *uy = wrk + 3 * _np;

    _elmt -> sideGrad (_side, u + _doffset, ux, uy, wrk);

    for (i = 0; i < _np; i++) {
      Force.x += (2.0*ux[i]*_nx[i] + uy[i]*_ny[i]) * _area[i];
      Force.y +=                     uy[i]*_nx[i]  * _area[i];
    }

    _elmt -> sideGrad (_side, v + _doffset, ux, uy, wrk);

    for (i = 0; i < _np; i++) {
      Force.x +=                     ux[i]*_ny[i]  * _area[i];
      Force.y += (2.0*uy[i]*_ny[i] + ux[i]*_nx[i]) * _area[i];
    }
  }

  return Force;
}


real Edge::normalFlux (const char* grp,
		       const real* u  ,
		       const real* v  ,
		       real*       wrk) const
// ---------------------------------------------------------------------------
// Compute edge-normal flux, with u being x-component velocity and v
// being y-component, if this edge lies in group grp. Work vector wrk
// is 2*_np long.
// ---------------------------------------------------------------------------
{
  integer i;
  real    flux = 0.0;
  real    *U = wrk, *V = wrk + _np;
  
  if (strcmp (grp, _group) == 0) {

    _elmt -> sideGet (_side, u, U);
    _elmt -> sideGet (_side, v, V);

    for (i = 0; i < _np; i++)
      flux += (U[i]*_nx[i] + V[i]*_ny[i]) * _area[i];
  }

  return flux;
}


real Edge::gradientFlux (const char* grp,
			 const real* src,
			 real*       wrk) const
// ---------------------------------------------------------------------------
// Compute wall-normal gradient flux of field src on this boundary
// segment, if it lies in group grp.  Wrk is a work vector, 4 *
// elmt_np_max long.  NB: n is a unit outward normal, with no
// component in Fourier direction.  NB: For cylindrical coords, it is
// assumed we are dealing with a scalar!
// ---------------------------------------------------------------------------
{
  register real dcdn = 0.0;
  
  if (strcmp (grp, _group) == 0) {
    register integer  i;
    register real *cx = wrk, *cy = wrk + _np, *r = wrk + _np + _np;

    _elmt -> sideGrad (_side, src + _doffset, cx, cy, r);
    for (i = 0; i < _np; i++)
      dcdn += (cx[i]*_nx[i] + cy[i]*_ny[i]) * _area[i];
  }

  return dcdn;
}


void Edge::addForGroup (const char* grp,
			const real  val,
			real*       tgt) const
// ---------------------------------------------------------------------------
// Add val to tgt if this Edge falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, _group) == 0) Veclib::sadd (_np, val, tgt, 1, tgt, 1);
}


void Edge::setForGroup (const char* grp,
			const real  val,
			real*       tgt) const
// ---------------------------------------------------------------------------
// Set tgt to val if this Edge falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, _group) == 0) Veclib::fill (_np, val, tgt, 1);
}


void Edge::get (const real* src,
		real*       tgt) const
// ---------------------------------------------------------------------------
// Load np-long tgt (representing storage along edge of element) from
// element-wise data storage src.
// ---------------------------------------------------------------------------
{
  Veclib::copy (_np, src + _doffset, _dskip, tgt, 1);
}


void Edge::mulY (real* tgt) const
// ---------------------------------------------------------------------------
// Multiply tgt by y (i.e. radius) along this edge.
// ---------------------------------------------------------------------------
{
  Veclib::vmul (_np, tgt, 1, _y, 1, tgt, 1);
}


void Edge::divY (real* tgt) const
// ---------------------------------------------------------------------------
// Divide tgt by y (typically, radius) along this edge.
// ---------------------------------------------------------------------------
{
  register integer i;
  real             invr;

  for (i = 0; i < _np; i++) {
    invr = (_y[i] > EPSDP) ? 1.0/_y[i] : 0.0;
    tgt[i] *= invr;
  }
}
