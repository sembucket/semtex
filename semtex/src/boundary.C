///////////////////////////////////////////////////////////////////////////////
// boundary.C: implement Boundary class functions.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// SYNOPSIS
// --------
// Boundaries correspond to domain edges that have boundary conditions
// applied (as opposed to periodic edges).  The ordering of internal
// storage for condition values and geometric factors corresponds to
// CCW traverse of 2D element edges.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


Boundary::Boundary (const integer    Ident ,
		    const char*      Bgroup,
		    const Condition* Bcondn,
		    const Element*   Elmt  ,
		    const integer    Side  ) :
// ---------------------------------------------------------------------------
// Constructor.  Allocate new memory for value & geometric factors.
// ---------------------------------------------------------------------------
  id      (Ident ),
  np      (Geometry::nP()),
  bgroup  (Bgroup),
  bcondn  (Bcondn),
  elmt    (Elmt  ),
  side    (Side  )
{
  const char    routine[] = "Boundary::Boundary";
  const integer np        = Geometry::nP();
  const integer npnp      = Geometry::nTotElmt();
  char          err[StrMax];

  x    = new real [(size_t) 5 * np];
  y    = x  + np;
  nx   = y  + np;
  ny   = nx + np;
  area = ny + np;

  doffset = elmt -> ID() * npnp;
  switch (side) {
  case 0: doffset += 0;             dskip = 1;   break;
  case 1: doffset += (np - 1);      dskip = np;  break;
  case 2: doffset += np * (np - 1); dskip = -1;  break;
  case 3: doffset += 0;             dskip = -np; break;
  default:
    sprintf (err, "cannot construct side %1d", side + 1);
    message (routine, err, ERROR);
  }
  elmt -> sideGeom (side, x, y, nx, ny, area);
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
  Veclib::copy (np, x, 1, X, 1);
  Veclib::copy (np, y, 1, Y, 1);
  if (Nx)   Veclib::copy (np, nx,   1, Nx,   1);
  if (Ny)   Veclib::copy (np, ny,   1, Ny,   1);
  if (Area) Veclib::copy (np, area, 1, Area, 1);
}


void Boundary::evaluate (const integer plane,
			 const integer step ,
			 real*         tgt  ) const
// ---------------------------------------------------------------------------
// Load boundary condition storage area with numeric values.
// ---------------------------------------------------------------------------
{
  const integer np = Geometry::nP();

  bcondn -> evaluate (np, id, plane, elmt, side, step, nx, ny, tgt);
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
  bcondn -> set (side, b2g, src, tgt);
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
  bcondn -> sum (side, b2g, src, area, wrk, tgt);
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
  bcondn -> augmentSC (side, nband, nsolve, b2g + bOff(), area, work, H);
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
  bcondn -> augmentOp (side, b2g + bOff(), area, src, tgt);
}


void Boundary::augmentDg (const integer* b2g,
			  real*          tgt) const
// ---------------------------------------------------------------------------
// This operation is used to augment the element-wise construction of
// the diagonal of the global Helmholtz matrix where there are mixed
// BCs.  Add in diagonal terms <K, w> to globally-numbered tgt.
// ---------------------------------------------------------------------------
{
  bcondn -> augmentDg (side, b2g + bOff(), area, tgt);
}


void Boundary::print () const
// ---------------------------------------------------------------------------
// (Debugging) utility to print internal information.
// ---------------------------------------------------------------------------
{
  const integer np = Geometry::nP();
  char          info[StrMax];

  cout << "** Boundary id: " << id + 1 << " -> ";
  cout <<     elmt ->  ID() + 1 << "." << side + 1;
  cout << " (Element id.side)" << endl;
  
  bcondn -> describe (info);

  cout << info << endl;

  cout << "  " << np << " (number of points along edge)" << endl;
  cout << "         nx             ny             area";
  cout << endl;
  
  printVector (cout, "rrr", np, nx, ny, area);
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
  const Geometry::CoordSys space = Geometry::system();

  const integer np   = Geometry::nP();
  const integer ntot = Geometry::nTotElmt();
  const integer doff = elmt -> ID() * ntot;

  vector<real> work (5 * ntot + np);
  real* gw = work();
  real* w  = gw + ntot + ntot;
  real* vx = w  + ntot;
  real* uy = vx + ntot;
  real* t  = uy + ntot;
  const real** DV;
  const real** DT;

  Femlib::quad (LL, np, np, 0, 0, 0, 0, 0, &DV, &DT);
      
  if (k == 0) {			// -- Zeroth mode / 2D.

    Veclib::copy (ntot, Ur + doff, 1, uy, 1);
    Veclib::copy (ntot, Vr + doff, 1, vx, 1);

    elmt -> grad (vx, uy, DV, DT, gw);

    // -- (Z-component of) vorticity, w = dv/dx - du/dy.

    Veclib::vsub (ntot, vx, 1, uy, 1, w, 1);

    // -- Find dw/dx & dw/dy on appropriate edge.

    elmt -> sideGrad (side, w, yr, xr);
    
    // -- Add in cylindrical space modification to complete x-component.

    if (space == Geometry::Cylindrical) {
      elmt -> sideDivR (side, w, t);
      Veclib::vadd     (np, xr, 1, t, 1, xr, 1);
    }

    // -- Sign change to complete y-component of curl curl u.
    
    Veclib::neg (np, yr, 1);

  } else {			// -- 3D.

    const real betaK  = k * Femlib::value ("BETA");
    const real betaK2 = sqr (betaK);
    const integer  loff   = doffset - doff; // -- Side offset in workspace.

    // -- Make the equivalents of the 2D terms above.

    Veclib::copy     (ntot, Ur + doff, 1, uy, 1);
    Veclib::copy     (ntot, Vr + doff, 1, vx, 1);
    elmt -> grad     (vx, uy, DV, DT, gw);
    Veclib::vsub     (ntot, vx, 1, uy, 1, w, 1);
    elmt -> sideGrad (side, w, yr, xr);
    Veclib::neg      (np, yr, 1);

    if (space == Geometry::Cylindrical) {
      elmt -> sideDivR (side, w, t);
      Veclib::vadd     (np, xr, 1, t, 1, xr, 1);
    }

    Veclib::copy     (ntot, Ui + doff, 1, uy, 1);
    Veclib::copy     (ntot, Vi + doff, 1, vx, 1);
    elmt -> grad     (vx, uy, DV, DT, gw);
    Veclib::vsub     (ntot, vx, 1, uy, 1, w, 1);
    elmt -> sideGrad (side, w, yi, xi);
    Veclib::neg      (np, yi, 1);

    if (space == Geometry::Cylindrical) {
      elmt -> sideDivR (side, w, t);
      Veclib::vadd     (np, xi, 1, t, 1, xi, 1);
    }

    // -- Semi-Fourier terms based on Wr.

    Veclib::copy      (ntot, Wr + doff, 1, vx, 1);
    Veclib::copy      (ntot, Wr + doff, 1, uy, 1);
    elmt -> grad      (vx, uy, DV, DT, gw);
    if (space == Geometry::Cylindrical) {
      elmt -> sideDivR  (side, vx,  t);
      Blas::axpy        (np, betaK, t,         1,     xi, 1);
      elmt -> sideDivR  (side, uy,  t);
      Blas::axpy        (np, betaK, t,         1,     yi, 1);
      elmt -> sideDivR2 (side, Wr,  t);
      Blas::axpy        (np, betaK, t,         1,     yi, 1);
    } else {
      Blas::axpy        (np, betaK, vx + loff, dskip, xi, 1);
      Blas::axpy        (np, betaK, uy + loff, dskip, yi, 1);
    }

    // -- Semi-Fourier terms based on Wi.

    Veclib::copy      (ntot, Wi + doff, 1, vx, 1);
    Veclib::copy      (ntot, Wi + doff, 1, uy, 1);
    elmt -> grad      (vx, uy, DV, DT, gw);
    if (space == Geometry::Cylindrical) {
      elmt -> sideDivR  (side, vx,   t);
      Blas::axpy        (np, -betaK, t,         1,     xr, 1);
      elmt -> sideDivR  (side, uy,   t);
      Blas::axpy        (np, -betaK, t,         1,     yr, 1);
      elmt -> sideDivR2 (side, Wi,   t);
      Blas::axpy        (np, -betaK, t,         1,     yr, 1);
    } else {
      Blas::axpy        (np, -betaK, vx + loff, dskip, xr, 1);
      Blas::axpy        (np, -betaK, uy + loff, dskip, yr, 1);
    }

    // -- Fourier second derivatives in the third direction.

    if (space == Geometry::Cylindrical) {
      elmt -> sideDivR  (side, Ur,   t);
      Blas::axpy        (np, betaK2, t,         1,     xr, 1);
      elmt -> sideDivR  (side, Ui,   t);
      Blas::axpy        (np, betaK2, t,         1,     xi, 1);
      elmt -> sideDivR2 (side, Vr,   t);
      Blas::axpy        (np, betaK2, t,         1,     yr, 1);
      elmt -> sideDivR2 (side, Vi,   t);
      Blas::axpy        (np, betaK2, t,         1,     yi, 1);
    } else {
      Blas::axpy        (np, betaK2, Ur + loff, dskip, xr, 1);
      Blas::axpy        (np, betaK2, Ui + loff, dskip, xi, 1);
      Blas::axpy        (np, betaK2, Vr + loff, dskip, yr, 1);
      Blas::axpy        (np, betaK2, Vi + loff, dskip, yi, 1);
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

  if (strcmp (grp, bcondn -> group()) == 0) {
    register integer i;
    const integer    np = Geometry::nP();

    Veclib::copy (np, p + doffset, dskip, wrk, 1);

    for (i = 0; i < np; i++) {
      Force.x += nx[i] * wrk[i] * area[i];
      Force.y += ny[i] * wrk[i] * area[i];
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
  Vector Force = {0.0, 0.0, 0.0};

  if (strcmp (grp, bcondn -> group()) == 0) {
    const integer    np     = Geometry::nP();
    const integer    offset = elmt -> ID() * Geometry::nTotElmt();
    register integer i;

    elmt -> sideGrad (side, u + offset, ux, uy);

    for (i = 0; i < np; i++) {
      Force.x += (2.0*ux[i]*nx[i] + uy[i]*ny[i]) * area[i];
      Force.y +=                    uy[i]*nx[i]  * area[i];
    }

    elmt -> sideGrad (side, v + offset, ux, uy);

    for (i = 0; i < np; i++) {
      Force.x +=                    ux[i]*ny[i]  * area[i];
      Force.y += (2.0*uy[i]*ny[i] + ux[i]*nx[i]) * area[i];
    }
  }

  return Force;
}


real Boundary::flux (const char* grp,
		     const real* src,
		     real*       wrk) const
// ---------------------------------------------------------------------------
// Compute wall-normal flux of field src on this boundary segment,
// if it lies in group grp.  Wrk is a work vector, 3 * elmt_np_max long.
// NB: n is a unit outward normal, with no component in Fourier direction.
// ---------------------------------------------------------------------------
{
  register real dcdn = 0.0;
  
  if (strcmp (grp, bcondn -> group()) == 0) {
    const integer    np = Geometry::nP();
    const real*      data = src + elmt -> ID() * Geometry::nTotElmt();
    register integer i;
    register real    *cx = wrk, *cy = wrk + np, *r = wrk + np + np;

    elmt -> sideGrad (side, data, cx, cy);

    if (Geometry::system() == Geometry::Cylindrical) {
      elmt -> sideGetR (side, r);
      for (i = 0; i < np; i++)
	dcdn += (cx[i]*nx[i] + cy[i]*ny[i]) * area[i] * r[i];
      elmt -> sideGet  (side, data, cy);
      for (i = 0; i < np; i++)
	dcdn -= cy[i]*ny[i] * area[i];

    } else {			// -- Cartesian.
      for (i = 0; i < np; i++)
	dcdn += (cx[i]*nx[i] + cy[i]*ny[i]) * area[i];
    }
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
  if (strcmp (grp, bcondn -> group()) == 0)
    Veclib::sadd (Geometry::nP(), val, tgt, 1, tgt, 1);
}


void Boundary::setForGroup (const char* grp,
			    const real  val,
			    real*       tgt) const
// ---------------------------------------------------------------------------
// Set tgt to val if this Boundary falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, bcondn -> group()) == 0)
    Veclib::fill (Geometry::nP(), val, tgt, 1);
}


void Boundary::get (const real* src,
		    real*       tgt) const
// ---------------------------------------------------------------------------
// Load np-long tgt (representing storage along edge of element) from
// element-wise data storage src.
// ---------------------------------------------------------------------------
{
  Veclib::copy (np, src + doffset, dskip, tgt, 1);
}


const char* Boundary::group () const
// ---------------------------------------------------------------------------
// Return group of underlying boundary Condition.
// ---------------------------------------------------------------------------
{
  return bcondn -> group();
}
