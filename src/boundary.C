///////////////////////////////////////////////////////////////////////////////
// boundary.C
//
// Boundaries correspond to domain edges that have boundary conditions
// applied (as opposed to periodic edges).  The ordering of internal
// storage for condition values and geometric factors corresponds to
// CCW traverse of 2D element edges.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>


Boundary::Boundary (const integer Ident ,
		    const integer Voffst,
		    const char*   Bgroup,
		    Condition*    Bcondn,
		    Element*      Elmt  ,
		    const integer Side  ) :

		    id           (Ident ),
		    bgroup       (Bgroup),
		    bcondn       (Bcondn),
		    elmt         (Elmt  ),
		    side         (Side  ),
		    voffset      (Voffst)
// ---------------------------------------------------------------------------
// Constructor.  Allocate new memory for value & geometric factors.
// ---------------------------------------------------------------------------
{
  const integer np = elmt -> nKnot();

  nx   = new real [np];
  ny   = new real [np];
  area = new real [np];

  elmt -> sideOffset (side, doffset, dskip);
  elmt -> sideGeom   (side, nx, ny, area);
}


void Boundary::evaluate (const integer plane,
			 const integer step ,
			 real*         tgt  ) const
// ---------------------------------------------------------------------------
// Load boundary condition storage area with numeric values.
// ---------------------------------------------------------------------------
{
  const integer np = elmt -> nKnot();

  bcondn -> evaluate (np, id, plane, elmt, side, step, nx, ny, tgt);
}


void Boundary::set (const real*    src,
		    const integer* b2g,
		    real*          tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to over-ride (set) values
// in globally-numbered tgt.  This will only take place on essential BCs.
// ---------------------------------------------------------------------------
{
  bcondn -> set (elmt, side, b2g, src, tgt);
}


void Boundary::sum (const real*    src,
		    const integer* b2g,
		    real*          tgt) const
// ---------------------------------------------------------------------------
// Use (boundary condition) values in src to add in the boundary-integral
// terms generated in constructing the weak form of the MWR into globally-
// numbered tgt.  This will only take place on natural BCs.
// ---------------------------------------------------------------------------
{
  bcondn -> sum (elmt, side, b2g, src, area, tgt);
}


void Boundary::print () const
// ---------------------------------------------------------------------------
// (Debugging) utility to print internal information.
// ---------------------------------------------------------------------------
{
  char info[StrMax];

  cout << "** Boundary id: " << id + 1 << " -> ";
  cout <<     elmt ->  ID() + 1 << "." << side + 1;
  cout << " (Element id.side)" << endl;
  
  bcondn -> describe (info);

  cout << info << endl;

  cout << "  " << elmt -> nKnot() << " (number of points along edge)" << endl;
  cout << "         nx             ny             area";
  cout << endl;
  
  printVector (cout, "rrr", elmt -> nKnot(), nx, ny, area);
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

  const integer np   = elmt -> nKnot();
  const integer ntot = elmt -> nTot();
  const integer doff = elmt -> dOff();

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
    const integer    np = nKnot();

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
				  const real  mu ,
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
    register integer i;
    const integer    np = nKnot(), offset = elmt -> dOff();

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

    Force.x *= -mu;
    Force.y *= -mu;
  }

  return Force;
}


void Boundary::addForGroup (const char* grp,
			    const real  val,
			    real*       tgt) const
// ---------------------------------------------------------------------------
// Add val to tgt if this Boundary falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, bcondn -> group()) == 0)
    Veclib::sadd (nKnot(), val, tgt, 1, tgt, 1);
}


void Boundary::setForGroup (const char* grp,
			    const real  val,
			    real*       tgt) const
// ---------------------------------------------------------------------------
// Set tgt to val if this Boundary falls in group.
// ---------------------------------------------------------------------------
{
  if (strcmp (grp, bcondn -> group()) == 0)
    Veclib::fill (nKnot(), val, tgt, 1);
}
