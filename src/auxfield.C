///////////////////////////////////////////////////////////////////////////////
// auxfield.C: routines for AuxField class, including Fourier expansions.
//
// For 2D problems, the data storage is organized by 2D Elements.
//
// For 3D problems, Fourier expansions are used in the 3rd direction, and
// each Fourier mode can be thought of as a separate 2D problem (with real
// and imaginary parts, or planes, of 2D data).  The data are then organized
// plane-by-plane, with each plane being a 2D AuxField; if in physical space
// there are nz planes of data, then there are nz/2 Fourier modes.
// Data for the Nyquist mode are stored as the imaginary part of the zeroth
// Fourier mode, but are kept zero and never evolve.
//
// The data are transformed to physical space for storage in restart files.
///////////////////////////////////////////////////////////////////////////////

static char 
RCSid[] = "$Id$";

#include <Sem.h>


AuxField::AuxField (vector<Element*>& Elts,
		    const char        name) :
		    
		    field_name       (name),
		    Elmt             (Elts)
// ---------------------------------------------------------------------------
// Allocate field storage area and integer size records.
// ---------------------------------------------------------------------------
{
  char         routine[] = "AuxField::AuxField";
  register int k;
  const int    nT = Geometry::nTot  ();
  const int    nZ = Geometry::nZ    ();
  const int    nP = Geometry::nPlane();

  if (Geometry::nElmt() != Elmt.getSize())
    message (routine, "conflicting number of elements in input data", ERROR);

  data  = new real  [nT];
  plane = new real* [nZ];

  for (k = 0; k < nZ; k++) plane[k] = data + k * nP;
}


AuxField& AuxField::operator = (const real val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (Geometry::nTot(),      data, 1);
  else              Veclib::fill (Geometry::nTot(), val, data, 1);

  return *this;
}


AuxField& AuxField::operator += (const real val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (Geometry::nTot(), val, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const real val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (Geometry::nTot(), -val, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator *= (const real val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (Geometry::nTot(),      data, 1);
  else              Blas::scal   (Geometry::nTot(), val, data, 1);

  return *this;
}


AuxField& AuxField::operator /= (const real val)
// ---------------------------------------------------------------------------
// Divide field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) message ("AuxField::op /= real", "divide by zero", ERROR);
  else              Blas::scal (Geometry::nTot(), 1.0 / val, data, 1);

  return *this;
}


AuxField& AuxField::operator = (const AuxField& f)
// ---------------------------------------------------------------------------
// This presumes the two fields conform, and copies f's value storage to LHS.
// ---------------------------------------------------------------------------
{
  register int k;
  const int    nZ = Geometry::nZ();
  const int    nP = Geometry::nPlane();
  
  for (k = 0; k < nZ; k++)
    Veclib::copy (nP, f.plane[k], 1, plane[k], 1);
  
  return *this;
}


AuxField& AuxField::operator += (const AuxField& f)
// ---------------------------------------------------------------------------
// Add f's value to this AuxField's.
// ---------------------------------------------------------------------------
{
  register int k;
  const int    nZ = Geometry::nZ();
  const int    nP = Geometry::nPlane();
  
  for (k = 0; k < nZ; k++)
    Veclib::vadd (nP, plane[k], 1, f.plane[k], 1, plane[k], 1);

  return *this;
}


AuxField& AuxField::operator -= (const AuxField& f)
// ---------------------------------------------------------------------------
// Subtract f's value from this AuxField's.
// ---------------------------------------------------------------------------
{
  register int k;
  const int    nZ = Geometry::nZ();
  const int    nP = Geometry::nPlane();
  
  for (k = 0; k < nZ; k++)  
    Veclib::vsub (nP, plane[k], 1, f.plane[k], 1, plane[k], 1);

  return *this;
}


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      i, k, offset = 0;
  const int         nE = Geometry::nElmt();
  const int         nZ = Geometry::nZ();
  const real        dz = Femlib::value ("TWOPI / BETA") / nZ;
  real              z;

  for (k = 0; k < nZ; k++) {
    z = k * dz;
    Femlib::value ("z", z);
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();

      E -> evaluate (function, plane[k] + offset);
    }
  }
  
  return *this;
}


AuxField& AuxField::product (const AuxField& a,
			     const AuxField& b)
// ---------------------------------------------------------------------------
// Set this AuxField's value as the product of a & b.
// Pseudospectrally in z, dealiased using 3/2 rule.
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim() == 2)
    Veclib::vmul (Geometry::nTot(), a.plane[0], 1, b.plane[0], 1, plane[0], 1);

  else {
    register int i, k;
    const int    nZ    = Geometry::nZ();
    const int    nP    = Geometry::nPlane();
    const int    nm    = nZ - 1;
    const int    nz32  = 3 * (nZ >> 1);
    const int    nzero = nz32 - nm;

    vector<real> work (4 * nz32 + 15);
    real         *at, *bt, *Wtab;

    at   = work();
    bt   = at + nz32;
    Wtab = bt + nz32;
    
    Femlib::rffti (nz32, Wtab);
    
    for (i = 0; i < nP; i++) {

      // -- Load dealiasing buffers at & bt, with zero padding at ends.
      
      at[0] = a.plane[0][i];
      bt[0] = b.plane[0][i];
      for (k = 1; k < nm; k++) {
	at[k] = a.plane[k + 1][i];
	bt[k] = b.plane[k + 1][i];
      }
      Veclib::zero (nzero, at + nm, 1);
      Veclib::zero (nzero, bt + nm, 1);

      // -- Multiply in physical space, transform back.
      
      Femlib::rfftb (nz32, at, Wtab);
      Femlib::rfftb (nz32, bt, Wtab);
      Veclib::vmul  (nz32, at, 1, bt, 1, at, 1);
      Femlib::rfftf (nz32, at, Wtab);

      // -- Copy transformed product back to current data area.

      plane[0][i] = at[0];
      plane[1][i] = 0.0;
      for (k = 1; k < nm; k++) plane[k + 1][i] = at[k];
    }

    // -- Normalize.
    
    Blas::scal (Geometry::nTot(), 1.0 / nz32, data, 1);
  }

  return *this;
}


AuxField& AuxField::addprod (const AuxField& a,
			     const AuxField& b)
// ---------------------------------------------------------------------------
// Add the product of a & b to this AuxField.
// Pseudospectrally in z, dealiased using 3/2 rule.
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim() == 2)
    Veclib::vvtvp (Geometry::nTot(), a.data, 1, b.data, 1, data, 1, data, 1);

  else {
    register int i, k;
    const int    nZ    = Geometry::nZ();
    const int    nP    = Geometry::nPlane();
    const int    nm    = nZ - 1;
    const int    nz32  = 3 * (nZ >> 1);
    const int    nzero = nz32 - nm;

    vector<real> work (4 * nz32 + 15);
    real         *at, *bt, *Wtab;

    at   = work();
    bt   = at + nz32;
    Wtab = bt + nz32;

    Femlib::rffti (nz32, Wtab);

    for (i = 0; i < nP; i++) {

      // -- Load dealiasing buffers at & bt, with zero padding at ends.

      at[0] = a.plane[0][i];
      bt[0] = b.plane[0][i];
      for (k = 1; k < nm; k++) {
	at[k] = a.plane[k + 1][i];
	bt[k] = b.plane[k + 1][i];
      }
      Veclib::zero (nzero, at + nm, 1);
      Veclib::zero (nzero, bt + nm, 1);

      // -- Multiply in physical space, transform back.

      Femlib::rfftb (nz32, at, Wtab);
      Femlib::rfftb (nz32, bt, Wtab);
      Veclib::vmul  (nz32, at, 1, bt, 1, at, 1);
      Femlib::rfftf (nz32, at, Wtab);
      
      // -- Normalize.
    
      Blas::scal (nZ, 1.0 / nz32, at, 1);

      // -- Add transformed product to current data area.

      plane[0][i] += at[0];
      for (k = 1; k < nm; k++) plane[k + 1][i] += at[k];
    }
  }

  return *this;
}


AuxField& AuxField::axpy (const real      alpha,
			  const AuxField& x    )
// ---------------------------------------------------------------------------
// Add alpha * x to this AuxField, plane-by-plane.
// ---------------------------------------------------------------------------
{
  register int k;
  const int    nZ = Geometry::nZ();
  const int    nP = Geometry::nPlane();
  
  for (k = 0; k < nZ; k++)
    Blas::axpy (nP, alpha, x.plane[k], 1, plane[k], 1);

  return *this;
}


AuxField& AuxField::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (Geometry::nTot(), data, 1, data, 1);

  return *this;
}


AuxField& AuxField::gradient (const int dir)
// ---------------------------------------------------------------------------
// Operate on AuxField to produce the nominated index of the gradient.
// dir == 0 ==> gradient in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
// ---------------------------------------------------------------------------
{
  char routine[] = "AuxField::gradient";

  register Element* E;
  register int      i, k, offset;
  const int         nE = Geometry::nElmt();
  const int         nZ = Geometry::nZ();
  const int         nP = Geometry::nPlane();

  switch (dir) {

  case 0:
    for (k = 0; k < nZ; k++) 
      for (i = 0; i < nE; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (plane[k] + offset, 0);
      }
    break;

  case 1:
    for (k = 0; k < nZ; k++) 
      for (i = 0; i < nE; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (0, plane[k] + offset);
      }
    break;

  case 2: {
    const int nmodes = Geometry::nMode();
    real      beta   = Femlib::value ("BETA");
    real*     tmp;

    for (k = 0; k < nmodes; k++) {
      tmp              = plane[2 * k];
      plane[2 * k]     = plane[2 * k + 1];
      plane[2 * k + 1] = tmp;
      Blas::scal (nP, -beta * k, plane[2 * k],     1);
      Blas::scal (nP,  beta * k, plane[2 * k + 1], 1);
    }
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }

  return *this;
}


AuxField& AuxField::errors (const Mesh& mesh    ,
			    const char* function)
// ---------------------------------------------------------------------------
// Compare F with function, print the infinity-norm Li, the 2-norm L2
// and the Sobolev 1-norm H1.
//
// The norms are found element-by-element, using projection onto higher-order
// elements and high-order quadrature.
//
// Warning: these routines only work in 2D at the moment.
// ---------------------------------------------------------------------------
{
  char routine[] = "AuxField::errors";

  if (!function) {
    message (routine, "empty function string", WARNING);
    return *this;
  }
  
  const int NQUAD = 15;

  Element*     E;
  Element*     P;
  real         area = 0.0;
  real         Li   = 0.0;
  real         L2   = 0.0;
  real         H1   = 0.0;
  vector<real> sol;
  vector<real> v;
  vector<real> tmp;
  const real   *z, **IN, **IT;
  real         *u;
  int          k, np, npp, ntot;
  const int    nE = Geometry::nElmt();

  Femlib::mesh (GLL, GLL, NQUAD, NQUAD, &z, 0, 0, 0, 0);

  for (k = 0; k < nE; k++) {
    E = Elmt[k];
    P = new Element (E -> ID(), mesh, z, NQUAD);

    u = data + E -> dOff();

    np   = E -> nKnot();
    npp  = P -> nKnot();
    ntot = P -> nTot ();

    tmp.setSize (np * npp);
    v  .setSize (ntot);
    sol.setSize (ntot);

    Femlib::mesh (GLL, GLL, np, npp, 0, &IN, &IT, 0, 0);

    Blas::mxm (*IN,   npp,  u,  np, tmp(), np );
    Blas::mxm (tmp(), npp, *IT, np, v(),   npp);

    P -> evaluate (function, sol());
    Veclib::vsub  (ntot, v(), 1, sol(), 1, v(), 1);

    area += P -> area ();
    Li    = max (Li, P -> norm_inf (v()));
    L2   += P -> norm_L2 (v());
    H1   += P -> norm_H1 (v());

    delete (P);
  }
  
  L2 /= area;
  H1 /= area;

  char  s[StrMax];
  ostrstream (s, StrMax) << "AuxField '"
                         << name()
                         << "' error norms (inf, L2, H1): "
			 << Li << "  " << L2 << "  " << H1 << ends;
  message ("", s, REMARK);

  return *this;
}


real AuxField::norm_inf () const
// ---------------------------------------------------------------------------
// Return infinity-norm (absolute max value) of AuxField data area.
// ---------------------------------------------------------------------------
{
  return fabs (data[Blas::iamax (Geometry::nTot(), data, 1)]);
}


real AuxField::mode_L2 (const int mode) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
// Mode numbers run 0 -- n_z/2 - 1.
// ---------------------------------------------------------------------------
{
  char         routine[] = "AuxField::mode_L2";
  real         area = 0.0, Ek = 0.0;
  register int j, offset;
  const int    nE = Geometry::nElmt();
  const int    nZ = Geometry::nZ();
  const int    kr = 2 * mode;
  const int    ki = 2 * mode + 1;
  Element*     E;
  
  if (mode < 0 ) message (routine, "negative mode number",        ERROR);
  if (ki   > nZ) message (routine, "mode number exceeds maximum", ERROR);

  for (j = 0; j < nE; j++) {
    E      = Elmt[j];
    offset = E -> dOff();
    area  += E -> area();
    Ek    += E -> norm_L2 (plane[kr] + offset);
    if (kr)
      Ek  += E -> norm_L2 (plane[ki] + offset);
  }

  return Ek / (2.0 * area);
}


ostream& operator << (ostream&  strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  strm.write ((char*) F.data, Geometry::nTot() * sizeof (real));

  return strm;
}


istream& operator >> (istream&  strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  strm.read ((char*) F.data, Geometry::nTot() * sizeof (real));

  return strm;
}


void AuxField::describe (char* s)  const
// ---------------------------------------------------------------------------
// Load s with a (prism-compatible) description of field geometry:
// NR NS NZ NEL.
// ---------------------------------------------------------------------------
{
  ostrstream (s, StrMax) << Geometry::nP()    << " "
                         << Geometry::nP()    << " "
                         << Geometry::nZ()    << " "
                         << Geometry::nElmt() << ends;
}


AuxField& AuxField::transform (const int sign)
// ---------------------------------------------------------------------------
// Discrete Fourier transform in homogeneous direction.  Number of points
// in that direction must be even, but is otherwise unrestricted. 
// Use sign = 1 for forward transform, -1 for inverse.
//
// Normalization is carried out on forward transform, so that the zeroth
// mode's real data are the average over the homogeneous direction of the
// physical space values.
// ---------------------------------------------------------------------------
{
  Femlib::DFTr (data              ,
		Geometry::nZ()    ,
		Geometry::nPlane(),
		1                 ,
		Geometry::nPlane(),
		sign              );
  
  return *this;
}


AuxField& AuxField::addToPlane (const int  k    ,
				const real alpha)
// ---------------------------------------------------------------------------
// Add in a constant to the values on nominated plane (if it exists),
// starting at plane zero.  You could call this a HACK.
// ---------------------------------------------------------------------------
{
  char routine[] = "AuxField::addToPlane";

  if (k < 0 || k >= Geometry::nZ())
    message (routine, "nominated plane doesn't exist", ERROR);
  else
    Veclib::sadd (Geometry::nPlane(), alpha, plane[k], 1, plane[k], 1);

  return *this;
}


void AuxField::swapData (AuxField* x,
			 AuxField* y)
// ---------------------------------------------------------------------------
// Swap data areas of two fields.
// ---------------------------------------------------------------------------
{
  register int   k;
  const int      nZ = Geometry::nZ();
  register real* tmp = x -> data;

  x -> data = y -> data;
  y -> data = tmp;

  for (k = 0; k < nZ; k++) {
    tmp           = x -> plane[k];
    x -> plane[k] = y -> plane[k];
    y -> plane[k] = tmp;
  }
}


void AuxField::couple (AuxField* v  ,
		       AuxField* w  ,
		       const int dir)
// ---------------------------------------------------------------------------
// Couples/uncouple field data for the radial and azimuthal
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
  register int k;
  const int    nZ    = Geometry::nZ(),
               nP    = Geometry::nPlane(),
               nMode = nZ >> 1;
  vector<real> work (nP);
  real         *Vr, *Vi, *Wr, *Wi, *tp;
  
  if (dir == 1) {

    for (k = 1; k < nMode; k++) {
      
      tp = work();

      Vr = v -> plane[2 * k];
      Vi = v -> plane[2 * k + 1];
      Wr = w -> plane[2 * k];
      Wi = w -> plane[2 * k + 1];

      Veclib::copy (nP, Vr, 1, tp, 1);
      Veclib::vsub (nP, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nP, Vi, 1, tp, 1);
      Veclib::vadd (nP, Vi, 1, Wr, 1, Vi, 1);
      Veclib::neg  (nP, Wr, 1);
      Veclib::vadd (nP, Wr, 1, tp, 1, Wr, 1);

      tp                    = w -> plane[2 * k];
      w -> plane[2 * k]     = w -> plane[2 * k + 1];
      w -> plane[2 * k + 1] = tp;
    }

  } else if (dir == -1) {

    for (k = 1; k < nMode; k++) {

      tp = work();

      Vr = v -> plane[2 * k];
      Vi = v -> plane[2 * k + 1];
      Wr = w -> plane[2 * k];
      Wi = w -> plane[2 * k + 1];

      Veclib::copy (nP, Vr, 1, tp, 1);
      Veclib::vadd (nP, Vr, 1, Wr, 1, Vr, 1);
      Veclib::vsub (nP, Wr, 1, tp, 1, Wr, 1);
      Veclib::copy (nP, Vi, 1, tp, 1);
      Veclib::vadd (nP, Vi, 1, Wi, 1, Vi, 1);
      Veclib::neg  (nP, Wi, 1);
      Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);

      tp                    = w -> plane[2 * k];
      w -> plane[2 * k]     = w -> plane[2 * k + 1];
      w -> plane[2 * k + 1] = tp;
    }

    Blas::scal (nP * (nZ - 2), 0.5, v -> data + 2 * nP, 1);
    Blas::scal (nP * (nZ - 2), 0.5, w -> data + 2 * nP, 1);

  } else
    message (routine, "unknown direction given", ERROR);
}


AuxField& AuxField::divR ()
// ---------------------------------------------------------------------------
// Divide data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      i, k, offset = 0;
  register real*    p;
  const int         nE = Geometry::nElmt();
  const int         nZ = Geometry::nZ();

  for (k = 0; k < nZ; k++) {
    p = plane[k];
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();
      E -> divR (p + offset);
    }
  }
  
  return *this;
}


AuxField& AuxField::mulR ()
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      i, k, offset = 0;
  register real*    p;
  const int         nE = Geometry::nElmt();
  const int         nZ = Geometry::nZ();

  for (k = 0; k < nZ; k++) {
    p = plane[k];
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();
      E -> mulR (p + offset);
    }
  }
  
  return *this;
}
