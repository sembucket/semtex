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
// Fourier mode, but are kept zero and never evolve.  The planes always
// point to the same storage locations within the data area.
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
  const char       routine[] = "AuxField::AuxField";
  register integer k;
  const integer    nZ = Geometry::nZ();
  const integer    nP = Geometry::planeSize();
  const integer    nT = Geometry::nTotal();

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
  if   (val == 0.0) Veclib::zero (Geometry::nTotal(),      data, 1);
  else              Veclib::fill (Geometry::nTotal(), val, data, 1);

  return *this;
}


AuxField& AuxField::operator += (const real val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (Geometry::nTotal(), val, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const real val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (Geometry::nTotal(), -val, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator *= (const real val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (Geometry::nTotal(),      data, 1);
  else              Blas  ::scal (Geometry::nTotal(), val, data, 1);

  return *this;
}


AuxField& AuxField::operator /= (const real val)
// ---------------------------------------------------------------------------
// Divide field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) message ("AuxField::op /= real", "divide by zero", ERROR);
  else              Blas::scal (Geometry::nTotal(), 1.0 / val, data, 1);

  return *this;
}


AuxField& AuxField::operator = (const AuxField& f)
// ---------------------------------------------------------------------------
// This presumes the two fields conform, and copies f's value storage to LHS.
// ---------------------------------------------------------------------------
{
  Veclib::copy (Geometry::nTotal(), f.data, 1, data, 1);
  
  return *this;
}


AuxField& AuxField::operator += (const AuxField& f)
// ---------------------------------------------------------------------------
// Add f's value to this AuxField's.
// ---------------------------------------------------------------------------
{
  Veclib::vadd (Geometry::nTotal(), data, 1, f.data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const AuxField& f)
// ---------------------------------------------------------------------------
// Subtract f's value from this AuxField's.
// ---------------------------------------------------------------------------
{
  Veclib::vsub (Geometry::nTotal(), data, 1, f.data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.  Physical space.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register integer  i, k, offset = 0;
  const integer     nE = Geometry::nElmt();
  const integer     nZ = Geometry::nZ();
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
    Veclib::vmul (Geometry::nTotal(), a.data, 1, b.data, 1, data, 1);

  else {
    register integer i, k;
    const integer    nZ    = Geometry::nZ();
    const integer    nP    = Geometry::nPlane();
    const integer    nm    = nZ - 1;
    const integer    nz32  = 3 * (nZ >> 1);
    const integer    nzero = nz32 - nm;

    vector<real>     work (4 * nz32 + 15);
    real             *at, *bt, *Wtab;

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
    
    Blas::scal (Geometry::nTotal(), 1.0 / nz32, data, 1);
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
    Veclib::vvtvp (Geometry::nTotal(), a.data, 1, b.data, 1, data, 1, data, 1);

  else {
    register integer i, k;
    const integer    nZ    = Geometry::nZ();
    const integer    nP    = Geometry::nPlane();
    const integer    nm    = nZ - 1;
    const integer    nz32  = 3 * (nZ >> 1);
    const integer    nzero = nz32 - nm;

    vector<real>     work (4 * nz32 + 15);
    real             *at, *bt, *Wtab;

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


AuxField& AuxField::innerProduct (const vector <AuxField*>& a,
                                  const vector <AuxField*>& b)
// ---------------------------------------------------------------------------
// Set this AuxField's value as the inner product of a & b
// in physical space --- don't worry about dealiasing.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer ntot = Geometry::nTotal();
  const integer ndim = Geometry::nDim();
  
  Veclib::zero (ntot, data, 1);

  for (i = 0; i < ndim; i++)
    Veclib::vvtvp (ntot, a[i] -> data, 1, b[i] -> data, 1, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::times (const AuxField& a,
			   const AuxField& b)
// ---------------------------------------------------------------------------
// Set this AuxField equal to the product of a & b (in physical space).
// No dealiasing.
// ---------------------------------------------------------------------------
{
  const int ntot = Geometry::nTotal();

  Veclib::vmul (ntot, a.data, 1, b.data, 1, data, 1);

  return *this;
}


AuxField& AuxField::timesPlus (const AuxField& a,
			       const AuxField& b)
// ---------------------------------------------------------------------------
// Add the product of a & b to this AuxField (in physical space).
// No dealiasing.
// ---------------------------------------------------------------------------
{
  Veclib::vvtvp (Geometry::nTotal(), a.data, 1, b.data, 1, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::timesMinus (const AuxField& a,
			        const AuxField& b)
// ---------------------------------------------------------------------------
// Subtract the product of a & b from this AuxField (in physical space).
// No dealiasing.
// ---------------------------------------------------------------------------
{
  Veclib::vvvtm (Geometry::nTotal(), data, 1, a.data, 1, b.data, 1, data, 1);

  return *this;
}


AuxField& AuxField::axpy (const real      alpha,
			  const AuxField& x    )
// ---------------------------------------------------------------------------
// Add alpha * x to this AuxField.
// ---------------------------------------------------------------------------
{
  Blas::axpy (Geometry::nTotal(), alpha, x.data, 1, data, 1);

  return *this;
}


AuxField& AuxField::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (Geometry::nTotal(), data, 1, data, 1);

  return *this;
}


AuxField& AuxField::gradient (const integer dir)
// ---------------------------------------------------------------------------
// Operate on AuxField to produce the nominated index of the gradient.
// dir == 0 ==> gradient in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
// ---------------------------------------------------------------------------
{
  const char        routine[] = "AuxField::gradient";
  register Element* E;
  register integer  i, k, offset;
  const integer     nE = Geometry::nElmt();
  const integer     nZ = Geometry::nZ();
  const integer     nP = Geometry::planeSize();
  vector<real>      tmp (2*Geometry::nTotElmt());
  register real*    work = tmp();
  const real        **DV, **DT;

  Femlib::quad (LL, Geometry::nP(), Geometry::nP(), 0, 0, 0, 0, 0, &DV, &DT);

  switch (dir) {

  case 0:
    for (k = 0; k < nZ; k++) 
      for (i = 0; i < nE; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (plane[k] + offset, 0, DV, DT, work);
      }
    break;

  case 1:
    for (k = 0; k < nZ; k++) 
      for (i = 0; i < nE; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (0, plane[k] + offset, DV, DT, work);
      }
    break;

  case 2: {
    integer       Re, Im;
    const integer nmodes = Geometry::nMode();
    const real    beta   = Femlib::value ("BETA");
    vector<real>  work (nP);
    real*         tmp = work();

    Veclib::zero (2 * nP, data, 1); // -- Zero real & Nyquist planes.
    for (k = 1; k < nmodes; k++) {
      Re = k  + k;
      Im = Re + 1;
      Veclib::copy (nP,             plane[Re], 1, tmp,       1);
      Veclib::smul (nP, -beta * k,  plane[Im], 1, plane[Re], 1);
      Veclib::smul (nP,  beta * k,  tmp,       1, plane[Im], 1);
    }
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }

  return *this;
}


void AuxField::gradient (const integer nZ ,
			 real*         src,
			 const integer dir) const
// ---------------------------------------------------------------------------
// Use Field structure to perform gradient operations on data area src,
// according to nominated direction.  Input value nZ is the number of planes
// to operate on.
// ---------------------------------------------------------------------------
{
  const char        routine[] = "AuxField::gradient";
  register Element* E;
  register integer  i, k, offset;
  const integer     nE = Geometry::nElmt();
  const integer     nP = Geometry::planeSize();
  vector<real>      tmp (2*Geometry::nTotElmt());
  register real*    work = tmp();
  const real        **DV, **DT;

  Femlib::quad (LL, Geometry::nP(), Geometry::nP(), 0, 0, 0, 0, 0, &DV, &DT);

  switch (dir) {

  case 0:
    for (k = 0; k < nZ; k++) 
      for (i = 0; i < nE; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (src + k * nP + offset, 0, DV, DT, work);
      }
    break;

  case 1:
    for (k = 0; k < nZ; k++) 
      for (i = 0; i < nE; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (0, src + k * nP + offset, DV, DT, work);
      }
    break;

  case 2: {
    const integer nmodes = Geometry::nMode();
    const real    beta   = Femlib::value ("BETA");
    vector<real>  work (nP);
    real          *Re, *Im, *tmp = work();

    Veclib::zero (2 * nP, src, 1); // -- Zero real & Nyquist planes.
    for (k = 1; k < nmodes; k++) {
      Re = src + 2 * k * nP;
      Im = Re  + nP;
      Veclib::copy (nP,             Re,  1, tmp, 1);
      Veclib::smul (nP, -beta * k,  Im,  1, Re,  1);
      Veclib::smul (nP,  beta * k,  tmp, 1, Im,  1);
    }
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;

  }
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
  const char routine[] = "AuxField::errors";

  if (!function) {
    message (routine, "empty function string", WARNING);
    return *this;
  }
  
  const integer NQUAD = 15;

  Element*      E;
  Element*      P;
  real          area = 0.0;
  real          Li   = 0.0;
  real          L2   = 0.0;
  real          H1   = 0.0;
  vector<real>  sol;
  vector<real>  v;
  vector<real>  tmp;
  const real    *z, **IN, **IT;
  real          *u;
  integer       k, np, npp, ntot;
  const integer nE = Geometry::nElmt();

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
  return fabs (data[Blas::iamax (Geometry::nTotal(), data, 1)]);
}


real AuxField::mode_L2 (const integer mode) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
// Mode numbers run 0 -- n_z/2 - 1.
// ---------------------------------------------------------------------------
{
  const char        routine[] = "AuxField::mode_L2";
  real              area = 0.0, Ek = 0.0;
  register integer  j, offset;
  const integer     nE = Geometry::nElmt();
  const integer     nZ = Geometry::nZ();
  const integer     kr = 2 * mode;
  const integer     ki = 2 * mode + 1;
  register Element* E;
  
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
  register integer i;
  const integer    nZ = Geometry::nZ();
  const integer    nP = Geometry::nPlane();
  
  for (i = 0; i < nZ; i++) strm.write ((char*) F.plane[i], nP * sizeof (real));

  return strm;
}


istream& operator >> (istream&  strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.  Zero any blank storage areas.
// ---------------------------------------------------------------------------
{
  register integer i;
  const integer    nZ = Geometry::nZ();
  const integer    nP = Geometry::nPlane();
  const integer    NP = Geometry::planeSize();
  
  for (i = 0; i < nZ; i++) strm.read ((char*) F.plane[i], nP * sizeof (real));

  if (NP > nP) Veclib::zero (nZ, F.data + nP, NP);

  return strm;
}


AuxField& AuxField::zeroNyquist ()
// ---------------------------------------------------------------------------
// Set storage for highest frequency mode to zero.  This mode is carried
// but never evolves.  This mode is stored as the second data plane.
// ---------------------------------------------------------------------------
{
  const integer nZ = Geometry::nZ();
  const integer nP = Geometry::planeSize();

  if (nZ > 1) Veclib::zero (nP, plane[1], 1);

  return *this;
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


AuxField& AuxField::transform (const integer sign)
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
  const integer nZ = Geometry::nZ();
  const integer nP = Geometry::planeSize();

  if (nZ > 1)
    if (nZ == 2)
      if   (sign == +1) Veclib::zero (nP, plane[1], 1);
      else              Veclib::copy (nP, plane[0], 1, plane[1], 1);
    else
      Femlib::DFTr (data, nZ, nP, sign);
  
  return *this;
}


AuxField& AuxField::transform32 (real*         phys,
				 const integer sign)
// ---------------------------------------------------------------------------
// Discrete Fourier transform in homogeneous direction, extended for
// dealiasing.  Input pointer phys points to data in physical space, 
// which acts as input area if sign == +1, output area if sign == -1.
// So transform is from phys to internal storage if sign == +1 and
// vice versa.
// ---------------------------------------------------------------------------
{
  const integer nZ     = Geometry::nZ();
  const integer nP     = Geometry::planeSize();
  const integer nTot   = nZ * nP;
  const integer nZ32   = (3 * nZ) >> 1;
  const integer nTot32 = nZ32 * nP;

  if (nZ <= 2) {
    if (sign == +1)
      Veclib::copy (nTot, phys, 1, data, 1);
    else
      Veclib::copy (nTot, data, 1, phys, 1);

  } else {
    if (sign == +1) {
      Femlib::DFTr (phys, nZ32, nP, +1);
      Veclib::copy (nTot, phys, 1, data, 1);
    } else {
      Veclib::copy (nTot, data, 1, phys, 1);
      Veclib::zero (nTot32 - nTot, phys + nTot, 1);
      Femlib::DFTr (phys, nZ32, nP, -1);
    }
  }

  return *this;
}


AuxField& AuxField::addToPlane (const integer k    ,
				const real    alpha)
// ---------------------------------------------------------------------------
// Add in a constant to the values on nominated plane (if it exists),
// starting at plane zero.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::addToPlane";

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
  register integer k;
  const integer    nZ = Geometry::nZ();
  register real*   tmp;

  tmp       = x -> data;
  x -> data = y -> data;
  y -> data = tmp;

  for (k = 0; k < nZ; k++) {
    tmp           = x -> plane[k];
    x -> plane[k] = y -> plane[k];
    y -> plane[k] = tmp;
  }
}


void AuxField::couple (AuxField*     v  ,
		       AuxField*     w  ,
		       const integer dir)
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

  const char       routine[] = "Field::couple";
  register integer k, Re, Im;
  const integer    nZ    = Geometry::nZ(),
                   nP    = Geometry::planeSize(),
                   nMode = nZ >> 1;
  vector<real>     work (nP);
  real             *Vr, *Vi, *Wr, *Wi, *tp = work();
  
  if (dir == 1) {

    for (k = 1; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> plane[Re];
      Vi = v -> plane[Im];
      Wr = w -> plane[Re];
      Wi = w -> plane[Im];

      Veclib::copy (nP, Vr, 1, tp, 1);
      Veclib::vsub (nP, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nP, Wr, 1, tp, 1);
      Veclib::copy (nP, Wi, 1, Wr, 1);
      Veclib::vsub (nP, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nP, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == -1) {

    for (k = 1; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> plane[Re];
      Vi = v -> plane[Im];
      Wr = w -> plane[Re];
      Wi = w -> plane[Im];

      Veclib::copy  (nP,      Vr, 1, tp, 1);
      Veclib::svvpt (nP, 0.5, Vr, 1, Wr, 1, Vr, 1);
      Veclib::svvmt (nP, 0.5, Wr, 1, tp, 1, Wr, 1);
      Veclib::copy  (nP,      Wi, 1, tp, 1);
      Veclib::copy  (nP,      Wr, 1, Wi, 1);
      Veclib::svvmt (nP, 0.5, Vi, 1, tp, 1, Wr, 1);
      Veclib::svvpt (nP, 0.5, Vi, 1, tp, 1, Vi, 1);
    }

  } else
    message (routine, "unknown direction given", ERROR);
}


AuxField& AuxField::divR ()
// ---------------------------------------------------------------------------
// Divide data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  register Element* E;
  register integer  i, k, offset = 0;
  register real*    p;
  const integer     nE = Geometry::nElmt();
  const integer     nZ = Geometry::nZ();

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


void AuxField::divR (const integer nZ ,
		     real*         src) const
// ---------------------------------------------------------------------------
// Divide src by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  register Element* E;
  register integer  i, k, offset = 0;
  const integer     nE = Geometry::nElmt();
  const integer     nP = Geometry::planeSize();

  for (k = 0; k < nZ; k++) {
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();
      E -> divR (src + k * nP + offset);
    }
  }
}


AuxField& AuxField::mulR ()
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  register Element* E;
  register integer  i, k, offset = 0;
  register real*    p;
  const integer     nE = Geometry::nElmt();
  const integer     nZ = Geometry::nZ();

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


void AuxField::mulR (const integer nZ ,
		     real*         src) const
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  register Element* E;
  register integer  i, k, offset = 0;
  const integer     nE = Geometry::nElmt();
  const integer     nP = Geometry::planeSize();

  for (k = 0; k < nZ; k++) {
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();
      E -> mulR (src + k * nP + offset);
    }
  }
}


real AuxField::probe (const Element* E,
		      const real     r,
		      const real     s,
		      const integer  k) const
// ---------------------------------------------------------------------------
// Return the value of data on plane k, in Element E, location r, s.
// ---------------------------------------------------------------------------
{
  const integer offset = E -> dOff();
  
  return E -> probe (r, s, plane[k] + offset);
}


real AuxField::probe (const Element* E,
		      const real     r,
		      const real     s,
		      const real     z) const
// ---------------------------------------------------------------------------
// Return the value of data, in Element E, location r, s, z.
//
// NB: interpolation assumes that AuxField is Fourier transformed.
// ---------------------------------------------------------------------------
{
  register integer k, Re, Im;
  register real    value, phase;
  const integer    NZ     = Geometry::nZ();
  const integer    NZH    = NZ >> 1;
  const integer    NHM    = NZH - 1;
  const integer    offset = E -> dOff();
  const real       betaZ  = z * Femlib::value ("BETA");
  vector<real>     work (NZ);
  register real*   fbuf = work();

  if (NZ < 3)
    value = E -> probe (r, s, plane[0] + offset);
  
  else {
    for (k = 0; k < NZ; k++)
      fbuf[k] = E -> probe (r, s, plane[k] + offset);

    Blas::scal (NZ - 2, 2.0, fbuf + 2, 1);

    value  = fbuf[0];
    value += fbuf[1] * cos (NZH * betaZ);
    for (k = 1; k < NHM; k++) {
      Re     = k  + k;
      Im     = Re + 1;
      phase  = k * betaZ;
      value += fbuf[Re] * cos (phase) - fbuf[Im] * sin (phase);
    }
  }
   
  return value;
}


AuxField& AuxField::Smagorinsky ()
// ---------------------------------------------------------------------------
// Given that *this is the strain-rate magnitude S, compute the Smagorinsky
// eddy-viscosity
//                      \nu_T = (Cs delta)^2 S.
// ---------------------------------------------------------------------------
{
  register integer  k, offset;
  register Element* E;
  const integer     nZ = Geometry::nZ();
  const integer     nE = Geometry::nElmt();
  const integer     nP = Geometry::nPlane();
  const real        Cs = Femlib::value ("C_SMAG");
  vector<real>      work (nP);

  for (k = 0; k < nE; k++) {
    E      = Elmt[k];
    offset = E -> dOff();
    E -> lengthScale (work() + offset);
  }
  
  Blas::scal   (nP, Cs, work(), 1);
  Veclib::vmul (nP, work(), 1, work(), 1, work(), 1);

  for (k = 0; k < nZ; k++) {
    if (k == 1) continue;
    Veclib::vmul (nP, work(), 1, plane[k], 1, plane[k], 1);
  }

  return *this;
}


real AuxField::CFL (const integer dir) const
// ---------------------------------------------------------------------------
// Return the inverse CFL timescale using this AuxField as a velocity 
// component in the nominated direction.  Computations only occur on the
// zeroth Fourier mode.
// dir == 0 ==> CFL estimate in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
// ---------------------------------------------------------------------------
{
  const char        routine[] = "AuxField::CFL";
  register Element* E;
  register integer  i, offset;
  const integer     nE = Geometry::nElmt();
  real              dxy, CFL = 0.0;
  
  {
    const integer nP = Geometry::nP();
    const real*   z;
    Femlib::quad (LL, nP, nP, &z, 0, 0, 0, 0, 0, 0);
    dxy = z[1] - z[0];
  }

  switch (dir) {
  case 0:
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();
      CFL    = max (CFL, E -> CFL (dxy, data + offset, 0));
      }
    break;
  case 1:
    for (i = 0; i < nE; i++) {
      E      = Elmt[i];
      offset = E -> dOff();
      CFL    = max (CFL, E -> CFL (dxy, 0, data + offset));
    }
    break;
  case 2: {
    const integer nP = Geometry::nPlane();
    const real    dz = Femlib::value ("TWOPI / (BETA * N_Z)");
    for (i = 0; i < nP; i++) CFL = max (CFL, fabs (data[i]));
    CFL /= dz;
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }

  return CFL;
}
