///////////////////////////////////////////////////////////////////////////////
// auxfield.C:  routines for AuxField class, including Fourier expansions.
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
		    const int         nz  ,
		    const char        name) :
		    
		    field_name       (name), 
		    Elmt             (Elts),
		    n_z              (nz  ),
		    n_plane          (0   ),
		    n_elmt_bnodes    (0   ),
		    n_elmt_inodes    (0   ),
		    elmt_np_max      (0   ),
		    elmt_nt_max      (0   ),
		    elmt_ne_max      (0   ),
		    elmt_ni_max      (0   )
// ---------------------------------------------------------------------------
// Allocate field storage area and integer size records.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      k;
  const int         N = Elmt.getSize();

  k = 0;
  for (k = 0; k < N; k++) {
    E = Elmt[k];

    n_plane       += E -> nTot();
    n_elmt_bnodes += E -> nExt();

    elmt_np_max = max (elmt_np_max, E -> nKnot());
    elmt_nt_max = max (elmt_nt_max, E -> nTot ());
    elmt_ne_max = max (elmt_ne_max, E -> nExt ());
    elmt_ni_max = max (elmt_ne_max, E -> nInt ());
  }
  n_elmt_inodes = n_plane - n_elmt_bnodes;

  data  = new real  [n_z * n_plane];
  plane = new real* [n_z];

  for (k = 0; k < n_z; k++) plane[k] = data + k * n_plane;
}


AuxField& AuxField::operator = (const real val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (nTot(),      data, 1);
  else              Veclib::fill (nTot(), val, data, 1);

  return *this;
}


AuxField& AuxField::operator += (const real val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (nTot(), val, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const real val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (nTot(), -val, data, 1, data, 1);

  return *this;
}


AuxField& AuxField::operator *= (const real val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (nTot(),      data, 1);
  else              Blas::scal   (nTot(), val, data, 1);

  return *this;
}


AuxField& AuxField::operator /= (const real val)
// ---------------------------------------------------------------------------
// Divide field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) message ("AuxField::op /= real", "divide by zero", ERROR);
  else              Blas::scal (nTot(), 1.0 / val, data, 1);

  return *this;
}


AuxField& AuxField::operator = (const AuxField& f)
// ---------------------------------------------------------------------------
// This presumes the two fields conform, and copies f's value storage to LHS.
// ---------------------------------------------------------------------------
{
  if (f.nTot() != nTot())
    message ("AuxField::operator = const AuxField&", "size mismatch", ERROR);

  register int k;
  
  for (k = 0; k < n_z; k++)
    Veclib::copy (nPlane(), f.plane[k], 1, plane[k], 1);
  
  return *this;
}


AuxField& AuxField::operator += (const AuxField& f)
// ---------------------------------------------------------------------------
// Add f's value to this AuxField's.
// ---------------------------------------------------------------------------
{
  if (f.nTot() != nTot())
    message ("AuxField::operator += const AuxField&", "size mismatch", ERROR);

  register int k;
  
  for (k = 0; k < n_z; k++)
    Veclib::vadd (nPlane(), plane[k], 1, f.plane[k], 1, plane[k], 1);

  return *this;
}


AuxField& AuxField::operator -= (const AuxField& f)
// ---------------------------------------------------------------------------
// Subtract f's value from this AuxField's.
// ---------------------------------------------------------------------------
{
  if (f.nTot() != nTot())
    message ("AuxField::operator -= const AuxField&", "size mismatch", ERROR);

  register int k;
  
  for (k = 0; k < n_z; k++)  
    Veclib::vsub (nPlane(), plane[k], 1, f.plane[k], 1, plane[k], 1);

  return *this;
}


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.
// ---------------------------------------------------------------------------
{
  register Element* E;
  register int      i, k, offset = 0;
  const int         ne = nEl();
  const real        dz = Femlib::value ("TWOPI/BETA") / n_z;
  real              z;

  for (k = 0; k < n_z; k++) {
    z = k * dz;
    Femlib::value ("z", z);
    for (i = 0; i < ne; i++) {
      E = Elmt[i];
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
  if (a.nTot() != nTot() || b.nTot() != nTot())
    message ("AuxField::product", "size mismatch", ERROR);

  if (n_z == 1)
    Veclib::vmul (nTot(), a.plane[0], 1, b.plane[0], 1, plane[0], 1);

  else {
    register int i, k;
    const int    nm    = n_z - 1;
    const int    nz32  = 3 * (n_z >> 1);
    const int    nzero = nz32 - nm;

    vector<real> work (4 * nz32 + 15);
    real         *at, *bt, *Wtab;

    at   = work();
    bt   = at + nz32;
    Wtab = bt + nz32;
    
    Femlib::rffti (nz32, Wtab);
    
    for (i = 0; i < n_plane; i++) {

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
    
    Blas::scal (nTot(), 1.0 / nz32, data, 1);
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
  if (a.nTot() != nTot() || b.nTot() != nTot())
    message ("AuxField::addprod", "size mismatch", ERROR);

  if (n_z == 1)
    Veclib::vvtvp (nTot(), a.data, 1, b.data, 1, data, 1, data, 1);

  else {
    register int i, k;
    const int    nm    = n_z - 1;
    const int    nz32  = 3 * (n_z >> 1);
    const int    nzero = nz32 - nm;

    vector<real> work (4 * nz32 + 15);
    real         *at, *bt, *Wtab;

    at   = work();
    bt   = at + nz32;
    Wtab = bt + nz32;

    Femlib::rffti (nz32, Wtab);

    for (i = 0; i < n_plane; i++) {

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
    
      Blas::scal (n_z, 1.0 / nz32, at, 1);

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
  if (x.nTot() != nTot())
    message ("AuxField::axpy", "size mismatch", ERROR);
  
  register int k;
  
  for (k = 0; k < n_z; k++)
    Blas::axpy (n_plane, alpha, x.plane[k], 1, plane[k], 1);

  return *this;
}


AuxField& AuxField::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (nTot(), data, 1, data, 1);

  return *this;
}


void AuxField::swapData (AuxField* x,
			 AuxField* y)
// ---------------------------------------------------------------------------
// Swap data areas of two fields.
// ---------------------------------------------------------------------------
{
  if (x -> nTot() != y -> nTot())
    message ("AuxField::swapData", "size mismatch", ERROR);

  register int   k;
  register real* tmp = x -> data;

  x -> data = y -> data;
  y -> data = tmp;

  for (k = 0; k < x -> n_z; k++) {
    tmp           = x -> plane[k];
    x -> plane[k] = y -> plane[k];
    y -> plane[k] = tmp;
  }
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
  const int         N  = nEl();

  switch (dir) {

  case 0:
    for (k = 0; k < n_z; k++) 
      for (i = 0; i < N; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (plane[k] + offset, 0);
      }
    break;

  case 1:
    for (k = 0; k < n_z; k++) 
      for (i = 0; i < N; i++) {
	E      = Elmt[i];
	offset = E -> dOff();

	E -> grad (0, plane[k] + offset);
      }
    break;

  case 2: {
    const int nmodes = n_z >> 1;
    real      beta   = Femlib::value ("BETA");
    real*     tmp;

    for (k = 0; k < nmodes; k++) {
      tmp              = plane[2 * k];
      plane[2 * k]     = plane[2 * k + 1];
      plane[2 * k + 1] = tmp;
      Blas::scal (n_plane, -beta * k, plane[2 * k],     1);
      Blas::scal (n_plane,  beta * k, plane[2 * k + 1], 1);
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
  real         *u, *z, **IN, **IT;
  int          k, np, npp, ntot;
  const int    N = nEl();

  Femlib::mesh (GLL, GLL, NQUAD, NQUAD, &z, 0, 0, 0, 0);

  for (k = 0; k < N; k++) {
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


ostream& operator << (ostream&  strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  strm.write ((char*) F.data, F.nTot() * sizeof (real));

  return strm;
}


istream& operator >> (istream&  strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  strm.read ((char*) F.data, F.nTot() * sizeof (real));

  return strm;
}


void AuxField::describe (char* s)  const
// ---------------------------------------------------------------------------
// Load s with a (prism-compatible) description of field geometry:
// NR NS NZ NEL.
// ---------------------------------------------------------------------------
{
  ostrstream (s, StrMax) << elmt_np_max << " "
                         << elmt_np_max << " "
                         << n_z         << " "
                         << nEl()       << ends;
}


AuxField& AuxField::transform (const int sign)
// ---------------------------------------------------------------------------
// Discrete Fourier transform in homogeneous direction.  Number of points
// in that direction must be even, but is otherwise unrestricted. 
// Use sign = 1 for forward transform, -1 for inverse.
//
// Much of the complication is caused by the "non-standard" data
// format produced by FFTPACK's rfftf.  Here the (Nyquist, real) datum gets
// stored as the imaginary part of the zeroth mode, while for rfftf it
// gets stored as the last point.  The Nyquist datum is always set to
// zero below.
//
// Normalization is carried out on forward transform, so that the zeroth
// mode's real data are the average over the homogeneous direction of the
// physical space values.
// ---------------------------------------------------------------------------
{
  if (n_z < 2) return *this;

  register int i;
  vector<real> work (3 * n_z + 15);
  real*        tmp  = work();
  real*        Wtab = tmp + n_z;
  real*        ptr  = data;

  Femlib::rffti (n_z, Wtab);

  switch (sign) {

  case 1:
    for (i = 0; i < n_plane; i++, ptr++) {
      Veclib::copy  (n_z, ptr, n_plane, tmp, 1);
      Femlib::rfftf (n_z, tmp, Wtab);
      Veclib::copy  (n_z - 2, tmp + 1, 1, ptr + 2 * n_plane, n_plane);
      ptr[0]       = tmp[0];
      ptr[n_plane] = 0.0;
    }
    Blas::scal (nTot(), 1.0 / n_z, data, 1);
    break;

  case -1:
    for (i = 0; i < n_plane; i++, ptr++) {
      tmp[n_z - 1] = 0.0;
      tmp[0]       = ptr[0];
      Veclib::copy  (n_z - 2, ptr + 2 * n_plane, n_plane, tmp + 1, 1);
      Femlib::rfftb (n_z, tmp, Wtab);
      Veclib::copy  (n_z, tmp, 1, ptr, n_plane);
    }
    break;

  default:
    message ("AuxField::transform", "illegal direction flag", ERROR);
    break;
  }
  
  return *this;
}


AuxField& AuxField::addToPlane (const int  k    ,
				const real alpha)
// ---------------------------------------------------------------------------
// Add in a constant to the values on nominated plane (if it exists),
// starting at plane zero.
// ---------------------------------------------------------------------------
{
  char routine[] = "AuxField::addToPlane";

  if (k < 0 || k >= n_z)
    message (routine, "nominated plane doesn't exist", ERROR);
  else
    Veclib::sadd (n_plane, alpha, plane[k], 1, plane[k], 1);

  return *this;
}

