//////////////////////////////////////////////////////////////////////////////
// auxfield.C: routines for AuxField class, including Fourier expansions.
//
// Copyright (c) <--> $Date$, Hugh Blackburn
//
// NB: this version is modified to suit dual.  For dual, data are NOT
// transformed to physical space for storage in restart files.  There
// are 3 planes of data: the first corresponds to mode 0 (which is
// real only), the second and third to the real an imaginary parts of
// the selected Fourier mode.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "sem.h"


AuxField::AuxField (real_t*           alloc,
		    const int_t       nz   ,
		    vector<Element*>& elmt ,
		    const char        name ) :
// ---------------------------------------------------------------------------
// Install field storage area and integer size records.
// ---------------------------------------------------------------------------
  _name (name),
  _elmt (elmt),
  _nz   (nz),
  _size (nz * Geometry::planeSize()),
  _data (alloc)
{
  const char     routine[] = "AuxField::AuxField";
  const int_t    nP = Geometry::planeSize();
  register int_t k;

  if (Geometry::nElmt() != _elmt.size())
    message (routine, "conflicting number of elements in input data", ERROR);

  _plane = new real_t* [static_cast<size_t>(_nz)];

  for (k = 0; k < _nz; k++) _plane[k] = _data + k * nP;
}


AuxField& AuxField::setInput (real_t*     alloc,
			      const int_t nz   )
// ---------------------------------------------------------------------------
// Install a new lot of field storage, with nominated number of planes.
// It is assumed that alloc is at least nZ * geometry::planeSize() long.
// ---------------------------------------------------------------------------
{
  register int_t k;
  const int_t    nP = Geometry::planeSize();

  _nz   = nz;
  _size = _nz * nP;
  _data = alloc;

  delete [] _plane;
  _plane = new real_t* [(size_t) _nz];

  for (k = 0; k < _nz; k++) _plane[k] = _data + k * nP;

  return *this;
}


AuxField& AuxField::operator = (const real_t val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (_size,      _data, 1);
  else              Veclib::fill (_size, val, _data, 1);

  return *this;
}


AuxField& AuxField::operator += (const real_t val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (_size, val, _data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const real_t val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (_size, -val, _data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator *= (const real_t val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (_size,      _data, 1);
  else              Blas  ::scal (_size, val, _data, 1);

  return *this;
}


AuxField& AuxField::operator /= (const real_t val)
// ---------------------------------------------------------------------------
// Divide field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) message ("AuxField::op /= real_t", "divide by zero", ERROR);
  else              Blas::scal (_size, 1.0 / val, _data, 1);

  return *this;
}


AuxField& AuxField::operator = (const AuxField& f)
// ---------------------------------------------------------------------------
// This presumes the two fields conform, and copies f's value storage to LHS.
// ---------------------------------------------------------------------------
{
  Veclib::copy (_size, f._data, 1, _data, 1);
  
  return *this;
}


AuxField& AuxField::operator += (const AuxField& f)
// ---------------------------------------------------------------------------
// Add f's value to this AuxField's.
// ---------------------------------------------------------------------------
{
  Veclib::vadd (_size, _data, 1, f._data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const AuxField& f)
// ---------------------------------------------------------------------------
// Subtract f's value from this AuxField's.
// ---------------------------------------------------------------------------
{
  Veclib::vsub (_size, _data, 1, f._data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator *= (const AuxField& f)
// ---------------------------------------------------------------------------
// Multiply *this storage vectorwise with f's.  You sort out which
// space you're in!
// ---------------------------------------------------------------------------
{
  Veclib::vmul (_size, _data, 1, f._data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.  Physical space.
//
// NB: function should not have z as a variable for dual: enforce z=0.
// ---------------------------------------------------------------------------
{
  const int_t    nel = Geometry::nElmt();
  const int_t    np2 = Geometry::nTotElmt();
  register int_t i, k;
  real_t*            p;

  Femlib::value ("z", 0.0);

  for (k = 0; k < _nz; k++)
    for (p = _plane[k], i = 0; i < nel; i++, p += np2)
      _elmt[i] -> evaluate (function, p);
  
  return *this;
}


AuxField& AuxField::times (const AuxField& a,
			   const AuxField& b)
// ---------------------------------------------------------------------------
// Set this AuxField equal to the product of a & b (in physical space).
// No dealiasing.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::times";

  if (_size != a._size || _size != b._size)
    message (routine, "non-congruent inputs", ERROR);
  
  Veclib::vmul (_size, a._data, 1, b._data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::timesPlus (const AuxField& a,
			       const AuxField& b)
// ---------------------------------------------------------------------------
// Add the product of a & b to this AuxField (in physical space).
// No dealiasing.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::timesPlus";

  if (_size != a._size || _size != b._size)
    message (routine, "non-congruent inputs", ERROR);

  Veclib::vvtvp (_size, a._data, 1, b._data, 1, _data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::timesMinus (const AuxField& a,
			        const AuxField& b)
// ---------------------------------------------------------------------------
// Subtract the product of a & b from this AuxField (in physical space).
// No dealiasing.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::timesMinus";

  if (_size != a._size || _size != b._size)
    message (routine, "non-congruent inputs", ERROR);

  Veclib::vvvtm (_size, _data, 1, a._data, 1, b._data, 1, _data, 1);

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
// ---------------------------------------------------------------------------
{
  const int_t nP = Geometry::planeSize();

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


AuxField& AuxField::axpy (const real_t    alpha,
			  const AuxField& x    )
// ---------------------------------------------------------------------------
// Add alpha * x to this AuxField.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::axpy";

  if (_size != x._size) message (routine, "non-congruent inputs", ERROR);

  Blas::axpy (_size, alpha, x._data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (_size, _data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::gradient (const int_t dir)
// ---------------------------------------------------------------------------
// Operate on AuxField to produce the nominated index of the gradient.
// dir == 0 ==> gradient in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
//
// NB: modfied for dual, only 1 non-zero mode, nominated as "K_FUND"
// in session.
// ---------------------------------------------------------------------------
{
  const char      routine[] = "AuxField::gradient";
  const int_t     nel  = Geometry::nElmt();
  const int_t     np   = Geometry::nP();
  const int_t     npnp = np  * np;
  const int_t     ntot = nel * npnp;
  const int_t     nP   = Geometry::planeSize();
  vector<real_t>  work;
  register real_t *xr, *xs, *tmp;
  register int_t  i, k;
  const real_t    *DV, *DT;

  Femlib::quadrature (0, 0, &DV, 0  , np, LL, 0.0, 0.0);
  Femlib::quadrature (0, 0, 0  , &DT, np, LL, 0.0, 0.0);

  switch (dir) {

  case 0:
    work.resize (2 * nP);
    xr = &work[0];
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      tmp = _plane[k];

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (tmp, tmp, xr, xs, DV, DT, np, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, tmp += npnp)
	_elmt[i] -> gradX (xr, xs, tmp);
      
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 1:
    work.resize (2 * nP);
    xr = &work[0];
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      tmp = _plane[k];

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (tmp, tmp, xr, xs, DV, DT, np, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, tmp += npnp)
	_elmt[i] -> gradY (xr, xs, tmp);
      
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 2: {
    const real_t betak = Femlib::value ("BETA * K_FUND");

    work.resize (nP);
    xr = &work[0];

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


void AuxField::gradient (const int_t nP ,
			 real_t*     src,
			 const int_t dir) const
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
// ---------------------------------------------------------------------------
{
  const char      routine[] = "AuxField::gradient";
  const int_t     nel  = Geometry::nElmt();
  const int_t     np   = Geometry::nP();
  const int_t     npnp = np  * np;
  const int_t     ntot = nel * npnp;
  register int_t  i, k;
  vector<real_t>  work;
  register real_t *plane, *xr, *xs, *Re, *Im;
  const real_t    *DV, *DT;

  Femlib::quadrature (0, 0, &DV, 0  , np, LL, 0.0, 0.0);
  Femlib::quadrature (0, 0, 0  , &DT, np, LL, 0.0, 0.0);

  switch (dir) {

  case 0:
    work.resize (2 * nP);
    xr = &work[0];
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      plane = src + k * nP;

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (plane, plane, xr, xs, DV, DT, np, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, plane += npnp)
	_elmt[i] -> gradX (xr, xs, plane);
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 1:
    work.resize (2 * nP);
    xr = &work[0];
    xs = xr + nP;

    for (k = 0; k < _nz; k++) {
      plane = src + k * nP;

      Veclib::zero  (2 * nP, xr, 1);
      Femlib::grad2 (plane, plane, xr, xs, DV, DT, np, np, nel);

      for (i = 0; i < nel; i++, xr += npnp, xs += npnp, plane += npnp)
	_elmt[i] -> gradY (xr, xs, plane);
      xr -= ntot;
      xs -= ntot;
    }
    break;

  case 2: {
    const real_t betak = Femlib::value ("BETA * K_FUND");

    work.resize (nP);
    xr = &work[0];

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


real_t AuxField::norm_inf () const
// ---------------------------------------------------------------------------
// Return infinity-norm (absolute max value) of AuxField data area.
// ---------------------------------------------------------------------------
{
  return fabs (_data[Blas::iamax (_size, _data, 1)]);
}


real_t AuxField::mode_L2 (const int_t mode) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
// Mode numbers run 0 -- n_z/2 - 1.
// ---------------------------------------------------------------------------
{
  const int_t       nel  = Geometry::nElmt();
  const int_t       npnp = Geometry::nTotElmt();
  register real_t   area = 0.0, Ek = 0.0, *Re, *Im;
  register int_t    i;
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


void AuxField::mode_en (const int_t mode  ,
			real_t&     RePart,
			real_t&     ImPart) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
//
// NB: Modified for dual: only 2 modes exist, 0 & 1.  We return
// mode_L2 split into components produced by real_t and imaginary parts.
// ---------------------------------------------------------------------------
{
  const int_t       nel  = Geometry::nElmt();
  const int_t       npnp = Geometry::nTotElmt();
  register real_t   area = 0.0, *Re, *Im;
  register int_t    i;
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


real_t AuxField::integral () const
// ---------------------------------------------------------------------------
// Return the total amount of scalar, integrated over spatial volume.
// It is assumed that the AuxField is in the Fourier-transformed state,
// so that the integration takes place over the zeroth Fourier mode
// only, then is scaled for Fourier normalisation.
// ---------------------------------------------------------------------------
{
  const int_t    nel  = Geometry::nElmt();
  const int_t    npnp = Geometry::nTotElmt();
  const real_t   Lz   = (Geometry::nDim()>2) ? Femlib::value("TWOPI/BETA"):1.;
  register int_t i;
  vector<real_t> work (npnp);
  real_t         total = 0.0, *p, *tmp = &work[0];

  ROOTONLY
    for (p = _plane[0], i = 0; i < nel; i++, p += npnp)
      total += _elmt[i] -> integral (p, tmp);

  return Lz * total;
}


real_t AuxField::integral (const int_t k) const
// ---------------------------------------------------------------------------
// Return the total amount of scalar, integrated over plane k.
// ---------------------------------------------------------------------------
{
  const int_t    nel  = Geometry::nElmt();
  const int_t    npnp = Geometry::nTotElmt();
  register int_t i;
  vector<real_t> work (npnp);
  real_t         total = 0.0, *p, *tmp = &work[0];

  for (p = _plane[k], i = 0; i < nel; i++, p += npnp)
    total += _elmt[i] -> integral (p, tmp);

  return total;
}


ofstream& operator << (ofstream& strm,
		       AuxField& F  )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "ofstream<<AuxField";
  const int_t    NP    = Geometry::planeSize();
  const int_t    nP    = Geometry::nPlane();
  register int_t i, k, n;

  for (i = 0; i < F._nz; i++) {
    strm.write((char*) F._plane[i], (int) (nP * sizeof (real_t))); 
    if (strm.bad())
      message (routine, "unable to write binary output", ERROR);
  }

  return strm;
}


ifstream& operator >> (ifstream&  strm,
		       AuxField&  F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.  Zero any unused storage areas.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "ifstream>>AuxField";
  const int_t    nP    = Geometry::nPlane();
  const int_t    NP    = Geometry::planeSize();
  const int_t    nProc = Geometry::nProc();
  register int_t i, k, n;

  for (i = 0; i < F._nz; i++) {
    strm.read ((char*) F._plane[i], (int) (nP * sizeof (real_t))); 
    if (strm.bad()) 
      message (routine, "unable to read binary input", ERROR);
    Veclib::zero (NP - nP, F._plane[i] + nP, 1);
  }

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


AuxField& AuxField::addToPlane (const int_t  k    ,
				const real_t alpha)
// ---------------------------------------------------------------------------
// Add in a constant to the values on nominated plane (if it exists),
// starting at plane zero.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::addToPlane";

  if (k < 0 || k >= _nz)
    message (routine, "nominated plane doesn't exist", ERROR);
  else
    Veclib::sadd (Geometry::nPlane(), alpha, _plane[k], 1, _plane[k], 1);

  return *this;
}


AuxField& AuxField::getPlane (const int_t k  ,
			      real_t*     tgt)
// ---------------------------------------------------------------------------
// Copy nominated plane to tgt.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::getPlane";

  if (k < 0 || k >= Geometry::nZProc())
    message (routine, "nominated plane doesn't exist", ERROR);
  else
    Veclib::copy (Geometry::nPlane(), _plane[k], 1, tgt, 1);

  return *this;
}


AuxField& AuxField::setPlane (const int_t   k  ,
			      const real_t* src)
// ---------------------------------------------------------------------------
// Copy copy src to nominated plane.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::setPlane";

  if (k < 0 || k >= _nz)
    message (routine, "nominated plane doesn't exist", ERROR);
  else
    Veclib::copy (Geometry::nPlane(), src, 1, _plane[k], 1);

  return *this;
}


AuxField& AuxField::setPlane (const int_t  k    ,
			      const real_t alpha)
// ---------------------------------------------------------------------------
// Set nominated plane to scalar alpha.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::setPlane";

  if (k < 0 || k >= _nz)
    message (routine, "nominated plane doesn't exist", ERROR);
  else {
    if (alpha == 0.0)
      Veclib::zero (Geometry::nPlane(), _plane[k], 1);
    else
      Veclib::fill (Geometry::nPlane(), alpha, _plane[k], 1);
  }

  return *this;
}


void AuxField::swapData (AuxField* x,
			 AuxField* y)
// ---------------------------------------------------------------------------
// Swap data areas of two fields.
// ---------------------------------------------------------------------------
{
  const char       routine[] = "AuxField::swapData";
  register int_t   k;
  register real_t* tmp;

  if (x -> _size != y -> _size)
    message (routine, "non-congruent inputs", ERROR);
 
  tmp        = x -> _data;
  x -> _data = y -> _data;
  y -> _data = tmp;

  for (k = 0; k < x -> _nz; k++) {
    tmp            = x -> _plane[k];
    x -> _plane[k] = y -> _plane[k];
    y -> _plane[k] = tmp;
  }
}


void AuxField::couple (AuxField*   v  ,
		       AuxField*   w  ,
		       const int_t dir)
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
// NB: modified for use with dual.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "Field::couple";
  const int_t    nP    =  Geometry::planeSize();
  vector<real_t> work (nP);
  real_t         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
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


AuxField& AuxField::divY ()
// ---------------------------------------------------------------------------
// Divide data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      npnp = Geometry::nTotElmt();
  register int_t   i, k;
  register real_t* p;

  for (k = 0; k < _nz; k++)
    for (p = _plane[k], i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> divY (p);
  
  return *this;
}


void AuxField::divY (const int_t nZ ,
		     real_t*     src) const
// ---------------------------------------------------------------------------
// Divide src by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      npnp = Geometry::nTotElmt();
  const int_t      ntot = Geometry::planeSize();
  register int_t   i, k;
  register real_t* p;   

  for (k = 0; k < nZ; k++)
    for (p = src + k * ntot, i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> divY (p);
}


AuxField& AuxField::mulY ()
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      npnp = Geometry::nTotElmt();
  register int_t   i, k;
  register real_t* p;

  for (k = 0; k < _nz; k++)
    for (p = _plane[k], i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> mulY (p);
  
  return *this;
}


void AuxField::mulY (const int_t nZ ,
		     real_t*     src) const
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      npnp = Geometry::nTotElmt();
  const int_t      ntot = Geometry::planeSize();
  register int_t   i, k;
  register real_t* p;   

  for (k = 0; k < nZ; k++)
    for (p = src + k * ntot, i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> mulY (p);
}


real_t AuxField::probe (const Element* E,
			const real_t   r,
			const real_t   s,
			const int_t    k) const
// ---------------------------------------------------------------------------
// Return the value of data on plane k, in Element E, location r, s.
// ---------------------------------------------------------------------------
{
  const int_t           offset = E -> ID() * Geometry::nTotElmt();
  static vector<real_t> work (3 * Geometry::nP());
  
  return E -> probe (r, s, _plane[k] + offset, &work[0]);
}


void AuxField::lengthScale (real_t* tgt) const
// ---------------------------------------------------------------------------
// Load tgt with data that represent the mesh resolution lengthscale
// at each planar location.
// ---------------------------------------------------------------------------
{
  const int_t      nel  = Geometry::nElmt();
  const int_t      npnp = Geometry::nTotElmt();
  register int_t   i;
  register real_t* p;

  for (p = tgt, i = 0; i < nel; i++, p += npnp)
    _elmt[i] -> lengthScale (p);
}


real_t AuxField::CFL (const int_t dir) const
// ---------------------------------------------------------------------------
// Return the inverse CFL timescale using this AuxField as a velocity 
// component in the nominated direction.  Computations only occur on the
// zeroth Fourier mode.
// dir == 0 ==> CFL estimate in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
// ---------------------------------------------------------------------------
{
  const char            routine[] = "AuxField::CFL";
  const int_t           nel  = Geometry::nElmt();
  const int_t           npnp = Geometry::nTotElmt();
  register int_t        i;
  register real_t*      p;
  real_t                dxy, CFL = 0.0;
  static vector<real_t> work (npnp);
  
  {
    const int_t   nP = Geometry::nP();
    const real_t* z;
    Femlib::quadrature (&z, 0, 0, 0, nP, LL, 0.0, 0.0);
    dxy = z[1] - z[0];
  }

  switch (dir) {
  case 0:
    for (p = _data, i = 0; i < nel; i++, p += npnp)
      CFL = max (CFL, _elmt[i] -> CFL (dxy, p, 0, &work[0]));
    break;
  case 1:
    for (p = _data, i = 0; i < nel; i++, p += npnp)
      CFL = max (CFL, _elmt[i] -> CFL (dxy, 0, p, &work[0]));
    break;
  case 2: {
    const int_t  nP = Geometry::nPlane();
    const real_t dz = Femlib::value ("TWOPI / BETA / N_Z");
    for (i = 0; i < nP; i++)
      CFL = max (CFL, fabs (_data[i]));
    CFL /= dz;
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }

  return CFL;
}


AuxField& AuxField::sqroot()
// ---------------------------------------------------------------------------
// Take sqrt of all data points.
// ---------------------------------------------------------------------------
{
  Veclib::vsqrt (_size, _data, 1, _data, 1);

  return *this;
}
