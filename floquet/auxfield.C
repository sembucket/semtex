//////////////////////////////////////////////////////////////////////////////
// auxfield.C: routines for AuxField class, including Fourier expansions.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn.
//
// Modified for stability analysis.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

AuxField::AuxField (real_t*           alloc,
		    const int_t       nz   ,
		    vector<Element*>& elmt ,
		    const char        name ) :
// ---------------------------------------------------------------------------
// Install field storage area and int size records.
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
  if   (val == 0.0) message ("AuxField::op /= real", "divide by zero", ERROR);
  else              Blas::scal (_size, 1.0 / val, _data, 1);

  return *this;
}


AuxField& AuxField::operator = (const AuxField& f)
// ---------------------------------------------------------------------------
// Modified for stability code. *This can have either 1 or 2 data
// planes and f has the same or less.
// ---------------------------------------------------------------------------
{
  const int_t nP = Geometry::planeSize();

  for (int_t i = 0; i < f._nz; i++)
    Veclib::copy (nP, f._plane[i], 1, _plane[i], 1);
  
  return *this;
}


AuxField& AuxField::operator += (const AuxField& f)
// ---------------------------------------------------------------------------
// Modified for stability code. *This can have either 1 or 2 data
// planes and f has the same or less.
// ---------------------------------------------------------------------------
{
  const int_t nP = Geometry::planeSize();

  for (int_t i = 0; i < f._nz; i++)
    Veclib::vadd (nP, _plane[i], 1, f._plane[i], 1, _plane[i], 1);

  return *this;
}


AuxField& AuxField::operator -= (const AuxField& f)
// ---------------------------------------------------------------------------
// Modified for stability code. *This can have either 1 or 2 data
// planes and f has the same or less.
// ---------------------------------------------------------------------------
{
  const int_t nP = Geometry::planeSize();

  for (int_t i = 0; i < f._nz; i++)
    Veclib::vsub (nP, _plane[i], 1, f._plane[i], 1, _plane[i], 1);

  return *this;
}


AuxField& AuxField::operator *= (const AuxField& f)
// ---------------------------------------------------------------------------
// Modified for stability code. *This can have either 1 or 2 data
// planes and f has the same or less.
// ---------------------------------------------------------------------------
{
  const int_t nP = Geometry::planeSize();

  for (int_t i = 0; i < f._nz; i++)
    Veclib::vmul (nP, _plane[i], 1, f._plane[i], 1, _plane[i], 1);

  return *this;
}


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.  Physical space.
// ---------------------------------------------------------------------------
{
  const int_t    nel = Geometry::nElmt();
  const int_t    np2 = Geometry::nTotElmt();
  const int_t    kb  = Geometry::basePlane();
  const real_t   dz  = Femlib::value ("TWOPI / BETA / N_Z");
  register int_t i, k;
  real_t*        p;

  for (k = 0; k < _nz; k++) {
    Femlib::value ("z", (kb + k) * dz);
    for (p = _plane[k], i = 0; i < nel; i++, p += np2)
      _elmt[i] -> evaluate (function, p);
  }
  
  return *this;
}


AuxField& AuxField::times (const AuxField& a,
			   const AuxField& b)
// ---------------------------------------------------------------------------
// Set this AuxField equal to the product of a & b (in physical space).
//
// This version is specialised for stability calculations: a can have
// 1 or 2 data planes, but b has one. *this conforms with a.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::times";
  
  if (b._nz > 1) message (routine, "second operand must be real", ERROR);

  Veclib::vmul (b._size, a._plane[0], 1, b._plane[0], 1, _plane[0], 1);
  if (a._nz == 2)
  Veclib::vmul (b._size, a._plane[1], 1, b._plane[0], 1, _plane[1], 1);

  return *this;
}


AuxField& AuxField::timesPlus (const AuxField& a,
			       const AuxField& b)
// ---------------------------------------------------------------------------
// Add the product of a & b to this AuxField.
//
// This version is specialised for stability calculations: a can have
// 1 or 2 data planes, but b has one. *this conforms with a.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::timesPlus";

  if (b._nz > 1) message (routine, "second operand must be real", ERROR);

  Veclib::vvtvp
    (b._size, a._plane[0], 1, b._plane[0], 1, _plane[0], 1, _plane[0], 1);
  if (a._nz == 2) Veclib::vvtvp
    (b._size, a._plane[1], 1, b._plane[0], 1, _plane[1], 1, _plane[1], 1);

  return *this;
}


AuxField& AuxField::timesMinus (const AuxField& a,
			        const AuxField& b)
// ---------------------------------------------------------------------------
// Subtract the product of a & b from this AuxField.
//
// This version is specialised for stability calculations: a can have
// 1 or 2 data planes, but b has one. *this conforms with a.
// ---------------------------------------------------------------------------
{
  const char routine[] = "AuxField::timesMinus";

  if (b._nz > 1) message (routine, "second operand must be real", ERROR);

  Veclib::vvvtm
    (b._size, a._plane[0], 1, b._plane[0], 1, _plane[0], 1, _plane[0], 1);
  if (a._nz == 2) Veclib::vvvtm
    (b._size, a._plane[1], 1, b._plane[0], 1, _plane[1], 1, _plane[1], 1);

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
// We check the how many planes of data exist, when taking the
// gradients in z direction; if only one, then the AuxField is
// half-complex (and the sign of z-gradient may need changing).
// ---------------------------------------------------------------------------
{
  const char     routine[] = "AuxField::gradient";
  const int_t    nel  = Geometry::nElmt();
  const int_t    np   = Geometry::nP();
  const int_t    npnp = np  * np;
  const int_t    ntot = nel * npnp;
  const int_t    nP   = Geometry::planeSize();
  int_t          i, k;
  vector<real_t> work;
  real_t         *tmp;

  switch (dir) {

  case 0:
    work.resize (2 * npnp);
    for (k = 0; k < _nz; k++) {
     tmp = _plane[k];
     for (i = 0; i < nel; i++, tmp += npnp)
       _elmt[i] -> grad (tmp, 0, &work[0]);
    }
    break;

  case 1:
    work.resize (2 * npnp);
    for (k = 0; k < _nz; k++) {
      tmp = _plane[k];
      for (i = 0; i < nel; i++, tmp += npnp)
	_elmt[i] -> grad (0, tmp, &work[0]);
    }
    break;

  case 2: {
    const real_t beta = Femlib::value ("BETA");

    if (_nz == 1)          // -- Half-complex (may need sign change).
      Blas::scal (nP, beta, _plane[0], 1);
    else {			// -- Full-complex.
      work.resize  (nP);
      Veclib::copy (nP,        _plane[0], 1,  &work[0], 1);
      Veclib::smul (nP, -beta, _plane[1], 1, _plane[0], 1);
      Veclib::smul (nP,  beta,  &work[0], 1, _plane[1], 1);
    }
  } break;
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }

  return *this;
}


real_t AuxField::mode_L2 (const int_t mode) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
// Mode numbers run 0 -- n_z/2 - 1.
// ---------------------------------------------------------------------------
{
  const char  routine[] = "AuxField::mode_L2";
  const int_t nel  = Geometry::nElmt();
  const int_t kr   = 2 * mode;
  const int_t ki   = kr + 1;
  const int_t npnp = Geometry::nTotElmt();
  real_t      area = 0.0, Ek = 0.0;
  int_t       i;
  Element*    E;
  
  if (kr < 0  ) message (routine, "negative mode number",        ERROR);
  if (ki > _nz) message (routine, "mode number exceeds maximum", ERROR);

  for (i = 0; i < nel; i++) {
    E      = _elmt[i];
    area  += E -> area();
    Ek    += sqr (E -> norm_L2 (_plane[kr] + i*npnp));
    if (_nz > 1) 
      Ek  += sqr (E -> norm_L2 (_plane[ki] + i*npnp));
  }

  return Ek / (2.0 * area);
}


real_t AuxField::integral () const
// ---------------------------------------------------------------------------
// Return the total amount of scalar, integrated over spatial volume.
// It is assumed that the AuxField is in the Fourier-transformed state,
// so that the integration takes place over the zeroth Fourier mode
// only, then is scaled for Fourier normalisation.
// ---------------------------------------------------------------------------
{
  const int_t  nel  = Geometry::nElmt();
  const int_t  npnp = Geometry::nTotElmt();
  const real_t Lz   = (Geometry::nPert()>2) ? Femlib::value("TWOPI/BETA"):1.0;
  int_t        i;
  vector<real_t> work (npnp);
  real_t         total = 0.0, *p;

  ROOTONLY
    for (p = _plane[0], i = 0; i < nel; i++, p += npnp)
      total += _elmt[i] -> integral (p, &work[0]);

  return Lz * total;
}


ostream& operator << (ostream& strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
//
// For multiple-processor jobs, only the root processor does output,
// receiving data from other processors.  This ensures that the data
// are written out in the correct order, and that only one processor
// needs access to the output stream.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "ostream<<AuxField";
  const int_t    NP    = Geometry::planeSize();
  const int_t    nP    = Geometry::nPlane();
  const int_t    nProc = Geometry::nProc();
  register int_t i, k;

  if (nProc > 1) {

    ROOTONLY {
      vector<real_t> buffer (NP);

      for (i = 0; i < F._nz; i++)
        strm.write((char*) F._plane[i], static_cast<int>(nP * sizeof(real_t)));
        if (strm.bad())
	  message (routine, "unable to write binary output", ERROR);

      for (k = 1; k < nProc; k++)
	for (i = 0; i < F._nz; i++) {
	  Femlib::recv (&buffer[0], NP, k);
	  strm.write((char*)&buffer[0], static_cast<int>(nP * sizeof(real_t)));
          if (strm.bad()) 
	    message (routine, "unable to write binary output", ERROR);
	}

    } else for (i = 0; i < F._nz; i++) Femlib::send (F._plane[i], NP, 0);

  } else {

    for (i = 0; i < F._nz; i++) {
      strm.write((char*) F._plane[i], static_cast<int> (nP * sizeof(real_t)));
      if (strm.bad())
	message (routine, "unable to write binary output", ERROR);
    }
  }

  return strm;
}


istream& operator >> (istream& strm,
		      AuxField& F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.  Zero any unused storage areas.
//
// As for the write operator, only the root processor accesses strm.
// This precaution is possibly unnecessary for input.
// ---------------------------------------------------------------------------
{
  const char     routine[] = "istream>>AuxField";
  const int_t    nP    = Geometry::nPlane();
  const int_t    NP    = Geometry::planeSize();
  const int_t    nProc = Geometry::nProc();
  register int_t i, k;

  if (nProc > 1) {

    ROOTONLY {
      vector<real_t> buffer (NP);

      for (i = 0; i < F._nz; i++) {
	strm.read ((char*) F._plane[i], static_cast<int>(nP * sizeof(real_t)));
        if (strm.bad()) 
	  message (routine, "unable to read binary input", ERROR);
	Veclib::zero (NP - nP, F._plane[i] + nP, 1);
      }

      for (k = 1; k < nProc; k++) {
	for (i = 0; i < F._nz; i++) {
	  strm.read ((char*)&buffer[0], static_cast<int>(nP * sizeof(real_t)));
          if (strm.bad()) 
	    message (routine, "unable to read binary input", ERROR);
	  Veclib::zero (NP - nP, &buffer[0] + nP, 1);
	  Femlib::send (&buffer[0], NP, k);
	}
      }
    } else for (i = 0; i < F._nz; i++) Femlib::recv (F._plane[i], NP, 0);

  } else {

    for (i = 0; i < F._nz; i++) {
      strm.read ((char*) F._plane[i], static_cast<int>(nP * sizeof(real_t))); 
      if (strm.bad()) 
	message (routine, "unable to read binary input", ERROR);
      Veclib::zero (NP - nP, F._plane[i] + nP, 1);
    }
  }

  return strm;
}


void AuxField::describe (char* s)  const
// ---------------------------------------------------------------------------
// Load s with a (prism-compatible) description of field geometry:
// NR NS NZ NEL.
// ---------------------------------------------------------------------------
{
#if 1
  ostringstream sf;
  sf << Geometry::nP()    << " "
     << Geometry::nP()    << " "
     << Geometry::nZ()    << " "
     << Geometry::nElmt() << ends;
  strcpy (s, sf.str().c_str());
#else
  ostrstream (s, StrMax) << Geometry::nP()    << " "
                         << Geometry::nP()    << " "
                         << Geometry::nZ()    << " "
                         << Geometry::nElmt() << ends;
#endif
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
// dir == FORWARD
// --------------
//           v~ <-- v + i w
//           w~ <-- v - i w
// dir == INVERSE
// --------------
//           v  <-- 0.5   * (v~ + w~)
//           w  <-- 0.5 i * (w~ - v~)
//
// Since there is no coupling for the viscous terms in the 2D equation,
// do nothing for the zeroth Fourier mode.
//
// For stability calc, we have to take account of possible
// half-complex case, where v~, w~ are both pure real_t.
// ---------------------------------------------------------------------------
{
  if (Geometry::problem() == Geometry::O2_2D ||
      Geometry::problem() == Geometry::SO2_2D ) return;

  const char     routine[] = "Field::couple";
  const int_t    nP        =  Geometry::planeSize();
  vector<real_t> work (nP);
  real_t         *Vr, *Vi, *Wr, *Wi, *tp = &work[0];
  
  if (dir == FORWARD) {
    if (v->_nz == 1) {		// -- Half-complex.
      Vr = v -> _plane[0];
      Wi = w -> _plane[0];

      Veclib::copy (nP, Vr, 1, tp, 1);
      Veclib::vsub (nP, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);
    } else {			// -- Full complex.
      Vr = v -> _plane[0];
      Vi = v -> _plane[1];
      Wr = w -> _plane[0];
      Wi = w -> _plane[1];

      Veclib::copy (nP, Vr, 1, tp, 1);
      Veclib::vsub (nP, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nP, Wr, 1, tp, 1);
      Veclib::copy (nP, Wi, 1, Wr, 1);
      Veclib::vsub (nP, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nP, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == INVERSE) {

    if (v->_nz == 1) {		// -- Half-complex.
      Vr = v -> _plane[0];
      Wr = w -> _plane[0];

      Veclib::copy  (nP,      Vr, 1, tp, 1);
      Veclib::svvpt (nP, 0.5, Vr, 1, Wr, 1, Vr, 1);
      Veclib::svvmt (nP, 0.5, Wr, 1, tp, 1, Wr, 1);
    } else {			// -- Full complex.
      Vr = v -> _plane[0];
      Vi = v -> _plane[1];
      Wr = w -> _plane[0];
      Wi = w -> _plane[1];

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


real_t AuxField::probe (const Element* E,
			const real_t   r,
			const real_t   s,
			const int_t    k) const
// ---------------------------------------------------------------------------
// Return the value of data on plane k, in Element E, location r, s.
// ---------------------------------------------------------------------------
{
  const int_t    offset = E -> ID() * Geometry::nTotElmt();
  vector<real_t> work (3 * Geometry::nP());
  
  return E -> probe (r, s, _plane[k] + offset, &work[0]);
}


real_t AuxField::probe (const Element* E,
			const real_t   r,
			const real_t   s,
			const real_t   z) const
// ---------------------------------------------------------------------------
// Return the value of data, in Element E, location r, s, z.
//
// NB: interpolation assumes that AuxField is Fourier transformed.
//
// For multiprocessor operation, the Fourier interpolation is done on
// the root processor, and the return value is only valid on that
// processor.  The approach taken here is inefficient for
// multiprocessor work, and it would be more rational to redesign to
// make the message buffers as big as possible, i.e. to collect all
// the data for each history point on each processor before passing it
// to the root processor for interpolation.
// ---------------------------------------------------------------------------
{
  const int_t      nZ     = Geometry::nZ();
  const int_t      nP     = Geometry::nProc();
  const int_t      np     = Geometry::nP();
  const int_t      NZH    = nZ >> 1;
  const int_t      NHM    = NZH - 1;
  const int_t      offset = E -> ID() * Geometry::nTotElmt();
  const real_t     betaZ  = z * Femlib::value ("BETA");

  register int_t   k, Re, Im;
  register real_t  value, phase;
  vector<real_t>   work (nZ + _nz + 3 * np);
  register real_t* fbuf = &work[0];
  register real_t* lbuf = fbuf + nZ;
  real_t*          ewrk = lbuf + _nz;

  if (nP > 1) {
    for (k = 0; k < _nz; k++)
      lbuf[k] = E -> probe (r, s, _plane[k] + offset, ewrk);
    
    ROOTONLY {
      Veclib::copy (_nz, lbuf, 1, fbuf, 1);
      for (k = 1; k < nP; k++) 
	Femlib::recv (fbuf + k * _nz, _nz, k);
    } else
      Femlib::send (lbuf, _nz, 0);

  } else {
    if (_nz < 3)			// -- Hey!  This is 2D!
      return value = E -> probe (r, s, _plane[0] + offset, ewrk);
  
    else {
      for (k = 0; k < _nz; k++)
	fbuf[k] = E -> probe (r, s, _plane[k] + offset, ewrk);
    }
  }

  ROOTONLY {
    Blas::scal (nZ - 2, 2.0, fbuf + 2, 1);
    
    value  = fbuf[0];
    value += fbuf[1] * cos (NZH * betaZ);
    for (k = 1; k < NHM; k++) {
      Re     = k  + k;
      Im     = Re + 1;
      phase  = k * betaZ;
      value += fbuf[Re] * cos (phase) - fbuf[Im] * sin (phase);
    }
  } else
    value = 0.0;
   
  return value;
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
  const char       routine[] = "AuxField::CFL";
  const int_t      nel  = Geometry::nElmt();
  const int_t      npnp = Geometry::nTotElmt();
  register int_t   i;
  register real_t* p;
  real_t           dxy, CFL = 0.0;
  vector<real_t>   work (npnp);
  
  {
    const int_t   nP = Geometry::nP();
    const real_t* z;
    Femlib::quadrature (&z, 0, 0, 0, nP, 'L', 0.0, 0.0);
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
    const int_t nP = Geometry::nPlane();
    const real_t    dz = Femlib::value ("TWOPI / BETA / N_Z");
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


AuxField& AuxField::update (const int_t   nSlice ,
			    const real_t* FFTdata,
			    const real_t  time   ,
			    const real_t  period )
// ---------------------------------------------------------------------------
// Fourier series reconstruction of base flow, if it is periodic in time.
//
// Fourier transformed fields have been pre-scaled.
// ---------------------------------------------------------------------------
{
  if (nSlice < 2) return *this;

  const int_t  nPlane = Geometry::planeSize();
  const real_t BetaT  = TWOPI * fmod (time, period) / period;
  real_t       phase;
  int_t        i;
  
  // -- For each point in plane do Fourier interpolation in time.

  Veclib::copy (nPlane, FFTdata, 1, _data, 1);
  Blas::axpy   (nPlane, cos(0.5*nSlice*BetaT), FFTdata + nPlane, 1, _data, 1);

  for (i = 2; i < nSlice; i += 2) {
    phase = (i>>1) * BetaT;
    Blas::axpy (nPlane,  cos(phase), FFTdata +  i   *nPlane, 1, _data, 1);
    Blas::axpy (nPlane, -sin(phase), FFTdata + (i+1)*nPlane, 1, _data, 1);
  }

  return *this;
}
