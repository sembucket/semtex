//////////////////////////////////////////////////////////////////////////////
// auxfield.C: routines for AuxField class, including Fourier expansions.
//
// For 2D problems, the data storage is organized by 2D Elements.
//
// For 3D problems, Fourier expansions are used in the 3rd direction,
// and each Fourier mode can be thought of as a separate 2D problem
// (with real and imaginary parts, or planes, of 2D data).  The data
// are then organized plane-by-plane, with each plane being a 2D
// AuxField; if in physical space there are nz planes of data, then
// there are nz/2 Fourier modes.  Data for the Nyquist mode are stored
// as the imaginary part of the zeroth Fourier mode, but are kept zero
// and never evolve.  The planes always point to the same storage
// locations within the data area.
//
// The data are transformed to physical space for storage in restart
// files.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>


AuxField::AuxField (real*             alloc,
		    const integer     nz   ,
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
  const char       routine[] = "AuxField::AuxField";
  const integer    nP = Geometry::planeSize();
  register integer k;

  if (Geometry::nElmt() != _elmt.getSize())
    message (routine, "conflicting number of elements in input data", ERROR);

  _plane = new real* [(size_t) _nz];

  for (k = 0; k < _nz; k++) _plane[k] = _data + k * nP;
}


AuxField& AuxField::setInput (real*         alloc,
			      const integer nz   )
// ---------------------------------------------------------------------------
// Install a new lot of field storage, with nominated number of planes.
// It is assumed that alloc is at least nZ * geometry::planeSize() long.
// ---------------------------------------------------------------------------
{
  register integer k;
  const integer    nP = Geometry::planeSize();

  _nz   = nz;
  _size = _nz * nP;
  _data = alloc;

  delete [] _plane;
  _plane = new real* [(size_t) _nz];

  for (k = 0; k < _nz; k++) _plane[k] = _data + k * nP;

  return *this;
}


AuxField& AuxField::operator = (const real val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (_size,      _data, 1);
  else              Veclib::fill (_size, val, _data, 1);

  return *this;
}


AuxField& AuxField::operator += (const real val)
// ---------------------------------------------------------------------------
// Add val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (_size, val, _data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator -= (const real val)
// ---------------------------------------------------------------------------
// Add -val to field storage area.
// ---------------------------------------------------------------------------
{
  if (val != 0.0) Veclib::sadd (_size, -val, _data, 1, _data, 1);

  return *this;
}


AuxField& AuxField::operator *= (const real val)
// ---------------------------------------------------------------------------
// Multiply field storage area by val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (_size,      _data, 1);
  else              Blas  ::scal (_size, val, _data, 1);

  return *this;
}


AuxField& AuxField::operator /= (const real val)
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


AuxField& AuxField::operator = (const char* function)
// ---------------------------------------------------------------------------
// Set AuxField's value to temporo-spatially varying function.  Physical space.
// ---------------------------------------------------------------------------
{
  const integer     nel = Geometry::nElmt();
  const integer     np2 = Geometry::nTotElmt();
  const integer     kb  = Geometry::basePlane();
  const real        dz  = Femlib::value ("TWOPI / BETA / N_Z");
  register integer  i, k;
  real              *p;

  for (k = 0; k < _nz; k++) {
    Femlib::value ("z", (kb + k) * dz);
    for (p = _plane[k], i = 0; i < nel; i++, p += np2)
      _elmt[i] -> evaluate (function, p);
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
  const char    routine[] = "AuxField::innerProduct";
  const integer ndim = Geometry::nDim();
  integer       i;

  if (_size != a[0]->_size || _size != b[0]->_size)
    message (routine, "non-congruent inputs", ERROR);
  
  Veclib::zero (_size, _data, 1);

  for (i = 0; i < ndim; i++)
    Veclib::vvtvp (_size, a[i]->_data, 1, b[i]->_data, 1, _data, 1, _data, 1);

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


AuxField& AuxField::axpy (const real      alpha,
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


AuxField& AuxField::gradient (const integer dir)
// ---------------------------------------------------------------------------
// Operate on AuxField to produce the nominated index of the gradient.
// dir == 0 ==> gradient in first direction, 1 ==> 2nd, 2 ==> 3rd.
// AuxField is presumed to have been Fourier transformed in 3rd direction.
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
    const integer nmodes = Geometry::nModeProc();
    const integer base   = Geometry::baseMode();
    const real    beta   = Femlib::value ("BETA");
    integer       Re, Im, klo;

    work.setSize (nP);
    xr = work();

    if (base == 0) { // -- We have real & Nyquist planes, to be set zero.
      klo = 1; Veclib::zero (2 * nP, _data, 1);
    } else
      klo = 0;

    for (k = klo; k < nmodes; k++) {
      Re = k  + k;
      Im = Re + 1;
      Veclib::copy (nP,                     _plane[Re], 1, xr,         1);
      Veclib::smul (nP, -beta * (k + base), _plane[Im], 1, _plane[Re], 1);
      Veclib::smul (nP,  beta * (k + base), xr,         1, _plane[Im], 1);
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
			 const integer nP ,
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
// NB: the Fourier mode index is assumed to start at zero for all processes.
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

    for (k = 0; k < nZ; k++) {
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

    for (k = 0; k < nZ; k++) {
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
    if (nZ == 1) break;

    const integer nmodes = nZ >> 1;
    const real    beta   = Femlib::value ("BETA");

    work.setSize (nP);
    xr = work();

    Veclib::zero (2 * nP, src, 1);

    for (k = 1; k < nmodes; k++) {
      Re = src + 2 * k * nP;
      Im = Re  + nP;
      Veclib::copy (nP,             Re, 1, xr, 1);
      Veclib::smul (nP, -beta * k,  Im, 1, Re, 1);
      Veclib::smul (nP,  beta * k,  xr, 1, Im, 1);
    }
    break;
  }
  default:
    message (routine, "nominated direction out of range [0--2]", ERROR);
    break;
  }
}


void AuxField::errors (const Mesh* mesh    ,
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
    return;
  }
  
  const integer nq   = 15;
  const integer nqnq = nq * nq;
  const integer np   = Geometry::nP();
  const integer npnp = Geometry::nTotElmt();
  const integer nel  = Geometry::nElmt();

  Element       *E, *P;
  real          area = 0.0, Li = 0.0, L2 = 0.0, H1 = 0.0;
  vector<real>  work (np * nq + 2 * nqnq);
  real          *u, *err = work(), *sol = err + nqnq, *tmp = sol + nqnq;
  const real    *z, **IN, **IT;
  integer       k;

  Femlib::mesh (GLL, GLL, nq, nq, &z, 0,   0,  0, 0);
  Femlib::mesh (GLL, GLL, np, nq, 0, &IN, &IT, 0, 0);

  for (u = _plane[0], k = 0; k < nel; k++, u += npnp) {

    E = _elmt[k];
    P = new Element (E -> ID(), mesh, z, nq);

    Blas::mxm (*IN, nq,  u,  np, tmp, np);
    Blas::mxm (tmp, nq, *IT, np, err, nq);

    P -> evaluate (function, sol);
    Veclib::vsub  (nqnq, err, 1, sol, 1, err, 1);

    area += P -> area ();
    Li    = max (Li, P -> norm_inf (err));
    L2   += P -> norm_L2 (err);
    H1   += P -> norm_H1 (err);

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
}


real AuxField::norm_inf () const
// ---------------------------------------------------------------------------
// Return infinity-norm (absolute max value) of AuxField data area.
// ---------------------------------------------------------------------------
{
  return fabs (_data[Blas::iamax (_size, _data, 1)]);
}


real AuxField::mode_L2 (const integer mode) const
// ---------------------------------------------------------------------------
// Return energy norm per unit area for indicated mode = 1/(2*A) \int u.u dA.
// Mode numbers run 0 -- n_z/2 - 1.
// ---------------------------------------------------------------------------
{
  const char        routine[] = "AuxField::mode_L2";
  const integer     nel  = Geometry::nElmt();
  const integer     kr   = 2 * mode;
  const integer     ki   = kr + 1;
  const integer     npnp = Geometry::nTotElmt();
  register real     area = 0.0, Ek = 0.0, *Re = _plane[kr], *Im = _plane[ki];
  register integer  i;
  register Element* E;
  
  if (kr < 0  ) message (routine, "negative mode number",        ERROR);
  if (ki > _nz) message (routine, "mode number exceeds maximum", ERROR);

  for (i = 0; i < nel; i++, Re += npnp, Im += npnp) {
    E      = _elmt[i];
    area  += E -> area();
    Ek    += sqr (E -> norm_L2 (Re));
    if (_nz > 1) 
      Ek  += sqr (E -> norm_L2 (Im));
  }

  return Ek / (2.0 * area);
}


real AuxField::integral () const
// ---------------------------------------------------------------------------
// Return the total amount of scalar, integrated over spatial volume.
// It is assumed that the AuxField is in the Fourier-transformed state,
// so that the integration takes place over the zeroth Fourier mode
// only, then is scaled for Fourier normalisation.
// ---------------------------------------------------------------------------
{
  const integer nel  = Geometry::nElmt();
  const integer npnp = Geometry::nTotElmt();
  const real    Lz   = (Geometry::nDim() > 2) ? Femlib::value("TWOPI/BETA"):1.;
  register integer i;
  vector<real>     work (npnp);
  real             total = 0.0, *p, *tmp = work();

  ROOTONLY
    for (p = _plane[0], i = 0; i < nel; i++, p += npnp)
      total += _elmt[i] -> integral (p, tmp);

  return Lz * total;
}


real AuxField::integral (const integer k) const
// ---------------------------------------------------------------------------
// Return the total amount of scalar, integrated over plane k.
// ---------------------------------------------------------------------------
{
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  register integer i;
  vector<real>     work (npnp);
  real             total = 0.0, *p, *tmp = work();

  for (p = _plane[k], i = 0; i < nel; i++, p += npnp)
    total += _elmt[i] -> integral (p, tmp);

  return total;
}


ofstream& operator << (ofstream& strm,
		       AuxField& F  )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
//
// For multiple-processor jobs, only the root processor does output,
// receiving data from other processors.  This ensures that the data
// are written out in the correct order, and that only one processor
// needs access to the output stream.
//
// UNIX interface used for IO in order to avoid buffering done by C++.
// For this to work strm has to be of ofstream rather than ostream class.
// ---------------------------------------------------------------------------
{
  const char       routine[] = "ofstream<<AuxField";
  const integer    NP    = Geometry::planeSize();
  const integer    nP    = Geometry::nPlane();
  const integer    nProc = Geometry::nProc();
  const int        fd    = strm.rdbuf() -> fd();
  register integer i, k, n;

  if (nProc > 1) {

    ROOTONLY {
      vector<real> buffer (NP);

      strm.rdbuf() -> sync();

      for (i = 0; i < F._nz; i++)
	n = write (fd, F._plane[i], nP * sizeof (real));
	if (n != nP * sizeof (real))
	  message (routine, "unable to write binary input", ERROR);

      for (k = 1; k < nProc; k++)
	for (i = 0; i < F._nz; i++) {
	  Femlib::recv (buffer(), NP, k);
	  n = write (fd, buffer(), nP * sizeof (real));
	  if (n != nP * sizeof (real))
	    message (routine, "unable to write binary input", ERROR);
	}

      strm.rdbuf() -> sync();      
      strm.flush();

    } else for (i = 0; i < F._nz; i++) Femlib::send (F._plane[i], NP, 0);

  } else {

    strm.rdbuf() -> sync();

    for (i = 0; i < F._nz; i++) {
      n = write (fd, F._plane[i], nP * sizeof (real));
      if (n != nP * sizeof (real))
	message (routine, "unable to write binary input", ERROR);
    }
  }

  strm.rdbuf() -> sync();
  return strm;
}


ifstream& operator >> (ifstream&  strm,
		       AuxField&  F   )
// ---------------------------------------------------------------------------
// Binary read of F's data area.  Zero any unused storage areas.
//
// As for the write operator, only the root processor accesses strm.
// This precaution is possibly unnecessary for input.
//
// UNIX interface used for IO in order to avoid buffering done by C++.
// ---------------------------------------------------------------------------
{
  const char       routine[] = "ifstream>>AuxField";
  const integer    nP    = Geometry::nPlane();
  const integer    NP    = Geometry::planeSize();
  const integer    nProc = Geometry::nProc();
  const int        fd    = strm.rdbuf() -> fd();
  register integer i, k, n;

  if (nProc > 1) {

    ROOTONLY {
      vector<real> buffer (NP);

      strm.rdbuf() -> sync();     

      for (i = 0; i < F._nz; i++) {
	n = read (fd, F._plane[i], nP * sizeof (real));
	if (n != nP * sizeof (real))
	  message (routine, "unable to read binary input", ERROR);
	Veclib::zero (NP - nP, F._plane[i] + nP, 1);
      }

      for (k = 1; k < nProc; k++) {
	for (i = 0; i < F._nz; i++) {
	  n = read (fd, buffer(), nP*sizeof (real));
	  if (n != nP * sizeof (real))
	    message (routine, "unable to read binary input", ERROR);
	  Veclib::zero (NP - nP, buffer() + nP, 1);
	  Femlib::send (buffer(), NP, k);
	}
      }
    } else for (i = 0; i < F._nz; i++) Femlib::recv (F._plane[i], NP, 0);

  } else {

    strm.rdbuf() -> sync();

    for (i = 0; i < F._nz; i++) {
      n = read (fd, F._plane[i], nP * sizeof (real));
      if (n != nP * sizeof (real))
	message (routine, "unable to read binary input", ERROR);
      Veclib::zero (NP - nP, F._plane[i] + nP, 1);
    }
  }

  strm.rdbuf() -> sync();    

  return strm;
}


AuxField& AuxField::zeroNyquist ()
// ---------------------------------------------------------------------------
// Set storage for highest frequency mode to zero.  This mode is
// carried but never evolves, and is stored as the second data plane
// on the lowest-numbered processor.
// ---------------------------------------------------------------------------
{
  ROOTONLY if (_nz > 1) Veclib::zero (Geometry::planeSize(), _plane[1], 1);

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
// Discrete Fourier transform in homogeneous direction.  Number of
// points in that direction must be even, but is otherwise
// unrestricted.  Use sign = 1 for forward transform, -1 for inverse.
//
// Normalization is carried out on forward transform, so that the zeroth
// mode's real data are the average over the homogeneous direction of the
// physical space values.
//
// For multiple-processor execution, data must be gathered across processors
// prior to Fourier transform, then scattered back.
// ---------------------------------------------------------------------------
{
  const integer nzt = Geometry::nZ();
  const integer nP  = Geometry::planeSize();
  const integer nPR = Geometry::nProc();
  const integer nPP = Geometry::nBlock();

  if (nPR == 1) {
    if (nzt > 1)
      if (nzt == 2)
	if   (sign == +1) Veclib::zero (nP, _plane[1], 1);
	else              Veclib::copy (nP, _plane[0], 1, _plane[1], 1);
      else
	Femlib::DFTr  (_data, nzt, nP, sign);

  } else {
    Femlib::exchange (_data, _nz,  nP,   +1);
    Femlib::DFTr     (_data, nzt, nPP, sign);
    Femlib::exchange (_data, _nz,  nP,   -1);

  }

  return *this;
}


AuxField& AuxField::transform32 (const integer sign,
				 real*         phys)
// ---------------------------------------------------------------------------
// Discrete Fourier transform in homogeneous direction, extended for
// dealiasing.  Input pointer phys points to data in physical space,
// which acts as input area if sign == +1, output area if sign == -1.
// So transform is from phys to internal storage if sign == +1 and
// vice versa.  After transform of either type, the data have normal
// planar configuration.
//
// NB: dealiasing does not occur in multiple-processor execution, so phys
// has the same number of data as *this.
// ---------------------------------------------------------------------------
{
  const integer nZ   = Geometry::nZ();
  const integer nP   = Geometry::planeSize();
  const integer nZ32 = Geometry::nZ32();

  if (Geometry::nProc() == 1) {	 // -- Single processor.

    const integer nTot32 = nZ32 * nP;
    const integer nPad   = nTot32 - _size;

    if (nZ <= 2) {
      if   (sign == +1) Veclib::copy (_size,  phys, 1, _data, 1);
      else              Veclib::copy (_size, _data, 1,  phys, 1);
    } else {
      if (sign == +1) {
	Femlib::DFTr (phys, nZ32, nP, +1);
	Veclib::copy (_size, phys, 1, _data, 1);
      } else {
	Veclib::copy (_size, _data, 1, phys, 1);
	Veclib::zero (nPad, phys + _size, 1);
	Femlib::DFTr (phys, nZ32, nP, -1);
      }
    }

  } else {			// -- Multiple processor.
    
    const integer nPP = Geometry::nBlock();

    if (sign == +1) {
      Femlib::exchange (phys, nZ32, nP,  +1);
      Femlib::DFTr     (phys, nZ,   nPP, +1);
      Femlib::exchange (phys, nZ32, nP,  -1);
      Veclib::copy     (_size, phys, 1, _data, 1);
    } else {
      Veclib::copy     (_size, _data, 1, phys, 1);
      Femlib::exchange (phys, nZ32, nP,  +1);
      Femlib::DFTr     (phys, nZ,   nPP, -1);
      Femlib::exchange (phys, nZ32, nP,  -1);
    }
  }

  return *this;
}


AuxField& AuxField::DLT2D (const integer sign,
			   const real*   filt)
// ---------------------------------------------------------------------------
// Carry out 2D discrete Legendre transform on an element-by-element
// basis, for each plane of data.  Sign = 1 ==> forward transform.
// Optionally apply a set of (1D) filter weights during forward or inverse
// transformation.
// ---------------------------------------------------------------------------
{
  const integer    np  = Geometry::nP();
  const integer    nel = Geometry::nElmt();
  const integer    np2 = sqr(np);

  register integer p, q, pq, r, s, rs;
  register real    cr, cs, P, Q;
  integer          i, k;
  vector<real>     work (np2);
  real             *src, *tmp = work();
  const real       *w, *legtab;

  Femlib::legCoef (np, &legtab);
  Femlib::quad    (LL, np, np, 0, 0, &w, 0, 0, 0, 0);

  if (filt) {			// -- Apply filter within transform.

    if (sign == 1) {		// -- Forward transform.
      for (k = 0; k < _nz; k++) {
	for (src = _plane[k], i = 0; i < nel; i++, src += np2) {
	  Veclib::zero (np2, tmp, 1);
	  for (rs = 0, r = 0; r < np; r++) {
	    cr = filt[r] * legtab[Veclib::row_major (np, r, np)];
	    for (s = 0; s < np; s++, rs++) {
	      cs = filt[s] * legtab[Veclib::row_major (np, s, np)];
	      for (pq = 0, p = 0; p < np; p++) {
		P = legtab[Veclib::row_major (r, p, np)];
		for (q = 0; q < np; q++, pq++) {
		  Q = legtab[Veclib::row_major (s, q, np)];
		  tmp[rs] += w[p] * w[q] * P * Q * src[pq];
		}
	      }
	      tmp[rs] *= cr * cs;
	    }
	  }
	  Veclib::copy (np2, tmp, 1, src, 1);
	}
      }
    } else {			// -- Inverse transform.
      for (k = 0; k < _nz; k++) {
	for (src = _plane[k], i = 0; i < nel; i++, src += np2) {
	  Veclib::zero (np2, tmp, 1);
	  for (rs = 0, r = 0; r < np; r++) {
	    for (s = 0; s < np; s++, rs++) {
	      for (pq = 0, p = 0; p < np; p++) {
		P = filt[p] * legtab[Veclib::row_major (p, r, np)];
		for (q = 0; q < np; q++, pq++) {
		  Q = filt[q] * legtab[Veclib::row_major (q, s, np)];
		  tmp[rs] += P * Q * src[pq];
		}
	      }
	    }
	  }
	  Veclib::copy (np2, tmp, 1, src, 1);
	}
      }
    }

  } else {			// -- Plain transform. no filter.

    if (sign == 1) {		// -- Forward transform.
      for (k = 0; k < _nz; k++) {
	for (src = _plane[k], i = 0; i < nel; i++, src += np2) {
	  Veclib::zero (np2, tmp, 1);
	  for (rs = 0, r = 0; r < np; r++) {
	    cr = legtab[Veclib::row_major (np, r, np)];
	    for (s = 0; s < np; s++, rs++) {
	      cs = legtab[Veclib::row_major (np, s, np)];
	      for (pq = 0, p = 0; p < np; p++) {
		P = legtab[Veclib::row_major (r, p, np)];
		for (q = 0; q < np; q++, pq++) {
		  Q = legtab[Veclib::row_major (s, q, np)];
		  tmp[rs] += w[p] * w[q] * P * Q * src[pq];
		}
	      }
	      tmp[rs] *= cr * cs;
	    }
	  }
	  Veclib::copy (np2, tmp, 1, src, 1);
	}
      }
    } else {			// -- Inverse transform.
      for (k = 0; k < _nz; k++) {
	for (src = _plane[k], i = 0; i < nel; i++, src += np2) {
	  Veclib::zero (np2, tmp, 1);
	  for (rs = 0, r = 0; r < np; r++) {
	    for (s = 0; s < np; s++, rs++) {
	      for (pq = 0, p = 0; p < np; p++) {
		P = legtab[Veclib::row_major (p, r, np)];
		for (q = 0; q < np; q++, pq++) {
		  Q = legtab[Veclib::row_major (q, s, np)];
		  tmp[rs] += P * Q * src[pq];
		}
	      }
	    }
	  }
	  Veclib::copy (np2, tmp, 1, src, 1);
	}
      }
    }
  }

  return *this;
}


AuxField& AuxField::DFfilt (const real* filter)
// ---------------------------------------------------------------------------
// Presuming data are already in Fourier-transformed state, apply a set
// of filter weights to the data on the present processor.  Filter
// coefficients are real only.
// ---------------------------------------------------------------------------
{
  const integer    nM  = Geometry::nModeProc();
  const integer    nP  = Geometry::planeSize();
  const integer    pID = Geometry::procID();
  register integer i, k, Re, Im, offset;

  ROOTONLY {
    Blas::scal (nP, filter[0], _plane[0], 1);
    for (k = 1; k < nM; k++) {
      Re = k  + k;
      Im = Re + 1;
      Blas::scal (nP, filter[k], _plane[Re], 1);
      Blas::scal (nP, filter[k], _plane[Im], 1);
    }
  } else {
    offset = pID * nM;
    for (k = 0; k < nM; k++) {
      i  = k  + offset;
      Re = k  + k;
      Im = Re + 1;
      Blas::scal (nP, filter[i], _plane[Re], 1);
      Blas::scal (nP, filter[i], _plane[Im], 1);
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

  if (k < 0 || k >= _nz)
    message (routine, "nominated plane doesn't exist", ERROR);
  else
    Veclib::sadd (Geometry::nPlane(), alpha, _plane[k], 1, _plane[k], 1);

  return *this;
}


AuxField& AuxField::getPlane (const integer k  ,
			      real*         tgt)
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


AuxField& AuxField::setPlane (const integer k  ,
			      const real*   src)
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


AuxField& AuxField::setPlane (const integer k    ,
			      const real    alpha)
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
  register integer k;
  register real*   tmp;

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
// ---------------------------------------------------------------------------
{
  if (Geometry::nDim() < 3) return;

  const char       routine[] = "Field::couple";
  const integer    nP    =  Geometry::planeSize();
  const integer    nMode =  Geometry::nModeProc();
  const integer    kLo   = (Geometry::procID() == 0) ? 1 : 0;
  register integer k, Re, Im;
  vector<real>     work (nP);
  real             *Vr, *Vi, *Wr, *Wi, *tp = work();
  
  if (dir == 1) {

    for (k = kLo; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> _plane[Re];
      Vi = v -> _plane[Im];
      Wr = w -> _plane[Re];
      Wi = w -> _plane[Im];

      Veclib::copy (nP, Vr, 1, tp, 1);
      Veclib::vsub (nP, Vr, 1, Wi, 1, Vr, 1);
      Veclib::vadd (nP, Wi, 1, tp, 1, Wi, 1);
      Veclib::copy (nP, Wr, 1, tp, 1);
      Veclib::copy (nP, Wi, 1, Wr, 1);
      Veclib::vsub (nP, Vi, 1, tp, 1, Wi, 1);
      Veclib::vadd (nP, Vi, 1, tp, 1, Vi, 1);
    }

  } else if (dir == -1) {

    for (k = kLo; k < nMode; k++) {
      Re = k  + k;
      Im = Re + 1;

      Vr = v -> _plane[Re];
      Vi = v -> _plane[Im];
      Wr = w -> _plane[Re];
      Wi = w -> _plane[Im];

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
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  register integer i, k;
  register real*   p;

  for (k = 0; k < _nz; k++)
    for (p = _plane[k], i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> divR (p);
  
  return *this;
}


void AuxField::divR (const integer nZ ,
		     real*         src) const
// ---------------------------------------------------------------------------
// Divide src by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  const integer    ntot = Geometry::planeSize();
  register integer i, k;
  register real*   p;   

  for (k = 0; k < nZ; k++)
    for (p = src + k * ntot, i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> divR (p);
}


AuxField& AuxField::mulR ()
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  register integer i, k;
  register real*   p;

  for (k = 0; k < _nz; k++)
    for (p = _plane[k], i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> mulR (p);
  
  return *this;
}


void AuxField::mulR (const integer nZ ,
		     real*         src) const
// ---------------------------------------------------------------------------
// Multiply data values by radius (i.e. y in cylindrical coords).
// ---------------------------------------------------------------------------
{
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  const integer    ntot = Geometry::planeSize();
  register integer i, k;
  register real*   p;   

  for (k = 0; k < nZ; k++)
    for (p = src + k * ntot, i = 0; i < nel; i++, p += npnp)
      _elmt[i] -> mulR (p);
}


real AuxField::probe (const Element* E,
		      const real     r,
		      const real     s,
		      const integer  k) const
// ---------------------------------------------------------------------------
// Return the value of data on plane k, in Element E, location r, s.
// ---------------------------------------------------------------------------
{
  const integer offset = E -> ID() * Geometry::nTotElmt();
  
  return E -> probe (r, s, _plane[k] + offset);
}


real AuxField::probe (const Element* E,
		      const real     r,
		      const real     s,
		      const real     z) const
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
  const integer    nZ     = Geometry::nZ();
  const integer    nP     = Geometry::nProc();
  const integer    NZH    = nZ >> 1;
  const integer    NHM    = NZH - 1;
  const integer    offset = E -> ID() * Geometry::nTotElmt();
  const real       betaZ  = z * Femlib::value ("BETA");

  register integer k, Re, Im;
  register real    value, phase;
  vector<real>     work (nZ + _nz);
  register real*   fbuf = work();
  register real*   lbuf = fbuf + nZ;

  if (nP > 1) {
    for (k = 0; k < _nz; k++)
      lbuf[k] = E -> probe (r, s, _plane[k] + offset);
    
    ROOTONLY {
      Veclib::copy (_nz, lbuf, 1, fbuf, 1);
      for (k = 1; k < nP; k++) 
	Femlib::recv (fbuf + k * _nz, _nz, k);
    } else
      Femlib::send (lbuf, _nz, 0);

  } else {
    if (_nz < 3)			// -- Hey!  This is 2D!
      return value = E -> probe (r, s, _plane[0] + offset);
  
    else {
      for (k = 0; k < _nz; k++)
	fbuf[k] = E -> probe (r, s, _plane[k] + offset);
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


void AuxField::lengthScale (real* tgt) const
// ---------------------------------------------------------------------------
// Load tgt with data that represent the mesh resolution lengthscale
// at each planar location.
// ---------------------------------------------------------------------------
{
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  register integer i;
  register real*   p;

  for (p = tgt, i = 0; i < nel; i++, p += npnp)
    _elmt[i] -> lengthScale (p);
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
  const char       routine[] = "AuxField::CFL";
  const integer    nel  = Geometry::nElmt();
  const integer    npnp = Geometry::nTotElmt();
  register integer i;
  register real*   p;
  real             dxy, CFL = 0.0;
  
  {
    const integer nP = Geometry::nP();
    const real*   z;
    Femlib::quad (LL, nP, nP, &z, 0, 0, 0, 0, 0, 0);
    dxy = z[1] - z[0];
  }

  switch (dir) {
  case 0:
    for (p = _data, i = 0; i < nel; i++, p += npnp)
      CFL = max (CFL, _elmt[i] -> CFL (dxy, p, 0));
    break;
  case 1:
    for (p = _data, i = 0; i < nel; i++, p += npnp)
      CFL = max (CFL, _elmt[i] -> CFL (dxy, 0, p));
    break;
  case 2: {
    const integer nP = Geometry::nPlane();
    const real    dz = Femlib::value ("TWOPI / BETA / N_Z");
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
