///////////////////////////////////////////////////////////////////////////////
// data2df.C: simple 2D x Fourier data class.
//
// Copyright (c) 2004 <--> $Date$, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>
#include <data2df.h>


Data2DF::Data2DF (const int_t nP  ,
		  const int_t nZ  ,
		  const int_t nEl ,
		  const char  Name) :
  _name (Name),
  _np   (nP  ),
  _nz   (nZ  ),
  _nel  (nEl ),
  _np2  (nP * nP)
// ---------------------------------------------------------------------------
// Data2DF constructor. 
// ---------------------------------------------------------------------------
{
  int_t i;
  
  _nplane = _np * _np * _nel;
  if (_nplane & 1) _nplane++;
  _ntot   = _nplane * _nz;

  _data  = new real_t  [_ntot];
  _plane = new real_t* [_nz];

  for (i = 0; i < _nz; i++) _plane[i] = _data + i*_nplane;
  Veclib::zero (_ntot, _data, 1);
}


Data2DF& Data2DF::DFT1D (const int_t sign)
// ---------------------------------------------------------------------------
// Carry out discrete Fourier transformation in z direction.
// ---------------------------------------------------------------------------
{
  if (_nz > 2) Femlib::DFTr (_data, _nz, _nplane, sign);

  return *this;
}


Data2DF& Data2DF::DPT2D (const int_t sign, 
			 const char  basis)
// ---------------------------------------------------------------------------
// Carry out 2D discrete polynomial transform (element-by-element) on planes.
//
// Use basis = 'l' to transform to Legendre polynomial basis, 'm' for
// modal polynomial basis.
// ---------------------------------------------------------------------------
{
  int_t          i;
  vector<real_t> work (_nplane);
  const real_t   *Fu, *Ft, *Bu, *Bt;

  if (basis == 'l')
    Femlib::legTran (_np, &Fu, &Ft, &Bu, &Bt, 0, 0);
  else
    Femlib::modTran (_np, &Fu, &Ft, &Bu, &Bt, 0, 0);

  if (sign == FORWARD)
    for (i = 0; i < _nz; i++)
      Femlib::tpr2d (_plane[i], _plane[i], &work[0], Fu, Ft, _np, _np, _nel);
  else
    for (i = 0; i < _nz; i++)
      Femlib::tpr2d (_plane[i], _plane[i], &work[0], Bu, Bt, _np, _np, _nel);

  return *this;
}


Data2DF& Data2DF::conjugate (const bool zero)
// ---------------------------------------------------------------------------
// Take complex conjugate. If zero is true, assume that mode zero is
// complex instead of two real_t modes packed together.
// ---------------------------------------------------------------------------
{
  int_t i, first = (zero) ? 1 : 3;

  if (_nz > 1)
    for (i = first; i < _nz; i += 2)
      Veclib::neg (_nplane, _plane[i], 1);

  return *this;
} 


Data2DF& Data2DF::symmetrize (const bool zero)
// ---------------------------------------------------------------------------
// Symmetrize velocity field component. This enforces a reflection
// symmetry of the velocity and pressure fields, as follows.
// 
// 'u': mode_k.Im = 0, k > 0
// 'v': mode_k.Im = 0, k > 0
// 'w': mode_k.Re = 0, k > 0
// 'p': mode_k.Im = 0, k > 0
//
// If zero is true, assume that mode zero is complex instead of two
// real modes packed together.
// ---------------------------------------------------------------------------
{
  int_t i, first;

  switch (_name) {
  case 'u': case 'v': case 'p': first = (zero) ? 1 : 3; break;
  case 'w':                     first = (zero) ? 0 : 2; break;
  default:                      return *this;           break;
  }

  if (_nz > 1)
    for (i = first; i < _nz; i += 2)
      Veclib::zero (_nplane, _plane[i], 1);

  return *this;
} 


Data2DF& Data2DF::operator = (const Data2DF& rhs)
// ---------------------------------------------------------------------------
// If the two fields conform, copy rhs's data storage to lhs.
//
// Otherwise perform projection/interpolation of rhs's data area to lhs.
// Interpolation ASSUMES THAT FOURIER TRANSFORMATION HAS ALREADY OCCURRED
// in z direction if rhs is 3D.  Truncation of Fourier modes occurs if this
// Data2DF has less modes than rhs (to avoid aliasing).
// ---------------------------------------------------------------------------
{
  if (rhs._nel != _nel)
    message ("Data2DF::operator =", "fields can't conform", ERROR);

  if (rhs._np == _np && rhs._nz == _nz)
    Veclib::copy (_ntot, rhs._data, 1, _data, 1);

  else {			// -- Perform projection.

    register int_t  i, k;
    register real_t *LHS, *RHS;
    const real_t    *IN,  *IT;
    const int_t     nzm = min (rhs._nz, _nz);
    vector<real_t>  work (rhs._np * _np);
    real_t*         tmp = &work[0];

    Femlib::projection (&IN, &IT, rhs._np, 'L', 0.0, 0.0, _np, 'L', 0.0, 0.0);

    for (k = 0; k < nzm; k++) {	// -- 2D planar projections.
      LHS = _plane[k];
      RHS = rhs._plane[k];

      if (rhs._np == _np)
	Veclib::copy (_nplane, RHS, 1, LHS, 1);
      else
	for (i = 0; i < _nel; i++, LHS += _np2, RHS += rhs._np2) {
	  Blas::mxm (IN,  _np, RHS, rhs._np, tmp, rhs._np);
	  Blas::mxm (tmp, _np, IT,  rhs._np, LHS,     _np);
	}
    }

    if ((i = _nz - rhs._nz) > 0) // -- Zero pad for Fourier projections.
      Veclib::zero (i * _nplane, _data + rhs._ntot, 1);
  }

  return *this;
}


Data2DF& Data2DF::operator = (const real_t val)
// ---------------------------------------------------------------------------
// Set field storage area to val.
// ---------------------------------------------------------------------------
{
  if   (val == 0.0) Veclib::zero (_ntot,      _data, 1);
  else              Veclib::fill (_ntot, val, _data, 1);

  return *this;
}


Data2DF& Data2DF::reverse ()
// ---------------------------------------------------------------------------
// Reverse order of bytes within each word of data.
// ---------------------------------------------------------------------------
{
  Veclib::brev (_ntot, _data, 1, _data, 1);

  return *this;
}


ostream& operator << (ostream&  strm,
		      Data2DF& F    )
// ---------------------------------------------------------------------------
// Binary write of F's data area.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  for (i = 0; i < F._nz; i++)
    strm.write (reinterpret_cast<char*> (F._plane[i]),
		F._np * F._np * F._nel * sizeof (real_t));

  return strm;
}


istream& operator >> (istream&  strm,
		      Data2DF& F    )
// ---------------------------------------------------------------------------
// Binary read of F's data area.
// ---------------------------------------------------------------------------
{
  int_t i;
  
  for (i = 0; i < F._nz; i++)
    strm.read (reinterpret_cast<char*> (F._plane[i]),
	       F._np * F._np * F._nel * sizeof (real_t));

  return strm;
}


Data2DF& Data2DF::filter1D (const real_t roll ,
			    const int_t  order)
// ---------------------------------------------------------------------------
// Carry out filtering in the Fourier domain, with a Boyd--VanDeven
// filter (see Leven, Iskandarani and Haidvogel JCP 137). The input
// parameter roll [0, 1] specifies where the filter begins to roll
// over, while parameter order, 2 or larger, specifies the shape of the
// filter.
//
// It is assumed that the data are already Fourier transformed on input.
// ---------------------------------------------------------------------------
{
  int_t          i;
  const int_t    nh = _nz >> 1, ord = max (2, order);
  vector<real_t> filter (nh + 1), mask (_nz);
  const real_t   lag = clamp (roll, 0.0, 1.0);

  Femlib::erfcFilter (nh, ord, lag * nh, 1.0, &filter[0]);
  mask[0] = filter[ 0];
  mask[1] = filter[nh];
  for (i = 1; i < nh; i++) mask[2*i] = mask[2*i+1] = filter[i];

  for (i = 0; i < _nz; i++)
    Veclib::smul  (_nplane, mask[i], _plane[i], 1, _plane[i], 1);

  return *this;
}


Data2DF& Data2DF::filter2D (const real_t roll ,
			    const int_t  order)
// ---------------------------------------------------------------------------
// Carry out 2D filtering element-by-element in the modal polynomial
// domain, with a Boyd--VanDeven filter (see Leven, Iskandarani and
// Haidvogel JCP 137). The input parameter roll [0, 1] specifies where
// the filter begins to roll over, while parameter order, 2 or larger,
// specifies the shape of the filter.
//
// Doing the filtering in the modal polynomial space ensures that the field
// retains C0 continuity across element boundaries, see Blackburn & Schmidt
// JCP 186.
// ---------------------------------------------------------------------------
{
  int_t          i, j;
  const real_t   lag = clamp (roll, 0.0, 1.0);
  const int_t    ord = max (2, order);
  vector<real_t> filter (_np);
  vector<real_t> work (_nplane);
  vector<real_t> Iu (_np2), It (_np2);
  const real_t   *Du, *Dt, *dpt;
  real_t         *dataplane;

  Femlib::erfcFilter (_np - 1, ord, lag, 1.0, &filter[0]);

  // -- Polynomial transform+filter (Du, Dt) and inversion (Iu, It) matrices.

  Femlib::modTran (_np, &Du, &Dt, 0, &dpt, 0, 0);

  for (i = 0; i < _np; i++)
    Veclib::smul (_np, filter[i], dpt + i*_np, 1, &It[i*_np], 1);

  Blas::mxm    (Dt, _np, &It[0], _np, &Iu[0], _np);
  Veclib::copy (_np2,    &Iu[0], 1,   &It[0], 1);

  for (i = 0; i < _np; i++)
    for (j = 0; j < _np; j++)
      Iu[Veclib::row_major (j, i, _np)] = It[Veclib::row_major (i, j, _np)];

  for (i = 0; i < _nz; i++)
    Femlib::tpr2d (_plane[i],_plane[i],&work[0],&Iu[0],&It[0],_np,_np,_nel);

  return *this;
}
