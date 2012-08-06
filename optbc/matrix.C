///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines for direct solution of Helmholtz problems.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

static vector<MatrixSys*> MS;



MatrixSys::MatrixSys (const real_t            lambda2,
		      const real_t            betak2 ,
		      const int_t             mode   ,
		      const vector<Element*>& elmt   ,
		      const BoundarySys*      bsys   ,
		      const SolverKind        method ,
		      const vector<AuxField*>&  Baseflow,
		      bool  forwards) :
// ---------------------------------------------------------------------------
// Initialize and factorise matrices in this system.
//
// For method == DIRECT:
//   Matrices are assembled using LAPACK-compatible ordering systems.
//   Global Helmholtz matrix uses symmetric-banded format; elemental
//   Helmholtz matrices (hii & hbi) use column-major formats.
// For method == JACPCG:
//   Build and invert diagonal preconditioner matrix.
// ---------------------------------------------------------------------------
// NB: these get evaluated in the order they appear in the class
// definition!:
  _HelmholtzConstant (lambda2),
  _FourierConstant   (betak2 ),
  _BC                (bsys -> BCs  (mode)),
  _NS                (bsys -> Nsys (mode)),
  _nel               (Geometry::nElmt()),
  _nglobal           (_NS -> nGlobal()),
  _singular          ((_HelmholtzConstant + _FourierConstant) < EPSSP &&
		     !_NS -> fmask() && !bsys -> mixBC()),
  _nsolve            ((_singular) ? _NS -> nSolve() - 1 : _NS -> nSolve()),
  _method            (method),
  _nband             (_NS -> nBand()),
  _npack             (_nband * _nsolve),
  _H                 (0),
  _hbi               (0),
  _hii               (0),
  _bipack            (0),
  _iipack            (0),
  _npts              (_nglobal + Geometry::nInode()),
  _PC                (0)
{
  const char     routine[] = "MatrixSys::MatrixSys";
  const int_t    verbose   = Femlib::ivalue ("VERBOSE");
  const int_t    np        = Geometry::nP();
  const int_t    next      = Geometry::nExtElmt();// number of external nodes per element
  const int_t    nint      = Geometry::nIntElmt();
  const int_t    npnp      = Geometry::nTotElmt();
  const int_t*   bmap;
  register int_t i, j, k, m, n;
  switch (_method) {

  case DIRECT: {
    vector<real_t> work     (sqr (next) + sqr (np) + sqr (npnp));
    vector<int_t>  pivotmap (nint);
    real_t*        hbb  = &work[0];
    real_t*        rmat = hbb  + sqr (next);
    real_t*        rwrk = rmat + sqr (np);
    int_t*         ipiv = &pivotmap[0];
    int_t          info;

    _hbi    = new real_t*[static_cast<size_t>(_nel)];
    _hii    = new real_t*[static_cast<size_t>(_nel)];
    _bipack = new int_t  [static_cast<size_t>(_nel)];
    _iipack = new int_t  [static_cast<size_t>(_nel)];

    if (_nsolve) {
      _H = new real_t [static_cast<size_t>(_npack)];
      Veclib::zero (_npack, _H, 1);

      if (verbose > 1)
	cout << endl
	     << "Helmholtz constant (lambda2): " << setw(10) << lambda2
	     << ", Fourier constant (betak2): "  << setw(10) << betak2;
      if (verbose)
	cout << endl << "System matrix: " << _nsolve << "x" << _nband
	     << "\t(" << _npack << " words)";
    }

    // -- Loop over elements, creating & posting elemental Helmholtz matrices.

    for (bmap = _NS -> btog(), j = 0; j < _nel; j++, bmap += next) {
      _bipack[j] = next * nint;
      _iipack[j] = nint * nint;
      
      if (nint) {
	_hbi[j] = new real_t [static_cast<size_t>(_bipack[j])];
	_hii[j] = new real_t [static_cast<size_t>(_iipack[j])];
	Veclib::zero (_bipack[j], _hbi[j], 1);
	Veclib::zero (_iipack[j], _hii[j], 1);
      } else
	_hbi[j] = _hii[j] = 0;
      

      elmt[j]->HelmholtzSC(lambda2,betak2,hbb,_hbi[j],_hii[j],rmat,rwrk,ipiv);

      for (i = 0; i < next; i++)
	if ((m = bmap[i]) < _nsolve)
	  for (k = 0; k < next; k++)
	    if ((n = bmap[k]) < _nsolve && n >= m)
	      _H[Lapack::band_addr (m, n, _nband)] += // band_addr: return (n+1)*_nband-n+m-1;
		hbb[Veclib::row_major (i, k, next)];  // row_major: return k+i*next;
    
      Family::adopt (_bipack[j], _hbi + j);
      Family::adopt (_iipack[j], _hii + j);

    }
    if (_nsolve) {
      // -- Loop over BCs and add diagonal contribution from mixed BCs.

      if (bsys -> mixBC()) {
	const int_t  nbound = bsys -> nSurf();
	const int_t* bmap   = _NS  -> btog();
	for (i = 0; i < nbound; i++)
	  _BC[i] -> augmentSC (_nband, _nsolve, bmap, rwrk, _H);
      }

      // -- Loop over BCs and add diagonal contribution from Toutflow BCs.
      if (bsys -> ToutflowBC()) {
	// if (Geometry::nSlice()>1) message (routine, "Sorry, The DIRECT solver is only for LNS based on steady base flow. This solver does not support time-dependent base flow or nonlinear works", ERROR); 
	const int_t  nbound = bsys -> nSurf();
	const int_t* bmap   = _NS  -> btog();
	for (i = 0; i < nbound; i++) {
	  _BC[i] -> switchK (Baseflow[0]->plane(0), Baseflow[1]->plane(0), forwards);  
	  _BC[i] -> augmentSC (_nband, _nsolve, bmap, rwrk, _H);
	}
      }
	  

		
		
      // -- Cholesky factor global banded-symmetric Helmholtz matrix.
    
      Lapack::pbtrf ("U",_nsolve,_nband-1,_H,_nband,info);

      Family::adopt (_npack, &_H);

      if (info) message (routine, "failed to factor Helmholtz matrix", ERROR);

      if (verbose) {
	real_t cond;
	pivotmap.resize (_nsolve);  ipiv = &pivotmap[0];
	work.resize (3 * _nsolve);  rwrk = &work[0];

	Lapack::pbcon ("U",_nsolve,_nband-1,_H,_nband,1.0,cond,rwrk,ipiv,info);
	cout << ", condition number: " << cond << endl;
      }
    }
  } break;
    
  case JACPCG: {

    const int_t    nbound = _BC.size();   
    real_t*        PCi;
    vector<real_t> work (2 * npnp + np);
    real_t         *ed = &work[0], *ewrk = &work[0] + npnp;
    
    _PC = new real_t [static_cast<size_t>(_npts)];
    _PC_notoutflow = new real_t [static_cast<size_t>(_NS -> nGlobal())];

    Veclib::zero (_npts, _PC, 1);
    
    PCi  = _PC + _NS -> nGlobal();
    bmap = _NS -> btog();
    
    // -- Mixed BC contributions.

    if (bsys -> mixBC())
      for (i = 0; i < nbound; i++)
	_BC[i] -> augmentDg (bmap, _PC);


    // -- Toutflow BC contributions are moved out and calculated separately.

  
    // -- Element contributions.

    for (i = 0; i < _nel; i++, bmap += next, PCi += nint) {
      elmt[i] -> HelmholtzDiag (lambda2, betak2, ed, ewrk);
      Veclib::scatr_sum (next, ed,  bmap,    _PC);
      Veclib::copy      (nint, ed + next, 1, PCi, 1);
    }
    
    Veclib::copy      (_NS -> nGlobal(), _PC, 1, _PC_notoutflow, 1);
    Veclib::vrecp (_npts, _PC, 1, _PC, 1);
    
  } break;
    
  default:
    message (routine, "no solver of type requested -- never happen", ERROR);
    break;
  }
	
}


bool MatrixSys::match (const real_t     lambda2,
		       const real_t     betak2 ,
		       const NumberSys* nScheme,
		       const SolverKind method ) const
// ---------------------------------------------------------------------------
// The unique identifiers of a MatrixSys are presumed to be given
// by the constants and the numbering system used.  Other things that
// could be checked but aren't (yet) include geometric systems and
// quadrature schemes.
// ---------------------------------------------------------------------------
{
  if (fabs (_HelmholtzConstant - lambda2) < EPSDP                       &&
      fabs (_FourierConstant   - betak2 ) < EPSDP                       &&
      _NS -> nGlobal() == nScheme -> nGlobal()                          &&
      _NS -> nSolve()  == nScheme -> nSolve()                           &&
      Veclib::same (_NS->nGlobal(), _NS->btog(), 1, nScheme->btog(), 1) &&
      _method == method                                                  )
    return true;

  else
    return false;
}


MatrixSys::~MatrixSys()
// ---------------------------------------------------------------------------
// Destructor.  Because there may be aliases to the internal vector
// storage we use the family class routines.
// ---------------------------------------------------------------------------
{
  switch (_method) {
  case JACPCG:
    Family::abandon (&_PC);
    break;
  case DIRECT: {
    int_t i;
    for (i = 0; i < _nel; i++) {
      Family::abandon (_hbi + i);
      Family::abandon (_hii + i);
    }
    Family::abandon (&_H);
    delete[] _bipack;
    delete[] _iipack;
  } break;
  default:
    break;
  }
}



void MatrixSys::Contribution_toutflow (  vector<AuxField*>&  Baseflow) 
// ---------------------------------------------------------------------------
// contribution of controlbc and toutflow to _PC
// ---------------------------------------------------------------------------

{  real_t*    PC_toutflow    = new  real_t [static_cast<size_t>(_NS -> nGlobal())];
  Veclib::zero (_NS -> nGlobal(), PC_toutflow, 1);
  
  register int_t  i;
  const    int_t  nbound = _BC.size();
  const    int_t* bmap   = _NS  -> btog();
  const    int_t  np     = Geometry::nP();	

	
  for (i = 0; i < nbound; i++) {
    if (*_BC[i]->group() =='t') 
      _BC[i] -> switchK (Baseflow[0]->plane(0), Baseflow[1]->plane(0), 0);//the last parameter is forwards==0, since this function is called only for the backward integration.  
    _BC[i] -> augmentDg (bmap, PC_toutflow);	 
  }
  
  
  Veclib::vadd (_NS -> nGlobal(), PC_toutflow, 1, _PC_notoutflow, 1, PC_toutflow, 1);
  
  Veclib::vrecp (_NS -> nGlobal(), PC_toutflow, 1, PC_toutflow, 1);
  
  Veclib::copy      (_NS -> nGlobal(), PC_toutflow, 1, _PC, 1);
  
  delete [] PC_toutflow;
  
  
}
