///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines for direct solution of Helmholtz problems.
//
// Copyright (c) 1994,2003 Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem_h>

static vector<MatrixSys*> MS;


ModalMatrixSys::ModalMatrixSys (const real              lambda2 ,
				const real              beta    ,
				const integer           baseMode,
				const integer           numModes,
				const vector<Element*>& Elmt    ,
				const BoundarySys*      Bsys    ,
				const SolverKind        method  )
// ---------------------------------------------------------------------------
// Generate or retrieve from internal database MS the vector of
// MatrixSys's which will be used to solve all the Fourier-mode
// discrete Helmholtz problems for the associated scalar Fields
// (called out by names).
//
// Input variables:
//   lambda2 : Helmholtz constant for the problem,	
//   beta    : Fourier length scale = TWOPI / Lz,
//   nmodes  : number of Fourier modes which will be solved,
//   Elmt    : vector of Element*'s used to make local Helmholtz matrices,
//   Bsys    : boundary system for this field.
//   method  : specify the kind of solver we want (Cholesky, PCG ...).
// ---------------------------------------------------------------------------
{
  const char name = Bsys -> field();
  integer    mode;
  bool       found;

  MatrixSys* M;
  vector<MatrixSys*>::iterator m;

  _fields = new char [strlen (Bsys -> Nsys (0) -> fields()) + 1];
  strcpy (_fields, Bsys -> Nsys (0) -> fields());
  _Msys.resize (numModes);

  if (method == DIRECT) {
    ROOTONLY cout << "-- Installing matrices for field '" << name << "' [";
    cout.flush();
    Femlib::synchronize();
  }

  for (mode = baseMode; mode < baseMode + numModes; mode++) {
    const integer    modeIndex = mode * Geometry::kFund();
    const NumberSys* N  = Bsys -> Nsys (modeIndex);
    const real    betak2   = sqr (Field::modeConstant (name, modeIndex, beta));
    const integer localMode = mode - baseMode;

    for (found = false, m = MS.begin(); !found && m != MS.end(); m++) {
      M     = *m;
      found = M -> match (lambda2, betak2, N, method);
    }
    if (found) {
      _Msys[localMode] = M;
      if (method == DIRECT) { cout << '.'; cout.flush(); }
    } else {
      _Msys[localMode] =
	new MatrixSys (lambda2, betak2, modeIndex, Elmt, Bsys, method);
      MS.insert (MS.end(), _Msys[localMode]);
      if (method == DIRECT) { cout << '*'; cout.flush(); }
    }
  }

  if (method == DIRECT) {
    Femlib::synchronize();
    ROOTONLY cout << "]" << endl;
    cout.flush();
  }
}


ModalMatrixSys::~ModalMatrixSys ()
// ---------------------------------------------------------------------------
// Destructor hands off calls to MatrixSys::~MatrixSys.  Note there
// can be side effects here, owing to the multiple use of MatrixSys*'s
// in different ModalMatrixSys's.  If multiple related ModalMatrixSys's
// exist, it is advisable to delete and recreate all of them before
// attempting reuse.
// ---------------------------------------------------------------------------
{
  integer N = _Msys.size();
  while (N--) delete (_Msys[N]);
  MS.resize (0);
}


MatrixSys::MatrixSys (const real              lambda2,
		      const real              betak2 ,
		      const integer           mode   ,
		      const vector<Element*>& elmt   ,
		      const BoundarySys*      bsys   ,
		      const SolverKind        method ) :
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
  const char       routine[] = "MatrixSys::MatrixSys";
  const integer    verbose   = Femlib::ivalue ("VERBOSE");
  const integer    np        = Geometry::nP();
  const integer    next      = Geometry::nExtElmt();
  const integer    nint      = Geometry::nIntElmt();
  const integer    npnp      = Geometry::nTotElmt();
  const integer*   bmap;
  register integer i, j, k, m, n;

  switch (_method) {

  case DIRECT: {
    vector<real>    work     (sqr (next) + sqr (np) + sqr (npnp));
    vector<integer> pivotmap (nint);
    real*           hbb  = &work[0];
    real*           rmat = hbb  + sqr (next);
    real*           rwrk = rmat + sqr (np);
    integer*        ipiv = &pivotmap[0];
    integer         info;

    _hbi    = new real*   [static_cast<size_t>(_nel)];
    _hii    = new real*   [static_cast<size_t>(_nel)];
    _bipack = new integer [static_cast<size_t>(_nel)];
    _iipack = new integer [static_cast<size_t>(_nel)];

    if (_nsolve) {
      _H = new real [static_cast<size_t>(_npack)];
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
	_hbi[j] = new real [static_cast<size_t>(_bipack[j])];
	_hii[j] = new real [static_cast<size_t>(_iipack[j])];
	Veclib::zero (_bipack[j], _hbi[j], 1);
	Veclib::zero (_iipack[j], _hii[j], 1);
      } else
	_hbi[j] = _hii[j] = 0;
      

      elmt[j]->HelmholtzSC(lambda2,betak2,hbb,_hbi[j],_hii[j],rmat,rwrk,ipiv);

      for (i = 0; i < next; i++)
	if ((m = bmap[i]) < _nsolve)
	  for (k = 0; k < next; k++)
	    if ((n = bmap[k]) < _nsolve && n >= m)
	      _H[Lapack::band_addr (m, n, _nband)] +=
		hbb[Veclib::row_major (i, k, next)];

      Family::adopt (_bipack[j], _hbi + j);
      Family::adopt (_iipack[j], _hii + j);

    }
    if (_nsolve) {
      // -- Loop over BCs and add diagonal contribution from mixed BCs.

      if (bsys -> mixBC()) {
	const integer  nbound = bsys -> nSurf();
	const integer* bmap   = _NS  -> btog();
	for (i = 0; i < nbound; i++)
	  _BC[i] -> augmentSC (_nband, _nsolve, bmap, rwrk, _H);
      }

      // -- Cholesky factor global banded-symmetric Helmholtz matrix.
    
      Lapack::pbtrf ("U",_nsolve,_nband-1,_H,_nband,info);

      Family::adopt (_npack, &_H);

      if (info) message (routine, "failed to factor Helmholtz matrix", ERROR);

      if (verbose) {
	real cond;
	pivotmap.resize (_nsolve);  ipiv = &pivotmap[0];
	work.resize (3 * _nsolve);  rwrk = &work[0];

	Lapack::pbcon ("U",_nsolve,_nband-1,_H,_nband,1.0,cond,rwrk,ipiv,info);
	cout << ", condition number: " << cond << endl;
      }
    }
  } break;

  case JACPCG: {
    const integer nbound = _BC.size();   
    real*         PCi;
    vector<real>  work (2 * npnp + np);
    real          *ed = &work[0], *ewrk = &work[0] + npnp;

    _PC = new real [static_cast<size_t>(_npts)];

    Veclib::zero (_npts, _PC, 1);

    PCi  = _PC + _NS -> nGlobal();
    bmap = _NS -> btog();

    // -- Mixed BC contributions.

    if (bsys -> mixBC())
      for (i = 0; i < nbound; i++)
	_BC[i] -> augmentDg (bmap, _PC);

    // -- Element contributions.

    for (i = 0; i < _nel; i++, bmap += next, PCi += nint) {
      elmt[i] -> HelmholtzDiag (lambda2, betak2, ed, ewrk);
      Veclib::scatr_sum (next, ed,  bmap,    _PC);
      Veclib::copy      (nint, ed + next, 1, PCi, 1);
    }

#if 1
    Veclib::vrecp (_npts, _PC, 1, _PC, 1);
#else  // -- Turn off preconditioner for testing.
    Veclib::fill  (_npts, 1.0, _PC, 1);
#endif

    Family::adopt (_npts, &_PC);

  } break;

  default:
    message (routine, "no solver of type requested -- never happen", ERROR);
    break;
  }
}


bool MatrixSys::match (const real       lambda2,
		       const real       betak2 ,
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
// storage we use the Femlib family routines.
// ---------------------------------------------------------------------------
{
  switch (_method) {
  case JACPCG:
    Family::abandon (&_PC);
    break;
  case DIRECT: {
    integer i;
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
