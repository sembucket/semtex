///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines for direct solution of Helmholtz problems.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

static List<MatrixSys*> MS;


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
  const char               name = Bsys -> field();
  integer                  found, mode;
  ListIterator<MatrixSys*> m (MS);
  MatrixSys*               M;

  _fields = new char [strlen (Bsys -> Nsys (0) -> fields()) + 1];
  strcpy (_fields, Bsys -> Nsys (0) -> fields());
  _Msys.setSize (numModes);

  if (method == DIRECT) {
    ROOTONLY cout << "-- Installing matrices for field '" << name << "' [";
    cout.flush();
    Femlib::synchronize();
  }

  for (mode = baseMode; mode < baseMode + numModes; mode++) {
    const NumberSys* N         = Bsys -> Nsys (mode);
    const real       betak2    = sqr (Field::modeConstant (name, mode, beta));
    const integer    localMode = mode - baseMode;

    for (found = 0, m.reset(); !found && m.more(); m.next()) {
      M     = m.current();
      found = M -> match (lambda2, betak2, N, method);
    }
    if (found) {
      _Msys[localMode] = M;
      if (method == DIRECT) { cout << '.'; cout.flush(); }
    } else {
      _Msys[localMode] =
	new MatrixSys (lambda2, betak2, mode, Elmt, Bsys, method);
      MS.add (_Msys[localMode]);
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
  integer       i;
  const integer N = _Msys.getSize();
  MatrixSys*    M;

  for (i = 0; i < N; i++) delete (M = MS.remove (_Msys[i]));
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
  const integer    verbose   = (integer) Femlib::value ("VERBOSE");
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
    real*           hbb  = work();
    real*           rmat = hbb  + sqr (next);
    real*           rwrk = rmat + sqr (np);
    integer*        ipiv = pivotmap();
    integer         info;

    _hbi    = new real*   [(size_t) _nel];
    _hii    = new real*   [(size_t) _nel];
    _bipack = new integer [(size_t) _nel];
    _iipack = new integer [(size_t) _nel];

    if (_nsolve) {
      _H = new real [(size_t) _npack];
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

      _hbi[j] = (nint) ? new real [(size_t) _bipack[j]] : 0;
      _hii[j] = (nint) ? new real [(size_t) _iipack[j]] : 0;

      elmt[j]->HelmholtzSC(lambda2,betak2,hbb,_hbi[j],_hii[j],rmat,rwrk,ipiv);

      for (i = 0; i < next; i++)
	if ((m = bmap[i]) < _nsolve)
	  for (k = 0; k < next; k++)
	    if ((n = bmap[k]) < _nsolve && n >= m)
	      _H[Lapack::band_addr (m, n, _nband)] +=
		hbb[Veclib::row_major (i, k, next)];

      Femlib::adopt (_bipack[j], _hbi + j);  
      Femlib::adopt (_iipack[j], _hii + j);
    }
    if (_nsolve) {
      // -- Loop over BCs and add diagonal contribution from mixed BCs.

      if (bsys -> mixBC()) {
	const integer  nbound = bsys -> nSurf();
	const integer* bmap   = _NS   -> btog();
	for (i = 0; i < nbound; i++)
	  _BC[i] -> augmentSC (_nband, _nsolve, bmap, rwrk, _H);
      }

      // -- Cholesky factor global banded-symmetric Helmholtz matrix.
    
      Lapack::pbtrf ("U",_nsolve,_nband-1,_H,_nband,info);
      Femlib::adopt (_npack, &_H);

      if (info) message (routine, "failed to factor Helmholtz matrix", ERROR);

      if (verbose) {
	real cond;
	pivotmap.setSize (_nsolve);  ipiv = pivotmap();
	work.setSize (3 * _nsolve);  rwrk = work();

	Lapack::pbcon ("U",_nsolve,_nband-1,_H,_nband,1.0,cond,rwrk,ipiv,info);
	cout << ", condition number: " << cond << endl;
      }
    }
  } break;

  case JACPCG: {
    const integer nbound = _BC.getSize();   
    real*         PCi;
    vector<real>  work (2 * npnp + np);
    real          *ed = work(), *ewrk = work() + npnp;

    _PC = new real [(size_t) _npts];
    Veclib::zero (_npts, _PC, 1);

    PCi  = _PC + _NS -> nGlobal();
    bmap = _NS -> btog();

    // -- Mixed BC contributions.

    if (bsys -> mixBC())
      for (i = 0; i < nbound; i++)
	_BC[i] -> augmentDg (bmap, _PC);

    // -- Element contributions.

    for (i = 0; i < _nel; i++, bmap += next, PCi += nint) {
      elmt[i] -> HelmholtzDg (lambda2, betak2, ed, ewrk);
      Veclib::scatr_sum (next, ed,  bmap,    _PC);
      Veclib::copy      (nint, ed + next, 1, PCi, 1);
    }

#if 1
    Veclib::vrecp (_npts, _PC, 1, _PC, 1);
#else  // -- Turn off preconditioner for testing.
    Veclib::fill  (_npts, 1.0, _PC, 1);
#endif

//    Femlib::adopt (_npts, &_PC);

  } break;

  default:
    message (routine, "no solver of type requested -- never happen", ERROR);
    break;
  }
}


integer MatrixSys::match (const real       lambda2,
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
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  if (fabs (_HelmholtzConstant - lambda2) < EPS                         &&
      fabs (_FourierConstant   - betak2 ) < EPS                         &&
      _NS -> nGlobal() == nScheme -> nGlobal()                          &&
      _NS -> nSolve()  == nScheme -> nSolve()                           &&
      Veclib::same (_NS->nGlobal(), _NS->btog(), 1, nScheme->btog(), 1) &&
      _method == method                                                  )
    return 1;

  else
    return 0;
}


MatrixSys::~MatrixSys()
// ---------------------------------------------------------------------------
// Destructor.  Because there may be aliases to the internal vector
// storage we use the Femlib family routines.
// ---------------------------------------------------------------------------
{
  switch (_method) {
  case JACPCG:
    Femlib::abandon (&_PC);
    break;
  case DIRECT: {
    integer i;
    for (i = 0; i < _nel; i++) {
      Femlib::abandon (_hbi + i);
      Femlib::abandon (_hii + i);
    }
    Femlib::abandon (&_H);
  } break;
  default:
    break;
  }
}


#if 0
ostream& operator << (ostream&      str,
		      MatrixSys& M  )
// ---------------------------------------------------------------------------
// Output a MatrixSys to file...one day we'll actually use this!
// ---------------------------------------------------------------------------
{
  char *hdr_fmt[] = {
    "-- Helmholtz MatrixSys Storage File --",
    "%-25d "    "Elements",
    "%-25d "    "Global matrix size (words)",
    "%-25d "    "Global matrix bandwidth",
    "%-25d "    "Global matrix singularity flag",
    "%-25.17e " "Helmholtz constant",
    "%-25d "    "Azimuthal constant",
    "%-25d "    "Word size (bytes)",
    "%-25s "    "Format"  
  };
  char    bufr[StrMax], fmt[StrMax];
  integer i, n;

  Veclib::describeFormat (fmt);
  
  sprintf (bufr, hdr_fmt[0]);                      str << bufr << endl;
  sprintf (bufr, hdr_fmt[1], M.nel);               str << bufr << endl;
  sprintf (bufr, hdr_fmt[2], M.npack);             str << bufr << endl;
  sprintf (bufr, hdr_fmt[3], M.nband);             str << bufr << endl;
  sprintf (bufr, hdr_fmt[4], M.singular);          str << bufr << endl;
  sprintf (bufr, hdr_fmt[5], M.HelmholtzConstant); str << bufr << endl;
  sprintf (bufr, hdr_fmt[6], sizeof (real));       str << bufr << endl;
  sprintf (bufr, hdr_fmt[7], fmt);                 str << bufr << endl;
  
  // -- Global Helmholtz matrix.

  str.write ((char*) M.H, M.npack * sizeof (real));

  // -- Elemental interior / exterior coupling matrices hbi.

  for (i = 0; i < M.nel; i++) {
    n = M.bipack[i];
    str.write ((char*) &n, sizeof (integer));
    str.write ((char*) M.hbi[i], n * sizeof (real));
  }
  
  // -- Elemental interior matrices hii.

  for (i = 0; i < M.nel; i++) {
    n = M.iipack[i];
    str.write ((char*) &n, sizeof (integer));
    str.write ((char*) M.hii[i], n * sizeof (real));
  }
  return str;
}


istream& operator >> (istream&      str,
		      MatrixSys& sys)
// ---------------------------------------------------------------------------
// Input a MatrixSys from file, with binary format conversion if required.
// ---------------------------------------------------------------------------
{
  char       routine[] = "MatrixSys::operator >>";
  char       bufr[StrMax], fmt[StrMax];
  istrstream s (bufr, strlen (bufr));
  integer    n, swab;
  real       f;
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  Veclib::describeFormat (fmt);
  
  str.getline (bufr);
  if (strcmp (bufr, "-- Helmholtz MatrixSys Storage File --"))
    message (routine, "input file lacks valid header",                  ERROR);

  str.getline (bufr);
  s >> n;
  if (n != hbi.getSize ())
    message (routine, "mismatch: number of elements in file & system",  ERROR);

  str.getline (bufr);
  s >> n;
  if (n != H.getSize ())
    message (routine, "mismatch: size of H in file & system",           ERROR);
  
  str.getline (bufr);
  s >> n;
  if (n != n_band)
    message (routine, "mismatch: bandwidth of H in file & system",      ERROR);
  
  str.getline (bufr);
  s >> n;
  if (n != singular)
    message (routine, "mismatch: singularity of H in file & system",    ERROR);

  str.getline (bufr);
  s >> f;
  if (fabs (f - HelmholtzConstant) > EPS)
    message (routine, "mismatch: Helmholtz constant in file & system",  ERROR);

  str.getline (bufr);
  s >> n;
  if (n != sizeof (real))
    message (routine, "mismatch: word size/precision in file & system", ERROR);

  str.getline (bufr);
  if (!strstr (bufr, "IEEE"))
  	message (routine, "unrecognized binary format", ERROR);
  	swab = (strstr (bufr, "little") && strstr (fmt, "big")  ||
		strstr (bufr, "big")    && strstr (fmt, "little"));

  str.read ((char*) H (), H.getSize () * sizeof (real));

  if (swab)   Veclib::brev (H.getSize(), H(), 1, H(), 1);
  
  register integer i;
  integer          n;
  const integer    N = hbi.getSize ();

  for (i = 0; i < N; i++) {
    str.read ((char*) &n, sizeof (integer));
    if (n != hbi[i].getSize ())
      message (routine, "mismatch: size of hbi in file & system", ERROR);
    else  {
      str.read ((char*) hbi[i], n * sizeof (real));
      if (swab)   Veclib::brev (n, hbi[i], 1, hbi[i], 1);
    }
    
    hbi[k] = FamilyMgr::insert (n, hbi[k]);  
  }
  
  for (i = 0; i < N; i++) {
    str.read ((char*) &n, sizeof (integer));
    if (n != hii[i].getSize ())
      message (routine, "mismatch: size of hii in file & system", ERROR);
    else  {
      str.read ((char*) hii[i], n * sizeof (real));
      if (swab)   Veclib::brev (n, hii[i], 1, hii[i], 1);
    }

    hii[k] = FamilyMgr::insert (n, hii[k]);
  }

  return str;
}
#endif

