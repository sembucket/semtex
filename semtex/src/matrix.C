///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines for direct solution of Helmholtz problems.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>

static List<MatrixSys*> MS;


ModalMatrixSys::ModalMatrixSys (const real              lambda2 ,
				const real              beta    ,
				const integer           baseMode,
				const integer           numModes,
				const vector<Element*>& Elmt    ,
				const BoundarySys*      Bsys    )
// ---------------------------------------------------------------------------
// Generate or retrieve from internal database MS the vector of
// MatrixSyss which will be used to solve all the Fourier-mode
// discrete Helmholtz problems for the associated scalar Fields
// (called out by names).
//
// Input variables:
//   lambda2: Helmholtz constant for the problem,	
//   beta   : Fourier length scale = TWOPI / Lz,
//   nmodes : number of Fourier modes which will be solved,
//   Elmt   : vector of Element*'s used to make local Helmholtz matrices,
//   Bsys   : boundary system for this field.
// ---------------------------------------------------------------------------
{
  const char               name = Bsys -> field();
  integer                  found, mode;
  real                     betak2;
  ListIterator<MatrixSys*> m (MS);
  MatrixSys*               M;

  fields = new char [strlen (Bsys -> Nsys (0) -> fields()) + 1];
  strcpy (fields, Bsys -> Nsys (0) -> fields());
  Msys.setSize (numModes);

  ROOTONLY cout << "-- Installing matrices for field '" << name << "' [";
  cout.flush();

  Femlib::synchronize();

  for (mode = baseMode; mode < baseMode + numModes; mode++) {
    const NumberSys* N      = Bsys -> Nsys (mode);
    const real       betak2 = sqr (Field::modeConstant (name, mode, beta));
    found  = 0;
    for (m.reset(); !found && m.more(); m.next()) {
      M     = m.current();
      found = M -> match (lambda2, betak2, N);
    }
    if (found) {
      Msys[mode] = M;
      cout << '.';  cout.flush();
    } else {
      Msys[mode] = new MatrixSys (lambda2, betak2, mode, Elmt, Bsys);
      MS.add (Msys[mode]);
      cout << '*'; cout.flush();
    }
  }

  Femlib::synchronize();

  ROOTONLY cout << "]" << endl;
  cout.flush();
}


MatrixSys::MatrixSys (const real               lambda2,
		      const real               betak2 ,
		      const integer            mode   ,
		      const vector<Element*>&  elmt   ,
		      const BoundarySys*       bsys   ) :
// ---------------------------------------------------------------------------
// Initialize and factorise matrices in this system.
//
// Matrices are assembled using LAPACK-compatible ordering systems.
// Global Helmholtz matrix uses symmetric-banded format; elemental
// Helmholtz matrices (hii & hbi) use column-major formats.
// ---------------------------------------------------------------------------
// NB: these get evaluated in the order they appear in the class
// definition!:
  HelmholtzConstant (lambda2),
  FourierConstant   (betak2 ),
  BC                (bsys -> BCs  (mode)),
  NS                (bsys -> Nsys (mode)),
  nel               (Geometry::nElmt()),
  nband             (NS -> nBand()),
  singular          ((HelmholtzConstant + FourierConstant) < EPSSP &&
		     !NS -> fmask() && !bsys -> mixBC()),
  nsolve            ((singular) ? NS -> nSolve() - 1 : NS -> nSolve()),
  npack             (nband * nsolve)
{
  const char        routine[] = "MatrixSys::MatrixSys";
  const integer     verbose = (integer) Femlib::value ("VERBOSE");
  const integer     np      = Geometry::nP();
  const integer     next    = Geometry::nExtElmt();
  const integer     nint    = Geometry::nIntElmt();
  const integer     npnp    = Geometry::nTotElmt();
  const integer*    bmap;
  vector<real>      work     (sqr (next) + sqr (np) + sqr (npnp));
  vector<integer>   pivotmap (nint);
  real*             hbb  = work();
  real*             rmat = hbb  + sqr (next);
  real*             rwrk = rmat + sqr (np);
  integer*          ipiv = pivotmap();
  register integer  i, j, k, m, n;
  integer           info;

  hbi    = new real*   [(size_t) nel];
  hii    = new real*   [(size_t) nel];
  bipack = new integer [(size_t) nel];
  iipack = new integer [(size_t) nel];

  if (nsolve) {
    H = new real [(size_t) npack];
    Veclib::zero (npack, H, 1);

    if (verbose > 1)
      cout << endl
	   << "Helmholtz constant (lambda2): " << setw(10) << lambda2
	   << ", Fourier constant (betak2): "  << setw(10) << betak2;
    if (verbose)
      cout << endl << "System matrix: " << nsolve << "x" << nband
	   << "\t(" << npack << " words)";
  }

  // -- Loop over elements, creating & posting elemental Helmholtz matrices.

  for (bmap = NS -> btog(), j = 0; j < nel; j++, bmap += next) {
    bipack[j] = next * nint;
    iipack[j] = nint * nint;

    hbi[j] = (nint) ? new real [(size_t) bipack[j]] : 0;
    hii[j] = (nint) ? new real [(size_t) iipack[j]] : 0;

    elmt[j] -> HelmholtzSC (lambda2,betak2, hbb,hbi[j],hii[j], rmat,rwrk,ipiv);

    for (i = 0; i < next; i++)
      if ((m = bmap[i]) < nsolve)
	for (k = 0; k < next; k++)
	  if ((n = bmap[k]) < nsolve && n >= m)
	    H[Lapack::band_addr (m, n, nband)] +=
	      hbb[Veclib::row_major (i, k, next)];

    Femlib::adopt (bipack[j], hbi + j);  
    Femlib::adopt (iipack[j], hii + j);
  }

  // -- Loop over BCs and add diagonal contribution from mixed BCs.
#if 1
  if (bsys -> mixBC()) {
    const integer      nbound = bsys -> nSurf();
    const integer*     bmap   = NS   -> btog();
    for (i = 0; i < nbound; i++)
      BC[i] -> augmentSC (nband, nsolve, bmap, rwrk, H);
  }
#endif
  // -- Cholesky factor global banded-symmetric Helmholtz matrix.

  Lapack::pbtrf ("U", nsolve, nband - 1, H, nband, info);

  if (info) message (routine, "failed to factor Helmholtz matrix", ERROR);

  if (verbose) {
    real cond;
    pivotmap.setSize (nsolve);  ipiv = pivotmap();
    work.setSize (3 * nsolve);  rwrk = work();

    Lapack::pbcon ("U", nsolve, nband-1, H, nband, 1.0,cond, rwrk,ipiv,info);
    cout << ", condition number: " << cond << endl;
  }
}


integer MatrixSys::match (const real       lambda2,
			  const real       betak2 ,
			  const NumberSys* nScheme) const
// ---------------------------------------------------------------------------
// The unique identifiers of a MatrixSys are presumed to be given
// by the constants and the numbering system used.  Other things that
// could be checked but aren't (yet) include geometric systems and
// quadrature schemes.
// ---------------------------------------------------------------------------
{
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  if (fabs (HelmholtzConstant - lambda2) < EPS  &&
      fabs (FourierConstant   - betak2 ) < EPS  &&
      NS -> nGlobal() == nScheme -> nGlobal() &&
      NS -> nSolve()  == nScheme -> nSolve()  &&
      Veclib::same (NS->nGlobal(), NS->btog(), 1, nScheme->btog(), 1))
    return 1;

  else
    return 0;
}


ostream& operator << (ostream&      str,
		      MatrixSys& M  )
// ---------------------------------------------------------------------------
// Output a MatrixSys to file.
// ---------------------------------------------------------------------------
{
#if 0
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
#endif
  return str;
}


istream& operator >> (istream&      str,
		      MatrixSys& sys)
// ---------------------------------------------------------------------------
// Input a MatrixSys from file, with binary format conversion if required.
// ---------------------------------------------------------------------------
{
#if 0
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

#endif

  return str;
}
