///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines for direct solution of Helmholtz problems.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>

static List<MatrixSystem*> MS;


ModalMatrixSystem::ModalMatrixSystem (const real              lambda2 ,
				      const real              beta    ,
				      const char              name    ,
				      const int               baseMode,
				      const int               numModes,
				      const vector<Element*>& Elmt    ,
				      const NumberSystem**    Nsys    )
// ---------------------------------------------------------------------------
// Generate or retrieve from internal database MS the vector of MatrixSystems
// which will be used to solve all the Fourier-mode discrete Helmholtz
// problems for the associated scalar Fields (called out by names).
//
// Input variables:
//   lambda2: Helmholtz constant for the problem,	
//   beta   : Fourier length scale = TWOPI / Lz,
//   nmodes : number of Fourier modes which will be solved,
//   Elmt   : vector of Element*'s used to make local Helmholtz matrices,
//   Nsys   : truncated modal vector of numbering systems.
// ---------------------------------------------------------------------------
{
  int                         k, trunc, found;
  real                        betak2;
  ListIterator<MatrixSystem*> m (MS);
  MatrixSystem*               M;

  fields = strdup (Nsys[0] -> fields());
  Msys.setSize (numModes);
  cout << "-- Building matrices for Fields \"" << fields << "\"\t[";

  for (k = baseMode; k < numModes; k++) {
    trunc   = min (k, 2);
    betak2  = sqr (Field::modeConstant (name, k) * beta);
    found   = 0;
    for (m.reset(); !found && m.more(); m.next()) {
      M     = m.current();
      found = M -> match (lambda2, betak2, Nsys[trunc]);
    }
    if (found) {
      cout << '.';
      Msys[k] = M;
    } else {
      cout << '*';
      Msys[k] = new MatrixSystem (lambda2, betak2, Elmt, Nsys[trunc]);
      MS.add (Msys[k]);
    }
  }

  cout << "]" << endl;
}


MatrixSystem::MatrixSystem (const real              lambda2,
			    const real              betak2 ,
			    const vector<Element*>& Elmt   ,
			    const NumberSystem*     Nsys   ) :
			    
			    HelmholtzConstant      (lambda2),
			    FourierConstant        (betak2 ),
			    nel                    (Geometry::nElmt()),
			    nband                  (Nsys -> nBand()),
			    N                      (Nsys)
// ---------------------------------------------------------------------------
// Initialize and factorize matrices in this system.
//
// Matrices are assembled using LAPACK-compatible ordering systems.
// Global Helmholtz matrix uses     symmetric-banded format,
// elemental Helmholtz matrices use symmetric-packed (hii)
//                              and column-major     (hbi) formats.
// ---------------------------------------------------------------------------
{
  char         routine[] = "MatrixSystem::MatrixSystem";
  const int    verbose = (int) Femlib::value ("VERBOSE");
  const real   EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  
  register int i, j, k, m, n;
  const int*   bmap;
  int          next, nint, info;
  real         *hbb, *rmat, *rwrk;
  Element*     E;

  vector<real> work (sqr (Geometry::nExtElmt()) + sqr (Geometry::nP()) +
		     Geometry::nExtElmt() * Geometry::nTotElmt());

  hbb      = work();
  rmat     = hbb  + sqr (Geometry::nExtElmt());
  rwrk     = rmat + sqr (Geometry::nP());

  singular = fabs (HelmholtzConstant+FourierConstant) < EPS && !N-> fmask();
  nsolve   = (singular) ? N -> nSolve() - 1 : N -> nSolve();

  if (nsolve) {
    npack = nband * nsolve; // -- Size for LAPACK banded format.
    H     = new real [npack];
    Veclib::zero (npack, H, 1);

    if (verbose)
      cout << endl
	   << "   System matrix: "
	   << nsolve
	   << "x"
	   << nband
	   << "\t("
	   << npack
	   << " words)";
  }

  // -- Loop over elements, creating & posting elemental Helmholtz matrices.

  hbi    = new real* [nel];
  hii    = new real* [nel];
  bipack = new int   [nel];
  iipack = new int   [nel];

  for (j = 0; j < nel; j++) {
    E = Elmt[j];

    next      = E -> nExt();
    nint      = E -> nInt();
    bipack[j] = next * nint;
    iipack[j] = ((nint + 1) * nint) >> 1; // -- Size for LAPACK packed format.

    hbi[j] = (nint) ? new real [bipack[j]] : 0;
    hii[j] = (nint) ? new real [iipack[j]] : 0;

    E -> HelmholtzSC (lambda2, betak2, hbb, hbi[j], hii[j], rmat, rwrk);

    bmap = N -> btog() + E -> bOff();

    for (i = 0; i < next; i++)
      if ((m = bmap[i]) < nsolve)
	for (k = 0; k < next; k++)
	  if ((n = bmap[k]) < nsolve && n >= m)
	    H[Lapack::band_addr (m, n, nband)] +=
	      hbb[Veclib::row_major (i, k, next)];

    Femlib::adopt (bipack[j], hbi + j);  
    Femlib::adopt (iipack[j], hii + j);
  }

  // -- Cholesky factor global banded-symmetric Helmholtz matrix.

  Lapack::pbtrf ("U", nsolve, nband - 1, H, nband, info);

  if (info) message (routine, "failed to factor Helmholtz matrix", ERROR);

  if (verbose) {
    real        cond;
    vector<int> iwrk (nsolve);
    work.setSize (3 * nsolve);
    rwrk = work();

    Lapack::pbcon ("U",nsolve,nband-1,H,nband,1.0,cond, rwrk,iwrk(),info);
    cout << ", condition number: " << cond << endl;
  }
}


int MatrixSystem::match (const real          lambda2,
			 const real          betak2 ,
			 const NumberSystem* nScheme) const
// ---------------------------------------------------------------------------
// The unique identifiers of a MatrixSystem are presumed to be given by the
// constants and the numbering system used.  Other things that could be
// checked but aren't (yet) include geometric systems and quadrature schemes.
// ---------------------------------------------------------------------------
{
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  if (fabs (HelmholtzConstant - lambda2) < EPS &&
      fabs (FourierConstant   - betak2 ) < EPS &&
      N -> nGlobal() == nScheme -> nGlobal()   &&
      N -> nSolve()  == nScheme -> nSolve()    &&
      Veclib::same (N -> nGlobal(), N -> btog(), 1, nScheme -> btog(), 1))
    return 1;

  else
    return 0;
}


ostream& operator << (ostream&      str,
		      MatrixSystem& M  )
// ---------------------------------------------------------------------------
// Output a MatrixSystem to file.
// ---------------------------------------------------------------------------
{
#if 0
  char *hdr_fmt[] = {
    "-- Helmholtz MatrixSystem Storage File --",
    "%-25d "    "Elements",
    "%-25d "    "Global matrix size (words)",
    "%-25d "    "Global matrix bandwidth",
    "%-25d "    "Global matrix singularity flag",
    "%-25.17e " "Helmholtz constant",
    "%-25d "    "Azimuthal constant",
    "%-25d "    "Word size (bytes)",
    "%-25s "    "Format"  
  };
  char bufr[StrMax], fmt[StrMax];
  int  i, n;

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
    str.write ((char*) &n, sizeof (int));
    str.write ((char*) M.hbi[i], n * sizeof (real));
  }
  
  // -- Elemental interior matrices hii.

  for (i = 0; i < M.nel; i++) {
    n = M.iipack[i];
    str.write ((char*) &n, sizeof (int));
    str.write ((char*) M.hii[i], n * sizeof (real));
  }
#endif
  return str;
}


istream& operator >> (istream&      str,
		      MatrixSystem& sys)
// ---------------------------------------------------------------------------
// Input a MatrixSystem from file, with binary format conversion if required.
// ---------------------------------------------------------------------------
{
#if 0
  char       routine[] = "MatrixSystem::operator >>";
  char       bufr[StrMax], fmt[StrMax];
  istrstream s (bufr, strlen (bufr));
  int        n, swab;
  real       f;
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  Veclib::describeFormat (fmt);
  
  str.getline (bufr);
  if (strcmp (bufr, "-- Helmholtz MatrixSystem Storage File --"))
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
  
  register int i;
  int          n;
  const int    N = hbi.getSize ();

  for (i = 0; i < N; i++) {
    str.read ((char*) &n, sizeof (int));
    if (n != hbi[i].getSize ())
      message (routine, "mismatch: size of hbi in file & system", ERROR);
    else  {
      str.read ((char*) hbi[i], n * sizeof (real));
      if (swab)   Veclib::brev (n, hbi[i], 1, hbi[i], 1);
    }
    
    hbi[k] = FamilyMgr::insert (n, hbi[k]);  
  }
  
  for (i = 0; i < N; i++) {
    str.read ((char*) &n, sizeof (int));
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











