///////////////////////////////////////////////////////////////////////////////
// matrix.C: routines for direct solution of Helmholtz problems.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <Sem.h>


ModalMatrixSystem::ModalMatrixSystem (const real              lambda2,
				      const real              beta   ,
				      const int               nmodes ,
				      const vector<Element*>& Elmt   ,
				      const NumberSystem*     Nsys   ,
				      const MatrixGenerator   Mgen   )
// ---------------------------------------------------------------------------
// Generate the vector of MatrixSystems which will be used to solve all
// the Fourier-mode discrete Helmholtz problems for the associated scalar
// Fields (called out by names).
// ---------------------------------------------------------------------------
{
  char routine[] = "ModalMatrixSystem::ModalMatrixSystem";
  int  k;
  real HelmholtzConstant;

  fields = strdup (Nsys -> fields());
  Msys.setSize (nmodes);

  cout << routine << ": building matrices for Fields \"" << fields << "\" [";
  for (k = 0; k < nmodes; k++) {
    HelmholtzConstant = lambda2 + sqr (k * beta);
    Msys[k] = new MatrixSystem (HelmholtzConstant, Elmt, Nsys, Mgen);
  }

  cout << "]" << endl;
}


MatrixSystem::MatrixSystem (const real              HCon,
			    const vector<Element*>& Elmt,
			    const NumberSystem*     Nsys,
			    const MatrixGenerator   MGen) :
			    
			    HelmholtzConstant      (HCon),
			    nel                    (Nsys -> nEl()),
			    nband                  (Nsys -> nBand())
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

  vector<real> work (Nsys -> neMax() * Nsys -> neMax() +
		     Nsys -> npMax() * Nsys -> npMax() +
		     Nsys -> neMax() * Nsys -> ntMax() );

  singular = fabs (HelmholtzConstant) < EPS && !Nsys -> fmask();
  nsolve   = (singular) ? Nsys -> nSolve() - 1 : Nsys -> nSolve();

  if (nsolve) {
    npack = nband * nsolve; // -- Size for LAPACK banded format.
    H     = new real [npack];
    Veclib::zero (npack, H, 1);

    if (verbose)
      cout << routine
	   << " : Banded system matrix: "
	   << nsolve
	   << "x"
	   << nband
	   << "\t("
	   << npack
	   << " words)"
	   << endl;
  }

  // -- Loop over elements, creating & posting elemental Helmholtz matrices.

  hbi    = new real* [nel];
  hii    = new real* [nel];
  bipack = new int   [nel];
  iipack = new int   [nel];

  hbb  = work();
  rmat = hbb  + Nsys -> neMax() * Nsys -> neMax();
  rwrk = rmat + Nsys -> npMax() * Nsys -> npMax();

  for (j = 0; j < nel; j++) {
    E = Elmt[j];

    next      = E -> nExt();
    nint      = E -> nInt();
    bipack[j] = next * nint;
    iipack[j] = ((nint + 1) * nint) >> 1; // -- Size for LAPACK packed format.

    hbi[j] = (nint) ? new real [bipack[j]] : 0;
    hii[j] = (nint) ? new real [iipack[j]] : 0;

    (E ->* MGen) (HelmholtzConstant, hbb, hbi[j], hii[j], rmat, rwrk);

    bmap = Nsys -> btog() + E -> bOff();

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
    cout << "System condition number:                   " << cond << endl;
  }

  cout << '*';
}


ostream& operator << (ostream&      str,
		      MatrixSystem& M  )
// ---------------------------------------------------------------------------
// Output a MatrixSystem to file.
// ---------------------------------------------------------------------------
{
  char *hdr_fmt[] = {
    "-- Helmholtz MatrixSystem Storage File --",
    "%-25d "    "Elements",
    "%-25d "    "Global matrix size (words)",
    "%-25d "    "Global matrix bandwidth",
    "%-25d "    "Global matrix singularity flag",
    "%-25.17e " "Helmholtz constant",
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
  sprintf (bufr, hdr_fmt[6], M.HelmholtzConstant); str << bufr << endl;
  sprintf (bufr, hdr_fmt[5], sizeof (real));       str << bufr << endl;
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
  char       bufr[StrMax];
  istrstream s (bufr, strlen (bufr));
  int        n;
  real       f;
  const real EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;

  strm.getline (bufr);
  if (strcmp (bufr, "-- Helmholtz MatrixSystem Storage File --"))
    message (routine, "input file lacks valid header",                  ERROR);

  strm.getline (bufr);
  s >> n;
  if (n != hbi.getSize ())
    message (routine, "mismatch: number of elements in file & system",  ERROR);

  strm.getline (bufr);
  s >> n;
  if (n != H.getSize ())
    message (routine, "mismatch: size of H in file & system",           ERROR);
  
  strm.getline (bufr);
  s >> n;
  if (n != n_band)
    message (routine, "mismatch: bandwidth of H in file & system",      ERROR);
  
  strm.getline (bufr);
  s >> n;
  if (n != singular)
    message (routine, "mismatch: singularity of H in file & system",    ERROR);

  strm.getline (bufr);
  s >> n;
  if (n != sizeof (real))
    message (routine, "mismatch: word size/precision in file & system", ERROR);

  strm.getline (bufr);
  s >> f;
  if (fabs (f - HelmholtzConstant) > EPS)
    message (routine, "mismatch: Helmholtz constant in file & system",  ERROR);


  str.read ((char*) H (), H.getSize () * sizeof (real));

  register int i;
  int          n;
  const int    N = hbi.getSize ();

  for (i = 0; i < N; i++) {
    str.read ((char*) &n, sizeof (int));
    if (n != hbi[i].getSize ())
      message (routine, "mismatch: size of hbi in file & system", ERROR);
    else
      str.read ((char*) hbi[i], n * sizeof (real));
    hbi[k] = FamilyMgr::insert (n, hbi[k]);  
  }
  
  for (i = 0; i < N; i++) {
    str.read ((char*) &n, sizeof (int));
    if (n != hii[i].getSize ())
      message (routine, "mismatch: size of hii in file & system", ERROR);
    else    
      str.read ((char*) hii[i], n * sizeof (real));
    hii[k] = FamilyMgr::insert (n, hii[k]);
  }
#endif
  
  return str;
}

