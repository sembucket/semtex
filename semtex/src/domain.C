///////////////////////////////////////////////////////////////////////////////
// domain.C
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";


#include <Sem.h>
#include <time.h>


Domain::Domain (FEML&       F   ,
		Mesh&       M   ,
		BCmgr&      B   ,
		const char* flds,
		const char* sess)
// ---------------------------------------------------------------------------
// Construct a new Domain with all user Fields, and NumberSystems.
//
// No initialization of Field MatrixSystems.
// ---------------------------------------------------------------------------
{
  int       j, k, doff, boff, found;
  const int np   = (int) Femlib::value ("N_POLY");
  const int nz   = (int) Femlib::value ("N_Z");
  const int verb = (int) Femlib::value ("VERBOSE");
  const int NE   = M.nEl();
  const int NF   = strlen (flds);
  real*     z;

  cout << "-- " << NE << " elements, ";
  cout << np << "x" << np << "x" << nz << endl;

  name   = strdup (sess);
  fields = strdup (flds);
  step   = 0;
  time   = 0.0;

  Femlib::value ("t", time);

  if (verb) cout << "   Building Elements ... ";

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  Esys.setSize (NE);
  doff = 0;
  boff = 0;
  for (k = 0; k < NE; k++) {
    Esys[k] = new Element (k, M, z, np, doff, boff);
    doff += Esys[k] -> nTot();
    boff += Esys[k] -> nExt();
  }

  if (verb) cout << "done" << endl;
  
  if (verb) cout << "   Retrieving prebuilt numbering systems ... ";
  
  number();
  const int NS = Nsys.getSize();
  
  if (verb) cout << "done" << endl;
  
  if (verb) cout << "   Building Fields ... ";

  u.setSize (NF);

  for (k = 0; k < NF; k++) {
    found = 0;
    for (j = 0; !found && j < NS; j++)
      if (strchr (Nsys[j] -> fields(), flds[k])) {
	found = 1;
	u[k]  = new Field (F, B, Esys, Nsys[j], nz, flds[k]);
      }
  }

  if (verb) cout << " done" << endl;
}


void Domain::number ()
// ---------------------------------------------------------------------------
// Attempt to retreive numbering schemes (btog and bmsk values) from
// file "name.num".  If this doesn't exist, first try to create it by
// running "enumerate" utility.
// ---------------------------------------------------------------------------
{
  char         routine[] = "Domain::number";
  char         buf[StrMax], err[StrMax];
  register int i, j, nset;
  ifstream     num;

  strcat   (strcpy (buf, name), ".num");
  num.open (buf);

  if (!num) {
    sprintf (buf, "enumerate -O1 %s > %s.num", name, name);
    if (system (buf)) {
      sprintf (err, "couldn't open session file %s or %s.num", name, name);
      message (routine, err, ERROR);
    }
    
    strcat    (strcpy (buf, name), ".num");
    num.clear ();
    num.open  (buf);

    if (!num) {
      sprintf (err, "couldn't find or create number file %s", buf);
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "FIELDS")) {
    sprintf (err, "expected \"FIELDS\", found %s in %s.num", buf, name);
    message (routine, err, ERROR);
  }

  num >> buf >> buf;
  if (strcmp (buf, fields)) {
    sprintf (err, "Fields nominated in %s.num (\"%s\") don't match \"%s\"",
	     name, buf, fields);
    message (routine, err, ERROR);
  }

  num . getline (buf, StrMax) . getline (buf, StrMax);

  num >> buf >> nset >> buf >> buf >> buf;
    
  Nsys.setSize (nset);
  for (i = 0; i < nset; i++) {
    num >> buf;
    Nsys[i] = new NumberSystem;
    Nsys[i] -> ns_fields = strdup (buf);
  }

  for (i = 0; i < nset; i++)
    for (j = 0; j < nset; j++) {
      if (i == j) continue;
      if (strpbrk (Nsys[i] -> ns_fields, Nsys[j] -> ns_fields)) {
	sprintf (err, "Field name duplication: %s <--> %s",
		 Nsys[i] -> ns_fields, Nsys[j] -> ns_fields);
	message (routine, err, ERROR);
      }
    }

  num >> buf >> buf;
  if (strcmp (buf, "NEL")) {
    sprintf (err, "expected \"NEL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_nel;

  num >> buf >> buf;
  if (strcmp (buf, "NP_MAX")) {
    sprintf (err, "expected \"NP_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_np_max;

  num >> buf >> buf;
  if (strcmp (buf, "NEXT_MAX")) {
    sprintf (err, "expected \"NEXT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_ex_max;

  num >> buf >> buf;
  if (strcmp (buf, "NINT_MAX")) {
    sprintf (err, "expected \"NINT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_in_max;

  num >> buf >> buf;
  if (strcmp (buf, "NTOTAL")) {
    sprintf (err, "expected \"NTOTAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_ntotal;

  num >> buf >> buf;
  if (strcmp (buf, "NBOUNDARY")) {
    sprintf (err, "expected \"NBOUNDARY\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> Nsys[i] -> ns_nbndry;
    Nsys[i] -> ns_btog  = new int [Nsys[i] -> ns_nbndry];
    Nsys[i] -> ns_bmask = new int [Nsys[i] -> ns_nbndry];
  }

  for (i = 1; i < nset; i++)
    if (Nsys[i] -> ns_nbndry != Nsys[0] -> ns_nbndry) {
      sprintf (err, "NBOUNDARY values for systems 1 & %1d don't match", i + 1);
      message (routine, err, ERROR);
    }

  num >> buf >> buf;
  if (strcmp (buf, "NGLOBAL")) {
    sprintf (err, "expected \"NGLOBAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_nglobal;

  num >> buf >> buf;
  if (strcmp (buf, "NSOLVE")) {
    sprintf (err, "expected \"NSOLVE\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_nsolve;

  num . getline (buf, StrMax) . getline (buf, StrMax);

  num >> buf >> buf;
  if (strcmp (buf, "BANDWIDTH")) {
    sprintf (err, "expected \"BANDWIDTH\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++)
    num >> Nsys[i] -> ns_nbandw;

  num . getline (buf, StrMax) . getline (buf, StrMax). getline (buf, StrMax);

  for (i = 0; i < Nsys[0] -> ns_nbndry; i++) {
    num >> buf >> buf >> buf;
    for (j = 0; j < nset; j++) 
      num >> Nsys[j] -> ns_btog[i] >> Nsys[j] -> ns_bmask [i];
  }

  // -- At this point, all external data have been read in.
  //    Now create element-boundary mass smoothing vectors.

  Element*      E;
  NumberSystem* N;
  vector<real>  work (Nsys[0] -> ns_ex_max + Nsys[0] -> ns_in_max);
  real*         unity = work();

  for (j = 0; j < nset; j++) {
    N = Nsys[j];
    N -> ns_inv_mass = new real [N -> ns_nglobal];
    Veclib::zero (N -> ns_nglobal, N -> ns_inv_mass, 1);

    for (i = 0; i < N -> ns_nel; i++) {
      E = Esys[i];
      Veclib::fill (E -> nTot(), 1.0, unity, 1);
      E -> bndryDsSum (N -> ns_btog + E -> bOff(), unity, N -> ns_inv_mass);
    }

    Veclib::vrecp (N -> ns_nglobal, N -> ns_inv_mass, 1, N -> ns_inv_mass, 1);
  }

  // -- And emask vectors.

  for (j = 0; j < nset; j++) {
    N = Nsys[j];
    N -> ns_emask = new int [N -> ns_nel];

    for (i = 0; i < N -> ns_nel; i++) {
      E = Esys[i];
      N -> ns_emask[i] =
	Veclib::any (E -> nExt(), N -> ns_bmask + E -> bOff(), 1);
    }
  }
}


void Domain::initialize (const char* src)
// ---------------------------------------------------------------------------
// If a file name is given, attempt to use it to load the Domain Field data
// areas.  Otherwise, if a restart file "name".rst can be found, use it for
// input.  If these options fail, initialize all Fields to zero ICs.
//
// Carry out forwards Fourier transformation.
// ---------------------------------------------------------------------------
{
  char restartfile[StrMax];

  if   (src)         strcpy (restartfile, src);
  else       strcat (strcpy (restartfile, name), ".rst");

  ifstream file (restartfile);

  if (file) {
    cout << "-- Restarting from file:  " << restartfile;
    file >> *this;
    transform (+1);
  } else {
    cout << "-- Initializing solution with zero IC";
    for (int i(0); i < nField (); i++)
      *u[i] = 0.0;
  }

  cout << endl;

  real dt = Femlib::value ("D_T"   );
  int  ns = Femlib::value ("N_STEP");
  int  nt = Femlib::value ("N_TIME");

  real t  = time;
  Femlib::value ("t", t);

  cout << "   Start time       : " << t           << endl;
  cout << "   Time step        : " << dt          << endl;
  cout << "   Number of steps  : " << ns          << endl;
  cout << "   End time         : " << t + ns * dt << endl;
  cout << "   Integration order: " << nt          << endl;
}


void  Domain::dump (ofstream& output)
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
//
// Fields are inverse Fouier transformed prior to dumping in order to
// provide physical space values.
// ---------------------------------------------------------------------------
{
  int periodic = !(step %  (int) Femlib::value ("IO_FLD"));
  int final    =   step == (int) Femlib::value ("N_STEP");

  if ( ! (periodic || final) ) return;

  char  routine[] = "Domain::dump";
  int   verbose   = (int) Femlib::value ("VERBOSE");
  int   chkpoint  = (int) Femlib::value ("CHKPOINT");

  transform (-1);

  if (chkpoint) {
    char s[StrMax], b[StrMax], c[StrMax];

    output.close ();

    if (final)
      strcat (strcpy (s, name), ".fld");
    else {
      strcat (strcpy (s, name), ".chk");
      strcat (strcpy (b, name), ".chk.bak");
    
      sprintf (c, "mv ./%s ./%s", s, b);
      system  (c);
    }

    output.open (s);
    if (!output) message (routine, "can't open checkpoint file", ERROR);
  }

  if (verbose) message (routine, ": writing field dump", REMARK);

  output << *this;

  transform (+1);
}


void Domain::transform (const int sign)
// ---------------------------------------------------------------------------
// Fourier transform all Fields according to sign.
// ---------------------------------------------------------------------------
{
  int       i;
  const int N = nField ();
  
  for (i = 0; i < N; i++)
    u[i] -> transform (sign);
}


ostream& operator << (ostream& strm, Domain& D)
// ---------------------------------------------------------------------------
// Output all Domain field variables on ostream in prism-compatible form.
// Binary output only.
// ---------------------------------------------------------------------------
{
  char *hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields written\n",
    "%-25s "    "Format\n"
  };

  char   routine[] = "ostream << Domain";
  char   s1[StrMax], s2[StrMax];
  time_t tp   (::time (0));

  sprintf (s1, hdr_fmt[0], D.name);
  strm << s1;

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  strm << s1;

  D.u[0] -> describe (s2);
  sprintf (s1, hdr_fmt[2], s2);
  strm << s1;

  sprintf (s1, hdr_fmt[3], D.step);
  strm << s1;

  sprintf (s1, hdr_fmt[4], D.time);
  strm << s1;

  sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
  strm << s1;

  sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
  strm << s1;

  sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
  strm << s1;

  int k;
  for (k = 0; k < D.nField(); k++) s2[k] = D.u[k] -> name();
  s2[k] = '\0';
  sprintf (s1, hdr_fmt[8], s2);
  strm << s1;

  sprintf (s2, "binary ");
  Veclib::describeFormat (s2 + strlen (s2));
  sprintf (s1, hdr_fmt[9], s2);
  strm << s1;

  for (int n (0); n < D.nField (); n++) strm << *D.u[n];

  if (!strm) message (routine, "failed writing field file", ERROR);
  strm << flush;

  return strm;
}


istream& operator >> (istream& strm,
		      Domain&  D   )
// ---------------------------------------------------------------------------
// Input all Domain field variables from prism-compatible istream.
//
// Only binary storage format is allowed.  Check if conversion to native
// format (IEEE little/big-endian) is required.
//
// Ordering of fields in file is allowed to differ from that in D.
// ---------------------------------------------------------------------------
{
  char routine[] = "strm>>Domain";
  int  i, j, np, nz, nel, ntot, nfields;
  int  npchk,  nzchk, nelchk, swap = 0;
  char s[StrMax], f[StrMax], err[StrMax], fields[StrMax];

  if (strm.getline(s, StrMax).eof()) return strm;

  strm.getline(s,StrMax).getline(s,StrMax);

  D.u[0] -> describe (f);
  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != D.u[0] -> nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline(s,StrMax).getline(s,StrMax);

  istrstream (s, strlen (s)) >> D.time;
  Femlib::value ("t", D.time);
  
  strm.getline(s,StrMax).getline(s,StrMax).getline(s,StrMax).getline(s,StrMax);

  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';
  if (nfields != strlen (D.fields)) {
    sprintf (err, "file: %1d fields, Domain: %1d", nfields, strlen (D.fields));
    message (routine, err, ERROR);
  }
  for (i = 0; i < nfields; i++) 
    if (!strchr (D.fields, fields[i])) {
      sprintf (err, "field %c not present in Domain (%s)", fields[i],D.fields);
      message (routine, err, ERROR);
    }

  strm.getline (s, StrMax);
  Veclib::describeFormat (f);

  if (!strstr (s, "binary"))
    message (routine, "input field file not in binary format", ERROR);
  
  if (!strstr (s, "endian"))
    message (routine, "input field file in unknown binary format", WARNING);
  else {
    swap = (   (strstr (s, "big") && strstr (f, "little"))
	    || (strstr (f, "big") && strstr (s, "little")) );
    if (swap) cout << " (byte-swapping input fields)";
  }

  for (j = 0; j < nfields; j++) {
    for (i = 0; i < nfields; i++)
      if (D.fields[i] == fields[j]) break;
    strm >> *D.u[i];
    if (swap) D.u[i] -> reverse();
  }

  if (strm.bad()) message (routine, "failed reading field file", ERROR);
    
  return strm;
}
