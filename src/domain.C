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
// By convention, all Fields stored in the Domain have
// single-character lower-case names.  On input, the names of the
// Fields to be created are stored in the string "flds".  See the file
// field.C for significance of the names.
//
// No initialization of Field MatrixSystems.
// ---------------------------------------------------------------------------
{
  integer       k, doff, boff;
  const integer np   = Geometry::nP();
  const integer nz   = Geometry::nZ();
  const integer verb = (integer) Femlib::value ("VERBOSE");
  const integer NE   = M.nEl();
  const integer NF   = strlen (flds);
  const real*   z;

  ROOTONLY if (verb) {
    cout << NE << " elements, ";
    cout << np << "x" << np << "x" << nz << endl;
  }

  name   = strdup (sess);
  fields = strdup (flds);
  step   = 0;
  time   = 0.0;

  Femlib::value ("t", time);

  ROOTONLY if (verb) cout << "Building Elements ... ";

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  Esys.setSize (NE);
  doff = 0;
  boff = 0;
  for (k = 0; k < NE; k++) {
    Esys[k] = new Element (k, M, z, np, doff, boff);
    doff += Esys[k] -> nTot();
    boff += Esys[k] -> nExt();
  }

  ROOTONLY if (verb) cout << "done" << endl;
  ROOTONLY if (verb) cout << "Retrieving prebuilt numbering systems ... ";

  // -- Load all the available numbering systems created by enumerate.

  getNumber();

  // -- Form each Field with a Fourier-mode ordered vector of systems.

  const integer        NS      = Nsys.getSize();
  const NumberSystem** systems = new const NumberSystem* [3];

  ROOTONLY if (verb) cout << "done" << endl << "Building Fields ... ";

  u.setSize (NF);
  for (k = 0; k < NF; k++) {
    setNumber (fields[k], systems);
    u[k] = new Field (F, B, Esys, fields[k], systems);
  }

  ROOTONLY if (verb) cout << " done" << endl;
}


void Domain::getNumber ()
// ---------------------------------------------------------------------------
// Attempt to retrieve numbering schemes (btog and bmsk values) from
// file "name.num".  If this doesn't exist, first try to create it by
// running "enumerate" utility on root processor.
//
// The names of fields and their numbering schemes are significant.
// The convention employed is that the fields have lower-case
// single-character names.  Numbering schemes have the same names,
// *except* in the case of cylindrical coordinate systems where the
// domain includes the symmetry axis.  See file field.C for
// mode-related significance for upper-cased names of numbering
// schemes.
// ---------------------------------------------------------------------------
{
  const char       routine[] = "Domain::getNumber";
  char             buf[StrMax], err[StrMax];
  register integer i, j, nset;
  ifstream         num;

  strcat   (strcpy (buf, name), ".num");
  num.open (buf);

  if (!num) {

    ROOTONLY {
      sprintf (buf, "enumerate -O1 %s > %s.num", name, name);
      if (system (buf)) {
	sprintf (err, "couldn't open session file %s, or %s.num", name, name);
	message (routine, err, ERROR);
      }
    }
    
    Femlib::synchronize();

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
    sprintf (err, "couldn't locate \"FIELDS\" in %s.num", name);
    message (routine, err, ERROR);
  }

  num >> buf >> buf;
  j = strlen (buf);
  for (i = 0; i < j; i++)
    if (!strchr (fields, tolower (buf[j]))) {
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
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != Geometry::nElmt()) {
      sprintf (err, "mismatch in number of elements: %1d vs. %1d",
	       j, Geometry::nElmt());
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NP_MAX")) {
    sprintf (err, "expected \"NP_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != Geometry::nP()) {
      sprintf (err, "mismatch in element order: %1d vs. %1d",
	       j, Geometry::nP());
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NEXT_MAX")) {
    sprintf (err, "expected \"NEXT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> j;

  num >> buf >> buf;
  if (strcmp (buf, "NINT_MAX")) {
    sprintf (err, "expected \"NINT_MAX\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) num >> j;

  num >> buf >> buf;
  if (strcmp (buf, "NTOTAL")) {
    sprintf (err, "expected \"NTOTAL\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != Geometry::nPlane()) {
      sprintf (err, "mismatch in Field storage requirements: %1d vs %1d",
	       j, Geometry::nPlane());
      message (routine, err, ERROR);
    }
  }

  num >> buf >> buf;
  if (strcmp (buf, "NBOUNDARY")) {
    sprintf (err, "expected \"NBOUNDARY\", found %s", buf);
    message (routine, err, ERROR);
  }
  num >> buf;
  for (i = 0; i < nset; i++) {
    num >> j;
    if (j != Geometry::nBnode()) {
      sprintf (err, "mismatch in number of boundary nodes: %1d vs. %1d",
	       j, Geometry::nBnode());
      message (routine, err, ERROR);
    }
    Nsys[i] -> ns_btog  = new integer [Geometry::nBnode()];
    Nsys[i] -> ns_bmask = new integer [Geometry::nBnode()];
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

  for (i = 0; i < Geometry::nBnode(); i++) {
    num >> buf >> buf >> buf;
    for (j = 0; j < nset; j++) 
      num >> Nsys[j] -> ns_btog[i] >> Nsys[j] -> ns_bmask [i];
  }

  if (num.bad())
    message (routine, "failed reading to end of node-number file", ERROR);

  // -- At this point, all external data have been read in.
  //    Now create element-boundary mass smoothing vectors.
  //    Need to take care for cases (in cylindrical coords) where
  //    mass weighting will be zero where r = 0.

  Element*      E;
  NumberSystem* N;
  vector<real>  work (Geometry::nTotElmt());
  const real    EPS = (sizeof (real) == sizeof (double)) ? EPSDP : EPSSP;
  register real m;
  real*         unity = work();

  for (j = 0; j < nset; j++) {
    N = Nsys[j];
    N -> ns_inv_mass = new real [N -> ns_nglobal];
    Veclib::zero (N -> ns_nglobal, N -> ns_inv_mass, 1);

    for (i = 0; i < Geometry::nElmt(); i++) {
      E = Esys[i];
      Veclib::fill (E -> nTot(), 1.0, unity, 1);
      E -> bndryDsSum (N -> ns_btog + E -> bOff(), unity, N -> ns_inv_mass);
    }

    for (i = 0; i < N -> ns_nglobal; i++) {
      m = N -> ns_inv_mass[i];
      m = (m > EPS) ? 1.0 / m : 0.0;
      N -> ns_inv_mass[i] = m;
    }
  }

  // -- And emask vectors.

  for (j = 0; j < nset; j++) {
    N = Nsys[j];
    N -> ns_emask = new integer [Geometry::nElmt()];

    for (i = 0; i < Geometry::nElmt(); i++) {
      E = Esys[i];
      N -> ns_emask[i] =
	Veclib::any (E -> nExt(), N -> ns_bmask + E -> bOff(), 1);
    }
  }
}


void Domain::setNumber (const char           field ,
			const NumberSystem** system) const
// ---------------------------------------------------------------------------
// Set up the vector of NumberSystems that will be required to
// construct the named field.  This is rather messy since it encodes
// the rules for usage of NumberSystems both for Cartesian and
// cylindrical Fields, as discussed in field.C.
//
// For Cartesian (and for 2D cylindrical) geometries, only the first
// element of the systems vector is used, so the selection of systems
// is quite straightforward, based only on a match of input variable
// "field" and the names encoded in the Domain's internal vector of
// NumberSystems.
//
// For 3D Cylindrical geometries, there will be three entries in
// systems, corresponding to the numbering schemes for the 0th, 1st
// 2nd (and higher) Fourier modes.  The selection of appropriate
// NumberSystems depends on the input variable "field".  See remarks
// at head of file field.C.
// ---------------------------------------------------------------------------
{
  const char    routine[] = "Domain::setNumber";
  char          err[StrMax];
  integer       i;
  const integer N     = Nsys.getSize();
  const integer cyl3D = Geometry::system() == Geometry::Cylindrical &&
                        Geometry::nDim()   == 3;

  system[0] = system[1] = system[2] = 0;

  for (i = 0; i < N; i++)
    if (strchr (Nsys[i] -> fields(), field))
      system[0] = system[1] = system[2] = Nsys[i];

  if (cyl3D)
    for (i = 0; i < N; i++)
    switch (field) {
    case 'u':			// -- Axial velocity.
      if (strchr (Nsys[i] -> fields(), 'U')) system[1] = system[2] = Nsys[i];
      break;
    case 'v':			// -- Radial velocity (do nothing).
      break;
    case 'w':			// -- Azimuthal veolcity (coupled).
      if (strchr (Nsys[i] -> fields(), 'W')) system[1] = Nsys[i];
      break;
    case 'c':			// -- Scalar (optional).
      if (strchr (Nsys[i] -> fields(), 'C')) system[1] = system[2] = Nsys[i];
      break;
    case 'p':			// -- Pressure.
      if (strchr (Nsys[i] -> fields(), 'P')) system[1] = system[2] = Nsys[i];
      break;
    default:
      sprintf (err, "unrecognized Field name: %c", field);
      message (routine, err, ERROR);
      break;
    }

  if (!(system[0] && system[1] && system[2])) {
    sprintf (err, "unsuccessful with field %c", field);
    message (routine, err, ERROR);
  }
}


void Domain::report ()
// ---------------------------------------------------------------------------
// Print a run-time summary of domain & timestep information on cout.
// ---------------------------------------------------------------------------
{
  const real    t   =           time;
  const real    dt  =           Femlib::value ("D_T");
  const real    lz  =           Femlib::value ("TWOPI / BETA");
  const integer ns  = (integer) Femlib::value ("N_STEP");
  const integer nt  = (integer) Femlib::value ("N_TIME");
  const integer chk = (integer) Femlib::value ("CHKPOINT");
  const integer per = (integer) Femlib::value ("IO_FLD");

  cout << "-- Coordinate system       : ";
  if (Geometry::system() == Geometry::Cylindrical)
    cout << "cylindrical" << endl;
  else
    cout << "Cartesian" << endl;

  cout << "   Solution fields         : " << fields             << endl;

  cout << "   Number of elements      : " << Geometry::nElmt()  << endl;
  cout << "   Number of planes        : " << Geometry::nZ()     << endl;
  cout << "   Number of processors    : " << Geometry::nProc()  << endl;
  if (Geometry::nZ() > 1) cout << "   Periodic length         : " << lz<< endl;
  cout << "   Polynomial order (np-1) : " << Geometry::nP() - 1 << endl;

  cout << "   Time integration order  : " << nt                 << endl;
  cout << "   Start time              : " << t                  << endl;
  cout << "   Finish time             : " << t + ns * dt        << endl;
  cout << "   Time step               : " << dt                 << endl;
  cout << "   Number of steps         : " << ns                 << endl;
  cout << "   Dump interval (steps)   : " << per;
  if (chk) cout << " (checkpoint)";  
  cout << endl;
}


void Domain::initialize ()
// ---------------------------------------------------------------------------
// If a restart file "name".rst can be found, use it for input.  If
// this fails, initialize all Fields to zero ICs.
//
// Carry out forwards Fourier transformation, zero Nyquist data.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer nF = nField();
  char          restartfile[StrMax];
  
  ROOTONLY cout << "-- Initial condition       : ";
  ifstream file (strcat (strcpy (restartfile, name), ".rst"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << restartfile;
      cout.flush();
    }
    file >> *this;
    file.close();
    transform (+1);
    ROOTONLY for (i = 0; i < nF; i++) u[i] -> zeroNyquist();
  } else {
    ROOTONLY cout << "set to zero";
    for (i = 0; i < nF; i++) *u[i] = 0.0;
  }

  ROOTONLY cout << endl;
  
  Femlib::value ("t", time);
  step = 0;
}


void Domain::dump ()
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
//
// Fields are inverse Fourier transformed prior to dumping in order to
// provide physical space values.
// ---------------------------------------------------------------------------
{
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;
  ofstream output;
  
  ROOTONLY {
    const char    routine[] = "Domain::dump";
    char          dumpfl[StrMax], backup[StrMax], command[StrMax];
    const integer verbose   = (integer) Femlib::value ("VERBOSE");
    const integer chkpoint  = (integer) Femlib::value ("CHKPOINT");

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, name), ".fld");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, name), ".chk");
	if (!initial) {
	  strcat  (strcpy (backup, name), ".chk.bak");
	  sprintf (command, "mv ./%s ./%s", dumpfl, backup);
	  system  (command);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, name), ".fld");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }
    
    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }

  transform (-1);
  output << *this;
  transform (+1);

  ROOTONLY output.close();
}


void Domain::transform (const integer sign)
// ---------------------------------------------------------------------------
// Fourier transform all Fields according to sign.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer N = nField ();
  
  for (i = 0; i < N; i++) u[i] -> transform (sign);
}


ofstream& operator << (ofstream& strm,
		       Domain&   D   )
// ---------------------------------------------------------------------------
// Output all Domain field variables on ostream in prism-compatible
// form.  Binary output only.  Note that output is only done on root
// processor.
// ---------------------------------------------------------------------------
{
  const char *hdr_fmt[] = { 
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

  const char routine[] = "ostream << Domain";
  char       s1[StrMax], s2[StrMax];
  time_t     tp (::time (0));
  integer    k;

  ROOTONLY {
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

    for (k = 0; k < D.nField(); k++) s2[k] = D.u[k] -> name();
    s2[k] = '\0';
    sprintf (s1, hdr_fmt[8], s2);
    strm << s1;

    sprintf (s2, "binary ");
    Veclib::describeFormat (s2 + strlen (s2));
    sprintf (s1, hdr_fmt[9], s2);
    strm << s1;
  }

  for (k = 0; k < D.nField(); k++) strm << *D.u[k];

  ROOTONLY {
    if (!strm) message (routine, "failed writing field file", ERROR);
    strm << flush;
  }

  return strm;
}


ifstream& operator >> (ifstream& strm,
		       Domain&   D   )
// ---------------------------------------------------------------------------
// Input all Domain field variables from prism-compatible istream.
//
// Only binary storage format is allowed.  Check if conversion to
// native format (IEEE little/big-endian) is required.
//
// Ordering of fields in file is allowed to differ from that in D.
// ---------------------------------------------------------------------------
{
  const char routine[] = "strm>>Domain";
  integer    i, j, np, nz, nel, ntot, nfields;
  integer    npchk,  nzchk, nelchk;
  integer    swap = 0, verb = (integer) Femlib::value ("VERBOSE");
  char       s[StrMax], f[StrMax], err[StrMax], fields[StrMax];

  if (strm.getline(s, StrMax).eof()) return strm;

  strm.getline(s,StrMax).getline(s,StrMax);
  
  D.u[0] -> describe (f);
  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline(s,StrMax);
  istrstream (s, strlen (s)) >> D.step;

  strm.getline(s,StrMax);
  istrstream (s, strlen (s)) >> D.time;
  Femlib::value ("t", D.time);
  
  strm.getline(s,StrMax).getline(s,StrMax);
  strm.getline(s,StrMax).getline(s,StrMax);

  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';
  if (nfields != strlen (D.fields)) {
    sprintf (err, "file: %1d fields, Domain: %1d", nfields,strlen(D.fields));
    message (routine, err, ERROR);
  }
  for (i = 0; i < nfields; i++) 
    if (!strchr (D.fields, fields[i])) {
      sprintf (err,"field %c not present in Domain (%s)",fields[i],D.fields);
      message (routine, err, ERROR);
    }

  strm.getline (s, StrMax);
  Veclib::describeFormat (f);

  if (!strstr (s, "binary"))
    message (routine, "input field file not in binary format", ERROR);
  
  if (!strstr (s, "endian"))
    message (routine, "input field file in unknown binary format", WARNING);
  else {
    swap = ((strstr (s, "big") && strstr (f, "little")) ||
	    (strstr (f, "big") && strstr (s, "little")) );
    ROOTONLY {
      if (swap) cout << " (byte-swapping)";
      cout.flush();
    }
  }

  for (j = 0; j < nfields; j++) {
    for (i = 0; i < nfields; i++)
      if (fields[j] == D.fields[i]) {
	strm >>  *D.u[i];
	if (swap) D.u[i] -> reverse();
	break;
      }
  }
    
  ROOTONLY if (strm.bad())
    message (routine, "failed reading field file", ERROR);
    
  return strm;
}




