///////////////////////////////////////////////////////////////////////////////
// domain.C: implement domain class functions.
//
// Copyright (c) 1994,2003 Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <Sem.h>


Domain::Domain (FEML*             F,
		vector<Element*>& E,
		BCmgr*            B) :
// ---------------------------------------------------------------------------
// Construct a new Domain with all user Fields.
//
// By convention, all Fields stored in the Domain have
// single-character lower-case names.  On input, the names of the
// Fields to be created are stored in the string "flds".  See the file
// field.C for significance of the names.
//
// No initialization of Field MatrixSystems.
// ---------------------------------------------------------------------------
  elmt (E)
{
  const integer verbose = (integer) Femlib::value ("VERBOSE");
  integer       i, nfield;
  const integer nz   = Geometry::nZProc();
  const integer ntot = Geometry::nTotProc();
  real*         alloc;

  strcpy ((name = new char [strlen (F -> root()) + 1]), F -> root());
  Femlib::value ("t", time = 0.0);
  step = 0;

  strcpy ((field = new char [strlen (B -> field()) + 1]), B -> field());
  nfield = strlen (field);
  
  VERBOSE cout << "  Domain will contain fields: " << field << endl;

  // -- Build boundary system and field for each variable.
  
  VERBOSE cout << "  Building domain boundary systems ... ";

  b.resize (nfield);
  for (i = 0; i < nfield; i++) b[i] = new BoundarySys (B, E, field[i]);

  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Building domain fields ... ";

  u   .resize (nfield);
  udat.resize (nfield);

  alloc = new real [(size_t) nfield * ntot];
  for (i = 0; i < nfield; i++) {
    udat[i] = alloc + i * ntot;
    u[i]    = new Field (b[i], udat[i], nz, E, field[i]);
  }

  VERBOSE cout << "done" << endl;
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

  cout << "   Solution fields         : " << field              << endl;

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


void Domain::restart ()
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
    transform (FORWARD);
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

  Femlib::synchronize();

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

  Femlib::synchronize();
  this -> transform (INVERSE);
  Femlib::synchronize();

  output << *this;

  Femlib::synchronize();
  this -> transform (FORWARD);
  Femlib::synchronize();

  ROOTONLY output.close();
}


void Domain::transform (const integer sign)
// ---------------------------------------------------------------------------
// Fourier transform all Fields according to sign.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer N = this -> nField ();
  
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
  int               i;
  const int         N = D.u.size();
  vector<AuxField*> field (N);

  for (i = 0; i < N; i++) field[i] = D.u[i];

  writeField (strm, D.name, D.step, D.time, field);

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
  if (nfields != strlen (D.field)) {
    sprintf (err, "file: %1d fields, Domain: %1d", nfields, strlen(D.field));
    message (routine, err, ERROR);
  }
  for (i = 0; i < nfields; i++) 
    if (!strchr (D.field, fields[i])) {
      sprintf (err,"field %c not present in Domain (%s)",fields[i],D.field);
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
      if (fields[j] == D.field[i]) {
	strm >>  *D.u[i];
	if (swap) D.u[i] -> reverse();
	break;
      }
  }
    
  ROOTONLY if (strm.bad())
    message (routine, "failed reading field file", ERROR);
    
  return strm;
}
