///////////////////////////////////////////////////////////////////////////////
// statistics.C: routines for statistical analysis of AuxFields.
//
// At present, this is limited to running averages, but could be
// modified for computing Reynolds stresses, etc.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include <Sem.h>
#include <time.h>

Statistics::Statistics (Domain&            D    ,
			vector<AuxField*>& extra) : 
                        name (D.name),
			base (D)
// ---------------------------------------------------------------------------
// Store averages for all Domain Fields, and any extra AuxFields
// supplied.
//
// Try to initialize from file session.avg, failing that set all
// buffers to zero.  Number of fields in file should be same as
// Statistics::avg buffer.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer NF = base.u.getSize();
  const integer NE = extra.getSize();
  const integer NT = NF + NE;

  ROOTONLY cout << "-- Initializing averaging  : ";  

  // -- Set pointers, allocate storage.

  src.setSize (NT);
  avg.setSize (NT);
  
  for (i = 0; i < NF; i++) src[     i] = (AuxField*) base.u[i];
  for (i = 0; i < NE; i++) src[NF + i] = extra[i];

  for (i = 0; i < NT; i++) avg[i] = new AuxField (base.Esys, src[i] -> name());

  // -- Initialize averages, either from file or zero.
  //    This is much the same as Domain input routine.

  char     s[StrMax];
  ifstream file (strcat (strcpy (s, name), ".avg"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << s;
      cout.flush();
    }
    file >> *this;
    file.close();
    for (i = 0; i < NT; i++) avg[i] -> transform (+1);
  
  } else {			// -- No file, set to zero.
    ROOTONLY cout << "set to zero";
    for (i = 0; i < NT; i++) *avg[i] = 0.0;
    navg = 0;
  }

  ROOTONLY cout << endl;
}


void Statistics::update (AuxField*** work)
// ---------------------------------------------------------------------------
// Update running averages.  Zeroth time level of work is available as
// workspace, but not used yet.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer N = avg.getSize();

  for (i = 0; i < N; i++) {
    *avg[i] *= (real)  navg;
    *avg[i] += *src[i];
    *avg[i] /= (real) (navg + 1);
  }
  navg++;
}


void Statistics::dump ()
// ---------------------------------------------------------------------------
// Similar to Domain::dump
// ---------------------------------------------------------------------------
{
  const integer step     = base.step;
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;

  integer   i;
  const int N = avg.getSize();
  ofstream  output;

  ROOTONLY {
    const char    routine[] = "Statistics::dump";
    const integer verbose   = (integer) Femlib::value ("VERBOSE");
    const integer chkpoint  = (integer) Femlib::value ("CHKPOINT");
    char          dumpfl[StrMax], backup[StrMax], command[StrMax];

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, name), ".avg");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, name), ".ave");
	if (!initial) {
	  strcat  (strcpy (backup, name), ".ave.bak");
	  sprintf (command, "mv ./%s ./%s", dumpfl, backup);
	  system  (command);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, name), ".avg");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }

    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }

  for (i = 0; i < N; i++) avg[i] -> transform (-1);
  output << *this;
  for (i = 0; i < N; i++) avg[i] -> transform (+1);

  ROOTONLY output.close();
}


ostream& operator << (ostream&    strm,
		      Statistics& src )
// ---------------------------------------------------------------------------
// Output Statistics class to file.  Like similar Domain routine.
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

  const char    routine[] = "ostream<<Statistics";
  const integer N    = src.avg.getSize();
  const real    time = src.base.time;
  integer       k;
  char          s1[StrMax], s2[StrMax];
  time_t        tp (::time (0));

  ROOTONLY {
    sprintf (s1, hdr_fmt[0], src.name);
    strm << s1;

    strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
    sprintf  (s1, hdr_fmt[1], s2);
    strm << s1;

    src.avg[0] -> describe (s2);
    sprintf (s1, hdr_fmt[2], s2);
    strm << s1;

    sprintf (s1, hdr_fmt[3], src.navg);
    strm << s1;

    sprintf (s1, hdr_fmt[4], time);
    strm << s1;

    sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
    strm << s1;

    sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
    strm << s1;
  
    sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
    strm << s1;

    for (k = 0; k < N; k++) s2[k] = src.avg[k] -> name();
    s2[k] = '\0';
    sprintf (s1, hdr_fmt[8], s2);
    strm << s1;

    sprintf (s2, "binary ");
    Veclib::describeFormat (s2 + strlen (s2));
    sprintf (s1, hdr_fmt[9], s2);
    strm << s1;
  }

  for (k = 0; k < N; k++) strm << *src.avg[k];

  ROOTONLY {
    if (!strm) message (routine, "failed writing average file", ERROR);
    strm << flush;
  }

  return strm;
}


istream& operator >> (istream&    strm,
		      Statistics& tgt )
// ---------------------------------------------------------------------------
// Input Statistics class from file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  const char routine[] = "strm>>Statistics";
  integer    i, j, np, nz, nel, ntot, nfields;
  integer    npchk,  nzchk, nelchk, swap = 0;
  char       s[StrMax], f[StrMax], err[StrMax], fields[StrMax];

  if (strm.getline(s, StrMax).eof()) return strm;
  
  strm.getline (s, StrMax) . getline (s, StrMax);

  tgt.avg[0] -> describe (f);
  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline (s,StrMax);
  istrstream (s, strlen (s)) >> tgt.navg;
    
  strm.getline (s, StrMax) . getline (s, StrMax);
  strm.getline (s, StrMax) . getline (s, StrMax) . getline (s, StrMax);
    
  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';
  if (nfields != tgt.avg.getSize()) {
    sprintf (err, "strm: %1d fields, avg: %1d", nfields, tgt.avg.getSize());
    message (routine, err, ERROR);
  }
  for (i = 0; i < nfields; i++) 
    if (!strchr (fields, tgt.avg[i] -> name())) {
      sprintf (err, "field %c not present in avg", fields[i]);
      message (routine, err, ERROR);
    }

  strm.getline (s, StrMax);
  Veclib::describeFormat (f);

  if (!strstr (s, "binary"))
    message (routine, "input field strm not in binary format", ERROR);
  
  if (!strstr (s, "endian"))
    message (routine, "input field strm in unknown binary format", WARNING);
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
      if (tgt.avg[i] -> name() == fields[j]) break;
    strm >> *tgt.avg[i];
    if (swap) tgt.avg[i] -> reverse();
  }
  
  ROOTONLY if (strm.bad())
    message (routine, "failed reading average file", ERROR);

  return strm;
}



