///////////////////////////////////////////////////////////////////////////////
// domain.C:  implement Domain class.
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include "Fem.h"


Domain::Domain (Mesh& M, const char* session, const int& np)
// ---------------------------------------------------------------------------
// Construct a new Domain with a single Field, using geometry from M.
//
// On return, all information in Element and Boundary lists is set for u[0].
// ---------------------------------------------------------------------------
{
  cout << "-- " << M.nEl () << " elements, ";
  cout << np << "x" << np << endl;
  
  domain_name = new char[StrMax];
  strcpy (domain_name, session);

  domain_step = 0;
  domain_time = 0.0;
  nfield      = 1;

  u    = new SystemField* [nfield];
  u[0] = new SystemField (M, np);
}


void Domain::addField (SystemField* F)
// ---------------------------------------------------------------------------
// Add a new Field pointer to Domain's array.
// ---------------------------------------------------------------------------
{
  SystemField** X = new SystemField* [nfield + 1];

  for (int i (0); i < nfield; i++) X[i] = u[i];
  X[nfield] = F;

  delete [] u;
  u = X;

  nfield++;
}


void Domain::restart ()
// ---------------------------------------------------------------------------
// If a restart file "name".rst can be found, use it to load
// velocity fields.  Otherwise initialize velocity fields to zero.
// ---------------------------------------------------------------------------
{
  char  restartfile[StrMax];

  strcat (strcpy (restartfile, domain_name), ".rst");

  ifstream file (restartfile);

  if (file) {
    cout << "-- Restarting from file:  " << restartfile;
    file >> *this;
  } else {
    cout << "-- Initializing solution at zero";
    for (int i(0); i < nField (); i++)
      *u[i] = 0.0;
  }

  cout << endl;

  real t  = time   ();
  real dt = dparam ("DELTAT");
  int  ns = iparam ("N_STEP");
  int  nt = iparam ("N_TIME");

  setDparam ("t", t);

  cout << "   Start time:            " << t           << endl;
  cout << "   Time step:             " << dt          << endl;
  cout << "   Number of steps:       " << ns          << endl;
  cout << "   End time:              " << t + ns * dt << endl;
  cout << "   Integration order:     " << nt          << endl;
}


void  Domain::dump (ofstream& output)
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
// ---------------------------------------------------------------------------
{
  int periodic = !(domain_step %  iparam ("IO_FLD"));
  int final    =   domain_step == iparam ("N_STEP");

  if ( ! (periodic || final) ) return;

  char  routine[] = "Domain::dump";
  int   verbose   = option ("VERBOSE");
  int   chkpoint  = option ("CHKPOINT");

  if (chkpoint) {
    char s[StrMax], b[StrMax], c[StrMax];

    output.close ();

    if (final)
      strcat (strcpy (s, domain_name), ".fld");
    else {
      strcat (strcpy (s, domain_name), ".chk");
      strcat (strcpy (b, domain_name), ".chk.bak");
    
      sprintf (c, "mv ./%s ./%s", s, b);
      system  (c);
    }

    output.open (s);
    if (!output) message (routine, "can't open checkpoint file", ERROR);
  }

  if (verbose) message (routine, ": writing field dump", REMARK);

  output << *this;

}


#include <time.h>


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

  char      routine[] = "ostream << Domain";
  char      s1[StrMax], s2[StrMax];
  time_t    tp   (::time (0));

  sprintf (s1, hdr_fmt[0], D.name());
  strm << s1;

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  strm << s1;

  D.u[0] -> describe (s2);
  sprintf (s1, hdr_fmt[2], s2);
  strm << s1;

  sprintf (s1, hdr_fmt[3], D.step());
  strm << s1;

  sprintf (s1, hdr_fmt[4], D.time());
  strm << s1;

  sprintf (s1, hdr_fmt[5], dparam ("DELTAT"));
  strm << s1;

  sprintf (s1, hdr_fmt[6], dparam ("KINVIS"));
  strm << s1;

  sprintf (s1, hdr_fmt[7], dparam ("BETA"));
  strm << s1;

  int k;
  for (k = 0; k < D.nField(); k++) s2[k] = D.u[k] -> getName ();
  s2[k] = '\0';
  sprintf (s1, hdr_fmt[8], s2);
  strm << s1;

  sprintf (s2, "binary, ");
  Veclib::describeFormat (s2 + strlen (s2));
  sprintf (s1, hdr_fmt[9], s2);
  strm << s1;

  for (int n (0); n < D.nField (); n++) strm << *D.u[n];

  if (!strm) message (routine, "failed writing field file", ERROR);
  strm << flush;

  return strm;
}


istream& operator >> (istream& strm, Domain& D)
// ---------------------------------------------------------------------------
// Input all Domain field variables from prism-compatible istream.
//
// Only binary storage format is allowed.  Check if conversion to native
// format (IEEE little/big-endian) is required.
// ---------------------------------------------------------------------------
{
  char  routine[] = "strm>>Domain";
  int   np, nz, nel, ntot, nfields;
  int   npchk,  nzchk, nelchk, swap = 0;
  char  s[StrMax], f[StrMax];

  if (strm.getline(s, StrMax).eof()) return strm;

  strm . getline (s, StrMax) . getline (s, StrMax);

  D.u[0] -> describe (f);
  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != D.u[0] -> nTot ())
    message (routine, "declared sizes mismatch",     ERROR);

  strm . getline (s, StrMax) . getline(s, StrMax);

  istrstream (s, strlen (s)) >> D.time ();
  setDparam  ("t", D.time ());
  
  strm.getline(s,StrMax).getline(s,StrMax).getline(s,StrMax).getline(s,StrMax);
  
  nfields = 0;
  while (isalpha (s[nfields])) nfields++;
  if (nfields > D.nField ())
    message (routine, "number of fields in file > number in domain", ERROR);
  nfields = 0;
  while (isalpha (s[nfields])) {
    D.u[nfields] -> setName (s[nfields]);
    nfields++;
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

  for (int n = 0; n < nfields; n++) {
    strm >> *D.u[n];
    if (swap) D.u[n] -> reverse ();
  }

  if (strm.bad ()) message (routine, "failed reading field file", ERROR);
    
  return strm;
}
