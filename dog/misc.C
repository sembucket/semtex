///////////////////////////////////////////////////////////////////////////////
// misc.C: miscellaneous routines for I/O, memory management, service
// routines that don't fit class structures.
//
// Copyright (C) 1994, 2001 Hugh BLackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>
#include <time.h>


ostream& printVector (ostream&      strm,
		      const char*   fmt , 
		      const integer ntot,
		                    ... )
// ---------------------------------------------------------------------------
// Print up a variable number of numeric vectors on strm, in columns.
//
// The format specifier is gives the number and type of the vectors, with
// type specified by the first character.  The vectors must all be of the
// same type & length.
//
// Allowed types are: int ("i"), real ("r").
// Examples: four integer vectors ==> fmt is "iiii".  Two reals ==> "rr".
// 
// Vectors are printed in a fixed field width of 15, regardless of type.
// ---------------------------------------------------------------------------
{
  char    routine[] = "printVector";
  integer nvect;
  va_list ap;

  nvect = strlen (fmt);
  if (! nvect   ) message (routine, "empty format string",   ERROR  );
  if (nvect > 10) message (routine, "more than 10 vectors?", WARNING);
  if (ntot  <  0) message (routine, "ntot < 0",              ERROR  );

  switch (fmt[0]) {

  case 'i': {
    integer** u = new integer* [nvect];
    va_start (ap, ntot);
    for (integer k = 0; k < nvect; k++) u[k] = va_arg (ap, integer*);
    va_end (ap);
    for (register integer l = 0; l < ntot; l++) {
      for (register integer j = 0; j < nvect; j++)
	strm << setw(15) << u[j][l];
      strm << endl;
    }
    delete [] u;
    break;
  }
  case 'r': {
    real** u = new real* [nvect];
    va_start (ap, ntot);
    for (integer k = 0; k < nvect; k++) u[k] = va_arg (ap, real*);
    va_end (ap);
    for (register integer l = 0; l < ntot; l++) {
      for (register integer j = 0; j < nvect; j++)
	strm << setw(15) << u[j][l];
      strm << endl;
    }
    delete [] u;
    break;
  }
  default:
    message (routine, fmt, ERROR);
    break;
  }

  if (!strm) message (routine, "output failed", ERROR);

  return strm;
}


char* upperCase (char *s)
// ---------------------------------------------------------------------------
// Uppercase characters in string.
// ---------------------------------------------------------------------------
{
  char *z(s); while (*z = toupper (*z)) z++; return s;
}


void writeField (ofstream&          file   ,
		 const char*        session,
		 const int          runstep,
		 const real         runtime,
		 vector<AuxField*>& field  )
// ---------------------------------------------------------------------------
// Write fields out to an opened file, binary prism format.  Output is
// only done by the root processor.
// ---------------------------------------------------------------------------
{
  const char routine [] = "writeField";
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

  char      s1[StrMax], s2[StrMax];
  time_t    tp (time (0));
  int       i;
  const int N = field.getSize();

  if (N < 1) return;

  ROOTONLY {
    sprintf (s1, hdr_fmt[0], session);
    file << s1;

    strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
    sprintf  (s1, hdr_fmt[1], s2);
    file << s1;

    field[0] -> describe (s2);
    sprintf (s1, hdr_fmt[2], s2);
    file << s1;

    sprintf (s1, hdr_fmt[3], runstep);
    file << s1;

    sprintf (s1, hdr_fmt[4], runtime);
    file << s1;

    sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
    file << s1;

    sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
    file << s1;

    sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
    file << s1;

    for (i = 0; i < N; i++) s2[i] = field[i] -> name();
    s2[i] = '\0';
    sprintf (s1, hdr_fmt[8], s2);
    file << s1;

    sprintf (s2, "binary ");
    Veclib::describeFormat (s2 + strlen (s2));
    sprintf (s1, hdr_fmt[9], s2);
    file << s1;
  }

  for (i = 0; i < N; i++) file << *field[i];

  ROOTONLY {
    if (!file) message (routine, "failed writing field file", ERROR);
    file << flush;
  }
}


Header::Header ()
// ---------------------------------------------------------------------------
// Allocate storage for Header.
// ---------------------------------------------------------------------------
{
  sess = new char [StrMax];
  sesd = new char [StrMax];
  flds = new char [StrMax];
  frmt = new char [StrMax];

  sess[0] = sesd[0] = flds[0] = '\0';
  sprintf (frmt, "binary "); Veclib::describeFormat (frmt + strlen (frmt));
  nr = ns = nz = nel = step = 0;
  time = dt = visc = beta = 0.0;
}


ifstream& operator >> (ifstream& file,
		       Header&   hdr )
// ---------------------------------------------------------------------------
// Insert data into Header struct from file.
// ---------------------------------------------------------------------------
{
  char routine[] = "operator: ifstream >> Header";
  char s[StrMax];

  if (file.get(hdr.sess, 25).eof()) return file; file.getline(s, StrMax);
  file.get(hdr.sesd, 25);                        file.getline(s, StrMax);
  file >> hdr.nr >> hdr.ns >> hdr.nz >> hdr.nel; file.getline(s, StrMax);
  file >> hdr.step;                              file.getline(s, StrMax);
  file >> hdr.time;                              file.getline(s, StrMax);
  file >> hdr.dt;                                file.getline(s, StrMax);
  file >> hdr.visc;                              file.getline(s, StrMax);
  file >> hdr.beta;                              file.getline(s, StrMax);
  file >> hdr.flds;                              file.getline(s, StrMax);
  file.get(hdr.frmt, 25);                        file.getline(s, StrMax);

  if (!file) message (routine, "failed reading header information", ERROR);
  return file;
}


ofstream& operator << (ofstream& file,
		       Header&   hdr )
// ---------------------------------------------------------------------------
// Put data from Header struct onto file. Use current time info.
// ---------------------------------------------------------------------------
{
  const char routine [] = "operator: ofstream << Header";
  const char *hdr_fmt[] = { 
    "%-25s "                "Session\n",
    "%-25s "                "Created\n",
    "%-5d %-5d %-5d %-5d "  "Nr, Ns, Nz, Elements\n",
    "%-25d "                "Step\n",
    "%-25.6g "              "Time\n",
    "%-25.6g "              "Time step\n",
    "%-25.6g "              "Kinvis\n",
    "%-25.6g "              "Beta\n",
    "%-25s "                "Fields written\n",
    "%-25s "                "Format\n"
  };

  char   s1[StrMax], s2[StrMax];
  time_t tp (time (0));

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));

  sprintf  (s1, hdr_fmt[0], hdr.sess);                        file << s1;
  sprintf  (s1, hdr_fmt[1], s2);                              file << s1;
  sprintf  (s1, hdr_fmt[2], hdr.nr, hdr.ns, hdr.nz, hdr.nel); file << s1;
  sprintf  (s1, hdr_fmt[3], hdr.step);                        file << s1;
  sprintf  (s1, hdr_fmt[4], Femlib::value ("t"));             file << s1;
  sprintf  (s1, hdr_fmt[5], Femlib::value ("D_T"));           file << s1;
  sprintf  (s1, hdr_fmt[6], Femlib::value ("KINVIS"));        file << s1;
  sprintf  (s1, hdr_fmt[7], Femlib::value ("BETA"));          file << s1;
  sprintf  (s1, hdr_fmt[8], hdr.flds);                        file << s1;
  sprintf  (s1, hdr_fmt[9], hdr.frmt);                        file << s1;

  if (!file) message (routine, "failed writing field file header", ERROR);
  file << flush;

  return file;
}
