///////////////////////////////////////////////////////////////////////////////
// dual2sem.C: convert a dual 3D file to a semtex/nekton 3D file.
//
// Copyright (c) 2000 <--> $Date$, Hugh Blackburn
//
// USAGE
// -----
// dual2sem [options] [file]
// options:
// -h       ... print this message.
// -f <num> ... the mode number of the 3D mode is <num>. [Default: 1]
// -i       ... invert (sem --> dual).
// 
// If file is not present, read from standard input.  Write to
// standard output.
//
// A dual field file describes a 3D data set, but carries only two
// Fourier modes, the zeroth (just real data) and one more mode
// (complex), with a nominated mode number.  So only three planes of
// data are held, and in addition the single complex mode is held in
// complex form (its first plane is real, its second is imaginary).
//
// Dual2sem projects a dual field file onto a nekton type field file,
// with the minimal number of modes required to contain it (subject to
// the 2, 3, 5 prime factor requirements for FFT).  Dual2sem cannot
// determine from the dual field file the number of the fundamental 3D
// mode (the header for the file is standard semtex format), so it
// needs a command-line argument to set it.  Likewise the mode number
// is needed in the inverse operation (in order to select the
// fundamental mode number from the semtex/nekton file).
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include "sem.h"

static char prog[] = "dual2sem";

static char* hdr_fmt[] = { 
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

typedef struct hdr_data {
  char   session[StrMax];
  char   created[StrMax];
  int_t  nr, ns, nz, nel;
  int_t  step;
  real_t time;
  real_t timestep;
  real_t kinvis;
  real_t beta;
  char   fields[StrMax];
  char   format[StrMax];
} hdr_info;

static void getargs   (int_t, char**, int_t&, int_t&, istream*&);
static void roundup   (const int_t, int_t&);
static void getheader (istream&, hdr_info&);
static void allocate  (hdr_info&, const int_t, vector<real_t*>&);
static void readdata  (hdr_info&, istream&, const int_t, vector<real_t*>&);
static void packdata  (hdr_info&, const int_t, const int_t, vector<real_t*>&,
		       const int_t);
static void transform (hdr_info&, const int_t, vector<real_t*>&, const int_t);
static void writedata (hdr_info&, ostream&, const int_t, vector<real_t*>&);
static bool doswap    (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  hdr_info        header;
  int_t           dir = FORWARD, mode = 1, nz;
  istream*        input;
  vector<real_t*> u;

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, dir, mode, input);

  if (dir == FORWARD) {		// -- Project a dual field file.

    roundup   (mode, nz);
    getheader (*input, header);
    allocate  (header, nz, u);
    readdata  (header, *input, 3, u);
    packdata  (header, mode, nz, u, INVERSE);
    transform (header, nz, u, INVERSE);
    writedata (header, cout, nz, u);

  } else {			// -- Restrict a sem field file.

    getheader (*input, header);
    allocate  (header, header.nz, u);
    readdata  (header, *input, header.nz, u);
    transform (header, header.nz, u, FORWARD);
    packdata  (header, mode, header.nz, u, FORWARD);
    writedata (header, cout, 3, u);

  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     int_t&    dir  ,
		     int_t&    mode ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: dual2sem [options] [file]\n"
    "options:\n"
    "-h       ... print this message\n"
    "-f <num> ... the mode number of the 3D mode is <num>. [Default: 1]\n"
    "-i       ... invert (sem --> dual).\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'i':
      dir = INVERSE;
      break;
    case 'f':
      if (*++argv[0]) mode = atoi (*argv);
      else {mode = atoi (*++argv); argc--;}
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> bad()) message (prog, "unable to open input file", ERROR);
  } else input = &cin;

  if (mode < 1) {
    cerr << prog << ": mode number must be a positive integer" << endl;
    exit (EXIT_FAILURE);
  }
}


static void roundup (const int_t mode, 
		     int_t&      nz  )
// ---------------------------------------------------------------------------
// Calculate nz from mode number, then increment until it suits FFT.
// ---------------------------------------------------------------------------
{
  int_t n, ip, iq, ir, ipqr2;

  nz = 2 * (mode + 1);

  do {
    n = nz;
    Femlib::primes235 (n, ip, iq, ir, ipqr2);
    nz += 2;
  } while (n == 0);
  nz = n;
}


static void getheader (istream&  file  ,
		       hdr_info& header)
// ---------------------------------------------------------------------------
// Load data structure from file header info.
// ---------------------------------------------------------------------------
{
  char  buf[StrMax];
  int_t i, j; 

  file.get (header.session, 25); file.getline (buf, StrMax);
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);

  file.get (header.created, 25); file.ignore (StrMax, '\n');

  file >> header.nr >> header.ns >> header.nz >> header.nel;
  file.ignore (StrMax, '\n');

  file >> header.step;     file.ignore (StrMax, '\n');
  file >> header.time;     file.ignore (StrMax, '\n');
  file >> header.timestep; file.ignore (StrMax, '\n');
  file >> header.kinvis;   file.ignore (StrMax, '\n');
  file >> header.beta;     file.ignore (StrMax, '\n');

  file.get (buf, 25);
  for (i = 0, j = 0; i < 25; i++)
    if (isalpha (buf[i])) header.fields[j++] = buf[i];
  file.ignore (StrMax, '\n');
  header.fields[j] = '\0';

  file.get (header.format, 25); file.getline (buf, StrMax);

  if (!strstr (header.format, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (header.format, "-endia"))
    message (prog, "input field file in unknown binary format", WARNING);
}


static void allocate (hdr_info&        header,
		      const int_t      nz    ,
		      vector<real_t*>& u     )
// ---------------------------------------------------------------------------
// Allocate enough storage to hold all the data fields (sem format).
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  const int_t ntot      = planesize * nz;
  int_t       i;

  u.resize (nfield);
  
  for (i = 0; i < nfield; i++) {
    u[i] = new real_t [ntot];
    Veclib::zero (ntot, u[i], 1);
  }
}


static void readdata (hdr_info&        header,
		      istream&         file  ,
		      const int_t      nplane,
		      vector<real_t*>& u     )
// ---------------------------------------------------------------------------
// Binary read of data area, followed by byte-swapping if required.
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  const int_t ntot      = planesize * nplane;
  const bool  swab      = doswap (header.format);
  int_t       i, j;
  
  for (i = 0; i < nfield; i++) {
    for (j = 0; j < nplane; j++) {
      file.read (reinterpret_cast<char*>(u[i]+j*planesize),
		 ntotelmt*sizeof(real_t));
    }
    if (swab) Veclib::brev (ntot, u[i], 1, u[i], 1);
  }
}


static void packdata (hdr_info&        header,
		      const int_t      mode  ,
		      const int_t      nz    ,
		      vector<real_t*>& u     ,
		      const int_t      dir   )
// ---------------------------------------------------------------------------
// Data in u are considered to be in Fourier-transformed state.  Move
// planes around to place them in required spots.  Zero other data.
// ---------------------------------------------------------------------------
{
  const int_t    nfield    = strlen (header.fields);
  const int_t    ntotelmt  = header.nr * header.ns * header.nel;
  const int_t    planesize = ntotelmt + (ntotelmt % 2);
  const int_t    nblock    = 2 * planesize;
  vector<real_t> tmp (2 * planesize);
  int_t          i;

  if (dir == FORWARD)		// -- Pack down to 3 data planes.
    for (i = 0; i < nfield; i++) {
      Veclib::copy (nblock, u[i] + mode * nblock, 1, &tmp[0], 1);
      Veclib::zero ((nz - 1) * planesize, u[i] + planesize, 1);
      Veclib::copy (nblock, &tmp[0], 1, u[i] + planesize, 1);
    }
  else				// -- Reposition nominated mode.
    for (i = 0; i < nfield; i++) {
      Veclib::copy (nblock, u[i] + planesize, 1, &tmp[0], 1);
      Veclib::zero ((nz - 1) * planesize, u[i] + planesize, 1);
      Veclib::copy (nblock, &tmp[0], 1, u[i] + mode * nblock, 1);
    }
}


static void transform (hdr_info&      header,
		       const int_t      nz    , 
		       vector<real_t*>& u     , 
		       const int_t      dir   )
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  int_t       i;

  for (i = 0; i < nfield; i++) Femlib::DFTr (u[i], nz, planesize, dir);
}


static void writedata (hdr_info&        header,
		       ostream&         file  ,
		       const int_t      nz    ,
		       vector<real_t*>& u     )
// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
{
  const int_t nfield    = strlen (header.fields);
  const int_t ntotelmt  = header.nr * header.ns * header.nel;
  const int_t planesize = ntotelmt + (ntotelmt % 2);
  char        buf[StrMax], tmp[StrMax];
  int_t       i, j;

  sprintf (buf, hdr_fmt[0], header.session);
  file << buf;
  sprintf (buf, hdr_fmt[1], header.created);
  file << buf;

  sprintf (tmp, "%1d %1d %1d %1d", header.nr, header.ns, nz, header.nel);
  sprintf (buf, hdr_fmt[2], tmp);
  file << buf;

  sprintf (buf, hdr_fmt[3], header.step);     file << buf;
  sprintf (buf, hdr_fmt[4], header.time);     file << buf;
  sprintf (buf, hdr_fmt[5], header.timestep); file << buf;
  sprintf (buf, hdr_fmt[6], header.kinvis);   file << buf;
  sprintf (buf, hdr_fmt[7], header.beta);     file << buf;

  sprintf (buf, hdr_fmt[8], header.fields);
  file << buf;  

  sprintf (tmp, "binary ");
  Veclib::describeFormat (tmp + strlen (tmp));
  sprintf (buf, hdr_fmt[9], tmp);
  file << buf;

  for (i = 0; i < nfield; i++) 
    for (j = 0; j < nz; j++)
      file.write(reinterpret_cast<char*>(u[i]+j*planesize),
		 ntotelmt*sizeof(real_t));
}


static bool doswap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping of input is required to make sense of input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "-endia"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}
