///////////////////////////////////////////////////////////////////////////////
// SYNOPSIS
// refreq.C: compute new body motion file, for cosine oscillation.
//
// USAGE
// refreq [options] session
//   options:
//   -k <num> ... factor original frequencies by num, or ...
//   -s <num> ... sets new frequencies to num.
//   -z       ... computes new phase angle for time zero.
//
// FILES
// Input files required are session.bdy, session.rst.
//
// Body motion file format: {entry not allowed in file}/(treated as comment).
// -- BOF --
// Title string.                          {Just a header, not used.}
// Another string or blank line.          (Cosmetic padding.)
// x-axis cosine     ampl freq phase      (3 real values in ASCII.)
// y-axis feedback   mass freq zeta       (3 real values in ASCII.)
// x-state           pos  vel             (2 real values, default to zero.)
// x-state           pos  vel             (Repetitions allowed, last used.)
// y-state           pos  vel             (X & Y states can be intermixed.)
// -- EOF --
//
// Restart file format: standard SEM-compatible format.  Header:
// -- BOF --
// chan1                     Session
// Wed Oct 02 22:49:27 1996  Created
// 9 9 1 4                   Nr, Ns, Nz, Elements
// 200                       Step
// 2                         Time
// 0.01                      Time step
// 0.5                       Kinvis
// 1                         Beta
// uvp                       Fields written
// binary IEEE little-endian Format
//
// Output is a new body motion file, written to standard output stream.
//
// NOTES
// 1.  The only information used from the restart file is the Time.
// 2.  Frequencies in body motion file are frequencies, not radiancies.
// 3.  No action takes place on lines that do not include "axis cosine".
// 4.  In the conversion to a new frequency, the original frequency and
//     phase angles are assumed to be constant with time.
//
///////////////////////////////////////////////////////////////////////////////

static char
RCSid[] = "$Id$";

#include "Sem.h"

static char prog[] = "refreq";
static void getargs (int, char**, char*&, double&, double&, int&);


int main (int argc, char** argv)
// ---------------------------------------------------------------------------
// Driver routine.
// ---------------------------------------------------------------------------
{
  char     err[StrMax], buf[StrMax];
  ifstream file;
  char*    session = 0;
  double   factor  = 0.0, freq = 0.0, time;
  int      i, zero = 0;

  getargs (argc, argv, session, factor, freq, zero);
  Femlib::prep();

  // -- Open restart file and extract the simulation time.

  strcat (strcpy (buf, session), ".rst");
  file.open (buf);

  if (!file) {
    sprintf (err, "couldn't open restart file %s", buf);
    message (prog, err, ERROR);
  }

  for (i = 0; i < 4; i++) file.getline (buf, StrMax);
  file >> time;
  file.getline (buf, StrMax);
  
  if (!strstr (buf, "Time")) {
    sprintf (err, "couldn't extract Time from 5th line of restart file");
    message (prog, err, ERROR);
  }

  file.close ();

  // -- Open & process body motion file.

  strcat (strcpy (buf, session), ".bdy");
  file.open (buf);

  if (!file) {
    sprintf (err, "couldn't open body motion file %s", buf);
    message (prog, err, ERROR);
  }

  while (file.getline (buf, StrMax)) {

    if (strstr (buf, "cosine")) { // -- This is what we were waiting for...
      char       axis[StrMax], motion[StrMax], amplitude[StrMax];
      char       frequency[StrMax], phaseangle[StrMax];
      istrstream line (buf, strlen (buf));

      line >> axis >> motion >> amplitude >> frequency >> phaseangle;
      if (line.bad ()) {
	sprintf (err, "couldn't parse all values from: %s", buf);
	message (prog, err, ERROR);
      } 

      double oldfreq  = Femlib::value (frequency),  newfreq,
	     oldphase = Femlib::value (phaseangle), newphase;

      if (factor != 0.0)
	newfreq  = factor * oldfreq;
      else if (freq != 0.0)
	newfreq  = freq;

      if (zero)
	newphase = 2.0 * M_PI *  oldfreq            * time + oldphase;
      else
	newphase = 2.0 * M_PI * (oldfreq - newfreq) * time + oldphase;

      newphase = fmod (newphase, 2.0* M_PI);

      cout << axis << " " << motion   << " " << amplitude << " "
	<< newfreq << " " << newphase << endl;
	
    } else
      cout << buf << endl;
  }

  file.close ();

  return EXIT_SUCCESS;
}


static void getargs (int     argc   ,
		     char**  argv   ,
		     char*&  session,
		     double& factor ,
		     double& freq   , 
		     int&    zero   )
// ---------------------------------------------------------------------------
// Parse command-line args.
// ---------------------------------------------------------------------------
{
  char err[StrMax], buf[StrMax];
  char usage[] =
    "Usage: %s [options] session\n"
    "options are:\n"
    "-h       ... print this message\n"
    "-k <num> ... factor original frequencies by num, or ...\n"
    "-s <num> ... sets new frequencies to num\n"
    "-z       ... computes new phase angle at time zero\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cerr << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'k':
      if   (*++argv[0])           factor = atof (*argv);
      else              { --argc; factor = atof (*++argv); }
      break;
    case 's':
      if   (*++argv[0])           freq   = atof (*argv);
      else              { --argc; freq   = atof (*++argv); }
      break;
    case 'z':
      zero = 1;
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }
   
  if (argc != 1) {
    sprintf (err, usage, prog);
    cerr << buf;
    exit (EXIT_FAILURE);
  } else
    session = *argv; 

  if (freq == 0.0 && factor == 0.0 && zero == 0)
    message (prog, "must either set factor or new frequency", ERROR);
}

