///////////////////////////////////////////////////////////////////////////////
// statistics.C: routines for statistical analysis of AuxFields.
//
// At present, this is limited to running averages of primitive
// variables and product terms for Reynolds stresses.
//
// Naming scheme for "Reynolds stresses":
//
//                      / uu uv uw \     /  A  B  D \
//                      | .  vv vw |  =  |  .  C  E |
//                      \ .  .  ww /     \  .  .  F /
//
// Reynolds stresses are selected for computation if AVERAGE > 1.  In this
// case it is assumed that there will be at least DIM fields to be multiplied
// together as indicated.  NB: what is computed are the running average of the 
// products uu, uv, etc, which are NOT the actual Reynolds stresses: they need
// to have the products of the mean values UU, UV etc subtracted, assumed
// to occur in postprocessing.
//
// Running averages are kept semi-Fourier, Reynolds stress products in
// physical space to minimize number of transforms needed.
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
//
// NR = Number of Reynolds stress averaging buffers, set if AVERAGE > 1.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer ND = Geometry::nDim();
  const integer NF = base.u.getSize();
  const integer NE = extra.getSize();
  const integer NR = ((integer) Femlib::value ("AVERAGE") > 1) ? 
                            ((ND + 1) * ND) >> 1 : 0;
  const integer NT = NF + NE + NR;

  ROOTONLY cout << "-- Initializing averaging  : ";  

  // -- Set pointers, allocate storage.

  src.setSize (NF + NE);	// -- Straight running average of these.
  avg.setSize (NT);		// -- Additional are computed from src.
  
  for (i = 0; i < NF; i++) src[     i] = (AuxField*) base.u[i];
  for (i = 0; i < NE; i++) src[NF + i] = extra[i];

  for (i = 0; i < NF + NE; i++)
    avg[i] = new AuxField (base.Esys, src[i] -> name());
  for (i = 0; i < NT - NF - NE; i++)
    avg[i + NF + NE] = new AuxField (base.Esys, 'A' + i);

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
    for (i = 0; i < NT - NR; i++) avg[i] -> transform (+1);
  
  } else {			// -- No file, set to zero.
    ROOTONLY cout << "set to zero";
    for (i = 0; i < NT; i++) *avg[i] = 0.0;
    navg = 0;
  }

  ROOTONLY cout << endl;
}


void Statistics::update (AuxField*** work)
// ---------------------------------------------------------------------------
// Update running averages, using zeroth time level of work as
// workspace.  Reynolds stress terms are calculated without
// dealiasing, and are held in physical space.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer NT = avg.getSize();
  const integer ND = Geometry::nDim();
  const integer NR = ((integer) Femlib::value ("AVERAGE") > 1) ? 
                            ((ND + 1) * ND) >> 1 : 0;
  const integer NA = NT - NR;

  if (NR) {
    
    // -- Running averages and Reynolds stresses.

    for (i = 0; i < ND; i++) {
      *work[i][0] = *src[i];
       work[i][0] -> transform (-1);
    }

    for (i = 0; i < NA; i++) {
      *avg[i] *= (real) navg;
      *avg[i] += *src[i];
    }

    avg[NA + 0] -> timesPlus (*work[0][0], *work[0][0]);
    avg[NA + 1] -> timesPlus (*work[0][0], *work[1][0]);
    avg[NA + 2] -> timesPlus (*work[1][0], *work[1][0]);
    
    if (ND > 2) {
      avg[NA + 3] -> timesPlus (*work[0][0], *work[2][0]);
      avg[NA + 4] -> timesPlus (*work[1][0], *work[2][0]);
      avg[NA + 5] -> timesPlus (*work[2][0], *work[2][0]);
    }

    for (i = 0; i < NT; i++) *avg[i] /= (real) (navg + 1);

  } else {

    // -- Running averages only.

    for (i = 0; i < NA; i++) {
      *avg[i] *= (real)  navg;
      *avg[i] += *src[i];
      *avg[i] /= (real) (navg + 1);
    }
  }

  navg++;
}


void Statistics::dump ()
// ---------------------------------------------------------------------------
// Similar to Domain::dump.
// ---------------------------------------------------------------------------
{
  const integer step     = base.step;
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;

  integer       i;
  ofstream      output;
  const integer NT = avg.getSize();
  const integer ND = Geometry::nDim();
  const integer NR = ((integer) Femlib::value ("AVERAGE") > 1) ? 
                            ((ND + 1) * ND) >> 1 : 0;
  const integer NA = NT - NR;

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

  for (i = 0; i < NA; i++) avg[i] -> transform (-1);
  output << *this;
  for (i = 0; i < NA; i++) avg[i] -> transform (+1);

  ROOTONLY output.close();
}


ofstream& operator << (ofstream&   strm,
		       Statistics& src )
// ---------------------------------------------------------------------------
// Output Statistics class to file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  int               i;
  const int         N = src.avg.getSize();
  vector<AuxField*> field (N);

  for (i = 0; i < N; i++) field[i] = src.avg[i];

  writeField (strm, src.name, src.navg, src.base.time, field);

  return strm;
}


ifstream& operator >> (ifstream&   strm,
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



