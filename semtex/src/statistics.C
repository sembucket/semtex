///////////////////////////////////////////////////////////////////////////////
// statistics.C: routines for statistical analysis of AuxFields.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// What statistics are collected is controlled by the setting of the
// AVERAGE token. Legal values are 0 (default), 1, 2, 3. The routines
// here do not control how often statistics are updated: that happens
// in Analyser class methods.
//
// All collected statistics are given 1-character names. Potentially
// these may clash with quantities derived/used by addfield utility
// for example. Caveat emptor. In the longer term we should move to
// CSV for field names.
//
// AVERAGE = 0.  No statistics.
//
// AVERAGE = 1.  Running averages of variables held by the Domain used
// for initialisation. The data are held in semi-Fourier state.
//
// AVERAGE = 2.  Additionally, product terms for computation of
// Reynolds stresses. (All product terms are computed and held in
// physical space.)
// 
// Naming scheme for components of the symmetric "Reynolds stresses" tensor:
//
//                      / uu uv uw \     /  A  B  D \
//                      | .  vv vw |  =  |  .  C  E |
//                      \ .  .  ww /     \  .  .  F /
//
// What is computed are the running average of the products uu, uv,
// etc, which are NOT the actual Reynolds stresses: they need to have
// the products of the mean values UU, UV etc subtracted, assumed to
// occur in postprocessing. If there are only two velocity components
// present in the initialising Domain, only averages for A, B & C are
// made.
//
// AVERAGE = 3. Additional products are kept for computation of energy
// equation terms. Again, the correct terms need to be made in
// post-processing.
//
// a) Scalar: 
//    i) q = 0.5 [u^2 + v^2 (+ w^2)]
//   ii) d = SijSij
//
// b) Vector: Naming:
//    i) p u_i
//                      / pu \   / m \
//                      | pv | = | n |
//                      \ pw /   \ o /
//   ii) q u_i
//                      / qu \   / r \
//                      | qv | = | s |
//                      \ qw /   \ t /
// 
// c) Tensor: symmetric rate-of-strain tensor S_ij. Naming:
//
//                      / xx xy xz \     /  G  H  J \
//                      | .  yy yz |  =  |  .  I  K |
//                      \ .  .  zz /     \  .  .  L /
//
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <ctime>
#include <sem.h>


Statistics::Statistics (Domain* D) :
// ---------------------------------------------------------------------------
// Store averages for all Domain Fields.
//
// Try to initialize from file session.avg, failing that set all
// buffers to zero.  Number of fields in file should be same as
// Statistics::avg buffer.
//
// NR = Number of Reynolds stress averaging buffers, set if AVERAGE > 1.
// ---------------------------------------------------------------------------
  _name (D -> name),
  _base (D),
  _iavg (Femlib::ivalue ("AVERAGE"))
{
  if (_iavg == 0) return;
  if ((_iavg  < 0) || (_iavg > 2))
    message ("Statistics::Statistics", "AVERAGE token out of [0,3]", ERROR);
					 
  int_t       i, j;
  const int_t NF    = _base -> nField();
  const int_t NC    = NF - 1;	// -- Number of velocity components.
  const int_t NR    = (_iavg > 1) ? ((NC+1)*NC)/2   : 0;
  const int_t NE    = (_iavg > 2) ? ((NC+5)*NC+4)/2 : 0;
  const int_t NT    = NF + NR + NE;
  const int_t nz    = Geometry::nZProc();
  const int_t ntot  = Geometry::nTotProc();
  real_t*     alloc = new real_t [static_cast<size_t> (NT * ntot)];

  // -- Set pointers, allocate storage.

  _src.resize (NF);	      // -- Straight running average of these.
  _avg.resize (NT);	      // -- Additional, computed from _src.
  
  for (i = 0; i < NF; i++) _src[i] = (AuxField*) _base -> u[i];

  if (_iavg > 0)
    for (j = 0, i = 0; i < NF; i++, j++)
      _avg[i] = new AuxField (alloc+j*ntot, nz, _base->elmt, _src[i]->name());
  if (_iavg > 1)
    for (i = 0; i < NT - NF; i++, j++)
      _avg[i + NF] = new AuxField (alloc+j*ntot, nz, _base->elmt, 'A' + i);
  if (_iavg > 2) {
    _avg[NF + NR]     = new AuxField (alloc+(j++)*ntot, nz, _base->elmt, 'q');
    _avg[NF + NR + 1] = new AuxField (alloc+(j++)*ntot, nz, _base->elmt, 'd');
  }
}


void Statistics::initialise ()
// ---------------------------------------------------------------------------
// This is for standard running averages. Try to initialize from file
// session.avg, failing that set all buffers to zero.  Number of
// fields in file should be same as Statistics::avg buffer.
//
// NR = Number of Reynolds stress averaging buffers, set if AVERAGE > 1.
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t NF = _base -> nField();
  const int_t NC = NF - 1;	// -- Number of velocity components.
  const int_t NR = (_iavg > 1) ? ((NC+1)*NC)>>1 : 0;
  const int_t NT = NF + NR;

  ROOTONLY cout << "-- Initialising averaging  : ";  

  // -- Initialise averages, either from file or zero.
  //    This is much the same as Domain input routine.

  char     s[StrMax];
  ifstream file (strcat (strcpy (s, _name), ".avg"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << s;
      cout.flush();
    }
    file >> *this;
    file.close();
    for (i = 0; i < NT - NR; i++) _avg[i] -> transform (FORWARD);
  
  } else {			// -- No file, set to zero.
    ROOTONLY cout << "set to zero";
    for (i = 0; i < NT; i++) *_avg[i] = 0.0;
    _navg = 0;
  }

  ROOTONLY cout << endl;
}


void Statistics::update (AuxField** work)
// ---------------------------------------------------------------------------
// Update running averages, using zeroth time level of work as
// workspace.  Reynolds stress terms are calculated without
// dealiasing, and are held in physical space.
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t NT = _avg.size();
  const int_t NC = _base -> nField() - 1;
  const int_t NR = (_iavg > 1) ? ((NC+1)*NC)>>1 : 0;
  const int_t NA = NT - NR;

  if (NR) {
    
    // -- Running averages and Reynolds stresses.

    for (i = 0; i < NC; i++) {
      *work[i] = *_src[i];
       work[i] -> transform (INVERSE);
    }

    for (i = 0; i < NT; i++) *_avg[i] *= static_cast<real_t>(_navg);

    for (i = 0; i < NA; i++) *_avg[i] += *_src[i];

    _avg[NA + 0] -> timesPlus (*work[0], *work[0]);
    _avg[NA + 1] -> timesPlus (*work[0], *work[1]);
    _avg[NA + 2] -> timesPlus (*work[1], *work[1]);
    
    if (NC > 2) {
      _avg[NA + 3] -> timesPlus (*work[0], *work[2]);
      _avg[NA + 4] -> timesPlus (*work[1], *work[2]);
      _avg[NA + 5] -> timesPlus (*work[2], *work[2]);
    }

    for (i = 0; i < NT; i++) *_avg[i] /= static_cast<real_t>(_navg + 1);

  } else {

    // -- Running averages only.

    for (i = 0; i < NA; i++) {
      *_avg[i] *= static_cast<real_t>(_navg);
      *_avg[i] += *_src[i];
      *_avg[i] /= static_cast<real_t>(_navg + 1);
    }
  }

  _navg++;
}


void Statistics::dump ()
// ---------------------------------------------------------------------------
// Similar to Domain::dump.
// ---------------------------------------------------------------------------
{
  const int_t step     = _base -> step;
  const bool  periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool  initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool  final    =   step == Femlib::ivalue ("N_STEP");

  if (!(periodic || final)) return;

  int_t       i;
  ofstream    output;
  const int_t NT = _avg.size();
  const int_t NC = Geometry::nDim();
  const int_t NR = (_iavg > 1) ? ((NC+1)*NC)>> 1 : 0;
  const int_t NA = NT - NR;

  Femlib::synchronize();

  ROOTONLY {
    const char  routine[] = "Statistics::dump";
    const int_t verbose   = Femlib::ivalue ("VERBOSE");
    const int_t chkpoint  = Femlib::ivalue ("CHKPOINT");
    char        dumpfl[StrMax], backup[StrMax], command[StrMax];

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, _name), ".avg");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, _name), ".ave");
	if (!initial) {
	  strcat  (strcpy (backup, _name), ".ave.bak");
	  sprintf (command, "mv ./%s ./%s", dumpfl, backup);
	  system  (command);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, _name), ".avg");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }

    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }

  Femlib::synchronize();
  for (i = 0; i < NA; i++) _avg[i] -> transform (INVERSE);
  Femlib::synchronize();

  output << *this;

  Femlib::synchronize();
  for (i = 0; i < NA; i++) _avg[i] -> transform (FORWARD);
  Femlib::synchronize();

  ROOTONLY output.close();
}


ofstream& operator << (ofstream&   strm,
		       Statistics& src )
// ---------------------------------------------------------------------------
// Output Statistics class to file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  int_t             i;
  const int_t       N = src._avg.size();
  vector<AuxField*> field (N);

  for (i = 0; i < N; i++) field[i] = src._avg[i];

  writeField (strm, src._name, src._navg, src._base -> time, field);

  return strm;
}


ifstream& operator >> (ifstream&   strm,
		       Statistics& tgt )
// ---------------------------------------------------------------------------
// Input Statistics class from file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  const char routine[] = "strm>>Statistics";
  int_t i, j, np, nz, nel, ntot, nfields;
  int_t npchk,  nzchk, nelchk;
  char  s[StrMax], f[StrMax], err[StrMax], fields[StrMax];
  bool  swap = false;

  if (strm.getline(s, StrMax).eof()) return strm;
  
  strm.getline (s, StrMax) . getline (s, StrMax);

  tgt._avg[0] -> describe (f);
  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline (s,StrMax);
  istrstream (s, strlen (s)) >> tgt._navg;
    
  strm.getline (s, StrMax) . getline (s, StrMax);
  strm.getline (s, StrMax) . getline (s, StrMax) . getline (s, StrMax);
    
  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';
  if (nfields != tgt._avg.size()) {
    sprintf (err, "strm: %1d fields, avg: %1d", nfields, tgt._avg.size());
    message (routine, err, ERROR);
  }
  for (i = 0; i < nfields; i++) 
    if (!strchr (fields, tgt._avg[i] -> name())) {
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
      if (tgt._avg[i] -> name() == fields[j]) break;
    strm >> *tgt._avg[i];
    if (swap) tgt._avg[i] -> reverse();
  }
  
  ROOTONLY if (strm.bad())
    message (routine, "failed reading average file", ERROR);

  return strm;
}


void Statistics::phaseUpdate (const int_t j   ,
			      AuxField**  work)
// ---------------------------------------------------------------------------
// Phase updates are running updates, like those for standard
// statistics, but they are only computed at (a presumed very limited)
// number of instants per run. And so there would be a number of
// averaging files, to be called e.g. session.0.phs, session.1.phs
// ... session.(N-1).phs. Instead of reserving enough memory for all
// these buffers, we just keep workspace reserved, and at every phase
// point the appropriate file (number j) is read in, updated, and
// written back out. If the file does not exist, it is created.
//
// NB: the number of averages computed (needed for the running
// averaging) is only updated when j == 0, which will happen as the
// last of a set of phase averages.
// ---------------------------------------------------------------------------
{
  int_t       i;
  const int_t NF = _base -> nField();
  const int_t NC = NF - 1;	// -- Number of velocity components.
  const int_t NR = (_iavg > 1) ? ((NC+1)*NC)>>1 : 0;
  const int_t NT = NF + NR;
  const bool  verbose = static_cast<bool> (Femlib::ivalue ("VERBOSE"));

  // -- Initialise averages, either from file or zero.
  //    This is much the same as Domain input routine.

  char     s[StrMax];
  sprintf  (s, "%s.%1d.phs", _name, j);
  ifstream ifile (s);
  
  VERBOSE cout << "-- Updating phase average " << j << ": ";

  if (ifile) {
    VERBOSE {
      cout << "read from file " << s;
      cout.flush();
    }
    ifile >> *this;
    ifile.close();
    for (i = 0; i < NT - NR; i++) _avg[i] -> transform (FORWARD);
  
  } else {			// -- No file, set to zero.
    VERBOSE cout << "set to zero";
    for (i = 0; i < NT; i++) *_avg[i] = 0.0;
    _navg = 0;
  }

  VERBOSE cout << endl;

  // -- Stuff has been read in to local buffers. Now do running average.

  if (NR) {
    
    // -- Running averages and Reynolds stresses.

    for (i = 0; i < NC; i++) {
      *work[i] = *_src[i];
       work[i] -> transform (INVERSE);
    }

    for (i = 0; i < NF; i++) *_avg[i] *= static_cast<real_t>(_navg);

    for (i = 0; i < NF; i++) *_avg[i] += *_src[i];

    _avg[NF + 0] -> timesPlus (*work[0], *work[0]);
    _avg[NF + 1] -> timesPlus (*work[0], *work[1]);
    _avg[NF + 2] -> timesPlus (*work[1], *work[1]);
    
    if (NC > 2) {
      _avg[NF + 3] -> timesPlus (*work[0], *work[2]);
      _avg[NF + 4] -> timesPlus (*work[1], *work[2]);
      _avg[NF + 5] -> timesPlus (*work[2], *work[2]);
    }

    for (i = 0; i < NT; i++) *_avg[i] /= static_cast<real_t>(_navg + 1);

  } else {

    // -- Running averages only.

    for (i = 0; i < NF; i++) {
      *_avg[i] *= static_cast<real_t>(_navg);
      *_avg[i] += *_src[i];
      *_avg[i] /= static_cast<real_t>(_navg + 1);
    }
  }

  // -- Increment the number of averages for output.
  //    NB: This value is re-read from file each time we do an update,
  //    so this increment does not get to corrupt other phase points.

  _navg++;

  // -- Write the updated averages back out to file.

  ofstream ofile (s);

  for (i = 0; i < NF; i++) _avg[i] -> transform (INVERSE);

  ofile << *this;

  for (i = 0; i < NF; i++) _avg[i] -> transform (FORWARD);

  ROOTONLY ofile.close();
}
