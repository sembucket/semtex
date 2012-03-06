///////////////////////////////////////////////////////////////////////////////
// statistics.C: routines for statistical analysis of AuxFields,
// specialised for use with dual -- so there are no Fourier transforms.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// The collection of statistics is controlled by the setting of the
// AVERAGE token. Legal values are 0 (default), 1, 2. The routines
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
// AVERAGE = 2.  Additionally, correlation terms for computation of
// Reynolds stresses. (Correlations, based on products of variables,
// are computed and held in physical space.)
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
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>


Statistics::Statistics (Domain* D) :
// ---------------------------------------------------------------------------
// Store averages for all Domain Fields, and correlations.
// ---------------------------------------------------------------------------
  _name (D -> name),
  _base (D),
  _iavg (Femlib::ivalue ("AVERAGE")),
  _nraw (_base -> nField()),
  _nvel (_nraw - 1),
  _nrey ((_iavg > 1) ? ((_nvel+1)*_nvel)/2 : 0),
  _neng (0)
{
  if (_iavg == 0) return;
  if ((_iavg  < 0) || (_iavg > 2))
    message ("Statistics::Statistics", "AVERAGE token out of [0,2]", ERROR);
					 
  int_t       i;
  const int_t nz   = Geometry::nZProc();
  const int_t ntot = Geometry::nTotProc();

  // -- Set pointers, allocate storage.

  for (i = 0; i < _nraw; i++)	// -- Local pointers to raw variables.
    _raw[_base -> u[i] -> name()] = (AuxField*) _base -> u[i];

  if (_iavg > 0) // -- Set up buffers for averages of raw variables.
    for (i = 0; i < _nraw; i++)
      _avg[_base -> u[i] -> name()] =
	new AuxField (new real_t[ntot],nz,_base->elmt,_base->u[i]->name());

  if (_iavg > 1) // -- Set up buffers for Reynolds stress correlations.
    for (i = 0; i < _nrey; i++)
      _avg['A' + i] = new AuxField (new real_t[ntot],nz,_base->elmt,'A'+i);
}


void Statistics::initialise (const char* filename)
// ---------------------------------------------------------------------------
// This is for standard running averages. Try to initialize from file
// filename (e.g. "session.avg"), failing that set all buffers to
// zero.  Number of fields in file should be same as Statistics::_avg
// buffer.
// ---------------------------------------------------------------------------
{
  cout << "-- Initialising averaging  : ";  

  ifstream file (filename);
  map<char, AuxField*>::iterator k;

  if (file) {
    cout << "read from file " << filename;
    cout.flush();

    file >> *this;
    file.close();
  
  } else {			// -- No file, set to zero.
    cout << "set to zero";
    for (k = _avg.begin(); k != _avg.end(); k++) *(k -> second) = 0.0;
    _navg = 0;
  }

  cout << endl;
}


void Statistics::update (AuxField** wrka,
			 AuxField** wrkb)
// ---------------------------------------------------------------------------
// Update running averages, using arrays wrka & wrkb as workspace.
// 
// NB: Reynolds stress products are done by convolution in dual.
// ---------------------------------------------------------------------------
{
  if (_iavg < 1) return;
  
  char   key;
  int_t  i, j;
  Field* master = _base -> u[0];
  map<char, AuxField*>::iterator k;

  // -- Weight old running averages.

  for (k = _avg.begin(); k != _avg.end(); k++)
    *(k -> second) *= static_cast<real_t>(_navg);
    
  // -- Always do running averages of raw data.

  for (k = _raw.begin(); k != _raw.end(); k++)
    *_avg[k -> second-> name()] += *(k -> second);

  // -- Reynolds stress correlations.

  if (_iavg > 1) {
    *wrka[0] = *_raw['u'];
    *wrka[1] = *_raw['v'];
    if (_nvel == 3) *wrka[2] = *_raw['w']; 

    wrkb[0] -> convolve (*wrka[0], *wrka[0]); *_avg['A'] += *wrkb[0];
    wrkb[0] -> convolve (*wrka[0], *wrka[1]); *_avg['B'] += *wrkb[0];
    wrkb[0] -> convolve (*wrka[1], *wrka[1]); *_avg['C'] += *wrkb[0];
    
    if (_nvel == 3) {
      wrkb[0] -> convolve (*wrka[0], *wrka[2]); *_avg['D'] += *wrkb[0];
      wrkb[0] -> convolve (*wrka[1], *wrka[2]); *_avg['E'] += *wrkb[0];
      wrkb[0] -> convolve (*wrka[2], *wrka[2]); *_avg['F'] += *wrkb[0];
    }
  }

  // -- Normalise and smooth running averages.

  for (k = _avg.begin(); k != _avg.end(); k++) {
    master -> smooth (k -> second);
    *(k -> second) /= static_cast<real_t>(_navg + 1);
  }

  _navg++;
}


void Statistics::dump (const char* filename)
// ---------------------------------------------------------------------------
// Similar to Domain::dump.
//
// As of 24/11/2004, we deleted the checkpointing that used to happen:
// all dumping now happens to file named on input.
//
// We also smooth all the outputs with the mass matrix.
// ---------------------------------------------------------------------------
{
  const int_t step     = _base -> step;
  const bool  periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool  initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool  final    =   step == Femlib::ivalue ("N_STEP");

  if (!(periodic || final)) return;

  ofstream    output;
  int_t       i;
  map<char, AuxField*>::iterator k;

  const char routine[] = "Statistics::dump";
  const bool verbose   = static_cast<bool> (Femlib::ivalue ("VERBOSE"));

  output.open (filename);
  if (!output) message (routine, "can't open dump file", ERROR);
  if (verbose) message (routine, ": writing field dump", REMARK);
  
  // -- All terms are written out in Fourier space.

  output << *this;

  output.close();
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
  map<char, AuxField*>::iterator k;

  for (i = 0, k = src._avg.begin(); i < N; i++, k++) field[i] = k -> second;

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
  map<char, AuxField*>::iterator k;

  if (strm.getline(s, StrMax).eof()) return strm;
  
  strm.getline (s, StrMax) . getline (s, StrMax);

  string ss(s);
  istringstream sss (ss);
  sss >> np >> np >> nz >> nel;
 
  tgt._avg.begin()->second->describe (f);
  sss.clear ();
  sss.str   (ss = f);
  sss >> npchk >> npchk >> nzchk >> nelchk;

  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline (s, StrMax);

  sss.clear();
  sss.str  (ss = s);
  sss >> tgt._navg;
    
  strm.getline (s, StrMax) . getline (s, StrMax);
  strm.getline (s, StrMax) . getline (s, StrMax) . getline (s, StrMax);
    
  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';
  if (nfields != tgt._avg.size()) {
    sprintf (err, "strm: %1d fields, avg: %1d", 
	     nfields, static_cast<int_t>(tgt._avg.size()));
    message (routine, err, ERROR);
  }

  for (i = 0, k = tgt._avg.begin(); k != tgt._avg.end(); k++, i++)
    if (!strchr (fields, k -> second -> name())) {
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
    if (swap) cout << " (byte-swapping)";
    cout.flush();
  }

  for (j = 0; j < nfields; j++) {
    strm >> *tgt._avg[fields[j]];
    if (swap) tgt._avg[fields[j]] -> reverse();
  }
  
  if (strm.bad())
    message (routine, "failed reading average file", ERROR);

  return strm;
}


void Statistics::phaseUpdate (const int_t j   ,
			      AuxField**  wrka,
			      AuxField**  wrkb)
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
  char filename[StrMax];
  sprintf (filename, "%s.%1d.phs", _name, j);

  this -> initialise (filename);
  this -> update     (wrka, wrkb);
  this -> dump       (filename);
}
