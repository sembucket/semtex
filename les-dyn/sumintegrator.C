///////////////////////////////////////////////////////////////////////////////
// sumintegrator.C: implement SumIntegrator class.
//
// Copyright (C) 2000 Hugh Blackburn
//
// This class exists to carry out first-order system temporal
// smoothing of the dynamic estimate of Smagorinsky constant.
// 
// Given token "TIME_CONST", which is the characteristic time for the
// first order system, this class initialises internal storage and
// provides a method for the update
//
//   Cs^2[n] = A*Cs^2_dyn + B*Cs^2[n-1],
//
// where A+B=1 and B = exp (-TIME_CONST * D_T).  This carries out a
// first-order system relaxation towards the ensemble-average dynamic
// estimate.
//
// The class relies on files session.smg (which contains the last
// finalised estimate of Cs^2) and session.sma (which contains the
// checkpoint estimate).  If it can't find session.smg, the
// constructor will open it, fill it with C_SMAG*C_SMAG.  The file
// update and rollover techniques are similar to what happens in the
// Statistics class.
//
// The update method is designed to be called every timestep.
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include "les.h"


SumIntegrator::SumIntegrator (Domain* info) :
// ---------------------------------------------------------------------------
// Open files, initialise storage and constants.
// ---------------------------------------------------------------------------
  _domain (info),
  _BB     (Femlib::value ("exp (-TIME_CONST * D_T)")),
  _AA     (1.0 - _BB),
  _ntot   (Geometry::nTotProc()),
  _nz     (Geometry::nZProc())
{
  _work = new real [(size_t) _ntot];
  _Cs2 = new AuxField (_work, _nz, _domain -> elmt, 'x');

  ROOTONLY cout << "-- Initialising Smagorinsky constant : ";

  char s[StrMax];
  ifstream file (strcat (strcpy (s, _domain -> name), ".smg"));

  if (file) {
    ROOTONLY { cout << "read from file " << s; cout.flush(); }
    file >> *this;
    file.close();
  } else {
    const real cs2 = Femlib::value ("C_SMAG*C_SMAG");
    ROOTONLY cout << " set to " << cs2;
    *_Cs2 = cs2;
  }

  ROOTONLY cout << endl;
}


void SumIntegrator::update (real* src)
// ---------------------------------------------------------------------------
// Take the new dynamic estimate and do a relaxation update.
//
// Actually the input is -2*Cs^2.
// ---------------------------------------------------------------------------
{ 
  Veclib::smul  (_ntot, _BB,      _work, 1, _work, 1);
  Veclib::svtvp (_ntot, -0.5*_AA, src,   1, _work, 1, _work, 1);
  Veclib::smul  (_ntot, -2.0,     _work, 1, src,   1);
}


void SumIntegrator::dump ()
// ---------------------------------------------------------------------------
// Dump out internal storage to file.
// ---------------------------------------------------------------------------
{
  const integer step     = _domain -> step;
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;

  ofstream output;

  ROOTONLY {
    const char    routine[] = "SumIntegrator::dump";
    const integer verbose   = (integer) Femlib::value ("VERBOSE");
    const integer chkpoint  = (integer) Femlib::value ("CHKPOINT");
    char          dumpfl[StrMax], backup[StrMax], command[StrMax];

    if (chkpoint) {
      if (final) {
	strcat (strcpy (dumpfl, _domain -> name), ".smg");
	output.open (dumpfl, ios::out);
      } else {
	strcat (strcpy (dumpfl, _domain -> name), ".sma");
	if (!initial) {
	  strcat  (strcpy (backup, _domain -> name), ".sma.bak");
	  sprintf (command, "mv ./%s ./%s", dumpfl, backup);
	  system  (command);
	}
	output.open (dumpfl, ios::out);
      }
    } else {
      strcat (strcpy (dumpfl, _domain -> name), ".smg");
      if   (initial) output.open (dumpfl, ios::out);
      else           output.open (dumpfl, ios::app);
    }

    if (!output) message (routine, "can't open dump file", ERROR);
    if (verbose) message (routine, ": writing field dump", REMARK);
  }
  
  output << *this;

  ROOTONLY output.close();
}


ofstream& operator << (ofstream&      strm,
		       SumIntegrator& src )
// ---------------------------------------------------------------------------
// Output SumIntegrator class to file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  const Domain* D = src._domain;
  vector<AuxField*> tmp(1); tmp[0] = src._Cs2;

  writeField (strm, D -> name, D -> step, D -> time, tmp);

  return strm;
}


ifstream& operator >> (ifstream&      strm,
		       SumIntegrator& tgt )
// ---------------------------------------------------------------------------
// Input SumIntegrator class from file.  Like similar Domain routine.
// ---------------------------------------------------------------------------
{
  const char routine[] = "strm>>SumIntegrator";
  integer    np, nz, nel, ntot;
  integer    npchk,  nzchk, nelchk, swap = 0;
  char       s[StrMax], f[StrMax];

  if (strm.getline(s, StrMax).eof()) return strm;
  
  strm.getline (s, StrMax) . getline (s, StrMax);

  tgt._Cs2 -> describe (f);
  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline (s,StrMax);
  strm.getline (s, StrMax) . getline (s, StrMax);
  strm.getline (s, StrMax) . getline (s, StrMax) . getline (s, StrMax);
    
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
    
  strm >> *tgt._Cs2;
  if (swap) tgt._Cs2 -> reverse();
  
  ROOTONLY if (strm.bad())
    message (routine, "failed reading Smagorinsky file", ERROR);

  return strm;
}

