//////////////////////////////////////////////////////////////////////////////
// domain.C: implement domain class functions.
//
// Copyright (c) 1994 <--> $Date: 2015/04/29 01:48:00 $, Hugh Blackburn
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: domain.cpp,v 1.2 2015/04/29 01:48:00 hmb Exp $";

#include <sem.h>
#include <stab.h>

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
// NB: we allocate one more base flow AuxField than the problem
// requires: this extra is used as temporary storage in constructing
// gradients (in linAdvect).
// ---------------------------------------------------------------------------
  elmt (E)
{
  const char  routine[] = "Domain::Domain";
  const int_t verbose = Femlib::ivalue ("VERBOSE");
  const int_t nbase   = Geometry::nBase();
  const int_t nz      = Geometry::nZProc();
  const int_t ntot    = Geometry::nTotProc();
  const int_t nplane  = Geometry::planeSize();
  int_t       i, nfield;
  real_t*     alloc;

  strcpy ((name = new char [strlen (F -> root()) + 1]), F -> root());
  Femlib::value ("t", time = 0.0);
  step = 0;

  strcpy ((field = new char [strlen (B -> field()) + 1]), B -> field());
  nfield = strlen (field);
  
  VERBOSE cout << routine << ": Domain will contain fields: " << field << endl;

  // -- Build boundary system and field for each variable.
  
  VERBOSE cout << "  Building domain boundary systems ... ";

  b.resize (nfield);
  for (i = 0; i < nfield; i++) b[i] = new BoundarySys (B, E, field[i]);

  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Building domain perturbation fields ... ";

  u   .resize (nfield);
  udat.resize (nfield);

  alloc = new real_t [static_cast<size_t> (nfield * ntot)];
  for (i = 0; i < nfield; i++) {
    udat[i] = alloc + i * ntot;
    u[i]    = new Field (b[i], udat[i], nz, E, field[i]);
  }

  VERBOSE cout << "done" << endl;

  if   (nbase == 2) strcpy ((baseField = new char [3]), "UV" );
  else              strcpy ((baseField = new char [4]), "UVW");

  VERBOSE cout << "  Building domain base flow fields: "
	       << baseField <<" ... ";

  U   .resize (nbase + 1);
  Udat.resize (nbase + 1);

  alloc = new real_t [static_cast<size_t> ((nbase + 1) * nplane)];
  for (i = 0; i < nbase + 1; i++) {
    Udat[i] = alloc + i * nplane;
    U[i]    = new AuxField (Udat[i], 1, E, baseField[i]);
  }

  VERBOSE cout << "done" << endl;
}


void Domain::report (ostream& file)
// ---------------------------------------------------------------------------
// Print a run-time summary of domain & timestep information on file.
// ---------------------------------------------------------------------------
{
  const real_t dt   = Femlib::value  ("D_T");
  const real_t lz   = Femlib::value  ("TWOPI / BETA");
  const real_t Re   = Femlib::value  ("1.0   / KINVIS");
  const int_t  ns   = Femlib::ivalue ("N_STEP");
  const int_t  nt   = Femlib::ivalue ("N_TIME");
  const int_t  chk  = Femlib::ivalue ("CHKPOINT");
  const int_t  per  = Femlib::ivalue ("IO_FLD");
  const int_t  kdim = Femlib::ivalue ("KRYLOV_KDIM");
  const int_t  nits = Femlib::ivalue ("KRYLOV_NITS");
  const int_t  nvec = Femlib::ivalue ("KRYLOV_NVEC");
  const real_t eps  = Femlib::value  ("KRYLOV_KTOL");
  const int_t  task = Femlib::ivalue ("TASK");

  file << "-- Coordinate system       : ";
  if (Geometry::cylindrical())
    file << "cylindrical" << endl;
  else
    file << "Cartesian"   << endl;

  file << "   Spatial symmetry        : " << Geometry::symmetry() << endl;

#ifdef FLIP
  file << "   with half-period-flip"                              << endl;
#endif

  file << "-- Task for solution       : ";
  if      (task == 1) file << "adjoint ";
  else if (task == 2) file << "growth ";
  else if (task == 3) file << "shrink ";
  else                file << "primal "; 
  file << "problem" << endl;

  file << "-- Krylov dimension        : " << kdim                 << endl;
  file << "   Convergence dimension   : " << nvec                 << endl;
  file << "   Convergence tolerance   : " << eps                  << endl;
  file << "   Maximum iterations      : " << nits                 << endl;

  file << "-- Solution fields         : " << field                << endl;
  file << "   Base flow fields        : " << baseField            << endl;
  file << "   Number of base slices   : " << Geometry::nSlice()   << endl;
  file << "   Number of elements      : " << Geometry::nElmt()    << endl;
  file << "   Number of planes        : " << Geometry::nZ()       << endl;
  file << "   Number of processors    : " << Geometry::nProc()    << endl;
  if (Geometry::nPert() > 2)
    file<<"   Periodic length         : " << lz                   << endl;
  file << "   Polynomial order (np-1) : " << Geometry::nP() - 1   << endl;
  file << "   Reynolds number         : " << Re                   << endl;
  file << "   Time integration order  : " << nt                   << endl;
  file << "   Start time              : " << time                 << endl;
  file << "   Finish time             : " << time + ns * dt       << endl;
  file << "   Time step               : " << dt                   << endl;
  if (Geometry::nSlice() > 1) 
    file<<"   Base flow period        : " << period               << endl;
  file << "   Number of steps         : " << ns                   << endl;
  file << "   Dump interval (steps)   : " << per;
  if (chk) file << " (checkpoint)";  
  file << endl;
}


bool Domain::restart()
// ---------------------------------------------------------------------------
// If a restart file "name".rst can be found, use it for input, and
// return true.  If this fails, initialise all fields to random noise,
// and return false.
// ---------------------------------------------------------------------------
{
  const int_t nF   = nField();
  const int_t ntot = Geometry::nTotProc();
  int_t       i;
  char        restartfile[StrMax];
  bool        restarted = false;
  
  ROOTONLY cout << "-- Initial condition       : ";
  ifstream file (strcat (strcpy (restartfile, name), ".rst"));

  if (file) {
    ROOTONLY cout << "read from file " << restartfile << flush;
    file >> *this;
    file.close();
    restarted = true;
  } else {
    ROOTONLY cout << "randomised";
    for (i = 0; i < nF; i++) Veclib::vnormal (ntot, 0.0, 1.0, udat[i], 1);
  }

  ROOTONLY cout << endl;
  
  Femlib::value ("t", time);
  step = 0;
  
  return restarted;
}


void Domain::dump ()
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
// ---------------------------------------------------------------------------
{
  const bool periodic = !(step %  Femlib::ivalue ("IO_FLD"));
  const bool initial  =   step == Femlib::ivalue ("IO_FLD");
  const bool final    =   step == Femlib::ivalue ("N_STEP");

  if (!(periodic || final)) return;
  ofstream output;

  Femlib::synchronize();

  ROOTONLY {
    const char  routine[] = "Domain::dump";
    char        dumpfl[StrMax], backup[StrMax], command[StrMax];
    const int_t verbose   = Femlib::ivalue ("VERBOSE");
    const int_t chkpoint  = Femlib::ivalue ("CHKPOINT");

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
    if (verbose) message (routine, "writing field dump",  REMARK);
  }

  output << *this;

  ROOTONLY output.close();
}


ofstream& operator << (ofstream& strm,
		       Domain&   D   )
// ---------------------------------------------------------------------------
// Output all Domain field variables on ostream in prism-compatible
// form.  Binary output only.  Note that output is only done on root
// processor.
// ---------------------------------------------------------------------------
{
  const int_t       N = D.u.size();
  vector<AuxField*> field (N);

  for (int_t i = 0; i < N; i++) field[i] = D.u[i];

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
  int_t      i, j, np, nz, nel, ntot, nfields;
  int_t      npchk,  nzchk, nelchk;
  int_t      verb = Femlib::ivalue ("VERBOSE");
  char       s[StrMax], f[StrMax], err[StrMax], fields[StrMax];
  bool       swap;

  if (strm.getline(s, StrMax).eof()) return strm;

  strm.getline(s,StrMax).getline(s,StrMax);
  
  string ss(s);
  istringstream sss (ss);
  sss >> np >> np >> nz >> nel;

  D.u[0] -> describe (f);
  sss.clear();
  sss.str  (ss = f);
  sss >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline(s,StrMax);
  sss.clear();
  sss.str  (ss = s);
  sss >> D.step;

  strm.getline(s,StrMax);
  sss.clear();
  sss.str  (ss = s);
  sss >> D.time;
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
    sprintf (err,"file: %1d fields, Domain: %1d",nfields,(int)strlen(D.field));
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


void Domain::loadBase()
// ---------------------------------------------------------------------------
// Open file name.bse, assumed to contain field dumps for base flow.
// Attempt to load velocity fields (ignoring pressure), then Fourier
// transform if appropriate, and prepare for Fourier reconstruction.
//
// Base flow dumps should contain only velocity and pressure fields,
// have same number of elements, np as perturbation.  Only the first
// plane of data in each dump is used, higher values are ignored
// (i.e. the base flow is treated as 2D, no matter what).
// ---------------------------------------------------------------------------
{
  const char  routine[] = "Domain::loadBase()";
  const int_t nP     = Geometry::nP();
  const int_t nEl    = Geometry::nElmt();
  const int_t nBase  = Geometry::nBase();
  const int_t nPlane = Geometry::nPlane();
  const int_t nTot   = Geometry::planeSize();
  const int_t nSlice = Geometry::nSlice();

  Header   H;
  int_t    i, j;
  real_t*  addr;
  real_t   t0, dt = 0;
  int_t    len;
  char     filename[StrMax];
  ifstream file (strcat (strcpy (filename, name), ".bse"));

  // -- Allocate storage.

  baseFlow.resize (nBase);

  if (nSlice > 1)
    for (i = 0; i < nBase; i++) baseFlow[i] = new real_t [nTot * nSlice];
  else
    for (i = 0; i < nBase; i++) baseFlow[i] = Udat[i];

  // -- Read base velocity fields, ignore pressure fields.

  ROOTONLY cout << "-- Base flow               : " << flush; 

  for (i = 0; file >> H; i++) {

    if (H.nr != nP || H.nel != nEl)
      message (routine, "base flow and perturbation do not conform", ERROR);
    if ((nBase == 2 && strcmp (H.flds, "uvp" )) ||
	(nBase == 3 && strcmp (H.flds, "uvwp")))
      message (routine, "mismatch: No. of base components/declaration", ERROR);

    for (j = 0; j < nBase; j++) {

      // -- Note that we wrap around so that the last dump goes first:
      addr = baseFlow[j] + ((i+1) % nSlice) * nTot;
      
      len = nPlane * sizeof(real_t);
      file.read (reinterpret_cast<char*>(addr), len);
      if (H.swab()) Veclib::brev (nTot, addr, 1, addr, 1);
    
      len = (H.nz - 1) * nPlane * sizeof (real_t);
      file.ignore (len); // -- Ignore higher planes.

      if (file.bad()) message (routine, "unable to read binary input", ERROR);
      Veclib::zero (nTot - nPlane, addr + nPlane, 1);
    }
    
    len = H.nz * nPlane * sizeof (real_t);
    file.ignore (len); // -- Ignore pressure field.

    if (i == 0) t0 = H.time;
    dt = H.time - t0;
  }

  file.close();

  if (i != nSlice)
    message (routine, "mismatch: No. of base slices/declaration", ERROR);

  if (nSlice > 1) {		// -- Prepare for base flow
                                //    reconstruction.  Default is
                                //    Fourier, but now we can also
                                //    slect 4-point (cubic) Lagrange
                                //    interpolation.

    period = Femlib::value ("BASE_PERIOD"); // -- Use this if installed.
    if (period < EPSDP)
      Femlib::value ("BASE_PERIOD", period = dt * i / (i - 1.0));

    if (! (Femlib::ivalue ("LAGRANGE_INT"))) {

	// -- Fourier transform in time, scale for reconstruction.
	for (i = 0; i < nBase; i++) {
	  Femlib::DFTr (baseFlow[i], nSlice, nTot, FORWARD);
	  Blas::scal   ((nSlice-2)*nTot, 2.0, baseFlow[i] + 2*nTot, 1);
	}
      }
  } else period = 0.0;

  ROOTONLY {
    cout << "read from file " << filename;
    if (H.swab()) cout << " (byte swapping)";
    cout << endl;
  }
}


void Domain::updateBase()
// ---------------------------------------------------------------------------
// Update base velocity fields by interpolation.
// ---------------------------------------------------------------------------
{
  const int_t nBase   = Geometry::nBase();
  const int_t nSlice  = Geometry::nSlice();
  int_t       i;
  
  if (nSlice < 2) return;

  for (i = 0; i < nBase; i++)
    U[i] -> update (nSlice, baseFlow[i], time, period);
}
