///////////////////////////////////////////////////////////////////////////////
// domain.C: implement domain class functions.
//
// Copyright (C) 1994, 1999 Hugh Blackburn
//
// $Id$
///////////////////////////////////////////////////////////////////////////////

#include <Sem.h>

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
// No initialization of Field MatrixSystems.
//
// New Base Fields will have Uppercase letters to distinquish from
// perturbed fields (JE - 26/5/99)
// ---------------------------------------------------------------------------
  elmt (E)
{
  const integer verbose = (integer) Femlib::value ("VERBOSE");
  integer       i, nfield;
  const integer nz   = Geometry::nZProc();
  const integer ntot = Geometry::nTotProc();
  real*         alloc;

  strcpy ((name = new char [strlen(F->root())+1]), F -> root());
  Femlib::value ("t", time = 0.0);
  step = 0;

  strcpy ((field = new char [strlen (B -> field()) + 1]), B -> field());
  nfield = strlen (field);

#ifdef STABILITY

  // make same length for Domain::loadfield reasons.
  strcpy ((basefield = new char [nfield]),"UV");

  // -- Output base field

  VERBOSE cout << "  Domain will contain base fields: " << basefield << endl;

  U   .setSize(2);
  Udat.setSize(2);
  // memory and auxfields are created as required in loadbase.

#endif 
  
  VERBOSE cout << "  Domain will contain fields: " << field << endl;

  // -- Build boundary system and field for each variable.
  
  VERBOSE cout << "  Building domain boundary systems ... ";

  b.setSize (nfield);
  for (i = 0; i < nfield; i++) b[i] = new BoundarySys (B, E, field[i]);

  VERBOSE cout << "done" << endl;

  VERBOSE cout << "  Building domain fields ... ";

  u   .setSize (nfield);
  udat.setSize (nfield);

  alloc = new real [(size_t) nfield * ntot];
  for (i = 0; i < nfield; i++) {
    udat[i] = alloc + i * ntot;
    u[i]    = new Field (b[i], udat[i], nz, E, field[i]);
  }

  VERBOSE cout << "done" << endl;
}


void Domain::report ()
// ---------------------------------------------------------------------------
// Print a run-time summary of domain & timestep information on cout.
// ---------------------------------------------------------------------------
{
  const real    t   =           time;
  const real    dt  =           Femlib::value ("D_T");
  const real    lz  =           Femlib::value ("TWOPI / BETA");
  const integer ns  = (integer) Femlib::value ("N_STEP");
  const integer nt  = (integer) Femlib::value ("N_TIME");
  const integer chk = (integer) Femlib::value ("CHKPOINT");
  const integer per = (integer) Femlib::value ("IO_FLD");

  cout << "-- Coordinate system       : ";
  if (Geometry::system() == Geometry::Cylindrical)
    cout << "cylindrical" << endl;
  else
    cout << "Cartesian" << endl;

  cout << "   Solution fields         : " << field              << endl;

#ifdef STABILITY
  
  cout << "   Base Solution fields    : " << basefield          << endl;
  cout << "   Base flow Dimension     : 2" << endl;
  cout << "   Perturbed Dimension     : " << nField()-1          << endl;

#endif

  cout << "   Number of elements      : " << Geometry::nElmt()  << endl;
  cout << "   Number of planes        : " << Geometry::nZ()     << endl;
  cout << "   Number of processors    : " << Geometry::nProc()  << endl;
  if (Geometry::nZ() > 1) cout << "   Periodic length         : " << lz<< endl;
  cout << "   Polynomial order (np-1) : " << Geometry::nP() - 1 << endl;

  cout << "   Time integration order  : " << nt                 << endl;
  cout << "   Start time              : " << t                  << endl;
  cout << "   Finish time             : " << t + ns * dt        << endl;
  cout << "   Time step               : " << dt                 << endl;
  cout << "   Number of steps         : " << ns                 << endl;
  cout << "   Dump interval (steps)   : " << per;
  if (chk) cout << " (checkpoint)";  
  cout << endl;
}


void Domain::restart ()
// ---------------------------------------------------------------------------
// If a restart file "name".rst can be found, use it for input.  If
// this fails, initialize all Fields to zero ICs.
//
// Carry out forwards Fourier transformation, zero Nyquist data.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer nF = nField();
  char          restartfile[StrMax];
  
  ROOTONLY cout << "-- Initial condition       : ";
  ifstream file (strcat (strcpy (restartfile, name), ".rst"));

  if (file) {
    ROOTONLY {
      cout << "read from file " << restartfile;
      cout.flush();
    }
    file >> *this;
    file.close();
    transform (FORWARD);
    ROOTONLY for (i = 0; i < nF; i++) u[i] -> zeroNyquist();
  } else {

#ifdef STABILITY
    cout << "generating random noise";
    for (i = 0; i < nF; i++) u[i] -> randomise();
#else
    ROOTONLY cout << "set to zero";
    for (i = 0; i < nF; i++) *u[i] = 0.0;
#endif
  }

  ROOTONLY cout << endl;
  
  Femlib::value ("t", time);
  step = 0;
}


void Domain::loadbase ()
// ---------------------------------------------------------------------------
// If a base file "base.bse" can be found, use it for input.  If
// this fails, initialize all Fields to random noise.
//
// Loads up base fields (U+V) -> adapts for multiple base fields.
// Sets -> period, d_t, PBF, U, Udat
  //
  // very clumsy at moment because of difference between passing auxfield and
  // field.. in loadfield
// ---------------------------------------------------------------------------
{
  const char     routine[]  = "Domain::loadbase()";
  const integer  nT         = Geometry::nTotProc();
  const integer  planeSize  = Geometry::planeSize();
  const integer  MAX        = (int) Femlib::value ("MAX_BASE_F");
  char       basefile[StrMax];
  integer    dump_step, i, x, y;
  real       first_time = 0.0, dump_time, last_time = -1.0;
  vector<AuxField*>     PBF ;    // Periodic base file fields
  vector<real*>         PBFdat;  // data storage area for periodic basefields.
  vector<Field*>        U_temp;  // tempory storage for dump reads.
  vector<real*>         U_alloc; // memory for U_temp
  const integer ntot = Geometry::nTotProc();


  U_temp.setSize(2);
  U_alloc.setSize(2);

  for(i = 0; i < 2; i++){
    U_alloc[i] = new real [(size_t) ntot];
    U_temp[i] = new Field(b[i],U_alloc[i],1,elmt,'U');
  }

  PBF.setSize(2*MAX);
  PBFdat.setSize(2*MAX);
  d_t = 0.0;


  cout << "-- Base Field condition    : ";
  ifstream file (strcat (strcpy (basefile, name), ".bse"));

  // read in data from basefile until EOF
  // base file format is sequence of 2D dns field dumps
  // store fields in temp variables PBF,
 
  if (basefile) {
    cout << "read from file " << basefile << endl;    

    n_basefiles = 0;

    cout << tab << "File" << tab   << "Memory" << tab << "Time" 
	 << tab << tab    << "d_t" << endl;

    while (loadfield (file, &dump_step, &dump_time, 
		      basefield, U_temp, BASE_LOAD)) {     

      if (n_basefiles == MAX)
	message(routine, "Limit of base file storage reached", ERROR);

      cout << tab << n_basefiles << tab;

      x = 0 + 2*n_basefiles;
      y = 1 + 2*n_basefiles;

      // PBF memory allocation.
      PBFdat[x] = new real[(size_t) nT];
      PBFdat[y] = new real[(size_t) nT];

      PBF[x] = new AuxField(PBFdat[x], 1, elmt, 'U');
      PBF[y] = new AuxField(PBFdat[y], 1, elmt, 'V');

      cout << ".." << tab;

      // copy auxfields to PBF

      *PBF[x] = *U_temp[0];
      *PBF[y] = *U_temp[1];
		    
      // evaluate d_t
      if  (!n_basefiles) first_time = dump_time;   
      cout << dump_time-first_time << tab;

      if (n_basefiles > 0) {
	d_t = dump_time - last_time;
      }
      last_time = dump_time;

      cout << tab << d_t << endl;
  
      // increment number of base files.
      n_basefiles++;
    
    }

    // accuracy of time is too low -> read from session tokens.
    d_t = (real) Femlib::value ("BASE_DT");
    cout << endl << tab << "BASE_DT: d_t = " << d_t << endl;

    file.close();

    for (i=0;i < 2;i++){
      delete U_temp[i];
      delete U_alloc[i];
    }

    // if only one base file -> no periodic files.
    if (n_basefiles == 1) {
      // only 1 base field - set U[0,1] = PBF[0,1]
      for(i = 0 ; i < 2; i++) { 
      	U[i] = PBF[i];
      }     

    } else { // more than one base field set.

      // check that an even number of fields was given -> FFT limitation.
      
      if (n_basefiles% 2) 
	message(routine, "require even # of base fields", ERROR);
      
      // allocate memory for fourier transformed data
      BaseUV.setSize(2);       // U & V real*
      BaseUV[0] = new real[n_basefiles*planeSize];
      BaseUV[1] = new real[n_basefiles*planeSize];

      // zero new arrays -> if data odd, last entry = 0.0
      Veclib::zero (planeSize*n_basefiles, BaseUV[0],  1);
      Veclib::zero (planeSize*n_basefiles, BaseUV[1],  1);

      // allocate memory for U[0] and U[1] & copy PBF[0,1] across
      Udat.setSize(2);
      Udat[0] = new real[(size_t) nT];
      Udat[1] = new real[(size_t) nT];

      U[0] = new AuxField(Udat[0], 1, elmt, 'U');
      U[1] = new AuxField(Udat[1], 1, elmt, 'V');
      *U[0] = *PBF[0];
      *U[1] = *PBF[1];

      // copy across PBF auxfields to BaseUV & delete PBF
      for(int n = 0; n < n_basefiles; n++) { // was n++
	PBF[2*n]   -> getPlane(0, BaseUV[0] + n*planeSize);
	PBF[2*n+1] -> getPlane(0, BaseUV[1] + n*planeSize);
		 
	// PBF[n] -> getVel( (n/2)*planeSize, BaseUV[n%2] );
	// delete PBF[n];
      }

      // transform data from real to complex.
      Femlib::DFTr (BaseUV[0], n_basefiles, planeSize, +1);
      Femlib::DFTr (BaseUV[1], n_basefiles, planeSize, +1);      


      // scale fields > 0,1 now to avoid duplication.
      
      for (int i = 2; i < n_basefiles; i++ ) {
	Blas::scal (planeSize, 2.0, BaseUV[0] + i*planeSize, 1);
	Blas::scal (planeSize, 2.0, BaseUV[1] + i*planeSize, 1);	
      }
      
      // set period time.
      period = d_t * n_basefiles;
      cout << "   - " << "Period time = "<< n_basefiles << " x " << d_t
	   << " = " << period << endl;

    }
    
  } else {
    message (routine, "Error opening base field file", ERROR);
  }

  cout << endl;
  
}

void Domain::Base_update()
  // -------------------------------------------------------------------
  // update fields U and V with probe data from BaseUV at time t.
  // -------------------------------------------------------------------

{
  
  U[0] -> update(n_basefiles, BaseUV[0], time, period);
  U[1] -> update(n_basefiles, BaseUV[1], time, period);

}


void Domain::dump ()
// ---------------------------------------------------------------------------
// Check if a field-file write is required, carry out.
//
// Fields are inverse Fourier transformed prior to dumping in order to
// provide physical space values.
// ---------------------------------------------------------------------------
{
  const integer periodic = !(step %  (integer) Femlib::value ("IO_FLD"));
  const integer initial  =   step == (integer) Femlib::value ("IO_FLD");
  const integer final    =   step == (integer) Femlib::value ("N_STEP");

  if (!(periodic || final)) return;
  ofstream output;
  
  ROOTONLY {
    const char    routine[] = "Domain::dump";
    char          dumpfl[StrMax], backup[StrMax], command[StrMax];
    const integer verbose   = (integer) Femlib::value ("VERBOSE");
    const integer chkpoint  = (integer) Femlib::value ("CHKPOINT");

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
    if (verbose) message (routine, ": writing field dump", REMARK);
  }

  this -> transform (INVERSE);
  output << *this;
  this -> transform (FORWARD);

  ROOTONLY output.close();
}


void Domain::transform (const integer sign)
// ---------------------------------------------------------------------------
// Fourier transform all Fields according to sign.
// ---------------------------------------------------------------------------
{
  integer       i;
  const integer N = this -> nField ();
  
  for (i = 0; i < N; i++) u[i] -> transform (sign);
}


ofstream& operator << (ofstream& strm,
		       Domain&   D   )
// ---------------------------------------------------------------------------
// Output all Domain field variables on ostream in prism-compatible
// form.  Binary output only.  Note that output is only done on root
// processor.
// ---------------------------------------------------------------------------
{
  int               i;
  const int         N = D.u.getSize();
  vector<AuxField*> field (N);

  for (i = 0; i < N; i++) field[i] = D.u[i];

  writeField (strm, D.name, D.step, D.time, field);

  return strm;
}

integer Domain::loadfield(ifstream&       strm      ,
			  integer*        strm_step ,
			  real*           strm_time ,
			  char*           strm_field,
			  vector<Field*>  strm_u    ,
			  const LoadKind  method    )
{
  const char     routine[] = "Domain::loadField";
  char           s[StrMax], f[StrMax], fields[StrMax];
  integer        npchk,  nzchk, nelchk;
  integer        i, j, np, nz, nel, ntot, nfields;
  integer        swap = 0;


  if (strm.getline(s, StrMax).eof()) return 0;

  strm.getline(s,StrMax).getline(s,StrMax);
  
  u[0] -> describe (f);

  istrstream (s, strlen (s)) >> np    >> np    >> nz    >> nel;
  istrstream (f, strlen (f)) >> npchk >> npchk >> nzchk >> nelchk;
  
  if (np  != npchk ) message (routine, "element size mismatch",       ERROR);
  if (nz  != nzchk ) message (routine, "number of z planes mismatch", ERROR);
  if (nel != nelchk) message (routine, "number of elements mismatch", ERROR);
  
  ntot = np * np * nz * nel;
  if (ntot != Geometry::nTot())
    message (routine, "declared sizes mismatch", ERROR);

  strm.getline(s,StrMax);
  istrstream (s, strlen (s)) >> *strm_step;

  strm.getline(s,StrMax);
  istrstream (s, strlen (s)) >> *strm_time;
  
  strm.getline(s,StrMax).getline(s,StrMax);
  strm.getline(s,StrMax).getline(s,StrMax);

  nfields = 0;
  while (isalpha (s[nfields])) {
    fields[nfields] = s[nfields];
    nfields++;
  }
  fields[nfields] = '\0';

#ifndef STABILITY
  char err[StrMax];
  if (nfields != (int) strlen (strm_field)) {
    sprintf (err, "file: %1d fields, Domain: %1d", nfields,
	     strlen(strm_field));
    message (routine, err, ERROR);
  }

  for (i = 0; i < nfields; i++) 
    if (!strchr (strm_field, fields[i])) {
      sprintf (err,"field %c not present in Domain (%s)",fields[i],field);
      message (routine, err, ERROR);
    }
#endif

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

  if (method == BASE_LOAD) {
    // check that field is only 2D
    if (nfields!=3) message (routine, "Base Field number != 3", ERROR);
    
    // allocate auxfield for temporary pressure
    real*      Pmem = new real[(size_t) ntot];
    AuxField*  P = new AuxField(Pmem, 1, elmt, 'P');

    for (i = 0; i < 2; i++) { // read U and V
      strm >> *strm_u[i];
      if (swap) strm_u[i] -> reverse();
    }
    strm >> *P;           // P -> discarded
    
    delete Pmem;
    delete P;
  }
  else { // method = STD_LOAD

    for (j = 0; j < nfields; j++) {
      for (i = 0; i < nfields; i++)
	if (fields[j] == tolower(strm_field[i])) {
	  strm >>  *strm_u[i];
	  if (swap) strm_u[i] -> reverse();
	  break;
	}
    }
  }

  ROOTONLY if (strm.bad())
    message (routine, "failed reading field file", ERROR);

  // succesful field read 
  return 1;

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

  D.loadfield (strm, &D.step, &D.time, D.field, D.u, STD_LOAD);
 
  Femlib::value ("t", D.time);

  return strm;
}

