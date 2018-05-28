///////////////////////////////////////////////////////////////////////////////
// nonlin.cpp: given a field file containing velocity data u, compute
// and output the nonlinear terms u.grad(u), in physical space.
//
// Copyright (c) 2016 <--> $Date$,
//   Hugh Blackburn
//
// USAGE
// -----
// nonlin [options] session [file]
// options:
// -h       ... print this message.
//
// The nonlinear terms are output as a standard semtex field file,
// where the names uv(w) now stand for components of the nonlinear
// terms associated with the xy(z) directions.  Any extra fields in
// the input file are ignored (and discarded).
//
// The computation of nonlinear terms u.grad(u) is as stated: in
// non-conservative (a.k.a. convective) form.  This differs from the
// default used by dns, which is the (alternating) skew symmetric
// form.
//
// i.e., in Cartesian component form
//
//           n  =   u  d(u ) / dx 
//            i      j    i      j
//
// in cylindrical coordinates
//
//           nx = {ud(u)/dx + vd(u)/dy + 1/y [wd(u)/dz]}
//           ny = {ud(v)/dx + vd(v)/dy + 1/y [wd(v)/dz - ww]}
//           nz = {ud(w)/dx + vd(w)/dy + 1/y [wd(w)/dz + wv]}
//
// The nonlinear terms are NOT mass-matrix smoothed (averaged) along
// element boundaries, as is done in dns.
// 
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA
///////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>

static char  prog[] = "nonlin";
static void  getargs  (int, char**, char*&, istream*&);
static void  getMesh  (const char*, vector<Element*>&);
static void  getDump  (istream&, vector<Element*>&,
		       vector<AuxField*>&, vector<AuxField*>&,
		       AuxField*&);
static bool  doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys         system;
  char*                      session;
  istream*                   file;
  vector<Element*>           elmt;
  vector<AuxField*>          u, n;     // -- Velocity fields, nonlinear terms.
  AuxField*                  work;
  int_t                      i, j, NDIM, NCOM;

  Femlib::initialize (&argc, &argv);

  // -- Read command line.

  getargs (argc, argv, session, file);

  // -- Set up geometrical information.

  getMesh (session, elmt);

  NDIM = Geometry::nDim();

  getDump (*file, elmt, u, n, work);

  NCOM = u.size();

  // -- Now manufacture the nonlinear terms.

  if (Geometry::cylindrical()) {

    for (i = 0; i < NCOM; i++) {

      // -- Terms involving azimuthal derivatives and frame components.

      if (NCOM == 3) {
	if (i == 1) n[1] -> timesMinus (*u[2], *u[2]);
	if (i == 2) n[2] -> timesPlus  (*u[2], *u[1]);

	if (NDIM > 2) {
	  (*work = *u[i]) . 
	    transform (FORWARD) . gradient (2) . transform (INVERSE);
	  n[i] -> timesPlus (*u[2], *work);
	}
      }

      n[i] -> divY ();

      // -- 2D convective derivatives.

      for (j = 0; j < 2; j++) {
	(*work = *u[i]) . gradient (j);
	n[i] -> timesPlus (*u[j], *work);
      }

//      work -> smooth (n[i]);
    }

  } else {			// -- Cartesian coordinates.

    for (i = 0; i < NCOM; i++) {
      for (j = 0; j < NDIM; j++) { // -- Perform n_i += u_j d(u_i) / dx_j.

	if (j == 2) (*work = *u[i]) . 
		      transform (FORWARD). gradient (j) . transform (INVERSE);
	else    (*work = *u[i]) . gradient (j);
	n[i] -> timesPlus (*u[j], *work);
      }
//      work -> smooth (n[i]);
    }
  }

  writeField (cout, session, 0, 0.0, n);
    
  return EXIT_SUCCESS;
}


static void getargs (int        argc   ,
		     char**     argv   ,
		     char*&     session,
		     istream*&  file   )
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session [file]"
    "  [options]:\n"
    "  -h       ... print this message\n";

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  switch (argc) {
  case 1:
    session = argv[0];
    file = &cin;
    break;
  case 2:
    session = argv[0];
    file = new ifstream (argv[1]);
    if (file -> fail()) {
      cerr << usage;
      sprintf (buf, "unable to open field file: %s", argv[1]);
      message (prog, buf, ERROR);
    }
    break;
  default:
    cerr << usage;
    message (prog, "session file not supplied", ERROR);
    break;
  }  
}


static void getMesh (const char*       session,
		     vector<Element*>& elmt   )
// ---------------------------------------------------------------------------
// Set up 2D mesh information. Note that parser tokens and Geometry
// are set here, too.
// ---------------------------------------------------------------------------
{
  FEML* F = new FEML (session);
  Mesh* M = new Mesh (F);
  
  const int_t nel = M -> nEl();  
  const int_t np  = Femlib::ivalue ("N_P");
  const int_t nz  = Femlib::ivalue ("N_Z");

  Geometry::CoordSys space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);
  elmt.resize   (nel);

  for (int_t k = 0; k < nel; k++) elmt[k] = new Element (k, np, M);
}


static void getDump (istream&           file,
		     vector<Element*>&  elmt,
		     vector<AuxField*>& u   ,
		     vector<AuxField*>& n   ,
		     AuxField*&         work)
// ---------------------------------------------------------------------------
// Load velocity data from field dump into 'u', with byte-swapping if
// required.
// 
// Also create (and clear) matching array of nonlinear data areas 'n',
// and an AuxField to be used for workspace.
// ---------------------------------------------------------------------------
{
  char    buf[StrMax], fields[StrMax];
  int_t   i, np, nz, nel, nf, ncom;
  bool    swab;
  real_t* alloc;

  file.getline (buf, StrMax);

  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> np >> nz >> nz >> nel;
  file.getline (buf, StrMax);
  
  if (np != Geometry::nP() || nz != Geometry::nZ() || nel != Geometry::nElmt())
    message (prog, "input file mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  if      (strstr (fields, "uvw")) ncom = 3;
  else if (strstr (fields, "uv"))  ncom = 2;
  else message (prog, "input file lacks velocity components", ERROR);

  u.resize (ncom);
  n.resize (ncom);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields, read in velocity data.

  for (i = 0; i < ncom; i++) {
    alloc = new real_t [Geometry::nTotProc()];
    u[i]  = new AuxField (alloc, nz, elmt, 'u' + i);

    alloc = new real_t [Geometry::nTotProc()];
    n[i]  = new AuxField (alloc, nz, elmt, 'u' + i);
  }
  alloc = new real_t [Geometry::nTotProc()];
  work  = new AuxField (alloc, nz, elmt, 0);

  // -- Read binary field data. Use work as a temporary staging location.

  for (i = 0; i < nf; i++) {
    file >> *work;
    if (fields[i] == 'u') { 
      if (swab) work -> reverse(); *u[0] = *work;
    } else if (fields[i] == 'v') {
      if (swab) work -> reverse(); *u[1] = *work;
    } else if (fields[i] == 'w') {
      if (swab) work -> reverse(); *u[2] = *work;
    }
  }

  if (file.bad()) message (prog, "problem reading input data", ERROR);

  // -- Clear all nonlinear storage areas.

  for (i = 0; i < ncom; i++) *n[i]  = 0.0;
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping is required to make sense of binary input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    message (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "endian"))
    message (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}
