//////////////////////////////////////////////////////////////////////////////
// addfield.cpp: process semtex/NEKTON-type field files, computing and
// adding vorticity and divergence, rate of strain magnitude, velocity
// gradient discriminant, etc.
//
// Copyright (c) 1998 <--> $Date: 2019/05/30 06:36:12 $, 
//   Hugh Blackburn, Murray Rudman, Jagmohan Singh
//
//
// Usage:
// -----
// addfield [options] -s session session.fld
//   options:
//   -h        ... print this message
//   -q        ... add kinetic energy per unit mass 0.5(u.u) (default)
//   -d        ... add divergence div(u)
//   -v        ... add vorticity w=curl(u)
//   -e        ... add enstrophy 0.5(w.w)
//   -H        ... add helicity 0.5(u.w)
//   -g        ... add strain rate magnitude sqrt(2SijSji)
//   -D        ... add discriminant of velocity gradient tensor
//                 NB: divergence is assumed to be zero.
//   -J        ... add vortex core measure of Jeong & Hussain. (3D only)
//   -a        ... add all fields derived from velocity (above)
//   -f <func> ... add a computed function <func> of x, y, z, t, etc.
//
// Reserved field names used/assumed:
// ---------------------------------
//
// u -- x velocity component (cylindrical: axial)
// v -- y velocity component (cylindrical: radial)
// w -- z velocity component (cylindrical: azimuthal)
// p -- pressure/density
// c -- scalar
//
// The following are reserved names, not used by addfield.
// ------------------------------------------------------
//
// A -- uu covariance
// B -- uv covariance
// C -- vv covariance
// D -- uw covariance
// E -- vw covariance
// F -- ww covariance
// 
// Computed variables:
// ------------------
//
// d -- divergence
// e -- enstrophy 0.5*(r^2 + s^2 + t^2) = 0.5 (omega . omega)
// f -- a computed function of spatial variables
// g -- strain rate magnitude sqrt(2SijSij)
// q -- kinetic energy per unit mass 0.5*(u^2 + v^2 + w^2) = 0.5 (u . u)
// r -- x component vorticity
// s -- y component vorticity
// t -- z component vorticity
// H -- helicity  0.5*(u*r + v*s + w*t) = 0.5 (u . omega) .
// J -- vortex core identification measure, see [2]. 3D only.
// D -- discriminant of velocity gradient tensor, see [1].

//
// NB: product terms -- such as are used to calculate enstrophy,
// helicity, the invariants and discriminant of the velocity gradient
// tensor, and the strain rate magnitude, all computed in physical
// space -- are not dealiased.  Therefore it is advisable to project
// the original field to a greater number of planes (3/2 rule) before
// these terms are calculated, otherwise the products can be quite
// different from what would be expected (especially if N_Z is small,
// say 4). If this is done you need to edit a matching session file
// with the appropriate value of N_Z.
//
// References
//-----------
//
// [1] Chong, Perry & Cantwell (1990) A general classification of
// three-dimensional flow fields, PF(A) 2:765--777; see also Blackburn
// et al. (1996) JFM 310:269--292
//
// [2] Jeong & Hussain (1995) On the identification of a vortex, JFM
// 285:69--94
//
//
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
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: addfield.cpp,v 9.1 2019/05/30 06:36:12 hmb Exp $";

#include <sem.h>
#include <tensorcalcs.h>
#define FLDS_MAX 64 // -- More than we'll ever want.
#define FLAG_MAX 10 // -- NB: FLAG_MAX should tally with the following enum:
enum {
  ENERGY      ,     // -- NB: the placing of ENERGY and FUNCTION in the first
  FUNCTION    ,	    //    two positions is significant: don't break this.
  DIVERGENCE  ,
  ENSTROPHY   ,
  DISCRIMINANT,
  HELICITY    ,
  STRAINRATE  ,
  VORTICITY   ,
  VORTEXCORE
};

static char  prog[] = "addfield";


static void  getargs  (int, char**, char*&, char*&, char*&, bool[]);
static void  getMesh (const char*,vector<Element*>&);
static bool  getDump (ifstream&,map<char, AuxField*>&,vector<Element*>&,char*&);
static bool  doSwap  (const char*);
static char* fieldNames(ifstream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver -- adapted from probe.C.
// ---------------------------------------------------------------------------
{
  char                       err[StrMax]; 
  Geometry::CoordSys         system;
  char                       *session, *dump, *func, *fields;
  ifstream                   file;
  map<char, AuxField*>       input,  addfield;
  vector<AuxField*>          addbuf, outbuf;
  int_t                      i , j, k, p, q, nComponent, nFields;
  int_t                      np, nz, nel, allocSize, NCOM, NDIM, outbuf_len;
  bool                       add[FLAG_MAX], need[FLAG_MAX], gradient;
  FEML*                      F;
  Mesh*                      M;
  BCmgr*                     B;
  Domain*                    D;
  vector<Element*>           elmt;
  AuxField                   *Func, *Vtx, *DivL, *Nrg, *work, *Disc, *Div;
  vector<AuxField*>          velocity;
  vector<vector<AuxField*> > Vij;     // -- Usually computed, for internal use.
  vector<vector<real_t*> >   VijData; // -- For pointwise access in Vij.
  vector<real_t*>            VorData; // -- Ditto in vorticity.
  map<char, AuxField*>::iterator  ki, ko;

  real_t          *DisData, *DivData, *StrData, *VtxData, *HelData, *EnsData;
  real_t          vel[3], vort[3], tensor[9];

  Femlib::initialize (&argc, &argv);
  for (i = 0; i < FLAG_MAX; i++) add [i] = need [i] = false;

  getargs (argc, argv, session, dump, func, add);

  file.open (dump, ios::in);
  if (!file) message (prog, "no field file", ERROR);
  
  // -- Set up domain.

  F      = new FEML (session);
  M      = new Mesh (F);
  nel    = M -> nEl ();  
  np     =  Femlib::ivalue ("N_P");
  nz     =  Femlib::ivalue ("N_Z");
  system = (Femlib::ivalue ("CYLINDRICAL") ) ?
    Geometry::Cylindrical : Geometry::Cartesian;
  Geometry::set (np, nz, nel, system);

  allocSize = Geometry::nTotal();

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, M);
  NDIM = Geometry::nDim();

  if (NDIM == 2) add[VORTEXCORE] = false;

  for (p = 0, i = 0; i < FLAG_MAX; i++) p += (add[i]) ? 1 : 0;
  if  (p == 0) message (prog, "nothing to be done", ERROR);

  // -- Check if we just have the (first two) cases not requiring derivatives.

  for (p = 0, i = 0; i < FLAG_MAX; i++) p += (add[i]) ? (i + 1) : 0;
  if (p <= 3) gradient = false; else gradient = true;

  for (i = 0; i < FLAG_MAX; i++) need[i] = add[i];  
 
  if (gradient) {
    Vij    .resize (3);
    VijData.resize (3);
    for (i = 0; i < 3; i++) {
      Vij    [i].resize (3);
      VijData[i].resize (3);
      for (j = 0; j < 3; j++) {
	VijData[i][j] = new real_t [allocSize];
	Vij    [i][j] = new AuxField (VijData[i][j], nz, elmt);
	*Vij   [i][j] = 0.0;
      }
    }
  }
  
  B = new BCmgr  (F, elmt);
  D = new Domain (F, elmt, B);

  // -- From the requested fields, flag dependencies.
  // -- First, only allow the "coherent structures" measures for flows
  //    that are 3D.
  // -- Need to link addbuf to addfield so we can use writeField.
  //    Maybe in the longer term we should overload writeField.

  addbuf.resize (addfield.size()); i = 0;
  for (map<char,AuxField*>::iterator k = addfield.begin();
       k != addfield.end(); k++, i++) addbuf[i] = k -> second;
    
  // -- Read fieldnames from dump and allocate input.
 
  fields = fieldNames(file);
  if      (strstr (fields, "uvw")) NCOM = 3;
  else if (strstr (fields, "uv"))  NCOM = 2;
  else message (prog, "lacking velocity components: is session valid?", ERROR);
  
  if (nFields == 1)		// -- Scalar original field.
    nComponent = 1;
  else				// -- Original field was vector.
    nComponent = (nFields == 3) ? 2 : 3;
    
  velocity.resize (NCOM);
  
  for (i = 0; i < strlen(fields); i++)
    input[fields[i]] = new AuxField
      (new real_t[allocSize], nz, elmt, fields[i]);
   
  // -- Allocate space for addfield variables.
  
  // -- Fields without dependants.

  if (need[ENERGY])
    addfield['q'] = new AuxField (new real_t[allocSize], nz, elmt, 'q');

  if (need[FUNCTION])
    addfield['f'] = new AuxField (new real_t[allocSize], nz, elmt, 'f');
  
  if (need[DIVERGENCE]) {
    DivData = new real_t [allocSize];
    *(Div =  new AuxField (DivData, nz, elmt, 'd')) = 0.0;
    addfield['d'] = Div;
  }

  if (need[DISCRIMINANT]) {
    DisData = new real_t [allocSize];
    *(Disc = new AuxField (DisData, nz, elmt, 'D')) = 0.0;
    addfield['D'] = Disc;
  }

  if (need[STRAINRATE]) {
    StrData       = new real_t [allocSize];
    addfield['g'] = new AuxField (StrData, nz, elmt, 'g');
  }

  if (need[VORTEXCORE]) {
    VtxData       = new real_t [allocSize];
    addfield['J'] = new AuxField (VtxData, nz, elmt, 'J');
  }

  if (need[VORTICITY])
    if (NCOM == 2) {
      VorData.resize (1);
      VorData[0]    = new real_t [allocSize];
      addfield['t'] = new AuxField (VorData[0], nz, elmt, 't');
    } else {
      VorData  .resize (3);
      for (i = 0; i < 3; i++) {
	VorData[i]      = new real_t [allocSize];
	addfield['r'+i] = new AuxField (VorData[i], nz, elmt, 'r'+i);
      }
    }

  if (need[ENSTROPHY]) {
    EnsData       = new real_t[allocSize];
    addfield['e'] = new AuxField (EnsData, nz, elmt, 'e');
  }
  
  if (need[HELICITY]) {
    HelData       = new real_t [allocSize];
    addfield['H'] = new AuxField (HelData, nz, elmt, 'H');
  }  

  // -- Read from input and calculate addfield.
  
  while (getDump (file, input, elmt, fields)) {
    
    for (i = 0; i < NCOM; i++) velocity[i] = input['u'+i];
    if (need[FUNCTION]) (*addfield['f']) = func;
    if (need[ENERGY]) ((*addfield['q']) .
		       innerProduct (velocity, velocity)) *= 0.5;
   
    if (gradient) {		// -- All other things.

      // -- First make all VG components.
   
      for (i = 0; i < NDIM ; i++){
	for (j = 0; j < NCOM ; j++) {
	  (*Vij[i][j] = *velocity[j]).gradient(i);
	  if (i == 2) (*Vij[i][j] = *velocity[j]) .
			transform(FORWARD).gradient(i).transform(INVERSE); 
	}	
      }

      if (Geometry::cylindrical()) {
	work = new AuxField (new real_t[allocSize],  nz, elmt);
	if (NDIM == 3) for (j = 0; j < NCOM; j++) Vij[2][j] -> divY();
	(*work = *velocity[1]) . divY(); *Vij[2][2] += *work;
#if 1
	if (NCOM == 3) { (*work = *velocity[2]) . divY(); *Vij[1][2] += *work; }
#else
	if (NCOM == 3) { (*work = *velocity[2]) . divY(); *Vij[1][2] -= *work; }
#endif
      }
  
      // -- Loop over every point in the mesh and compute everything
      //    from Vij.  Quite likely this could be made more efficient
      //    but for now simplicity is the aim.

      for (i = 0; i < allocSize; i++) {
	
	for (k = 0, p = 0; p < 3; p++) {
	  for (q = 0; q < 3; q++, k++)
	    tensor [k] = VijData [p][q][i];
	}

	// -- These operations produce a simple scalar result from Vij.

	if (need[DIVERGENCE])   DivData[i] = tensor3::trace      (tensor);
	if (need[ENSTROPHY])    EnsData[i] = tensor3::enstrophy  (tensor);
	if (need[DISCRIMINANT]) DisData[i] = tensor3::discrimi   (tensor);
	if (need[STRAINRATE])   StrData[i] = tensor3::strainrate (tensor);
	if (need[VORTEXCORE])   VtxData[i] = tensor3::lambda2    (tensor);
	
	// -- Vorticity could be considered scalar in 2D.

	if (need[VORTICITY]) {
	  tensor3::vorticity (tensor, vort);
	  if (NCOM == 2) 
	    VorData[0][i] = vort[2];
	  else { 
	    VorData[0][i] = vort[0]; 
	    VorData[1][i] = vort[1]; 
	    VorData[2][i] = vort[2];
	  }
	}

	// -- Helicity requies velocity too.

	if (need[HELICITY]) {
	  vel[0] = velocity[0] -> data()[i];
	  vel[1] = velocity[1] -> data()[i];
	  if (NCOM ==3) {vel[2] = velocity[2] -> data()[i];}
	  else {vel[2] =0.0;}

	  HelData[i] = tensor3::helicity (tensor, vel);
	}
      }
    }
  
    for (map<char,AuxField*>::iterator k = addfield.begin();
         k != addfield.end(); k++, i++)
      D -> u[0] -> smooth(addfield[k-> first]);
  
    outbuf_len = input.size()+addfield.size();

    // -- Find if input already contains some of the variables that are
    //    being calculated accordingly figure out the required size of
    //    outbuf.

    for (map<char,AuxField*>::iterator k = addfield.begin();
	 k != addfield.end(); k++) {
      ki = input.find (k -> first);
      if (ki != input.end()) outbuf_len--;
    }
    outbuf.resize (outbuf_len); i = 0;         
          
    // -- Set outbuf to input, skip the entries already available in
    //    addfield.
  
    for (k = 0; k < strlen(fields); k++) {
      ki = input.find(fields[k]);
      ko = addfield.find(fields[k]);
      if (ko != addfield.end()) {
	sprintf (err,"found field %c in input, overwriting", fields[k]);
	message (prog, err, WARNING);  
      } else
	outbuf[i++] = ki -> second;
    }
  
    // -- Append addfield data into outbuf.
  
    for (map<char,AuxField*>::iterator k = addfield.begin();
	 k != addfield.end(); k++, i++) outbuf[i] = k -> second;

    writeField (cout, session, 0, 0.0, outbuf);
  }
  
  Femlib::finalize();
  return EXIT_SUCCESS;
}

static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session,
		     char*& dump   ,
		     char*& func   ,
		     bool*  flag   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] =
    "Usage: %s [options] -s session dump.fld\n"
    "options:\n"
    "  -h        ... print this message \n"
    "  -q        ... add kinetic energy per unit mass 0.5(u.u) (default)\n"
    "  -d        ... add divergence div(u)\n"
    "  -v        ... add vorticity w=curl(u)\n"
    "  -e        ... add enstrophy 0.5(w.w)\n"
    "  -H        ... add helicity 0.5(u.w) \n"
    "  -g        ... add strain rate magnitude sqrt(2SijSji)\n"
    "  -D        ... add discriminant of velocity gradient tensor\n"
    "                NB: divergence is assumed to be zero. \n"
    "  -J        ... add vortex core measure of Jeong & Hussain (3D only)\n"
    "  -a        ... add all fields derived from velocity (above)\n"
    "  -f <func> ... add a computed function <func> of x, y, z, t, etc.\n";
              
  int_t i, sum = 0;
  char  buf[StrMax];
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 's':
      if (*++argv[0]) session = *argv; else { --argc; session = *++argv; }
      break;
    case 'v': flag[VORTICITY]    = true; break;
    case 'e': flag[ENSTROPHY]    = true; break;
    case 'H': flag[HELICITY]     = true; break;
    case 'd': flag[DIVERGENCE]   = true; break;
    case 'g': flag[STRAINRATE]   = true; break;
    case 'D': flag[DISCRIMINANT] = true; break;
    case 'J': flag[VORTEXCORE]   = true; break;
    case 'q': flag[ENERGY]       = true; break;
    case 'a': flag[0]=true;
      for (i = 2; i < FLAG_MAX ; i++) flag[i]=true;
      break;
    case 'f':
      if (*++argv[0]) func = *argv; else { --argc; func = *++argv; }
      flag[FUNCTION] = true;
      break;
    default: sprintf (buf, usage, prog); cout<<buf; exit(EXIT_FAILURE); break;
    }

  for (i = 0; i < FLAG_MAX; i++) sum += (flag[i]) ? 1 : 0;
  if (!sum) flag[ENERGY] = true;

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
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


static char* fieldNames (ifstream& file)
// ---------------------------------------------------------------------------
// Return string containing single-character names of fields.
// ---------------------------------------------------------------------------
{
  static char fields[StrMax];
  char        buf[StrMax];
  int_t       i;
  
  for (i = 0; i < 8; i++) file.getline (buf, StrMax);
  file >> fields;
  file.clear();
  file.seekg(0);
  return fields;
}

static bool getDump (ifstream&             file  ,
		     map<char, AuxField*>& u     ,
		     vector<Element*>&     elmt  ,
		     char*&                fieldn)
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  const int_t ntot = Geometry::nTotal();
  char        buf[StrMax], fields[StrMax];
  int_t       i, nf, np, nz, nel;
  bool        swab;
  char         err[StrMax];
  map<char, AuxField*>::iterator  k;
  
  if (file.getline(buf, StrMax).eof()) return false;
  
  if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> np >> nz >> nz >> nel;
  file.getline (buf, StrMax);
  
  if (np != Geometry::nP() || nz != Geometry::nZ() || nel != Geometry::nElmt())
    message (prog, "size of dump mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);
  
  if (u.size() != 0) {
    if (strcmp (fieldn, fields) != 0)
      message (prog, "fields mismatch with first dump in file", ERROR);
  }
  
  for (i = 0; i < nf; i++) {
    k = u.find(fields[i]);
    file >> *(k -> second);
    if (swab) k -> second -> reverse();
  }

  return file.good();
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


