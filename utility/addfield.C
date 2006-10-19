//////////////////////////////////////////////////////////////////////////////
// addfield.C: process sem field files, computing and adding vorticity and
// divergence, rate of strain tensor, velocity gradient invariants, etc.
//
// Copyright (c) 1998 <--> $Date$, 
//   Hugh Blackburn, Murray Rudman
//
// NB: the input field file is assumed to contain only velocity and
// pressure data.
//
// Usage:
// -----
// addfield [options] -s session session.fld
//   options:
//   -h        ... print this message
//   -q        ... add kinetic energy per unit mass (default)
//   -v        ... add vorticity
//   -e        ... add enstrophy
//   -B        ... add elliptic instability growth rate
//   -b        ... add ellipticity
//   -h        ... add helicity (3D only)
//   -d        ... add divergence
//   -t        ... add components of rate of strain tensor Sij
//   -g        ... add strain rate magnitude sqrt(2SijSij)
//   -i        ... add discriminant of velocity gradient tensor, 3D only
//   -l        ... add Lamb vector (vorticityxu) and its divergence, 3D only
//   -a        ... add all fields derived from velocity
//   -f <func> ... add a computed function (of x, y, z, t, etc.)
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
// The following are reserved names, not used:
// ------------------------------------------
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
// a -- elliptic instability growth rate
// b -- ellipticity parameter \gamma/\omega = 0.5*G/sqrt(e)
// d -- divergence
// e -- enstrophy 0.5*(r^2 + s^2 + t^2)
// f -- a computed function of spatial variables
// G -- Strain rate magnitude sqrt(2SijSij)
// h -- helicity  0.5*(u*r + v*s + w*t)
// J -- The vortex core identification measure of Jeong and Hussain, JFM 285.
// L -- Discriminant of velocity gradient tensor 27/4 R^2 + Q^3
// m -- x component Lamb vector
// n -- y component Lamb vector
// o -- z component Lamb vector
// Q -- divergence  Lamb vector
// q -- kinetic energy per unit mass 0.5*(u^2 + v^2 + w^2)
// r -- x component vorticity
// s -- y component vorticity
// t -- z component vorticity
//
// NB: product terms -- such as are used to calculate enstrophy,
// helicity, the invariants and discriminant of the velocity gradient
// tensor, and the strain rate magnitude -- are not dealiased.
// Therefore it is advisable to project the original field to a
// greater number of planes (3/2 rule) before these terms are
// calculated, otherwise the products can be quite different from what
// would be expected (especially if N_Z is small, say 4).
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <sem.h>
#include <tensorcalcs.h>

#define FLDS_MAX 64 // -- More than we'll ever want.
#define FLAG_MAX 13 // -- NB: FLAG_MAX should tally with the following enum:
enum {
  ENERGY      ,     // -- NB: the placing of ENERGY and FUNCTION in the first
  FUNCTION    ,	    //    two positions is significant: don't break this.
  EGROWTH     ,
  ELLIPTICITY ,
  ENSTROPHY   ,
  DISCRIMINANT,
  DIVERGENCE  ,
  HELICITY    ,
  LAMBVECTOR  ,
  STRAINRATE  ,
  STRAINTENSOR,
  VORTICITY   ,
  VORTEXCORE
};

static char  prog[] = "addfield";
static void  getargs  (int, char**, char*&, char*&, char*&, bool[]);
static int_t getDump  (Domain*, ifstream&);
static void  putDump  (Domain*, vector<AuxField*>&, int_t, ostream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys         system;
  char                       *session, *dump, *func, fields[StrMax];
  int_t                      i, j, k, p, q;
  int_t                      np, nz, nel, allocSize, nComponent;
  int_t                      iAdd = 0;
  bool                       add[FLAG_MAX], need[FLAG_MAX], gradient;
  ifstream                   file;
  FEML*                      F;
  Mesh*                      M;
  BCmgr*                     B;
  Domain*                    D;
  vector<Element*>           elmt;
  AuxField                   *Ens, *Hel, *Div, *InvQ, *InvR, *Disc, *Strain;
  AuxField                   *Func, *Ell, *Egr, *Vtx, *work, *DivL, *Nrg;
  vector<AuxField*>          lamb, velocity, vorticity, addField(FLDS_MAX);
  vector<vector<AuxField*> > Sij;
  vector<vector<AuxField*> > Vij;     // -- Usually computed, for internal use.
  vector<vector<real_t*> >   VijData; // -- For pointwise access.
  real_t                     *egrow, *VtxData; // -- Ditto.
  real_t                     tensor[9];

  Femlib::initialize (&argc, &argv);
  for (i = 0; i < FLAG_MAX; i++) add [i] = need [i] = false;

  // -- Read command line.

  getargs (argc, argv, session, dump, func, add);

  file.open (dump);
  if (!file) message (prog, "can't open input field file", ERROR);

  // -- Set up domain.

  F      = new FEML (session);
  M      = new Mesh (F);
  nel    = M -> nEl ();  
  np     =  Femlib::ivalue ("N_P");
  nz     =  Femlib::ivalue ("N_Z");
  system = (Femlib::ivalue ("CYLINDRICAL") ) ?
                       Geometry::Cylindrical : Geometry::Cartesian;
  Geometry::set (np, nz, nel, system);
  if   (nz > 1) strcpy (fields, "uvwp");
  else          strcpy (fields, "uvp");
  nComponent = Geometry::nDim();
  allocSize  = Geometry::nTotal();

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, M);

  B = new BCmgr  (F, elmt);
  D = new Domain (F, elmt, B);

  velocity.resize   (nComponent);
  vorticity.resize ((nComponent == 2) ? 1 : 3);
  for (i = 0; i < nComponent; i++) velocity[i] = D -> u[i];

  // -- From the requested fields, flag dependencies.

  if  (nComponent < 3)
    add[HELICITY] = add[DISCRIMINANT] = add[VORTEXCORE] = false;

  for (p = 0, i = 0; i < FLAG_MAX; i++) p += (add[i]) ? 1 : 0;
  if  (p == 0) message (prog, "nothing to be done", ERROR);

  for (p = 0, i = 0; i < FLAG_MAX; i++)
    p += (add[i]) ? (i + 1) : 0;
  if (p <= 3) gradient = false; else gradient = true;

  for (i = 0; i < FLAG_MAX; i++) need[i] = add[i];

  if (add[ENSTROPHY]) {
    need[VORTICITY]    = true;
    need[ENSTROPHY]    = true;
  }
  if (add[HELICITY]) {
    need[VORTICITY]    = true;
  }
  if (add[ELLIPTICITY]) {
    need[VORTICITY]    = true;
    need[ENSTROPHY]    = true;
    need[STRAINRATE]   = true;
    need[STRAINTENSOR] = true;
  }
  if (add[STRAINRATE]) {
    need[STRAINTENSOR] = true;
  }
  if (add[EGROWTH]) {
    need[VORTICITY]    = true;
    need[ENSTROPHY]    = true;
    need[STRAINRATE]   = true;
    need[STRAINTENSOR] = true;
    need[ELLIPTICITY]  = true;
  }
  if (add[LAMBVECTOR]) {
    need[VORTICITY]    = true;
    need[LAMBVECTOR]   = true;
  }

  // -- If any gradients need to be computed, we make all of the.
  
  if (gradient) {
    Vij    .resize (nComponent);
    VijData.resize (nComponent);
    for (i = 0; i < nComponent; i++) {
      Vij    [i].resize (nComponent);
      VijData[i].resize (nComponent);
      for (j = 0; j < nComponent; j++) {
	VijData[i][j] = new real_t [allocSize];
	Vij   [i][j] = new AuxField (VijData[i][j], nz, elmt);
	*Vij  [i][j] = 0.0;
      }
    }
  }
  
  // -- Fields without dependants.

  work = new AuxField (new real_t[allocSize], nz, elmt);

  if (need[ENERGY]) {
    Nrg = new AuxField (new real_t[allocSize], nz, elmt, 'q');
    addField[iAdd++] = Nrg;
  }

  if (need[FUNCTION]) {
    Func = new AuxField (new real_t[allocSize], nz, elmt, 'f');
    addField[iAdd++] = Func;
  }

  if (need[DIVERGENCE]) {
    Div  =  new AuxField (new real_t[allocSize], nz, elmt, 'd');
    *Div = 0.0;
    addField[iAdd++] = Div;
  }

  if (need[DISCRIMINANT]) {
    *(InvQ = new AuxField (new real_t[allocSize], nz, elmt, 'Q')) = 0.0;
    *(InvR = new AuxField (new real_t[allocSize], nz, elmt, 'R')) = 0.0;
    *(Disc = new AuxField (new real_t[allocSize], nz, elmt, 'L')) = 0.0;
    addField[iAdd++] = Disc;
  }

  if (need[VORTEXCORE]) {
    VtxData = new real_t [allocSize];
    Vtx = new AuxField (VtxData, nz, elmt, 'J');
    addField[iAdd++] = Vtx;
  }

  // -- Vorticity and its dependants.

  if (need[VORTICITY])
    if (nComponent == 2) {
      vorticity[0] = new AuxField (new real_t[allocSize], nz, elmt, 't');
      if (add[VORTICITY])
	addField[iAdd++] = vorticity[0];
    } else {
      for (i = 0; i < nComponent; i++)
	vorticity[i] = new AuxField (new real_t[allocSize], nz, elmt, 'r' + i);
      if (add[VORTICITY])
	for (i = 0; i < 3; i++) addField[iAdd++] = vorticity[i];
    }

  if (add[LAMBVECTOR]) {
    if (nComponent == 2) {
      lamb.resize(2);
      for (i = 0; i < nComponent; i++)
	lamb[i] = new AuxField (new real_t[allocSize], nz, elmt, 'm' + i);
      for (i = 0; i < 2; i++) addField[iAdd++] = lamb[i];
    } else {
      lamb.resize(3);
      for (i = 0; i < nComponent; i++)
	lamb[i] = new AuxField (new real_t[allocSize], nz, elmt, 'm' + i);
      for (i = 0; i < 3; i++) addField[iAdd++] = lamb[i];
    }
    DivL = new AuxField (new real_t[allocSize], nz, elmt, 'Q');
    addField[iAdd++] = DivL;
  }

  if (need[ENSTROPHY]) {
    Ens = new AuxField (new real_t[allocSize], nz, elmt, 'e');
    if (add[ENSTROPHY]) addField[iAdd++] = Ens;
  }
  
  if (need[HELICITY]) {
    Hel = new AuxField (new real_t[allocSize], nz, elmt, 'h');
    if (add[HELICITY]) addField[iAdd++] = Hel;
  }

  if (need[ELLIPTICITY]) {
    Ell = new AuxField (new real_t[allocSize], nz, elmt, 'b');
    if (add[ELLIPTICITY]) addField[iAdd++] = Ell;
  }

  // -- Rate of strain tensor and its dependants.

  if (need[STRAINTENSOR]) {
    Sij.resize (nComponent);
    for (k = 0, i = 0; i < nComponent; i++) {
      Sij[i].resize (nComponent);
      for (j = 0; j < nComponent; j++) {
	if (j >= i) {
	  Sij[i][j] = new AuxField (new real_t[allocSize], nz, elmt, 'i'+k++);
	} else
	  Sij[i][j] = Sij[j][i];
      }
    }
  }

  if (need[STRAINRATE]) {
    *(Strain = new AuxField (new real_t[allocSize], nz, elmt, 'G')) = 0.0;
    if (add[STRAINRATE]) addField[iAdd++] = Strain;
  }

  if (need[EGROWTH]) {
    egrow = new real_t[allocSize]; // -- A handle for direct access to data.
    Egr = new AuxField (egrow, nz, elmt, 'a');
    if (add[EGROWTH]) addField[iAdd++] = Egr;
  }
  
  // -- Cycle through field dump, first computing the velocity
  //    gradient tensor and then the needed quantities -- vorticity
  //    etc.  The order of computation is determined by dependencies.
  //    Then write requested output, listed in addField.
  
  while (getDump (D, file)) {
        
    if (need[FUNCTION]) *Func = func;

    // -- Velocity gradient tensor; transform required for z (2) gradient.  

    if (nComponent > 2) D -> transform (FORWARD);

    if (gradient)
      for (i = 0; i < nComponent ; i++)
	for (j = 0; j < nComponent ; j++) {
	  (*Vij[i][j] = *velocity[i]) . gradient (j);
	  Vij[i][j] -> transform (INVERSE);
	}

    if (need[ENERGY]) {
      *Nrg = 0.0;
      for (i = 0; i < nComponent ; i++)
	Nrg -> timesPlus (*D -> u[i], *D -> u[i]);
      *Nrg *= 0.5;
    }

    if (nComponent > 2) D -> transform (INVERSE);

    if (Geometry::cylindrical() && nComponent == 3) {
      for (i = 0; i < nComponent; i++) Vij[i][2] -> divY();
      (*work = *D -> u[2]) . divY();  *Vij[1][2] -= *work;
      (*work = *D -> u[1]) . divY();  *Vij[2][2] += *work;
    }
    
    if (need[DISCRIMINANT]) {

      // -- 2nd invariant (Q from Chong et al.).

      InvQ -> times      (*Vij[0][0], *Vij[1][1]);
      InvQ -> timesMinus (*Vij[0][1], *Vij[1][0]);
      InvQ -> timesPlus  (*Vij[0][0], *Vij[2][2]);
      InvQ -> timesMinus (*Vij[0][2], *Vij[2][0]);
      InvQ -> timesPlus  (*Vij[1][1], *Vij[2][2]);
      InvQ -> timesMinus (*Vij[1][2], *Vij[2][1]);

      // -- 3rd invariant: determinant of Vij (R from Chong et al.).

      work -> times      (*Vij[1][1], *Vij[2][2]);
      work -> timesMinus (*Vij[2][1], *Vij[1][2]);
      InvR -> times      (*work,      *Vij[0][0]);

      work -> times      (*Vij[1][2], *Vij[2][0]);
      work -> timesMinus (*Vij[2][2], *Vij[1][0]);
      InvR -> timesPlus  (*work,      *Vij[0][1]);

      work -> times      (*Vij[2][1], *Vij[1][0]);
      work -> timesMinus (*Vij[1][1], *Vij[2][0]);
      InvR -> timesPlus  (*work,      *Vij[0][2]);

      // -- Discriminant L of Vij.
      //    NB: DIVERGENCE (P from Chong et al.) ASSUMED = 0.

      work -> times (*InvQ, *InvQ);
      Disc -> times (*work, *InvQ);
      work -> times (*InvR, *InvR);
      *work *= 6.75;
      *Disc += *work;
    }
	
    if (need[DIVERGENCE])
      for (i = 0; i < nComponent; i++) *Div += *Vij[i][i];
    
    if (need[VORTICITY]) {
      if (nComponent == 2) {
	*vorticity[0]  = *Vij[1][0];
	*vorticity[0] -= *Vij[0][1];
      } else {
	*vorticity[0]  = *Vij[2][1];
	*vorticity[0] -= *Vij[1][2];
	*vorticity[1]  = *Vij[0][2];
	*vorticity[1] -= *Vij[2][0];
	*vorticity[2]  = *Vij[1][0];
 	*vorticity[2] -= *Vij[0][1];
      }
    }

    if (need[LAMBVECTOR]) {
      if (nComponent == 2) {
	*lamb[0] = 0.0;
	lamb[0] -> timesMinus (*D -> u[1], *vorticity[0]);
	lamb[1] -> times      (*D -> u[0], *vorticity[0]);
	*DivL  = (*work = *lamb[0]) . gradient (0);
	*DivL += (*work = *lamb[1]) . gradient (1);
	if (Geometry::cylindrical())
	  *DivL += (*work = *lamb[1]) . divY();
      } else {
	lamb[0] -> times      (*D -> u[2], *vorticity[1]);
	lamb[0] -> timesMinus (*D -> u[1], *vorticity[2]);
	lamb[1] -> times      (*D -> u[0], *vorticity[2]);
	lamb[1] -> timesMinus (*D -> u[2], *vorticity[0]);
	lamb[2] -> times      (*D -> u[1], *vorticity[0]);
	lamb[2] -> timesMinus (*D -> u[0], *vorticity[1]);
	*DivL  = (*work = *lamb[0]) . gradient (0);
	*DivL += (*work = *lamb[1]) . gradient (1);
	(*work = *lamb[2]).transform(FORWARD).gradient(2).transform(INVERSE);
	if (Geometry::cylindrical()) {
	  *DivL += work -> divY();
	  *DivL += (*work = *lamb[1]) . divY();
	} else
	  *DivL += *work;
      }
    }
    
    if (need[ENSTROPHY])
      Ens -> innerProduct (vorticity, vorticity) *= 0.5;

    if (need[HELICITY])
      Hel -> innerProduct (vorticity, velocity)  *= 0.5;

    if (need[STRAINTENSOR]) {
      for (i = 0; i < nComponent; i++)
	for (j = i; j < nComponent; j++) {
	  *Sij[i][j]  = *Vij[i][j];
	  *Sij[i][j] += *Vij[j][i];
	  *Sij[i][j] *= 0.5;
	}
    }

    if (need[STRAINRATE]) {
      *Strain = 0.0;
      for (i = 0; i < nComponent; i++)
	for (j = 0; j < nComponent; j++)
	  Strain -> timesPlus (*Sij[i][j], *Sij[j][i]);
      (*Strain *= 2.0) . sqroot();
    }

    if (need[ELLIPTICITY]) {
      ((*work = *Ens) . sqroot() *= 2.0) += 1.0e-1;
      Ell -> divide (*Strain, *work);
    }

    if (need[EGROWTH]) {
      *Egr = *Ell;
      Veclib::clip (allocSize, 0.0, 1.0, egrow, 1, egrow, 1);
      Veclib::spow (allocSize, 2.811,    egrow, 1, egrow, 1);
      Veclib::vneg (allocSize,           egrow, 1, egrow, 1);
      Veclib::sadd (allocSize, 1.0,      egrow, 1, egrow, 1);
      Veclib::spow (allocSize, 0.3914,   egrow, 1, egrow, 1);
      Veclib::smul (allocSize, 9.0/16.0, egrow, 1, egrow, 1);
      (*work = *Strain) *= sqrt (1.0/8.0);
      Egr -> times (*Egr, *work);
    }

    if (need[VORTEXCORE])	// -- Only done for 3-component fields.
      for (i = 0; i < allocSize; i++)
	for (k = 0, p = 0; p < 3; p++)
	  for (q = 0; q < 3; q++, k++) {
	    tensor [k] = VijData[p][q][i];
	    VtxData[i] = lambda2 (tensor);
	  }

    // -- Finally, add mass-projection smoothing on everything.

    for (i = 0; i < iAdd; i++) D -> u[0] -> smooth (addField[i]);

    putDump (D, addField, iAdd, cout);
  }
  
  file.close();
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
    "  -q        ... add kinetic energy per unit mass (default)\n"
    "  -v        ... add vorticity\n"
    "  -e        ... add enstrophy\n"
    "  -H        ... add helicity (3D only)\n"
    "  -b        ... add ellipticity\n"
    "  -B        ... add elliptic instability growth rate\n"
    "  -d        ... add divergence\n"
    "  -l        ... add Lamb vector (uxw) and its divergence\n"
    "  -g        ... add strain rate magnitude sqrt(2SijSij)\n"
    "  -i        ... add discriminant of velocity gradient tensor\n"
    "                NB: divergence is assumed to be zero. (3D only)\n"
    "  -j        ... add vortex core measure of Jeong & Hussain\n"
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
    case 'B': flag[EGROWTH]      = true; break;
    case 'b': flag[ELLIPTICITY]  = true; break;
    case 'd': flag[DIVERGENCE]   = true; break;
    case 'g': flag[STRAINRATE]   = true; break;
    case 'i': flag[DISCRIMINANT] = true; break;
    case 'j': flag[VORTEXCORE]   = true; break;
    case 'l': flag[LAMBVECTOR]   = true; break;
    case 'q': flag[ENERGY]       = true; break;
    case 'a': for (i = 0; i < FLAG_MAX; i++) flag[i] = true; break;
    case 'f':
      if (*++argv[0]) func = *argv; else { --argc; func = *++argv; }
      flag[FUNCTION] = true; break;
    default: sprintf (buf, usage, prog); cout<<buf; exit(EXIT_FAILURE); break;
    }

  for (i = 0; i < FLAG_MAX; i++) sum += (flag[i]) ? 1 : 0;
  if (!sum) flag[ENERGY] = true;

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static int_t getDump (Domain*   D   ,
		      ifstream& dump)
// ---------------------------------------------------------------------------
// Read next set of field dumps from file.
// ---------------------------------------------------------------------------
{
  dump >> *D;
  return dump.good ();
}


static void putDump  (Domain*            D       ,
		      vector<AuxField*>& addField,
		      int_t              nOut    ,
		      ostream&           strm    )
// ---------------------------------------------------------------------------
// This is a version of the normal Domain dump that adds extra AuxFields.
// ---------------------------------------------------------------------------
{
  char *hdr_fmt[] = { 
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

  int_t  i, nComponent;
  char   routine[] = "putDump";
  char   s1[StrMax], s2[StrMax];
  time_t tp (::time (0));

  if (D -> nField() == 1)	// -- Scalar original field.
    nComponent = 1;
  else				// -- Original field was vector.
    nComponent = (D -> nField() == 3) ? 2 : 3;

  sprintf (s1, hdr_fmt[0], D -> name);
  strm << s1;

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  strm << s1;

  D -> u[0] -> describe (s2);
  sprintf (s1, hdr_fmt[2], s2);
  strm << s1;

  sprintf (s1, hdr_fmt[3], D -> step);
  strm << s1;

  sprintf (s1, hdr_fmt[4], D -> time);
  strm << s1;

  sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
  strm << s1;

  sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
  strm << s1;

  sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
  strm << s1;

  for (i = 0; i <= nComponent; i++) s2[i] = D -> u[i] -> name();
  for (i = 0; i <  nOut; i++)
    s2[nComponent + i + 1] = addField[i] -> name();
  s2[nComponent + nOut + 1] = '\0';

  sprintf (s1, hdr_fmt[8], s2);
  strm << s1;

  sprintf (s2, "binary ");
  Veclib::describeFormat (s2 + strlen (s2));
  sprintf (s1, hdr_fmt[9], s2);
  strm << s1;

  for (i = 0; i <= nComponent; i++) strm << *D -> u[i];
  for (i = 0; i < nOut; i++) strm << *addField[i];

  if (!strm) message (routine, "failed writing field file", ERROR);
  strm << flush;
}
