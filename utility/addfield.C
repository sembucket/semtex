//////////////////////////////////////////////////////////////////////////////
// addfield.C: process sem field files, computing and adding vorticity and
// divergence, rate of strain tensor, velocity gradient invariants, etc.
//
// Copyright (c) 1998,2004 Hugh Blackburn, Murray Rudman
//
// NB: the input field file is assumed to contain only velocity and
// pressure data.
//
// Usage:
// -----
// addfield [options] -s session session.fld
//   options:
//   -h        ... print this message
//   -v        ... add vorticity (default)
//   -e        ... add enstrophy
//   -B        ... add elliptic instability growth rate
//   -b        ... add ellipticity
//   -h        ... add helicity (3D only)
//   -d        ... add divergence
//   -t        ... add components of rate of strain tensor Sij
//   -g        ... add strain rate magnitude sqrt(2SijSij)
//   -i        ... add invariants and discriminant of Vij (NOT Vij itself)
//                 (NB: Divergence is ASSUMED equal to zero.) 3D only
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
// i -- uu strain rate component
// j -- uv strain rate component
// k -- vv strain rate component
// L -- Discriminant  of velocity gradient tensor 27/4 R^2 + Q^3
// l -- uw strain rate component
// m -- vw strain rate component
// n -- ww strain rate component
// Q -- 2nd invariant of velocity gradient tensor
// R -- 3rd invariant of velocity gradient tensor
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
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <sem.h>
#include <ctime>

#define FLDS_MAX 64 // -- More than we'll ever want.
#define FLAG_MAX 10 // -- NB: FLAG_MAX should tally with the following enum:
enum {
  EGROWTH     ,
  ELLIPTICITY ,
  ENSTROPHY   ,
  DIVERGENCE  ,
  FUNCTION    ,
  HELICITY    ,
  INVARIANTS  ,
  STRAINRATE  ,
  STRAINTENSOR,
  VORTICITY
};

static char prog[] = "addfield";
static void getargs  (int, char**, char*&, char*&, char*&, int[]);
static int  getDump  (Domain*, ifstream&);
static void putDump  (Domain*, vector<AuxField*>&, int, ostream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys system;
  char               *session, *dump, *func, fields[StrMax];
  int                i, j, k, np, nz, nel, allocSize, nComponent, iAdd = 0;
  int                add[FLAG_MAX], need[FLAG_MAX];
  ifstream           file;
  FEML*              F;
  Mesh*              M;
  BCmgr*             B;
  Domain*            D;
  vector<Element*>   elmt;
  real*              egrow;
  AuxField           *Ens, *Hel, *Div, *InvQ, *InvR, *Disc, *Strain;
  AuxField           *Func, *Ell, *Egr, *work;
  vector<AuxField*>  velocity, vorticity, addField(FLDS_MAX);
  AuxField***        Sij;
  AuxField***        Vij; // -- Always computed, for internal use.

  Femlib::initialize (&argc, &argv);
  Veclib::zero (FLAG_MAX, add,  1);  // -- Requested fields.
  Veclib::zero (FLAG_MAX, need, 1);  // -- Requests + dependencies.

  // -- Read command line.

  getargs (argc, argv, session, dump, func, add);

  file.open (dump);
  if (!file) message (prog, "can't open input field file", ERROR);

  // -- Set up domain.

  F      = new FEML (session);
  M      = new Mesh (F);
  nel    = M -> nEl ();  
  np     =  Femlib::ivalue ("N_POLY");
  nz     =  Femlib::ivalue ("N_Z"   );
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
    add[HELICITY] = add[INVARIANTS] = 0;
  if (Veclib::sum (FLAG_MAX, add, 1) == 0)
    message (prog, "nothing to be done", ERROR);
  for (i = 0; i < FLAG_MAX; i++) need[i] = add[i];

  if (add[ENSTROPHY]) {
    need[VORTICITY]    = 1;
    need[ENSTROPHY]    = 1;
  }
  if (add[HELICITY]) {
    need[VORTICITY]    = 1;
  }
  if (add[ELLIPTICITY]) {
    need[VORTICITY]    = 1;
    need[ENSTROPHY]    = 1;
    need[STRAINRATE]   = 1;
    need[STRAINTENSOR] = 1;
  }
  if (add[STRAINRATE]) {
    need[STRAINTENSOR] = 1;
  }
  if (add[EGROWTH]) {
    need[VORTICITY]    = 1;
    need[ENSTROPHY]    = 1;
    need[STRAINRATE]   = 1;
    need[STRAINTENSOR] = 1;
    need[ELLIPTICITY]  = 1;
  }

  // -- Always compute the velocity gradient tensor for internal use.

  Vij = new AuxField** [nComponent];
  for (i = 0; i < nComponent; i++) {
    Vij[i] = new AuxField* [nComponent];
    for (j = 0; j < nComponent; j++) {
       Vij[i][j] = new AuxField (new real[allocSize], nz, elmt);
      *Vij[i][j] = 0.0;
    }
  }
  
  // -- Fields without dependants.

  work = new AuxField (new real[allocSize], nz, elmt);

  if (need[FUNCTION]) {
    Func = new AuxField (new real[allocSize], nz, elmt, 'f');
    addField[iAdd++] = Func;
  }

  if (need[DIVERGENCE]) {
    Div  =  new AuxField (new real[allocSize], nz, elmt, 'd');
    *Div = 0.0;
    addField[iAdd++] = Div;
  }

  if (need[INVARIANTS]) {
    *(InvQ = new AuxField (new real[allocSize], nz, elmt, 'Q')) = 0.0;
    *(InvR = new AuxField (new real[allocSize], nz, elmt, 'R')) = 0.0;
    *(Disc = new AuxField (new real[allocSize], nz, elmt, 'L')) = 0.0;
    addField[iAdd++] = InvQ;
    addField[iAdd++] = InvR;
    addField[iAdd++] = Disc;
  }

  // -- Vorticity and its dependants.

  if (need[VORTICITY])
    if (nComponent == 2) {
      vorticity[0] = new AuxField (new real[allocSize], nz, elmt, 't');
      if (add[VORTICITY])
	addField[iAdd++] = vorticity[0];
    } else {
      for (i = 0; i < nComponent; i++)
	vorticity[i] = new AuxField (new real[allocSize], nz, elmt, 'r' + i);
      if (add[VORTICITY])
	for (i = 0; i < 3; i++) addField[iAdd++] = vorticity[i];
    }

  if (need[ENSTROPHY]) {
    Ens = new AuxField (new real[allocSize], nz, elmt, 'e');
    if (add[ENSTROPHY]) addField[iAdd++] = Ens;
  }
  
  if (need[HELICITY]) {
    Hel = new AuxField (new real[allocSize], nz, elmt, 'h');
    if (add[HELICITY]) addField[iAdd++] = Hel;
  }

  if (need[ELLIPTICITY]) {
    Ell = new AuxField (new real[allocSize], nz, elmt, 'b');
    if (add[ELLIPTICITY]) addField[iAdd++] = Ell;
  }

  // -- Rate of strain tensor and its dependants.

  if (need[STRAINTENSOR]) {
    Sij = new AuxField** [nComponent];
    for (k = 0, i = 0; i < nComponent; i++) {
      Sij[i] = new AuxField* [nComponent];
      for (j = 0; j < nComponent; j++) {
	if (j >= i) {
	  Sij[i][j] = new AuxField (new real[allocSize], nz, elmt, 'i' + k++);
	} else
	  Sij[i][j] = Sij[j][i];
      }
    }
    if (add[STRAINTENSOR])
      for (i = 0; i < nComponent; i++)
	for (j = i; j < nComponent; j++)
	  addField[iAdd++] = Sij[i][j];
  }

  if (need[STRAINRATE]) {
    *(Strain = new AuxField (new real[allocSize], nz, elmt, 'G')) = 0.0;
    if (add[STRAINRATE]) addField[iAdd++] = Strain;
  }

  if (need[EGROWTH]) {
    egrow = new real[allocSize]; // -- A handle for direct access to data.
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

    for (i = 0; i < nComponent ; i++)
      for (j = 0; j < nComponent ; j++) {
	(*Vij[i][j] = *velocity[i]) . gradient (j);
	D -> u[0] -> smooth (Vij[i][j]);
	Vij[i][j] -> transform (INVERSE);
      }

    if (nComponent > 2) D -> transform (INVERSE);

    if (system == Geometry::Cylindrical && nComponent == 3) {
      for (i = 0; i < nComponent; i++) Vij[i][2] -> divY();
      (*work = *D -> u[2]) . divY();  *Vij[1][2] -= *work;
      (*work = *D -> u[1]) . divY();  *Vij[2][2] += *work;
    }
    
    if (need[INVARIANTS]) {

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
		     int*   flag   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] =
    "Usage: %s [options] -s session dump.fld\n"
    "options:\n"
    "  -h        ... print this message \n"
    "  -v        ... add vorticity (default)\n"
    "  -e        ... add enstrophy\n"
    "  -H        ... add helicity (3D only)\n"
    "  -b        ... add ellipticity\n"
    "  -B        ... add elliptic instability growth rate\n"
    "  -d        ... add divergence\n"
    "  -t        ... add components of rate of strain tensor Sij\n"
    "  -g        ... add strain rate magnitude sqrt(2SijSij)\n"
    "  -i        ... add invariants and discriminant of"
                     " velocity gradient tensor\n"
    "                NB: divergence is assumed to be zero. (3D only)\n"
    "  -a        ... add all fields derived from velocity (above)\n"
    "  -f <func> ... add a computed function <func> of x, y, z, t, etc.\n";
              
  int  i, sum;
  char buf[StrMax];
 
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
    case 'v': flag[VORTICITY]    = 1; break;
    case 'e': flag[ENSTROPHY]    = 1; break;
    case 'H': flag[HELICITY]     = 1; break;
    case 'B': flag[EGROWTH]      = 1; break;
    case 'b': flag[ELLIPTICITY]  = 1; break;
    case 'd': flag[DIVERGENCE]   = 1; break;
    case 'g': flag[STRAINRATE]   = 1; break;
    case 't': flag[STRAINTENSOR] = 1; break;
    case 'i': flag[INVARIANTS]   = 1; break;
    case 'a': Veclib::fill(FLAG_MAX, 1, flag, 1); break;
    case 'f':
      if (*++argv[0]) func = *argv; else { --argc; func = *++argv; }
      flag[FUNCTION] = 1; break;
    default: sprintf (buf, usage, prog); cout<<buf; exit(EXIT_FAILURE); break;
    }

  sum = Veclib::sum (FLAG_MAX, flag, 1);
  if (!sum) flag[VORTICITY] = 1;

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static int getDump (Domain*   D   ,
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
		      int                nOut    ,
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

  int    i, nComponent;
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
