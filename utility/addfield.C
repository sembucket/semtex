//////////////////////////////////////////////////////////////////////////////
// addfield.C: process sem field files, computing and adding vorticity and
//             div, rate of strain tensor and invariants.
//
// Copyright (c) 1998--1999 Hugh Blackburn, Murray Rudman
//
// Usage:
// -----
// addvort [options] -s session session.fld
//   options:
//   -h   ... print this message
//   -v   ... add vorticity
//   -e   ... add enstrophy and helicity (3D only)
//   -d   ... add divergence
//   -t   ... add rate of strain tensor Sij
//   -g   ... add gamma = total strain rate = sqrt(2 Sij Sji)
//   -i   ... add invariants of Vij (but NOT Vij itself)
//            (NB: Divergence is ASSUMED equal to zero)
//   -a   ... add them all
//
// Field names used/assumed here:
//
// u -- x velocity component (cylindrical: axial)
// v -- y velocity component (cylindrical: radial)
// w -- z velocity component (cylindrical: azimuthal)
// p -- pressure / density
// A -- uu covariance
// B -- uv covariance
// C -- vv covariance
// D -- uw covariance
// E -- vw covariance
// F -- ww covariance
// r -- x component vorticiy
// s -- y component vorticity
// t -- z component vorticity
// d -- divergence
// Q -- 2nd invariant of velocity gradient tensor
// R -- 3rd invariant of velocity gradient tensor
// L -- Discriminant  of velocity gradient tensor 27/4 R^2 + Q^3
// G -- Strain rate sqrt (2 Sij Sji)
// e -- enstrophy = r^2 + s^2 + t^2
// h -- helicity  = ur  + vs  + wt
// i -- uu strain rate component
// j -- uv strain rate component
// k -- vv strain rate component
// l -- uw strain rate component
// m -- vw strain rate component
// n -- ww strain rate component
//
// $Id$
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>
#include <new.h>
#include <time.h>

static char  prog[] = "addfield";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }

static void    getargs  (integer, char**, char*&, char*&,
			 integer&, integer&, integer&,
			 integer&, integer&, integer&);
static integer getDump  (Domain*, ifstream&);
static void    putDump  (Domain*, vector<AuxField*>&, integer, ofstream&);


integer main (integer argc,
	      char**  argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);

  Geometry::CoordSys system;
  char               *session, *dump, fields[StrMax];
  integer            i, j, np, nz, nel, DIM;
  ifstream           file;
  ofstream           outp (1);
  FEML*              F;
  Mesh*              M;
  BCmgr*             B;
  Domain*            D;
  AuxField*          Ens;
  AuxField*          Hel;
  AuxField*          Div;
  AuxField*          InvQ;
  AuxField*          InvR;
  AuxField*          Disc;
  AuxField*          Strain;
  AuxField*          work;
  integer            nData = 0;
  const real*        z;
  vector <Element*>  elmt;
  vector <AuxField*> AuxPoint(3);
  vector<AuxField*>  vorticity;
  vector <AuxField*> dataField(12);
  AuxField***        Sij;

  // -- Set command line defaults

  integer add_vort   = 0,
          add_enst   = 0,
          add_div    = 0,
          add_sij    = 0,
          add_strain = 0,
          add_inv    = 0;
  
  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, session, dump,
	   add_vort, add_enst, add_div, add_sij, add_inv, add_strain);

  // -- Set up domain.

  F = new FEML  (session);
  M = new Mesh  (*F);


  nel    = M -> nEl();  
  np     =  (integer) Femlib::value ("N_POLY");
  nz     =  (integer) Femlib::value ("N_Z"   );
  system = ((integer) Femlib::value ("CYLINDRICAL") ) ?
                      Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, system);
  if   (nz > 1) strcpy (fields, "uvwp");
  else          strcpy (fields, "uvp");
  DIM = Geometry::nDim();

  vorticity.setSize ((DIM == 2) ? 1 : 3);
  if (DIM < 3) add_enst = 0;

  Femlib::mesh (GLL, GLL, np, np, &z, 0, 0, 0, 0);

  elmt.setSize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, M, z, np);

  B = new BCmgr (F, elmt);
  D = new Domain (F, elmt, B);

  file.open (dump);

  // -- Compute the velocity gradient tensor irrespective of
  //    what fields we want to output.

  AuxField*** Vij = new AuxField** [DIM];
  for (i = 0; i < DIM; i++) {
    Vij[i] = new AuxField* [DIM];
    for (j = 0; j < DIM; j++) {
       Vij[i][j] = new AuxField (elmt);
      *Vij[i][j] = 0.0;
    }
  }
  
  // -- Create workspace and dataField storage areas.

  work = new AuxField (elmt);

  if (add_vort + add_enst) {
    if (DIM == 2)
      vorticity[0] = new AuxField (elmt, 't');
    else
      for (i = 0; i < DIM; i++)
	vorticity [i] = new AuxField (elmt, 'r' + i);
    if (add_vort == 1) {
      for (i = 0 ; i < 2*DIM - 3 ; i++ ) dataField(i) = vorticity[i];
      nData += 2*DIM - 3;
    }
    if (add_enst == 1) {
      dataField(nData)   = new AuxField (elmt, 'e');
      dataField(nData+1) = new AuxField (elmt, 'h');
      Ens = dataField(nData);
      Hel = dataField(nData+1);
      nData += 2;
    }
  }

  if(add_div+add_inv > 0) {
    dataField(nData) = new AuxField (elmt, 'd'); Div = dataField(nData);
    nData += 1;
  }
  
  integer iadd = 0;

  if (add_sij + add_inv + add_strain) {
    Sij = new AuxField** [DIM];
    for (i = 0; i < DIM; i++) {
      Sij[i] = new AuxField* [DIM];
      for ( j = 0 ; j < DIM ; j++ ) {
	if (i <= j) {
	  Sij[i][j] = new AuxField (elmt, 'i' + iadd);
	  iadd++;
	} else
	  Sij[i][j] = Sij[j][i];
      }
    }
    if (add_sij == 1)
      for ( i = 0 ; i < DIM ; i++)
	for ( j = i ; j < DIM ; j++ ) {
	  dataField(nData) = Sij[i][j];
	  nData++;
	}
    if (add_inv == 1) {
      dataField(nData)   = new AuxField (elmt, 'Q');
      dataField(nData+1) = new AuxField (elmt, 'R');
      dataField(nData+2) = new AuxField (elmt, 'L');

      *(InvQ = dataField(nData))   = 0.0;
      *(InvR = dataField(nData+1)) = 0.0;
      *(Disc = dataField(nData+2)) = 0.0;
      nData += 3;
    }
    if (add_strain == 1) {
      dataField(nData) = new AuxField (elmt, 'G');
      *(Strain = dataField(nData)) = 0.0;
      nData++;
    }
  }
  
  // -- Cycle through field dump, first computing the velocity gradient 
  //    tensor and then the other quantities - vorticity etc. Then write
  //    output defined in the dataField AuxFields. musn't change the order
  //    below because we are relying on the inverse transform being done
  //    BEFORE dumping enstrophy/helicity because inner products are done
  //    in physical space.
  
  
  while (getDump (D, file)) {
    
    if( DIM > 2 ) D -> transform (+1);

    // -- Velocity gradient tensor, calculated in Fourier, transformed back.

    for (i = 0; i < DIM ; i++)
      for (j = 0; j < DIM ; j++) {
	(*Vij[i][j] = *D -> u[i]) . gradient (j);
	D -> u[0] -> smooth (Vij[i][j]);
	Vij[i][j] -> transform(-1);
      }

    if (DIM > 2) D -> transform (-1);

    if (system == Geometry::Cylindrical && DIM == 3) {
      
      for (i = 0; i < DIM; i++) Vij[i][2] -> divR();
      (*work = *D -> u[2]) . divR(); *Vij[1][2] -= *work;
      (*work = *D -> u[1]) . divR(); *Vij[2][2] += *work;

    }
	
    if (add_div + add_inv) {
      *Div = 0.0;
      for (i = 0; i < DIM; i++) *Div += *Vij[i][i];
    }
    
    if (add_vort + add_enst) {
      
      if (DIM == 2) {
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

    if (add_sij + add_inv + add_strain) {
      
      for (i = 0; i < DIM; i++)
	for (j = i; j < DIM; j++) {
	  *Sij[i][j]  = *Vij[i][j];
	  *Sij[i][j] += *Vij[j][i];
	  *Sij[i][j] *= 0.5;
	}

      if (add_inv == 1) {

	// -- 2nd invariant (Q from Chong et al.).

	InvQ -> times      (*Vij[0][0], *Vij[1][1]);
	InvQ -> timesMinus (*Vij[0][1], *Vij[1][0]);
	InvQ -> timesPlus  (*Vij[0][0], *Vij[2][2]);
	InvQ -> timesMinus (*Vij[0][2], *Vij[2][0]);
	InvQ -> timesPlus  (*Vij[1][1], *Vij[2][2]);
	InvQ -> timesMinus (*Vij[1][2], *Vij[2][1]);

	// -- 3rd invariant - determinant of Vij (R from Chong et al.).

	work -> times      (*Vij[1][1], *Vij[2][2]);
	work -> timesMinus (*Vij[2][1], *Vij[1][2]);
	InvR -> times      (*work,      *Vij[0][0]);

	work -> times      (*Vij[1][2], *Vij[2][0]);
	work -> timesMinus (*Vij[2][2], *Vij[1][0]);
	InvR -> timesPlus  (*work,      *Vij[0][1]);

	work -> times      (*Vij[2][1], *Vij[1][0]);
	work -> timesMinus (*Vij[1][1], *Vij[2][0]);
	InvR -> timesPlus  (*work,      *Vij[0][2]);


	// -- Discriminant of Vij
	//    NB: DIVERGENCE (P from Chong et al.) ASSUMED = 0.

	work -> times (*InvQ, *InvQ);
	Disc -> times (*work, *InvQ);
	work -> times (*InvR, *InvR);
	*work *= 6.75;
	*Disc += *work;

      }
    }    

    if (add_strain == 1) {
      *Strain = 0.0;
      for (i = 0; i < DIM; i++)
	for (j = 0; j < DIM; j++)
	  Strain -> timesPlus ( *Sij[i][j], *Sij[j][i] );
      *Strain *= 2.0;
    }
    
    if (add_enst == 1) {
      Ens -> innerProduct (vorticity, vorticity);
      for (i = 0; i < DIM; i++) AuxPoint[i] = D -> u[i];
      Hel -> innerProduct (vorticity, AuxPoint);
    }
      
    putDump (D, dataField, nData, outp);
    
  }
  
  file.close();
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (integer  argc    ,
		     char**   argv    ,
		     char*&   session ,
		     char*&   dump    ,
		     integer& add_vort,
		     integer& add_enst,
		     integer& add_div ,
		     integer& add_sij ,
		     integer& add_inv ,
                     integer& add_strain )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: %s [options] -s session dump.fld\n"
                 "options:\n"
                 "  -h   ... print this message \n"
                 "  -v   ... add vorticity\n"
		 "  -e   ... add enstrophy and helicity (3D only) \n"
		 "  -d   ... add divergence\n"
		 "  -t   ... add rate of strain tensor Sij\n"
		 "  -g   ... add gamma = total strain rate = sqrt(2 Sij Sji)\n"
		 "  -i   ... add invariants of Vij (but NOT Vij itself) \n"
		 "           (NB: Divergence is ASSUMED equal to zero) \n"
                 "  -a   ... add them all \n";

  char buf[StrMax];
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 's':
      if (*++argv[0])
	session = *argv;
      else {
	--argc;
	session = *++argv;
      }
      break;
    case 'v':
      add_vort = 1;
      break;
    case 'e':
      add_enst = 1;
      break;
    case 'd':
      add_div = 1;
      break;
    case 'g':
      add_strain = 1;
      break;
    case 't':
      add_sij = 1;
      break;
    case 'i':
      add_inv = 1;
      add_div = 1;
      break;
    case 'a':
      add_enst = 1;
      add_vort = 1;
      add_div = 1;
      add_sij = 1;
      add_inv = 1;
      break;

    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (add_vort+add_enst+add_div+add_inv+add_sij+add_strain == 0) add_vort=1;
  // message (prog, "No fields to add", ERROR);

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static integer getDump (Domain*   D   ,
			ifstream& dump)
// ---------------------------------------------------------------------------
// Read next set of field dumps from file.
// ---------------------------------------------------------------------------
{
  dump >> *D;
  return dump.good ();
}


static void putDump  (Domain*            D       ,
		      vector<AuxField*>& outField,
		      integer            nOut    ,
		      ofstream&          strm    )
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

  integer       i;
  const integer DIM = (D -> nField() == 3) ? 2 : 3;
  char      routine[] = "putDump";
  char      s1[StrMax], s2[StrMax];
  time_t    tp (::time (0));

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

  for (i = 0; i <= DIM; i++) s2[i] = D -> u[i] -> name();
  for (i = 0; i <  nOut; i++)
    s2[DIM + i + 1] = outField(i) -> name();
  s2[DIM + nOut + 1] = '\0';

  sprintf (s1, hdr_fmt[8], s2);
  strm << s1;

  sprintf (s2, "binary ");
  Veclib::describeFormat (s2 + strlen (s2));
  sprintf (s1, hdr_fmt[9], s2);
  strm << s1;

  for (i = 0; i <= DIM; i++) strm << *D -> u[i];
  for (i = 0; i < nOut; i++) strm << *outField(i);

  if (!strm) message (routine, "failed writing field file", ERROR);
  strm << flush;
}
