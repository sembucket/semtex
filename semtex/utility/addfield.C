//////////////////////////////////////////////////////////////////////////////
// addfield.C: process sem field files, computing and adding vorticity,
//             optionally divergence, helicity and enstrophy.
//
// Usage:
// -----
// addfield [options] -s session session.fld
//   options:
//   -h     ... print this message
//   -v     ... add vorticity (default)
//   -e     ... add enstrophy and helicity (3D only)
//   -d     ... add divergence
//   -a     ... add vorticity/enstrophy/helicity/divergence
//////////////////////////////////////////////////////////////////////////////

#include <Sem.h>
#include <new.h>
#include <time.h>

static char
RCSid[] = "$Id$";

static char prog[] = "addfield";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }

static void    getargs  (integer, char**, char*&, char*&,
			 integer&, integer&, integer&);
static integer getDump  (Domain*, istream&);
static void    putDump  (Domain*, vector<AuxField*>&, integer, ostream&);


integer main (integer argc,
	      char**  argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  set_new_handler (&memExhaust);

  Geometry::CoordSys system;
  char               *session, *dump, fields[StrMax];
  integer            i, np, nz, nel, DIM;
  ifstream           file;
  FEML*              F;
  Mesh*              M;
  BCmgr*             B;
  Domain*            D;
  AuxField*          Ens;
  AuxField*          Hel;
  AuxField*          Div;
  AuxField*          work;
  integer            nData = 0;
  vector <AuxField*> AuxPoint(3);
  vector <AuxField*> dataField(10);

  // -- Set command line defaults.

  integer add_vort=0,
          add_enst=0,
          add_div =0;
  
  Femlib::prep ();
  getargs      (argc, argv, session, dump, add_vort, add_enst, add_div);

  // -- Set up domain.

  F = new FEML (session);
  M = new Mesh   (*F);
  B = new BCmgr  (*F);

  nel    = M -> nEl();  
  np     =  (integer) Femlib::value ("N_POLY");
  nz     =  (integer) Femlib::value ("N_Z"   );
  system = ((integer) Femlib::value ("CYLINDRICAL") ) ?
                     Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, system);
  if   (nz > 1) strcpy (fields, "uvwp");
  else          strcpy (fields, "uvp");
  DIM = Geometry::nDim();

  if ( DIM < 3 ) add_enst = 0;

  vector<AuxField*> vorticity ((DIM == 2) ? 1 : 3);

  D = new Domain (*F, *M, *B, fields, session);

  file.open (dump);

  // -- Create workspace and dataField storage areas.

  work = new AuxField (D -> Esys);

  if(add_vort + add_enst > 0) {
    if (DIM == 2)
      vorticity[0] = new AuxField (D -> Esys, 't');
    else
      for (i = 0; i < DIM; i++)
	vorticity [i] = new AuxField (D -> Esys, 'r' + i);
    if (add_vort == 1) {
      for (i = 0; i < 2*DIM - 3; i++) dataField(i) = vorticity[i];
      nData += 2*DIM - 3;
    }
    if (add_enst == 1) {
      dataField(nData)   = new AuxField (D -> Esys, 'e');
      Ens = dataField(nData);

      dataField(nData+1) = new AuxField (D -> Esys, 'h');
      Hel = dataField(nData+1);

      nData += 2;
    }
  }
  
  if (add_div == 1) {
    dataField(nData) = new AuxField (D -> Esys, 'd');
    Div = dataField(nData);

    nData += 1;
  }
  
  // -- Cycle through field dump, computing vorticity, writing output.
  // -- Musn't change the order below because we are relying on the
  // -- inverse transform being done BEFORE dumping enstrophy/helicity.
  
  while (getDump (D, file)) {
    
    if (DIM > 2 ) D -> transform (+1);
    
    if (add_div == 1)  {
      
      *Div = 0.0;
      
      if(system == Geometry::Cartesian)
	for( i = 0 ; i< DIM ; i++) {
	  (*work = *D -> u[i]) . gradient (i); *Div += *work;
	}
      else {
	(*work = *D -> u[0]) . gradient (0);                   *Div += *work;
	(*work = *D -> u[1]) . mulR() . gradient (1) . divR(); *Div += *work;
	if( DIM == 3 ) {
	  (*work = *D -> u[2]) . gradient (2) . divR();        *Div += *work;
	}	  
	D -> u[0] -> smooth (Div);
      }
      if ( DIM > 2) Div -> transform (-1);
    }
    
    if (add_vort+add_enst > 0) {
      
      if (DIM == 2) {
	(*work = *D -> u[1]) . gradient (0); *vorticity[0]  = *work;
	(*work = *D -> u[0]) . gradient (1); *vorticity[0] -= *work;
	D -> u[0] -> smooth (vorticity [0]);
      } else {
	if(system == Geometry::Cartesian) {
	  (*work = *D -> u[2]) . gradient (1); *vorticity[0]  = *work;
	  (*work = *D -> u[1]) . gradient (2); *vorticity[0] -= *work;
	  (*work = *D -> u[0]) . gradient (2); *vorticity[1]  = *work;
	  (*work = *D -> u[2]) . gradient (0); *vorticity[1] -= *work;
	}
	else {
	  (*work = *D -> u[2]) . mulR() . gradient (1); *vorticity[0]  = *work;
	  (*work = *D -> u[1]) . gradient (2);          *vorticity[0] -= *work;
          vorticity[0] -> divR();
	  
	  (*work = *D -> u[0]) . gradient (2) . divR(); *vorticity[1]  = *work;
	  (*work = *D -> u[2]) . gradient (0);          *vorticity[1] -= *work;
	}
	(*work = *D -> u[1]) . gradient (0); *vorticity[2]  = *work;
	(*work = *D -> u[0]) . gradient (1); *vorticity[2] -= *work;
	for (i = 0; i < DIM; i++) D -> u[i] -> smooth (vorticity [i]);
	vorticity [0] -> transform (-1);
	vorticity [1] -> transform (-1);
	vorticity [2] -> transform (-1);
      }
      
    }
    
    D -> transform (-1);
    
    if(add_enst == 1) {
      Ens -> innerProduct(vorticity,vorticity);
      for ( i = 0 ; i < DIM ; i++ ) AuxPoint[i] = D -> u[i];
      Hel -> innerProduct(vorticity,AuxPoint);
    }
      
    putDump (D, dataField, nData, cout);
    
  }
  
  file.close();
  
  return EXIT_SUCCESS;
}


static void getargs (integer  argc    ,
		     char**   argv    ,
		     char*&   session ,
		     char*&   dump    ,
		     integer& add_vort,
		     integer& add_enst,
		     integer& add_div )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[]   = "Usage: %s [options] -s session dump.fld\n"
                   "options:\n"
                   "  -h   ... print this message \n"
                   "  -v   ... add vorticity (default) \n"
		   "  -e   ... add enstrophy and helicity (3D only) \n"
		   "  -d   ... add divergence\n"
		   "  -a   ... add vorticity/enstrophy/helicity/divergence\n";

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
    case 'a':
      add_enst = 1;
      add_vort = 1;
      add_div = 1;
      break;

    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (add_vort+add_enst+add_div == 0) add_vort=1;

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static integer getDump (Domain*  D   ,
			istream& dump)
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



