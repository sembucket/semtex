//////////////////////////////////////////////////////////////////////////////
// addfield.C: process sem field files, computing and adding vorticity and
//            div, rate of strain tensor and invariants
//
// Usage:
// -----
// addvort -s session session.fld
//////////////////////////////////////////////////////////////////////////////

#include "Sem.h"
#include <new.h>
#include <time.h>

static char
RCSid[] = "$Id$";

static char  prog[] = "addfield";
static void  memExhaust () { message ("new", "free store exhausted", ERROR); }

static void  getargs  (int, char**, char*&, char*&, int&, int&, int&, int&, int&);
static int   getDump  (Domain*, istream&);
static void  putDump  (Domain*, vector<AuxField*>&, int, ostream&);


int main (int argc, char *argv[])
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{

  set_new_handler (&memExhaust);

  Geometry::CoordSys system;
  char      *session, *dump, fields[StrMax];
  int       i, j, np, nz, nel, DIM;
  ifstream  file;
  FEML*     F;
  Mesh*     M;
  BCmgr*    B;
  Domain*   D;
  AuxField* Ens;
  AuxField* Hel;
  AuxField* Div;
  AuxField* InvQ;
  AuxField* InvR;
  AuxField* Disc;
  AuxField* work;
  int       nData = 0;
  vector <AuxField*> AuxPoint(3);
  vector <AuxField*> dataField(12);
  AuxField*** Sij;

//  Geometry::CoordSys space=Geometry::system();

  // -- Set command line defaults

  int add_vort=0,
      add_enst=0,
      add_div =0,
      add_sij =0,
      add_inv =0;
  
  Femlib::prep ();
  getargs      (argc, argv, session, dump, add_vort, add_enst, add_div, add_sij, add_inv);

  // -- Set up domain.

  F = new FEML (session);
  M = new Mesh   (*F);
  B = new BCmgr  (*F);

  nel    = M -> nEl();  
  np     =  (int) Femlib::value ("N_POLY");
  nz     =  (int) Femlib::value ("N_Z"   );
  system = ((int) Femlib::value ("CYLINDRICAL") ) ?
                     Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, system);
  if   (nz > 1) strcpy (fields, "uvwp");
  else          strcpy (fields, "uvp");
  DIM = Geometry::nDim();

  if ( DIM < 3 ) add_enst = 0;

  vector<AuxField*> vorticity ((DIM == 2) ? 1 : 3);

  D = new Domain (*F, *M, *B, fields, session);

  file.open (dump);

  // we want the velocity gradient tensor irrespective of
  // what fields we want to output

  AuxField*** Vij = new AuxField** [DIM];
  for ( i = 0 ; i < DIM ; i++) {
    Vij[i] = new AuxField* [DIM];
    for ( j = 0 ; j < DIM ; j++ ) {
      Vij[i][j] = new AuxField (D -> Esys);
      *Vij[i][j] = 0.0;
    }
  }
  
  // -- Create workspace and dataField storage areas.

  work = new AuxField (D -> Esys);

  if(add_vort + add_enst > 0) {
    if (DIM == 2)
      vorticity[0] = new AuxField (D -> Esys, 't');
    else
      for (i = 0; i < DIM; i++)
	vorticity [i] = new AuxField (D -> Esys, 'r' + i);
    if(add_vort == 1) {
      for( i = 0 ; i < 2*DIM - 3 ; i++ ) dataField(i) = vorticity[i];
      nData += 2*DIM - 3;
    }
    if(add_enst == 1) {
      dataField(nData)   = new AuxField (D -> Esys, 'e'); Ens = dataField(nData);
      dataField(nData+1) = new AuxField (D -> Esys, 'h'); Hel = dataField(nData+1);
      nData += 2;
    }
  }

  if(add_div+add_inv > 0) {
    dataField(nData) = new AuxField (D -> Esys, 'd'); Div = dataField(nData);
    nData += 1;
  }
  
  int iadd = 0;

  if( add_sij+add_inv > 0 ) {
    Sij = new AuxField** [DIM];
    for ( i = 0 ; i < DIM ; i++) {
      Sij[i] = new AuxField* [DIM];
      for ( j = 0 ; j < DIM ; j++ ) {
	if( i <= j ) {
	  Sij[i][j] = new AuxField (D -> Esys, 'i' + iadd);
	  iadd ++;
	}
	else
	  Sij[i][j] = Sij[j][i];
      }
    }
    if( add_sij == 1 ){
      for ( i = 0 ; i < DIM ; i++) {
	for ( j = i ; j < DIM ; j++ ) {
	  dataField(nData) = Sij[i][j];
	  nData += 1;
	}
      }
    }
    if( add_inv == 1 ) {
      dataField(nData)   = new AuxField (D -> Esys, 'Q'); InvQ = dataField(nData);   	*InvQ = 0.0;
      dataField(nData+1) = new AuxField (D -> Esys, 'R'); InvR = dataField(nData+1); 	*InvR = 0.0;
      dataField(nData+2) = new AuxField (D -> Esys, 'L'); Disc = dataField(nData+2); 	*Disc = 0.0;
      nData += 3;
    }
  }
  
  // -- Cycle through field dump, first computing the velocity gradient 
  // -- tensor and then the other quantities - vorticity etc. Then write
  // -- output defined in the dataField AuxFields. musn't change the order
  // -- below because we are relying on the inverse transform being done
  // -- BEFORE dumping enstrophy/helicity because inner products are done
  // -- in physical space
  
  
  while (getDump (D, file)) {
    
    if( DIM > 2 ) D -> transform (+1);

    // Velocity gradient tensor - calculated in Fourier, then transformed back

    for ( i = 0 ; i < DIM ; i++ ) {
      for ( j = 0 ; j < DIM ; j++ ) {
	(*Vij[i][j] = *D -> u[i]) . gradient (j);
	D -> u[0] -> smooth (Vij[i][j]);
	Vij[i][j] -> transform(-1);
      }
    }

    if( DIM > 2 ) D -> transform (-1);

    if(system == Geometry::Cylindrical && DIM == 3) {
      
      for( i=0; i<DIM; i++ ) Vij[i][2] -> divR();
      (*work = *D -> u[2]) . divR(); *Vij[1][2] -= *work;
      (*work = *D -> u[1]) . divR(); *Vij[2][2] += *work;

    }
	
    if(add_div + add_inv > 0)  {
      *Div = 0.0;
      for ( i = 0 ; i < DIM ; i++ ) *Div += *Vij[i][i];
    }
    
    if(add_vort+add_enst > 0) {
      
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

    if( add_sij+add_inv > 0 ) {
      
      for( i=0; i<DIM; i++ ) {
	for( j=1; j<DIM; j++ ) {
	  *Sij[i][j]  = *Vij[i][j];
	  *Sij[i][j] += *Vij[j][i];
	  *Sij[i][j] *= 0.5;
	}
      }

      if( add_inv == 1 ) {

	// 2nd invariant (Q from Chong et al.)

	InvQ -> times      (*Vij[0][0],*Vij[1][1]);
	InvQ -> timesMinus (*Vij[0][1],*Vij[1][0]);
	InvQ -> timesPlus  (*Vij[0][0],*Vij[2][2]);
	InvQ -> timesMinus (*Vij[0][2],*Vij[2][0]);
	InvQ -> timesPlus  (*Vij[1][1],*Vij[2][2]);
	InvQ -> timesMinus (*Vij[1][2],*Vij[2][1]);

	// 3rd invariant - determinant of Vij (R from Chong et al.)

	work -> times      (*Vij[1][1],*Vij[2][2]);
	work -> timesMinus (*Vij[2][1],*Vij[1][2]);
	InvR -> times      (*work,     *Vij[0][0]);

	work -> times      (*Vij[1][2],*Vij[2][0]);
	work -> timesMinus (*Vij[2][2],*Vij[1][0]);
	InvR -> timesPlus  (*work,     *Vij[0][1]);

	work -> times      (*Vij[2][1],*Vij[1][0]);
	work -> timesMinus (*Vij[1][1],*Vij[2][0]);
	InvR -> timesPlus  (*work,     *Vij[0][2]);


	// Discriminant of Vij (N.B. DIVERGENCE (P from Chong et al.) ASSUMED = 0)

	work -> times (*InvQ,*InvQ);
	Disc -> times (*work,*InvQ);
	work -> times (*InvR,*InvR);
	*work *= 6.75;
	*Disc += *work;

      }
    }    
    
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


static void getargs (int argc       ,
		     char** argv    ,
		     char*& session ,
		     char*& dump    ,
		     int&   add_vort,
		     int&   add_enst,
		     int&   add_div ,
		     int&   add_sij ,
                     int&   add_inv )

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
		   "  -t   ... add rate of strain tensor Sij\n"
		   "  -i   ... add invariants of Sij (but NOT Sij itself) \n"
		   "           (NB - Divergence is ASSUMED equal to zero) \n"
		   "  -a   ... add them all \n";

  char buf[StrMax], c;
 
  while (--argc  && **++argv == '-')
    switch (c = *++argv[0]) {
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

  if(add_vort+add_enst+add_div+add_inv+add_sij == 0) message (prog, "No fields to add", ERROR);

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static int getDump (Domain*  D   ,
		    istream& dump)
// ---------------------------------------------------------------------------
// Read next set of field dumps from file.
// ---------------------------------------------------------------------------
{
  dump >> *D;
  return dump.good ();
}


static void putDump  (Domain*            D    ,
		      vector<AuxField*>& outField,
		      int                nOut,
		      ostream&           strm)
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

  int       i;
  const int DIM = (D -> nField() == 3) ? 2 : 3;
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



