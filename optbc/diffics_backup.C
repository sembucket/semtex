#include <stab.h>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cstring>

#include <iostream>
#include <fstream>

#include <iomanip>
#include <vector>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <blas.h>
#include <lapack.h>
#include <veclib.h>
#include <femlib.h>
#include <math.h>

static char* hdr_fmt[] = { 
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

typedef struct hdr_data {
	char   session[StrMax];
	char   created[StrMax];
	int_t  nr, ns, nz, nel;
	int_t  step;
	real_t time;
	real_t timestep;
	real_t kinvis;
	real_t beta;
	char   fields[StrMax];
	char   format[StrMax];
} hdr_info;

static char             prog[] = "diffics";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;


static void  gethead   (istream&, hdr_info&);
static int_t preprocess (const char*, bool&);
static real_t L2norm (); 
static real_t L2norm_mixed (vector<real_t*>); 

void installdomainu (vector<real_t*>);
static void getargs (int, char**, char*&, ifstream&, ifstream&);

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{  	ifstream        u1File, u2File;	// -- b ==> base, p ==> perturbation.
	hdr_info        u1Head, u2Head;
	 bool      restart = false;
  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session, u1File, u2File);
	gethead (u1File, u1Head);
	gethead (u2File, u2Head);

 
//load base and restart files  
 int_t ntot = preprocess (session, restart);
	int_t i;
 // Allocate memory for adjoint_integration, uc and lagragian gradient.
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  ntot = NZ * NPERT * NP;

	real_t normv1, normv2, normv_mixed;


  vector<real_t*>  u1;
  vector<real_t*>  u2;
  u1.resize(NPERT);
  u2.resize(NPERT); 
  real_t*     alloc_adjoint = new real_t [static_cast<size_t>(2*ntot)];
  for (i=0; i<NPERT; i++)        u1[i] = alloc_adjoint+ i          * NZ * NP;
  for (i=0; i<NPERT; i++)        u2[i] = alloc_adjoint+ (NPERT+i)  * NZ * NP;


	
	u1File.read (reinterpret_cast<char*>(u1[0]), ntot);
	u2File.read (reinterpret_cast<char*>(u2[0]), ntot);
	u1File.close(); u2File.close();
	
	
	installdomainu(u1);
	normv1=L2norm();
	installdomainu(u2);
	normv2=L2norm();
	normv_mixed=L2norm_mixed(u1);
	printf("%e", normv_mixed/sqrt(normv1*normv2));

}



static void getargs (int       argc ,
					 char**    argv ,
					 char*&     session,
					 ifstream& u1file,
					 ifstream& u2file)
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{  if   (--argc != 3) message (prog, "three input files are required",   ERROR);
else         
{  session = argv[1];
	u1file.open (argv[2], ios::in);
	u2file.open (argv[3], ios::in);		
}
}




static int_t preprocess (const char* session,
			 bool&       restart)
// ---------------------------------------------------------------------------
// Create objects needed for semtex execution, given the session file name.
//
// Return length of an eigenvector: the amount of storage required for
// a velocity component * number of components.
// ---------------------------------------------------------------------------
{
  const real_t* z;
  int_t         i, np, nel, npert;

  // -- Set default additional tokens.

  Femlib::value  ("BASE_PERIOD", 0.0);
  Femlib::value  ("T_OFFSET",    0.0);
  Femlib::ivalue ("BIG_RESIDS",  0);

  // -- Start up dealing with session file.

  file  = new FEML (session);
  mesh  = new Mesh (file);

  np    = Femlib::ivalue ("N_P");
  nel   = mesh -> nEl();
  npert = file -> attribute ("FIELDS", "NUMBER") - 1;
  Geometry::set (nel, npert);

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);

  bman   = new BCmgr  (file, elmt);
  domain = new Domain (file, elmt, bman);

  // -- Load restart and base flow data.

  restart = domain -> restart ();
  domain -> loadBase();
  domain -> report  ();

  analyst = new StabAnalyser (domain, file);

  // -- Over-ride any CHKPOINT flag in session file.

  Femlib::value ("CHKPOINT", 1);

  return Geometry::nPert() * Geometry::planeSize() * Geometry::nZ();
}


static real_t L2norm ()
// ---------------------------------------------------------------------------
// get the L2norm of the velocity vectors. Not norm but norm^2
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
  real_t      norm = 0.0;
  for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> mode_L2 (0);
   return norm;
}


static real_t L2norm_mixed (vector<real_t*> anotheru)
// ---------------------------------------------------------------------------
// get the mixed L2norm. Not norm but norm^2
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
	real_t      norm = 0.0;
	for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> mode_L2_mixed (domain -> u[i] -> plane (0),anotheru[i]);
	return norm;
}




void installdomainu(vector<real_t*> utoinstall)
 //install utoinstall to domain
 {  const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();
    const int_t NP = Geometry::planeSize();
	int_t i,j;
	
  for (i=0;i<NPERT;i++)
   for(j=0;j<NZ;j++)
    Veclib::copy (NP, utoinstall[i]+j*NP, 1, domain -> u[i] -> plane (j), 1);
}


static void gethead (istream&  file  ,
					 hdr_info& header)
// ---------------------------------------------------------------------------
// Load data structure from file header info. Note that the field
// names are packed into a string without spaces.
// ---------------------------------------------------------------------------
{
	char  buf[StrMax];
	int_t i, j; 
	
	file.get (header.session, 25); file.getline (buf, StrMax);
	
	if (!strstr (buf, "Session")) message (prog, "not a field file", ERROR);
	
	file.get (header.created, 25); file.ignore (StrMax, '\n');
	
	file >> header.nr >> header.ns >> header.nz >> header.nel;
	file.ignore (StrMax, '\n');
	
	file >> header.step;     file.ignore (StrMax, '\n');
	file >> header.time;     file.ignore (StrMax, '\n');
	file >> header.timestep; file.ignore (StrMax, '\n');
	file >> header.kinvis;   file.ignore (StrMax, '\n');
	file >> header.beta;     file.ignore (StrMax, '\n');
	
	file.get (buf, 25);
	for (i = 0, j = 0; i < 25; i++)
		if (isalpha (buf[i])) header.fields[j++] = buf[i];
	file.ignore (StrMax, '\n');
	header.fields[j] = '\0';
	
	file.get (header.format, 25); file.getline (buf, StrMax);
	
	if (!strstr (header.format, "binary"))
		message (prog, "input field file not in binary format", ERROR);
	else if (!strstr (header.format, "-endia"))
		message (prog, "input field file in unknown binary format", WARNING);
}