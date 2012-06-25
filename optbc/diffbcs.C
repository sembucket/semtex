#include <stab.h>
static char             prog[] = "diffbcs";
static char*            session;
static char*      fileuc1;
static char*    fileuc2;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void  getargs    (int, char**, char*&, char*&, char*&);
static int_t preprocess (const char*, bool&);
static real_t norm_control (vector<real_t*>); 
static real_t norm_control_mixed (vector<real_t*>,vector<real_t*>); 
void initializebc (int_t, real_t*, char*); 
int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  real_t   normc1, normc2,normc_mixed;
  int_t     i, j, k;
  

  bool      restart = false;


  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, session, fileuc1, fileuc2);
  
	//load base and restart files  
  const int_t ntot = preprocess (session, restart);
 
 // Allocate memory for adjoint_integration, uc and lagragian gradient.
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP = Geometry::planeSize();


  
  int_t sizecontrolbc =  domain-> u[0] -> size_controlbc();
   
  vector<real_t*> uc1;
  vector<real_t*> uc2; 
  uc1.resize(NPERT);
  uc2.resize(NPERT);  
  real_t*     alloc_adjoint = new real_t [static_cast<size_t>(NPERT*NZ*sizecontrolbc*2)];
  for (i=0; i<NPERT; i++)        uc1[i] = alloc_adjoint+ i          * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)        uc2[i] = alloc_adjoint+(i+NPERT)   * NZ * sizecontrolbc;

	initializebc(sizecontrolbc,uc1[0],fileuc1);
initializebc(sizecontrolbc,uc2[0],fileuc2);

	normc1=norm_control(uc1);
	normc2=norm_control(uc2);
	normc_mixed=norm_control_mixed(uc1,uc2);
	printf("       %e\n\n", normc_mixed/sqrt(normc1*normc2));
  Femlib::finalize();
  return (EXIT_SUCCESS);
}

static void getargs (int        argc   ,
					 char**     argv   ,
					 char*&     session,
					 char*&     fileuc1,
					 char*&     fileuc2)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "calculate the difference between uc1 and uc2."
    "options:\n"
    "-h       ... print this message\n"
    "sessionfile, uc1 file, uc2 file   ";

  if   (--argc != 3) message (prog, "three input files are required",   ERROR);
  else         
  {  session = argv[1];
	  fileuc1=argv[2];	
	  fileuc2=argv[3];	
	
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




static real_t norm_control (vector<real_t*> uc)
// ---------------------------------------------------------------------------
// get the time-indepdent energy of the control force (from control bc). it's not norm, but norm^2
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
  real_t      norm = 0.0;
  for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> normc (uc[i]);

    return norm;
}

static real_t norm_control_mixed (vector<real_t*> uc, vector<real_t*> lag)
// ---------------------------------------------------------------------------
// get the time-indepdent energy of the control force (from control bc and Lag). it's not norm, but norm^2:[uc,lag].
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
	real_t      norm = 0.0;
	for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> normc_mixed (uc[i],lag[i]);
	
    return norm;
}

 
 


void initializebc(int_t sizeplane,
       			   real_t* uc0,
			   char* initialbc)
// randomly initially the controlbc if restart file session_initialbc does not exist.
{  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();  
   const int_t       np     = Geometry::nP();
   const int_t nbound=sizeplane/np;
   int i, j,k;
  const int_t NBASE = Geometry::nBase();
ifstream file (initialbc);
	real_t value;
    if (file)  {
     while (! file.eof() ){
       file >> *(uc0++);
	 }
     file.close();
    }
    else
	Veclib::fill (sizeplane*NPERT*NZ, 1 , uc0, 1);

	
	// if nz=2 and npert=3, wr =-wi.
	if (NBASE==2 && NPERT==3 && NZ==2) 
		Veclib::fill (sizeplane, -1 , uc0+sizeplane*4, 1);

//the same nodes have the same value
	for (i=0;i<NPERT;i++)
	 for(j=0;j<NZ;j++)
	  for(k=0;k<nbound-1;k++)
	   *(uc0+k*np+ sizeplane*(i*NZ+j))= *(uc0+k*np+2*np-1+ sizeplane*(i*NZ+j));
	  
	
//fit boundary conditions
	for (i=0;i<NPERT;i++)
	 for(j=0;j<NZ;j++)
	 {
//	 *(uc0+np-1 + sizeplane*(i*NZ+j))=0;
  	 *(uc0 + sizeplane*(i*NZ+j))=0; // TESTING::: for stenotic flow
	 //   *(uc0+sizeplane-np+ sizeplane*(i*NZ+j))=0;
	  *(uc0+sizeplane-1+ sizeplane*(i*NZ+j))=0; // TESTING::: for stenotic flow
	 }
  }
  