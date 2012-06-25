#include <stab.h>
static char             prog[] = "dog";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;
static real_t L2norm (); 

static void  getargs    (int, char**, char*&, problem_t&, real_t& );
static int_t preprocess (const char*, bool&);
static void Zerodomain ();
static real_t norm_control (vector<real_t*>); 
void dumpbc (int_t, real_t*, char*);
void initializebc (int_t, real_t*, char*); 
int main (int    argc,
		  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
	int_t     i, j, k;
	real_t  Ec=1,normic;
	char      buf[StrMax];
	char      initialbc[StrMax],finalbc[StrMax],growthbc[StrMax];
	ofstream  growthfile;
	bool      restart = true;
    problem_t        task = PRIMAL;
	Femlib::ivalue ("conjugate_gradient", 0); //default is steepest method. using conjugate_gradient= 1 to switch to CG.
	Femlib::ivalue ("TIMEDECAY", 0);   //default is the double Gaussian time decay factor.
	Femlib::initialize (&argc, &argv);
	
	getargs (argc, argv, session, task, Ec);
	
	strcat (strcpy (initialbc, session), ".rstbc");
	strcat (strcpy (finalbc, session), ".fldbc");
	strcat (strcpy (growthbc, session), ".growthbc");

	
	//load base and restart files  
	const int_t ntot = preprocess (session, restart);
	analyst = new StabAnalyser (domain, file);
	// Allocate memory for adjoint_integration, uc and lagragian gradient.
	const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();
	const int_t NP = Geometry::planeSize();
	int_t sizecontrolbc =  domain-> u[0] -> size_controlbc();
	
	vector<real_t*>  adjoint_integration;
	vector<real_t*>  uc;
	
	adjoint_integration.resize(NPERT);
	uc.resize(NPERT);

	
	
	
	real_t*     alloc_adjoint = new real_t [static_cast<size_t>(NPERT*NZ*sizecontrolbc*2)];
	for (i=0; i<NPERT; i++)      adjoint_integration[i] = alloc_adjoint+ i          * NZ * sizecontrolbc;
	for (i=0; i<NPERT; i++)                       uc[i] = alloc_adjoint+(i+NPERT)   * NZ * sizecontrolbc;
    
	
	//set the control force for v and w zero to avoid singularity in cylindrical condition??
	initializebc(sizecontrolbc,uc[0],initialbc);
	//normalize uc
	//Veclib::smul (NPERT*NZ*sizecontrolbc, sqrt(Ec/norm_control(uc)), uc[0], 1, uc[0], 1);
			
			
	if (task == PRIMAL)
	{  	if (restart ==0)  Zerodomain ();
		//normalize the initial condition
	else {
		normic=L2norm();
		for (i=0; i<NPERT; i++) 
		Veclib::smul (Geometry::planeSize()*NZ, 1.0/sqrt(normic), domain->u[i]->plane(0), 1,  domain->u[i]->plane(0), 1);
	}
	
		integrate (linAdvect, domain, analyst, adjoint_integration,uc);
	}
	else
	{   if (restart ==0)   message (prog, "initial condition is required",   ERROR);
		Veclib::fill (sizecontrolbc*NPERT*NZ, 0 , uc[0], 1);
		Veclib::fill (sizecontrolbc*NPERT*NZ, 0 , adjoint_integration[0], 1);
		integrate (linAdvectT, domain, analyst, adjoint_integration,uc);
		 dumpbc(sizecontrolbc,adjoint_integration[0],finalbc);
	}
			
	

	Femlib::finalize();
	return (EXIT_SUCCESS);
}

static void getargs (int        argc   ,
					 char**     argv   ,
					  char*&     session,
					  problem_t& task ,
					  real_t&     Ec  )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
	char usage[] = "dsa(-H) [options] session\n TG only."
    "options:\n"
	" -a      ... ADJOINT evolution\n"
    "-h       ... print this message\n"
	"-e       ... set energy norm of boundary pert, default value is Ec=1 \n";
	
	while (--argc && **++argv == '-')
		switch (*++argv[0]) {
			case 'a':
				task = ADJOINT;
				break;
		    case 'l':	
				do
					Femlib::ivalue ("FLinterpolation", Femlib::ivalue ("FLinterpolation") + 1);
				while (*++argv[0] == 'l');
				break;
			case 'i':
				do
					Femlib::ivalue ("ITERATIVE", Femlib::ivalue ("ITERATIVE") + 1);
				while (*++argv[0] == 'i');
				break;
			case 'e':
				if (*++argv[0]) Ec  = atof (  *argv);
				else { --argc;  Ec  = atof (*++argv); }
				break;	
			default:
				cerr << usage;
				exit (EXIT_FAILURE);
				break;
		}
	
	if   (argc != 1) message (prog, "no session file",   ERROR);
	else             session = *argv;
	
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





static void Zerodomain ()
// ---------------------------------------------------------------------------
// Set zero initial condition.
// ---------------------------------------------------------------------------
{
	int_t       i, k;
	const int_t ND = Geometry::nPert();
	const int_t NP = Geometry::planeSize();
	const int_t NZ = Geometry::nZ();
	real_t*  zerou        = new real_t [NP]; 
	Veclib::zero (NP,     zerou, 1);  
	for (i = 0; i < ND; i++)
		for (k = 0; k < NZ; k++)
			domain -> u[i] -> setPlane (k,zerou);
	
	delete [] zerou;
}










void initializebc(int_t sizeplane,
				  real_t* uc0,
				  char* initialbc)
// zero initially the controlbc if restart file session_initialbc does not exist.
{  const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();  
	const int_t       np     = Geometry::nP();
	const int_t nbound=sizeplane/np;
	int i, j,k;
	ifstream file (initialbc);
	real_t value;
    if (file)  {
		while (! file.eof() ){
			file >> *(uc0++);
		}
		file.close();
    }
    else{
		//Veclib::zero (sizeplane*NPERT*NZ, uc0, 1);
		//Veclib::vnormal (sizeplane*NPERT*NZ, 0.0, 1.0, uc0, 1);
		Veclib::fill (sizeplane*NPERT*NZ, 0 , uc0, 1);
    }
	

}




void  dumpbc(int_t size,
			 real_t* uc0,
			 char* finalbc)
//dump the adjoint control bc to session_finalbc.
{const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();
	char      name[StrMax];
	strcpy (name,finalbc);
	strcat(name,"optimalbc");
	ofstream file;
	ofstream file1;
	file.open(name);
	file1.open(finalbc);
	
	int_t i,j,k;
	real_t* meshx;
	real_t* meshy;
	meshx= new real_t [size];
	meshy= new real_t [size];
	
	domain -> u[0] -> controlmesh(meshx,meshy);
	file << "x, y, realu, imagu, realv,imagv,realw,imagw"<<endl;
	for (i=0; i < size ; i++, uc0++)
	{file<<setprecision(20)<< *(meshx+i)  
		<< setw( 25) <<setprecision(20)<< *(meshy+i) <<" "; 
		for (k=0; k < NPERT; k++)
			for (j=0; j < NZ   ; j++)
				file<< setw( 30)<<setprecision(20) <<*(uc0+size*(j+k*NZ)); 
		file<< endl;
	}
	file.close();
	
	uc0 -=size;
	for (i=0;i<size*NZ*NPERT;i++,uc0++)
		file1<<setprecision(20)<<*uc0  <<endl;
	file1.close();
	
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


static real_t L2norm ()
// ---------------------------------------------------------------------------
// get the L2norm of the velocity vectors. Not norm but norm^2
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
	real_t      norm = 0.0;
	for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> mode_L2 (0);
	return norm;
}


