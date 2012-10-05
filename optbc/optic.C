#include <stab.h>
static char             prog[] = "optic";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void  getargs    (int, char**, real_t&, int_t&, real_t&, char*&);
static int_t preprocess (const char*, bool&);
static real_t L2norm (); 
static real_t L2norm_mixed (vector<real_t*>); 
void choose_step_length(vector<real_t*>, vector<real_t*>,vector<real_t*>,vector<real_t*>, real_t,real_t, real_t&);
void getu(vector<real_t*>);
void installdomainu (vector<real_t*>);
real_t innerproduct (int_t, real_t*, real_t*);

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  int_t     converged = 0,iterations=0;
  real_t   steplengthcommand=0, steplength, resnorm, evtol = 1.0e-16, ro=0.1,c1const=0.000000000000001;
  int_t     i, j, k;
  char      buf[StrMax];
  char      initialbc[StrMax],finalbc[StrMax],growthbc[StrMax];
  ofstream  growthfile;
  bool      restart;
  
  Femlib::ivalue ("TIMEDECAY", 0);   //default is the double Gaussian time decay factor.
  Femlib::initialize (&argc, &argv);

  getargs (argc, argv, steplengthcommand, iterations, evtol, session);

  strcat (strcpy (growthbc, session), ".growthbc");

  growthfile.open(growthbc);

//load base and restart files  
  const int_t ntot = preprocess (session, restart);

 // Allocate memory for adjoint_integration, uc and lagragian gradient.
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  real_t energygrowth[iterations];

  real_t normu0[iterations];
  real_t normu0_mixed[iterations];
  real_t normutau[iterations];
  real_t norml[iterations];
  real_t norma[iterations];


  int_t sizecontrolbc =  domain-> u[0] -> size_controlbc();
  if (sizecontrolbc >0)  message (prog, "this is an optimal ic problem, therefore no control bc",   ERROR);
  sizecontrolbc=1;
	  
  
	  
  vector<real_t*>  adjoint_integration;
  vector<real_t*>  uc;
  vector<real_t*>  uu0;  //field u from uc
  vector<real_t*>  udu0;  //field u from \delta uc
  vector<real_t*>  u0;
  vector<real_t*>  direction;
  vector<real_t*>  lagrangian_gradient;

  
  adjoint_integration.resize(NPERT);
  uc.resize(NPERT);
  uu0.resize(NPERT);
  udu0.resize(NPERT);
  u0.resize(NPERT);  
  direction.resize(NPERT);
  lagrangian_gradient.resize(NPERT);	
  
  real_t*     alloc_adjoint = new real_t [static_cast<size_t>(NPERT*NZ*sizecontrolbc*2+5*NZ*NPERT*NP)];
  for (i=0; i<NPERT; i++)        adjoint_integration[i] = alloc_adjoint+ i          * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                         uc[i] = alloc_adjoint+(i+1*NPERT) * NZ * sizecontrolbc;

  for (i=0; i<NPERT; i++)                        uu0[i] = alloc_adjoint+   2*NPERT  * NZ * sizecontrolbc+ i         *NZ*NP;
  for (i=0; i<NPERT; i++)                       udu0[i] = alloc_adjoint+   2*NPERT  * NZ * sizecontrolbc+(i+1*NPERT)*NZ*NP;
  for (i=0; i<NPERT; i++)                         u0[i] = alloc_adjoint+   2*NPERT  * NZ * sizecontrolbc+(i+2*NPERT)*NZ*NP;
  for (i=0; i<NPERT; i++)                  direction[i] = alloc_adjoint+   2*NPERT  * NZ * sizecontrolbc+(i+3*NPERT)*NZ*NP;	
  for (i=0; i<NPERT; i++)        lagrangian_gradient[i] = alloc_adjoint+   2*NPERT  * NZ * sizecontrolbc+(i+4*NPERT)*NZ*NP;	

 	

  
 //controlbc always zero.	
  Veclib::zero (sizecontrolbc*NZ*NPERT,  uc[0], 1);  
  getu(u0);
  
  for (j = 0; !converged && j < iterations; j++) {
    if (j>0){
      Veclib::copy (NZ*NP*NPERT, lagrangian_gradient[0], 1, direction[0], 1);
      
      //use duc as uc. first get the direction There are only gradient and conjugate methods here.
      if (j==1 || Femlib::ivalue  ("gradient_method")==0)
	Veclib::copy (NZ*NP*NPERT, lagrangian_gradient[0], 1, direction[0], 1);
      else if (Femlib::ivalue  ("gradient_method")==1)
	Veclib::svtvp (NZ*NP*NPERT, norml[j-1]/norml[j-2], direction[0], 1, lagrangian_gradient[0], 1, direction[0], 1); 
      
      installdomainu(direction);
      
      
      integrate (linAdvect, domain, analyst, adjoint_integration,uc);
      getu(udu0);
      //steplength is determined here.
      choose_step_length(uu0,udu0,u0, direction, normu0[j-1],normutau[j-1],steplength);
      
      
      Veclib::svtvp (NZ*NP*NPERT, steplength, direction[0], 1,u0[0], 1, u0[0], 1); 
      installdomainu(u0);
      normu0[j]=L2norm();		
      
      Veclib::svtvp (NZ*NP*NPERT, steplength,  udu0[0], 1,uu0[0], 1, uu0[0], 1); 
      installdomainu(uu0);
      normutau[j]=L2norm();	
      
      energygrowth[j]=normutau[j]/normu0[j];	
    } else {
      
      normu0[0]=L2norm();	
      integrate (linAdvect, domain, analyst, adjoint_integration,uc);
      normutau[0]=L2norm();
      getu(uu0);
      energygrowth[0]=normutau[0]/normu0[0];
      
      
    }
    //Scale u(T) as input of backward integration.
    for (i = 0; i < NPERT; i++)  	 *(domain -> u[i]) *=2/normu0[j];
    
    //backward integration to obtain the adjoint gradient term.
    integrate (linAdvectT, domain, analyst,adjoint_integration,uc); 
    norma[j]=L2norm();
    normu0_mixed[j]=L2norm_mixed(u0);
    getu(lagrangian_gradient);
    
    //Obtain the gradient of Lagrangian function with respect to the control force.
    Veclib::svtvp (NZ*NP*NPERT,  -2 * normutau[j] / normu0[j] / normu0[j], u0[0], 1, lagrangian_gradient[0], 1, lagrangian_gradient[0], 1); 
    
    
    // convergence criterions.
    installdomainu(lagrangian_gradient);
    norml[j]= L2norm();
    
    growthfile<< j<<setw(15)<<setprecision(8)<<steplength<< setw(15)<<energygrowth[j]<<setw(15)<<normu0_mixed[j]/sqrt(normu0[j]*norma[j])<<endl;
    if (steplength<EPSDP*EPSDP*EPSDP*EPSDP && j>0) {converged=1; printf("converged !\n"); growthfile<<"converged!"<<endl;} 
    installdomainu(u0);
    char          msg[StrMax], nom[StrMax];
    ofstream      file;
    sprintf   (msg, ".eig.%d", j);
    strcat    (strcpy (nom, domain -> name), msg);
    //scale u0 so that (u0,u0)=1.
    printf("%e     a\n", normu0[j]);
    for (i = 0; i < NPERT; i++)  	 *(domain -> u[i]) /=sqrt(normu0[j]);
    file.open (nom, ios::out); file << *domain; file.close();
  }
  
  growthfile.close();
  Femlib::finalize();
  return (EXIT_SUCCESS);
}

static void getargs (int        argc   ,
		     char**     argv   ,
		     real_t&     steplength   , 
			 int_t&     iterations,
		     real_t&    evtol  ,
		     char*&     session)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "dsa(-H) [options] session\n TG only."
    "options:\n"
    "-h       ... print this message\n"
    "-s       ... set step length\n"
	"-m       ... set number of iterations\n"
    "-t <num> ... set eigenvalue tolerance to num [Default 1e-6]\n";

  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 's':
      if (*++argv[0]) steplength  = strtod (  *argv,NULL);
      else { --argc;  steplength  = strtod (*++argv,NULL); }
      break;
	case 'm':
      if (*++argv[0]) iterations  = atoi (  *argv);
      else { --argc;  iterations  = atoi (*++argv); }
      break;
    case 't':
      if (*++argv[0]) evtol = atof (  *argv);
      else { --argc;  evtol = atof (*++argv); }
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


void getu(vector<real_t*> uu0)
//take u from domain and saved in uu0.
{  const int_t NZ = Geometry::nZ();
   const int_t NPERT = Geometry::nPert();
   const int_t NP = Geometry::planeSize();
	
	int_t i,j;
	for (i=0;i<NPERT;i++)
		for(j=0;j<NZ;j++)
			Veclib::copy (NP, domain -> u[i] -> plane (j), 1, uu0[i]+j*NP, 1);
}




void choose_step_length(vector<real_t*>uu0,
						vector<real_t*>udu0,
						vector<real_t*>u0, 
						vector<real_t*>direction, 
                        real_t  normu0_old,
                        real_t  normutau_old, 
						real_t&  steplength)
//choose a best steplength
{ const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();
	const int_t NP = Geometry::planeSize();
	int_t i,j;
	real_t  normu0_new, normutau_new,normutau_mixed, normu0_mixed,steplength1, steplength2,steplength_singular, a1, a2, a3, a4, G1, G2, G3, G4, G_singular;
	
	installdomainu    (direction);
	normu0_new =L2norm();
	normu0_mixed = L2norm_mixed (u0); 	
	
	installdomainu    (udu0);
	normutau_new =L2norm();
	normutau_mixed=L2norm_mixed(uu0);
	
  
    
	a1=2*normu0_mixed/normu0_new;
	a2=normu0_old/normu0_new;
	a3=2*normutau_mixed/normu0_new-2*normu0_mixed*normutau_new/normu0_new/normu0_new;
	a4=normutau_old/normu0_new-normutau_new*normu0_old/normu0_new/normu0_new;
	
	if (a4*a4/a3/a3+a2-a1*a4/a3<0){
		steplength1=-1;//remove it from competition.
		steplength2=-1;//remove it from competition.
	}
	else {
		steplength1=-a4/a3+sqrt(a4*a4/a3/a3+a2-a1*a4/a3);
		steplength2=-a4/a3-sqrt(a4*a4/a3/a3+a2-a1*a4/a3);		
	}

		

    if (steplength1 <0)
	  G1=-10; //remove it from competition.
	else
	  G1=(steplength1*steplength1*normutau_new+2*steplength1*normutau_mixed+normutau_old)/(normu0_old+steplength1*steplength1*normu0_new+2*steplength1*normu0_mixed);
	

	if (steplength2 <0)
		G2=-10; //remove it from competition.
	else
		G2=(steplength2*steplength2*normutau_new+2*steplength2*normutau_mixed+normutau_old)/(normu0_old+steplength2*steplength2*normu0_new+2*steplength2*normu0_mixed);
	

	G3=normutau_old/normu0_old;
	G4=normutau_new/normu0_new;
	
	//a3=0 is a singular point and has to be considered separately. 
    if (a3<EPSDP && a1<0)  
		steplength_singular=-a1/2;
	else 
		steplength_singular=0;
	G_singular = (steplength_singular*steplength_singular*normutau_new+2*steplength_singular*normutau_mixed+normutau_old)/(normu0_old+steplength_singular*steplength_singular*normu0_new+2*steplength_singular*normu0_mixed);
	
	if (G1 >= G2 && G1 >= G3 && G1 >= G4  &&  G1 >= G_singular) 
			steplength=steplength1; 
	else  if (G2 >= G1 && G2 >= G3 && G2 >= G4 && G2 >= G_singular) 
			steplength=steplength2;
	else  if (G_singular >= G1 && G_singular >= G2 &&G_singular >= G3 && G_singular >= G4)
			steplength=steplength_singular; 
	else  if (G4 >= G1 && G4 >= G2 &&G4 >= G3 && G4 >= G_singular)
		    steplength=1.0e50; //big enough when G4 is the largest
	else 
		steplength=0;
	
	 	printf("%e  %e  %e  %e  %e", G1, G2, G3, G4, G_singular);
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


real_t innerproduct (int_t size, 
				   real_t* v1,
				   real_t* v2)
// scalar value=v1*v2.
{int_t i; 
 real_t	value=0;
	for (i=0; i<size; i++)
		value+=*(v1+i) * *(v2+i);	
	
	return value;
}
