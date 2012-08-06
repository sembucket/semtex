#include <stab.h>
static char             prog[] = "optbc";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void   getargs            (int, char**, real_t&, int_t&, real_t&, char*&, bool*);
static int_t  preprocess         (const char*, bool&);
static void   Zerodomain         ();
static real_t L2norm             (); 
static real_t L2norm_mixed       (vector<real_t*>); 
static real_t norm_control       (vector<real_t*>); 
static real_t norm_control_mixed (vector<real_t*>,vector<real_t*>); 
       void   dumpbc             (int_t, real_t*, char*, real_t);
       void   initializebc       (int_t, real_t*, char*); 
       void   choose_step_length (vector<real_t*>, vector<real_t*>,
				  vector<real_t*>,vector<real_t*>,
				  int_t, real_t,real_t, real_t&);
       void   getu               (vector<real_t*>);
       void   bfgs_Hessian       (int_t, real_t*, real_t*, real_t*);
       void   vvtom              (int_t, real_t*, real_t*, real_t*);
       real_t innerproduct       (int_t, real_t*, real_t*);
       void   controlmesh        (real_t*, real_t*);


       void   bccheck            (int_t, real_t*, char*);
       void   check_adjoint_integration (int_t, real_t*, char*);

int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  int_t     converged = 0,iterations=0;
  real_t    steplengthcommand = 0, steplength, resnorm, normlbefore, evtol = 1.0e-16, ro=0.1, c1const=0.000000000000001;
  int_t     i, j, k;
  char      buf[StrMax];
  char      initialbc[StrMax],finalbc[StrMax],growthbc[StrMax], checkbc[StrMax], checkadj[StrMax];
  ofstream  growthfile;
  bool      restart = false, LNS = false;
  bool*     LNS_P = &LNS;
  
  

  Femlib::ivalue ("gradient_method", 0);  //default is steepest method. using conjugate_gradient = 1 to switch to CG.
  Femlib::ivalue ("TIMEDECAY", 0);        //default is the double Gaussian time decay factor.
  Femlib::initialize (&argc, &argv);
	
  getargs (argc, argv, steplengthcommand, iterations, evtol, session, LNS_P);
  
  strcat (strcpy (initialbc, session), ".rstbc");
  strcat (strcpy (finalbc, session), ".fldbc");
  strcat (strcpy (growthbc, session), ".growthbc");
  strcpy (checkbc, session);
  strcpy (checkadj, session);
  
  if (LNS==0){
  growthfile.open(growthbc);
  }

  //load base and restart files  
  const int_t ntot = preprocess (session, restart);
 
  //Allocate memory for adjoint_integration, uc and lagragian gradient.
  const int_t NZ    = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP    = Geometry::planeSize();
  real_t      energygrowth[iterations];
  real_t      normc[iterations];
  real_t      normc_mixed[iterations];
  real_t      normv[iterations];
  real_t      norml[iterations];
  real_t      normg[iterations];

  //  int_t cnodes = domain -> u[0] -> size_controlbc();       //Alternative name for sizecontrolbc
  int_t sizecontrolbc = bman -> sizecontrolbc();               //number of control boundary nodes for one field in one plane
  int_t cnodetot      = NPERT*NZ*sizecontrolbc;                //total number of control boundary nodes in problem
   
  vector<real_t*>  adjoint_integration;
  vector<real_t*>  lagrangian_gradient;
  vector<real_t*>  lagrangian_gradient_old;
  vector<real_t*>  direction;
  vector<real_t*>  lagrangian_gradient_change;
  vector<real_t*>  uc;
  vector<real_t*>  uc_next;
  vector<real_t*>  uc_change; 
  vector<real_t*>  uuc;                         //field u from uc
  vector<real_t*>  uduc;                        //field u from \delta uc

  adjoint_integration.resize(NPERT);
  lagrangian_gradient.resize(NPERT);
  lagrangian_gradient_old.resize(NPERT);
  direction.resize(NPERT);
  lagrangian_gradient_change.resize(NPERT);
  uc.resize(NPERT);
  uc_next.resize(NPERT);
  uc_change.resize(NPERT);
  uuc.resize(NPERT);
  uduc.resize(NPERT);
  
  //NOTE:Allocate a block of memory to store all field and gradient data? [SEAN]
  real_t*  alloc_adjoint = new real_t [static_cast<size_t>(cnodetot*8+2*NZ*NPERT*NP)];

  for (i=0; i<NPERT; i++)        adjoint_integration[i] = alloc_adjoint+ i          * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)        lagrangian_gradient[i] = alloc_adjoint+(i+NPERT)   * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)    lagrangian_gradient_old[i] = alloc_adjoint+(i+2*NPERT) * NZ * sizecontrolbc; 
  for (i=0; i<NPERT; i++)                  direction[i] = alloc_adjoint+(i+3*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++) lagrangian_gradient_change[i] = alloc_adjoint+(i+4*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                         uc[i] = alloc_adjoint+(i+5*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                    uc_next[i] = alloc_adjoint+(i+6*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                  uc_change[i] = alloc_adjoint+(i+7*NPERT) * NZ * sizecontrolbc;	
  for (i=0; i<NPERT; i++)                        uuc[i] = alloc_adjoint+8*cnodetot+i*NZ*NP;
  for (i=0; i<NPERT; i++)                       uduc[i] = alloc_adjoint+8*cnodetot+(i+NPERT)*NZ*NP;
  
  //initialize the hessian matrix to be an identity matrix
  real_t*  Hessian = new  real_t [static_cast<size_t>(cnodetot*cnodetot)];
  Veclib::zero (cnodetot*cnodetot, Hessian, 1);  
  for (j=0; j< cnodetot;j++)
    *(Hessian+ j*cnodetot + j)=1;
 

  //set the control force for v and w zero to avoid singularity in cylindrical condition??
  initializebc(sizecontrolbc,uc[0],initialbc);  

  //Sean Loh edit: Check original BC is valid
  bccheck(sizecontrolbc, uc[0],checkbc);

  for (j = 0; !converged && j < iterations; j++) {
    if (j>0){
      //use duc as uc. first get the direction
      if (j==1 || Femlib::ivalue("gradient_method")==0)
	Veclib::copy (cnodetot, lagrangian_gradient[0], 1, direction[0], 1);
      
      else if (Femlib::ivalue("gradient_method")==1)
	Veclib::svtvp (cnodetot, norml[j-1]/norml[j-2], direction[0], 1, lagrangian_gradient[0], 1, direction[0], 1); 
      
      else if (Femlib::ivalue("gradient_method")==2){
	Veclib::vsub (cnodetot, lagrangian_gradient[0] , 1, lagrangian_gradient_old[0], 1, lagrangian_gradient_change[0], 1); 
	Veclib::smul (cnodetot, steplength, direction[0], 1, uc_change[0], 1); 
	//scale the original Hessian.
	if (j==2 && innerproduct (cnodetot, uc_change[0], lagrangian_gradient_change[0]) <0)
	  Veclib::smul (cnodetot*cnodetot, -innerproduct (cnodetot, uc_change[0], lagrangian_gradient_change[0])/innerproduct (cnodetot, lagrangian_gradient_change[0], lagrangian_gradient_change[0]) , Hessian, 1, Hessian, 1); 
	if (j==2 && innerproduct (cnodetot, uc_change[0], lagrangian_gradient_change[0]) >0)
	  Veclib::smul (cnodetot*cnodetot,  innerproduct (cnodetot, uc_change[0], lagrangian_gradient_change[0])/innerproduct (cnodetot, lagrangian_gradient_change[0], lagrangian_gradient_change[0]) , Hessian, 1, Hessian, 1); 
	
	
	bfgs_Hessian (cnodetot, uc_change[0],  lagrangian_gradient_change[0], Hessian);
	Blas::mxv    (Hessian, cnodetot, lagrangian_gradient[0], cnodetot, direction[0]);
      }
		
      //forward integration. Start from zero initial condition.
      // NOTE: I think this forward integration gets u field induced by the change in uc from the last iteration, that's why "direction" is parsed to integrate not "uc". May also explain why there are steps to update uc, the "plane" data for each field and "uuc" folowing this integration. This seems to be a strange way of doing things which begs the question, why?!
      Zerodomain ();
      integrate (linAdvect, domain, analyst, adjoint_integration, direction);
      getu(uduc);

      //steplength is determined here.
      choose_step_length(uuc,uduc,uc, direction, sizecontrolbc, normc[j-1],normv[j-1],steplength);
      if (steplengthcommand > EPSDP) steplength = steplengthcommand;                                 //if steplength is given from command line, use it.
      
      for (i=0; i<NPERT; i++) 
	for(k=0;k<NZ;k++){
	  Veclib::svtvp (sizecontrolbc, steplength, direction[i]+sizecontrolbc*k, 1, uc[i]+sizecontrolbc*k, 1, uc[i]+sizecontrolbc*k, 1);      //update uc.
	  Veclib::svtvp (NP, steplength, uduc[i]+NP*k, 1, uuc[i]+NP*k, 1, domain -> u[i] -> plane (k), 1);                                     //update domain -> u
	  Veclib::svtvp (NP, steplength, uduc[i]+NP*k, 1, uuc[i]+NP*k, 1, uuc[i]+NP*k, 1);                                                     //update uuc.
	}
      
      normc[j]=norm_control(uc);
      normv[j]=L2norm();
      energygrowth[j]=normv[j]/normc[j];
      
    } 
    
    else {
      //forward integration. Start from zero initial condition.
      Zerodomain ();
      
      integrate (linAdvect, domain, analyst, adjoint_integration,uc);
      getu(uuc);
      
      normc[0]=norm_control(uc);
      normv[0]=L2norm();
      energygrowth[0]=normv[0]/normc[0];
      
      // If using LNS to verify gain of optimal perturbation output energy gain and end iteration by setting converged=1.      
      if (LNS){
	converged=1;
	printf("LNS evolution complete!\n");
	printf("Kinetic energy gain at final time is: ");
	cout << energygrowth[0] << endl;
	break;
      }
      
      if (steplengthcommand > EPSDP) steplength= steplengthcommand; // does not affect any results. only for clarification.

    }

    //Scale u(T) as input of backward integration.
    for (i = 0; i < NPERT; i++)
      *(domain -> u[i]) *=2/normc[j];
    
    //backward integration to obtain the adjoint gradient term.
    Veclib::zero (cnodetot, adjoint_integration[0], 1);
    Veclib::zero (cnodetot, uc_next[0], 1);
    integrate (linAdvectT, domain, analyst,adjoint_integration,uc_next); 

    //Sean's checking routine
    check_adjoint_integration (sizecontrolbc, adjoint_integration[0], checkadj);

    //Obtain the gradient of the Lagrangian functional with respect to the control force.
    //firstly save the previous gradient
    if (j>0)
      Veclib::copy (cnodetot, lagrangian_gradient[0], 1, lagrangian_gradient_old[0], 1);
    
    Veclib::smul  (cnodetot, -2 * normv[j] / normc[j] / normc[j], uc[0], 1, lagrangian_gradient[0], 1);
    Veclib::svtvp (cnodetot, 1.000, adjoint_integration[0], 1, lagrangian_gradient[0], 1, lagrangian_gradient[0], 1); 
        
    //project lagrangian to the uc normal (tangent?) direction.
    normc_mixed[j] = norm_control_mixed (uc, adjoint_integration); 
    normg[j]       = norm_control(adjoint_integration);

    normlbefore = norm_control(lagrangian_gradient);	

    // NOTE: What does this step do? What does the norm_control_mixed/normc term physically represent? Does this just remove any component of the gradient that is not due to a change in the uc direction?

    Veclib::svtvp (cnodetot, - norm_control_mixed(uc, lagrangian_gradient)/normc[j], uc[0], 1, lagrangian_gradient[0], 1, lagrangian_gradient[0], 1); 
    norml[j] = norm_control(lagrangian_gradient);
    
    growthfile << j << setw(15)
	       << steplength << setw(15)
	       << energygrowth[j] << setw(15)
	       << normc_mixed[j]/2 << setw(15)
	       << energygrowth[j]-(j>0? energygrowth[j-1] : 0) << setw(15)
	       << norml[j]*steplength << setw(15)
	       << normc_mixed[j]/ sqrt(normg[j]*normc[j]) << endl;
    
    if (steplength < EPSDP*EPSDP*EPSDP*EPSDP && j>0) {
      converged = 1;
      printf("converged!\n");
      growthfile << "converged!" << endl;
    } 
    
    //output the control velocities.
    dumpbc(sizecontrolbc,uc[0],finalbc, normc[j]);
    
  
  }


  growthfile.close();
  Femlib::finalize();
  delete[] alloc_adjoint;
  delete[] Hessian;
  return (EXIT_SUCCESS);
}

static void getargs (int        argc   ,
		     char**     argv   ,
		     real_t&    steplength   , 
		     int_t&     iterations,
		     real_t&    evtol  ,
		     char*&     session,
		     bool*      LNS)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "optbc [options] session\n TG only."
    "options:\n"
    "-h       ... print this message\n"
    "-s       ... set step length\n"
    "-m       ... set number of iterations\n"
    "-t <num> ... set eigenvalue tolerance to num [Default 1e-6]\n"
    "-g       ... evolve perturbation to final time using LNS and compute gain\n";

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
    case 'g':
      *LNS = true;
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
{ const int_t NC   = Geometry::nPert();
  real_t      norm = 0.0;
  for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> mode_L2_mixed (domain -> u[i] -> plane (0),anotheru[i]);

  return norm;
}



static real_t norm_control (vector<real_t*> uc)
// ---------------------------------------------------------------------------
// get the time-independent energy of the control force (from control bc). it's not norm, but norm^2
// ---------------------------------------------------------------------------
{ const int_t NC   = Geometry::nPert();
  real_t      norm = 0.0;
  for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> normc (uc[i]);

  return norm;
}

static real_t norm_control_mixed (vector<real_t*> uc, vector<real_t*> lag)
// ---------------------------------------------------------------------------
// get the time-independent energy of the control force (from control bc and Lag). it's not norm, but norm^2:[uc,lag].
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
  real_t      norm = 0.0;
  for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> normc_mixed (uc[i],lag[i]);
  
  return norm;
}

 
 
void  dumpbc( int_t   size,
              real_t* uc0,
	      char*   finalbc,
              real_t  normc)
//dump the control bc to session_finalbc.
{
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();

  char     name[StrMax];
  ofstream file;
  ofstream file1;

  strcpy (name,finalbc);
  strcat (name,"optimalbc");
  file.open(name);
  file1.open(finalbc);
  
  int_t   i,j,k;
  real_t* meshx;
  real_t* meshy;
 
  meshx = new real_t [size];
  meshy = new real_t [size];

  //  domain -> u[0] -> controlmesh(meshx,meshy);
  controlmesh(meshx, meshy);

  file << "x, y, realu, imagu, realv,imagv,realw,imagw"<<endl;

  for (i = 0; i < size ; i++, uc0++)
    {
      file << setprecision(20) << *(meshx+i) << setw( 25) << setprecision(20) << *(meshy+i) << " "; 
      
      for (k=0; k < NPERT; k++)
	for (j=0; j < NZ   ; j++)
	  file << setw(30) << setprecision(20) << *(uc0+size*(j+k*NZ)) / sqrt(normc); 
      
      file << endl;
    }
  
  file.close();
  
  uc0 -=size;
  for (i=0;i<size*NZ*NPERT;i++,uc0++)
    {
      file1 << setprecision(20) << *uc0 / sqrt(normc) << endl;
    }
  file1.close();
}



void initializebc(int_t   sizeplane,
		  real_t* uc0,
		  char*   initialbc)
// Randomly initialize the controlbc if restart file session.rstbc does not exist.
{
   const int_t NZ      = Geometry::nZ();
   const int_t NPERT   = Geometry::nPert();  
   const int_t np      = Geometry::nP();
   const int_t nel     = sizeplane/np;
   const int_t NBASE   = Geometry::nBase();
   const int_t step    = 0;
   int_t       problem = Geometry::problem();

   vector<real_t> work;
   real_t         *uc, *gvector;
   int            i,j,k,nbound;

   ifstream file (initialbc);
   
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
       for(k=0;k<nel-1;k++)
	 *(uc0+k*np+ sizeplane*(i*NZ+j))= *(uc0+k*np+2*np-1+ sizeplane*(i*NZ+j));

   // -- Modify uc to ensure agreement with other BCs.
   for (i = 0; i < NPERT; i++){         
     // Write BC's to _line storage.
     domain -> u[i] -> evaluateControl (step, uc0+i*sizeplane*NZ);
     domain -> u[i] -> evaluateBoundaries(step);
   }

   // Couple v-w BC data in _line storage if required.
   if(Geometry::cylindrical() && 
      (problem != Geometry::O2_2D || problem != Geometry::SO2_2D)){
     Field::coupleBCs (domain -> u[1], domain -> u[2], FORWARD);
   }

   // Modify control BC data in line storage to accommodate other BCs.
   for (i = 0; i < NPERT; i++){         
     domain -> u[i] -> fixControl ();
   }
  
   // Uncouple v-w BC data in _line storage if required.
   if(Geometry::cylindrical() && 
      (problem != Geometry::O2_2D || problem != Geometry::SO2_2D)){
     Field::coupleBCs (domain -> u[1], domain -> u[2], INVERSE);
   }

   // Retrieve control velocity data from _line storage.
   for (i = 0; i < NPERT; i++){         
     domain -> u[i] -> getControl (uc0+i*sizeplane*NZ);
   }   
}



void getu(vector<real_t*> uuc)
//take u from domain and save in uuc.
{
  const int_t NZ    = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP    = Geometry::planeSize();
  
  int_t i,j;
  for (i=0;i<NPERT;i++)
    for(j=0;j<NZ;j++)
      Veclib::copy (NP, domain -> u[i] -> plane (j), 1, uuc[i]+j*NP, 1);
}




void choose_step_length(vector<real_t*>uuc,
			vector<real_t*>uduc,
			vector<real_t*>uc, 
			vector<real_t*>direction, 
                        int_t          sizecontrolbc,
                        real_t         normc_old,
                        real_t         normv_old, 
			real_t&        steplength)
//choose a best steplength
{ const int_t  NZ    = Geometry::nZ();
  const int_t  NPERT = Geometry::nPert();
  const int_t  NP    = Geometry::planeSize();
	int_t  i,j;
	real_t normc_new, normv_new,normv_mixed, normc_mixed,steplength1, steplength2,steplength_singular, a1, a2, a3, a4, G1, G2, G3, G4, G_singular;
	
	normc_new = norm_control(direction);
	normv_new = L2norm();
	normv_mixed = L2norm_mixed(uuc);
	normc_mixed = norm_control_mixed (uc, direction); 	
  
    
	a1=2*normc_mixed/normc_new;
	a2=normc_old/normc_new;
	a3=2*normv_mixed/normc_new-2*normc_mixed*normv_new/normc_new/normc_new;
	a4=normv_old/normc_new-normv_new*normc_old/normc_new/normc_new;
	
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
	  G1=(steplength1*steplength1*normv_new+2*steplength1*normv_mixed+normv_old)/(normc_old+steplength1*steplength1*normc_new+2*steplength1*normc_mixed);
	
	
	if (steplength2 <0)
	  G2=-10; //remove it from competition.
	else
	  G2=(steplength2*steplength2*normv_new+2*steplength2*normv_mixed+normv_old)/(normc_old+steplength2*steplength2*normc_new+2*steplength2*normc_mixed);
	

	G3=normv_old/normc_old;
	G4=normv_new/normc_new;
	
	//a3=0 is a singular point and has to be considered separately. 
	if (a3<EPSDP && a1<0)  
	  steplength_singular=-a1/2;
	else 
	  steplength_singular=0;
	G_singular = (steplength_singular*steplength_singular*normv_new+2*steplength_singular*normv_mixed+normv_old)/(normc_old+steplength_singular*steplength_singular*normc_new+2*steplength_singular*normc_mixed);
	
	if (G1 >= G2 && G1 >= G3 && G1 >= G4  &&  G1 >= G_singular) 
	  steplength=steplength1; 
	else  if (G2 >= G1 && G2 >= G3 && G2 >= G4 && G2 >= G_singular) 
	  steplength=steplength2;
	else  if (G_singular >= G1 && G_singular >= G2 &&G_singular >= G3 && G_singular >= G4)
	  steplength=steplength_singular; 
	else  if (G4 >= G1 && G4 >= G2 &&G4 >= G3 && G4 >= G_singular)
	  steplength=0.1; // not sure how to set it
	else 
	  steplength=0;
	
	//printf("%e  %e  %e  %e  %e", G1, G2, G3, G4, G_singular);
}


void bfgs_Hessian (int_t size, 
				   real_t* duc, 
				   real_t* dgradient, 
				   real_t* Hessian)
//calculate direction for BFGS .
{ int_t i,j,k;
	real_t ro=0, a;
	real_t* m1= new  real_t [static_cast<size_t>(size*size)];
	real_t* m2= new  real_t [static_cast<size_t>(size*size)];
	real_t* m3= new  real_t [static_cast<size_t>(size*size)];
	real_t* minter= new  real_t [static_cast<size_t>(size*size)]; //intermidate 
	real_t* vinter= new  real_t [static_cast<size_t>(size)]; //intermidate 
	real_t*      I= new  real_t [static_cast<size_t>(size*size)];  	
	//identity matrix
	Veclib::zero (size*size, I, 1);  
	for (j=0; j< size;j++)
		*(I+ j*size + j)=1;
	
	//get ro
	ro= innerproduct (size, duc, dgradient);
	ro=1.0/ro;// printf("      ro=       %e\n\n", ro);
	if (ro>0 || fabs(ro)< EPSDP ) return;
	
	vvtom (size, duc, dgradient, minter);
	Veclib::svtvp (size*size, -ro, minter, 1, I, 1, m1, 1); 
	
	vvtom (size, dgradient, duc, minter);
	Veclib::svtvp (size*size, -ro, minter, 1, I, 1, m2, 1); 
	
	vvtom (size, duc, duc, m3);
	
	Blas::mxm (m1, size, Hessian, size, minter, size);
	Blas::mxm (minter, size, m2, size, Hessian, size);
	Veclib::svtvp (size*size, -ro, m3, 1, Hessian, 1, Hessian, 1);
	
	
	delete [] m1;
	delete [] m2;
	delete [] m3;
	delete [] minter;
	delete [] I;
}


void vvtom (int_t size, 
            real_t* v1,
			real_t* v2,
            real_t* m)
// vector v1^T * vector v2 =matrix m
{int_t i;
	for (i=0; i<size; i++)
		Veclib::smul (size, *(v1+i), v2, 1, m+i*size, 1); 	
}

real_t innerproduct (int_t size, 
				real_t* v1,
				real_t* v2)
	// scalar value=v1*v2.
	{int_t i; 
	 real_t value=0;
		for (i=0; i<size; i++)
			value+=*(v1+i) * *(v2+i);	
		return value;	
	}

//Sean Loh edit
void bccheck (int_t  size,
              real_t* uc0,
	      char*   checkbc)
//dump the control bc to session.checkbc
{
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();

  char     name[StrMax];
  ofstream file;
  ofstream file1;

  strcpy (name,checkbc);
  strcat (name,".checkbc");
  file.open(name);
  
  int_t   i,j,k;
  real_t* meshx;
  real_t* meshy;
 
  meshx = new real_t [size];
  meshy = new real_t [size];

  controlmesh(meshx,meshy);

  file << "x, y, realu, imagu, realv,imagv,realw,imagw"<<endl;

  for (i = 0; i < size ; i++, uc0++)
    {
      file << setprecision(20) << *(meshx+i) << setw( 25) << setprecision(20) << *(meshy+i) << " "; 
      
      for (k=0; k < NPERT; k++)
	for (j=0; j < NZ   ; j++)
	  file << setw(30) << setprecision(20) << *(uc0+size*(j+k*NZ)); 

      file << endl;
    }

  file.close();
  
}

void check_adjoint_integration(int_t   size,
			       real_t* adj0,
			       char*   checkadj)
//dump adjoint_integration to session.checkadj
{
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();

  char name[StrMax];
  ofstream file;

  strcpy (name, checkadj);
  strcat (name, ".checkadj");
  file.open(name);

  for (int i = 0; i < NPERT*NZ*size; i++)
    {
      file << *(adj0+i) << endl;
    }
  file.close();
}

void controlmesh(real_t* controlx,real_t* controly)
// Install the control force in uci.                                                          
{
  const int_t       np = Geometry::nP();
  int_t       i,k, nbound;
  vector<Boundary*> BC;

  if (Geometry::nPert() == 2) BC = domain -> b[0] -> BCs (0);
  else                        BC = domain -> b[0] -> BCs (Femlib::ivalue ("BETA"));

  nbound = domain -> b[0] -> nSurf();

  for (i = 0; i < nbound; i++ )
    if(*BC[i] -> group()=='c') {
      BC[i] -> controlbcmesh(controlx, controly);
      controlx += np;
      controly += np;
    }
}
