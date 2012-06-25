#include <stab.h>
static char             prog[] = "icbc";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void  getargs    (int, char**, real_t&, int_t&, real_t&, char*&, real_t&);
static int_t preprocess (const char*, bool&);
static void Zerodomain ();
static real_t L2norm (); 
static real_t L2norm_mixed (vector<real_t*>); 
static real_t norm_control (vector<real_t*>); 
static real_t norm_control_mixed (vector<real_t*>,vector<real_t*>); 
void dumpbc (int_t, real_t*, real_t*, char*);
void initializebc (int_t, real_t*, char*, real_t*, char*); 
void choose_step_length(vector<real_t*>, vector<real_t*>, vector<real_t*>,vector<real_t*>,vector<real_t*>,int_t, real_t,real_t, real_t&);
void getu(vector<real_t*>);
void vvtom (int_t, real_t*, real_t*, real_t*);
real_t innerproduct (int_t, real_t*, real_t*);
void normal_to_uvw (int_t, real_t*, real_t*, real_t*, real_t*);
int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver routine for stability analysis code.
// ---------------------------------------------------------------------------
{
  int_t     converged = 0,iterations=0;
  real_t   steplengthcommand=0, steplength, resnorm, evtol = 1.0e-16, Ec=1,scale_factor;
  int_t     i, j, k;
  char      buf[StrMax];
  char      initialbc[StrMax],finalbc[StrMax],growthbc[StrMax],tocontrol[StrMax];
  ofstream  growthfile;
  bool      restart = false;
  Femlib::ivalue ("conjugate_gradient", 0); //default is steepest method. using conjugate_gradient= 1 to switch to CG.
  Femlib::ivalue ("TIMEDECAY", 0);   //default is the double Gaussian time decay factor.
 Femlib::ivalue ("MAXMIZATION ", 0); // 1 for maximization problem and 0 for minimization problem only for icbc.

	Femlib::initialize (&argc, &argv);

  getargs (argc, argv, steplengthcommand, iterations, evtol, session, Ec);
  
  strcat (strcpy (initialbc, session), ".rstbc");
  strcat (strcpy (finalbc, session), ".fldbc");
  strcat (strcpy (growthbc, session), ".growthbc");
  strcat (strcpy (tocontrol, session), ".tocontrol");
  
  growthfile.open(growthbc);

//load base and restart files  
  const int_t ntot = preprocess (session, restart);
 
 // Allocate memory for adjoint_integration, uc and lagragian gradient.
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  real_t energygrowth[iterations];
  real_t normc[iterations];
  real_t normc_mixed[iterations];
  real_t normv[iterations];
  real_t normg[iterations];
  real_t norml[iterations];

  
  int_t sizecontrolbc =  domain-> u[0] -> size_controlbc();
   
  vector<real_t*>  adjoint_integration;
  vector<real_t*>  lagrangian_gradient;
  vector<real_t*>  uc;
  vector<real_t*>  direction;
  vector<real_t*>  uc_next;
  vector<real_t*>  uuc;  //field u from uc
  vector<real_t*>  uduc;  //field u from \delta uc
  vector<real_t*>  uctocontrol;  //u0 to control
  
  adjoint_integration.resize(NPERT);
  lagrangian_gradient.resize(NPERT);
  uc.resize(NPERT);
  direction.resize(NPERT);
  uc_next.resize(NPERT);
  uuc.resize(NPERT);
  uduc.resize(NPERT);
  uctocontrol.resize(NPERT);  
   real_t*     alloc_adjoint = new real_t [static_cast<size_t>(NPERT*NZ*sizecontrolbc*6+2*NZ*NPERT*NP)];
  for (i=0; i<NPERT; i++)        adjoint_integration[i] = alloc_adjoint+ i          * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)        lagrangian_gradient[i] = alloc_adjoint+(i+  NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                         uc[i] = alloc_adjoint+(i+2*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                    uc_next[i] = alloc_adjoint+(i+3*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                uctocontrol[i] = alloc_adjoint+(i+4*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                  direction[i] = alloc_adjoint+(i+5*NPERT) * NZ * sizecontrolbc;	
															  
  for (i=0; i<NPERT; i++)                        uuc[i] = alloc_adjoint+   6*NPERT  * NZ * sizecontrolbc+ i       *NZ*NP;
  for (i=0; i<NPERT; i++)                       uduc[i] = alloc_adjoint+   6*NPERT  * NZ * sizecontrolbc+(i+NPERT)*NZ*NP;
	
  real_t*  uc_n                 = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  adjoint_integration_n= new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  controlnx            = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  controlny            = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  uctocontrol_n        = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
//readin control boundary normals
	
	domain -> u[0]->control_normal_direction(controlnx,controlny);
	
 // readin uc_n, uctocontrol, fill uc from uc_n and fill uctocontrol from uctocontrol_n
 initializebc(sizecontrolbc, uc[0], initialbc, uctocontrol[0], tocontrol);
 Veclib::vvtvvtp (NZ*sizecontrolbc,          uc[0], 1, controlnx, 1,          uc[1], 1, controlny, 1,          uc_n, 1);	//make sure u, v,w in uc is linear dependent.
 Veclib::vvtvvtp (NZ*sizecontrolbc, uctocontrol[0], 1, controlnx, 1, uctocontrol[1], 1, controlny, 1, uctocontrol_n, 1);	//squash and then expand.
 normal_to_uvw (sizecontrolbc, controlnx, controlny,          uc_n,          uc[0]);
 normal_to_uvw (sizecontrolbc, controlnx, controlny, uctocontrol_n, uctocontrol[0]);

	//normalize uc
  Veclib::smul (NPERT*NZ*sizecontrolbc, sqrt(Ec/norm_control(uc)), uc[0], 1, uc[0], 1);

    for (j = 0; !converged && j < iterations; j++) {
		if (j>0){
	 // get the direction from the gradient
		  if (j==1 || Femlib::ivalue  ("conjugate_gradient")==0)
			Veclib::copy (NZ*sizecontrolbc*NPERT, lagrangian_gradient[0], 1, direction[0], 1);
		  else if (Femlib::ivalue  ("conjugate_gradient")==1) {
			Veclib::svtvp (NZ*sizecontrolbc*NPERT, norml[j-1]/norml[j-2], direction[0], 1, lagrangian_gradient[0], 1, direction[0], 1); 
            Veclib::svtvp (NZ*sizecontrolbc*NPERT, - norm_control_mixed (uc, direction) /Ec, uc[0], 1, direction[0], 1, direction[0], 1); // remove the uc component in direction.
		  }
	
        //forward integration. Start from zero initial condition.
	    Zerodomain ();
		integrate (linAdvect, domain, analyst, adjoint_integration,direction);
	    getu(uduc);
		//steplength is determined here.
			
		choose_step_length(uctocontrol, uuc, uduc, uc, direction, sizecontrolbc, Ec,normv[j-1], steplength);
		if (steplengthcommand > EPSDP) steplength= steplengthcommand; // if steplength is given from command line, use it.		

		
        for (i=0; i<NPERT; i++) 
			for(k=0;k<NZ;k++){
				Veclib::svtvp (sizecontrolbc, steplength, direction[i]+sizecontrolbc*k, 1,  uc[i]+sizecontrolbc*k, 1, uc[i]+sizecontrolbc*k      , 1);  //update uc.
				Veclib::svtvp (           NP, steplength,      uduc[i]+           NP*k, 1, uuc[i]+           NP*k, 1, domain -> u[i] -> plane (k), 1); //update domain ->u
				Veclib::svtvp (           NP, steplength,      uduc[i]+           NP*k, 1, uuc[i]+           NP*k, 1,      uuc[i]+           NP*k, 1);  //update uuc.
			}
      //scale uc domain->u and uuc:
		scale_factor = sqrt(Ec/norm_control(uc));
		for (i=0; i<NPERT; i++) 
			for(k=0;k<NZ;k++){
				Veclib::smul (sizecontrolbc, scale_factor, uc[i]+sizecontrolbc*k      , 1, uc[i]+sizecontrolbc*k      , 1);  //scale uc.
				Veclib::smul (           NP, scale_factor, domain -> u[i] -> plane (k), 1, domain -> u[i] -> plane (k), 1);  //scale domain -> u.
				Veclib::smul (           NP, scale_factor, uuc[i] + NP*k              , 1, uuc[i] + NP*k              , 1);  //scale uuc.
			}
		normv[j]=L2norm();
		energygrowth[j]=2*norm_control_mixed(uctocontrol, uc)+ normv[j];
		
    } else {
	//forward integration. Start from zero initial condition.
	Zerodomain ();

	integrate (linAdvect, domain, analyst, adjoint_integration,uc);
	getu(uuc);

	normv[0]=L2norm();
	energygrowth[0]=2*norm_control_mixed(uctocontrol, uc)+ normv[0];
		
	}

	
	//backward integration to obtain the adjoint gradient term.
	Veclib::zero (sizecontrolbc*NZ*NPERT, adjoint_integration[0], 1);
    Veclib::zero (sizecontrolbc*NZ*NPERT,             uc_next[0], 1);
	integrate (linAdvectT, domain, analyst,adjoint_integration,uc_next); 
	
	// from uvw form to normal form.
	Veclib::vvtvvtp (NZ*sizecontrolbc, adjoint_integration[0], 1, controlnx, 1, adjoint_integration[1], 1, controlny, 1, adjoint_integration_n, 1);
	// spread adjoint_integration_n from normal version to uvw version.
	normal_to_uvw (sizecontrolbc, controlnx, controlny, adjoint_integration_n, adjoint_integration[0]);	
		
		
    //Obtain the gradient of Lagrangian function with respect to the control force.
	Veclib::svtvp (NZ*sizecontrolbc*NPERT, 1, uctocontrol[0], 1, adjoint_integration[0], 1, lagrangian_gradient[0], 1); 	
	Veclib::svtvp (NZ*sizecontrolbc*NPERT, - norm_control_mixed (uc, lagrangian_gradient) /Ec, uc[0], 1, lagrangian_gradient[0], 1, lagrangian_gradient[0], 1); 
	Veclib::svtvp (NZ*sizecontrolbc*NPERT, - norm_control_mixed (uc, lagrangian_gradient) /Ec, uc[0], 1, lagrangian_gradient[0], 1, lagrangian_gradient[0], 1); 
	//testing::: to scale gradient so as to reduce the divergence free energy
	Veclib::smul (NZ*sizecontrolbc*NPERT, 1.0/sqrt(norm_control(lagrangian_gradient)), lagrangian_gradient[0], 1, lagrangian_gradient[0], 1);
		
	norml[j]= norm_control(lagrangian_gradient);
	
			
	growthfile<< j<<setw(15)<<steplength<< setw(15)<<energygrowth[j]<<setw(15)<<energygrowth[j]-(j>0? energygrowth[j-1] : 0)<<endl;
		if (abs(steplength)<EPSDP*EPSDP*EPSDP*EPSDP*EPSDP*EPSDP  && j>0) {converged=1; printf("converged !\n"); growthfile<<"converged!"<<endl;} 

    //output the control velocities.
Veclib::vvtvvtp (NZ*sizecontrolbc,                  uc[0], 1, controlnx, 1,                  uc[1], 1, controlny, 1,                  uc_n, 1);
		dumpbc(sizecontrolbc,uc_n, uc[0], finalbc);
	   
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
		     char*&     session ,
			 real_t&     Ec  )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "dsa(-H) [options] session\n TG only."
    "options:\n"
    "-h       ... print this message\n"
    "-s       ... set step length\n"
	"-m       ... set number of iterations\n"
	"-e       ... Ec=[uc,uc]."
    "-t <num> ... set tolerance (not used in the code) to num [Default 1e-6]\n";

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
	case 'e':
			if (*++argv[0]) Ec  = atof (  *argv);
			else { --argc;  Ec  = atof (*++argv); }
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
{ const int_t NC = Geometry::nPert();
	real_t      norm = 0.0;
	for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> mode_L2_mixed (domain -> u[i] -> plane (0),anotheru[i]);
	return norm;
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

 
 
 void  dumpbc(int_t size,
              real_t* uc_n,
			   real_t* uc0,
			  char* finalbc)
 //dump the control bc to session_finalbc.
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
  file << "x, y, real (un) imag (un)"<<endl;
  for (i=0; i < size ; i++, uc_n++)
   {file<<setprecision(20)<< *(meshx+i)  
         << setw( 25) <<setprecision(20)<< *(meshy+i) <<" "; 
        for (j=0; j < NZ   ; j++)
	file<< setw( 30)<<setprecision(20) <<*(uc_n+size*j); 
    file<< endl;
   }
  file.close();

for (i=0;i<size*NZ*NPERT;i++,uc0++)
file1<<setprecision(20)<<*uc0 <<endl;
file1.close();

}


void initializebc(int_t sizecontrolbc,
       			  real_t* uc0,
			      char* initialbc,
				  real_t* uctocontrol,
				  char* tocontrol)
// randomly initially the controlbc if restart file session_initialbc does not exist.
{  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();  
   const int_t       np     = Geometry::nP();
   const int_t nbound=sizecontrolbc/np;
   int i, j,k;
  const int_t NBASE = Geometry::nBase();
//read in to control files
	ifstream filetocontrol (tocontrol);
    if (filetocontrol)  {
	
		while (! filetocontrol.eof() ){
			filetocontrol >> *(uctocontrol++); //which is the output of lnsbc -a ...(in u,v,w form)
		}
		filetocontrol.close();
    }
	else
		Veclib::fill (sizecontrolbc*NZ*NPERT, 0 , uctocontrol, 1);
	
// read in initial guess of optimal bc.
	ifstream fileuc (initialbc);
    if (fileuc)  {
     while (! fileuc.eof() ){
       fileuc >> *(uc0++);
	 }
     fileuc.close();
    }
    else
	Veclib::fill (sizecontrolbc*NZ*NPERT, 1 , uc0, 1);
	

	
  }
  


void getu(vector<real_t*> uuc)
//take u from domain and saved in uuc.
{  const int_t NZ = Geometry::nZ();
   const int_t NPERT = Geometry::nPert();
   const int_t NP = Geometry::planeSize();
	
	int_t i,j;
	for (i=0;i<NPERT;i++)
		for(j=0;j<NZ;j++)
			Veclib::copy (NP, domain -> u[i] -> plane (j), 1, uuc[i]+j*NP, 1);
}




void choose_step_length(vector<real_t*>u0,
						vector<real_t*>uuc,
						vector<real_t*>uduc,
						vector<real_t*>uc, 
						vector<real_t*>direction, 
                        int_t sizecontrolbc,
                        real_t  Ec,
                        real_t  normv_old, 
						real_t&  steplength)
//choose a best steplength
{   const char  routine[] = "choose_step_length";
	const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();
	const int_t NP = Geometry::planeSize();
	int_t N, i,j, info, ILO, IHI;
	real_t   a1, a2, a3, a4, a5, a6, alpha, L, maxL;
	real_t*  c = new  real_t [static_cast<size_t>(5)];
	bool problem=Femlib::ivalue ("MAXMIZATION");

	a1 = norm_control(direction)/Ec;
	a2 = 2*norm_control_mixed(u0,direction);
	a3 = normv_old;
	a4 = 2*L2norm_mixed(uuc);
	a5 = L2norm();
	a6 = 2* norm_control_mixed(u0,uc);
	
	*(c+0)=a4*a4-a2*a2;
	*(c+1)=-4*a1*a3*a4+4*a4*a5+2*a1*a2*a6;
	*(c+2)=(2*a1*a3-2*a5)*(2*a1*a3-2*a5) - 2*a1*a4*a4-a2*a2*a1 - a1*a1*a6*a6;
	*(c+3)=2*a1*(2*a1*a3*a4-2*a4*a5+a1*a2*a6);
	*(c+4)=a1*a1*(a4*a4-a1*a6*a6);


	N=4;
	//solve the roots of polynomial as eigenvalues
	real_t*     H = new real_t [static_cast<size_t>(N*N)];
	real_t*     WR = new real_t [static_cast<size_t>(N)];
	real_t*     WI = new real_t [static_cast<size_t>(N)];
	real_t*     Z = new real_t [static_cast<size_t>(N*N)];
	real_t*    WORK = new real_t [static_cast<size_t>(N)];
    real_t*    SCALE = new real_t [static_cast<size_t>(N)];
	
	//build matrix H.
	Veclib::fill (N*N, 0 , H, 1);
    for (i=0;i<N-1;i++)		*(H+i*N+i+1)=1;
	for (i=0;i<N;i++)		 *(H+(N-1)*N+i)=-*(c+i)/ *(c+N);

     //solve eigenvalues		
    Lapack::gebal("B", N, H, N, ILO, IHI, SCALE, info);
	Lapack::hseqr("E","N", N, ILO, IHI, H, N, WR, WI, Z, N, WORK, N, info);
  
	//choose the optimal alpha among WR which maximize L, 6 possibilities for steplength, 0, WR[0]..WR[5] and 10^50.
	maxL=a6+a3;
	steplength =0;
	for (i=0; i<N; i++)
		{ alpha =  *(WR+i); 
		  L = ((a6+a2*alpha)*sqrt(1+a1*alpha*alpha) +a3 + a4*alpha + a5*alpha*alpha)/ (1+a1*alpha*alpha);
			printf("%e      %e \n", alpha,   L);
	      if ((L > maxL && problem && alpha >0) || (L < maxL && !problem && alpha<0)) {
		    maxL = L; steplength = alpha;
		  }
		}
	
	if ((maxL > (a2*sqrt(a1) +a5) /a1  &&  !problem)  ||  (maxL < (a2*sqrt(a1) +a5) /a1  &&  problem)) steplength =  1e50; 
    if ((maxL > (-a2*sqrt(a1)+a5) /a1  &&  !problem)  ||  (maxL < (-a2*sqrt(a1)+a5) /a1  &&  problem)) steplength = -1e50;
//	printf("%e   %e \n", steplength, 1/sqrt(1+norm_control(direction)/Ec*steplength*steplength));
printf("%e      %e  %e  \n", steplength,   (-a2*sqrt(a1)+a5) /a1, ((a6+a2*steplength)*sqrt(1+a1*steplength*steplength) +a3 + a4*steplength + a5*steplength*steplength)/ (1+a1*steplength*steplength));
	
 delete [] H; delete [] WR; delete [] WI; delete [] Z;  delete [] WORK;  
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




void normal_to_uvw ( int_t   sizecontrolbc,
					 real_t* controlnx, 
					 real_t* controlny, 
					 real_t* uc_n, 
					 real_t*  uc0)
//converge the normal velocity component to u, v, w forms.
{  const int_t NZ = Geometry::nZ();
	const int_t NPERT = Geometry::nPert();
	Veclib::fill (sizecontrolbc*NZ*NPERT, 0 , uc0, 1);
	Veclib::vmul (sizecontrolbc*NZ, uc_n, 1, controlnx, 1, uc0, 1);
	Veclib::vmul (sizecontrolbc*NZ, uc_n, 1, controlny, 1, uc0+sizecontrolbc*NZ, 1);
}

