#include <stab.h>
static char             prog[] = "unsteady";
static char*            session;
static Domain*          domain;
static StabAnalyser*    analyst;
static vector<Element*> elmt;
static FEML*            file;
static Mesh*            mesh;
static BCmgr*           bman;

static void  getargs    (int, char**, real_t&, int_t&, real_t&, char*&, real_t&);
static int_t preprocess (const char*, bool&);
static real_t L2norm (); 
static real_t L2norm_mixed (vector<real_t*>, vector<real_t*>); 
static real_t norm_control (vector<real_t*>); 
static real_t norm_control_mixed (vector<real_t*>,vector<real_t*>); 
void dumpbc (int_t, real_t*, real_t*, char*);
void initializebc (int_t, real_t*, char*); 
void vvtom (int_t, real_t*, real_t*, real_t*);
real_t innerproduct (int_t, real_t*, real_t*);
void average_base(real_t,  real_t*);	
void normal_to_uvw (int_t, real_t*, real_t*, real_t*, real_t*);
void average_base (vector<real_t*>, real_t*, real_t*, real_t*);
void average_fields (vector<real_t*>, vector<real_t*>, vector<real_t*>,  real_t*, real_t*);
void choose_step_length (real_t*, real_t*, real_t&, const real_t, const real_t);
void readinmesh(real_t*, real_t*, real_t*);
void smooth_uc_gradient (int_t, vector<real_t*>);

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
  char      initialbc[StrMax],finalbc[StrMax],growthbc[StrMax];
  ofstream  growthfile;
  bool      restart = false;

  Femlib::initialize (&argc, &argv);

  Femlib::ivalue ("FLinterpolation", 2); //Lagrangian interpolation. switch off the Fourier interpolation.	

 getargs (argc, argv, steplengthcommand, iterations, evtol, session, Ec);
  
  strcat (strcpy (initialbc, session), ".rstbc");
  strcat (strcpy (finalbc, session), ".fldbc");
  strcat (strcpy (growthbc, session), ".growthbc");
  
  growthfile.open(growthbc);

//load base and restart files  
  const int_t ntot = preprocess (session, restart);
  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
  const int_t NP = Geometry::planeSize();
  const int_t  nStep    = Femlib::ivalue ("N_STEP");	
  if (NZ != 1) 	message (prog, "three-dimentional control is meaningless",   ERROR);

	
	
	//allocate memory for velocity fields
	vector<real_t*>  av_base;
	vector<real_t*>  av_pert;
	vector<real_t*>  av_grad;
	av_base.resize(NPERT);
	av_pert.resize(NPERT);
	av_grad.resize(NPERT);
	
   real_t*   alloc_average = new real_t [static_cast<size_t>(ntot*3)];
	for (i=0; i<NPERT; i++)        av_base[i] = alloc_average+ i          * NZ * NP;
	for (i=0; i<NPERT; i++)        av_pert[i] = alloc_average+(i+  NPERT) * NZ * NP;
    for (i=0; i<NPERT; i++)        av_grad[i] = alloc_average+(i+2*NPERT) * NZ * NP;
	
	real_t E_domain[8];
	real_t*  E_base     =E_domain; 
	real_t*  E_pert     =E_base+1; 
	real_t*  E_grad     =E_base+2;  
	real_t*  E_base_pert=E_base+3;  
	real_t*  E_base_grad=E_base+4;  
	real_t*  E_pert_grad=E_base+5;
    real_t*  E_steady   =E_base+6;
	real_t*  E_total    =E_base+7;

	real_t   Eav_domain[6];
	real_t*  Eav_base     =Eav_domain; 
	real_t*  Eav_pert     =Eav_base+1; 
	real_t*  Eav_grad     =Eav_base+2;  
	real_t*  Eav_base_pert=Eav_base+3;  
	real_t*  Eav_base_grad=Eav_base+4;  
	real_t*  Eav_pert_grad=Eav_base+5;

	
	
 // Allocate memory for adjoint_integration, uc and lagragian gradient.


  real_t unsteadiness[iterations];
  real_t normlag[iterations];

  
  int_t sizecontrolbc =  domain-> u[0] -> size_controlbc();
   
  vector<real_t*>  adjoint_integration;
  vector<real_t*>  lagrangian_gradient;
  vector<real_t*>  uc;
  vector<real_t*>  direction;
  vector<real_t*>  uc_next;
  vector<real_t*>  uuc;  //field u from uc
  vector<real_t*>  ugrad;  //field u from \delta uc
  
  adjoint_integration.resize(NPERT);
  lagrangian_gradient.resize(NPERT);
  uc.resize(NPERT);
  direction.resize(NPERT);
  uc_next.resize(NPERT);
  uuc.resize(NPERT);
  ugrad.resize(NPERT);

real_t*     alloc_adjoint = new real_t [static_cast<size_t>(NPERT*NZ*sizecontrolbc*5+2*NZ*NPERT*NP)];
  for (i=0; i<NPERT; i++)        adjoint_integration[i] = alloc_adjoint+ i          * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)        lagrangian_gradient[i] = alloc_adjoint+(i+  NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                         uc[i] = alloc_adjoint+(i+2*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                    uc_next[i] = alloc_adjoint+(i+3*NPERT) * NZ * sizecontrolbc;
  for (i=0; i<NPERT; i++)                  direction[i] = alloc_adjoint+(i+4*NPERT) * NZ * sizecontrolbc;	
															  
  for (i=0; i<NPERT; i++)                        uuc[i] = alloc_adjoint+   5*NPERT  * NZ * sizecontrolbc+ i       *NZ*NP;
  for (i=0; i<NPERT; i++)                       ugrad[i] = alloc_adjoint+   5*NPERT  * NZ * sizecontrolbc+(i+NPERT)*NZ*NP;
	
  real_t*  uc_n                 = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  adjoint_integration_n= new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  controlnx            = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];
  real_t*  controlny            = new  real_t [static_cast<size_t>(NZ*sizecontrolbc)];

  real_t*  meshx            = new  real_t [static_cast<size_t>(NP)];
  real_t*  meshy            = new  real_t [static_cast<size_t>(NP)];
  real_t*  massweight         = new  real_t [static_cast<size_t>(NP)];
	//readin control boundary normals
  domain -> u[0]->control_normal_direction(controlnx,controlny);
	
 //read in the whole mesh and mass matrix.
	readinmesh(meshx, meshy, massweight);
 // reading uc and normalize it 
 initializebc(sizecontrolbc, uc[0], initialbc);
//smooth uc
 smooth_uc_gradient (sizecontrolbc, uc);
 
if (Femlib::ivalue  ("NORMAL_CONTROL")){
    Veclib::vvtvvtp (NZ*sizecontrolbc,          uc[0], 1, controlnx, 1,          uc[1], 1, controlny, 1,          uc_n, 1);	//make sure u, v,w in uc is linear dependent.
    normal_to_uvw (sizecontrolbc, controlnx, controlny,          uc_n,          uc[0]);
 }

Veclib::smul (NPERT*NZ*sizecontrolbc, sqrt(Ec/norm_control(uc)), uc[0], 1, uc[0], 1); //normalize uc

 
//average the base flow to obtain av_base and the energy of base flow
average_base(av_base,E_base,Eav_base, massweight);	
	
	for (j = 0; !converged && j < iterations; j++) {
		if (j>0){
			// get the direction from the gradient
			if (j==1 || Femlib::ivalue  ("conjugate_gradient")==0)
				Veclib::copy (NZ*sizecontrolbc*NPERT, lagrangian_gradient[0], 1, direction[0], 1);
			else if (Femlib::ivalue  ("conjugate_gradient")==1) 
				Veclib::svtvp (NZ*sizecontrolbc*NPERT, normlag[j-1]/normlag[j-2], direction[0], 1, lagrangian_gradient[0], 1, direction[0], 1); 
			Veclib::svtvp (NZ*sizecontrolbc*NPERT, - norm_control_mixed (uc, direction) /Ec, uc[0], 1, direction[0], 1, direction[0], 1); // remove the uc component in direction.
            Veclib::smul (NZ*sizecontrolbc*NPERT, sqrt(Ec/norm_control(direction)), direction[0], 1, direction[0], 1); //normalize direction so as to simplify the step length function.
			
			//smooth the gradient
			smooth_uc_gradient (sizecontrolbc, direction);
			
			
			//integrate forwards the direction.
			integrate (linAdvect_unsteady, domain, analyst, adjoint_integration, direction, av_grad, E_domain, massweight, 2);
		    average_fields (av_base, av_pert, av_grad, Eav_domain, massweight);
			//choose steplength
			choose_step_length(E_domain, Eav_domain, steplength, norm_control_mixed(uc, direction)/Ec, norm_control_mixed(direction, direction)/Ec);
			if (steplengthcommand > EPSDP) steplength= steplengthcommand; // if steplength is given from command line, use it.		
        	// the nonlinear evolution may require extra forward integrations. The linear optimal steplength is gradually reduced.
			steplength=2*steplength;
			do{
			     steplength=steplength/2;
			     //update uc.
				Veclib::svtvp (sizecontrolbc*NZ*NPERT, steplength, direction[0], 1,  uc[0], 1, uc[0], 1);  
				smooth_uc_gradient (sizecontrolbc, uc); //smooth again
				Veclib::smul (NZ*sizecontrolbc*NPERT, sqrt(Ec/norm_control(uc)), uc[0], 1, uc[0], 1);
			
				//integrate forwards the updated uc.
				integrate (linAdvect_unsteady, domain, analyst, adjoint_integration,uc, av_pert, E_domain, massweight, 1);	
				average_fields (av_base, av_pert, av_grad, Eav_domain, massweight);
			
				*E_steady = *Eav_base + 2* *Eav_base_pert + *Eav_pert;
				*E_total  = *E_base   + 2* *E_base_pert   + *E_pert;
				unsteadiness[j]= *E_total -  *E_steady;
				growthfile<< j<<setw(15)<<steplength<< setprecision(15) << setw(25)<<unsteadiness[j]<<setw(25)<<unsteadiness[j]-unsteadiness[j-1] <<endl;
				
 }
			while(unsteadiness[j]>unsteadiness[j-1] && ! (steplengthcommand > EPSDP)  && steplength > 1.0e-10);
			
		} else {
			//forward integration. Start from zero initial condition.	
			integrate (linAdvect_unsteady, domain, analyst, adjoint_integration,uc, av_pert, E_domain, massweight, 1);	
		    average_fields (av_base, av_pert, av_grad,Eav_domain, massweight);
		    *E_steady = *Eav_base + 2* *Eav_base_pert + *Eav_pert;
		    *E_total  = *E_base   + 2* *E_base_pert   + *E_pert;
		    unsteadiness[j]= *E_total -  *E_steady;
			growthfile<< j<<setw(15)<<steplength<< setprecision(15) << setw(25)<<unsteadiness[j]<<setprecision(15) << setw(25)<<*E_base-*Eav_base<<endl;
		  }

		
		//backward integration to obtain the adjoint gradient term.
		Veclib::zero (sizecontrolbc*NZ*NPERT, adjoint_integration[0], 1);
		Veclib::zero (sizecontrolbc*NZ*NPERT,             uc_next[0], 1);
		//use av_grad to save av_pert + av_base.
		Veclib::svtvp (ntot, 1, av_pert[0], 1,  av_base[0], 1, av_grad[0], 1);  	
		integrate (linAdvectT_unsteady, domain, analyst,adjoint_integration,uc_next, av_grad, E_domain, massweight, 3); 
		
		if (Femlib::ivalue  ("NORMAL_CONTROL")){
		  // from uvw form to normal form.
		   Veclib::vvtvvtp (NZ*sizecontrolbc, adjoint_integration[0], 1, controlnx, 1, adjoint_integration[1], 1, controlny, 1, adjoint_integration_n, 1);
		  // spread adjoint_integration_n from normal version to uvw version.
		   normal_to_uvw (sizecontrolbc, controlnx, controlny, adjoint_integration_n, adjoint_integration[0]);	
		}
		
		//Obtain the gradient of Lagrangian function with respect to the control force.
		Veclib::smul (NZ*sizecontrolbc*NPERT, 1.0, adjoint_integration[0], 1, lagrangian_gradient[0], 1);
			
		normlag[j]= norm_control(lagrangian_gradient);
		
		//growthfile<< j<<setw(15)<<steplength<< setprecision(15) << setw(25)<<unsteadiness[j]<<setw(25)<<(j>0? unsteadiness[j]-unsteadiness[j-1] : unsteadiness[0])<<setw(15)<< setprecision(6) <<*E_base-*Eav_base<<setw(15)<<*E_pert-*Eav_pert<<endl;
		if (abs(steplength)<1.0e-10  && j>0) {converged=1; printf("converged !\n"); growthfile<<"converged!"<<endl;} 
		
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


static real_t L2norm_mixed (vector<real_t*> u1, vector<real_t*> u2)
// ---------------------------------------------------------------------------
// get the mixed L2norm. Not norm but norm^2
// ---------------------------------------------------------------------------
{ const int_t NC = Geometry::nPert();
	real_t      norm = 0.0;
	for (int_t i = 0; i < NC; i++) norm += domain -> u[i] -> mode_L2_mixed (u1[i],u2[i]);
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
  file << "x, y, un, u, v, real_imag"<<endl;
  for (i=0; i < size ; i++, uc_n++, uc0++)
   {file<< setw( 35) <<setprecision(20)<< *(meshx+i)  
    << setw(35) <<setprecision(20)<< *(meshy+i) <<" "; 
    for (j=0; j < NZ   ; j++)
	    file<< setw( 35)<<setprecision(20) <<*(uc_n+size*j); 
    for (k=0; k < NPERT; k++)
		for (j=0; j < NZ; j++)
		   file<< setw( 30)<<setprecision(20) <<*(uc0+size*(j+k*NZ)); 
	   file<< endl;
   }
  file.close();

uc0 -=size;
for (i=0;i<size*NZ*NPERT;i++,uc0++)
file1<<setprecision(20)<<*uc0 <<endl;
file1.close();
	 delete [] meshx;
	 delete [] meshy;
}



void initializebc(int_t sizecontrolbc,
       			  real_t* uc0,
			      char* initialbc)
//  initially the controlbc with zero if restart file session_initialbc does not exist.
{  const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();  
   const int_t       np     = Geometry::nP();
   const int_t nbound=sizecontrolbc/np;
   int i, j,k;
  const int_t NBASE = Geometry::nBase();
	
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




void	average_base(vector<real_t*>  av_base,
					        real_t*   E_base,
					        real_t*   Eav_base,
					        real_t*   massweight)
//calculate <U,U> and save in E_base. Average U and save in av_base. Note the size of av_base=NPERT*NZ*NP
{  const int_t  nStep    = Femlib::ivalue ("N_STEP");
   const int_t NP = Geometry::planeSize();
   const int_t NZ = Geometry::nZ();
   const int_t NPERT = Geometry::nPert();
   const real_t dt = Femlib::value ("D_T");
   const int_t NBASE = Geometry::nBase();
   int i, k, ntot = NP*NZ*NPERT;
	*E_base   = 0;
	*Eav_base = 0;
	Veclib::fill (ntot, 0 , av_base[0], 1);
	
	domain -> step = 0;
	domain -> time = 0.0;
		while (domain -> step < nStep) {
			domain -> step += 1; 
			domain -> time += dt;
			Femlib::value ("t", domain -> time);
			domain -> updateBase();
			for (i = 0; i < NBASE ; i++) {//note NBASE is always \le  NPERT
				Veclib::svtvp (NP, 1.0/nStep, domain -> U[i] -> plane (0), 1, av_base[i], 1, av_base[i], 1); //base flow always has only one plane.
				*E_base  += domain -> U[i] -> mode_L2_mixed_weight (domain -> U[i] -> plane (0), domain -> U[i] -> plane (0), massweight)/nStep;
			}
		}
   for (i = 0; i < NBASE ; i++)
	  *Eav_base += domain -> U[i] ->mode_L2_mixed_weight (av_base[i],av_base[i], massweight);
	//for (i=0; i<NPERT; i++) 			Veclib::smul (  NP*NZ, 1, av_base[i], 1, domain -> u[i] -> plane (0), 1); 
   //domain->dumppersecond();
		printf("unsteadiness of base flow is   %e  %e  %e\n ", *E_base,  *Eav_base, *E_base- *Eav_base);
}





void    average_fields(vector<real_t*> av_base, 
					   vector<real_t*> av_pert, 
					   vector<real_t*> av_grad, 
					          real_t*  Eav_domain,
                              real_t*  massweight)
//calculate the averaged values.
{  	const int_t NPERT = Geometry::nPert();
	int_t i;
	*(Eav_domain+1)=0;
	*(Eav_domain+2)=0;
	*(Eav_domain+3)=0;
	*(Eav_domain+4)=0;
	*(Eav_domain+5)=0;	
		
	for (i = 0; i < NPERT; i++) {
		*(Eav_domain+1) += domain -> u[i] -> mode_L2_mixed_weight (av_pert[i],av_pert[i], massweight);
		*(Eav_domain+2) += domain -> u[i] -> mode_L2_mixed_weight (av_grad[i],av_grad[i], massweight);
		*(Eav_domain+3) += domain -> u[i] -> mode_L2_mixed_weight (av_base[i],av_pert[i], massweight);
		*(Eav_domain+4) += domain -> u[i] -> mode_L2_mixed_weight (av_base[i],av_grad[i], massweight);	
		*(Eav_domain+5) += domain -> u[i] -> mode_L2_mixed_weight (av_grad[i],av_pert[i], massweight);
	}	
	
}







void 	choose_step_length(real_t*  E_domain, 
						   real_t*  Eav_domain, 
						   real_t&  steplength,
				const      real_t   norm_direction_uc,
				const	   real_t    norm_direction)
//choose a best steplength
{   const char  routine[] = "choose_step_length";
	int_t i, j, steps1=1000000;
	real_t scale, alpha, steps2=100000.0;
	real_t  E_base , E_pert, E_base_pert, Eav_base , Eav_pert, Eav_base_pert, E_steady, E_total; 
	real_t unsteadiness[steps1], min;

	for (i = 0; i < steps1; i++){
	  alpha = -i/steps2;
      scale = 1.0 / sqrt(1+alpha*alpha*norm_direction + 2*norm_direction_uc*alpha);
		E_base      =        *E_domain;
		E_pert      = scale*scale* (*(E_domain+1) + 2*alpha* *(E_domain+5) + alpha*alpha * *(E_domain+2));
		E_base_pert = scale*     *(E_domain+3) + scale*alpha**(E_domain+4);
		Eav_base = *Eav_domain;
		Eav_pert = scale*scale* (*(Eav_domain+1) + 2*alpha* *(Eav_domain+5) + alpha*alpha * *(Eav_domain+2));
		Eav_base_pert = scale* *(Eav_domain+3) + scale*alpha**(Eav_domain+4);
		
		E_steady = Eav_base + 2* Eav_base_pert + Eav_pert;
		E_total  = E_base   + 2* E_base_pert   + E_pert;
		unsteadiness[i]=  E_total - E_steady;
        
		if (i==0)  {min=unsteadiness[0]; steplength=0;}
		
		if (unsteadiness[i]<min) {min= unsteadiness[i]; 		steplength=alpha;}

	}
		
	if (abs(steplength)<EPSDP)    
		for (i = 0; i < steps1; i++){
			alpha = -i/steps2/steps1;
			scale = 1.0 / sqrt(1+alpha*alpha*norm_direction + 2*norm_direction_uc*alpha);
			E_base      =        *E_domain;
			E_pert      = scale*scale* (*(E_domain+1) + 2*alpha* *(E_domain+5) + alpha*alpha * *(E_domain+2));
			E_base_pert = scale*     *(E_domain+3) + scale*alpha**(E_domain+4);
			Eav_base = *Eav_domain;
			Eav_pert = scale*scale* (*(Eav_domain+1) + 2*alpha* *(Eav_domain+5) + alpha*alpha * *(Eav_domain+2));
			Eav_base_pert = scale* *(Eav_domain+3) + scale*alpha**(Eav_domain+4);
			
			E_steady = Eav_base + 2* Eav_base_pert + Eav_pert;
			E_total  = E_base   + 2* E_base_pert   + E_pert;
			unsteadiness[i]=  E_total - E_steady;
			
			if ( unsteadiness[i]<min) {min= unsteadiness[i]; 		steplength=alpha;}
		}	
		
			
	printf("(in linear framework) steplength =%e ,  unsteadiness[i-1]= %e  , normdirection = %e \n ", steplength, min, norm_direction);
}




void readinmesh(real_t* meshx, 
				real_t* meshy, 
				real_t* massweight)
// readin the whole mesh in meshx and meshy, the mass matrix in mesh mass, and filter this matrix.
{ const int_t NP = Geometry::planeSize();
	domain -> u[0] -> meshxymass(meshx, meshy);
	int_t i; 
	real_t radius, criterion=10; 
	//massweight has to be real and positive. It is a square of another vector.
		//Veclib::fill (NP, 1 , massweight, 1);
	for (i=0;i<NP;i++){
		radius= sqrt(*(meshx+i) * *(meshx+i)  + *(meshy+i) * *(meshy+i));
		if ( *(meshx+i) < criterion) *(massweight+i)=1;
	   else                     *(massweight+i)=exp(-(*(meshx+i)-criterion) * (*(meshx+i)-criterion));
     }
}



void smooth_uc_gradient (int_t sizecontrolbc,
						 vector<real_t*>  tosmooth)
//to smooth uc or the gradient so that the same node has the same value.
{ const int_t NZ = Geometry::nZ();
  const int_t NPERT = Geometry::nPert();
	int_t i,j,k,h;
	real_t* meshx;
	real_t* meshy;
	meshx= new real_t [sizecontrolbc];
	meshy= new real_t [sizecontrolbc];
	
	domain -> u[0] -> controlmesh(meshx,meshy);
	
	for (i = 0; i < NPERT; i++)
		for (j = 0; j < NZ; j++)
			for (k=0; k<sizecontrolbc; k++)
				for (h=k+1; h<sizecontrolbc; h++)
					if (abs(*(meshx+k) - *(meshx+h))<EPSDP && abs(*(meshy+k) - *(meshy+h))<EPSDP) {
					 *(tosmooth[i]+j * sizecontrolbc + h) = *(tosmooth[i]+j * sizecontrolbc + k);
						//printf("%d  %d \n ", k,h); 
					}
//	printf("%e   %e\n ", *(tosmooth[1]+5), *(tosmooth[1]+6)); 
	
	delete [] meshx;
	delete [] meshy;	
	
}