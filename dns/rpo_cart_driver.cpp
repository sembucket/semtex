//////////////////////////////////////////////////////////////////////////////
// drive.C: control spectral element DNS for incompressible flows.
//
// Copyright (c) 1994 <--> $Date$, Hugh Blackburn
//
// USAGE:
// -----
// dns [options] session
//   options:
//   -h       ... print usage prompt
//   -i       ... use iterative solver for viscous [and pressure] steps
//   -t[t]    ... select time-varying BCs for mode 0 [or all modes]
//   -v[v...] ... increase verbosity level
//   -chk     ... turn off checkpoint field dumps [default: selected]
//   -S|C|N   ... regular skew-symm || convective || Stokes advection
//
// AUTHOR:
// ------
// Hugh Blackburn
// Department of Mechanical & Aerospace Engineering
// Monash University
// Vic 3800
// Australia
// hugh.blackburn@monash.edu
//
// --
// This file is part of Semtex.
//
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA.
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id$";

#include <dns.h>

#include "rpo_cart_utils.h"

#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>

static char prog[] = "rpo";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
void integrate (void (*)(Domain*, BCmgr*, AuxField**, AuxField**, FieldForce*),
		Domain*, BCmgr*, DNSAnalyser*, FieldForce*);

#define PIPE
#define TRAVELLING_WAVE

#define NFIELD 3
#define YMIN (0.0)
#define YMAX (1.0)
#define ZMAX (0.7391982714328925) // L = 2.pi/beta
#define NELS_X 12
#define NELS_Y 9

#define XMIN (-1.0)
#define XMAX (+1.0)

void build_constraints(Context* context, Vec x_delta, double* f_phi, double* f_tau) {
  int          elOrd       = Geometry::nP() - 1;
  int          nNodesX     = NELS_X*elOrd;
  int          index;
  int          nDofsCube_l = Geometry::nZProc() * context->nDofsPlane;
  int          nl          = 3 * nDofsCube_l;
  int          plane_j;
  int          pt_j, el_j;
  double       p_y;
  double       k_z;
  double       f_phi_l, f_tau_l;
  double*      rz          = new double[nl];
  double*      rt          = new double[nl];
  double*      data_r      = new double[context->nDofsPlane];
  double*      data_i      = new double[context->nDofsPlane];
  PetscScalar* xArray;
  Vec          xl;
  int          nStep;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->global_to_semtex, x_delta, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArray(xl, &xArray);

  if(!context->travelling_wave) {
    for(int field_i = 0; field_i < 3; field_i++) {
      *context->domain->u[field_i] = *context->u0[field_i];
    }
    if(!Geometry::procID()) cout << "time step in constraints evaluation: " << scientific << Femlib::value("D_T") << endl;

    delete context->analyst;
    context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);

    nStep = Femlib::ivalue("N_STEP");
    Femlib::ivalue("N_STEP", 1);
    integrate(convective, context->domain, context->bman, context->analyst, context->ff);
    Femlib::ivalue("N_STEP", nStep);

    for(int field_i = 0; field_i < 3; field_i++) {
      *context->domain->u[field_i] -= *context->u0[field_i];
      *context->domain->u[field_i] *= 1.0 / Femlib::value("D_T");
    }
  } 

  for(int dof_i = 0; dof_i < nl; dof_i++) { rz[dof_i] = rt[dof_i] = 0.0; }

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
      elements_to_logical(context->u0[field_i]->plane(plane_i+0), data_r);
      elements_to_logical(context->u0[field_i]->plane(plane_i+1), data_i);
      for(int dof_i = 0; dof_i < NELS_Y * elOrd * nNodesX; dof_i++) {
//#ifdef PIPE
//        el_j = dof_i / (nNodesX * elOrd);
//        pt_j = (dof_i / nNodesX) % elOrd;
//        p_y  = context->elmt[el_j]->_ymesh[pt_j*(elOrd+1)];
//        if(fabs(p_y) < 1.0e-6) p_y = 1.0;
        // assume a radius of 1, so scale by (2*pi) / (2*pi*r)
        //k_z  = (1.0 / p_y) * context->phi_i * (plane_j / 2);
//        k_z  = 1.0 * (plane_j / 2);
//#else
        //k_z  = (2.0 * M_PI / ZMAX) * context->phi_i * (plane_j / 2);
        //k_z  = (1.0 / ZMAX) * (plane_j / 2);
//        k_z  = (2.0 * M_PI / ZMAX) * (plane_j / 2);
//#endif
        k_z  = (2.0 * M_PI / ZMAX) * (plane_j / 2);

        index = field_i * nDofsCube_l + (plane_i+0) * context->nDofsPlane + dof_i;
        rz[index] = -k_z * data_i[dof_i];
        index = field_i * nDofsCube_l + (plane_i+1) * context->nDofsPlane + dof_i;
        rz[index] = +k_z * data_r[dof_i];
      }
    }
    if(!context->travelling_wave) {
      for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
        elements_to_logical(context->domain->u[field_i]->plane(plane_i), data_r);
        for(int dof_i = 0; dof_i < context->nDofsPlane; dof_i++) {
          index = field_i * nDofsCube_l + plane_i * context->nDofsPlane + dof_i;
          rt[index] = data_r[dof_i];
        }
      }
    }
  }

  f_phi_l = f_tau_l = 0.0;
  for(int dof_i = 0; dof_i < nl; dof_i++) {
    f_phi_l -= rz[dof_i] * xArray[dof_i];
    if(!context->travelling_wave) f_tau_l -= rt[dof_i] * xArray[dof_i];
  }
  MPI_Allreduce(&f_phi_l, f_phi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&f_tau_l, f_tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecRestoreArray(xl, &xArray);

  delete[] rz;
  delete[] rt;
  delete[] data_r;
  delete[] data_i;
  VecDestroy(&xl);
}

PetscErrorCode _snes_jacobian(SNES snes, Vec x, Mat J, Mat P, void* ctx) {
  Context* context = (Context*)ctx;

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  J, MAT_FINAL_ASSEMBLY);

  return 0;
}

PetscErrorCode _snes_function(SNES snes, Vec x, Vec f, void* ctx) {
  Context* context = (Context*)ctx;
  real_t dt;
  int field_i, mode_i, mode_j, dof_i, nStep;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  real_t f_norm, x_norm, dx_norm, dx_norm_2;
  double f_phi, f_tau, dummy[2];
  Vec dx_test;
  double runTime;
  bool update = false;
  char filename[100];
  PetscViewer viewer;

  // create the \delta x vector for use in the assembly of the constraints into the residual vector
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &dx_test);
  VecZeroEntries(dx_test);
  VecAXPY(dx_test, +1.0, x);
  VecAXPY(dx_test, -1.0, context->x_prev);
  VecNorm(dx_test, NORM_2, &dx_norm);
  VecNorm(x, NORM_2, &x_norm);

  if(context->iteration == 0) update = true;
  if(dx_norm/x_norm > 1.0e-6) update = true;

  if(!Geometry::procID()) cout << "\tupdate test; |dx|: " << dx_norm << "\t|x|: " << x_norm << "\t|dx|/|x|: " << dx_norm/x_norm << endl;
  VecNorm(context->x_delta, NORM_2, &dx_norm_2);
  if(!Geometry::procID()) cout << "\tactual |dx|: " << dx_norm_2 << endl;
  if(update) {
    if(!Geometry::procID()) cout << "\tunpacking solution vector" << endl;
    VecCopy(dx_test, context->x_delta);

    UnpackX(context, context->u0, &dummy[0], &dummy[1], x);

    // write the current state vector
    sprintf(filename, "x_curr_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(x, viewer);
    PetscViewerDestroy(&viewer);
    // write the previous state vector
    sprintf(filename, "x_prev_%.4u.vec", context->iteration);
    PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    VecView(context->x_prev, viewer);
    PetscViewerDestroy(&viewer);
  } else {
    if(!Geometry::procID()) cout << "\t existing solution vector" << endl;
  }
  UnpackX(context, context->ui, &context->phi_i, &context->tau_i, x);

  // update the starting time
  Femlib::value("t", 0.0);
  context->domain->step = 0;

  // initialise the flow map fields with the solution fields
  for(field_i = 0; field_i < 3; field_i++) {
    *context->domain->u[field_i] = *context->ui[field_i];
  }

  // solve the flow map for time tau_i
  MPI_Bcast(&context->tau_i, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&context->phi_i, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  dt = context->tau_i / Femlib::ivalue("N_STEP");
  Femlib::value("D_T", dt);
  runTime = Femlib::ivalue("N_STEP")*dt;

  if(!Geometry::procID()) {
    cout << "run time: " << runTime << "\tnstep: " << Femlib::ivalue("N_STEP") << "\tdt: " << Femlib::value("D_T") << endl;
    cout << scientific << "\ttau:   " << context->tau_i << "\tphi:   " << context->phi_i * (2.0*M_PI/ZMAX) << endl;
  }

  // don't want to call the dns analysis, use custom integrate routine instead
  delete context->analyst;
  context->analyst = new DNSAnalyser (context->domain, context->bman, context->file);
  integrate(convective, context->domain, context->bman, context->analyst, context->ff);

  // phase shift in fourier direction
  phase_shift_z(context, context->phi_i * (2.0*M_PI/ZMAX), 1.0, context->domain->u);

  // set f
  for(field_i = 0; field_i < 3; field_i++) {
    *context->fi[field_i]  = *context->domain->u[field_i];
    *context->fi[field_i] -= *context->ui[field_i];
  }

  build_constraints(context, context->x_delta, &f_phi, &f_tau);

  RepackX(context, context->fi, f_phi, f_tau, f);

  for(field_i = 0; field_i < 3; field_i++) *context->write->u[field_i] = *context->ui[field_i];
  *context->write->u[3] = *context->domain->u[3]; //dump the pressure as a sanity check
  context->write->dump();

  VecNorm(x, NORM_2, &x_norm);
  VecNorm(f, NORM_2, &f_norm);
  context->iteration++;
  if(!Geometry::procID()) cout << context->iteration << ":\tevaluating function, |x|: " << x_norm << "\t|f|: " << f_norm << endl;

  VecCopy(x, context->x_prev);

  VecDestroy(&dx_test);

  return 0;
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<Field*> ui, vector<Field*> fi, vector<Field*> u0)
{
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = NELS_X*elOrd;
  Vec x, f;
  Mat P;
  KSP ksp;
  SNES snes;
  FEML* file_i;
  char session_i[100];

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

  context->mesh     = mesh;
  context->elmt     = elmt;
  context->domain   = domain;
  context->bman     = bman;
  context->file     = file;
  context->analyst  = analyst;
  context->ff       = FF;
  context->ui       = ui;
  context->fi       = fi;
  context->u0       = u0;
  context->xmin     = XMIN;
  context->xmax     = XMAX;
  context->phi_i    = 0.0;
  context->tau_i    = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  context->iteration = 0;
#ifdef TRAVELLING_WAVE
  context->travelling_wave = true;
#else
  context->travelling_wave = false;
#endif

  sprintf(session_i, "cav01_rpo_1");
  file_i = new FEML(session_i);
  context->write = new Domain(file_i, context->elmt, context->bman);
  delete file_i;

  // add dofs for phi and tau for each time slice
  context->nDofsPlane = NELS_X*elOrd*NELS_Y*elOrd;
  context->nDofsSlice = 3 * Geometry::nZ() * context->nDofsPlane + 2;
  context->localSize  = 3 * Geometry::nZProc() * context->nDofsPlane;
  if(!Geometry::procID()) context->localSize += 2;

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &f);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &context->x_prev);

  MatCreate(MPI_COMM_WORLD, &P);
  MatSetType(P, MATMPIAIJ);
  MatSetSizes(P, context->localSize, context->localSize, context->nDofsSlice, context->nDofsSlice);
  //MatMPIAIJSetPreallocation(P, 2*nNodesX, PETSC_NULL, 2*nNodesX, PETSC_NULL);
  MatMPIAIJSetPreallocation(P, 1, PETSC_NULL, 1, PETSC_NULL);
  MatSetOptionsPrefix(P, "P_");
  MatSetFromOptions(P);

  SNESCreate(MPI_COMM_WORLD, &snes);
  SNESSetFunction(snes, f,    _snes_function, (void*)context);
  SNESSetJacobian(snes, P, P, _snes_jacobian, (void*)context);
  SNESGetKSP(snes, &ksp);
  KSPSetType(ksp, KSPGMRES);
  SNESSetType(snes, SNESNEWTONTR);
  SNESSetNPCSide(snes, PC_LEFT);
  SNESSetFromOptions(snes);

  context->snes = snes;
  context->x_norm = 1.0;//1000.0;
  RepackX(context, context->ui, context->phi_i, context->tau_i, x);
  VecNorm(x, NORM_2, &context->x_norm);
  if(!Geometry::procID()) cout << "|x_0|: " << context->x_norm << endl;
  if(context->x_norm < 1.0e-4) {
    if(!Geometry::procID()) cout << "ERROR: initial state vector norm is SMALL! "
                                 << "Are you sure you loaded the initial condition correctly??\n";
    abort();
  }
  context->x_norm = 1.0;//1000.0;
  VecCopy(x, context->x_prev);
  SNESSolve(snes, NULL, x);
  UnpackX(context, context->ui, &context->phi_i, &context->tau_i, x);

  if(!Geometry::procID()) cout << "rpo solve complete.\n";
  if(!Geometry::procID()) cout << "\tshift phi:   " << context->phi_i * (2.0*M_PI/ZMAX) << endl;
  if(!Geometry::procID()) cout << "\tshift tau:   " << context->tau_i << endl;

  VecDestroy(&x);
  VecDestroy(&f);
  MatDestroy(&P);
  VecDestroy(&context->x_prev);
}

int main (int argc, char** argv) {
#ifdef _GNU_SOURCE
  feenableexcept (FE_OVERFLOW);    // -- Force SIG8 crash on FP overflow.
#endif
  char*            session;
  bool             freeze = false;
  vector<Element*> elmt;
  FEML*            file;
  Mesh*            mesh;
  BCmgr*           bman;
  Domain*          domain;
  DNSAnalyser*     analyst;
  FieldForce*      FF;
  static char      help[] = "petsc";
  vector<Field*>   ui;  // Solution fields for velocities, pressure at the i time slices
  vector<Field*>   fi;  // Solution fields for flow maps at the i time slices
  vector<Field*>   u0;  // Additional fields for building the constraints
  char             session_i[100];
  FEML*            file_i;
  char*            fname;
  BoundarySys*     bndry;

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  //domain -> restart ();
  //ROOTONLY domain -> report ();
  
  // load in the time slices
  ui.resize(3);
  fi.resize(3);
  u0.resize(3);
  delete file;
  delete domain;
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();
  for(int field_i = 0; field_i < NFIELD; field_i++) {
    ui[field_i] = domain->u[field_i];

    strcpy ((fname = new char [strlen (bman -> field()) + 1]), bman -> field());

    bndry = new BoundarySys(bman, elmt, fname[0]);
    fi[field_i] = new Field(bndry, new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, fname[0]);

    bndry = new BoundarySys(bman, elmt, fname[0]);
    u0[field_i] = new Field(bndry, new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, fname[0]);
  }
  delete file_i;
  delete domain;

  // allocate the temporary fields
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0);
  delete file_i;
  delete domain;

  // dump the output
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  for(int field_i = 0; field_i < 3; field_i++) {
    *domain->u[field_i] = *ui[field_i];
  }
  domain->dump();
  delete file_i;
  delete domain;

  PetscFinalize();

  return EXIT_SUCCESS;
}

static void getargs (int    argc   ,
		     char** argv   ,
		     bool&  freeze ,
		     char*& session)
// ---------------------------------------------------------------------------
// Install default parameters and options, parse command-line for optional
// arguments.  Last argument is name of a session file, not dealt with here.
// ---------------------------------------------------------------------------
{
  char       buf[StrMax];
  const char routine[] = "getargs";
  const char usage[]   = "Usage: %s [options] session-file\n"
    "  [options]:\n"
    "  -h       ... print this message\n"
    "  -f       ... freeze velocity field (scalar advection/diffusion only)\n"
    "  -i       ... use iterative solver for viscous steps\n"
    "  -v[v...] ... increase verbosity level\n"
    "  -chk     ... turn off checkpoint field dumps [default: selected]\n"
    "  -S|C|N   ... regular skew-symm || convective || Stokes advection\n";

  Femlib::ivalue ("ADVECTION", 1); // -- Default is alternating skew symmetric.

  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 'S': Femlib::ivalue ("ADVECTION", 0); break;
    case 'C': Femlib::ivalue ("ADVECTION", 2); break;
    case 'N': Femlib::ivalue ("ADVECTION", 3); break;
    case 'f':
      freeze = true;
      break;
    case 'i':
      do			// -- Only allowing ITERATIVE=1 (Viscous).
	Femlib::ivalue ("ITERATIVE", 1);
      while (*++argv[0] == 'i');
      break;
    case 'v':
      do
	Femlib::ivalue ("VERBOSE",   Femlib::ivalue ("VERBOSE")   + 1);
      while (*++argv[0] == 'v');
      break;
    case 'c':
      if (strstr ("chk", *argv))     Femlib::ivalue ("CHKPOINT",    0);
      else { fprintf (stdout, usage, prog); exit (EXIT_FAILURE); }
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  //if   (argc != 1) message (routine, "no session definition file", ERROR);
  //else             session = *argv;
  session = *argv;

  Femlib::value ("DTBDX", 0.0);
}

static void preprocess (const char*       session,
			FEML*&            file   ,
			Mesh*&            mesh   ,
			vector<Element*>& elmt   ,
			BCmgr*&           bman   ,
			Domain*&          domain ,
			FieldForce*&      FF     )
// ---------------------------------------------------------------------------
// Create objects needed for execution, given the session file name.
// They are listed above in order of creation.
// ---------------------------------------------------------------------------
{

  const char routine[] = "preprocess";
  const int_t        verbose = Femlib::ivalue ("VERBOSE");
  Geometry::CoordSys space;
  int_t              i, np, nz, nel, procid, seed;

  // -- Initialise problem and set up mesh geometry.
  VERBOSE cout << "Building mesh ..." << endl;
  file = new FEML (session);
  mesh = new Mesh (file);
  VERBOSE cout << "done" << endl;

  // -- Set up global geometry variables.
  VERBOSE cout << "Setting geometry ... ";
  nel   =  mesh -> nEl();
  np    =  Femlib::ivalue ("N_P");
  nz    =  Femlib::ivalue ("N_Z");
  space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (np, nz, nel, space);
  VERBOSE cout << "done" << endl;

  // -- If token RANSEED > 0 then initialize the random number
  //    generator based on wall clock time and process ID (i.e. a "truly"
  //    pseudo-random number).  NB: it is important to have done this
  //    before any other possible call to random number routines.

  if (Femlib::ivalue("RANSEED") > 0) {
    procid = Geometry::procID();
    seed   = -abs((procid + 1) * (char) time(NULL));
  } else seed = -1;
  Veclib::ranInit (seed);

  // -- Build all the elements.
  VERBOSE cout << "Building elements ... ";
  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, mesh);
  VERBOSE cout << "done" << endl;

  // -- Build all the boundary condition applicators.
  VERBOSE cout << "Building boundary condition manager ..." << endl;
  bman = new BCmgr (file, elmt);
  VERBOSE cout << "done" << endl;

  // -- Build the solution domain.
  VERBOSE cout << "Building domain ..." << endl;
  domain = new Domain (file, elmt, bman);
  VERBOSE cout << "done" << endl;

  // -- Build field force.
  VERBOSE cout << "Building field force ..." << endl;
  FF = new FieldForce (domain, file);
  VERBOSE cout << "done" << endl;
}
