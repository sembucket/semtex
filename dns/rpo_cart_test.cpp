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

#define NFIELD 3

static PetscErrorCode RPOVecNormL2_Hookstep(void* ctx,Vec v,PetscScalar* norm) {
  Context* context = (Context*)ctx;
  PetscInt nDofs_l = Geometry::nZProc() * context->n_mesh_sum;
  PetscInt ind_i;
  double norm_sq, norm_l_sq;
  PetscScalar* vArray;
  Vec vl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl);

  VecScatterBegin(context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v, vl, INSERT_VALUES, SCATTER_FORWARD);

  norm_l_sq = 0.0;
  VecGetArray(vl, &vArray);
  for(ind_i=0; ind_i<nDofs_l; ind_i++) {
    norm_l_sq += vArray[ind_i]*vArray[ind_i];
  }
  VecRestoreArray(vl, &vArray);

  norm_sq = 0.0;
  MPI_Allreduce(&norm_l_sq, &norm_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *norm = sqrt(norm_sq);

  VecDestroy(&vl);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDot_Hookstep(void* ctx,Vec v1,Vec v2,PetscScalar* dot) {
  Context* context = (Context*)ctx;
  PetscInt nDofs_l = Geometry::nZProc() * context->n_mesh_sum;
  PetscInt ind_i;
  double dot_l;
  PetscScalar *v1Array, *v2Array;
  Vec vl1, vl2;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl1);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &vl2);

  VecScatterBegin(context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v1, vl1, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, v2, vl2, INSERT_VALUES, SCATTER_FORWARD);

  dot_l = 0.0;
  VecGetArray(vl1, &v1Array);
  VecGetArray(vl2, &v2Array);
  for(ind_i=0; ind_i<nDofs_l; ind_i++) {
    dot_l += v1Array[ind_i]*v2Array[ind_i];
  }
  VecRestoreArray(vl1, &v1Array);
  VecRestoreArray(vl2, &v2Array);

  *dot = 0.0;
  MPI_Allreduce(&dot_l, dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  VecDestroy(&vl1);
  VecDestroy(&vl2);

  PetscFunctionReturn(0);
}

static PetscErrorCode RPOVecDiff_Hookstep(void* ctx,Vec y,Vec F,PetscScalar h) {
  Context* context = (Context*)ctx;
  PetscInt nDofs_l = Geometry::nZProc() * context->n_mesh_sum;
  PetscInt ind_i;
  PetscScalar *yArray, *FArray;
  Vec yl, Fl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &yl);
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &Fl);

  VecScatterBegin(context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, y, yl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterBegin(context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd  (context->global_to_semtex, F, Fl, INSERT_VALUES, SCATTER_FORWARD);

  VecGetArray(yl, &yArray);
  VecGetArray(Fl, &FArray);
  for(ind_i=0; ind_i<nDofs_l; ind_i++) {
    yArray[ind_i] = (yArray[ind_i] - FArray[ind_i])/h;
  }
  VecRestoreArray(yl, &yArray);
  VecRestoreArray(Fl, &FArray);

  VecScatterBegin(context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd  (context->global_to_semtex, yl, y, INSERT_VALUES, SCATTER_REVERSE);

  VecDestroy(&yl);
  VecDestroy(&Fl);

  PetscFunctionReturn(0);
}

void rpo_solve(Mesh* mesh, vector<Element*> elmt, BCmgr* bman, FEML* file, Domain* domain, DNSAnalyser* analyst, FieldForce* FF, 
               vector<AuxField*> ui, vector<AuxField*> fi, vector<AuxField*> u0, char* session)
{
  Context* context = new Context;
  Vec x, f;
  double norm;

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
  context->phi_i    = 0.0;
  context->tau_i    = Femlib::ivalue("N_STEP") * Femlib::value("D_T");
  context->iteration = Femlib::ivalue("RPO_LOAD_VEC");
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->build_dx = false;

  // add dofs for phi and tau for each time slice
  build_addToVector(context, context->domain->u);
  context->nDofsSlice = Geometry::nZ() * context->n_mesh_sum + 1;
  if(!Femlib::ivalue("TRAV_WAVE")) context->nDofsSlice++;
  context->localSize  = Geometry::nZProc() * context->n_mesh_sum;
  if(!Geometry::procID()) {
    context->localSize++;
    if(!Femlib::ivalue("TRAV_WAVE")) context->localSize++;
  }
  
  // assign the coordinate weights
  build_coordWeights(context);

  assign_scatter_semtex(context);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &x);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &f);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &context->x_delta);
  VecZeroEntries(context->x_delta);

  context->uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'b');
  base_profile(context, context->domain->u[2], Femlib::value("BASE_PROFILE_SCALE"), context->uBar);
  *context->ui[2] -= *context->uBar;
  velocity_scales(context);
  *context->ui[2] += *context->uBar;

  _RepackX(context, context->ui, context->phi_i, context->tau_i, x);

  VecNorm(x, NORM_2, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
  if(norm < 1.0e-6) {
    if(!Geometry::procID()) cout << "ERROR: initial state vector norm is SMALL! "
                                 << "Are you sure you loaded the initial condition correctly??\n";
    abort();
  }

  context->c_scale = context->tau_i / norm;
  if(!Geometry::procID()) printf("c scale: %g\n", context->c_scale);

  _UnpackX(context, context->ui, &context->phi_i, &context->tau_i, x);

  if(!Geometry::procID()) cout << "rpo solve complete.\n";

  VecDestroy(&x);
  VecDestroy(&f);
  VecDestroy(&context->x_delta);
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
  vector<AuxField*>   ui;  // Solution fields for velocities, pressure at the i time slices
  vector<AuxField*>   fi;  // Solution fields for flow maps at the i time slices
  vector<AuxField*>   u0;  // Additional fields for building the constraints
  char             session_i[100];
  char*            fname;
  BoundarySys*     bndry;

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  domain -> restart ();
  //ROOTONLY domain -> report ();
  
  // load in the time slices
  ui.resize(3);
  fi.resize(3);
  u0.resize(3);
  for(int field_i = 0; field_i < NFIELD; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }

  // solve the newton-rapheson problem
  rpo_solve(mesh, elmt, bman, file, domain, analyst, FF, ui, fi, u0, session);

  domain->dump();

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
