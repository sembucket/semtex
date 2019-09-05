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

#include "rpo_utils.h"
#include "rpo_preconditioner.h"

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

#define X_FOURIER
#define NFIELD 3
#define XMIN 0.0
#define YMIN 0.0
#define YMAX 1.0
#define NSLICE 1
//#define RM_2FOLD_SYM

void int_base_flow_ke(Context* context, Field* ux) {
  int elOrd = Geometry::nP() - 1;
  int nModesX = Femlib::ivalue("NELS_X")*elOrd;
  double rh, dr, vol_k, u_bar, ke = 0.0, vol = 0.0;
  real_t* data_r = new real_t[context->nDofsPlane];
  real_t* data_i = new real_t[context->nDofsPlane];

  SEM_to_Fourier(0, context, ux, data_r, data_i);

  if(!Geometry::procID()) {
    for(int pt_y = 0; pt_y < Femlib::ivalue("NELS_Y")*elOrd; pt_y++) {
      rh = context->rad_coords[pt_y];
      if(!pt_y) rh = 0.5*(context->rad_coords[pt_y+1] + context->rad_coords[pt_y]);
      dr = context->rad_weights[pt_y];
      vol_k = context->xmax * dr * 2.0*M_PI*rh;
      u_bar = data_r[pt_y * nModesX] / nModesX;

      ke += 0.5 * vol_k * u_bar * u_bar;

      vol += vol_k;
    }

    cout << "base flow ke: " << ke << ", total volume " << vol << endl;
  }

  delete[] data_r;
  delete[] data_i;
}

void plane_norms(vector<Field*> ui) {
  int plane_j;
  int elOrd = Geometry::nP() - 1;
  real_t* plane;
  double norm_2;
  real_t* data = new real_t[Femlib::ivalue("NELS_X")*Femlib::ivalue("NELS_Y")*elOrd*elOrd];

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int proc_i = 0; proc_i < Geometry::nProc(); proc_i++) {
      if(Geometry::procID() == proc_i) {
        for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i++) {
          plane_j = Geometry::procID() * Geometry::nZProc() + plane_i;
          norm_2 = 0.0;

          plane = ui[field_i]->plane(plane_i);
          for(int dof_i = 0; dof_i < Femlib::ivalue("NELS_X")*Femlib::ivalue("NELS_Y")*(elOrd+1)*(elOrd+1); dof_i++) {
            norm_2 += plane[dof_i]*plane[dof_i];
          }
          if(plane_i%2 == 0) cout << field_i << ",\t" << plane_j/2 << ",\t";
          cout << norm_2;
          if(plane_i%2 == 0) cout << ",\t";
          if(plane_i%2 == 1) cout << endl;
        }
      }
    }
  }

  delete[] data;
}

void rpo_solve(int nSlice, Mesh* mesh, vector<Element*> elmt, BCmgr* bman, Domain* domain, DNSAnalyser* analyst, vector<AuxField*> ui, vector<AuxField*> fi, vector<AuxField*> u0)
{
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  int nModesX;
  real_t dx, er, es, ex, ey;
  real_t ckt, skt, cTmp, rTmp;
  const real_t *qx, *wx;
  double dy, dy_sum, norm;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j, io;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  PetscScalar* xArray;
  Vec x, xl;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

  context->nField   = NFIELD;
  context->mesh     = mesh;
  context->elmt     = elmt;
  context->domain   = domain;
  context->bman     = bman;
  context->analyst  = analyst;
  context->ui       = ui;
  context->fi       = fi;
  context->u0       = u0;
  context->build_PC = true;
  context->x_fourier = true;
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->nElsX    = Femlib::ivalue("NELS_X");
  context->nElsY    = Femlib::ivalue("NELS_Y");

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  //context->nModesX = nNodesX/2;
  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");

  nModesX = context->nModesX;

  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    context->theta_i[slice_i] = 0.0;
    context->phi_i[slice_i] = 0.0;
    context->tau_i[slice_i] = 0.0;
  }
  context->domain->step++;

  // setup the fourier mapping data
  context->el = new int_t[Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd];
  context->r  = new real_t[Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd];
  context->s  = new real_t[Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd];

  context->rad_weights = new double[context->nElsY*elOrd+1];
  context->rad_coords  = new double[context->nElsY*elOrd+1];

  dx = (Femlib::value("XMAX")/* - XMIN*/)/nNodesX;
  Femlib::quadrature(&qx, &wx, 0, 0, elOrd+1, GLJ, 0.0, 0.0);

  for(int pt_i = 0; pt_i < Femlib::ivalue("NELS_X")*elOrd*Femlib::ivalue("NELS_Y")*elOrd; pt_i++) {
    pt_x = pt_i%(Femlib::ivalue("NELS_X")*elOrd);
    pt_y = pt_i/(Femlib::ivalue("NELS_X")*elOrd);
    el_x = pt_x/elOrd;
    el_y = pt_y/elOrd;
    el_i = el_y*Femlib::ivalue("NELS_X") + el_x;
    ex = XMIN + pt_x*dx;
    // element y size increases with distance from the boundary
    ey = elmt[el_i]->_ymesh[(pt_y%elOrd)*(elOrd+1)];

    found = false;  
    for(el_j = 0; el_j < mesh->nEl(); el_j++) {
      if(elmt[el_j]->locate(ex, ey, er, es, &work[0], true) && !found) {
        if(fabs(er) < 1.0000000001) {
          context->el[pt_i] = el_j;

          if(er > +0.99999999) er = +0.99999999;
          if(er < -0.99999999) er = -0.99999999;
          if(es > +0.99999999) es = +0.99999999;
          if(es < -0.99999999) es = -0.99999999;

          context->r[pt_i] = er;
          context->s[pt_i] = es;
          found = true;
          //break;
        }
      }
      if(found) break;
    }
    if(!found || el_j==mesh->nEl()) {
      cout << "ERROR! element does not contain point: " << pt_i << "\tx: " << ex << "\ty: " << ey << endl;
      cout << "       pt x: " << pt_x << "\tpt y: " << pt_y << endl;
      abort();
    }
  }

  // compute the radial weights
  for(int pt_i = 0; pt_i <= context->nElsY*elOrd; pt_i++) context->rad_weights[pt_i] = 0.0;

  dy_sum = 0.0;
  for(int el_y = 0; el_y < context->nElsY; el_y++) {
    el_i = el_y*context->nElsX;
    dy = fabs(elmt[el_i]->_ymesh[elOrd*(elOrd+1)] - elmt[el_i]->_ymesh[0]);
    if(!Geometry::procID())printf("%d\tdy: %g\n",el_y,dy);
    for(int qp_i = 0; qp_i <= elOrd; qp_i++) {
      pt_y = el_y*elOrd + qp_i;
      context->rad_weights[pt_y] += 0.5 * dy * wx[qp_i];
      context->rad_coords[pt_y]   = dy_sum + 0.5 * dy * (qx[qp_i]+1.0);
    }
    dy_sum += dy;
  }
  context->u_scale[0] = Femlib::value("U_SCALE");
  context->u_scale[1] = Femlib::value("V_SCALE");
  context->u_scale[2] = Femlib::value("W_SCALE");
  if(!Geometry::procID()) printf("u scales: %g, %g, %g\n", context->u_scale[0], context->u_scale[1], context->u_scale[2]);

  context->rdr = Int_rdr(context);
  for(int ii = 0; ii < context->nElsY*elOrd; ii++) context->rdr[ii] = 0.0;
  for(int el_y = 0; el_y < context->nElsY; el_y++) {
    double ir[99];
    double is[99];
    double dr[99];
    double ds[99];
    double drg;
    double rad;

    el_i = el_y*context->nElsX;
    dy = fabs(elmt[el_i]->_ymesh[elOrd*(elOrd+1)] - elmt[el_i]->_ymesh[0]);

    for(int qp_i = 0; qp_i <= elOrd; qp_i++) {
      Femlib::interpolation (ir,is,dr,ds,elOrd+1,GLJ,0.0,0.0,elOrd+1,GLJ,0.0,0.0,qx[0],qx[qp_i]);
      pt_y = el_y*elOrd + qp_i;
      if(pt_y == context->nElsY*elOrd) continue;

      //drg = 0.0;
      //for(int qp_j = 0; qp_j <= elOrd; qp_j++) drg += context->rad_coords[el_y*elOrd+qp_j]*ds[qp_j];
      drg = 1.0/(elmt[el_i]->_dsdy[qp_i*(elOrd+1)]);

      //context->rdr[pt_y] += 0.5*dy*wx[qp_i]*context->rad_coords[pt_y]*drg;
      rad = 0.5*(context->rad_coords[pt_y] + context->rad_coords[pt_y+1]);
      context->rdr[pt_y] += 0.5 * dy * wx[qp_i] * rad * drg;
    }
  }
  if(!Geometry::procID()) {for(int ii = 0; ii < context->nElsY*elOrd; ii++)cout << "\t" << context->rdr[ii]; cout << endl;}

  // add dofs for theta and tau for each time slice
  context->nSlice     = NSLICE;
  context->nField     = NFIELD;
  context->nDofsPlane = nModesX*Femlib::ivalue("NELS_Y")*elOrd;
  context->nDofsSlice = NFIELD * Geometry::nZ() * context->nDofsPlane + 2;
  context->localSize  = nSlice * NFIELD * Geometry::nZProc() * context->nDofsPlane;
  if(!Femlib::ivalue("TRAV_WAVE")) context->nDofsSlice++;
  if(!Geometry::procID()) {
    context->localSize += (nSlice*2);
    if(!Femlib::ivalue("TRAV_WAVE")) context->localSize += nSlice;
  }
#ifdef RM_2FOLD_SYM
  if(Geometry::procID()%2==1) context->localSize = 0;
  MPI_Allreduce(&context->localSize, &context->nDofsSlice, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

if(!Geometry::procID())cout<<"assigning scatter...\n";
  assign_scatter_semtex(context);
if(!Geometry::procID())cout<<"....done.           \n";

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice * context->nDofsSlice, &x);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);
  VecNorm(x, NORM_2, &norm);
  if(!Geometry::procID()) cout << "|x_0|: " << norm << endl;
/*
  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  for(int field_i = 0; field_i < NFIELD; field_i++) *context->domain->u[field_i] = *context->ui[field_i];
  if(!Geometry::procID()) cout << "applying phase shifts...\n";
  phase_shift_x(context, -0.3*M_PI, +1.0, context->domain->u);
  for(int field_i = 0; field_i < NFIELD; field_i++) *context->ui[field_i] = *context->domain->u[field_i];

  RepackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);
*/
  UnpackX(context, context->ui, context->theta_i, context->phi_i, context->tau_i, x);

  VecDestroy(&x);
  VecDestroy(&xl);
  VecScatterDestroy(&context->global_to_semtex);
  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
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
  vector<AuxField*>   fi;  // Solution fields for velocities, pressure at the i time slices
  vector<AuxField*>   u0;  // Solution fields for velocities, pressure at the i time slices
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
  ui.resize(NSLICE * NFIELD);
  fi.resize(NSLICE * NFIELD);
  u0.resize(NSLICE * NFIELD);
  delete file;
  delete domain;
  sprintf(session_i, "%s", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  domain->restart();

  for(int field_i = 0; field_i < NFIELD; field_i++) {
    ui[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'U'+field_i);
    fi[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'F'+field_i);
    u0[field_i] = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), elmt, 'A'+field_i);

    *ui[field_i] = *domain->u[field_i];
  }

  // solve the newton-rapheson problem
  rpo_solve(NSLICE, mesh, elmt, bman, domain, analyst, ui, fi, u0);
  delete file_i;
  delete domain;

  // dump the output
  sprintf(session_i, "%s_shift", session);
  file_i = new FEML(session_i);
  domain = new Domain(file_i, elmt, bman);
  for(int field_i = 0; field_i < NFIELD; field_i++) {
    *domain->u[field_i] = *ui[field_i];
  }
  domain->dump();
  delete file_i;
  delete domain;

  if(!Geometry::procID()) cout << "...done.\n";

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
