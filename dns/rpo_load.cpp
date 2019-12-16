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

#define NFIELD 3
#define XMIN 0.0
#define YMIN 0.0
#define YMAX 1.0
#define NSLICE 1

void rpo_solve(int nSlice, Mesh* mesh, vector<Element*> elmt, BCmgr* bman, Domain* domain, DNSAnalyser* analyst) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int nNodesX = Femlib::ivalue("NELS_X")*elOrd;
  int nModesX = nNodesX;
  real_t dx, er, es, ex, ey;
  real_t ckt, skt, cTmp, rTmp;
  const real_t *qx, *wx;
  double dy, dy_sum, norm;
  int_t pt_x, pt_y, el_x, el_y, el_i, el_j, ii, mode_l;
  bool found;
  vector<real_t> work(static_cast<size_t>(max (2*Geometry::nTotElmt(), 5*Geometry::nP()+6)));
  real_t *data_r, *data_i;
  ifstream file;
  double value;
  string line;
  Vector du;
  char filename[100];
  PetscScalar* xArray;
  Vec x, xl;
  PetscViewer viewer;
  int np2 = Geometry::nP() * Geometry::nP();
  AuxField* uBar;

  if(!Geometry::procID()) cout << "time step: " << Femlib::value("D_T") << endl;

  context->nField   = NFIELD;
  context->mesh     = mesh;
  context->elmt     = elmt;
  context->domain   = domain;
  context->bman     = bman;
  context->analyst  = analyst;
  context->travelling_wave = Femlib::ivalue("TRAV_WAVE");
  context->nElsX    = Femlib::ivalue("NELS_X");
  context->nElsY    = Femlib::ivalue("NELS_Y");

  context->theta_i = new real_t[NSLICE];
  context->phi_i   = new real_t[NSLICE];
  context->tau_i   = new real_t[NSLICE];

  context->nModesX = nNodesX;
  context->xmax    = Femlib::value("XMAX");

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

  data_r = new real_t[Femlib::ivalue("NELS_Y")*elOrd*context->nModesX];
  data_i = new real_t[Femlib::ivalue("NELS_Y")*elOrd*context->nModesX];

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
if(!Geometry::procID()){
cout << "radial coords:\n";
for(int pt_i=0; pt_i<context->nElsY*elOrd; pt_i++)cout<<context->rad_coords[pt_i]<<", ";
cout << "\n";
}

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

  assign_scatter_semtex(context);

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecCreateMPI(MPI_COMM_WORLD, context->localSize, nSlice * context->nDofsSlice, &x);

  // setup the complex fft in the axial direction
  context->data_s = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->data_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(context->nModesX));
  context->trans_fwd = fftw_plan_dft_1d(context->nModesX, context->data_s, context->data_f, FFTW_FORWARD,  FFTW_ESTIMATE);
  context->trans_bck = fftw_plan_dft_1d(context->nModesX, context->data_f, context->data_s, FFTW_BACKWARD, FFTW_ESTIMATE);

  if(!Geometry::procID()) {
    cout << "nModesX:    " << context->nModesX << endl;
    cout << "nDofsPlane: " << context->nDofsPlane << endl;
  }

  // load from files (one for each proc)
  VecZeroEntries(xl);
  VecGetArray(xl, &xArray);
  sprintf(filename, "u.%.3d.rpo", Geometry::procID());
  cout << "loading file: " << filename << endl;
  file.open(filename);
  ii = 0;
  while (std::getline(file, line)) {
    stringstream ss(line);
    ss >> value;
    xArray[ii++] = value;
  }
  file.close();

  for(int field_i = 0; field_i < 3; field_i++) {
    for(int plane_i = 0; plane_i < Geometry::nZProc(); plane_i += 2) {
      for(int point_y = 0; point_y < context->nElsY*elOrd; point_y++) {
        for(int point_x = 0; point_x < context->nModesX; point_x++) {
          mode_l = (point_x <= context->nModesX/2) ? point_x : point_x - context->nModesX; // fftw ordering of complex data

          ii = LocalIndex(context, field_i, plane_i+0, point_x, point_y);
          data_r[point_y*context->nModesX+point_x] = xArray[ii];
          // divergence free: mean component of radial velocity is 0
          // note that we are in \tilde{} variables, and the nyquist frequency is also 0
          if(field_i  > 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) data_r[point_y*context->nModesX+point_x] = 0.0;

if(field_i==0 && !Geometry::procID() && !plane_i && !mode_l) cout << data_r[point_y*context->nModesX+point_x] << ", ";

          ii = LocalIndex(context, field_i, plane_i+1, point_x, point_y);
          data_i[point_y*context->nModesX+point_x] = xArray[ii];
          // don't include the nyquist frequency
          if(field_i == 0 && Geometry::procID() == 0 && plane_i == 0 && mode_l == 0) data_i[point_y*context->nModesX+point_x] = 0.0;
        }
      }

      Fourier_to_SEM(plane_i, context, domain->u[field_i], data_r, data_i, field_i);
    }
  }
  VecRestoreArray(xl, &xArray);

  if(!Geometry::procID())cout<<"\n";
  VecNorm(xl,NORM_2,&norm);
  cout << Geometry::procID() << "\t" << norm << endl;

  VecScatterBegin(context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->global_to_semtex, xl, x, INSERT_VALUES, SCATTER_REVERSE);

  sprintf(filename, "x_curr_%.4u.vec", 1);
  PetscViewerBinaryOpen(MPI_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
  VecView(x, viewer);
  PetscViewerDestroy(&viewer);

  for(int field_i = 0; field_i < 3; field_i++) {
    *domain->u[3] = *domain->u[field_i];
    domain->u[3]->gradient(1);
    du = Field::normTraction(domain->u[3]);
    if(!Geometry::procID() && !field_i) cout << "dudy:      " << du.y/Femlib::value("XMAX") << endl;
    domain->u[3]->transform(INVERSE);
    *domain->u[3] *= *domain->u[3];
    domain->u[3]->transform(FORWARD);
    norm = domain->u[3]->integral();

    if(!Geometry::procID()) cout << field_i << " integral: " << norm << endl;
  }

  // add in the base flow
  uBar = new AuxField(new real_t[(size_t)Geometry::nTotProc()], Geometry::nZProc(), domain->elmt, 'b');
  for(int pl_i = 0; pl_i < Geometry::nZProc(); pl_i++) {
    for(int el_i = 0; el_i < domain->elmt.size(); el_i++) {
      for(int pt_i = 0; pt_i < np2; pt_i++) {
        uBar->plane(pl_i)[el_i*np2+pt_i] = (1.0 - domain->elmt[el_i]->_ymesh[pt_i]*domain->elmt[el_i]->_ymesh[pt_i]);
      }
    }
  }
  uBar->transform(FORWARD);
  *domain->u[0] += *uBar;

  if(!Geometry::procID()) cout << "dumping fields...\n";
  domain->dump();

  VecDestroy(&x);
  VecDestroy(&xl);
  VecScatterDestroy(&context->global_to_semtex);
  delete[] context->el;
  delete[] context->r;
  delete[] context->s;
  delete[] data_r;
  delete[] data_i;
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

  PetscInitialize(&argc, &argv, (char*)0, help);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);
  preprocess (session, file, mesh, elmt, bman, domain, FF);

  analyst = new DNSAnalyser (domain, bman, file);
  //domain -> restart ();
  //ROOTONLY domain -> report ();
  
  // solve the newton-rapheson problem
  rpo_solve(NSLICE, mesh, elmt, bman, domain, analyst);
  //delete domain;

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
