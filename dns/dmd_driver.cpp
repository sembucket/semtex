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
#include <petsc.h>
#include <petscis.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsnes.h>
#include <slepcsvd.h>

static char prog[] = "dmd";
static void getargs    (int, char**, bool&, char*&);
static void preprocess (const char*, FEML*&, Mesh*&, vector<Element*>&,
			BCmgr*&, Domain*&, FieldForce*&);
static int myRank;

struct Context {
    int              nSlice;
    int              nDofsSlice;
    int              nDofsPlane;
    Mesh*            mesh;
    Domain*          domain;
    vector<Field*>   ui;
    // parallel vector scattering data
    int              localSize;
    IS               isl;
    IS               isg;
    VecScatter       ltog;
};

#define XMIN 0.0
#define XMAX (2.0*M_PI)
#define YMIN 0.0
#define YMAX 0.5
#define NELS_X 30
#define NELS_Y 7
#define NSLICE 16

void data_transpose(real_t* data, int nx, int ny) {
  real_t* temp = new real_t[nx*ny];

  for(int iy = 0; iy < ny; iy++) {
    for(int ix = 0; ix < nx; ix++) {
      temp[ix*ny + iy] = data[iy*nx + ix];
    }
  }
  for(int ii = 0; ii < nx*ny; ii++) {
    data[ii] = temp[ii];
  }

  delete[] temp;
}

void elements_to_logical(real_t* data_els, real_t* data_log) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int pt_x, pt_y;
  int index = -1;

  for(int el_y = 0; el_y < NELS_Y; el_y++) {
    for(int el_x = 0; el_x < NELS_X; el_x++) {
      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        index++;

        // skip over right and top edges for each element, as these are redundant
        if(pt_i%(elOrd+1) == elOrd || pt_i/(elOrd+1) == elOrd) continue;

        pt_x = el_x*elOrd + pt_i%(elOrd+1);
        pt_y = el_y*elOrd + pt_i/(elOrd+1);

        data_log[pt_y*NELS_X*elOrd + pt_x] = data_els[index];
      }
    }
  }
}

void logical_to_elements(real_t* data_log, real_t* data_els) {
  int elOrd = Geometry::nP() - 1;
  int nodes_per_el = (elOrd + 1)*(elOrd + 1);
  int shift_els, pt_r, pt_s, pt_x, pt_y;
  
  for(int el_y = 0; el_y < NELS_Y; el_y++) {
    for(int el_x = 0; el_x < NELS_X; el_x++) {
      shift_els = (el_y*NELS_X + el_x)*nodes_per_el;

      for(int pt_i = 0; pt_i < nodes_per_el; pt_i++) {
        pt_r = pt_i%(elOrd+1);
        pt_s = pt_i/(elOrd+1);

        pt_x = el_x*elOrd + pt_r;
        pt_y = el_y*elOrd + pt_s;
        // asseume periodic in x
        if(pt_x == NELS_X*elOrd) pt_x = 0;
        // don't do axis for now
        if(pt_y == NELS_Y*elOrd) continue;

        data_els[shift_els+pt_i] = data_log[pt_y*NELS_X*elOrd + pt_x];
      }
    }
  }
}

void UnpackX(Context* context, vector<Field*> fields, int slice_i, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, field_i, index;
  int nZ = Geometry::nZProc();
  int nNodesX = NELS_X*elOrd;
  AuxField* field;
  const PetscScalar *xArray;
  real_t* data = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecScatterBegin(context->ltog, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(  context->ltog, x, xl, INSERT_VALUES, SCATTER_FORWARD);
  VecGetArrayRead(xl, &xArray);

  index = 0;
  for(field_i = 0; field_i < context->domain->nField(); field_i++) {
    field = fields[slice_i * context->domain->nField() + field_i];

    for(kk = 0; kk < nZ; kk++) {
      // skip over redundant real dofs
      for(jj = 0; jj < NELS_Y*elOrd; jj++) {
        for(ii = 0; ii < nNodesX; ii++) {
          data[jj*nNodesX + ii] = xArray[index++];
        }
      }
      logical_to_elements(data, field->plane(kk));
    }
  }
  VecRestoreArrayRead(xl, &xArray);
  VecDestroy(&xl);

  delete[] data;
}

void RepackX(Context* context, vector<Field*> fields, int slice_i, Vec x) {
  int elOrd = Geometry::nP() - 1;
  int ii, jj, kk, field_i, index;
  int nZ = Geometry::nZProc();
  int nNodesX = NELS_X*elOrd;
  AuxField* field;
  PetscScalar *xArray;
  real_t* data = new real_t[NELS_X*elOrd*NELS_Y*elOrd];
  Vec xl;

  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  VecGetArray(xl, &xArray);

  index = 0;
  for(field_i = 0; field_i < context->domain->nField(); field_i++) {
    field = fields[slice_i * context->domain->nField() + field_i];

    for(kk = 0; kk < nZ; kk++) {
      elements_to_logical(field->plane(kk), data);

      // skip over redundant real dofs
      for(jj = 0; jj < NELS_Y*elOrd; jj++) {
        for(ii = 0; ii < nNodesX; ii++) {
          xArray[index] = data[jj*nNodesX + ii];
          index++;
        }
      }
    }
  }
  VecRestoreArray(xl, &xArray);

  VecScatterBegin(context->ltog, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecScatterEnd(  context->ltog, xl, x, INSERT_VALUES, SCATTER_REVERSE);
  VecDestroy(&xl);

  delete[] data;
}

void dmd_solve(int nSlice, Domain* domain, vector<Field*> ui) {
  Context* context = new Context;
  int elOrd = Geometry::nP() - 1;
  int localShift;
  int slice_i, nEig, eig_i, ind_i, ind_f, ii;
  int *localInds;
  double sigma, vSigmaInv, lambda_r, lambda_i, freq;
  Vec *x, xl, lVec, rVec, vr, vi;
  Mat X, U, VSI, XT, XTVSI, Atilde;
  SVD svd;
  EPS eps;
  PetscScalar* xArray;

  context->nSlice  = nSlice;
  context->domain  = domain;
  context->ui      = ui;

  // add dofs for theta and tau for each time slice
  context->nDofsPlane = NELS_X*elOrd*NELS_Y*elOrd;
  context->nDofsSlice = context->domain->nField() * Geometry::nZ() * context->nDofsPlane;
  context->localSize = context->domain->nField() * Geometry::nZProc() * context->nDofsPlane;
  localShift = myRank * context->localSize;

  x = new Vec[nSlice];
  for(slice_i = 0; slice_i < nSlice; slice_i++) {
    VecCreateMPI(MPI_COMM_WORLD, context->localSize, context->nDofsSlice, &x[slice_i]);
  }

  // create the local to global scatter object
  VecCreateSeq(MPI_COMM_SELF, context->localSize, &xl);
  ISCreateStride(MPI_COMM_SELF, context->localSize, 0, 1, &context->isl);
  ISCreateStride(MPI_COMM_WORLD, context->localSize, localShift, 1, &context->isg);
  VecScatterCreate(x[0], context->isg, xl, context->isl, &context->ltog);

  localInds = new int[context->localSize];
  for(ind_i = 0; ind_i < context->localSize; ind_i++) {
    localInds[ind_i] = localShift + ind_i;
  }

  // pack the data into vectors
  for(slice_i = 0; slice_i < nSlice; slice_i++) {
    RepackX(context, context->ui, slice_i, x[slice_i]);
  }

  // initialise the matrices
  MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nSlice, context->nDofsSlice, NULL, &X);

  MatCreate(MPI_COMM_WORLD, &U);
  MatSetSizes(U, PETSC_DECIDE, PETSC_DECIDE, nSlice, context->nDofsSlice);
  MatSetType(U, MATMPIAIJ);
  MatMPIAIJSetPreallocation(U, context->nDofsSlice/Geometry::nProc()+1, PETSC_NULL, context->nDofsSlice, PETSC_NULL);
  MatZeroEntries(U);

  MatCreate(MPI_COMM_WORLD, &VSI);
  MatSetSizes(VSI, PETSC_DECIDE, PETSC_DECIDE, nSlice, nSlice);
  MatSetType(VSI, MATMPIAIJ);
  MatMPIAIJSetPreallocation(VSI, nSlice, PETSC_NULL, nSlice, PETSC_NULL);
  MatZeroEntries(VSI);
 
  // populate the operator X with the state data 
  for(slice_i = 0; slice_i < nSlice; slice_i++) {
    VecScatterBegin(context->ltog, x[slice_i], xl, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(  context->ltog, x[slice_i], xl, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(xl, &xArray);
    MatSetValues(X, 1, &slice_i, context->localSize, localInds, xArray, INSERT_VALUES);
    VecRestoreArray(xl, &xArray);
  }
  MatAssemblyBegin(X, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  X, MAT_FINAL_ASSEMBLY);
  MatCreateVecs(X, &rVec, &lVec);

  // solve the svd
  SVDCreate(MPI_COMM_WORLD, &svd);
  SVDSetOperator(svd, X);
  SVDSetFromOptions(svd);
  SVDSolve(svd);
  SVDGetConverged(svd, &nEig);

  if(!myRank) cout << "number of svd eigenvalues: " << nEig << endl;
  for(eig_i = 0; eig_i < nEig; eig_i++) {
    SVDGetSingularTriplet(svd, eig_i, &sigma, lVec, rVec);
    if(!myRank) cout << eig_i << "\tsigma: " << sigma << endl;

    VecScatterBegin(context->ltog, rVec, xl, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(  context->ltog, rVec, xl, INSERT_VALUES, SCATTER_FORWARD); 

    // populate the left eigenvectors matrix
    VecGetArray(xl, &xArray);
    MatSetValues(U, 1, &eig_i, context->localSize, localInds, xArray, INSERT_VALUES);
    VecRestoreArray(xl, &xArray);

    // populate the right eigenvectors / sigma matrix
    VecGetOwnershipRange(lVec, &ind_i, &ind_f);
    VecGetArray(lVec, &xArray);
    for(ii = ind_i; ii < ind_f; ii++) {
      vSigmaInv = xArray[ii-ind_i]/sigma;
      MatSetValues(VSI, 1, &eig_i, 1, &ii, &vSigmaInv, INSERT_VALUES);
    }
    VecRestoreArray(lVec, &xArray);
  }
  MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  U, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(VSI, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  VSI, MAT_FINAL_ASSEMBLY);

  // Assemble the approximation of X.X^{-1}
  MatTranspose(X, MAT_INITIAL_MATRIX, &XT);
  MatMatMult(XT, VSI, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XTVSI);
  MatMatMult(U, XTVSI, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Atilde);

  MatCreateVecs(Atilde, NULL, &vr);
  MatCreateVecs(Atilde, NULL, &vi);

  if(!myRank) cout << "solving the eigenvalues....................\n";
  EPSCreate(MPI_COMM_WORLD, &eps);
  EPSSetOperators(eps, Atilde, NULL);
  EPSSetFromOptions(eps);
  EPSSolve(eps);

  // write to file
  for(slice_i = 0; slice_i < nSlice; slice_i++) {
    EPSGetEigenpair(eps, ii, &lambda_r, &lambda_i, vr, vi);
    freq = log(fabs(lambda_r))/1.0;

    if(!myRank) cout << slice_i << "\tlambda:   " << lambda_r << " + " << lambda_i << endl;
    if(!myRank) cout << slice_i << "\tfreqency: " << freq << endl;

    MatMult(XTVSI, vr, x[slice_i]);
    UnpackX(context, context->ui, slice_i, x[slice_i]);
  }

  // clean up
  delete[] localInds;
  VecDestroy(&xl);
  for(slice_i = 0; slice_i < nSlice; slice_i++) {
    VecDestroy(&x[slice_i]);
  }
  delete[] x;
  VecScatterDestroy(&context->ltog);
  VecDestroy(&lVec);
  VecDestroy(&rVec);
  VecDestroy(&vr);
  VecDestroy(&vi);
  ISDestroy(&context->isl);
  ISDestroy(&context->isg);
  MatDestroy(&X);
  MatDestroy(&U);
  MatDestroy(&VSI);
  MatDestroy(&XT);
  MatDestroy(&XTVSI);
  SVDDestroy(&svd);
  EPSDestroy(&eps);
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
  FieldForce*      FF;
  static char      help[] = "petsc";
  vector<Field*>   ui;  // Solution fields for velocities, pressure at the i time slices
  char             session_i[100];

  SlepcInitialize(&argc, &argv, (char*)0, help);

  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  Femlib::initialize (&argc, &argv);
  getargs (argc, argv, freeze, session);

  preprocess (session, file, mesh, elmt, bman, domain, FF);

  //domain -> restart ();
  //ROOTONLY domain -> report ();
  
  // load in the time slices
  ui.resize(NSLICE * domain->nField());
  delete file;
  delete domain;
  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    //sprintf(session_i, "%s.%u", session, slice_i + 1);
    sprintf(session_i, "%s.%u", session, 4*slice_i);
    FEML* file_i = new FEML(session_i);
    domain = new Domain(file_i, elmt, bman);
    domain->restart();
    for(int field_i = 0; field_i < domain->nField(); field_i++) {
      ui[slice_i*domain->nField()+field_i] = domain->u[field_i];
      {
        char* field;
        real_t* data;
        BoundarySys* bndry;

        strcpy ((field = new char [strlen (bman -> field()) + 1]), bman -> field());
        data = new real_t[static_cast<size_t>(Geometry::nTotProc())];
        bndry = new BoundarySys(bman, elmt, field[0]);
      }
    }
    delete file_i;
    delete domain;
  }

  //Femlib::initialize (&argc, &argv);
  file = new FEML(session_i);
  domain = new Domain(file, elmt, bman);

  // solve the newton-rapheson problem
  dmd_solve(NSLICE, domain, ui);

  // dump the output
  for(int slice_i = 0; slice_i < NSLICE; slice_i++) {
    sprintf(session_i, "%s_dmd_%u", session, slice_i + 1);
    FEML* file_i = new FEML(session_i);
    domain = new Domain(file_i, elmt, bman);
    for(int field_i = 0; field_i < domain->nField(); field_i++) {
      domain->u[field_i] = ui[slice_i*domain->nField()+field_i];
    }
    domain->dump();
    delete file_i;
    delete domain;
  } 

  SlepcFinalize();

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
