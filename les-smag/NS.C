/*****************************************************************************
 * NS.C:  Unsteady Navier--Stokes solver, using "stiffly-stable" integration.
 *
 * Reference: Karniadakis, Israeli & Orszag 1991.  "High-order splitting
 * methods for the incompressible Navier--Stokes equations", JCP 9(2).
 *****************************************************************************/

static char RCSid[] = "$Id$";

#include <Fem.h>

#ifdef __DECCXX
  #pragma define_template roll<Field*>
#endif


static void  setUForce (Domain*, Field***);
static void  setPForce (Domain*, Field***, Field***);
static void  project   (Domain*, Field***, Field***);
static void  waveProp  (Domain*, Field***, Field***, Vector, Vector);
static real  gamma0    (int);


void  NavierStokes (Domain* D, Mesh* M, char** forcing)
// ---------------------------------------------------------------------------
// On entry, D contains storage for velocity Field 'u'.
//
// On exit, D is set up so that it contains DIM velocity Fields
// and a pressure Field (in that order).
//
// Us is multi-level auxillary Field storage for velocities and 
// Uf is multi-level auxillary Field storage for nonlinear forcing terms.
// ---------------------------------------------------------------------------
{
  char        routine[] = "NavierStokes";

  const real  dt        = dparam ("DELTAT");
  const int   nOrder    = iparam ("N_TIME");
  const int   nStep     = iparam ("N_STEP");
  const int   DIM       = iparam ("N_VAR" );
  
  Field*      Ftmp;
  Field*      Pressure;

  Field***    Us = new Field** [DIM];
  Field***    Uf = new Field** [DIM];

  Vector      a  = {0.0, 0.0, 0.0};   // Frame acceleration for N--S.
  Vector      f  = {0.0, 0.0, 0.0};   // Spatially-constant forcing for N--S.
  Vector      fv = {0.0, 0.0, 0.0};
  Vector      fp = {0.0, 0.0, 0.0};

  // -- Set up to solve velocity viscous step, build other velocity fields.

  const real lambda2 = gamma0 (nOrder) / (dt * dparam ("KINVIS"));

  message (routine, ": -- Building velocity matrices", REMARK);
  D -> u[0] -> buildSys (lambda2);

  for (int i = 1; i < DIM; i++) {
    Ftmp = new Field (*D -> u[0], *M, 'u' + i);
    D -> addField (Ftmp);
  }

  // -- Set up multi-level storage for velocities and forcing.

  for (i = 0; i < DIM; i++) {
    Us[i] = new Field* [nOrder];
    Uf[i] = new Field* [nOrder];
     for (int j = 0; j < nOrder; j++) {
      Us[i][j] = new Field (*D -> u[i], *M);
      Uf[i][j] = new Field (*D -> u[i], *M);
    }
  }

  // -- Set up pressure field and solution matrices.

  Pressure = new Field (*D -> u[0], *M, 'p');
  PBCmanager::build    (*Pressure);
  Pressure -> connect  (*M, iparam ("N_POLY"));
  D -> addField (Pressure);

  message (routine, ": -- Building pressure matrices", REMARK);
  Pressure  -> buildSys (0.0);
 
  // -- Initialize velocity fields and associated boundary conditions.

  setDparam ("t", D -> time ());

  D -> restart ();

  for (i = 0; i < DIM; i++)
    D -> u[i] -> evaluateBoundaries (0);

  // -- Timestepping loop.

  while (D -> step () < nStep) {
 
    D -> step () += 1; 
    D -> time () += dt;
    setDparam ("t", D -> time ());

    // -- Unconstrained forcing step.
    
    waveProp (D, Us, Uf, a, f);

    // -- Pressure projection step.

    PBCmanager::maintain (D -> step (), Pressure, Us, Uf);
    Pressure -> evaluateBoundaries (D -> step ());
    for (i = 0; i < DIM; i++) {
      roll (Us[i], nOrder);
      roll (Uf[i], nOrder);
    }
    setPForce (D, Us, Uf);
    Pressure -> solveSys (Uf[0][0]);
    project   (D, Us, Uf);

    // -- Viscous correction step.

    setUForce (D, Us);
    for (i = 0; i < DIM; i++)
      D -> u[i] -> evaluateBoundaries (D -> step ());

    D -> u[0] -> solveSys (Us[0][0]);
    D -> u[1] -> solveSys (Us[1][0]);

    // -- Analysis.

    fp = Field::normalTraction  (*D -> u[2], *D -> u[0]);
    fv = Field::tangentTraction (*D -> u[0], *D -> u[1], *Us[0][0], *Us[1][0]);

    char     s[StrMax];
    sprintf (s, ": Step: %d  Time: %f"
	        "  Fvx: %f  Fpx: %f  Fx: %f  Fvy: %f  Fpy: %f  Fy: %f",
	     D -> step (), D -> time (),
	     fv.x, fp.x, fv.x + fp.x, fv.y, fp.y, fv.y + fp.y);
    message ("NS", s, REMARK);

    D -> dump ();
  }
}


static void  waveProp (Domain*   D ,   // Using,

		       Field***  Us,   // update.
		       Field***  Uf,

		       Vector    a ,
		       Vector    f )
// ---------------------------------------------------------------------------
// Compute the first substep of stiffly-stable timestepping scheme.
//
// The next stage of nonlinear forcing terms N(u) - a + f are computed from
// velocity fields held in D and left in the first level of Uf, then the
// intermediate velocity field u^ is computed and left in Domain velocity
// storage.  The Domain pressure field is used as scratch.
//
// Nonlinear terms are computed in the skew-symmetric form (Zang 1991)
//                    -0.5 ( u . grad u + div uu ).
//                           ~        ~       ~~
// ---------------------------------------------------------------------------
{
  // -- Put old velocities in Us, use H to construct u^, which goes in D.

  *Us[0][0] = *D -> u[0];
  *Us[1][0] = *D -> u[1];

  Field& Ux = *Us[0][0];
  Field& Uy = *Us[1][0];
  Field& Tp = *D -> u[D -> nField () - 1];

  Field& Hx = *D -> u[0];
  Field& Hy = *D -> u[1];

  Field& Nx = *Uf[0][0];
  Field& Ny = *Uf[1][0];

#ifdef STOKES  /* -- No nonlinear terms. */
  
  Nx = 0.0;
  Ny = 0.0;

#else          /* -- Build skew-symmetric nonlinear terms. */

  // -- Conservative NL terms.

  Nx.prod (Ux, Uy);
  Ny = Nx;
  Nx.grad (0, 1);
  Ny.grad (1, 0);

  Tp.prod (Ux, Ux);
  Tp.grad (1, 0);
  Nx += Tp;

  Tp.prod (Uy, Uy);
  Tp.grad (0, 1);
  Ny += Tp;

  // -- Nonconservative NL terms.

  Tp = Ux;
  Nx.addprod (Ux, Tp.grad (1, 0));

  Tp = Ux;
  Nx.addprod (Uy, Tp.grad (0, 1));

  Tp = Uy;
  Ny.addprod (Ux, Tp.grad (1, 0));

  Tp = Uy;
  Ny.addprod (Uy, Tp.grad (0, 1));

  // -- Smooth result & scale.

  Nx.smooth ();
  Ny.smooth ();

  Nx *= -0.5;
  Ny *= -0.5;
  
#endif

  // -- Add in distributed forcing terms.

  if (fabs (f.x -= a.x) > EPSDP) Nx += f.x;
  if (fabs (f.y -= a.y) > EPSDP) Ny += f.y;

  // -- Construct u^ in Hx & Hy.

  int  Je = iparam ("N_TIME");
  Je = min (D -> step (), Je);

  real  alpha[TIME_ORDER_MAX], beta [TIME_ORDER_MAX];

  Veclib::copy (Je, Icoef[Je - 1].alpha, 1, alpha, 1);
  Veclib::copy (Je, Icoef[Je - 1].beta,  1, beta,  1);
  Blas::scal   (Je, dparam ("DELTAT"),      beta,  1);

  Hx = 0.0;
  Hy = 0.0;

  for (int q = 0; q < Je; q++) {
    Hx.axpy (alpha[q], *Us[0][q]).axpy (beta[q], *Uf[0][q]);
    Hy.axpy (alpha[q], *Us[1][q]).axpy (beta[q], *Uf[1][q]);
  }
}


static real gamma0 (int order)
// ---------------------------------------------------------------------------
// Return coefficient gamma0 for time-stepping scheme. 
// ---------------------------------------------------------------------------
{
  return Icoef[order - 1].gamma;
}


static void setUForce (Domain*  D, Field***  Us)
// ---------------------------------------------------------------------------
// On entry, intermediate velocity storage u^^ is in lowest levels of Us.
// Multiply by -1.0 / (DELTAT * KINVIS) to create forcing for viscous step.
// ---------------------------------------------------------------------------
{
  const int   DIM   = D -> nField () - 1;
  const real  alpha = -1.0 / (dparam ("DELTAT") * dparam ("KINVIS"));

  for (int i = 0; i < DIM; i++)
    *Us[i][0] *= alpha;
}


static void setPForce (Domain*  D, Field***  Us, Field***  Uf)
// ---------------------------------------------------------------------------
// On input, intermediate velocity storage u^ is in D.  Transfer this to
// the first levels of Us and create div u^ / DELTAT in the first dimension,
// first level storage of Uf as a forcing field for discrete PPE.
// ---------------------------------------------------------------------------
{
  const int DIM = D -> nField () - 1;

  for (int i = 0; i < DIM; i++) {
    *Us[i][0] = *D -> u[i];
    *Uf[i][0] = *D -> u[i];
  }
  
  Uf[0][0] -> grad (1, 0);
  Uf[1][0] -> grad (0, 1);

  *Uf[0][0] += *Uf[1][0];
  *Uf[0][0] /= dparam ("DELTAT");

  Uf[0][0] -> smooth ();
}


static void  project (Domain*  D, Field***  Us, Field***  Uf)
// ---------------------------------------------------------------------------
// On input, new pressure field is stored in D and intermediate velocity
// level u^ is stored in lowest level of Us.  Constrain velocity field:
//                    u^^ = u^ - DELTAT * grad P;
// u^^ is left in lowest level of Us.
// ---------------------------------------------------------------------------
{
  const int DIM = D -> nField () - 1;
  
  for (int i = 0; i < DIM; i++) *Uf[i][0] = *D -> u[DIM];

  Uf[0][0] -> grad (1, 0);
  Uf[1][0] -> grad (0, 1);

  const real  dt = dparam ("DELTAT");
  for (i = 0; i< DIM; i++) {
    Us[i][0] -> axpy (-dt, *Uf[i][0]);
    Us[i][0] -> smooth ();
  }
}

