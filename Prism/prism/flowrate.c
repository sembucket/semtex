/*
 * FLOWRATE --- Drive the velocity field at a fixed flowrate
 *                                                                           
 * An alternative method for driving the flow is to specifiy the flowrate.   
 * The following function solves for a Green's function velocity field cor-  
 * responding to a unit force in the predominant flow direction (x or z).    
 * Subsequent calls to flowrate() will add a correction to the unforced      
 * velocity fields to provide the proper mass flow in the appropriate dir-   
 * ection.                                                                   
 *                                                                           
 * NOTE: For flowrate in the x-direction, Q(t) is evaluated along the        
 *       periodic edges.                                                     
 *                                                                           
 * RCS Information                                                           
 * -----------------------------                                             
 * $Author$  
 * $Date$ 
 * $Source$
 * $Revision$                                         
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "veclib/veclib.h"
#include "prism/prism.h"

#if DIM == 2
#    define VelocitySetup(Omega,V,Vs,Vf,Vbc) \
  V  [0] = Omega->U  ; V  [1] = Omega->V  ;  \
  Vbc[0] = Omega->Ubc; Vbc[1] = Omega->Vbc;  \
  Vf [0] = Omega->Uf ; Vf [1] = Omega->Vf ;  \
  Vs [0] = Omega->Us ; Vs [1] = Omega->Vs ; 
#
#else
#    define VelocitySetup(Omega,V,Vs,Vf,Vbc) \
  V  [0] = Omega->U  ; V  [1] = Omega->V  ; V  [2] = Omega->W  ;  \
  Vbc[0] = Omega->Ubc; Vbc[1] = Omega->Vbc; Vbc[2] = Omega->Wbc;  \
  Vf [0] = Omega->Uf ; Vf [1] = Omega->Vf ; Vf [2] = Omega->Wf ;  \
  Vs [0] = Omega->Us ; Vs [1] = Omega->Vs ; Vs [2] = Omega->Ws ; 
#
# endif

/*
 * Compute the flowrate for a given velocity field
 *
 * 2D: Compute the flow across a given set of edges (assumed to be both
 *     inflow and outflow) and divide by 2.
 *
 * 3D: Compute the flow across the entire (x,y)-cross section.
 */

static double flow_calc (GreensF *G, Field *U)
{
  register double sum = 0.;
  register int    np  = 0 , k;
#
# if DIM == 2
#
  Bedge   *bc   = G->Gbc;
  Edge    *edg;
  double  *w;

  double tmp [_MAX_NORDER];

  while (bc) {
    k    = bc->elmt->id;
    edg  = bc->edge;
    
    if (np != edg->np) getops(np = edg->np, 0, &w, 0, 0);
    
    dvmul (np, edg->area, 1, *U[k].base + edg->start, edg->skip, tmp, 1);
    sum += ddot(np, w, 1, tmp, 1);
    bc   = bc->next;
  }

  sum *= .5;
#
# else
#
  int  nrns  = U->nr * U->ns;
  int  nel   = Field_count(U);
  
  for(k = 0; k < nel; k++)
    sum += ddot(nrns, *U[k].mass, 1, *U[k].base, 1);
#
# endif
#    
  return sum;
}


/*
 * Initialization for a flowrate-driven simulation
 */

void flowrate_init (Domain *omega)
{
  int      ntot, ntotz;
  int      Jmax  = iparam("TORDER");
  double   delt  = dparam("DT");
  Field    *V  [DIM], **Vs[DIM], **Vf[DIM], *P;
  Bedge    *Vbc[DIM];
  double   *Vic[DIM];
  GreensF  *G;
  register int i, n;

  G     = (GreensF *) calloc(1, sizeof(GreensF));
  P     = omega->P;
  ntot  = P->nr * P->ns * Field_count(omega->U);
  ntotz = P->nz * ntot;

  VelocitySetup (omega, V, Vs, Vf, Vbc);

  /* Create the Green's function velocity fields */

  for (i = 0; i < DIM; i++)
    for (n = 0; n < Jmax; n++)
      (*G).Gv[i][n] = Field_aux (P, P->nr, P->ns, 1, 's');
  for (n = 0; n < Jmax; n++)
    (*G).Gp[n] = Field_aux (P, P->nr, P->ns, 1, 's');

  /* Save the initial conditions */

  for (i = 0; i < DIM; i++) {
    dcopy(ntotz,  *V[i]->base, 1, Vic[i] = dvector(0, ntotz), 1);
    dzero(ntotz,  *V[i]->base, 1);
  }

  /* Compute Stokes Flow with "delta" Forcing */
#
# if DIM == 2
#
  dparam_set("FFX", 1.);
  dparam_set("FFY", 0.);
  dparam_set("FFZ", 0.);
#
# else
#
  dparam_set("FFX", 0.);
  dparam_set("FFY", 0.);
  dparam_set("FFZ", 1.);
#
# endif
# 
  for (n = 1; n <= Jmax; n++) {

    set_order (n);

    Vorticity (omega);
    MakeF     (omega, Stokes, n, 0.);
    for (i = 0; i < DIM; i++)
      Integrate (V[i], Vs[i], Vf[i], n, delt);

    MakeF      (omega, Pressure, n, 0.);
    MultiSolve (omega, &P, Vf, &(omega->Pbc), omega->Pressure, 1, n);

    MakeF      (omega, Viscous, n, 0.);
    MultiSolve (omega, V, Vf, Vbc, omega->Velocity, DIM, n);
    
    dcopy (ntot, *P->base, 1, *(*G).Gp[n-1]->base, 1);
    for (i = 0; i < DIM; i++) {
      dcopy (ntot , *V[i]->base, 1, *(*G).Gv[i][n-1]->base, 1);
      dzero (ntotz, *V[i]->base, 1);
    }
  }
#
# if DIM == 2
#
  G->order  =  1;
  G->basis  =  V[0];
  G->Gbc    =  BC_get(Vbc[0], 'P');

  for (n = 0; n < Jmax; n++) (G->Fg)[n] = flow_calc (G, (*G).Gv[0][n]);
#
# else
#
  G->order  =  1;
  G->basis  =  V[2];

  for (n = 0; n < Jmax; n++) (G->Fg)[n] = flow_calc (G, (*G).Gv[2][n]);
#
# endif
#
  omega->G  = G;

  dparam_set("FFX", 0.);                     /* reset the forcing */
  dparam_set("FFY", 0.);
  dparam_set("FFZ", 0.);

  /* Restore initial conditions */

  for (i = 0; i < DIM; i++) {
    dcopy (ntotz, Vic[i], 1, *V[i]->base, 1);
    free  (Vic[i]);
  }

  return;
}

/*
 * Apply a flowrate-correction to the velocity field
 */

void flowrate(Domain *omega)
{
  GreensF *Gf = omega -> G;
  int      Je = iparam("TORDER");
  Field   *V[DIM], *G[DIM], *Gp;
  double   Fa, Fg, Fr, dP;
  register int i, n, ntot;

# if DIM == 2
#
  V[0] = omega->U;
  V[1] = omega->V;
#
# else
#
  V[0] = omega->U;
  V[1] = omega->V;
  V[2] = omega->W;
#
# endif

  Je   = MIN (Je, Gf->order);
  n    = Je-1;
  ntot = V[0]->nr * V[0]->ns * Field_count(V[0]);

  for (i = 0; i < DIM; i++) G[i] = (*Gf).Gv[i][n];
  Gp   = (*Gf).Gp[n];

  Fr   = dparam("FLOWRATE");              /* Specified flowrate        */
  Fa   = flow_calc (Gf, Gf->basis);       /* Actual flowrate           */
  Fg   = (Gf->Fg)[n];                     /* Green's function flowrate */
  dP   = (Fr - Fa) / Fg;                  /* Required pressure drop    */

  for (i = 0; i < DIM; i++)
    daxpy (ntot, dP, *G[i]->base, 1, *V[i]->base, 1);
  daxpy (ntot, dP, *Gp->base, 1, *omega->P->base, 1);

  dparam_set("PDROP", dP);
  dparam_set("FLOWR", Fa);

  Gf->order++;

  return;
}
