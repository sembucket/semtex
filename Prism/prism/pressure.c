/*
 * Pressure boundary conditions
 *
 * $Id$
 * -------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "prism/prism.h"
#include "veclib/veclib.h"

/*
 * Build boundary conditions for the pressure
 */

Bedge *BuildPBCs (Field *P, Bedge *template)
{
  int    nbcs = 0;
  Bedge *Pbc  = NULL;

  const int nb  = (FIELD_NR(P) + FIELD_NS(P) - 2) * 2;
  const int nz  = template->elmt->nz;
  const int Je  = iparam("TORDER");
  const int nel = Field_count(P);
  int i, k;

  /* Count the number of boundaries to create */

  while (template[nbcs++].next)
    /* empty while loop */;

  if (nbcs) {

    Pbc = (Bedge*) malloc (nbcs*sizeof(Bedge));

    for (i = 0; i < nbcs; i++) {
      Pbc[i].id          = i;
      Pbc[i].type        = template[i].type;
      Pbc[i].elmt        = Element_get  (P, template[i].elmt->id);
      Pbc[i].edge        = Element_edge (Pbc[i].elmt, template[i].edge->id);
      Pbc[i].next        = Pbc + i + 1;
    }
    Pbc[nbcs-1].next = (Bedge*) NULL;

    /* Translate the boundary conditions */

    for (i = 0; i < nbcs; i++) {
      switch (Pbc[i].type) {

      case 'D': 
      case 'N': 
	Pbc[i].bc.value = 
	  (double*) calloc (Pbc[i].edge->np*nz,sizeof(double));
	break;     	  
	
      case 'O': {
	int *solve;
	Pbc[i].type = 'o';     
	Pbc[i].bc.value = 
	  (double*) calloc (Pbc[i].edge->np*nz*Je,sizeof(double));
	break;
      }
      
      case 'V': case 'v':          /* allocate additional multi-step storage */
      case 'W': 
	Pbc[i].type = 'F';
	Pbc[i].bc.value = 
	  (double*) calloc (Pbc[i].edge->np*nz*Je,sizeof(double));
	break;
	
      case 'M': case 'S':                /* make a copy of the boundary info */
	memcpy 
	  (Pbc[i].bc.value = 
	   (double*) calloc (Pbc[i].edge->np,sizeof(double)),
	   template[i].bc.value, Pbc[i].edge->np * sizeof(double));
	break;
	
      default:
	Pbc[i].bc.value = 
	  (double*) calloc (Pbc[i].edge->np, sizeof(double));
	break;
      }
    }
  }

  return Pbc;
}

/* ------------------------------------------------------------------------- *
 * ComputePBCs() -- Compute the high-order Pressure BC's                     *
 *                                                                           *
 * This function computes the (x,y)-components of the curl of the vorticity  *
 * field for use in the pressure boundary conditions.  It loops over the     *
 * velocity boundary condition array, and each one that represents a solid   *
 * wall boundary triggers the curl computation for that element.  A flag ar- *
 * ray is used to prevent elements from being repeated.                      *
 *                                                                           *
 * All of the input (<), output (>) and workspace (*) for this function      *
 * comes from the domain data structure.  The following members are used:    *
 *                                                                           *
 * struct domain {                                                           *
 *                                                                           *
 *   Pbc     >   pressure boundary conditions                                *
 *   Q_i     <   vorticity                                                   *
 *   Us_i    *   curl of the vorticity                                       *
 *   Uf_i    *   workspace (3D-only)                                         *
 *                                                                           *
 * }                                                                         *
 *                                                                           *
 * Results from this function are stored in the pressure boundary condition  *
 * array.                                                                    *
 *                                                                           *
 * NOTE: This function overwrites level 0 of the velocity multi-step array   *
 * ------------------------------------------------------------------------- */

void ComputePBCs (Domain *omega)
{
  Bedge   *Pbc    = omega->Pbc;
  BSystem *M      = omega->Pressure;
  const int nel   = M->elements;
  const int nrns  = omega->U->nr * omega->U->ns;
  const int nz    = omega->U->nz;
  const int ntot  = nrns * nel;

  const double viscos = scalar("1./Re");

  Vector Q, Work, BC;
  int flag[_MAX_NEL];
  int i, k;

  Q.x = omega->Qx; BC.x = omega->Us[0]; Work.x = omega->Uf[0];  
  Q.y = omega->Qy; BC.y = omega->Vs[0]; Work.y = omega->Vf[0];
  Q.z = omega->Qz;

  /* .......... High-Order Pressure Boundary Conditions .......... */

  izero (nel, flag, 1);
  while (Pbc) {
    if (Pbc->type == 'F') {
      const int id    = Pbc->elmt->id;
      const int np    = Pbc->edge->np;
      const int start = Pbc->edge->start;
      const int skip  = Pbc->edge->skip;

      double *nx     = Pbc->edge->unx;
      double *ny     = Pbc->edge->uny;
      double *dpdn   = Pbc->bc  .value;
      double *curl_x = *BC.x[id].base + start;
      double *curl_y = *BC.y[id].base + start;
      
      /* Check to see if we need to compute the B.C. for this element */

      if (!flag[id]) {
#if DIM == 2
	Element_grad 
	         (Q.z+id, BC.y+id, BC.x+id);
	dneg     (nrns,  *BC.y[id].base, 1);
#else
	Element_grad_3D
	         (Q.z+id, Work.x+id, Work.y+id, ntot);
	Element_gradz
	         (Q.x+id, Q.z+id, ntot);
	Element_gradz
                 (Q.y+id, Q.y+id, ntot);

	for (k = 0; k < nz; k++) {
	  dvsub (nrns, *Work.y[id].base + k*ntot, 1, 
		          *Q.y[id].base + k*ntot, 1,
		         *BC.x[id].base + k*ntot, 1);
	  dvsub (nrns,    *Q.z[id].base + k*ntot, 1,
		       *Work.x[id].base + k*ntot, 1,
		         *BC.y[id].base + k*ntot, 1);
	}
#endif
	flag[id] = 1;
      }

      /* Now store the boundary condition */

#if DIM == 3
      for (k = 0; k < nz; k++, dpdn += np, curl_x += ntot, curl_y += ntot) 
#endif
	for (i = 0; i < np; i++)
	  dpdn[i] = -viscos * 
	    (nx[i] * curl_x[i*skip] + ny[i] * curl_y[i*skip]);
    }
    Pbc = Pbc->next;
  }
}

/* ------------------------------------------------------------------------- *
 * SetPBCs() -  Set boundary conditions for the pressure                     *
 *                                                                           *
 *                 Je                     1                                  *
 *      dP/dn = SUM    beta[q] * [ N(u) - - curl ( curl u ) ] * n            *
 *                 q=0                    R                                  *
 *                                                                           *
 * where n is the unit outward normal along the edge, u is the velocity      *
 * field, and N(u) are the nonlinear terms in the momentum equation.  This   *
 * routine computes the RHS of the above equation at the <new> time level    *
 * and saves it in the corresponding Bedge for Pressure.                     *
 * ------------------------------------------------------------------------- */

static void outflow        (Bedge*, Domain*);
static void high_order_pbc (Bedge*, Domain*);
static void rotate         (Bedge*, double*, const, const int);

void SetPBCs (Domain *omega)
{
  Bedge *Pbc = omega->Pbc;

  const int Je = iparam("TORDER");
  const int EQ = iparam("EQTYPE");

  double beta[3];
  get_beta (beta);          /* Get the integration coefficients */

  while (Pbc) {
    switch (Pbc->type) {
    case 'D': case 'N': case 'P': case 'M': case 'S':
      break;

    case 'o':
      if (EQ == Rotational) {
	outflow (Pbc, omega);
	rotate  (Pbc, beta, Je, Pbc->edge->np * Pbc->elmt->nz);
      }
      break;


    case 'F':
      high_order_pbc 
	      (Pbc, omega);
      rotate  (Pbc, beta, Je, Pbc->edge->np * Pbc->elmt->nz);
      break;

    default:
      Prism_error("Prism: uknown bc in SetPBCs\n");
      break;
    }
    
    Pbc = Pbc->next;
  }
}

/* -------------------------- Private Functions ---------------------------- */

/* Rotate time dependent boundary conditions */

static void rotate (Bedge *Pbc, double *beta, const int Je, const int nplevel)
{
  double tmp [_MAX_NORDER * _MAX_NZ * _MAX_TORDER];
  register int q;


  /* Integrate the boundary conditions */

  dzero (nplevel, tmp, 1);
  for (q = 0; q < Je; q++)
    daxpy (nplevel, beta[q], Pbc->bc.value + q*nplevel, 1, tmp, 1);
  

  /* Rotate the time levels */

  for (q--; 0 < q; q--)
    dcopy (nplevel, Pbc->bc.value + (q-1)*nplevel, 1,
	            Pbc->bc.value +  q   *nplevel, 1);


  /* Install the new boundary condition */

  dcopy (nplevel, tmp, 1, Pbc->bc.value, 1);

  return;
}



/*
 * PI = p + 1/2 U.U   (outflow boundary conditions)
 */

static void outflow (Bedge *bc, Domain *omega)
{
  const int    np    = bc->edge->np,
           start = bc->edge->start,
           skip  = bc->edge->skip,
           id    = bc->elmt->id;
  double   half  = .5, **u, **v, **w;
  register int i, j, m;

# if DIM == 2
#
  u = omega->U[id].base;
  v = omega->V[id].base;
  for (i = 0; i < np; i++)
    (bc->bc.value)[i] = half * 
      ((*u+start)[i*skip] * (*u+start)[i*skip] + 
       (*v+start)[i*skip] * (*v+start)[i*skip] );
#      
# else
#
  const int pid    = option("procid");
  const int nprocs = option("nprocs");
  const int nzl    = bc->elmt->nz;
  const int nz     = nzl*nprocs;
  const int ntot   = (omega->U->nr * omega->U->ns) * Field_count(omega->U);

  u = dmatrix (0, np-1, 0, nz-1);   dzero (np*nz, *u, 1);
  v = dmatrix (0, np-1, 0, nz-1);   dzero (np*nz, *v, 1);
  w = dmatrix (0, np-1, 0, nz-1);   dzero (np*nz, *w, 1);
  
  for (m = 0; m < nzl; m = m || pid ? m+1 : m+2) {
    dcopy (np, *(omega->U[id]).base + start + m*ntot, skip, 
	       *u + m + pid*nzl, nz);
    dcopy (np, *(omega->V[id]).base + start + m*ntot, skip, 
	       *v + m + pid*nzl, nz);
    dcopy (np, *(omega->W[id]).base + start + m*ntot, skip, 
	       *w + m + pid*nzl, nz);
  }

#ifdef PARALLEL
  comm_dsum (np*nz, *u, *omega->P->base); /* merge these arrays across the */
  comm_dsum (np*nz, *v, *omega->P->base); /* processors.                   */
  comm_dsum (np*nz, *w, *omega->P->base);
# endif  

  for (i = 0; i < np; i++) {
    realft (nz/2, u[i], 1);
    realft (nz/2, v[i], 1);
    realft (nz/2, w[i], 1);
  }

  for (i = 0; i < np; i++) {          /* Replace u with u.u/2 */
    for (m = 0; m < nz; m++)
      u[i][m] = half * (u[i][m]*u[i][m] + v[i][m]*v[i][m] + w[i][m]*w[i][m]);
    realft (nz/2, u[i], -1);
    dcopy  (nzl,  u[i] + pid*nzl, 1, bc->bc.value + i, np);
  }

  free_dmatrix (u, 0, 0);
  free_dmatrix (v, 0, 0);
  free_dmatrix (w, 0, 0);
#
#endif
}

/*
 * dP/dn = n * [ N(u) - 1/Re curl (curl U) ]
 */

static void high_order_pbc (Bedge *bc, Domain *omega)
{
  const int nz       = bc -> elmt -> nz;
  const int id       = bc -> elmt -> id;
  const int np       = bc -> edge -> np;
  const int skip     = bc -> edge -> skip;
  const int start    = bc -> edge -> start;

  double *nx       = bc -> edge -> unx,
         *ny       = bc -> edge -> uny,
         *dpdn     = bc -> bc   .  value;
          
  double *Nu       = *(*omega->Uf)[id].base + start,
         *Nv       = *(*omega->Vf)[id].base + start;

  const int   ntot     = omega->U->nr * omega->U->ns * Field_count(omega->U);
  
  register int i;


  /* NOTE: In the following loop, n.N(u) is added to the current pressure  *
   * boundary condition.  ComputePBCs() MUST be called first to initialize *
   * the bc.value array with the curl of the vorticity.                    */

#
# if DIM == 3
#
  register int m;
  for (m = 0; m < nz; m++, dpdn += np, Nu += ntot, Nv += ntot)
# endif
#
    for (i = 0; i < np; i++)
      dpdn[i] += nx[i] * Nu[i*skip] + ny[i] * Nv[i*skip];
  
  return;
}

