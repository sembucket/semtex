/* Anaylze() -- flow analysis routines
 * 
 * Copyright (c) 1997 Ronald D. Henderson and Caltech
 *
 * $Id$
 * ------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "veclib/veclib.h"
#include "prism/prism.h"
#include "prism/measure.h"

/* Private functions */

static gatherPts (HisPoint *hp, Field *V[DIM+1], double *vbuf[_MAX_HP]);
static hisHeader (Domain *omega);

/*
 * Run-time Analyzer ... called every time step
 */

void Analyzer (Domain *d, double time, int step)
{
  const int verbose = option("verbose");
  const int measure = d->step.measure;
  const int cfl     = d->step.cfl;
  const int history = d->step.history;
  const int field   = d->step.field;
  const int stats   = d->step.stats;
  
  if (step == 0) hisHeader(d);
  
  /* .......... General Output ......... */

  Prism_log(1,"Time step = %d, Time = %-14.8g\n", step, time);

  if (verbose > 1) {        /* ..... Fluid State ..... */
    Divergence (d);     /* Total divergence        */
    Momentum   (d);     /* Total momentum          */
  }

  ROOTONLY {
    if (cfl && (step % cfl == 0))    
      estimate_CFL(d);
  }

  if (history && (step % history == 0))
    History (d, time);

  if (measure && (step % measure == 0))
    measure_analyze(d->mea);

  /* ........... User Defined .......... */

#if defined(USER)
  USER (d);
#endif

  if (step == 0)  /* Everything else is for step > 0 */
    return;

  if (stats && (step % stats == 0))
    stat_update (d->stats);

  /* Is checkpointing turned on?  If so, backup the checkpoint file  *
   * and open a new one.  There should always be at least two check- *
   * point files available for restarts.                             */

  if (step % field == 0) {
    if (option ("checkpt") && step < iparam("NSTEPS"))
      Domain_checkpt(d);
    else
      Domain_save(d);
  }
}

/* ------------------------------------------------------------------------- *
 * History() -- Process history points                                       *
 *                                                                           *
 * History points can be tracked either in Physical or Fourier space.  This  *
 * function processes the history point list and transforms data points if   *
 * necessary to physical space.                                              *
 * ------------------------------------------------------------------------- */

int History (Domain *d, double time)
{
  FILE     *fp = d->his_file;
  HisPoint *hp = d->his_list;
  Field    *V[DIM+1];
  double   *hbuf[_MAX_HP];
  int      i, n, npts, cnt;

  if (!hp) return 0;

  V[0]   = d->U; 
  V[1]   = d->V; 
  V[2]   = d->W;
  V[DIM] = d->P;
  npts   = gatherPts(hp, V, hbuf);
  cnt    = 0;
  
  ROOTONLY {
    do { 
      fprintf (fp, "%#13.8g ", time);
      for (n = 0; n < strlen(hp->fields); n++)
	fprintf (fp, "%#13.6g ", hbuf[cnt][n]);
      fprintf (fp, ":%d\n", cnt+1);
      free (hbuf[cnt++]);    } 
    while 
      (hp = hp->next);
  }

  return cnt;
}

/* Collect history points from the processors (?) */

static gatherPts (HisPoint *hp, Field *V[DIM+1], double *hbuf[_MAX_HP])
{
  const int pid    = option("procid");
  const int nprocs = option("nprocs");
  const int nel    = Field_count (V[0]);
  const int ntot   = FIELD_NR(V[0]) * FIELD_NS(V[0]) * nel;
  const int nz     = FIELD_NZ(V[0]);

  int i, k, n, nflds, pos, p;
  tempVector (tmp, nz*nprocs);

  for (i = 0; hp != NULL; i++, hp = hp->next) {
    ROOTONLY {
      nflds   = strlen(hp->fields);
      hbuf[i] = (double*) calloc(nflds,sizeof(double));
    }

    for  (n = 0, pos = 0; n <= DIM; n++) {
      if (strchr (hp->fields, V[n]->type)) {

	/* Sample the field and store the values in tmp[] */
	
	for (k = 0; k < nz; k++) {
	  Frame_set_one (k, V[n]);
	  tmp[pid*nz + k] = Probe_eval(hp->locator, V[n]);
	}
	  
#ifdef PARALLEL
#  define MSGTYPE(i) (1000+i)
	if (pid)
	  comm_send (MSGTYPE(pid), tmp+pid*nz, nz*sizeof(double), 0);
	else
	  for (p = 1; p < nprocs; p++)
	    comm_recv (MSGTYPE(p), tmp+p*nz, nz*sizeof(double));
#  undef MSGTYPE
# endif	

	ROOTONLY
	  switch (hp->mode) {
	  case Fourier:
	    hbuf[i][pos++] = tmp[hp->frame];
	    break;
	  case Physical:
	    realft(nz*nprocs/2, tmp, 1);
	    hbuf[i][pos++] = tmp[hp->frame];
	    break;
	  default:
	    Prism_error("history: unknown history point mode");
	    break;
	  }
      }
    }
  }
  
  freeVector(tmp);
  return i;
}

/* Write the header for the history point file */

static hisHeader (Domain *d)
{
  FILE     *fp = d->his_file;
  HisPoint *hp = d->his_list;
  Field    *U  = d->U;
  int       n  = 1;

  if (!fp) return 0;

  fputs ("# Prism history point file\n"
	 "# \n"
	 "# ID   Fields @    x       y    Element  Frame\n"
	 "# ----------------------------------------------------\n", fp);
  
  do {  
    fprintf (fp, "# %2d   %6s   %#7.4lf %#7.4lf   %3d    ",
	     hp->id, 
	     hp->fields, 
	     hp->locator->x, 
	     hp->locator->y, 
	     hp->locator->elmt->id+1);

    if (hp->mode == Physical)
      fprintf (fp, "%3d", hp->frame);
    else
#if DIM==2
      fprintf (fp, "%3d", hp->frame);
#else
      fprintf (fp, "%3d%c", hp->frame >> 1,
	       (hp->frame & 1 ? 'i' : 'r'));
#endif
    fputc ('\n', fp);
  } while
    (hp = hp->next);

  fputs ("# ----------------------------------------------------\n", fp);
  return 0;
}
    
/* ------------------------------------------------------------------------- *
 *                   F L O W F I E L D    A N A L Y S I S                    *
 * ------------------------------------------------------------------------- *

 * ------------------------------------------------------------------------- *
 * Divergence() -- Velocity field divergence                                 *
 *                                                                           *
 * Integrate the divergence of the velocity field over the entire computa-   *
 * tional domain.  Prints the total divergence and the maximum divergence,   *
 * returning the total.                                                      *
 * ------------------------------------------------------------------------- */

double Divergence (Domain *d)
{
  int     pid  = option("procid");
  int     nel  = Field_count (d->U), 
          nrns = d->U->nr * d->U->ns;
  double  div  = 0.,
          div2 = 0.;
  Field  *U, *Ux, *V, *Vy, *Div;
  register int k;
  
  if (pid) return 0.;

  U   =  d->U;  Ux  =  *d->Uf;
  V   =  d->V;  Vy  =  *d->Vf;
  Div = *d->Uf;

# if DIM == 3
  Frame_set (0, 4, U, V, Ux, Vy);
# endif

  Field_grad (U, Ux, 0);
  Field_grad (V, 0, Vy);
  dvadd (nrns * nel, *Ux->base, 1, *Vy->base, 1, *Div->base, 1);

  for (k = 0; k < nel; k++)
    div  += ddot (nrns, *Div[k].mass, 1, *Div[k].base, 1);

  k = idamax (nrns*nel, *Div->base, 1);
  Prism_log (2,"\tDivergence: %#.2g, L2 = %#.2g, max = %#.2g @ (%.2g,%.2g)\n", 
	     div, Field_L2(Div), (*Div->base)[k], (*Div->xmesh)[k], 
	     (*Div->ymesh)[k]);

  return div;
}

/* ------------------------------------------------------------------------- *
 * Momentum() -- Fluid Momentum calculation                                  *
 *                                                                           *
 * Integrates the fluid momentum in each direction over the entire computa-  *
 * tional domain.  Calculated in Fourier space.                              *
 *                                                                           *
 * Returns the total fluid momentum.                                         *
 * ------------------------------------------------------------------------- */

double Momentum (Domain *d)
{
  const int pid = option("procid");

  int    nrns, nel;
  double mom, f[3];
  Field  *Ux, *Uy, *Uz;
  int k;

  if (pid) return 0.;

  Ux = d->U;
  Uy = d->V;
  Uz = d->W;

  nrns  = Ux->nr * Ux->ns;
  nel   = Field_count(Ux);

  dzero (DIM, f, 1);
  for (k = 0; k < nel; k++) {
    f[0] += ddot (nrns, *Ux[k].mass, 1, *Ux[k].base, 1);
    f[1] += ddot (nrns, *Uy[k].mass, 1, *Uy[k].base, 1);
#if DIM  == 3
    f[2] += ddot (nrns, *Uz[k].mass, 1, *Uz[k].base, 1);
#endif
  }

  Prism_log(2, "\tMomentum  : ");
  for (k = 0; k < DIM; k++)
    Prism_log(2, "%-#10.4g ", f[k]);
  Prism_log(2, "\n");

  return sqrt(ddot(DIM,f,1,f,1));
}

/* ------------------------------------------------------------------------- *
 * Vorticity() - Fluid Vorticity calculation                                 *
 *                                                                           *
 * Compute the vorticity, and store it in level 0 of the velocity multi-step *
 * arrays.  This function also uses level 0 of the force multi-step arrays   *
 * as workspace, and it resets the values D.{Qx,Qy,Qz} to point to the   *
 * correct location for accessing the vorticity.                             *
 * ------------------------------------------------------------------------- */

void Vorticity (Domain *d)
{
  BSystem *B   =  d->Velocity;
  Field   *Ux  =  d->U,  *Uf = *(d->Uf),  *Qx;
  Field   *Uy  =  d->V,  *Vf = *(d->Vf),  *Qy;
  Field   *Uz  =  d->W,  *Wf = *(d->Wf),  *Qz;
  const int ntot =  Ux->nr * Ux->ns * B->elements;
  const int nz   =  Ux->nz;

#if DIM == 2
  Qz = d->Qz = *(d->Us);   /* Reset the pointer for Qz */

  Field_grad (Ux, 0, Uf);
  Field_grad (Uy, Vf, 0);

  dvsub (ntot, *Vf->base, 1, *Uf->base, 1, *Qz->base, 1);
  Field_davg(Qz, B);
#else
  Qx = d->Qx = *(d->Us);
  Qy = d->Qy = *(d->Vs);
  Qz = d->Qz = *(d->Ws);

  Field_grad_3D (Ux, 0,  Qz, *Vf->base, *Wf->base);
  Field_grad_3D (Uy, Uf,  0, *Vf->base, *Wf->base);
  Field_grad_3D (Uz, Qy, Qx, *Vf->base, *Wf->base);

  dvsub   (ntot*nz, *Uf->base, 1, *Qz->base, 1, *Qz->base, 1);
  Field_gradz   (Ux, Uf);
  dvsub   (ntot*nz, *Uf->base, 1, *Qy->base, 1, *Qy->base, 1);
  Field_gradz   (Uy, Uf);
  dvsub   (ntot*nz, *Qx->base, 1, *Uf->base, 1, *Qx->base, 1);

  Field_davg_3D (Qx, B);
  Field_davg_3D (Qy, B);
  Field_davg_3D (Qz, B);
#endif

  return;
}
