/* ------------------------------------------------------------------------- *
 * StokesBC() - Calculate the high-order boundary conditions for Stokes flow *
 *                                                                           *
 * This routine simply computes the curl of the vorticity to be used in the  *
 * high-order pressure boundary conditions and sets the non-linear terms to  *
 * zero.                                                                     *
 *                                                                           *
 * RCS Information                                  
 * ---------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 * ------------------------------------------------------------------------- */

#include "prism/prism.h"
#include "veclib/veclib.h"

#if DIM == 3

void StokesBC (Domain *omega)
{
  const int    nz    = omega->U->nz,
           ntot  = omega->U->nr * omega->U->ns * Field_count(omega->U),
           ntotz = ntot * nz;
  register int k;
  
  if (option("need.vorticity")) 
    Vorticity(omega);
  ComputePBCs(omega);


  /* Set the non-linear terms to zero */

  dzero (ntotz, *(omega->Uf[0])->base, 1);
  dzero (ntotz, *(omega->Vf[0])->base, 1);
  dzero (ntotz, *(omega->Wf[0])->base, 1);

  return;
}

#else  /* -----------------  2-D Boundary Terms  ------------------- */

void StokesBC (Domain *omega)
{
  const int ntot = omega->U->nr * omega->U->ns * Field_count(omega->U);
  
  if (option("need.vorticity")) 
    Vorticity(omega);
  ComputePBCs(omega);
  

  /* Set the non-linear terms to zero */

  dzero (ntot, *(omega->Uf[0])->base, 1);
  dzero (ntot, *(omega->Vf[0])->base, 1);

  return;
}

#endif


