/*===========================================================================
 * RCS Information:
 * ----------------
 * $Author$
 * $Date$
 * $Source$
 * $Revision$
 *===========================================================================*/

#include "iso.h"


#define SIZE 8


main()
/*===========================================================================*/
/* Test real DFT of SIZE*SIZE*SIZE data: do DFT/IDFT on first element of IC. */
/*===========================================================================*/
{
  int                   c, i, j, k, Npts, TabLen;
  CVF  IC;
  real**      head;
  int*               Dimension;
  complex*               Wtab;
  real                 DFTF;


  Dimension = ivect(1, 3);
  Dimension[1] = (Dimension[2] = SIZE);
  Dimension[3] = SIZE/2;
  Npts = Dimension[1]*Dimension[2]*Dimension[3];
  printf("Dimensions: %d, %d, %d\n", Dimension[1], Dimension[2], Dimension[3]);
  
  head = cfield(Dimension, &IC);

  TabLen = SIZE/2;
  Wtab = cvect(0, TabLen-1);
  preFFT(TabLen, Wtab);

  for (i=0; i<Npts*2; i++)
    for (c=1; c<=3; c++)
      head[c][i] = ((c-1)*512 + 1) + i;

  DFTF = 1.0 / Npts;
  for (i=0; i<Npts*2; i++)
    head[1][i] *= DFTF;
  rc3DFT(IC[1], Dimension, Wtab, FORWARD);
  rc3DFT(IC[1], Dimension, Wtab, INVERSE);

  for (i=0; i<Npts*2; i++)
    printf("%f\t%f\t%f\n", head[1][i], head[2][i], head[3][i]);

  exit(0);
}


  
