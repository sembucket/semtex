/*****************************************************************************
 * fourier.c
 *
 * Fourier transform routines for real data fields based on FFTPACK.
 *****************************************************************************/

static char
RCSid_fourier[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <alplib.h>
#include <femdef.h>
#include <femlib.h>


void sDFTr (float*    data ,
	    const int len  ,
	    const int nopr ,
	    const int skips,
	    const int skipd,
	    const int sign )
/* ------------------------------------------------------------------------- *
 * Carry out sequential 1D single--complex Fourier transforms of data.
 *
 * Input parameters:
 * ----------------
 * len:   number of real data in each transform, product of prime numbers.
 * nopr:  number of transforms to perform.
 * skips: skip of starting point of successive transforms.
 * skipd: skip within data for each point in a single transform.
 * sign:  transform direction: +1 ==> r-->c, -1 ==> c-->r.
 *
 * Notes:
 * -----
 * (1) Data must contain at least len*skips*skipd points.
 * (2) All data are scaled by 1/len when sign is +1.
 * (3) Len must be an even number.
 * (4) After forward (r-->c) transform, data are ordered so that within
 *     each transform, the zeroth mode datum comes first, followed by the
 *     real datum from the maximum frequency (Nyquist) mode, after which
 *     the real and complex data for each mode alternate.
 * (5) The maximum frequency datum is always set to zero.
 * ------------------------------------------------------------------------- */
{
  char           routine[] = "sDFTr";
  register int   i;
  register float *work, *Wtab, *ptr;

  if (len < 2) return;
  if (len & 1) message (routine, "input transform length must be even", ERROR);

  work = svector (0, 3 * len + 14);
  Wtab = work + len;
  ptr  = data;

  srffti (len, Wtab);

  switch (sign) {

  case 1:
    for (i = 0; i < nopr; i++, ptr += skips) {
      scopy  (len, ptr, skipd, work, 1);
      srfftf (len, work, Wtab);
      scopy  (len - 2, work + 1, 1, ptr + 2 * skipd, skipd);
      ptr[0]     = work[0];
      ptr[skipd] = work[len - 1];
    }
    sscal (len * skipd * skips, 1.0 / len, data, 1);
    break;

  case -1:
    for (i = 0; i < nopr; i++, ptr += skips) {
      work[len - 1] = ptr[nopr];
      work[0]       = ptr[0];
      scopy  (len - 2, ptr + 2 * skipd, skipd, work + 1, 1);
      srfftb (len, work, Wtab);
      scopy  (len, work, 1, ptr, skipd);
    }
    break;

  default:
    message (routine, "illegal direction flag", ERROR);
    break;
  }
  
  freeSvector (work, 0);
}


void dDFTr (double*   data ,
	    const int len  ,
	    const int nopr ,
	    const int skips,
	    const int skipd,
	    const int sign )
/* ------------------------------------------------------------------------- *
 * Carry out sequential 1D double--zomplex Fourier transforms of data.
 * 
 * See remarks for sDFTr().
 * ------------------------------------------------------------------------- */
{
  char            routine[] = "dDFTr";
  register int    i;
  register double *work, *Wtab, *ptr;

  if (len < 2) return;
  if (len & 1) message (routine, "input transform length must be even", ERROR);

  work = dvector (0, 3 * len + 14);
  Wtab = work + len;
  ptr  = data;

  drffti (len, Wtab);

  switch (sign) {

  case 1:
    for (i = 0; i < nopr; i++, ptr += skips) {
      dcopy  (len, ptr, skipd, work, 1);
      drfftf (len, work, Wtab);
      dcopy  (len - 2, work + 1, 1, ptr + 2 * skipd, skipd);
      ptr[0]     = work[0];
      ptr[skipd] = work[len - 1];
    }
    dscal (len * skipd * skips, 1.0 / len, data, 1);
    break;

  case -1:
    for (i = 0; i < nopr; i++, ptr += skips) {
      work[len - 1] = ptr[nopr];
      work[0]       = ptr[0];
      dcopy  (len - 2, ptr + 2 * skipd, skipd, work + 1, 1);
      drfftb (len, work, Wtab);
      dcopy  (len, work, 1, ptr, skipd);
    }
    break;

  default:
    message (routine, "illegal direction flag", ERROR);
    break;
  }
  
  freeDvector (work, 0);
}

