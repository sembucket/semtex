/*****************************************************************************
 * fourier.c
 *
 * 1D Fourier transform routines for real data fields based on FFTPACK
 * or Temperton FFT routines, or vendor-supplied alternatives.
 * NB: different restrictions may apply to input args depending on
 * selected routine. 
 *****************************************************************************/

static char
RCSid_fourier[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <femdef.h>
#include <alplib.h>
#include <femlib.h>


void sDFTr (float*        data,
	    const integer tlen,
	    const integer ntrn,
	    const integer sign)
/* ------------------------------------------------------------------------- *
 * Carry out multiple 1D single--complex Fourier transforms of data.
 * Data is to be Fourier transformed in the direction normal to the most
 * rapid traverse through memory, with sucessive points in the transform
 * separated by ntrn.  Data has a total of tlen * ntrn real points.
 *
 * Input parameters:
 * ----------------
 * tlen: number of real data in each transform, product of prime numbers.
 * ntrn: number of transforms to perform (also skip in data).
 * sign: transform direction: +1 ==> r-->c, -1 ==> c-->r.
 *
 * Notes:
 * -----
 * (1) Data are scaled/normalized with 1/tlen when sign is +1, so that
 *     the zeroth Fourier mode contains the spatial average value.
 * (2) After forward (r-->c) transform, data are ordered so that within
 *     each transform, the zeroth mode datum comes first.  Then ordering
 *     depends on tlen (even or odd): if even, the zeroth mode is followed
 *     by the real datum from the maximum frequency mode, after which
 *     the real and imaginary parts for each mode alternate (both cases).
 * ------------------------------------------------------------------------- */
{
  char             routine[] = "sDFTr";
  const integer    ntot = tlen * ntrn;
  register integer i;
  integer          dum, ip, iq, ir, ipqr2;
  register float  *work, *Wtab, *ptr;

  if (tlen < 2 || !ntrn) return;

#if 1
  
  /* -- Temperton FFT routine is default. */

  dum = tlen;
  prf235 (&dum, &ip, &iq, &ir, &ipqr2);
  
  if (!dum    ) message (routine, "tlen needs prime factors 2, 3, 5", ERROR);
  if (ntrn & 1) message (routine, "ntrn must be even",                ERROR);

  work = svector (0, ntot + ipqr2 - 1);
  Wtab = work + ntot;

  ssetpf (Wtab, tlen, ip, iq, ir);
  smpfft (data, work, ntrn, tlen, ip, iq, ir, Wtab, sign);

  if (sign == +1) sscal (ntot, 1.0 / tlen, data, 1);

  freeSvector (work, 0);

#else
  work = svector (0, 3 * tlen + 14);
  Wtab = work + tlen;
  ptr  = data;

  srffti (tlen, Wtab);

  switch (sign) {

  case +1:
    if (tlen & 1) {
      for (i = 0; i < ntrn; i++, ptr++) {
	scopy  (tlen, ptr, ntrn, work, 1);
	srfftf (tlen, work, Wtab);
	scopy  (tlen, work, 1, ptr, ntrn);
      }
    } else {
      for (i = 0; i < ntrn; i++, ptr++) {
	scopy  (tlen, ptr, ntrn, work, 1);
	srfftf (tlen, work, Wtab);
	scopy  (tlen - 2, work + 1, 1, ptr + 2 * ntrn, ntrn);
	ptr[0]    = work[0];
	ptr[ntrn] = work[tlen - 1];
      }
    }
    sscal (ntot, 1.0 / tlen, data, 1);
    break;

  case -1:
    if (tlen & 1) {
      for (i = 0; i < ntrn; i++, ptr++) {
	scopy  (tlen, ptr, ntrn, work, 1);
	srfftb (tlen, work, Wtab);
	scopy  (tlen, work, 1, ptr, ntrn);
      }
    } else {
      for (i = 0; i < ntrn; i++, ptr++) {
	work[tlen - 1] = ptr[ntrn];
	work[0]        = ptr[0];
	scopy  (tlen - 2, ptr + 2 * ntrn, ntrn, work + 1, 1);
	srfftb (tlen, work, Wtab);
	scopy  (tlen, work, 1, ptr, ntrn);
      }
    }
    break;

  default:
    message (routine, "illegal direction flag", ERROR);
    break;
  }
  
  freeSvector (work, 0);
#endif
}


void dDFTr (double*       data,
	    const integer tlen,
	    const integer ntrn,
	    const integer sign)
/* ------------------------------------------------------------------------- *
 * Double-precision version of sDFTr.
 * ------------------------------------------------------------------------- */
{
  char             routine[] = "dDFTr";
  const integer    ntot = tlen * ntrn;
  register integer i;
  integer          dum, ip, iq, ir, ipqr2;
  register double  *work, *Wtab, *ptr;

  if (tlen < 2 || !ntrn) return;

#if 1

  /* -- Temperton FFT routine is default. */

  dum = tlen;
  prf235 (&dum, &ip, &iq, &ir, &ipqr2);
  
  if (!dum    ) message (routine, "tlen needs prime factors 2, 3, 5", ERROR);
  if (ntrn & 1) message (routine, "ntrn must be even",                ERROR);

  work = dvector (0, ntot + ipqr2 - 1);
  Wtab = work + ntot;

  dsetpf (Wtab, tlen, ip, iq, ir);
  dmpfft (data, work, ntrn, tlen, ip, iq, ir, Wtab, sign);
  if (sign == +1) dscal (ntot, 1.0 / tlen, data, 1);

  freeDvector (work, 0);

#else
  work = dvector (0, 3 * tlen + 14);
  Wtab = work + tlen;
  ptr  = data;

  drffti (tlen, Wtab);

  switch (sign) {

  case +1:
    if (tlen & 1) {
      for (i = 0; i < ntrn; i++, ptr++) {
	dcopy  (tlen, ptr, ntrn, work, 1);
	drfftf (tlen, work, Wtab);
	dcopy  (tlen, work, 1, ptr, ntrn);
      }
    } else {
      for (i = 0; i < ntrn; i++, ptr++) {
	dcopy  (tlen, ptr, ntrn, work, 1);
	drfftf (tlen, work, Wtab);
	dcopy  (tlen - 2, work + 1, 1, ptr + 2 * ntrn, ntrn);
	ptr[0]    = work[0];
	ptr[ntrn] = work[tlen - 1];
      }
    }
    dscal (ntot, 1.0 / tlen, data, 1);
    break;

  case -1:
    if (tlen & 1) {
      for (i = 0; i < ntrn; i++, ptr++) {
	dcopy  (tlen, ptr, ntrn, work, 1);
	drfftb (tlen, work, Wtab);
	dcopy  (tlen, work, 1, ptr, ntrn);
      }
    } else {
      for (i = 0; i < ntrn; i++, ptr++) {
	work[tlen - 1] = ptr[ntrn];
	work[0]        = ptr[0];
	dcopy  (tlen - 2, ptr + 2 * ntrn, ntrn, work + 1, 1);
	drfftb (tlen, work, Wtab);
	dcopy  (tlen, work, 1, ptr, ntrn);
      }
    }
    break;

  default:
    message (routine, "illegal direction flag", ERROR);
    break;
  }
  
  freeDvector (work, 0);
#endif
}
