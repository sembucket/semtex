/*****************************************************************************
 * FFT.c: routines for real--complex 3D FFTs
 *
 * $Id$
 *****************************************************************************/

#include "iso.h"


static void  rcZFFT (CF,  const int, const int, const int,
		     const complex*, const int, const int);

static void  pcYFFT (CF, const int, const int, const int);
static void  pcXFFT (CF, const int, const int, const int);

static void  cZFFT  (CF,  const int, const int, const int,
		     const complex*, const int, const int);
static void  cYFFT  (CF,  const int, const int, const int,
		     const complex*, const int, const int);
static void  cXFFT  (CF,  const int, const int, const int,
		     const complex*, const int, const int);

static int N0 = 0, N1 = 0, N2 = 0; /* -- Globals used in 3D indexing. */

#define SWAP(a, b)  W = (a); (a) = (b); (b) = W
#define rm(i,j,k)   ((k) + N2 * ((j) + N1 * (i)))


void preFFT (complex*  Wtab,
	     const int K   )
/* ------------------------------------------------------------------------- *
 * Make angular factors for FFT, put in Wtab[0..N-1].
 * N is the HALF length of FFT; i.e., half the number of complex values in
 * any direction.
 * SIGN is -1 for exp(-ikx) on forward transform (the opposite to Numerical
 * Recipes; define it as 1 for their definition).
 * ------------------------------------------------------------------------- */
{
#define SIGN -1

  int     i;
  double  theta;

  Wtab[0].Re = 1.0;   Wtab[0].Im = 0.0;

  for (i = 1; i < K; i++) {
    theta = i * M_PI / (double) K;
    Wtab[i].Re =        cos(theta);
    Wtab[i].Im = SIGN * sin(theta);
  }

#undef SIGN
}


void preShift (complex*  Stab,
	       const int N   )
/* ------------------------------------------------------------------------- *
 * Make angular factors for shift in convolution.
 * N is the length of storage in any direction (2*K).
 * ------------------------------------------------------------------------- */
{
  int    i;
  double theta;

  Stab[0].Re = 1.0;   Stab[0].Im = 0.0;

  for (i = 1; i < N; i++) {
    theta = i * M_PI / (double) N;
    Stab[i].Re = cos(theta);   Stab[-i].Re =  Stab[i].Re;
    Stab[i].Im = sin(theta);   Stab[-i].Im = -Stab[i].Im;
  }
}


void rc3DFT(/* update */ CF              Data   ,
	    /* using  */ const int*      Dim    ,
	                 const complex*  Wtab   ,
                         const int       Forward)
/* ------------------------------------------------------------------------- *
 * Perform 3-dimensional FFT of real-stored-as-complex data, or its inverse.
 *                                                                          
 * Data is a 3-dimensional complex matrix, for which all the indices are    
 *   assumed to begin at ZERO.                                                
 * Dim is an int vector[1..3], in which the three values give the     
 *   dimensions of Data[0..Dim[1]-1,0..Dim[2]-1,0..Dim[3]-1].
 * Note that the number of complex data in the k-direction is Dim[3]  
 * (typically HALF the numbers in the i & j directions for a Fourier cube). 
 * Forward = 1 ==> forward transform, 0 ==> inverse.                        
 *                                                                          
 * Since the transforms of real data are conjugate-symmetric, half the full-
 * complex data would be redundant, and we take advantage by using ideas    
 * developed for the transforms of real data.  (See the discussions for that
 * go with the procedures twofft() and realft() in Chapter 12 of Numerical  
 * Recipes.)  To be able to easily implement the ideas involved, we exploit 
 * the interpretation of a 3-D DFT as a sequence of 1-D DFTs applied in the 
 * three directions in turn (see The Fast Fourier Transform by E.O. Brigham).
 *                                                                          
 * The storage scheme we use (contiguous storage allocated by cbox()) has   
 * "storage by row", that is, as we move through locations in memory, it is 
 * the third (k) index which varies most rapidly, then the second (j), and  
 * finally the first (i).  Therefore it is most natural to consider the     
 * initial real-to-complex packing as occuring in the k-direction, and so do
 * the real-data DFTs first, in the k-direction, at the start of our forward
 * transform.  Subsequently, we have to take account of the fact that the   
 * data on the k=0 face have come from the real parts of k-direction trans- 
 * forms at two different frequencies: k-0 and k-Nyquist (the data resulting
 * from the k-transforms are both purely real at these two frequencies, and 
 * thus we can pack the real-only k-Nyquist data into the imaginary parts of
 * the zero-k-frequency storage).                                           
 *                                                                          
 * When we then go ahead to do the j & i direction transforms, these can    
 * proceed as transforms of complex data everywhere except on the k=0 face, 
 * where special account has to be taken of the packing just described.     
 * On the k=0 face, the j-direction transforms of the complex data can be   
 * unpacked as the simultaneous transforms of two sets of real data (see the
 * discussion of twofft() provided in Press et al.).  We do things just a   
 * little differently than they do, choosing to pack the positive-frequency 
 * parts of the two transforms into the single complex buffer.              
 * With this scheme, the j-Nyquist part of the first transform (data which  
 * came from the k=0 face after k-direction DFT) gets packed into the j-0.Im
 * storage location.  This leaves the entire j-upper-half of the storage    
 * free, and into that we put the positive-frequency part of the transform  
 * of the data which came from the k-Nyquist frequencies of the k-direction 
 * transforms, but all of this is reflected, so that along the j-maximum    
 * line are packed the j-0 and j-Nyquist real parts of the j-direction      
 * transforms.  Got all that?  Next the i-direction transforms are done, and
 * we are free to do them all as full-complex DFTs except along the j-0 and 
 * j-max lines, where we must again interpret what is stored there as the   
 * packed parts of two sets of real data, transform and unpack with the same
 * reflective scheme as was just described in the j-direction, k=0 face.    
 *                                                                          
 * The inverse transform (Forwards == 0) inverts the procedure above.       
 * Note that as written, no normalization is implemented, so that it can be 
 * done either after the forward or inverse transform, as preferred.        
 * ------------------------------------------------------------------------- */
{
  register int i, j, k;
  const int    N      = Dim[1];
  const int    Non2   = Dim[3];
  const int    Nm     = N - 1;
  const int    TabLen = Non2;

  N0 = N1 = N;			/* -- Set global vars used for indexing. */
  N2 = Non2;

  if (Forward) {

    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
	rcZFFT (Data, i, j, Non2, Wtab, TabLen, FORWARD);
    
    for (i = 0; i < N; i++) {
      cYFFT  (Data, i, 0, N, Wtab, TabLen,      FORWARD);
      pcYFFT (Data, i, N,                       FORWARD);
      for (k = 1; k < Non2; k++)
	cYFFT (Data, i, k, N, Wtab, TabLen,     FORWARD);
    }
    
    cXFFT  (Data, 0, 0, N, Wtab, TabLen,        FORWARD);
    pcXFFT (Data, 0, N,                         FORWARD);
    cXFFT  (Data, Nm, 0, N, Wtab, TabLen,       FORWARD);
    pcXFFT (Data, Nm, N,                        FORWARD);
    for (j = 1; j < Nm; j++)
      cXFFT (Data, j, 0, N, Wtab, TabLen,       FORWARD);
    for (j = 0; j < N; j++)
      for (k = 1; k < Non2; k++)
	cXFFT (Data, j, k, N, Wtab, TabLen,     FORWARD);
    
  } else {
    
    for (j = 0; j < N; j++)
      for (k = 1; k < Non2; k++)
	cXFFT (Data, j, k, N, Wtab, TabLen,     INVERSE);
    for (j = 1; j < Nm; j++)
      cXFFT (Data, j, 0, N, Wtab, TabLen,       INVERSE);       
    pcXFFT (Data, 0,  N,                        INVERSE);
    cXFFT  (Data, 0, 0, N, Wtab, TabLen,        INVERSE);
    pcXFFT (Data, Nm, N,                        INVERSE);
    cXFFT  (Data, Nm, 0, N, Wtab, TabLen,       INVERSE);
    
    for (i = 0; i < N; i++) {
      pcYFFT (Data, i, N,                       INVERSE);
      cYFFT  (Data, i, 0, N, Wtab, TabLen,      INVERSE);
      for (k = 1; k < Non2; k++)
	cYFFT (Data, i, k, N, Wtab, TabLen,     INVERSE);
    }
    
    for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
	rcZFFT (Data, i, j, Non2, Wtab, TabLen, INVERSE);
  }
}


void scaleFT (/* update */ CF         Data,
	      /* using  */ const int* Dim )
/* ------------------------------------------------------------------------- *
 * Scale data for FFT: the convention adopted is that physical
 * variables are the inverse transform without scaling, but the forward
 * transform incorporates a factor of 1 / Npts, where Npts is the number
 * of complex data carried into the transform (1/2 the number of real data).
 * ------------------------------------------------------------------------- */
{
  const    int   Npts     = Dim[1] * Dim[2] * Dim[3];
  const    int   NpN      = Npts + Npts;
  register real  DFT_Fact = 1.0  / Npts;
  register real *u        = &Data[0][0][0].Re;
  register int   i;

  for (i = 0; i < NpN; i++) u[i] *= DFT_Fact;
}


int ispow2 (int num)
/* ------------------------------------------------------------------------- *
 * Is num a strictly positive integer power of two?  Yes ==> 1, No ==> 0.
 * ------------------------------------------------------------------------- */
{
  while (num % 2 == 0 && num > 2) num /= 2;

  if (num != 2) return 0; else return 1;
}


/* ################# Remaining routines have file scope. ################## */


void rcZFFT (CF             Data,
	     const int      i,
	     const int      j,
	     const int      Nz,
	     const complex* Wtab,
	     const int      TabLen,
	     const int      Forward)
/* ------------------------------------------------------------------------- *
 * Perform FFT along the third dimension of a 3-D array of real-stored-as-
 * complex data.  Forward = 1 ==> forward FFT, 0 ==> inverse.
 *
 * The array Data is assumed to be indexed starting at zero in each
 * direction.  Values i, j, supply indices (invariant during transform) in
 * the first and second directions in Data.  Nz specfies the number of
 * complex values in the third direction.  The storage scheme used for
 * complex values has the real and imaginary parts alternating as the third
 * direction in Data is traversed, so the actual number of real values in
 * the third direction is 2*Nz.
 *
 * The code is a straight-forward adaption of procedure realft() given by
 * Press et al. in Numerical Recipes.
 *
 * No normalization is carried out; divide by Nz as appropriate.
 * The unsymmetrical scaling values on forward & inverse transforms
 * are needed to achieve the expected magnitudes in the PHYSICAL and
 * FOURIER domains.
 * ------------------------------------------------------------------------- */
{
  register int	    _i, _j, revi, k, revk;
  register real     h1r, h1i, h2r, h2i;
  register complex* data = &Data[0][0][0];
  const    int      Non2 = Nz >> 1;
  
  if (Forward) {
    cZFFT (Data, i, j, Nz, Wtab, TabLen, FORWARD);
    for (k = 1; k < Non2; k++) {
      revk = Nz - k;

      _i = rm (i, j, k); _j = rm (i, j, revk);
      h1r =  0.25 * (data[_i].Re + data[_j].Re);
      h1i =  0.25 * (data[_i].Im - data[_j].Im);
      h2r =  0.25 * (data[_i].Im + data[_j].Im);
      h2i = -0.25 * (data[_i].Re - data[_j].Re);
    
      data[_i].Re =  h1r + Wtab[k].Re * h2r - Wtab[k].Im * h2i;
      data[_i].Im =  h1i + Wtab[k].Re * h2i + Wtab[k].Im * h2r;
      data[_j].Re =  h1r - Wtab[k].Re * h2r + Wtab[k].Im * h2i;
      data[_j].Im = -h1i + Wtab[k].Re * h2i + Wtab[k].Im * h2r;
    }

    _i = rm (i, j, Non2);
    data[_i].Re =  0.5 * data[_i].Re;
    data[_i].Im = -0.5 * data[_i].Im;

    _i = rm (i, j, 0);
    h1r = data[_i].Re;
    data[_i].Re =  0.5 * (h1r + data[_i].Im);
    data[_i].Im =  0.5 * (h1r - data[_i].Im);

  } else {
    for (k = 1; k < Non2; k++) {
      revk = Nz - k;

      _i = rm (i, j, k); _j = rm (i, j, revk);
      h1r =  (data[_i].Re + data[_j].Re);
      h1i =  (data[_i].Im - data[_j].Im);
      h2r = -(data[_i].Im + data[_j].Im);
      h2i =  (data[_i].Re - data[_j].Re);
    
      data[_i].Re =  h1r + Wtab[k].Re * h2r + Wtab[k].Im * h2i;
      data[_i].Im =  h1i + Wtab[k].Re * h2i - Wtab[k].Im * h2r;
      data[_j].Re =  h1r - Wtab[k].Re * h2r - Wtab[k].Im * h2i;
      data[_j].Im = -h1i + Wtab[k].Re * h2i - Wtab[k].Im * h2r;
    }

    _i = rm (i, j, Non2);
    data[_i].Re =  2.0 * data[_i].Re;
    data[_i].Im = -2.0 * data[_i].Im;

    _i = rm (i, j, 0);
    h1r = data[_i].Re;
    data[_i].Re = h1r + data[_i].Im;
    data[_i].Im = h1r - data[_i].Im;
    cZFFT (Data, i, j, Nz, Wtab, TabLen, INVERSE);

  }
}

 
void cZFFT (CF             Data   ,
	    const int      i      ,
	    const int      j      ,
	    const int      Nz     ,
	    const complex* Wtab   ,
	    const int      TabLen ,
	    const int      Forward)
/* ------------------------------------------------------------------------- *
 * Perform FFT along the third dimension of a 3-D array of complex data.
 * Forward = 1 ==> forward FFT, 0 ==> inverse.
 *
 * The array Data is assumed to be indexed starting at zero in each
 * direction.  Values i, j, supply indices (invariant during transform) in
 * the first and second directions in Data.  Nz specfies the number of
 * complex values in the third direction.
 *
 * The code is an adaption of procedure four1() given by
 * Press et al. in Numerical Recipes.  Some slight modifications change the
 * sign of the angles used in forward or reverse transform (to recover their
 * definition, define SIGN as 1), and start the indexing of data at 0,
 * rather than 1.
 *
 * No normalization is carried out; divide by Nz as appropriate.
 * ------------------------------------------------------------------------- */
{
  register int	   _i, _j, mmax, m, p, q, s, t, tstep;
  register real    tempr, tempi;
  register complex W, *data = &Data[0][0][0];
  const int        Non2 = Nz >> 1;
  const int        Nm   = Nz  - 1;
  const real       s1   = (Forward) ? -1.0 :  1.0;
  const real       s2   = (Forward) ?  1.0 : -1.0;

  /* -- Bit reversal. */

  s = 0;
  for (t = 0; t < Nm; t++) {
    if (s > t) {
      _i = rm (i, j, s); _j = rm (i, j, t);
      SWAP (data[_i], data[_j]);
    }
    m = Non2;
    while (m <= s) {
      s -= m;
      m >>= 1;
    }
    s += m;
  }

  /* -- Butterfly. */

  mmax = 1;
  q = TabLen;
  while (Nz > mmax) {
    p     = 0;
    tstep = mmax << 1;
    for (m = 0; m < mmax; m++) {
      for (t = m; t < Nz; t += tstep) {
	s = t + mmax;

	_i = rm (i, j, s); _j = rm (i, j, t);
	tempr = Wtab[p].Re*data[_i].Re + s1 * Wtab[p].Im*data[_i].Im;
	tempi = Wtab[p].Re*data[_i].Im + s2 * Wtab[p].Im*data[_i].Re;
	data[_i].Re = data[_j].Re - tempr;
	data[_i].Im = data[_j].Im - tempi;
	data[_j].Re += tempr;
	data[_j].Im += tempi;
      }
      p += q;
    }
    mmax = tstep;
    q >>= 1;
  }
}

 
void pcYFFT (CF        Zbuf   ,
	     const int i      ,
	     const int Ny     ,
	     const int Forward)
/* ------------------------------------------------------------------------- *
 * This is a modification of the procedure packedFFTs() (which dealt with   
 * two conjugate-symmetric DFTs of real data packed into one complex buffer)
 * so that it does the same thing on the k=0 face of a 3-D DFT box.  In this
 * case (talking about forwards transform here), the real parts of the data 
 * would have derived from the zero-frequency part of a k-direction DFT,    
 * while the imaginary parts would have derived from the Nyquist frequencies
 * of the same transforms.  Since the transforms are on the k=0 face of the 
 * DFT box, only the i-coordinate will vary from one tranform to the next.  
 *                                                                          
 * Parameter Ny is the number of complex data which would exist in each of  
 * the two DFTs if they were not packed.  The +ve-j-frequency part of the   
 * first DFT is packed into the low j-half of Zbuf, with the real values    
 * from the j-Nyquist frequency stored in the imaginary part of the  zero-j-
 * frequency allocation.  The +ve-j-frequency part of the second j-direct-  
 * ion-DFT is stored in the second j-half of Zbuf, but reflected, so that   
 * the zero and Nyquist j-frequency data are located at the high j-end of   
 * Zbuf, with progressively higher j-frequencies towards the centre of Zbuf.
 * Zbuf must have j-length Ny, allocated as Zbuf[i][0..Ny-1][k]. The real   
 * and imaginary parts are packed in the k-direction.                       
 *                                                                          
 * To carry out the forward FFT of two lots of real data, the first lot of  
 * data are stored in the real locations in Zbuf[i][j][0], while the        
 * second lot of data are stored in the imaginary locations.                
 * This will already be the case if a k-direction transform of real data has
 * been carried out.  Forward YFFT Zbuf[i][j][0], then call pcYFFT(Zbuf,    
 * i, Ny, 1).  Unscrambling of the transform is carried out until the +ve-  
 * frequency parts of the DFTs are separated as described above.            
 *                                                                          
 * When using the inverse process, the two DFTs are assumed to be packed as 
 * described.  First call pcYFFT(Zbuf, i, Ny, 0), which does the           
 * mixing-up of the two DFTs, then do inverse FFT.  The two lots of real    
 * data are in the real and imaginary locations on the 0-k face of Zbuf.    
 *                                                                          
 * References: Numerical Recipes sect 12.3, Bendat & Piersol 1971 sect 9.84.
 * ------------------------------------------------------------------------- */
{
  register int     _i, _j, _k, j, revj;
  const int        Nyo2 = Ny >> 1;
  register complex A, B, *data = &Zbuf[0][0][0];
  real             s1, s2, s3;

  if (Forward) {	/* -- Take mixed-up DFTs & unscramble. */
    _i = rm (i, Nyo2, 0); _j = rm (i, 0, 0);
    s1 = data[_i].Re;
    s2 = data[_i].Im;
    s3 = data[_j].Im;

    for (j = Nyo2 - 1; j > 0; j--) {
      revj = Ny - j;
      _i = rm (i, j, 0); _j = rm (i, revj, 0); _k = rm (i, revj-1, 0);
      A = data[_i];
      B = data[_j];
      data[_i].Re = 0.5*(A.Re + B.Re);
      data[_i].Im = 0.5*(A.Im - B.Im);
      data[_k].Re = 0.5*(A.Im + B.Im);
      data[_k].Im = 0.5*(B.Re - A.Re);
    }

    _i = rm (i, 0, 0); _j = rm (i, Ny-1, 0);
    data[_i].Im = s1;
    data[_j].Re = s3;
    data[_j].Im = s2;

  } else {			/* -- Take DFTs of real data & scramble up. */

    _i = rm (i, Ny-1, 0); _j = rm (i, 0, 0);
    s1 = data[_i].Re;
    s2 = data[_i].Im;
    s3 = data[_j].Im;

    for (j = 1; j < Nyo2; j++) {
      revj = Ny - j;

      _i = rm (i, j, 0); _j = rm (i, revj-1, 0); _k = rm (i, revj, 0);
      A = data[_i];
      B = data[_j];
      data[_i].Re = A.Re - B.Im;
      data[_i].Im = A.Im + B.Re;
      data[_k].Re = A.Re + B.Im;
      data[_k].Im = B.Re - A.Im;
    }
    
    _i = rm (i, 0, 0); _j = rm (i, Nyo2, 0);
    data[_i].Im = s1;
    data[_j].Re = s3;
    data[_j].Im = s2;
  }
}
      
 
void cYFFT (CF             Data   ,
	    const int      i      ,
	    const int      k      ,
	    const int      Ny     ,
	    const complex* Wtab   ,
	    const int      TabLen ,
	    const int      Forward)
/* ------------------------------------------------------------------------- *
 * Perform FFT along the second dimension of a 3-D array of complex data.   
 * Forward = 1 ==> forward FFT, 0 ==> inverse.                              
 *                                                                          
 * The array Data is assumed to be indexed starting at zero in each         
 * direction.  Values i, k, supply indices (invariant during transform) in  
 * the firat and third directions in Data.                                  
 * Ny specfies the number of complex values in the second direction.        
 *                                                                          
 * No normalization is carried out; divide by Ny as appropriate.            
 *-------------------------------------------------------------------------- */
{
  register int	   _i, _j, mmax, m, p, q, s, t, tstep;
  register real    tempr, tempi;
  register complex W, *data = &Data[0][0][0];
  const int        Non2 = Ny >> 1;
  const int        Nym  = Ny  - 1;
  const real       s1   = (Forward) ? -1.0 :  1.0;
  const real       s2   = (Forward) ?  1.0 : -1.0;

  /* -- Bit reversal. */

  s = 0;
  for (t = 0; t < Nym; t++) {
    if (s > t) {
      _i = rm (i, s, k); _j = rm (i, t, k);
      SWAP (data[_i], data[_j]);
    }
    m = Non2;
    while (m <= s) {
      s -= m;
      m >>= 1;
    }
    s += m;
  }
  
  /* -- Butterfly. */

  mmax = 1;
  q = TabLen;
  while (Ny > mmax) {		
    p     = 0;
    tstep = mmax << 1;
    for (m = 0; m < mmax; m++) {
      for (t = m; t < Ny; t+=tstep) {
	s = t + mmax;
	_i = rm (i, s, k); _j = rm (i, t, k);
	tempr = Wtab[p].Re*data[_i].Re + s1 * Wtab[p].Im*data[_i].Im;
	tempi = Wtab[p].Re*data[_i].Im + s2 * Wtab[p].Im*data[_i].Re;
	data[_i].Re = data[_j].Re - tempr;
	data[_i].Im = data[_j].Im - tempi;
	data[_j].Re += tempr;
	data[_j].Im += tempi;
      }
      p += q;
    }
    mmax = tstep;
    q >>= 1;
  }
}
 

void pcXFFT (CF        Zbuf   ,
	     const int j      ,
	     const int Nx     ,
	     const int Forward)
/* ------------------------------------------------------------------------- *
 * This is a modification of the procedure packedFFTs() (which dealt with   
 * two conjugate-symmetric DFTs of real data packed into one complex buffer)
 * so that it does the same thing on the k=0 face of a 3-D DFT box.  In this
 * case (talking about forwards transform here), we will be doing the DFTs  
 * of data on the two j-outermost k=0 lines in the i-direction (i.e. the    
 * parameter j should only be given values 0 of Ny-1).  The data would have 
 * derived from packed j-direction transforms of real data, so that on the  
 * very outer-most edges of the DFT-box (which we're dealing with here), the
 * data running in the i-direction are the packed values from the j-0 and   
 * j-Nyquist frequencies.                                                   
 *                                                                          
 * Parameter Nx is the number of complex data which would exist in each of  
 * the two DFTs if they were not packed.  The +ve-i-frequency part of the   
 * first DFT is packed into the low i-half of Zbuf, with the real values    
 * from the i-Nyquist frequency stored in the imaginary part of the  zero-i-
 * frequency allocation.  The +ve-i-frequency part of the second i-direct-  
 * ion-DFT is stored in the second j-half of Zbuf, but reflected, so that   
 * the zero and Nyquist i-frequency data are located at the high i-end of   
 * Zbuf, with progressively higher i-frequencies towards the centre of Zbuf.
 * Zbuf must have i-length Nx, allocated as Zbuf[0..Nx-1][j][k]. The real   
 * and imaginary parts are packed in the k-direction.                       
 *                                                                          
 * To carry out the forward FFT of two lots of real data, the first lot of  
 * data are stored in the real locations in Zbuf[i][j][0], while the        
 * second lot of data are stored in the imaginary locations Zbuf[i][j][0].  
 * This will already be the case if a k-direction transform of real data has
 * been carried out.  Forward pcYFFT Zbuf[i][j][0], then call               
 * pcXFFT(Zbuf, j, Nx, 1).  Unscrambling of the transform is carried        
 * out until the +ve-frequency parts of the DFTs are separated as described 
 * above.                                                                   
 *                                                                          
 * When using the inverse process, the two DFTs are assumed to be packed as 
 * described.  First call pcXFFT(Zbuf, j, Nx, 0), which does the           
 * mixing-up of the two DFTs, then do inverse FFT.  The two lots of real    
 * data are in the real and imaginary locations on the 0-k face of Zbuf.    
 *                                                                          
 * References: Numerical Recipes sect 12.3, Bendat & Piersol 1971 sect 9.84.
 * ------------------------------------------------------------------------- */
{
  register int      _i, _j, _k, i, revi;
  const int         Nxo2 = Nx >> 1;
  register complex  A, B, *data = &Zbuf[0][0][0];
  real              s1, s2, s3;

  if (Forward) {	/* -- Take mixed-up DFTs & unscramble. */

    _i = rm (Nxo2, j, 0); _j = rm (0, j, 0);
    s1 = data[_i].Re;
    s2 = data[_i].Im;
    s3 = data[_j].Im;

    for (i = Nxo2-1; i > 0; i--) {
      revi = Nx - i;
      
      _i = rm (i, j, 0); _j = rm (revi, j, 0); _k = rm (revi-1, j, 0);
      A = data[_i];
      B = data[_j];
      data[_i].Re = 0.5*(A.Re + B.Re);
      data[_i].Im = 0.5*(A.Im - B.Im);
      data[_k].Re = 0.5*(A.Im + B.Im);
      data[_k].Im = 0.5*(B.Re - A.Re);
    }

    _i = rm (0, j, 0); _j = rm (Nx-1, j, 0);
    data[_i].Im = s1;
    data[_j].Re = s3;
    data[_j].Im = s2;

  } else {			/* -- Take DFTs of real data & scramble up. */

    _i = rm (Nx-1, j, 0); _j = rm (0, j, 0);
    s1 = data[_i].Re;
    s2 = data[_i].Im;
    s3 = data[_j].Im;

    for (i = 1; i < Nxo2; i++) {
      revi = Nx - i;

      _i = rm (i, j, 0); _j = rm (revi-1, j, 0); _k = rm (revi, j, 0);
      A = data[_i];
      B = data[_j];
      data[_i].Re = A.Re - B.Im;
      data[_i].Im = A.Im + B.Re;
      data[_k].Re = A.Re + B.Im;
      data[_k].Im = B.Re - A.Im;
    }
    
    _i = rm (0, j, 0); _j = rm (Nxo2, j, 0);
    data[_i].Im = s1;
    data[_j].Re = s3;
    data[_j].Im = s2;
  }
}
      
 
void cXFFT (CF             Data   ,
	    const int      j      ,
	    const int      k      ,
	    const int      Nx     ,
	    const complex* Wtab   ,
	    const int      TabLen ,
	    const int      Forward)
/* -------------------------------------------------------------------------- *
 * Perform FFT along the first dimension of a 3-D array of real-stored-as-  
 * complex data.  Forward = 1 specifies a forward FFT, 0 ==>  inverse.      
 *                                                                          
 * The array Data is assumed to be indexed starting at zero in each         
 * direction.  Values i, j, supply indices (invariant during transform) in  
 * the first and second directions in Data.  Nx specfies the number of      
 * complex values in the first direction.  The storage scheme used for      
 * complex values has the real and imaginary parts alternating as the third 
 * direction in Data is traversed.                                          
 *                                                                          
 * The code is an adaption of procedure four1() given by Press et al.
 * in Numerical Recipes.  A table of precomputed complex       
 * angular factors is used, and array indices start from 0, not 1.          
 *                                                                          
 * No normalization is carried out; divide by Nx as appropriate.            
 * ------------------------------------------------------------------------- */
{
  register int	   _i, _j, mmax, m, p, q, s, t, tstep;
  register real    tempr, tempi;
  register complex W, *data = &Data[0][0][0];
  const int        Non2 = Nx >> 1;
  const int        Nm   = Nx  - 1;
  const real       s1   = (Forward) ? -1.0 :  1.0;
  const real       s2   = (Forward) ?  1.0 : -1.0;

  /* -- Bit reversal. */

  s = 0;
  for (t = 0; t < Nm; t++) {
    if (s > t) {
      _i = rm(s,j,k); _j = rm(t,j,k);
      SWAP (data[_i], data[_j]);
    }
    m = Non2;
    while (m <= s) {
      s -= m;
      m >>= 1;
    }
    s += m;
  }

  /* -- Butterfly. */

  mmax = 1;
  q    = TabLen;
  while (Nx > mmax) {
    p     = 0;
    tstep = mmax << 1;
    for (m = 0; m < mmax; m++) {
      for (t = m; t < Nx; t += tstep) {
	s = t + mmax;
	_i = rm (s, j, k); _j = rm (t, j, k);
	tempr = Wtab[p].Re*data[_i].Re + s1 * Wtab[p].Im*data[_i].Im;
	tempi = Wtab[p].Re*data[_i].Im + s2 * Wtab[p].Im*data[_i].Re;
	data[_i].Re = data[_j].Re - tempr;
	data[_i].Im = data[_j].Im - tempi;
	data[_j].Re += tempr;
	data[_j].Im += tempi;
      }
      p += q;
    }
    mmax = tstep;
    q >>= 1;
  }
}
