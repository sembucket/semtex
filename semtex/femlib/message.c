/*****************************************************************************
 * message.c: message-passing routines, currently MPI-specific.
 *****************************************************************************/

static char
RCSid[] = "$Id$";

#include <stdio.h>
#include <malloc.h>
#include <femdef.h>
#include <femlib.h>

#if defined(MPI)
#include <mpi.h>
#endif


void message_init (int*    argc,
		   char*** argv)
/* ------------------------------------------------------------------------- *
 * Do whatever is required to initialize message-passing.  Set up global 
 * variables.  The parser must have been initialized already.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  int  n;
  char s[STR_MAX];

  MPI_Init      (argc, argv);

  MPI_Comm_rank (MPI_COMM_WORLD,   &n);
  sprintf       (s, "I_PROC = %1d", n);
  yy_interpret  (s);

  MPI_Comm_size (MPI_COMM_WORLD,   &n);
  sprintf       (s, "N_PROC = %1d", n);
  yy_interpret  (s);

#endif
}


void message_stop ()
/* ------------------------------------------------------------------------- *
 * Shut down message passing.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Barrier  (MPI_COMM_WORLD);
  MPI_Finalize ();

#endif
}


void message_sync ()
/* ------------------------------------------------------------------------- *
 * Block until all processes have entered.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Barrier (MPI_COMM_WORLD);

#endif
}


static int first (int n, const int* x)
{ 
  register int i;
  for (i = 0; i < n; i++) if (x[i]) return i;
  return 0;
}


void message_dtranspose (double*       data,
			 const integer nZ  ,
			 const integer nP  ,
			 const integer sign)
/* ------------------------------------------------------------------------- *
 * Transpose blocks of data across processors.  Data is a double-precision 
 * vector, total length nP*nZ on each processor.
 *
 * The amount of data held by each processor is nP*nZ.  Each nP-sized plane
 * can be split into nB = nP/nProc sized blocks (it is assumed that nP is an
 * integer multiple of nProc, i.e. that nB is a whole number).  Initially
 * the data are ordered by as nZ nP-sized planes, with memory traversed 
 * fastest moving over the planes, and the block indices vary more rapidly
 * than the z-indices.
 *
 * The aim of the the transpose is to re-order the data so that each
 * processor holds all the nZ data for one of the blocks, i.e. a gather of
 * all the z-information for a block onto each a single processor, which e.g.
 * can be followed by multiple 1D Fourier transformations over each block.
 *
 * First the data are transposed within a single processor so that (in terms
 * of blocks) the block (rather than the z) indices vary slowest as memory 
 * is traversed.  Then a block-transpose of data across processors is
 * carried out using message passing.
 *
 * NB: order of inter- and intra-processor transposes must be reversed 
 * in order to invert the transpose with a second transpose: this is the use
 * of input variable sign.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  register int i, j;
  const int    ip = (int) yy_interpret ("I_PROC");
  const int    np = (int) yy_interpret ("N_PROC");
  const int    nB = nP / np;	   /* Size of intra-processor block.     */
  const int    NB = nP / nB;	   /* Number of these blocks in a plane. */
  const int    NM = nP * nZ / np;  /* Size of message block.             */
  const int    dsize = sizeof (double);
  double*      tmp;
  MPI_Status   status;

  if (np == 1) return;

  tmp = (double*) malloc (nB * dsize);

  if (sign == 1) {		/* -- "Forwards" transpose. */

    /* -- Intra-processor transpose. */

    if (NB == nZ) {		/* -- Symmetric transpose. */
      for (i = 0; i < nZ; i++)
	for (j = i; j < nZ; j++) {
	  if (i != j) {
	    memcpy (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	    memcpy (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	    memcpy (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	  }
	}

    } else {			/* -- Asymmetric transpose. */

      int       k, knext, kconf, *kmove, *kpost;
      const int NBnZm = NB * nZ - 1;

      kmove = (int*) malloc (2*nZ*NB * sizeof (int));
      kpost = kmove + nZ*NB;

      /* -- Build scatter indices. */
    
      for (i = 0; i < NB; i++)
	for (j = 0; j < nZ; j++)
	  kpost[j * NB + i] = i * nZ + j;

      for (i = 1; i < NBnZm; i++) kmove[i] = 1;
      kmove[0] = kmove[NBnZm] = 0;

      /* -- Do "in-place" scatter. */

      while (k = first (nZ*NB, kmove)) {
	knext = kconf = kpost[k];
	memcpy (tmp, data + kconf * nB, nB * dsize);
	do {
	  memcpy (data + knext * nB, data + k * nB, nB * dsize);
	  kmove[k] = 0;
	  for (i = 1; i < NBnZm; i++)
	    if (kpost[i] == k) {
	      k     = i;
	      knext = kpost[k];
	      break;
	    }
	} while (k != kconf);
	memcpy (data + knext * nB, tmp, nB * dsize);
	kmove[k] = 0;
      }
      
      free (kmove);
    }

    /* -- Inter-processor transpose, with NB blocks size nZ*nB / processor. */

    for (i = 0; i < np; i++)
      if (i != ip) {
	MPI_Send (data + i*NM, NM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
	MPI_Recv (data + i*NM, NM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      }

  } else {			/* -- "Backwards" transpose. */

    for (i = 0; i < np; i++)
      if (i != ip) {
	MPI_Send (data + i*NM, NM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
	MPI_Recv (data + i*NM, NM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
      }

    if (NB == nZ) {
      for (i = 0; i < nZ; i++)
	for (j = i; j < nZ; j++) {
	  if (i != j) {
	    memcpy (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	    memcpy (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	    memcpy (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	  }
	}

    } else {
      int       k, knext, kconf, *kmove, *kpost;
      const int NBnZm = NB * nZ - 1;

      kmove = (int*) malloc (2*nZ*NB * sizeof (int));
      kpost = kmove + nZ*NB;
  
      for (i = 0; i < NB; i++)
	for (j = 0; j < nZ; j++)
	  kpost[i * nZ + j] = j * NB + i;

      for (i = 1; i < NBnZm; i++) kmove[i] = 1;
      kmove[0] = kmove[NBnZm] = 0;

      while (k = first (nZ*NB, kmove)) {
	knext = kconf = kpost[k];
	memcpy (tmp, data + kconf * nB, nB * dsize);
	do {
	  memcpy (data + knext * nB, data + k * nB, nB * dsize);
	  kmove[k] = 0;
	  for (i = 1; i < NBnZm; i++)
	    if (kpost[i] == k) {
	      k     = i;
	      knext = kpost[k];
	      break;
	    }
	} while (k != kconf);
	memcpy (data + knext * nB, tmp, nB * dsize);
	kmove[k] = 0;
      }
      
      free (kmove);
      
    }
  }

  free (tmp);
#endif
}


void message_stranspose (float*        data,
			 const integer nZ  ,
			 const integer nP  ,
			 const integer sign)
/* ------------------------------------------------------------------------- *
 * Single-precision version of message_dtranspose.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  register int i, j;
  const int    ip = (int) yy_interpret ("I_PROC");
  const int    np = (int) yy_interpret ("N_PROC");
  const int    nB = nP / np;
  const int    NB = nP / nB;
  const int    NM = nP * nZ / np;
  const int    dsize = sizeof (float);
  float*       tmp;
  MPI_Status   status;

  if (np == 1) return;

  tmp = (float*) malloc (nB * dsize);

  if (sign == 1) {		/* -- "Forward" transpose. */

    /* -- Intra-processor transpose. */
    
    if (NB == nZ) {		/* -- Symmetric transpose. */
      for (i = 0; i < nZ; i++)
	for (j = i; j < nZ; j++) {
	  if (i != j) {
	    memcpy (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	    memcpy (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	    memcpy (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	  }
	}

    } else {			/* -- Asymmetric transpose. */
      int       k, knext, kconf, *kmove, *kpost;
      const int NBnZm = NB * nZ - 1;

      kmove = (int*) malloc (2*nZ*NB * sizeof (int));
      kpost = kmove + nZ*NB;

      /* -- Build scatter indices. */
    
      for (i = 0; i < NB; i++)
	for (j = 0; j < nZ; j++)
	  kpost[j * NB + i] = i * nZ + j;

      for (i = 1; i < NBnZm; i++) kmove[i] = 1;
      kmove[0] = kmove[NBnZm] = 0;

      /* -- Do "in-place" scatter. */

      while (k = first (nZ*NB, kmove)) {
	knext = kconf = kpost[k];
	memcpy (tmp, data + kconf * nB, nB * dsize);
	do {
	  memcpy (data + knext * nB, data + k * nB, nB * dsize);
	  kmove[k] = 0;
	  for (i = 1; i < NBnZm; i++)
	    if (kpost[i] == k) {
	      k     = i;
	      knext = kpost[k];
	      break;
	    }
	} while (k != kconf);
	memcpy (data + knext * nB, tmp, nB * dsize);
	kmove[k] = 0;
      }
      
      free (kmove);
    }

    /* -- Inter-processor transpose, with NB blocks size nZ*nB / processor. */

    for (i = 0; i < np; i++)
      if (i != ip) {
	MPI_Send (data + i*NM, NM, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
	MPI_Recv (data + i*NM, NM, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
      }

  } else {			/* -- "Backwards" transpose. */

    for (i = 0; i < np; i++)
      if (i != ip) {
	MPI_Send (data + i*NM, NM, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
	MPI_Recv (data + i*NM, NM, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status);
      }

    if (NB == nZ) {
      for (i = 0; i < nZ; i++)
	for (j = i; j < nZ; j++) {
	  if (i != j) {
	    memcpy (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	    memcpy (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	    memcpy (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	  }
	}

    } else {
      int       k, knext, kconf, *kmove, *kpost;
      const int NBnZm = NB * nZ - 1;

      kmove = (int*) malloc (2*nZ*NB * sizeof (int));
      kpost = kmove + nZ*NB;
  
      for (i = 0; i < NB; i++)
	for (j = 0; j < nZ; j++)
	  kpost[i * nZ + j] = j * NB + i;

      for (i = 1; i < NBnZm; i++) kmove[i] = 1;
      kmove[0] = kmove[NBnZm] = 0;

      while (k = first (nZ*NB, kmove)) {
	knext = kconf = kpost[k];
	memcpy (tmp, data + kconf * nB, nB * dsize);
	do {
	  memcpy (data + knext * nB, data + k * nB, nB * dsize);
	  kmove[k] = 0;
	  for (i = 1; i < NBnZm; i++)
	    if (kpost[i] == k) {
	      k     = i;
	      knext = kpost[k];
	      break;
	    }
	} while (k != kconf);
	memcpy (data + knext * nB, tmp, nB * dsize);
	kmove[k] = 0;
      }
      
      free (kmove);
    }
  }

  free (tmp);
#endif
}


void message_itranspose (integer*      data,
			 const integer nZ  ,
			 const integer nP  ,
			 const integer sign)
/* ------------------------------------------------------------------------- *
 * Integer transpose.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  register int i, j;
  const int    ip = (int) yy_interpret ("I_PROC");
  const int    np = (int) yy_interpret ("N_PROC");
  const int    nB = nP / np;	   /* Size of intra-processor block.     */
  const int    NB = nP / nB;	   /* Number of these blocks in a plane. */
  const int    NM = nP * nZ / np;  /* Size of message block.             */
  const int    dsize = sizeof (integer);
  integer*     tmp;
  MPI_Status   status;

  if (np == 1) return;

  tmp = (integer*) malloc (nB * dsize);

  if (sign == 1) {		/* -- "Forward" transpose. */

    /* -- Intra-processor transpose. */

    if (NB == nZ) {		/* -- Symmetric transpose. */
      for (i = 0; i < nZ; i++)
	for (j = i; j < nZ; j++) {
	  if (i != j) {
	    memcpy (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	    memcpy (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	    memcpy (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	  }
	}

    } else {			/* -- Asymmetric transpose. */
      int       k, knext, kconf, *kmove, *kpost;
      const int NBnZm = NB * nZ - 1;

      kmove = (int*) malloc (2*nZ*NB * sizeof (int));
      kpost = kmove + nZ*NB;

      /* -- Build scatter indices. */
    
      for (i = 0; i < NB; i++)
	for (j = 0; j < nZ; j++)
	  kpost[j * NB + i] = i * nZ + j;

      for (i = 1; i < NBnZm; i++) kmove[i] = 1;
      kmove[0] = kmove[NBnZm] = 0;

      /* -- Do "in-place" scatter. */

      while (k = first (nZ*NB, kmove)) {
	knext = kconf = kpost[k];
	memcpy (tmp, data + kconf * nB, nB * dsize);
	do {
	  memcpy (data + knext * nB, data + k * nB, nB * dsize);
	  kmove[k] = 0;
	  for (i = 1; i < NBnZm; i++)
	    if (kpost[i] == k) {
	      k     = i;
	      knext = kpost[k];
	      break;
	    }
	} while (k != kconf);
	memcpy (data + knext * nB, tmp, nB * dsize);
	kmove[k] = 0;
      }
      
      free (kmove);
    }    
    
    /* -- Inter-processor transpose, with NB blocks size nZ*nB / processor. */

    if (sizeof (integer) == sizeof (int)) {
      for (i = 0; i < np; i++)
	if (i != ip) {
	  MPI_Send (data + i*NM, NM, MPI_INT, i, 0, MPI_COMM_WORLD);
	  MPI_Recv (data + i*NM, NM, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
	}
    } else {
      for (i = 0; i < np; i++)
	if (i != ip) {
	  MPI_Send (data + i*NM, NM, MPI_LONG, i, 0, MPI_COMM_WORLD);
	  MPI_Recv (data + i*NM, NM, MPI_LONG, i, 0, MPI_COMM_WORLD, &status);
	}
    }

  } else {			/* -- "Backwards" transpose. */

    if (sizeof (integer) == sizeof (int)) {
      for (i = 0; i < np; i++)
	if (i != ip) {
	  MPI_Send (data + i*NM, NM, MPI_INT, i, 0, MPI_COMM_WORLD);
	  MPI_Recv (data + i*NM, NM, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
	}
    } else {
      for (i = 0; i < np; i++)
	if (i != ip) {
	  MPI_Send (data + i*NM, NM, MPI_LONG, i, 0, MPI_COMM_WORLD);
	  MPI_Recv (data + i*NM, NM, MPI_LONG, i, 0, MPI_COMM_WORLD, &status);
	}
    }

    if (NB == nZ) {
      for (i = 0; i < nZ; i++)
	for (j = i; j < nZ; j++) {
	  if (i != j) {
	    memcpy (tmp,                data + (i*nZ+j)*nB, nB * dsize);
	    memcpy (data + (i*nZ+j)*nB, data + (j*nZ+i)*nB, nB * dsize);
	    memcpy (data + (j*nZ+i)*nB, tmp,                nB * dsize);
	  }
	}

    } else {
      int       k, knext, kconf, *kmove, *kpost;
      const int NBnZm = NB * nZ - 1;

      kmove = (int*) malloc (2*nZ*NB * sizeof (int));
      kpost = kmove + nZ*NB;
  
      for (i = 0; i < NB; i++)
	for (j = 0; j < nZ; j++)
	  kpost[i * nZ + j] = j * NB + i;

      for (i = 1; i < NBnZm; i++) kmove[i] = 1;
      kmove[0] = kmove[NBnZm] = 0;

      while (k = first (nZ*NB, kmove)) {
	knext = kconf = kpost[k];
	memcpy (tmp, data + kconf * nB, nB * dsize);
	do {
	  memcpy (data + knext * nB, data + k * nB, nB * dsize);
	  kmove[k] = 0;
	  for (i = 1; i < NBnZm; i++)
	    if (kpost[i] == k) {
	      k     = i;
	      knext = kpost[k];
	      break;
	    }
	} while (k != kconf);
	memcpy (data + knext * nB, tmp, nB * dsize);
	kmove[k] = 0;
      }
      
      free (kmove);
    }
  }

  free (tmp);
#endif
}


void message_dsend (double*       data,
		    const integer N   ,
		    const integer tgt )
/* ------------------------------------------------------------------------- *
 * Send data to processor number tgt.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Send (data, (int) N, MPI_DOUBLE, (int) tgt, 0, MPI_COMM_WORLD);

#endif
}


void message_drecv (double*       data,
		    const integer N   ,
		    const integer src )
/* ------------------------------------------------------------------------- *
 * Receive data from processor number src.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Status status;

  MPI_Recv (data, (int) N, MPI_DOUBLE, (int) src, 0, MPI_COMM_WORLD, &status);

#endif
}


void message_ssend (float*        data,
		    const integer N   ,
		    const integer tgt )
/* ------------------------------------------------------------------------- *
 * Send data to processor number tgt.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Send (data, (int) N, MPI_FLOAT, (int) tgt, 0, MPI_COMM_WORLD);

#endif
}


void message_srecv (float*        data,
		    const integer N   ,
		    const integer src )
/* ------------------------------------------------------------------------- *
 * Receive data from processor number src.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Status status;

  MPI_Recv (data, (int) N, MPI_FLOAT, (int) src, 0, MPI_COMM_WORLD, &status);

#endif
}


void message_isend (integer*      data,
		    const integer N   ,
		    const integer tgt )
/* ------------------------------------------------------------------------- *
 * Send data to processor number tgt.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  if (sizeof (integer) == sizeof (int))
    MPI_Send (data, (int) N, MPI_INT,  (int) tgt, 0, MPI_COMM_WORLD);
  else
    MPI_Send (data, (int) N, MPI_LONG, (int) tgt, 0, MPI_COMM_WORLD);

#endif
}


void message_irecv (integer*      data,
		    const integer N   ,
		    const integer src )
/* ------------------------------------------------------------------------- *
 * Receive data from processor number src.
 * ------------------------------------------------------------------------- */
{
#if defined(MPI)

  MPI_Status status;

  if (sizeof (integer) == sizeof (int))
    MPI_Recv (data, (int) N, MPI_INT,  (int) src, 0, MPI_COMM_WORLD, &status);
  else
    MPI_Recv (data, (int) N, MPI_LONG, (int) src, 0, MPI_COMM_WORLD, &status);

#endif
}

