/*****************************************************************************
 * Convert an ISO field file into a TECPLOT file, or a binary semtex
 * file.  Vorticity is computed and added as well.  Optionally leave
 * data in Fourier space.  If semtex output is selected, also create a
 * matching prism/semtex mesh file (input.msh).
 *
 * usage: iso2tec [options] input.fld
 * options:
 * -h ... print this message
 * -s ... create semtex/prism format output [Default: TECPLOT]
 * -f ... suppress inverse Fourier transformation
 *
 * : iso2tec.c,v 2.1 1995/11/08 00:42:04 hmb Exp hmb $
 *****************************************************************************/

#include "iso.h"

int N, K, FourKon3;

static const char prog[] = "iso2tec";

static void getargs   (int, char**, const char**, int* , int*);
static void enstrophy (CF, const CVF, const int*);
static void tecout    (const char*, const CVF, const CVF, const CF,
		       const int*, const int);
static void semout    (const char*, const CVF, const CVF, const CF,
		       const int*, const int);


int main (int    argc,
	  char** argv)
/* ------------------------------------------------------------------------- *
 * Driver.
 * ------------------------------------------------------------------------- */
{
  CVF          U;
  CVF          Q;
  CF           Work;
  Param*       Info = (Param*) calloc (1, sizeof (Param));
  complex*     Wtab;
  int          N, Npts;
  FILE*        fp;
  const char*  session;
  int          i, semtex = 0, transform = 1;

  /* -- Process command-line arguments. */

  getargs (argc, argv, &session, &semtex, &transform);

  fp = efopen (session, "r");
  readParam   (fp, Info);

  /* -- Set up the problem size. */

  N        = Info -> ngrid;
  K        = N / 2;
  FourKon3 = (4 * K) / 3;

  /* -- Allocate storage for the solution. */

  Wtab = cvector (0, K-1);
  cfield (&U);
  cfield (&Q);
  cbox   (0, N-1, 0, N-1, 0, K-1, &Work);

  /* -- Input Fourier-transformed velocity field. */

  readCVF (fp, U);
  fclose  (fp);

  /* -- Compute vorticity. */

  curl (U, Q, Work);

  /* -- Transform to physical space. */

  if (transform) {
    preFFT (Wtab, K);
    for (i = 1; i <= 3; i++) {
      rc3DFT (U[i], Wtab, INVERSE);
      rc3DFT (Q[i], Wtab, INVERSE);
    }
    enstrophy (Work, Q);
  }

  /* -- Create flagged output. */

  if   (semtex) semout (session, U, Q, Work, transform);
  else          tecout (session, U, Q, Work, transform);

  return EXIT_SUCCESS;
}


static void getargs (int          argc   ,
		     char**       argv   ,
		     const char** session,
		     int*         semtex ,
		     int*         trans  )
/* ------------------------------------------------------------------------- *
 * Process command-line arguments.
 * ------------------------------------------------------------------------- */
{
  const char usage[] =
    "Usage: %s [options] input.fld\n"
    "options:\n"
    "-f ... suppress inverse Fourier transformation\n"
    "-h ... print this message\n"
    "-s ... create semtex/prism format output [Default: TECPLOT]\n";
  char c;
 
  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'f':
      *trans  = FALSE;
      break;
    case 'h':
      fprintf (stdout, usage, prog);
      exit (EXIT_SUCCESS);
      break;
    case 's':
      *semtex = TRUE;
      break;
    default:
      break;
    }

  if   (argc != 1)  message (prog, "no session definition file", ERROR);
  else             *session = *argv;
}
 

static void enstrophy (CF         Q,
		       const CVF  U)
/* ------------------------------------------------------------------------- *
 * Make enstrophy out of physical-space vorticity field.
 * ------------------------------------------------------------------------- */
{
  const int     ntot = N * N * N;
  register real *u = &U[1][0][0][0].Re;
  register real *v = &U[2][0][0][0].Re;
  register real *w = &U[3][0][0][0].Re;
  register real *q = &Q   [0][0][0].Re;
  register int  i;

  for (i = 0; i < ntot; i++) q[i] = u[i]*u[i] + v[i]*v[i] + w[i]*w[i];
}


static void tecout (const char* session,
		    const CVF   U      ,
		    const CVF   Q      ,
		    const CF    S      ,
		    const int   trans  )
/* ------------------------------------------------------------------------- *
 * Amtec TECPLOT format ASCII output to stdout.
 * ------------------------------------------------------------------------- */
{
  const real   dxyz = 2.0 * M_PI / N;
  real         x, y, z;
  real         *u, *v, *w, *r, *s, *t, *q;
  register int i, j, k;
  
  u = &U[1][0][0][0].Re;
  v = &U[2][0][0][0].Re;
  w = &U[3][0][0][0].Re;

  r = &Q[1][0][0][0].Re;
  s = &Q[2][0][0][0].Re;
  t = &Q[3][0][0][0].Re;
  
  q = &S   [0][0][0].Re;

  /* -- TECPLOT header. */

  printf ("TITLE = ISO FIELD FILE\n");
  if (trans) {
    printf ("VARIABLES = x y z u v w r s t q\n");
    printf ("ZONE T = \"BOX\", I=%d, J=%d, K=%d\n", N, N, N);
    for (k = 0; k < N; k++) {
      const real z = k * dxyz;
      for (j = 0; j < N; j++) {
	const real y = j * dxyz;
	for (i = 0; i < N; i++) {
	  const real x = i * dxyz;

	  printf ("%g %g %g ", x, y, z);
	  printf ("%g %g %g %g %g %g %g\n",
		  u[ k + N * (j + i * N) ],
		  v[ k + N * (j + i * N) ],
		  w[ k + N * (j + i * N) ],
		  r[ k + N * (j + i * N) ],
		  s[ k + N * (j + i * N) ],
		  t[ k + N * (j + i * N) ],
		  q[ k + N * (j + i * N) ]);
	}
      }
    }
  } else {
    printf ("VARIABLES = x y z u v w r s t\n");
    printf ("ZONE T = \"BOX\", I=%d, J=%d, K=%d\n", N, N, N);
    printf ("VARIABLES = x y z u v w r s t q\n");
    printf ("ZONE T = \"BOX\", I=%d, J=%d, K=%d\n", N, N, N);
    for (k = 0; k < N; k++) {
      const real z = k * dxyz;
      for (j = 0; j < N; j++) {
	const real y = j * dxyz;
	for (i = 0; i < N; i++) {
	  const real x = i * dxyz;

	  printf ("%g %g %g ", x, y, z);
	  printf ("%g %g %g %g %g %g\n",
		  u[ k + N * (j + i * N) ],
		  v[ k + N * (j + i * N) ],
		  w[ k + N * (j + i * N) ],
		  r[ k + N * (j + i * N) ],
		  s[ k + N * (j + i * N) ],
		  t[ k + N * (j + i * N) ]);
	}
      }
    }
  }
}


static void semout (const char* session,
		    const CVF   U      ,
		    const CVF   Q      ,
		    const CF    S      ,
		    const int   trans  )
/* ------------------------------------------------------------------------- *
 * Semtex/prism: deal with this as a single element.  Output a mesh
 * file, tmp.msh, then output the field file on stdout.
 * ------------------------------------------------------------------------- */
{
  const char *hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields written\n",
    "%-25s "    "Format\n"
  };

  const int    ntot = N * N * N;
  const real   dxyz = 2.0 * M_PI / N;
  real         x, y, z;
  real         *u, *v, *w, *r, *s, *t, *q;
  char         buf[STR_MAX];
  FILE*        mesh;
  register int i, j, k;

  /* -- Output mesh file. */

  mesh = efopen ("tmp.msh", "w");

  fprintf (mesh, "%1d %1d %1d 1 NR NS NZ NEL\n", N, N, N);
  for (j = 0; j < N; j++) {
    const real y = j * dxyz;
    for (i = 0; i < N; i++) {
      const real x = i * dxyz;
      fprintf (mesh, "%g %g\n", x, y);
    }
  }
  for (k = 0; k <= N; k++) {
    const real z = k * dxyz;
    fprintf (mesh, "%g\n", z);
  }
  fclose (mesh);

  /* -- Output header. */

  printf  (hdr_fmt[0], session);
  printf  (hdr_fmt[1], "Y2K bug");
  sprintf (buf, "%1d %1d %1d 1", N, N, N);
  printf  (hdr_fmt[2], buf);
  printf  (hdr_fmt[3], 0);
  printf  (hdr_fmt[4], 0.0);
  printf  (hdr_fmt[5], 0.0);
  printf  (hdr_fmt[6], 0.0);
  printf  (hdr_fmt[7], 1.0);
  if   (trans) printf  (hdr_fmt[8], "uvwrstq");
  else         printf  (hdr_fmt[8], "uvwrst");
  sprintf (buf, "binary ");
  format  (buf + strlen (buf));
  printf  (hdr_fmt[9], buf);

  /* -- Output each field, plane-by-plane...the order it's in anyway. */

  u = &U[1][0][0][0].Re;
  v = &U[2][0][0][0].Re;
  w = &U[3][0][0][0].Re;

  r = &Q[1][0][0][0].Re;
  s = &Q[2][0][0][0].Re;
  t = &Q[3][0][0][0].Re;
  
  q = &S   [0][0][0].Re;

  fwrite (u, sizeof (real), ntot, stdout);
  fwrite (v, sizeof (real), ntot, stdout);
  fwrite (w, sizeof (real), ntot, stdout);
  fwrite (r, sizeof (real), ntot, stdout);
  fwrite (s, sizeof (real), ntot, stdout);
  fwrite (t, sizeof (real), ntot, stdout);
  if (trans) fwrite (q, sizeof (real), ntot, stdout);
}
