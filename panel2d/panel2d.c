/*****************************************************************************
 * PANEL2D: compute panel-method solution to an arbitrary lifting airfoil.
 *
 * USAGE: panel2d -i incidenceangle panelfile
 *
 * SYNOPSIS
 * --------
 * Panel2d computes inviscid, irrotational flows past 2D bodies of
 * arbitrary shape, set at an arbitrary angle of incidence to the
 * oncoming flow.  A LINEAR VORTICITY distribution is used on each
 * element, and a Kutta condition is enforced on the start of the
 * first panel (typically the trailing edge point).  The body is made
 * up of flat panels.  A Neumann or normal-velocity boundary condition
 * is used with collocation points at panel midpoints.  The arbitrary
 * angle of incidence is set on the command line (deg.).
 *
 * FILES
 * -----
 * panelfile contains the body information: the first line is treated
 * as a comment, the second line gives the number of panels, the
 * remaining lines supply the x--z coordinates at the start of each
 * panel.  The end of the last panel is assumed to coincide with the
 * start of the first, with this point also being the point at which
 * the Kutta condition is met.  NB: Shape is assumed to be traversed
 * CW.
 *
 * OUTPUT
 * ------
 * The pressure coefficient at each collocation point is printed on
 * standard output.
 *
 * REFERENCES
 * ----------
 * 1. Katz & Plotkin 1991, Low-Speed Aerodynamics. McGraw-Hill.
 * 2. Keuthe & Chow 1986, Foundations of Aerodynamics. Wiley.
 *
 * VERIFICATION
 * ------------
 * Code gives the same values as Keuthe & Chow, NACA 2412 at 8 deg, p 135:
 *    panel2d -i 8 NACA2412
 *
 * AUTHOR
 * ------
 * Hugh Blackburn
 * Department of Mechanical Engineering
 * Monash University
 * Vic 3800
 * Australia
 * hugh.blackburn@eng.monash.edu.au
 *****************************************************************************/

static char RCSid[] = "$Id$";


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <math.h>
#include "panel.h"


static void    getargs  (int     , char *[], char ** , double *);
static Point  *getgeom  (FILE *  , int *                       );
static Point  *midpnts  (Point * , int                         );
static Point  *norm2D   (Point * , int                         );
static Point  *tang2D   (Point * , int                         );
static double  dot2D    (Point   , Point                       );
static Point   sum2D    (Point   , Point                       ); 
static double *press2D  (int     , double *, Point   , Point  *);
static void    load2D   (int     , double *, double  , Point  *,
			 Point  *, double *, double *, double *);
static void    printup  (FILE *  , char   *, int     , Point  *,
			 Point  *, double *, double *, double *,
			 double  , double  , double            );

static char prog[] = "panel2d";





int main(int argc, char *argv[])
/* ========================================================================= *
 * Driver.                                                                   *
 * ========================================================================= */
{
  FILE     *fp;
  char     *session;

  double    incidence = 0.0;
  int       i, j, npanel;
  Point    *panel, *colloc, *normal, *tangent;
  Point     u_a, u_b, v_a, v_b, U_inf;
  double   *Cp, Cl, Cm, CentP;

  double  **ICn, **ICt, *vel, *gamma, det_sign;
  int      *permut;


  /* ----------------------------------------------------------------------- *
   * Deal with command-line arguments.                                       *
   * ----------------------------------------------------------------------- */

  getargs(argc, argv, &session, &incidence);
  U_inf.x = cos(incidence*PI_180);
  U_inf.z = sin(incidence*PI_180);

  /* ----------------------------------------------------------------------- *
   * Read panel information from file, set normals & collocation points.     *
   * ----------------------------------------------------------------------- */

  fp      = efopen  (session, "r");
  panel   = getgeom (fp,   &npanel);
  colloc  = midpnts (panel, npanel);
  normal  = norm2D  (panel, npanel);
  tangent = tang2D  (panel, npanel);

  /* ----------------------------------------------------------------------- *
   * Allocate influence coefficient matrix and BCs (RHS vector).             *
   * ----------------------------------------------------------------------- */

  ICn    = dmatrix(1, npanel+1, 1, npanel+1);
  ICt    = dmatrix(1, npanel,   1, npanel+1);
  gamma  = dvector(1, npanel+1);
  vel    = dvector(1, npanel  );
  permut = ivector(1, npanel+1);

  /* ----------------------------------------------------------------------- *
   * Fill out influence coefficient matrices & RHS.  A row for each colloc-  *
   * ation point.  K&P eqs (11.102--104)].                                   *
   * ----------------------------------------------------------------------- */

  for (i=1; i<=npanel; i++) {

    vor2Dl(colloc[i], panel[1], panel[2], normal[1], &u_a, &u_b, i == 1);
    ICn[i][1] = dot2D(u_a, normal [i]);
    ICt[i][1] = dot2D(u_a, tangent[i]);

    for (j=2; j<=npanel; j++) {

      vor2Dl(colloc[i], panel[j], panel[j+1], normal[j], &v_a, &v_b, i == j);
      u_a       = sum2D(u_b, v_a);
      ICn[i][j] = dot2D(u_a, normal [i]);
      ICt[i][j] = dot2D(u_a, tangent[i]);
      u_b       = v_b;
    }

    ICn  [i][npanel+1] =  dot2D(u_b,   normal [i]);
    ICt  [i][npanel+1] =  dot2D(u_b,   tangent[i]);
    gamma[i]           = -dot2D(U_inf, normal [i]);
  }

  /* ----------------------------------------------------------------------- *
   * Set up equation to get Kutta condition [K&P eq (11.105)].               *
   * ----------------------------------------------------------------------- */

  dzero(npanel-1, ICn[npanel+1]+2, 1); 
  ICn  [npanel+1][       1] = 1.0;
  ICn  [npanel+1][npanel+1] = 1.0;
  gamma[npanel+1]           = 0.0;

  /* ----------------------------------------------------------------------- *
   * Solve matrix system.                                                    *
   * ----------------------------------------------------------------------- */

  ludcmp(ICn, npanel+1, permut, &det_sign);
  lubksb(ICn, npanel+1, permut, gamma);

  /* ----------------------------------------------------------------------- *
   * Compute collocation point panel-induced tangent velocities.             *
   * ----------------------------------------------------------------------- */

  dmxv(ICt[1]+1, npanel, gamma+1, npanel+1, vel+1);

  /* ----------------------------------------------------------------------- *
   * Compute pressure coefficients and loads.                                *
   * ----------------------------------------------------------------------- */

  Cp = press2D(npanel, vel, U_inf, tangent);
  load2D(npanel, gamma, incidence, panel, colloc, &Cl, &Cm, &CentP);

  /* ----------------------------------------------------------------------- *
   * Print up & exit.                                                        *
   * ----------------------------------------------------------------------- */

  printup(stdout, session, npanel, panel, colloc, vel, gamma,
	  Cp, incidence, Cl, CentP);


  exit(0);
}




static void getargs(int argc, char *argv[], char **session, double *incidence)
/* ========================================================================= *
 * Process command-line arguments.                                           *
 * ========================================================================= */
{
  char c;
  char usage[] = "usage: panel2d [options] panelfile\n"
                 "  options:\n"
		 "  -h            ... print this message\n"
		 "  -i incidence  ... incidence angle (deg.)\n";


  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      fprintf(stderr, usage);
      exit(0);
      break;
    case 'i':
      if (*++argv[0])
	*incidence = atof(*argv);
      else {
	--argc;
	*incidence = atof(*++argv);
      }
      break;
    default:
      sprintf(buf, "illegal option: %c", c);
      message(prog, buf, WARNING);
      break;
    }

  if (argc == 1)
    *session = *argv;
  else {
    fprintf(stderr, usage);
    exit(1);
  }
}

static Point *getgeom(FILE *fp, int *npanel)
/* ========================================================================= *
 * Return panel end-points from file.                                        *
 * ========================================================================= */
{
  Point  *p;
  int     i;

  fgets (buf, STR_MAX, fp);
  fgets (buf, STR_MAX, fp);
  sscanf(buf, "%d", npanel);
  
  p = (Point *) malloc((*npanel+1)*sizeof(Point));
  p--;
  for (i=1; i<=*npanel; i++) {
    fgets (buf, STR_MAX, fp);
    sscanf(buf, "%lf %lf", &p[i].x, &p[i].z);
  }
  p[*npanel+1]= p[1];

  return p;
}


static Point *midpnts(Point *panel, int npanel)
/* ========================================================================= *
 * Return collocation points at panel midpoints.                             *
 * ========================================================================= */
{
  int     i;
  Point  *c;

   
  c = (Point *) malloc(npanel*sizeof(Point));
  c--;
  for (i=1; i<=npanel; i++) {
    c[i].x = 0.5*(panel[i + 1].x + panel[i].x);
    c[i].z = 0.5*(panel[i + 1].z + panel[i].z);
  }

  return c;
}


static Point *norm2D(Point *panel, int npanel)
/* ========================================================================= *
 * Return panel unit outward normals.  CW traverse assumed.                  *
 * ========================================================================= */
{
  int     i;
  Point  *n;
  double  dx, dz, l;

   
  n = (Point *) malloc(npanel*sizeof(Point));
  n--;
  for (i=1; i<=npanel; i++) {
    dx     =  panel[i+1].x - panel[i].x;
    dz     =  panel[i+1].z - panel[i].z;
    l      =  hypot(dx, dz);
    n[i].x = -dz / l;
    n[i].z =  dx / l;
  }

  return n;
}


static Point *tang2D(Point *panel, int npanel)
/* ========================================================================= *
 * Return panel unit outward tangents.  CW traverse assumed.                 *
 * ========================================================================= */
{
  int     i;
  Point  *t;
  double  dx, dz, l;

   
  t = (Point *) malloc(npanel*sizeof(Point));
  t--;
  for (i=1; i<=npanel; i++) {
    dx     =  panel[i+1].x - panel[i].x;
    dz     =  panel[i+1].z - panel[i].z;
    l      =  hypot(dx, dz);
    t[i].x =  dx / l;
    t[i].z =  dz / l;
  }

  return t;
}


static double dot2D(Point a, Point b)
/* ========================================================================= *
 * Dot product of two 2D vectors.                                            *
 * ========================================================================= */
{
  return a.x*b.x + a.z*b.z;
}


static Point sum2D(Point a, Point b)
/* ========================================================================= *
 * Sum of two 2D vectors.                                                    *
 * ========================================================================= */
{
  Point sum;


  sum.x = a.x + b.x;
  sum.z = a.z + b.z;
  
  return sum;
}


static double *press2D(int npanel, double *vel, Point U_inf, Point *tangent)
/* ========================================================================= *
 * Compute panel pressure coefficients.                                      *
 * ========================================================================= */
{
  int     j;
  double *cp = dvector(1, npanel);
  double  Qti;

  for (j=1; j<=npanel; j++) {
    Qti     = dot2D(tangent[j], U_inf);
    vel[j] += Qti;
    cp [j]  = 1.0 - SQR(vel[j]);
  }

  return cp;
}


static void load2D(int     npanel , double *gamma,   double inc   ,
		   Point  *panel  , Point  *colloc ,
		   double *Cl     , double *Cm     , double *CentP)
/* ========================================================================= *
 * Compute moment and lift coefficients, centre of pressure.                 *
 * ========================================================================= */
{
  int     j;
  double *panlen;
  double  x0, xC, cosA, dl;

  
  cosA = cos(inc * PI_180);
  x0 = xC = *Cl = *Cm = 0.0;
  panlen = dvector(1, npanel);
 
  for (j=1; j<=npanel; j++) {
    panlen[j] = hypot(panel[j+1].x - panel[j].x, panel[j+1].z - panel[j].z);
    x0 = MIN(x0, panel[j].x);
    xC = MAX(xC, panel[j].x);
  }

  for (j=1; j<=npanel; j++) {
    dl   = 0.5 * (gamma[j+1] + gamma[j]) * panlen[j];
    *Cl += dl;
    *Cm += dl * (colloc[j].x - x0) * cosA;
  }

  *Cl   *= 2.0 / (xC - x0);
  *Cm   *= 2.0 / SQR(xC - x0);
  *CentP = *Cm / *Cl;

  free_dvector(panlen, 1);
}


static void printup(FILE   *fp    , char  *session,
		    int    npanel ,
		    Point  *panel , Point *colloc ,
		    double *vel   , double *gamma ,
		    double *Cp    , double inc    ,
		    double  Cl    , double CentP  )
/* ========================================================================= *
 * Print solution information on fp.                                         *
 * ========================================================================= */
{
  int j;


  fprintf(fp, "# %-20s" 
	  "Panel method solution with linear vorticity distribution\n",
	  session);
  fprintf(fp, "# %-20d"   "Number of panels\n"           , npanel );
  fprintf(fp, "# %-20.5g" "Incidence angle (deg.)\n"     , inc    );
  fprintf(fp, "# %-20.5g" "Coefficient of lift\n"        , Cl     );
  fprintf(fp, "# %-20.5g" "Centre of pressure / chord\n" , CentP  );
  fprintf(fp, "#\n");
  fprintf(fp, "# j    x[j]      z[j]      cx[j]     cz[j]"
	  "    gamma[j]   vel[j]    Cp[j]\n");
  fprintf(fp, "#\n");

  for (j=1; j<=npanel; j++)
    fprintf(fp, "%3d%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
	    j, panel[j].x, panel[j].z, colloc[j].x, colloc[j].z,
	    gamma[j], vel[j], Cp[j]);

  fprintf(fp, "%3d%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
	  j, panel[j].x, panel[j].z, 0.0, 0.0, gamma[j], 0.0, 0.0);
}


void message (const char *routine, const char *text, int level)
/* ------------------------------------------------------------------------- *
 * A simple error handler.
 * ------------------------------------------------------------------------- */
{
  switch (level) {
  case WARNING:
    fprintf (stderr, "WARNING: %s: %s\n", routine, text); 
    break;
  case ERROR:
    fprintf (stderr, "ERROR: %s: %s\n", routine, text); 
    break;
  case REMARK:
    fprintf (stderr, "%s: %s\n", routine, text);
    break;
  default:
    fprintf (stderr, "bad error level in message: %d\n", level);
    exit (EXIT_FAILURE);
    break;
  }

  if (level == ERROR) exit (EXIT_FAILURE);
}


FILE *efopen (const char *file, const char *mode)
/* ------------------------------------------------------------------------- *
 * fopen file, die if can't.
 * ------------------------------------------------------------------------- */
{
  FILE *fp;

  if (fp = fopen (file, mode)) return fp;

  sprintf (buf, "can't open %s mode %s", file, mode);
  message ("efopen", buf, ERROR);
  
  return (FILE*) 0;
}


double* dvector (int nl, int nh)
/* ------------------------------------------------------------------------- *
 * Allocates a double vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  double* v;

  v = (double*) malloc((size_t) ((nh-nl+1)*sizeof(double)));
  if (v) return v-nl;

  message("dvector()", "allocation failure", WARNING);
  return NULL;
}


int* ivector (int nl, int nh)
/* ------------------------------------------------------------------------- *
 * Allocates an int vector with subscript range [nl..nh].
 * ------------------------------------------------------------------------- */
{
  int* v;
  
  v = (int*) malloc((size_t) ((nh-nl+1)*sizeof(int)));
  if (v) return v-nl;

  message("ivector()", "allocation failure", WARNING);
  return NULL;
}


double **dmatrix (int nrl, int nrh,
		  int ncl, int nch)
/* ------------------------------------------------------------------------- *
 * Allocate a double  matrix with subscript ranges [nrl..nrh][ncl..nch].
 * ------------------------------------------------------------------------- */
{
  int i, nrow = nrh-nrl+1, ncol = nch-ncl+1;
  double  **m;

  m = (double**) malloc((size_t) (nrow*sizeof(double*)));
  if (!m) {
    message("dmatrix()", "allocation failure 1", WARNING);
    return NULL;
  }
  m -= nrl;

  m[nrl] = (double*) malloc((size_t) (nrow*ncol*sizeof(double)));
  if (!m[nrl]) {
    message("dmatrix()", "allocation failure 2", WARNING);
    return NULL;
  }
  m[nrl] -= ncl;

  for (i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  return m;
}


void free_dvector (double *v, int nl)
/* ------------------------------------------------------------------------- *
 * Frees a double vector allocated by dvector().
 * ------------------------------------------------------------------------- */
{
  free((void*) (v+nl));
}


void dzero (int n, double* x, int incx)
{
  if (incx == 1)
    memset (x, '\0', n * sizeof (double));

  else {
    register int i;

    x += (incx < 0) ? (-n+1)*incx : 0;

    for (i = 0; i < n; i++) x[i*incx] = 0.0;
  }
}


void dmxv (double* A, int nra, double* B, int nca, double* C)
{
  register double  *a = A,
                   *c = C;
  register double  sum;
  register int     i, j;

  for (i = 0; i < nra; i++) {
    sum  = 0.0;
    for (j = 0; j < nca; j++) sum += (*a++) * B[j]; 
    *c++ = sum;
  }
}
