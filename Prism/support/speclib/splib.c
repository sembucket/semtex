/*

  LIBRARY ROUTINES FOR SPECTRAL METHODS 

  C Version (double precision)

  Fortran version by:  Einar Malvin Ronquist
                       Room 3-243
                       Department of Mechanical Engineering
		       Massachusetts Institute of Technology
		       77 Massachusette Avenue
		       Cambridge, MA 02139

  C translation:       Ron Henderson
                       Department of Mechanical and Aerospace Engineering
                       Engineering Quadrangle D-209B
		       Princeton University
		       Princeton, NJ 08541

  There are additional routines which did not appear in the original
  FORTRAN library.


  Abbreviations

  M      -    Set of mesh points
  Z      -    Set of collocation/quadrature points
  W      -    Set of quadrature weights
  H      -    Lagrangian interpolant
  D      -    Derivative operator
  I      -    Interpolation operator
  GL     -    Gauss Legendre
  GLL    -    Gauss-Lobatto Legendre
  GJ     -    Gauss Jacobi
  GLJ    -    Gauss-Lobatto Jacobi
  GRL    -    Gauss-Radau Legendre

  -----------------------------------------------------------------------
                         M A I N     R O U T I N E S
  -----------------------------------------------------------------------

  Points and Weights:

  zwgl        Compute Gauss Legendre points and weights
  zwgll       Compute Gauss-Lobatto Legendre points and weights
  zwgrl       Compute Gauss-Radau Legendre points and weights

  Lagrangian Interpolants:

  hgll        Compute Gauss-Lobatto Legendre Lagrangian interpolants

  Derivative Operators:

  dgll        Compute Gauss-Lobatto Legendre derivative matrix


  Interpolation Operators:

  iglm        Compute interpolation operator GL->M
  igllm       Compute interpolation operator GLL->M

  Polynomial Evaluation:

  pnleg       Compute Legendre polynomial of degree N
  pndleg      Compute derivative of Legendre polynomial of degree N
  pnd2leg     Compute 2nd derivative of Legendre polynomial of degree N


  ------------------------------------------------------------------------


  Useful references:

  [1] Gabor Szego: Orthogonal Polynomials, American Mathematical Society,
      Providence, Rhode Island, 1939.
  [2] Abramowitz & Stegun: Handbook of Mathematical Functions,
      Dover, New York, 1972.
  [3] Canuto, Hussaini, Quarteroni & Zang: Spectral Methods in Fluid
      Dynamics, Springer-Verlag, 1988.



  NOTES
  -----
  (1) All routines are double precision.  
  (2) All array subscripts start from zero, i.e. vector[0..N-1] 
  (3) Matrices should be allocated as true 2-dimensional arrays with
      row and column indices starting from 0.


  RCS Information
  -------------------------
  $Author$
  $Source$
  $Date$
  $Revision$
  -------------------------
*/

#include <sys/types.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "speclib/speclib.h"
#include "speclib/splib.h"

#define STOP  100
#define EPS   DBL_EPSILON

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

/* External variables */

static double splib_alp;
static double splib_bet;

/*-----------------------------------------------------------------------*
 * zwgll() - Gauss-Lobatto Legendre points and weights                   *
 *                                                                       *
 * Generate NP Gauss-Lobatto Legendre points (z) and weights (w) as-     *
 * sociated with the Jacobi polynomial P(n) (alpha = 0, beta = 0).  The  *
 * polynomial degree n = np - 1                                          *
 *                                                                       *
 *-----------------------------------------------------------------------*/


void zwgll (double *z, double *w, int np)
{
  if( np == 1 ) {
    z[0] = w[0] = 0.0;
  }
  else {
    double *zd = (double*) malloc( np * sizeof(double) );
    double poly,pder,polym1,pderm1,polym2,pderm2;
    double alpha,beta,con;
    int i, n;

    alpha = beta = 0;
    n     = np - 1;
    con   = 2. / ( n * (n + 1) );
    jacgl (n,alpha,beta,zd);
    
    for(i = 0;i < np;i++) {
      jacobf(n,&poly,&pder,&polym1,&pderm1,&polym2,&pderm2,zd[i]);
      z[i] = zd[i];
      w[i] = con / ( poly * poly );
    }
    free (zd);
  }
  return;
}


/*-----------------------------------------------------------------------*
 * jacgl() - Gauss-Lobatto points and weights                            *
 *                                                                       *
 * Compute Gauss-Lobatto points xjac for a polynomial of degree n. Alpha *
 * and beta determine the specific type of Gauss-Lobatto points.         *
 *                                                                       *
 * For example...     alpha = beta = 0        -> Legendre points         *
 *                    alpha = beta = -1/2     -> Chebychev points        *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void jacgl (int n, double alpha, double beta, double *xjac)
{
  int    np,nh,i,j,k;
  double pnp1p,pdn1p,pnp,pdnp,pnm1p,pdnm1p;
  double pnp1m,pdn1m,pnm,pdnm,pnm1m,pdnm1m;
  double det,rp,rm,a,b,con;

  splib_alp = alpha;
  splib_bet = beta;
  np  = n + 1;
  
  jacobf(np,&pnp1p,&pdn1p,&pnp,&pdnp,&pnm1p,&pdnm1p, 1.0);
  jacobf(np,&pnp1m,&pdn1m,&pnm,&pdnm,&pnm1m,&pdnm1m,-1.0);

  det = pnp * pnm1m - pnm * pnm1p;
  rp  = -pnp1p;
  rm  = -pnp1m;
  a   = (rp * pnm1m - rm * pnm1p) / det;
  b   = (rm * pnp - rp * pnm) / det;

  xjac[np - 1] = 1.0;
  nh  = (n + 1) >> 1;
  con = M_PI / (double) n;

  for (j = 1; j < nh; j++)
    {
      double x = cos (con*j), dx = 1.;
      
      for (k = 0; k < STOP && fabs(dx) > EPS; k++)
	{
	  double poly,pder,recsum;
	  double pnp1,pdnp1,pn,pdn,pnm1,pdnm1;

	  jacobf(np,&pnp1,&pdnp1,&pn,&pdn,&pnm1,&pdnm1,x);
	  poly   = pnp1 + a * pn + b * pnm1;
	  pder   = pdnp1 + a * pdn + b * pdnm1;
	  recsum = 0.;

	  for (i = 0; i < j-1; i++) 
	    recsum += 1. / (x - xjac [np-i-1]);
	  
	  x += (dx = -poly / (pder - recsum * poly));
	}

      if (k == STOP)
	speclib_warning
	  ("splib: jacgl did NOT converge, resid = %g\n", fabs(dx));

      xjac[np - j - 1] = x;
    }

  xjac[0] = -1.0;
  for(i = 1;i < nh;i++) xjac[i] = -xjac[np - i - 1];

  if (np & 1) xjac[ nh ] = 0.0;
  return;
}

/*-----------------------------------------------------------------------*
 * jacobf() - Jacobi Polynomial                                          *
 *                                                                       *
 * Computes the Jacobi polynomial (POLY) and its derivative (PDER) of    *
 * degree N at the point X.                                              *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void jacobf (int n, double *poly, double *pder, double *polym1, 
             double *pderm1, double *polym2, double *pderm2, double x)
{
  *poly = 1.0;
  *pder = 0.0;

  if (n > 0) {
    double rv    = splib_alp + 1.;
    double apb   = splib_alp + splib_bet;
    double polyl = *poly;
    double pderl = *pder;
    
    *poly = rv * x;
    *pder = rv;
    
    if (n == 1) 
      return;
    else {
      int k;
      double a0,a1,a2,a3,a4,a5;
      double polyn,pdern,psave,pdsave;
      double one = 1.0, two = 2.0;
      
      for (k = 2;k <= n;k++) {
	a0 = two * k + apb - two;
	a1 = two * k * (k + apb) * a0;
	a2 = (a0 + one) * (splib_alp * splib_alp - splib_bet 
			   * splib_bet);
	a3 = a0 * (a0 + one) * (a0 + two);
	a4 = two * (k + splib_alp - one) * (k + splib_bet - one) 
	  * (two * k + apb);
	a5 = a2 + a3 * x;
	
	polyn  = (a5 * (*poly) - a4 * polyl) / a1;
	pdern  = (a5 * (*pder) - a4 * pderl + a3 * (*poly) ) / a1;
	psave  = polyl;
	pdsave = pderl;
	polyl  = *poly;
	*poly  = polyn;
	pderl  = *pder;
	*pder  = pdern;
      }
      *polym1 = polyl;
      *pderm1 = pderl;
      *polym2 = psave;
      *pderm2 = pdsave;
    }
  }
  
  return;
}

/*-----------------------------------------------------------------------*
 * zwgrl() - Gauss-Radau Legendre Points and Weights                     *
 *                                                                       *
 * Generate NP Gauss-Radau Legendre points (Z) and weights (W) associat- *
 * ed with the Jacobi polynomial P(N) (alpha = 0,beta = 0).  The poly-   *
 * nomial degree N = NP-1                                                *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void zwgrl (double *z, double *w, int np)
{
  if (np == 1)
    {
      z[0] = 0.;
      w[0] = 2.;
    }
  else
    {
      int n = np - 1;
      double *zd = (double*) malloc(np * sizeof(double));
      double poly,pder,polym1,pderm1,polym2,pderm2;
      double alpha,beta,con;
      int i;

      alpha = beta = 0.;
      con   = 1. / ( (n + 1) * (n + 1) );

      jacgr (n,alpha,beta,zd);
      jacobf(n,&poly,&pder,&polym1,&pderm1,&polym2,&pderm2,zd[0]);

      z[0]  = zd[0];
      w[0]  = 2. * con;
      for (i = 1;i < np;i++)
	{
	  jacobf(n,&poly,&pder,&polym1,&pderm1,&polym2,&pderm2,zd[i]);
	  z[i] = zd[i]; 
	  w[i] = con * ( 1. - zd[i] ) / ( poly * poly );
	}
      free (zd);
    }
  return;
}

/*-----------------------------------------------------------------------*
 * jacgr() - Gauss-Radau points                                          *
 *                                                                       *
 * Compute Gauss-Radau points XJAC for a polynomial of degree N.  ALPHA  *
 * and BETA determine the specific type of Gauss-Radau points.  For ex-  *
 * ample...                                                              *
 *                  alpha = beta =  0.0  -> Legendre points              *
 *                  alpha = beta = -0.5  -> Chebchev points              *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void jacgr (int n, double alpha, double beta, double *xjac)
{
  int np = n + 1;
  double con = 2.0 * M_PI / ((n<<1) + 1);
  double pnp1,pdnp1,pn,pdn,pnm1,pdnm1;
  double func,funcd;
  int i,j,k;

  splib_alp = alpha;
  splib_bet = beta;

  for(j = 0;j < np;j++)
    {
      double x = -cos (con*j), dx = 1.;

      for (k = 1; k < STOP && fabs(dx) > EPS; k++)
	{
	  double recsum = 0.;

	  jacobf(np,&pnp1,&pdnp1,&pn,&pdn,&pnm1,&pdnm1,x);
	  func  = pn + pnp1;
	  funcd = pdn + pdnp1;

	  for(i = 0;i < j-1;i++) recsum += 1.0 / (x - xjac[i]);

	  x += (dx = - func / (funcd - recsum * func));
	}

      if (k == STOP) 
	speclib_warning
	  ("splib: jacgr did NOT converge, resid = %g\n", fabs(dx));

      xjac[j] = x;
    }
}

/*-----------------------------------------------------------------------*
 * zwgl() - Gauss Legendre Points and Weights                            *
 *                                                                       *
 * Generate NP Gauss Legendre points (Z) and weights (W) associated with *
 * the Jacobi polynomial P(N) (alpha = 0,beta = 0).  The polynomial deg- *
 * ree N = NP-1                                                          *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void zwgl (double *z, double *w, int np)
{
  if( np == 1 )
    {
      z[0] = 0.0;
      w[0] = 2.0;
    }
  else
    {
      int n = np - 1;
      double *zd = (double *) malloc((unsigned) np * sizeof(double));
      double poly,pder,polym1,pderm1,polym2,pderm2;
      double alpha,beta;
      double one = 1.0, two = 2.0;
      int i;

      alpha = beta = 0.0;
      jacg(n,alpha,beta,zd);

      for(i = 0;i < np;i++)
	{
	  jacobf(n,&poly,&pder,&polym1,&pderm1,&polym2,&pderm2,zd[i]);
	  z[i] = zd[i];
	  w[i] = two / ( ( one - zd[i] * zd[i]) * pder * pder );
	}
      free( zd );
    }
  return;
}

/*-----------------------------------------------------------------------*
 * jacg() - Gauss points                                                 *
 *                                                                       *
 * Compute Gauss points XJAC for a polynomial of degree N.  ALPHA and    *
 * BETA determine the specific type of Gauss points.  For example,       *
 *                                                                       *
 *                 alpha = beta =  0.0  -> Legendre points               *
 *                 alpha = beta = -0.5  -> Chebchev points               *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void jacg (int n, double alpha, double beta, double *xjac)
{
  int np = n + 1;
  int nh = (n + 1) >> 1;
  double dth = M_PI / ((n<<1) + 2 );
  double poly,pder,polym1,pderm1,polym2,pderm2;
  int i,j,k;

  splib_alp = alpha;
  splib_bet = beta;

  for (j = 0;j < nh;j++)
    {
      double recsum = 0.;
      double x = cos(((j<<1) + 1) * dth), dx = 1.;
      
      for (k = 1; k < STOP && fabs(dx) > EPS; k++)
	{
	  jacobf(np,&poly,&pder,&polym1,&pderm1,&polym2,&pderm2,x);
	  
	  for(i = 0;i < j-1;i++) 
	    recsum += 1. / (x - xjac[np - i - 1]);
	  
	  x += (dx = -poly / (pder - recsum * poly));
	}

      if (k == STOP) 
	speclib_warning
	  ("splib: jacg did NOT converge, resid = %g\n", fabs(dx));
      
      xjac [np-j-1] = x;
    }

  for(i = 0;i < nh;i++) xjac[i] = -xjac[np - i - 1];

  if( np & 1 ) xjac[ nh ] = 0.0;
}

/*-----------------------------------------------------------------------*
 * dgll() - Compute the Derivative Matrix                                *
 *                                                                       *
 * Compute the derivative matrix D and its transpose DT associated with  *
 * the Nth order Lagrangian interpolants through the NZ Gauss-Lobatto    *
 * Legendre points Z.                                                    *
 *                                                                       *
 *              dU                                                       *
 *              --  = Dij * Uj evaluated at z = zi                       *
 *              dz                                                       *
 *                                                                       *
 * NOTE: D and DT are both square matrices.                              *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void dgll (double **d, double **dt, double *z, int nz)
{
  if( nz == 1 )
    {
      d[0][0] = dt[0][0] = 0.0;
    }
  else
    {
      int n = nz - 1;
      double d0 = n * (n + 1) * 0.25;
      int i,j;

      d[0][0] = -d0;
      for(i = 1;i < n;i++) d[i][i] = 0.0;
      d[n][n] =  d0;

      for(i = 0;i < nz;i++) 
	for(j = 0;j < nz;j++) {
	  if( i != j ) d[i][j] = pnleg(z[i],n) / 
                               ( pnleg(z[j],n) * (z[i]-z[j]) );
	  dt[j][i] = d[i][j];
	}
    }
  return;
}

/*-----------------------------------------------------------------------*
 * hgl() - Compute the GL Lagrangian Interpolant                         *
 *                                                                       *
 * Compute the value of the Lagrangian interpolant HGL through the NZ    *
 * Gauss-Lobatto points ZGL at the point Z.                              *
 *                                                                       *
 *-----------------------------------------------------------------------*/
 
double hgl (int i, double z, double *zgl, int nz)
{
  double dz = z - zgl[i];

  if (fabs(dz) < EPS) return 1.;

  return pnleg(z,nz) / (pndleg(zgl[i],nz) * dz);
}

/*-----------------------------------------------------------------------*
 * hgll() - Compute the GLL Lagrangian Interpolant                       * 
 *                                                                       *
 * Compute the value of the Lagrangian interpolant HGLL through the NZ   *
 * Gauss-Lobatto Legendre points ZGLL at the point Z.                    *
 *                                                                       *
 *-----------------------------------------------------------------------*/

double hgll (int i, double z, double *zgll, int nz)
{
  int    n     = nz - 1;
  double dz    = z - zgll[i];
  double alpha = n * (n + 1);

  if (fabs(dz) < EPS) return 1.;

  return  -(1. - z*z) * pndleg(z,n) / (alpha * pnleg(zgll[i],n) * dz);
}

/*-----------------------------------------------------------------------*
 * pnleg() - Legendre Polynomial Evaluation                              *
 *                                                                       *
 * Compute the value of the Nth order Legendre polynomial at z.  Based   *
 * on the recursion formula for the Legendre polynomial.                 *
 *-----------------------------------------------------------------------*/

double pnleg (double z, int n)
{
  register double p0 = 1., 
                  p1 = z;
  register int    k  = 2;
  register double p2, p, dk;

  switch (n) {
  case 0:
    p = p0;
    break;
  case 1:
    p = p1;
    break;    
  default:
    do {  
      dk = (double) k;
      p2 = ((2.*dk - 1.) * z * p1 - (dk-1.) * p0) / dk;
      p0 = p1;
      p1 = p2;  
    } while (++k <= n);  
    p = p2;
    break;
  }

  return p;
}


/*-----------------------------------------------------------------------*
 * pndleg() - Legendre Polynomial Derivative Evaluation                  *
 *                                                                       *
 * Compute the derivative of the Nth Order Legendre polynomial at Z.     *
 * Based on the recursion formula for the Legendre polynomial.           * 
 *                                                                       *
 *-----------------------------------------------------------------------*/

double pndleg (double z, int n)
{
  register double p0  = 1.,
                  p0d = 0.,
                  p1  = z,
                  p1d = 1.;
  register int    k   = 2;
  register double p2, p2d, p, dk, fac;

  switch (n) {
  case 0:
    p = p0d;
    break;
  case 1:
    p = p1d;
    break;
  default:
    do {
      dk  = (double) k;
      fac = 2.* dk - 1.;
      p2  = (fac * z * p1 - (dk-1.) * p0) / dk;
      p2d = (fac * (p1 + z * p1d) - (dk-1.) * p0d) / dk;
      
      p0  = p1;     p1  = p2;
      p0d = p1d;    p1d = p2d;

    } while (++k <= n);  
    p = p2d;
  }
  return p;
}


/*-----------------------------------------------------------------------*
 * pnd2leg() - Legendre Polynomial 2nd Derivative Evaluation             *
 *                                                                       *
 * Compute the 2nd derivative of the Nth-order Legendre Polynomial at z. *
 * Based on Legendre's equation,                                         *
 *                                                                       *
 *                    2                                                  *
 *              (1 - z ) L'' - 2 z L' + n (n+1) L = 0                    *
 *                                                                       *
 *-----------------------------------------------------------------------*/

double pnd2leg (double z, int n)
{
  double d = 1.0 - z*z;
  return (2. * z * pndleg(z,n) - n*(n+1) * pnleg(z,n)) / d;
}

/*-----------------------------------------------------------------------*
 * igllm() - Interpolation Operator GLL -> M                             *
 *                                                                       *
 * Compute the one-dimensional interpolation operator (matrix) I12 and   *
 * its transpose IT12 for interpolating a function from a Gauss-Labatto  *
 * Legendre mesh (1) to another mesh M (2).                              *
 *                                                                       *
 * If the argument "im12t" is passed as NULL, it will not be computed.   *
 *                                                                       *
 *-----------------------------------------------------------------------*/

void igllm(double **im12, double **im12t, 
	   double *zgll, double *zm, int nz, int mz)
{
  double zp;
  register int i, j;

  if (nz == 1)
    **im12 = 1.;
  else
    for (i = 0; i < mz; i++) {
      zp = zm[i];
      for (j = 0; j < nz; j++)
	im12 [i][j] = hgll(j, zp, zgll, nz);
    }

  if (im12t)
    for (i = 0; i < nz; i++)
      for (j = 0; j < mz; j++)
	im12t[i][j] = im12[j][i];
  
  return;
}
