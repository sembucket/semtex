/*****************************************************************************
 * POLYOPS.C:  Routines for manipulating polynomials.
 *
 * Summary of routines:
 * --------------------
 * dermat_g: Derivative operator for Lagrange interpolant, arbitrary points.
 * dermat_k: Derivative operator for Lagrange interpolant, at nodes (knots).
 * intmat_g: Interpolation operator for Lagrange interpolant, arbitrary.
 * jacobf  : Jacobi polynomial and its derivative.
 * jacg    : Points for Gauss-Jacobi quadrature.
 * jacgr   : Points for Gauss-Radau-Jacobi quadrature.
 * jacgl   : Points for Gauss-Lobatto-Jacobi quadrature.
 * zwgl    : Points and weights for Gauss-Legendre quadrature.
 * zwgrl   : Points and weights for Gauss-Radau-Legendre quadrature.
 * zwgll   : Points and weights for Gauss-Lobatto-Legendre quadrature.
 * pnleg   : Evaluate Legendre polynomial.
 * pndleg  : Evaluate derivative of Legendre polynomial.
 * pnd2leg : Evaluate second derivative of Legendre polynomial.
 * legtr2d : Compute 2D Legendre transform based on GLL grid.
 * dgll    : Derivative operator for Gauss-Lobatto-Legendre interpolant.
 * uniknot : Points uniformly distributed on [-1, 1].
 *
 * Routines that deal specifically with orthogonal polynomials come from
 * a library of spectral routines written in FORTRAN by Einar Ronquist, MIT.
 * Many of the formulae used may be found in Canuto, Hussaini, Quarteroni &
 * Zang, "Spectral Methods in Fluid Dynamics", Springer, 1988.
 *
 * Everything here is double precision.
 *
 * $Id$
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <femdef.h>
#include <alplib.h>

#define STOP 16


void dermat_g (const integer K   ,
	       const double* zero, 
	       const integer I   ,
	       const double* x   ,
	       double**      D   ,
	       double**      DT  )
/* ------------------------------------------------------------------------- *
 * Given a set of zeros for a Lagrange interpolant order K-1 and a set of I
 * points x, return the derivative operator matrix D and its transpose DT.
 *
 * The matrix D is IxK and DT is KxI, both are zero-offset.  D is defined
 * so that
 *                           [D]{y} = {y'},
 * that is, when it premultiplies a vector of points (ordinates) {y},
 * located at (abscissae) {x}, it returns (ordinates) {y'}, the derivative
 * of the Lagrange interpolant with knots {zero}.
 *
 * The set of points {x} may be identical to the knot points {zero} in which
 * case the matrix [D] maps the values {y} onto the derivative of their
 * interpolating polynomial evaluated at the same locations.   In this case,
 * however, the routine dermat_k(), which should have lower rounding errors,
 * may be used instead.
 *
 * Reference: Abramowitz & Stegun 25.3.2.
 * ------------------------------------------------------------------------- */
{
  register integer i, j, k, l;
  register double* a;
  register double  sum, prod;

  a = dvector (0, K-1);
  dfill (K, 1.0, a, 1);

  for (k=0; k<K; k++)
    for (l=0; l<K; l++)
      if (l != k) a[k] *= zero[k] - zero[l];

   for (j=0; j<I; j++)
     for (k=0; k<K; k++) {
       sum = 0.0;
       for (l=0; l<K; l++) {
	 if (l != k) {
	   prod = 1.0;
	   for (i=0; i<K; i++)
	     if (i != k && i != l) prod *= x[j] - zero[i];
	   sum += prod;
	 }
       }
       D[j][k]  = DT[k][j] = sum / a[k];
     }

  freeDvector (a, 0);
}


void dermat_k (const integer K   ,
	       const double* zero,
	       double**      D   ,
	       double**      DT  )
/* ------------------------------------------------------------------------- *
 * Given a set of zeros for a Lagrange interpolant order K-1, return the
 * (collocation) derivative operator matrix D and its transpose DT.
 *
 * The matrix D is KxK and is zero-offset.  D is defined so that
 *                           [D]{y} = {y'},
 * that is, when it premultiplies a vector of points (ordinates) {y},
 * located at (abscissae) {zero}, it returns (ordinates) {y'}, the deriv-
 * ative of the Lagrange interpolant with those knot points (roots).
 *
 * The fact that the evaluation points lie at the zeros of the Lagrange
 * functions simplifies the operations that must be performed when compared
 * to the more general routine dermat_g().
 *
 * It may be more accurate to use closed-form results for any given class of
 * e.g. orthogonal Lagrange interpolants if available.
 *
 * Reference: Solomonoff, A. & Turkel, E.  1989.  "Global properties of
 *   pseudospectral methods".  JCP 81, 239--276
 * ------------------------------------------------------------------------- */
{
  register integer j, k, l;
  register double* a;
  register double  sum, prod;

  a = dvector (0, K-1);
  dfill (K, 1.0, a, 1);

  for (k = 0; k < K; k++)
    for (l = 0; l < K; l++)
      if (l != k) a[k] *= zero[k] - zero[l];

   for (j = 0; j < K; j++)
     for (k = 0; k < K; k++)
       if (j == k) {		    /* -- Diagonal term:     use eq (2.5). */
	 sum = 0.0;
	 for (l = 0; l < K; l++)
	   if (l != k) sum += 1.0 / (zero[k] - zero[l]);
	 D[k][k] = DT[k][k] = sum;
       } else {			    /* -- Off-diagonal term: use eq (2.7). */
	 prod = 1.0;
	 for (l = 0; l < K; l++) 
	   if (l != j) prod *= zero[j] - zero[l];
	 D[j][k] = DT[k][j] = prod / ((zero[j] - zero[k]) * a[k]);
       }

  freeDvector (a, 0);
}


void intmat_g (const integer K   ,
	       const double* zero,
	       const integer I   ,
	       const double* x   ,
	       double**      IN  ,
	       double**      IT  )
/* ------------------------------------------------------------------------- *
 * Given a set of zeros for a Lagrange interpolant order K-1 and a set of I
 * points x, return the interpolation matrix IN and its transpose IT.
 *
 * The matrix IN is IxK and IT is KxI, both are zero-offset.  IN is defined
 * so that
 *                           [IN]{y} = {z},
 * that is, it maps a set of K (ordinates) {y}, given at the Lagrange knot
 * points (abscissae, roots) onto a set of I ordinates {z}, given at arbitr-
 * ary locations {x}
 *
 * Reference: Abramowitz & Stegun 25.2.2.
 * ------------------------------------------------------------------------- */
{
  register integer i, j, k, l;
  register double* a;
  register double  prod;

  a = dvector (0, K-1);
  dfill (K, 1.0, a, 1);

  for (k = 0; k < K; k++)
    for (l = 0; l < K; l++)
      if (l != k) a[k] *= zero[k] - zero[l];

  for (j = 0; j < I; j++)
    for (k = 0; k < K; k++) {
      prod = 1.0;
      for (i = 0; i < K; i++)
	if (i != k) prod *= x[j] - zero[i];
      IN[j][k] = IT[k][j] = prod / a[k];
    }

  freeDvector (a, 0);
}  


static void jacobf (const integer n     ,
		    const double  x     ,
		    const double  alpha ,
		    const double  beta  ,
		    double*       poly  ,
		    double*       pder  ,
		    double*       polym1,
		    double*       pderm1,
		    double*       polym2,
		    double*       pderm2)
/* ------------------------------------------------------------------------- *
 * Computes the Jacobi polynomial (poly) of degree n, and its derivative
 * (pder), at location x.  Values for lower degree are also returned.
 *
 * n          :  degree of approximation,
 * alpha, beta:  parameters in Jacobi weight function
 *                                   alpha           beta
 *                    w(x) = (1 - x)^      * (1 + x)^
 *               where special cases alpha = beta =  0.0 <==> G-L-Legendre
 *                                   alpha = beta = -0.5 <==> G-L-Chebyshev.
 *
 * References:
 *   LIBRARY ROUTINES FOR SPECTRAL METHODS, Einar Ronquist, 1988.
 *   Canuto et al.,  "Spectral Methods in Fluid Dynamics".  1988. App C.
 * ------------------------------------------------------------------------- */
{
  register integer k;
  register double  apb, a1, a2, a3, a4, b3, dk, kk;
  register double  polylst, pderlst, polyn, pdern, psave, pdsave;

  *poly = 1.0;
  *pder = 0.0;

  if (n < 1)  return;

  apb     = alpha + beta;
  polylst = *poly;
  pderlst = *pder;
  *pder   = (1.0 + alpha);
  *poly   = *pder * x;

  if (n == 1) return;

  for (k = 2; k <= n; k++) {
    dk      = (double) k;
    kk      = (double) (k << 1);
    a1      = kk * (dk + apb) * (kk + apb - 2.0);
    a2      = (kk + apb - 1.0) * (SQR(alpha) - SQR(beta));
    b3      = (kk + apb - 2.0);
    a3      = b3 * (b3 + 1.0) * (b3 + 2.0);
    a4      = 2.0 * (dk + alpha - 1.0) * (dk + beta - 1.0) * (kk + apb);
    polyn   = ((a2 + a3*x)*(*poly) - a4*polylst) / a1;
    pdern   = ((a2 + a3*x)*(*pder) - a4*pderlst + a3*(*poly)) / a1;
    psave   = polylst;
    pdsave  = pderlst;
    polylst = *poly;
    *poly   = polyn;
    pderlst = *pder;
    *pder   = pdern;
  }

  *polym1 = polylst;
  *pderm1 = pderlst;
  *polym2 = psave;
  *pderm2 = pdsave;
}


void jacg (const integer n    ,
	   const double  alpha,
	   const double  beta ,
	   double*       xjac )
/* ------------------------------------------------------------------------- *
 * Compute Gauss points xjac for a polynomial of degree n on range [-1, 1].
 * Points are returned in ascending order.  N + 1 points are found.
 *
 * Alpha, beta = 0.0 => Gauss-Legendre; = -0.5 => Chebyshev.
 * ------------------------------------------------------------------------- */
{
  register integer i, j, k, np, nh;
  register double  dth, recsum, x, delx;
  double           poly, pder, polym1, pderm1, polym2, pderm2;

  np  = n + 1;
  nh  = np >> 1;
  dth = M_PI / (np << 1);

  for (j = 0; j < nh; j++) {
    x = cos (((j<<1) + 1) * dth);

    k = 0;
    do {
      jacobf (np, x, alpha, beta, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);

      recsum = 0.0;
      for (i = 0; i < j-1; i++) recsum += 1.0 / (x - xjac[n - i]);

      delx = -poly / (pder - recsum*poly);
      x   += delx;
    } while (fabs (delx) > EPSDP && k++ < STOP);
    xjac [n - j] = x;
  }

  for (i = 0; i < nh; i++) xjac[i] = -xjac[n - i];

  if (np & 1) xjac[nh] = 0.0;
}


void jacgr (const integer n    ,
	    const double  alpha,
	    const double  beta ,
	    double*       xjac )
/* ------------------------------------------------------------------------- *
 * Computes the Gauss-Radau points xjac for polynomial of degree n on
 * [-1, 1].  Points are returned in ascending order.
 *
 * Alpha, beta = 0.0 => Gauss-Legendre; = -0.5 => Chebyshev.
 * ------------------------------------------------------------------------- */
{
  register integer np, i, j, k;
  register double  x, delx, con, recsum;
  double           pn, pdn, pnp1, pdnp1, pnm1, pdnm1, func, funcd;

  np  = n + 1;
  con = 2.0 * M_PI / (n<<1 + 1);

  for (j = 0; j < np; j++) {
    x = -cos (con * j);

    k = 0;
    do {
      jacobf (np, x, alpha, beta, &pn, &pdn, &pnp1, &pdnp1, &pnm1, &pdnm1);
      func  = pn  + pnp1;
      funcd = pdn + pdnp1;

      recsum = 0.0;
      for (i = 0; i < j-1; i++) recsum += 1.0 / (x - xjac[i]);

      delx  = -func  / (funcd - recsum*func);
      x    += delx;
    } while (fabs (delx) > EPSDP && k++ < STOP);

    xjac[j] = x;
  }
}


void jacgl (const integer n    ,
	    const double  alpha,
	    const double  beta ,
	    double*       xjac )
/* ------------------------------------------------------------------------- *
 * Computes the Gauss-Lobatto points xjac for polynomial of degree n on
 * [-1, 1].  Points are returned in ascending order
 *
 * Alpha, beta = 0.0 => Gauss-Legendre, = -0.5 => Chebyshev.
 * ------------------------------------------------------------------------- */
{
  register integer i, j, k, np, jm, nh;
  register double  a, b, det, con, rp, rm, x, delx, recsum;
  double           poly, pder, pnp1p, pdnp1p, pnp, pdnp, pnm1p, pdnm1;
  double           pnp1m, pdnp1m, pnm, pdnm, pnm1m;

  np      = n + 1;
  nh      = np >> 1;
  xjac[n] = 1.0;
  con     = M_PI / (double) n;

  jacobf (np,  1.0, alpha, beta, &pnp1p, &pdnp1p, &pnp, &pdnp, &pnm1p, &pdnm1);
  jacobf (np, -1.0, alpha, beta, &pnp1m, &pdnp1m, &pnm, &pdnm, &pnm1m, &pdnm1);
  det = pnp*pnm1m - pnm*pnm1p;
  rp  = -pnp1p;
  rm  = -pnp1m;
  a   = (rp*pnm1m - rm*pnm1p) / det;
  b   = (rm*pnp   -   rp*pnm) / det;

  for (j = 1; j < nh; j++) {
    jm = j - 1;
    x  = cos (con * j);

    k  = 0;

    do {
      jacobf (np, x, alpha,beta, &pnp1p, &pdnp1p, &pnp, &pdnp, &pnm1p, &pdnm1);
      poly = pnp1p  + a*pnp  + b*pnm1p;
      pder = pdnp1p + a*pdnp + b*pdnm1;
      
      recsum = 0.0;
      for (i = 0; i < jm; i++) recsum += 1.0 / (x - xjac[n - i]);

      delx = -poly / (pder - recsum*poly);
      x   += delx;
    } while (fabs (delx) > EPSDP && k++ < STOP);

    xjac[n - j] = x;
  }

  xjac[0] = -1.0;
  for (i = 1; i < nh; i++) xjac[i] = -xjac[n - i];

  if (np & 1) xjac[nh] = 0.0;
}


void zwgl (double*       z ,
	   double*       w ,
	   const integer np)
/* ------------------------------------------------------------------------- *
 * Gauss-Legendre points and weights.
 *
 * Generate np G-L points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.10).
 * ------------------------------------------------------------------------- */
{
  register integer i, n;
  double           poly, pder, polym1, pderm1, polym2, pderm2;

  if (np < 2) {
    z[0] = 0.0;  w[0] = 2.0;
    return;
  }
  
  n  = np - 1;
  jacg (n, 0.0, 0.0, z);

  for (i = 0; i < np; i++) {
    jacobf (np, z[i], 0.0, 0.0, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = 2.0 / ( (1.0 - SQR(z[i])) * SQR(pder) );
  }
}
  

void zwgrl (double*       z ,
	    double*       w ,
	    const integer np)
/* ------------------------------------------------------------------------- *
 * Gauss-Radau-Legendre points and weights.
 *
 * Generate np G-R-L points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.11).
 * ------------------------------------------------------------------------- */
{
  register integer i, n;
  double           poly, pder, polym1, pderm1, polym2, pderm2, con;

  if (np < 2) {
    z[0] = 0.0; w[0] = 2.0;
    return;
  }

  n   = np - 1;
  con = 1.0 / SQR(np);

  jacgr (n, 0.0, 0.0, z);
  
  w[0] = 2.0 * con;

  for (i = 1; i < np; i++) {
    jacobf (n, z[i], 0.0, 0.0, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = con * (1.0 * z[i]) / SQR(poly);
  }
}
   

void zwgll (double*       z ,
	    double*       w ,
	    const integer np)
/* ------------------------------------------------------------------------- *
 * Gauss-Lobatto-Legendre points and weights.
 *
 * Generate np G-L-L points (z) and weights (w) for integration over the
 * range [-1, 1] for a polynomial of degree n = np - 1.
 *
 * Reference: Canuto et al., eq (2.3.12).
 * ------------------------------------------------------------------------- */
{
  register integer i, n;
  double           poly, pder, polym1, pderm1, polym2, pderm2, con;

  if (np < 2) {
    z[0] = w[0] = 0.0;
    return;
  }

  if (np == 2) {
    z[0] = -(z[1] = w[0] = w[1] = 1.0);
    return;
  }

  n   = np - 1;
  con = 2.0 / (double) (n * np);
  
  jacgl (n, 0.0, 0.0, z);

  for (i = 0; i < np; i++) {
    jacobf (n, z[i], 0.0, 0.0, &poly,&pder,&polym1,&pderm1,&polym2,&pderm2);
    w[i] = con / (SQR(poly));
  }
}


double pnleg (const double  z,
	      const integer n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the nth order Legendre polynomial at z, based on the
 * recursion formula for Legendre polynomials.
 * ------------------------------------------------------------------------- */
{
  register integer k;
  register double  dk, p1, p2, p3;
 
  if (n == 0) return 1.0;

  p1 = 1.0;
  p3 = p2 = z;

  for (k = 1; k < n; k++) {
    dk = (double) k;
    p3 = ((2.0*dk + 1.0)*z*p2 - dk*p1) / (dk + 1.0);
    p1 = p2;
    p2 = p3;
  }
 
  return p3;
}


double pndleg (const double  z,
	       const integer n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the derivative of the nth order Legendre polynomial
 * at z, based on the recursion formula for Legendre polynomials.
 * ------------------------------------------------------------------------- */
{
  register integer k;
  register double  dk, p1, p2, p3, p1d, p2d, p3d;

  if (n == 0) return 0.0;

  p2  = z;
  p1d = 0.0;
  p1  = p2d = p3d = 1.0;

  for (k = 1; k < n; k++) {
    dk  = (double) k;
    p3  = ((2.0*dk + 1.0)*z*p2 - dk*p1) / (dk + 1.0);
    p3d = ((2.0*dk + 1.0)*p2 + (2.0*dk + 1.0)*z*p2d - dk*p1d) / (dk + 1.0);
    p1  = p2;
    p2  = p3;
    p1d = p2d;
    p2d = p3d;
  }

  return p3d;
}


double pnd2leg (const double  z,
		const integer n)
/* ------------------------------------------------------------------------- *
 * Compute the value of the second derivative of the nth order Legendre
 * polynomial at z, based on the definition of the singular Sturm-Liouville
 * problem that generates them (Canuto et al. eq 2.3.1):
 *
 *               (1 - z*z)L(z)'' - 2 z L' + n (n+1) L = 0.
 * ------------------------------------------------------------------------- */
{
  return (2.0 * z * pndleg (z, n) - n * (n+1) * pnleg (z, n) / (1.0 - SQR(z)));
}


#define U(i, j) (u[i*n + j])
#define A(i, j) (a[i*n + j])

void legtr2d (const integer n,
	      const double* u,
	      double*       a)
/* ------------------------------------------------------------------------- *
 * Compute 2D discrete Legendre transform of values in array u[n * n],
 * return in array a[n * n].  The spatial locations of u lie on the GLL grid.
 *
 * See Canuto et al. eqs 2.2.21--3, 2.3.13.
 * ------------------------------------------------------------------------- */
{
  const integer N = n - 1;
  register int  i, j, p, q;
  double        *z, *w;

  dQuadOps (LL, n, n, &z, 0, &w, 0, 0, 0, 0); 

  for (i = 0; i < n; i++) {
    const double ci = (i < N) ? 0.5 * (i + i + 1) : 0.5 * N;
    for (j = 0; j < n; j++) {
      const double cj = (j < N) ? 0.5 * (j + j + 1) : 0.5 * N;
      A(i, j) = 0.0;
      for (p = 0; p < n; p++) {
	 const double P = pnleg (z[p], i);
	 for (q = 0; q < n; q++) 
	   A(i, j) += w[p] * w[q] * U(p, q) * P * pnleg (z[q], j);
      }
      A(i, j) *= ci * cj;
    }
  }
}


void legtr2i (const integer n,
	      const double* u,
	      double*       a)
/* ------------------------------------------------------------------------- *
 * Compute 2D inverse discrete Legendre transform of values in array u[n * n],
 * return in array a[n * n].  The spatial locations of u lie on the GLL grid.
 *
 * See Canuto et al. eqs 2.2.21--3, 2.3.13.
 * ------------------------------------------------------------------------- */
{
  const integer N = n - 1;
  register int  i, j, p, q;
  double        *z;

  dQuadOps (LL, n, n, &z, 0, 0, 0, 0, 0, 0);

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      A(i, j) = 0.0;
      for (p = 0; p < n; p++) {
	 const double P = pnleg (z[i], p);
	 for (q = 0; q < n; q++) 
	   A(i, j) += U(p, q) * P * pnleg (z[j], q);
      }
    }
  }
}

#undef U
#undef A


void dgll (const integer  nz,
	   const double*  z ,
	   double**       D ,
	   double**       DT)
/* ------------------------------------------------------------------------- *
 * Compute the derivative operator matrix D and its transpose DT associated
 * with the nth order Lagrangian interpolants through the nz Gauss-Lobatto-
 * Legendre points z.
 *                     dU
 *                     --   = D_ij U_j evaluated at z = z_i
 *                     dz
 *
 * NB: D & DT are both nz x nz.  Canuto et al. eq (2.3.25).
 * ------------------------------------------------------------------------- */
{
  register integer i, j, n;
  register double  d0;

  if (nz < 2) {
    D[0][0] = DT[0][0] = 0.0;
    return;
  }
  
  n  = nz - 1;
  d0 = n * (n + 1) * 0.25;

  D[0][0] = -d0;
  for (i = 1; i < n; i++) D[i][i] = 0.0;
  D[n][n] =  d0;

  for (i = 0; i < nz; i++)
    for (j=0; j<nz; j++) {
      if (i != j) D[i][j] = pnleg (z[i], n) / (pnleg(z[j], n) * (z[i] - z[j]));
      DT[j][i] = D[i][j];
    }
}

					     
void uniknot (const integer nk,
	      double*       k )
/* ------------------------------------------------------------------------- *
 * Return nk knot points with uniform spacing on [-1, 1].
 * ------------------------------------------------------------------------- */
{
  register integer i, nh;
  register double  dx;

  if (nk < 2) {
    *k = 0.0;
    return;
  }

  nh = nk >> 1;

  dx      =  2.0 / (double) (nk - 1);
  k[0]    = -1.0;
  k[nk-1] =  1.0;
  for (i = 1; i < nh; i++) {
    k[i]      =  k[i-1] + dx;
    k[nk-1-i] = -k[i];
  }
  if (nk & 1) k[i] = 0.0;
}


integer quadComplete (const integer dim,
		      const integer np )
/* ------------------------------------------------------------------------- *
 * Return the number of Gauss-Legendre quadrature points sufficient to
 * achieve the full rate of convergence for tensor-product element bases.
 *
 * Dim is the number of space dimensions, np the number of points defining
 * basis polynomials.
 *
 * References: Hughes \S 4.1, Strang & Fix \S 4.3.
 * ------------------------------------------------------------------------- */
{
  register integer  n, ktot;

  ktot = (dim + 1)*(np - 1) - 2;
  n = (ktot & 0) ? ktot + 2 : ktot + 1;
  n >>= 1;

  return MAX (n, 2);
}
