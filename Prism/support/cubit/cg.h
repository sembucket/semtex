#ifndef CG_H
#define CG_H

/* CG
 *
 * $Revision$
 *
 * Author:    R. D. Henderson
 *
 * CG solves the symmetric positive-definite system A*x = b using the
 * conjugate gradient method.
 *
 * Arguments:
 *
 *   n            (input) size of the solution vector
 *
 *   x            (input/output) initial guess and final solution
 *
 *   b            (input) right-hand-side
 *
 *   max_iter     (input/output) maximum number of iterations to take; on
 *                return it gives the iterations required for convergence
 *
 *   tol          (input/output) desiered tolerance for the solution vector,
 *                ||b - A*x|| < tol*||b|| using the supplied norm; on return
 *                it gives the final residual
 *
 * 
 * Call-Back Functions
 *
 * You need to pass pointers to three functions for CG to do its work.
 * These functions provide a norm for the solution vector, a preconditioner
 * for the iterations, and a routine to apply the operator to a given 
 * vector.  Here are the prototypes for each one:
 *
 *   double (*norm) (int n, const double *x)
 *        ...returns the norm of x
 *
 *   void (*M_solve) (int n, const double *x, double *y)
 *        ...solve for y such that M y = x; M is the preconditioner
 *
 *   void (*A_apply) (int n, const double *x, double *y)
 *        ...compute the product y = A*x; A is the operator you're inverting
 *
 * The return value of CG is 0 if the iteration converges, and 1 if it does 
 * not.
 *
 * Based on CG, pg. 15 of the SIAM Templates book.
 * ------------------------------------------------------------------------- */
 
int CG (int n, double *x, const double *b, int *max_iter, double *tol,
	double (*norm)(), double (*dot)(), void (*M_solve)(), 
	void (*A_apply)());

#endif

