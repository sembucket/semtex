      SUBROUTINE DSPMVC( TRANS, N, M, A, COLPTR, ROWIND, X, LDX, 
     $                   Y, LDY )
*     ..
*     .. Scalar Arguments .. 
      INTEGER            LDX, LDY, M, N, TRANS 
*     ..
*     .. Array Arguments .. 
      INTEGER            ROWIND( * ), COLPTR( * )       
      DOUBLE PRECISION   A( * ), X( LDX, * ), Y( LDY, * ) 
*
*  Purpose
*  =======
*
*  computes the product of a sparse matrix and (dense block) vector(s)
*
*                   Y = A * X   if trans = 0 
*                   Y = A'* X   otherwise 
*
*  where A is stored in the compressed column format.
*
*  Arguments
*  =========
*
*  TRANS   (input) INTEGER
*          If TRANS = 0, compute Y = A*X
*          If TRANS = 1, compute Y = A'*X
*
*  N       (input) INTEGER
*          The order of the square matrix A
*
*  M       (input) INTEGER
*          the number of columns in the (block) vectors X. 
*
*  A       (input) DOUBLE PRECISION array, dimension (NZ) 
*          the numerical values of the nonzero elements in matrix A.
*          NZ is the total number of nonzeros. 
*
*  COLPTR  (input) INTEGER array, dimension ( N+1 )
*          Column pointer of matrix A. 
*
*  ROWIND  (input) INTEGER array, dimension ( NZ ) 
*          Row indices of matrix A. 
*
*  X       (input) DOUBLE PRECISION, dimension ( N, M ) 
*          the block vectors X to be multiplied. 
*
*  LDX     (input) INTEGER
*          The leading dimension of array X, LDX >= max( 1, N ).
*
*  Y       (output) DOUBLE PRECISION, dimension ( N, M )
*          the product of the matrix - block vectors.
*
*  LDY     (input) INTEGER
*          The leading dimension of array Y, LDY >= max( 1, N ).
* 
*  ==============================================================
*
*     .. Parameter .. 
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*
*     .. Local Scalars .. 
      INTEGER            I, J, K, L
*
*     Initialization 
*
      DO 20 L = 1, M 
         DO 10 I = 1, N
	    Y( I, L ) = ZERO
  10     CONTINUE
  20  CONTINUE
*
      IF( TRANS.EQ.0 )THEN
*
*        Compute y = A*x
*
         DO 50 L = 1, M  
   	    DO 40 J = 1,N
	       DO 30 K = COLPTR( J ), COLPTR( J+1 ) - 1
                  I = ROWIND( K )
*
*                 Compute y(i) = y(i) + a(i,j) * x(j)
*
	          Y( I, L ) = Y( I, L ) + A( K )*X( J, L )
*
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE 
*
      ELSE
*
*        Compute y = A'*x 
*
         DO 80 L = 1, M 
   	    DO 70 J = 1,N
	       DO 60 K = COLPTR( J ), COLPTR( J+1 ) - 1 
                  I = ROWIND( K )
*
*                 Compute y(j) = y(j) + a(i,j) * x(i)
*
	          Y( J, L ) = Y( J, L ) + A( K )* X( I, L )
*
   60          CONTINUE
   70  	    CONTINUE
   80    CONTINUE 
*
      ENDIF
*
      RETURN
*
*     End of DSPMVC
*
      END
