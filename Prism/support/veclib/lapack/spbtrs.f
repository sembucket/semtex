      SUBROUTINE SPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
*
*  -- LAPACK routine (version 1.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL               AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  SPBTRS solves a system of linear equations A*X = B with a symmetric
*  positive definite band matrix A using the Cholesky factorization
*  A = U'*U or A = L*L' computed by SPBTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the factor stored in AB is upper or lower
*          triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of super-diagonals of the matrix A if UPLO = 'U',
*          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) REAL array, dimension (LDAB,N)
*          The triangular factor U or L from the Cholesky factorization
*          A = U'*U or A = L*L' of the band matrix A, stored in the
*          first KD+1 rows of the array.  The j-th column of U or L is
*          stored in the array AB as follows:
*          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors B for the system of
*          linear equations.
*          On exit, the solution vectors, X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           STBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B where A = U'*U.
*
         DO 10 J = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
*
*           Solve U*X = B, overwriting B with X.
*
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Solve A*X = B where A = L*L'.
*
         DO 20 J = 1, NRHS
*
*           Solve L*X = B, overwriting B with X.
*
            CALL STBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
*
*           Solve L'*X = B, overwriting B with X.
*
            CALL STBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB,
     $                  LDAB, B( 1, J ), 1 )
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of SPBTRS
*
      END
