*
      SUBROUTINE BICGSTAB( N, B, X, WORK, LDW, ITER, RESID, MATVEC,
     $                     PSOLVE, INFO )
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), B( * ), WORK( * )
*     ..
*     .. Function Arguments ..
      EXTERNAL           MATVEC, PSOLVE
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Convergence test: ( norm( b - A*x ) / norm( b ) ) < TOL.
*  For other measures, see the above reference.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  MATVEC  (external subroutine)
*          The user must provide a subroutine to perform the
*          matrix-vector product
*
*               y := alpha*A*x + beta*y,
*
*          where alpha and beta are scalars, x and y are vectors,
*          and A is a matrix. Vector x must remain unchanged.
*          The solution is over-written on vector y.
*
*          The call is:
*
*             CALL MATVEC( ALPHA, X, BETA, Y )
*
*          The matrix is passed into the routine in a common block.
*
*  PSOLVE  (external subroutine)
*          The user must provide a subroutine to perform the
*          preconditioner solve routine for the linear system
*
*               M*x = b,
*
*          where x and b are vectors, and M a matrix. Vector b must 
*          remain unchanged.
*          The solution is over-written on vector x.
*
*          The call is:
*
*             CALL PSOLVE( X, B )
*
*          The preconditioner is passed into the routine in a common block.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA < BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  ==============================================================
*     .. Local Scalars ..
*
*This variable used to communicate requests between BiCGSTAB() 
*  and BiCGSTABREVCOM()
*BiCGSTAB -> BiCGSTABREVCOM: 1 = init, 
*                            2 = use saved state to resume flow.
*BiCGSTABREVCOM -> BiCGSTAB: -1 = done, return to main, 
*                             1 = matvec using SCLR1/2, NDX1/2 
*                             2 = solve using NDX1/2
      INTEGER          IJOB
      LOGICAL          FTFLG
*
*     Arg/Result indices into WORK[].
      INTEGER NDX1, NDX2
*     Scalars passed from BiCGSTABREVCOM to BiCGSTAB.
      DOUBLE PRECISION SCLR1, SCLR2
*     Vars reqd for STOPTEST2
      DOUBLE PRECISION TOL, BNRM2
*     ..
*     .. External subroutines ..
      EXTERNAL         BICSTABGREVCOM, STOPTEST2
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Test the input parameters.
*
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT.MAX( 1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE.0 ) RETURN
*
*     Stop test may need some indexing info from REVCOM
*     use the init call to send the request across. REVCOM
*     will note these requests, and everytime it asks for
*     stop test to be done, it will provide the indexing info.
*
*     1 == R; 2 == RTLD; 3 == P; 4 == V; 5 == T; 6 == PHAT; 7 == SHAT;
*     8 == S; -1 == ignore; any other == error
      NDX1 = 1
      NDX2 = -1
      TOL = RESID
      FTFLG = .TRUE.
*
*     First time call always init.
*
      IJOB = 1

 1    CONTINUE

      CALL BICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO, 
     $                    NDX1, NDX2, SCLR1, SCLR2, IJOB)

*     On a return from REVCOM() we use the table (REVCOM -> BiCGSTAB)
*     to decode IJOB.
      IF (IJOB .eq. -1) THEN
*        revcom wants to terminate, so do it.
         GOTO 2
      ELSEIF (IJOB .eq. 1) THEN
*        call matvec.
         CALL MATVEC(SCLR1, WORK(NDX1), SCLR2, WORK(NDX2))
      ELSEIF (IJOB .eq. 2) THEN
*        call solve.
         CALL PSOLVE(WORK(NDX1), WORK(NDX2))
      ELSEIF (IJOB .eq. 3) THEN
*        call matvec with X.
         CALL MATVEC(SCLR1, X, SCLR2, WORK(NDX2))
      ELSEIF (IJOB .EQ. 4) THEN
*        do stopping test 2
*        if first time, set INFO so that BNRM2 is computed.
         IF( FTFLG ) INFO = -1
         CALL STOPTEST2(N, WORK(NDX1), B, BNRM2, RESID, TOL, INFO)
         FTFLG = .FALSE.
      ENDIF
*
*         done what revcom asked, set IJOB & go back to it.
          IJOB = 2
          GOTO 1
*
*     come here to terminate
*
 2    CONTINUE
*
      RETURN
*
*     End of BICGSTAB
*
      END
*
C 23456789012345678901234567890123456789012345678901234567890123456789012
C =======================================================================
C
      SUBROUTINE BICGSTABREVCOM(N, B, X, WORK, LDW, ITER, RESID, INFO,
     $                    NDX1, NDX2, SCLR1, SCLR2, IJOB)
*
*  -- Iterative template routine --
*     Univ. of Tennessee and Oak Ridge National Laboratory
*     October 1, 1993
*     Details of this algorithm are described in "Templates for the 
*     Solution of Linear Systems: Building Blocks for Iterative 
*     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra, 
*     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
*     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
*
*     .. Scalar Arguments ..
      INTEGER            N, LDW, ITER, INFO 
      DOUBLE PRECISION   RESID
      INTEGER            NDX1, NDX2
      DOUBLE PRECISION   SCLR1, SCLR2
      INTEGER            IJOB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), B( * ), WORK( LDW,* )
*     ..
*
*  Purpose
*  =======
*
*  BICGSTAB solves the linear system A*x = b using the
*  BiConjugate Gradient Stabilized iterative method with 
*  preconditioning.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER. 
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
* 
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  X       (input/output) DOUBLE PRECISION array, dimension N.
*          On input, the initial guess. This is commonly set to 
*          the zero vector. 
*          On exit, if INFO = 0, the iterated approximate solution.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDW,7)
*          Workspace for residual, direction vector, etc.
*          Note that vectors R and S shared the same workspace.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array WORK. LDW >= max(1,N).
*
*  ITER    (input/output) INTEGER
*          On input, the maximum iterations to be performed.
*          On output, actual number of iterations performed.
*
*  RESID   (input/output) DOUBLE PRECISION
*          On input, the allowable convergence measure for
*          norm( b - A*x ) / norm( b ).
*          On output, the final value of this measure.
*
*  INFO    (output) INTEGER
*
*          =  0: Successful exit. Iterated approximate solution returned.
*
*          >  0: Convergence to tolerance not achieved. This will be 
*                set to the number of iterations performed.
*
*          <  0: Illegal input parameter, or breakdown occurred
*                during iteration.
*
*                Illegal parameter:
*
*                   -1: matrix dimension N < 0
*                   -2: LDW < N
*                   -3: Maximum number of iterations ITER <= 0.
*                   -5: Erroneous NDX1/NDX2 in INIT call.
*                   -6: Erroneous RLBL.
*
*                BREAKDOWN: If parameters RHO or OMEGA become smaller
*                   than some tolerance, the program will terminate.
*                   Here we check against tolerance BREAKTOL.
*
*                  -10: RHO < BREAKTOL: RHO and RTLD have become 
*                                       orthogonal.
*                  -11: OMEGA < BREAKTOL: S and T have become 
*                                         orthogonal relative to T'*T.
*
*                  BREAKTOL is set in function GETBREAK.
*
*  NDX1    (input/output) INTEGER. 
*  NDX2    On entry in INIT call contain indices required by interface
*          level for stopping test.
*          All other times, used as output, to indicate indices into
*          WORK[] for the MATVEC, PSOLVE done by the interface level.
*
*  SCLR1   (output) DOUBLE PRECISION.
*  SCLR2   Used to pass the scalars used in MATVEC. Scalars are reqd because
*          original routines use dgemv.
*
*  IJOB    (input/output) INTEGER. 
*          Used to communicate job code between the two levels.
*
*  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*  ==============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT,
     $                   NEED1, NEED2
      DOUBLE PRECISION   TOL, ALPHA, BETA, RHO, RHO1, BNRM2, OMEGA, 
     $                   RHOTOL, OMEGATOL, GETBREAK, DDOT, DNRM2
*     indicates where to resume from. Only valid when IJOB = 2!
      INTEGER RLBL
*
*     saving all.
      SAVE
*     ..
*     .. External Functions ..
      EXTERNAL           GETBREAK, DAXPY, DCOPY, DDOT, DNRM2, DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Entry point, so test IJOB
      IF (IJOB .eq. 1) THEN
         GOTO 1
      ELSEIF (IJOB .eq. 2) THEN
*        here we do resumption handling
         IF (RLBL .eq. 2) GOTO 2
         IF (RLBL .eq. 3) GOTO 3
         IF (RLBL .eq. 4) GOTO 4
         IF (RLBL .eq. 5) GOTO 5
         IF (RLBL .eq. 6) GOTO 6
         IF (RLBL .eq. 7) GOTO 7
*        if neither of these, then error
         INFO = -6
         GOTO 20
      ENDIF
*
*
*****************
 1    CONTINUE
*****************
*
      INFO = 0
      MAXIT = ITER
      TOL   = RESID
*
*     Alias workspace columns.
*
      R    = 1
      RTLD = 2
      P    = 3
      V    = 4
      T    = 5
      PHAT = 6
      SHAT = 7
      S    = 1
*
*     Check if caller will need indexing info.
*
      IF( NDX1.NE.-1 ) THEN
         IF( NDX1.EQ.1 ) THEN
            NEED1 = ((R - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.2 ) THEN
            NEED1 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.3 ) THEN
            NEED1 = ((P - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.4 ) THEN
            NEED1 = ((V - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.5 ) THEN
            NEED1 = ((T - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.6 ) THEN
            NEED1 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.7 ) THEN
            NEED1 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX1.EQ.8 ) THEN
            NEED1 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED1 = NDX1
      ENDIF
*
      IF( NDX2.NE.-1 ) THEN
         IF( NDX2.EQ.1 ) THEN
            NEED2 = ((R - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.2 ) THEN
            NEED2 = ((RTLD - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.3 ) THEN
            NEED2 = ((P - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.4 ) THEN
            NEED2 = ((V - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.5 ) THEN
            NEED2 = ((T - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.6 ) THEN
            NEED2 = ((PHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.7 ) THEN
            NEED2 = ((SHAT - 1) * LDW) + 1
         ELSEIF( NDX2.EQ.8 ) THEN
            NEED2 = ((S - 1) * LDW) + 1
         ELSE
*           report error
            INFO = -5
            GO TO 20
         ENDIF
      ELSE
         NEED2 = NDX2
      ENDIF
*
*     Set parameter tolerances.
*
      RHOTOL = GETBREAK()
      OMEGATOL = GETBREAK()
*
*     Set initial residual.
*
      CALL DCOPY( N, B, 1, WORK(1,R), 1 )
      IF ( DNRM2( N, X, 1 ).NE.ZERO ) THEN
*********CALL MATVEC( -ONE, X, ONE, WORK(1,R) )
*        Note: using RTLD[] as temp. storage.
*********CALL DCOPY(N, X, 1, WORK(1,RTLD), 1)
         SCLR1 = -ONE
         SCLR2 = ONE
         NDX1 = -1
         NDX2 = ((R - 1) * LDW) + 1
*
*        Prepare for resumption & return
         RLBL = 2
         IJOB = 3
         RETURN
*
*****************
 2       CONTINUE
*****************
*
         IF ( DNRM2( N, WORK(1,R), 1 ).LE.TOL ) GO TO 30
      ENDIF
      CALL DCOPY( N, WORK(1,R), 1, WORK(1,RTLD), 1 )
*
      BNRM2 = DNRM2( N, B, 1 )
      IF ( BNRM2 .EQ. ZERO ) BNRM2 = ONE
*
      ITER = 0
*
   10 CONTINUE
*
*     Perform BiConjugate Gradient Stabilized iteration.
*
         ITER = ITER + 1
*
         RHO = DDOT( N, WORK(1,RTLD), 1, WORK(1,R), 1 )
         IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
*
*        Compute vector P.
*
         IF ( ITER.GT.1 ) THEN
            BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
            CALL DAXPY( N, -OMEGA, WORK(1,V), 1, WORK(1,P), 1 )
            CALL DSCAL( N, BETA, WORK(1,P), 1 )
            CALL DAXPY( N, ONE, WORK(1,R), 1, WORK(1,P), 1 )
         ELSE
            CALL DCOPY( N, WORK(1,R), 1, WORK(1,P), 1 )
         ENDIF
*
*        Compute direction adjusting vector PHAT and scalar ALPHA.
*
*********CALL PSOLVE( WORK(1,PHAT), WORK(1,P) )
*
         NDX1 = ((PHAT - 1) * LDW) + 1
         NDX2 = ((P    - 1) * LDW) + 1
*        Prepare for return & return
         RLBL = 3
         IJOB = 2
         RETURN
*
*****************
 3       CONTINUE
*****************
*
*********CALL MATVEC( ONE, WORK(1,PHAT), ZERO, WORK(1,V) )
*
         NDX1 = ((PHAT - 1) * LDW) + 1
         NDX2 = ((V    - 1) * LDW) + 1
*        Prepare for return & return
         SCLR1 = ONE
         SCLR2 = ZERO
         RLBL = 4
         IJOB = 1
         RETURN
*
*****************
 4       CONTINUE
*****************
*
         ALPHA = RHO / DDOT( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
*
*        Early check for tolerance.
*
         CALL DAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 )
         CALL DCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
         IF ( DNRM2( N, WORK(1,S), 1 ).LE.TOL ) THEN
            CALL DAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
            RESID = DNRM2( N, WORK(1,S), 1 ) / BNRM2
            GO TO 30
         ELSE
*
*           Compute stabilizer vector SHAT and scalar OMEGA.
*
************CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
*
            NDX1 = ((SHAT - 1) * LDW) + 1
            NDX2 = ((S    - 1) * LDW) + 1
*           Prepare for return & return
            RLBL = 5
            IJOB = 2
            RETURN
*
*****************
 5          CONTINUE
*****************
*
************CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
*
            NDX1 = ((SHAT - 1) * LDW) + 1
            NDX2 = ((T    - 1) * LDW) + 1
*           Prepare for return & return
            SCLR1 = ONE
            SCLR2 = ZERO
            RLBL = 6
            IJOB = 1
            RETURN
*
*****************
 6          CONTINUE
*****************
*
            OMEGA = DDOT( N, WORK(1,T), 1, WORK(1,S), 1 ) / 
     $              DDOT( N, WORK(1,T), 1, WORK(1,T), 1 )
*
*           Compute new solution approximation vector X.
*
            CALL DAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
            CALL DAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
*
*           Compute residual R, check for tolerance.
*
            CALL DAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
*
************RESID = DNRM2( N, WORK(1,R), 1 ) / BNRM2 
************IF ( RESID.LE.TOL  ) GO TO 30
*
            NDX1 = NEED1
            NDX2 = NEED2
*           Prepare for resumption & return
            RLBL = 7
            IJOB = 4
            RETURN
*
*****************
 7          CONTINUE
*****************
            IF( INFO.EQ.1 ) GO TO 30
*
            IF ( ITER.EQ.MAXIT ) THEN
               INFO = 1
               GO TO 20
            ENDIF
*
         ENDIF
*
      IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         GO TO 25
      ELSE
         RHO1 = RHO
         GO TO 10
      ENDIF
*
   20 CONTINUE
*
*     Iteration fails.
*
      RLBL = -1
      IJOB = -1
      RETURN
*
   25 CONTINUE
*
*     Set breakdown flag.
*
      IF ( ABS( RHO ).LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ).LT.OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RLBL = -1
      IJOB = -1
      RETURN
*
   30 CONTINUE
*
*     Iteration successful; return.
*
      INFO = 0
      RLBL = -1
      IJOB = -1
      RETURN
*
*     End of BICGSTABREVCOM
*
      END
*
C 23456789012345678901234567890123456789012345678901234567890123456789012
C =======================================================================
C
      SUBROUTINE STOPTEST2( N, R, B, BNRM2, RESID, TOL, INFO )
*
*     .. Scalar Arguments ..
      INTEGER            N, INFO
      DOUBLE PRECISION   RESID, TOL, BNRM2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   R( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  Computes the stopping criterion 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER.
*          On entry, the dimension of the matrix.
*          Unchanged on exit.
*
*  INFO    (output) INTEGER
*          On exit, 1/0 depending on whether stopping criterion
*          was met or not.
*
*  BNRM2   (input/output) DOUBLE PRECISION.
*          On first time entry, will be -1.0.
*          On first time exit will contain norm2(B)
*          On all subsequent entry/exit's unchanged.
*
*  RESID   (output) DOUBLE PRECISION.
*          On exit, the computed stopping measure.
*
*  TOL     (input) DOUBLE PRECISION.
*          On input, the allowable convergence measure.
*
*  R       (input) DOUBLE PRECISION array, dimension N.
*          On entry, the residual.
*          Unchanged on exit.
*
*  B       (input) DOUBLE PRECISION array, dimension N.
*          On entry, right hand side vector B.
*          Unchanged on exit.
*
*  BLAS CALLS:   DNRM2
*  ============================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER        ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. External Routines ..
      DOUBLE PRECISION DNRM2
      EXTERNAL         DNRM2
*     ..
*     .. Executable Statements ..
*
      IF( INFO.EQ.-1 ) THEN
         BNRM2 = DNRM2( N, B, 1 )
         IF ( BNRM2.EQ.ZERO ) BNRM2 = ONE
      ENDIF
*
      RESID = DNRM2( N, R, 1 ) / BNRM2
*
      INFO = 0
      IF ( RESID.LE.TOL )
     $     INFO = 1
*
      RETURN
*
*     End of STOPTEST2
*
      END
*     ================================================================
      DOUBLE PRECISION FUNCTION GETBREAK()
*
*     Get breakdown parameter tolerance; for the test routine,
*     set to machine precision.
*
      DOUBLE PRECISION EPS, DLAMCH
*
      EPS = DLAMCH('EPS')
      GETBREAK = EPS**2
*
      RETURN
*
      END
