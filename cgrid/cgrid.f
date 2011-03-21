      PROGRAM HYPGG
C     ************************************************************************
C     HYPERBOLIC GRID GENERATION.
C     
C     BY Z. ALSALIHI, 1987, VKI INSTITUTE TN 162.
C
C     $Id$
C     ************************************************************************
      PARAMETER (ND1=420,ND2=90)
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*4    CPU,CPU2
      REAL*8    X(ND1,ND2),Y(ND1,ND2)
      DIMENSION R(ND1,ND2,2),RJ(ND1,ND2),XI(ND1,ND2,2),ETA(ND1,ND2,2)
     1     ,BA(2,2),BF(2),AA(2,2,ND1),BB(2,2,ND1),CC(2,2,ND1),YY(2,ND1),
     2     D(2,ND1),IP(2,ND1)

C     -- INPUT MODULE.

      OPEN(UNIT=3, FILE='cgrid.in',  STATUS='old')     
      OPEN(UNIT=2, FILE='cgrid.out', STATUS='unknown')     

      READ (3,*) IMAX, JMAX, EPS, NTYPE
      WRITE(6,*) IMAX, JMAX, EPS, NTYPE
      READ (3,*) (R(I,1,1),R(I,1,2), I = 1, IMAX)
      DO I = 1, IMAX
         READ(3,*) (RJ(I, J), J = 1, JMAX)
      ENDDO

C     -- START SOLUTION.

      IM = (IMAX+1) / 2
C     -- START ON THE BODY.
      J  = 1

C     -- CALCULATION OF THE XI AND ETA DERIVATIVES ON THE KNOWN J LINE.

 700  DO L = 1, 2
         DO I = 2, IMAX - 1
            XI(I,J,L) = (R(I+1,J,L) - R(I-1,J,L)) * 0.5D00
         ENDDO
         XI(1,J,L) = (-3.0D00*R(1,J,L)+4.0D00*R(2,J,L)-R(3,J,L))*0.5D00
         XI(IMAX,J,L) = (3.0D00*R(IMAX,J,L)-4.0D00*R(IMAX-1,J,L)+
     1        R(IMAX-2,J,L))*0.5D00
      ENDDO
      DO I = 1, IMAX
         SQ         =  RJ(I,J)/(XI(I,J,1)**2+XI(I,J,2)**2)
         ETA(I,J,1) = -XI(I,J,2)*SQ
         ETA(I,J,2) =  XI(I,J,1)*SQ
      ENDDO

C     --  CALCULATION OF MATRICES.

      DO KK = 2, IMAX - 1
         K  = KK - 1
         A1 = ETA(KK,J,1)
         A2 = ETA(KK,J,2)
         B1 = XI(KK,J,1)
         B2 = XI(KK,J,2)

         DET = B1*B1+B2*B2
         B1  = B1/DET
         B2  = B2/DET

         BA(1,1) = (B1*A1-B2*A2)*0.500D00
         BA(1,2) = (B1*A2+B2*A1)*0.500D00
         BA(2,1) = (B2*A1+B1*A2)*0.500D00
         BA(2,2) = (B2*A2-B1*A1)*0.500D00

         F     =  RJ(KK,J)+RJ(KK,J+1)
         BF(1) = -B2*F+R(KK,J,1)
         BF(2) =  B1*F+R(KK,J,2)

         DO L = 1, 2
            DO N = 1, 2
               AA(L,N,K) =  0.00D00
               BB(L,N,K) =  BA(L,N)
               CC(L,N,K) = -BA(L,N)
            ENDDO
            AA(L,L,K) = 1.00D00
            YY(L,K)   = BF(L)
         ENDDO
      ENDDO
      
      AA(1,1,1) =  1.0
      AA(1,2,1) =  CC(1,2,1)
      AA(2,1,1) =  0.0
      AA(2,2,1) =  1.0+CC(2,2,1)
      YY(1,1)   = -CC(1,1,1)*R(1,J,1)+YY(1,1)
      YY(2,1)   = -CC(2,1,1)*R(1,J,1)+YY(2,1)

C     -- AT I=IMAX-1.

      KZ=IMAX-2
      AA(1,1,KZ)=1.0
      AA(1,2,KZ)=BB(1,2,KZ)
      AA(2,1,KZ)=0.0
      AA(2,2,KZ)=1.0+BB(2,2,KZ)
      DFE=-BB(1,1,KZ)*R(IMAX,J,1)+YY(1,KZ)
      YY(1,IMAX-2)=DFE
      DFE=-BB(2,1,KZ)*R(IMAX,J,1)+YY(2,KZ)
      YY(2,IMAX-2)=DFE

C     -- END OF COMPUTATION OF MATRICES.

      DO L=1,2
         DO M=1,2
            CC(L,M,1)=0.0
            BB(L,M,KZ)=0.0
         ENDDO
      ENDDO

C     -- CALCULATION OF NUMERICAL DISSIPATION TERMS.

      IF(EPS.EQ.0.OR.J.LT.4) GO TO 9
      DO L=1,2
         DO K=IM-20,IM+20
            D(L,K)=R(K+2,J,L)-4.0D00*R(K+1,J,L)+6.0D00*R(K,J,L)-
     1           4.0D00*R(K-1,J,L)+R(K-2,J,L)
            YY(L,K-1)=YY(L,K-1)-D(L,K)*EPS
         ENDDO
      ENDDO

C     -- SOLVE THE SYSTEM OF EQUATIONS.

 9    CALL DECBT(2,IMAX-2,AA,BB,CC,IP,IER)
      IF(IER.NE.0) STOP 9999
      CALL SOLBT(2,IMAX-2,AA,BB,CC,YY,IP)
      DO L=1,2
         DO K=1,IMAX-2
            R(K+1,J+1,L)=YY(L,K)
         ENDDO
      ENDDO
      R(1,J+1,1)=R(1,J,1)
      R(IMAX,J+1,1)=R(IMAX,J,1)
      R(1,J+1,2)=R(2,J+1,2)
      R(IMAX,J+1,2)=R(IMAX-1,J+1,2)
      IF(R(IM,J+1,1).GT.R(IM-1,J+1,1).AND.
     1     R(IM,J+1,1).GT.R(IM+1,J+1,1)) THEN
         R(IM,J+1,1)=R(IM,J+1,1)-
     1        2.0D00*(DABS(R(IM-1,J+1,1)-R(IM,J+1,1))
     1        +DABS(R(IM+1,J+1,1)-R(IM,J+1,1)))*.25D00
      ENDIF
C     -- NTYPE = 2 CORRESPONDS TO "O-GRID".
      IF(NTYPE.EQ.2) THEN
         IF(R(1,J+1,1).LT.R(2,J+1,1)) THEN
            R(1,J+1,1)=(R(IMAX-1,J+1,1)+R(2,J+1,1))*0.5D00
         ENDIF
         R(1,J+1,2)=0.0D00
         R(IMAX,J+1,1)=R(1,J+1,1)
         R(IMAX,J+1,2)=R(1,J+1,2)
      END IF

C     -- CONTINUE FOR THE NEXT J LINE IF J IS LT JMAX.

      IF(J.LT.JMAX) THEN
         J=J+1
         GO TO 700
      END IF

C     -- OUTPUT MODULE.
      
C      WRITE(2,199)
C 199  FORMAT ('1',20X, 'HYPERBOLIC GRID GENERATION'//)
C      WRITE(2,2000) IMAX,JMAX,EPS
C 2000 FORMAT(5X,'IMAX = ',I3,2X,'JMAX = ',I3,2X,
C     1     ' NUM. DIS. COF. = ',F10.5/)
C      DO I = 1, IMAX
C         DO J = 1, JMAX
C            WRITE(2,*)I,J,R(I,J,1),R(I,J,2)
C         ENDDO
C      ENDDO

      WRITE(2,2000) (IMAX-1)*(JMAX-1)
 2000 FORMAT('    2    2    1',1X,I5,'   NR   NS   NZ  NEL')
      DO I = 1, IMAX-1
         DO J = 1, JMAX-1
            WRITE(2,*) R(I,   J,   1), R(I,   J,   2)
            WRITE(2,*) R(I+1, J,   1), R(I+1, J,   2)
            WRITE(2,*) R(I,   J+1, 1), R(I,   J+1, 2)
            WRITE(2,*) R(I+1, J+1, 1), R(I+1, J+1, 2)
         ENDDO
      ENDDO

      STOP 001
      END





      SUBROUTINE DECBT (M,N,A,B,C,IP,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(M,M,N),B(M,M,N),C(M,M,N),IP(M,N)
C-----------------------------------------------------------------------------
C     BLOCK-TRIDIAGONAL MATRIX DECOMPOSITION ROUTINE.
C     WRITTEN BY A.C. HINDMARSH,
C     THE INPUT MATRIX CONTAINS THREE BLOCKS OF ELEMENTS IN EACH BLOCK ROW,
C     INCLUDING BLOCKS IN THE (1.3) AND (N,N-2) BLOCK POSITIONS.
C     DECBT USES BLOCK GAUSS ELIMINATION AND S
C     FOR SOLUTION OF BLOCKS, PARTIAL PIVOTING IS DONE WITHIN
C     BLOCK ROWS ONLY.CINPUT
C     M=ORDER OF EACH BLOCK.
C     N=NUMBER OF BLOCKS IN EACH DIRECTION OF THE MATRIX.
C     N MUST BE 4 OR MORE.  THE COMPLETE MATRIX HAS ORDER M*N.
C     A=M BY M BY N ARRAYCONTAINING DIAGONAL BLOCKS.
C     A(I,J,K) CONTAINS THE (I.J) ELEMENT OF THE K-TH BLOCK.
C     B=M BY M BY N ARRAY CONTAINING THE SUPER-DIAGONAL BLOCKS
C     (IN B(,,K) FOR K=1,....N-1) AND THE BLOCK IN THE (N,N-2)
C     C=M BY M BY N ARRAY CONTAINING THE SUBDIAGONAL BLOCKS
C     (IN C(,,K) FOR K=2,3,....N) AND THE BLOCK IN THE 
C     (1.3) BLOCK POSITION (IN C(,,1)).
C     IP=INTEGER ARRAY OF LENGTH M+N FOR WORKING STORAGE.
C     OUTPUT..
C     A,B,C=M BY M BYN ARRAYS CONTAINING THE BLOCK LU DECOMPOSITION
C     OF THE INPUT MATRIX.
C     IP=M BY N ARRAY OF PIVOT INFORMATION.IP9,,K0 CONTAINS
C     INFORMATION FOR THE K-TH DIAGONAL BLOCK.
C     IER=O IF NO TROUBLE OCCURRED, OR
C     =-1 IF THE INPUTVALUE OF M OR N WAS ILLEGAL, OR 
C     =K IF A SINGULAR MATRIX WAS FOUND IN THE K-TH DIAGONAL BLOCK.
C     USE SOLBT TO SOLVE THE ASSOCIATED LINEAR SYSTEM
C     DECBT CALLS SUBROUTINES DEC (M,MO,A,IP,IER0 AND SOL9M,MO,A,Y,IP)
C     FOR SOLUTION OF M BY M LINEAR SYSTEMS.
C-----------------------------------------------------------------------------
      IF (M.LT.1.OR.N.LT.4) GO TO 210
      NM1=N-1
      NM2=N-2
C     PROCESS THE FIRST BLOCK-ROW.--------------------------------------------
      CALL DEC (M,M,A,IP,IER)
      K=1
      IF (IER .NE.0) GO TO 200
      DO J=1,M
         CALL SOL (M,M,A,B(1,J,1),IP)
         CALL SOL (M,M,A,C(1,J,1),IP)
      ENDDO

C     ADJUST B(,,2).----------------------------------------------------------
      DO J=1,M
         DO I=1,M
            DP=0.
            DO L=1,M
               DP = DP+C(I,L,2)*C(L,J,1)
            ENDDO
          B(I,J,2)= B(I,J,2)-DP
         ENDDO
      ENDDO

C     MAIN LOOP.  PROCESS BLOCK-ROWS 2 TO N-1. -------------------------------
      DO K=2,NM1
         KM1=K-1
         DO J=1,M
            DO I=1,M
               DP=0.
               DO L=1,M
                  DP=DP+C(I,L,K)*B(L,J,KM1)
               ENDDO
               A(I,J,K)=A(I,J,K)-DP
            ENDDO
         ENDDO

         CALL DEC (M,M,A(1,1,K), IP(1,K),IER)
         IF (IER.NE.0) GO TO 200
         DO J=1,M
            CALL SOL (M,M,A(1,1,K),B(1,J,K),IP(1,K))
         ENDDO
      ENDDO
      
C     PROCESS LAST BLOCK-ROW AND RETURN --------------------------------------
      DO J=1,M
         DO I=1,M
            DP=0
            DO L=1,M
               DP= DP+B(I,L,N)*B(L,J,NM2)
            ENDDO
            C(I,J,N)=C(I,J,N)-DP
         ENDDO
      ENDDO

      DO J=1,M
         DO I=1,M
            DP=0.
            DO L=1,M
               DP=DP+C(I,L,N)*B(L,J,NM1)
            ENDDO
            A(I,J,N)=A(I,J,N)-DP
         ENDDO
      ENDDO

      CALL DEC (M,M,A(1,1,N),IP(1,N),IER)
      K=N
      IF (IER.NE.0) GO TO 200
      RETURN
C     ERROR RETURNS-.-------------------------------------------------
 200  IER=K
      RETURN
 210  IER=-1
      RETURN
      END
C
C
C
      SUBROUTINE SOLBT (M,N,A,B,C,Y,IP)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------------
C     SOLUTION OF BLOCK-TRIDIAGONAL LINEAR SYSTEM
C     COEFFICIENT MATRIX MUST HAVE BEEN PREVIOUSLY PROCESSEDBY DECBT.
C     INPUT
C     M=ORDER OF EACH BLOCK.
C     N=NUMBER OF BLOCKS IN EACH DIRECTION OF MATRIX
C     A,B,C=M BY M BY N ARRAYS CONTAINING BLOCK LU DECOMPOSITION
C     OF COEFFICIENT MATRIX FROM DECBT.
C     P=M BY N INTEGER ARRAY OF PIVOT INFORMATION FROM DECBT .
C     Y=ARRAY OF LENGHT M*N CONTAINING THE RIGHT-HAND SIDE VECTOR
C     (TREATED AS AN M BY N ARRAY HERE).
C     OUTPUT..
C     Y=SOLUTION VECTOR OF LENGTH M*N.
C     SOLBT MAKES CALLS TO SUBROUTINE SOL(M,MO,A,Y,IP)
C     FOR SOLUTION OF M BY M LINEAR SYSTEMS.
C-----------------------------------------------------------------------------
      DIMENSION A(M,M,N),B(M,M,N),C(M,M,N),Y(M,N),IP(M,N)
C
      NM1=N-1
      NM2=N-2
C     FORWARD SOLUTION SWEEP
      CALL SOL(M,M,A,Y,IP)
      DO K=2,NM1
         KM1=K-1
         DO I=1,M
            DP=0.
            DO J=1,M
               DP=DP+C(I,J,K)*Y(J,KM1)
            ENDDO
            Y(I,K)=Y(I,K)-DP
         ENDDO
         CALL SOL (M,M,A(1,1,K),Y(1,K),IP(1,K))
      ENDDO

      DO I=1,M
         DP=0
         DO J=1,M
            DP=DP+C(I,J,N)*Y(J,NM1)+B(I,J,N)*Y(J,NM2)
         ENDDO
         Y(I,N)=Y(I,N)-DP
      ENDDO

      CALL SOL (M,M,A(1,1,N),Y(1,N),IP(1,N))

C     BACKWARD SOLUTION SWEEP. -----------------------------------------------

      DO KB=1,NM1
         K=N-KB
         KP1=K+1
         DO I=1,M
            DP=0.
            DO J=1,M
               DP=DP+B(I,J,K)*Y(J,KP1)
            ENDDO
            Y(I,K)=Y(I,K)-DP
         ENDDO
      ENDDO

      DO I=1,M
         DP=0.
         DO J=1,M
            DP=DP+C(I,J,1)*Y(J,3)
         ENDDO
         Y(I,1)=Y(I,1)-DP
      ENDDO
      RETURN
      END
C
C
C
C
C
C
      SUBROUTINE DEC  (N,NDIM,A,IP,IER)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------------
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------------
C     MATRIX TRIANGULARISATION BY GAUSS ELIMINATION WITH PARTIAL PIVOTING 
C     INPUT
C     N=ORDER OF MATRIX
C     NDIM=DECLARED FIRST DIMENSION OF ARRAY A.
C     A=MATRIX TO BE TRIANGULARIZED.
C     OUTPUT..
C     A(I,J),I.LE.J=UPPER TRIANGULAR FACTOR,U.
C     A(I,J),I.GT.J=MULTIPLIERS=LOWER TRIANGULAR FACTOR, I-L.
C     IP(K),K.LT.N=INDEX OF  K-TH PIVOT ROW.
C     IER=0 IF MATRIX A IS NONSINGULAR,OR  K IF FOUND TO BE
C     SINGULAR AT STAGE K.
C     ROW INTERCHANGES ARE FINISHED IN U.ONLY PARTLY IN L.
C     USE SOL TO OBTAIN SOLUTION OF LINEAR SYSTEM
C     IF IER .NE. O,A IS SINGULAR,SOL WILL DIVIDE BY ZERO
C-----------------------------------------------------------------------------
      IER=0
      IF (N.EQ.1) GO TO 70
      NM1=N-1
      DO K=1,NM1
         KP1=K+1
C     FIND THE PIVOT IN COLUMN K.SEARCH ROWS K TO N.------------------------
         M=K
         DO I=KP1,N
            IF (ABS(A(I,K)) .GT.ABS(A(M,K))) M=I
         ENDDO
         IP(K)=M
C     INTERCHANGE ELEMENTS IN RWOS K AND M.---------------------------------
         T=A(M,K)
         IF (M.EQ.K) GO TO 20
         A(M,K)=A(K,K)
         A(K,K)=T
 20      IF (T.EQ.0.) GO TO 80
C     STORE MULTIPLIERS IN A (I,K), I=K+1,....N.----------------------------
         T=1./T
         DO I=KP1,N
            A(I,K)=-A(I,K)*T
         ENDDO
C     APPLY MULTIPLIERS TO OTHER COLUMNS OFF A.-----------------------------

         DO 50 J=KP1,N
            T=A(M,J)
            A(M,J)=A(K,J)
            A(K,J)=T
            IF (T.EQ.0) GO TO 50
            DO I=KP1,N
               A(I,J)=A(I,J) +A(I,K)*T
            ENDDO
 50         CONTINUE
         ENDDO

 70      K=N
         IF (A(N,N).EQ.0.) GO TO 80
         RETURN
 80      IER=K
         RETURN
         END
C
C
C
C
      SUBROUTINE SOL  (N,NDIM,A,B,IP)
      IMPLICIT REAL*8(A-H,O-Z)
C-----------------------------------------------------------------------------
      DIMENSION A(NDIM,N),B(N),IP(N)
C-----------------------------------------------------------------------------
C     SOLUTION OF LINEAR SYSTEME A=X =B USING OUTPUT OF DEC
C     INPUT..
C     NDIM=DECLARED FIRST DIMENSION OF ARRAY A.
C     N=ORDER OF MATRIX.
C     A=TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C     B=RIGHT HAND SIDE VECTOR.
C     IP=PIVOT INFORMATION VECTOR OBTAINED FROM DEC.
C     DO NOT USE IF DEC HAS SET IER.NE.O.
C     OUTPUT..
C     B=SOLUTION VECTOR,X.
C-----------------------------------------------------------------------------
      IF (N.EQ.1) GO TO 50
      NM1=N-1
C     APPLY ROW PERMUTATIONS AND MULTIPLIERS TO B.--------------------------
      DO K=1,NM1
         KP1=K+1
         M=IP(K)
         T=B(M)
         B(M)=B(K)
         B(K)=T
         DO I=KP1,N
            B(I)=B(I)+A(I,K)*T
         ENDDO
      ENDDO

C     BACK SOLVE.-----------------------------------------------------------
      DO KB=1,NM1
         KM1=N-KB
         K= KM1+1
         B(K)=B(K)/A(K,K)
         T=-B(K)
         DO  I=1,KM1
            B(I)=B(I)+A(I,K)*T
         ENDDO
      ENDDO

 50   B(1)=B(1)/A(1,1)
      RETURN
      END
