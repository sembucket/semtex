C12345678901234567890123456789012345678901234567890123456789012345678901
C
C     -- Matrix operators for gradients, etc.
C
C     $Id$
C
C
      SUBROUTINE DGRADR (U, UR, DV, NP, NEL)
C
C     -- Ur = dU/dr, over all NEL elements each NP*NP.
C
      IMPLICIT NONE
C
      INTEGER          NP, NEL, I, J, K, L
      DOUBLE PRECISION U(NP,NP,NEL), UR(NP,NP,NEL), DV(NP,NP)
C
#if   defined (_SX)
      DO J = 1, NP
         DO K = 1, NP
            DO L = 1, NP
               DO I = 1, NEL
                  UR(K,J,I) = UR(K,J,I) + U(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#elif defined (__uxp__)
C     -- NB restricted to 512 elements.
      DO J = 1, NP
         DO L = 1, NP
!OCL UNROLL(4)
            DO K = 1, NP
!OCL NOVREC
!OCL REPEAT(512)
               DO I = 1, NEL
                  UR(K,J,I) = UR(K,J,I) + U(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#else
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  UR(K,J,I) = UR(K,J,I) + U(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#endif
      RETURN
      END
C
C
      SUBROUTINE DGRADS (U, US, DT, NP, NEL)
C
C     -- Us = dU/ds over all NEL elements each NP*NP.
C
      IMPLICIT NONE
C
      INTEGER          NP, NEL, I, J, K, L
      DOUBLE PRECISION U(NP,NP,NEL), US(NP,NP,NEL), DV(NP,NP)
C
#if   defined (_SX)
      DO J = 1, NP
         DO L = 1, NP
            DO K = 1, NP
*VDIR NODEP, UNROLL=2
               DO I = 1, NEL
                  US(K,J,I) = US(K,J,I) + U(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#elif defined (__uxp__)
C     -- NB restricted to 512 elements.
      DO J = 1, NP
         DO L = 1, NP
!OCL UNROLL(4)
            DO K = 1, NP
!OCL NOVREC
!OCL REPEAT(512)
               DO I = 1, NEL
                  US(K,J,I) = US(K,J,I) + U(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#else
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  US(K,J,I) = US(K,J,I) + U(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#endif
      RETURN
      END
C
C
      SUBROUTINE SGRADR (U, UR, DV, NP, NEL)
C
C     -- Ur = dU/dr, over all NEL elements each NP*NP.
C
      IMPLICIT NONE
C
      INTEGER NP, NEL, I, J, K, L
      REAL    U(NP,NP,NEL), UR(NP,NP,NEL), DV(NP,NP)
C
#if   defined (_SX)
      DO J = 1, NP
         DO K = 1, NP
            DO L = 1, NP
               DO I = 1, NEL
                  UR(K,J,I) = UR(K,J,I) + U(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#elif defined (__uxp__)
C     -- NB restricted to 512 elements.
      DO J = 1, NP
         DO L = 1, NP
!OCL UNROLL(4)
            DO K = 1, NP
!OCL NOVREC
!OCL REPEAT(512)
               DO I = 1, NEL
                  UR(K,J,I) = UR(K,J,I) + U(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#else
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  UR(K,J,I) = UR(K,J,I) + U(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#endif
      RETURN
      END
C
C
      SUBROUTINE SGRADS (U, US, DT, NP, NEL)
C
C     -- Us = dU/ds over all NEL elements each NP*NP.
C
      IMPLICIT NONE
C
      INTEGER NP, NEL, I, J, K, L
      REAL    U(NP,NP,NEL), US(NP,NP,NEL), DV(NP,NP)
C
#if   defined (_SX)
      DO J = 1, NP
         DO L = 1, NP
            DO K = 1, NP
*VDIR NODEP, UNROLL=2
               DO I = 1, NEL
                  US(K,J,I) = US(K,J,I) + U(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#elif defined (__uxp__)
C     -- NB restricted to 512 elements.
      DO J = 1, NP
         DO L = 1, NP
!OCL UNROLL(4)
            DO K = 1, NP
!OCL NOVREC
!OCL REPEAT(512)
               DO I = 1, NEL
                  US(K,J,I) = US(K,J,I) + U(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#else
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  US(K,J,I) = US(K,J,I) + U(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#endif
      RETURN
      END
C
C
      SUBROUTINE DMXMC (A, NRA, B, NCA, C, NCB)
C
C     -- C += A * B, but to be called from C so compute B * A.
C
      IMPLICIT NONE
C
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA,NRA), B(NCB,NCA), C(NCB,NRA)
C
      DO I = 1, NRA
         DO K = 1, NCA
            DO J = 1, NCB
               C(J,I) = C(J,I) + B(J,K) * A(K,I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C
      SUBROUTINE SMXMC (A, NRA, B, NCA, C, NCB)
C
C     -- C += A * B, but to be called from C so compute B * A.
C
      IMPLICIT NONE
C
      INTEGER NRA, NCA, NCB, I, J, K
      REAL    A(NCA,NRA), B(NCB,NCA), C(NCB,NRA)
C
      DO I = 1, NRA
         DO K = 1, NCA
            DO J = 1, NCB
               C(J,I) = C(J,I) + B(J,K) * A(K,I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
