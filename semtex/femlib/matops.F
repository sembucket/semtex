C12345678901234567890123456789012345678901234567890123456789012345678901
C
C     -- Matrix operators for gradients, etc.
C
C     $Id$
C
C
      SUBROUTINE DGRAD2 (X, Y, UR, US, DV, DT, NP, NEL)
C
C     Ur = dX/dr, Us = dY/ds over all nel elements each np*np.
C     Ur & Us are updated by addition, so should be zero on entry.
C     Conceivably Ur & Us can be aliased to one another, so
C     the operations are broken out into two loops.
C
      IMPLICIT NONE
C
      INTEGER          NP, NEL, I, J, K, L
      DOUBLE PRECISION  X(NP,NP,NEL),  Y(NP,NP,NEL)
      DOUBLE PRECISION UR(NP,NP,NEL), US(NP,NP,NEL)
      DOUBLE PRECISION DV(NP,NP),     DT(NP,NP)
C
#if defined (_SX) || defined (__uxp__)
      DO K = 1, NP
         DO J = 1, NP
            DO L = 1, NP
               DO I = 1, NEL
                  UR(K,J,I) = UR(K,J,I) + X(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      DO J = 1, NP
         DO L = 1, NP
            DO K = 1, NP
               DO I = 1, NEL
                  US(K,J,I) = US(K,J,I) + Y(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#elif defined (__sgi) || defined (__alpha)
      CALL DGEMM ('N','N',NP,NP*NEL,NP,1.0D0,DT,NP,X,NP,1.0D0,UR,NP)
      DO I = 1, NEL
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0,Y(1,1,I),
     +        NP,DV,NP,1.0D0,US(1,1,I),NP)
      ENDDO
#else
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  UR(K,J,I) = UR(K,J,I) + X(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  US(K,J,I) = US(K,J,I) + Y(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#endif
      RETURN
      END
C
C
      SUBROUTINE SGRAD2 (X, Y, UR, US, DV, DT, NP, NEL)
C
C     -- Ur = dX/dr, Us = dY/ds over all nel elements each np*np.
C
      IMPLICIT NONE
C
      INTEGER  NP, NEL, I, J, K, L
      REAL      X(NP,NP,NEL),  Y(NP,NP,NEL)
      REAL     UR(NP,NP,NEL), US(NP,NP,NEL)
      REAL     DV(NP,NP),     DT(NP,NP)
C
#if defined (_SX) || defined (__uxp__)
      DO K = 1, NP
         DO J = 1, NP
            DO L = 1, NP
               DO I = 1, NEL
                  UR(K,J,I) = UR(K,J,I) + X(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      DO J = 1, NP
         DO L = 1, NP
            DO K = 1, NP
               DO I = 1, NEL
                  US(K,J,I) = US(K,J,I) + Y(K,L,I) * DT(J,L)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
#elif defined (__sgi) || defined (__alpha)
      CALL SGEMM ('N','N',NP,NP*NEL,NP,1.0D0,DT,NP,X,NP,1.0D0,UR,NP)
      DO I = 1, NEL
         CALL SGEMM ('N','N',NP,NP,NP,1.0D0,Y(1,1,I),
     +        NP,DV,NP,1.0D0,US(1,1,I),NP)
      ENDDO
#else
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  UR(K,J,I) = UR(K,J,I) + X(L,J,I) * DV(L,K)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
      DO I = 1, NEL
         DO J = 1, NP
            DO L = 1, NP
               DO K = 1, NP
                  US(K,J,I) = US(K,J,I) + Y(K,L,I) * DT(J,L)
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
