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
C     This operator can be used for any field operation that can be 
C     described by the tensor-product forms
C
C     UR = [DV]  X
C     US =  Y  [DT]
C
C     Typically for tensor-product forms, DT = transpose (DV).
C     This routine is designed for use with row-major (C-like) matrices.
C
      IMPLICIT NONE
C
      INTEGER          NP, NEL, I, J, K, L
      DOUBLE PRECISION  X(NP,NP,NEL),  Y(NP,NP,NEL)
      DOUBLE PRECISION UR(NP,NP,NEL), US(NP,NP,NEL)
      DOUBLE PRECISION DV(NP,NP),     DT(NP,NP)
C
#if defined (_SX) || defined (__uxp__)
      IF (NEL.EQ.1) THEN
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0,DT,NP, X,NP,1.0D0,UR,NP)
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0, Y,NP,DV,NP,1.0D0,US,NP)
      ELSE
         DO K = 1, NP
            DO J = 1, NP
               DO L = 1, NP
                  DO I = 1, NEL
                     UR(K,J,I) = UR(K,J,I) + X(L,J,I) * DV(L,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#if defined (_SX)
!CDIR NOVECTOR
#endif
         DO J = 1, NP
#if defined (_SX)
!CDIR NOVECTOR
#endif
            DO L = 1, NP
#if defined (_SX)
!CDIR NOVECTOR
#endif
               DO K = 1, NP
                  DO I = 1, NEL
                     US(K,J,I) = US(K,J,I) + Y(K,L,I) * DT(J,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
#else
      CALL DGEMM ('N','N',NP,NP*NEL,NP,1.0D0,DT,NP,X,NP,1.0D0,UR,NP)
      DO I = 1, NEL
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0,Y(1,1,I),
     +        NP,DV,NP,1.0D0,US(1,1,I),NP)
      ENDDO
#endif
      RETURN
      END
C
C
      SUBROUTINE DTPR2D (X, Y, T, DV, DT, NP, NEL)
C
C     This operator can be used for any field operation that can be 
C     described by the tensor-product form
C
C     Y  = [DV] X [DT]  (here X and Y may be aliases).
C
C     Typically for tensor-product forms, DT = transpose (DV).
C     This routine is designed for use with row-major (C-like) matrices.
C
C     Input T is workspace, the same size as X & Y.
C
      IMPLICIT NONE
C
      INTEGER          NTOT, NP, NEL, I, J, K, L
      DOUBLE PRECISION  X(NP,NP,NEL),  Y(NP,NP,NEL), T(NP,NP,NEL)
      DOUBLE PRECISION DV(NP,NP),     DT(NP,NP)
C
#if defined (_SX) || defined (__uxp__)
      IF (NEL.EQ.1) THEN
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0, X,NP,DV,NP,0.0D0,T,NP)
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0,DT,NP, T,NP,0.0D0,Y,NP)
      ELSE
         NTOT = NP*NP*NEL
         DO I = 1, NTOT
            T(I,1,1) = 0.0D0
         ENDDO
#if defined (_SX)
!CDIR NOVECTOR
#endif
         DO J = 1, NP
#if defined (_SX)
!CDIR NOVECTOR
#endif
            DO L = 1, NP
#if defined (_SX)
!CDIR NOVECTOR
#endif
               DO K = 1, NP
                  DO I = 1, NEL
                     T(K,J,I) = T(K,J,I) + X(K,L,I) * DT(J,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO I = 1, NTOT
            Y(I,1,1) = 0.0D0
         ENDDO
         DO K = 1, NP
            DO J = 1, NP
               DO L = 1, NP
                  DO I = 1, NEL
                     Y(K,J,I) = Y(K,J,I) + T(L,J,I) * DV(L,K)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
#else
      CALL DGEMM ('N','N',NP,NP*NEL,NP,1.0D0,DT,NP,X,NP,0.0D0,T,NP)
      DO I = 1, NEL
         CALL DGEMM ('N','N',NP,NP,NP,1.0D0,T(1,1,I),
     +        NP,DV ,NP,0.0D0,Y(1,1,I),NP)
      ENDDO
#endif
      RETURN
      END
