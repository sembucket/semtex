C12345678901234567890123456789012345678901234567890123456789012345678901
C
C     $Id$
C
C     Matrix-matrix, matrix-vector multiply routines,
C     designed to be called from C.
C     E.g. where C = A * B; the FORTRAN equivalent is C' = B' * A'.
C
C
C
C
      SUBROUTINE DMXM (A, NRA, B, NCA, C, NCB)
C
C     C = A * B.
C
      IMPLICIT NONE
C
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
C
      DO J = 1, NCB
         DO I = 1, NRA
            C(J, I) = 0.0D0
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C
      SUBROUTINE SMXM (A, NRA, B, NCA, C, NCB)
C
C     C = A * B.
C
      IMPLICIT NONE
C
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
C
      DO J = 1, NCB
         DO I = 1, NRA
            C(J, I) = 0.0
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C
C
C
      SUBROUTINE DMXMA (A, NRA, B, NCA, C, NCB)
C
C     C += A * B.
C
      IMPLICIT NONE
C
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
C
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C
      SUBROUTINE SMXMA (A, NRA, B, NCA, C, NCB)
C
C     C += A * B.
C
      IMPLICIT NONE
C
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
C
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

C
C
C
C
      SUBROUTINE DMXMS (A, NRA, B, NCA, C, NCB)
C
C     C -= A * B.
C
      IMPLICIT NONE
C
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
C
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C
      SUBROUTINE SMXMS (A, NRA, B, NCA, C, NCB)
C
C     C -= A * B.
C
      IMPLICIT NONE
C
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
C
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

C
C
C
C
      SUBROUTINE DMXV (A, NRA, X, NCA, Y)
C
C     Y = A * X.
C
      IMPLICIT NONE
C
      INTEGER          NRA, NCA, I, J
      DOUBLE PRECISION A(NCA, NRA), X(NRA), Y(NCA)
C     
C     -- ALTERNATIVE ORDERING FOR TESTING:
C
C      DO I = 1, NCA
C         Y(I) = 0.0D0
C      ENDDO
C      
C      DO J = 1, NRA
C         DO I = 1, NCA
C            Y(I) = Y(I) + A(J, I) * X(J)
C         ENDDO
C      ENDDO
      DO I = 1, NCA
         Y(I) = 0.0D0
         DO J = 1, NRA
            Y(I) = Y(I) + A(J, I) * X(J)
         ENDDO
      ENDDO
      RETURN
      END
C
C
      SUBROUTINE SMXV (A, NRA, X, NCA, Y)
C
C     Y = A * X.
C
      IMPLICIT NONE
C
      INTEGER NRA, NCA, I, J
      REAL    A(NCA, NRA), X(NRA), Y(NCA)
C
      DO I = 1, NCA
         Y(I) = 0.0
         DO J = 1, NRA
            Y(I) = Y(I) + A(J, I) * X(J)
         ENDDO
      ENDDO
      RETURN
      END
