C12345678901234567890123456789012345678901234567890123456789012345678901
C
C     $Id$
C
C     Matrix-matrix, matrix-vector multiply routines,
C     designed to be called from C.
C     E.g. where C = A * B; the FORTRAN equivalent is C' = B' * A'.
C
C     ------------------------------------------------------------------
C     C = A * B.  (As written in C, with row-major storage.)
C
      SUBROUTINE DMXM (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
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
      SUBROUTINE SMXM (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
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
C     ------------------------------------------------------------------
C     C += A * B.
C
      SUBROUTINE DMXMA (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXMA (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
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
C     ------------------------------------------------------------------
C     C -= A * B.
C
      SUBROUTINE DMXMS (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXMS (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
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
C     ------------------------------------------------------------------
C               t
C     C -= A * B.
C
      SUBROUTINE DMXMTS (A, NRA, B, NCA, C, NCBT)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCBT, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCA, NCBT), C(NCBT, NRA)
      DO J = 1, NCBT
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(K, J) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXMTS (A, NRA, B, NCA, C, NCBT)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCBT, I, J, K
      REAL     A(NCA, NRA), B(NCA, NCBT), C(NCBT, NRA)
      DO J = 1, NCBT
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(K, J) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C     ------------------------------------------------------------------
C     Y = A * X.  Matrix-vector multiply.
C
      SUBROUTINE DMXV (A, NRA, X, NCA, Y)
      IMPLICIT NONE
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
      SUBROUTINE SMXV (A, NRA, X, NCA, Y)
      IMPLICIT NONE
      INTEGER NRA, NCA, I, J
      REAL    A(NCA, NRA), X(NRA), Y(NCA)
      DO I = 1, NCA
         Y(I) = 0.0
         DO J = 1, NRA
            Y(I) = Y(I) + A(J, I) * X(J)
         ENDDO
      ENDDO
      RETURN
      END
