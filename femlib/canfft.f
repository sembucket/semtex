C123456789012345678901234567890123456789012345678901234567890123456789012
C     ===================================================================
C     Reference: Spectral Methods in Fluid Dynamics, App. B.
C     Canuto, Hussaini, Quarteroni & Zang, Springer-Verlag.  1988.
C
C     -- User routines:
C     FACTOR:  compute prime factors.
C     xPREFT:  set up prime factors, trigonometric data.
C     xMRCFT:  multiple real--complex 1D FFTs.
C
C     -- The following are not intended to be called by the user:
C     FACTOR:  compute prime factors.
C     xFFT1:   multiple complex--complex 1D FFT.
C     xPASS1:  FFT kernel.
C
C     $Id$
C     ===================================================================
C
C
      SUBROUTINE SPREFT (N,NFAX,IFAX,TRIG)
C     -------------------------------------------------------------------
C     Set up trigonometric array, prime factor array prior to calling
C     other routines.
C
C     NB: note sign change on the SINes compared to Canuto.
C     -------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          K, N, NFAX, IFAX(*)
      REAL             TRIG(2,0:N-1)
      DOUBLE PRECISION ARG, TWOPI
      PARAMETER        ( TWOPI = 6.28318530717958647688D0 )
C     
      CALL FACTOR (N,NFAX,IFAX)
      DO 10 K = 0, N-1
         ARG = TWOPI * K / N
         TRIG (1,K) =  DCOS (ARG)
         TRIG (2,K) = -DSIN (ARG)
 10   CONTINUE
C     
      RETURN
      END
C
C
      SUBROUTINE FACTOR (N,NFAX,IFAX)
C     -------------------------------------------------------------------
C     Compute the (2, 3) prime factors of N.
C     If there is an error, NFAX is returned as zero.
C     -------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  N, NFAX, IFAX(*), II, NN
C
      NFAX = 0
      NN   = N
      DO 10 II = 1 ,20
         IF (NN .EQ. 3*(NN/3)) THEN
            NFAX = NFAX + 1
            IFAX(NFAX) = 3
            NN = NN / 3
         ELSE
            GO TO 20
         END IF
 10   CONTINUE
 20   CONTINUE
      DO 30 II = NFAX + 1, 20
         IF (NN .EQ. 2*(NN/2)) THEN
            NFAX = NFAX + 1
            IFAX(NFAX) = 2
            NN = NN / 2
         ELSE
            GO TO 40
         END IF
 30   CONTINUE
 40   CONTINUE
      IF (NN .NE. 1) THEN
         NFAX = 0
      END IF
C     
      RETURN
      END
C
C
      SUBROUTINE SMRCFT (V,NP,NZ,W,NFAX,IFAX,TRIG,ISIGN)
C     -------------------------------------------------------------------
C     Compute multiple 1D real--complex transforms of v, return in place.
C     
C     For the real-->complex transform (ISIGN=+1), V is interpreted on 
C     input as a 2D array with each row as a distinct data set: DFTs are
C     to be carried out along each of the NP rows, each of which has NZ
C     points in it.  NP must be even, NZ must have prime factors
C     of 2 & 3.
C
C     After the transform, the real and imaginary parts of the V~ are
C     stored alternating in the rows of V.  For NZ even, the Nyquist
C     frequency datum is stored in the imaginary storage location of
C     the zeroth mode Fourier component.
C
C     W is used as workspace, total amount of storage is also NP*NZ.
C     -------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  NP, NPH, NZ, NZH, NZHM, ISIGN, NFAX, IFAX(*)
      INTEGER  I, J, IP, JE, JO, MJ
      REAL     V    (NP,0:NZ-1)
      REAL     W    (NP/2,2,0:NZ-1)
      REAL     TRIG (2,0:NZ-1)
C 
      NPH  = NP / 2
      NZH  = NZ / 2
      NZHM = NZH - 1
C         
      IF (ISIGN .EQ. +1) THEN
         CALL SFFT1 (V,W,NZ,NFAX,IFAX,+1,TRIG,NPH)
         DO 10 J = 1, NZHM
            JE = J  + J
            JO = JE + 1
            MJ = NZ - J
            DO 20 I = 1, NPH
               IP = I + NPH
               V (I,   0) = W (I, 1,   0)
               V (I,   1) = W (I, 1, NZH)
               V (IP,  0) = W (I, 2,   0)
               V (IP,  1) = W (I, 2, NZH)
               V (I,  JE) = 0.5 * (W (I, 1,  J) + W (I, 1, MJ))
               V (I,  JO) = 0.5 * (W (I, 2,  J) - W (I, 2, MJ))
               V (IP, JE) = 0.5 * (W (I, 2,  J) + W (I, 2, MJ))
               V (IP, JO) = 0.5 * (W (I, 1, MJ) - W (I, 1,  J))
 20         CONTINUE
 10      CONTINUE
C     
      ELSE
C     
         DO 30 J = 1, NZHM
            JE = J  + J
            JO = JE + 1
            MJ = NZ - J
            DO 40 I = 1, NPH
               IP = I + NPH
               W (I, 1,   0) = V (I,  0)
               W (I, 2,   0) = V (IP, 0)
               W (I, 1, NZH) = V (I,  1)
               W (I, 2, NZH) = V (IP, 1)
               W (I, 1,   J) = V (I,  JE) - V (IP, JO)
               W (I, 2,   J) = V (I,  JO) + V (IP, JE)
               W (I, 1,  MJ) = V (I,  JE) + V (IP, JO)
               W (I, 2,  MJ) = V (IP, JE) - V (I,  JO)
 40         CONTINUE
 30      CONTINUE
         CALL SFFT1 (W,V,NZ,NFAX,IFAX,-1,TRIG,NPH)
C     
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE SFFT1 (A,C,N,NFAX,IFAX,ISIGN,TRIG,LEN)
C     -------------------------------------------------------------------
C     Performs complex FFT on multiple data in A, return the result in C.
C     No normalization.
C     
C     A:     Input array (destroyed during calculation):
C              A is dimensioned A (LEN, 2, 0:N-1);
C              First index labels distinct data to be transformed,
C              Second index labels real (1) or imaginary (2) parts,
C              Third index labels the N points to be transformed.
C              (Transform is along the third dimension).
C     C:     Output array (must be distinct from the input array):
C              C is dimensioned the same as A.
C     N:     Number of points in transform direction, must have
C              only prime factors of 2 and 3.
C     NFAX:  Number of prime factors of N.
C     IFAX:  Integer array containing prime factors of N.
C     ISIGN: Transform direction.
C              +1 to compute Fourier coefficients,
C              -1 to compute real space values.
C     TRIG:  Array containing trigonometric factors:
C              TRIG is dimensioned (2, 0:N-1);
C              TRIG (1, J) =    COS (2 * PI * J / N),
C              TRIG (2, J) =  - SIN (2 * PI * J / N).
C     LEN:   Number of distinct transforms to be performed.
C     -------------------------------------------------------------------
      IMPLICIT NONE
      REAL     A(LEN,2,0:N-1)
      REAL     C(LEN,2,0:N-1)
      REAL     TRIG(2,0:N-1)
      INTEGER  ISIGN, N, LEN, NFAX, IFAX(*), IFAC, LA, I, IJ
      LOGICAL  ODD
C
      LA  = 1
      ODD = .TRUE.
C     
      DO 10 I = 1, NFAX
         IFAC = IFAX (I)
         IF (ODD) THEN
            CALL SPASS1 (A,C,N,ISIGN,IFAC,LA,TRIG,LEN)
         ELSE
            CALL SPASS1 (C,A,N,ISIGN,IFAC,LA,TRIG,LEN)
         END IF
         ODD = .NOT. ODD
         LA  = LA * IFAC
 10   CONTINUE
C
      IF (ODD) THEN
         DO 30 I = 0, N-1
            DO 20 IJ = 1, LEN
               C (IJ,1,I) = A (IJ,1,I)
               C (IJ,2,I) = A (IJ,2,I)
 20         CONTINUE
 30      CONTINUE
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE SPASS1 (A,C,N,ISIGN,IFAC,LA,TRIG,LEN)
C     -------------------------------------------------------------------
C     Performs one pass of FFT.
C
C     This routine is never called directly by the user.
C     The arguments are similar to those of FFT1.
C     
C     The following constants presume a 64-bit machine.  They are the
C     sines of 45 degrees and 60 degrees.
C     -------------------------------------------------------------------
      IMPLICIT  NONE
      INTEGER   IND(0:20),JND(0:20)
      INTEGER   N, LEN, ISIGN, IFAC
      INTEGER   M, I, J, K, L, LA, I0, I1, I2, J0, J1, J2, JUMP, IJ
      REAL      A(LEN,2,0:N-1)
      REAL      C(LEN,2,0:N-1)
      REAL      TRIG(2,0:N-1)
      REAL      SN60, CC, SS, C1, C2, T1, T2
      REAL      AM1, AM2, TA1, TA2, AP1, AP2
      REAL      S1, S2, ROOT, ASN60
      PARAMETER ( ROOT  = 0.707106781186548 )
      PARAMETER ( ASN60 = 0.866025403784439 )
C
      SN60 = ISIGN * ASN60
      M    = N / IFAC
      DO 10 K = 0, IFAC-1
         IND(K) = K * M
         JND(K) = K * LA
 10   CONTINUE
C
      I    = 0
      J    = 0
      JUMP = (IFAC-1) * LA
      DO 130 K = 0, M-LA, LA
         DO 120 L = 1, LA
            IF (IFAC.EQ.2) THEN
               I0 = IND(0) + I
               I1 = IND(1) + I
               J0 = JND(0) + J
               J1 = JND(1) + J
               CC =         TRIG(1,K)
               SS = ISIGN * TRIG(2,K)
               IF (K.EQ.0) THEN
                  DO 20 IJ = 1, LEN
                     C(IJ,1,J0) = A(IJ,1,I0) + A(IJ,1,I1)
                     C(IJ,2,J0) = A(IJ,2,I0) + A(IJ,2,I1)
                     C(IJ,1,J1) = A(IJ,1,I0) - A(IJ,1,I1)
                     C(IJ,2,J1) = A(IJ,2,I0) - A(IJ,2,I1)
 20               CONTINUE
               ELSE
                  DO 50 IJ = 1, LEN
                     C(IJ,1,J0) = A(IJ,1,I0) + A(IJ,1,I1)
                     C(IJ,2,J0) = A(IJ,2,I0) + A(IJ,2,I1)
                     AM1 = A(IJ,1,I0) - A(IJ,1,I1)
                     AM2 = A(IJ,2,I0) - A(IJ,2,I1)
                     C(IJ,1,J1) = CC * AM1 - SS * AM2
                     C(IJ,2,J1) = SS * AM1 + CC * AM2
 50               CONTINUE
               END IF
            ELSE IF (IFAC.EQ.3) THEN
               I0 =IND(0) + I
               I1 =IND(1) + I
               I2 =IND(2) + I
               J0 =JND(0) + J
               J1 =JND(1) + J
               J2 =JND(2) + J
               IF (K.EQ.0) THEN
                  DO 60 IJ = 1, LEN
                     AP1 = A(IJ,1,I1) + A(IJ,1,I2)
                     AP2 = A(IJ,2,I1) + A(IJ,2,I2)
                     C(IJ,1,J0) = A(IJ,1,I0) + AP1
                     C(IJ,2,J0) = A(IJ,2,I0) + AP2
                     TA1 = A(IJ,1,I0) - 0.5 * AP1
                     TA2 = A(IJ,2,I0) - 0.5 * AP2
                     AM1 = SN60 * (A(IJ,1,I1) - A(IJ,1,I2))
                     AM2 = SN60 * (A(IJ,2,I1) - A(IJ,2,I2))
                     C(IJ,1,J1) = TA1 - AM2
                     C(IJ,2,J1) = TA2 + AM1
                     C(IJ,1,J2) = TA1 + AM2
                     C(IJ,2,J2) = TA2 - AM1
 60               CONTINUE
               ELSE
                  C1 =         TRIG(1,K)
                  C2 =         TRIG(1,2*K)
                  S1 = ISIGN * TRIG(2,K)
                  S2 = ISIGN * TRIG(2,2*K)
                  DO 70 IJ = 1, LEN
                     AP1 = A(IJ,1,I1) + A(IJ,1,I2)
                     AP2 = A(IJ,2,I1) + A(IJ,2,I2)
                     C(IJ,1,J0) = A(IJ,1,I0) + AP1
                     C(IJ,2,J0) = A(IJ,2,I0) + AP2
                     TA1 = A(IJ,1,I0) - 0.5 * AP1
                     TA2 = A(IJ,2,I0) - 0.5 * AP2
                     AM1 = SN60 * (A(IJ,1,I1) - A(IJ,1,I2))
                     AM2 = SN60 * (A(IJ,2,I1) - A(IJ,2,I2))
                     T1 = TA1 - AM2
                     T2 = TA2 + AM1
                     C(IJ,1,J1) = C1 * T1 - S1 * T2
                     C(IJ,2,J1) = S1 * T1 + C1 * T2
                     T1 = TA1 + AM2
                     T2 = TA2 - AM1
                     C(IJ,1,J2) = C2 * T1 - S2 * T2
                     C(IJ,2,J2) = S2 * T1 + C2 * T2
 70               CONTINUE
               END IF
            END IF
            I = I + 1
            J = J + 1
 120     CONTINUE
         J = J + JUMP
 130  CONTINUE
C     
      RETURN
      END
C
C
      SUBROUTINE DPREFT (N,NFAX,IFAX,TRIG)
C     -------------------------------------------------------------------
C     See comments for SPREFT.
C     -------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          K, N, NFAX, IFAX(*)
      DOUBLE PRECISION TRIG(2,0:N-1), ARG, TWOPI
      PARAMETER        ( TWOPI = 6.28318530717958647688D0 )
C     
      CALL FACTOR (N,NFAX,IFAX)
      DO 10 K = 0, N-1
         ARG = TWOPI * K / N
         TRIG (1,K) =  DCOS (ARG)
         TRIG (2,K) = -DSIN (ARG)
 10   CONTINUE
C     
      RETURN
      END
C
C
      SUBROUTINE DMRCFT (V,NP,NZ,W,NFAX,IFAX,TRIG,ISIGN)
C     -------------------------------------------------------------------
C     See comments for SMRCFT.
C     -------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          NP, NPH, NZ, NZH, NZHM, ISIGN, NFAX, IFAX(*)
      INTEGER          I, J, IP, JE, JO, MJ
      DOUBLE PRECISION V    (NP,0:NZ-1)
      DOUBLE PRECISION W    (NP/2,2,0:NZ-1)
      DOUBLE PRECISION TRIG (2,0:NZ-1)
C 
      NPH  = NP / 2
      NZH  = NZ / 2
      NZHM = NZH - 1
C         
      IF (ISIGN .EQ. +1) THEN
         CALL DFFT1 (V,W,NZ,NFAX,IFAX,+1,TRIG,NPH)
         DO 10 J = 1, NZHM
            JE = J  + J
            JO = JE + 1
            MJ = NZ - J
            DO 20 I = 1, NPH
               IP = I + NPH
               V (I,   0) = W (I, 1,   0)
               V (I,   1) = W (I, 1, NZH)
               V (IP,  0) = W (I, 2,   0)
               V (IP,  1) = W (I, 2, NZH)
               V (I,  JE) = 0.5 * (W (I, 1,  J) + W (I, 1, MJ))
               V (I,  JO) = 0.5 * (W (I, 2,  J) - W (I, 2, MJ))
               V (IP, JE) = 0.5 * (W (I, 2,  J) + W (I, 2, MJ))
               V (IP, JO) = 0.5 * (W (I, 1, MJ) - W (I, 1,  J))
 20         CONTINUE
 10      CONTINUE
C
      ELSE
C     
         DO 30 J = 1, NZHM
            JE = J  + J
            JO = JE + 1
            MJ = NZ - J
            DO 40 I = 1, NPH
               IP = I + NPH
               W (I, 1,   0) = V (I,  0)
               W (I, 2,   0) = V (IP, 0)
               W (I, 1, NZH) = V (I,  1)
               W (I, 2, NZH) = V (IP, 1)
               W (I, 1,   J) = V (I,  JE) - V (IP, JO)
               W (I, 2,   J) = V (I,  JO) + V (IP, JE)
               W (I, 1,  MJ) = V (I,  JE) + V (IP, JO)
               W (I, 2,  MJ) = V (IP, JE) - V (I,  JO)
 40         CONTINUE
 30      CONTINUE
         CALL DFFT1 (W,V,NZ,NFAX,IFAX,-1,TRIG,NPH)
C
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE DFFT1 (A,C,N,NFAX,IFAX,ISIGN,TRIG,LEN)
C     -------------------------------------------------------------------
C     See comments for SFFT1.
C     -------------------------------------------------------------------
      IMPLICIT         NONE
      DOUBLE PRECISION A(LEN,2,0:N-1)
      DOUBLE PRECISION C(LEN,2,0:N-1)
      DOUBLE PRECISION TRIG(2,0:N-1)
      INTEGER          ISIGN, N, LEN, NFAX, IFAX(*), IFAC, LA, I, IJ
      LOGICAL          ODD
C
      LA  = 1
      ODD = .TRUE.
C     
      DO 10 I = 1, NFAX
         IFAC = IFAX (I)
         IF (ODD) THEN
            CALL DPASS1 (A,C,N,ISIGN,IFAC,LA,TRIG,LEN)
         ELSE
            CALL DPASS1 (C,A,N,ISIGN,IFAC,LA,TRIG,LEN)
         END IF
         ODD = .NOT. ODD
         LA  = LA * IFAC
 10   CONTINUE
C
      IF (ODD) THEN
         DO 30 I = 0, N-1
            DO 20 IJ = 1, LEN
               C (IJ,1,I) = A (IJ,1,I)
               C (IJ,2,I) = A (IJ,2,I)
 20         CONTINUE
 30      CONTINUE
      END IF
C
      RETURN
      END
C
C
      SUBROUTINE DPASS1 (A,C,N,ISIGN,IFAC,LA,TRIG,LEN)
C     -------------------------------------------------------------------
C     See comments for SPASS1.
C     -------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          IND(0:20),JND(0:20)
      INTEGER          N, LEN, ISIGN, IFAC
      INTEGER          M, I, J, K, L, LA, I0, I1, I2
      INTEGER          J0, J1, J2, JUMP, IJ
      DOUBLE PRECISION A(LEN,2,0:N-1)
      DOUBLE PRECISION C(LEN,2,0:N-1)
      DOUBLE PRECISION TRIG(2,0:N-1)
      DOUBLE PRECISION SN60, CC, SS, C1, C2, T1, T2
      DOUBLE PRECISION AM1, AM2, TA1, TA2, AP1, AP2
      DOUBLE PRECISION S1, S2, ROOT, ASN60
      PARAMETER        ( ROOT  = 0.707106781186548D0 )
      PARAMETER        ( ASN60 = 0.866025403784439D0 )
C
      SN60 = ISIGN * ASN60
      M    = N / IFAC
      DO 10 K = 0, IFAC-1
         IND(K) = K * M
         JND(K) = K * LA
 10   CONTINUE
C
      I    = 0
      J    = 0
      JUMP = (IFAC-1) * LA
      DO 130 K = 0, M-LA, LA
         DO 120 L = 1, LA
            IF (IFAC.EQ.2) THEN
               I0 = IND(0) + I
               I1 = IND(1) + I
               J0 = JND(0) + J
               J1 = JND(1) + J
               CC =         TRIG(1,K)
               SS = ISIGN * TRIG(2,K)
               IF (K.EQ.0) THEN
                  DO 20 IJ = 1, LEN
                     C(IJ,1,J0) = A(IJ,1,I0) + A(IJ,1,I1)
                     C(IJ,2,J0) = A(IJ,2,I0) + A(IJ,2,I1)
                     C(IJ,1,J1) = A(IJ,1,I0) - A(IJ,1,I1)
                     C(IJ,2,J1) = A(IJ,2,I0) - A(IJ,2,I1)
 20               CONTINUE
               ELSE
                  DO 50 IJ = 1, LEN
                     C(IJ,1,J0) = A(IJ,1,I0) + A(IJ,1,I1)
                     C(IJ,2,J0) = A(IJ,2,I0) + A(IJ,2,I1)
                     AM1 = A(IJ,1,I0) - A(IJ,1,I1)
                     AM2 = A(IJ,2,I0) - A(IJ,2,I1)
                     C(IJ,1,J1) = CC * AM1 - SS * AM2
                     C(IJ,2,J1) = SS * AM1 + CC * AM2
 50               CONTINUE
               END IF
            ELSE IF (IFAC.EQ.3) THEN
               I0 =IND(0) + I
               I1 =IND(1) + I
               I2 =IND(2) + I
               J0 =JND(0) + J
               J1 =JND(1) + J
               J2 =JND(2) + J
               IF (K.EQ.0) THEN
                  DO 60 IJ = 1, LEN
                     AP1 = A(IJ,1,I1) + A(IJ,1,I2)
                     AP2 = A(IJ,2,I1) + A(IJ,2,I2)
                     C(IJ,1,J0) = A(IJ,1,I0) + AP1
                     C(IJ,2,J0) = A(IJ,2,I0) + AP2
                     TA1 = A(IJ,1,I0) - 0.5 * AP1
                     TA2 = A(IJ,2,I0) - 0.5 * AP2
                     AM1 = SN60 * (A(IJ,1,I1) - A(IJ,1,I2))
                     AM2 = SN60 * (A(IJ,2,I1) - A(IJ,2,I2))
                     C(IJ,1,J1) = TA1 - AM2
                     C(IJ,2,J1) = TA2 + AM1
                     C(IJ,1,J2) = TA1 + AM2
                     C(IJ,2,J2) = TA2 - AM1
 60               CONTINUE
               ELSE
                  C1 =         TRIG(1,K)
                  C2 =         TRIG(1,2*K)
                  S1 = ISIGN * TRIG(2,K)
                  S2 = ISIGN * TRIG(2,2*K)
                  DO 70 IJ = 1, LEN
                     AP1 = A(IJ,1,I1) + A(IJ,1,I2)
                     AP2 = A(IJ,2,I1) + A(IJ,2,I2)
                     C(IJ,1,J0) = A(IJ,1,I0) + AP1
                     C(IJ,2,J0) = A(IJ,2,I0) + AP2
                     TA1 = A(IJ,1,I0) - 0.5 * AP1
                     TA2 = A(IJ,2,I0) - 0.5 * AP2
                     AM1 = SN60 * (A(IJ,1,I1) - A(IJ,1,I2))
                     AM2 = SN60 * (A(IJ,2,I1) - A(IJ,2,I2))
                     T1 = TA1 - AM2
                     T2 = TA2 + AM1
                     C(IJ,1,J1) = C1 * T1 - S1 * T2
                     C(IJ,2,J1) = S1 * T1 + C1 * T2
                     T1 = TA1 + AM2
                     T2 = TA2 - AM1
                     C(IJ,1,J2) = C2 * T1 - S2 * T2
                     C(IJ,2,J2) = S2 * T1 + C2 * T2
 70               CONTINUE
               END IF
            END IF
            I = I + 1
            J = J + 1
 120     CONTINUE
         J = J + JUMP
 130  CONTINUE
C     
      RETURN
      END
