      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
C  MANY MACHINES YET.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS.
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
*    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* Standard C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return DBL_MIN;
*	  case 2: return DBL_MAX;
*	  case 3: return DBL_EPSILON/FLT_RADIX;
*	  case 4: return DBL_EPSILON;
*	  case 5: return log10((double)FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCRY(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
      INTEGER FUNCTION I1MACH(I)
      INTEGER I
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
C    INTEGERS HAVE FORM SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C    I1MACH( 7) = A, THE BASE.
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C    FLOATS HAVE FORM  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C               WHERE  EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, THE BASE.
C  SINGLE-PRECISION
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C  DOUBLE-PRECISION
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
      INTEGER IMACH(16), OUTPUT, SC, SMALL(2)
      SAVE IMACH, SC
      REAL RMACH
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
      INTEGER I3, J, K, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
C  INCLUDING AUTO-DOUBLE COMPILERS.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SC/987/
C
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532
     *           .AND. SMALL(2) .EQ. -448790528)
     *       .OR.     (SMALL(2) .EQ. 1117925532
     *           .AND. SMALL(1) .EQ. -448790528)) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SC = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -128
               IMACH(13) = 127
               IMACH(14) = 53
               IMACH(15) = -1024
               IMACH(16) = 1023
            ELSE
               DO 10 I3 = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(I3)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
*              *** CRAY T3E ***
               IMACH( 1) = 5
               IMACH( 2) = 6
               IMACH( 3) = 0
               IMACH( 4) = 0
               IMACH( 5) = 64
               IMACH( 6) = 8
               IMACH( 7) = 2
               IMACH( 8) = 63
               CALL I1MCR1(IMACH(9), K, 32767, 16777215, 16777215)
               IMACH(10) = 2
               IMACH(11) = 53
               IMACH(12) = -1021
               IMACH(13) = 1024
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               GO TO 35
 20            CALL I1MCR1(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
               IMACH(1) = 5
               IMACH(2) = 6
               IMACH(3) = 102
               IMACH(4) = 6
               IMACH(5) = 46
               IMACH(6) = 8
               IMACH(7) = 2
               IMACH(8) = 45
               CALL I1MCR1(IMACH(9), K, 0, 4194303, 16777215)
               IMACH(10) = 2
               IMACH(11) = 47
               IMACH(12) = -8188
               IMACH(13) = 8189
               IMACH(14) = 94
               IMACH(15) = -8141
               IMACH(16) = 8189
               GO TO 35
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
 35      SC = 987
         END IF
 9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/
     * ' statements appropriate for your machine and setting'/
     * ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
 9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/
     * ' appropriate for your machine.')
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 40
      I1MACH = IMACH(I)
      RETURN
 40   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
      STOP
* /* C source for I1MACH -- remove the * in column 1 */
* /* Note that some values may need changing. */
*#include <stdio.h>
*#include <float.h>
*#include <limits.h>
*#include <math.h>
*
*long i1mach_(long *i)
*{
*	switch(*i){
*	  case 1:  return 5;	/* standard input */
*	  case 2:  return 6;	/* standard output */
*	  case 3:  return 7;	/* standard punch */
*	  case 4:  return 0;	/* standard error */
*	  case 5:  return 32;	/* bits per integer */
*	  case 6:  return sizeof(int);
*	  case 7:  return 2;	/* base for integers */
*	  case 8:  return 31;	/* digits of integer base */
*	  case 9:  return LONG_MAX;
*	  case 10: return FLT_RADIX;
*	  case 11: return FLT_MANT_DIG;
*	  case 12: return FLT_MIN_EXP;
*	  case 13: return FLT_MAX_EXP;
*	  case 14: return DBL_MANT_DIG;
*	  case 15: return DBL_MIN_EXP;
*	  case 16: return DBL_MAX_EXP;
*	  }
*	fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
*	exit(1);return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCR1(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
