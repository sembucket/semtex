c     Routines from SPARSPAK that carry out Reverse Cuthill--McKee 
c     node ordering.  Main routine is GENRCM, at and of this file.


C----- SUBROUTINE ROOTLS
C***************************************************************           1.
C***************************************************************           2.
C********     ROOTLS ..... ROOTED LEVEL STRUCTURE      *********           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C     PURPOSE - ROOTLS GENERATES THE LEVEL STRUCTURE ROOTED                7.
C        AT THE INPUT NODE CALLED ROOT. ONLY THOSE NODES FOR               8.
C        WHICH MASK IS NONZERO WILL BE CONSIDERED.                         9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        ROOT - THE NODE AT WHICH THE LEVEL STRUCTURE IS TO               12.
C               BE ROOTED.                                                13.
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE                14.
C               GIVEN GRAPH.                                              15.
C        MASK - IS USED TO SPECIFY A SECTION SUBGRAPH. NODES              16.
C               WITH MASK(I)=0 ARE IGNORED.                               17.
C                                                                         18.
C     OUTPUT PARAMETERS -                                                 19.
C        NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE.           20.
C        (XLS, LS) - ARRAY PAIR FOR THE ROOTED LEVEL STRUCTURE.           21.
C                                                                         22.
C***************************************************************          23.
C                                                                         24.
      SUBROUTINE  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )      25.
C                                                                         26.
C***************************************************************          27.
C                                                                         28.
         INTEGER ADJNCY(1), LS(1), MASK(1), XLS(1)                        29.
         INTEGER XADJ(1), I, J, JSTOP, JSTRT, LBEGIN,                     30.
     1           CCSIZE, LVLEND, LVSIZE, NBR, NLVL,                       31.
     1           NODE, ROOT                                               32.
C                                                                         33.
C***************************************************************          34.
C                                                                         35.
C        ------------------                                               36.
C        INITIALIZATION ...                                               37.
C        ------------------                                               38.
         MASK(ROOT) = 0                                                   39.
         LS(1) = ROOT                                                     40.
         NLVL = 0                                                         41.
         LVLEND = 0                                                       42.
         CCSIZE = 1                                                       43.
C        -----------------------------------------------------            44.
C        LBEGIN IS THE POINTER TO THE BEGINNING OF THE CURRENT            45.
C        LEVEL, AND LVLEND POINTS TO THE END OF THIS LEVEL.               46.
C        -----------------------------------------------------            47.
  200    LBEGIN = LVLEND + 1                                              48.
         LVLEND = CCSIZE                                                  49.
         NLVL = NLVL + 1                                                  50.
         XLS(NLVL) = LBEGIN                                               51.
C        -------------------------------------------------                52.
C        GENERATE THE NEXT LEVEL BY FINDING ALL THE MASKED                53.
C        NEIGHBORS OF NODES IN THE CURRENT LEVEL.                         54.
C        -------------------------------------------------                55.
         DO 400 I = LBEGIN, LVLEND                                        56.
            NODE = LS(I)                                                  57.
            JSTRT = XADJ(NODE)                                            58.
            JSTOP = XADJ(NODE + 1) - 1                                    59.
            IF ( JSTOP .LT. JSTRT )  GO TO 400                            60.
               DO 300 J = JSTRT, JSTOP                                    61.
                  NBR = ADJNCY(J)                                         62.
                  IF (MASK(NBR) .EQ. 0) GO TO 300                         63.
                     CCSIZE = CCSIZE + 1                                  64.
                     LS(CCSIZE) = NBR                                     65.
                     MASK(NBR) = 0                                        66.
  300          CONTINUE                                                   67.
  400    CONTINUE                                                         68.
C        ------------------------------------------                       69.
C        COMPUTE THE CURRENT LEVEL WIDTH.                                 70.
C        IF IT IS NONZERO, GENERATE THE NEXT LEVEL.                       71.
C        ------------------------------------------                       72.
         LVSIZE = CCSIZE - LVLEND                                         73.
         IF (LVSIZE .GT. 0 ) GO TO 200                                    74.
C        -------------------------------------------------------          75.
C        RESET MASK TO ONE FOR THE NODES IN THE LEVEL STRUCTURE.          76.
C        -------------------------------------------------------          77.
         XLS(NLVL+1) = LVLEND + 1                                         78.
         DO 500 I = 1, CCSIZE                                             79.
            NODE = LS(I)                                                  80.
            MASK(NODE) = 1                                                81.
  500    CONTINUE                                                         82.
         RETURN                                                           83.
      END                                                                 84.





C----- SUBROUTINE FNROOT
C***************************************************************           1.
C***************************************************************           2.
C*******     FNROOT ..... FIND PSEUDO-PERIPHERAL NODE    *******           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C    PURPOSE - FNROOT IMPLEMENTS A MODIFIED VERSION OF THE                 7.
C       SCHEME BY GIBBS, POOLE, AND STOCKMEYER TO FIND PSEUDO-             8.
C       PERIPHERAL NODES.  IT DETERMINES SUCH A NODE FOR THE               9.
C       SECTION SUBGRAPH SPECIFIED BY MASK AND ROOT.                      10.
C                                                                         11.
C    INPUT PARAMETERS -                                                   12.
C       (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE GRAPH.          13.
C       MASK - SPECIFIES A SECTION SUBGRAPH. NODES FOR WHICH              14.
C              MASK IS ZERO ARE IGNORED BY FNROOT.                        15.
C                                                                         16.
C    UPDATED PARAMETER -                                                  17.
C       ROOT - ON INPUT, IT (ALONG WITH MASK) DEFINES THE                 18.
C              COMPONENT FOR WHICH A PSEUDO-PERIPHERAL NODE IS            19.
C              TO BE FOUND. ON OUTPUT, IT IS THE NODE OBTAINED.           20.
C                                                                         21.
C    OUTPUT PARAMETERS -                                                  22.
C       NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE             23.
C              ROOTED AT THE NODE ROOT.                                   24.
C       (XLS,LS) - THE LEVEL STRUCTURE ARRAY PAIR CONTAINING              25.
C                  THE LEVEL STRUCTURE FOUND.                             26.
C                                                                         27.
C    PROGRAM SUBROUTINES -                                                28.
C       ROOTLS.                                                           29.
C                                                                         30.
C***************************************************************          31.
C                                                                         32.
      SUBROUTINE  FNROOT ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )      33.
C                                                                         34.
C***************************************************************          35.
C                                                                         36.
         INTEGER ADJNCY(1), LS(1), MASK(1), XLS(1)                        37.
         INTEGER XADJ(1), CCSIZE, J, JSTRT, K, KSTOP, KSTRT,              38.
     1           MINDEG, NABOR, NDEG, NLVL, NODE, NUNLVL,                 39.
     1           ROOT                                                     40.
C                                                                         41.
C***************************************************************          42.
C                                                                         43.
C        ---------------------------------------------                    44.
C        DETERMINE THE LEVEL STRUCTURE ROOTED AT ROOT.                    45.
C        ---------------------------------------------                    46.
         CALL  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )         47.
         CCSIZE = XLS(NLVL+1) - 1                                         48.
         IF ( NLVL .EQ. 1 .OR. NLVL .EQ. CCSIZE ) RETURN                  49.
C        ----------------------------------------------------             50.
C        PICK A NODE WITH MINIMUM DEGREE FROM THE LAST LEVEL.             51.
C        ----------------------------------------------------             52.
  100    JSTRT = XLS(NLVL)                                                53.
         MINDEG = CCSIZE                                                  54.
         ROOT = LS(JSTRT)                                                 55.
         IF ( CCSIZE .EQ. JSTRT )  GO TO 400                              56.
            DO 300 J = JSTRT, CCSIZE                                      57.
               NODE = LS(J)                                               58.
               NDEG = 0                                                   59.
               KSTRT = XADJ(NODE)                                         60.
               KSTOP = XADJ(NODE+1) - 1                                   61.
               DO 200 K = KSTRT, KSTOP                                    62.
                  NABOR = ADJNCY(K)                                       63.
                  IF ( MASK(NABOR) .GT. 0 )  NDEG = NDEG + 1              64.
  200          CONTINUE                                                   65.
               IF ( NDEG .GE. MINDEG ) GO TO 300                          66.
                  ROOT = NODE                                             67.
                  MINDEG = NDEG                                           68.
  300       CONTINUE                                                      69.
C        ----------------------------------------                         70.
C        AND GENERATE ITS ROOTED LEVEL STRUCTURE.                         71.
C        ----------------------------------------                         72.
  400    CALL  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NUNLVL, XLS, LS )       73.
         IF (NUNLVL .LE. NLVL)  RETURN                                    74.
            NLVL = NUNLVL                                                 75.
            IF ( NLVL .LT. CCSIZE )  GO TO 100                            76.
            RETURN                                                        77.
      END                                                                 78.





C----- SUBROUTINE DEGREE
C***************************************************************           1.
C***************************************************************           2.
C********     DEGREE ..... DEGREE IN MASKED COMPONENT   ********           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C     PURPOSE - THIS ROUTINE COMPUTES THE DEGREES OF THE NODES             7.
C        IN THE CONNECTED COMPONENT SPECIFIED BY MASK AND ROOT.            8.
C        NODES FOR WHICH MASK IS ZERO ARE IGNORED.                         9.
C                                                                         10.
C     INPUT PARAMETER -                                                   11.
C        ROOT - IS THE INPUT NODE THAT DEFINES THE COMPONENT.             12.
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR.                       13.
C        MASK - SPECIFIES A SECTION SUBGRAPH.                             14.
C                                                                         15.
C     OUTPUT PARAMETERS -                                                 16.
C        DEG - ARRAY CONTAINING THE DEGREES OF THE NODES IN               17.
C              THE COMPONENT.                                             18.
C        CCSIZE-SIZE OF THE COMPONENT SPECIFED BY MASK AND ROOT           19.
C                                                                         20.
C     WORKING PARAMETER -                                                 21.
C        LS - A TEMPORARY VECTOR USED TO STORE THE NODES OF THE           22.
C               COMPONENT LEVEL BY LEVEL.                                 23.
C                                                                         24.
C***************************************************************          25.
C                                                                         26.
      SUBROUTINE  DEGREE ( ROOT, XADJ, ADJNCY, MASK,                      27.
     1                     DEG, CCSIZE, LS )                              28.
C                                                                         29.
C***************************************************************          30.
C                                                                         31.
         INTEGER ADJNCY(1), DEG(1), LS(1), MASK(1)                        32.
         INTEGER XADJ(1), CCSIZE, I, IDEG, J, JSTOP, JSTRT,               33.
     1           LBEGIN, LVLEND, LVSIZE, NBR, NODE, ROOT                  34.
C                                                                         35.
C***************************************************************          36.
C                                                                         37.
C        -------------------------------------------------                38.
C        INITIALIZATION ...                                               39.
C        THE ARRAY XADJ IS USED AS A TEMPORARY MARKER TO                  40.
C        INDICATE WHICH NODES HAVE BEEN CONSIDERED SO FAR.                41.
C        -------------------------------------------------                42.
         LS(1) = ROOT                                                     43.
         XADJ(ROOT) = -XADJ(ROOT)                                         44.
         LVLEND = 0                                                       45.
         CCSIZE = 1                                                       46.
C        -----------------------------------------------------            47.
C        LBEGIN IS THE POINTER TO THE BEGINNING OF THE CURRENT            48.
C        LEVEL, AND LVLEND POINTS TO THE END OF THIS LEVEL.               49.
C        -----------------------------------------------------            50.
  100    LBEGIN = LVLEND + 1                                              51.
         LVLEND = CCSIZE                                                  52.
C        -----------------------------------------------                  53.
C        FIND THE DEGREES OF NODES IN THE CURRENT LEVEL,                  54.
C        AND AT THE SAME TIME, GENERATE THE NEXT LEVEL.                   55.
C        -----------------------------------------------                  56.
         DO 400 I = LBEGIN, LVLEND                                        57.
            NODE = LS(I)                                                  58.
            JSTRT = -XADJ(NODE)                                           59.
            JSTOP = IABS(XADJ(NODE + 1)) - 1                              60.
            IDEG = 0                                                      61.
            IF ( JSTOP .LT. JSTRT ) GO TO 300                             62.
               DO 200 J = JSTRT, JSTOP                                    63.
                  NBR = ADJNCY(J)                                         64.
                  IF ( MASK(NBR) .EQ. 0 )  GO TO  200                     65.
                     IDEG = IDEG + 1                                      66.
                     IF ( XADJ(NBR) .LT. 0 ) GO TO 200                    67.
                        XADJ(NBR) = -XADJ(NBR)                            68.
                        CCSIZE = CCSIZE + 1                               69.
                        LS(CCSIZE) = NBR                                  70.
  200          CONTINUE                                                   71.
  300       DEG(NODE) = IDEG                                              72.
  400    CONTINUE                                                         73.
C        ------------------------------------------                       74.
C        COMPUTE THE CURRENT LEVEL WIDTH.                                 75.
C        IF IT IS NONZERO , GENERATE ANOTHER LEVEL.                       76.
C        ------------------------------------------                       77.
         LVSIZE = CCSIZE - LVLEND                                         78.
         IF ( LVSIZE .GT. 0 ) GO TO 100                                   79.
C        ------------------------------------------                       80.
C        RESET XADJ TO ITS CORRECT SIGN AND RETURN.                       81.
C        ------------------------------------------                       82.
         DO 500 I = 1, CCSIZE                                             83.
            NODE = LS(I)                                                  84.
            XADJ(NODE) = -XADJ(NODE)                                      85.
  500    CONTINUE                                                         86.
         RETURN                                                           87.
      END                                                                 88.





C----- SUBROUTINE RCM
C***************************************************************           1.
C***************************************************************           2.
C********     RCM ..... REVERSE CUTHILL-MCKEE ORDERING   *******           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C     PURPOSE - RCM NUMBERS A CONNECTED COMPONENT SPECIFIED BY             7.
C        MASK AND ROOT, USING THE RCM ALGORITHM.                           8.
C        THE NUMBERING IS TO BE STARTED AT THE NODE ROOT.                  9.
C                                                                         10.
C     INPUT PARAMETERS -                                                  11.
C        ROOT - IS THE NODE THAT DEFINES THE CONNECTED                    12.
C               COMPONENT AND IT IS USED AS THE STARTING                  13.
C               NODE FOR THE RCM ORDERING.                                14.
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR                    15.
C               THE GRAPH.                                                16.
C                                                                         17.
C     UPDATED PARAMETERS -                                                18.
C        MASK - ONLY THOSE NODES WITH NONZERO INPUT MASK                  19.
C               VALUES ARE CONSIDERED BY THE ROUTINE.  THE                20.
C               NODES NUMBERED BY RCM WILL HAVE THEIR                     21.
C               MASK VALUES SET TO ZERO.                                  22.
C                                                                         23.
C     OUTPUT PARAMETERS -                                                 24.
C        PERM - WILL CONTAIN THE RCM ORDERING.                            25.
C        CCSIZE - IS THE SIZE OF THE CONNECTED COMPONENT                  26.
C               THAT HAS BEEN NUMBERED BY RCM.                            27.
C                                                                         28.
C     WORKING PARAMETER -                                                 29.
C        DEG - IS A TEMPORARY VECTOR USED TO HOLD THE DEGREE              30.
C               OF THE NODES IN THE SECTION GRAPH SPECIFIED               31.
C               BY MASK AND ROOT.                                         32.
C                                                                         33.
C     PROGRAM SUBROUTINES -                                               34.
C        DEGREE.                                                          35.
C                                                                         36.
C***************************************************************          37.
C                                                                         38.
      SUBROUTINE  RCM ( ROOT, XADJ, ADJNCY, MASK,                         39.
     1                  PERM, CCSIZE, DEG )                               40.
C                                                                         41.
C***************************************************************          42.
C                                                                         43.
         INTEGER ADJNCY(1), DEG(1), MASK(1), PERM(1)                      44.
         INTEGER XADJ(1), CCSIZE, FNBR, I, J, JSTOP,                      45.
     1           JSTRT, K, L, LBEGIN, LNBR, LPERM,                        46.
     1           LVLEND, NBR, NODE, ROOT                                  47.
C                                                                         48.
C***************************************************************          49.
C                                                                         50.
C        -------------------------------------                            51.
C        FIND THE DEGREES OF THE NODES IN THE                             52.
C        COMPONENT SPECIFIED BY MASK AND ROOT.                            53.
C        -------------------------------------                            54.
         CALL  DEGREE ( ROOT, XADJ, ADJNCY, MASK, DEG,                    55.
     1                  CCSIZE, PERM )                                    56.
         MASK(ROOT) = 0                                                   57.
         IF ( CCSIZE .LE. 1 ) RETURN                                      58.
         LVLEND = 0                                                       59.
         LNBR = 1                                                         60.
C        --------------------------------------------                     61.
C        LBEGIN AND LVLEND POINT TO THE BEGINNING AND                     62.
C        THE END OF THE CURRENT LEVEL RESPECTIVELY.                       63.
C        --------------------------------------------                     64.
  100    LBEGIN = LVLEND + 1                                              65.
         LVLEND = LNBR                                                    66.
         DO 600 I = LBEGIN, LVLEND                                        67.
C           ----------------------------------                            68.
C           FOR EACH NODE IN CURRENT LEVEL ...                            69.
C           ----------------------------------                            70.
            NODE = PERM(I)                                                71.
            JSTRT = XADJ(NODE)                                            72.
            JSTOP = XADJ(NODE+1) - 1                                      73.
C           ------------------------------------------------              74.
C           FIND THE UNNUMBERED NEIGHBORS OF NODE.                        75.
C           FNBR AND LNBR POINT TO THE FIRST AND LAST                     76.
C           UNNUMBERED NEIGHBORS RESPECTIVELY OF THE CURRENT              77.
C           NODE IN PERM.                                                 78.
C           ------------------------------------------------              79.
            FNBR = LNBR + 1                                               80.
            DO 200 J = JSTRT, JSTOP                                       81.
               NBR = ADJNCY(J)                                            82.
               IF ( MASK(NBR) .EQ. 0 )  GO TO 200                         83.
                  LNBR = LNBR + 1                                         84.
                  MASK(NBR) = 0                                           85.
                  PERM(LNBR) = NBR                                        86.
  200       CONTINUE                                                      87.
            IF ( FNBR .GE. LNBR )  GO TO 600                              88.
C              ------------------------------------------                 89.
C              SORT THE NEIGHBORS OF NODE IN INCREASING                   90.
C              ORDER BY DEGREE. LINEAR INSERTION IS USED.                 91.
C              ------------------------------------------                 92.
               K = FNBR                                                   93.
  300          L = K                                                      94.
                  K = K + 1                                               95.
                  NBR = PERM(K)                                           96.
  400             IF ( L .LT. FNBR )  GO TO 500                           97.
                     LPERM = PERM(L)                                      98.
                     IF ( DEG(LPERM) .LE. DEG(NBR) )  GO TO 500           99.
                        PERM(L+1) = LPERM                                100.
                        L = L - 1                                        101.
                        GO TO 400                                        102.
  500             PERM(L+1) = NBR                                        103.
                  IF ( K .LT. LNBR )  GO TO 300                          104.
  600    CONTINUE                                                        105.
         IF (LNBR .GT. LVLEND) GO TO 100                                 106.
C        ---------------------------------------                         107.
C        WE NOW HAVE THE CUTHILL MCKEE ORDERING.                         108.
C        REVERSE IT BELOW ...                                            109.
C        ---------------------------------------                         110.
         K = CCSIZE/2                                                    111.
         L = CCSIZE                                                      112.
         DO 700 I = 1, K                                                 113.
            LPERM = PERM(L)                                              114.
            PERM(L) = PERM(I)                                            115.
            PERM(I) = LPERM                                              116.
            L = L - 1                                                    117.
  700    CONTINUE                                                        118.
         RETURN                                                          119.
      END                                                                120.





C----- SUBROUTINE GENRCM
C***************************************************************           1.
C***************************************************************           2.
C********   GENRCM ..... GENERAL REVERSE CUTHILL MCKEE   *******           3.
C***************************************************************           4.
C***************************************************************           5.
C                                                                          6.
C     PURPOSE - GENRCM FINDS THE REVERSE CUTHILL-MCKEE                     7.
C        ORDERING FOR A GENERAL GRAPH. FOR EACH CONNECTED                  8.
C        COMPONENT IN THE GRAPH, GENRCM OBTAINS THE ORDERING               9.
C        BY CALLING THE SUBROUTINE RCM.                                   10.
C                                                                         11.
C     INPUT PARAMETERS -                                                  12.
C        NEQNS - NUMBER OF EQUATIONS                                      13.
C        (XADJ, ADJNCY) - ARRAY PAIR CONTAINING THE ADJACENCY             14.
C               STRUCTURE OF THE GRAPH OF THE MATRIX.                     15.
C                                                                         16.
C     OUTPUT PARAMETER -                                                  17.
C        PERM - VECTOR THAT CONTAINS THE RCM ORDERING.                    18.
C                                                                         19.
C     WORKING PARAMETERS -                                                20.
C        MASK - IS USED TO MARK VARIABLES THAT HAVE BEEN                  21.
C               NUMBERED DURING THE ORDERING PROCESS. IT IS               22.
C               INITIALIZED TO 1, AND SET TO ZERO AS EACH NODE            23.
C               IS NUMBERED.                                              24.
C        XLS - THE INDEX VECTOR FOR A LEVEL STRUCTURE.  THE               25.
C               LEVEL STRUCTURE IS STORED IN THE CURRENTLY                26.
C               UNUSED SPACES IN THE PERMUTATION VECTOR PERM.             27.
C                                                                         28.
C     PROGRAM SUBROUTINES -                                               29.
C        FNROOT, RCM.                                                     30.
C                                                                         31.
C***************************************************************          32.
C                                                                         33.
      SUBROUTINE  GENRCM ( NEQNS, XADJ, ADJNCY, PERM, MASK, XLS )         34.
C                                                                         35.
C***************************************************************          36.
C                                                                         37.
         INTEGER ADJNCY(1), MASK(1), PERM(1), XLS(1)                      38.
         INTEGER XADJ(1), CCSIZE, I, NEQNS, NLVL,                         39.
     1           NUM, ROOT                                                40.
C                                                                         41.
C***************************************************************          42.
C                                                                         43.
         DO 100 I = 1, NEQNS                                              44.
            MASK(I) = 1                                                   45.
  100    CONTINUE                                                         46.
         NUM = 1                                                          47.
         DO 200 I = 1, NEQNS                                              48.
C           ---------------------------------------                       49.
C           FOR EACH MASKED CONNECTED COMPONENT ...                       50.
C           ---------------------------------------                       51.
            IF (MASK(I) .EQ. 0) GO TO 200                                 52.
               ROOT = I                                                   53.
C              -----------------------------------------                  54.
C              FIRST FIND A PSEUDO-PERIPHERAL NODE ROOT.                  55.
C              NOTE THAT THE LEVEL STRUCTURE FOUND BY                     56.
C              FNROOT IS STORED STARTING AT PERM(NUM).                    57.
C              THEN RCM IS CALLED TO ORDER THE COMPONENT                  58.
C              USING ROOT AS THE STARTING NODE.                           59.
C              -----------------------------------------                  60.
               CALL  FNROOT ( ROOT, XADJ, ADJNCY, MASK,                   61.
     1                        NLVL, XLS, PERM(NUM) )                      62.
               CALL     RCM ( ROOT, XADJ, ADJNCY, MASK,                   63.
     1                        PERM(NUM), CCSIZE, XLS )                    64.
               NUM = NUM + CCSIZE                                         65.
               IF (NUM .GT. NEQNS) RETURN                                 66.
  200    CONTINUE                                                         67.
         RETURN                                                           68.
      END                                                                 69.
