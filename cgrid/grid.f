     PROGRAM BOFITS

C A PROGRAM FOR THE NUMBERICAL GENERATION OF NON-ORTHOGONAL
C COORDINATES IN A SIMPLY-CONNECTED DOMAIN
C THIS PROGRAM HAS BEEN DEVELOPED BY DR. Z.U.A. WARSI, VISITING PROFESSOR
C AT VKI FROM MISSISSIPPI STATE UNIVERSITY DEPT. AEROSPACE ENG.
C MISS STATE, MS 39762, USA
C NOTE: BEFORE RUNNING THE PROG. BE SURE TO CHECK THE VALUES OF
C KMAS, IM BOTH IN THE MAIN PROG. AND THE SUBROUTINE
C ALSO SET THE VALUE OF LMAS, KM, TYPE, CK, AND MOD.
C $Id$

C *************************************************************************

          PARAMETER (KMAX=85, IM=25)
          COMMON /XY/ XB(KMAX,IM), YB(KMAX,IM), X(KMAX,IM), Y(KMAX,IM)
     1,KM,KMM,LMAX
          COMMON /ZUWA/G11(KMAX,IM), G13(KMAX,IM), G33(KMAX,IM)
     1,GG(KMAX,IM), CJ2(KMAX,IM)
          COMMON /ZUWB/ AA(KMAX,IM), BB(KMAX,IM), CC(KMAX,IM)

          DIMENSION A(KMAX,IM), B(KMAX,IM), C(KMAX,IM), RT(KMAX, IM)
     1, ST(KMAX,IM), XP(KMAX,IM), XC(KMAX,IM), YP(KMAX,IM),
     2YC(KMAX,IM), D(KMAX,IM), H(KMAX,IM), E(KMAX,IM),
     3T(KMAX,IM), P1(KMAX,IM), P2(KMAX,IM), P3(KMAX,IM),
     4Q1(KMAX,IM), Q2(KMAX,IM),Q3(KMAX,IM), PBAR(KMAX,IM),
     5QBAR(KMAX,IM), DIFX(KMAX,IM), DIFY(KMAX,IM)

C MOD STANDS FOR VARIOUS MODIFICATIONS IN THE INPUT DATA.

C MOD=1,2,3 ARE MERELY EXAMPLES.

C CK IS THE COORDINATE CONTROL PARAMETER FOR THE I-COORDINATE

C CK=1.0 IMPLIES NO CONTRACTION.

C SET CK=1.0 IF MOD=1 OR 2.

C TYPE =1.  STANDS FOR THOSE BOUNDARY COORDINATES IN WHICH EACH

C GIVEN BOUNDARY BELONGS TO A DISTINCT FAMILY OF COORDINATE LINES.

C TYPE=2.  STANDS FOR THOSE BOUNDARY COORDINATES IN WHICH ONE GIVEN







C BOUNDARY IS COMPOSED OF TWO DIFFERENT COORDINATE SPECIES.

         TYPE=1.
         MOD=4
         CK=1.1
         PI=3.14159264
         LMAX=21
         KM=81
         KMM=KM-1
         KPM=KM+1
         LMM=LMAX-1
         LPM=LMAX+1

C SET LM=LMAX FOR MOD=1 OR 2, AND ALSO FOR THE CASE WHEN

C CONTRACTION IS DESIRED ONLY NEAR I=1.  FOR MOD=3 OR 4, IF

C CONTRACTION IS DESIRED BOTH NEAR I=1 AND ALSO NEAR I=LMAX, THEN

C CHOOSE A VALUE OF LM LESS THAN LMAX.

         LM=LPM/2
         LM=LMAX
         ITM=40
         ELM=.1
         W=1.5
         ANG=PI
         IF(MOD.EQ.1) GO TO 6500
         IF(MOD.EQ.2) GO TO 7000
         IF(MOD.EQ.3) GO TO 8500
         IF(MOD.EQ.4) GO TO 9500
  6500   CONTINUE
         K1=9
         K2=26
         DA=PI/(K2-K1)
         DO 108 K=1, K1
         XB(K,1)=(K+1. -2.*K1)/(K1-1.)
         YB(K,1)=0
  108    CONTINUE
         DO 111 K=K1, K2
         XB(K,1)=COS(ANG)
         YB(K,1)=SIN(ANG)
         ANG=ANG-DA
  111    CONTINUE
         DO 112 K=K2, KM
         XB(K,1)=(K+KM-2.*K2)/(KM-K2)
         YB(K,1)=0
  112    CONTINUE
         DO 113 K=1, KM
         XB(K,LMAX)=(4.*K-2.*KM-2.)/KMM

         YB(K,LMAX)=2.
  113    CONTINUE
         DO 114 I=1,LMAX






         XB(1,I)=-2
         YB(1,I)=2.*(I-1)/LMM
         XB(KM,I)=2.
         YB(KM,I)=2.*(I-1)/LMM
  114    CONTINUE
         IF(MOD.EQ.1) GO TO 8000
  7000   CONTINUE
         K1=9
         K2=26
         AL=2.
         BT=1.
         EP=.1
         GM=1.
         DO 7001 K=1,K1-1
         XB(K,1)=((AL-BT-EP)*K-AL*(K1-1.)+BT+EP)/(K1-2.)
         YB(K,1)=0.
  7001   CONTINUE
         XB(K1,1)=-BT
         YB(K1,1)=GM
         DO 7002 K=K+1,K2-1
         XB(K,1)=(BT-EP)*(2.*K-K1-K2)/(K2-K1-2.)
         YB(K,1)=GM
  7002   CONTINUE
         XB(K2,1)=BT
         YB(K2,1)=GM
         DO 7003 K=K2+1,KM
         XB(K,1)=((AL-BT-EP)*K+(BT+EP)*KM-(K2+1.)*AL)
     1/  (KM-K2-1.)
         YB((K,1)=0.
  7003   CONTINUE
         DO 7004 K=1, KM
         XB(K,LMAX)=(2.*AL*K-AL*KM-AL)/KMM
         YB(K,LMAX)=AL
  7004   CONTINUE
         DO 7005 I=1,LMAX
         XB(1,I)=-AL
         XB(KM,I)=AL
         YB(1,I)=AL*(I-1.)/LMM
         YB(KM,I)=AL*(I-1.)/LMM
  7005   CONTINUE
         IF(MOD.EQ.2) GO TO 8000
  8500   CONTINUE
         DA=PI/KMM
         AI=1.
         BI=1.
         AO=3.
         BO=2.
         AS=AO**2
         BS=BO**2
         DO 351 K=1,KM
         XB(K,1)=AI*COS(ANG)
         YB(K,1)=BI*SIN(ANG)
         WD=SQRT(AS*(SIN(ANG)**2)+BS*(COS(ANG)**2))
         XB(K,LMAX)=AO*BO*COS(ANG)/WD






         YB(K,LMAX)=AO*BO*SIN(ANG)/WD
C        XB(K,LMAX)=AO*COS(ANG)
C        YB(K,LMAX)=BO*SIN(ANG)
         ANG=ANG-DA
  351    CONTINUE
         F2=FLOAT(LM-1)
         F4=FLOAT(LMAX-LM)
         F6=FLOAT(LMM)
         IF(LM.EQ.LMAX) GO TO 358
         DO 349 I=1,LMAX
         F1=FLOAT(I-1)
         F3=FLOAT(LMAX-I)
         G1=F1/F2
         G2=F3/F4
         IF(I.LE.LM) GO TO 352
         IF(I.GT.LM) GO TO 353
  352    CONTINUE
         XB(1,I)=-AI-G1*((CK)**(I-LM))
         GO TO 354
  353    CONTINUE
         XB(1,I)=-AO+G2*((CK)**(LM-I))
  354    CONTINUE
         YB(1,I)=0.
  349    CONTINUE
         DO 350 I=1,LMAX
         F1=FLOAT(I-1)
         F3=FLOAT(LMAX-I)
         G1=F1/F2
         G2=F3/F4
         IF(I.LE.LM) GO TO 355
         IF(I.GT.LM) TO TO 356
  355    CONTINUE
         XB(KM,I)=AI+G1*((CK)**(I-LM))
         GO TO 357
  356    CONTINUE
         XB(KM,I)=AO-G2*((CK)**(LM-I))
  357    CONTINUE
         YB(KM,I)=0.
  350    CONTINUE
         GO TO 8000
  358    CONTINUE
         DO 359 I=1,LMAX
         F1=FLOAT(I-1)
         G1=F1/F6
         XB(1,I)=-AI-(AO-AI)*G1*((CK)**(I-LMAX))
         YB(1,I)=0.
  359    CONTINUE
         DO 360 I=1,LMAX
         F1=FLOAT(I-1)
         G1=F1/F6
         XB(KM,I)=AI+(AO-AI)*G1*((CK)**(I-LMAX))
         YB(KM,I)=0.
  360    CONTINUE
         IF(MOD.EQ.3) GO TO 8000






  9500   CONTINUE
         READ(20,8001)((XB(K,1),YB(K,1)),K=1,KM)
         READ(20,8001)((XB(K,LMAX),YB(K,LMAX)),K=1,KM)
         READ(20,8001)((XB(1,I),YB(1,I)),I=1.LMAX)
         READ(20,8001)((XB(KM,I),YB(KM,I)),I=1,LMAX)
  8001   FORMAT(2F10.5)
         IF(TYPE.EQ.1.) GO TO 8000
         F3=FLOAT(LM-1)
         F4=FLOAT(LMAX-LM)
         F5=FLOAT(LMM)
         IF(LM.EQ.LMAX) GO TO 9910
         DO 9900 I=1,LMAX
         F1=FLOAT(I-1)
         F2=FLOAT(LMAX-I)
         G1=F1/F3
         G2=F2/F4
         IF(I.LE.LM) GO TO 9901
         IF(I.GT.LM) GO TO 9902
  9901   CONTINUE
         XB(1,I)=XB(1,1)+(XB(1,LM)-XB(1,1))*
     1   G1*(CK**(I-LM))
         YB(1,I)=YB(1,1)+(YB(1,LM)-YB(1,1))*
     1   G1*(CK**(I-LM))
         XB(KM,I)=XB(KM,1)+(XB(KM,LM)-XB(KM,1))*
         GI*(CK**(I-LM))
         YB(KM,I)=YB(KM,1)+(YB(KM,LM)-YB(KM,1))*
     1   G1*(CK**(I-LM))
         GO TO 9905
  9902   CONTINUE
         XB(1,I)=XB(1,LMAX)-(XB(1,LMAX)-XB(1,LM))*
     1   G2*(CK**(LM-I))
         YB(1,I)=YB(1,LMAX)-(YB(1,LMAX)-YB(1,LM))*
     1   G2*(CK**(LM-I))
         XB(KM,I)=XB(KM,LMAX)-(XB(KM,LMAX)-XB(KM,LM))*
     1   G2*(CK**(LM-I))
         YB(KM,I)=YB(KM,LMAX)-(YB(KM,LMAX)-YB(KM,LM))*
     1   G2*(CK**(LM-I))
  9905   CONTINUE
  9900   CONTINUE
         GO TO 8000
  9910   CONTINUE
         DO 9909 I=1,LMAX
         F1=FLOAT(I-1)
         G1=F1/F5
         XB(1,I0=XB(1,1)+(XB(1,LMAX)-XB(1,1))*
     1   G1*(CK**)(I-LMAX))
         YB(1,I)=YB(1,1)+(YB(1,LMAX)-YB(1,1))*
     1   G1*(CK**)(I-LMAX))
         XB(KM,I)=XB(KM,1)+(XB(KM,LMAX)-XB(KM,1))*
     1   G1*(CK**(1-LMAX))
         YB(KM,I)=YB(KM,1)+(YB(KM,LMAX)-YB(KM,1))*
     1   G1*(CK**(I-LMAX))
  9909   CONTINUE
  8000   CONTINUE






         DO 28 K=1, KM
         X(K,1)=XB(K,1)
         Y(K,1)=YB(K,1)
         X(K,LMAX)=XB(K,LMAX)
         Y(K,LMAX)=YB(K,LMAX)
  28     CONTINUE
         DO 29 I=1,LMAX
         X(1,I)=XB(1,I)
         Y(1,I0=YB(1,I)
         X(KM,I)=XB(KM,I)
         Y(KM,I)=YB(KM,I)
  29     CONTINUE
         DO 116 K=1,KM
         WRITE(1,109) K,XB(K,1),YB(K,1),XB(K,LMAX),YB(K,LMAX)
  109    FORMAT(  INPUT DATA  ,14,4F12.4)
  116    CONTINUE
         DO 117 I=1,LMAX
         WRITE(1,118) I,XB(1,I),YB(1,I),XB(KM,I),YB(KM,I)
  118    FORMAT(  INPUT DATA  ,14,4F12.4)
  117    CONTINUE

C INITIAL GUESS:
         C5=FLOAT(LMM)
         C6=FLOAT(KMM)
         D5=FLOAT(LMM*KMM)
         DO 4 K=2,KMM
         DO 3 I=2,LMM
         C1=FLOAT(LMAX-I)
         C2=FLOAT(I-1)
         C3=FLOAT(KM-K)
         C4=FLOAT(K-1)
         D1=FLOAT((LMAX-I)*(KM-K))
         D2=FLOAT((LMAX-I)*(K-1))
         D3=FLOAT((I-1)*(KM-K))
         D4=FLOAT((I-1)*(K-1))
         E1=(C1*X(K,1)+C2*X(K,LMAX))/C5+
     1   (C3*X(1,I)+C4*X(KM,I))/C6
         E2=(D1*X(1,1)+D2*X(KM,1)+D3*X(1,LMAX)+
     1   D4*X(KM,LMAX))/D5
         X(K,I)=E1-E2
         H1=(C1*Y(K,1)+C2*Y(K,LMAX))/C5+
     1   (C3*Y(1,I)+C4*Y(KM,I))/C6
         H2=(D1*Y(1,1)+D2*Y(KM,1)+D3*Y(1,LMAX)+
     1   DF*Y(KM,LMAX))/D5
         Y(K,I)=H1-H2
  3      CONTINUE
  4      CONTINUE
         DO 110 IJK=1, ITM
         CHECK=0.0
         DO 100 II=2,LMAX-1
         CALL GCOEF(II)

C    INSERT THE EXPRESSIONS FOR P1,P2,P3,Q1,Q2,Q3 BELOW.







         DO 777 K=1,KM
         P1(K,II)=0.
         P2(K,II)=0.
         P3(K,II)=0.
         Q1(K,II)=0.
         Q2(K,II)=0
         IF(II.LE.LM) GO TO 9000
         IF(II.GT.LM) GO TO 9001
  9000   CONTINUE
         Q3(K,II)=-2.+(II-1)*ALOG(CK))*ALOG(CK)
     1   /(1.+(II-1)*ALOG(CK))
         GO TO 9003
  9001   CONTINUE
         Q3(K,II)=(2.+(LMAX-II)*ALOG(CK))*ALOG(CK)
     1   /(1.+(LMAX-II)*ALOG(CK))
  9003   CONTINUE
  777    CONTINUE
         DO 504 K=1,KM
         PBAR(K,II)-P1(K,II)*AA(K,II)-2.0*P2(K,II)*
     18B(K,II)+P3(K,II)*CC(K,II)
         QBAR(K,II)=Q1(K,11*AA(K,II)-2.0*Q2(K,II)*
     18B(K,II)+Q3(K,II)*CC(K,II)
  504    CONTINUE
  620    CONTINUE
         DO 6 K=1,KM
         A(K,II)=AA(K,II)+0.5*PBAR(K,II)
         B(K,II)=-2.*(AA(K,II)+CC(K,II))
         C(K,II)=AA(K,II)-0.5*PBAR(K,II)
  6      CONTINUE
         RT(1,II)=-(CC(1,II)+0.5*QBAR(1,II))*X(1,II+1)-(CC(1,II)-
     10.5*QBAR(1,II*X(1,II-1)+0.5*BB(1,II)*(X(2,II+1)-
     1X(2,II-1)-X(KMM,II+1)+X(KMM,II-1))
         ST(1,II)=-(CC(1,II)+0.5*QBAR(1,II))*Y(1,II+1)-(CC(1,II
     1)-0.5*QBAR(1,II))*Y(1,II-1)+BB(1,II)*(Y(2,II+1)-Y(2,II-1)-
     2Y(KMM,II+1)+Y(KMM,II-1))*0.5
         DO 8 K=2,KMM
         RT(K,II)=-(CC(K,II)+0.5*QBAR(K,II))*X(K,II+1)-(CC(K,II)
     1-0.5*QBAR(1,II))*X(K,II-1)+0.5*BB(K,II)*(X(K+1,II+1)-
     2X(K+1,II-1)-X(K-1,II+1)+X(K-1,II-1))
         ST(K,II)=-(CC(K,II)+0.5*QBAR(K,II))*Y(K,II+1)-(CC(K,II)
     1-0.5*QBAR(K,II))*Y(K,II-1)+BB(K,II)*(Y(K+1,II+1)-Y(K+1,II-1)-
     2Y(K-1,II+1)+Y(K-1,II-1))*0.5
  8  CONTINUE
C XP,YP  DENOTE THE PRESENTLY AVAILABLE ITERATIVE VALUES
C XC,YC  DENOTE THE PREVIOUS ITERATIVE VALUES
         DO 5 K=1,KM
         XP(K,II)=X(K,II)
         YP(K,II)=Y(K,II)
         XC(K,II)-X(K,II)
         YC(K,II)=Y(K,II)
  5      CONTINUE
         D(1,II)=0
         H(1,II)=0
         E(1,II)=XP(1,II)






         T(1,II)=YP(1,II)
         DO 10 K=2,KMM
         D(K,II)=-A(K,II)/(B(K,II)+C(K,II)*D(K-1,II))
         H(K,II)=-A(K,II)/(B(K,II)+C(K,II)*H(K-1,II))
         E(K,II)=(RT(K,II)-C(K,II)*E(K-1,II))/(B(K,II)+C(K,II)*D(K-1,II))
         T(K,II)=(ST(K,II)-C(K,II)*T(K-1,II))/(B(K,II)+C(K,II)*H(K-1,II))
  10     CONTINUE
         X(1,II)=W*XP(1,II)+(1.-W)*XC(1,II)
         Y(1,II)=W*YP(1,II)+(1.-W)*YC(1,II)
         DO 20 K=KMM,2,-1
         XP(K,II)=D(K,II)*XP(K+1,II)+E(K,II)
         YP(K,II)=H(K,II)*YP(K+1,II)+T(K,II)
         X(K,II)=W*XP(K,II)+(1.-W)*XC(K,II)
         Y(K,II)=W*YP(K,II)+(1.-W)*YC(K,II)
  20     CONTINUE
  100    CONTINUE
         ERX=0.
         ERY=0
         DO 30 I=2,LMAX-1
         DO 30 K=1,KMM
         DIFX(K,I)=X(K,I)-XC(K,I)
         DIFY(K,I)=Y(K,I)-YC(K,I)
         CHECK-CHECK+ABS(DIFX(K,I))+ABS(DIFY(K,I))
         ERX=AMAX1(ERX,DIFX(K,I))
         ERY=AMAX1(ERY,DIFY(K,I))
  30     CONTINUE
         KKK=IJK
         WRITE(1,609)KKK,CHECK,ERX,ERY
         IF(CHECK.LT.ELM) GO TO 120
  110    CONTINUE
  120    WRITE(1,609)KKK,CHECK,ERX,ERY
  609    FORMAT('ITERATION NUMBER' ,15,' TOTAL AREAS ARE'
     1,F10.6,/, MAX POINT ERROR EX= ,F10.7,2X'EY=',F10.7)
         DO 1002 I=1,LMAX
         WRITE(1,1001)I,G13(3,I)
  1001   FORMAT('G13 DATA', 15,F12.5)
  1002   CONTINUE
         WRITE(1,900)
         DO 150 I=1,LMAX
         DO 160 K=1,KM
         WRITE(1,910)K,I,X(K,I),Y(K,I),AA(K,I),CC(K,I),
     1   G13(K,I),CJ2(K,I)
  160    CONTINUE
  150    CONTINUE
  900    FORMAT(1H, ' K I ',4X,'X','Y',FX,'AA',5X,'CC'
     1   ,5X,'G13',5X,'CJ2')
  910    FORMAT('  ',213,6F7.3)
         WRITE(10,*) ((X(K,I),Y(K,I),K=1,KM),I=1,LMAX)
         STOP

         END
         SUBROUTINE GCOEF(II)
         PARAMETER (KMAX=85,IM=25)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)






     1,KM,KMM,LMAX
         COMMON /ZUWA/ G11(KMAX,IM),G13(KMAX,IM),G33(KMAX,IM)
     1,GG(KMAX,IM),CJ2(KMAX,IM)
         COMMON /ZUWB/ AA(KMAX,IM),BB(KMAX,IM),CC(KMAX,IM)
                 DIMENSION DXXI(KMAX,IM),DYXI(KMAX,IM),
     1DXZT(KMAX,IM),DYZT(KMAX,IM),PG(KMAX,IM),PH(KMAX,IM)
         DO 70 I=II-1,II+1
         DXZT(1,I)=(-3.*X(1,I)+4.*X(2,I)-X(3,I))/2.
         DYZT(1,I)=(-3.*Y(1,I)+4.*Y(2,I)-Y(3,I))/2.
         DXZT(KM,I)=(3.*X(KM,I)-4.*X(KMM,I)+X(KMM-1,I))/2.
         DYZT(KM,I)=(3.*Y(KM,I)-4.*Y(KMM,I)+Y(KMM-1,I))/2.
         DO 5 K=2,KMM
         DXZT(K,I)=(X(K+1,I)-X(K-1,I))/2.
         DYZT(K,I)=(Y(K+1,I)-Y(K-1,I))/2.
  5      CONTINUE
         DO 35 K=1,KM
         IF(I.EQ.1) THEN
         DXXI(K,I)=(-3.*X(K,I)+4.*X(K,I+1)-X(K,I+2))/2.
         DYXI(K,I)=(-3.*Y(K,I)+4.*Y(K,I+1)-Y(K,I+2))/2.
         GO TO 34
         END IF
         IF(I.EQ.LMAX) THEN
         DXXI(K,I)=(3.*(K,I)-4.*X(K,I-1)+X(K,I-2))/2.
         DYXI(K,I)=(3.*Y(K,I)-4.*(K,I-1)+Y(K,I-2))/2.
         GO TO 34
         END IF
         DXXI(K,I)=(X(K,I+1)-X(K,I-1))/2.
         DYXI(K,I)=(Y(K,I+1)-Y(K,I-1))/2.
  34     CONTINUE
  35     CONTINUE
  70     CONTINUE
         DO 80 I=II-1,II+1
         DO 82 K=1,KM
         AA(K,I)=DXXI(K,I)**2+DYXI(K,I)**2
         G13(K,I)=DXXI(K,I)*DXZT(K,I)+DYXI(K,I)*DYZT(K,I)
         BB(K,I)=G13(K,I)
  630    CONTINUE
         CC(K,I)=DXZT(K,I)**2+DYZT(K,I)**2
         CJ2(K,I)=SQRT(AA(K,I)*CC(K,I)-BB(K,I)**2)
  82     CONTINUE
  80     CONTINUE
  635    CONTINUE
         RETURN
         END







                __________ - _______________________________
                   _______________________________________
                  (Generally Orthogonal.  Based on Ref. 13)


     PROGRAMq BOFITSOL

C A PROGRAM FOR THE NUMERICAL GENERATION OF ORTHOGONAL
C COORDINATES IN A SIMPLY-CONNECTED DOMAIN
C THIS PROGRAM HAS BEEN DEVELOPED BY DR.Z.U.A.WARSI,VISITING PROFESSOR
C AT VKI FROM MISSISSIPPI STATE UNIVERSITY DEPT. AEROSPACE ENG.
C MISS. STATE, MS 39762, USA.
         PARAMETER(KMAX=40,IM=30)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)
     1,KM,KMM,LMAX,LM,LMM,J1,CK
         COMMON /ZUWA/ G13(KMAX,IM),CJ2(KMAX,IM)
         COMMON /ZUWB/ AA(KMAX,IM),BB(KMAX,IM),CC(KMAX,IM)
         COMMON /ZUWC/ DXXI(KMAX,IM),DYXI(KMAX,IM)
     1,DXZT(KMAX,IM),DYZT(KMAX,IM),F(KMAX,IM),FINV(KMAX,IM)
         DIMENSION A(KMAX,IM),B(KMAX,IM),C(KMAX,IM),RT(KMAX,IM)
     1,ST(KMAX,IM),XP(KMAX,IM),XC(KMAX,IM),YP(KMAX,IM),
     2YC(KMAX,IM),D(KMAX,IM),H(KMAX,IM),E(KMAX,IM),
     3T(KMAX,IM),PBAR(KMAX,IM),QBAR(KMAX,IM),
     4DIFX(KMAX,IM),DIFY(KMAX,IM)

C ****************************************************************************

C NOTE:  THIS PROGRAM CAN GENERATE OROTHOGONAL COORDINATES ONLY

C WHEN THE GIVEN BOUNDARY SEGMENTS OF DIFFERENT SPECIES ARE ORTHOGONAL

C NOTE:  BEFORE RUNNING THE PROG. BE SURE TO CHECK THE VALUES OF

C KMAX,IM BOTH IN THE MAIN PROG. AND THE SUBROUTINE.

C ALSO SET THE VALUES OF LMAX,KM,TYPE,CK, AND MOD.

C MOD=1 STANDS FOR ANALYTIC INPUT DATA.

C MOD=2 STANDS FOR THE INPUT DATA THROUGH A DATA FILE.

C CK IS THE COORDINATE CONTROL PARAMETER FOR THE I-COORDINATES

C TYPE=1.  STANDS FOR THOSE COORDINATES IN WHICH EACH GIVEN BOUNDARY

C BELONGS TO A DISTINCT FAMILY OF COORINDATE LINES.

C TYPE=2 STANDS FOR THOSE BOUNDARY COORDINATES IN WHICH ONE GIVEN

C BOUNDARY IS COMPOSED OF TWO DIFFERENT COORDINATE SPECIES

         TYPE=2.
         CK=1.1
         MOD=1






         PI=3.14159264
         W=1 5
     MAKE LMAX AN ODD INTEGER
         LMAX=25
         KM=34
         KMM=KM-1
         KPM=KM+1
         LMM=LMAX 1
         LPM=LMAX+1

C SET LM-LMAX WHEN EITHER NO COORDINATE CONCENTRATION IS NEEDED

C OR WHEN CONCENTRATION IS NEEDED ONLY NEAR I=1

C SET LM LESS THAN LMAX WHEN CONCENTRATION IS NEEDED BOTH NEAR I=1

C AND I=LMAX.

C J1 IS THE INITIAL ITERATION INDEX.  MAY BE SET=2 TO 4

         LM=LPM/2
         LM=LMAX
         ITM=40
         01=4
         ELM=.002
         ANG=PI
         IF(MOD.EQ.1) GO TO 8500
         IF(MOD.EQ.2) GO TO 9500
  8500   CONTINUE
         DA=PI/KMM
         DO 351 K=1,KM
         XB(K,1)=COS(ANG)
         YB(K,1)=SIN(ANG)
         WD-SQRT(9.*(SIN(ANG)**2)+4.*(COS(ANG)**2))
         XB(K,LMAX)=6.*COS(ANG)/WD
         YB(K,LMAX)=6.*SIN(ANG)/WD
C        XB(K,LMAX)=3.*COS(ANG)
C        YB(K,LMAX)=3.*SIN(ANG)
         ANG=ANG-DA
  351    CONTINUE
         F2=FLOAT(LMM)
         DO 349 I=1,LMAX
         F1=FLOAT(I-1)
         XB(1,I)=1.-2.*F1/F2
         XB(KM,I)=1.+2.*F1/F2
         YB(1,I)=0.
         YB(KM,I)=0.
  349    CONTINUE
         IF(MOD.EQ.1) GO TO 360
  9500   CONTINUE
         READ(16,8001)((XB(K,1),YB(K,1)),K=1,KM)
         READ(16,8001)((XB(K,LMAX)YB(K,LMAX)),K=1,KM)
         READ(16,8001)((XB(1,I),YB(1,I)),I=1,LMAX)
         READ(16,8001)((XB(KM,I),YB(KM,I)),I=1,LMAX)






  8001   FORMAT(2F10.5)
  360    CONTINUE
         IF(TYPE.EQ.1.) GO TO 8000
         F3=FLOAT(LM-1)
         F4=FLOAT(LMAX-LM)
         IF(LM.EQ.LMAX) GO TO 9910
         DO 9900 I 1.LMAX
         F1=FLOAT(I,1)
         F2=FLOAT(LMAX I)
         G1=F1/F3
         G2=F2/F4
         IF(I.LE.LM) GO TO 9901
         IF(I.GT.LM) GO TO 9902
  9901   CONTINUE
         XB(1,I)=XB(1,I),(XB(1,LM) XB(1,1))*
     1   GK*(CK**(I-LM))
         YB(1,I)=YB(1,1),(YB(1,LM) YB(1,1))*
     1   G1*(CK**(I LM))
         XB(KM,I)-XB(KM,1),(XB(KM,LM)-XB(KM,1))*
     1   G1*(CK**(I-LM))
         YB(KM,I)-YB(KM,1),(YB(KM,LM)-YB(KM,1))*
     1   G1*(CK**(I-LM))
         GO TO 9905
  9902   CONTINUE
         XB(1,I)=XB(1,LMAX) (XB(1,LMAX) XB(1,LM)*
     1   G2*(CK**(LM I))
         YB(1,I)=YB(1,LMAX) (YB(1,LMAX) YB(1,LM))*
     1   G2*(CK**(LM-1))
         XB(KM,I)=XB(KM,LMAX) (XB(KM,LMAX) XB(KM,LM))*
     1   G2*(CK**(LM I))
         YB(KM,I)=YB(KM,LMAX) (YB(KM,LMAX) YB(KM,LM))*
     1   G2*(CK**(LM I))
  9905   CONTINUE
  9900   CONTINUE
         GO TO 8000
  9910   CONTINUE
         FS=FLCAT(LMM)
         DO 9909 I = 1,LMAX
         F1=FLOAT(I, 1)
         G1=F1/F5
         XB(1,I)=XB(1,1)+(XB(1,LMAX) XB(1,1))*
     1   G1*(CK**(I LMAX))
         YB(1,I)=YB(1,1)+(YB(1,LMAX) YB(1,1))*
     1   C1*(CK**(I LMAX))
         XB(KM,I)=XB(KM,1)+(XB(KM,LMAX) XB(KM,1))*
     1   G1*(CK**(I-LMAX))
         YB(KM,I)-YB(KM,1),(YB(KM,LMAX) YB(KM,1))*
     1   G1*(CK**(I-LMAX))
  9909   CONTINUE
  8000   CONTINUE
         DO 28 K 1,KM
         X(K,1)=XB(K,1)
         Y(K,1)=YB(K,1)
         X(K,LMAX)=XB(K,LMAX)






         Y(K,LMAX)=YB(K,LMAX)
  28     CONTINUE
         DO 29 I=1,LMAX
         X(1,I)=XB(1,I)
         Y(1,I)=YB(1,I)
         X(KM,I)=XB(KM,I)
         Y(KM,I)=YB(KM,I)
  29     CONTINUE
         DO 116 K=1,KM
         WRITE(2,109) K,XB(K,1), YB(K,1),XB(K,LMAX),YB(K,LMAX)
  109    FORMAT( 'INPUT DATA',I4,4F12.4)
  118    CONTINUE
         DO 117 I=1,LMAX
         WRITE(2,118) I,XB(1,I),YB(1,I),XB(KM,I),YB(KM,I)
  118    FORMAT( 'INPUT DATA ',I4,4F12.4)
  117    CONTINUE
C INITIAL GUESS:
         C5=FLOAT(LMM)
         C6=FLOAT(KMM)
         D5=FLOAT(LMM*KMM)
         DO 4 K=2,KMM
         DO 3 I=2,LMM
         C1=FLOAT(LMAX-I)
         C2=FLOAT(I-1)
         C3=FLOAT(KM-K)
         C4=FLOAT(K-1)
         D1=FLOAT(LMAX-I)*(KM-K))
         D2=FLOAT(LMAX-I)*(K-1))
         D3=FLOAT(I-1)*(KM-K))
         D4=FLOAT((I-1)*(K-1))
         E1=(C1*X(K,1)+C2*X(K,LMAX))/C5+
     1   (C3*X91,I)+C4*X(KM,I))/C6
         E2=(D1*X(1,1)+D2*X(KM,1)+D3*X(1,LMAX)+
     1   D4*X(KM,LMAX))/D5
         X(K,I)=E1-E2
         H1=(C1*Y(K,1)+C2*(K,LMAX))/C5+
     1   (C3*Y(1,I)+C4*Y(KM,I))/C6
         H2=(D1*Y(1,1)+D2*(KM,1)+D3*Y(1,LMAX)+
     1   D4*Y(KM,LMAX))/D5
         Y(K,I)=H1-H2
  3      CONTINUE
  4      CONTINUE
         DO 110 1JK=1,ITM
         CHECK=0.0
         DO 100 II=2,LMAX-1
         CALL GCOEF(IJK,II)
         PBAR(1,II)=.5*(-3.*F(1,II)+4.*F(2,II)-F(3,II))
         PBAR(KM,II)=.5*(3.*F(KM,II)-4.*F(KMM,II)+F(KMM-1,II))
         DO 204 K=2,KMM
         PBAR(K,II)=.5*(F(K+1,II)-F(K-1,II))
  204    CONTINUE
         DO 205 K=1,KM
         QBAR(K,II)=.5*(FINV(K,II+1)-FINV(K,II-1))
  205    CONTINUE






  620    CONTINUE
         DO 6 K=1,KM
         A(K,II)=F(K,II)+0.5*PBAR(K,II)
         B(K,II)=-2.*(F(K,II)+FINV(K,II)
         C(K,II)=F(K,II)-0.5*PBAR(K,II)
  6      CONTINUE
         DO 8 K=1,KM
         RT(K,II)=-(FINV(K,II)+0.5*QBAR(K,II))*X(K,II+1)
     1-(FINV(K,II)-.5*QBAR(K,II))*X(K,II-1)
         ST(K,II)=-(FINV(K,II)+0.5*QBAR(K,II))*Y(K,II+1)
     1-(FINV(K,II)-.5*QBAR(K,II))*Y(K,II-1)
  8      CONTINUE
C XP,YP DENOTE THE PRESENTLY AVAILABLE ITERATIVE VALUES
C XC,YC DENOTE THE PREVIOUS ITERATIVE VALUES
         DO 5 K=1,KM
         XP(K,II)=X(K,II)
         YP(K,II)=Y(K,II)
         XC(K,II)=X(K,II)
         YC(K,II)=Y(K,II)
  5      CONTINUE
         D(1,II)=0.
         H(1,II)=0.
         E(1,II)=XP(1,II)
         T(1,II)=YP(1,II)
         DO 10 K=2,KMM
         D(K,II)=-A(K,II)/(B(K,II)+C(K,II)*D(K-1,II))
         H(K,II)=-A(K,II)/(B(K,II)+C(K,II)*H(K-1,II))
         E(K,II)=(RT(K,II)-C(K,II)*E(K-1,II))/(B(K,II)+C(K,II)*D(K-1,II))
         T(K,II)=(ST(K,II)-C(K,II)*T(K-1,II))/(B(K,II)_C(K,II)*H(K-1,II))
  10     CONTINUE
         X(1,II)=W*XP(1,II)+(1.-W)*XC(1,II)
         Y(1,II)=W*YP(1,II)+(1.-W)*YC(1,II)
         DO 20 K=KMM,2,-1
         XP(K,II)=D(K,II)*XP(K+1,II)+E(K,II)
         YP(K,II)=H(K,II)*YP(K+1,II)+T(K,II)
         X(K,II)=W*XP(K,II)+(1.-W)*XC(K,II)
         Y(K,II)=W*YP(K,II)+(1.-W)*YC(K,II)
  20     CONTINUE
  100    CONTINUE
         ERX=0
         ERY=0
         DO 30 I=2,LMAX-1
         DO 30 K=1,KMM
         DIFX(K,I)=X(K,I)-XC(K,I)
         DIFY(K,I)=Y(K,I)-YC(K,I)
         CHECK=CHECK+ABS(DIFX(K,I))+ABS(DIFY(K,I))
         ERX=AMAX19ERX,DIFX(K,I))
         ERY=AMAX1(ERY,DIFY(K,I))
  30     CONTINUE
         KKK=IJK
         WRITE(2,609)KKK,CHECK,ERX,ERY
         IF(CHECK.LT.ELM) GO TO 120
  110    CONTINUE
  120    WRITE(2,609)KKK,CHECK,ERX,ERY






  609    FORMAT(' ITERATION NUMBER' ,15,' TOTAL ERRORS ARE'
     1,F10.6,/,' MAX POINT ERROR EX=',F10.7,2X,'EY=',F10.7)
         DO 1002 I=1,LMAX
         WRITE(2,1001)I,G13(1,I)
  1001   FORMAT('G13 DATA', I5,F12.5)
  1002   CONTINUE
         WRITE(2,900)
         DO 150 I=1,LMAX
         DO 160 K=1,KM
         WRITE(2,910)K,I,X(K,I),Y(K,I),AA(K,I),CC(K,I)
     1   ,G13(K,I),CJ2(K,I),F(K,I)
  160    CONTINUE
  150    CONTINUE
  900    FORMAT(1H,' K I ',FX, 'X',5X,'Y',5X,'AA',5X
     1   ,'CC',5X,'G13',5X,'CJ2',5X,'F')
  910    FORMAT(' ',2I3,7F7.3)
         WRITE(12,*) ((X(K,I),Y(K,I),K=1,KM),I=1,LMAX)
         STOP

         END
         SUBROUTINE GCOEF(IJK,II)
         PARAMETER (KMAX=40,IM=30)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)
     1,KM,KMM,LMAX,LM,LMM,J1,CK
         COMMON /ZUWA/ G13(KMAX,IM),CJ2(KMAXIM)
         COMMON /ZUWB/ AA(KMAXIM),BB(KMAX,IM),CC(KMAX,IM)
         COMMON /ZUWC/ DXXI(KMAX,IM),DYXI(KMAX,IM)
     1,DXZT(KMAX,IM),DYZT(KMAX,IM),F(KMAX,IM),FINV(KMAX,IM)
         DO 70 I=II-1,II+1
         DXZT(1,I)=(-3.*X(1,I)+4.*X(2,I)-X(3,I))/2.
         DYZT(1,I)=(-3.*Y(1,I)+4.*Y(2,I)-Y(3,I))/2.
         DXZT(KM,I)=(3.*X(KM,I)-4.*X(KMM,I)+X(KMM-1,I))/2.
         DYZT(KM,I)=(3.*Y(KM,I)-4.*Y(KMM,I)+Y(KMM-1,I))/2.
         DO 5 K=2,KMM
         DXZT(K,I)=(X(K+1,I)-X(K-1,I))/2.
         DYZT(K,I)=(Y(K+1,I)-Y(K-1,I))/2.
  5      CONTINUE
         DO 35 K=1,KM
         IF(I.EQ.1) THEN
         DXXI(K,I)=(-3.*X(K,I)+4.*X(K,I+1)-X(K,I+2))/2.
         DYXI(K,I)=(-3.*Y(K,I)+4.*Y(K,I+1)-Y(K,I+2))/2.
         GO TO 34
         END IF
         IF(I.EQ.LMAX) THEN
         DXXI(K,I)=(3.*X(K,I)-4.*X(K,I-1)+X(K,I-2))/2.
         DYXI(K,I)=(3.*Y(K,I)-4.*Y(K,I-1)+Y(K,I-2))/2.
         GO TO 34
         END IF
         DXXI(K,I)=(X(K,I+1)-X(K,I-1))/2.
         DYXI(K,I)=(Y(K,I+1)-Y(K,I-1))/2.
  34     CONTINUE
  35     CONTINUE
  70     CONTINUE







         IF(IJK.LT.J1) GO TO 558
         IF(IJK.GE.J1) GO TO 559
  559    CONTINUE
         DO 560 K=1,KM
         DO 561 I=II-1,II+1
         DXXI(K,I)=-F(K,I)*DYZT(K,I)
         DYXI(K,I)=F(K,I)*DXZT(K,I)
  561    CONTINUE
  560    CONTINUE
  558    CONTINUE
         DO 80 I=II-1,II+1
         DO 82 K=1,KM
         AA(K,I)=DXXI(K,I)**2+DYXI(K,1)**2
         G13(K,I)=DXXI(K,I)*DXZT(K,I)+DYXI(K,I)*DYZT(K,I)
         BB(K,I)=0.
         CC(K,I)=DXZT(K,I)**2+DYZT(K,I)**2
         CJ2(K,I)=SQRT(AA(K,I)/CC(K,I))
         FINV(K,I)=1./F(K,I)
  82     CONTINUE
  80     CONTINUE
         IF(IJK.LT.J1) GO TO 9013
         IF(IJK.GE.J1) GO TO 555
  555    CONTINUE
         IF(LM.EQ.LMAX) GO TO 777
         F4=FLOAT(LM-1)
         F5=FLOAT(LMAX-LM)
         G3=FLOAT(KMM)
         DO 778 I=II,II+1
         DO 779 K=2,KMM
         F1=FLOAT(I-1)
         F2=FLOAT(LM-1)
         F3=FLOAT(LMAX-I)
         G1=FLOAT(K-1)
         G2=FLOAT(KM-K)
         IF(I.LE.LM) GO TO 780
         IF(I.GT.LM) GO TO 781
  780    CONTINUE
         A1=(1.-(F1/F4)*(CK**(I-LM))*F(K,1)+(F1/F4)
     1   *(CK**(I-LM))*F(K,LM)+(G2/G3)*F(1,I)
     2   +(G1/G3)*F(KM,I)
         A2=((G2/G3)*F(1,1)+(G1/G3)*F(KM,1))*
     1   (1.-(F1/F4)*(CK**(I-LM)))+((G2/G3)*
     2   F(1,LM)+(G1/G3)*F(KM,LM))*(F1/F4)
     3   *(CK**(I-LM))
         F(K,I)=A1-A2
         FINV(K,I)=1./F(K,I)
         GO TO 782
  781    CONTINUE
         A1=(F3/F5)*(CK**(LM-I))*F(K,LM)+(1.-(F3/F5)*
     1   (CK**(LM-I)))*F(K,LMAX)+(G2/G3)*F(1,I)
     2   +(G1/G3)*F(KM,I)
         A2=((G2/G3)*F(1,LM)+(G1/G3*F(KM,LM))*(F3/F5)
     1   *(CK**(LM-1))+((G2/G3)*F(1,LMAX)+G1/G3)
     2   *F(KM,LMAX))*(1.-(F3/F5)*(CK**(LM-1)))






         F(K,I)=A1-A2
         FINV(K,I)=1./F(K,I)
  782    CONTINUE
  779    CONTINUE
  778    CONTINUE
         GO TO 9013
  777    CONTINUE
         F3=FLOAT(LMM)
         G3=FLOAT(KMM)
         DO 9111 I=II,II+1
         DO 9112 K=2,KMM
         F1=FLOAT(I-1)
         F2=FLOAT(K-1)
         G2=FLOAT(KM-K)
         A1=(1.-(F1/F3)*(CK**I-LMAX)))*F(K,1)
     1   +(F1/F3)*(CK**(I-LMAX))*F(K,LMAX)+(G2/G3)
     2   *F(1,I)+(G1/G3)*F(KM,I)
         A2=((G2/G3)*F(1,1)+(G1/G3)*F(KM,1))*(1.-
     1   (F1/F3)*(CK**(I-LMAX)))+((G2/G3)*F(1,LMAX)
     2   +(G1/G3)*F(KM,LMAX))*(F1/F3)*(CK**(I-LMAX))
         F(K,I)=A1-A2
         FINV(K,I)=1./F(K,I)
  9112   CONTINUE
  9111   CONTINUE
  9013   CONTINUE
         RETURN
         END







                __________ - _______________________________
                   _______________________________________
                         (Generally Non-Orthogonal)



The coordinates will be of the 0-type.


  PROGRAM BOFITO

C A PROGRAM FOR THE NUMERICAL GENERATION OF GENERAL COORDINATES

C IN 2D DOUBLY-CONNECTED REGIONS.

C THE GENERATED COORDINATES WILL BE OF THE 0-TYPE.

C THIS PROGRAM HAS BEEN DEVELOPED BY DR. Z.U.A. WARSI, VISITNG PROFESSOR

C AT VKI FROM MISSISSIPPI STATE UNIVERSITY DEPT. AEROSPACE ENG.

C MISS.STATE. MS39762, USA.

C THE COORDINATE K FOLLOWS THE GIVEN INNER AND OUTER BOUNDARIES

C IN THE CLOCKWISE SENSE

C ****************************************************************************

         PARAMEER(KMAX=40,IM=30)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)
     1,KM,KMM,LMAX
         COMMON /ZUWA/ G11(KMAX,IM),G13(KMAX,IM),G33(KMAX,IM)
     1,GG(KMAX,IM),CJ2(KMAX,IM)
         COMMON /ZUWB/ AA(KMAX,IM),BB(KMAX,IM),CC(KMAX,IM)
         COMMON /ZUWC/ DXXI(KMAX,IM),DYXI(KMAX,IM)
     1,DXZT(KMAX,IM),DYZT(KMAX,IM)

         DIMENSION A(KMAX,IM),B(KMAX,IM),C(KMAX,IM),RT(KMAX,IM)
     1,ST(KMAX,IM),XP(KMAX,IM),XC(KMAX,IM),YP(KMAX,IM),
     2YC(KMAX,IM),D(KMAX,IM),H(KMAX,IM),E(KMAX,IM)
     3T(KMAX,IM),P1(KMAX,IM),P2(KMAX,IM),PE(KMAX,IM),
     4QT(KMAX,IM),Q2(KMAX,IM),Q3(KMAX,IM),PBAR(KMAX,IM),
     5QBAR(KMAX,IM),DIFX(KMAX,IM),DIFY(KMAX,IM)

C MOD=1 STANDS FOR THE ANALYTIC INPUT DATA BOTH FOR THE INNER
C AND THE OUTER BOUNDARIES.
C MOD=2 STANDS FOR READING THE DATA THROUGH AN INPUT FILE.

         MOD=1
         W=1.5
         PI=3.14159264
         LMAX=25
         KM=37






         KMM=KM-1
         ELM=0.1
         ITM=40
         CK=1.1
         IF(MOD.EQ.1) GO TO 570
         IF(MOD.EQ.2) GO TO 571
  570    CONTINUE
         ANG=2.*PI
         DA=2.*PI/KMM
         AI=2.
         BI=1.5
         AO=3.
         BO=3.
         DO 108 K=1,KM
         XB(K,1)=AI*COS(ANG)
         YB(K,1)=BI*SIN(ANG)
         XB(K,LMAX)=AO*COS(ANG)
         YB(K,LMAX)=BO*SIN(ANG)
         DEG=ANG*180./PI
         WRITE(6,109) DEG,XB(K,1),YB(K,1),XB(K,LMAX),YB(K,LMAX)
  109    FORMAT(' INPUT DATA ',F10.4,4F12.4)
         ANG=ANG-DA
  108    CONTINUE
         IF(MOD.EQ.1) GO TO 572
  571    CONTINUE
         READ(18,8001)((XB(K,1),YB(K,1)),K=1,KM)
         READ(18,8001)((XB(K,LMAX),YB(K,LMAX)),K=1,KM)
  8001   FORMAT(2F10.5)
  572    CONTINUE
         DO 28 K=1,KM
         X(K,1)=XB(K,1)
         Y(K,1)=YB(K,1)
         X(K,LMAX)=XB(K,LMAX)
         Y(K,LMAX)=YB(K,LMAX)
  28     CONTINUE
C INITIAL GUESS: LINEAR INTERPOLATION
         DO 4K=1,KMM
         DX=(X(K,LMAX)-X(K,1))/(LMAX-1)
         DY=(Y(K,LMAX)-Y(K,1))/(LMAX-1)
         DO 3 I=2,LMAX-1
         X(K,I)=X(K,1)+(I-1)*DX
         Y(K,I)=Y(K,1)+(I-1)*DY
  3      CONTINUE
  4      CONTINUE
         DO 110 IJK=1,ITM
         CHECK=0.0
         DO 560 I=2,LMAX-1
         X(KM,I)=X(1,I)
         Y(KM,I)=Y(1,I)
  560    CONTINUE
  555    CONTINUE
         DO 100 II=2,LMAX 1
         CALL GCOEF(II)







C INSERT EXPRESSIONS FOR P1,P2,P3,Q1,Q2,Q3 BELOW.

C THE PARAMETER CK IS THE COORDINATE CONTROL PARAMETER

         DO 777 K=1,KMM
         P1(K,II)=0
         P2(K,II)=0
         P3(K,II)=0
         Q1(K,II)=0
         Q2(K,II)=0
         Q3(K,II)=-(2.+(II-1)*ALOG(CK))*ALOG(CK)
     1   /(1.*(II-1)*ALOG(CK))
  777    CONTINUE
         DO 504 K=1,KMM
         PBAR(K,II)=P1(K,II)*AA(K,II)-2.0*P2(K,II)*
     1BB(K,II)+P3(K,II)*AA(K,II)-2.0*Q2(K,II)*
         QBAR(K,II)=Q1(K,II)*AA(K,II)-2.0*Q2(K,II)*
     1BB(K,II)+Q3(K,II)*CC(K,II)
  504    CONTINUE
  620    CONTINUE
         DO 6 K=1,KMM
         A(K,II)=AA(K,II)+0.5*PBAR(K,II)
         B(K,II)=-2.*(AA(K,II)+CC(K,II))
         C(K,II)=AA(K,II)-0.5*PBAR(K,II)
  6      CONTINUE
         RT(1,II)=-(CC(1,II)+0.5*QBAR(1,II))*X(1,II+1)-(CC(1,II)-
     10.5*QBAR(1,II)*X(1,II-1)+0.5*BB(1,II)*(X(2,II+1)-
     1X(2,II-1)-X(KMM,II+1)+X(KMM,II-1))
         ST(1,II)=-(CC(1,II)+0.5*QBAR(1,II))*Y(1,II+1)-(CC(1,II
     1)-0.5*QBAR(1,II)*Y(KMM,II-1))*0.5
         DO 8 K=2,KMM
         RT(K,II)=-(CC(K,II)+0*QBAR(K,II))*X(K,II+1)-(CC(K,II)
     1-0.5*QBAR(K,II))*X(K,II-1)+0.5*BB(K,II)*(X(K+1,II+1)-
     2X(K+1,II-1)-X(K-1,II+1)+X(K-1,II-1))
         ST(K,II)=-(CC(K,II)+0.5*QBAR(K,II))*Y(K,II+1)-(CC(K,II)
     1-0.5*QBAR(K,II)*Y(K,II-1)+BB(K,II)*(Y(K+1,II+1)-Y(K+1,II-1)-
     2Y(K-1,II+1)+Y(K-1,II-1))*0.5
  8      CONTINUE
C XP,YP ARE THE 'PRESENTLY AVAILABLE' SOLUTIONS
C XC,YC ARE THE 'PREVIOUSLY AVAILABLE' SOLUTIONS
         A(KM,II)=A(1,II)
         B(KM,II)=B(1,II)
         C(KM,II)=C(1,II)
         RT(KM,II)=RT(1,II)
         ST(KM,II)=ST(1,II)
         DO 5 K=1,KM
         XP(K,II)=X(K,II)
         YP(K,II)=Y(K,II)
         XC(K,II)=X(K,II)
         YC(K,II)=Y(K,II)
  5      CONTINUE
         D(1,II)=0.
         H(1,II)=0.
         E(1,II)=XP(1,II)






         T(1,II)=YP(1,II)
         DO 10 K=2,KM
         D(K,II)=-A(K,II)/(B(K,II)+C(K,II)*D(K-1,II))
         H(K,II)=-A(K,II)/(B(K,II)+C(K,II)*H(K-1,II))
         E(K,II)=(RT(K,II)-C(K,II)*E(K-1,II))/(B(K,II)+C(K,II)*D(K-1,II))
         T(K,II)=(ST(K,II)-C(K,II)*T(K-1,II))/(B(K,II)+C(K,II)*H(K-1,II))
  10     CONTINUE
         XP(1,II)=D(KM,II)*XP(2,II)+E(KM,II)
         YP(1,II)=H(KM,II)*YP(2,II)+T(KM,II)
         X(1,II)=W*XP(1,II)+(1.-W)*XC(1,II)
         Y(1,II)=W*YP(1,II)+(1.-W)*YC(1,II)
         DO 20 K=KMM,2,-1
         XP(K,II)=D(K,II)*XP(K+1,II)+E(K,II)
         YP(K,II)=H(K,II)*YP(K+1,II)+T(K,II)
         X(K,II)=W*XP(K,II)+(1.-W)*XC(K,II)
         Y(K,II)=W*YP(K,II)+(1.-W)*YC(K,II)
  20     CONTINUE
  100    CONTINUE
         ERX=0.
         ERY=0.
         DO 30 I=2,LMAX-1
         DO 30 K=1,KMM
         DIFX(K,I)=X(K,I)-XC(K,I)
         DIFY(K,I)=Y(K,I)-YC(K,I)
         CHECK=CHECK+ABS(DIFX(K,I))+ABS(DIFY(K,I))
         ERX=AMAX1(ERX,DIFX(K,I))
         ERY=AMAX1(ERY,DIFY(K,I))
  30     CONTINUE
         KKK=IJK
         WRITE(6,609)KKK,CHECK,ERX,ERY
         IF(CHECK.LT.ELM) GO TO 120
  110    CONTINUE
  120    WRITE(6,609)KKK,CHECK,ERX,ERY
  609    FORMAT(' ITERATION NUMBER' ,I5,' TOTAL ERRORS ARE'
     1,F10.6,/,' MAX POINT ERROR EX=' ,F10.7,2X,'EY=',F10.7)
C THE FOLLOWING STATEMENTS ARE MEANT TO CHECK THE VALUE
C OF G13 AT SELECTED K-VALUES.
         DO 1002 I=1,LMAX
         WRITE(6,1001)I,G13(3,I)
  1001   FORMAT(' G13 DATA',I5,F12.5)
  1002   CONTINUE
         WRITE(6,900)
         DO 150 I=1,LMAX
         DO 160 K=1,KMM
         WRITE(6,910)K,I,X(K,I),Y(K,I),AA(K,I),CC(K,I)
     1   ,G13(K,I),CJ2(K,I)
  160    CONTINUE
  150    CONTINUE
  900    FORMAT(1H,' K I ',5X, 'X',5X,'Y',5X,'AA',5X,'CC'
     1   ,5X,'G13',5X,'CJ2')
  910    FORMAT(' ',213,6F7.3)
         WRITE(10,*) ((X(K,I,Y(K,I),K=1,KM),I=1,LMAX)
         STOP
         END







         SUBROUTINE GCOEF(II)
         PARAMETER (KMAX=40,IM=30)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)
     1,KM,KMM,LMAX
         COMMON /ZUWA/ G11(KMAX,IM),G13(KMAX,IM),G33(KMAX,IM)
     1,GG(KMAX,IM),CJ2(KMAX,IM)
         COMMON /ZUWB/ AA(KMAX,IM),BB(KMAX,IM),CC(KMAX,IM)
         COMMON /ZUWC/ DXXI(KMAX,IM),DYXI(KMAX,IM)
     1,DXZT(KMAX,IM),DYZT(KMAX,IM)

         DO 1 I=1,LMAX
         X(KM,I)=X(1,I)
         Y(KM,I)=Y(1,I)
     1   CONTINUE
         DO 70 I=II-1,II+1
         DXZT(1,I0=(X92,I)-X(KMM,I))/2.
         DYZT(1,I)=(Y(2,I)-Y(KMM,I))/2.
         DO 5 K=2,KMM
         DXZT(K,I)=(X(K+1,I)-X(K-1,I))/2.
         DYZT(K,I)=(Y(K+1,I)-Y(K-1,I))/2.
  5      CONTINUE
         DO 35 K=1,KMM
         IF(I.EQ.1) THEN
         DXXI(K,I)=(-3.*X(K,I)+4.*X(K,I+1)-X(K,I+2))/2.
         DYXI(K,I)=(-3.*Y(K,I)+4.*Y(K,I+1)-Y(K,I+2))/2.
         GO TO 34
         END IF
         IF(I.EQ.LMAX) THEN
         DXXI(K,I)=(3.*X(K,I)-4.*X(K,I-1)+X(K,I-2))/2.
         DYXI(K,I)=(3.*Y(K,I)-4.*Y(K,I-1)+Y(K,I-2))/2.
         GO TO 34
         END IF
         DXXI(K,I)=(X(K,I+1)-X(K,I-1))/2.
         DYXI(K,I)=(Y(K,I+1)-Y(K,I-1))/2.
  34     CONTINUE
  35     CONTINUE
  70     CONTINUE
         DO 80 I=II,II+1
         DO 82 K=1,KMM
         AA(K,I)=DXXI(K,I)**2+DYXI(K,I)**2
         G13(K,I)=DXXI(K,I)*DXZT(K,I)+DYXI(K,I)*DYZT(K,I)
         BB(K,I)=G13(K,I)
         CC(K,I)=DXZT(K,I)**2+DYZT(K,I)**2
         CJ2(K,I)=SQRT(AA(K,I)*CC(K,I)-BB(K,I)**2)
  82     CONTINUE
  80     CONTINUE
         DO 1111 I=II,II+1
         AA(KM,I)=AA(1,I)
         BB(KM,I)=BB(1,I)
         CC(KM,I)=CC(1,I)
         CJ2(KM,I)=CJ2(1,I)
  1111   CONTINUE
         RETURN






         END









                __________ - _______________________________
                   _______________________________________
                           (Generally Orthogonal)


The coordinates will be of the 0-type.


     PROGRAM BOFITOOR

C A PROGRAM FOR THE NUMERICAL GENERATION OF ORTHOGONAL COORDINATES

C IN 2D DOUBLY-CONNECTED REGIONS.

C THE RESULTING COORDINATES WILL BE OF THE O-TYPE.

C AT VKI FROM MISSISSIPPI STATE UNIVERSITY DEPT.AEROSPACE ENG.

C MISS.STATE. MS39762,USA.

C THE COORDINATE K FOLLOWS THE GIVEN INNER AND OUTER BOUNDARIES

C IN THE CLOCKWISE SENSE.

C ****************************************************************************

         PARAMETER(KMAX=40,IM=30)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)
     1,KM,KMM,LMAX,LMM,J1,CK
         COMMON /ZUWA/ G13(KMAX,IM),CJ2(KMAX,IM)
         COMMON /ZUWB/ AA(KMAX,IM),BB(KMAX,IM),CC(KMAX,IM)
         COMMON /ZUWC/ DXXI(KMAX,IM),DYXI(KMAX,IM)
     1,DXZT(KMAX,IM),DYZT(KMAX,IM),F(KMAX,IM),FINV(KMAX,IM)

         DIMENSION A(KMAX,IM),B(KMAX,IM),C(KMAX,IM),RT(KMAX,IM)
     1,ST(KMAX,IM),XP(KMAX,IM),XC(KMAX,IM),YP(KMAX,IM),
     2YC(KMAX,IM),D(KMAX,IM),H(KMAX,IM),E(KMAX,IM),
     3T(KMAX,IM),PBAR(KMAX,IM),
     4QBAR(KMAX,IM),DIFX(KMAX,IM),DIFY(KMAX,IM)

C MOD=1 STANDS FOR THE ANALYTIC INPUT DATA BOTH FOR THE INNER
C AND THE OUTER BOUNDARIES.
C MOD=2 STANDS FOR THE INTPUT DATA TO BE READ THROUGH A FILE.

         MOD=1
         CK=1.1
         W=1.4
         PI=3.14159264
         LMAX=25
         KM=37
         KMM=KM-1
         01=4
         ELM=.1
         ITM=40






         IF(MOD.EQ.1) GO TO 570
         IF(MOD.EQ.2) GO TO 571
  570    CONTINUE
         ANG=2.*PI
         DA=2.*PI/KMM
         AI=1.
         BI=1.
         AO=3.
         BO=2.
         AS=AO**2
         BS=BO**2
         DO 108 K=1,KM
         XB(K,1)=AI*COS(ANG)
         YB(K,1)=BI*SIN(ANG)
         WD=SQRT(AS*(SIN(ANG)**2)+BS*(COS(ANG)**2))
         XB(K,LMAX)=AO*BO*COS(ANG)/WD
         YB(K,LMAX)=AO*BO*SIN(ANG)/WD
         DEG=ANG*180./PI
         WRITES(6,109) DEG,XB(K,1),YB(K,1),XB(K,LMAX),YB(K,LMAX)
  109    FORMAT(' INPUT DATA ',F10.4,4F12.4)
         ANG=ANG-DA
  108    CONTINUE
         IF(MOD.EQ.1) GO TO 572
  571    CONTINUE
         READ(19,8001)((XB(K,1),YB(K,1)),K=1,KM)
         READ(19,8001)((XB(K,LMAX),YB(K,LMAX)),K=1,KM)
  8001   FORMAT(2F10.5)
  572    CONTINUE
         DO 28 K=1,KM
         X(K,1)=XB(K,1)
         Y(K,1)=YB(K,1)
         X(K,LMAX)=XB(K,LMAX)
         Y(K,LMAX)=YB(K,LMAX)
  28     CONTINUE
C INITIAL GUESS: LINEAR INTERPOLATION
         DO 4 K=1,KMM
         DX=(X(K,LMAX)-X(K,1))/(LMAX-1)
         DY=(Y(K,LMAX)-Y(K,1))/(LMAX-1)
         DO 3 I=2,LMAX-1
         X(K,I)=X(K,1)+(I-1)*DX
         Y(K,I)=Y(K,1)+(I-1)*DY
  3      CONTINUE
  4      CONTINUE
         DO 100 IJK=1,ITM
         CHECK=0.0
         DO 560 I=2,LMAX-1
         X(KM,I)=X(1,I)
         Y(KM,I)=Y(1,I)
  560    CONTINUE
         DO 100 II=2,LMAX-1
         CALL GCOEF(IJK,II)
         PBAR(1,II)=.5*(F(2,II)-F(KMM,II))
         PBAR(KM,II)=PBAR(1,II)
         DO 204 K=2,KMM






         PBAR(K,II)=.5*(F(K+1,II)-F(K-1,II))
  204    CONTINUE
         DO 206 K=1,KMM
         QBAR(K,II)=.5*(FINV(K,II+1)-FINV(K,II-1))
  206    CONTINUE
         QBAR(KM,II)=QBAR(1,II)
         DO 6 K=1,KMM
         A(K,II)=F(K,II)+0.5*PBAR(K,II)
         B(K,II)=-2.*(F(K,II)+FINV(K,II))
         C(K,II)=F(K,II)-0.5*PBAR(K,II)
  6      CONTINUE
         RT(1,II)=-(FINV(1,II)+0.5*QBAR(1,II))*X(1,II+1)-(FINV(1,II)-
     1.5*QBAR(1,II))*X(1,II-1)
         ST(1,II)=-(FINV(1,II)+0.5*QBAR(1,II))*Y(1,II+1)-(FINV(1,II
     1)-0.5*QBAR(1,II))*Y(1,II-1)
         DO 8 K=2,KMM
         RT(K,II)=-(FINV(K,II)+0.5*QBAR(K,II))*X(K,II+1)-(FINV(K,II)
     1-0.5*QBAR(K,II))*X(K,II-1)
         ST(K,II)=-(FINV(K,II)+0.5*QBAR(K,II))*Y(K,II+1)-(FINV(K,II)
     1-0.5*QBAR(K,II))*Y(K,II-1)
  8      CONTINUE
C XP,YP ARE THE 'PRESENTLT AVAILABLE' SOLUTIONS
C XC,YC ARE THE 'PRVIOUSLY AVAILABLE' SOLUTIONS
         A(KM,II)=A(1,II)
         B(KM,II)=B(1,II)
         C(KM,II)=C(1,II)
         RT(KM,II)=RT(1,II)
         ST(KM,II)=ST(1,II)
         DO 5 K=1,KM
         XP(K,II)=X(K,II)
         YP(K,II)=Y(K,II)
         XC(K,II)=X(K,II)
         YC(K,II)=Y(K,II)
  5      CONTINUE
         D(1,II)=0.
         H(1,II)=0.
         E(1,II)=XP(1,II)
         T(1,II)=YP(1,II)
         DO 10 K=2,KM
         D(K,II)=-A(K,II)/(B(K,II)+C(K,II)*D(K-1,II))
         H(K,II)=-A(K,II)/(B(K,II)+C(K,II)*H(K-1,II))
         E(K,II)=(RT(K,II)-C(K,II)*E(K-1,II))/(B(K,II)+C(K,II)*D(K-1,II))
         T(K,II)=(ST(K,II)-C(K,II)*T(K-1,II))/(B(K,II)+C(K,II)*H(K-1,II))
  10     CONTINUE
         XP(1,II)=D(KM,II)*XP(2,II)+E(KM,II)
         YP(1,II)=H(KM,II)*YP(2,II)+T(KM,II)
         X(1,II)=W*(XP(1,II)+(1.-W)*XC(1,II)
         Y(1,II)=W*(YP(1,II)+(1.-W)*YC(1,II)
         DO 20 K=KMM,2,-1
         XP(K,II)=D(K,II)*XP(K+1,II)+E(K,II)
         YP(K,II)=H(K,II)*YP(K+1,II)+T(K,II)
         X(K,II)=W*XP(K,II)+(1.-W)*XC(K,II)
         Y(K,II)=W*YP(K,II)+(1.-W)*YC(K,II)
  20     CONTINUE






  100    CONTINUE
         ERX=0.
         ERY=0.
         DO 30 I=2,LMAX-1
         DO 30 K=1,KMM
         DIFX(K,I)=X(K,I)-XC(K,I)
         DIFY(K,I)=Y(K,I)-YC(K,I)
         CHECK=CHECK+ABS(DIFX(K,I))+ABS(DIFY(K,I))
         ERX=AMAX1(ERX,DIFX(K,I))
         ERY=AMAX1(ERY,DIFY(K,I))
  30     CONTINUE
         KKK=IJK
         WRITE(6,609)KKK,CHECK,ERX,ERY
         IF(CHECK.LT.ELM) GO TO 120
  110    CONTINUE
  120    WRITE(6,609)KKK,CHECK,ERX,ERY
  609    FORMAT('ITERATION NUMBER',15,' TOTAL ERRORS ARE'
     1,F10.6,/,'MAX POINT ERROR EX=',F10.7,2X,'EY=',F10.7)
C THE FOLLOWING STATEMENTS ARE MEANT TO CHECK THE VALUES OF
C G13 AT SELECTED K-VALUES.
         DO 1002 I=1,LMAX
         WRITE(6,1001)I,G13(7,I)
  1001   FORMAT('G13 DATA',I5,F12.5)
  1002   CONTINUE
         WRITE(6,900)
         DO 150 I=1,LMAX
         DO 160 K=1,KMM
         WRITE(6,910)K,I,X(K,I),Y(K,I),AA(K,I),CC(K,I)
     1   ,G13(K,I),CJ2(K,I),F(K,I)
  160    CONTINUE
  150    CONTINUE
  900    FORMAT(1H,'K I ',5X, 'X',5X,'Y',5X,'AA',5X,'CC'
     1   ,5X,'G13',5X,'CJ2',5X,'F')
  910    FORMAT(' ',2I3,7F7.3)
         WRITE(10,*) ((X(K,I),Y(K,I),K=1,KM),I=1,LMAX)
         STOP
         END

         SUBROUTINE GCOEF(IJK,II)
         PARAMETER (KMAX=40,IM=30)
         COMMON /XY/ XB(KMAX,IM),YB(KMAX,IM),X(KMAX,IM),Y(KMAX,IM)
     1,KM,KMM,LMAX,LMM,J1,CK
         COMMON /ZUWA/ G13(KMAX,IM),CJ2(KMAX,IM)
         COMMON /ZUWB/ AA(KMAX,IM),BB(KMAX,IM),CC(KMAX,IM)
         COMMON /ZUWC/ DXXI(KMAX,IM),DYXI(KMAX,IM)
     1,DXZT(KMAX,IM),DYZT(KMAX,IM),F(KMAX,IM),FINV(KMAX,IM)
         DO 1 I=1,LMAX
         X(KM,I)=X(1,I)
         Y(KM,I)=Y(1,I)
  1      CONTINUE
         DO 70 I=II-1,II+1
         DXZT(1,I)=(X(2,I)-X(KMM,I))/2.
         DYZT(1,I)=(Y(2,I)-Y(KMM,I))/2.
         DO 5 K=2,KMM






         DXZT(K,I)=(X(K+1,I)-X(K-1,I))/2.
         DYZT(K,I)=(Y(K+1,I)-Y(K-1,I))/2.
  5      CONTINUE
         DO 35 K=1,KMM
         IF(I.EQ.1) THEN
         DXXI(K,I)=(-3.*X(K,I)+4.*X(K,I+1)-X(K,I+2))/2.
         DYXI(K,I)=(-3.*Y(K,I)+4.*Y(K,I+1)-Y(K,I+2))/2.
         GO TO 34
         END IF
         IF(I.EQ.LMAX) THEN
         DXXI(K,I)=(3.*X(K,I)-4.*X(K,I-1)+X(K,I-2))/2.
         DYXI(K,I)=(3.*Y(K,I)-4.*Y(K,I-1)+Y(K,I-2))/2.
         GO TO 34
         END IF
         DXXI(K,I)=(X(K,I+1)-X(K,I-1))/2.
         DYXI(K,I)=(Y(K,I+1)-Y(K,I-1))/2.
  34     CONTINUE
  35     CONTINUE
  70     CONTINUE
         IF(IJK.LT.J1) GO TO 558
         IF(IJK.GE.J1) GO TO 559
  559    CONTINUE
         DO 560 K=1,KM
         DO 561 I=II-1,II+1
         DXXI(K,I)=-F(K,I)*DYZT(K,I)
         DYXI(K,I)=F(K,I)*DXZT(K,I)
  561    CONTINUE
  560    CONTINUE
  558    CONTINUE
         DO 80 I=II,II+1
         DO 82 K=1,KMM
         AA(K,I)=DXXI(K,I)**2+DYXI(K,I)**2
         G13(K,I)=DXXI(K,I)*DXZT(K,I)+DYXI(K,I)*DYZT(K,I)
         BB(K,I)=G13(K,I)
         BB(K,I)=0.
         CC(K,I)=DXZT(K,I)**2+DYZT(K,I)**2
         CJ2(K,I)=SQRT(AA(K,I)*CC(K,I)-BB(K,I)**2)
         F(K,I)=SQRT(AA(K,I)/CC(K,I))
         FINV(K,I)=1./F(K,I)
  82     CONTINUE
  80     CONTINUE
         DO 1111 I=II,II+1
         AA(KM,I)=AA(1,I)
         BB(KM,I)=BB(1,I)
         CC(KM,I)=CC(1,I)
         CJ2(KM,I)=CJ2(1,I)
         F(KM,I)=F(1,I)
         FINV(KM,I)=FINV(1,I)
  1111   CONTINUE
         IF(IJK.LT.J1) GO TO 635
         IF(IJK.GE.J1) GO TO 555
  555    CONTINUE
         DO 556 K=1,KMM
         F(K,1)=SQRT(AA(K,1)/CC(K,1))






         F(K,LMAX)=SQRT(AA(K,LMAX)/CC(K,LMAX))
  556    CONTINUE
         F(KM,1)=F(1,1)
         F(KM,LMAX)=F(1,LMAX)
         DO 557 I=II-1,II+1
         F(1,I)=SQRT(AA(1,I)/CC(1,I))
         F(KM,I)=F(1,I)
  557    CONTINUE
         F2=FLOAT(LMAX-1)
         DO 9111 I=II-1,II+1
         DO 9112 K=2,KMM
         F1=FLOAT(I-1)
         A1=(1.-(F1/F2)*(CK**(I-LMAX)))*(F(K,1)-F(1,1))
         A2=(F1/F2)*(CK**(I-LMAX))*(F(K,LMAX)-F(1,LMAX))
         F(K,I)=A1+A2+F(1,I)
         FINV(K,I)=1./F(K,I)
  9112   CONTINUE
  9111   CONTINUE
  635    CONTINUE
         RETURN
         END

