c     Routines from SPARSPAK that carry out Reverse Cuthill--McKee 
c     node ordering.  Main routine is GENRCM, at and of this file.


C----- SUBROUTINE ROOTLS
C***************************************************************        
C***************************************************************        
C********     ROOTLS ..... ROOTED LEVEL STRUCTURE      *********        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - ROOTLS GENERATES THE LEVEL STRUCTURE ROOTED             
C        AT THE INPUT NODE CALLED ROOT. ONLY THOSE NODES FOR            
C        WHICH MASK IS NONZERO WILL BE CONSIDERED.                      
C                                                                       
C     INPUT PARAMETERS -                                                
C        ROOT - THE NODE AT WHICH THE LEVEL STRUCTURE IS TO             
C               BE ROOTED.                                              
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE              
C               GIVEN GRAPH.                                            
C        MASK - IS USED TO SPECIFY A SECTION SUBGRAPH. NODES            
C               WITH MASK(I)=0 ARE IGNORED.                             
C                                                                       
C     OUTPUT PARAMETERS -                                               
C        NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE.         
C        (XLS, LS) - ARRAY PAIR FOR THE ROOTED LEVEL STRUCTURE.         
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )    
C                                                                       
C***************************************************************        
C                                                                       
         INTEGER ADJNCY(1), LS(1), MASK(1), XLS(1)                      
         INTEGER XADJ(1), I, J, JSTOP, JSTRT, LBEGIN,                   
     1           CCSIZE, LVLEND, LVSIZE, NBR, NLVL,                     
     1           NODE, ROOT                                             
C                                                                       
C***************************************************************        
C                                                                       
C        ------------------                                             
C        INITIALIZATION ...                                             
C        ------------------                                             
         MASK(ROOT) = 0                                                 
         LS(1) = ROOT                                                   
         NLVL = 0                                                       
         LVLEND = 0                                                     
         CCSIZE = 1                                                     
C        -----------------------------------------------------          
C        LBEGIN IS THE POINTER TO THE BEGINNING OF THE CURRENT          
C        LEVEL, AND LVLEND POINTS TO THE END OF THIS LEVEL.             
C        -----------------------------------------------------          
  200    LBEGIN = LVLEND + 1                                            
         LVLEND = CCSIZE                                                
         NLVL = NLVL + 1                                                
         XLS(NLVL) = LBEGIN                                             
C        -------------------------------------------------              
C        GENERATE THE NEXT LEVEL BY FINDING ALL THE MASKED              
C        NEIGHBORS OF NODES IN THE CURRENT LEVEL.                       
C        -------------------------------------------------              
         DO 400 I = LBEGIN, LVLEND                                      
            NODE = LS(I)                                                
            JSTRT = XADJ(NODE)                                          
            JSTOP = XADJ(NODE + 1) - 1                                  
            IF ( JSTOP .LT. JSTRT )  GO TO 400                          
               DO 300 J = JSTRT, JSTOP                                  
                  NBR = ADJNCY(J)                                       
                  IF (MASK(NBR) .EQ. 0) GO TO 300                       
                     CCSIZE = CCSIZE + 1                                
                     LS(CCSIZE) = NBR                                   
                     MASK(NBR) = 0                                      
  300          CONTINUE                                                 
  400    CONTINUE                                                       
C        ------------------------------------------                     
C        COMPUTE THE CURRENT LEVEL WIDTH.                               
C        IF IT IS NONZERO, GENERATE THE NEXT LEVEL.                     
C        ------------------------------------------                     
         LVSIZE = CCSIZE - LVLEND                                       
         IF (LVSIZE .GT. 0 ) GO TO 200                                  
C        -------------------------------------------------------        
C        RESET MASK TO ONE FOR THE NODES IN THE LEVEL STRUCTURE.        
C        -------------------------------------------------------        
         XLS(NLVL+1) = LVLEND + 1                                       
         DO 500 I = 1, CCSIZE                                           
            NODE = LS(I)                                                
            MASK(NODE) = 1                                              
  500    CONTINUE                                                       
         RETURN                                                         
      END                                                               
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
C----- SUBROUTINE FNROOT                                                
C***************************************************************        
C***************************************************************        
C*******     FNROOT ..... FIND PSEUDO-PERIPHERAL NODE    *******        
C***************************************************************        
C***************************************************************        
C                                                                       
C    PURPOSE - FNROOT IMPLEMENTS A MODIFIED VERSION OF THE              
C       SCHEME BY GIBBS, POOLE, AND STOCKMEYER TO FIND PSEUDO-          
C       PERIPHERAL NODES.  IT DETERMINES SUCH A NODE FOR THE            
C       SECTION SUBGRAPH SPECIFIED BY MASK AND ROOT.                    
C                                                                       
C    INPUT PARAMETERS -                                                 
C       (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR THE GRAPH.        
C       MASK - SPECIFIES A SECTION SUBGRAPH. NODES FOR WHICH            
C              MASK IS ZERO ARE IGNORED BY FNROOT.                      
C                                                                       
C    UPDATED PARAMETER -                                                
C       ROOT - ON INPUT, IT (ALONG WITH MASK) DEFINES THE               
C              COMPONENT FOR WHICH A PSEUDO-PERIPHERAL NODE IS          
C              TO BE FOUND. ON OUTPUT, IT IS THE NODE OBTAINED.         
C                                                                       
C    OUTPUT PARAMETERS -                                                
C       NLVL - IS THE NUMBER OF LEVELS IN THE LEVEL STRUCTURE           
C              ROOTED AT THE NODE ROOT.                                 
C       (XLS,LS) - THE LEVEL STRUCTURE ARRAY PAIR CONTAINING            
C                  THE LEVEL STRUCTURE FOUND.                           
C                                                                       
C    PROGRAM SUBROUTINES -                                              
C       ROOTLS.                                                         
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  FNROOT ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )    
C                                                                       
C***************************************************************        
C                                                                       
         INTEGER ADJNCY(1), LS(1), MASK(1), XLS(1)                      
         INTEGER XADJ(1), CCSIZE, J, JSTRT, K, KSTOP, KSTRT,            
     1           MINDEG, NABOR, NDEG, NLVL, NODE, NUNLVL,               
     1           ROOT                                                   
C                                                                       
C***************************************************************        
C                                                                       
C        ---------------------------------------------                  
C        DETERMINE THE LEVEL STRUCTURE ROOTED AT ROOT.                  
C        ---------------------------------------------                  
         CALL  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NLVL, XLS, LS )       
         CCSIZE = XLS(NLVL+1) - 1                                       
         IF ( NLVL .EQ. 1 .OR. NLVL .EQ. CCSIZE ) RETURN                
C        ----------------------------------------------------           
C        PICK A NODE WITH MINIMUM DEGREE FROM THE LAST LEVEL.           
C        ----------------------------------------------------           
  100    JSTRT = XLS(NLVL)                                              
         MINDEG = CCSIZE                                                
         ROOT = LS(JSTRT)                                               
         IF ( CCSIZE .EQ. JSTRT )  GO TO 400                            
            DO 300 J = JSTRT, CCSIZE                                    
               NODE = LS(J)                                             
               NDEG = 0                                                 
               KSTRT = XADJ(NODE)                                       
               KSTOP = XADJ(NODE+1) - 1                                 
               DO 200 K = KSTRT, KSTOP                                  
                  NABOR = ADJNCY(K)                                     
                  IF ( MASK(NABOR) .GT. 0 )  NDEG = NDEG + 1            
  200          CONTINUE                                                 
               IF ( NDEG .GE. MINDEG ) GO TO 300                        
                  ROOT = NODE                                           
                  MINDEG = NDEG                                         
  300       CONTINUE                                                    
C        ----------------------------------------                       
C        AND GENERATE ITS ROOTED LEVEL STRUCTURE.                       
C        ----------------------------------------                       
  400    CALL  ROOTLS ( ROOT, XADJ, ADJNCY, MASK, NUNLVL, XLS, LS )     
         IF (NUNLVL .LE. NLVL)  RETURN                                  
            NLVL = NUNLVL                                               
            IF ( NLVL .LT. CCSIZE )  GO TO 100                          
            RETURN                                                      
      END                                                               
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
C----- SUBROUTINE DEGREE                                                
C***************************************************************        
C***************************************************************        
C********     DEGREE ..... DEGREE IN MASKED COMPONENT   ********        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - THIS ROUTINE COMPUTES THE DEGREES OF THE NODES          
C        IN THE CONNECTED COMPONENT SPECIFIED BY MASK AND ROOT.         
C        NODES FOR WHICH MASK IS ZERO ARE IGNORED.                      
C                                                                       
C     INPUT PARAMETER -                                                 
C        ROOT - IS THE INPUT NODE THAT DEFINES THE COMPONENT.           
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR.                     
C        MASK - SPECIFIES A SECTION SUBGRAPH.                           
C                                                                       
C     OUTPUT PARAMETERS -                                               
C        DEG - ARRAY CONTAINING THE DEGREES OF THE NODES IN             
C              THE COMPONENT.                                           
C        CCSIZE-SIZE OF THE COMPONENT SPECIFED BY MASK AND ROOT         
C                                                                       
C     WORKING PARAMETER -                                               
C        LS - A TEMPORARY VECTOR USED TO STORE THE NODES OF THE         
C               COMPONENT LEVEL BY LEVEL.                               
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  DEGREE ( ROOT, XADJ, ADJNCY, MASK,                    
     1                     DEG, CCSIZE, LS )                            
C                                                                       
C***************************************************************        
C                                                                       
         INTEGER ADJNCY(1), DEG(1), LS(1), MASK(1)                      
         INTEGER XADJ(1), CCSIZE, I, IDEG, J, JSTOP, JSTRT,             
     1           LBEGIN, LVLEND, LVSIZE, NBR, NODE, ROOT                
C                                                                       
C***************************************************************        
C                                                                       
C        -------------------------------------------------              
C        INITIALIZATION ...                                             
C        THE ARRAY XADJ IS USED AS A TEMPORARY MARKER TO                
C        INDICATE WHICH NODES HAVE BEEN CONSIDERED SO FAR.              
C        -------------------------------------------------              
         LS(1) = ROOT                                                   
         XADJ(ROOT) = -XADJ(ROOT)                                       
         LVLEND = 0                                                     
         CCSIZE = 1                                                     
C        -----------------------------------------------------          
C        LBEGIN IS THE POINTER TO THE BEGINNING OF THE CURRENT          
C        LEVEL, AND LVLEND POINTS TO THE END OF THIS LEVEL.             
C        -----------------------------------------------------          
  100    LBEGIN = LVLEND + 1                                            
         LVLEND = CCSIZE                                                
C        -----------------------------------------------                
C        FIND THE DEGREES OF NODES IN THE CURRENT LEVEL,                
C        AND AT THE SAME TIME, GENERATE THE NEXT LEVEL.                 
C        -----------------------------------------------                
         DO 400 I = LBEGIN, LVLEND                                      
            NODE = LS(I)                                                
            JSTRT = -XADJ(NODE)                                         
            JSTOP = IABS(XADJ(NODE + 1)) - 1                            
            IDEG = 0                                                    
            IF ( JSTOP .LT. JSTRT ) GO TO 300                           
               DO 200 J = JSTRT, JSTOP                                  
                  NBR = ADJNCY(J)                                       
                  IF ( MASK(NBR) .EQ. 0 )  GO TO  200                   
                     IDEG = IDEG + 1                                    
                     IF ( XADJ(NBR) .LT. 0 ) GO TO 200                  
                        XADJ(NBR) = -XADJ(NBR)                          
                        CCSIZE = CCSIZE + 1                             
                        LS(CCSIZE) = NBR                                
  200          CONTINUE                                                 
  300       DEG(NODE) = IDEG                                            
  400    CONTINUE                                                       
C        ------------------------------------------                     
C        COMPUTE THE CURRENT LEVEL WIDTH.                               
C        IF IT IS NONZERO , GENERATE ANOTHER LEVEL.                     
C        ------------------------------------------                     
         LVSIZE = CCSIZE - LVLEND                                       
         IF ( LVSIZE .GT. 0 ) GO TO 100                                 
C        ------------------------------------------                     
C        RESET XADJ TO ITS CORRECT SIGN AND RETURN.                     
C        ------------------------------------------                     
         DO 500 I = 1, CCSIZE                                           
            NODE = LS(I)                                                
            XADJ(NODE) = -XADJ(NODE)                                    
  500    CONTINUE                                                       
         RETURN                                                         
      END                                                               
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
C----- SUBROUTINE RCM                                                   
C***************************************************************        
C***************************************************************        
C********     RCM ..... REVERSE CUTHILL-MCKEE ORDERING   *******        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - RCM NUMBERS A CONNECTED COMPONENT SPECIFIED BY          
C        MASK AND ROOT, USING THE RCM ALGORITHM.                        
C        THE NUMBERING IS TO BE STARTED AT THE NODE ROOT.               
C                                                                       
C     INPUT PARAMETERS -                                                
C        ROOT - IS THE NODE THAT DEFINES THE CONNECTED                  
C               COMPONENT AND IT IS USED AS THE STARTING                
C               NODE FOR THE RCM ORDERING.                              
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR FOR                  
C               THE GRAPH.                                              
C                                                                       
C     UPDATED PARAMETERS -                                              
C        MASK - ONLY THOSE NODES WITH NONZERO INPUT MASK                
C               VALUES ARE CONSIDERED BY THE ROUTINE.  THE              
C               NODES NUMBERED BY RCM WILL HAVE THEIR                   
C               MASK VALUES SET TO ZERO.                                
C                                                                       
C     OUTPUT PARAMETERS -                                               
C        PERM - WILL CONTAIN THE RCM ORDERING.                          
C        CCSIZE - IS THE SIZE OF THE CONNECTED COMPONENT                
C               THAT HAS BEEN NUMBERED BY RCM.                          
C                                                                       
C     WORKING PARAMETER -                                               
C        DEG - IS A TEMPORARY VECTOR USED TO HOLD THE DEGREE            
C               OF THE NODES IN THE SECTION GRAPH SPECIFIED             
C               BY MASK AND ROOT.                                       
C                                                                       
C     PROGRAM SUBROUTINES -                                             
C        DEGREE.                                                        
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  RCM ( ROOT, XADJ, ADJNCY, MASK,                       
     1                  PERM, CCSIZE, DEG )                             
C                                                                       
C***************************************************************        
C                                                                       
         INTEGER ADJNCY(1), DEG(1), MASK(1), PERM(1)                    
         INTEGER XADJ(1), CCSIZE, FNBR, I, J, JSTOP,                    
     1           JSTRT, K, L, LBEGIN, LNBR, LPERM,                      
     1           LVLEND, NBR, NODE, ROOT                                
C                                                                       
C***************************************************************        
C                                                                       
C        -------------------------------------                          
C        FIND THE DEGREES OF THE NODES IN THE                           
C        COMPONENT SPECIFIED BY MASK AND ROOT.                          
C        -------------------------------------                          
         CALL  DEGREE ( ROOT, XADJ, ADJNCY, MASK, DEG,                  
     1                  CCSIZE, PERM )                                  
         MASK(ROOT) = 0                                                 
         IF ( CCSIZE .LE. 1 ) RETURN                                    
         LVLEND = 0                                                     
         LNBR = 1                                                       
C        --------------------------------------------                   
C        LBEGIN AND LVLEND POINT TO THE BEGINNING AND                   
C        THE END OF THE CURRENT LEVEL RESPECTIVELY.                     
C        --------------------------------------------                   
  100    LBEGIN = LVLEND + 1                                            
         LVLEND = LNBR                                                  
         DO 600 I = LBEGIN, LVLEND                                      
C           ----------------------------------                          
C           FOR EACH NODE IN CURRENT LEVEL ...                          
C           ----------------------------------                          
            NODE = PERM(I)                                              
            JSTRT = XADJ(NODE)                                          
            JSTOP = XADJ(NODE+1) - 1                                    
C           ------------------------------------------------            
C           FIND THE UNNUMBERED NEIGHBORS OF NODE.                      
C           FNBR AND LNBR POINT TO THE FIRST AND LAST                   
C           UNNUMBERED NEIGHBORS RESPECTIVELY OF THE CURRENT            
C           NODE IN PERM.                                               
C           ------------------------------------------------            
            FNBR = LNBR + 1                                             
            DO 200 J = JSTRT, JSTOP                                     
               NBR = ADJNCY(J)                                          
               IF ( MASK(NBR) .EQ. 0 )  GO TO 200                       
                  LNBR = LNBR + 1                                       
                  MASK(NBR) = 0                                         
                  PERM(LNBR) = NBR                                      
  200       CONTINUE                                                    
            IF ( FNBR .GE. LNBR )  GO TO 600                            
C              ------------------------------------------               
C              SORT THE NEIGHBORS OF NODE IN INCREASING                 
C              ORDER BY DEGREE. LINEAR INSERTION IS USED.               
C              ------------------------------------------               
               K = FNBR                                                 
  300          L = K                                                    
                  K = K + 1                                             
                  NBR = PERM(K)                                         
  400             IF ( L .LT. FNBR )  GO TO 500                         
                     LPERM = PERM(L)                                    
                     IF ( DEG(LPERM) .LE. DEG(NBR) )  GO TO 500         
                        PERM(L+1) = LPERM                               
                        L = L - 1                                       
                        GO TO 400                                       
  500             PERM(L+1) = NBR                                       
                  IF ( K .LT. LNBR )  GO TO 300                         
  600    CONTINUE                                                       
         IF (LNBR .GT. LVLEND) GO TO 100                                
C        ---------------------------------------                        
C        WE NOW HAVE THE CUTHILL MCKEE ORDERING.                        
C        REVERSE IT BELOW ...                                           
C        ---------------------------------------                        
         K = CCSIZE/2                                                   
         L = CCSIZE                                                     
         DO 700 I = 1, K                                                
            LPERM = PERM(L)                                             
            PERM(L) = PERM(I)                                           
            PERM(I) = LPERM                                             
            L = L - 1                                                   
  700    CONTINUE                                                       
         RETURN                                                         
      END                                                               
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
C----- SUBROUTINE GENRCM                                                
C***************************************************************        
C***************************************************************        
C********   GENRCM ..... GENERAL REVERSE CUTHILL MCKEE   *******        
C***************************************************************        
C***************************************************************        
C                                                                       
C     PURPOSE - GENRCM FINDS THE REVERSE CUTHILL-MCKEE                  
C        ORDERING FOR A GENERAL GRAPH. FOR EACH CONNECTED               
C        COMPONENT IN THE GRAPH, GENRCM OBTAINS THE ORDERING            
C        BY CALLING THE SUBROUTINE RCM.                                 
C                                                                       
C     INPUT PARAMETERS -                                                
C        NEQNS - NUMBER OF EQUATIONS                                    
C        (XADJ, ADJNCY) - ARRAY PAIR CONTAINING THE ADJACENCY           
C               STRUCTURE OF THE GRAPH OF THE MATRIX.                   
C                                                                       
C     OUTPUT PARAMETER -                                                
C        PERM - VECTOR THAT CONTAINS THE RCM ORDERING.                  
C                                                                       
C     WORKING PARAMETERS -                                              
C        MASK - IS USED TO MARK VARIABLES THAT HAVE BEEN                
C               NUMBERED DURING THE ORDERING PROCESS. IT IS             
C               INITIALIZED TO 1, AND SET TO ZERO AS EACH NODE          
C               IS NUMBERED.                                            
C        XLS - THE INDEX VECTOR FOR A LEVEL STRUCTURE.  THE             
C               LEVEL STRUCTURE IS STORED IN THE CURRENTLY              
C               UNUSED SPACES IN THE PERMUTATION VECTOR PERM.           
C                                                                       
C     PROGRAM SUBROUTINES -                                             
C        FNROOT, RCM.                                                   
C                                                                       
C***************************************************************        
C                                                                       
      SUBROUTINE  GENRCM ( NEQNS, XADJ, ADJNCY, PERM, MASK, XLS )       
C                                                                       
C***************************************************************        
C                                                                       
         INTEGER ADJNCY(1), MASK(1), PERM(1), XLS(1)                    
         INTEGER XADJ(1), CCSIZE, I, NEQNS, NLVL,                       
     1           NUM, ROOT                                              
C                                                                       
C***************************************************************        
C                                                                       
         DO 100 I = 1, NEQNS                                            
            MASK(I) = 1                                                 
  100    CONTINUE                                                       
         NUM = 1                                                        
         DO 200 I = 1, NEQNS                                            
C           ---------------------------------------                     
C           FOR EACH MASKED CONNECTED COMPONENT ...                     
C           ---------------------------------------                     
            IF (MASK(I) .EQ. 0) GO TO 200                               
               ROOT = I                                                 
C              -----------------------------------------                
C              FIRST FIND A PSEUDO-PERIPHERAL NODE ROOT.                
C              NOTE THAT THE LEVEL STRUCTURE FOUND BY                   
C              FNROOT IS STORED STARTING AT PERM(NUM).                  
C              THEN RCM IS CALLED TO ORDER THE COMPONENT                
C              USING ROOT AS THE STARTING NODE.                         
C              -----------------------------------------                
               CALL  FNROOT ( ROOT, XADJ, ADJNCY, MASK,                 
     1                        NLVL, XLS, PERM(NUM) )                    
               CALL     RCM ( ROOT, XADJ, ADJNCY, MASK,                 
     1                        PERM(NUM), CCSIZE, XLS )                  
               NUM = NUM + CCSIZE                                       
               IF (NUM .GT. NEQNS) RETURN                               
  200    CONTINUE                                                       
         RETURN                                                         
      END                                                               
