      SUBROUTINE DRV000(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
* MODULE TO FIND A ZERO USING BRENT'S METHOD WITH
*                             GIVEN FUNCTION VALUES.
*
* INPUT/OUTPUT PARAMETERS:
*  NENTRY : NUMBER OF LINKED LISTS AND FILES USED BY THE MODULE.
*  HENTRY : CHARACTER*12 NAME OF EACH LINKED LIST OR FILE.
*  IENTRY : =0 CLE-2000 VARIABLE; =1 LINKED LIST; =2 XSM FILE;
*           =3 SEQUENTIAL BINARY FILE; =4 SEQUENTIAL ASCII FILE;
*           =5 DIRECT ACCESS FILE.
*  JENTRY : =0 THE LINKED LIST OR FILE IS CREATED.
*           =1 THE LINKED LIST OR FILE IS OPEN FOR MODIFICATIONS;
*           =2 THE LINKED LIST OR FILE IS OPEN IN READ-ONLY MODE.
*  KENTRY : =FILE UNIT NUMBER; =LINKED LIST ADDRESS OTHERWISE.
*           DIMENSION HENTRY(NENTRY),IENTRY(NENTRY),JENTRY(NENTRY),
*           KENTRY(NENTRY)
*
*-------------------------------------- AUTHOR: R. ROY ;   29/11/94 ---
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NENTRY
      CHARACTER        HENTRY(NENTRY)*12
      INTEGER          IENTRY(NENTRY), JENTRY(NENTRY)
      TYPE(C_PTR)      KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      INTEGER          NITMA, ITYP, I
      CHARACTER        TEXT12*12, SGNTUR*12
      LOGICAL          LSTART, LCONV
      DOUBLE PRECISION DFLOTT
      INTEGER          ITER,ITMAX,IPRT, ISGNTR(3),ITMD
      INTEGER          ITERV(3), ICONV
      REAL             A,B,C, D,E,   FA,FB,FC, P,Q,R,S, TOL1,XM,TOL
      REAL             X(3),  DE(2), Y(3),     PQRS(4), ATOL(3)
      REAL             FLOTT
      REAL             EPM,TOLD,Z0,ZH,Z1,Z2,Z3,ZBESTM
      TYPE(C_PTR)      IPL0
      PARAMETER       (EPM=3.E-8,TOLD=1.E-5,ITMD=100)
      PARAMETER       (Z0=0.0,ZH=0.5,Z1=1.0,Z2=2.0,Z3=3.0)
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.1) CALL XABORT('DRV000: ONLY ONE ENTRY EXPECTED.')
      TEXT12=HENTRY(1)
      IF(IENTRY(1).NE.1) CALL XABORT('DRV000: LHS L_0 OBJECT EXPECTED ('
     > //TEXT12//').')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('DRV000: LH'
     > //'S L_0 OBJECT IN CREATE OR MODIFICATION MODE EXPECTED.')
*
      LSTART= JENTRY(1).EQ.0
      IPL0= KENTRY(1)
      IF( LSTART )THEN
*        DEFINE ALL TEMP VARIABLES
         D= 0.0
         E= 0.0
         P= 0.0
         Q= 0.0
         R= 0.0
         S= 0.0
*
         LCONV= .FALSE.
         TOL=   TOLD
         ITMAX= ITMD
         ITER=  0
         IPRT=  0
  10     CONTINUE
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF( ITYP.NE.3 )
     >    CALL XABORT('DRV000: KEYWORDS *TOL*, *POINT*... EXPECTED.')
         IF( TEXT12.EQ.'TOL' )THEN
           CALL REDGET(ITYP,NITMA,TOL,TEXT12,DFLOTT)
           IF( ITYP.NE.2 )
     >      CALL XABORT('DRV000: A REAL TOLERANCE *TOL* IS EXPECTED.')
           IF( TOL.LT.1.E-7 )
     >      CALL XABORT('DRV000: TOLERANCE .LT. 1.E-7.')
           GO TO 10
         ELSEIF( TEXT12.EQ.'ITMAX' )THEN
           CALL REDGET(ITYP,ITMAX,FLOTT,TEXT12,DFLOTT)
           IF( ITYP.NE.1 )
     >      CALL XABORT('DRV000: AN INTEGER *ITMAX* IS EXPECTED.')
           GO TO 10
         ELSEIF( TEXT12.EQ.'DEBUG' )THEN
           IPRT= 1
           GO TO 10
         ELSEIF( TEXT12.EQ.'POINT' )THEN
           CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
           IF(ITYP.NE.3.OR.TEXT12.NE.'X')
     >       CALL XABORT('DRV000: *X* KEYWORD EXPECTED.')
           CALL REDGET(ITYP,NITMA,A    ,TEXT12,DFLOTT)
           IF(ITYP.NE.2 ) CALL XABORT('DRV000: AFTER *X*,'
     >      // ' A REAL IS EXPECTED.')
           CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
           IF(ITYP.NE.3.OR.TEXT12.NE.'Y')
     >       CALL XABORT('DRV000: *Y* KEYWORD EXPECTED.')
           CALL REDGET(ITYP,NITMA,FA   ,TEXT12,DFLOTT)
           IF(ITYP.NE.2 ) CALL XABORT('DRV000: AFTER *Y*,'
     >      // ' A REAL IS EXPECTED.')
         ELSE
          CALL XABORT('DRV000: KEYWORDS *TOL* OR *POINT* EXPECTED.')
         ENDIF
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.3.OR.TEXT12.NE.'POINT')
     >     CALL XABORT('DRV000: ONCE MORE, *POINT* KEYWORD EXPECTED.')
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.3.OR.TEXT12.NE.'X')
     >     CALL XABORT('DRV000: *X* KEYWORD EXPECTED.')
         CALL REDGET(ITYP,NITMA,B    ,TEXT12,DFLOTT)
         IF(ITYP.NE.2 ) CALL XABORT('DRV000: AFTER *X*,'
     >    // ' A REAL IS EXPECTED.')
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.3.OR.TEXT12.NE.'Y')
     >     CALL XABORT('DRV000: *Y* KEYWORD EXPECTED.')
         CALL REDGET(ITYP,NITMA,FB   ,TEXT12,DFLOTT)
         IF(ITYP.NE.2 ) CALL XABORT('DRV000: AFTER *Y*,'
     >    // ' A REAL IS EXPECTED.')
         SGNTUR='L_0'
         READ(SGNTUR,'(3A4)') ISGNTR
         CALL LCMSIX(IPL0,' ',0)
*
*        PUT SIGNATURE
         CALL LCMPUT(IPL0,'SIGNATURE',3,3,ISGNTR)
*
*        PUT CONVERGENCE FLAG
         ICONV=-1
         CALL LCMPUT(IPL0,'ICONV',1,1,ICONV)
      ELSE
*
         CALL LCMSIX(IPL0,' ',0)
*
*        VERIFY SIGNATURE
         CALL LCMGET(IPL0,'SIGNATURE',ISGNTR)
         WRITE(SGNTUR,'(3A4)') (ISGNTR(I),I=1,3)
         IF(SGNTUR.NE.'L_0')
     >     CALL XABORT('DRV000: L_0 OBJECT IS EXPECTED')
*
         CALL LCMGET(IPL0,'ICONV',ICONV)
         LCONV= ICONV.EQ.+1
*
*        NOTIFY USER IF ALREADY CONVERVED
         IF( LCONV )
     >     CALL XABORT('DRV000: PROCESS IS ALREADY CONVERGED')
*
*        GET L_0 OBJECT VALUES
         CALL LCMGET(IPL0,'X',X)
         A=X(1)
         B=X(2)
         C=X(3)
         CALL LCMGET(IPL0,'DE',DE)
         D=DE(1)
         E=DE(2)
         CALL LCMGET(IPL0,'Y',Y)
         FA=Y(1)
         FB=Y(2)
         FC=Y(3)
         CALL LCMGET(IPL0,'PQRS',PQRS)
         P=PQRS(1)
         Q=PQRS(2)
         R=PQRS(3)
         S=PQRS(4)
         CALL LCMGET(IPL0,'ATOL',ATOL)
         TOL=ATOL(1)
         XM=ATOL(2)
         TOL1=ATOL(3)
         CALL LCMGET(IPL0,'ITERV',ITERV)
         ITER=ITERV(1)
         ITMAX=ITERV(2)
         IPRT=ITERV(3)
*
*        GET NEW *Y* VALUE
         CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.3.OR.TEXT12.NE.'Y')
     >     CALL XABORT('DRV000: *Y* KEYWORD EXPECTED.')
         CALL REDGET(ITYP,NITMA,FB   ,TEXT12,DFLOTT)
         IF(ITYP.NE.2 )
     >     CALL XABORT('DRV000: AFTER *Y*, A REAL IS EXPECTED.')
      ENDIF
*
*       METHOD:  BRENT'S METHOD FOR FINDING ZEROS.
*                FIRST INTERVAL MUST BE BRACKETED: FA*FB < 0.
*
*        INPUT:  ITER=   NUMBER OF ITERATIONS (0 AT START)
*                TOL=    TOLERANCE FOR ZERO FINDING
*                (A,FA)= FIRST  POINT
*                (B,FB)= SECOND POINT
*
*       OUTPUT:  ZBESTM= ESTIMATION OF NEXT ZERO
*
      IF( ITER.EQ.0 )THEN
         IF((FA.GT.Z0.AND.FB.GT.Z0).OR.(FA.LT.Z0.AND.FB.LT.Z0))
     >     CALL XABORT(' DRV000: ROOT MUST BE BRACKETED')
         C=B
         FC=FB
      ENDIF
      IF((FB.GT.Z0.AND.FC.GT.Z0).OR.(FB.LT.Z0.AND.FC.LT.Z0))THEN
         C=A
         FC=FA
         D=B-A
         E=D
      ENDIF
      IF(ABS(FC).LT.ABS(FB)) THEN
         A=B
         B=C
         C=A
         FA=FB
         FB=FC
         FC=FA
      ENDIF
      TOL1=Z2*EPM*ABS(B)+ ZH*TOL
      XM=ZH*(C-B)
      IF(ABS(XM).LE.TOL1 .OR. FB.EQ.Z0)THEN
         ZBESTM=B
         LCONV= .TRUE.
         GO TO 20
      ENDIF
      IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
         S=FB/FA
         IF(A.EQ.C) THEN
            P=Z2*XM*S
            Q=Z1-S
         ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(Z2*XM*Q*(Q-R)-(B-A)*(R-Z1))
            Q=(Q-Z1)*(R-Z1)*(S-Z1)
         ENDIF
         IF(P.GT.Z0) Q=-Q
         P=ABS(P)
         IF(Z2*P .LT. MIN(Z3*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
         ELSE
            D=XM
            E=D
         ENDIF
      ELSE
         D=XM
         E=D
      ENDIF
      A=B
      FA=FB
      IF(ABS(D) .GT. TOL1) THEN
         B=B+D
      ELSE
         B=B+SIGN(TOL1,XM)
      ENDIF
      ZBESTM=B
      ITER= ITER + 1
      IF( ITER.GT.ITMAX )
     >   CALL XABORT('DRV000: MAX NUMBER OF ITERATIONS REACHED.')
*
   20 CONTINUE
*
*     PUT L_0 OBJECT VALUES
      X(1)=A
      X(2)=B
      X(3)=C
      CALL LCMPUT(IPL0,'X',3,2,X)
      DE(1)=D
      DE(2)=E
      CALL LCMPUT(IPL0,'DE',2,2,DE)
      Y(1)=FA
      Y(2)=FB
      Y(3)=FC
      CALL LCMPUT(IPL0,'Y',3,2,Y)
      PQRS(1)=P
      PQRS(2)=Q
      PQRS(3)=R
      PQRS(4)=S
      CALL LCMPUT(IPL0,'PQRS',4,2,PQRS)
      ATOL(1)=TOL
      ATOL(2)=XM
      ATOL(3)=TOL1
      CALL LCMPUT(IPL0,'ATOL',3,2,ATOL)
      ITERV(1)=ITER
      ITERV(2)=ITMAX
      ITERV(3)=IPRT
      CALL LCMPUT(IPL0,'ITERV',3,1,ITERV)
      IF( LCONV )THEN
         ICONV=+1
      ELSE
         ICONV=-1
      ENDIF
*
*     SAVE CONVERGENCE FLAG
      CALL LCMPUT(IPL0,'ICONV',1,1,ICONV)
      CALL LCMSIX(IPL0,' ',0)
      IF( IPRT.EQ.1 )THEN
        WRITE(6,*) 'DEBUG:  A=', A,'  B=', B,'  C=', C
        WRITE(6,*) 'DEBUG: FA=',FA,' FB=',FB,' FC=',FC
        WRITE(6,*) 'DEBUG:  D=', D,'  E=', E
        WRITE(6,*) 'DEBUG:  P=', P,'  Q=', Q,'  R=',R,'  S=',S
        WRITE(6,*) 'DEBUG: TOL1=',TOL1,' XM=',XM,' TOL=',TOL
        WRITE(6,*) 'DEBUG: ITER=',ITER,' ITMAX=',ITMAX
        WRITE(6,*) 'DEBUG: ICONV=',ICONV
      ENDIF
*
*     NOW, RETURN BACK LOGICAL VALUE AND ZERO ESTIMATE
      CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
      ITYP= -ITYP
      IF( ITYP.NE.5 )THEN
         CALL XABORT('DRV000: MUST WRITE LOGICAL FLAG INTO >>.<<')
      ENDIF
      CALL REDPUT(ITYP,ICONV,FLOTT,TEXT12,DFLOTT)
      CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
      ITYP= -ITYP
      IF( ITYP.NE.2 )THEN
         CALL XABORT('DRV000: MUST WRITE REAL ZERO INTO >>.<<')
      ENDIF
      CALL REDPUT(ITYP,NITMA,ZBESTM,TEXT12,DFLOTT)
*
      CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(ITYP.NE.3.OR.TEXT12.NE.';')
     >  CALL XABORT('DRV000: *;* IS EXPECTED FOR ENDING THE SENTENCE.')
      RETURN
      END
