c============================================================================================
c    DDFSA - A Distributed Derivative-free Simulated Annealing Method for
c    bound constrained global optimization
c    Copyright (C) 2011  G.Liuzzi, S.Lucidi, V.Piccialli
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
c    G. Liuzzi, S. Lucidi, V. Piccialli, A. Sotgiu. A magnetic resonance device designed via 
c    global optimization techniques, Mathematical Programming, 101: 339-364 (2004)
c
c============================================================================================
c************************************************************
c*
c***********************************************************
      SUBROUTINE MAINBOX1(N,X,D,D1,Z,Z1,XOLD,NUM_ITER,
     *DOLDALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,FFSTOP,BL,BU)

      IMPLICIT NONE
      INTEGER N,I,J,I_CORR,NUM_FUNCT,NUM_ITER
      INTEGER NUM_FAL,NUM_SUCC,ISTOP
      INTEGER IPRINT
      INTEGER IMIN,IMAX,IMINALFA,IMAXALFA

      REAL*8 DCONV(N),DNR,FFM
      REAL*8 X(N),Z(N),Z1(N),D(N),D1(N),XOLD(N)
      REAL*8 DOLDALFA(N),ALFA,BL(N),BU(N)
      REAL*8 F,FZ 
      REAL*8 DOLDALFAMEDIO,DALFAMAX
      REAL*8 FMIN,FMAX,ALFA0,DOLDALFAMIN,DOLDALFAMAX
      REAL*8 RAPALFA,FINIT(N,2)

C     VETTORE DEI VALORI DI F SUI PUNTI DI UN SIMPLESSO N+1 DIM.

      REAL*8 FSTOP(N+1),XFSTOP(N,N+1),FFSTOP,FMM

c     num_fal rappresenta il numero di fallimenti consecutivi

c     i_corr rappresenta l'indice della direzione corrente

      COMMON /NUM/F
      COMMON /NUMNEW/NUM_FUNCT

c     inizializzazione

      NUM_SUCC=0

      NUM_FUNCT = 0
      NUM_ITER = 0 
      NUM_FAL=0
      ISTOP = 0

      I_CORR=1

      ALFA0=1.D0

      DO I=1,N

C      
        IF(IPRINT.GE.2) THEN
          WRITE(*,*) ' ALFAiniz(',I,')=',DOLDALFA(I)
          WRITE(1,*) ' ALFAiniz(',I,')=',DOLDALFA(I)
        ENDIF
C
      END DO

c     scelta iniziale delle direzioni

      
      CALL FUNCT(X,N,F)
!	WRITE(*,*)'MAINBOX1: ',X,F
!	WRITE(1,*)'MAINBOX1: ',X,F
	NUM_FUNCT=NUM_FUNCT+1

      FSTOP(1)=F

      DO I=1,N
        XFSTOP(I,1)=X(I)
	  Z(I)=X(I)
      END DO

      IF(IPRINT.GE.2) THEN
        WRITE(*,*) ' ----------------------------------'
        WRITE(1,*) ' ----------------------------------'
        WRITE(*,*) ' Finiz =',F
        WRITE(1,*) ' Finiz =',F
        DO I=1,N
          WRITE(*,*) ' Xiniz(',I,')=',X(I)
          WRITE(1,*) ' Xiniz(',I,')=',X(I)
        ENDDO
      ENDIF
 
      NUM_FAL=0


C---------------------------
C     CICLO PRINCIPALE
C---------------------------

  1   CONTINUE

      IF(I_CORR.EQ.1) THEN
           DO I=1,N
                DCONV(I)=-D(I)
           END DO
      ENDIF

	ISTOP = 0

      DALFAMAX=MAXVAL(DOLDALFA)

      IF(IPRINT.GE.1) THEN
        WRITE(*,*) '----------------------------------------------'
        WRITE(1,*) '----------------------------------------------'
        WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
      ENDIF

      IF (ISTOP.EQ.1) THEN
!        WRITE(2,50)'&',N,'&',NUM_FUNCT,'&',F,'\\'  
   50   FORMAT(1X,a3,i5,a3,i5,a3,d13.5,a6)
!        WRITE(*,*) '----------------------------------------------'
!        WRITE(1,*) '----------------------------------------------'
!        WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
!        WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        RETURN
      END IF

C-------------------------------------
C    CAMPIONAMENTO LUNGO ASSE I_CORR
C-------------------------------------
 
C      CALL LINESEARCH_INVERTI(N,X,F,D,ALFA,DOLDALFA,Z,FZ,
C     *I_CORR,NUM_FAL,DALFAMAX,IPRINT)
      CALL linesearchbox_cont(n,x,f,d,alfa,doldalfa,z,fz,
     *i_corr,num_fal,dalfamax,iprint,bl,bu)

      IF(DABS(ALFA).GE.1.D-24) THEN

          X(I_CORR) = X(I_CORR)+ALFA*D(I_CORR)
          Z(I_CORR) = X(I_CORR)

          F=FZ
          DO I=N+1,2,-1
             FSTOP(I)=FSTOP(I-1)
             DO J=1,N
                XFSTOP(J,I)=XFSTOP(J,I-1)
             END DO
          ENDDO
          FSTOP(1)=F
          DO J=1,N
             XFSTOP(J,1)=X(J)
          END DO            
          NUM_FAL=0
          NUM_ITER=NUM_ITER+1
          NUM_SUCC=NUM_SUCC+1
    
          IF(IPRINT.GE.1) THEN
             WRITE(*,*) ' F =',F
             WRITE(1,*) ' F =',F
          ENDIF

          IF(IPRINT.GE.2) THEN
	       DO I=1,N
                WRITE(*,*) ' X(',I,')=',X(I)
                WRITE(1,*) ' X(',I,')=',X(I)
             ENDDO
          ENDIF
      
	ELSE
      
	    DO I=N+1,2,-1
             FSTOP(I)=FSTOP(I-1)
             DO J=1,N
                XFSTOP(J,I)=XFSTOP(J,I-1)
             END DO
          ENDDO
          FSTOP(1)=FZ
          DO J=1,N
             XFSTOP(J,1)=Z(J)  
          END DO       
          NUM_FAL=NUM_FAL+1
          NUM_ITER=NUM_ITER+1
          Z(I_CORR) = X(I_CORR)
		      
	END IF

      IF(I_CORR.LT.N) THEN

          I_CORR=I_CORR+1

      ELSE

          I_CORR=1
	  ISTOP = 1

	    IF (ISTOP.EQ.1) THEN
!             WRITE(2,50)'&',N,'&',NUM_FUNCT,'&',F,'\\'  
!             WRITE(*,*) '----------------------------------------------'
!             WRITE(1,*) '----------------------------------------------'
!             WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
!             WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
              FFM=0.D0
		    DO I=1,N+1
			  FFM=FFM+FSTOP(I)
			ENDDO
			FFM=FFM/DFLOAT((N+1))

			FFSTOP=0.D0
			DO I=1,N+1
			  FFSTOP=FFSTOP+(FSTOP(I)-FFM)*(FSTOP(I)-FFM)
			ENDDO
		!	write(*,*) ' ffstop =',ffstop,  ' dfloat =',dfloat(n+1),' ffm =',ffm
			FFSTOP=DSQRT(FFSTOP/DFLOAT(N+1))

             RETURN
          END IF

      END IF 

      GO TO 1

      END
        
c************************************************************
c*
c***********************************************************
      SUBROUTINE MAINBOX2(N,X,D,D1,Z,Z1,XOLD,NUM_ITER,DOLDALFA,
     *ALFA_STOP,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,FFSTOP,
     *FFSTOPTOL,FOTT,ALFAOTT,XOTT,WKS,IDIMWKS,ISTOP,MASCH,MAXNF,fglob,
     *BL,BU)
	
      IMPLICIT NONE
	INCLUDE 'TYPEDECL.FI'

      INTEGER N,I,J,I_CORR,NUM_FUNCT,NUM_ITER,IDIMWKS,MAXNF
      INTEGER NUM_FAL,NUM_SUCC,ISTOP
      INTEGER IPRINT
	INTEGER IMIN,IMAX,IMINALFA,IMAXALFA

      REAL*8 DCONV(N),DNR,ALFA_STOP,fglob,BL(N),BU(N)
      REAL*8 X(N),Z(N),Z1(N),D(N),D1(N),XOLD(N),XOTT(N)
      REAL*8 DOLDALFA(N),ALFA,FOTT,ALFAOTT
      REAL*8 F,FZ 
      REAL*8 DOLDALFAMEDIO,DALFAMAX
      REAL*8 FMIN,FMAX,ALFA0,DOLDALFAMIN,DOLDALFAMAX
      REAL*8 RAPALFA,FINIT(N,2)
	TYPE(COMP_WKS_TYP)	:: WKS(IDIMWKS)

C     VETTORE DEI VALORI DI F SUI PUNTI DI UN SIMPLESSO N+1 DIM.

      REAL*8 FSTOP(N+1),XFSTOP(N,N+1),FFSTOP,FFSTOPTOL
	LOGICAL MASCH(IDIMWKS)

c     num_fal rappresenta il numero di fallimenti consecutivi

c     i_corr rappresenta l'indice della direzione corrente

      COMMON /NUM/F
      COMMON /NUMNEW/NUM_FUNCT

c     inizializzazione

      NUM_SUCC=0

      NUM_FUNCT = 0
      NUM_ITER = 0 
      NUM_FAL=0
      ISTOP = 0

      I_CORR=1

      ALFA0=1.D0

c      DO I=1,N

c        DOLDALFA(I)=1.D0

c        DOLDALFA(I)=DMAX1(5.D-1,DMIN1(2.D+0,5.D-1*DABS(X(I))))
c        DOLDALFA(I)=1.0D0*DMAX1(1.D-3,DMIN1(1.D+3,0.5D0*DABS(X(I))))
      
c        IF(IPRINT.GE.2) THEN
c          WRITE(*,*) ' ALFAiniz(',I,')=',DOLDALFA(I)
c          WRITE(1,*) ' ALFAiniz(',I,')=',DOLDALFA(I)
c        ENDIF

C      END DO

c     scelta iniziale delle direzioni

	DALFAMAX = MAXVAL(DOLDALFA)

      DO I=1,N      
        D(I)=1.D0 
      END DO
       
      CALL FUNCT(X,N,F)
	NUM_FUNCT=NUM_FUNCT+1

      FSTOP(1)=F

      DO I=1,N
        XFSTOP(I,1)=X(I)
	  Z(I)=X(I)
      END DO

      IF(IPRINT.GE.2) THEN
        WRITE(*,*) ' ----------------------------------'
        WRITE(1,*) ' ----------------------------------'
        WRITE(*,*) ' Finiz =',F
        WRITE(1,*) ' Finiz =',F
        DO I=1,N
          WRITE(*,*) ' Xiniz(',I,')=',X(I)
          WRITE(1,*) ' Xiniz(',I,')=',X(I)
        ENDDO
      ENDIF
 
      NUM_FAL=0

C---------------------------
C     CICLO PRINCIPALE
C---------------------------

  1   CONTINUE

      IF(I_CORR.EQ.1) THEN
           DO I=1,N
                DCONV(I)=-D(I)
           END DO
      ENDIF

      IF(NUM_ITER.GE.N) THEN
        CALL STOP(N,DOLDALFA,ISTOP,DALFAMAX,NUM_FUNCT,FSTOP,
     *			    ALFA_STOP,FFSTOPTOL,FFSTOP,FOTT,X,XOTT,ALFAOTT,
     *				WKS,IDIMWKS,MASCH,MAXNF,fglob)
      ENDIF

      IF(IPRINT.GE.1) THEN
        WRITE(*,*) '----------------------------------------------'
        WRITE(1,*) '----------------------------------------------'
        WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
      ENDIF

      IF (ISTOP.GE.1) THEN
!        WRITE(2,50)'&',N,'&',NUM_FUNCT,'&',F,'\\'  
   50   FORMAT(1X,a3,i5,a3,i5,a3,d13.5,a6)
!        WRITE(*,*) '----------------------------------------------'
!        WRITE(1,*) '----------------------------------------------'
!        WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
!        WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        RETURN
      END IF

C-------------------------------------
C    CAMPIONAMENTO LUNGO ASSE I_CORR
C-------------------------------------
 

C      CALL LINESEARCH_INVERTI(N,X,F,D,ALFA,DOLDALFA,Z,FZ,
C     *I_CORR,NUM_FAL,DALFAMAX,IPRINT)
      CALL linesearchbox_cont(n,x,f,d,alfa,doldalfa,z,fz,
     *i_corr,num_fal,dalfamax,iprint,bl,bu)

      IF(DABS(ALFA).GE.1.D-24) THEN

          X(I_CORR) = X(I_CORR)+ALFA*D(I_CORR)
          Z(I_CORR) = X(I_CORR)

          F=FZ
          DO I=N+1,2,-1
             FSTOP(I)=FSTOP(I-1)
             DO J=1,N
                XFSTOP(J,I)=XFSTOP(J,I-1)
             END DO
          ENDDO
          FSTOP(1)=F
          DO J=1,N
             XFSTOP(J,1)=X(J)
          END DO            
          NUM_FAL=0
          NUM_ITER=NUM_ITER+1
          NUM_SUCC=NUM_SUCC+1
    
          IF(IPRINT.GE.1) THEN
             WRITE(*,*) ' F =',F
             WRITE(1,*) ' F =',F
          ENDIF

          IF(IPRINT.GE.2) THEN
	       DO I=1,N
                WRITE(*,*) ' X(',I,')=',X(I)
                WRITE(1,*) ' X(',I,')=',X(I)
             ENDDO
          ENDIF
      
	ELSE
      
	    DO I=N+1,2,-1
             FSTOP(I)=FSTOP(I-1)
             DO J=1,N
                XFSTOP(J,I)=XFSTOP(J,I-1)
             END DO
          ENDDO
          FSTOP(1)=FZ
          DO J=1,N
             XFSTOP(J,1)=Z(J)  
          END DO       
          NUM_FAL=NUM_FAL+1
          NUM_ITER=NUM_ITER+1
          Z(I_CORR) = X(I_CORR)
		      
	END IF

      IF(I_CORR.LT.N) THEN

          I_CORR=I_CORR+1

      ELSE

          I_CORR=1

      END IF 

      GO TO 1

      END


C     #######################################################

      SUBROUTINE STOP(N,DOLDALFA,ISTOP,DALFAMAX,NUM_FUNCT,FSTOP,
     *			    ALFA_STOP,FFSTOPTOL,FFSTOP,FOTT,X,
     *		      XOTT,ALFAOTT,WKS,IDIMWKS,MASCH,MAXNF,fglob)
C	USE IMSLF90

      IMPLICIT NONE

	INCLUDE 'TYPEDECL.FI'
      
      INTEGER N,ISTOP,I,NUM_FUNCT,IDIMWKS,J,IMAX,MAXNF
      REAL*8 DOLDALFA(N),DALFAMAX,FSTOP(N+1),FFSTOP,FFM,ALFA_STOP,F
     *	,XOTT(N),X(N)
	REAL*8 FFSTOPTOL,FOTT,ALFAOTT,fglob
	TYPE(COMP_WKS_TYP)	:: WKS(IDIMWKS)

	LOGICAL NORMADUE,L_F_DISTANTI,NON_MINIMIZZ,L_F_VICINI
        LOGICAL MASCH(IDIMWKS)

	COMMON/CALFANR2/NORMADUE
      COMMON/NUM/F

      ISTOP=0
	NON_MINIMIZZ = .FALSE.

!	WRITE(*,*)'STO DENTRO STOP ! ALFA_STOP = ',ALFA_STOP
!	WRITE(1,*)'STO DENTRO STOP ! ALFA_STOP = ',ALFA_STOP
CC	PAUSE
      DALFAMAX=DOLDALFA(1)
      DO I=1,N
        IF(DOLDALFA(I).GT.DALFAMAX) THEN
          DALFAMAX=DOLDALFA(I)
        END IF
      END DO
     
C      FFM=1.D+30
C      DO I=1,N+1
C        IF(FSTOP(I).LT.FFM) FFM=FSTOP(I)
C      ENDDO

      FFM=0.D0
      DO I=1,N+1
        FFM=FFM+FSTOP(I)
      ENDDO
      FFM=FFM/DFLOAT((N+1))

      FFSTOP=0.D0
      DO I=1,N+1
        FFSTOP=FFSTOP+(FSTOP(I)-FFM)*(FSTOP(I)-FFM)
      ENDDO
!	write(*,*) ' ffstop =',ffstop,  ' dfloat =',dfloat(n+1),' ffm =',ffm
      FFSTOP=DSQRT(FFSTOP/DFLOAT(N+1))

	IF(NORMADUE) THEN
		IF(DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA)).LE.ALFA_STOP) THEN
		  ISTOP = 1
		END IF
	ELSE
		IF(DALFAMAX.LE.ALFA_STOP) THEN
		  ISTOP = 1
		END IF
	ENDIF

!--------------------------------------
!		calcola l'indice imax della
!		componente di wks contenente
!		il valore di fob piu' alto
!--------------------------------------
		IMAX = 1
		DO I = 2,IDIMWKS
			IF((WKS(I)%FOB > WKS(IMAX)%FOB).AND.(.NOT.MASCH(I))) THEN
				IMAX = I
			ENDIF
		ENDDO

	L_F_DISTANTI = ((F-WKS(IMAX)%FOB)/DMAX1(1.D-8,
     *					DMIN1(DABS(F),DABS(WKS(IMAX)%FOB))) ).GE.1.D-6
	L_F_DISTANTI = L_F_DISTANTI.AND.
     *			   (DALFAMAX.LE.1.D+0*WKS(IMAX)%ALFAMAX)

	IF (L_F_DISTANTI) THEN
		ISTOP = 1
		RETURN
	ENDIF

	DO J = 1,IDIMWKS

		IF(F > FOTT) THEN

			L_F_VICINI = (DABS(F-WKS(J)%FOB)/(DMAX1(1.D-8,
     *			DABS(DMAX1(F,WKS(J)%FOB) ))) ).LE.1.D-1
			L_F_VICINI = L_F_VICINI.AND.((DSQRT(													
     *			DOT_PRODUCT(X-WKS(J)%X,X-WKS(J)%X)) <= 1.D-1)
     *			.OR.(DALFAMAX.LE.1.D+1*WKS(J)%ALFAMAX))
!			L_F_VICINI = L_F_VICINI.AND.(DSQRT(													
!     *			DOT_PRODUCT(X-WKS(J)%X,X-WKS(J)%X)) <= 1.D-1)

			IF (L_F_VICINI) THEN 
			
!				ISTOP=1
!				RETURN

			ENDIF
		ENDIF
	ENDDO

	L_F_DISTANTI = ((F-FOTT)/DMAX1(1.D-3,DMIN1(DABS(F),DABS(FOTT))
     *										  ) ).GE.1.D+1
	L_F_DISTANTI = L_F_DISTANTI.AND.
     *			   (DALFAMAX.LE.1.D+1*ALFA_STOP)

	L_F_VICINI = ((F-FOTT)/DMAX1(1.D-8,DABS(DMAX1(F,FOTT)) 
     *										  ) ).LE.1.D-1

C	L_F_VICINI = L_F_VICINI.AND.((DSQRT(DOT_PRODUCT(X-XOTT,X-XOTT))
C     *			<= 1.D-1))

	L_F_VICINI = L_F_VICINI.AND.((DSQRT(DOT_PRODUCT(X-XOTT,X-XOTT))
     *			<= 1.D-1) .OR. (DALFAMAX <= 0.5D+0*ALFAOTT))

c	L_F_DISTANTI = ((F-FOTT)/DMAX1(1.D-16,DMAX1( DABS(FOTT),
c     *										DABS(F) ) ) ).GE.5.D-1
c	NON_MINIMIZZ = L_F_DISTANTI.AND.
c     *			   (DALFAMAX.LE.DMIN1(1.0D-0,1.D+2*ALFA_STOP))

!	NON_MINIMIZZ = L_F_DISTANTI.AND.(DALFAMAX.LE.1.D-3)
	
	NON_MINIMIZZ = L_F_DISTANTI .OR. L_F_VICINI

	IF(NON_MINIMIZZ) ISTOP = 1

C      IF(FFSTOP.LE.FFSTOPTOL) THEN
C        ISTOP = 1
C      END IF
      IF(NUM_FUNCT.GT.MAXNF) THEN
        ISTOP = 1
      END IF
	
	if((f-fglob)/max(1.0d0,abs(fglob))<1.0d-4) then
		istop = 1
	endif
	if((fott-fglob)/max(1.0d0,abs(fglob))<1.0d-4) then
		istop = 1
	endif

      RETURN

      END

C *****************************************************************
C *****************************************************************
C *****************************************************************

      SUBROUTINE LINESEARCH_INVERTI(N,X,F,D,ALFA,DOLDALFA,Z,FZ,
     *I_CORR,NUM_FAL,DALFAMAX,IPRINT)
     
      IMPLICIT NONE

      INTEGER N,I_CORR,NUM_FUNCT
      INTEGER I,J,L,LL
      INTEGER NUM_FAL
      INTEGER IPRINT
      REAL*8 X(N),D(N),DOLDALFA(N),Z(N),Z1(N)
      REAL*8 F,ALFA,FZ,GAMMA,DNR
      REAL*8 DELTA,DELTA1,FPAR,FZDELTA
      REAL*8 ALFAEX,FMIN,ALFAMIN 
      REAL*8 DALFAMAX,FCOMMON
      REAL*8 ALFAOLD,FZOLD
	LOGICAL PRIMA_VOLTA

      COMMON /NUM/FCOMMON
	COMMON /NUMNEW/NUM_FUNCT

	PRIMA_VOLTA = .TRUE.
      GAMMA = 1.d-6
      DELTA =0.5D0
      DELTA1 =0.5D0
      DNR=1.D0


c     indice della direzione corrente

      J=I_CORR

      IF(IPRINT.GE.1) THEN
         WRITE(*,*) ' J=',J,'  D(J)=',D(J),'  DOLDALFA=',DOLDALFA(J)
         WRITE(1,*) ' J=',J,'  D(J)=',D(J),'  DOLDALFA=',DOLDALFA(J)
      ENDIF

 10   CONTINUE

      ALFA=DOLDALFA(J)

      IF(DABS(ALFA).LE.1.D-3*DALFAMAX) THEN
	   
	   ALFA=0.D0
         RETURN
      
      END IF

20	CONTINUE

      Z(J) = X(J)+ALFA*D(J)

      ALFAEX=ALFA
       
      CALL FUNCT(Z,N,FZ)
      NUM_FUNCT=NUM_FUNCT+1

      IF(IPRINT.GE.1) THEN
         WRITE(*,*) ' FZ =',FZ,'   ALFA =',ALFA
         WRITE(1,*) ' FZ =',FZ,'   ALFA =',ALFA
      ENDIF
      IF(IPRINT.GE.2) THEN
         DO I=1,N
            WRITE(*,*) ' Z(',I,')=',Z(I)
            WRITE(1,*) ' Z(',I,')=',Z(I)
         ENDDO
      ENDIF

        
      FPAR= F-GAMMA*ALFA*ALFA

      IF(FZ.LE.FPAR) THEN
         FMIN=FZ
         ALFAMIN=ALFA 

C       controllo sull'espansione

   11    CONTINUE

         ALFAEX=ALFA/DELTA1

         Z(J) = X(J)+ALFAEX*D(J)
               
         CALL FUNCT(Z,N,FZDELTA)
         NUM_FUNCT=NUM_FUNCT+1

         IF(IPRINT.GE.1) THEN
            WRITE(*,*) ' FZex=',FZDELTA,'  ALFAEX=',ALFAEX  
            WRITE(1,*) ' FZex=',FZDELTA,'  ALFAEX=',ALFAEX
         ENDIF
         IF(IPRINT.GE.2) THEN
             DO I=1,N
                WRITE(*,*) ' Z(',I,')=',Z(I)
                WRITE(1,*) ' Z(',I,')=',Z(I)
             ENDDO
         ENDIF

         FPAR= F-GAMMA*ALFAEX*ALFAEX
         IF(FZDELTA.LE.FPAR) THEN
                FZ=FZDELTA
                ALFA=ALFAEX
                GO TO 11
         ELSE
                  
              DOLDALFA(J)=DELTA*ALFA

              RETURN
         END IF

      ELSE 

          D(J)=-D(J)
		IF(PRIMA_VOLTA) THEN
			PRIMA_VOLTA = .FALSE.
			GO TO 20
		ENDIF

          DOLDALFA(J)=DELTA*DOLDALFA(J)


          IF(IPRINT.GE.1) THEN
              WRITE(*,*) ' direzione opposta'
              WRITE(1,*) ' direzione opposta'
          ENDIF

          ALFA=0.D0
          RETURN

      END IF

      END

