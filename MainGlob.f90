!============================================================================================
!    DDFSA - A Distributed Derivative-free Simulated Annealing Method for
!    bound constrained global optimization
!    Copyright (C) 2011  G.Liuzzi, S.Lucidi, V.Piccialli
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G. Liuzzi, S. Lucidi, V. Piccialli, A. Sotgiu. A magnetic resonance device designed via 
!    global optimization techniques, Mathematical Programming, 101: 339-364 (2004)
!
!============================================================================================
PROGRAM MAIN_GLOB

	INTEGER, PARAMETER	:: TABELLA = 14
	INTEGER			:: iok
	INTEGER			:: I, MYDATE(8)
	INTEGER			:: NF, NFAL, NFTOT, NFALTOT
	INTEGER			:: ICHECK		

	OPEN (TABELLA,FILE='risultato.tex',STATUS='REPLACE')
	CALL DATE_AND_TIME(VALUES = MYDATE)
	WRITE(TABELLA,1000) MYDATE(3),MYDATE(2),MYDATE(1)

	CALL SENZADERGLOB(TABELLA,iok)
	
	WRITE(TABELLA,1010)
	CLOSE(TABELLA)

1020 FORMAT(I6,'  ',I3)
1010 FORMAT(//,'\end{tabular}				',/,&
'\end{center}				')

1000 FORMAT('\begin{center}\par							',/,&
'{\bf (',I2,'/',I2,'/',I4,')}\par			',/,&
'\begin{tabular}{l|c||c|c|c|c|c|c}			',/,&
'PROBLEM & n & nF & nFott & F. min & F. medio & F. min & F. medio\\ \hline',///)

END PROGRAM

SUBROUTINE SENZADERGLOB(TABELLA,iok)
	!USE IMSLF90

	IMPLICIT NONE

	INTEGER, INTENT(IN)				:: TABELLA 

	INTEGER, PARAMETER				:: IRUNS      = 10

	INCLUDE 'TYPEDECL.FI'

	INTEGER						:: N, NF, NFAILS, IDIMWKS, NUMFAL, IOTT, ISTOP, iok
	INTEGER						:: VECTNF(IRUNS), VECTNFOTT(IRUNS)
	INTEGER						:: LSTOT, NFTOT, NPCTOT, NLMTOT, IT, IMAX
	INTEGER						:: LS, NPC, NLM, ICONTFAL, IDFAL, NUMITER
	INTEGER						:: NUMCAS, K, I, J, MAXLS, MODLS
	INTEGER						:: MAXITER, MAXNF, PRINT_LEVEL, nfott, SUMNF
	INTEGER*4					:: IVAR, IVAR2
	
	TYPE(COMP_WKS_TYP),ALLOCATABLE	:: WKS(:)
	DOUBLE PRECISION,  ALLOCATABLE	:: X(:), XOTT(:), DOLDALFA(:)
	DOUBLE PRECISION,  ALLOCATABLE	:: Z(:), VETT(:), XPROP(:), D(:)
	DOUBLE PRECISION,  ALLOCATABLE	:: XINF(:), XSUP(:)
	DOUBLE PRECISION,  ALLOCATABLE  :: RA(:), RB(:)
	INTEGER			,  ALLOCATABLE	:: IPERM(:)
	DOUBLE PRECISION				:: VECTFVAL(IRUNS), VECTF_OTT(IRUNS)
	DOUBLE PRECISION				:: TOL, FMAX_EFF, fglob
	DOUBLE PRECISION				:: RFAL, RNUMFAL, RMED, RMED2, f_ott, FOTT, F, ALFAMAX, ALFAOTT
	DOUBLE PRECISION				:: FPROP, VAR2, VAR1, ALFA, ZETA, TRES, DIFF
	DOUBLE PRECISION				:: FSTOP, FSTOPTOL, FFSTOP
	DOUBLE PRECISION				:: TEMPERATURA, SUMITER, SUMTIME
	DOUBLE PRECISION				:: F_MIN, F_MIN1, F_MAX, F_MINTOT, F_MINOTT, TCOEFF, FMIN, FMAX

	REAL							:: RRAND
	LOGICAL							:: L_F_DISTANTI, L_F_VICINI, NON_MINIMIZZ, SOLOSULPEGG, NORMADUE
	LOGICAL, ALLOCATABLE			:: MASCH(:)

	CHARACTER ( LEN = 40 )			:: nomefun
	CHARACTER ( LEN = 2  )			:: CSTEP

	COMMON/NUMNEW2/		NF
	COMMON/CALFANR2/    NORMADUE

	NFOTT  = 0

	CALL     NINIT(N)
	ALLOCATE(XINF(N),XSUP(N))
	CALL FUNCTINIT(N,XINF,XSUP,nomefun,fglob)

!-----------------------------------------------
!	SETTA LE TOLLERANZE E I MASSIMI
!-----------------------------------------------
	TOL         = 1.0D-4
	FSTOPTOL    = 1.0D-15
	MAXITER     = 100000000
	MAXNF       = 1000*n
	NUMFAL		= 0
	PRINT_LEVEL = 0
	TCOEFF      = 1.0D0
!--------------------------------------------------------
!	modificato il 23/9/2003
!--------------------------------------------------------
!	IDIMWKS     = 20 !MIN(10*N,10)
!	IDIMWKS     = MAX(10,N)
	IDIMWKS     = MIN(20,MAX(10,N))

!--------------------------------------------------------
	SOLOSULPEGG = .FALSE.
	NORMADUE    = .FALSE.
!-----------------------------------------------
	ALLOCATE(X(N),XOTT(N),Z(N),VETT(N),XPROP(N))
	ALLOCATE(DOLDALFA(N),D(N))
	ALLOCATE(RA(IDIMWKS),RB(IDIMWKS),IPERM(IDIMWKS))
!
!	-------
!	 Iniz.
!	-------
!
	LSTOT  = 0
	NFTOT  = 0
	NFOTT  = 0
	NPCTOT = 0
	NLMTOT = 0

	CALL RANDOM_SEED()
	
	DO 171 IT=1,IRUNS

		ALLOCATE(WKS(IDIMWKS),MASCH(IDIMWKS))
		DO I = 1,IDIMWKS
			ALLOCATE(WKS(I)%X(N))
			ALLOCATE(WKS(I)%D(N))
			ALLOCATE(WKS(I)%DOLDALFA(N))
		ENDDO

		MASCH    = .FALSE.
		LS       = 0
		NF       = 0
		NFOTT    = 0
		NPC      = 0
		NLM      = 0
		RFAL     = 0.0D0
		RNUMFAL  = 0.0D0
		ICONTFAL = 0
		IDFAL    = 0  	
		NUMITER  = 0

		ivar     = 2*it + 3*(it+7)
		ivar2    = 3*it + 2*(it+9) -5

		RMED     = 0.D0
		RMED2    = 0.D0

!	-------------------------------------
!	  Generazione del campione originale
!	-------------------------------------
		FOTT     = 1.D+32
!--------------------------------------------------------
!	modificato il 23/9/2003
!--------------------------------------------------------
		NUMCAS   = 30
!		NUMCAS   = 1
!--------------------------------------------------------
		IF(PRINT_LEVEL >= 0) THEN
			WRITE(*,2070)
			WRITE(*,2080) 
			WRITE(1,2070)
			WRITE(1,2080) 
		ENDIF

		DO K=1,NUMCAS

3000		CONTINUE
			DO I=1,N
				CALL RANDOM_NUMBER(RRAND)
				VETT(I)=DBLE(RRAND) !RAN(IVAR2)
				!VETT(I)=RAN(IVAR2)
				XPROP(I)=XINF(I)+VETT(I)*(XSUP(I)-XINF(I))
			ENDDO

			CALL FUNCT(XPROP,N,FPROP)
			
			NPC=NPC+1
			NF=NF+1
			RMED=((NPC-1.0D0)*RMED+FPROP)/NPC

			IF(NPC.LE.2) THEN
				RMED2=RMED2+FPROP**2
			ELSE
				RMED2=((NPC-2.)*RMED2+FPROP**2)/(NPC-1.0D0)
			ENDIF

			IF(FPROP.LT.FOTT) THEN
				CALL AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
			ENDIF

			IF(PRINT_LEVEL >= 0) THEN
				WRITE(*,2130) NUMITER, NF, FPROP,'--'
				WRITE(1,2130) NUMITER, NF, FPROP,'--'
			ENDIF
		ENDDO

		IF(PRINT_LEVEL >= 2) THEN
			write(*,*)'fine generazione random iniziale nf = ', nf
			write(1,*)'fine generazione random iniziale nf = ', nf
		ENDIF

		VAR2=RMED2-(RMED**2)*(NPC/(NPC-1.))

!---------------------------------------------------------
!	RIEMPIO IL WKS DI IDIMWKS PUNTI
!	CHE SODDISFANO IL TEST DEL S.A.
!---------------------------------------------------------
!	alfa     = 1.0d0
!	DOLDALFA = ALFA
    DO I=1,N
		DOLDALFA(I)=1.0D0*DMAX1(1.D-3,DMIN1(1.D+3,.5D0*DABS(XPROP(I))))
!		DOLDALFA(I)=1.D0
    END DO
	D        = 1.0D0
	xprop    = xott
	fprop    = fott
!	fmax_eff = fott
	if(print_level >= 2) then
		write(*,*) 'ricerca lungo assi'
	endif
	CALL RICERCA_LUNGO_ASSI(N,XPROP,FPROP,ALFA,DOLDALFA,D,FSTOP,XINF,XSUP,Fmax_eff)
	if(print_level >= 2) then
		write(*,*) 'fine ricerca lungo assi'
	endif

	WKS(1)%X       = XPROP 
	WKS(1)%FOB     = FPROP 
	WKS(1)%ALFAMAX = ALFA
	WKS(1)%DOLDALFA= DOLDALFA			
	if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
		write(*,*) '---- 1 ----'
	endif
	WKS(1)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
	WKS(1)%D       = D
	WKS(1)%FSTOP   = FSTOP
	WKS(1)%NUMCOST = 0

	IF(FPROP.LT.FOTT) THEN
		CALL AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
	ENDIF

	IF(PRINT_LEVEL >= 0) THEN
		WRITE(*,2130) NUMITER, NF, FOTT,'IN'
		WRITE(1,2130) NUMITER, NF, FOTT,'IN'
	ENDIF

	DO I = 2,IDIMWKS
		ZETA = 1.0D0
		TRES = 0.0D0

		DO WHILE (ZETA.GT.TRES)
			CALL RANDOM_NUMBER(RRAND)
			ZETA=DBLE(RRAND) !RAN(IVAR)
			!ZETA=RAN(IVAR)

3010		CONTINUE
			DO J = 1,N
				CALL RANDOM_NUMBER(RRAND)
				VETT(J)=DBLE(RRAND) !RAN(IVAR2)
				!VETT(J)=RAN(IVAR2)
				XPROP(J)=XINF(J)+VETT(J)*(XSUP(J)-XINF(J))
			ENDDO

			CALL FUNCT(XPROP,N,FPROP)

			NPC = NPC+1
			NF  = NF+1

			RMED =((NPC-1.)*RMED+FPROP)/NPC
			RMED2=((NPC-2.)*RMED2+FPROP**2)/(NPC-1.)
			VAR2 =RMED2-(RMED**2)*(NPC/(NPC-1.))

			DIFF=RMED-FOTT

			TRES=DEXP(-TCOEFF*(RFAL/DIFF)*(FPROP-FOTT))
!--------------------------------------------------------
!	modificato il 23/9/2003
!--------------------------------------------------------
			IF(FPROP.LE.FOTT) THEN
				TRES = 1.0D0
			ENDIF
!			TRES = 1.0D0
!--------------------------------------------------------
		ENDDO
				
!		DO J=1,N
!			DOLDALFA(J)=1.0D0*DMAX1(1.D-3,DMIN1(1.D+3,.5D0*DABS(XPROP(J))))
!			DOLDALFA(J)=1.0D0
!		END DO

		ALFA      =  WKS(1)%ALFAMAX
		DOLDALFA  =  WKS(1)%DOLDALFA
		DO J = 2,I-1
		   IF(ALFA.LT.WKS(J)%ALFAMAX) THEN
				ALFA      =  WKS(J)%ALFAMAX
				DOLDALFA  =  WKS(J)%DOLDALFA
		   ENDIF
		ENDDO

!		DO K = 1,N
!			ALFA=WKS(1)%DOLDALFA(K)
!			DO J = 2,I-1
!			   IF(ALFA.LT.WKS(J)%DOLDALFA(K)) THEN
!				   ALFA      =  WKS(J)%DOLDALFA(K)
!			   ENDIF
!			ENDDO
!			DOLDALFA(K) = ALFA
!		ENDDO

		D        = 1.0D0
!		fmax_eff = fprop
		CALL RICERCA_LUNGO_ASSI(N,Xprop,Fprop,ALFA,DOLDALFA,D,FSTOP,XINF,XSUP,Fmax_eff)

		WKS(I)%X       = XPROP
		WKS(I)%FOB     = FPROP
		WKS(I)%ALFAMAX = ALFA
		WKS(I)%DOLDALFA= DOLDALFA
		if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
			write(*,*) '---- 2 ----'
		endif
		WKS(I)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
		WKS(I)%D       = D
		WKS(I)%FSTOP   = FSTOP
		WKS(I)%NUMCOST = 0

		IF(FPROP.LT.FOTT) THEN
			CALL AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
		ENDIF
		IF(PRINT_LEVEL >= 0) THEN
			WRITE(*,2130) NUMITER, NF, FOTT,'IN'
			WRITE(1,2130) NUMITER, NF, FOTT,'IN'
		ENDIF

	ENDDO 

!--------------------------------------------------------
!	modificato il 23/9/2003
!--------------------------------------------------------
!	OPEN(99,FILE='MIN_LOCALI.txt',STATUS='REPLACE')
!	DO I = 1,IDIMWKS
!		WRITE(99,*) WKS(I)%X
!		WRITE(99,*) 'FOB = ',WKS(I)%FOB
!		WRITE(99,*) '--------------------------------'
!	ENDDO
!	CLOSE(99)
!	STOP
!--------------------------------------------------------

	DO I=1,N
		X(I)=XOTT(I)
	ENDDO

	MAXLS=100000000      
	MODLS=1

!	--------------------
!	  Ciclo principale
!	--------------------
	IF(PRINT_LEVEL >= 0) THEN
		WRITE(*,2070)
		WRITE(*,2080) 
		WRITE(1,2070)
		WRITE(1,2080) 
	ENDIF
2070 FORMAT( '    ITER      NFTOT       FOTT      SALTO       ALFA       FMIN        FMAX')
2080 FORMAT( '-----------------------------------------------------------------------------')
!              123456  123456789012  +1.2345E-02    OK     +1.2345E-02 +1.2345E-02 +1.2345E-02 
2090 FORMAT(1X,  I6, 2X,  I12,  2X,  ES11.4, 4X, A2, 5X,  ES11.4, 1X,   ES11.4, 1X, ES11.4 )
2130 FORMAT(1X,  I6, 2X,  I12,  2X,  ES11.4, 4X, A2 )

100	FORMAT(/,1X,' GENERAZ.  N = ', I9,/)
101	FORMAT(/,1X,'  Attuale valore ottimo = ',D13.6,/)
102	FORMAT(3(:,1X,'XOTT(',I1,') =',D13.6,1X))
109	FORMAT(3(:,1X,'XOTT(',I2,') =',D13.6,1X))
6	FORMAT(/,1X,' NF=',I15,4X,' NPC=',I7,4X,' NLM=',I7,4X,' LS=',I7,/)
96	FORMAT(/,1X,'   MED = ',D13.6,'    VAR = ',D13.6,/)

	CSTEP ='  '

7	CONTINUE
8	CONTINUE

!--------------------------------------
!	stopping criterion
!--------------------------------------
!	calcola il max tra gli alfamax
!   o alfanr2 a seconda che normadue 
!	sia falso o veroo
!--------------------------------------
	ALFA     = 0.0D0
	FSTOP    = 0.0D0
	FMIN     = 1.0D+30
	FMAX     =-1.0D+30
	DO I = 1,IDIMWKS
		IF(.NOT.MASCH(I)) THEN
			IF(WKS(I)%FSTOP > FSTOP) FSTOP = WKS(I)%FSTOP
			IF(NORMADUE) THEN
				IF(ALFA.LT.WKS(I)%ALFANR2) THEN
					ALFA     = WKS(I)%ALFANR2
				ENDIF
			ELSE
				IF(ALFA.LT.WKS(I)%ALFAMAX) THEN
					ALFA     = WKS(I)%ALFAMAX
				ENDIF
			ENDIF
!			IF(WKS(I)%FOB > FMAX) FMAX = WKS(I)%FOB
		ENDIF
		IF(WKS(I)%FOB > FMAX) FMAX = WKS(I)%FOB
		IF(WKS(I)%FOB < FMIN) THEN
			FMIN = WKS(I)%FOB
			IF(NORMADUE) THEN
				ALFAOTT = WKS(I)%ALFANR2
			ELSE
				ALFAOTT = WKS(I)%ALFAMAX
			ENDIF
		ENDIF
	ENDDO

	ALFAMAX = ALFA

	IF(PRINT_LEVEL >= 1) THEN
		WRITE(*,2070)
		WRITE(*,2080) 
		WRITE(1,2070)
		WRITE(1,2080) 
	ENDIF
	IF(PRINT_LEVEL >= 0) THEN
		WRITE(*,2090) NUMITER, NF, FOTT,CSTEP,ALFA,FMIN,FMAX
		WRITE(1,2090) NUMITER, NF, FOTT,CSTEP,ALFA,FMIN,FMAX
	ENDIF

	IF((ALFAMAX <= TOL).OR.(FSTOP <= FSTOPTOL).OR. &
	   ((FOTT-FGLOB)/max(1.0d0,abs(fglob)) < tol)) THEN
		IF(NFOTT == 0) THEN
			NFOTT = NF
			F_OTT = FOTT
		ENDIF
	ENDIF

	!IF(NF >= MAXNF) THEN
	IF((ALFAMAX <= TOL).OR.(FSTOP <= FSTOPTOL).OR.(NUMITER >= MAXITER).OR.(NF >= MAXNF)) THEN
		IF(PRINT_LEVEL >= 0) THEN	 
			write(*,*) 'alfamax = ',alfamax	
			write(1,*) 'alfamax = ',alfamax	
			write(*,*) '  fstop = ',fstop	
			write(1,*) '  fstop = ',fstop	
			WRITE(1,100) NUMITER
			WRITE(*,100) NUMITER
			WRITE(1,*) FOTT
			WRITE(*,*) FOTT
			WRITE(1,101) FOTT
			WRITE(*,101) FOTT

			IF(N.LT.10) WRITE(1,102) (J,XOTT(J),J=1,N)
			IF(N.LT.10) WRITE(*,102) (J,XOTT(J),J=1,N)

			IF((N.GT.9).AND.(N.LT.100)) WRITE(1,109) (J,XOTT(J),J=1,N)
			IF((N.GT.9).AND.(N.LT.100)) WRITE(*,109) (J,XOTT(J),J=1,N)

			WRITE(1,6) NF,NPC,NLM,LS
			WRITE(6,6) NF,NPC,NLM,LS

			if(var2 < 0.0d0) then
				write(*,*) '---- 3 ----'
			endif
			VAR1=DSQRT(VAR2)
			WRITE(6,96) RMED,VAR1
			WRITE(1,96) RMED,VAR1
		ENDIF
		GO TO 170
	ENDIF

	!IF (FMIN /= FOTT) PAUSE

	IF(PRINT_LEVEL >= 1) THEN


		WRITE(*,*)'/---------------------------------\'
		WRITE(1,*)'/---------------------------------\'
		WRITE(*,*)'ALFAOTT = ',ALFAOTT,' FOTT = ',FOTT
		WRITE(1,*)'ALFAOTT = ',ALFAOTT,' FOTT = ',FOTT
!		WRITE(*,*) XOTT
!		WRITE(1,*) XOTT
		WRITE(*,*)'\---------------------------------/'
		WRITE(1,*)'\---------------------------------/'
		WRITE(*,*)
		WRITE(1,*)
	ENDIF


!--------------------------------------

	NUMITER = NUMITER + 1
	IF(PRINT_LEVEL >= 0) THEN
		IF(MOD(NUMITER,30)==0) THEN
			WRITE(*,2070)
			WRITE(*,2080) 
			WRITE(1,2070)
			WRITE(1,2080) 
		ENDIF
	ENDIF

!--------------------------------------
!	genera zeta random tra 0 e 1
!	e un punto xprop su cui calcolare
!	la f.ob.
!--------------------------------------
	CALL RANDOM_NUMBER(RRAND)
	ZETA=DBLE(RRAND) !RAN(IVAR)
	!ZETA=RAN(IVAR)


3020 CONTINUE

	DO I=1,N
		CALL RANDOM_NUMBER(RRAND)
		VETT(I)=DBLE(RRAND) !RAN(IVAR2)
		!VETT(I)=RAN(IVAR2)
		XPROP(I)=XINF(I)+VETT(I)*(XSUP(I)-XINF(I))
	ENDDO	

	CALL FUNCT(XPROP,N,FPROP)

	NPC=NPC+1
	NF=NF+1

103	FORMAT(/,1X,'      Proposta di salto')
105	FORMAT(/,1X,'   Zeta =',D13.6,'   Tres =',D13.6,/)

!	----------------------------
!	  Calcolo media, varianza e
!	  la temperatura del SA
!	----------------------------

	RMED  = ((NPC-1.)*RMED+FPROP)/NPC
	RMED2 = ((NPC-2.)*RMED2+FPROP**2)/(NPC-1.)
	VAR2  = RMED2-(RMED**2)*(NPC/(NPC-1.))

	DIFF  = RMED-FOTT

	TRES  = DEXP(-TCOEFF*(RFAL/DIFF)*(FPROP-FOTT))
	temperatura = diff/rfal

	IF(FPROP.LE.FOTT) TRES=1.0D0

!--------------------------------------
!	se zeta e' minore o uguale a tres
!	allora accetta il punto e fa
!	partire una min. locale
!--------------------------------------
!	IF(.FALSE.) THEN
	IF(ZETA.LE.TRES) THEN

		NUMFAL = 0
		IF(PRINT_LEVEL >= 2) THEN
			WRITE(1,405) 
			WRITE(*,405) 
		ENDIF
405		FORMAT(1X,'    Proposta di salto accettata ',/)
		DO I=1,N
			X(I)=XPROP(I)
		ENDDO

		NLM     = NLM+1

!--------------------------------------
!		calcola l'indice imax della
!		componente di wks contenente
!		il valore di fob piu' alto
!--------------------------------------
		FMAX    =-1.0D+30
		DO I = 1,IDIMWKS
!			IF((WKS(I)%FOB > FMAX)) THEN
			IF((WKS(I)%FOB > FMAX) .AND. (.NOT.MASCH(I))) THEN
				FMAX = WKS(I)%FOB
				IMAX = I
			ENDIF
		ENDDO

		ALFA  = WKS(1)%ALFAMAX
		FSTOP = WKS(1)%FSTOP
		DO J = 2,IDIMWKS
			IF(FSTOP.LT.WKS(J)%FSTOP) FSTOP = WKS(J)%FSTOP
			IF(ALFA.LT.WKS(J)%ALFAMAX) THEN
				ALFA      =  WKS(J)%ALFAMAX
				DOLDALFA  =  WKS(J)%DOLDALFA
			ENDIF
		ENDDO

!		DO K = 1,N
!			ALFA=WKS(1)%DOLDALFA(K)
!			DO J = 1,IDIMWKS
!			   IF(ALFA.LT.WKS(J)%DOLDALFA(K)) THEN
!				   ALFA      =  WKS(J)%DOLDALFA(K)
!			   ENDIF
!			ENDDO
!			DOLDALFA(K) = ALFA
!		ENDDO

		ALFA = ALFAMAX
!		IF(NORMADUE) THEN
!			ALFA = WKS(IMAX)%ALFANR2
!		ELSETABELLA
!			ALFA = WKS(IMAX)%ALFAMAX
!		ENDIF

		CALL INTERFACCIA(N,X,F,ALFA,DOLDALFA,D,FFSTOP,FSTOP,XINF,XSUP,FOTT,ALFAMAX,ALFAOTT,XOTT,WKS,IDIMWKS,ISTOP,MASCH,MAXNF,fglob)

		IF(F.LT.FOTT) THEN
			CALL AGGIORNO_OTTIMO(N,X,F,XOTT,FOTT)
			ISTOP = 1
		ENDIF

		IF((F.LT.WKS(IMAX)%FOB).AND.(ISTOP.NE.2)) THEN
!--------------------------------------
!			se la min. produce un punto
!			migliore del peggiore
!			lo sostituisce
!--------------------------------------
			WKS(IMAX)%X       = X
			WKS(IMAX)%FOB     = F
			WKS(IMAX)%ALFAMAX = ALFA
			WKS(IMAX)%DOLDALFA= DOLDALFA
			WKS(IMAX)%D		  = D
			WKS(IMAX)%FSTOP	  = FFSTOP
			if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
				write(*,*) '---- 4 ----'
			endif
			WKS(IMAX)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
			WKS(IMAX)%NUMCOST = 0
			CSTEP = 'OK'
			MASCH(IMAX) = .FALSE.
			GO TO 7
		ELSE
!--------------------------------------
!			non accetta il punto
!			prodotto dalla minimizzazione
!--------------------------------------
			RNUMFAL  = RNUMFAL  + 1.0D0
			ICONTFAL = ICONTFAL + 1
			RFAL     = RFAL     + 1.0D0/RNUMFAL
			IDFAL    = 1

			CSTEP = 'NO'
	
		ENDIF

	ELSE
!--------------------------------------
!	altrimenti (zeta > tres) scarta
!	il punto proposto e aggiorna
!	il cont. numfal
!--------------------------------------
		NUMFAL   = NUMFAL   + 1

!		RNUMFAL  = RNUMFAL  + 1.0D0
!		ICONTFAL = ICONTFAL + 1
!		RFAL     = RFAL     + 1.0D0/RNUMFAL
!		IDFAL    = 1
	ENDIF

!====================================================
!   QUI DEVO FARE UNA RICERCA LUNGO GLI ASSI
!	SU I PUNTI DEL WKS. 
!   SOLOSULPEGG = .TRUE.  --> LA FA SOLO SU QUELLO
!						 	  CON ALFA PEGGIORE
!				= .FALSE. --> LA FA SU TUTTI
!====================================================	

	if(.NOT.SOLOSULPEGG) then !fai la ricerca lungo gli assi per tutti
							  !i punti dell'array
		IF(NUMFAL >= 0) THEN
			NUMFAL  = 0

			FMAX_EFF    =-1.0D+30
			DO J = 1,IDIMWKS
!				IF(.NOT.MASCH(J)) THEN
					IF(WKS(J)%FOB > FMAX_EFF) THEN
						FMAX_EFF = WKS(J)%FOB
					ENDIF
!				ENDIF
				IF(WKS(J)%FOB <= FOTT) IOTT = J
				RA(J)    = WKS(J)%FOB
				IPERM(J) = J
			ENDDO

!-------------------------------------------------
! ORDINO RA E RESTITUISCE RB ORDINATO
! IN IPERM CI SONO LE POSIZIONI DEGLI ELEMENTI
! DI RB IN RA
!-------------------------------------------------
			call qsortd(RA,IPERM,IDIMWKS)
			DO I = 1,IDIMWKS
				RB(I) = RA(IPERM(I))
			ENDDO
			!CALL DSVRGP(IDIMWKS,RA,RB,IPERM)
!			MASCH = .FALSE.

			I = 1
			J = 2
			DO WHILE ( J <= IDIMWKS )

				!CONFRONTO J,I --> MASCH
				IF (.NOT.MASCH(IPERM(J))) THEN


					L_F_VICINI = (DABS(RB(I)-RB(J))/(DMAX1(1.D-8,DABS(RB(J)) )) ).LE.1.D-1 
!					L_F_VICINI = L_F_VICINI.AND.((DSQRT(													&
!						DOT_PRODUCT(WKS(IPERM(I))%X-WKS(IPERM(J))%X,WKS(IPERM(I))%X-WKS(IPERM(J))%X)) <= 1.D-2)& !ERA 1.D-1
!						.OR. (WKS(IPERM(J))%ALFAMAX.LE.1.D-1*WKS(IPERM(I))%ALFAMAX) )
					L_F_VICINI = L_F_VICINI.AND.((DSQRT(													&
						DOT_PRODUCT(WKS(IPERM(I))%X-WKS(IPERM(J))%X,WKS(IPERM(I))%X-WKS(IPERM(J))%X)) <= 1.D-1))
	!				MASCH(IPERM(J)) = MASCH(IPERM(J)).OR.L_F_VICINI

					IF (L_F_VICINI) THEN
						
						IF ( (WKS(IPERM(J))%ALFAMAX <= WKS(IPERM(I))%ALFAMAX) .OR. (I == 1) ) THEN
							
							MASCH(IPERM(J))=.TRUE.
							J = J+1 

						ELSE
							
							MASCH(IPERM(I)) =.TRUE.
							I = J
							J = J+1 								
								
						ENDIF

					ELSE

						I = J   
						J = J+1

					ENDIF

				ELSE
					
					J = J+1

				ENDIF

			ENDDO

			

			DO I = 1,IDIMWKS

				IF(NORMADUE) THEN
					ALFA     = WKS(I)%ALFANR2
				ELSE
					ALFA     = WKS(I)%ALFAMAX
				ENDIF
				DOLDALFA = WKS(I)%DOLDALFA
				D        = WKS(I)%D
				X        = WKS(I)%X
				F        = WKS(I)%FOB


				IF(I.NE.IOTT) THEN

					RA(I) = (F-FOTT)/DMAX1(1.D-3,DMIN1(DABS(F),DABS(FOTT)) )
					RB(I) = DSQRT(DOT_PRODUCT(X-XOTT,X-XOTT))
					L_F_DISTANTI = ((F-FOTT)/DMAX1(1.D-3,DMIN1(DABS(F),DABS(FOTT)) ) ).GE.1.D+1
!					L_F_DISTANTI = L_F_DISTANTI.AND.(ALFA <= 1.D-4)
					L_F_DISTANTI = L_F_DISTANTI.AND.(ALFA <= 1.D+1*ALFAOTT)
					L_F_VICINI = ((F-FOTT)/DMAX1(1.D-3,DABS(FOTT) ) ).LE.1.D-1
!					L_F_VICINI = L_F_VICINI.AND.( (DSQRT(DOT_PRODUCT(X-XOTT,X-XOTT)) <= 1.D-1) &
!								.OR. (ALFA.LE.1.D-1*ALFAOTT))
					L_F_VICINI = L_F_VICINI.AND.( (DSQRT(DOT_PRODUCT(X-XOTT,X-XOTT)) <= 1.D-1))
					NON_MINIMIZZ = L_F_DISTANTI .OR. L_F_VICINI
!					NON_MINIMIZZ = L_F_DISTANTI.AND.(ALFA.LE.1.D-4)
				ELSE
					NON_MINIMIZZ = .FALSE.
					RA(I) = -1.D0
					RB(I) = -1.D0
				ENDIF
				
				MASCH(I) = MASCH(I) .OR. NON_MINIMIZZ
				  
!			ENDDO



!--------------------------------------------
!	apro il file per fare i controlli
!	sull'alfamax
!--------------------------------------------
			IF (PRINT_LEVEL >= 1) THEN

				open(99,FILE='alfamax.txt',STATUS='REPLACE')
				DO J = 1,IDIMWKS
					write(99,1999) J,wks(J)%alfamax,wks(J)%fob,masch(J), RA(J), RB(J)
					write(*,1999) J,wks(J)%alfamax,wks(J)%fob,masch(J), RA(J), RB(J)
				enddo

			1999 format(1x,i3,1x,es18.10,1x,es18.10,1x,l1,1x, es9.2, 1x, es9.2)
			!--------------------------------------------
			!	apro il file per fare i controlli
			!	sull'alfamax
			!--------------------------------------------
				close(99)
				!pause

			ENDIF

!			DO I = 1,IDIMWKS

				IF(NORMADUE) THEN
					ALFA     = WKS(I)%ALFANR2
				ELSE
					ALFA     = WKS(I)%ALFAMAX
				ENDIF
				DOLDALFA = WKS(I)%DOLDALFA
				D        = WKS(I)%D
				X        = WKS(I)%X
				F        = WKS(I)%FOB


!				IF(.TRUE.) THEN
!				IF(	  (ALFA.GT.TOL).AND..NOT.MASCH(I)	)	THEN
!				IF(	  (ALFA.GT.TOL).AND..NOT.MASCH(I).AND. (ALFA.GT.1.D-3*ALFAMAX)	)	THEN
				IF(	 ((ALFA.GT.TOL).AND..NOT.MASCH(I).AND. (ALFA.GT.1.D-1*ALFAMAX)).OR. &
				   	 ((ALFA.GT.TOL).AND..NOT.MASCH(I).AND. (I == IOTT)))	THEN

					MASCH(I) = .FALSE.
!					FMAX_EFF = F
!					WRITE(*,*)'X=',X
!					WRITE(*,*)'F=',F
					CALL RICERCA_LUNGO_ASSI(N,X,F,ALFA,DOLDALFA,D,FSTOP,XINF,XSUP,FMAX_EFF)
					if(f > fmax_eff) then
						write(*,*) '---> 1 <--- diff = ',f-fmax_eff
!						pause
					endif
!					WRITE(*,*)'X=',X
!					WRITE(*,*)'F=',F
					WKS(I)%X       = X
					WKS(I)%FOB     = F
!					PAUSE
					!IF(WKS(I)%ALFAMAX > ALFA) THEN
					!	WRITE(*,*) '\/\/\/\/\/ ',F,' \/\/\/\/\/\ ',FMAX_EFF !/\/\/\/\/', I
					!ENDIF
					IF(ALFA >= WKS(I)%ALFAMAX) THEN
						WKS(I)%NUMCOST = WKS(I)%NUMCOST + 1
					ELSE
						WKS(I)%NUMCOST = 0
					ENDIF

					WKS(I)%ALFAMAX = ALFA
					WKS(I)%DOLDALFA= DOLDALFA
					WKS(I)%D       = D
					WKS(I)%FSTOP   = FSTOP
					if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
						write(*,*) '---- 5 ----'
					endif
					WKS(I)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
					IF(F.LT.FOTT) THEN	
						CALL AGGIORNO_OTTIMO(N,X,F,XOTT,FOTT)
					ENDIF
					IF(WKS(I)%NUMCOST <= -1) THEN
						DO WHILE (ALFA >= WKS(I)%ALFAMAX)
							CALL RICERCA_LUNGO_ASSI(N,X,F,ALFA,DOLDALFA,D,FSTOP,XINF,XSUP,FMAX_EFF)
							if(f > fmax_eff) then
								write(*,*) '---> 2 <---'
!								pause
							endif
							WKS(I)%X       = X
							WKS(I)%FOB     = F
							WKS(I)%ALFAMAX = ALFA
							WKS(I)%DOLDALFA= DOLDALFA
							WKS(I)%D       = D
							WKS(I)%FSTOP   = FSTOP
							if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
								write(*,*) '---- 6 ----'
							endif
							WKS(I)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
							WRITE(*,*) ALFA,F
						ENDDO
						WKS(I)%NUMCOST = 0
					ENDIF
				ELSE
					IF(.FALSE.) THEN
!					IF(NON_MINIMIZZ) THEN
						IF(PRINT_LEVEL >= 10) THEN
							WRITE(1,*) 'SOSTITUISCO IL PUNTO ',I
						ENDIF
						CALL GENERA_PUNTO_SA(N,XPROP,FPROP,XINF,XSUP,FOTT,TCOEFF,IVAR,IVAR2, &
								     VAR2,RMED,RMED2,NLM,NUMFAL,RFAL,RNUMFAL,ICONTFAL,IDFAL,NPC,NF)
						ALFA     = 0.0D0
						FSTOP    = 0.0D0
						DO J = 1,IDIMWKS
							IF(.NOT.MASCH(J)) THEN
								IF(WKS(J)%FSTOP > FSTOP) FSTOP = WKS(J)%FSTOP								
								IF(NORMADUE) THEN
									IF(ALFA.LT.WKS(J)%ALFANR2) THEN
										ALFA     = WKS(J)%ALFANR2
									ENDIF
								ELSE
									IF(ALFA.LT.WKS(J)%ALFAMAX) THEN
										ALFA     = WKS(J)%ALFAMAX
									ENDIF
								ENDIF
							ENDIF
						ENDDO
!						ALFAMAX = ALFA

						ALFA = ALFAMAX
						CALL INTERFACCIA(N,XPROP,FPROP,ALFA,DOLDALFA,D,FFSTOP,FSTOP,XINF,	&
							XSUP,FOTT,ALFAMAX,ALFAOTT,XOTT,WKS,IDIMWKS,ISTOP,MASCH,MAXNF,fglob)
						WKS(I)%X       = XPROP
						WKS(I)%FOB     = FPROP
						WKS(I)%ALFAMAX = ALFA
						WKS(I)%DOLDALFA= DOLDALFA
						WKS(I)%D       = D
						WKS(I)%FSTOP   = FFSTOP
						if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
							write(*,*) '---- 7 ----'
						endif
						WKS(I)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
						IF(FPROP.LT.FOTT) THEN
							CALL AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
						ENDIF
						MASCH(I) = .FALSE.
					ELSE
						MASCH(I) = .TRUE.
					ENDIF
				ENDIF

			ENDDO
		ENDIF
	else !fai la ricerca lungo gli assi solo per
		 !il punto con l'alfa maggiore

		ALFA     = 0.0D0
		FMIN     = 1.0D+30
		FMAX     =-1.0D+30
		DO I = 1,IDIMWKS
			IF(.NOT.MASCH(I)) THEN
				IF(NORMADUE) THEN
					IF(ALFA.LT.WKS(I)%ALFANR2) THEN
						ALFA     = WKS(I)%ALFANR2
						DOLDALFA = WKS(I)%DOLDALFA
						IMAX     = I
					ENDIF
				ELSE
					IF(ALFA.LT.WKS(I)%ALFAMAX) THEN
						ALFA     = WKS(I)%ALFAMAX
						DOLDALFA = WKS(I)%DOLDALFA
						IMAX     = I
					ENDIF
				ENDIF
				IF(WKS(I)%FOB < FMIN) FMIN = WKS(I)%FOB
				IF(WKS(I)%FOB > FMAX) FMAX = WKS(I)%FOB
			ENDIF
		ENDDO

		D    = WKS(IMAX)%D
		X    = WKS(IMAX)%X
		F    = WKS(IMAX)%FOB
		CALL RICERCA_LUNGO_ASSI(N,X,F,ALFA,DOLDALFA,D,FSTOP,XINF,XSUP,FMAX_EFF)
		if(f > fmax_eff) then
			write(*,*) '---> 3 <---'
!			pause
		endif
		WKS(IMAX)%X       = X
		WKS(IMAX)%FOB     = F
		WKS(IMAX)%ALFAMAX = ALFA
		WKS(IMAX)%DOLDALFA= DOLDALFA
		WKS(IMAX)%D       = D
		WKS(IMAX)%FSTOP   = FSTOP
		if(DOT_PRODUCT(DOLDALFA,DOLDALFA) < 0.0d0) then
			write(*,*) '---- 8 ----'
		endif
		WKS(IMAX)%ALFANR2 = DSQRT(DOT_PRODUCT(DOLDALFA,DOLDALFA))
		IF(F.LT.FOTT) THEN
			CALL AGGIORNO_OTTIMO(N,X,F,XOTT,FOTT)
		ENDIF

	endif

	IF(CSTEP.NE.'NO') CSTEP = 'KO'
		
	GO TO 8

170 CONTINUE

	IF(NFOTT == 0) THEN
		NFOTT = NF
		F_OTT = FOTT
	ENDIF

	LSTOT  = LSTOT  + LS
	NFTOT  = NFTOT  + NF
	NPCTOT = NPCTOT + NPC
	NLMTOT = NLMTOT + NLM

	VECTNF(IT)   = NF
	VECTFVAL(IT) = FOTT
	VECTNFOTT(IT)= NFOTT
	VECTF_OTT(IT)= F_OTT

	DO I = 1,IDIMWKS
		DEALLOCATE(WKS(I)%X)
		DEALLOCATE(WKS(I)%DOLDALFA)
		DEALLOCATE(WKS(I)%D)
	ENDDO
	DEALLOCATE(WKS,MASCH)

171 CONTINUE

	sumiter  = 0.D0 ; sumnf = 0 ; sumtime = 0.0D0 
	f_min1   = 1.D+30
	f_minott = 0.D0
	f_min    = 1.D+30
	f_mintot = 0.D0
	f_max    = -1.D+30

	DO IT=1,Iruns
		SUMNF = SUMNF + VECTNFOTT(IT)
		f_mintot = f_mintot + vectfVAL(IT)
		f_minott = f_minott + vectf_ott(IT)
		IF (f_min.GT.vectfVAL(IT)) f_min = vectfVAL(IT)
		IF (f_min1.GT.vectf_ott(IT)) f_min1 = vectf_ott(IT)
		IF (f_max.LT.vectfVAL(IT)) f_max = vectfVAL(IT)
	END DO

	write(TABELLA,2030) nomefun, n, NFTOT/IRUNS, SUMNF/IRUNS
	
	F_MINTOT = F_MINTOT/IRUNS
	F_MINOTT = F_MINOTT/IRUNS

	IF (DABS(F_MIN).LT.1.D-3) THEN
		WRITE(TABELLA,2110) F_MIN
	ELSE
		WRITE(TABELLA,2100) F_MIN
	ENDIF
	IF (DABS(F_MINTOT).LT.1.D-3) THEN
		WRITE(TABELLA,2110) F_MINTOT
	ELSE
		WRITE(TABELLA,2100) F_MINTOT
	ENDIF
	IF (DABS(F_MIN1).LT.1.D-3) THEN
		WRITE(TABELLA,2110) F_MIN1
	ELSE
		WRITE(TABELLA,2100) F_MIN1
	ENDIF
	IF (DABS(F_MINOTT).LT.1.D-3) THEN
		WRITE(TABELLA,2110) F_MINOTT
	ELSE
		WRITE(TABELLA,2100) F_MINOTT
	ENDIF
	WRITE(TABELLA,2110) (F_MINOTT-fglob)/max(1.0d0,abs(fglob))
	
	if((F_MINOTT-fglob)/max(1.0d0,abs(fglob))<1.d-4) then
		iok = 1
	else
		iok = 0
	endif

	NFAILS = 0
	DO IT=1,IRUNS
		IF(DABS(VECTFVAL(IT)-F_MIN).GT.1.0d-2) NFAILS = NFAILS+1
	ENDDO
	!WRITE(TABELLA,2120) NFAILS
	WRITE(TABELLA,2180)

	DEALLOCATE(XINF,XSUP)
	DEALLOCATE(X,XOTT,Z,VETT,XPROP,DOLDALFA,D)
	DEALLOCATE(RA,RB,IPERM)

2100 FORMAT(' & ',ES12.5,$)
2110 FORMAT(' & ',D9.3,$)

2030 FORMAT(A40,' & ',I2,' & ',I6,' & ',I6, $)
!2100 FORMAT(' & ',F13.4,$)
!2110 FORMAT(' & ',D9.3,$)
2120 FORMAT(' & ',I4,'\\')
2180 FORMAT(' \\ ')
2140 FORMAT(I3,'  ',$)
2150 FORMAT(ES11.4,'  ',$)
2160 FORMAT(I6,'  ',$)
2170 FORMAT(I6,'  ',I3)

691 FORMAT(/,1X,' NFTOT=',I7,4X,' NPCTOT=',I7,4X,' NLMTOT=',I7,4X,' LSTOT=',I7,/)

END SUBROUTINE SENZADERGLOB

!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE AGGIORNO_OTTIMO(N,XPROP,FPROP,XOTT,FOTT)
	IMPLICIT NONE

	INTEGER,	  INTENT(IN)		:: N
	DOUBLE PRECISION, INTENT(IN)		:: XPROP(N), FPROP
	DOUBLE PRECISION, INTENT(OUT)		:: XOTT(N),  FOTT
	INTEGER					:: I

	XOTT = XPROP
	FOTT = FPROP

	RETURN

	OPEN(99,FILE='OttimoCorrente.txt',STATUS='REPLACE')
	DO I = 1,N
		WRITE(99,1000) I, XOTT(I)
	ENDDO
	WRITE(99,1010) FOTT
	CLOSE(99)

1000 FORMAT(1X,'XOTT(',I3,') = ',ES11.4)
1010 FORMAT(1X,'FOTT      = ',ES11.4)

	RETURN

END SUBROUTINE AGGIORNO_OTTIMO

!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE RICERCA_LUNGO_ASSI(N,X,FF,ALFA,DOLDALFA,D,FFSTOP,XINF,XSUP,FMAX_EFF)
	IMPLICIT NONE
	INTEGER N, NN , nf  

	INTEGER I,J,IDIR_CORR,NUM_FUNCT,NUM_ITER,NUM_NOTFAL(N)
	INTEGER MAX_NOTFAL,NUM_FAL,ISTOP
	INTEGER IPUNT,NUM_PARAM
	INTEGER IPRINT
	INTEGER	ENV, LINK

	REAL*8 FMAX_EFF
	REAL*8 X(N),DIREZIONI(N),Z(N),Z1(N),D(N),XOLD(N)
	REAL*8 GAMMA,RHO,ALFA,F,FF,FZ,FTAR,F0,DSTOP,FFSTOP   
	REAL*8 DOLDALFA(N),DCONV(N),FINIT(N,2) 
	REAL*8 FSTOP(N+1),XFSTOP(N,N+1),XINF(N),XSUP(N),ALFA_MIN(N)

	COMMON /NUM/F
	COMMON /NUMNEW2/NUM_FUNCT
	COMMON /NUMNEW/NF


	F = FF


	IPRINT=-1

	call funct(x,n,ff)
!	write(*,*) '/--------------------------------------\'
!	write(*,*) 'ric.assi: ff = ',ff,' f = ',f

!	CALL MAINBOX1(N,X,D,DIREZIONI,Z,Z1,XOLD,NUM_ITER,DOLDALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,XINF,XSUP,ALFA_MIN,FFSTOP,FMAX_EFF)
	CALL MAINBOX1(N,X,D,DIREZIONI,Z,Z1,XOLD,NUM_ITER,DOLDALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,FFSTOP,XINF,XSUP)

	ALFA = MAXVAL(DOLDALFA)

	num_funct = num_funct + nf

	call funct(x,n,ff)

	if(f.ne.ff) then
		write(*,*) 'ric.assi: ff = ',ff,' f = ',f
		write(*,*) x
	    write(*,*) '\--------------------------------------/'
!		pause
	endif

	FF = F
	

	RETURN

END SUBROUTINE RICERCA_LUNGO_ASSI

!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE INTERFACCIA(N,X,FF,ALFA,DOLDALFA,D,FFSTOP,FFSTOPTOL,XINF,XSUP,FOTT,ALFAMAX,ALFAOTT,	&
		       XOTT,WKS,IDIMWKS,ISTOP,MASCH,MAXNF,fglob)

	IMPLICIT NONE

	INCLUDE 'TYPEDECL.FI'

	INTEGER N, NN , nf  

	INTEGER I,J,IDIR_CORR,NUM_FUNCT,NUM_ITER,NUM_NOTFAL(N),IDIMWKS,MAXNF
	INTEGER MAX_NOTFAL,NUM_FAL,ISTOP
	INTEGER IPUNT,NUM_PARAM
	INTEGER IPRINT
	INTEGER	ENV, LINK

	REAL*8 X(N),DIREZIONI(N),Z(N),Z1(N),D(N),XOLD(N),XOTT(N),fglob
	REAL*8 GAMMA,RHO,ALFA,F,FF,FZ,FTAR,F0,DSTOP,ALFAOTT    
	REAL*8 FFSTOP,FFSTOPTOL,FOTT,ALFAMAX
	REAL*8 DOLDALFA(N),DCONV(N),FINIT(N,2) 
	REAL*8 FSTOP(N+1),XFSTOP(N,N+1),XINF(N),XSUP(N),ALFA_MIN(N)
	TYPE(COMP_WKS_TYP)	:: WKS(IDIMWKS)
	LOGICAL MASCH(IDIMWKS)


	COMMON /NUM/F
	COMMON /NUMNEW2/NUM_FUNCT
	COMMON /NUMNEW/NF


	F = FF

	IPRINT=-1

!	CALL MAINBOX2(N,X,D,DIREZIONI,Z,Z1,XOLD,NUM_ITER,DOLDALFA,ALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,XINF,XSUP,ALFA_MIN,FFSTOP,FFSTOPTOL,FOTT,ALFAMAX)
	CALL MAINBOX2(N,X,D,DIREZIONI,Z,Z1,XOLD,NUM_ITER,DOLDALFA,ALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,FFSTOP,FFSTOPTOL,&
					FOTT,ALFAOTT,XOTT,WKS,IDIMWKS,ISTOP,MASCH,MAXNF,fglob,XINF,XSUP)
	ALFA = MAXVAL(DOLDALFA)

	num_funct = num_funct + nf

	call funct(x,n,ff)
!	write(*,*) 'ric.alfa: ff = ',ff,' f = ',f
	
	FF = F

	IPRINT=-1
	IF(IPRINT>=1) THEN
		WRITE(*,*) X
	ENDIF
		

	RETURN

END SUBROUTINE INTERFACCIA

!==============================================================================================
!==============================================================================================
!==============================================================================================
SUBROUTINE GENERA_PUNTO_SA(N,XPROP,FPROP,XINF,XSUP,FOTT,TCOEFF,IVAR,IVAR2, &
			   VAR2,RMED,RMED2,NLM,NUMFAL,RFAL,RNUMFAL,ICONTFAL,IDFAL,NPC,NF)
	IMPLICIT NONE
	
	INTEGER					:: N, NLM, NUMFAL, ICONTFAL, IDFAL
	INTEGER					:: NPC, NF
	INTEGER*4				:: IVAR, IVAR2
	DOUBLE PRECISION			:: XPROP(N), XINF(N), XSUP(N), FPROP, FOTT, TCOEFF
	DOUBLE PRECISION			:: VAR2, RMED, RMED2, RFAL, RNUMFAL
	DOUBLE PRECISION			:: ZETA, TRES, DIFF
	DOUBLE PRECISION			:: VETT(N)
	LOGICAL					:: TROVATO
	INTEGER					:: I
	REAL					:: RRAND

	TROVATO = .FALSE.

	DO WHILE (.NOT.TROVATO)
		CALL RANDOM_NUMBER(RRAND)
		ZETA=DBLE(RRAND) !RAN(IVAR)
		!ZETA=RAN(IVAR)

100		CONTINUE
		DO I=1,N
			CALL RANDOM_NUMBER(RRAND)
			VETT(I)=DBLE(RRAND) !RAN(IVAR2)
			!VETT(I)=RAN(IVAR2)
			XPROP(I)=XINF(I)+VETT(I)*(XSUP(I)-XINF(I))
		ENDDO	

		CALL FUNCT(XPROP,N,FPROP)

		NPC = NPC + 1
		NF  = NF  + 1

!	----------------------------
!	  Calcolo media e varianza
!	----------------------------

		RMED  = ((NPC-1.)*RMED+FPROP)/NPC
		RMED2 = ((NPC-2.)*RMED2+FPROP**2)/(NPC-1.)
		VAR2  = RMED2-(RMED**2)*(NPC/(NPC-1.))

		DIFF  = RMED-FOTT

		TRES  = DEXP(-TCOEFF*(RFAL/DIFF)*(FPROP-FOTT))

		IF(FPROP.LE.FOTT) TRES=1.

		IF(ZETA.LE.TRES) THEN

			NUMFAL = 0
			NLM    = NLM+1
			TROVATO = .TRUE.

		ELSE
			NUMFAL = NUMFAL + 1
			RNUMFAL  = RNUMFAL+1.
			ICONTFAL = ICONTFAL+1
			RFAL     = RFAL+1./RNUMFAL
			IDFAL    = 1
		ENDIF

	ENDDO

	RETURN

END SUBROUTINE GENERA_PUNTO_SA
