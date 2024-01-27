!===============================================================
!
! Funzione Cosine Mixture  (n=2,4)
!
!===============================================================
!			VARIABLE BOUNDS
!===============================================================
!       -1.D0                1.D0          Limits of variable 1
!       -1.D0                1.D0             "        "      2
!       -1.D0                1.D0             "        "      3
!       -1.D0                1.D0             "        "      4
!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE NINIT(N)
	IMPLICIT NONE
	
	INTEGER				:: N
	DOUBLE PRECISION, POINTER	:: LB(:), UB(:)
	CHARACTER(LEN=40)	:: nomefun
	INTEGER				:: I

	N = 4

	RETURN
END SUBROUTINE NINIT 

!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE FUNCTINIT(N,LB,UB,nomefun,fglob)
	IMPLICIT NONE
	
	INTEGER				:: INDFUN
	INTEGER				:: N
	DOUBLE PRECISION	:: LB(N), UB(N), fglob
	CHARACTER(LEN=40)	:: nomefun
	INTEGER				:: I

	DO I = 1,N
		LB(I) = -1.D0
		UB(I) =  1.D0
	ENDDO

	LB(1:n) = LB(1:n) + 4.0d0*ATAN(1.0d0)/10.0d0
	UB(1:n) = UB(1:n) + 4.0d0*ATAN(1.0d0)/10.0d0

	nomefun = 'Cosine mixture'

	fglob = -0.4d0

	RETURN
END SUBROUTINE FUNCTINIT 
          
!==============================================================================================
!==============================================================================================
!==============================================================================================

SUBROUTINE FUNCT(X,N,F)
	IMPLICIT NONE
	INTEGER			:: N
	DOUBLE PRECISION	:: X(N), F

	CALL COSMIX(X,N,F)

	RETURN

END SUBROUTINE FUNCT 

SUBROUTINE COSMIX(X,N,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION, PARAMETER :: PI = 3.1415923

      INTEGER          :: I

      F = 0.0D0

      DO I=1,N
         F = F + 0.1D0*DCOS(5.0D0*PI*X(I))-X(I)**2
      END DO
      
      F = - F

	  !IF(N == 2) F = F + 0.2d0
	  !IF(N == 4) F = F + 0.4d0

      RETURN

END
