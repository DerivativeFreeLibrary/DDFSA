      subroutine linesearchbox_cont(n,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                                 alfa_max,iprint,bl,bu)
      
!      SUBROUTINE LINESEARCH_INVERTI(N,X,F,D,ALFA,DOLDALFA,Z,FZ,
!     *I_CORR,NUM_FAL,DALFAMAX,IPRINT)
     
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int
      real*8 :: delta,delta1,fpar,fzdelta

      COMMON /NUMNEW/NF
	  
	  gamma=1.d-6    

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0


      j=i_corr

	  if(iprint.ge.1) then
			write(*,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
			write(1,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif


	  if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
			alfa=0.d0
			if(iprint.ge.1) then
				 write(*,*) '  alfa piccolo'
				 write(1,*) '  alfa piccolo'
				 write(*,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
				 write(1,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
			endif
			return
	  endif
      

	  do ielle=1,2

		 if(d(j).gt.0.d0) then

		     if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
   			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
			    alfa=bu(j)-x(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *'
					   write(1,*) ' punto espan. sulla front. *'
				endif
			 endif

		  else

			 if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
				alfa=x(j)-bl(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *'
					   write(1,*) ' punto espan. sulla front. *'
				endif
			 endif

		  endif

		  if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
  
			 d(j)=-d(j)
			 i_corr_fall=i_corr_fall+1
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta per alfa piccolo'
				   write(1,*) ' direzione opposta per alfa piccolo'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
				   write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
				   write(1,*) ' alfa=',alfa,'    alfamax=',alfa_max
			  endif
			  alfa=0.d0
			  cycle

		  endif

		  alfaex=alfa

		  z(j) = x(j)+alfa*d(j)
	  
		  call funct(z,n,fz)

		  nf=nf+1

		  if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
				write(1,*) ' fz =',fz,'   alfa =',alfa
		  endif
		  if(iprint.ge.2) then
			  do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
				  write(1,*) ' z(',i,')=',z(i)
			  enddo
		  endif

		  fpar= f-gamma*alfa*alfa


		  if(fz.lt.fpar) then


			 do

				  if((ifront.eq.1)) then

			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif
				     alfa_d(j)=delta*alfa

				     return

				 end if

				 if(d(j).gt.0.d0) then
							
					 if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=bu(j)-x(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 else

					 if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=x(j)-bl(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 endif
						 
				 z(j) = x(j)+alfaex*d(j) 
				   
     
				 
			     call funct(z,n,fzdelta)
							      
				
				 nf=nf+1

				 if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
					  write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
				 endif
				 if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
						 write(1,*) ' z(',i,')=',z(i)
					  enddo
				 endif

				 fpar= f-gamma*alfaex*alfaex

				 if(fzdelta.lt.fpar) then

					 fz=fzdelta
					 alfa=alfaex

				 else               

					 alfa_d(j)=delta*alfa
!					 alfa_d(j)=alfa
			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto fz =',fz,'   alfa =',alfa
				         write(1,*) ' accetta punto fz =',fz,'   alfa =',alfa
			         endif
					 return
				 end if

		     enddo
		  else      

			 d(j)=-d(j)
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   write(1,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   write(1,*) ' j =',j,'    d(j) =',d(j)
			 endif

		  endif      
			  
	  enddo     

	  if(i_corr_fall.eq.2) then
			 alfa_d(j)=alfa_d(j)
	  else
			 alfa_d(j)=delta*alfa_d(j)
	  end if

	  alfa=0.d0

	  if(iprint.ge.1) then
			write(*,*) ' fallimento direzione'
			write(1,*) ' fallimento direzione'
	  endif

	  return      
	  
      end

