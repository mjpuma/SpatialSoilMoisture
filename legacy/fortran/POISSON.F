      subroutine poisson(xnu,nx) 
        implicit none
!***********************************************************************
!     Routine:  poisson.f                                              *
!     Function: Generates poisson distributed random numbers           *
!     Written:  August 2004                                            *
!     By:       M.J. Puma                                              *
!***********************************************************************  
      integer,parameter :: prec10 = selected_real_kind(10)
        
      real(kind = prec10) t,xnu,tau 
      integer*4 nx        
           
      nx=0
      t=0
   2  call expo(1./xnu,tau)
      t=t+tau
        
      if(t.lt.1)then
         nx = nx + 1
         goto 2
      endif
          
      return
      end