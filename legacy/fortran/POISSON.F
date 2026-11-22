      subroutine poisson(xnu,nx) 
        implicit none
!***********************************************************************
!     Routine:  poisson.f                                              *
!     Function: Poisson count with mean xnu.                            *
!     nx is the count. Used for the number of woody plants in the      *
!     grid block (mean treelam * area) and the number of rain cells    *
!     in one storm (mean celllam * area of the extended domain).       *
!     The count is the number of exponential waits that fit in [0, 1]  *
!     when the wait has mean 1/xnu, which is a Poisson(xnu) draw.      *
!     Written:  August 2004, M.J. Puma                                 *
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