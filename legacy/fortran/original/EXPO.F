      subroutine expo(xnu,tau)   
        implicit none
!***********************************************************************
!     Routine:  expo.f                                                 *
!     Function: Generates exponentially distributed random numbers     *
!               with mean xnu                                          *
!     Written:  August 2004                                            *
!     By:       M.J. Puma                                              *
!***********************************************************************  
      integer*4,parameter :: prec10 = selected_real_kind(10)
       
      real(kind = prec10) xnu, dum, random, tau 

    1 dum = random()
      if(abs(dum).lt.1.e-9) goto 1
!   tau is the time (distance) btw occurrences
      tau = -xnu*log(1-dum)                        

      return
      end