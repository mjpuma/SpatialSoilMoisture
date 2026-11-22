      subroutine cunif(lowA,upB,valU)
        implicit none
!***********************************************************************
!     Routine:  cunif.f                                                *
!     Function: Generates continuously uniformly distributed random    *
!               numbers                                                *
!     Written:  August                                                 *
!     By:       M.J. Puma                                              *
!***********************************************************************  
      integer*4,parameter :: prec10 = selected_real_kind(10)
      
      real(kind = prec10) valU, random, dum1 
      integer*4 lowA,upB  
      
      dum1=random()
      valU=(upB-lowA)*dum1 + lowA
        
      return
      end
 