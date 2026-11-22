      subroutine storm(stormlam,numstorm,tmax,ideltat,nextstorm)
        implicit none 
!***********************************************************************
!     Routine:  storm.f                                                 *
!     Function: Generates storm arrival times.  That is, the time       *
!               between storms and individual storm velocity.           *
!     Written:  August 2004                                             *
!     By:       M.J. Puma                                               *
!***********************************************************************  
      integer*4,parameter :: prec10 = selected_real_kind(10)
      integer*4,parameter :: ntdim = 35000
                    
      real(kind = prec10) stormlam, deltat, tmax
      integer*4 ideltat,j,nextstorm(ntdim),numstorm
        
!    Calculates the time btw storms from an exponential distribution        
      j = 1

  20  call expo(1./stormlam,deltat)   
      ideltat = anint(deltat)  
      if(ideltat.eq.0) goto 20 
      
      if (j.ne.1)then   
         nextstorm(j) = ideltat + nextstorm(j-1)
      else
         nextstorm(j) = ideltat
      endif
      
      if(nextstorm(j).lt.tmax)then
         numstorm = numstorm + 1
         j = j+1
         goto 20
      endif          
                      
      return
      end