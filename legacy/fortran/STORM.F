      subroutine storm(stormlam,numstorm,tmax,ideltat,nextstorm)
        implicit none 
!***********************************************************************
!     Routine:  storm.f                                                 *
!     Function: Storm arrival days over [1, tmax].                       *
!     Interarrival times are exponential with mean 1/stormlam (days).  *
!     Each interarrival is rounded to the nearest integer day. A gap   *
!     that rounds to 0 is drawn again, so two storms do not fall on    *
!     the same day. nextstorm(j) is the day of storm j. numstorm is    *
!     how many storms fall before tmax.                                 *
!     Written:  August 2004, M.J. Puma                                  *
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