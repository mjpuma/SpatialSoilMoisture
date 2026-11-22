      subroutine locate(xx,n,x,j)
        implicit none
!*********************************************************************
!   Routine:  locate.f                                               *
!   Purpose:  Given a monotonic array xx(1:n) and a value x, return  *
!             j such that xx(j) <= x <= xx(j+1). j = 0 or j = n      *
!             means x is outside the array. bucket2d uses this to    *
!             find which sample point a woody-plant center falls     *
!             nearest, then checks the crown radius in meters.       *
!   Written:  August 2004, M.J. Puma. Algorithm from Numerical       *
!             Recipes in FORTRAN, 2nd edition.                       *
!********************************************************************* 
      integer*4,parameter :: prec10 = selected_real_kind(10)
      integer*4 j,n,jl,jm,ju
      real(kind = prec10) x,xx(n)
   
!  Intitialize lower and upper limits     
      jl = 0
      ju = n+1
  
!  Repeat until test condition 10 is satisfied     
  10  if(ju-jl.gt.1)then
!      Compute a midpoint       
         jm = (ju+jl)/2
!      Replace either a lower or upper limit as appropriate       
         if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
            jl=jm
         else
            ju=jm
         endif 
      
         goto 10
      
      endif
  
!  Set the output    
      if(x.eq.xx(1))then
         j=1
      elseif(x.eq.xx(n))then
         j=n-1
      else
         j=jl
      endif
      
      return
      end   