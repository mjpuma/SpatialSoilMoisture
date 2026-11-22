      subroutine expo(xnu,tau)
        implicit none
!***********************************************************************
!     Exponential deviate with mean xnu.
!
!     Used for the time between storms (mean 1/stormlam, days), the
!     rainfall depth at a cell center (mean rain, centimeters), and the
!     radius of a woody crown (mean treerad, meters).
!
!     tau = -xnu * log(1-U), with U uniform on (0, 1). U is redrawn if
!     it is so close to 0 or 1 that the log would overflow.
!
!     Written: August 2004, M.J. Puma
!***********************************************************************
      integer*4,parameter :: prec10 = selected_real_kind(10)

      real(kind = prec10) xnu, dum, tau

    1 call random_number(dum)
      if(abs(dum).lt.1.e-9) goto 1
      if(abs(1.0-dum).lt.1.e-9) goto 1
      tau = -xnu*log(1.0-dum)

      return
      end
