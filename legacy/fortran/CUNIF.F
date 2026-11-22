      subroutine cunif(lowA,upB,valU)
        implicit none
!***********************************************************************
!     Uniform deviate on the interval [lowA, upB].
!
!     lowA and upB are the ends of the interval. In bucket2d they are
!     positions in meters. valU is the returned deviate, in the same
!     units.
!
!     The draw uses the standard Fortran generator, random_number, so
!     the program builds with gfortran. Without a call to random_seed
!     the sequence changes from run to run. A check of the program does
!     not depend on the seed. It depends on the files in check/ being
!     written and on saturations staying inside [0, 1].
!
!     Written: August 2004, M.J. Puma
!***********************************************************************
      integer*4,parameter :: prec10 = selected_real_kind(10)

      real(kind = prec10) valU, dum1
      integer*4 lowA,upB

      call random_number(dum1)
      valU=(upB-lowA)*dum1 + lowA

      return
      end
