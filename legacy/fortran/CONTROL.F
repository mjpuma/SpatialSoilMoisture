      subroutine control(flname,tmax, stormlam, celllam, aa,
     +                    meanrain,treelam, treerad,dt,ngridx,ngridy,
     +                    nsubgridx, nsubgridy)
        implicit none
!**********************************************************************
!   Routine:  control.for                                             *
!   Purpose:  Read CONTROL.IN. One value per line, in this order:     *
!     dt            time step, days (the program steps once a day)    *
!     tmax          last day                                          *
!     ngridx,ngridy grid-block size, meters                           *
!     nsubgridx,    spacing of the underlying 5 m plant grid, meters  *
!     nsubgridy     The water balance is solved every 10th point.     *
!     stormlam      storms per day                                    *
!     celllam       rain cells per square meter                       *
!     aa            1/a', with a' the rain-cell scale in meters       *
!     meanrain      mean cell-center depth, centimeters               *
!     treelam       woody plants per square meter                     *
!     treerad       mean crown radius, meters                         *
!   Written:  August 2004, M.J. Puma                                  *
!**********************************************************************
      integer*4, parameter :: prec10 = selected_real_kind(10)                    

      character*12 flname
      integer*4, parameter :: nutemp = 15
                  
      real(kind = prec10) tmax, stormlam, celllam, aa, meanrain,
     +                    treelam, treerad,dt 
      integer*4 ngridx, ngridy, nsubgridx, nsubgridy


!   Open file
      open(unit=nutemp,file=flname)

!   Read Initial Time Step Size
      read(nutemp,*) dt

!   Read Max. simulation time
      read(nutemp,*) tmax

!   Read dimensions of grid block [m]
      read(nutemp,*) ngridx
      read(nutemp,*) ngridy

!   Read spacing between grid points [m]
      read(nutemp,*) nsubgridx
      read(nutemp,*) nsubgridy 
      
!   Read number of storms per day
      read(nutemp,*) stormlam
      
!   Read number of cells per area
      read(nutemp,*) celllam

!   Read storm cell parameter aa
      read(nutemp,*) aa
                                                      
!!   Read mean radius of storm cells
!      read(nutemp,*) meanrad

!   Read mean depth of storm cells
      read(nutemp,*) meanrain 
      
!   Read number of trees per area
      read(nutemp,*) treelam
                           
!   Read mean radius of trees
      read(nutemp,*) treerad      
            
!   This ends control parameter input.

      close(unit=nutemp)

      return
      end
