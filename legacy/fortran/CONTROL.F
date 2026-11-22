      subroutine control(flname,tmax, stormlam, celllam, aa,
     +                    meanrain,treelam, treerad,dt,ngridx,ngridy,
     +                    nsubgridx, nsubgridy)
        implicit none
!**********************************************************************
!   Routine:  control.for                                             *
!   Purpose:  To read control parameters for the spatial bucket model *
!   Written:  August 2004                                             *
!   By:       M.J. Puma                                               *
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
