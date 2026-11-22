      subroutine soilprop(flname,s_init,s_zero,s_fc,por,Ksat,b)
        implicit none
!*********************************************************************
!   Routine:  soilprop.f                                             *
!   Purpose:  Input soil properties                                  *
!   Written:  August 2004                                            *
!   By:       M.J. Puma                                              *
!*********************************************************************
      integer*4,parameter :: prec10 = selected_real_kind(10)
      
      character*12 flname
      integer*4,parameter :: nutemp = 15  
                            
      real(kind = prec10) s_init,s_zero,s_fc,por,Ksat,b

!   Open file
      open(unit=nutemp,file=flname)

!   Read Initial Saturation
      read(nutemp,*) s_init

!   Read hygroscopic saturation
      read(nutemp,*) s_zero

!   Read saturation at field capacity
      read(nutemp,*) s_fc

!   Read porosity
      read(nutemp,*) por
      
!   Read saturated hydraulic conductivity
      read(nutemp,*) Ksat
      
!   Read soil water retention curve parameter b (Clapp and Hornberger)
      read(nutemp,*) b   
            
!   This ends control parameter input.

      close(unit=nutemp)
    

      return
      end
