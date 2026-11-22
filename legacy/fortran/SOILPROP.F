      subroutine soilprop(flname,s_init,s_zero,s_fc,por,Ksat,b)
        implicit none
!*********************************************************************
!   Routine:  soilprop.f                                             *
!   Purpose:  Read soil properties from flname (SOIL.IN).            *
!   One number per line, in this order:                              *
!     s_init   initial relative saturation, 0 to 1                    *
!     s_zero   hygroscopic saturation Sh                             *
!     s_fc     field capacity                                        *
!     por      porosity, volume of voids per volume of soil          *
!     Ksat     saturated conductivity, centimeters per day           *
!     b        Clapp-Hornberger pore-size index                      *
!   These match Table 3 of the dissertation chapter when centimeters *
!   are converted to millimeters (82.2 cm/d = 822 mm/d).             *
!   Written:  August 2004, M.J. Puma                                 *
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
