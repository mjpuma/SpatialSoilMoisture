      subroutine vegprop(flname,s_wilt,s_star,Zr,Emax,Ewilt,deltarain)
        implicit none 
!*********************************************************************
!   Routine:  vegprop.f                                              *
!   Purpose:  Input vegetation properties                            *
!   Written:  August 2004                                            *
!   By:       M.J. Puma                                              *
!*********************************************************************
      character*12 flname
      integer*4,parameter :: nutemp = 15       
      integer*4,parameter :: prec10 = selected_real_kind(10)                    
      
      real(kind = prec10) s_wilt(2),s_star(2),Zr(2),Emax(2),Ewilt(2),
     +                    deltarain(2) 

!   Open file
      open(unit=nutemp,file=flname)

! GRASS PARAMETERS
!   Read wilting point saturation
      read(nutemp,*) s_wilt(1)

!   Read saturation when stomata closure begins
      read(nutemp,*) s_star(1)

!   Read depth of root zone
      read(nutemp,*) Zr(1)

!   Read maximum daily evapotranspiration
      read(nutemp,*) Emax(1)
      
!   Read evaporation rate at the wilting point 
      read(nutemp,*) Ewilt(1) 
      
!   Read canopy interception depth for grasses 
      read(nutemp,*) deltarain(1)
      
! TREE PARAMETERS
!   Read wilting point saturation
      read(nutemp,*) s_wilt(2)

!   Read saturation when stomata closure begins
      read(nutemp,*) s_star(2)

!   Read depth of root zone
      read(nutemp,*) Zr(2)

!   Read maximum daily evapotranspiration
      read(nutemp,*) Emax(2)
      
!   Read evaporation rate at the wilting point 
      read(nutemp,*) Ewilt(2) 
     
!   Read canopy interception depth for trees 
      read(nutemp,*) deltarain(2)
            
!   This ends control parameter input.

      close(unit=nutemp)
    
      return
      end      
            