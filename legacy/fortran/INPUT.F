      subroutine input(tmax,stormlam,celllam,aa,meanrain,
     +                 treelam,treerad,dt,ngridx,ngridy, 
     +                 nsubgridx,nsubgridy,s_init,s_zero,s_fc,
     +                 por,Ksat,b,s_wilt,s_star,Zr,Emax,Ewilt,
     +                 deltarain)
        implicit none
!*********************************************************************
!   Routine:  input.f                                                *
!   Purpose:  Input variables for grid and rainfall generator        *     
!   Written:  August 2004                                            *
!   By:       M.J. Puma                                              *
!*********************************************************************                    
      integer*4, parameter :: prec10 = selected_real_kind(10)
      
      character*12 flname
      integer*4, parameter :: nuin   = 11                 
                                                   
      real(kind = prec10) tmax, stormlam, celllam, aa, meanrain,
     +             treelam, treerad,dt,s_init,s_zero,s_fc,por,Ksat,b,
     +             s_wilt(2),s_star(2),Zr(2),Emax(2),Ewilt(2),
     +             deltarain(2) 
      integer*4 ngridx,ngridy,nsubgridx,nsubgridy 
                                                                                                                                                                                                                                                                                                                                                                                                     
!   Open Main Input File
      open(unit=nuin, file='bucket.in')

!   Read Name of File for Control Parameter Input
      read(nuin,10) flname

!   Control Parameter Input
      call control(flname,tmax, stormlam, celllam,aa, meanrain,
     +                    treelam, treerad, dt, ngridx, ngridy,
     +                    nsubgridx, nsubgridy) 

!   Read Name of File for Soil Input
      read(nuin,10) flname

!   Soil Properties Input
      call soilprop(flname,s_init,s_zero,s_fc,por,Ksat,b) 

!   Read Name of File for Vegetation Properties
      read(nuin,10) flname

!   Materials Properties Input
      call vegprop(flname,s_wilt,s_star,Zr,Emax,Ewilt,deltarain)

!   End of Input
      close(unit=nuin)

   10 format (a)

      return
      end