      Program bucket2d
        implicit none
!***********************************************************************
!     Routine:  bucket2d.f                                             *
!     Function: Main routine for the spatial bucket model              *
!     Written:  August 2004                                            *
!     By:       M.J. Puma                                              *
!***********************************************************************
!     What this program does
!     ----------------------
!     One grid block of a land-surface model. Each sample point is grass
!     (type 1) or woody (type 2) and has its own root zone. Storms are a
!     Poisson process in time. Each storm is a set of Gaussian rain cells.
!     Between storms the root zone dries by the analytical solution of
!     Laio et al. (2001). Rain is added at the end of the storm day, after
!     interception. Saturation is held at 1; the excess is runoff.
!
!     How to run
!     ----------
!     From this directory, with gfortran:
!         gfortran -std=legacy -O2 -o bucket2d *.F
!     The La Copita inputs (CONTROL.IN, 30 km, 20200 days) are the
!     dissertation experiment and need a long run. To check that the
!     program works, use the short inputs:
!         cd check && ../bucket2d
!     See README.md in this directory for the file list, the units, and
!     what to look for in the output. Depths in the input files are
!     centimeters. The Python model in src/ uses millimeters.
!
!     Daily order, for every sample point
!     ------------------------------------
!     1. Record evapotranspiration and leakage from today's saturation.
!        These are instantaneous rates, not the millimeter total for the
!        day. et.out and leak.out are those rates, converted back to a
!        depth per day (cm/d).
!     2. Dry the soil for one day.
!     3. On a storm day, add today's rainfall after subtracting canopy
!        interception, and hold saturation at 1.
!***********************************************************************
      integer*4,parameter :: prec10 = selected_real_kind(10)
       
      integer*4,parameter :: nuin=11 
      integer*4,parameter :: nuout=12 
      integer*4,parameter :: nusout=13
      integer*4,parameter :: nuetout=14
      integer*4,parameter :: nuleakout=21
      integer*4,parameter :: nurainout=22
      
      integer*4,parameter :: nsparrain = 23
      integer*4,parameter :: nsparsat = 24
      integer*4,parameter :: nsparet = 25
      integer*4,parameter :: nsparetn = 26

      integer*4,parameter :: nvarsout=33
      integer*4,parameter :: nvaretout=34
      integer*4,parameter :: nvarrainout=37
      integer*4,parameter :: nvarleakout=39
      
	  integer*4,parameter :: nuwetout=38
      
      integer*4,parameter :: nxydim = 3600
      integer*4,parameter :: ntdim = 35000
      
      integer*4 i,j,k,m,p,ideltat,numstorm,
     +        ncell,ngridx,ngridy,nsubgridx,nsubgridy,t,nx,ny, 
     +        iveggrid(nxydim,nxydim),nextstorm(ntdim),
     +        nsubxspar, nsubyspar                      

      real(kind = prec10) totrain(nxydim,nxydim),
     +     rainspar(nxydim,nxydim),S(nxydim,nxydim),
     +     scell(nxydim,nxydim),et(nxydim,nxydim),
     +     leak(nxydim,nxydim),  
     +     actrain(nxydim),etnout(nxydim,nxydim),
     +     rgcorx(nxydim,nxydim),rgcory(nxydim,nxydim),
     +     countree,countgras
                 
      real(kind = prec10) sin,Snew,s_init,s_zero,s_fc,por,Ksat,b,
     +     etn,celEwilt,celEmax,celZr,sstar,swilt,
     +     tmax,stormlam,celllam,rintrcpt,
     +     aa, meanrain,treelam,treerad,dt,half 
          
      real(kind = prec10) s_wilt(2),s_star(2),Zr(2),Emax(2),Ewilt(2),
     +                    deltarain(2) 
      
      real(kind = prec10) eta,eta_w,kk,beta,mm,
     +                          temp,tsfc,tsstar,tswilt,BB  
      
      real(kind = prec10) xyavgrain(nxydim),xyavgsnew(nxydim),
     + 		xyavget(nxydim),xyavgleak(nxydim),xyvarrain(nxydim),
     + 		xyvarsnew(nxydim),xyvarleak(nxydim),
     +      xyvaret(nxydim),fracgras(nxydim),fractree(nxydim)
     
      real(kind = prec10)sumw_S,m_S,tot_S,q_S,temp_S,r_S,
     + 		sumw_et,m_et,tot_et,q_et,temp_et,r_et,
     + 		sumw_leak,m_leak,tot_leak,q_leak,temp_leak,r_leak,
     + 		sumw_rain,m_rain,tot_rain,q_rain,temp_rain,r_rain
          
      real(kind = prec10)  cellx(nxydim),celly(nxydim),raindepth
                           
      real(kind = prec10) treex,treey,arrayx(nxydim),
     +                     arrayy(nxydim),deltax,deltay,canopy,
     +                     canx,cany,dist,distx,disty,
     +                     clx,cly        
      
      integer*4 numtree,xloc,yloc,xlow,xhi,ylow,yhi
      
      integer*4 numavgblcks,new_nrg,lowx,hix,lowy,hiy,nspars, 
     +          nxspar,nyspar,iprtflag,iprtflag2,iprtflag3
     
      real(kind = prec10) wetfrac,sfccount,area
                       
! -------------------------------------------------------------------

!   Open read file and files to save specific pieces of output
      open(nuin,file='bucket.in')
      open(nuout,file='vegmap.out')
      
      open(nusout,file='sat.out')
      open(nuetout,file='et.out') 
      open(nuleakout,file='leak.out') 
      open(nurainout,file='rain.out')
      open(nuwetout,file='wetf.out')
      
      open(nvarsout,file='varsat.out')
      open(nvaretout,file='varet.out')
      open(nvarrainout,file='varrain.out')
      open(nvarleakout,file='varleak.out')
            
      open(nsparrain,file='sparrain.out')
      open(nsparsat,file='sparsat.out')
      open(nsparet,file='sparet.out')
      open(nsparetn,file='sparetn.out')
            
!   Input data
      call input(tmax, stormlam,celllam,aa,meanrain,
     +                 treelam,treerad,dt,ngridx,ngridy, 
     +                 nsubgridx,nsubgridy,s_init,s_zero,s_fc,
     +                 por,Ksat,b,s_wilt,s_star,Zr,Emax,Ewilt,
     +                 deltarain)
      
      celllam = celllam*ngridx*ngridy*9
      treelam = treelam*ngridx*ngridy
      nspars = 10
      
!   Initialize variables 
!    ngridX is the length in the x-direction
!    nsubgridx is the length of cell in x-direction
!    nx is the number of cells in the x-direction
!    nxspar is the number of cells sampled in the x-direction
      nx = ngridX/nsubgridx
      ny = ngridY/nsubgridy  
      nxspar = nx/nspars
      nyspar = ny/nspars
      nsubxspar = nsubgridx*nspars
      nsubyspar = nsubgridy*nspars
      half = nxspar/2
      numavgblcks = ceiling(half)
 
!    Initialize grid
      do i=1,nxspar
         do j=1,nyspar 
            rgcorx(i,j)=0.0
            rgcory(i,j)=0.0 
         enddo
      enddo 

!    Initialize variables
      t = 0   
      numstorm=0 
      ncell = 0 
      iprtflag = 0 
      iprtflag2 = 0
      iprtflag3 = 0   
      
	  do i=1,nxspar         
         do j=1,nyspar 
            totrain(i,j)= 0.0
         enddo
      enddo    
      
	  do i = 1,1000
         nextstorm(i) = 0
      enddo
 
!   Get coordinates of the subgrid cells (raingages - rg)   
!    Set first coordinate points
      do i = 1,nyspar
         rgcorx(i,1)=nsubxspar/2.
      enddo      
      do i = 1,nxspar 
         rgcory(1,i)=nsubyspar/2.
      enddo

!    Create grid              
      do i=1,nxspar 
         do j=2,nyspar     
            rgcorx(i,j)=rgcorx(i,j-1)+nsubxspar            
         enddo
      enddo             
      do j=1,nyspar
         do i=2,nxspar
            rgcory(i,j)=rgcory(i-1,j)+nsubyspar
         enddo
      enddo
                             
!   Assign soil and vegetation properties to horizontal grid               
      do i=1,nxspar 
         do j=1,nyspar
            iveggrid(i,j) = 1
         enddo
      enddo        
      
      do i=1,nxspar 
         do j=1,nyspar
            scell(i,j) = s_init 
         enddo
      enddo       
                                                
      call poisson(treelam,numtree)
              
      do i=1,nxspar
         arrayx(i) = rgcorx(1,i)
      enddo
      
      do i=1,nyspar
         arrayy(i) = rgcory(i,1)
      enddo   
      
      do j=1,numtree
!       Find tree location within the X and Y limits        
         call cunif(0,ngridX,treex)  
         call cunif(0,ngridY,treey) 
         call expo(treerad,canopy)
          
         call locate(arrayx,nxspar,treex,xloc)
         call locate(arrayy,nyspar,treey,yloc)
         
         canx=canopy/nsubgridx
         cany=canopy/nsubgridy
         if(xloc.eq.0)then
            xlow =  1
            xhi = 1+ceiling(canx) 
            xloc = 1                  
         elseif (xloc.eq.nxspar)then
            xlow =  nxspar-1
            xhi = nxspar                  
         else
            xlow = xloc - ceiling(canx)
            xhi =  xloc + ceiling(canx)
            if (xlow.le.0) xlow = 1
            if (xhi.gt.nxspar) xhi = nxspar
         endif   
         

         if(yloc.eq.0)then
            ylow =  1
            yhi = 1+ceiling(cany) 
            yloc = 1           
         elseif (yloc.eq.nyspar)then
            ylow =  nyspar-1
            yhi = nyspar
         else
            ylow =  yloc - ceiling(cany) 
            yhi =   yloc + ceiling(cany) 
            
            if (ylow.le.0) ylow = 1
            if (yhi.gt.nyspar) yhi = nyspar
         endif
                        
         do i=xlow,xhi
            do k=ylow,yhi
               deltax = rgcorx(1,i)-treex
               deltay = rgcory(k,1)-treey
               dist=sqrt(deltax**2+deltay**2)           
               if(dist.lt.canopy)then
!                 Type 2 is the woody plant. Points start as grass (type 1).
                  iveggrid(i,k) = 2
               endif
            enddo
         enddo   
         
      enddo 

!   Print vegetation spatial distribution to file
      do i = 1,nxspar
         write(nuout,5) (iveggrid(i,j),j=1,nyspar)
      enddo
  
  5   format(3000i10)         
                   
      if(nxspar/2.eq.numavgblcks)then
         lowx=half
         hix=half+1
         lowy=half
         hiy=half+1 
         new_nrg = 4
      else
         lowx=half
         hix=half
         lowy=half
         hiy=half 
         new_nrg = 1
      endif      
                                             
      do k=1,numavgblcks 
         countgras = 0
         countree = 0
         do i=lowx,hix
            do j=lowy,hiy 
               if (iveggrid(i,j).eq.1)then
                  countgras = countgras + 1
               else
                  countree = countree + 1
               endif    
            enddo
         enddo 
         fracgras(k)=countgras/(countgras+countree)
         fractree(k)=countree/(countgras+countree)
         lowx=lowx-1
         hix=hix+1
         lowy=lowy-1
         hiy=hiy+1      
      enddo
           
!   Compute the arrival times of storms
      call storm(stormlam,numstorm,tmax,ideltat,nextstorm)

!*************************
!     Begin Time Loop
!*************************

!   New Time Step                             
      m = 1
  
  100 if(t < tmax)then  
         t = t+dt
         sfccount = 0. 
         if (t >= nextstorm(m))then
            m = m+1    
            
            do i = 1,nxydim
               actrain(i) = 0.
               cellx(i) = 0.
               celly(i)=0.
            enddo
            
            raindepth = 0.      
!!            Compute a spatially uniform rainfield  
!            call expo(meanrain,actrain) ! Find rainfall depth    
!            do i=1,nxspar 
!               do j=1,nyspar
!                  totrain(i,j) = actrain(1)
!               enddo
!             enddo

!         Compute a spatially heterogeneous rainfield  
!         Compute the number of cells for the current storm
            call poisson(celllam,ncell) 
                            
!         Compute the rainfall depth and coordinates for each cell            
            do p=1,ncell 
!              Find cell location within the X and Y limits
               call cunif(0,ngridX*3,cellx(p)) 
               call cunif(0,ngridY*3,celly(p))
!              Find rainfall depth at the cell center               
               call expo(meanrain,actrain(p))  
            enddo              

!         Rainfall depth for this storm only, in centimeters.
!         Each cell adds h * exp(-2 a^2 r^2). Cells are placed on a
!         domain three times as wide as the grid, and the grid sits in
!         the middle, so the edges are not systematically dry.
            do i=1,nxspar
               do j=1,nyspar
                  totrain(i,j) = 0.0
               enddo
            enddo

!         Compute rainfall at each grid point   
            do k = 1,ncell         
               do i=1,nxspar 
                  do j=1,nyspar
                     clx = cellx(k)
                     cly = celly(k)
                     distx=(rgcorx(i,j)+ ngridX)-clx
                     disty=(rgcory(i,j)+ ngridY)-cly
                  
                     dist = sqrt(distx**2 + disty**2)               
                     raindepth = actrain(k)*exp(-2*(aa**2)*(dist**2))
                     
                     if (raindepth.ge.0.001)then
                         totrain(i,j) = totrain(i,j) + raindepth
                     endif
                  
                  enddo
               enddo
            enddo

!         Calculate saturation, transpiration, evaporation and leakage
!         for a sparse matrix of grid points
            do i=1,nxspar 
               do j=1,nyspar
                  S(i,j) = scell(i,j) 
                  swilt = s_wilt(iveggrid(i,j))
                  sstar = s_star(iveggrid(i,j))
                  celZr = Zr(iveggrid(i,j)) 
                  celEmax = Emax(iveggrid(i,j))
                  celEwilt = Ewilt(iveggrid(i,j))
                  rintrcpt = deltarain(iveggrid(i,j))
                  sin = scell(i,j)

!               Normalization of parameters by nZr
                  eta = celEmax/(por*celZr)
                  eta_w = celEwilt/(por*celZr)
                  kk = Ksat/(por*celZr)
                  beta = 2*b+4                   
                  mm = kk/(exp(beta*(1-s_fc))-1)

                 if (sin <= s_zero)then
                     Snew = sin
                     sin = 0
                  endif

                  if (sin > s_fc)then
                     temp = mm*exp(beta*(sin-s_fc))
                     tsfc = (1/(beta*(mm-eta)))*( beta*(s_fc-sin)+ 
     +               log((eta-mm+temp)/eta) )
           
                     if (dt <= tsfc) then
                        Snew = sin - (1/beta)*(log( 
     +                        ((eta-mm+temp)*exp(beta*(eta-mm)
     +                        *dt)-temp )/(eta-mm) ))
                        sin = 0
                     else
                        dt = dt - tsfc
                        sin = s_fc
                     endif 
           
                  endif

                  if (sin > sstar)then
                     tsstar = (sin-sstar)/eta
           
                     if (dt <= tsstar)then
                        Snew = sin-eta*dt
                        sin = 0 
                     else
                        dt = dt-tsstar
                        sin = sstar
                     endif
           
                  endif


                  if (sin > swilt)then                    
                     tswilt = -((sstar-swilt)/(eta-eta_w))*log(
     +                        eta_w/(eta_w+(eta-eta_w)*(sin-swilt)
     +                        /(sstar-swilt)))    
     
                
                     if (dt <= tswilt)then
                        BB = (eta-eta_w)/(sstar-swilt)
                        Snew = (swilt-(eta_w/BB))*(1-exp(-BB*dt))+
     +                          sin*exp(-BB*dt) 
                        sin = 0
                     else
                        dt = dt-tswilt
                        sin = swilt
                     endif
        
                  endif
        
                  if (sin > s_zero)then
                     Snew = s_zero + (sin-s_zero)*
     +                      exp(-(eta_w/(swilt-s_zero))*dt)
                     
                  endif
                   
                  if (S(i,j) > sstar)then 
                     etn = eta
                  elseif( S(i,j) > swilt)then
                     etn=(eta-eta_w)*(S(i,j)-swilt)/
     +                         (sstar-swilt)+eta_w 
                  elseif(S(i,j) > s_zero)then 
                     etn = eta_w*((S(i,j)-s_zero)/(swilt-s_zero)) 
                  else
                     etn = 0 
                  endif

                  etnout(i,j) = etn                  
                  et(i,j) = etn*por*celZr 
                  
                  if(S(i,j) > s_fc)then
                     leak(i,j) = (Ksat/(exp(beta*(1-s_fc))-1))
     +                            *(exp(beta*(S(i,j)-s_fc))-1)
                     sfccount = sfccount+1    
                  else
                     leak(i,j) = 0
                  endif
                     
                  dt = 1
                  rainspar(i,j) =  totrain(i,j)-rintrcpt
                  if (rainspar(i,j).lt.1e-5)then 
                     rainspar(i,j) = 0.
                  endif
                  scell(i,j) = min(1.,(rainspar(i,j)/(por*celZr)+Snew))
    
               enddo
              
            enddo
         
         else 
            ncell=0
            do i=1,nxspar 
               do j=1,nyspar      
                  S(i,j) = scell(i,j) 
                  swilt = s_wilt(iveggrid(i,j))
                  sstar = s_star(iveggrid(i,j))
                  celZr = Zr(iveggrid(i,j)) 
                  celEmax = Emax(iveggrid(i,j))
                  celEwilt = Ewilt(iveggrid(i,j))
                  sin = scell(i,j)  

!               Normalization of parameters by nZr
                  eta = celEmax/(por*celZr)
                  eta_w = celEwilt/(por*celZr)
                  kk = Ksat/(por*celZr)
                  beta = 2*b+4                   
                  mm = kk/(exp(beta*(1-s_fc))-1)
                  

                 if (sin <= s_zero)then
                     Snew = sin
                     sin = 0
                  endif

                  if (sin > s_fc)then
                     temp = mm*exp(beta*(sin-s_fc))
                     tsfc = (1/(beta*(mm-eta)))*( beta*(s_fc-sin)+ 
     +               log((eta-mm+temp)/eta) )
           
                     if (dt <= tsfc) then
                        Snew = sin - (1/beta)*(log( 
     +                        ((eta-mm+temp)*exp(beta*(eta-mm)
     +                        *dt)-temp )/(eta-mm) ))
                        sin = 0
                     else
                        dt = dt - tsfc
                        sin = s_fc
                     endif 
           
                  endif

                  if (sin > sstar)then
                     tsstar = (sin-sstar)/eta
           
                     if (dt <= tsstar)then
                        Snew = sin-eta*dt
                        sin = 0 
                     else
                        dt = dt-tsstar
                        sin = sstar
                     endif
           
                  endif


                  if (sin > swilt)then                    
                     tswilt = -((sstar-swilt)/(eta-eta_w))*log(
     +                        eta_w/(eta_w+(eta-eta_w)*(sin-swilt)
     +                        /(sstar-swilt)))    
     
                
                     if (dt <= tswilt)then
                        BB = (eta-eta_w)/(sstar-swilt)
                        Snew = (swilt-(eta_w/BB))*(1-exp(-BB*dt))+
     +                          sin*exp(-BB*dt) 
                        sin = 0
                     else
                        dt = dt-tswilt
                        sin = swilt
                     endif
        
                  endif
        
                  if (sin > s_zero)then
                     Snew = s_zero + (sin-s_zero)*
     +                      exp(-(eta_w/(swilt-s_zero))*dt)
                     
                  endif
                   
                  if (S(i,j) > sstar)then 
                     etn = eta
                  elseif( S(i,j) > swilt)then
                     etn=(eta-eta_w)*(S(i,j)-swilt)/
     +                         (sstar-swilt)+eta_w 
                  elseif(S(i,j) > s_zero)then 
                     etn = eta_w*((S(i,j)-s_zero)/(swilt-s_zero)) 
                  else
                     etn = 0 
                  endif

                  etnout(i,j) = etn                  
                  et(i,j) = etn*por*celZr

                  if(S(i,j) > s_fc)then
                     leak(i,j) = (Ksat/(exp(beta*(1-s_fc))-1))
     +                            *(exp(beta*(S(i,j)-s_fc))-1)
                     sfccount = sfccount+1  
                  else
                     leak(i,j) = 0 
                  endif                  
                  
                  dt = 1
                  rainspar(i,j) = 0.0
                  totrain(i,j) = 0.0 
                  scell(i,j) = Snew
                  
               enddo   
            enddo
         endif
!!      Computer the mean and variance in space for the different 
!!      averaging subregions                                                                                   
!         if (nxspar.eq.nyspar)then  
!            half = nxspar/2
!            numavgblcks = ceiling(half)
         
!         First pass to get the mean
!            do k = 1,numavgblcks
!               xyavgrain(k) = 0.0
!               xyavgsnew(k) =  0.0
!               xyavget(k) = 0.0
!               xyavgleak(k) = 0.0
!            enddo              
!            if(nxspar/2.eq.numavgblcks)then
!               lowx=half
!               hix=half+1
!               lowy=half
!               hiy=half+1 
!               new_nrg = 4
!            else
!               lowx=half
!               hix=half
!               lowy=half
!               hiy=half 
!               new_nrg = 1
!            endif      
                                             
!            do k=1,numavgblcks
!               do i=lowx,hix
!                  do j=lowy,hiy 
!                     xyavgrain(k) = xyavgrain(k) + rainspar(i,j)
!                     xyavgsnew(k) =  xyavgsnew(k) + S(i,j)
!                     xyavget(k) = xyavget(k) + et(i,j)
!                     xyavgleak(k) = xyavgleak(k) + leak(i,j)
!                  enddo
!               enddo
!               lowx=lowx-1
!               hix=hix+1
!               lowy=lowy-1
!               hiy=hiy+1      
!            enddo
           
!            do i=1,numavgblcks      
!               xyavgrain(i) = xyavgrain(i)/new_nrg
!               xyavgsnew(i) =  xyavgsnew(i)/new_nrg
!               xyavget(i) = xyavget(i)/new_nrg
!               xyavgleak(i) = xyavgleak(i)/new_nrg
!               new_nrg = ((new_nrg**.5)+2)**2
!            enddo

!      Computer the mean and variance in space for the different 
!      averaging subregions of S, ET, L, and R using algorithm 
!      of West (Updating mean and variance estimates: 
!      an improved method, 1979)                                
  	     if (nxspar.eq.nyspar)then           
            do k = 1,numavgblcks
               xyavgsnew(k) = 0.0
               xyvarsnew(k) = 0.0
               
			   xyavget(k) = 0.0               
               xyvaret(k) = 0.0                

			   xyavgrain(k) = 0.0
               xyvarrain(k) = 0.0

               xyavgleak(k)=0.0              
			   xyvarleak(k) = 0.0  
            enddo
               
            if(nxspar/2.eq.numavgblcks)then
               lowx=half
               hix=half+1
               lowy=half
               hiy=half+1 
               new_nrg = 4
            else
               lowx=half
               hix=half
               lowy=half
               hiy=half 
               new_nrg = 1
            endif      
                                             
!           Means on the window [lowx:hix, lowy:hiy]. The window grows
!           outward from the center. Saturation is weighted by porosity
!           times root depth (a volume). Fluxes are an area mean.
!           lowx is the first row of this window.
            do k=1,numavgblcks
               sumw_S = por*Zr(iveggrid(lowx,lowy))
               m_S = S(lowx,lowy)
               tot_S = 0
               
               sumw_et = 1
               m_et = et(lowx,lowy)
               tot_et = 0
               
               sumw_rain = 1
               m_rain = rainspar(lowx,lowy)
               tot_rain = 0
               
               sumw_leak = 1
               m_leak = leak(lowx,lowy)
               tot_leak = 0
               
			   do j=lowy+1,hiy
         	      q_S = S(lowx,j)-m_S
         		  temp_S = sumw_S + por*Zr(iveggrid(lowx,j))
         		  r_S = q_S*por*Zr(iveggrid(lowx,j))/temp_S
         		  m_S = m_S + r_S
         		  tot_S = tot_S + r_S*sumw_S*q_S
         		  sumw_S = temp_S
                    
         		  q_et = et(lowx,j)-m_et
         		  temp_et = sumw_et + 1
         		  r_et = q_et*1/temp_et
         		  m_et = m_et + r_et
         		  tot_et = tot_et + r_et*sumw_et*q_et
         		  sumw_et = temp_et
         			
         		  q_rain = rainspar(lowx,j)-m_rain
         		  temp_rain = sumw_rain + 1
         		  r_rain = q_rain*1/temp_rain
         		  m_rain = m_rain + r_rain
         		  tot_rain = tot_rain + 
     +                        r_rain*sumw_rain*q_rain
         		  sumw_rain = temp_rain

         		  q_leak = leak(lowx,j)-m_leak
         		  temp_leak = sumw_leak + 1
        		  r_leak = q_leak*1/temp_leak
         		  m_leak = m_leak + r_leak
         		  tot_leak = tot_leak + 
     +              		   r_leak*sumw_leak*q_leak
         		  sumw_leak = temp_leak
               enddo

			   do i=lowx+1,hix
                  do j=lowy,hiy
         			q_S = S(i,j)-m_S
         			temp_S = sumw_S + por*Zr(iveggrid(i,j))
         			r_S = q_S*por*Zr(iveggrid(i,j))/temp_S
         			m_S = m_S + r_S
         			tot_S = tot_S + r_S*sumw_S*q_S
         			sumw_S = temp_S
                    
         			q_et = et(i,j)-m_et
         			temp_et = sumw_et + 1
         			r_et = q_et*1/temp_et
         			m_et = m_et + r_et
         			tot_et = tot_et + r_et*sumw_et*q_et
         			sumw_et = temp_et
         			
         			q_rain = rainspar(i,j)-m_rain
         			temp_rain = sumw_rain + 1
         			r_rain = q_rain*1/temp_rain
         			m_rain = m_rain + r_rain
         			tot_rain = tot_rain + 
     +                         r_rain*sumw_rain*q_rain
         			sumw_rain = temp_rain

         			q_leak = leak(i,j)-m_leak
         			temp_leak = sumw_leak + 1
         			r_leak = q_leak*1/temp_leak
         			m_leak = m_leak + r_leak
         			tot_leak = tot_leak + 
     +               		   r_leak*sumw_leak*q_leak
         			sumw_leak = temp_leak
                  enddo
               enddo
               
               xyavgsnew(k) = m_S
      		   xyvarsnew(k) = tot_S*new_nrg/((new_nrg-1)*sumw_S) 
			   
			   xyavget(k) = m_et
      		   xyvaret(k) = tot_et*new_nrg/((new_nrg-1)*sumw_et) 

			   xyavgrain(k) = m_rain
     		   xyvarrain(k)=tot_rain*new_nrg/((new_nrg-1)*sumw_rain)
			   
			   xyavgleak(k) = m_leak
      		   xyvarleak(k)=tot_leak*new_nrg/((new_nrg-1)*sumw_leak)
                              			   
			   new_nrg = ((new_nrg**.5)+2)**2
              					                    
			   lowx=lowx-1
               hix=hix+1
               lowy=lowy-1
               hiy=hiy+1      
            enddo          
                               
!         Write statistics of the various areas to file
            write(nusout,20) (xyavgsnew(i),i=1,numavgblcks)
            write(nuetout,20) (xyavget(i),i=1,numavgblcks)
            write(nurainout,20) (xyavgrain(i),i=1,numavgblcks)
            write(nuleakout,20) (xyavgleak(i),i=1,numavgblcks)            

            write(nvarsout,20) (xyvarsnew(i),i=1,numavgblcks)
            write(nvaretout,20) (xyvaret(i),i=1,numavgblcks)
            write(nvarrainout,20) (xyvarrain(i),i=1,numavgblcks)
            write(nvarleakout,20) (xyvarleak(i),i=1,numavgblcks)  

  20        format(600f15.10)
         
         endif

!   Write sparse matrix of output to file       

         if ((xyavgsnew(numavgblcks)>.7).and.(iprtflag==0)
     +      		.and.(t>200))then
            do i=1,nxspar
                write(nsparrain,50) (rainspar(i,j),j=1,nyspar)
                write(nsparsat,50) (S(i,j),j=1,nyspar)  
                write(nsparet,50) (et(i,j),j=1,nyspar)
                write(nsparetn,50) (etnout(i,j),j=1,nyspar) 
            enddo     
            iprtflag = 1                     
  50        format(600f15.10) 
         endif     
      
         if ((xyavgsnew(numavgblcks)<.4).and.(iprtflag2==0)
     +   		.and.(t>200))then
            do i=1,nxspar
               write(nsparrain,51) (rainspar(i,j),j=1,nyspar)
               write(nsparsat,51) (S(i,j),j=1,nyspar)  
               write(nsparet,51) (et(i,j),j=1,nyspar)
               write(nsparetn,51) (etnout(i,j),j=1,nyspar) 
            enddo     
            iprtflag2 = 1                        
  51        format(600f15.10) 
         endif           

         if ((xyavgsnew(numavgblcks)>.4).and.(xyavgsnew(numavgblcks)<.5)
     +      		.and.(iprtflag3==0).and.(t>200))then
            do i=1,nxspar
               write(nsparrain,52) (rainspar(i,j),j=1,nyspar)
               write(nsparsat,52) (S(i,j),j=1,nyspar)  
               write(nsparet,52) (et(i,j),j=1,nyspar)
               write(nsparetn,52) (etnout(i,j),j=1,nyspar) 
            enddo     
            iprtflag3 = 1                        
  52        format(600f15.10) 
         endif   
    
!   Fraction of the land surface above field capacity
         area = nxspar*nyspar
         wetfrac = sfccount/area
         write(nuwetout,70) wetfrac

!    Write information to screen
         write(*,65) t,fracgras(numavgblcks),nextstorm(m),
     +    			xyavgsnew(numavgblcks),wetfrac
  65     format(1x,'t= ',i5,', %grass= ',f10.4,', next t= ',i5,
     +                ', <S>=',f10.4,',wetfrac=',f10.4)
  70     format(f15.10)   

!      Individual storm finished                  
         goto 100       

!    Simulation endtime has been reached      
      endif                
      
      stop
      end
