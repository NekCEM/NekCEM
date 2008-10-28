                                
                                
c-----------------------------------------------------------------------
c FIXME  1. include mstep 
c FIXME  2. sdomain shifting 
c FIXME  3. free postion of source
c-----------------------------------------------------------------------
      subroutine cem_wakepotential
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

      integer  e
      integer  icalld
      save     icalld
      data     icalld /0/
   
      !if (ifmovw) return ! FIXME

      if (param(35).eq.7) return

      if (.not.IFWAKE) return

c     if (icalld.eq.0) then

         if (nid.eq.0) then
             write(6,*)        
             write(6,*) 'Setup Wakepotential Calculation >>>>>> '
         endif
         icalld = 1

         call total_charge         
         call interpolation_init         
         call interpolation_init_glo     
         call parallel_zframe_index           
          
         return  ! FIXME : temporary for moving window
         call extract_integral_path_glo

         !if (IFNAPY) call get_napoly_path      

c     endif


      if (mod(istep,mstep).eq.0) then

c ... longitudinal wakepoten 
      call perimeter_integral 
     $ (Ez_int,lwght,nzlev,indxw,izlvi,swght,n_wp_grid,Ez,workwp,zind)
      call interpolation_path_gll_to_s(Ez_wp,Ez_int)           
      !if (IFNAPY) call napoly_interpolation          

c ... transverse wakepoten
      if  (IFWAKET) then
          call perimeter_integral 
     $(Ez_int1,lwght1,nzlev,indxw1,izlvi1,swght1,n_wp_grid1,Ez,workwp1)
          call interpolation_path_gll_to_s(Ez_wp1,Ez_int1)           
      endif

      endif

      call wake_longitudinal                      
      
      return
      end
c-----------------------------------------------------------------------
c ... parallel
      subroutine napoly_interpolation                                  
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'
      integer  e

      nxyz =nx1*ny1*nz1

      do iz=1,nz1           
      do i =1,napmax(1)  
         lel =0
         l2=0
         l5=0
         tmp2=0
         tmp5=0
         rr=rnp(1,i)
         do k= 1,nwnpel(1)  
            e= wnpel(1,k)
            rmin= wnprad1(1,e)
            rmax= wnprad1(2,e)
            tmin= wnptht1(1,e)
            tmax= wnptht1(2,e)
            if ((rmin.le.rr).and.(rr.le.rmax)) then
              lel=lel+1

              if     (eltdir(e).EQ.'RADIS') then  
                !grids starting in radial direction
                !interpolation in radial first
                 rint = 2.0/(rmax-rmin)*(rr-rmin)-1.0 
                 call igllm (ju,jut,zi,rint,nz1,1,nz1,1)    ! define jut ??? 
                 do iy= 1,ny1
                    call mxm (Ez(1,iy,iz,e),1,jut,ny1,Eznp,1)
                    Ez_np1(iy)=Eznp
                 enddo 
                 do itheta= 1,nx1              
                    dtheta=(tmax-tmin)/nx1
                    theta = tmin+dtheta/2.0+dtheta*(itheta-1)
                    tint  = 2.0/(tmax-tmin)*(theta-tmin)-1.0
                    call igllm (ju,jut,zi,tint,nz1,1,nz1,1) ! define jut ??? 
                    call mxm (Ez_np1,1,jut,ny1,Eznp,1)
                    Ez_np2(itheta)=Eznp
                    l2=l2+1
                    tmp2=tmp2+Eznp
                 enddo      
                   
              elseif (eltdir(e).EQ.'THETA') then 
                 rint = 2.0/(rmax-rmin)*(rr-rmin)-1.0 
                 do ix=1,nx1
                    do iy=1,ny1
                       Ez_np3(iy)=Ez(ix,iy,iz,e)
                    enddo
                    call igllm (ju,jut,zi,rint,nz1,1,nz1,1)    ! define jut ??? 
                    call mxm (Ez_np3,1,jut,ny1,Eznp,1)
                    Ez_np4(ix)=Eznp
                 enddo 
                 do itheta= 1,nx1        
                    dtheta=(tmax-tmin)/nx1
                    theta = tmin+dtheta/2.0+dtheta*(itheta-1)
                    tint  = 2.0/(tmax-tmin)*(theta-tmin)-1.0
                    call igllm (ju,jut,zi,tint,nz1,1,nz1,1) ! define jut ??? 
                    call mxm (Ez_np4,1,jut,ny1,Eznp,1)
                    Ez_np5(itheta)=Eznp
                    l5=l5+1
                    tmp5=tmp5+Eznp
                 enddo      
               else
                   if (nid.eq.0) write(6,*) 'ERR: check eltdir'
                   call exitt

              endif  ! radis or theta

            endif    ! rmin< rr
         enddo
      enddo
         tmp2=  glsum(tmp2,1)
         tmp5=  glsum(tmp5,1)
         lel = iglsum( lel,1)
         !write(6,*) lel,lelt,napmax(1),l2+l5,nwnpel(1)
         Ez_znp(1,iz)=tmp5+tmp2
         Ez_znp(1,iz)=Ez_znp(1,iz)*2*pi/(lel*nx1) 
      enddo

      do iz = 1,nz1           
      do i  = 1,napmax(2)  
         lel= 0
         rr = rnp(2,i)
         l2 = 0
         l5 = 0
         tmp2= 0
         tmp5= 0

         do k= 1,nwnpel(2) 
            e= wnpel(2,k)
            rmin= wnprad2(1,e)
            rmax= wnprad2(2,e)
            tmin= wnptht2(1,e)
            tmax= wnptht2(2,e)

            if ((rmin.le.rr).and.(rr.le.rmax)) then
               lel=lel+1
               if     (eltdir(e).EQ.'RADIS') then  
                !grids starting in radial direction
                !interpolation in radial first
                 rint = 2.0/(rmax-rmin)*(rr-rmin)-1.0 
                 call igllm (ju,jut,zi,rint,nz1,1,nz1,1)    ! define jut ??? 
                 do iy= 1,ny1
                    call mxm (Ez(1,iy,iz,e),1,jut,ny1,Eznp,1)
                    Ez_np1(iy)=Eznp
                 enddo 
                 do itheta= 1,nuniform
                    dtheta=(tmax-tmin)/nx1
                    theta = tmin+dtheta/2.0+dtheta*(itheta-1)
                    tint  = 2.0/(tmax-tmin)*(theta-tmin)-1.0
                    call igllm (ju,jut,zi,tint,nz1,1,nz1,1) ! define jut ??? 
                    call mxm (Ez_np1,1,jut,ny1,Eznp,1)
                    Ez_np2(itheta)=Eznp
                    l2=l2+1
                    tmp2=tmp2+Eznp
                 enddo      
                   
               elseif (eltdir(e).EQ.'THETA') then 
                 rint = 2.0/(rmax-rmin)*(rr-rmin)-1.0 
                 do ix=1,nx1
                    do iy=1,ny1
                       Ez_np3(iy)=Ez(ix,iy,iz,e)
                    enddo
                    call igllm (ju,jut,zi,rint,nz1,1,nz1,1)    ! define jut ??? 
                    call mxm (Ez_np3,1,jut,ny1,Eznp,1)
                    Ez_np4(ix)=Eznp
                 enddo 
                 do itheta= 1,nuniform
                    dtheta=(tmax-tmin)/nx1
                    theta = tmin+dtheta/2.0+dtheta*(itheta-1)
                    tint  = 2.0/(tmax-tmin)*(theta-tmin)-1.0
                    call igllm (ju,jut,zi,tint,nz1,1,nz1,1) ! define jut ??? 
                    call mxm (Ez_np4,1,jut,ny1,Eznp,1)
                    Ez_np5(itheta)=Eznp
                    l5=l5+1
                    tmp5=tmp5+Eznp
                 enddo      
               else
                   if (nid.eq.0) write(6,*) 'ERR: check eltdir'
                   call exitt

              endif  ! radis or theta

            endif    ! rmin< rr
         enddo
      enddo
         !write(6,*) lel,l2+l5,nwnpel(2)
         tmp2=  glsum(tmp2,1)
         tmp5=  glsum(tmp5,1)
         lel = iglsum( lel,1)
         Ez_znp(2,iz)=tmp5+tmp2
         Ez_znp(2,iz)=Ez_znp(1,iz)*2*pi/(lel*nx1) 
      enddo

      ! interpolate in znp 
      call igllm  (ju,jut,zi,zu,nz1,nuniform,nz1,nuniform)

      do i=1,2
         call mxm (Ez_znp(i,1),1,jut,nz1,Ez_np(i,1),nuniform)
      enddo

      do i=1,nuniform
         Ez_wp(npindx(1,i))=Ez_np(1,i)  ! overwite
         Ez_wp(npindx(2,i))=Ez_np(2,i)
      enddo

      ! get global napoly index for Ez_r into ez_wp

      return
      end
c-----------------------------------------------------------------------
c ... parallel
      subroutine get_napoly_path                                        
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'
      integer  e
c
c     + 
c      . r1
c        . r2
c          . r3
c            .
c              .
c                . r_nuniform
c     +-----------+              
c     wnapy1   wzone1                       
c
c     find rmax, rmin
c     define znp1,znp2,rnp1,rnp2

      !if (.NOT.IFNAPY) return

      if (nid.eq.0)  write(6,*)         

      nxyz = nx1*ny1*nz1
      rnapmax = -1.e+20
      rnapmin =  1.e+20
      do iel = 1,nwnpel(1)
         e = wnpel(1,iel)
         do i=1,nxyz
            xx= xm1(i,1,1,e)
            yy= ym1(i,1,1,e)
            rr = xx**2+yy**2
            theta=atan2(yy,xx)
            if (rr.gt.0) rr=sqrt(rr)
            rnapmax=max(rr,rnapmax)
            rnapmin=min(rr,rnapmin)
         enddo
      enddo
      rnapmax=glmax(rnapmax,1)  ! global rmax in napoly zone
      rnapmin=glmin(rnapmin,1)  ! global rmin in napoly zone

      if (nid.eq.0)  write(6,10) nwnpel(1),rnapmax,rnapmin 

      do i= 1,napmax(1) 
         znp(1,i)= zu0(npindx(1,i))
         rnp(1,i)=(rnapmin-rnapmax)/(wzone1-wnapy1)
     $           *(znp(1,i)-wnapy1)+rnapmax
      enddo

      rnapmax = -1.e+20
      rnapmin =  1.e+20
      do iel = 1,nwnpel(2)
         e = wnpel(2,iel)
         do i=1,nxyz
            xx= xm1(i,1,1,e)
            yy= ym1(i,1,1,e)
            rr = xx**2+yy**2
            if (rr.gt.0) rr=sqrt(rr)
            rnapmax=max(rr,rnapmax)
            rnapmin=min(rr,rnapmin)
         enddo
      enddo

      rnapmax=glmax(rnapmax,1)  ! this is equal to the one obtained above
      rnapmin=glmin(rnapmin,1)

      if (nid.eq.0)  write(6,20) nwnpel(2),rnapmin,rnapmax 

      do i= 1,napmax(2)  
         znp(2,i)= zu0(npindx(2,i))
         rnp(2,i)=(rnapmax-rnapmin)/(wnapy2-wzone2)
     $           *(znp(2,i)-wzone2)+rnapmin
      enddo

      if (nid.eq.0)  write(6,*)          

  10  format(' elts:rmax:rmin in napy zone1=',i8,2e18.10)          
  20  format(' elts:rmin:rmax in napy zone2=',i8,2e18.10)          

      return
      end
c-----------------------------------------------------------------------
c parallel
      subroutine interpolation_init_glo                                
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

c     FIXME : implementation needed for parallel  4/24/07  
c     
c     +---------------+---------------+ .......   +----------------+
c     z_1            z_2                                       z_nzarry0 
c
c     +--|----|----|--+--|----|----|--+ .......   +----------------+
c        <-dz->
c        <--------->
c         nuniform


      dzz = (zw0(2)-zw0(1))            
      dzt = (zw0(2)-zw0(1))/cspeed
      

      return
      end
c-----------------------------------------------------------------------
      subroutine interpolation_init                                     
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

      dz = 2.0                                                         
      dz = dz/nuniform       
      dz2= 0.5*dz

      do i=1,nuniform
         zu(i)=-1.0+dz*(i-1)  
      enddo

      l = 0
      do k= 1,nzarry0-1
      do i= 1,nuniform
         dd2= zarry0(k+1)-zarry0(k)
         dd1= 2.0
         l  = i+nuniform*(k-1)
         zu0(l)= zarry0(k)+dd2/dd1*(zu(i)+1)  
      enddo
      enddo
      dzt = (zu0(2)-zu0(1))/cspeed

      return
      end   
c-----------------------------------------------------------------------
c  parallel
      subroutine parallel_zframe_index            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BEAM'
      include 'NEKUSE'
      include 'WAKEPOT'
      
      const   = param(72)                
      smax    = const*zsigma                   
      slast   = zposition+smax          ! s frame max at beginning  
      sfrst   = zposition-smax          ! s frame min at beginning
      scenter = sfrst+0.5*(slast-sfrst) 
      tol     = 1.e-12

      ! [sindx ]= index array for the first frame s
      ! [sgrid ]= set a z array for the first frame s    
      ! [sbnch ]= bunch at the first frame s
      !  nsindx = length(sindx)

      l=0
      nsindx= 0
      do i = 1,nzw0                          
         zz= zw0(i)   
         zznew= -(cspeed*time+zposition-zz)                    
         if ((zz.ge.sfrst-tol).and.(zz.le.slast+tol)) then
              l=l+1
              sindx(l)= i     ! initial s-frame global index         
              sgrid(l)= zznew ! fixed s values for s-frame           
              sbnch(l)= bnch_zdensity(zz) ! chg in s-frame  
              nsindx  = l     ! total number of s-frame index 
         endif
      enddo

      !if (nid.eq.0) write(6,*) 'sindx',node,nsindx,slast        
      !call exitt
      !do i=1,nsindx
      !write(6,*) 'sindx',node,i,sindx(i)
      !enddo
      !call exitt

      if (nsindx.eq.0) then
      ! do nothing
      else
        call rzero(wakez,nsindx)
        call rzero(waker,nsindx)
      endif

      ! [dindx]= index within wakepoten calculation zone
      ! wnapy1,wzone1,wzone2,wnapy2: globally defined                              
      call rzero( dindx,nzw0) ! wake zone
      call rzero(pindx1,nzw0) ! napoly zone 1
      call rzero(pindx2,nzw0) ! napoly zone 2

      do i = 1,nzw0                      
         zz= zw0(i)
         if ((zz.ge.wnapy1-tol).and.(zz.le.wnapy2+tol)) then
         ! zz never be equal to wzone1 or wzone2
              dindx(i)= 1.0
         ! if (nid.eq.0) write(6,*) 'tt',sindx(1),sindx(nsindx),i      
         endif
      enddo

      k1=0
      k2=0
      napmax(1)=0
      napmax(2)=0

      do i = 1,nzw0                      
         zz= zw0(i)
         if ((zz.ge.wnapy1).and.(zz.lt.wzone1)) then
            ! zz never be equal to wnapy1 or wzone1
              pindx1(i)=1.0
              k1=k1+1
              npindx(1,k1)=i  ! global index
              napmax(1) = k1  ! max index              
         endif
         if ((zz.ge.wzone2).and.(zz.lt.wnapy2)) then
            ! zz never be equal to wnapy2 or wzone2
              pindx2(i)=1.0
              k2=k2+1
              npindx(2,k2)=i  ! global index      
              napmax(2) = k2    ! max index
         endif
      enddo

      !self-checking
      if (.not.(napmax(1).eq.nuniform)) then
           write(6,*) 'npg1 is not equal nuniform',wzone1,wzone2
           write(6,*) 'npg1 is not equal nuniform',wnapy1,wnapy2
           write(6,*) 'npg1 is not equal nuniform',npg1,nuniform
           call exitt
      endif

      if (.not.(napmax(2).eq.nuniform)) then       
           write(6,*) 'npg2 is not equal nuniform',npg2,nuniform
           call exitt
      endif

c     timechk=zw0(sindx(1))+cspeed*dzt*nsteps
c     if (timechk.lt.wnapy2) then
c        if (nid.eq.0) write(6,*) 'timechk:wnapy2',timechk,wnapy2
c        if (nid.eq.0) write(6,*) 'ERR: increase nsteps',nsteps
c        call exitt
c     endif

      do nsteps0=1,1000000
         timechk=zw0(sindx(1))+cspeed*dzt/mstep*nsteps0*100
         if (timechk.gt.wnapy2) then
             go to 10
         endif
      enddo

  10  nsteps=nsteps0*100
      do i=1,nsteps
         timechk=zw0(sindx(1))+cspeed*dzt/mstep*i
         if (timechk.gt.wnapy2) then
            nsteps=i  !new nsteps 
            !if (nid.eq.0) write(6,*) 'Maxtime Steps=',i                   
            return
         endif
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine wake_longitudinal            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BEAM'
      include 'NEKUSE'
      include 'WAKEPOT'
      
      integer    wnum,wnum1
      parameter (wnum = 331)
      parameter (wnum1= 332)

      integer    e
      integer    icalld
      save       icalld
      data       icalld /0/
      integer    nsubstep
      save       nsubstep
      data       nsubstep /0/

      character*20 filename
      character*20 filename1
      data         filename /'./wake/wk-t00000.dat'/
      data         filename1/'./wake/wk-r00000.dat'/
      CHARACTER*1  fnarray (20)
      CHARACTER*1  fnarray1(20)
      equivalence  (filename ,fnarray )
      equivalence  (filename1,fnarray1)

      real tim1,tim2,tim3,tim0
      save tim1,tim2,tim3,tim0
      data tim1,tim2,tim3,tim0 /4*0./

      real stim1,stim2,stim3,stim0
      save stim1,stim2,stim3,stim0
      data stim1,stim2,stim3,stim0 /4*0./

      integer    dumpno
      save       dumpno
      data       dumpno /0/

      real       zend1, zend2
      save       zend1, zend2

      character*1  numrl(0:9)
      data numrl /'0','1','2','3','4','5','6','7','8','9'/

      fnarray(16)= numrl(mod(dumpno,10))
      fnarray(15)= numrl(mod(dumpno,100)/10)
      fnarray(14)= numrl(mod(dumpno,1000)/100)
      fnarray(13)= numrl(mod(dumpno,10000)/1000)
      fnarray(12)= numrl(mod(dumpno,100000)/10000)

      fnarray1(16)= numrl(mod(dumpno,10))
      fnarray1(15)= numrl(mod(dumpno,100)/10)
      fnarray1(14)= numrl(mod(dumpno,1000)/100)
      fnarray1(13)= numrl(mod(dumpno,10000)/1000)
      fnarray1(12)= numrl(mod(dumpno,100000)/10000)

      nlev   = nuniform*nlevg           
      nwstep = nuniform*mstep

      if (icalld.eq.0) then
          n=nx1*ny1*nz1*nelt
          zend1=  glmin(zm1,n)
          zend2=  glmax(zm1,n)
          icalld = 1
      endif
 

      timechk=zw0(sindx(1))+cspeed*dzt/mstep*istep

      if (timechk.gt.wnapy2) then
c         if (nid.eq.0) write(6,*) 'Bunch has passed Wake Zone !'
C         call exitt
      else
          !if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) then
          !write(6,*) 'maxtime steps=',nsteps
          !endif
      endif

c ... Longitudinal wake potential calculation
      if (mod(istep,mstep).eq.0) then

      do i= 1,nsindx                    
         ind     = sindx(i)+nsubstep        ! nsubstep=istep/mstep: global index - shifting by nsubstep 
         wake1   = ez_wp (ind)* dindx  (ind)! wake mid zone                
         wake2   = ez_wp (ind)* pindx1 (ind)! napoly curve left zone           
         wake3   = ez_wp (ind)* pindx2 (ind)! napoly curve right zone 
         wake    = wake1+wake2+wake3      
         wake    = wake1                  
         wakez(i)= wakez(i)+dzz*wake*1./totchg 
      enddo


c ... self-checking 
      if (nid.eq.0) then
          write(6,20) istep,zsigma,zend1,zend2,wnapy1,wnapy2,
     $    zw0(sindx(1)+nsubstep),zw0(sindx(nsindx)+nsubstep)
      endif

      if ((mod(istep,nwstep).eq.0).and.(nid.eq.0)) then
          open  (wnum, file=filename)
          do i= 1,nsindx                    
             write(wnum,*) sgrid(i),sbnch(i),wakez(i)
          enddo
          close (wnum)
      endif

c ... IFWAKET=1: transverse wake potential calculation
      if (IFWAKET) then     

         do i= 1,nsindx                    
          ind     = sindx(i)+nsubstep        ! index shifting 
          waker1  = ez_wp1(ind)*dindx (ind)  ! wake mid zone                  
          waker2  = ez_wp1(ind)*pindx1(ind)  ! napoly curve left zone     
          waker3  = ez_wp1(ind)*pindx2(ind)  ! napoly curve right zone   
          wake    = waker1+waker2+waker3
          wakez1(i)= wakez1(i)+dzz*wake      ! write(6,*) 'ww',istep,dz,i,wakez(i)
          waker (i)=(wakez(i)-wakez1(i))/drt !
         enddo

c ... write wake transverse         
        if ((mod(istep,nwstep).eq.0).and.(nid.eq.0)) then
          open  (wnum1, file=filename1)
          do i= 1,nsindx                    
             write(wnum1,*) sgrid(i),sbnch(i),waker(i)
          enddo
          close (wnum1)
        endif

      endif  ! wake transverse

      nsubstep= nsubstep +1

      endif  ! mod(istep,mstep)==0

    5 format(1p3e19.12)       
   10 format(/,' no contribution for wakepoten',i5)       
   20 format(/,' calculating wakepoten: istep=',i7,'; zsigma=',f9.4,
     $ '; zframe=', '[',f9.4,',',f9.4,']','; wake zone=',
     $ '[',f9.4,',',f9.4,']','; s=',
     $ '[',f9.4,',',f9.4,']')
   30 format(' index error wakepotet:istep= ',i5) 
   40 format(3e18.10)

      dumpno  = dumpno + 1

      return
      end
c-----------------------------------------------------------------------
      subroutine extract_integral_path_glo
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BEAM'
      include 'NEKUSE'
      include 'WAKEPOT'

      !integer n_wp_surf,indx_wp_srf(1),e_wp_surf(1),f_wp_surf(1)

      rcirc = wakerad1
      nzlev = nz1*nlevg       

      call extract_cylindrical_surface
     $         (n_wp_surf,e_wp_surf,f_wp_surf,rcirc)
      call extract_wp_surface_indx
     $         (n_wp_grid,n_wp_surf,e_wp_surf,f_wp_surf)
      call set_perimeter_integral_weights
     $     (lwght,nzlev,indxw,izlvi,swght,n_wp_grid,wwrk,zind,node)

      rcirc1 = wakerad2                       
      if (IFWAKET) then
          drt = abs(rcirc1-rcirc)
          call extract_cylindrical_surface
     $         (n_wp_surf1,e_wp_surf1,f_wp_surf1,rcirc1)
          call extract_wp_surface_indx
     $         (n_wp_grid1,n_wp_surf1,e_wp_surf1,f_wp_surf1)
          call set_perimeter_integral_weights
     $    (lwght1,nzlev,indxw1,izlvi1,swght1,n_wp_grid1,wwrk,zind,node)
      endif

      return
      end
c-----------------------------------------------------------------------
c parallel checked
      subroutine extract_cylindrical_surface(n_surf,e_surf,f_surf,rcirc)
c     Find elments and faces for the specified surfaces ar r=rcirc
c     integer  e_surf (n_surf)
c     integer  f_surf (n_surf)
c
      include 'SIZE'
      include 'TOTAL'

      integer n_surf,e_surf(1),f_surf(1)
      integer e,f

c     Initially, we support circles
      rcirc2= rcirc*rcirc
      eps2  = rcirc2*1.e-3

      !write(6,*) 'rcirc',rcirc2,eps2
      !call exitt
 
      mode=0
c ... TEMP
      rad1  = param(75)**2
      rad2  = param(76)**2
      if (mode.eq.1) then
      if (abs(rad1).le.1e-12) then
          if (nid.eq.0) write(6,*) 'ERR: wakerad param (75)',param(75)
          call exitt
      endif
      if (abs(rad2).le.1e-12) then
          if (nid.eq.0) write(6,*) 'ERR: wakerad param (76)',param(76)
          call exitt
      endif
      endif

c ... TEMP

      n_surf = 0

      nfaces=2*ndim
      do e=1,nelt
      do f=1,nfaces
         call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
         kxm = (kx1+kx2)/2
         kym = (ky1+ky2)/2
         kzm = (kz1+kz2)/2
         xmd = xm1(kxm,kym,kzm,e)             ! Face midpoint
         ymd = ym1(kxm,kym,kzm,e)
         rd2 = xmd*xmd + ymd*ymd

         if (abs(rd2-rcirc2) .le. eps2) then  ! face potentially on circle.
         !if (nid.eq.0) write(6,*) 'ee ff',e,f,sqrt(rd2),rd2-rcirc2,eps2
c        if ((rd2.le.rad2).and.(rd2.ge.rad1)) then  ! face potentially on circle.
            !FIXME this doesn't always work                    
            call facind2 (i0,i1,iskip,j0,j1,jskip,f)

            l = 0                  !  Check normals
            dmx = 0
            do j=j0,j1,jskip
            do i=i0,i1,iskip
               l = l+1
               xmd = xm1(i,j,1,e)
               ymd = ym1(i,j,1,e)
               rad = xmd*xmd+ymd*ymd
               if (rad.gt.0) rad = sqrt(rad)
               xnm = xmd/rad
               ynm = ymd/rad
               dot = xnm*unx(l,1,f,e) + ynm*uny(l,1,f,e)
               dm1 = dot - 1.
               dmx = max(dmx,abs(dml))
            enddo
            enddo
            if (dmx.le.1.e-3) then ! found face on inside surf of rcirc
                n_surf = n_surf+1
                e_surf(n_surf) = e
                f_surf(n_surf) = f
               !write(6,*) 'surf',n_surf,e_surf(n_surf),f_surf(n_surf)
            endif

         endif

      enddo
      enddo

c     if (n_surf.eq.0) then
c        if (nid.eq.0) write(6,*) 
c    $      'ERR: integral path not recognized, surf=',n_surf
c        call exitt
c     endif
 
      return
      end
c-----------------------------------------------------------------------
c parallel checked
      subroutine extract_wp_surface_indx(n_indx,n_surf,e_surf,f_surf)

c     Find point index and  weights for the specified surfaces ar rcirc
c          indxw (n_indx) point index
c          swght (n_indx) weights            
c          izlvi (n_indx) indicate level at z on grids on the surface

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

      integer n_indx,n_surf,e_surf(1),f_surf(1)
      integer e,f,ef
c
      nxz  = nx1*nz1
      nxzf = nx1*nz1*2*ndim
      nxyz = nx1*ny1*nz1

      npts = 0
      do k = 1,n_surf

         e = e_surf(k)
         f = f_surf(k)
         ef= eface (f)     !
         lev = izelev (e)  ! global z-level for element e

         call facind2 (i0,i1,iskip,j0,j1,jskip,f)

         l = 0
         klev = nz1*(lev-1)

         do j = j0,j1,jskip     ! This is not the local z index
            klev= klev + 1      ! Level for this slice
            do i= i0,i1,iskip   ! March around "perimeter" at a fixed z
               l= l+1
               npts= npts+1
               indxw (npts)= i  + nx1*(j -1) + nxyz*(e-1)
               swght (npts)= area(l,1,f,e) ! FIXME : fixed 10/31/07 mmin            
               izlvi (npts)= klev           ! level at z grid Between 1 and nz1*nlevg     
               n_indx      = npts
            enddo

         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_perimeter_integral_weights
     $           (lwght,nlev,indx,izlvi,swght,npts,w,zind,node)
c     w() is a work array of length nlev = nz1*(nzarry0-1)
c     lwght(nlev)                                                   
 
      integer nlev
      integer indx(npts),izlvi(npts)
      real    p_int(nlev),lwght(nlev),zind(nlev),swght(npts),u(1),w(1)

      call rzero(lwght,nlev)  ! nlev= nz1*nlevg: global size
      do i=1,npts
         klev= izlvi(i)       ! global index of z slice: zni  
         lwght(klev)= lwght(klev)+swght(i)
         !if (klev.eq.190) write(6,*)  'klev-1-',klev,swght(i)
         !if (klev.eq.200) write(66,*) 'klev-1-',klev,swght(i)
      enddo

      ! FIXME nelv is diferent at each node
      call gop (lwght,w,'+  ',nlev)   ! sum across all processors
      !call invcol2(lwght,zind,nlev)  ! divide by redundant number
      call invcol1(lwght,nlev)

      !checked for np=1,..,np=64
      !ggg=glmax(lwght,nlev)
      !write(6,*) 'ggg--',npts,ggg                           
      !call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine perimeter_integral  ! PERIMETER INTEGRAL OF U
     $           (p_int,lwght,nlev,indx,izlvi,swght,npts,u,w,zind)
 
 
c     w() is a work array of length nlev = nz1*(nzarry0-1)

      integer indx (npts),izlvi(npts)
      real    p_int(nlev),lwght(nlev),zind(nlev),swght(npts),u(1),w(1)
 
      call rzero(p_int,nlev)
      do i=1,npts
         klev = izlvi(i)   ! global info for zslice
         p_int(klev) = p_int(klev) + u(indx(i))*swght(i)
         ! write(6,*) 't',i,indx(i),swght(i),p_int(klev)
         ! write(6,*) 't',i,indx(i),swght(i),klev
      enddo
      !call exitt
 
      call col2(p_int, lwght ,nlev)
      call gop (p_int,w,'+  ',nlev)    ! sum across all processors
      !call invcol2(p_int,zind,nlev)   ! divide by redundant number
 
      return
      end
c-----------------------------------------------------------------------
      subroutine interpolation_path_gll_to_s(uout,uin)           
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

      integer  e
      real     uin  (lz1,levg)
      real     uout (luniform,levg)

      nxy  = 1               
   
      call zwgll  (zi,wi,nz1)
      call igllm  (ju,jut,zi,zw,nz1,nuniform,nz1,nuniform)
      ! ju (k,i)= l_i(x_k) interpolation at x_k
      ! jut(i,k)= l_i(x_k) uout=[uin][jut]              
    
      do e=1,nlevg 
         call mxm (uin(1,e),nxy,jut,nz1,uout(1,e),nuniform)
      enddo
     
      return
      end
c-----------------------------------------------------------------------
      subroutine total_charge           
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

      n = nx1*ny1*nz1*nelt
      do i=1,n
         xx=xm1(i,1,1,1) 
         yy=ym1(i,1,1,1) 
         zz=zm1(i,1,1,1) 
         rr=xx**2+yy**2
         if (rr.gt.0) rr=sqrt(rr)
         !zznew   = -(cspeed*time+zposition-zz)                    
         rdensity= bnch_rdensity(rr)
         zdensity= bnch_zdensity(zz)
         workfld(i,1,1,1)= rdensity*zdensity
      enddo
c     write(6,*) 'vol',volvm1
      totchg = glsc2(workfld,bm1,n)/volvm1
c     write(6,*) 'total_charge',totchg
  
      return
      end
c-----------------------------------------------------------------------
