
c-----------------------------------------------------------------------
c
c  USER SPECIFIED ROUTINES:
c
c     - boundary conditions
c     - initial conditions
c     - variable properties
c     - forcing function for fluid (f)
c     - forcing function for passive scalar (q)
c     - general purpose routine for checking errors etc.
c
c-----------------------------------------------------------------------
      subroutine usersol(tt, myshx, myshy, myshz, mysex, mysey, mysez)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'

      real tt
      real myshx(lx1,ly1,lz1,lelt)
      real myshy(lx1,ly1,lz1,lelt)
      real myshz(lx1,ly1,lz1,lelt)
      real mysex(lx1,ly1,lz1,lelt)
      real mysey(lx1,ly1,lz1,lelt)
      real mysez(lx1,ly1,lz1,lelt)
      real xmaxsize
      save xmaxsize

      do  i= 1,nx1*ny1*nz1*nelt
         xx= xm1(i,1,1,1)
         yy= ym1(i,1,1,1)
         zz= zm1(i,1,1,1)

       ! myshx(i,1,1,1) = 0
       ! myshy(i,1,1,1) = 0
       ! myshz(i,1,1,1) = 0
       ! mysex(i,1,1,1) = 0
       ! mysey(i,1,1,1) = 0
       ! mysez(i,1,1,1) = 0
      enddo

      return
      end
c-------------------------------------------------------------------
      subroutine userinc

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'
      include 'SCRATCH'
      include 'RK5'

      ! Add source terms to the fields given in src{e,h}{x,y,z}.
      ! Read file doc/rea-file.txt for specification.
      ! Since we want to make it possible to have several sources,
      ! we specify a `baseidx', relative to which all parameters are
      ! read.

      ! Extent of domain
      real    kx , ky , kz 
      real    impedence_freespace                                     
      real    impedence                                               

c ... presetting
      twopi  = 2.0*pi
      impedence_freespace=sqrt(permit_0/permea_0)

c ... spatial mode set  
      wavelength = param(69)*nmscale
      refractive_index= 1.0 

      kx    = param(66)                  ! smodty=3: planewave
      ky    = param(67)                  
      kz    = param(68)              

      kx = refractive_index*kx
      ky = refractive_index*ky
      kz = refractive_index*kz

      tmp= permea_0/(permit_0*refractive_index**2)
      impedence=1./sqrt(tmp)

      n= nx1*ny1*nz1*nelt
 
      do i= 1,n   
         xx = xm1(i,1,1,1)
         yy = ym1(i,1,1,1)
         zz = zm1(i,1,1,1)

         if     (incflag.eq.1) then
           if ((cspeed*rktime-zz).ge.zminval) then
              tfac=outpml(i,1,1,1)
           else
              tfac=0.0
           endif
         elseif (incflag.eq.-1) then
           if ((cspeed*rktime+zz).ge.zmaxval) then
              tfac=outpml(i,1,1,1)
           else
              tfac=0.0
           endif
         endif

         incey(i,1,1,1)=  0             
         incex(i,1,1,1)= tfac*                                         
     $     cos(twopi/wavelength*(cspeed*rktime-(kx*xx+ky*yy+kz*zz)))
         incez(i,1,1,1)=  0             
         inchx(i,1,1,1)=  0             
         inchy(i,1,1,1)= tfac* impedence*
     $     cos(twopi/wavelength*(cspeed*rktime-(kx*xx+ky*yy+kz*zz)))
         inchz(i,1,1,1)=  0             
         workvec(i)=incex(i,1,1,1)*inchy(i,1,1,1)
      enddo
      !fff=glmax(workvec,n)
      !ggg=glmin(workvec,n)
      !write(6,*) 'fff---',istep,fff,ggg
  
      ! call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine usersrc(baseidx,srchx,srchy,srchz,srcex,srcey,srcez)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'
      include 'SCRATCH'

      ! Add source terms to the fields given in src{e,h}{x,y,z}.
      ! Read file doc/rea-file.txt for specification.
      ! Since we want to make it possible to have several sources,
      ! we specify a `baseidx', relative to which all parameters are
      ! read.

      integer baseidx

      real srchx(lx1,ly1,lz1,lelt),
     $     srchy(lx1,ly1,lz1,lelt),
     $     srchz(lx1,ly1,lz1,lelt),
     $     srcex(lx1,ly1,lz1,lelt),
     $     srcey(lx1,ly1,lz1,lelt),
     $     srcez(lx1,ly1,lz1,lelt)

      n= nx1*ny1*nz1*nelt

      do i = 1,nx1*ny1*nz1*nelt

         xx= xm1(i,1,1,1)
         yy= ym1(i,1,1,1)
         zz= zm1(i,1,1,1)
         tmp= bm1(i,1,1,1)*outpml(i,1,1,1)

         srcex(i,1,1,1)=srcex(i,1,1,1)+incex(i,1,1,1)*tmp           
         srchy(i,1,1,1)=srchy(i,1,1,1)+inchy(i,1,1,1)*tmp          

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine userini(tt, myhx, myhy, myhz, myex, myey, myez)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real tt
      real myhx(lx1,ly1,lz1,lelt)
      real myhy(lx1,ly1,lz1,lelt)
      real myhz(lx1,ly1,lz1,lelt)
      real myex(lx1,ly1,lz1,lelt)
      real myey(lx1,ly1,lz1,lelt)
      real myez(lx1,ly1,lz1,lelt)

      do i=1,nx1*ny1*nz1*nelt
        xx=XM1(i,1,1,1)
        yy=YM1(i,1,1,1)
        zz=ZM1(i,1,1,1)

        myhx(i,1,1,1) = 0
        myhy(i,1,1,1) = 0
        myhz(i,1,1,1) = 0
        myex(i,1,1,1) = 0
        myey(i,1,1,1) = 0
        myez(i,1,1,1) = 0
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ielg)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy, iz, ielg, iel
      real     pi, one, xcc, ycc, zcc, x1, y1, z1
      real     lambda, fc, wc                     

      iel=gllel(ielg)

      call usermedia   (ix,iy,iz,iel)

      if (IFDRUDE ) call userdrude (ix,iy,iz,iel)

      return
      end

c-----------------------------------------------------------------------

      subroutine usermedia(ix,iy,iz,iel) ! Vector form of userf
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy,iz, ieg
      real     pi, one, xcc, ycc, zcc, x1, y1, z1
      real     lambda, fc, wc                     

      udiff  = 0.
      utrans = 0.

      permea = permea_0             
      permit = permit_0                

      return
      end

c-----------------------------------------------------------------------
      subroutine userdrude(ix,iy,iz,iel) ! Vector form of userf

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy,iz, iel
      real     pi, one, xcc, ycc, zcc, x1, y1, z1
      real     rr, r2, rad, rad2                                      
      real     lambda, fc, wc                     
      real n_subs

      n_subs=1.59 ! index of refraction of substrate from Jason's .ppt


      if     (if_in_region(iel).eq.0) then ! air
             permit      = permit_0
             permea      = permea_0
             drude_alpha = 0.0
             drude_beta  = 0.0

      elseif (if_in_region(iel).eq.1) then ! substrate
             permit_r    = n_subs**2
             permit      = permit_0*permit_r
             permea      = permea_0
             drude_alpha = 0.0
             drude_beta  = 0.0

      elseif (if_in_region(iel).eq.2) then ! gold
            !eV -> SI: conversion from eV to inverse seconds
             tmp =  (27.2114*2.41888e-17)
            !permit_r    = 5.398334980
            !gamma_p     = 0.068017714/tmp
            !freq_p      = 9.200688/tmp
             permit_r    = 1.0000000000
             gamma_p     = 0.23670/tmp
             freq_p      = 7.32318/tmp
             permit      = permit_0*permit_r
             permea      = permea_0
             drude_alpha =-gamma_p
             drude_beta  = permit_0*freq_p**2
            ! write(6,*) 'drude---',iel,drude_alpha,drude_beta
      else
            if (nid.eq.0) write(6,*) 'dielectric or metal:not defiened'
            call exitt
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine userpml(ix,iy,iz,iel) ! Vector form of userf
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      character*5 pmlside   

      n= nx1*ny1*nz1*nelt

      pmlR0 = 10.0**(-15)
      pmlM  = 2.0

      return
      end

c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer ix, iy, iz, ieg
C
      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C
      qvol   = 0.0
      source = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine useric(ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C
      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux=0.0
      uy=0.0
      uz=0.0
      temp=0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      return
      end

c-----------------------------------------------------------------------
      subroutine usrdat2

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real pi
      real r,r0,z0,z1,z2,z3,z4,z5,dz

      n    = nx1*ny1*nz1*nelv

c     ifxyo= .true.

      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)
      zmin = glmin(zm1,n)
      zmax = glmax(zm1,n)

      if (IFNM) then
          nmscale = 1.e-9
          sx = nmscale            
          sy = nmscale                  
          sz = nmscale
          if (IF3D) sz = nmscale
      else
          nmscale = 1.0
          sx = nmscale
          sy = nmscale
          if (IF3D) sz = nmscale
      endif

      if (if3d) then

      do i=1,n
         xm1(i,1,1,1) = sx*xm1(i,1,1,1)
         ym1(i,1,1,1) = sy*ym1(i,1,1,1)
         zm1(i,1,1,1) = sz*zm1(i,1,1,1)
      enddo

      else

      do i=1,n
         xm1(i,1,1,1) = sx*xm1(i,1,1,1)       
         ym1(i,1,1,1) = sy*ym1(i,1,1,1)                
      enddo

      endif

      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)
      zmin = glmin(zm1,n)
      zmax = glmax(zm1,n)

      r0=257*nmscale ! EDITME
      dz=36*nmscale ! EDITME
      z0=0
      z1=944*nmscale ! EDITME 
      z2=z1+dz
      z3=1304*nmscale ! EDITME 
      z4=z3+dz
      z5=1780*nmscale ! EDITME 

      nxy2 = (nx1+1)/2

      do ie=1,nelv
         xmid = xm1(nxy2,nxy2,1,ie)
         ymid = ym1(nxy2,nxy2,1,ie)
         zmid = zm1(nxy2,nxy2,1,ie)
         rr = sqrt(xmid**2+ymid**2)
         
         if (zmid.le.z1) then
           if_in_region(ie)=1 ! substrate
         elseif (zmid.le.z2) then
           if (rr.le.r0) then
             if_in_region(ie)=2 ! gold
           else
             if_in_region(ie)=1 ! substrate
           endif
         elseif (zmid.le.z3) then
           if (rr.le.r0) then
             if_in_region(ie)=0 ! air
           else
             if_in_region(ie)=1 ! substrate
           endif
         elseif (zmid.le.z4) then
           if (rr.le.r0) then
             if_in_region(ie)=0 ! air
           else
             if_in_region(ie)=2 ! gold
           endif
         elseif (zmid.le.z5) then
           if_in_region(ie)=0 ! air
         else
           print *, 'zmid=', zmid, ' is invalid.'
           call exit
         endif
      enddo

      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)

      !write(6,*) 'xmin/xmax--',xmin,xmax,ymin,ymax
      !call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine surface_integration                 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer icalld
      save    icalld
      data    icalld  /0/

      if (icalld.eq.0) then
          call cem_set_fc_ptr
          icalld = 1
      endif

      n=nx1*ny1*nz1*nelt

      do i=1,n
         workfld(i,1,1,1)=ex(i,1,1,1)**2+hx(i,1,1,1)**2
      enddo

      srfint=0.0
      do j=1,ncemface
         i=cemface(j)
         faceval= 0.1*faceid(i,1,1,1)*workfld(i,1,1,1)
         srfint = srfint+faceval*area(j,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      common /ccpu/ cpu_t,cpu_t0

c ... compute errors

      if (istep.le.10.or.mod(istep,iocomm).eq.0) then

        ntott = nx1*ny1*nz1*nelt
        nsrf  = nx1*nz1*nelt*6

        if (ifsol) then
          call usersol(time,shx,shy,shz,sex,sey,sez)
        else
          call rzero(shx, ntott) 
          call rzero(shy, ntott) 
          call rzero(shz, ntott) 
          call rzero(sex, ntott) 
          call rzero(sey, ntott) 
          call rzero(sez, ntott) 
        endif

        call sub3(ErrHX,SHX,HX,NTOTT) 
        call sub3(ErrHY,SHY,HY,NTOTT) 
        call sub3(ErrHZ,SHZ,HZ,NTOTT) 
        call sub3(ErrEX,SEX,EX,NTOTT) 
        call sub3(ErrEY,SEY,EY,NTOTT) 
        call sub3(ErrEZ,SEZ,EZ,NTOTT) 

        HXmr = glamax(ErrHX,NTOTT) 
        HYmr = glamax(ErrHY,NTOTT) 
        HZmr = glamax(ErrHZ,NTOTT) 
        EXmr = glamax(ErrEX,NTOTT) 
        EYmr = glamax(ErrEY,NTOTT) 
        EZmr = glamax(ErrEZ,NTOTT) 

        HXer = glsc3(ErrHX,bm1,ErrHX,NTOTT) / volvm1
        HYer = glsc3(ErrHY,bm1,ErrHY,NTOTT) / volvm1
        HZer = glsc3(ErrHZ,bm1,ErrHZ,NTOTT) / volvm1
        EXer = glsc3(ErrEX,bm1,ErrEX,NTOTT) / volvm1
        EYer = glsc3(ErrEY,bm1,ErrEY,NTOTT) / volvm1
        EZer = glsc3(ErrEZ,bm1,ErrEZ,NTOTT) / volvm1

        if (hxer.gt.0) hxer = sqrt(hxer)
        if (hyer.gt.0) hyer = sqrt(hyer)
        if (hzer.gt.0) hzer = sqrt(hzer)
        if (exer.gt.0) exer = sqrt(exer)
        if (eyer.gt.0) eyer = sqrt(eyer)
        if (ezer.gt.0) ezer = sqrt(ezer)


c ... compute max gap size at interfaces between elmts.

        call CEM_GAP                     

        HXgp = glamax(gHX,NSRF) 
        HYgp = glamax(gHY,NSRF) 
        HZgp = glamax(gHZ,NSRF) 
        EXgp = glamax(gEX,NSRF) 
        EYgp = glamax(gEY,NSRF) 
        EZgp = glamax(gEZ,NSRF) 

c ... print statement

        snel = nelt + .5
        snel = sqrt(snel)
        nels = snel
        nnn  = nx1-1
        npt1 = nnn*nels
        npts = nx1*ny1*nz1*nelt ! npts = npt1**ndim

        if (istep.eq.0)       cpu_p_t= 0
        if (istep.gt.0)       cpu_p_t= cpu_t/(npts*istep)
        if (istep.eq.nsteps)  cpu_p_t= glsum(cpu_p_t,1)/np

        if (nid.eq.0) then

          write(6,91) ISTEP,nelt,nels,nnn,npt1,npts
     $         ,time,dt,HXer,HYer,HZer,EXer,EYer,EZer,cpu_t,cpu_p_t
          write(6,95) ISTEP,nelt,nels,nnn,npt1,npts
     $         ,time,dt,HXmr,HYmr,HZmr,EXmr,EYmr,EZmr,cpu_t,cpu_p_t

 91       format(/,i9,3i5,i5,i9,1p9e10.3,e9.2,' CPU: L2  '  )
 95       format(  i9,3i5,i5,i9,1p9e10.3,e9.2,' CPU: Linf'  )

        endif ! nid=0

      endif

      call process_poynting

      return
      end

c-----------------------------------------------------------------------
      subroutine process_poynting
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer icalld
      save    icalld
      data    icalld /0/

      integer e,f,nkey,iisum
      save    nkey,iisum
      common /cpoynt/ flag(2*ldim,lelt)
      integer         flag

      if (icalld.eq.0) then  ! customize your selection here

         do e=1,nelv
         do f=1,2*ndim
            flag(f,e) = 0
         enddo
         enddo

c        user supplies z0 and z1 that bracket the surface of interest

         n=nx1*ny1*nz1*nelt
         xmin = glmin(xm1,n)
         xmax = glmax(xm1,n)
         ymin = glmin(ym1,n)
         ymax = glmax(ym1,n)

c        Find faces #6 and with z0 < z < z1:
         if     (incflag.eq.1) then
              tsurface = 1580*nmscale
         elseif (incflag.eq.-1) then
              tsurface =  200*nmscale
         endif

         zz_0= tsurface -1.0e-14
         zz_1= tsurface +1.0e-14
         write(6,*) 'checking  =',tsurface,zz_0,zz_1 

         nkey=1
         isum = 0
         f = 6          ! assume it's an upward pointing z face
         do e=1,nelv
            call fcsum2(zsum,asum,zm1,e,f)  ! f in prex notation
            if (asum.gt.0) zavg = zsum/asum
            if (zz_0.le.zavg.and.zavg.le.zz_1) then
               flag(f,e) = nkey
               isum = isum+1
         !write(6,*) 'checking face numbers =',xmax,zavg,f,e
            endif
         enddo
         iisum=iglsum(isum,1)
         if (nid.eq.0) write(6,*)
     $       'transmission surface face numbers =',iisum

         icalld = 1

      endif

      if (iisum.eq.0) then
         if (nid.eq.0)  write(6,*)
     $      'integration surface:not correct',nkey,iisum
          call exitt
      endif
      !poynting vectors (avgSz,avgISz) computed: cem_energy 
        
      nkey = 1
      n = nx1*ny1*nz1*nelt

      ubar1 = 0
      ubar2 = 0
      ubar3 = 0
      ubar4 = 0
      atot1 = 0
      atot2 = 0
      atot3 = 0
      atot4 = 0

      nfaces= 2*ndim

      if (ifscat) then
          call sub3(totex,ex,incex,n)
          call sub3(totey,ey,incey,n)
          call sub3(totez,ez,incez,n)
          call sub3(tothx,hx,inchx,n)
          call sub3(tothy,hy,inchy,n)
          call sub3(tothz,hz,inchz,n)

      call vcross (poytSx,poytSy,poytSz,
     $             totex,totey,totez,tothx,tothy,tothz,n)
      call vcross (poySx,poySy,poySz,ex,ey,ez,hx,hy,hz,n)
      call vcross (poyISx,poyISy,poyISz,
     $             incex,incey,incez,inchx,inchy,inchz,n)
      elseif (ifsftf) then
      call vcross (poySx,poySy,poySz,ex,ey,ez,hx,hy,hz,n)
      call vcross (poyISx,poyISy,poyISz,
     $             incex,incey,incez,inchx,inchy,inchz,n)
      endif

      do e = 1,nelv
      do f = 1,nfaces
         if (flag(f,e).eq.nkey) then ! add this to this average
            call fcsum2(usum1,asum1,poytSz,e,f)  ! f in prex notation
            call fcsum2(usum2,asum2,poySz ,e,f)  ! f in prex notation
            call fcsum2(usum3,asum3,poyISz,e,f)  ! f in prex notation
            call fcsum2(usum4,asum4,workvec,e,f)  ! f in prex notation
            ubar1 = ubar1 + usum1
            ubar2 = ubar2 + usum2
            ubar3 = ubar3 + usum3
            ubar4 = ubar4 + usu4m
            atot1 = atot1 + asum1
            atot2 = atot2 + asum2
            atot3 = atot3 + asum3
            atot4 = atot4 + asum4
            isum  = isum  + 1
         endif
      enddo
      enddo

      ubar1= glsum(ubar1,1)
      ubar2= glsum(ubar2,1)
      ubar3= glsum(ubar3,1)
      ubar4= glsum(ubar4,1)

      atot1= glsum(atot1,1)
      atot2= glsum(atot2,1)
      atot3= glsum(atot3,1)
      atot4= glsum(atot4,1)

      if (istep.ge.1) then
          if (nid.eq.0) write(33,1) nperiod,istep,time,dt,ubar1,ubar2,
     $        ubar3,cspeed,param(69)*nmscale,nopstep,n_period,nx1-1
c         if (nid.eq.0) write(*,1) nperiod,istep,time,dt,ubar1,ubar2,
c    $        ubar3,cspeed,param(69)*nmscale,nopstep,n_period,nx1-1
      endif

      !if (istep.le.10.or.mod(istep,iocomm).eq.0) then 
      !   fff=glmax(workvec,n)
      !   write(6,*) 'ubar3', istep,ubar3,ubar4,fff
      !endif

      call trans_calculation(ubar1,ubar2,ubar3)

    1 format(2i9,1p7e25.15,3i5)

      return
      end

c-----------------------------------------------------------------------
      subroutine trans_calculation(ubar1,ubar2,ubar3)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      real     ubar1,ubar2,ubar3

      if (nsteps.eq.nopstep) then
          if (nid.eq.0) then
             write(6,*) 'one period? increase # of periods param(14)'
             call exitt
          endif
      endif

      if ((istep.ge.(n_period-1)*nopstep+1).and.
     $    (istep.le. n_period*nopstep)) then

       if (istep.eq.(n_period-1)*nopstep+1) then
          trans_tot =  0.0
          trans_scat=  0.0
          trans_inc =  0.0
       endif
       tmp = dt/cspeed/(param(69)*nmscale)
       trans_tot =trans_tot +ubar1*tmp
       trans_scat=trans_scat+ubar2*tmp 
       trans_inc =trans_inc +ubar3*tmp

       if (istep.eq.n_period*nopstep) then
        if ((nsteps.eq.nopstep).or.(istep.eq.nsteps)) then
          if (trans_inc.le.1.0e-14) then
             if (nid.eq.0)
     $       write(6,*) 'trans_inc=0; increase # of periods param(14)'
             call exitt
          endif
        endif
        if (nid.eq.0) write(6,*)  ' '
        if (nid.eq.0) then
           write(6,*) '==============================================='
           write(6,*)  n_period,'-th period:'
           write(6,*) '-----------------------------------------------'
           write(6,*) 'period=',(param(69)*nmscale)/cspeed
           write(6,*) 'finaltime=',dt*nsteps
           if (ifnm) then
              write(6,*) 'wavelength(nm)=',param(69)
           else
              write(6,*) 'wavelength= ',param(69)
           endif
           write(6,*) '-----------------------------------------------'
           transmit_t=trans_tot /trans_inc
           transmit_s=trans_scat/trans_inc
           write(6,*) 'Transmission_Tot =',transmit_t
           write(6,*) 'Transmission_Scat=',transmit_s
           write(6,*)  '=============================================='
        endif
       endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gl_fcavg(ubar,atot,u,flag,nkey,nx,nelv,ndim)
      real     u(nx,nx,nx,nelv)
      integer  e,f,nkey
      real     ubar,atot
      integer  flag(2*ndim,nelv)

      ubar = 0
      atot = 0

      nfaces= 2*ndim

      do e=1,nelv
      do f=1,nfaces
         if (flag(f,e).eq.nkey) then ! add this to this average
            call fcsum2(usum,asum,u,e,f)  ! f in prex notation
            ubar = ubar + usum
            atot = atot + asum
         endif
      enddo
      enddo

      ubar= glsum(ubar,1)
      atot= glsum(atot,1)
      if (atot.gt.0) ubar = ubar/atot
      return
      end
c-----------------------------------------------------------------------

