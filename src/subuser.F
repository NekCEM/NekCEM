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

      real tt
      real myshx(lx1,ly1,lz1,lelt)
      real myshy(lx1,ly1,lz1,lelt)
      real myshz(lx1,ly1,lz1,lelt)
      real mysex(lx1,ly1,lz1,lelt)
      real mysey(lx1,ly1,lz1,lelt)
      real mysez(lx1,ly1,lz1,lelt)

      do i=1,nx1*ny1*nz1*nelt
        xx=XM1(i,1,1,1)
        yy=YM1(i,1,1,1)
        zz=ZM1(i,1,1,1)

        myshx(i,1,1,1) = 0
        myshy(i,1,1,1) = 0
        myshz(i,1,1,1) = 0
        mysex(i,1,1,1) = 0
        mysey(i,1,1,1) = 0
        mysez(i,1,1,1) = 0
      enddo

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

      ! Extent of domain
      integer nxyzl
      save    nxyzl

      real    xsize, ysize, zsize
      save    xsize, ysize, zsize
      real    xmin, xmax, ymin, ymax, zmin, zmax
      save    xmin, xmax, ymin, ymax, zmin, zmax

      integer icalld
      save    icalld
      data    icalld /0/             

      real    srcscal
      real    tstart, tstop
      real    tfactor
      integer tmodtyp
      integer smodtyp
      integer fldchs

      real    xx, yy, zz
      real    xloc,yloc,zloc,xw,yw,zw,xwidth,ywidth,zwidth          
      real    xgauss, ygauss, zgauss
      real    kx0, ky0, kz0
      real    kx , ky , kz 
      real    ww, w, twopi
      real    fldhere
      real    wavelength
      real    impedence_freespace                                     

c ... presetting
      twopi  = 2.0*pi
      impedence_freespace=sqrt(permit_0/permea_0)
      nxyzl  = nx1*ny1*nz1*nelt

c ... define the source modes in time and space            

      tmodtyp = int(param(57))
      smodtyp = int(param(58))         
      srcscal =     param(59)
      fldchs  = int(param(60))
      tstart  =     param(62)          
      tstop   =     param(63)             

c ... temporal mode set  
      pulse_loc    = param(64)  ! tmodty=2: gaussian 
      pulse_width  = param(65)

c ... temporal mode set  
      temporal_freq= param(64)  ! tmodty=3: cosine pulse

c ... temporal mode set  
      duration     = param(65)  ! tmodty=4: rectangular pulse

c ... spatial mode set  
      xloc  = param(66)*nmscale         ! smodty=2: gaussian 
      yloc  = param(67)*nmscale  
      zloc  = param(68)*nmscale  
      xw    = param(69)  
      yw    = param(70)  
      zw    = param(71)  
      xwidth= param(72)*nmscale  
      ywidth= param(73)*nmscale  
      zwidth= param(74)*nmscale  

c ... spatial mode set  
      kx    = param(66)                  ! smodty=3: planewave
      ky    = param(67)                  
      kz    = param(68)              

      wavelength = param(69)*nmscale

c ... begin                   

      if ((rktime.lt.tstart).or.(rktime.gt.tstop)) then
        return
      endif

c ... Evaluate temporal modulation.

      if  (tmodtyp.eq.0)  tfactor =  0.0                  
      if  (tmodtyp.eq.1)  tfactor =  1.0
      if  (tmodtyp.eq.2)   
     $       call gussian_pulse(tfactor,rktime,pulse_loc,pulse_width)       
      if  (tmodtyp.eq.3)        
     $       call cosine_pulse (tfactor,rktime,twopi,temporal_freq)
      if  (tmodtyp.eq.4)         
     $       call rectangular_pulse(tfactor,rktime,duration)
      if ((tmodtyp.ge.5).or.(tmodty.lt.0)) then
             write (*,*) 'Invalid time modulation type', tmodtyp
             call exitt
      endif


c ... Evaluate spatial modulation.
      if      (smodtyp.eq.1) then

          do i=1,nxyzl
             fullary(i,1,1,1) = tfactor * srcscal
          enddo

      elseif  (smodtyp.eq.2) then    ! Gaussian

          xgauss = xloc+xw*rktime                                       
          ygauss = yloc+yw*rktime                                       
          zgauss = zloc+zw*rktime                                       

          do i=1,nxyzl

             xx = xm1(i,1,1,1)
             yy = ym1(i,1,1,1)
             zz = zm1(i,1,1,1)

             fldhere = srcscal

             if (xwidth.ne.0)             
     $           fldhere = fldhere*exp(-(2.0*(xx-xgauss)/xwidth)**2)   
             if (ywidth.ne.0) 
     $           fldhere = fldhere*exp(-(2.0*(yy-ygauss)/ywidth)**2)
             if (zwidth.ne.0) 
     $           fldhere = fldhere*exp(-(2.0*(zz-zgauss)/zwidth)**2)

             fullary(i,1,1,1) = tfactor * fldhere

          enddo

      elseif (smodtyp.eq.3) then

          fldhere = tfactor*srcscal

          do i=1,nxyzl

             xx = xm1(i,1,1,1)
             yy = ym1(i,1,1,1)
             zz = zm1(i,1,1,1)

             fullary(i,1,1,1)=fldhere*bm1(i,1,1,1)
     $       *sin(twopi/wavelength*(cspeed*rktime-(kx*xx+ky*yy+kz*zz))) 

             fullaryprime(i,1,1,1)=fldhere*bm1(i,1,1,1)                
     $       *cos(twopi/wavelength*(cspeed*rktime-(kx*xx+ky*yy+kz*zz)))
     $       *twopi/wavelength*cspeed

           enddo

      elseif (smodtyp.eq.4) then
        ! APS source term
        ! empty-> reference:: generic-source.inc
      else
          write (*,*) 'Invalid spatial modulation type', smodtyp
          call exitt
      endif

      ! Add fullary to selected field

      if (IFSCAT) then

          if (icalld.eq.0) then
             do i =1,nxyzl
                conscat1(i,1,1,1)=-(epsilo(i,1,1,1)-permit_0)
                conscat2(i,1,1,1)=-(    mu(i,1,1,1)-permea_0)
     $                            * impedence_freespace
             enddo
             icalld=1        
          endif

         if     (fldchs.eq.1) then

           call addcol3(srchy,fullaryprime,conscat2,nxyzl)
           call addcol3(srchz,fullaryprime,conscat2,nxyzl)
           call addcol3(srcey,fullaryprime,conscat1,nxyzl)
           call addcol3(srcez,fullaryprime,conscat1,nxyzl)

           call copy  (incEY,fullary,nxyzl)
           call copy  (incEZ,fullary,nxyzl)
           call copy  (incHY,fullary,nxyzl)
           call copy  (incHZ,fullary,nxyzl)

         endif

         if     (fldchs.eq.2) then

           call addcol3(srcex,fullaryprime,conscat1,nxyzl)
           call addcol3(srcez,fullaryprime,conscat1,nxyzl)
           call addcol3(srchx,fullaryprime,conscat2,nxyzl)
           call addcol3(srchz,fullaryprime,conscat2,nxyzl)

           call copy (incEX,fullary,nxyzl)
           call copy (incEZ,fullary,nxyzl)
           call copy (incHX,fullary,nxyzl)
           call copy (incHZ,fullary,nxyzl)

         endif

         if     (fldchs.eq.3) then

           call addcol3(srcex,fullaryprime,conscat1,nxyzl)
           call addcol3(srcey,fullaryprime,conscat1,nxyzl)
           call addcol3(srchx,fullaryprime,conscat2,nxyzl)
           call addcol3(srchy,fullaryprime,conscat2,nxyzl)

           call copy (incEX,fullary,nxyzl)
           call copy (incEY,fullary,nxyzl)
           call copy (incHX,fullary,nxyzl)
           call copy (incHY,fullary,nxyzl)

         endif

         if ((fldchs.le.0).or.(fldchs.ge.4)) then 
             write(*,*) 'Invalid value for fldchs in general source.'
             call exitt
         endif

      else
     
         if  (fldchs.eq.1) call fadd2(srchx, fullary, nxyzl)
         if  (fldchs.eq.2) call fadd2(srchy, fullary, nxyzl)
         if  (fldchs.eq.3) call fadd2(srchz, fullary, nxyzl)
         if  (fldchs.eq.4) call fadd2(srcex, fullary, nxyzl)
         if  (fldchs.eq.5) call fadd2(srcey, fullary, nxyzl)
         if  (fldchs.eq.6) call fadd2(srcez, fullary, nxyzl)
         if ((fldchs.le.0).or.(fldchs.ge.7)) then 
             write(*,*) 'Invalid value for fldchs in general source.'
             call exitt
         endif

      endif

      return
      end


      subroutine gussian_pulse(tfactor,tt,param1,param2)       
      real tfactor,tt,param

      tfactor = exp(-(2.0*(tt-param1)/param2)**2)

      return
      end


      subroutine cosine_pulse(tfactor,tt,twopi,param1)
      real tfacfor,tt,twopi,param1

      tfactor= cos(twopi*param1*tt)             

      return
      end


      subroutine rectangular_pulse(tfactor,tt,param1)
      real tfacfor,tt,param1

      if (mod(tt,param1).le.0.5*param1) then
          tfactor = 1.
      else
          tfactor = 0. 
      endif

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

      if (IFPML)    call userpml   (ix,iy,iz,iel)
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

      xx  = xm1(ix,iy,iz,iel)
      yy  = ym1(ix,iy,iz,iel)
      zz  = zm1(ix,iy,iz,iel)


      permit      = permit_0
      permea      = permea_0
      drude_alpha = 0.0               
      drude_beta  = 0.0                         

      radius = param(80)*nmscale            
      const  = param(81)
    
      r2 = radius*radius

      if  (IF3D) then

      if  (xx**2+yy**2+zz**2 .le.r2)    then 
           permit      = permit_0*permit_r
           permea      = permea_0
           drude_alpha =-gamma_p
           drude_beta  = permit_0*freq_p**2
      endif

      else

      if  (xx**2+yy**2 .le.r2)    then 
           permit      = permit_0*permit_r
           permea      = permea_0
           drude_alpha =-gamma_p
           drude_beta  = permit_0*freq_p**2
      endif

      endif

      face_value= 1.0

      if  (IF3D) then

           if  (xx**2+yy**2+zz**2 .le.r2)  face_value =-1

      else

           if  (xx**2+yy**2 .le.r2)        face_value =-1

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

      real     pi

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
          sx =nmscale
          sy =nmscale
          sz =nmscale
      else
          sx = 1.0/(xmax-xmin)
          sy = 1.0/(ymax-ymin)
          if (IF3D) sz = 1.0/(zmax-zmin)
          nmscale = 1.0/(xmax-xmin)
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

      if (nid.eq.0) write(6,*) 'xmin,ymin,zmin',xmin,ymin,zmin
      if (nid.eq.0) write(6,*) 'xmax,ymax,zmax',xmax,ymax,zmax

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

          !h2err = hxer**2 + hyer**2
          !if (h2err.gt.0) h2err = sqrt(h2err)
          !hmgp = max(hxgp,hygp)
          !hmer = max(hxmr,hymr)

          write(6,91) ISTEP,nelt,nels,nnn,npt1,npts
     $         ,time,dt,HXer,HYer,HZer,EXer,EYer,EZer,cpu_t,cpu_p_t
          write(6,95) ISTEP,nelt,nels,nnn,npt1,npts
     $         ,time,dt,HXmr,HYmr,HZmr,EXmr,EYmr,EZmr,cpu_t,cpu_p_t

 91       format(/,i5,3i5,i5,i9,1p9e10.3,e9.2,' CPU: L2  '  )
 95       format(  i5,3i5,i5,i9,1p9e10.3,e9.2,' CPU: Linf'  )

c93       format(i7,1pe10.3)

        endif ! nid=0

      endif

      return
      end

c-----------------------------------------------------------------------
