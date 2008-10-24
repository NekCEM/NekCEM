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
      temporal_freq= param(64)*cspeed  ! tmodty=3: cosine pulse

c ... temporal mode set  
      duration     = param(65)  ! tmodty=4: rectangular pulse

c ... spatial mode set  
      xloc  = param(66)                 ! smodty=2: gaussian 
      yloc  = param(67)          
      zloc  = param(68)             
      xw    = param(69)  
      yw    = param(70)  
      zw    = param(71)  
      xwidth= param(72)          
      ywidth= param(73)            
      zwidth= param(74)            

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
     $	         fldhere = fldhere*exp(-(2.0*(zz-zgauss)/zwidth)**2)

             fullary(i,1,1,1) = tfactor * fldhere*bm1(i,1,1,1)

          enddo

      elseif (smodtyp.eq.3) then

          fldhere = tfactor*srcscal

          do i=1,nxyzl

             xx = xm1(i,1,1,1)
             yy = ym1(i,1,1,1)
             zz = zm1(i,1,1,1)

             fullary(i,1,1,1)=fldhere
     $       *sin(twopi/wavelength*(cspeed*rktime-(kx*xx+ky*yy+kz*zz)))   

             fullaryprime(i,1,1,1)=fldhere                                                              
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
      subroutine usersol(tt, myshx, myshy, myshz, mysex, mysey, mysez)
c-----------------------------------------------------------------------
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
      subroutine userini(tt, myhx, myhy, myhz, myex, myey, myez)
c-----------------------------------------------------------------------
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

      n = nx1*ny1*nz1*nelt
      do i=1,n                 
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

      integer  n, ix, iy, iz, ieg
      real     pi, one, xcc, ycc, zcc, x1, y1, z1
      real     lambda, fc, wc                     

      iel= gllel(ielg)

      call usermedia   (ix,iy,iz,iel)

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

      permit   = permit_0
      permea   = permea_0
      drude_alpha = 0.0
      drude_beta  = 0.0

      return
      end

c-----------------------------------------------------------------------
      subroutine userdrude(ix,iy,iz,iel) ! Vector form of userf

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      return
      end

c-----------------------------------------------------------------------
      subroutine userpml(ix,iy,iz,iel) ! Vector form of userf
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

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
      subroutine useric (ix,iy,iz,ieg)
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

      nmscale = 1.0 !2*pi*1.0/(xmax-xmin)   ! nanoscale                                                  
      sx = 1.0*nmscale !2*pi*1.0/(xmax-xmin)
      sy = 1.0*nmscale !2*pi*1.0/(ymax-ymin)       
      if (if3d) sz =1.0*nmscale! 2.*pi/(zmax-zmin)            

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

      do ie=1,nelv

         xmid = xm1(nxy2,nxy2,nxy2,ie)
         ymid = ym1(nxy2,nxy2,nxy2,ie)
         zmid = zm1(nxy2,nxy2,nxy2,ie)
         rmd2 = xmid*xmid + ymid*ymid

         if_in_region(ie) = 0               ! default: free space

      enddo


      return
      end

c-----------------------------------------------------------------------
      subroutine userchk    
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      common /ccpu/ cpu_t,cpu_t0

c ... check maximum dt for every 10 steps. 
      n    = nx1*ny1*nz1*nelt
      nsrf = nx1*nz1*nelt*6
      idt  = 10

c ... regular checks for iocomm 
      if (istep.le.10.or.mod(istep,iocomm).eq.0) dochk=.true.    

c ... compute errors
      if (dochk) then

      call rzero(shx,n) 
      call rzero(shy,n) 
      call rzero(shz,n) 
      call rzero(sex,n) 
      call rzero(sey,n) 
      call rzero(sez,n) 

      if (ifsol) call usersol(time,shx,shy,shz,sex,sey,sez)

      call sub3(errhx,shx,hx,n)         
      call sub3(errhy,shy,hy,n)         
      call sub3(errhz,shz,hz,n)         
      call sub3(errex,sex,ex,n)         
      call sub3(errey,sey,ey,n)         
      call sub3(errez,sez,ez,n)         

      HXmr = glamax(errhx,n)        
      HYmr = glamax(errhy,n)        
      HZmr = glamax(errhz,n)        
      EXmr = glamax(errex,n)        
      EYmr = glamax(errey,n)        
      EZmr = glamax(errez,n)        

      HXer = glsc3(errhx,bm1,errhx,n)/volvm1            
      HYer = glsc3(errhy,bm1,errhy,n)/volvm1            
      HZer = glsc3(errhz,bm1,errhz,n)/volvm1            
      EXer = glsc3(errex,bm1,errex,n)/volvm1            
      EYer = glsc3(errey,bm1,errey,n)/volvm1            
      EZer = glsc3(errez,bm1,errez,n)/volvm1            

      if (hxer.gt.0) hxer = sqrt(hxer)
      if (hyer.gt.0) hyer = sqrt(hyer)
      if (hzer.gt.0) hzer = sqrt(hzer)
      if (exer.gt.0) exer = sqrt(exer)
      if (eyer.gt.0) eyer = sqrt(eyer)
      if (ezer.gt.0) ezer = sqrt(ezer)

c ... compute max gap size at interfaces between elmts.

      call cem_gap                     

      HXgp = glamax(ghx,nsrf) 
      HYgp = glamax(ghy,nsrf) 
      HZgp = glamax(ghz,nsrf) 
      EXgp = glamax(gex,nsrf) 
      EYgp = glamax(gey,nsrf) 
      EZgp = glamax(gez,nsrf) 

c ... print statement

      snel= nelt + .5
      snel= sqrt(snel)
      nels= snel
      nnn = nx1-1            ! degree 
      npt1= nnn*nels
      npts= nx1*ny1*nz1*nelt ! npts = npt1**ndim

      if (nid.eq.0) then

          h2err = hxer**2 + hyer**2
          if (h2err.gt.0) h2err = sqrt(h2err)

          hmgp = max(hxgp,hygp)
          hmer = max(hxmr,hymr)

          write(6,91) ISTEP,nelt,nels,nnn,npt1,npts
     $    ,time,dt,HXer,HYer,HZer,EXer,EYer,EZer,cpu_t       
          write(6,92) ISTEP,nelt,nels,nnn,npt1,npts
     $    ,time,dt,HXmr,HYmr,HZmr,EXmr,EYmr,EZmr,cpu_t      

          if (ifvardt) then            
              dtstd= HXmr+HYmr+HZmr+EXmr+EYmr+EZmr
              dtstd= abs(dtstd)
              cntdt= 1./idtchg
              write(6,93) idtchg,cntdt
              if (dtstd.gt.fldmax) dtchk=.true.
              if (dtstd.le.fldmax) dtchk=.false.
              if (dtstd.gt.fldmax) write(6,94) fldmax      
          endif
      endif ! nid=0

      endif ! ifdochk

 91   format(/,i5,3i5,i5,i9,1p9e10.3,' CPU: L2  ')
 92   format(/,i5,3i5,i5,i9,1p9e10.3,' CPU: Linf')
 93   format('  dxmin/(cspeed*',i2,')',', s=',f5.3)      
 94   format('  Reduce DT <- maximum norm is bigger than:',e12.7)      

      return
      end
c-----------------------------------------------------------------------
