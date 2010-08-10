
c-----------------------------------------------------------------------
c  USER SPECIFIED ROUTINES:
c
c     - boundary conditions
c     - initial conditions
c     - variable properties
c     - forcing function for fluid (f)
c     - forcing function for passive scalar (q)
c     - general purpose routine for checking errors etc.
c-----------------------------------------------------------------------

      subroutine userinc
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      return
      end

c-----------------------------------------------------------------------
      subroutine usersol(tt, myshx, myshy, myshz, mysex, mysey, mysez)
      ! cf. David A. de Wolf, Essentials of Electromagnetics for
      ! Engineering, Cambridge University Press, 2001
      ! Section 19.7
      ! This is the TM wave case.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'GEOMBOUND'

      real tt
      real myshx(lx1,ly1,lz1,lelt)
      real myshy(lx1,ly1,lz1,lelt)
      real myshz(lx1,ly1,lz1,lelt)
      real mysex(lx1,ly1,lz1,lelt)
      real mysey(lx1,ly1,lz1,lelt)
      real mysez(lx1,ly1,lz1,lelt)

      ! Mode selection
      integer m, rootidx
      real tmscale
      parameter (m = 0)         ! Index of the Bessel function
      parameter (rootidx = 1)   ! Index of Bessel J_nu root
      parameter (tmscale = 1.0) ! Scale factor of TM component

      ! Extent of domain
      real radius
      save radius

      ! Local coordinates
      real myx, myy, myz
      real rho, phi

      ! Following variables named as in reference.
      real kz, bigk, omega
      real erho, ephi, ezz
      real hrho, hphi, hzz

      ! Helper variables
      complex zfactor, tfactor, allfac

      integer zwaves

      ! Loop variables
      integer i

      ! Helper for initialize-once variables
      logical icalld
      data icalld/.false./
      save icalld

      complex iunit

      iunit = (0,1)

      n= nx1*ny1*nz1*nelt

      ! compute radius
      if (.not.icalld) then
        radius = geom_xyradius()
        xmax=glmax(xm1,n)
        ymax=glmax(ym1,n)
        zmax=glmax(zm1,n)
        xmin=glmin(xm1,n)
        ymin=glmin(ym1,n)
        zmin=glmin(zm1,n)
        icalld = .true.

        do i=1,nx1*ny1*nz1*nelt
           xx = xm1(i,1,1,1)
           yy = ym1(i,1,1,1)
           zz = zm1(i,1,1,1)
           rho = sqrt (xx**2+yy**2)
           if (rho.eq.0) then
             if (nid.eq.0) 
     $           write(6,*) 'change N: radius zero/solution infinity'
               call exitt
          ! NOTE:
          ! The analyic solution has a discontinuity at zero,
          ! so there really is no "correct" value there.
          ! (Consider approach to zero along X- and Y-axis.)
          ! The best, but still fishy, thing to do is to put 
          ! all fields to zero.

          ! Alternate is to set as below
            ezz  = 0
            erho = 0
            ephi = 0
            hrho = 0
            hphi = 0
           endif
        enddo
      endif

      zsize = zmax-zmin
      zwaves= 1.0  !   int(param(93))
      if (abs(zsize).le.1e-20) then
          kz = 0
      else
          kz = 2*pi*zwaves/zsize
      endif

      do i=1,nx1*ny1*nz1*nelt
        myx = xm1(i,1,1,1)
        myy = ym1(i,1,1,1)
        myz = zm1(i,1,1,1)

        ! find local coordinates
        rho = sqrt (myx**2 + myy**2)
        phi = atan2(myy, myx)

        ! compute solution
        bigk = bssljrt(m, rootidx)/radius
        omega = sqrt(kz**2 + bigk**2)*cspeed
            
        zfactor = cexp((0,-1)*kz*myz)
        tfactor = cexp((0,1)*omega*tt)
        allfac = tmscale * zfactor * tfactor
            
        if (rho.ne.0) then
        ezz  = real(allfac * bsslj(m, bigk*rho) * cos(m*phi))
        erho = real(allfac *(-iunit)*kz/bigk * bssljp(m, bigk*rho)
     $                     * cos(m*phi))
        ephi = real(allfac * iunit*m*kz/bigk**2/rho
     $                     * bsslj(m,bigk*rho) * sin(m*phi))
        hrho = real(allfac *(-iunit)*m*omega*permit_0/bigk**2/rho
     $                     * bsslj(m,bigk*rho) * sin(m*phi))
        hphi = real(allfac *(-iunit)*omega*permit_0/bigk
     $                     * bssljp(m,bigk*rho) * cos(m*phi))
        hzz = 0. ! TM case
           
        myshx(i,1,1,1) = cos(phi)*hrho - sin(phi)*hphi
        myshy(i,1,1,1) = sin(phi)*hrho + cos(phi)*hphi
        myshz(i,1,1,1) = hzz
        mysex(i,1,1,1) = cos(phi)*erho - sin(phi)*ephi
        mysey(i,1,1,1) = sin(phi)*erho + cos(phi)*ephi
        mysez(i,1,1,1) = ezz
        endif

      enddo

      return
      end
       
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real xsize, ysize, zsize
      real xmin,xmax,ymin,ymax,zmin,zmax
      real radius

      real e1

      integer n = 8*lelt

      radius = geom_xyradius()
      radius = int(10*radius+0.1)/10. 
      e1 = radius*radius - 1e-4

      do ie=1,nelt
      do i=1,8
        rr = xc(i,ie)*xc(i,ie) + yc(i,ie)*yc(i,ie)
        if (rr.gt.e1) then
           rr = sqrt(rr)
           rn = radius/rr
           xc(i,ie) = rn*xc(i,ie)
           yc(i,ie) = rn*yc(i,ie)
        endif
      enddo
      enddo
    1 format(2i5,1p3e18.9,a5)

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      n = nx1*ny1*nz1*nelv

      xmax=glmax(xm1,n)
      ymax=glmax(ym1,n)
      zmax=glmax(zm1,n)
      xmin=glmin(xm1,n)
      ymin=glmin(ym1,n)
      zmin=glmin(zm1,n)


      sz=2*pi*xmax/(zmax-zmin)

      if (if3d) then
        do i=1,n
           zm1(i,1,1,1) =sz*(zm1(i,1,1,1)-zmin)               
        enddo
      endif

      xmax=glmax(xm1,n)
      ymax=glmax(ym1,n)
      zmax=glmax(zm1,n)
      xmin=glmin(xm1,n)
      ymin=glmin(ym1,n)
      zmin=glmin(zm1,n)

      if (nid.eq.0) then
      write(6,*) 'xmax,ymax,zmax :: ', xmax,ymax,zmax
      write(6,*) 'xmin,ymin,zmin :: ', xmin,ymin,zmin
      endif

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
      subroutine usersrc(baseidx,srchx,srchy,srchz,srcex,srcey,srcez)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real srchx(lx1,ly1,lz1,lelt),
     $     srchy(lx1,ly1,lz1,lelt),
     $     srchz(lx1,ly1,lz1,lelt),
     $     srcex(lx1,ly1,lz1,lelt),
     $     srcey(lx1,ly1,lz1,lelt),
     $     srcez(lx1,ly1,lz1,lelt) 

      integer baseidx

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy, iz, ieg
   
      iel     = gllel(ieg)
      call usermedia   (ix,iy,iz,iel)

      return
      end

c-----------------------------------------------------------------------

      subroutine usermedia(ix,iy,iz,iel) ! Vector form of userf
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy,iz, iel

      permea = permea_0             
      permit = permit_0                
      if_in_region(iel)=0


      return
      end

c-----------------------------------------------------------------------
      subroutine userdrude(ix,iy,iz,iel) ! Vector form of userf

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy,iz, iel

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
     $             ,time,dt,HXer,HYer,HZer,EXer,EYer,EZer,cpu_t,cpu_p_t

          write(6,95) ISTEP,nelt,nels,nnn,npt1,npts
     $             ,time,dt,HXmr,HYmr,HZmr,EXmr,EYmr,EZmr,cpu_t,cpu_p_t

 91       format(/,i9,3i5,i5,i9,1p9e10.3,e9.2,' CPU: L2  ')
 95       format(/,i9,3i5,i5,i9,1p9e10.3,e9.2,' CPU: Linf')


        endif ! nid=0

      endif


      return
      end


