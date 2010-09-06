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
      subroutine userinc
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'


      return
      end
c-----------------------------------------------------------------------
      subroutine usersol(tt, myshx, myshy, myshz, mysex, mysey, mysez)
c-----------------------------------------------------------------------
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

      real     omega, tmpH, tmpE
      real     xx, yy, zz, u    
      integer  i, j, k, e, ii, ind                 


      integer nn, mm ! Mode selection, as in reference
      ! m >= 0, n>= 0
      real    tmscale
      logical ifcosine ! Mode selection, as in reference
      parameter (tmscale = 1) ! Scale factor of TM component

c     Local coordinates
      real x, y, z

c     Following variables named as in reference.
      real alpha, beta, cplus, cminus

      complex iunit

      iunit = (0,1)

      nn    = param(80)
      mm    = param(81)
      ifcosine = param(82).ne.0

      do i=1,nx1*ny1*nz1*nelv
        ! find local coordinates
        x=XM1(i,1,1,1)
        y=YM1(i,1,1,1)
        z=ZM1(i,1,1,1)

        ! compute helper quantities
        alpha = 2*pi*nn
        if (ifcosine) then
          cplus = 1
          cminus = 1
          beta = pi*(2*mm+1)
        else
          cplus = 1
          cminus = -1
          beta = 2*pi*mm
        endif

        omega = (alpha**2+beta**2)**0.5

        mysez(i,1,1,1) = tmscale * real(exp((0,1)*(alpha*x+omega*tt))
     $                * (cplus*exp((0,1)*beta*y)
     $                   +cminus*exp(-(0,1)*beta*y)))
        mysex(i,1,1,1) = 0
        mysey(i,1,1,1) = 0
        myshx(i,1,1,1) = tmscale * (-beta)/omega
     $                   *real(exp(iunit*(alpha*x+omega*tt))
     $                    *(cplus*exp(iunit*beta*y)
     $                     -cminus*exp((-iunit)*beta*y)))
        myshy(i,1,1,1) = tmscale * alpha/omega
     $                   *real(exp(iunit*(alpha*x+omega*tt))
     $                    *(cplus*exp(iunit*beta*y)
     $                     +cminus*exp((-iunit)*beta*y)))
        myshz(i,1,1,1) = 0
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

      subroutine uservp (ix,iy,iz,iel)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy, iz, ieg
      real     pi, one, xcc, ycc, zcc, x1, y1, z1
      real     lambda, fc, wc                     

      call usermedia   (ix,iy,iz,iel)
      if (IFDRUDE)  call userdrude (ix,iy,iz,iel)

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
      subroutine useric  (ix,iy,iz,ieg)
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

      n    = nx1*ny1*nz1*nelt

c     ifxyo= .true.

      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)
      zmin = glmin(zm1,n)
      zmax = glmax(zm1,n)

      !if (nid.eq.0) write(6,*) wavenumber,' This is wavenumber'


      if (IFNM) then
          sx = 1.0
          sy = 1.0
	  sz = 1.0
          nmscale = 1.0               ! nanoscale                                                  
          dt      = dt*10.0**(-18)                    
      else
          sx = 1.0/(xmax-xmin)
          sy = 1.0/(ymax-ymin)       
          if (if3d) sz = 1.0/(zmax-zmin)            
          nmscale = 1.0               ! nanoscale                                                  
          dt      = dt                                
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

      write(6,*) 'xmin,xmax',xmin,xmax,ymin,ymax

      return
      end
c-----------------------------------------------------------------------
      subroutine userft           
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      common /ccpu/ cpu_t,cpu_t0

c ... compute errors
      write(6,*) 'iocomm',iocomm
        call exitt

      if ((istep.le.10).or.(mod(istep,iocomm).eq.0).or.
     $    (istep.eq.nsteps)) then

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

          write(6,91) ISTEP,nelt,nels,nnn,npt1,npts
     $         ,time,dt,HXer,HYer,HZer,EXer,EYer,EZer,cpu_t,cpu_p_t
          write(6,95) ISTEP,nelt,nels,nnn,npt1,npts
     $         ,time,dt,HXmr,HYmr,HZmr,EXmr,EYmr,EZmr,cpu_t,cpu_p_t

 91       format(/,i9,3i5,i5,i9,1p2e10.3,1p6e18.8,2e12.4,' L2  '  )
 95       format(  i9,3i5,i5,i9,1p2e10.3,1p6e18.8,2e12.4,' Linf'  )

        endif ! nid=0

      endif

      return
      end





