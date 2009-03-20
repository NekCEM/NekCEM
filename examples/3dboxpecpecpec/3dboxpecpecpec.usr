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

      if (.not.ifpec) then
         write(6,*) 'ERROR:: this .usr file is for PEC boundary'
         call exitt
      endif

      if (ndim.eq.3) then

        ww    = pi                       

        omega = 2                                
        sqrt2 = sqrt(omega)
        omega = 3                                
        sqrt3 = sqrt(omega)
        omega = 6                                
        sqrt6 = sqrt(omega)

        tmpH=cos(ww*sqrt3*tt)/sqrt6
        tmpE=sin(ww*sqrt3*tt)/sqrt2   

        do i=1,nx1*ny1*nz1*nelt
          xx=XM1(i,1,1,1)
          yy=YM1(i,1,1,1)
          zz=ZM1(i,1,1,1)

          myshx(i,1,1,1) =  -sin(ww*xx)*cos(ww*yy)*cos(ww*zz)*tmpH
          myshy(i,1,1,1) =  -cos(ww*xx)*sin(ww*yy)*cos(ww*zz)*tmpH
          myshz(i,1,1,1) = 2*cos(ww*xx)*cos(ww*yy)*sin(ww*zz)*tmpH
          mysex(i,1,1,1) =  -cos(ww*xx)*sin(ww*yy)*sin(ww*zz)*tmpE
          mysey(i,1,1,1) =   sin(ww*xx)*cos(ww*yy)*sin(ww*zz)*tmpE
          mysez(i,1,1,1) =   0 
        enddo

      else 

        ww    = 1*pi                       

        omega = 2                                
        omega = sqrt(omega)

        tmpH=sin(ww*omega*tt)/omega
        tmpE=cos(ww*omega*tt)   

        if (iftm) then
          do i=1,nx1*ny1*nz1*nelt
            xx=XM1(i,1,1,1)
            yy=YM1(i,1,1,1)
            zz=ZM1(i,1,1,1)

            myshx(i,1,1,1) =  cos(ww*xx)*sin(ww*yy)*tmpH  
            myshy(i,1,1,1) = -sin(ww*xx)*cos(ww*yy)*tmpH  
            myshz(i,1,1,1) =  0.0
            mysex(i,1,1,1) =  0.0
            mysey(i,1,1,1) =  0.0
            mysez(i,1,1,1) =  cos(ww*xx)*cos(ww*yy)*tmpE  
          enddo

        elseif (ifte) then

          do i=1,nx1*ny1*nz1*nelt
            xx=XM1(i,1,1,1)
            yy=YM1(i,1,1,1)
            zz=ZM1(i,1,1,1)

            myshx(i,1,1,1) =   0.0
            myshy(i,1,1,1) =   0.0
            myshz(i,1,1,1) =  sin(ww*xx)*sin(ww*yy)*tmpE   
            mysex(i,1,1,1) =  sin(ww*xx)*cos(ww*yy)*tmpH  
            mysey(i,1,1,1) = -cos(ww*xx)*sin(ww*yy)*tmpH  
            mysez(i,1,1,1) =   0.0
          enddo

        else
          call exitt
        endif
      endif

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

      if (IFPML)    call userpml   (ix,iy,iz,iel)
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
      include 'EMWAVE'
 
      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
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

c     if (nid.eq.0) write(6,*) wavenumber,' This is wavenumber'

      sx = 2.0/(xmax-xmin)
      sy = 2.0/(ymax-ymin)       
      if (if3d) sz = 2.0/(zmax-zmin)            
      nmscale = 1.0               ! nanoscale                                                  

      if (if3d) then

      do i=1,n
         xm1(i,1,1,1) = sx*(xm1(i,1,1,1)-xmin)-1                
         ym1(i,1,1,1) = sy*(ym1(i,1,1,1)-ymin)-1                
         zm1(i,1,1,1) = sz*(zm1(i,1,1,1)-zmin)-1                
      enddo

      else

      do i=1,n
         xm1(i,1,1,1) = sx*xm1(i,1,1,1)               
         ym1(i,1,1,1) = sy*ym1(i,1,1,1)                        
      enddo
      endif

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
      
      integer errunit
      parameter (errunit=77)

c     Helper for initialize-once variables
      logical icalld
      data icalld/.false./
      save icalld


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

 91       format(/,i9,3i5,i5,i9,1p10e12.3,' CPU: L2  '  )
 95       format(  i9,3i5,i5,i9,1p10e12.3,' CPU: Linf'  )

        endif ! nid=0

      endif

      return
      end


