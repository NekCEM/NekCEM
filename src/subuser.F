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


      if (ndim.eq.3) then
        omega = 3             
        omega = sqrt(omega)
   
        tmpH=sin(omega*tt)/omega
        tmpE=cos(omega*tt)   

        do i=1,nx1*ny1*nz1*nelt
          xx=XM1(i,1,1,1)
          yy=YM1(i,1,1,1)
          zz=ZM1(i,1,1,1)

          myshx(i,1,1,1) = 2*cos(xx)*sin(yy)*cos(zz)*tmpH  
          myshy(i,1,1,1) = - sin(xx)*cos(yy)*cos(zz)*tmpH  
          myshz(i,1,1,1) =   sin(xx)*sin(yy)*sin(zz)*tmpH  
          mysex(i,1,1,1) =   0.0
          mysey(i,1,1,1) =   cos(xx)*sin(yy)*sin(zz)*tmpE  
          mysez(i,1,1,1) =   cos(xx)*cos(yy)*cos(zz)*tmpE  
        enddo

      else 

        omega = 2                                
        omega = sqrt(omega)

        if (iftm) then
          tmpH=sin(omega*tt)/omega
          tmpE=cos(omega*tt)   

          do i=1,nx1*ny1*nz1*nelt
            xx=XM1(i,1,1,1)
            yy=YM1(i,1,1,1)
            zz=ZM1(i,1,1,1)

            myshx(i,1,1,1) =   cos(xx)*sin(yy)*tmpH  
            myshy(i,1,1,1) = - sin(xx)*cos(yy)*tmpH  
            myshz(i,1,1,1) =   0.0
            mysex(i,1,1,1) =   0.0
            mysey(i,1,1,1) =   0.0
            mysez(i,1,1,1) =   cos(xx)*cos(yy)*tmpE  
          enddo

        elseif (ifte) then

          tmpH=cos(omega*tt)
          tmpE=sin(omega*tt)/omega

          do i=1,nx1*ny1*nz1*nelt
            xx=XM1(i,1,1,1)
            yy=YM1(i,1,1,1)
            zz=ZM1(i,1,1,1)

            myshx(i,1,1,1) =   0.0
            myshy(i,1,1,1) =   0.0
            myshz(i,1,1,1) =   sin(xx)*sin(yy)*tmpH  
            mysex(i,1,1,1) =   sin(xx)*cos(yy)*tmpE 
            mysey(i,1,1,1) =  -cos(xx)*sin(yy)*tmpE  
            mysez(i,1,1,1) =   0.0
          enddo

        else
          call exitt
        endif
      endif

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
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n, ix, iy, iz, ieg
      real     pi, one, xcc, ycc, zcc, x1, y1, z1
      real     lambda, fc, wc                     

      iel =gllel(ieg)
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

      if (if_in_region(iel).eq.0) then
          permea = permea_0             
          permit = permit_0                
      endif


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
      include 'GEOMBOUND'
      include 'SCRATCH'

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

      wavenum = 1.0
      sx = wavenum*2.0*pi/(xmax-xmin)
      sy = wavenum*2.0*pi/(ymax-ymin)       
      if (if3d) sz = wavenum*2.0*pi/(zmax-zmin)            
      nmscale = 1.0 !2.0*pi/(xmax-xmin)   ! nanoscale                                                  

      if (if3d) then

      do i=1,n
         xm1(i,1,1,1) = sx*(xm1(i,1,1,1)-xmin)                
         ym1(i,1,1,1) = sy*(ym1(i,1,1,1)-ymin)                
         zm1(i,1,1,1) = sz*(zm1(i,1,1,1)-zmin)                
      enddo

      else

      do i=1,n
         xm1(i,1,1,1) = sx*(xm1(i,1,1,1)-xmin)               
         ym1(i,1,1,1) = sy*(ym1(i,1,1,1)-ymin)                        
      enddo
      endif

      do ie=1,nelv
         if_in_region(ie)=0
      enddo

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


      if (istep.le.10.or.mod(istep,iocomm).eq.0.or.istep.eq.nsteps) then

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

        eshx= glamax(shx,NTOTT)
        eshy= glamax(shy,NTOTT)
        eshz= glamax(shz,NTOTT)
        esex= glamax(sex,NTOTT)
        esey= glamax(sey,NTOTT)
        esez= glamax(sez,NTOTT)
        
        if (eshx.eq.0)  eshx=1.e+300
        if (eshy.eq.0)  eshy=1.e+300
        if (eshz.eq.0)  eshz=1.e+300
        if (esex.eq.0)  esex=1.e+300
        if (esey.eq.0)  esey=1.e+300
        if (esez.eq.0)  esez=1.e+300
         
        HXmr = glamax(ErrHX,NTOTT) /eshx              
        HYmr = glamax(ErrHY,NTOTT) /eshy              
        HZmr = glamax(ErrHZ,NTOTT) /eshz                
        EXmr = glamax(ErrEX,NTOTT) /esex               
        EYmr = glamax(ErrEY,NTOTT) /esey                 
        EZmr = glamax(ErrEZ,NTOTT) /esez                      
 
        eshx= glsc3(SHX,bm1,SHX,NTOTT)
        eshy= glsc3(SHY,bm1,SHY,NTOTT)
        eshz= glsc3(SHZ,bm1,SHZ,NTOTT)
        esex= glsc3(SEX,bm1,SEX,NTOTT)
        esey= glsc3(SEY,bm1,SEY,NTOTT)
        esez= glsc3(SEZ,bm1,SEZ,NTOTT)
        
        if (eshx.eq.0)  eshx=1.e+300
        if (eshy.eq.0)  eshy=1.e+300
        if (eshz.eq.0)  eshz=1.e+300
        if (esex.eq.0)  esex=1.e+300
        if (esey.eq.0)  esey=1.e+300
        if (esez.eq.0)  esez=1.e+300

        HXer = glsc3(ErrHX,bm1,ErrHX,NTOTT) /eshx                     
        HYer = glsc3(ErrHY,bm1,ErrHY,NTOTT) /eshy                     
        HZer = glsc3(ErrHZ,bm1,ErrHZ,NTOTT) /eshz 
        EXer = glsc3(ErrEX,bm1,ErrEX,NTOTT) /esex
        EYer = glsc3(ErrEY,bm1,ErrEY,NTOTT) /esey
        EZer = glsc3(ErrEZ,bm1,ErrEZ,NTOTT) /esez

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

 91       format(/,i9,3i5,i5,i9,1p2e10.3,1p6e18.8,2e10.2,' L2  '  )
 95       format(  i9,3i5,i5,i9,1p2e10.3,1p6e18.8,2e10.2,' Linf'  )

c93       format(i7,1pe10.3)

        endif ! nid=0

      endif

      return
      end

