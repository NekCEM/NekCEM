c-----------------------------------------------------------------------
      subroutine q_filter(wght)
c
c     filter vx,vy,vz, and p by simple interpolation
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
c
c
c     These are the dimensions that we interpolate onto for v and p:
      parameter(lxv=lx1-1)
      parameter(lxp=lx2-1)
c
      real intdv(lx1,lx1)
      real intuv(lx1,lx1)
      real intdp(lx1,lx1)
      real intup(lx1,lx1)
      real intv(lx1,lx1)
      real intp(lx1,lx1)
c
      save intdv
      save intuv
      save intdp
      save intup
      save intv
      save intp

      common /ctmp0/ intw,intt
     $             , wk1,wk2
     $             , zgmv,wgtv,zgmp,wgtp,tmax(100),omax(103)

      real intw(lx1,lx1)
      real intt(lx1,lx1)
      real wk1  (lx1,lx1,lx1,lelt)
      real wk2  (lx1,lx1,lx1)
      real zgmv(lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)


      character*18 sfmt

      integer icalld
      save    icalld
      data    icalld /0/

      imax = nid
      imax = iglmax(imax,1)
      jmax = iglmax(imax,1)
      if (icalld.eq.0) then
         icalld = 1
         ncut = param(85) ! 1 or 2   
         call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)
      endif

      call filterq(ex,intv,nx1,nz1,wk1,wk2,intt,if3d,exmx)
      call filterq(ey,intv,nx1,nz1,wk1,wk2,intt,if3d,eymx)
      call filterq(ez,intv,nx1,nz1,wk1,wk2,intt,if3d,ezmx)

      call filterq(hx,intv,nx1,nz1,wk1,wk2,intt,if3d,hxmx)
      call filterq(hy,intv,nx1,nz1,wk1,wk2,intt,if3d,hymx)
      call filterq(hz,intv,nx1,nz1,wk1,wk2,intt,if3d,hzmx)

      write(6,1) istep,time,exmx,eymx,ezmx,hxmx,hymx,hzmx
    1 format(i8,1p7e10.3,' filt')


      return
      end
c-----------------------------------------------------------------------
      subroutine filterq(v,f,nx,nz,w1,w2,ft,if3d,dmax)
c
      include 'SIZE'
      include 'TSTEP'

      real v(nx*nx*nz,nelt),w1(1),w2(1)
      logical if3d
c
      real f(nx,nx),ft(nx,nx)
c
      integer e
c
      call transpose(ft,nx,f,nx)
c
      nxyz=nx*nx*nz
      dmax = 0.

      nel = nelv                  

      if (if3d) then
         do e=1,nel
c           Filter
            call copy(w2,v(1,e),nxyz)
            call mxm(f,nx,w2,nx,w1,nx*nx)
            i=1
            j=1
            do k=1,nx
               call mxm(w1(i),nx,ft,nx,w2(j),nx)
               i = i+nx*nx
               j = j+nx*nx
            enddo
            call mxm (w2,nx*nx,ft,nx,w1,nx)
            call sub3(w2,v(1,e),w1,nxyz)
            call copy(v(1,e),w1,nxyz)
            smax = vlamax(w2,nxyz)
            dmax = max(dmax,abs(smax))
         enddo
c
      else
         do e=1,nel
c           Filter
            call copy(w1,v(1,e),nxyz)
            call mxm(f ,nx,w1,nx,w2,nx)
            call mxm(w2,nx,ft,nx,w1,nx)
c
            call sub3(w2,v(1,e),w1,nxyz)
            call copy(v(1,e),w1,nxyz)
            smax = vlamax(w2,nxyz)
            dmax = max(dmax,abs(smax))
         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outmatx(a,m,n,io,name)
      real a(m*n)
      character*4 name
c
      open(unit=io,file=name)
      do i=1,m*n
         write(io,1) a(i)
      enddo
    1 format(1p1e22.13)
      close(unit=io)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine add2sxy(x,a,y,b,n)
      real x(1),y(1)
c
      do i=1,n
         x(i) = a*x(i) + b*y(i)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_new_filter(intv,zpts,nx,kut,wght,nid)
c
c     This routing builds a 1D filter with a transfer function that
c     looks like:
c
c
c        ^
c    d_k |
c        |                 |
c     1  |__________      _v_
c        |          -_     
c        |            \  wght
c        |             \  ___
c        |             |   ^
c     0  |-------------|---|>
c
c        0         c   N   k-->
c
c        Where c := N-kut is the point below which d_k = 1.
c
c
c
c      Here, nx = number of points
c
      real intv(nx,nx),zpts(nx)
c
      parameter (lm=40)
      parameter (lm2=lm*lm)
      real      phi(lm2),pht(lm2),diag(lm2),rmult(lm),Lj(lm)
      integer   indr(lm),indc(lm),ipiv(lm)
c
      if (nx.gt.lm) then
         write(6,*) 'ABORT in build_new_filter:',nx,lm
         call exitt
      endif
c
      kj = 0
      n  = nx-1
      do j=1,nx
         z = zpts(j)
         call legendre_poly(Lj,z,n)
         kj = kj+1
         pht(kj) = Lj(1)
         kj = kj+1
         pht(kj) = Lj(2)
         do k=3,nx
            kj = kj+1
            pht(kj) = Lj(k)-Lj(k-2)
         enddo
      enddo
      call transpose (phi,nx,pht,nx)
      call copy      (pht,phi,nx*nx)
      call gaujordf  (pht,nx,nx,indr,indc,ipiv,ierr,rmult)
c
c     Set up transfer function
c
      call ident   (diag,nx)
c
      k0 = nx-kut
      do k=k0+1,nx
         kk = k+nx*(k-1)
         amp = wght*(k-k0)*(k-k0)/(kut*kut)   ! quadratic growth
         diag(kk) = 1.-amp
      enddo
c
      call mxm  (diag,nx,pht,nx,intv,nx)      !          -1
      call mxm  (phi ,nx,intv,nx,pht,nx)      !     V D V
      call copy (intv,pht,nx*nx)
c
      do k=1,nx*nx
         pht(k) = 1.-diag(k)
      enddo
      np1 = nx+1
      if (nid.eq.0) then
         write(6,6) 'filt amp',(pht (k),k=1,nx*nx,np1)
         write(6,6) 'filt trn',(diag(k),k=1,nx*nx,np1)
   6     format(a8,16f7.4,6(/,8x,16f7.4))
      endif
c
      return
      end
c-----------------------------------------------------------------------
