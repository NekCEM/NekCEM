c-----------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,N,e,D,Dt)
c     Output: ur,us,ut         Input:u,N,e,D,Dt
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real u (0:N,0:N,0:N,1)
      real D (0:N,0:N),Dt(0:N,0:N)
      integer e
c
      m1 = N+1
      m2 = m1*m1
c
      call mxm(D ,m1,u(0,0,0,e),m1,ur,m2)
      do k=0,N
         call mxm(u(0,0,k,e),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u(0,0,0,e),m2,Dt,m1,ut,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad2(ur,us,u,N,e,D,Dt)
c     Output: ur,us         Input:u,N,e,D,Dt
      real ur(0:N,0:N),us(0:N,0:N)
      real u (0:N,0:N,1)
      real D (0:N,0:N),Dt(0:N,0:N)
      integer e
c
      m1 = N+1
c
      call mxm (D ,m1,u(0,0,e),m1,ur,m1)
      call mxm (u(0,0,e),m1,Dt,m1,us,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gradm1(ux,uy,uz,u)
c
c     Compute gradient of T -- mesh 1 to mesh 1 (vel. to vel.)
c
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
c
      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz,1),uy(lxyz,1),uz(lxyz,1),u(lxyz,1)
c
      common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz)
c
      common /dudxyj/ jacmi(lx1*ly1*lz1,lelt)
      real jacmi
c
      integer e
c
      integer icalld
      save    icalld
      data    icalld /0/
C
      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelt
C
      if (icalld.eq.0) then
         call invers2(jacmi,jacm1,ntot)
         icalld=1
      endif

c
      N = nx1-1
      do e=1,nelt
         if (if3d) then
            call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
            do i=1,lxyz
               ux(i,e) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                             + us(i)*sxm1(i,1,1,e)
     $                             + ut(i)*txm1(i,1,1,e) )
               uy(i,e) = jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                             + us(i)*sym1(i,1,1,e)
     $                             + ut(i)*tym1(i,1,1,e) )
               uz(i,e) = jacmi(i,e)*(ur(i)*rzm1(i,1,1,e)
     $                             + us(i)*szm1(i,1,1,e)
     $                             + ut(i)*tzm1(i,1,1,e) )
            enddo
         else
            if (ifaxis) call setaxdy (ifrzer(e))
            call local_grad2(ur,us,u,N,e,dxm1,dytm1)
            do i=1,lxyz
               ux(i,e) =jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)
     $                            + us(i)*sxm1(i,1,1,e) )
               uy(i,e) =jacmi(i,e)*(ur(i)*rym1(i,1,1,e)
     $                            + us(i)*sym1(i,1,1,e) )
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outpost(v1,v2,v3,vp,vt,name3)
c
      include 'SIZE'
      include 'SOLN'
      parameter(ltot1=lx1*ly1*lz1*lelv)
      parameter(ltot2=lx2*ly2*lz2*lelv)
      common /scrout1/ w1(ltot1),w2(ltot1),w3(ltot1),wp(ltot2)
      common /scrout2/ wt(ltot1)
c
      real v1(1),v2(1),v3(1),vp(1),vt(1)
      character*3 name3
c
      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv
c
      call copy(w1,vx,ntot1)
      call copy(w2,vy,ntot1)
      call copy(w3,vz,ntot1)
      call copy(wp,pr,ntot2)
      call copy(wt,t ,ntot1)
c
      call copy(vx,v1,ntot1)
      call copy(vy,v2,ntot1)
      call copy(vz,v3,ntot1)
      call copy(pr,vp,ntot2)
      call copy(t ,vt,ntot1)
c
c     if (nid.eq.0) write(6,*) 'calling prepost ',name3
      call prepost(.true.,name3)
c
      call copy(vx,w1,ntot1)
      call copy(vy,w2,ntot1)
      call copy(vz,w3,ntot1)
      call copy(pr,wp,ntot2)
      call copy(t ,wt,ntot1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine comp_vort3(vort,work1,work2,u,v,w)
c
      include 'SIZE'
      include 'TOTAL'
c
      parameter(lt=lx1*ly1*lz1*lelv)
      real vort(lt,3),work1(1),work2(1),u(1),v(1),w(1)
c
      ntot  = nx1*ny1*nz1*nelv
      if (if3d) then
c        work1=dw/dy ; work2=dv/dz
           call dudxyz(work1,w,rym1,sym1,tym1,jacm1,1,2)
           call dudxyz(work2,v,rzm1,szm1,tzm1,jacm1,1,3)
           call sub3(vort(1,1),work1,work2,ntot)
c        work1=du/dz ; work2=dw/dx
           call dudxyz(work1,u,rzm1,szm1,tzm1,jacm1,1,3)
           call dudxyz(work2,w,rxm1,sxm1,txm1,jacm1,1,1)
           call sub3(vort(1,2),work1,work2,ntot)
c        work1=dv/dx ; work2=du/dy
           call dudxyz(work1,v,rxm1,sxm1,txm1,jacm1,1,1)
           call dudxyz(work2,u,rym1,sym1,tym1,jacm1,1,2)
           call sub3(vort(1,3),work1,work2,ntot)
      else
c        work1=dv/dx ; work2=du/dy
           call dudxyz(work1,v,rxm1,sxm1,txm1,jacm1,1,1)
           call dudxyz(work2,u,rym1,sym1,tym1,jacm1,1,2)
           call sub3(vort,work1,work2,ntot)
      endif
c
c    Avg at bndry
c
      ifielt = ifield
      ifield = 1
      if (if3d) then
         do idim=1,ndim
            call col2  (vort(1,idim),bm1,ntot)
            call dssum (vort(1,idim),nx1,ny1,nz1)
            call col2  (vort(1,idim),binvm1,ntot)
         enddo
      else
         call col2  (vort,bm1,ntot)
         call dssum (vort,nx1,ny1,nz1)
         call col2  (vort,binvm1,ntot)
      endif
      ifield = ifielt
c
      return
      end
c-----------------------------------------------------------------------
      subroutine surface_int(sint,sarea,a,ie,iface1)
C
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'TOPOL'
      real a(lx1,ly1,lz1,1)
c
      integer icalld
      save    icalld
      data    icalld/0/
      logical ifpf
      save    ifpf
c
      if (icalld.eq.0) then
         icalld=icalld+1
         if (skpdat(1,2).eq.nx1) then
c           write(6,*) 'In surface_int, using pf version of skpdat.'
            ifpf = .true.
         else
c           write(6,*) 'In surface_int, using std version of skpdat.'
            ifpf = .false.
         endif
      endif
C
      sarea = 0.
      sint  = 0.
C
      call dsset(nx1,ny1,nz1)
      iface  = eface1(iface1)
c
c     Check skpdat (because of difference in pf vs. commercial version...arrghh)
c
      if (ifpf) then
c        pf version
         js1    = skpdat(1,iface)
         jf1    = skpdat(2,iface)
         jskip1 = skpdat(3,iface)
         js2    = skpdat(4,iface)
         jf2    = skpdat(5,iface)
         jskip2 = skpdat(6,iface)
      else
c        std version
         js1    = skpdat(iface,1)
         jf1    = skpdat(iface,2)
         jskip1 = skpdat(iface,3)
         js2    = skpdat(iface,4)
         jf2    = skpdat(iface,5)
         jskip2 = skpdat(iface,6)
      endif
C
      I = 0
      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         I = I+1
         sarea = sarea+area(i,1,iface1,ie)
         sint  = sint +area(i,1,iface1,ie)*a(j1,j2,1,ie)
  100 continue
C
      return
      end
c-----------------------------------------------------------------------
      subroutine surface_flux(dq,qx,qy,qz,ie,iface,w)
C
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      parameter (l=lx1*ly1*lz1)
      real w(lx1,ly1,lz1),qx(l,1),qy(l,1),qz(l,1)
c
      integer icalld
      save    icalld
      data    icalld/0/
      logical ifpf
      save    ifpf
c
      call dsset(nx1,ny1,nz1)
      if (icalld.eq.0) then
         icalld=icalld+1
         if (skpdat(1,2).eq.nx1) then
            write(6,*) 'In surface_flux, using pf version of skpdat.'
            ifpf = .true.
         else
            write(6,*) 'In surface_flux, using std version of skpdat.'
            ifpf = .false.
         endif
      endif
C
      ifacepf  = eface1(iface)
c
c     Check skpdat (because of difference in pf vs. commercial version...arrghh)
c
      if (ifpf) then
c        pf version
         js1    = skpdat(1,ifacepf)
         jf1    = skpdat(2,ifacepf)
         jskip1 = skpdat(3,ifacepf)
         js2    = skpdat(4,ifacepf)
         jf2    = skpdat(5,ifacepf)
         jskip2 = skpdat(6,ifacepf)
      else
c        std version
         js1    = skpdat(ifacepf,1)
         jf1    = skpdat(ifacepf,2)
         jskip1 = skpdat(ifacepf,3)
         js2    = skpdat(ifacepf,4)
         jf2    = skpdat(ifacepf,5)
         jskip2 = skpdat(ifacepf,6)
      endif
C
      call faccl3 (w,qx(1,ie),unx(1,1,iface,ie),iface)
      call faddcl3(w,qy(1,ie),uny(1,1,iface,ie),iface)
      if (if3d)
     $call faddcl3(w,qz(1,ie),unz(1,1,iface,ie),iface)
c
      dq = 0
      i  = 0
      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
         i = i+1
         dq    = dq   +area(i,1,iface,ie)*w(j1,j2,1)
  100 continue
C
      return
      end
c-----------------------------------------------------------------------
      subroutine legendre_poly(L,x,N)
c
c     Evaluate Legendre polynomials of degrees 0-N at point x
c
      real L(0:N)
c
      L(0) = 1.
      L(1) = x
c
      do j=2,N
         L(j) = ( (2*j-1) * x * L(j-1) - (j-1) * L(j-2) ) / j 
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dudxyz (du,u,rm1,sm1,tm1,jm1,imsh,isd)
C--------------------------------------------------------------
C
C     DU   - dU/dx or dU/dy or dU/dz
C     U    - a field variable defined on mesh 1
C     RM1  - dr/dx or dr/dy or dr/dz  
C     SM1  - ds/dx or ds/dy or ds/dz
C     TM1  - dt/dx or dt/dy or dt/dz
C     JM1  - the Jacobian   
C     IMESH - topology: velocity (1) or temperature (2) mesh
C
C--------------------------------------------------------------
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
C
      REAL  DU  (LX1,LY1,LZ1,1)
      REAL  U   (LX1,LY1,LZ1,1)
      REAL  RM1 (LX1,LY1,LZ1,1)
      REAL  SM1 (LX1,LY1,LZ1,1)
      REAL  TM1 (LX1,LY1,LZ1,1)
      REAL  JM1 (LX1,LY1,LZ1,1)
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
C     Store the inverse jacobian to speed this operation up
C
      COMMON /DUDXYJ/ JACMI(LX1,LY1,LZ1,LELT)
      REAL JACMI
C
      REAL  DRST(LX1,LY1,LZ1)
C
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /-1/
C
      IF (imsh.EQ.1) NEL = NELV
      IF (imsh.EQ.2) NEL = NELT
      NXY1  = NX1*NY1
      NYZ1  = NY1*NZ1
      NXYZ1 = NX1*NY1*NZ1
      NTOT  = NXYZ1*NEL
C
      if (istep.ne.icalld) then
         ntott = nxyz1*nelt
         call invers2(jacmi,jm1,ntott)
         icalld=istep
      endif
C
      DO 1000 IEL=1,NEL
C
      IF (IFAXIS) CALL SETAXDY (IFRZER(IEL) )
C
      IF (NDIM.EQ.2) THEN
            CALL MXM     (DXM1,NX1,U(1,1,1,IEL),NX1,DU(1,1,1,IEL),NYZ1)
            CALL COL2    (DU(1,1,1,IEL),RM1(1,1,1,IEL),NXYZ1)
            CALL MXM     (U(1,1,1,IEL),NX1,DYTM1,NY1,DRST,NY1)
            CALL ADDCOL3 (DU(1,1,1,IEL),DRST,SM1(1,1,1,IEL),NXYZ1)
      ELSE
            CALL MXM   (DXM1,NX1,U(1,1,1,IEL),NX1,DU(1,1,1,IEL),NYZ1)
            CALL COL2  (DU(1,1,1,IEL),RM1(1,1,1,IEL),NXYZ1)
            DO 20 IZ=1,NZ1
               CALL MXM  (U(1,1,IZ,IEL),NX1,DYTM1,NY1,DRST(1,1,IZ),NY1)
 20         CONTINUE
            CALL ADDCOL3 (DU(1,1,1,IEL),DRST,SM1(1,1,1,IEL),NXYZ1)
            CALL MXM     (U(1,1,1,IEL),NXY1,DZTM1,NZ1,DRST,NZ1)
            CALL ADDCOL3 (DU(1,1,1,IEL),DRST,TM1(1,1,1,IEL),NXYZ1)
      ENDIF
C
 1000 CONTINUE
C
c     CALL INVCOL2 (DU,JM1,NTOT)
      CALL COL2 (DU,JACMI,NTOT)
C
      return
      END
c-----------------------------------------------------------------------
