
c=======================================================================
      subroutine poissonq(myex,myey,myez,phi,n)
c=======================================================================
c     solve  nabla^2 (potent) = -charge/permittivity
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'

      common /ccpu/  cpu_t,cpu_dtime,cpu_chk

      integer  n
      real     myex(n)
      real     myey(n)
      real     myez(n)
      real     phi (n)
       
      cpu_dtime= dclock()          
      cpu_t    = 0.0

      if (nid.eq.0) write(6,*) '\n call poissonq'       
      npts =  nx1*ny1*nz1*nelt
      call rone    (mask,npts)      ! mask is one everywhere                 
      call setmaskp(mask,npts)      ! mask is zero for PEC, PML 

      if (nid.eq.0) write(6,*) 'done: setmaskp'       

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) 'done: dssum mult in poissonq'       

      call copy    (rhs,charge,npts)          ! rhs=charge
      call invcol2 (rhs,permittivity,npts)    ! rhs=rhs/permittivity
      call col2    (rhs,bm1,npts)             ! rhs=rhs*bm1       
      call dssum   (rhs,nx1,ny1,nz1)          ! rhs=dssum(rhs)
      call col2    (rhs,mask,npts)            ! rhs=0 for dirichlet
      if (nid.eq.0) write(6,*) 'done:: dssum rsh in poissonq'       

      call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call rzero   (h2,npts)
      if (nid.eq.0) write(6,*) 'done:: set h1 and h2 in poissonq'       

      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = param(22)
     
      if (nid.eq.0) write(6,*) 'start IFCG/IFGMRES',IFCG,IFGMRES

      if (IFCG)    then
         call cggo2 (phi,rhs,h1,h2,mult,mask,imsh,isd,npts,tol,maxit)
      endif
      if (IFGMRES) then
         call hmh_gmres2(phi,rhs,h1,h2,mult,mask,isd,imsh,npts,tol)
      endif
 
      if (nid.eq.0) write(6,*) 'done:: IFCG/IFGMRES',IFCG,IFGMRES
      if (nid.eq.0) write(6,*) 'done:: poissonq'
      cpu_dtime= dclock()-cpu_dtime
      cpu_t= cpu_t+cpu_dtime

      return
      end
c=======================================================================
      subroutine maxwell_init_poissonq(eix,eiy,eiz)
c=======================================================================
c     solve  nabla^2 phi = -q/epsilon,  q from userq()
c     set ex,ey,ez = -grad (phi)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      
      call rone    (mask,npts)           ! mask is one everywhere                 
      call setmaskp(mask,npts)           ! mask is zero for PEC, PML 

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)

c     call setqvol (charge)                ! rhs
      call copy    (rhs,charge,npts)       ! rhs= charge
      call col2    (rhs,bm1,npts)          ! rhs= B*rhs
      call invcol2 (rhs,permittivity,npts) ! rhs= rhs/permittivity
      call chsign  (rhs,npts)              ! rhs= -rhs
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,npts)

      call rone    (h1,npts)             ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call rzero   (h2,npts)

      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = 1.e-18 ! 1.e-10

c     call cggo 
c    $    (phi,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binv,'init')

c     call cgg2 (phi,rhs,mult,mask,r,w,p,z,imsh,isd,n,maxit)
c     if (nid.eq.0) write(6,*) 'POISSON Solver Done!'
c
c     call maxwell_grad(EN(1,1),EN(1,2),EN(1,3),potent)
      call chsign  (EN(1,1),npts)
      call chsign  (EN(1,2),npts)
      call chsign  (EN(1,3),npts)

      !call outpost (eix,eiy,eiz,p,phi,'   ')

      return
      end
c=======================================================================
      subroutine setmaskp(mask,n)
c=======================================================================
C     set mask for poisson problem
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer   e,f,i,j,n
      real      mask(n)
      character CB*3

      do e=1,nelt
      do f=1,nfaces

         CB =CBC(f,e,1)
c        Homogeneous Dirichlet boundary conditions
         if (CB.eq.'PEC') call facev (mask,e,f,0.0,nx1,ny1,nz1)
         if (CB.eq.'PML') call facev (mask,e,f,0.0,nx1,ny1,nz1) ! fischer/mmin: 1/23/07

      enddo
      enddo

      call dsop(mask,'MUL',nx1,ny1,nz1)

c     if (nid.eq.0) write(6,*) ' done: setmaskp'

      return
      end
c=======================================================================
      subroutine cggo2(x,f,h1,h2,mult,mask,imsh,isd,n,tol,niter)
c=======================================================================
      include 'SIZE'
      include 'GEOM'
      include 'MASS'
      include 'TSTEP'
      include 'INPUT'
      include 'FDMH1'
      common  /srcuz/ d(lpts),r(lpts),w(lpts),p(lpts),z(lpts)
      common  /srccg/ scalar(2)
      real     tol,r,w,p,z,scalar
      logical ifio

c     Solve Ax=f where A is SPD and is invoked by the routine ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided routine ax(w,z,n) returns  w := Az,  
c
c     User-provided routine solveM(z,r,n) ) returns  z := M^-1 r,  
c
c     User-provided array wght() is used to scale inner-products
c     of the form (p,z) = p'*wght*z
c
      integer  n,niter,imsh,isd
      real     eps, rhsnorm
      real     x(n),f(n),h1(n),h2(n),mult(n),mask(n)
      logical  ifmcor

      call rzero  (r,n)      ! init                                   
      call rzero  (w,n)      ! init                                   
      call rzero  (p,n)      ! init                                   
      call rzero  (z,n)      ! init      

c     set machine tolerances
C     one = 1.
C     eps = 1.e-20
C     rhsnorm = sqrt(glsc2(f,f,n))
C     if (one+eps .eq. one) eps = 1.e-14
C     if (one+eps .eq. one) eps = 1.e-7
C     eps= eps*rhsnorm


c     if (nid.eq.0) write(6,*) 'start: cggo2'
c     if (nid.eq.0) write(6,5) tol,n,niter
c  5  format(' tol=',e15.7,', n=',i6, ', iter=',i6)

      pap = 0.0
      vol = volvm1
      
      call rone   (bminv,n)
      call invcol2(bminv,bmn,n)

      kfldfdm=1

      if (param(23).eq.1) call setprec(D,h1,h2,imsh,isd)       !IFFDM=true
      if (param(23).eq.2) call set_fdm_prec_h1b(d,h1,h2,nel)          

      call copy   (r,f,n) ! f=b-Ax_0  
      call rzero  (x,n)
      call rzero  (p,n)

c...  check non-trivial null-space: "ifmcor := if mean correction is true"
c...  choose the non-trivial solution that has zero mean; this is equivalent to
c...  orthogonalization to one vector (1,....,1) of length n
c...
c...  eg), A= [ 2 0 ]  u=[1/2,0] or u[1/2,10] for Au=0
c...          [ 0 0 ]
c...       Our routine will give result of u=[1/2,1/2] 

      ifmcor= .false.
      h2max = glmax(h2  ,n)
      skmin = glmin(mask,n)                ! check dirichlet boundary
      if ((skmin.gt.0).and.(h2max.eq.0)) ifmcor = .true.  ! no dirichlet, no h2 term

      if (ifmcor) then
          smean = -1./glsum(bmn,n)         ! glsum(bmn)= volume integration    
          rmean = smean*glsc2(r,mult,n)    !    
          call copy  (x,bmn,n) 
          call dssum (x,nx1,ny1,nz1) 
          call add2s2(r,x,rmean,n)         ! r = r - rmean*x
          call rzero (x,n)
      endif

      rtz1= 1.0
      do iter=1,niter

         if     (param(23).eq.1) then
             call col3(z,r,d,n)
         elseif (param(23).eq.2) then
             call fdm_h1(z,r,d,mask,mult,nel,ktype(1,1,kfldfdm),w)
         else
             call copy(z,r,n)
         endif

         if (ifmcor) then
             rmean = smean*glsc2(z,bmn,n)  ! inner product could give different numbers on different # cores.   
             call cadd(z,rmean,n)          ! if do "sort" and then "glsc2" might give the same value all the time.
         endif

         rtz2      = rtz1
         scalar(1) = vlsc3 (r,mult,z,n)    ! tmp=tmp+r(i)*mult(i)*z(i)                      
         scalar(2) = vlsc32(r,mult,bminv,n)! tmp=tmp+r(i)*mult(i)*r(i)/bmn(i)
         call gop(scalar,w,'+  ',2)

         rtz1 = scalar(1)                   
         rbn2 = sqrt(scalar(2)/vol)        ! rbn2= norm(r,2) 

         if (iter.eq.1) rbn0 = rbn2
              
         if ((nid.eq.0).and.(istep.le.5))      
     $           write(6,1001) iter-1,ifmcor,rbn2,rbn0,tol,istep
         if (rbn2.le.tol) go to 1000

         beta = rtz1/rtz2

         if  (iter.eq.1) beta=0.0
         call add2s1  (p,z,beta,n)         ! p = z + beta*p
         call axhelm2 (w,p,h1,h2,imsh,isd) ! w = Ap
         call dssum   (w,nx1,ny1,nz1)      !
         call col2    (w,mask,n)           !

         pap  = glsc3 (w,mult,p,n)        ! pap  = (p^T)Ap
         alpha= rtz1/pap                  ! alpha= (r^T r)/pap
         alphm= -alpha

         call add2s2 (x,p,alpha,n)         ! x = x + alpha*p
         call add2s2 (r,w,alphm,n)         ! r = r + alphm*(Ap)

      enddo   
 1000 continue

      if (nid.eq.0) then
         if ((mod(istep,iocomm).eq.0).or.(istep.le.10))      
     $    write(6,1002) istep,iter,ifmcor,rbn2,rbn0,tol
      endif 

 1001 format(i5,l4,2x,1p3e12.4,i8,' cggo2')
 1002 format(i9,' cggo2: iteration#',i8,l4,2x,1p3e12.4)

      return
      end
c=======================================================================
      subroutine axhelm2(au,u,helm1,helm2,imsh,isd)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C     This routine solves  "A= D^T (B*h1) D" 
C
C------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
C
      COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)

      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)

      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)
C
      NEL = NELT
      NXY = NX1*NY1
      NYZ = NY1*NZ1
      NXZ = NX1*NZ1
      NXYZ= NX1*NY1*NZ1
      NTOT= NXYZ*NEL
C
#ifndef NOTIMER
      if (icalld.eq.0) taxhm=0.0
      icalld=icalld+1
      naxhm=icalld
      etime1=dnekclock()
#endif
 
      IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)

      CALL RZERO (AU,NTOT)
 
      DO 100 IEL=1,NEL
 
c       IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
 
        IF (NDIM.EQ.2) THEN ! 2-d case ...............
 
C          IF (IFFAST(IEL)) THEN 
c          Fast 2-d mode: constant properties and undeformed element
 
C          H1 = HELM1(1,1,1,IEL)
C          CALL MXM   (WDDX,NX1,U(1,1,1,IEL),NX1,TM1,NYZ)
C          CALL MXM   (U(1,1,1,IEL),NX1,WDDYT,NY1,TM2,NY1)
C          CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ)
C          CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)
C          CALL ADD3  (AU(1,1,1,IEL),TM1,TM2,NXYZ)
C          CALL CMULT (AU(1,1,1,IEL),H1,NXYZ)
C
C          ELSE
c
c          General case, speed-up for undeformed elements
c
           CALL MXM  (DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           CALL MXM  (U(1,1,1,IEL),NX1,DYTM1,NY1,DUDS,NY1)
           CALL COL3 (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3 (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)

C          IF (IFDFRM(IEL)) THEN
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
C          ENDIF

           CALL COL2 (TMP1,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,HELM1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           CALL MXM  (TMP2,NX1,DYM1,NY1,TM2,NY1)
           CALL ADD2 (AU(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM2,NXYZ)
 
          ! ENDIF
 
        ELSE

c          3-d case ...............
c          IF (IFFAST(IEL)) THEN
 
c          Fast 3-d mode: constant properties and undeformed element
 
c          H1 = HELM1(1,1,1,IEL)
c          CALL MXM   (WDDX,NX1,U(1,1,1,IEL),NX1,TM1,NYZ)
c          DO 5 IZ=1,NZ1
c          CALL MXM   (U(1,1,IZ,IEL),NX1,WDDYT,NY1,TM2(1,1,IZ),NY1)
c5         CONTINUE
c          CALL MXM   (U(1,1,1,IEL),NXY,WDDZT,NZ1,TM3,NZ1)
c          CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ)
c          CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)
c          CALL COL2  (TM3,G6M1(1,1,1,IEL),NXYZ)
c          CALL ADD3  (AU(1,1,1,IEL),TM1,TM2,NXYZ)
c          CALL ADD2  (AU(1,1,1,IEL),TM3,NXYZ)
c          CALL CMULT (AU(1,1,1,IEL),H1,NXYZ)
c
c          ELSE
 
c          General case, speed-up for undeformed elements

           CALL MXM(DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           DO 10 IZ=1,NZ1
              CALL MXM(U(1,1,IZ,IEL),NX1,DYTM1,NY1,DUDS(1,1,IZ),NY1)
   10      CONTINUE
           CALL MXM     (U(1,1,1,IEL),NXY,DZTM1,NZ1,DUDT,NZ1)
           CALL COL3    (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3    (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)
           CALL COL3    (TMP3,DUDT,G3M1(1,1,1,IEL),NXYZ)

C          IF (IFDFRM(IEL)) THEN
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
C          ENDIF

           CALL COL2 (TMP1,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP3,HELM1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           DO 20 IZ=1,NZ1
              CALL MXM(TMP2(1,1,IZ),NX1,DYM1,NY1,TM2(1,1,IZ),NY1)
   20      CONTINUE
           CALL MXM  (TMP3,NXY,DZM1,NZ1,TM3,NZ1)
           CALL ADD2 (AU(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM2,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM3,NXYZ)
 
c          ENDIF
 
        ENDIF
 
 100  CONTINUE
 
      CALL ADDCOL4 (AU,HELM2,BM1,U,NTOT) ! always called here: sacrify efficiency
 
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
      IF (IFAXIS.AND.(ISD.EQ.2)) THEN
         DO 200 IEL=1,NEL
C
            IF (IFRZER(IEL)) THEN
               CALL MXM(U  (1,1,1,IEL),NX1,DATM1,NY1,DUAX,1)
               CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
            ENDIF
C
            DO 190 J=1,NY1
            DO 190 I=1,NX1
C              IF (YM1(I,J,1,IEL).NE.0.) THEN
                  TERM1 = BM1(I,J,1,IEL)*U(I,J,1,IEL)/YM1(I,J,1,IEL)**2
                  IF (IFRZER(IEL)) THEN
                     TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J)*DUAX(I)
     $                       *JACM1(I,1,1,IEL)/YSM1(I)
                  ELSE
                     TERM2 = 0.
                  ENDIF
                  AU(I,J,1,IEL) = AU(I,J,1,IEL)
     $                          + HELM1(I,J,1,IEL)*(TERM1+TERM2)
C              ENDIF
  190       CONTINUE
  200    CONTINUE
      ENDIF
C
#ifndef NOTIMER
      taxhm=taxhm+(dnekclock()-etime1)
#endif
      return
      END
C
c=======================================================================
      subroutine hmh_gmres2(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,outer,isd,imsh
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp
      real*8   etime1,dnekclock

        common /ctmp0/   wk1(lgmres),wk2(lgmres)

      if (nid.eq.0) write(6,*) 'start: hmh_gmres2'

      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p= 0.
      divex  = 0.
      iter   = 0
      m      = lgmres

      !norm_fac = 1./sqrt(volvm1)

      tolps = tol
      tolpss= tolps
      iconv = 0
      call rzero(x,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r,res,n)                  ! r = res
         else
            !update residual
            call copy   (r,res,n)                  ! r = res
            call axhelm2(w,x,h1,h2,imsh,isd)       ! w = A x
            call dssum  (w,nx1,ny1,nz1)
            call col2   (w,mask,n)
            call add2s2 (r,w,-1.,n)                ! r = r - w
         endif
         
         gamma(1) = glsc3(r,r,wt,n)                ! gamma  = (r,r) 
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            !call ortho_new  (v(1,j)) 
            call axhelm2(w,v(1,j),h1,h2,imsh,isd)  ! w = A v
            call dssum  (w,nx1,ny1,nz1)
            call col2   (w,mask,n)

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
            enddo                                 !         i,j  i
                                                  ! i,j       i
            call gop(h(1,j),wk1,'+  ',j) 

            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1)) !*norm_fac
c            ratio = rnorm/div0
            if (nid.eq.0) write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' Divergence')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1            
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,v(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x,n)
      call copy(res,x,n)
c     call ortho_new   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1

      if (nid.eq.0) write(6,9999) istep,iter,tolpss,etime_p,etime1     
 9999 format(i9,' GMRES:',i5,1p3e12.4)

      return
      end
c=======================================================================
      subroutine ortho_new (respr)
      implicit none
C     Orthogonalize the residual in the pressure solver with respect
C     to (1,1,...,1)T  (only if all Dirichlet b.c. ).

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'TSTEP'
      real      respr(lx1,ly1,lz1,lelv)
      integer*8 ntot,ntotg,nxyz1
      real      rlam, glsum

      nxyz1 = nx1*ny1*nz1
      ntot  = nxyz1*nelv
      ntotg = nxyz1*nelgv

      rlam  = glsum (respr,ntot)/ntotg
      call cadd (respr,-rlam,ntot)

      return
      end
c=======================================================================

