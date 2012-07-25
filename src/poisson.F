
c=======================================================================
      subroutine poissonq(myex,myey,myez,phi,n)
c=======================================================================
c     solve  -nabla^2 (potent) = charge/permittivity
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

c     call setqvol (charge)         ! rhs
      call copy    (rhs,charge,npts)
      call invcol2 (rhs,permittivity,npts)   ! B*rhs
      call col2    (rhs,bm1,npts)   ! B*rhs
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,npts)
      if (nid.eq.0) write(6,*) 'done: dssum rsh in poissonq'       

      call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call rzero   (h2,npts)
      if (nid.eq.0) write(6,*) 'done: set h1 and h2 in poissonq'       

      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = 1.e-18 ! 1.e-10

      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES

      if (IFCG) then
      call cggo2 (phi,rhs,h1,h2,mult,mask,imsh,isd,npts,maxit)
      endif
      if (IFGMRES) then
      call  hmh_gmres2(phi,rhs,h1,h2,mult,mask,isd,imsh,npts,tol)
      endif
 
      call maxwell_grad (myex,myey,myez,phi)
      if (nid.eq.0) write(6,*) 'done: maxwell_grad'
      if (nid.eq.0) write(6,*) 'done: poissonq'
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
      
      npts =  nx1*ny1*nz1*nelt

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

      if (nid.eq.0) write(6,*) 'call setmaskp'
      nface = 2*ndim

      do e=1,nelt
      do f=1,nface

         CB =CBC(f,e,2)

c        Homogeneous Dirichlet boundary conditions
         if (CB.eq.'PEC') call facev (mask,e,f,0.0,nx1,ny1,nz1)
         if (CB.eq.'PML') call facev (mask,e,f,0.0,nx1,ny1,nz1) ! fischer/mmin: 1/23/07

      enddo
      enddo

      call dsop(mask,'MUL',nx1,ny1,nz1)
      if (nid.eq.0) write(6,*) 'done: mask dsop'

      return
      end
c=======================================================================
      subroutine cggo2(x,f,h1,h2,mult,mask,imsh,isd,n,niter)
c=======================================================================
      include 'SIZE'
      common /srcuz/ r(lpts1),w(lpts1),p(lpts1),z(lpts1)

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
      integer n, niter,imsh,isd
      real    eps, rhsnorm
      real    x(n),f(n),h1(n),h2(n),mult(n),mask(n)

      logical ifmcor

      call rzero   (r,n)      ! init                                   
      call rzero   (w,n)      ! init                                   
      call rzero   (p,n)      ! init                                   
      call rzero   (z,n)      ! init      

c     set machine tolerances

      one = 1.
      eps = 1.e-20
      rhsnorm = sqrt(glsc2(f,f,n))
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7
      eps= eps*rhsnorm

      !tube-9: just remain above: no additional setting for eps

      if (nid.eq.0) write(6,*) 'start: cggo2'
      if (nid.eq.0) write(6,5) eps,n,niter
   5  format(' eps=',e15.7,', n=',i6, ', iter=',i6)

      pap = 0.0

      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)

      ifmcor= .false.
      skmin =  glmin(mask,n)
      if (skmin.gt.0) ifmcor = .true.

      if (ifmcor) then
          smean = -1./glsum(mult,n)
          rmean = smean*glsc2(r,mult,n)
          call cadd(r,rmean,n)
      endif

      rnorm= sqrt(glsc3(r,mult,r,n)) ! inner product

      iter = 0
      if (nid.eq.0) write(6,*) 'iter;rnorm',iter,rnorm,ifmcor

      miter = niter

      rtz1=1.0
      do 1000 iter=1,miter

         call solveM(z,r,n)         ! Invert preconditioner here

         if (ifmcor) then
             rmean = smean*glsc2(z,mult,n)
             call cadd(z,rmean,n)
         endif

         rtz2 = rtz1
         rtz1 = glsc3(r,mult,z,n)   ! parallel inner product
         beta = rtz1/rtz2

         if (iter.eq.1) beta=0.0

         call add2s1 (p,z,beta,n)    ! p=p*beta+z
           
         call axhelm2(w,p,h1,h2,imsh,isd)
         call dssum  (w,nx1,ny1,nz1)
         call col2   (w,mask,n)

         pap = glsc3 (w,mult,p,n)

         alpha=  rtz1/pap
         alphm= -alpha

         call add2s2(x,p,alpha,n)
         call add2s2(r,w,alphm,n)

         rtr = glsc2(r,r,n)

         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)

         if (nid.eq.0)  write(6,6) iter,rnorm,alpha,beta,pap
c        if (iter.le.10.or.mod(iter,10).eq.0) write(6,6) iter,rnorm
    6    format(' cg:',i8,1p4e12.4)

         if (rtr.le.rlim2) goto 1001

 1000 continue
 1001 continue
      if (nid.eq.0) write(6,*) 'done: cggo2'

      return
      end



c=======================================================================
      subroutine hmholt2(name,u,rhs,h1,h2,mask,mult,imsh,tli,maxit,isd)
c=======================================================================
      INCLUDE 'SIZE'
      INCLUDE 'CTIMER'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      include 'FDMH1'
      include 'TSTEP'

      CHARACTER      NAME*4
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           RHS  (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           MASK (LX1,LY1,LZ1,1)
      REAL           MULT (LX1,LY1,LZ1,1)

      tol = abs(tli)

      thmhz = 0.1
      etime1= dnekclock()
 
      ntot = nx1*ny1*nz1
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,ntot)

      call cgg2(u,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binvm1,name)
      thmhz=thmhz+(dnekclock()-etime1)

      return
      END

c=======================================================================
      subroutine axhelm2(au,u,helm1,helm2,imesh,isd)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
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
 
c     IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)

      CALL RZERO (AU,NTOT)
 
      DO 100 IEL=1,NEL
        ie=iel
 
c       IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
 
        IF (NDIM.EQ.2) THEN ! 2-d case ...............
 
c          IF (IFFAST(IEL)) THEN 
c          Fast 2-d mode: constant properties and undeformed element
c
c          H1 = HELM1(1,1,1,IEL)
c          CALL MXM   (WDDX,NX1,U(1,1,1,IEL),NX1,TM1,NYZ)
c          CALL MXM   (U(1,1,1,IEL),NX1,WDDYT,NY1,TM2,NY1)
c          CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ)
c          CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)
c          CALL ADD3  (AU(1,1,1,IEL),TM1,TM2,NXYZ)
c          CALL CMULT (AU(1,1,1,IEL),H1,NXYZ)
c
c          ELSE
c
c          General case, speed-up for undeformed elements
c
           CALL MXM  (DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           CALL MXM  (U(1,1,1,IEL),NX1,DYTM1,NY1,DUDS,NY1)
           CALL COL3 (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3 (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)

c          IF (IFDFRM(IEL)) THEN
c             CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
c             CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
c          ENDIF

           CALL COL2 (TMP1,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,HELM1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           CALL MXM  (TMP2,NX1,DYM1,NY1,TM2,NY1)
           CALL ADD2 (AU(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM2,NXYZ)
 
c          ENDIF
 
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

c          IF (IFDFRM(IEL)) THEN
c             CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
c             CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
c             CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
c             CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
c             CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
c             CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
c          ENDIF

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
 
      return  !FIXME  temporarily use general case: misun 7/16/2012 


      IF (IFH2) CALL ADDCOL4 (AU,HELM2,BM1,U,NTOT)
 
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
               IF (YM1(I,J,1,IEL).NE.0.) THEN
                  TERM1 = BM1(I,J,1,IEL)*U(I,J,1,IEL)/YM1(I,J,1,IEL)**2
                  IF (IFRZER(IEL)) THEN
                     TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J)*DUAX(I)
     $                       *JACM1(I,1,1,IEL)/YSM1(I)
                  ELSE
                     TERM2 = 0.
                  ENDIF
                  AU(I,J,1,IEL) = AU(I,J,1,IEL)
     $                          + HELM1(I,J,1,IEL)*(TERM1+TERM2)
               ENDIF
  190       CONTINUE
  200    CONTINUE
      ENDIF
      if (nid.eq.0) write(6,*) 'done: axhelm2' 
C
#ifndef NOTIMER
      taxhm=taxhm+(dnekclock()-etime1)
#endif
      return
      END
C
c=======================================================================
      subroutine cgg2(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name)
C-----------------------------------------------------------------------
C
C     Solve the Helmholtz equation, H*U = RHS,
C     using preconditioned conjugate gradient iteration.
C     Preconditioner: diag(H).
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      include 'FDMH1'
      include 'GEOM'
c
      COMMON  /CPRINT/ IFPRINT, IFHZPC
      LOGICAL          IFPRINT, IFHZPC
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      logical ifmcor,ifprint_hmh
C
      real x(1),f(1),h1(1),h2(1),mask(1),mult(1),binv(1)
      parameter        (lg=lx1*ly1*lz1*lelt)
      COMMON /SCRCG/ d (lg) , scalar(2)
      common /SCRMG/ r (lg) , w (lg) , p (lg) , z (lg)
c
      parameter (maxcg=900)
      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /iterhm/ niterhm
      character*4 name
c
      if (ifsplit.and.name.eq.'PRES'.and.param(42).eq.0) then
         n = nx1*ny1*nz1*nelv
         call copy      (x,f,n)
         call hmh_gmres (x,h1,h2,mult,iter)
         niterhm = iter
         return
      endif
c      write(6,*) ifsplit,name,param(44),' P44 C'

c **  zero out stuff for Lanczos eigenvalue estimator
      call rzero(diagt,maxcg)
      call rzero(upper,maxcg)
C
C     Initialization
C
      NXYZ   = NX1*NY1*NZ1
      NEL    = NELV
      VOL    = VOLVM1
      IF (IMSH.EQ.2) NEL=NELT
      IF (IMSH.EQ.2) VOL=VOLTM1
      n      = NEL*NXYZ
c
      tol=abs(tin)
      if (param(22).ne.0) tol=abs(param(22))
      if (name.eq.'PRES'.and.param(21).ne.0) tol=abs(param(21))
      if (tin.lt.0)       tol=abs(tin)
      niter = min(maxit,maxcg)

C     Speed-up for undeformed elements and constant properties.
      if (.not.ifsolv) then
         call setfast(h1,h2,imesh)
         ifsolv = .true.
      endif
C
C     Set up diag preconditioner.
C
      if (kfldfdm.lt.0) then
         call setprec(D,h1,h2,imsh,isd)
      elseif(param(100).ne.2) then
         call set_fdm_prec_h1b(d,h1,h2,nel)
      endif
c
      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)
c
c     Check for non-trivial null-space
c
      ifmcor = .false.
      h2max = glmax(h2  ,n)
      skmin = glmin(mask,n)
      if (skmin.gt.0.and.h2max.eq.0) ifmcor = .true.
C
      if (name.eq.'PRES') then
c        call ortho (r)           ! Commented out March 15, 2011,pff
      elseif (ifmcor) then
         smean = -1./glsum(mult,n)
         rmean = smean*glsc2(r,mult,n)
         call cadd(r,rmean,n)
      endif
C
      krylov = 0
      rtz1=1.0
      niterhm = 0

      do iter=1,niter
C
         if (kfldfdm.lt.0) then  ! Jacobi Preconditioner
c           call copy(z,r,n)
            call col3(z,r,d,n)
         else                                       ! Schwarz Preconditioner
            if (name.eq.'PRES'.and.param(100).eq.2) then
               call h1_overlap_2(z,r,mask)
               call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
               call add2         (z,w,n)
            else   
               call fdm_h1(z,r,d,mask,mult,nel,ktype(1,1,kfldfdm),w)
               if (name.eq.'PRES') then 
                 call crs_solve_h1 (w,r)  ! Currently, crs grd only for P
                 call add2         (z,w,n)
               endif
            endif
         endif
c
         if (name.eq.'PRES') then
            call ortho (z)
         elseif (ifmcor) then
            rmean = smean*glsc2(z,mult,n)
            call cadd(z,rmean,n)
         endif
c
         rtz2=rtz1
         scalar(1)=vlsc3 (z,r,mult,n)
         scalar(2)=vlsc32(r,mult,binv,n)
         call gop(scalar,w,'+  ',2)
         rtz1=scalar(1)
         rbn2=sqrt(scalar(2)/vol)
         if (iter.eq.1) rbn0 = rbn2
         if (param(22).lt.0) tol=abs(param(22))*rbn0
         if (tin.lt.0)       tol=abs(tin)*rbn0

         ifprint_hmh = .false.
         if (nid.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
         if (nid.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.

         if (ifprint_hmh)
     $      write(6,3002) istep,iter,name,ifmcor,rbn2,tol,h1(1),h2(1)


c        Always take at least one iteration   (for projection) pff 11/23/98
#ifndef TST_WSCAL
         IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
#else
         iter_max = param(150)
         if (name.eq.'PRES') iter_max = param(151)
         if (iter.gt.iter_max) then
#endif
c        IF (rbn2.LE.TOL) THEN
            NITER = ITER-1
c           IF(NID.EQ.0.AND.((.NOT.IFHZPC).OR.IFPRINT))
            if (nid.eq.0)
     $         write(6,3000) istep,name,niter,rbn2,rbn0,tol
            goto 9999
         ENDIF
c
         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
         call add2s1 (p,z,beta,n)
         call axhelm2(w,p,h1,h2,imsh,isd)
         call dssum  (w,nx1,ny1,nz1)
         call col2   (w,mask,n)
c
         rho0 = rho
         rho  = glsc3(w,p,mult,n)
         alpha=rtz1/rho
         alphm=-alpha
         call add2s2(x,p ,alpha,n)
         call add2s2(r,w ,alphm,n)
c
c        Generate tridiagonal matrix for Lanczos scheme
         if (iter.eq.1) then
            krylov = krylov+1
            diagt(iter) = rho/rtz1
         elseif (iter.le.maxcg) then
            krylov = krylov+1
            diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
            upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
         endif
 1000 enddo
      niter = iter-1
c
      if (nid.eq.0) write (6,3001) istep,niter,name,rbn2,rbn0,tol
 3000 format(4x,i7,4x,'Hmholtz ',a4,': ',I6,1p6E13.4)
 3001 format(2i6,' **ERROR**: Failed in HMHOLTZ: ',a4,1p6E13.4)
 3002 format(i3,i6,' Helmholtz ',a4,1x,l4,':',1p6E13.4)
 9999 continue
      niterhm = niter
      ifsolv = .false.
c
c
c     Call eigenvalue routine for Lanczos scheme:
c          two work arrays are req'd if you want to save "diag & upper"
c
c     if (iter.ge.3) then
c        niter = iter-1
c        call calc (diagt,upper,w,z,krylov,dmax,dmin)
c        cond = dmax/dmin
c        if (nid.eq.0) write(6,6) istep,cond,dmin,dmax,' lambda'
c     endif
c   6 format(i9,1p3e12.4,4x,a7)
c
c     if (n.gt.0) write(6,*) 'quit in cggo'
c     if (n.gt.0) call exitt
c     call exitt
      return
      end
c=======================================================================
      subroutine solveM(z,r,n)
c=======================================================================
      real z(n),r(n)

      call copy(z,r,n)

      return
      end
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

      if (nid.eq.0) write(6,*) 'start: hmh_gmres2'

      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

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

            call axhelm2(w,v(1,j),h1,h2,imsh,isd)  ! w = A v
            call dssum  (w,nx1,ny1,nz1)
            call col2   (w,mask,n)

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i

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

            rnorm = abs(gamma(j+1))!*norm_fac
c            ratio = rnorm/div0
            if (ifprint.and.nid.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,istep
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

c      call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,tolpss,etime_p,etime1     
 9999 format(i9,' GMRES:',i5,1p3e12.4)

      return
      end
c=======================================================================

