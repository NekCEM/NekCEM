c=======================================================================
      subroutine hmh2DtN
c=======================================================================
c     solve  alpha nabla^2 (u) +k^2 u=f
c     with periodic in X direction
c          u(x,0)=0 and DtN du/dy(x,1)-\beta u(x,1)=g(x)
c     remark: sign(alpha)*sign(h1)=-1
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer  icalld
      save     icalld
      data     icalld /0/

      common /ccpu/  cpu_t,cpu_dtime,cpu_chk

      cpu_dtime= dclock()          
      cpu_t    = 0.0

      if (nid.eq.0) write(6,*) '\n call hmh2DtN'       
      npts2 =  npts*2
      if (nid.eq.0) write(6,*) 'before: setmaskp'       
      call rzero   (rhs ,npts2)      ! init                                   
      if (nid.eq.0) write(6,*) 'before: setmaskp'       
      call rone    (mask0,npts )      ! mask0 is one everywhere                 
      if (nid.eq.0) write(6,*) 'before: setmaskp'       
      call setmaskp(mask0,npts )      ! mask0 is zero for PEC, PML 
      if (nid.eq.0) write(6,*) 'done: setmaskp'       
      do i=1,npts
         mask(i)     =mask0(i)
         mask(npts+i)=mask0(i)
      enddo

      call rone    (mult0,npts)
      call dssum   (mult0,nx1,ny1,nz1)
      call invcol1 (mult0,npts)
      do i=1,npts
         mult(i)     =mult0(i)
         mult(npts+i)=mult0(i)
      enddo
      if (nid.eq.0) write(6,*) 'done: dssum mult in poissonq'       
      call rone    (h1,npts2)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call chsign  (h1,npts2)        ! negative for helmholtz, positive for poisson
      call rone    (h2,npts2)
      do i=1,npts2
         h2(i)=h2(i)*kwave(i)**2
      enddo
c      call rzero    (h2,npts2)     ! zero for poisson problem
      if (nid.eq.0) write(6,*) 'done: set h1 and h2 in acoustic'       
ccccccc below for force term
c      if (IFSRC) then
c      call copy    (rhs,charge,npts2)
c      call invcol2 (rhs,permittivity,npts)   ! B*rhs
c      call col2    (rhs,bm1,npts)   ! B*rhs
c      call col2    (rhs(npts+1),bm1,npts)   ! B*rhs
c      if (nid.eq.0) write(6,*) 'done: adding source term in acoustic'       
c      endif
      IF (IFDTN) THEN
ccccc RHS for DTN uinc
      maxdtn=0
      ncemface_dtn=0
      if (icalld.eq.0) call cem_set_dtn_fc_ptr
      icalld=1
      call cem_acoustic_dtn
      call copy(mult2,mult,npts2)
      call invcol1 (mult2,npts2)
      if (nid.eq.0) write(6,*) 'done:  percomputation for DtN'       
      call cem_acoustic_bd_dtn(rhs,charge)
      call rzero(charge,npts2)
      call cem_acoustic_bd(charge(1),potent(1))
      call cem_acoustic_bd(charge(npts+1),potent(npts+1))
      call add2s2  (rhs,charge,-1.,npts2)                ! rhs = rhs - charge
      endif
ccccccc below for neumann boundary condition
      if (IFNEU) THEN
      call rzero(charge,npts2)
      call cem_acoustic_bd(charge(1),potent(1))
      call cem_acoustic_bd(charge(npts+1),potent(npts+1))
      call add2s2  (rhs,charge,-1.,npts2)                ! rhs = rhs + charge
      if (nid.eq.0) write(6,*) 'done: dssum boundary in acoustic'
      endif
      
      call dssum   (rhs(1),     nx1,ny1,nz1)   
      call dssum   (rhs(npts+1),nx1,ny1,nz1)   
      call col2    (rhs,mask,npts2)   
      if (nid.eq.0) write(6,*) 'done: dssum rsh in acoustic'       

      isd   = 1
      imsh  = 1
      tol   = 1.e-10 ! 1.e-10
      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES
      if (IFGMRES) then
      call copy(potent,spotent,npts2)!passing initial guess other than zero if necessary
      call  hmh_gmres3(potent,rhs,h1,h2,mult,mask,isd,imsh,npts2,tol)
c      call  hmh_gmres4(potent,rhs,h1,h2,mult,mask,isd,imsh,npts2,tol)
      else 
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
      endif
      if (nid.eq.0) write(6,*) 'done: cgg2DtN'
      
      call copy(charge,potent,npts2)
      call add2s2  (charge,en(1,1),-1.,npts)                ! rhs = rhs + charge
      call add2s2  (charge(npts+1),en(1,2),-1.,npts)                ! rhs = rhs + charge
      call maxwell_grad (en(1,1),en(1,2),en(1,3),charge(1))
      call maxwell_grad (hn(1,1),hn(1,2),hn(1,3),charge(npts+1))
      if (nid.eq.0) write(6,*) 'done: maxwell_grad'
      if (nid.eq.0) write(6,*) 'done: acoustic'
      call col2   (hn(1,1),charge,npts)        ! dui/dy=dui/dy*ur
      call col2   (en(1,1),charge(npts+1),npts)!dur/dy=dur/dy*ui
      call add2s2  (hn(1,1),en(1,1),-1.,npts)                !  =dui/dy-dur/dy 
      call col2   (hn(1,1),mask,npts)
      
      call cem_acoustic_test(ssum,hn(1,1))
 
      if (nid.eq.0) write(6,*) 'test-energy error='
     $                        ,(ssum-kwave(1))/kwave(1)
      
      cpu_dtime= dclock()-cpu_dtime
      cpu_t= cpu_t+cpu_dtime

      return
      end
c=======================================================================
      subroutine hmh_gmres4(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,outer,isd,imsh
      real     x2(n),r2(n),wt2(n),w2(n),v2(n,lgmres)
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp
      real*8   etime1,dnekclock

      if (nid.eq.0) write(6,*) 'start: hmh_gmres3'

      n1 = npts
      n2 = npts*2

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      call rzero(h,m*m)
      call rzero(x2,n)!zero inital guess
c/block below if for zero inital guess
c      call copy(x2,phi,n) ! none zero inital guess
c      iter  = 1
c      call dssum(x2,nx1,ny1,nz1)
c      call dssum(x2(n1+1),nx1,ny1,nz1)
c      call col2 (x2,mask,n)
c/block above
      outer = 0
      call copy(wt2,wt,n)
      call rone(wt,n)
      do while (iconv.eq.0.and.iter.lt.3000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r2,res,n)                  ! r = res
            call col2  (r2,wt2,n)
         else
            !update residual
            call copy   (r2,res,n)                  ! r = res
            call col2   (r2,wt2,n)
            call axhelm_dtn(w2,x2,h1,h2,imsh,isd)       ! w = A x
            call dssum  (w2(1),nx1,ny1,nz1)
            call dssum  (w2(n1+1),nx1,ny1,nz1)
            call col2   (w2,wt2,n)
            call col2   (w2,mask,n)
            call add2s2 (r2,w2,-1.,n)                ! r = r - w
         endif
         
         gamma(1) = glsc3(r2,r2,wt,n)                ! gamma  = (r,r) 
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1

         !check for lucky convergence
c         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v2(1,1),r2,temp,n)             !  v  = r / gamma
                                                    !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call axhelm_dtn(w2,v2(1,j),h1,h2,imsh,isd)  ! w = A v
            call dssum  (w2(1),nx1,ny1,nz1)
            call dssum  (w2(n1+1),nx1,ny1,nz1)
            call col2   (w2,wt2,n)
            call col2   (w2,mask,n)

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w2,v2(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i

               call add2s2(w2,v2(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = glsc3(w2,w2,wt,n)       ! alpha =  \/ (w,w)
            alpha = sqrt(alpha)        ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))!*norm_fac
c      if (nid.eq.0) write(6,*) iter,'gmres-rnorm',rnorm
c            ratio = rnorm/div0
            if (nid.eq.0)write(6,66) iter,tolpss,rnorm,istep 
   66       format(i5,1p2e12.5,i8,' Divergence')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v2(1,j+1),w2,temp,n)   ! v    = w / alpha
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
            call add2s2(x2,v2(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

       call copy(phi,x2,n)

c      call ortho   (res) ! Orthogonalize wrt null space, if present
      if (nid.eq.0) write(6,*) 'gmres-rnorm',rnorm
      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,tolpss,etime_p,etime1     
 9999 format(i9,' GMRES:',i5,1p3e12.4)

      return
      end
c=======================================================================
      subroutine hmh_gmres3(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,outer,isd,imsh
      real     x2(n),r2(n),w2(n),v2(n,lgmres)
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp
      real*8   etime1,dnekclock

      if (nid.eq.0) write(6,*) 'start: hmh_gmres3'

      n1 = npts
      n2 = npts*2

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      call rzero(h,m*m)
      call rzero(x2,n)!zero inital guess
c/block below if for zero inital guess
c      call copy(x2,phi,n) ! none zero inital guess
      call rone (x2,n)!zero inital guess
      iter  = 1
      call dssum(x2,nx1,ny1,nz1)
      call dssum(x2(n1+1),nx1,ny1,nz1)
      call col2 (x2,mask,n)
c/block above
      outer = 0
      do while (iconv.eq.0.and.iter.lt.3000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r2,res,n)                  ! r = res
         else
            !update residual
            call copy   (r2,res,n)                  ! r = res
            call axhelm_dtn(w2,x2,h1,h2,imsh,isd)       ! w = A x
            call dssum  (w2(1),nx1,ny1,nz1)
            call dssum  (w2(n1+1),nx1,ny1,nz1)
            call col2   (w2,mask,n)
            call add2s2 (r2,w2,-1.,n)                ! r = r - w
         endif
         
         gamma(1) = glsc3(r2,r2,wt,n)                ! gamma  = (r,r) 
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1

         !check for lucky convergence
c         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v2(1,1),r2,temp,n)             !  v  = r / gamma
                                                    !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call axhelm_dtn(w2,v2(1,j),h1,h2,imsh,isd)  ! w = A v
            call dssum  (w2(1),nx1,ny1,nz1)
            call dssum  (w2(n1+1),nx1,ny1,nz1)
            call col2   (w2,mask,n)

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w2,v2(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i

               call add2s2(w2,v2(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = glsc3(w2,w2,wt,n)       ! alpha =  \/ (w,w)
            alpha = sqrt(alpha)        ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))!*norm_fac
c      if (nid.eq.0) write(6,*) iter,'gmres-rnorm',rnorm
c            ratio = rnorm/div0
            if (nid.eq.0)write(6,66) iter,tolpss,rnorm,istep 
   66       format(i5,1p2e12.5,i8,' Divergence')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v2(1,j+1),w2,temp,n)   ! v    = w / alpha
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
            call add2s2(x2,v2(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x2,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present
      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,tolpss,etime_p,etime1     
 9999 format(i9,' GMRES:',i5,1p3e12.4)

      return
      end
c=======================================================================
      subroutine axhelm_dtn(au,u,helm1,helm2,imesh,isd)
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      include 'EMWAVE'
      INCLUDE 'ACOUSTIC'
C
      real au(1),u(1),helm1(1),helm2(1)
      real tmp1(lx1*ly1*lz1*lelt*2)
      real tmp2(lx1*ly1*lz1*lelt*2)
      integer imesh,isd
      npts2=npts*2
      n2=npts+1
      call copy(tmp1,u,npts2)
      call axhelm3(au(1),tmp1(1),helm1(1),helm2(1),imesh,isd)
      call axhelm3(au(n2),tmp1(n2),helm1(n2),helm2(n2),imesh,isd)
      if (IFNEU) THEN
ccc      adding the Neumann term
c      call copy(tmp1,u,npts2)
c      call col2    (tmp1,kwave,npts2)   
c      call rzero(tmp2,npts2)
c      call cem_acoustic_bd(tmp2(1),tmp1(1))
c      call cem_acoustic_bd(tmp2(n2),tmp1(n2))
c      call add2s2  (au(1),tmp2(n2),-1.,npts)                ! Aur = Aur + ui
c      call add2s2  (au(n2),tmp2(1), 1.,npts)               ! Aui = Aui - ur
      ENDIF
      if (IFDTN) THEN
      call copy(tmp1,u,npts2)
      call cem_acoustic_bd_dtn(tmp2,tmp1)
      call add2s2  (au,tmp2,1.,npts2)                ! Aur = Aur + ui
      endif       
      return
      end
C------------------------------------------------------------------
      subroutine axhelm3(au,u,helm1,helm2,imesh,isd)
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
      include 'EMWAVE'
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
c      NXYZ= NX1*NY1*NZ1
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
 
      CALL ADDCOL4 (AU,HELM2,BM1,U,NTOT)
 
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
c      if (nid.eq.0) write(6,*) 'done: axhelm3' 
C
#ifndef NOTIMER
      taxhm=taxhm+(dnekclock()-etime1)
#endif
      return
      END

c-----------------------------------------------------------------------
      subroutine cem_acoustic_bd(bdint,bdfunc)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'

      real     aa,bdint(1),bdfunc(1)
      integer  i,i0


          do i=1,ncemface_dtn
             i0= cemface_dtn(i)
             aa=aream_dtn(i)
             bdint(i0)=bdfunc(i0)*aa
          enddo

      return 
      end

c-----------------------------------------------------------------------
      subroutine cem_acoustic_dtn
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'
 
      integer  nxyzf,nface,nxz,k,i,e,f,ef,p,pz
      character CB*3
      real     length,tmp,sdtn(lxzf1*2*np*lfp2*lfp2z)
      integer  p1r,p1i,p2r,p2i
      integer  p1rz,p1iz,p2rz,p2iz
      integer  maxdtn2
c
      integer  icalld
      save     icalld
      data     icalld /0/
      maxdtn2= maxdtn*2
      ndtn   = maxdtn2*np*lfp2*lfp2z
      k      = maxdtn
c      length=abs(xmax-xmin)
      call rzero(sdtn,ndtn)
      if (if3d) then
       length=2.0  ! FIXME LATER decide the exact value
      else
       length=2.0*pi  ! FIXME LATER decide the exact value
      endif
      tmp=sqrt(length)
c      tmp=tmp**(ndim-1)
      alfa=0.0  ! FIXME LATER decide the exact value
      gama0=0.0  ! FIXME LATER decide the exact value
      do pz= 1,lfp2z
         gamap(pz)=gama0+(pz-lfpz-1)*2*pi/length
      do p= 1,lfp2
         alphap(p)=alfa+(p-lfp-1)*2*pi/length
         npxz=p+(pz-1)*lfp2
         betap(npxz) =abs(kwave(1)**2-alphap(p)**2-gamap(pz)**2)
         betap(npxz)=sqrt(betap(npxz))
         do i= 1,ncemface_dtn
         i0= cemface_dtn(i)
         kk=i+(p-1)*maxdtn2+(pz-1)*lfp2*maxdtn2+maxdtn2*lfp2*lfp2z*nid
c simple case
             unx0=unxm_dtn(i)
             uny0=unym_dtn(i)
             unz0=unzm_dtn(i)
             if (if3d) then
c             xx=xmn(i0)
             xx=ymn(i0)
c             zz=zmn(i0)
             zz=ymn(i0)
             aa=aream_dtn(i)
             sdtn(kk)= (cos(alphap(p)*xx)*cos(gamap(pz)*zz)
     $                -sin(alphap(p)*xx)*sin(gamap(pz)*zz))*aa/length
             sdtn(maxdtn+kk)=-(sin(alphap(p)*xx)*cos(gamap(pz)*zz)
     $                +cos(alphap(p)*xx)*sin(gamap(pz)*zz))*aa/length
             else
             yy=ymn(i0)
             aa=aream_dtn(i)
             sdtn(kk)=cos(alphap(p)*yy)*aa/tmp
             sdtn(maxdtn+kk)=-sin(alphap(p)*yy)*aa/tmp
             endif
         enddo
      enddo
      enddo
      do i=1,ndtn
        sdtn(i)= glsum(sdtn(i),1)
      enddo
      k=2*maxdtn**2
      do nip=0,np-1
         do i2=1,maxdtn
         do i =1,maxdtn
             npp2=maxdtn2*lfp2*lfp2z*nip
             npp =maxdtn2*lfp2*lfp2z*nid
             ssumr=0.0
             ssumi=0.0
             do pz=1,lfp2z
             npz =lfp2z-pz+1
             p1rz=maxdtn2*lfp2*(pz-1)+npp
             p2rz=maxdtn2*lfp2*(npz-1)+npp2
             p1iz=maxdtn2*lfp2*(pz-1)+npp
             p2iz=maxdtn2*lfp2*(npz-1)+npp2
             do p=1,lfp2
             npx=lfp2-p+1
             p1r=i+maxdtn2*(p-1)+p1rz
             p2r=i2 +maxdtn2*(npx-1)+p2rz
             p1i=maxdtn+i+maxdtn2*(p-1)+p1iz
             p2i=maxdtn+i2 +maxdtn2*(npx-1)+p2iz
             ssumr0=sdtn(p1r)*sdtn(p2r)-sdtn(p1i)*sdtn(p2i)
             ssumi0=sdtn(p1r)*sdtn(p2i)+sdtn(p1i)*sdtn(p2r)
             npxz=p+(pz-1)*lfp2
             if ((kwave(1)**2-alphap(p)**2-gamap(pz)**2).gt.0) then
             ssumr=ssumr-betap(npxz)*ssumi0
             ssumi=ssumi+betap(npxz)*ssumr0
             else
             ssumr=ssumr-betap(npxz)*ssumr0
             ssumi=ssumi-betap(npxz)*ssumi0
             endif
             enddo
             enddo
             wrkdtn2(i+(i2-1)*maxdtn2+k*nip)=ssumr
             wrkdtn2(maxdtn+i+(i2-1)*maxdtn2+k*nip)=ssumi
         enddo
         enddo
      enddo
      return 
      end

c-----------------------------------------------------------------------
      subroutine cem_acoustic_bd_dtn(bout,bin)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'

      real     bout(1),bin(1)
      real     fbin(lxzf3*2)
      real     wrk(lxzf3*2*np)
      real     ssum,ssum2
      integer  nxz,k,i,e,f,p,ef
      integer  nelt0,nface0,nxz0
      integer  p1r,p1i,p2r,p2i
      integer  pz, p1rz,p1iz,p2rz,p2iz
      integer  i2,e2,f2
      integer  maxdtn2
      
     
      maxdtn2= maxdtn*2
      npts2=npts*2
      k    = maxdtn2
      call rzero(bout, npts2)
      call rzero(fbin, maxdtn2)
      call rzero(srflxc,maxdtn2*np)
      call rzero(wrk,maxdtn2*np)
      call cem_restrict_to_dtn_face(fbin(1),bin(1))
      call cem_restrict_to_dtn_face(fbin(maxdtn+1),bin(npts+1))
         do nip=0,np-1
            kk   =2*maxdtn**2*nip
         do i= 1,maxdtn
             ssum=0.0
             ssum2=0.0
             ki   =(i-1)*maxdtn2+kk
             do i2 =1,ncemface_dtn
             ssumr=wrkdtn2(i2+ki)
             ssumi=wrkdtn2(maxdtn+i2+ki)
             ssum= ssum+fbin(i2)*ssumr
     $                 -fbin(maxdtn+i2)*ssumi
             ssum2=ssum2+fbin(i2)*ssumi
     $                 +fbin(maxdtn+i2)*ssumr
             enddo
             srflxc(i+k*nip)=ssum
             srflxc(maxdtn+i+k*nip)=ssum2
         enddo
         enddo
c      do i=1,maxdtn2*np
c         srflxc(i)=glsum(srflxc(i),1) 
c      enddo
         call gop(srflxc,wrk,'+  ',maxdtn2*np) 
      do i=1,ncemface_dtn
         i0=cemface_dtn(i)
         bout(i0)=srflxc(i+k*nid)
         bout(npts+i0)=srflxc(maxdtn+i+k*nid)
      enddo
      return 
      end

c-----------------------------------------------------------------------
      subroutine cem_acoustic_test(ssum,bdfunc)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'

      real     bdfunc(1),ssum
      integer  nxyzf,nface,nxz,k,i,e,f,ef
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/
      nxyzf= nx1*nz1*2*ndim*nelt
      nface= 2*ndim
      nxz  = nx1*nz1
      k    = nxyzf
      
      ssum=0.0
      call cem_acoustic_bd(rhs,bdfunc)
c      call dssum   (rhs,nx1,ny1,nz1)
c      call col2    (rhs,mask0,npts)   
         do e= 1,nelt
         do f= 1,nfaces
         do j= 1,nxzf

         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,2)
         if (CB.eq.'DTN') then
             ssum=ssum+rhs(i0)
         endif
         enddo
         enddo
         enddo
      if (if3d)then
      ssum=ssum/4.0
      else
      ssum=ssum/(2.0*pi)
      endif
      ssum= glsum(ssum,1)
      return 
      end

c----------------------------------------------------------------------

      subroutine cem_set_dtn_fc_ptr
c
c     Set up pointer to restrict u to dtn faces ! NOTE: compact
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
c
      integer i,j,k,e,f,ef
      character CB*3
      nxz   = nx1*nz1
      nface = 2*ndim

      if (nid.eq.0) then
          write(6,*) 'start: cem_set_dtn_fc_ptr, 
     $                volume index j=cemface_dtn(i)'
      endif

      k= 0

      do e= 1,nelt
      do f= 1,nface
      do j= 1,nxz
         i = (e-1)*nxz*nface + nxz*(f-1) + j  ! face   numbering
         i0= cemface(i)                       ! global numbering
         CB = CBC(f,e,2)
         if (CB.eq.'DTN') then
            k=k+1                           ! face   numbering
            cemface_dtn(k) =i0               ! global numbering
            aream_dtn  (k) =aream(i)   
            unxm_dtn   (k) =unxm (i)   
            unym_dtn   (k) =unym (i)   
            unzm_dtn   (k) =unzm (i)   
         endif
      enddo
      enddo
      enddo
      ncemface_dtn = k
      maxdtn=iglmax(ncemface_dtn,1)
      mindtn=iglmin(ncemface_dtn,1)
      if (nid.eq.0) then
         write(6,*) 'done:cem_set_dtn_fc_ptr,ncemface_dtn=',ncemface_dtn
         write(6,*) 'done:cem_set_dtn_fc_ptr,max_dtn=',maxdtn
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine cem_restrict_to_dtn_face(fu,u)
c     Restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      real fu(1),u(1)

      do j=1,ncemface_dtn
          i= cemface_dtn(j)
           fU(j)   = U(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
