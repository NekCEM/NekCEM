c=======================================================================
      subroutine hmhPEC
c=======================================================================
c     solve  alpha nabla^2 (u) +k^2 u=f
c     with periodic in X direction
c          u(x,0)=0 and DtN du/dy(x,1)-\beta u(x,1)=g(x)
c     remark: sign(alpha)*sign(h1)=-1
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer  icalld
      save     icalld
      data     icalld /0/

      common /ccpu/  cpu_t,cpu_dtime,cpu_chk
      real*8         cpu_t,cpu_dtime,cpu_chk
      real*8         dclock,tol,ssum
      integer        i,isd,imsh,npts2
      real           constant
      integer        idir,l

      cpu_dtime= dclock()          
      cpu_t    = 0.0

      if (nid.eq.0) write(6,*) '\n call hmhDtN'       

      npts2 =  npts*2
      call rone    (mask, npts2)                                      
      call setmaskp(mask(     1),npts) ! mask is zero for PEC, PML 
      call setmaskp(mask(npts+1),npts) ! mask is zero for PEC, PML 
      if (nid.eq.0) write(6,*) ' done: setmaskp'       

      call rone    (mult, npts2)
      call dssum   (mult(     1),nx1,ny1,nz1)
      call dssum   (mult(npts+1),nx1,ny1,nz1)
      call invcol1 (mult,npts2)
      if (nid.eq.0) write(6,*) ' done: mult'       

      call rone    (h1, npts2)         ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call chsign  (h1, npts2)         ! negative for helmholtz, positive for poisson
      call rone    (h2, npts2)         ! call rzero   (h2,npts2)  if poisson 
      call col2    (h2, kwave2, npts2) ! FIXME kwave: to be redefined for multilayer
      if (nid.eq.0) write(6,*) ' done: set h1 and h2 in acoustic'       

     
      call rzero(rhs,npts2)            ! init
      call copy (rhs,potent,npts2)     ! force term: potent initially
      call col2 (rhs(     1),bm1,npts) ! B*f
      call col2 (rhs(npts+1),bm1,npts) ! B*f
      call dssum(rhs(     1),nx1,ny1,nz1)   
      call dssum(rhs(npts+1),nx1,ny1,nz1)   
      
      isd   = 1
      imsh  = 1
      tol   = 1.e-10 ! 1.e-10

      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES

      if (IFGMRES) then
          call copy      (potent,spotent,npts2)    !passing initial guess in usr (other than zero if necessary)
          call hmh_gmres5(potent,rhs,h1,h2,mult,mask,isd,imsh,npts2,tol)
         if (nid.eq.0) write(6,*) 'done: hmh_gmres3'
      else 
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
         call exitt
      endif
       

      if (nid.eq.0) write(6,*) '\n done hmhNEU'       

      cpu_dtime= dclock()-cpu_dtime
      cpu_t= cpu_t+cpu_dtime

      return
      end
c=======================================================================
      subroutine hmhNEU
c=======================================================================
c     solve  alpha nabla^2 (u) +k^2 u=f
c     with periodic in X direction
c          u(x,0)=0 and DtN du/dy(x,1)-\beta u(x,1)=g(x)
c     remark: sign(alpha)*sign(h1)=-1
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer  icalld
      save     icalld
      data     icalld /0/

      common /ccpu/  cpu_t,cpu_dtime,cpu_chk
      real*8         cpu_t,cpu_dtime,cpu_chk
      real*8         dclock,tol,ssum
      integer        i,isd,imsh,npts2
      real           constant
      integer        idir,l

      cpu_dtime= dclock()          
      cpu_t    = 0.0

      if (nid.eq.0) write(6,*) '\n call hmhDtN'       

      npts2 =  npts*2
      call rone    (mask, npts2)                                      
      call setmaskp(mask(     1),npts) ! mask is zero for PEC, PML 
      call setmaskp(mask(npts+1),npts) ! mask is zero for PEC, PML 
      if (nid.eq.0) write(6,*) ' done: setmaskp'       

      call rone    (mult, npts2)
      call dssum   (mult(     1),nx1,ny1,nz1)
      call dssum   (mult(npts+1),nx1,ny1,nz1)
      call invcol1 (mult,npts2)
      if (nid.eq.0) write(6,*) ' done: mult'       

      call rone    (h1, npts2)         ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call chsign  (h1, npts2)         ! negative for helmholtz, positive for poisson
      call rone    (h2, npts2)         ! call rzero   (h2,npts2)  if poisson 
      call col2    (h2, kwave2, npts2) ! FIXME kwave: to be redefined for multilayer
      if (nid.eq.0) write(6,*) ' done: set h1 and h2 in acoustic'       

     
      call rzero              (rhs,npts2)                   ! init
      
      if (IFNEU) then
      call cem_set_nmn_fc_ptr        
      call cem_acoustic_bd_nmn(rhs(     1),rhs_nmn(     1)) ! rhs  <-- neumann condition on DTN faces      
      call cem_acoustic_bd_nmn(rhs(npts+1),rhs_nmn(npts+1)) ! 
      call sub2               (rhs,rhs_dtn,npts2)           ! rhs =  rhs_nmn - rhs_dtn 
      call chsign             (rhs,npts2)                   ! rhs = -rhs
      if (nid.eq.0) write(6,*) ' done: neumann conditions'      
      endif
      
      call dssum   (rhs(     1),nx1,ny1,nz1)   
      call dssum   (rhs(npts+1),nx1,ny1,nz1)   
    
      isd   = 1
      imsh  = 1
      tol   = 1.e-10 ! 1.e-10

      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES

      if (IFGMRES) then
          call copy      (potent,spotent,npts2)    !passing initial guess in usr (other than zero if necessary)
          call hmh_gmres3(potent,rhs,h1,h2,mult,mask,isd,imsh,npts2,tol)
         if (nid.eq.0) write(6,*) 'done: hmh_gmres3'
      else 
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
         call exitt
      endif
       

      if (nid.eq.0) write(6,*) '\n done hmhNEU'       

      cpu_dtime= dclock()-cpu_dtime
      cpu_t= cpu_t+cpu_dtime

      return
      end
c=======================================================================
      subroutine hmhDtN
c=======================================================================
c     solve  alpha nabla^2 (u) +k^2 u=f
c     with periodic in X direction
c          u(x,0)=0 and DtN du/dy(x,1)-\beta u(x,1)=g(x)
c     remark: sign(alpha)*sign(h1)=-1
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer  icalld
      save     icalld
      data     icalld /0/

      common /ccpu/  cpu_t,cpu_dtime,cpu_chk
      real*8         cpu_t,cpu_dtime,cpu_chk
      real*8         dclock,tol,ssum
      integer        i,isd,imsh,npts2
      real           constant
      integer        idir,l

      cpu_dtime= dclock()          
      cpu_t    = 0.0

      if (nid.eq.0) write(6,*) '\n call hmhDtN'       

      npts2 =  npts*2
      call rone    (mask0,npts )      ! for grids                         
      call setmaskp(mask0,npts )      ! mask0 is zero for PEC, PML 
      call copy    (mask(1),mask0,npts)
      call copy    (mask(npts+1),mask0,npts)
      if (nid.eq.0) write(6,*) ' done: setmaskp'       

      call rone    (mult0,npts)
      call dssum   (mult0,nx1,ny1,nz1)
      call invcol1 (mult0,npts)
      call copy    (mult(1),mult0,npts)
      call copy    (mult(npts+1),mult0,npts)
      if (nid.eq.0) write(6,*) ' done: mult'       

      call rone    (h1,npts2)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call chsign  (h1,npts2)        ! negative for helmholtz, positive for poisson
      call rone    (h2,npts2)        ! call rzero   (h2,npts2)  if poisson 
      call col2    (h2,kwave2,npts2) ! FIXME kwave: to be redefined for multilayer
      if (nid.eq.0) write(6,*) ' done: set h1 and h2 in acoustic'       

      call cem_set_dtn_fc_ptr        

      if (param(75).eq.1)  call cem_acoustic_dtn2  ! quadrature
      if (param(75).eq.0)  call cem_acoustic_dtn   ! bessel

      constant = -1.0
 
      call rzero              (rhs_dtn,npts2)               ! init                                   
      call cem_acoustic_bd_dtn(rhs_dtn,rhs_inc)             ! rhs_dtn <-- incident field on DTN faces

      call rzero              (rhs,npts2)                   ! init
      call cem_acoustic_bd_nmn(rhs(     1),rhs_nmn(     1)) ! rhs  <-- neumann condition on DTN faces      
      call cem_acoustic_bd_nmn(rhs(npts+1),rhs_nmn(npts+1)) ! 
      call add2s2             (rhs,rhs_dtn,constant,npts2)  ! rhs =  rhs_nmn - rhs_dtn 
      call chsign             (rhs,npts2)                   ! rhs = -rhs
      if (nid.eq.0) write(6,*) ' done: dtn and neumann conditions'      

      call dssum   (rhs(1),     nx1,ny1,nz1)   
      call dssum   (rhs(npts+1),nx1,ny1,nz1)   

      isd   = 1
      imsh  = 1
      tol   = 1.e-10 ! 1.e-10

      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES

      if (IFGMRES) then
          call copy      (potent,spotent,npts2)    !passing initial guess in usr (other than zero if necessary)
          call hmh_gmres3(potent,rhs,h1,h2,mult,mask,isd,imsh,npts2,tol)
         if (nid.eq.0) write(6,*) 'done: hmh_gmres3'
      else 
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
         call exitt
      endif
      
       

      if (nid.eq.0) write(6,*) '\n done hmhDtN'       

      cpu_dtime= dclock()-cpu_dtime
      cpu_t= cpu_t+cpu_dtime

      return
      end
c=======================================================================
      subroutine hmh_gmres4(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,n1,n2,outer,isd,imsh,iter,m,i,j,k,iconv
      real     x2(n),r2(n),wt2(n),w2(n),v2(n,lgmres)
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp,tolpss,glsc3,rnorm
      real*8   etime1,etime_p,divex,dnekclock

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
c      call rone(wt,n)
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
c           if (nid.eq.0) write(6,*) iter,'gmres-rnorm',rnorm
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
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'

      common /wgmres/ x2(lpts3*2)
     $              , r2(lpts3*2)
     $              , w2(lpts3*2)
     $              , v2(lpts3*2,lgmres)  
      real            x2,r2,w2,v2


      integer  n,n1,n2,outer,isd,imsh,iter,iconv,i,j,k,m
c     real     x2(n),r2(n),w2(n),v2(n,lgmres)
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp,tolpss,glsc3,rnorm
      real*8   etime1,etime_p,divex,dnekclock

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
c     call copy(x2,phi,n) ! none zero inital guess
      call rone (x2,n)    !zero inital guess
      iter  = 1
      call dssum(x2(1),nx1,ny1,nz1)
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
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
c      INCLUDE 'CTIMER'
      include 'EMWAVE'
      INCLUDE 'ACOUSTIC'
C
      real au(1),u(1),helm1(1),helm2(1)
      real tmp1(lx1*ly1*lz1*lelt*2)
      real tmp2(lx1*ly1*lz1*lelt*2)
      real wrkdx(lx1*ly1*lz1*lelt*2)
      real wrkdy(lx1*ly1*lz1*lelt*2)
      real wrkdz(lx1*ly1*lz1*lelt*2)
      integer npts2,n2,imesh,isd
      real alp0
      npts2 =npts*2
      n2    =npts+1
      call copy(tmp1,u,npts2)
      call axhelm3(au(1),tmp1(1),helm1(1),helm2(1),imesh,isd)
      call axhelm3(au(n2),tmp1(n2),helm1(n2),helm2(n2),imesh,isd)
c      call acoustic_grad(wrkdx( 1),wrkdy( 1),wrkdz( 1),tmp1( 1))
c      call acoustic_grad(wrkdx(n2),wrkdy(n2),wrkdz(n2),tmp1(n2))
c      call col2    (wrkdx(1),bm1,npts)   ! B*dx
c      call col2    (wrkdx(n2),bm1,npts)   ! B*dx
      call acoustic_wght_grad(wrkdx( 1),wrkdy( 1),wrkdz( 1),tmp1( 1))
      call acoustic_wght_grad(wrkdx(n2),wrkdy(n2),wrkdz(n2),tmp1(n2))

      alp0=2.0*alpha0(1)
      call add2s2  (au( 1),wrkdx(n2),-alp0,npts)              ! Aur = Aur + ui
      call add2s2  (au(n2),wrkdx( 1), alp0,npts)              ! Aui = Aui - ur

      if (IFDTN) THEN
      call copy(tmp1,u,npts2)
      call rzero(tmp2,npts2)
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
      implicit none
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
      REAL WDDX,WDDYT,WDDZT
      INTEGER NEL,NXY,NYZ,NXZ,NTOT

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

      real            TMP1, TMP2,TMP3
      real            DUDR, DUDS,DUDT
      
      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      real           term1,term2,dnekclock
      integer  iel,ie,iz,isd,imesh,j,i
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)
C
      NEL = NELT
      NXY = NX1*NY1
      NYZ = NY1*NZ1
      NXZ = NX1*NZ1
      NTOT= NX1*NY1*NZ1*NEL
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
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
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
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
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
      subroutine cem_acoustic_bd_nmn(bdint,bdfunc)
      implicit none
      include 'SIZE'
      include 'TOTAL'  
      include 'ACOUSTIC'

      real     aa,bdint(1),bdfunc(1)
      integer  i,i0,idtn,idir

      do idir= 1,6
      if (l_dtnface(idir).ne.0) then 
         do i = 1,l_dtnface(idir)
            i0= dtn_face(i,idir)
            aa= area_dtn(i,idir)
            bdint(i0)=bdfunc(i0)*aa
         enddo
      endif
      enddo

      return 
      end

c-----------------------------------------------------------------------
      subroutine cem_acoustic_bd2(bdint,bdfunc)
      implicit none
      include 'SIZE'
      include 'TOTAL'  
      include 'ACOUSTIC'

      real     aa,bdint(1),bdfunc(1)
      integer  i,i0,idtn,idir


      do idir= 1,6
      if (l_dtnface(idir).ne.0) then
         do i = 1,l_dtnface(idir)
         i0= dtn_face(i,idir)
         aa= area_dtn(i,idir)
         bdint(i0)=bdfunc(i0)*aa*(-1)**idir
         enddo
      endif
      enddo

      return 
      end

c-----------------------------------------------------------------------
      subroutine cem_acoustic_dtn2        
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'

      common /DTNV1/ sdtn(lxzfl3*2*lfp2*lfp2z)
      real           sdtn
      real      length,tmp,tmp_area      
      integer   p,pz,p1r,p1i,p2r,p2i,ip,pn,k2
      integer   ndtn,npxz,i,i0,kk,k,maxdtn2,nip,i2,npp2,npp,npz,npx
      character CB*3
c
      real*8   alfa,gama0,unx0,uny0,unz0,xx,zz,aa,yy,glsum,ssumr,ssumi
      real*8   ssumr0,ssumi0
      integer  icalld
      save     icalld
      data     icalld /0/
      real     glmax,glmin,xmax,xmin,ymax,ymin,zmax,zmin
      integer  idtn,idir

      if (nid.eq.0) write(6,*) 'start: quadrature dtn'

      do idir= 1,6        !FIXME multilayer

      if (l_dtnface(idir).ne.0) then

         if     ((idir.eq.1).or.(idir.eq.2)) then ! DTN in x-direction
          xmax   = glmax(xm1,npts)
          xmin   = glmin(xm1,npts)
          length = sqrt(xmax-xmin)
         elseif ((idir.eq.3).or.(idir.eq.4)) then ! DTN in y-direction
          ymax   = glmax(ym1,npts)
          ymin   = glmin(ym1,npts)
          length = sqrt(ymax-ymin)
         elseif ((idir.eq.5).or.(idir.eq.6)) then ! DTN in y-direction
          zmax   = glmax(zm1,npts)
          zmin   = glmin(zm1,npts)
          length = sqrt(zmax-zmin)
         endif

         maxdtn = l_dtnface(idir)
         maxdtn2= maxdtn*2
         ndtn   = maxdtn2*np*lfp2*lfp2z
         k      = maxdtn
         pn     = lfp

         call rzero(sdtn,ndtn)

         do p = -pn,pn
            ip        = p+pn+1
            alphap(ip)= 2*p*pi/length
            betap(ip) = abs(kw(idir)**2-alphap(ip)**2)     !FIXME
            betap(ip) = sqrt(betap(ip))
          
         do i  = 1,l_dtnface(idir)
             i0 = dtn_face(i,idir)
             kk = i+(ip-1)*maxdtn2

             if ((idir.eq.1).or.(idir.eq.2)) tmp=xmn(i0) !FIXME
             if ((idir.eq.3).or.(idir.eq.4)) tmp=ymn(i0) !FIXME
             if ((idir.eq.5).or.(idir.eq.6)) tmp=zmn(i0) !FIXME

             tmp_area  = area_dtn(i,idir)

             sdtn(kk)  = cos(alphap(ip)*tmp)*tmp_area/length
             sdtn(k+kk)=-sin(alphap(ip)*tmp)*tmp_area/length

         enddo
         enddo

         do i=1,ndtn
            sdtn(i)= glsum(sdtn(i),1)
         enddo

         k = maxdtn
         k2= maxdtn*2

         do i2=1,maxdtn
         do i =1,maxdtn
         ssumr=0.0
         ssumi=0.0
             
         do p =-pn,pn
             ip = p+pn+1
             p1r= i   +k2*(ip-1)
             p2r= i2  +k2*(ip-1)
             p1i= k+i +k2*(ip-1)
             p2i= k+i2+k2*(ip-1)

             ssumr0=sdtn(p1r)*sdtn(p2r)-sdtn(p1i)*sdtn(p2i)
             ssumi0=sdtn(p1r)*sdtn(p2i)+sdtn(p1i)*sdtn(p2r)

             npxz=p+(pz-1)*lfp2

             if ((kw(idir)**2-alphap(ip)**2).gt.0) then
             ssumr=ssumr-betap(ip)*ssumi0
             ssumi=ssumi+betap(ip)*ssumr0
             else
             ssumr=ssumr-betap(ip)*ssumr0
             ssumi=ssumi-betap(ip)*ssumi0
             endif
         enddo
         wrkdtn2(i+(i2-1)*k2)=ssumr
         wrkdtn2(k+i+(i2-1)*k2)=ssumi
         enddo
         enddo

      endif
      enddo !FIXME multilayer

      if (nid.eq.0) write(6,*) ' done : quadrature dtn'

      return 
      end
c-----------------------------------------------------------------------
      subroutine cem_acoustic_dtn              
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'

      common /DTNV1/ sdtn(lxzfl3*2*lfp2*lfp2z)
      real           sdtn
      common /DTNV2/ xmaxl(lelt),xminl(lelt)
     $             , ymaxl(lelt),yminl(lelt)
     $             , zmaxl(lelt),zminl(lelt)
     $             , tmaxl(lelt),tminl(lelt)
     $             , tmph (lelt),tmpd (lelt)
      real    glsum,glmax,glmin,vlmax,vlmin
      real    xmaxl,xminl,ymaxl,yminl,zmaxl,zminl
      real    xmax,xmin,ymax,ymin,zmax,zmin,tmax,tmin
      real    tmaxl,tminl,tmph,tmpd
      real    length,ar,bi,wr,wi
      integer p,pz,p1r,p1i,p2r,p2i
      integer p1rz,p1iz,p2rz,p2iz
      integer ndtn,npxz,maxdtn2,nip,i2,npp2,npp,npz,npx
      integer i0,i,j,ip,k,k2,ksum,kn,pn,e     
c
      real*8   alpha,gama0,xx,zz,aa,yy
      real*8   ssumr0,ssumi0,ssumr,ssumi
      integer  icalld
      save     icalld
      data     icalld /0/
      
      real     PNLEG,rr(lx1),ww(lx1),wk(lx1,lx1),wk2(lx1*lx1)
      real     tmp,tmp0,tmp2,tmp1,alp,bsslj
      real     tmp3(lx1),tmp3r(lx1),tmp3i(lx1)
      real     dtnr (lx1,2*lfp+1),dtni (lx1,2*lfp+1),
     $         dtnr1(lx1,2*lfp+1),dtni1(lx1,2*lfp+1)
      integer  ipiv(lx1),info,ip1,ip2
      integer  idtn,idir,dir


      ksum=0
      do idir= 1,6        !FIXME multilayer
      if (l_dtnface(idir).ne.0) then

         if     ((idir.eq.1).or.(idir.eq.2)) then ! DTN in x-direction
          tmax  = glmax(ym1,npts)
          tmin  = glmin(ym1,npts)
          length= tmax-tmin       
          tmp   = sqrt(length)
          do e= 1,nelt
             tmaxl(e)= vlmax(ym1(1,1,1,e),nxyz) 
             tminl(e)= vlmin(ym1(1,1,1,e),nxyz)
             tmph (e)= tmaxl(e)-tminl(e)       ! FIXME
             tmpd (e)= tmaxl(e)+tminl(e)       ! FIXME
          enddo
         elseif ((idir.eq.3).or.(idir.eq.4)) then ! DTN in y-direction
          tmax  = glmax(xm1,npts)
          tmin  = glmin(xm1,npts)
          length= tmax-tmin       
          tmp   = sqrt(length)
    
          do e= 1,nelt
             tmaxl(e)= vlmax(xm1(1,1,1,e),nxyz)
             tminl(e)= vlmin(xm1(1,1,1,e),nxyz)
             tmph (e)= tmaxl(e)-tminl(e)
             tmpd (e)= tmaxl(e)+tminl(e)
          enddo
         elseif ((idir.eq.5).or.(idir.eq.6)) then ! FIXME               
          tmax  = glmax(zm1,npts)
          tmin  = glmin(zm1,npts)
          length= tmax-tmin       
          tmp   = sqrt(length)
          do e= 1,nelt
             tmaxl(e)= vlmax(zm1(1,1,1,e),nxyz) ! nxyz=nx1*ny1*nz1
             tminl(e)= vlmin(zm1(1,1,1,e),nxyz)
             tmph (e)= tmaxl(e)-tminl(e)
             tmpd (e)= tmaxl(e)+tminl(e)
          enddo
         endif

cdddd
         maxdtn = l_dtnface(idir)
         k      = maxdtn
         k2     = maxdtn*2
         pn     = lfp 
         ndtn   = k2*(pn*2+1)     

         call rzero(sdtn,ndtn)
      
      !obtain Legendre points, weights on [-1,1]
         call ZWGLL(rr,ww,nx1)
      
      !compute the value of the 1th order ledendre poly at rr(j)
         do j=1,nx1
         do i=1,nx1
            wk(i,j)=PNLEG(rr(j),i-1)
         enddo
         enddo 

      !get the inverse of matrix wk:= A
      !use lapack subroutine:DGETRF+DGETRI
         call DGETRF( nx1, nx1, wk, nx1, IPIV, INFO )
         call DGETRI( nx1, wk, nx1, IPIV, wk2, nx1*nx1, INFO )

      !compute the fourier of Ledendr using Bessel:= B
      !for p=0

         do i = 1,nx1
         info= mod(i-1,4)
         if      (info .eq. 0) then
                 tmp3r(i)= 1.0
                 tmp3i(i)= 0.0
         else if (info .eq. 1) then
                 tmp3r(i)= 0.0
                 tmp3i(i)= 1.0
         else if (info .eq. 2) then
                 tmp3r(i)=-1.0
                 tmp3i(i)= 0.0
         else if (info .eq. 3) then
                 tmp3r(i)= 0.0
                 tmp3i(i)=-1.0
         endif
         enddo

         pn = lfp
         dtnr1(1,pn+1)=2.0
         dtni1(1,pn+1)=0.0

         do i=2,nx1
         dtnr1(i,pn+1)=0.0
         dtni1(i,pn+1)=0.0
         enddo      
      
         do p=1,pn 
           ip1 = p+pn+1
           tmp0= 2*p*pi/length
           tmp1= tmp0*tmph(1)/2.0  ! size info  
           tmp2= sqrt(2.0*pi/tmp1)
           call rjbesl(tmp1,0.5,nx1,tmp3,info)
           do i=1,nx1
           dtnr1(i,ip1)=tmp2*tmp3(i)*tmp3r(i)*(-1)**(i-1)
           dtni1(i,ip1)=tmp2*tmp3(i)*tmp3i(i)*(-1)**(i-1)
           enddo
         enddo

         do p = pn,1,-1
           ip2=-p+pn+1
           ip1= p+pn+1
           do i=1,nx1
           dtni1(i,ip2)=dtni1(i,ip1)*(-1)**(i-1)
           dtnr1(i,ip2)=dtnr1(i,ip1)*(-1)**(i-1)
           enddo
         enddo

      ! compute: sdtn(nx1,p) instead of matlab
      ! sdtn=A*B
         call mxm(wk,nx1,dtnr1,nx1,dtnr,lfp2)
         call mxm(wk,nx1,dtni1,nx1,dtni,lfp2)

         do p = -pn,pn
         
           ip = p+pn+1
           alphap(ip)= 2*p*pi/length
           betap (ip)= abs (kw(idir)**2-(alpha0(1)+alphap(ip))**2) ! FIXME: kwave's dimension
           betap (ip)= sqrt(betap(ip))

           do i = 1,l_dtnface(idir)
              j = dtn_face(i,idir)
              kn= i+(ip-1)*k2      
              i0= mod(i,nx1)  ! FIXME   
              if (i0.eq.1) i2=0    
              i2= i2+1
              e         = dtn_e(i,idir)
              wr        = tmph(e)/(2.0*tmp)*dtnr(i2,ip)  !update with input from matlab 
              wi        = tmph(e)/(2.0*tmp)*dtni(i2,ip)  !update
              ar        = cos(alphap(ip)*tmpd(e)/2.0)
              bi        =-sin(alphap(ip)*tmpd(e)/2.0)
              sdtn(  kn)= ar*wr-bi*wi
              sdtn(k+kn)= ar*wi+bi*wr
           enddo

         enddo

c     k=2*maxdtn**2

         do i2=1,k     
         do i =1,k      
         ssumr=0.0
         ssumi=0.0

         do p =-pn,pn

            ip = p+pn+1
            p1r=   i  +k2*(ip-1)
            p2r=   i2 +k2*(ip-1)
            p1i= k+i  +k2*(ip-1)
            p2i= k+i2 +k2*(ip-1)

            ssumr0= sdtn(p1r)*sdtn(p2r)-sdtn(p1i)*sdtn(p2i)
            ssumi0= sdtn(p1r)*sdtn(p2i)+sdtn(p1i)*sdtn(p2r)

            if ((kw(idir)**2-alphap(ip)**2).gt.0) then ! FIXME

                ssumr=ssumr-betap(ip)*ssumi0
                ssumi=ssumi+betap(ip)*ssumr0
            else
                ssumr=ssumr-betap(ip)*ssumr0
                ssumi=ssumi-betap(ip)*ssumi0
            endif
         enddo
         
         wrkdtn2(  i+(i2-1)*k2+ksum)=ssumr
         wrkdtn2(k+i+(i2-1)*k2+ksum)=ssumi
         enddo
         enddo

         ksum=ksum+k2*k

      endif
      enddo !FIXME multilayer

      if (nid.eq.0) write(6,*) ' done : bessel dtn'

      return 
      end

c-----------------------------------------------------------------------
      subroutine cem_acoustic_bd_dtn(bout,bin)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'
      integer  idtn,idir
      real     bout(1),bin(1)
      real     fbin(lxzfl3*2)
      real     wrk (lxzfl3*2)
      real     wrk1(lxzfl3*2)
      real     wrk2(lxzfl3*2)
      real     ssum1,ssum2,ssumr,ssumi,glsum
      integer  maxdtn_tot,maxdtn2,mm,i,j     
      integer  ksum1,ksum2,jump,j1,j2
      
      maxdtn_tot= 0
      do idir= 1,6
         if (l_dtnface(idir).ne.0) then
             maxdtn_tot = maxdtn_tot + l_dtnface(idir)
         endif
      enddo

      mm = maxdtn_tot
      call rzero(fbin,mm*2)
      call cem_restrict_to_dtn_face(fbin(   1),bin(     1))
      call cem_restrict_to_dtn_face(fbin(mm+1),bin(npts+1))
      
      ksum1= 0
      ksum2= 0

      do idir= 1,6      

        if (l_dtnface(idir).ne.0) then

         maxdtn = l_dtnface(idir) 
         maxdtn2= maxdtn*2

         call rzero (wrk ,maxdtn2)
         call rzero (wrk1,maxdtn2)
         call rzero (wrk2,maxdtn2)

         do i= 1,maxdtn
            ssum1= 0.0
            ssum2= 0.0
            jump = (i-1)*maxdtn2
         do j= 1,maxdtn 
            j1= j+ksum1+jump 
            j2= j+ksum2 
            ssumr= wrkdtn2(       j1)
            ssumi= wrkdtn2(maxdtn+j1)
            ssum1= ssum1 + fbin(j2)*ssumr - fbin(maxdtn_tot+j2)*ssumi
            ssum2= ssum2 + fbin(j2)*ssumi + fbin(maxdtn_tot+j2)*ssumr
         enddo
            wrk(       i)=ssum1  
            wrk(maxdtn+i)=ssum2
         enddo
         ksum1= ksum1+maxdtn2*maxdtn
         ksum2= ksum2+maxdtn

         do i= 1,l_dtnface(idir)
            j= dtn_face(i,idir)
            bout(j)     = wrk(       i)
            bout(npts+j)= wrk(maxdtn+i)
         enddo

      endif
      enddo

      return 
      end
c-----------------------------------------------------------------------
      subroutine cem_acoustic_test(ssum,bdfunc)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'ACOUSTIC'

      integer  idtn,idir
      real     bdfunc(1),glsum,ssum
      integer  i,i0
      ssum=0.0

      call cem_acoustic_bd_nmn(rhs,bdfunc)

c      call dssum   (rhs,nx1,ny1,nz1)
c      call col2    (rhs,mask0,npts)   

      do idir= 1,6
      if (l_dtnface(idir).ne.0) then

         do i= 1,l_dtnface(idir)
         i0= dtn_face(i,idir)
         ssum=ssum+rhs(i0)
         enddo

      endif
      enddo

      if (if3d)then
          ssum=ssum/4.0
      else
          ssum=ssum/(2.0*pi)!Fixme here is case sensitive
      endif
      ssum= glsum(ssum,1)   !Fixme here is case senstitive

      return 
      end
c----------------------------------------------------------------------
      subroutine cem_set_dtn_direction  
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer   i,j,k,e,f,mx,mz
      real      tol,tmp,unx0,uny0,unz0
      character CB*3

c...  this routine returns type of direction for dtn boundary
c...  dtn_dir(f,e)=0 : not a DTN boundary                         
c...  dtn_dir(f,e)=1 : negative direction of DTN in x           
c...  dtn_dir(f,e)=2 : positive direction of DTN in x           
c...  dtn_dir(f,e)=3 : negative direction of DTN in y           
c...  dtn_dir(f,e)=4 : positive direction of DTN in y           
c...  dtn_dir(f,e)=5 : negative direction of DTN in z           
c...  dtn_dir(f,e)=6 : positive direction of DTN in z           


c...  setting for direction vectors: dir_v
c...  dtn_dir(f,e)= 1 if same direction with
      dir_v(1,1)=-1.0  
      dir_v(2,1)= 0.0
      dir_v(3,1)= 0.0

c...  dtn_dir(f,e)= 2 if same direction with
      dir_v(1,2)= 1.0
      dir_v(2,2)= 0.0
      dir_v(3,2)= 0.0

c...  dtn_dir(f,e)= 3 if same direction with
      dir_v(1,3)= 0.0
      dir_v(2,3)=-1.0
      dir_v(3,3)= 0.0

c...  dtn_dir(f,e)= 4 if same direction with
      dir_v(1,4)= 0.0
      dir_v(2,4)= 1.0
      dir_v(3,4)= 0.0

c...  dtn_dir(f,e)= 5 if same direction with
      dir_v(1,5)= 0.0
      dir_v(2,5)= 0.0
      dir_v(3,5)=-1.0

c...  dtn_dir(f,e)= 6 if same direction with
      dir_v(1,6)= 0.0
      dir_v(2,6)= 0.0
      dir_v(3,6)= 1.0

      if (nid.eq.0) then
          do i=1,6
             write(6,10) i,dir_v(1,i),dir_v(2,i),dir_v(3,i)
          enddo
      endif
 10   format('  direction:: ',i4,':: (nx,ny,nz)=',3f10.5)
      
c...  assign dtn directions into dtn_dir from 0 to 6:
      mx=nx1/2  ! get rough midpoint index: enough
      mz=1      ! default if not 3d
      if (if3d) mz=nz1/2

      tol=1.e-10 !

      do e = 1,nelt
      do f = 1,nfaces

         CB= CBC(f,e,2)  ! Ed's face numbering
         dtn_dir(f,e)= 0 ! defalut

         if (CB.eq.'DTN') then
            unx0= unx(mx,mz,f,e)
            uny0= uny(mx,mz,f,e)
            unz0= unz(mx,mz,f,e)
            do i= 1,6    !compute inner product: tmp= (dir_v,normal_vector)
               tmp= unx0*dir_v(1,i)+uny0*dir_v(2,i)+unz0*dir_v(3,i)
               if (abs(tmp-1.0).lt.tol) then ! same direction if (tmp.eq.1) 
                   dtn_dir(f,e)= i           ! save the i-direction info
                  if (nid.eq.0) write(6,20) f,e,dtn_dir(f,e)
               endif
            enddo
         endif
      enddo
      enddo

c     do e = 1,nelt
c     do f = 1,nfaces
c        CB= CBC(f,e,2)  ! Ed's face numbering
c        if (CB.eq.'PEC') then
c            if (nid.eq.0) write(6,30) f,e,dtn_dir(f,e)
c        endif
c     enddo
c     enddo

c     do e = 1,nelt
c     do f = 1,nfaces
c        CB= CBC(f,e,2)  ! Ed's face numbering
c        if (CB.eq.'P  ') then
c            if (nid.eq.0) write(6,40) f,e,dtn_dir(f,e)
c        endif
c     enddo
c     enddo

 20   format('DTN: f, e, dtn_dir(f,e)::',3i6)
 30   format('PEC: f, e, dtn_dir(f,e)::',3i6)
 40   format('P  : f, e, dtn_dir(f,e)::',3i6)

      if (nid.eq.0) write(6,*) ' done : cem_set_dtn_direction'

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_dtn_fc_ptr
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer  i,j,k,l,e,f,mx,mz,k1,k2,k3,k4,k5,k6
      integer  idir 
      integer  itmp_old,itmp_new
      integer  icalld
      save     icalld
      data     icalld/0/
      logical  ifnew

c...  set direction vectors and direction index on DTN faces
      call cem_set_dtn_direction


c...  count the index on dtn faces for each direction
      k1= 0
      k2= 0
      k3= 0
      k4= 0
      k5= 0
      k6= 0
      l = 0

      do e = 1,nelt
      do f = 1,nfaces
      do j = 1,nxzf
         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face   numbering
         k = cemface(i)                          ! global numbering
         if     (dtn_dir(f,e).eq.1) then
                 k1=k1+1
                 dtn_face(k1,1) = k 
                 dtn_e   (k1,1) = e 
                 area_dtn(k1,1) = aream(i)
         elseif (dtn_dir(f,e).eq.2) then
                 k2=k2+1
                 dtn_face(k2,2) = k  
                 dtn_e   (k2,2) = e  
                 area_dtn(k2,2) = aream(i)
         elseif (dtn_dir(f,e).eq.3) then
                 k3=k3+1
                 dtn_face(k3,3) = k  
                 dtn_e   (k3,3) = e  
                 area_dtn(k3,3) = aream(i)
         elseif (dtn_dir(f,e).eq.4) then
                 k4=k4+1
                 dtn_face(k4,4) = k  
                 dtn_e   (k4,4) = e  
                 area_dtn(k4,4) = aream(i)
         elseif (dtn_dir(f,e).eq.5) then
                 k5=k5+1
                 dtn_face(k5,5) = k  
                 dtn_e   (k5,5) = e  
                 area_dtn(k5,5) = aream(i)
         elseif (dtn_dir(f,e).eq.6) then
                 k6=k6+1
                 dtn_face(k6,6) = k  
                 dtn_e   (k6,6) = e  
                 area_dtn(k6,6) = aream(i)
         endif

      enddo
      enddo
      enddo
   
      l_dtnface(1) = k1
      l_dtnface(2) = k2
      l_dtnface(3) = k3
      l_dtnface(4) = k4
      l_dtnface(5) = k5
      l_dtnface(6) = k6

      call igop(l_dtnface,n_dtnface,'+  ',6)    
      
      if (nid.eq.0) write(6,10) n_dtnface(1),n_dtnface(2),n_dtnface(3)
     $                         ,n_dtnface(4),n_dtnface(5),n_dtnface(6)
 10   format('  global: n_dtnface(1:6)=',6i6)

      dtn_num= 0
      do idir=1,6
         if (l_dtnface(idir).ne.0) then
             dtn_num=dtn_num+1
         endif
      enddo

      if (nid.eq.0) then
         write(6,*) ' total number of dtn boundary:: dtn_num=',dtn_num
         write(6,*) ' done : cem_set_dtn_fc_ptr'
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine cem_restrict_to_dtn_face(fu,u)
c     Restrict u to faces
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      real     fu(1),u(1)
      integer  idir,i,j,jump

      jump= 0
      do idir= 1,6       
         if (l_dtnface(idir).ne.0) then
          do j= 1,l_dtnface(idir)
             i= dtn_face(j,idir)
             fu(j+jump)= u(i)
          enddo
          jump=jump+l_dtnface(idir)
         endif
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_dtn_quasi_sol(out1,alp,iflag)
c     Restrict u to faces
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      real     out1(1),alp,cos_tmp,sin_tmp
      real     xx,yy,zz,wr,wi
      integer  i,n2,iflag
      if (iflag .eq. 1) then
         do i=1,npts
         xx          =xm1(i,1,1,1)
         yy          =ym1(i,1,1,1)
         zz          =zm1(i,1,1,1)
         cos_tmp     =cos(alp*xx)
         sin_tmp     =sin(alp*xx)
         wr          =out1(i)
         wi          =out1(npts+i)
         out1(     i)=wr*cos_tmp-wi*sin_tmp
         out1(npts+i)=wr*sin_tmp+wi*cos_tmp
         enddo
      else
         do i=1,npts
         xx          =xm1(i,1,1,1)
         yy          =ym1(i,1,1,1)
         zz          =zm1(i,1,1,1)
         cos_tmp     =cos(-alp*xx)
         sin_tmp     =sin(-alp*xx)
         wr          =out1(i)
         wi          =out1(npts+i)
         out1(     i)=wr*cos_tmp-wi*sin_tmp
         out1(npts+i)=wr*sin_tmp+wi*cos_tmp
         enddo
      endif
      return
      end

c-----------------------------------------------------------------------

C------------------------------------------------------------------
      subroutine axhelm_quasi(au,bux,buy,u,helm1,helm2,imesh,isd)
c      subroutine axhelm3(au,u,helm1,helm2,imesh,isd)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
      implicit none
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
      REAL WDDX,WDDYT,WDDZT
      INTEGER NEL,NXY,NYZ,NXZ,NTOT

      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             BUX   (LX1,LY1,LZ1,1)
     $ ,             BUY   (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)

      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)

      real            TMP1, TMP2,TMP3
      real            DUDR, DUDS,DUDT
      
      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      real           term1,term2,dnekclock
      integer  iel,ie,iz,isd,imesh,j,i
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)
C
      NEL = NELT
      NXY = NX1*NY1
      NYZ = NY1*NZ1
      NXZ = NX1*NZ1
      NTOT= NX1*NY1*NZ1*NEL
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
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
c          ENDIF
           
           CALL copy (BUX(1,1,1,IEL),TMP1,NXYZ)
           CALL copy (BUY(1,1,1,IEL),TMP2,NXYZ)

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
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
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
      subroutine acoustic_wght_grad(w1,w2,w3,u1)
      implicit none
c     Weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
c
      integer mode, imode
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real    u1r,u1s,u1t
c
      nn    = nx1-1

      if (if3d) then

         do e=1,nelt
         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
c...        rym1,sym1,tym1,...jac are multiplied. w3m1 is purely weight.
c...        area(i,1,1,1) is multiplied with jacobina weight.
            u1rw = u1r(i)*w3mn(i)
            u1sw = u1s(i)*w3mn(i)
            u1tw = u1t(i)*w3mn(i)
            w1(k) =      u1rw*rxmn(k)
     $                 + u1sw*sxmn(k)
     $                 + u1tw*txmn(k)

            w2(k) =      u1rw*rymn(k)
     $                 + u1sw*symn(k)
     $                 + u1tw*tymn(k)

            w3(k) =      u1rw*rzmn(k)
     $                 + u1sw*szmn(k)
     $                 + u1tw*tzmn(k)

          enddo

          enddo

       else

c        ...2d grad 
          do e= 1,nelt

           j= nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j


             u1rw = u1r(i)*w3mn(i)
             u1sw = u1s(i)*w3mn(i)

c             u1rw = u1r(i) !*w3mn(i)
c             u1sw = u1s(i) !*w3mn(i)

             w1(k) = ( u1rw*rxmn(k)
     $               + u1sw*sxmn(k))!/jacm(k)

             w2(k) =(  u1rw*rymn(k)
     $               + u1sw*symn(k))!/jacm(k)
             w3(k)=w2(k)
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine acoustic_grad(w1,w2,w3,u1)
      implicit none
c     Non-weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
c
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real    u1r,u1s,u1t
c
      nn    = nx1-1

      if (if3d) then

         do e=1,nelt
         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (    u1r(i)*rxmn(k)
     $                 + u1s(i)*sxmn(k)
     $                 + u1t(i)*txmn(k))/jacm(k)

            w2(k) = (    u1r(i)*rymn(k)
     $                 + u1s(i)*symn(k)
     $                 + u1t(i)*tymn(k))/jacm(k)

            w3(k) = (    u1r(i)*rzmn(k)
     $                 + u1s(i)*szmn(k)
     $                 + u1t(i)*tzmn(k))/jacm(k)

          enddo

          enddo

       else

c        ...2d grad 
          do e= 1,nelt

           j= nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j

             w1(k) = ( u1r(i)*rxmn(k)
     $               + u1s(i)*sxmn(k))/jacm(k)

             w2(k) =(  u1r(i)*rymn(k)
     $               + u1s(i)*symn(k))/jacm(k)
             w3(k)=w2(k)
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_nmn_direction  
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer   i,j,k,e,f,mx,mz
      real      tol,tmp,unx0,uny0,unz0
      character CB*3

c...  this routine returns type of direction for dtn boundary
c...  dtn_dir(f,e)=0 : not a DTN boundary                         
c...  dtn_dir(f,e)=1 : negative direction of DTN in x           
c...  dtn_dir(f,e)=2 : positive direction of DTN in x           
c...  dtn_dir(f,e)=3 : negative direction of DTN in y           
c...  dtn_dir(f,e)=4 : positive direction of DTN in y           
c...  dtn_dir(f,e)=5 : negative direction of DTN in z           
c...  dtn_dir(f,e)=6 : positive direction of DTN in z           


c...  setting for direction vectors: dir_v
c...  dtn_dir(f,e)= 1 if same direction with
      dir_v(1,1)=-1.0  
      dir_v(2,1)= 0.0
      dir_v(3,1)= 0.0

c...  dtn_dir(f,e)= 2 if same direction with
      dir_v(1,2)= 1.0
      dir_v(2,2)= 0.0
      dir_v(3,2)= 0.0

c...  dtn_dir(f,e)= 3 if same direction with
      dir_v(1,3)= 0.0
      dir_v(2,3)=-1.0
      dir_v(3,3)= 0.0

c...  dtn_dir(f,e)= 4 if same direction with
      dir_v(1,4)= 0.0
      dir_v(2,4)= 1.0
      dir_v(3,4)= 0.0

c...  dtn_dir(f,e)= 5 if same direction with
      dir_v(1,5)= 0.0
      dir_v(2,5)= 0.0
      dir_v(3,5)=-1.0

c...  dtn_dir(f,e)= 6 if same direction with
      dir_v(1,6)= 0.0
      dir_v(2,6)= 0.0
      dir_v(3,6)= 1.0

      if (nid.eq.0) then
          do i=1,6
             write(6,10) i,dir_v(1,i),dir_v(2,i),dir_v(3,i)
          enddo
      endif
 10   format('  direction:: ',i4,':: (nx,ny,nz)=',3f10.5)
      
c...  assign dtn directions into dtn_dir from 0 to 6:
      mx=nx1/2  ! get rough midpoint index: enough
      mz=1      ! default if not 3d
      if (if3d) mz=nz1/2

      tol=1.e-10 !

      do e = 1,nelt
      do f = 1,nfaces

         CB= CBC(f,e,2)  ! Ed's face numbering
         dtn_dir(f,e)= 0 ! defalut

         if (CB.eq.'N  ') then
            unx0= unx(mx,mz,f,e)
            uny0= uny(mx,mz,f,e)
            unz0= unz(mx,mz,f,e)
            do i= 1,6 ! inner product: tmp= (dir_v,normal_vector)
               tmp= unx0*dir_v(1,i)+uny0*dir_v(2,i)+unz0*dir_v(3,i)
            if (abs(tmp-1.0).lt.tol) then          ! same direction if (tmp.eq.1) 
                dtn_dir(f,e)= i                    ! save the i-direction info
                if (nid.eq.1) write(6,20) f,e,dtn_dir(f,e)
            endif
            enddo
         endif
      enddo
      enddo

c     do e = 1,nelt
c     do f = 1,nfaces
c        CB= CBC(f,e,2)  ! Ed's face numbering
c        if (CB.eq.'PEC') then
c            if (nid.eq.0) write(6,30) f,e,dtn_dir(f,e)
c        endif
c     enddo
c     enddo

c     do e = 1,nelt
c     do f = 1,nfaces
c        CB= CBC(f,e,2)  ! Ed's face numbering
c        if (CB.eq.'P  ') then
c            if (nid.eq.0) write(6,40) f,e,dtn_dir(f,e)
c        endif
c     enddo
c     enddo

 20   format('  N  : f, e, dtn_dir(f,e) ::', 3i8)
 30   format('  PEC: f, e, dtn_dir(f,e) ::', 3i8)
 40   format('  P  : f, e, dtn_dir(f,e) ::', 3i8)

      if (nid.eq.0) write(6,*) ' done : cem_set_nmn_direction'

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_nmn_fc_ptr
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'ACOUSTIC'
      integer  i,j,k,l,e,f,mx,mz,k1,k2,k3,k4,k5,k6
      integer  idir 
      integer  itmp_old,itmp_new
      integer  icalld
      save     icalld
      data     icalld/0/
      logical  ifnew

c...  set direction vectors and direction index on DTN faces
      call cem_set_nmn_direction


c...  count the index on dtn faces for each direction
      k1= 0
      k2= 0
      k3= 0
      k4= 0
      k5= 0
      k6= 0
      l = 0

      do e = 1,nelt
      do f = 1,nfaces
      do j = 1,nxzf
         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face   numbering
         k = cemface(i)                          ! global numbering
         if     (dtn_dir(f,e).eq.1) then
                 k1=k1+1
                 dtn_face(k1,1) = k 
                 dtn_e   (k1,1) = e 
                 area_dtn(k1,1) = aream(i)
         elseif (dtn_dir(f,e).eq.2) then
                 k2=k2+1
                 dtn_face(k2,2) = k  
                 dtn_e   (k2,2) = e  
                 area_dtn(k2,2) = aream(i)
         elseif (dtn_dir(f,e).eq.3) then
                 k3=k3+1
                 dtn_face(k3,3) = k  
                 dtn_e   (k3,3) = e  
                 area_dtn(k3,3) = aream(i)
         elseif (dtn_dir(f,e).eq.4) then
                 k4=k4+1
                 dtn_face(k4,4) = k  
                 dtn_e   (k4,4) = e  
                 area_dtn(k4,4) = aream(i)
         elseif (dtn_dir(f,e).eq.5) then
                 k5=k5+1
                 dtn_face(k5,5) = k  
                 dtn_e   (k5,5) = e  
                 area_dtn(k5,5) = aream(i)
         elseif (dtn_dir(f,e).eq.6) then
                 k6=k6+1
                 dtn_face(k6,6) = k  
                 dtn_e   (k6,6) = e  
                 area_dtn(k6,6) = aream(i)
         endif

      enddo
      enddo
      enddo
   
      l_dtnface(1) = k1
      l_dtnface(2) = k2
      l_dtnface(3) = k3
      l_dtnface(4) = k4
      l_dtnface(5) = k5
      l_dtnface(6) = k6

      call icopy(i_dtnface,l_dtnface,6)
      call igop (i_dtnface,n_dtnface,'+  ',6)    
      
      if (nid.eq.0) write(6,10) n_dtnface(1),n_dtnface(2),n_dtnface(3)
     $                         ,n_dtnface(4),n_dtnface(5),n_dtnface(6)
 10   format('  global: n_dtnface(1:6)=',6i6)

      do idir=1,6
         if (l_dtnface(idir).ne.0) then
             dtn_nid(idir)= 1        !to count the node with DTNs
         else
             dtn_nid(idir)= 0        !to count the node with DTNs
         endif
      enddo

      dtn_num= 0
      do idir=1,6
         if (n_dtnface(idir).ne.0) then
             dtn_num=dtn_num+1
         endif
      enddo

      if (nid.eq.0) then
         write(6,*) ' total number of neumann boundary:: nmn_num=',dtn_num
         write(6,*) ' done : cem_set_nmn_fc_ptr'
      endif

      return
      end
c----------------------------------------------------------------------
c=======================================================================
      subroutine hmh_gmres5(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'

      common /wgmres/ x2(lpts3*2)
     $              , r2(lpts3*2)
     $              , w2(lpts3*2)
     $              , v2(lpts3*2,lgmres)  
      real            x2,r2,w2,v2


      integer  n,n1,n2,outer,isd,imsh,iter,iconv,i,j,k,m
c     real     x2(n),r2(n),w2(n),v2(n,lgmres)
      real     phi(n),res(n),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp,tolpss,glsc3,rnorm
      real*8   etime1,etime_p,divex,dnekclock

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
c     call copy(x2,phi,n) ! none zero inital guess
c      call rone (x2,n)    !zero inital guess
c      iter  = 1
c      call dssum(x2(1),nx1,ny1,nz1)
c      call dssum(x2(n1+1),nx1,ny1,nz1)
c      call col2 (x2,mask,n)
c/block above
      outer = 0
      do while (iconv.eq.0.and.iter.lt.3000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r2,res,n)                  ! r = res
         else
            !update residual
            call copy   (r2,res,n)                  ! r = res
            call axhelm_pec(w2,x2,h1,h2,imsh,isd)       ! w = A x
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

            call axhelm_pec(w2,v2(1,j),h1,h2,imsh,isd)  ! w = A v
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
      subroutine axhelm_pec(au,u,helm1,helm2,imesh,isd)
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
      INCLUDE 'DXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
c      INCLUDE 'CTIMER'
      include 'EMWAVE'
      INCLUDE 'ACOUSTIC'
C
      real au(1),u(1),helm1(1),helm2(1)
      real tmp1(lx1*ly1*lz1*lelt*2)
      real tmp2(lx1*ly1*lz1*lelt*2)
      real wrkdx(lx1*ly1*lz1*lelt*2)
      real wrkdy(lx1*ly1*lz1*lelt*2)
      real wrkdz(lx1*ly1*lz1*lelt*2)
      integer npts2,n2,imesh,isd
      real alp0
      npts2 =npts*2
      n2    =npts+1
      call copy(tmp1,u,npts2)
      call axhelm3(au(1),tmp1(1),helm1(1),helm2(1),imesh,isd)
      call axhelm3(au(n2),tmp1(n2),helm1(n2),helm2(n2),imesh,isd)
c      call acoustic_grad(wrkdx( 1),wrkdy( 1),wrkdz( 1),tmp1( 1))
c      call acoustic_grad(wrkdx(n2),wrkdy(n2),wrkdz(n2),tmp1(n2))
c      call col2    (wrkdx(1),bm1,npts)   ! B*dx
c      call col2    (wrkdx(n2),bm1,npts)   ! B*dx
      call acoustic_wght_grad(wrkdx( 1),wrkdy( 1),wrkdz( 1),tmp1( 1))
      call acoustic_wght_grad(wrkdx(n2),wrkdy(n2),wrkdz(n2),tmp1(n2))

      alp0=2.0*alpha0(1)
      call add2s2  (au( 1),wrkdx(n2),-alp0,npts)              ! Aur = Aur + ui
      call add2s2  (au(n2),wrkdx( 1), alp0,npts)              ! Aui = Aui - ur

      return
      end
