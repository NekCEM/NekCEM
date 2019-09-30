c-----------------------------------------------------------------------
c Below: This multigrid part is from the Version r1808
c-----------------------------------------------------------------------
      subroutine axh(w,x,h1,h2,mask,n)
      include 'SIZE'
      include 'TOTAL'

c
c     w = A*x for pressure iteration
c

      integer n
      real w(n),x(n),h1(n),h2(n),mask(n)

      imsh = 1
      isd  = 1
      call axhelm3(w,x,h1,h2,imsh,isd)
      call dssum  (w,nx1,ny1,nz1)
      call col2   (w,mask,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine hmh_gmres_mg(res,h1,h2,wt,mask,iter)

c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.


      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real    res(1),h1(1),h2(1),wt(1),mask(1)

      common /scrcg/   d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

      common /cgmres1/ y(lgmres)
      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      real alpha, l, temp
      integer outer,icalld
      save          icalld
      data          icalld /0/

      logical iflag
      save    iflag
      data    iflag  /.false. /
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock
      logical ifio

      ifmgrid = .true.
      ifield  = 1

      n = nx1*ny1*nz1*nelt

      if (icalld.eq.0) call setup_mg(h1,h2)
      icalld=1

      etime1 = dnekclock()
      etime_p= 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml,mu,bm1,binvm1,n)
         norm_fac = 1./sqrt(volvm1)
      endif

c     call chktcg1(tolps,res,h1,h2,mask,wt,1,1)
c     if (param(22).gt.0.and.tolps.gt.abs(param(22)))
c    $   tolps = abs(param(22))

      tolps = abs(param(22))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      iconv = 0
      call rzero(x,n)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.100)
         outer = outer+1

         if(iter.eq.0) then                       !      -1
            call copy(r,res,n)                 ! r = L  res
c           call col3(r,ml,res,n)                 ! r = L  res
         else                                     ! update residual
            call copy(r,res,n)                  ! r = res
            call axh(w,x,h1,h2,mask,n)         ! w = A x
            call sub2(r,w,n)                    ! r = r - w
c           call col2(r,ml,n)                     ! r = L   r
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc3(r,r,wt,n))         ! gamma  = \/ (r,r)
                                                  !      1
         if(iter.eq.0) then
            div0 = gamma(1)*norm_fac
            if (param(22).lt.0) tolpss=abs(param(22))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                  !       -1
            call col3(w,mu,v(1,j),n)              ! w  = U   v
                                                  !           j
            call h1mg_solve(z(1,j),w,wt)          ! z  = M   w
c           call copy(z(1,j),w,n)                 ! NO preconditioner

            call axh (w,z(1,j),h1,h2,mask,n)      ! w = A z
                                                  !      -1
            call col2(w,ml,n)                     ! w = L   w

            do i=1,j
               h(i,j)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
            enddo                                 !  i,j       i
            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !          i,j  i

            do i=1,j-1       ! apply Givens rotations to new column
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if ((nid.eq.0).and.(istep.le.5))
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' gmres_mg')


            if (rnorm .lt. tolpss) then
                goto 900  !converged
            endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1
         enddo

  900    iconv = 1

 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,z(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
      enddo
 9000 continue

      divex = rnorm
      call copy(res,x,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1

      if (nid.eq.0) then
        if ((mod(istep,iocomm).eq.0).or.(istep.le.10))
     $  write(6,9999) istep,iter,divex,tolpss,div0,etime_p,etime1
      endif
 9999 format(' ',i9,' gmres_mg: iteration#',i5,1p5e12.4)

      return
      end
c-----------------------------------------------------------------------
c
c  To do:
c
c  1)  Why does hsmg_schwarz_toext2d not zero out a, whereas 3d does??  DONE
c  2)  Convert all nelt refs to nelfld(ifield) or (nelmg?)  DONE
c  3)  Define mg_schwarz_wt for up to and including mg_h1_lmax   DONE
c  4)  MAKE CERTAIN common /hsmgw/ is LARGE enough in hsmg_tnsr and  DONE
c      elsewhere!
c  5)  Devise and implement UNIT tests, now, so that you can test
c      pieces of the setup code in stages.
c  6)  Start developing and testing, in a linear fashion, the SETUP driver.
c  7)  Make certain dssum flags declared for all levels  DONE
c  8)  Need TWO masks for each level:  one for A*x, and one for Schwarz!
c      NO -- one is fine.
c  9)  Modify axml so addition comes after dssum.  DONE
c
c-----------------------------------------------------------------------
c
c Some relevant parameters
c
c param(41):
c     0 - use additive SEMG
c     1 - use hybrid SEMG (not yet working... but coming soon!)
c
c param(42):   navier0.f, fasts.f
c     0 - use GMRES for iterative solver, use non-symmetric weighting
c     1 - use PCG for iterative solver, do not use weighting
c
c param(43):   uzawa_gmres.f, navier6.f
c     0 - use additive multilevel scheme (requires param(42).eq.0)
c     1 - use original 2 level scheme
c
c param(44):   fast3d.f, navier6.f
c     0 - base top-level additive Schwarz on restrictions of E
c     1 - base top-level additive Schwarz on restrictions of A
c
c----------------------------------------------------------------------
      subroutine hsmg_setup(h1,h2)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'SEMHAT'
      include 'TSTEP'

      real h1(1),h2(1)

      mg_fld = 1
      if (ifield.gt.1) mg_fld = 2
      if (ifield.gt.1) nz1=1/(nx1-ny1)

      if (ifield.eq.1) call hsmg_index_0 ! initialize index sets

      call hsmg_setup_mg_nx  ! set nx values for each level of multigrid
      call hsmg_setup_semhat ! set spectral element hat matrices
      call hsmg_setup_intp
      call hsmg_setup_dssum  ! set direct stiffness summation handles
      call hsmg_setup_wtmask ! set restriction weight matrices and bc masks
      call hsmg_setup_fdm    ! set up fast diagonalization method
      call hsmg_setup_schwarz_wt(.false.)
      call hsmg_setup_solve  ! set up the solver
c     call hsmg_setup_dbg

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_semhat
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      include 'SEMHAT'
      integer n,l
c     generate the SEM hat matrices for each level
c     top level
      n = mg_nx(mg_lmax)
      call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
      call copy(mg_zh(1,mg_lmax),zgl,n-1) !top level based on gl points
      mg_nh(mg_lmax)=n-1
      mg_nhz(mg_lmax)=n-1
      if(.not.if3d) mg_nhz(mg_lmax)=1
c     lower levels
      do l=1,mg_lmax-1
         n = mg_nx(l)
         if(n.gt.1) then
            call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
            call copy(mg_ah(1,l),ah,(n+1)*(n+1))
            call copy(mg_bh(1,l),bh,n+1)
            call copy(mg_dh(1,l),dh,(n+1)*(n+1))
            call transpose(mg_dht(1,l),n+1,dh,n+1)
            call copy(mg_zh(1,l),zh,n+1)
         else
            mg_zh(1,l) = -1.
            mg_zh(2,l) =  1.
         endif
         mg_nh(l)=n+1
         mg_nhz(l)=mg_nz(l)+1
      enddo
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_intp
      include 'SIZE'
      include 'HSMG'
      include 'SEMHAT'
      integer l,nf,nc

      do l=1,mg_lmax-1

         nf=mg_nh(l+1)
         nc=mg_nh(l)

!        Standard multigrid coarse-to-fine interpolation
         call hsmg_setup_intpm(
     $           mg_jh(1,l),mg_zh(1,l+1),mg_zh(1,l),nf,nc)
         call transpose(mg_jht(1,l),nc,mg_jh(1,l),nf)

!        Fine-to-coarse interpolation for variable-coefficient operators
         call hsmg_setup_intpm(
     $           mg_jhfc(1,l),mg_zh(1,l),mg_zh(1,l+1),nc,nf)
         call transpose(mg_jhfct(1,l),nf,mg_jhfc(1,l),nc)
c        call outmat(mg_jhfc(1,l),nc,nf,'MG_JHFC',l)

      enddo
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_intpm(jh,zf,zc,nf,nc)
      integer nf,nc
      real jh(nf,nc),zf(1),zc(1)
      include 'SIZE'
      real w(2*lx1+2)
      do i=1,nf
         call fd_weights_full(zf(i),zc,nc-1,1,w)
         do j=1,nc
            jh(i,j)=w(j)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_dssum
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'ZPER'
      parameter (lxyz=(lx1+2)*(ly1+2)*(lz1+2))
      common /c_is1/ glo_num(lxyz*lelv)
      common /ivrtx/ vertex ((2**ldim)*lelt)

      integer*8 glo_num
      integer vertex
      integer nx,ny,nz
      integer l

c     nz1=1/(nx1-ny1)

c     set up direct stiffness summation for each level
      ncrnr = 2**ndim

      if (nid.eq.0) write(6,*) 'hsmg_setup_dssum, ifgtp=',ifgtp
      if (ifgtp) then
          if (nid.eq.0) write(6,*) 'use regular mesh for multigrid'
          call gen_gtp_vertex    (vertex, ncrnr)
      else
          call get_vert
      endif

c++   write(6,*) mg_fld,' mgfld in hsmg_setup_dssum'

      do l=1,mg_lmax-1
         nx=mg_nh(l)
         ny=mg_nh(l)
         nz=mg_nhz(l)
         call setup_ds(mg_gsh_handle(l,mg_fld),nx,ny,nz
     $                ,nelt,nelgv,vertex,glo_num)
         nx=nx+2
         ny=ny+2
         nz=nz+2
         if(.not.if3d) nz=1
         call setup_ds(mg_gsh_schwarz_handle(l,mg_fld),nx,ny,nz
     $                ,nelt,nelgv,vertex,glo_num)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_wtmask
      include 'SIZE'
      include 'HSMG'
      integer i,l
      i = mg_mask_index(mg_lmax,mg_fld-1)
      do l=1,mg_lmax
         mg_rstr_wt_index(l,mg_fld)=i
         mg_mask_index   (l,mg_fld)=i
         i=i+mg_nh(l)*mg_nhz(l)*2*ndim*nelt
         if(i .gt. lmgs*lmg_rwt*2*ldim*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
            call exitt(1)
         endif
         call hsmg_setup_rstr_wt(
     $           mg_rstr_wt(mg_rstr_wt_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
         call hsmg_setup_mask(
     $           mg_mask(mg_mask_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
      enddo
      mg_mask_index(l,mg_fld)=i
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_wtmask
      include 'SIZE'
      include 'HSMG'
      integer i,l
      i = mg_mask_index(mg_lmax,mg_fld-1)
      do l=1,mg_lmax-1
         mg_rstr_wt_index(l,mg_fld)=i
         mg_mask_index   (l,mg_fld)=i
         i=i+mg_nh(l)*mg_nhz(l)*2*ndim*nelt
         if(i .gt. lmgs*lmg_rwt*2*ldim*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
            call exitt(1)
         endif
         call hsmg_setup_rstr_wt(
     $           mg_rstr_wt(mg_rstr_wt_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
         call hsmg_setup_mask(
     $           mg_mask(mg_mask_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
      enddo
      mg_mask_index(l,mg_fld)=i
      end
c----------------------------------------------------------------------
      subroutine hsmg_intp(uf,uc,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      call hsmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_rstr(uc,uf,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      if(l.ne.mg_lmax-1)
     $   call hsmg_do_wt(uf,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
      call hsmg_tnsr(uc,mg_nh(l),uf,mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))
      call hsmg_dssum(uc,l)
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_rstr_no_dssum(uc,uf,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      if(l.ne.mg_lmax-1)
     $   call hsmg_do_wt(uf,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
      call hsmg_tnsr(uc,mg_nh(l),uf,mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_tnsr(v,nv,u,nu,A,At)
c     computes
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u
      include 'SIZE'
      include 'INPUT'
      integer nv,nu
      real v(1),u(1),A(1),At(1)
      if (.not. if3d) then
         call hsmg_tnsr2d(v,nv,u,nu,A,At)
      else
         call hsmg_tnsr3d(v,nv,u,nu,A,At,At)
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_tnsr2d(v,nv,u,nu,A,Bt)
c     computes
c              T
c     v = A u B
      include 'SIZE'
      integer nv,nu
      real v(nv*nv,nelt),u(nu*nu,nelt),A(1),Bt(1)
      common /hsmgw/ work((lx1+2)*(lx1+2))
      integer ie
      do ie=1,nelt
         call mxm(A,nv,u(1,ie),nu,work,nu)
         call mxm(work,nv,Bt,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_tnsr3d(v,nv,u,nu,A,Bt,Ct)
c     computes
c
c     v = [C (x) B (x) A] u
      include 'SIZE'
      integer nv,nu
      real v(nv*nv*nv,nelt),u(nu*nu*nu,nelt),A(1),Bt(1),Ct(1)
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer ie, i
      do ie=1,nelt
         call mxm(A,nv,u(1,ie),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
c     computes
c              T
c     v = A u B
      subroutine hsmg_tnsr2d_el(v,nv,u,nu,A,Bt)
      integer nv,nu
      real v(nv*nv),u(nu*nu),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work((lx1+2)*(lx1+2))
c
      call mxm(A,nv,u,nu,work,nu)
      call mxm(work,nv,Bt,nu,v,nv)
c
      return
      end
c----------------------------------------------------------------------
c     computes
c
c     v = [C (x) B (x) A] u
      subroutine hsmg_tnsr3d_el(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real v(nv*nv*nv),u(nu*nu*nu),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer i
c
      call mxm(A,nv,u,nu,work,nu*nu)
      do i=0,nu-1
         call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
      enddo
      call mxm(work2,nv*nv,Ct,nu,v,nv)
c
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_dssum(u,l)
      include 'SIZE'
      include 'HSMG'
      include 'CTIMER'
      real u(l)
      if (ifsync) call nekgsync()
      etime1=dnekclock()
      call gs_op(mg_gsh_handle(l,mg_fld),u,1,1,0)
      tdadd =tdadd + dnekclock()-etime1

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_dsprod(u,l)
      include 'SIZE'
      include 'HSMG'
      include 'CTIMER'
      real u(l)

      if (ifsync) call nekgsync()

      call gs_op(mg_gsh_handle(l,mg_fld),u,1,2,0)
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_dssum(u,l)
      include 'SIZE'
      include 'HSMG'
      include 'CTIMER'
      real     u(l)

      if (ifsync) call nekgsync()
      etime1=dnekclock()
      i=mg_gsh_schwarz_handle(l,mg_fld)
      call gs_op(mg_gsh_schwarz_handle(l,mg_fld),u,1,1,0)
      tdadd =tdadd + dnekclock()-etime1

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_extrude(arr1,l1,f1,arr2,l2,f2,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer l1,l2,nx,ny,nz
      real arr1(nx,ny,nz,nelt),arr2(nx,ny,nz,nelt)
      real f1,f2

      integer i,j,k,ie,i0,i1
      i0=2
      i1=nx-1

      if(.not.if3d) then
         do ie=1,nelt
            do j=i0,i1
               arr1(l1+1 ,j,1,ie) = f1*arr1(l1+1 ,j,1,ie)
     $                             +f2*arr2(l2+1 ,j,1,ie)
               arr1(nx-l1,j,1,ie) = f1*arr1(nx-l1,j,1,ie)
     $                             +f2*arr2(nx-l2,j,1,ie)
            enddo
            do i=i0,i1
               arr1(i,l1+1 ,1,ie) = f1*arr1(i,l1+1 ,1,ie)
     $                             +f2*arr2(i,l2+1 ,1,ie)
               arr1(i,ny-l1,1,ie) = f1*arr1(i,ny-l1,1,ie)
     $                             +f2*arr2(i,nx-l2,1,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=i0,i1
            do j=i0,i1
               arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie)
     $                             +f2*arr2(l2+1 ,j,k,ie)
               arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie)
     $                             +f2*arr2(nx-l2,j,k,ie)
            enddo
            enddo
            do k=i0,i1
            do i=i0,i1
               arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie)
     $                             +f2*arr2(i,l2+1 ,k,ie)
               arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie)
     $                             +f2*arr2(i,nx-l2,k,ie)
            enddo
            enddo
            do j=i0,i1
            do i=i0,i1
               arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie)
     $                             +f2*arr2(i,j,l2+1 ,ie)
               arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie)
     $                             +f2*arr2(i,j,nx-l2,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_schwarz(e,r,sigma,l)
      include 'SIZE'
      include 'HSMG'

      real e(l),r(l)

      n = mg_h1_n(l,mg_fld)

      call h1mg_schwarz_part1 (e,r,l)
      call hsmg_schwarz_wt    (e,l)          ! e  := W e
      call cmult              (e,sigma,n)    !  l       l

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_schwarz_part1 (e,r,l)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield
      include 'HSMG'

      real e(l),r(l)

      integer enx,eny,enz,pm

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)
      pm = p_mg_msk(l,mg_fld)

      call h1mg_mask  (r,mg_imask(pm),nelfld(ifield))  ! Zero Dirichlet nodes

      if (if3d) then ! extended array
         call hsmg_schwarz_toext3d(mg_work,r,mg_nh(l))
      else
         call hsmg_schwarz_toext2d(mg_work,r,mg_nh(l))
      endif

      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      if(.not.if3d) enz=1
      i = enx*eny*enz*nelt+1

c     exchange interior nodes
      call hsmg_extrude(mg_work,0,zero,mg_work,2,one,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work,l)
      call hsmg_extrude(mg_work,0,one ,mg_work,2,onem,enx,eny,enz)

      call hsmg_fdm(mg_work(i),mg_work,l) ! Do the local solves

c     Sum overlap region (border excluded)
      call hsmg_extrude(mg_work,0,zero,mg_work(i),0,one ,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work(i),l)
      call hsmg_extrude(mg_work(i),0,one ,mg_work,0,onem,enx,eny,enz)
      call hsmg_extrude(mg_work(i),2,one,mg_work(i),0,one,enx,eny,enz)

      if(.not.if3d) then ! Go back to regular size array
         call hsmg_schwarz_toreg2d(e,mg_work(i),mg_nh(l))
      else
         call hsmg_schwarz_toreg3d(e,mg_work(i),mg_nh(l))
      endif

      call hsmg_dssum(e,l)                           ! sum border nodes
      call h1mg_mask (e,mg_imask(pm),nelfld(ifield)) ! apply mask

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz(e,r,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      real e(l),r(l)
      integer l
      integer enx,eny,enz
      integer i

      real zero,one,onem
      zero =  0
      one  =  1
      onem = -1

c     apply mask (zeros Dirichlet nodes)
      !!!!! uncommenting
      call hsmg_do_wt(r,mg_mask(mg_mask_index(l,mg_fld)),
     $                mg_nh(l),mg_nh(l),mg_nhz(l))

c     go to extended size array (room for overlap)
      if (if3d) then
         call hsmg_schwarz_toext3d(mg_work,r,mg_nh(l))
      else
         call hsmg_schwarz_toext2d(mg_work,r,mg_nh(l))
      endif

      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      if(.not.if3d) enz=1
      i = enx*eny*enz*nelt+1

c     exchange interior nodes
      call hsmg_extrude(mg_work,0,zero,mg_work,2,one,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work,l)
      call hsmg_extrude(mg_work,0,one ,mg_work,2,onem,enx,eny,enz)

c     do the local solves
      call hsmg_fdm(mg_work(i),mg_work,l)
c     sum overlap region (border excluded)
      call hsmg_extrude(mg_work,0,zero,mg_work(i),0,one ,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work(i),l)
      call hsmg_extrude(mg_work(i),0,one ,mg_work,0,onem,enx,eny,enz)
      call hsmg_extrude(mg_work(i),2,one,mg_work(i),0,one,enx,eny,enz)
c     go back to regular size array
      if(.not.if3d) then
         call hsmg_schwarz_toreg2d(e,mg_work(i),mg_nh(l))
      else
         call hsmg_schwarz_toreg3d(e,mg_work(i),mg_nh(l))
      endif
c     sum border nodes
      call hsmg_dssum(e,l)
c     apply mask (zeros Dirichlet nodes)
      !!!!!! changing r to e
      call hsmg_do_wt(e,mg_mask(mg_mask_index(l,mg_fld)),
     $                mg_nh(l),mg_nh(l),mg_nhz(l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toext2d(a,b,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,nelt),b(n,n,nelt)

      integer i,j,ie
c      call rzero(a,(n+2)*(n+2)*nelt)
      do ie=1,nelt
         do i=0,n+1
            a(i,0,ie)=0.
         enddo
         do j=1,n
            a(0  ,j,ie)=0.
            do i=1,n
               a(i,j,ie)=b(i,j,ie)
            enddo
            a(n+1,j,ie)=0.
         enddo
         do i=0,n+1
            a(i,n+1,ie)=0.
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toext3d(a,b,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelt),b(n,n,n,nelt)

      integer i,j,k,ie
      call rzero(a,(n+2)*(n+2)*(n+2)*nelt)
      do ie=1,nelt
      do k=1,n
      do j=1,n
      do i=1,n
         a(i,j,k,ie)=b(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toreg2d(b,a,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,nelt),b(n,n,nelt)

      integer i,j,ie
      do ie=1,nelt
      do j=1,n
      do i=1,n
         b(i,j,ie)=a(i,j,ie)
      enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toreg3d(b,a,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelt),b(n,n,n,nelt)

      integer i,j,k,ie
      do ie=1,nelt
      do k=1,n
      do j=1,n
      do i=1,n
         b(i,j,k,ie)=a(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_fdm()
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'

      integer l,i,j,nl
      i = mg_fast_s_index(mg_lmax,mg_fld-1)
      j = mg_fast_d_index(mg_lmax,mg_fld-1)
      do l=2,mg_lmax
         mg_fast_s_index(l,mg_fld)=i
         nl = mg_nh(l)+2
         i=i+nl*nl*2*ndim*nelt
         if(i .gt. lmg_fasts*2*ldim*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'lmg_fasts too small',i,itmp,lmg_fasts,l
            call exitt(1)
         endif
         mg_fast_d_index(l,mg_fld)=j
         j=j+(nl**ndim)*nelt
         if(j .gt. lmg_fastd*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'lmg_fastd too small',i,itmp,lmg_fastd,l
            call exitt(1)
         endif
         call hsmg_setup_fast(
     $             mg_fast_s(mg_fast_s_index(l,mg_fld))
     $            ,mg_fast_d(mg_fast_d_index(l,mg_fld))
     $            ,mg_nh(l)+2,mg_ah(1,l),mg_bh(1,l),mg_nx(l))
      enddo
      mg_fast_s_index(l,mg_fld)=i
      mg_fast_d_index(l,mg_fld)=j
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fdm()
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'

      integer l,i,j,nl
      i = mg_fast_s_index(mg_lmax,mg_fld-1)
      j = mg_fast_d_index(mg_lmax,mg_fld-1)
      do l=2,mg_lmax-1
         mg_fast_s_index(l,mg_fld)=i
         nl = mg_nh(l)+2
         i=i+nl*nl*2*ndim*nelt
         if(i .gt. lmg_fasts*2*ldim*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'lmg_fasts too small',i,itmp,lmg_fasts,l
            call exitt(1)
         endif
         mg_fast_d_index(l,mg_fld)=j
         j=j+(nl**ndim)*nelt
         if(j .gt. lmg_fastd*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'lmg_fastd too small',i,itmp,lmg_fastd,l
            call exitt(1)
         endif
         call hsmg_setup_fast(
     $             mg_fast_s(mg_fast_s_index(l,mg_fld))
     $            ,mg_fast_d(mg_fast_d_index(l,mg_fld))
     $            ,mg_nh(l)+2,mg_ah(1,l),mg_bh(1,l),mg_nx(l))
      enddo
      mg_fast_s_index(l,mg_fld)=i
      mg_fast_d_index(l,mg_fld)=j
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast(s,d,nl,ah,bh,n)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      real s(nl*nl,2,ndim,nelt)
      real d(nl**ndim,nelt)
      real ah(1),bh(1)
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt

      integer i,j,k
      integer ie,il,nr,ns,nt
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
      real eps,diag

      two  = 2
      ierr = 0
      do ie=1,nelt
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
         nr=nl
         ns=nl
         nt=nl
         call hsmg_setup_fast1d(s(1,1,1,ie),lr,nr,lbr,rbr
     $            ,llr(ie),lmr(ie),lrr(ie),ah,bh,n,ie)
         call hsmg_setup_fast1d(s(1,1,2,ie),ls,ns,lbs,rbs
     $            ,lls(ie),lms(ie),lrs(ie),ah,bh,n,ie)
         if(if3d) call hsmg_setup_fast1d(s(1,1,3,ie),lt,nt,lbt,rbt
     $                     ,llt(ie),lmt(ie),lrt(ie),ah,bh,n,ie)
         il=1
         if(.not.if3d) then
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            do j=1,ns
            do i=1,nr
               diag = lr(i)+ls(j)
               if (diag.gt.eps) then
                  d(il,ie) = 1.0/diag
               else
c                 write(6,2) ie,'Reset Eig in hsmg setup fast:',i,j,l
c    $                         ,eps,diag,lr(i),ls(j)
    2             format(i6,1x,a21,3i5,1p4e12.4)
                  d(il,ie) = 0.0
               endif
               il=il+1
            enddo
            enddo
         else
            eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $                  + vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            do k=1,nt
            do j=1,ns
            do i=1,nr
               diag = lr(i)+ls(j)+lt(k)
               if (diag.gt.eps) then
                  d(il,ie) = 1.0/diag
               else
c                 write(6,3) ie,'Reset Eig in hsmg setup fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
    3             format(i6,1x,a21,4i5,1p5e12.4)
                  d(il,ie) = 0.0
               endif
               il=il+1
            enddo
            enddo
            enddo
         endif
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('A INVALID BC FOUND in genfast$',ierrmx,1)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d(s,lam,nl,lbc,rbc,ll,lm,lr,ah,bh,n,ie)
      integer nl,lbc,rbc,n
      real s(nl,nl,2),lam(nl),ll,lm,lr
      real ah(0:n,0:n),bh(0:n)

      include 'SIZE'
      parameter(lxm=lx1+2)
      common /ctmp0/ b(2*lxm*lxm),w(2*lxm*lxm)

      call hsmg_setup_fast1d_a(s,lbc,rbc,ll,lm,lr,ah,n)
      call hsmg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)

c     if (nid.eq.0) write(6,*) 'THIS is generalev call',nl,lbc
      call generalev(s,b,lam,nl,w)
      if(lbc.gt.0) call row_zero(s,nl,nl,1)
      if(lbc.eq.1) call row_zero(s,nl,nl,2)
      if(rbc.gt.0) call row_zero(s,nl,nl,nl)
      if(rbc.eq.1) call row_zero(s,nl,nl,nl-1)

      call transpose(s(1,1,2),nl,s,nl)
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d_a(a,lbc,rbc,ll,lm,lr,ah,n)
      integer lbc,rbc,n
      real a(0:n+2,0:n+2),ll,lm,lr
      real ah(0:n,0:n)

      real fac
      integer i,j,i0,i1
      i0=0
      if(lbc.eq.1) i0=1
      i1=n
      if(rbc.eq.1) i1=n-1

      call rzero(a,(n+3)*(n+3))
      fac = 2.0/lm
      a(1,1)=1.0
      a(n+1,n+1)=1.0
      do j=i0,i1
         do i=i0,i1
            a(i+1,j+1)=fac*ah(i,j)
         enddo
      enddo
      if(lbc.eq.0) then
         fac = 2.0/ll
         a(0,0)=fac*ah(n-1,n-1)
         a(1,0)=fac*ah(n  ,n-1)
         a(0,1)=fac*ah(n-1,n  )
         a(1,1)=a(1,1)+fac*ah(n  ,n  )
      else
         a(0,0)=1.0
      endif
      if(rbc.eq.0) then
         fac = 2.0/lr
         a(n+1,n+1)=a(n+1,n+1)+fac*ah(0,0)
         a(n+2,n+1)=fac*ah(1,0)
         a(n+1,n+2)=fac*ah(0,1)
         a(n+2,n+2)=fac*ah(1,1)
      else
         a(n+2,n+2)=1.0
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
      integer lbc,rbc,n
      real b(0:n+2,0:n+2),ll,lm,lr
      real bh(0:n)

      real fac
      integer i,j,i0,i1
      i0=0
      if(lbc.eq.1) i0=1
      i1=n
      if(rbc.eq.1) i1=n-1

      call rzero(b,(n+3)*(n+3))
      fac = 0.5*lm
      b(1,1)=1.0
      b(n+1,n+1)=1.0
      do i=i0,i1
         b(i+1,i+1)=fac*bh(i)
      enddo
      if(lbc.eq.0) then
         fac = 0.5*ll
         b(0,0)=fac*bh(n-1)
         b(1,1)=b(1,1)+fac*bh(n  )
      else
         b(0,0)=1.0
      endif
      if(rbc.eq.0) then
         fac = 0.5*lr
         b(n+1,n+1)=b(n+1,n+1)+fac*bh(0)
         b(n+2,n+2)=fac*bh(1)
      else
         b(n+2,n+2)=1.0
      endif
      return
      end
c----------------------------------------------------------------------
c     clobbers r
      subroutine hsmg_fdm(e,r,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      call hsmg_do_fast(e,r,
     $      mg_fast_s(mg_fast_s_index(l,mg_fld)),
     $      mg_fast_d(mg_fast_d_index(l,mg_fld)),
     $      mg_nh(l)+2)
      return
      end
c----------------------------------------------------------------------
c     clobbers r
      subroutine hsmg_do_fast(e,r,s,d,nl)
      include 'SIZE'
      include 'INPUT'
      real e(nl**ndim,nelt)
      real r(nl**ndim,nelt)
      real s(nl*nl,2,ndim,nelt)
      real d(nl**ndim,nelt)

      integer ie,nn,i
      nn=nl**ndim
      if(.not.if3d) then
         do ie=1,nelt
            call hsmg_tnsr2d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,2,1,ie),s(1,1,2,ie))
            do i=1,nn
               r(i,ie)=d(i,ie)*e(i,ie)
            enddo
            call hsmg_tnsr2d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,1,1,ie),s(1,2,2,ie))
         enddo
      else
         do ie=1,nelt
            call hsmg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,2,1,ie),s(1,1,2,ie),s(1,1,3,ie))
            do i=1,nn
               r(i,ie)=d(i,ie)*e(i,ie)
            enddo
            call hsmg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,1,1,ie),s(1,2,2,ie),s(1,2,3,ie))
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
c     u = wt .* u
      subroutine hsmg_do_wt(u,wt,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz
      real u(nx,ny,nz,nelt)
      real wt(nx,nz,2,ndim,nelt)

      integer e

      if (.not. if3d) then
         do ie=1,nelt
            do j=1,ny
               u( 1,j,1,ie)=u( 1,j,1,ie)*wt(j,1,1,1,ie)
               u(nx,j,1,ie)=u(nx,j,1,ie)*wt(j,1,2,1,ie)
            enddo
            do i=2,nx-1
               u(i, 1,1,ie)=u(i, 1,1,ie)*wt(i,1,1,2,ie)
               u(i,ny,1,ie)=u(i,ny,1,ie)*wt(i,1,2,2,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=1,nz
            do j=1,ny
               u( 1,j,k,ie)=u( 1,j,k,ie)*wt(j,k,1,1,ie)
               u(nx,j,k,ie)=u(nx,j,k,ie)*wt(j,k,2,1,ie)
            enddo
            enddo
            do k=1,nz
            do i=2,nx-1
               u(i, 1,k,ie)=u(i, 1,k,ie)*wt(i,k,1,2,ie)
               u(i,ny,k,ie)=u(i,ny,k,ie)*wt(i,k,2,2,ie)
            enddo
            enddo
            do j=2,ny-1
            do i=2,nx-1
               u(i,j, 1,ie)=u(i,j, 1,ie)*wt(i,j,1,3,ie)
               u(i,j,nz,ie)=u(i,j,nz,ie)*wt(i,j,2,3,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_rstr_wt(wt,nx,ny,nz,l,w)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelt)
      real wt(nx,nz,2,ndim,nelt)

      integer ie
      !init border nodes to 1
      call rzero(w,nx*ny*nz*nelt)
c     print *, 'Setup rstr wt: ',nx,ny,nz,nelt
      if (.not.if3d) then
         do ie=1,nelt
            do i=1,nx
               w(i,1,1,ie)=1.0
               w(i,ny,1,ie)=1.0
            enddo
            do j=1,ny
               w(1,j,1,ie)=1.0
               w(nx,j,1,ie)=1.0
            enddo
         enddo
      else
         do ie=1,nelt
            do j=1,ny
            do i=1,nx
               w(i,j,1,ie)=1.0
               w(i,j,nz,ie)=1.0
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               w(i,1,k,ie)=1.0
               w(i,ny,k,ie)=1.0
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               w(1,j,k,ie)=1.0
               w(nx,j,k,ie)=1.0
            enddo
            enddo
         enddo
      endif
      call hsmg_dssum(w,l)
      !invert the count w to get the weight wt
      if (.not. if3d) then
         do ie=1,nelt
            do j=1,ny
               wt(j,1,1,1,ie)=1.0/w(1,j,1,ie)
               wt(j,1,2,1,ie)=1.0/w(nx,j,1,ie)
            enddo
            do i=1,nx
               wt(i,1,1,2,ie)=1.0/w(i,1,1,ie)
               wt(i,1,2,2,ie)=1.0/w(i,ny,1,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=1,nz
            do j=1,ny
               wt(j,k,1,1,ie)=1.0/w(1,j,k,ie)
               wt(j,k,2,1,ie)=1.0/w(nx,j,k,ie)
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               wt(i,k,1,2,ie)=1.0/w(i,1,k,ie)
               wt(i,k,2,2,ie)=1.0/w(i,ny,k,ie)
            enddo
            enddo
            do j=1,ny
            do i=1,nx
               wt(i,j,1,3,ie)=1.0/w(i,j,1,ie)
               wt(i,j,2,3,ie)=1.0/w(i,j,nz,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_mask(wt,nx,ny,nz,l,w)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelt)
      real wt(nx,nz,2,ndim,nelt)

      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt,two
c     init everything to 1

      n = nx*ny*nz*nelt
      call rone(w,n)

c     set dirichlet nodes to zero
      ierr = 0
      two  = 2
      do ie=1,nelt
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
         if (ierr.ne.0) then
            ierr = -1
            call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,two,ierr)
         endif

         if(lbr.eq.1) then
            do k=1,nz
            do j=1,ny
               w(1,j,k,ie)=0.0
            enddo
            enddo
         endif
         if(rbr.eq.1) then
            do k=1,nz
            do j=1,ny
               w(nx,j,k,ie)=0.0
            enddo
            enddo
         endif
         if(lbs.eq.1) then
            do k=1,nz
            do i=1,nx
               w(i,1,k,ie)=0.0
            enddo
            enddo
         endif
         if(rbs.eq.1) then
            do k=1,nz
            do i=1,nx
               w(i,ny,k,ie)=0.0
            enddo
            enddo
         endif
         if(if3d) then
            if(lbt.eq.1) then
               do j=1,ny
               do i=1,nx
                  w(i,j,1,ie)=0.0
               enddo
               enddo
            endif
            if(rbt.eq.1) then
               do j=1,ny
               do i=1,nx
                  w(i,j,nz,ie)=0.0
               enddo
               enddo
            endif
         endif
      enddo
c     do direct stiffness multiply

      call hsmg_dsprod(w,l)


c     store weight
      if (.not. if3d) then
         do ie=1,nelt
            do j=1,ny
               wt(j,1,1,1,ie)=w(1,j,1,ie)
               wt(j,1,2,1,ie)=w(nx,j,1,ie)
            enddo
            do i=1,nx
               wt(i,1,1,2,ie)=w(i,1,1,ie)
               wt(i,1,2,2,ie)=w(i,ny,1,ie)
            enddo
         enddo
      else
         do ie=1,nelt
            do k=1,nz
            do j=1,ny
               wt(j,k,1,1,ie)=w(1,j,k,ie)
               wt(j,k,2,1,ie)=w(nx,j,k,ie)
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               wt(i,k,1,2,ie)=w(i,1,k,ie)
               wt(i,k,2,2,ie)=w(i,ny,k,ie)
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               wt(j,k,1,3,ie)=w(i,j,1,ie)
               wt(j,k,2,3,ie)=w(i,j,nz,ie)
            enddo
            enddo
         enddo
      endif

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('B INVALID BC FOUND in genfast$',ierrmx,1)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_schwarz_wt(ifsqrt)
      logical ifsqrt
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'

      integer l,i,nl,nlz

      i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)
      do l=2,mg_lmax-1
         mg_schwarz_wt_index(l,mg_fld)=i
         nl = mg_nh(l)
         nlz = mg_nh(l)
         if(.not.if3d) nlz=1
         i=i+nl*nlz*4*ndim*nelt
         if(i .gt. lmg_swt*4*ldim*lelv) then
            itmp = i/(4*ldim*lelv)
            write(6,*) 'lmg_swt too small',i,itmp,lmg_swt,l
            call exitt(1)
         endif

         call h1mg_setup_schwarz_wt_1(
     $      mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),l,ifsqrt)
      enddo
      mg_schwarz_wt_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt(ifsqrt)
      logical ifsqrt
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'

      integer l,i,nl,nlz

      i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)
      do l=2,mg_lmax

         mg_schwarz_wt_index(l,mg_fld)=i
         nl  = mg_nh(l)
         nlz = mg_nhz(l)
         i   = i+nl*nlz*4*ndim*nelt

         if (i .gt. lmg_swt*4*ldim*lelv) then
            itmp = i/(4*ldim*lelv)
            write(6,*) 'lmg_swt too small',i,itmp,lmg_swt,l
            call exitt(1)
         endif

         call h1mg_setup_schwarz_wt_1(
     $      mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),l,ifsqrt)

      enddo

      mg_schwarz_wt_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_wt(e,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'

      if(.not.if3d) call hsmg_schwarz_wt2d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      if(if3d) call hsmg_schwarz_wt3d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_wt2d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,nelt)
      real wt(n,4,2,nelt)

      integer ie,i,j
      do ie=1,nelt
         do j=1,n
            e(1  ,j,ie)=e(1  ,j,ie)*wt(j,1,1,ie)
            e(2  ,j,ie)=e(2  ,j,ie)*wt(j,2,1,ie)
            e(n-1,j,ie)=e(n-1,j,ie)*wt(j,3,1,ie)
            e(n  ,j,ie)=e(n  ,j,ie)*wt(j,4,1,ie)
         enddo
         do i=3,n-2
            e(i,1  ,ie)=e(i,1  ,ie)*wt(i,1,2,ie)
            e(i,2  ,ie)=e(i,2  ,ie)*wt(i,2,2,ie)
            e(i,n-1,ie)=e(i,n-1,ie)*wt(i,3,2,ie)
            e(i,n  ,ie)=e(i,n  ,ie)*wt(i,4,2,ie)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_wt3d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,n,nelt)
      real wt(n,n,4,3,nelt)

      integer ie,i,j,k
      do ie=1,nelt
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*wt(j,k,1,1,ie)
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*wt(j,k,2,1,ie)
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*wt(j,k,3,1,ie)
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*wt(j,k,4,1,ie)
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*wt(i,k,1,2,ie)
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*wt(i,k,2,2,ie)
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*wt(i,k,3,2,ie)
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*wt(i,k,4,2,ie)
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*wt(i,j,1,3,ie)
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*wt(i,j,2,3,ie)
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*wt(i,j,3,3,ie)
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*wt(i,j,4,3,ie)
         enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_coarse_solve(e,r)
      include 'SIZE'
      include 'DOMAIN'
      include 'ESOLV'
      include 'GEOM'
      include 'SOLN'
      include 'PARALLEL'
      include 'HSMG'
      include 'CTIMER'
      include 'INPUT'
      include 'TSTEP'
      real e(1),r(1)
c
      integer n_crs_tot
      save    n_crs_tot
      data    n_crs_tot /0/
c
      if (icalld.eq.0) then ! timer info
         ncrsl=0
         tcrsl=0.0
      endif
      icalld = 1

      if (ifsync) call nekgsync()

      ncrsl  = ncrsl  + 1
      etime1=dnekclock()

      call crs_solve(xxth(ifield),e,r)

      nz=nelv*(2**ndim)*(2**ndim)
c     write(6,*) (r(k),k=1,nz),' CRS R',ifield
c     write(6,*) (e(k),k=1,nz),' CRS E',ifield

c     call rzero(e,nz)

      tcrsl=tcrsl+dnekclock()-etime1

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_solve
      include 'SIZE'
      include 'HSMG'

      integer l,i,nl,nlz
      i = mg_solve_index(mg_lmax+1,mg_fld-1)
      do l=1,mg_lmax
         mg_solve_index(l,mg_fld)=i
         i=i+mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelt
         if(i .gt. lmg_solve*lelv) then
            itmp = i/lelv
            write(6,*) 'lmg_solve too small',i,itmp,lmg_solve,l
            call exitt(1)
         endif
      enddo
      mg_solve_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_solve(e,r)
      real e(1),r(1)
      include 'SIZE'
      include 'HSMG'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'

      integer ilstep,iter
      save    ilstep,iter
      data    ilstep,iter /0,0/

      real    rhoavg,copt(2),copw(2)
      save    rhoavg,copt1,copt2
      data    rhoavg,copt1,copt2 /3*1./  ! Default copt = 1 for additive

      integer l,nt
      integer*8 ntotg,nxyz2


      mg_fld = 1
      if (ifield.gt.1) mg_fld = 2

      if (istep.ne.ilstep) then
         ilstep = istep
         ntot1  = nx1*ny1*nz1*nelt
         rhoavg = glsc2(vtrans,bm1,ntot1)/volvm1
      endif

      n = nx2*ny2*nz2*nelt
c     call copy(e,r,n)
c     return

      if (icalld.eq.0) then

         tddsl=0.0
         nddsl=0

         icalld = 1
         taaaa = 0
         tbbbb = 0
         tcccc = 0
         tdddd = 0
         teeee = 0
      endif

      nddsl  = nddsl  + 1
#ifndef NOTIMER
      etime1 = dnekclock()
#endif

c     n = nx2*ny2*nz2*nelt
c     rmax = glmax(r,n)
c     if (nid.eq.0) write(6,*) istep,n,rmax,' rmax1'

      iter = iter + 1

      l = mg_lmax
      nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelt
      ! e := W M        r
      !         Schwarz
      time_0 = dnekclock()
      call local_solves_fdm(e,r)

      time_1 = dnekclock()
      ! Additive
      ! w := r - w
      do i = 1,nt
         mg_work2(i) = r(i)
      enddo
      time_2 = dnekclock()

      do l = mg_lmax-1,2,-1

c        rmax = glmax(mg_work2,nt)
c        if (nid.eq.0) write(6,*) l,nt,rmax,' rmax2'

         nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelt
         !          T
         ! r   :=  J w
         !  l
         call hsmg_rstr(mg_solve_r(mg_solve_index(l,mg_fld)),mg_work2,l)

         ! w  := r
         !        l
         call copy(mg_work2,mg_solve_r(mg_solve_index(l,mg_fld)),nt)
         ! e  := M        w
         !  l     Schwarz
         call hsmg_schwarz(
     $          mg_solve_e(mg_solve_index(l,mg_fld)),mg_work2,l)

         ! e  := W e
         !  l       l
         call hsmg_schwarz_wt(mg_solve_e(mg_solve_index(l,mg_fld)),l)

c        call exitti('quit in mg$',l,1)

         ! w  := r  - w
         !        l
         do i = 0,nt-1
            mg_work2(i+1) = mg_solve_r(mg_solve_index(l,mg_fld)+i)
     $         !-alpha*mg_work2(i+1)
         enddo
      enddo

      call hsmg_rstr_no_dssum(
     $   mg_solve_r(mg_solve_index(1,mg_fld)),mg_work2,1)

      nzw = ndim-1

      call hsmg_do_wt(mg_solve_r(mg_solve_index(1,mg_fld)),
     $                mg_mask(mg_mask_index(1,mg_fld)),2,2,nzw)

      !        -1
      ! e  := A   r
      !  1         1
      call hsmg_coarse_solve(mg_solve_e(mg_solve_index(1,mg_fld)),
     $                       mg_solve_r(mg_solve_index(1,mg_fld)))

      call hsmg_do_wt(mg_solve_e(mg_solve_index(1,mg_fld)),
     $                mg_mask(mg_mask_index(1,mg_fld)),2,2,nzw)
      time_3 = dnekclock()
      do l = 2,mg_lmax-1
         nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelt
         ! w   :=  J e
         !            l-1
         call hsmg_intp
     $      (mg_work2,mg_solve_e(mg_solve_index(l-1,mg_fld)),l-1)

         ! e   :=  e  + w
         !  l       l
         do i = 0,nt-1
            mg_solve_e(mg_solve_index(l,mg_fld)+i) =
     $        + mg_solve_e(mg_solve_index(l,mg_fld)+i) + mg_work2(i+1)
         enddo
      enddo
      l = mg_lmax
      nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelt
      ! w   :=  J e
      !            m-1

      call hsmg_intp(mg_work2,
     $   mg_solve_e(mg_solve_index(l-1,mg_fld)),l-1)

      ! e := e + w

      do i = 1,nt
         e(i) = e(i) + copt2*mg_work2(i)
      enddo
      time_4 = dnekclock()
      print *, 'Did an MG iteration'

      taaaa = taaaa + (time_1 - time_0)
      tbbbb = tbbbb + (time_2 - time_1)
      tcccc = tcccc + (time_3 - time_2)
      tdddd = tdddd + (time_4 - time_3)
      teeee = teeee + (time_4 - time_0)
c
c     A typical time breakdown:
c
c  1.3540E+01  5.4390E+01  1.1440E+01  1.2199E+00  8.0590E+01 HSMG time
c
c  ==>  54/80 = 67 % of preconditioner time is in residual evaluation!
c
c     call ortho (e)

#ifndef NOTIMER
      tddsl  = tddsl + ( dnekclock()-etime1 )
#endif


      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_mg_nx()
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'SEMHAT'
      real w(lx1+2)
      integer nf,nc,nr
      integer nx,ny,nz

      integer mgn2(10)
      save    mgn2
      data    mgn2 / 1, 2, 2, 2, 2, 3, 3, 5, 5, 5/
c     data    mgn2 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 0

c     if (param(82).eq.0) param(82)=2  ! nek default
c     if (np.eq.1)        param(82)=2  ! single proc. too slow
      p82 = 2                          ! potentially variable nxc

c     mg_lmax = 4
c     mg_lmax = 3
c     mg_lmax = 2
      mg_lmax = 3

c     uuuu

      if (lx1.eq.4) mg_lmax = 2

      mgnx1    = p82-1 !1
      mg_nx(1) = mgnx1
      mg_ny(1) = mgnx1
      mg_nz(1) = mgnx1
      if (.not.if3d) mg_nz(1) = 0

      mgnx2 = 2*(lx2/4) + 1
      if (lx1.eq.5)  mgnx2 = 3
c     if (lx1.eq.6)  mgnx2 = 3
      if (lx1.le.10) mgnx2 = mgn2(nx1)
      if (lx1.eq.8)  mgnx2 = 4
      if (lx1.eq.8)  mgnx2 = 3

c     mgnx2 = min(3,mgnx2)


      mg_nx(2) = mgnx2
      mg_ny(2) = mgnx2
      mg_nz(2) = mgnx2
      if (.not.if3d) mg_nz(2) = 0

      mg_nx(3) = mgnx2+1
      mg_ny(3) = mgnx2+1
      mg_nz(3) = mgnx2+1
      if (.not.if3d) mg_nz(3) = 0

      mg_nx(mg_lmax) = lx1-1
      mg_ny(mg_lmax) = ly1-1
      mg_nz(mg_lmax) = lz1-1

      if (nid.eq.0) write(*,*) 'mg_nx:',(mg_nx(i),i=1,mg_lmax)
      if (nid.eq.0) write(*,*) 'mg_ny:',(mg_ny(i),i=1,mg_lmax)
      if (nid.eq.0) write(*,*) 'mg_nz:',(mg_nz(i),i=1,mg_lmax)

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_index_0 ! initialize index sets
      include 'SIZE'
      include 'HSMG'

      n = lmgn*(lmgs+1)

      call izero( mg_rstr_wt_index      , n )
      call izero( mg_mask_index         , n )
      call izero( mg_solve_index        , n )
      call izero( mg_fast_s_index       , n )
      call izero( mg_fast_d_index       , n )
      call izero( mg_schwarz_wt_index   , n )

      return
      end
c----------------------------------------------------------------------
      subroutine outfldn (x,n,txt10,ichk) ! writes into unit=40+ifiled
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(n,n,1,lelt)
      character*10 txt10
c
      integer idum,e
      save    idum
      data    idum /3/
      if (idum.lt.0)   return
      m = 40 + ifield                 ! unit #
c
C
      mtot = n*n*nelt
      if (n.gt.7.or.nelt.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelt
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelt-ne+1
      do ie=ne1,1,-ne
         l=ie-1
         do k=1,1
            if (ie.eq.ne1) write(m,116) txt10,k,ie,xmin,xmax,ichk,time
            write(m,117)
            do j=n,1,-1
              if (n.eq.2) write(m,102) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.3) write(m,103) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.4) write(m,104) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.5) write(m,105) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.6) write(m,106) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.7) write(m,107) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.8) write(m,108) ((x(i,j,k,e+l),i=1,n),e=1,ne)
            enddo
         enddo
      enddo

C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
c     if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldn0 (x,n,txt10,ichk)
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(n,n,1,lelt)
      character*10 txt10
c
      integer idum,e
      save idum
      data idum /3/
      if (idum.lt.0) return
c
C
      mtot = n*n*nelt
      if (n.gt.7.or.nelt.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelt
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelt-ne+1
      do ie=ne1,1,-ne
         l=ie-1
         do k=1,1
            if (ie.eq.ne1) write(6,116) txt10,k,ie,xmin,xmax,ichk,time
            write(6,117)
            do j=n,1,-1
              if (n.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.8) write(6,108) ((x(i,j,k,e+l),i=1,n),e=1,ne)
            enddo
         enddo
      enddo

C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
c     if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outflda (x,n,txt10,ichk) ! writes into unit=p130+ifiled
      INCLUDE 'SIZE'                      ! or into std. output for p130<9
      INCLUDE 'TSTEP'                     ! truncated below eps=p131
      INCLUDE 'INPUT'                     ! param(130)
      real x(1)
      character*10 txt10                  ! note: n is not used
c     parameter (eps=1.e-7)
C
      p130 = param(130)
      eps  = param(131)
      if (p130.le.0)    return
      m    = 6
      if (p130.gt.9)  m = p130 + ifield

      ntot = nx1*ny1*nz1*nelfld(ifield)

      xmin = glmin(x,ntot)
      xmax = glmax(x,ntot)
      xavg = glsum(x,ntot)/ntot

      if (abs(xavg).lt.eps) xavg = 0.     ! truncation

      if (nid.eq.0) write(m,10) txt10,ichk,ntot,xavg,xmin,xmax

   10 format(3X,a10,2i8,' pts, avg,min,max = ',1p3g14.6)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldan(x,n,txt10,ichk) ! writes x(1:n) into unit=p130+ifiled
      INCLUDE 'SIZE'                      ! or into std. output for 0<p130<9
      INCLUDE 'TSTEP'                     ! truncated below eps=p131
      INCLUDE 'INPUT'
      real x(1)
      character*10 txt10
c     parameter (eps=1.e-7)
C
      p130 = param(130)
      eps  = param(131)
      if (p130.le.0)    return
      m    = 6
      if (p130.gt.9)  m = p130 + ifield

      ntot = n

      xmin = glmin(x,ntot)
      xmax = glmax(x,ntot)
      xavg = glsum(x,ntot)/ntot

      if (abs(xavg).lt.eps) xavg = 0.     ! truncation

      if (nid.eq.0) write(m,10) txt10,ichk,ntot,xavg,xmin,xmax

   10 format(3X,a10,2i8,' pts, avg,min,max = ',1p3g11.3)
c  10 format(3X,a10,2i8,' pts, avg,min,max = ',1p3g14.6)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_solve(z,rhs,mlt)  !  Solve preconditioner: Mz=rhs
      real z(1),rhs(1),mlt(1)

c     Assumes that preprocessing has been completed via h1mg_setup()


      include 'SIZE'
      include 'HSMG'       ! Same array space as HSMG
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      integer p_msk,p_b


      nel   = nelfld(ifield)

      op    =  1.                                     ! Coefficients for h1mg_ax
      om    = -1.
      sigma =  1.

      l     = mg_h1_lmax
      n     = mg_h1_n(l,mg_fld)
      is    = 1                                       ! solve index

      call h1mg_schwarz(z,rhs,sigma,l)                ! z := sigma W M     rhs
                                                      !               Schwarz
      call copy(r,rhs,n)                              ! r  := rhs

      do l = mg_h1_lmax-1,2,-1                        ! DOWNWARD Leg of V-cycle
         is = is + n
         n  = mg_h1_n(l,mg_fld)
                                                      !          T
         call h1mg_rstr(r,l,.true.)                   ! r   :=  J r
                                                      !  l         l+1
!        OVERLAPPING Schwarz exchange and solve:
         call h1mg_schwarz(e(is),r,sigma,l)           ! e := sigma W M       r
      enddo                                           !  l            Schwarz l
      is = is+n
                                                      !         T
      call h1mg_rstr(r,1,.false.)                     ! r  :=  J  r
                                                      !  l         l+1
      p_msk = p_mg_msk(l,mg_fld)
      call h1mg_mask(r,mg_imask(p_msk),nel)           !        -1
      call hsmg_coarse_solve ( e(is) , r )            ! e  := A   r
      call h1mg_mask(e(is),mg_imask(p_msk),nel)       !  1     1   1


      do l = 2,mg_h1_lmax-1                           ! UNWIND.  No smoothing.
         im = is
         is = is - n
         n  = mg_h1_n(l,mg_fld)
         call hsmg_intp (w,e(im),l-1)                 ! w   :=  J e
         i1=is-1                                      !            l-1
         do i=1,n
            e(i1+i) = e(i1+i) + w(i)                  ! e   :=  e  + w
         enddo                                        !  l       l
      enddo

      l  = mg_h1_lmax
      n  = mg_h1_n(l,mg_fld)
      im = is  ! solve index
      call hsmg_intp(w,e(im),l-1)                     ! w   :=  J e
      do i = 1,n                                      !            l-1
         z(i) = z(i) + w(i)                           ! z := z + w
      enddo

c     call dsavg(z) ! Emergency hack --- to ensure continuous z!
      call dssum(z,nx1,ny1,nz1)
      call col2 (z,mlt,n)


      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_mask(w,mask,nel)
      include 'SIZE'

      real    w   (nel)
      integer mask(nel)        ! Pointer to Dirichlet BCs
      integer e

      do e=1,nel
         im = mask(e)
         call mg_mask_e(w,mask(im)) ! Zero out Dirichlet conditions
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine mg_mask_e(w,mask) ! Zero out Dirichlet conditions
      include 'SIZE'
      real w(1)
      integer mask(0:1)

      n=mask(0)
      do i=1,n
c        write(6,*) i,mask(i),n,' MG_MASK'
         w(mask(i)) = 0.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine hsmg_tnsr1(v,nv,nu,A,At)
c
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u
c
      integer nv,nu
      real v(1),A(1),At(1)
      include 'SIZE'
      include 'INPUT'
      if (.not. if3d) then
         call hsmg_tnsr1_2d(v,nv,nu,A,At)
      else
         call hsmg_tnsr1_3d(v,nv,nu,A,At,At)
      endif
      return
      end
c-------------------------------------------------------T--------------
      subroutine hsmg_tnsr1_2d(v,nv,nu,A,Bt) ! u = A u B
      integer nv,nu
      real v(1),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work(lx1*lx1)
      integer e

      nv2 = nv*nv
      nu2 = nu*nu

      if (nv.le.nu) then
         iv=1
         iu=1
         do e=1,nelt
            call mxm(A,nv,v(iu),nu,work,nu)
            call mxm(work,nv,Bt,nu,v(iv),nv)
            iv = iv + nv2
            iu = iu + nu2
         enddo
      else
         do e=nelt,1,-1
            iu=1+nu2*(e-1)
            iv=1+nv2*(e-1)
            call mxm(A,nv,v(iu),nu,work,nu)
            call mxm(work,nv,Bt,nu,v(iv),nv)
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_tnsr1_3d(v,nv,nu,A,Bt,Ct) ! v = [C (x) B (x) A] u
      integer nv,nu
      real v(1),A(1),Bt(1),Ct(1)
      include 'SIZE'
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
      integer e,e0,ee,es

      e0=1
      es=1
      ee=nelt

      if (nv.gt.nu) then
         e0=nelt
         es=-1
         ee=1
      endif

      nu3 = nu**3
      nv3 = nv**3

      do e=e0,ee,es
         iu = 1 + (e-1)*nu3
         iv = 1 + (e-1)*nv3
         call mxm(A,nv,v(iu),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(iv),nv)
      enddo

      return
      end
c------------------------------------------   T  -----------------------
      subroutine h1mg_rstr(r,l,ifdssum) ! r =J r,   l is coarse level
      include 'SIZE'
      include 'HSMG'
      logical ifdssum

      real r(1)
      integer l

      call hsmg_do_wt(r,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))

      call hsmg_tnsr1(r,mg_nh(l),mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))

      if (ifdssum) call hsmg_dssum(r,l)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup(h1,h2)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMG'

      real h1(1),h2(1)

      integer p_h1,p_h2,p_g,p_b,p_msk

      
      param(59) = 1
      call geom_reset(1)  ! Recompute g1m1 etc. with deformed only
      ifield=1

      n = nx1*ny1*nz1*nelt

      call h1mg_setup_mg_nx
      call h1mg_setup_semhat ! SEM hat matrices for each level
      call hsmg_setup_intp   ! Interpolation operators
      call h1mg_setup_dssum  ! set direct stiffness summation handles
      call h1mg_setup_wtmask ! set restriction weight matrices and bc masks
      call h1mg_setup_schwarz_wt(.false.)
      call hsmg_setup_solve  ! set up the solver

      l=mg_h1_lmax
      call mg_set_h1  (p_h1,h1,l)
      call mg_set_h2  (p_h2,h2,l)
      call mg_set_gb  (p_g,p_b,l)
      call mg_set_msk (p_msk,l)

      call h1mg_setup_fdm    ! set up fast diagonalization method
      call h1mg_setup_schwarz_wt(.false.)
      call hsmg_setup_solve  ! set up the solver

      l=mg_h1_lmax
      call mg_set_h1  (p_h1,h1,l)
      call mg_set_h2  (p_h2,h2,l)
      call mg_set_gb  (p_g,p_b,l)
      call mg_set_msk (p_msk,l)

      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_setup_mg_nx()
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'SEMHAT'
      include 'TSTEP'   ! nelfld
      real w(lx1+2)
      integer nf,nc,nr
      integer nx,ny,nz

      integer mgn2(10)
      save    mgn2
      data    mgn2 / 1, 2, 2, 2, 2, 3, 3, 5, 5, 5/
c     data    mgn2 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 0

      p82 = 2                          ! potentially variable nxc

c     mg_h1_lmax = 4
c     mg_h1_lmax = 3
c     mg_h1_lmax = 2
      mg_h1_lmax = 3

c     uuuu

      if (lx1.eq.4) mg_h1_lmax = 2


      mgnx1    = p82-1 !1
      mg_nx(1) = mgnx1
      mg_ny(1) = mgnx1
      mg_nz(1) = mgnx1
      if (.not.if3d) mg_nz(1) = 0

      mgnx2 = 2*(lx2/4) + 1
      if (lx1.eq.5)  mgnx2 = 3
c     if (lx1.eq.6)  mgnx2 = 3
      if (lx1.le.10) mgnx2 = mgn2(nx1)
      if (lx1.eq.8)  mgnx2 = 4
      if (lx1.eq.8)  mgnx2 = 3

      mgnx2 = min(3,mgnx2)  ! This choice seems best (9/24/12)

      mg_nx(2) = mgnx2
      mg_ny(2) = mgnx2
      mg_nz(2) = mgnx2
      if (.not.if3d) mg_nz(2) = 0

      mg_nx(3) = mgnx2+1
      mg_ny(3) = mgnx2+1
      mg_nz(3) = mgnx2+1
      if (.not.if3d) mg_nz(3) = 0

      mg_nx(mg_h1_lmax) = lx1-1
      mg_ny(mg_h1_lmax) = ly1-1
      mg_nz(mg_h1_lmax) = lz1-1

      if (nid.eq.0) write(*,*) 'h1_mg_nx:',(mg_nx(i),i=1,mg_h1_lmax)
      if (nid.eq.0) write(*,*) 'h1_mg_ny:',(mg_ny(i),i=1,mg_h1_lmax)
      if (nid.eq.0) write(*,*) 'h1_mg_nz:',(mg_nz(i),i=1,mg_h1_lmax)

      if (nid.eq.0) write(6,*) 'ldimt1',ldimt1
      do ifld=1,ldimt1
      do l=1,mg_lmax
         mg_h1_n(l,ifld)=(mg_nx(l)+1)
     $                  *(mg_ny(l)+1)
     $                  *(mg_nz(l)+1)*nelfld(ifld)
      enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_semhat ! SEM hat matrices for each level
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      include 'SEMHAT'

      do l=1,mg_h1_lmax
         n = mg_nx(l)     ! polynomial order
         call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
         call copy(mg_ah(1,l),ah,(n+1)*(n+1))
         call copy(mg_bh(1,l),bh,n+1)
         call copy(mg_dh(1,l),dh,(n+1)*(n+1))
         call transpose(mg_dht(1,l),n+1,dh,n+1)
         call copy(mg_zh(1,l),zh,n+1)

         mg_nh(l)=n+1
         mg_nhz(l)=mg_nz(l)+1

      enddo
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_dssum
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'ZPER'
      parameter (lxyz=(lx1+2)*(ly1+2)*(lz1+2))
      common /c_is1/ glo_num(lxyz*lelt)
      common /ivrtx/ vertex ((2**ldim)*lelt)

      integer*8 glo_num
      integer vertex
      integer nx,ny,nz
      integer l

      ncrnr = 2**ndim
      if (nid.eq.0) write(6,*) 'hsmg_setup_dssum, ifgtp=',ifgtp
      if (ifgtp) then
          if (nid.eq.0) write(6,*) 'use regular mesh for multigrid'
          call gen_gtp_vertex    (vertex, ncrnr)
      else
          call get_vert
      endif


      do l=1,mg_lmax  ! set up direct stiffness summation for each level
         nx=mg_nh(l)
         ny=mg_nh(l)
         nz=mg_nhz(l)
         call setup_ds(mg_gsh_handle(l,mg_fld),nx,ny,nz
     $                ,nelt,nelgv,vertex,glo_num)
         nx=nx+2
         ny=ny+2
         nz=nz+2
         if(.not.if3d) nz=1
         call setup_ds(mg_gsh_schwarz_handle(l,mg_fld),nx,ny,nz
     $                ,nelt,nelgv,vertex,glo_num)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine mg_set_msk(p_msk ,l0)
      include 'SIZE'
      include 'HSMG'
      include 'TSTEP'
      integer p_msk

      l                  = mg_h1_lmax
      p_mg_msk(l,mg_fld) = 0
      n                  = mg_h1_n(l,mg_fld)


      do l=mg_h1_lmax,1,-1
         nx = mg_nh  (l)
         ny = mg_nh  (l)
         nz = mg_nhz (l)

         p_msk = p_mg_msk(l,mg_fld)

         call h1mg_setup_mask
     $     (mg_imask(p_msk),nm,nx,ny,nz,nelfld(ifield),l,mg_work)

         if (l.gt.1) p_mg_msk(l-1,mg_fld)=p_mg_msk(l,mg_fld)+nm

      enddo

      p_msk = p_mg_msk(l0,mg_fld)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_mask(mask,nm,nx,ny,nz,nel,l,w)
      include 'SIZE'
      include 'INPUT'        ! if3d

      integer mask(1)        ! Pointer to Dirichlet BCs
      integer nx,ny,nz,l
      real w(nx,ny,nz,nel)

      integer e,count,ptr
      integer lbr,rbr,lbs,rbs,lbt,rbt,two

      zero = 0
      nxyz = nx*ny*nz
      n    = nx*ny*nz*nel

      call rone(w,n)   ! Init everything to 1

      ierrmx = 0       ! BC verification
      two    = 2
      do e=1,nel       ! Set dirichlet nodes to zero

         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,two,ierr)
c        write(6,6) e,lbr,rbr,lbs,rbs,ierr,nx
c   6    format(i5,2x,4i3,2x,i2,3x,i5,'  lbr,rbr,lbs')

         if (lbr.eq.1) call facev(w,e,4,zero,nx,ny,nz)
         if (rbr.eq.1) call facev(w,e,2,zero,nx,ny,nz)
         if (lbs.eq.1) call facev(w,e,1,zero,nx,ny,nz)
         if (rbs.eq.1) call facev(w,e,3,zero,nx,ny,nz)
         if (if3d) then
            if (lbt.eq.1) call facev(w,e,5,zero,nx,ny,nz)
            if (rbt.eq.1) call facev(w,e,6,zero,nx,ny,nz)
         endif
         ierrmx = max(ierrmx,ierr)
      enddo

      call hsmg_dsprod(w,l)    ! direct stiffness multiply

c
c     Prototypical mask layout, nel=5:
c
c    e=1 ...             10
c      1  2  3  4  5 ... 10 | 11 12 13 14 | 15 | 16 |
c     11 15 16 ...          |  3 p1 p2 p3 |  0 |  0 | ...
c                              ^
c                              |
c                              |_count for e=1
c

      nm  = 1                  ! store mask
      do e=1,nel

         mask(e) = nel+nm
         count   = 0          ! # Dirchlet points on element e
         ptr     = mask(e)

         do i=1,nxyz
            if (w(i,1,1,e).eq.0) then
               nm    = nm   +1
               count = count+1
               ptr   = ptr  +1
               mask(ptr) = i + nxyz*(e-1)   ! where I mask on element e
            endif
         enddo


         ptr       = mask(e)
         mask(ptr) = count

         nm        = nm+1     ! bump pointer to hold next count

      enddo

      nm = nel + nm-1 ! Return total number of mask pointers/counters

      ierrmx = iglmax(ierrmx,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL h1'
         call exitti('D INVALID BC FOUND in h1mg_setup_mask$',ierrmx,1)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine mg_set_h1  (p_h1,h1,l0)
      include 'SIZE'
      include 'HSMG'
      integer pf,pc

      real h1 (lx1,ly1,lz1,lelv)

      integer p_h1

      l                 = mg_h1_lmax
      p_mg_h1(l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)

      call copy (mg_h1,h1,n)   ! Fine grid is just original h1

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)

      do l=mg_h1_lmax-1,1,-1

         p_mg_h1(l,mg_fld) = p_mg_h1(l+1,mg_fld) + n
         n                 = mg_h1_n(l  ,mg_fld)

         pf                = p_mg_h1(l+1,mg_fld)
         pc                = p_mg_h1(l  ,mg_fld)

         call hsmg_intp_fc (mg_h1(pc),mg_h1(pf),l)

      enddo

      p_h1 = p_mg_h1(l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_set_h2  (p_h2,h2,l0)
      include 'SIZE'
      include 'HSMG'

      real h2 (lx1,ly1,lz1,lelv)

      integer p_h2,pf,pc

      l                 = mg_h1_lmax
      p_mg_h2(l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)

      call copy (mg_h2,h2,n)   ! Fine grid is just original h2

      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nhz(l)

      do l=mg_h1_lmax-1,1,-1

         p_mg_h2(l,mg_fld) = p_mg_h2(l+1,mg_fld) + n
         n                 = mg_h1_n(l  ,mg_fld)

         pf                = p_mg_h2(l+1,mg_fld)
         pc                = p_mg_h2(l  ,mg_fld)

         call hsmg_intp_fc (mg_h2(pc),mg_h2(pf),l)

      enddo

      p_h2 = p_mg_h2(l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine hsmg_intp_fc(uc,uf,l) ! l is coarse level

      include 'SIZE'
      include 'HSMG'

      real uc(1),uf(1)


      nc = mg_nh(l)
      nf = mg_nh(l+1)
      call hsmg_tnsr(uc,nc,uf,nf,mg_jhfc(1,l),mg_jhfct(1,l))

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_intp_fc_e(uc,uf,nxc,nyc,nzc,nxf,nyf,nzf,e,l,w)
      include 'SIZE'
      include 'INPUT'      ! if3d
      include 'HSMG'
      integer  e

      real uf(nxf,nyf,nzf),uc(nxc,nyc,nzc),w(1)

      if (if3d) then

         n1=nxf*nyf
         n2=nzf
         n3=nzc
         call mxm(uf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()
         lc=1 + n1*n3
         lc0=lc

         n1=nxf
         n2=nyf
         n3=nyc

         do k=1,nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,w(lc),n3)
            lf = lf + n1*n2
            lc = lc + n1*n3
         enddo

         lf=lc0  ! Rewind fine pointer to start of coarse data
         n1=nxc
         n2=nxf
         n3=nyc*nzc
         call mxm(mg_jhfc(1,l),n1,w(lf),n2,uc,n3)

      else ! 2D

         n1=nxf
         n2=nyf
         n3=nyc
         call mxm(uf,n1,mg_jhfct(1,l),n2,w,n3)

         n1=nxc
         n2=nxf
         n3=nyc
         call mxm(mg_jhfc(1,l),n1,w,n2,uc,n3)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_intp_gfc_e(gc,gf,ng,nxc,nyc,nzc,nxf,nyf,nzf,e,l,w)
      include 'SIZE'
      include 'INPUT'      ! if3d
      include 'HSMG'
      integer e
      real gf(ng,nxf,nyf,nzf),gc(ng,nxc,nyc,nzc),w(1)


      if (if3d) then

         n1=ng*nxf*nyf
         n2=nzf
         n3=nzc
         call mxm(gf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()
         lc=1 + n1*n3
         lc0=lc

         n1=ng*nxf
         n2=nyf
         n3=nyc

         do k=1,nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,w(lc),n3)
            lf = lf + n1*n2
            lc = lc + n1*n3
         enddo

         lf=lc0  ! Rewind fine pointer to start of coarse data
         n1=ng
         n2=nxf
         n3=nxc

         do k=1,nyc*nzc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,gc(1,1,k,1),n3)
            lf = lf + n1*n2
         enddo

      else ! 2D

         n1=ng*nxf
         n2=nyf
         n3=nyc
         call mxm(gf,n1,mg_jhfct(1,l),n2,w,n3)

         lf=1           ! pointers into work array w()

         n1=ng
         n2=nxf
         n3=nxc

         do k=1,nyc
            call mxm(w(lf),n1,mg_jhfct(1,l),n2,gc(1,1,k,1),n3)
            lf = lf + n1*n2
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_scale_mass (b,g,wt,ng,nx,ny,nz,wk,ifinv)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'HSMG'

      real b(1),g(ng,1),wt(1),wk(1)
      logical ifinv

      common /ctmp0/ wi(2*lx1+4)

      n = nx*ny*nz

      if (nx.le.2*lx1) then
         if (ifinv) then
            call invers2(wi,wt,nx)
         else
            call copy(wi,wt,nx)
         endif
      else
         call exitti('mg_scale_mass: wi too small$',nx,1)
      endif

      if (if3d) then
         l=0
         do k=1,nz
         do j=1,ny
            wjk=wi(j)*wi(k)
            do i=1,nx
               l=l+1
               wk(l) = wjk*wi(i)
            enddo
         enddo
         enddo

         do k=1,n
            b(k)   = wk(k)*b(k)
            g(1,k) = wk(k)*g(1,k)
            g(2,k) = wk(k)*g(2,k)
            g(3,k) = wk(k)*g(3,k)
            g(4,k) = wk(k)*g(4,k)
            g(5,k) = wk(k)*g(5,k)
            g(6,k) = wk(k)*g(6,k)
         enddo

      else      ! 2D
         l=0
         do j=1,ny
         do i=1,nx
            l=l+1
            wk(l) = wi(i)*wi(j)
         enddo
         enddo

         do k=1,n
            b(k)   = wk(k)*b(k)
            g(1,k) = wk(k)*g(1,k)
            g(2,k) = wk(k)*g(2,k)
            g(3,k) = wk(k)*g(3,k)
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mg_set_gb  (p_g,p_b,l0)
      include 'SIZE'
      include 'HSMG'
      include 'MASS'   ! bm1
      include 'TSTEP'  ! nelfld

      integer p_g,p_b,e
      common /ctmp1/ w(lx1*ly1*lz1*lelt*2)

      l                 = mg_h1_lmax
      p_mg_b (l,mg_fld) = 0
      p_mg_g (l,mg_fld) = 0
      n                 = mg_h1_n(l,mg_fld)


      ng = 3*(ndim-1)  ! 3 or 6 elements to symm dxd tensor

      do l=mg_h1_lmax-1,1,-1

         p_mg_b (l,mg_fld) = p_mg_b (l+1,mg_fld) + n
         p_mg_g (l,mg_fld) = p_mg_g (l+1,mg_fld) + n*ng
         n                 = mg_h1_n(l  ,mg_fld)

      enddo

      do e=1,nelfld(ifield) !FIXME misun 10/16/2014 ifield=1 ??
       do l=mg_h1_lmax,1,-1

         nx = mg_nh(l)
         ny = mg_nh(l)
         nz = mg_nhz(l)
         nxyz = nx*ny*nz

         p_g = p_mg_g (l,mg_fld) + ng*nx*ny*nz*(e-1)
         p_b = p_mg_b (l,mg_fld) +    nx*ny*nz*(e-1)

         if (l.eq.mg_h1_lmax) then
            call gxfer_e (mg_g(p_g) ,ng,e             ) ! Fine grid=original G
            call copy    (mg_b(p_b) ,bm1(1,1,1,e),nxyz) ! Fine grid=original B
            call mg_scale_mass                          ! Divide out Wghts
     $         (mg_b(p_b),mg_g(p_g),mg_bh(1,l),ng,nx,ny,nz,w,.true.)
         else

c           Generate G and B by interpolating their continous counterparts onto
c           the coarse grid and collocating with coarse-grid quadrature weights

            call mg_intp_gfc_e
     $            (mg_g(p_g),mg_g(l_g),ng,nx,ny,nz,nxl,nyl,nzl,e,l,w)

            call mg_intp_fc_e
     $            (mg_b(p_b),mg_b(l_b)   ,nx,ny,nz,nxl,nyl,nzl,e,l,w)

            call mg_scale_mass                         ! Reinstate weights
     $      (mg_b(l_b),mg_g(l_g),mg_bh(1,l+1),ng,nxl,nyl,nzl,w,.false.)

         endif

         l_b = p_b
         l_g = p_g

         nxl = nx
         nyl = ny
         nzl = nz

       enddo

       call mg_scale_mass                         ! Reinstate weights
     $      (mg_b(l_b),mg_g(l_g),mg_bh(1,1),ng,nxl,nyl,nzl,w,.false.)


      enddo

      p_b  = p_mg_b (l0,mg_fld)
      p_g  = p_mg_g (l0,mg_fld)

      return
      end
c-----------------------------------------------------------------------
      subroutine gxfer_e (g,ng,e)
      include 'SIZE'
      include 'TOTAL'

      real g(ng,1)
      integer e

      nxyz = nx1*ny1*nz1

c     ifdfrm(e) = .true.  ! TOO LATE

      if (if3d) then
         do i=1,nxyz
            g(1,i) = g1m1(i,1,1,e)
            g(2,i) = g2m1(i,1,1,e)
            g(3,i) = g3m1(i,1,1,e)
            g(4,i) = g4m1(i,1,1,e)
            g(5,i) = g5m1(i,1,1,e)
            g(6,i) = g6m1(i,1,1,e)
         enddo
      else
         do i=1,nxyz
            g(1,i) = g1m1(i,1,1,e)
            g(2,i) = g2m1(i,1,1,e)
            g(3,i) = g4m1(i,1,1,e)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine chkr(name3,ii)
      include 'SIZE'
      include 'TOTAL'
      include 'HSMG'
      character*3 name3

      write(6,*) mg_h1_lmax,ii,' ',name3,' CHKR'

      return
      end
c-----------------------------------------------------------------------
      subroutine outgmat(a,ng,nx,ny,name6,k,e)

      integer e
      real a(ng,nx,ny)
      common /ctmp0/ w(100000)
      character*6 name6

c     do i=1,ng
      do i=1,1
         sum = 0.
         do ii=1,nx*ny
            w(ii)=a(i,ii,1)
            sum = sum + a(i,ii,1)
         enddo

         write(6,1) name6,i,k,e,nx,ny,ng,sum
    1    format(a6,6i5,f12.5,'  outgmat')

         call outmatz(w,nx,ny,name6,i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine outmatz(a,m,n,name6,ie)
      real a(m,n)
      character*6 name6

      sum=0.
      sua=0.
      do i=1,m*n
         sum=sum+    a(i,1)
         sua=sua+abs(a(i,1))
      enddo
      sum=sum/(m*n)
      sua=sua/(m*n)

      write(6,*)
      write(6,1) ie,name6,m,n,sum,sua
    1 format(i8,' matrix: ',a6,2i5,1p2e12.4)

      n12 = min(m,12)
      do j=m,1,-1
         write(6,6) ie,name6,(a(i,j),i=1,n12)
      enddo
    6 format(i3,1x,a6,12f9.5)
c     write(6,*)
      return
      end
c-----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      real wt(1),work(1)
      logical ifsqrt

      if (ndim.eq.2) call h1mg_setup_schwarz_wt2d_2(wt,ie,n,work,ifsqrt)
      if (ndim.eq.3) call h1mg_setup_schwarz_wt3d_2(wt,ie,n,work,ifsqrt)

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt2d_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      logical ifsqrt
      integer n
      real wt(n,4,2,nelt)
      real work(n,n)

      integer ie,i,j
      do j=1,n
         wt(j,1,1,ie)=1.0/work(1,j)
         wt(j,2,1,ie)=1.0/work(2,j)
         wt(j,3,1,ie)=1.0/work(n-1,j)
         wt(j,4,1,ie)=1.0/work(n,j)
      enddo
      do i=1,n
         wt(i,1,2,ie)=1.0/work(i,1)
         wt(i,2,2,ie)=1.0/work(i,2)
         wt(i,3,2,ie)=1.0/work(i,n-1)
         wt(i,4,2,ie)=1.0/work(i,n)
      enddo
      if(ifsqrt) then
         do ii=1,2
         do j=1,4
         do i=1,n
            wt(i,j,ii,ie)=sqrt(wt(i,j,ii,ie))
         enddo
         enddo
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt3d_2(wt,ie,n,work,ifsqrt)
      include 'SIZE'
      logical ifsqrt
      integer n
      real wt(n,n,4,3,nelt)
      real work(n,n,n)

      integer ie,i,j,k
      integer lbr,rbr,lbs,rbs,lbt,rbt

      ierr = 0
      do k=1,n
      do j=1,n
         wt(j,k,1,1,ie)=1.0/work(1,j,k)
         wt(j,k,2,1,ie)=1.0/work(2,j,k)
         wt(j,k,3,1,ie)=1.0/work(n-1,j,k)
         wt(j,k,4,1,ie)=1.0/work(n,j,k)
      enddo
      enddo
      do k=1,n
      do i=1,n
         wt(i,k,1,2,ie)=1.0/work(i,1,k)
         wt(i,k,2,2,ie)=1.0/work(i,2,k)
         wt(i,k,3,2,ie)=1.0/work(i,n-1,k)
         wt(i,k,4,2,ie)=1.0/work(i,n,k)
      enddo
      enddo
      do j=1,n
      do i=1,n
         wt(i,j,1,3,ie)=1.0/work(i,j,1)
         wt(i,j,2,3,ie)=1.0/work(i,j,2)
         wt(i,j,3,3,ie)=1.0/work(i,j,n-1)
         wt(i,j,4,3,ie)=1.0/work(i,j,n)
      enddo
      enddo
      if(ifsqrt) then
         do ii=1,3
         do k=1,4
         do j=1,4
         do i=1,n
            wt(i,j,k,ii,ie)=sqrt(wt(i,j,k,ii,ie))
         enddo
         enddo
         enddo
         enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine h1mg_setup_schwarz_wt_1(wt,l,ifsqrt)
      include 'SIZE'
      include 'INPUT'  ! if3d
      include 'TSTEP'  ! ifield
      include 'HSMG'

      real wt(1),work(1)
      logical ifsqrt

      integer enx,eny,enz,pm

      zero =  0
      one  =  1
      onem = -1

      n  = mg_h1_n (l,mg_fld)
      pm = p_mg_msk(l,mg_fld)

      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      if(.not.if3d) enz=1
      ns = enx*eny*enz*nelfld(ifield)
      i  = ns+1

      call rone(mg_work(i),ns)

c     Sum overlap region (border excluded)
      call hsmg_extrude(mg_work,0,zero,mg_work(i),0,one ,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work(i),l)
      call hsmg_extrude(mg_work(i),0,one ,mg_work,0,onem,enx,eny,enz)
      call hsmg_extrude(mg_work(i),2,one,mg_work(i),0,one,enx,eny,enz)

      if(.not.if3d) then ! Go back to regular size array
         call hsmg_schwarz_toreg2d(mg_work,mg_work(i),mg_nh(l))
      else
         call hsmg_schwarz_toreg3d(mg_work,mg_work(i),mg_nh(l))
      endif

      call hsmg_dssum(mg_work,l)                           ! sum border nodes


      nx = mg_nh(l)
      ny = mg_nh(l)
      nz = mg_nh(l)
      if (.not.if3d) nz=1
      nxyz = nx*ny*nz
      k    = 1
      do ie=1,nelfld(ifield)
c        call outmat(mg_work(k),nx,ny,'NEW WT',ie)
         call h1mg_setup_schwarz_wt_2(wt,ie,nx,mg_work(k),ifsqrt)
         k = k+nxyz
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine gen_fast(df,sr,ss,st,x,y,z)
c
c     Generate fast diagonalization matrices for each element
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'WZ'
c
      parameter(lxx=lx1*lx1)
      real df(lx1*ly1*lz1,1),sr(lxx*2,1),ss(lxx*2,1),st(lxx*2,1)
c
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
c
      integer lbr,rbr,lbs,rbs,lbt,rbt,e
c
      real x(nx1,ny1,nz1,nelt)
      real y(nx1,ny1,nz1,nelt)
      real z(nx1,ny1,nz1,nelt)
      real axwt(lx2)

      ierr = 0

      do e=1,nelt
c
         if (param(44).eq.1) then ! FIXME misun 9/25/2014
           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,2,ierr)
         else
           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,3,ierr)
         endif
c
c        Set up matrices for each element.
c
         if (param(44).eq.1) then
           call set_up_fast_1D_fem( sr(1,e),lr,nr ,lbr,rbr
     $                      ,llr(e),lmr(e),lrr(e),zgm2(1,1),nx2,e)
         else
           call set_up_fast_1D_sem( sr(1,e),lr,nr ,lbr,rbr
     $                      ,llr(e),lmr(e),lrr(e),e)
         endif
         if (ifaxis) then
            xsum = vlsum(wxm2,nx2)
            do i=1,ny2
               yavg = vlsc2(y(1,i,1,e),wxm2,nx2)/xsum
               axwt(i) = yavg
            enddo
            call set_up_fast_1D_fem_ax( ss(1,e),ls,ns ,lbs,rbs
     $                 ,lls(e),lms(e),lrs(e),zgm2(1,2),axwt,ny2,e)
         else
            if (param(44).eq.1) then
               call set_up_fast_1D_fem( ss(1,e),ls,ns ,lbs,rbs
     $                      ,lls(e),lms(e),lrs(e),zgm2(1,2),ny2,e)
            else
               call set_up_fast_1D_sem( ss(1,e),ls,ns ,lbs,rbs
     $                      ,lls(e),lms(e),lrs(e),e)
            endif
         endif
         if (if3d) then
            if (param(44).eq.1) then
               call set_up_fast_1D_fem( st(1,e),lt,nt ,lbt,rbt
     $                      ,llt(e),lmt(e),lrt(e),zgm2(1,3),nz2,e)
            else
               call set_up_fast_1D_sem( st(1,e),lt,nt ,lbt,rbt
     $                      ,llt(e),lmt(e),lrt(e),e)
            endif
         endif
c
c        DIAGNOSTICS
c
c        n12 = min(9,nr)
c        write(6,1) e,'1D lr',llr(e),lmr(e),lrr(e),(lr(k),k=1,n12)
c        write(6,1) e,'1D ls',lls(e),lms(e),lrs(e),(ls(k),k=1,n12)
c        if (if3d)
c    $   write(6,1) e,'1D lt',llt(e),lmt(e),lrt(e),(lt(k),k=1,n12)
c   1    format(i6,1x,a5,1p12e12.4)
c
c
c        Set up diagonal inverse
c
         if (if3d) then
            eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $                  +  vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            l   = 1
            do k=1,nt
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j) + lt(k)
               if (diag.gt.eps) then
                  df(l,e) = 1.0/diag
               else
c                 write(6,3) e,'Reset Eig in gen fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
c   3             format(i6,1x,a21,4i5,1p5e12.4)
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
            enddo
         else
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            l   = 1
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j)
               if (diag.gt.eps) then
                  df(l,e) = 1.0/diag
               else
c                 write(6,2) e,'Reset Eig in gen fast:',i,j,l
c    $                         ,eps,diag,lr(i),ls(j)
c   2             format(i6,1x,a21,3i5,1p4e12.4)
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
         endif
c
c        Next element ....
c
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('E INVALID BC FOUND in genfast$',ierrmx,1)
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine gen_fast_spacing(x,y,z)
c
c     Generate fast diagonalization matrices for each element
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'WZ'
c
      parameter(lxx=lx1*lx1)
c
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
c
      integer lbr,rbr,lbs,rbs,lbt,rbt,e
c
      real x(nx1,ny1,nz1,nelt)
      real y(nx1,ny1,nz1,nelt)
      real z(nx1,ny1,nz1,nelt)
      real axwt(lx2)

      ierr = 0

      if (param(44).eq.1) then
c                                    __ __ __
c        Now, for each element, compute lr,ls,lt between specified planes
c
         n1 = nx2
         n2 = nx2+1
         nz0 = 1
         nzn = 1
         if (if3d) then
            nz0= 0
            nzn=n2
         endif
         eps = 1.e-7
         if (wdsize.eq.8)  eps = 1.e-14
c
c        Find mean spacing between "left-most" planes
         call plane_space2(llr,lls,llt, 0,wxm2,x,y,z,n1,n2,nz0,nzn)
c
c        Find mean spacing between "middle" planes
         call plane_space (lmr,lms,lmt, 1,n1,wxm2,x,y,z,n1,n2,nz0,nzn)
c
c        Find mean spacing between "right-most" planes
         call plane_space2(lrr,lrs,lrt,n2,wxm2,x,y,z,n1,n2,nz0,nzn)
c
      else
         call load_semhat_weighted    !   Fills the SEMHAT arrays
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine plane_space_std(lr,ls,lt,i1,i2,w,x,y,z,nx,nxn,nz0,nzn)
c
c     This routine now replaced by "plane_space()"
c
c     Here, spacing is based on arithmetic mean.
c     New verision uses harmonic mean.  pff 2/10/07
c
      include 'SIZE'
      include 'INPUT'
c
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
      j2 = i2
      k2 = i2
c
      do ie=1,nelt
c
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
               lr2  = lr2  + ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
     $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
     $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
               ls2  = ls2  + ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
     $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
     $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
               lt2  = lt2  + ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
     $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
     $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = sqrt(lt2)
c
         else
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
               lr2  = lr2  + ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
     $                     +   (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
               ls2  = ls2  + ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
     $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine plane_space(lr,ls,lt,i1,i2,w,x,y,z,nx,nxn,nz0,nzn)
c
c     Here, spacing is based on harmonic mean.  pff 2/10/07
c
c
      include 'SIZE'
      include 'INPUT'
c
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
      j2 = i2
      k2 = i2
c
      do ie=1,nelt
c
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
c              lr2  = lr2  + ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
c    $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
c    $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
c    $                     *   weight
               lr2  = lr2  +   weight /
     $                       ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
     $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
     $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = 1./sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
c              ls2  = ls2  + ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
c    $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
c    $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
c    $                     *   weight
               ls2  = ls2  +   weight /
     $                       ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
     $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
     $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = 1./sqrt(ls2)
c
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
c              lt2  = lt2  + ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
c    $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
c    $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
c    $                     *   weight
               lt2  = lt2  +   weight /
     $                       ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
     $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
     $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = 1./sqrt(lt2)
c
         else              ! 2D
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
c              lr2  = lr2  + ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
c    $                     +   (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
c    $                     *   weight
               lr2  = lr2  + weight /
     $                       ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
     $                       + (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = 1./sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
c              ls2  = ls2  + ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
c    $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
c    $                     *   weight
               ls2  = ls2  + weight /
     $                       ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
     $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = 1./sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine plane_space2(lr,ls,lt,i1,w,x,y,z,nx,nxn,nz0,nzn)
c
c     Here, the local spacing is already given in the surface term.
c     This addition made to simplify the periodic bdry treatment.
c
c
      include 'SIZE'
      include 'INPUT'
c
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
c
      do ie=1,nelt
c
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
               lr2  = lr2  + ( (x(i1,j,k,ie))**2
     $                     +   (y(i1,j,k,ie))**2
     $                     +   (z(i1,j,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
               ls2  = ls2  + ( (x(i,j1,k,ie))**2
     $                     +   (y(i,j1,k,ie))**2
     $                     +   (z(i,j1,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
               lt2  = lt2  + ( (x(i,j,k1,ie))**2
     $                     +   (y(i,j,k1,ie))**2
     $                     +   (z(i,j,k1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = sqrt(lt2)
c           write(6,1) 'lrlslt',ie,lr(ie),ls(ie),lt(ie)
    1       format(a6,i5,1p3e12.4)
c
         else
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
               lr2  = lr2  + ( (x(i1,j,1,ie))**2
     $                     +   (y(i1,j,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
               ls2  = ls2  + ( (x(i,j1,1,ie))**2
     $                     +   (y(i,j1,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie),lt(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_fem(s,lam,n,lbc,rbc,ll,lm,lr,z,nz,e)
      real s(1),lam(1),ll,lm,lr,z(1)
      integer lbc,rbc,e
c
      parameter (m=100)
      real dx(0:m)
      integer icalld
      save    icalld
      data    icalld/0/
c
      icalld=icalld+1
c
      if (nz.gt.m-3) then
         write(6,*) 'ABORT. Error in set_up_fast_1D_fem. Increase m to'
     $             , nz
         call exitt(1)
      endif
c
c     In the present scheme, each element is viewed as a d-fold
c     tensor of (1+nz+1) arrays, even if funky bc's are applied
c     on either end of the 1D array.
c
      n = nz+2
c
c     Compute spacing, dx()
c
      call set_up_1D_geom(dx,lbc,rbc,ll,lm,lr,z,nz)
c
      nn1 = n*n + 1
      call gen_eigs_A_fem(s,s(nn1),lam,n,dx,lbc,rbc)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_1D_geom(dx,lbc,rbc,ll,lm,lr,z,nz)
c
c     real dx(0:1),ll,lm,lr,z(2)
      real dx(0:1),ll,lm,lr,z(1) !FIXME misun 10/16/2014
      integer lbc,rbc
c
c
c     Set up the 1D geometry for the tensor-product based overlapping Schwarz
c
c     Upon return:
c
c       dx() contains the spacing required to set up the stiffness matrix.
c
c
c     Input:
c
c       lbc (rbc) is 0 if the left (right) BC is Dirichlet, 1 if Neumann.
c
c       ll is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the LEFT element
c
c       lm is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the MIDDLE element
c
c       lr is the space between the right-most Gauss point of the middle
c          element and the left-most Gauss point of the RIGHT element
c
c       --- if ll (lr) is very small (0), it indicates that there is no
c           left (right) spacing, and that the left (right) BC is Neumann.
c
c
c       z() is the array of nz Gauss points on the interval ]-1,1[.
c
c     Boundary conditions:
c
c     bc = 0  -- std. Dirichlet bc applied 2 points away from interior
c     bc = 1  -- Dirichlet bc applied 1 point away from interior (outflow)
c     bc = 2  -- Neumann bc applied on interior point (W,v,V,SYM,...)
c
c
c
c     Geometry:
c
c
c        dx0       dx1   dx2    dx3   dx5    dx5        dx6
c
c    bl        |<--ll-->|<------lm------>|<---lr--->|           br
c     0--------x-----|--x---x--------x---x--|-------x-----------0
c
c       left elem.         middle elem.         right elem.
c                   -1                     +1
c
c
c    "bl" = (extrapolated) location of Gauss point NX2-1 in left elem.
c
c    "br" = (extrapolated) location of Gauss point 2 in right elem.
c
c    Overlapping Schwarz applied with homogeneous Dirichlet boundary
c    conditions at "bl" and "br", and with a single d.o.f. extending
c    in to each adjacent domain.
c
      eps = 1.e-5
      call rone(dx,nz+3)
c
c     Middle
      scale = lm/(z(nz)-z(1))
      do i=1,nz-1
         dx(i+1) = (z(i+1)-z(i))*scale
      enddo

c     Left end
      if (lbc.eq.0) then
         dzm0   = z(1) + 1.
         dxm0   = scale*dzm0
         dxln   = ll - dxm0
         scalel = dxln/dzm0
         dx(0)  = scalel*(z(2)-z(1))
         dx(1)  = ll
      elseif (lbc.eq.1) then
         dx(1)  = ll
      endif
c
c     Right end
      if (rbc.eq.0) then
         dzm0      = z(1) + 1.
         dxm0      = scale*dzm0
         dxr0      = lr - dxm0
         scaler    = dxr0/dzm0
         dx(nz+1)  = lr
         dx(nz+2)  = scaler*(z(2)-z(1))
      elseif (rbc.eq.1) then
         dx(nz+1)  = lr
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_eigs_A_fem(sf,sft,atd,n,l,lbc,rbc)
c
c     Set up Fast diagonalization solver for FEM on mesh 2
c
      real sf(n,n),sft(n,n),atd(1),l(0:1)
      integer lbc,rbc
c
      parameter (m=100)
      real atu(m),ad(m),au(m),c(m),bh(m),li(0:m)
c
      if (n.gt.m) then
         write(6,*) 'ABORT. Error in gen_eigs_A_fem. Increase m to',n
         call exitt(1)
      endif
c
c     Get delta x's
c
      do i=0,n
         li(i) = 1.0/l(i)
      enddo
c                          ^   ^
c     Fill initial arrays, A & B:
c
      call rzero(ad,n)
      call rzero(au,n-1)
      call rzero(bh,n)
c
      ie1 = lbc
      ien = n-rbc
      do ie=ie1,ien
c
c        il,ir are the left and right endpts of element ie.
         il = ie
         ir = ie+1
c
         if (ie.gt.0) ad(il) = ad(il) +       li(ie)
         if (ie.lt.n) ad(ir) = ad(ir) +       li(ie)
         if (ie.gt.0) au(il) =        -       li(ie)
         if (ie.gt.0) bh(il) = bh(il) + 0.5 * l(ie)
         if (ie.lt.n) bh(ir) = bh(ir) + 0.5 * l(ie)
      enddo
c
c     Take care of bc's (using blasting)
      bhm = vlmax(bh(2),n-2)/(n-2)
      ahm = vlmax(ad(2),n-2)/(n-2)
c
      if (lbc.gt.0) then
         au(1) = 0.
         ad(1) = ahm
         bh(1) = bhm
      endif
c
      if (rbc.gt.0) then
         au(n-1) = 0.
         ad(n  ) = ahm
         bh(n  ) = bhm
      endif
c
c
      do i=1,n
         c(i) = sqrt(1.0/bh(i))
      enddo
c                                        ~
c     Scale rows and columns of A by C:  A = CAC
c
      do i=1,n
         atd(i) = c(i)*ad(i)*c(i)
      enddo
c
c     Scale upper diagonal
c
      atu(1) = 0.
      do i=1,n-1
         atu(i) = c(i)*au(i)*c(i+1)
      enddo
c                                             ~
c     Compute eigenvalues and eigenvectors of A
c
      call calcz(atd,atu,n,dmax,dmin,sf,ierr)
      if (ierr.eq.1) then
         nid = mynode()
         write(6,6) nid,' czfail:',(l(k),k=0,n)
    6    format(i5,a8,1p16e10.2)
         call exitt(1)
      endif
c
c     Sort eigenvalues and vectors
c
      call sort(atd,atu,n)
      call transpose(sft,n,sf,n)
      do j=1,n
         call swap(sft(1,j),atu,n,au)
      enddo
      call transpose(sf,n,sft,n)
c
c     Make "like" eigenvectors of same sign (for ease of diagnostics)
c
      do j=1,n
         avg = vlsum(sf(1,j),n)
         if (avg.lt.0) call chsign(sf(1,j),n)
      enddo
c
c     Clean up zero eigenvalue
c
      eps = 1.0e-6*dmax
      do i=1,n
c        if (atd(i).lt.eps)
c    $      write(6,*) 'Reset Eig in gen_Afem:',i,n,atd(i)
         if (atd(i).lt.eps) atd(i) = 0.0
      enddo
c
c     scale eigenvectors by C:
c
      do i=1,n
         do j=1,n
            sf(i,j) = sf(i,j)*c(i)
         enddo
      enddo
c                                        ^
c     Orthonormalize eigenvectors w.r.t. B inner-product
c
      do j=1,n
         alpha = vlsc3(bh,sf(1,j),sf(1,j),n)
         alpha = 1.0/sqrt(alpha)
         call cmult(sf(1,j),alpha,n)
      enddo
c
c     Diagnostics
c
c     do j=1,n
c        do i=1,n
c           sft(i,j) = vlsc3(bh,sf(1,i),sf(1,j),n)
c        enddo
c     enddo
c
c     n8 = min(n,8)
c     do i=1,n
c        write(6,2) (sft(i,j),j=1,n8)
c     enddo
c   2 format('Id:',1p8e12.4)
c
      call transpose(sft,n,sf,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ierr)
      integer                lbr,rbr,lbs,rbs,lbt,rbt,e,bsym
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'TSTEP'
c
      integer fbc(6)
c
c     ibc = 0  <==>  Dirichlet
c     ibc = 1  <==>  Dirichlet, outflow (no extension)
c     ibc = 2  <==>  Neumann,


      do iface=1,2*ndim
         ied = eface(iface)
         ibc = -1

         if (ifmhd) call mhd_bc_dn(ibc,iface,e) ! can be overwritten by 'mvn'

         if (cbc(ied,e,ifield).eq.'   ') ibc = 0
         if (cbc(ied,e,ifield).eq.'E  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'msi') ibc = 0
         if (cbc(ied,e,ifield).eq.'MSI') ibc = 0
         if (cbc(ied,e,ifield).eq.'P  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'p  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'O  ') ibc = 1
         if (cbc(ied,e,ifield).eq.'ON ') ibc = 1
         if (cbc(ied,e,ifield).eq.'o  ') ibc = 1
         if (cbc(ied,e,ifield).eq.'on ') ibc = 1
         if (cbc(ied,e,ifield).eq.'PEC') ibc = 1       ! Dirichlet
         if (cbc(ied,e,ifield).eq.'MS ') ibc = 1
         if (cbc(ied,e,ifield).eq.'ms ') ibc = 1
         if (cbc(ied,e,ifield).eq.'MM ') ibc = 1
         if (cbc(ied,e,ifield).eq.'mm ') ibc = 1
         if (cbc(ied,e,ifield).eq.'mv ') ibc = 2
         if (cbc(ied,e,ifield).eq.'mvn') ibc = 2
         if (cbc(ied,e,ifield).eq.'v  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'V  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'W  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'N  ') ibc = 2       ! Neumann
         if (cbc(ied,e,ifield).eq.'R  ') ibc = 2       ! Robin
         if (cbc(ied,e,ifield).eq.'SYM') ibc = bsym
         if (cbc(ied,e,ifield).eq.'SL ') ibc = 2
         if (cbc(ied,e,ifield).eq.'sl ') ibc = 2
         if (cbc(ied,e,ifield).eq.'SHL') ibc = 2
         if (cbc(ied,e,ifield).eq.'shl') ibc = 2
         if (cbc(ied,e,ifield).eq.'A  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'S  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'s  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'J  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'SP ') ibc = 0

         fbc(iface) = ibc

c        write(6,2) ibc,ied,e,ifield,cbc(ied,e,ifield)
c 2      format(2i3,i8,i3,2x,a3,'  get_fast_bc')


         if (ierr.eq.-1) write(6,1) ibc,ied,e,ifield,cbc(ied,e,ifield)
  1      format(2i3,i8,i3,2x,a3,'  get_fast_bc_error')

      enddo

      if (ierr.eq.-1) call exitti('Error A get_fast_bc$',e,1)

      lbr = fbc(1)
      rbr = fbc(2)
      lbs = fbc(3)
      rbs = fbc(4)
      lbt = fbc(5)
      rbt = fbc(6)

      ierr = 0
      if (ibc.lt.0) ierr = lglel(e)

c     write(6,6) e,lbr,rbr,lbs,rbs,(cbc(k,e,ifield),k=1,4)
c   6 format(i5,2x,4i3,3x,4(1x,a3),'  get_fast_bc')

      return
      end
c-----------------------------------------------------------------------
      subroutine outv(x,n,name3)
      character*3 name3
      real x(1)
c
      nn = min (n,10)
      write(6,6) name3,(x(i),i=1,nn)
    6 format(a3,10f12.6)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_fem_ax
     $   (s,lam,n,lbc,rbc,ll,lm,lr,z,y,nz,ie)
      real s(1),lam(1),ll,lm,lr,z(1),y(1)
      integer lbc,rbc
c
      parameter (m=100)
      real dx(0:m)
      integer icalld
      save    icalld
      data    icalld/0/
c
      icalld=icalld+1
c
      if (nz.gt.m-3) then
         write(6,*) 'ABORT. Error in set_up_fast_1D_fem. Increase m to'
     $             , nz
         call exitt(1)
      endif
c
c     In the present scheme, each element is viewed as a d-fold
c     tensor of (1+nz+1) arrays, even if funky bc's are applied
c     on either end of the 1D array.
c
      n = nz+2
c
c     Compute spacing, dx()
c
      call set_up_1D_geom_ax(dx,lbc,rbc,ll,lm,lr,z,y,nz)
c
      nn1 = n*n + 1
      call gen_eigs_A_fem_ax(s,s(nn1),lam,n,dx,y,lbc,rbc)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_1D_geom_ax(dx,lbc,rbc,ll,lm,lr,z,y,nz)
c
c     real dx(0:1),ll,lm,lr,z(2),y(1)
      real dx(0:1),ll,lm,lr,z(1),y(1) !FIXME misun 10/16/2014
      integer lbc,rbc
c
c
c     Set up the 1D geometry for the tensor-product based overlapping Schwarz
c
c     Upon return:
c
c       dx() contains the spacing required to set up the stiffness matrix.
c
c
c     Input:
c
c       lbc (rbc) is 0 if the left (right) BC is Dirichlet, 1 if Neumann.
c
c       ll is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the LEFT element
c
c       lm is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the MIDDLE element
c
c       lr is the space between the right-most Gauss point of the middle
c          element and the left-most Gauss point of the RIGHT element
c
c       --- if ll (lr) is very small (0), it indicates that there is no
c           left (right) spacing, and that the left (right) BC is Neumann.
c
c
c       z() is the array of nz Gauss points on the interval ]-1,1[.
c
c     Boundary conditions:
c
c     bc = 0  -- std. Dirichlet bc applied 2 points away from interior
c     bc = 1  -- Dirichlet bc applied 1 point away from interior (outflow)
c     bc = 2  -- Neumann bc applied on interior point (W,v,V,SYM,...)
c
c
c
c     Geometry:
c
c
c        dx0       dx1   dx2    dx3   dx5    dx5        dx6
c
c    bl        |<--ll-->|<------lm------>|<---lr--->|           br
c     0--------x-----|--x---x--------x---x--|-------x-----------0
c
c       left elem.         middle elem.         right elem.
c                   -1                     +1
c
c
c    "bl" = (extrapolated) location of Gauss point NX2-1 in left elem.
c
c    "br" = (extrapolated) location of Gauss point 2 in right elem.
c
c    Overlapping Schwarz applied with homogeneous Dirichlet boundary
c    conditions at "bl" and "br", and with a single d.o.f. extending
c    in to each adjacent domain.
c
      eps = 1.e-5
      call rone(dx,nz+3)
c
c     Middle
      scale = lm/(z(nz)-z(1))
      do i=1,nz-1
         dx(i+1) = (z(i+1)-z(i))*scale
      enddo

c     Left end
      if (lbc.eq.0) then
         dzm0   = z(1) + 1.
         dxm0   = scale*dzm0
         dxln   = ll - dxm0
         scalel = dxln/dzm0
         dx(0)  = scalel*(z(2)-z(1))
         dx(1)  = ll
      elseif (lbc.eq.1) then
         dx(1)  = ll
      endif
c
c     Right end
      if (rbc.eq.0) then
         dzm0      = z(1) + 1.
         dxm0      = scale*dzm0
         dxr0      = lr - dxm0
         scaler    = dxr0/dzm0
         dx(nz+1)  = lr
         dx(nz+2)  = scaler*(z(2)-z(1))
      elseif (rbc.eq.1) then
         dx(nz+1)  = lr
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_eigs_A_fem_ax(sf,sft,atd,n,l,y,lbc,rbc)
c
c     Set up Fast diagonalization solver for FEM on mesh 2
c
      real sf(n,n),sft(n,n),atd(1),l(0:1),y(1)
      integer lbc,rbc
c
      parameter (m=100)
      real atu(m),ad(m),au(m),c(m),bh(m),li(0:m)
c
      if (n.gt.m) then
         write(6,*) 'ABORT. Error in gen_eigs_A_fem. Increase m to',n
         call exitt(1)
      endif
c
c     Get delta x's
c
      do i=0,n
         li(i) = 1.0/l(i)
      enddo
c                          ^   ^
c     Fill initial arrays, A & B:
c
      call rzero(ad,n)
      call rzero(au,n-1)
      call rzero(bh,n)
c
      ie1 = lbc
      ien = n-rbc
      do ie=ie1,ien
c
c        il,ir are the left and right endpts of element ie.
         il = ie
         ir = ie+1
c
         if (ie.gt.0) ad(il) = ad(il) +       li(ie)
         if (ie.lt.n) ad(ir) = ad(ir) +       li(ie)
         if (ie.gt.0) au(il) =        -       li(ie)
         if (ie.gt.0) bh(il) = bh(il) + 0.5 * l(ie)
         if (ie.lt.n) bh(ir) = bh(ir) + 0.5 * l(ie)
      enddo
c
c     Take care of bc's (using blasting)
      bhm = vlmax(bh(2),n-2)/(n-2)
      ahm = vlmax(ad(2),n-2)/(n-2)
c
      if (lbc.gt.0) then
         au(1) = 0.
         ad(1) = ahm
         bh(1) = bhm
      endif
c
      if (rbc.gt.0) then
         au(n-1) = 0.
         ad(n  ) = ahm
         bh(n  ) = bhm
      endif
c
c
      do i=1,n
         c(i) = sqrt(1.0/bh(i))
      enddo
c                                        ~
c     Scale rows and columns of A by C:  A = CAC
c
      do i=1,n
         atd(i) = c(i)*ad(i)*c(i)
      enddo
c
c     Scale upper diagonal
c
      atu(1) = 0.
      do i=1,n-1
         atu(i) = c(i)*au(i)*c(i+1)
      enddo
c                                             ~
c     Compute eigenvalues and eigenvectors of A
c
      call calcz(atd,atu,n,dmax,dmin,sf,ierr)
      if (ierr.eq.1) then
         nid = mynode()
         write(6,6) nid,' czfail2:',(l(k),k=0,n)
    6    format(i5,a8,1p16e10.2)
         call exitt(1)
      endif
c
c     Sort eigenvalues and vectors
c
      call sort(atd,atu,n)
      call transpose(sft,n,sf,n)
      do j=1,n
         call swap(sft(1,j),atu,n,au)
      enddo
      call transpose(sf,n,sft,n)
c
c     Make "like" eigenvectors of same sign (for ease of diagnostics)
c
      do j=1,n
         avg = vlsum(sf(1,j),n)
         if (avg.lt.0) call chsign(sf(1,j),n)
      enddo
c
c     Clean up zero eigenvalue
c
      eps = 1.0e-6*dmax
      do i=1,n
c        if (atd(i).lt.eps)
c    $      write(6,*) 'Reset Eig in gen_Afm_ax:',i,n,atd(i)
         if (atd(i).lt.eps) atd(i) = 0.0
      enddo
c
c     scale eigenvectors by C:
c
      do i=1,n
         do j=1,n
            sf(i,j) = sf(i,j)*c(i)
         enddo
      enddo
c                                        ^
c     Orthonormalize eigenvectors w.r.t. B inner-product
c
      do j=1,n
         alpha = vlsc3(bh,sf(1,j),sf(1,j),n)
         alpha = 1.0/sqrt(alpha)
         call cmult(sf(1,j),alpha,n)
      enddo
c
c     Diagnostics
c
c     do j=1,n
c        do i=1,n
c           sft(i,j) = vlsc3(bh,sf(1,i),sf(1,j),n)
c        enddo
c     enddo
c
c     n8 = min(n,8)
c     do i=1,n
c        write(6,2) (sft(i,j),j=1,n8)
c     enddo
c   2 format('Id:',1p8e12.4)
c
      call transpose(sft,n,sf,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem(s,lam,n,lbc,rbc,ll,lm,lr,ie)
      include 'SIZE'
      include 'SEMHAT'
c
      common /fast1dsem/ g(lr2),w(lr2)
c
      real g,w
      real s(1),lam(1),ll,lm,lr
      integer lbc,rbc

      integer bb0,bb1,eb0,eb1,n,n1
      logical l,r

      n=nx1-1
      !bcs on E are from normal vel component
      if(lbc.eq.2 .or. lbc.eq.3) then !wall,sym - dirichlet velocity
         eb0=1
      else !outflow,element - neumann velocity
         eb0=0
      endif
      if(rbc.eq.2 .or. rbc.eq.3) then !wall,sym - dirichlet velocity
         eb1=n-1
      else !outflow,element - neumann velocity
         eb1=n
      endif
      !bcs on B are from tangent vel component
      if(lbc.eq.2) then !wall - dirichlet velocity
         bb0=1
      else !outflow,element,sym - neumann velocity
         bb0=0
      endif
      if(rbc.eq.2) then !wall - dirichlet velocity
         bb1=n-1
      else !outflow,element,sym - neumann velocity
         bb1=n
      endif
c
      l = (lbc.eq.0)
      r = (rbc.eq.0)
c
c     calculate E tilde operator
      call set_up_fast_1D_sem_op(s,eb0,eb1,l,r,ll,lm,lr,bh,dgl,0)
c     call outmat(s,n+1,n+1,'  Et  ',ie)
c     calculate B tilde operator
      call set_up_fast_1D_sem_op(g,bb0,bb1,l,r,ll,lm,lr,bh,jgl,1)
c     call outmat(g,n+1,n+1,'  Bt  ',ie)

      n=n+1
      call generalev(s,g,lam,n,w)
      if(.not.l) call row_zero(s,n,n,1)
      if(.not.r) call row_zero(s,n,n,n)
      call transpose(s(n*n+1),n,s,n) ! compute the transpose of s

      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem_op(g,b0,b1,l,r,ll,lm,lr,bh,jgl,jscl)
c            -1 T
c     G = J B  J
c
c     gives the inexact restriction of this matrix to
c     an element plus one node on either side
c
c     g - the output matrix
c     b0, b1 - the range for Bhat indices for the element
c              (enforces boundary conditions)
c     l, r - whether there is a left or right neighbor
c     ll,lm,lr - lengths of left, middle, and right elements
c     bh - hat matrix for B
c     jgl - hat matrix for J (should map vel to pressure)
c     jscl - how J scales
c            0: J = Jh
c            1: J = (L/2) Jh
c
c     result is inexact because:
c        neighbor's boundary condition at far end unknown
c        length of neighbor's neighbor unknown
c        (these contribs should be small for large N and
c         elements of nearly equal size)
c
      include 'SIZE'
      real g(0:lx1-1,0:lx1-1)
      real bh(0:lx1-1),jgl(1:lx2,0:lx1-1)
      real ll,lm,lr
      integer b0,b1
      logical l,r
      integer jscl
c
      real bl(0:lx1-1),bm(0:lx1-1),br(0:lx1-1)
      real gl,gm,gr,gll,glm,gmm,gmr,grr
      real fac
      integer n
      n=nx1-1
c
c
c     compute the scale factors for J
      if (jscl.eq.0) then
         gl=1.
         gm=1.
         gr=1.
      elseif (jscl.eq.1) then
         gl=0.5*ll
         gm=0.5*lm
         gr=0.5*lr
      endif
      gll = gl*gl
      glm = gl*gm
      gmm = gm*gm
      gmr = gm*gr
      grr = gr*gr
c
c     compute the summed inverse mass matrices for
c     the middle, left, and right elements
      do i=1,n-1
         bm(i)=2. /(lm*bh(i))
      enddo
      if (b0.eq.0) then
         bm(0)=0.5*lm*bh(0)
         if(l) bm(0)=bm(0)+0.5*ll*bh(n)
         bm(0)=1. /bm(0)
      endif
      if (b1.eq.n) then
         bm(n)=0.5*lm*bh(n)
         if(r) bm(n)=bm(n)+0.5*lr*bh(0)
         bm(n)=1. /bm(n)
      endif
c     note that in computing bl for the left element,
c     bl(0) is missing the contribution from its left neighbor
      if (l) then
         do i=0,n-1
            bl(i)=2. /(ll*bh(i))
         enddo
         bl(n)=bm(0)
      endif
c     note that in computing br for the right element,
c     br(n) is missing the contribution from its right neighbor
      if (r) then
         do i=1,n
            br(i)=2. /(lr*bh(i))
         enddo
         br(0)=bm(n)
      endif
c
      call rzero(g,(n+1)*(n+1))
      do j=1,n-1
         do i=1,n-1
            do k=b0,b1
               g(i,j) = g(i,j) + gmm*jgl(i,k)*bm(k)*jgl(j,k)
            enddo
         enddo
      enddo
c
      if (l) then
         do i=1,n-1
            g(i,0) = glm*jgl(i,0)*bm(0)*jgl(n-1,n)
            g(0,i) = g(i,0)
         enddo
c        the following is inexact
c        the neighbors bc's are ignored, and the contribution
c        from the neighbor's neighbor is left out
c        that is, bl(0) could be off as noted above
c        or maybe i should go from 1 to n
         do i=0,n
            g(0,0) = g(0,0) + gll*jgl(n-1,i)*bl(i)*jgl(n-1,i)
         enddo
      else
         g(0,0)=1.
      endif
c
      if (r) then
         do i=1,n-1
            g(i,n) = gmr*jgl(i,n)*bm(n)*jgl(1,0)
            g(n,i) = g(i,n)
         enddo
c        the following is inexact
c        the neighbors bc's are ignored, and the contribution
c        from the neighbor's neighbor is left out
c        that is, br(n) could be off as noted above
c        or maybe i should go from 0 to n-1
         do i=0,n
            g(n,n) = g(n,n) + grr*jgl(1,i)*br(i)*jgl(1,i)
         enddo
      else
         g(n,n)=1.
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine swap_lengths

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WZ'
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt

      real l,l2d
      integer e

      n2 = nx1-1
      nz0 = 1
      nzn = 1
      nx  = nx1-2
      if (if3d) then
         nz0 = 0
         nzn = n2
      endif
      call plane_space(lmr,lms,lmt,0,n2,wxm1,xm1,ym1,zm1,nx,n2,nz0,nzn)

      n=n2+1
      if (if3d) then
         do e=1,nelt
         do j=2,n2
         do k=2,n2
            l(1,k,j,e) = lmr(e)
            l(n,k,j,e) = lmr(e)
            l(k,1,j,e) = lms(e)
            l(k,n,j,e) = lms(e)
            l(k,j,1,e) = lmt(e)
            l(k,j,n,e) = lmt(e)
         enddo
         enddo
         enddo
         call dssum(l,n,n,n)
         do e=1,nelt
            llr(e) = l(1,2,2,e)-lmr(e)
            lrr(e) = l(n,2,2,e)-lmr(e)
            lls(e) = l(2,1,2,e)-lms(e)
            lrs(e) = l(2,n,2,e)-lms(e)
            llt(e) = l(2,2,1,e)-lmt(e)
            lrt(e) = l(2,2,n,e)-lmt(e)
         enddo
      else
         do e=1,nelt
         do j=2,n2
            l(1,j,1,e) = lmr(e)
            l(n,j,1,e) = lmr(e)
            l(j,1,1,e) = lms(e)
            l(j,n,1,e) = lms(e)
c           call outmat(l(1,1,1,e),n,n,' L    ',e)
         enddo
         enddo
c        call outmat(l(1,1,1,25),n,n,' L    ',25)
         call dssum(l,n,n,1)
c        call outmat(l(1,1,1,25),n,n,' L    ',25)
         do e=1,nelt
c           call outmat(l(1,1,1,e),n,n,' L    ',e)
            llr(e) = l(1,2,1,e)-lmr(e)
            lrr(e) = l(n,2,1,e)-lmr(e)
            lls(e) = l(2,1,1,e)-lms(e)
            lrs(e) = l(2,n,1,e)-lms(e)
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine row_zero(a,m,n,e)
      integer m,n,e
      real a(m,n)
      do j=1,n
         a(e,j)=0.
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mhd_bc_dn(ibc,face,e)
      integer                  face,e
c
c     sets Neumann BC on pressure (ibc=2) for face and e(lement) except
c     when ifield normal component has (homogeneous) Neumann
c     boundary condition setting ibc=1 (i.e. Direchlet BC on pressure)
c
c     Note: 'SYM' on a plane with r,s,t-normal is 'dnn','ndn','nnd'? bsym?
c
      include 'SIZE'
      include 'TOPOL'
      include 'INPUT'
      include 'TSTEP'

      ied = eface(face)	! symmetric -> preprocessor notation
      nfc = face+1
      nfc = nfc/2	! = 1,2,3 for face 1 & 2,3 & 4,5 & 6

      if (indx1(cbc(ied,e,ifield),'d',1).gt.0)   ibc=2

      if (indx1(cbc(ied,e,ifield),'n',1).gt.nfc) ibc=1 ! 'n' for V_n

      return
      end
c-----------------------------------------------------------------------
      subroutine local_solves_fdm(u,v)
c
c     Given an input vector v, this returns the additive Schwarz solution
c
c                       -1
c                  u = M   v
c
c                      T     -1
c     where M = sum ( R_i A_i    R_i )
c                i
c
c     The R_i's are simply index_set restriction operators.
c
c     The local solves are performed with the fast diagonalization method.
c
      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'PARALLEL'
c
      include 'TSTEP'
      include 'CTIMER'
c
      real u(lx2,ly2,lz2,lelv),v(lx2,ly2,lz2,lelv)
      common /scrpre/ v1(lx1,ly1,lz1,lelv)
     $               ,w1(lx1,ly1,lz1),w2(lx1,ly1,lz1)
      common /scrover/ ar(lelv)

      parameter(lxx=lx1*lx1, levb=lelv+lbelv)
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)
      integer e,eb,eoff

c
      if (icalld.eq.0) tsolv=0.0
      icalld=icalld+1
      nsolv=icalld
c
      ntot1 = nx1*ny1*nz1*nelt
      ntot2 = nx2*ny2*nz2*nelt
c
c     Fill interiors
      iz1 = 0
      if (if3d) iz1=1
      call rzero(v1,ntot1)
      do e=1,nelt
         do iz=1,nz2
         do iy=1,ny2
         do ix=1,nx2
            v1(ix+1,iy+1,iz+iz1,e) = v(ix,iy,iz,e)
         enddo
         enddo
         enddo
      enddo
      call dface_ext    (v1)
      call dssum        (v1,nx1,ny1,nz1)
      call dface_add1si (v1,-1.)
c
c     Now solve each subdomain problem:
c
      etime1=dnekclock()

      eoff  = 0
      if (ifield.gt.1) eoff  = nelt

      do e = 1,nelt
         eb = e + eoff
         call fastdm1(v1(1,1,1,e),df(1,eb)
     $                           ,sr(1,eb),ss(1,eb),st(1,eb),w1,w2)
      enddo
      tsolv=tsolv+dnekclock()-etime1
c
c     Exchange/add elemental solutions
c
      call s_face_to_int (v1,-1.)
      call dssum         (v1,nx1,ny1,nz1)
      call s_face_to_int (v1, 1.)
      if(param(42).eq.0) call do_weight_op(v1)
c
c     Map back to pressure grid (extract interior values)
c
      do e=1,nelt
         do iz=1,nz2
         do iy=1,ny2
         do ix=1,nx2
            u(ix,iy,iz,e) = v1(ix+1,iy+1,iz+iz1,e)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fastdm1(r,df,sr,ss,st,w1,w2)
c
c     Fast diagonalization solver for FEM on mesh 1
c
      include 'SIZE'
      parameter (lxx=lx1*lx1,lxyz=lx1*ly1*lz1)

      real r(1),df(1),sr(lxx,2),ss(lxx,2),st(lxx,2),w1(1),w2(1)
c
c
c      T
c     S  r
      call tensr3 (w1,nx1,r ,nx1,sr(1,2),ss(1,1),st(1,1),w2)
c
c
c      -1 T
c     D  S  r
c
      call col2   (w1,df,lxyz)
c
c
c        -1 T
c     S D  S  r
c
      call tensr3 (r ,nx1,w1,nx1,sr(1,1),ss(1,2),st(1,2),w2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine tensr3(v,nv,u,nu,A,Bt,Ct,w)
C
C     -  Tensor product application of v = (C x B x A) u
C        NOTE -- the transpose of B & C must be input, rather than B & C.
C
C     -  scratch arrays: w(nu*nu*nv)
C
C
      include 'SIZE'
      include 'INPUT'
      real v(nv,nv,nv),u(nu,nu,nu)
      real A(1),Bt(1),Ct(1)
      real w(1)

      if (nu.gt.nv) then
         write(6,*) nid,nu,nv,' ERROR in tensr3. Contact P.Fischer.'
         write(6,*) nid,nu,nv,' Memory problem.'
         call exitt(1)
      endif

      if (if3d) then
         nuv = nu*nv
         nvv = nv*nv
         call mxm(A,nv,u,nu,v,nu*nu)
         k=1
         l=1
         do iz=1,nu
            call mxm(v(k,1,1),nv,Bt,nu,w(l),nv)
            k=k+nuv
            l=l+nvv
         enddo
         call mxm(w,nvv,Ct,nu,v,nv)
      else
         call mxm(A,nv,u,nu,w,nu)
         call mxm(w,nv,Bt,nu,v,nv)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine s_face_to_int(x,c)
c
c     Scale face and add to interior of element
c
      include 'SIZE'
      include 'INPUT'
      real x(nx1,ny1,nz1,1)
c
      do ie=1,nelt
c
        if (if3d) then
c
          do iz=2,nz1-1
          do ix=2,nx1-1
            x(ix,2    ,iz,ie) = c*x(ix,1  ,iz,ie) + x(ix,2    ,iz,ie)
            x(ix,ny1-1,iz,ie) = c*x(ix,ny1,iz,ie) + x(ix,ny1-1,iz,ie)
          enddo
          enddo
c
          do iz=2,nz1-1
          do iy=2,ny1-1
            x(2    ,iy,iz,ie) = c*x(1  ,iy,iz,ie) + x(2    ,iy,iz,ie)
            x(nx1-1,iy,iz,ie) = c*x(nx1,iy,iz,ie) + x(nx1-1,iy,iz,ie)
          enddo
          enddo
c
          do iy=2,ny1-1
          do ix=2,nx1-1
            x(ix,iy,2    ,ie) = c*x(ix,iy,1  ,ie) + x(ix,iy,2    ,ie)
            x(ix,iy,nz1-1,ie) = c*x(ix,iy,nz1,ie) + x(ix,iy,nz1-1,ie)
          enddo
          enddo
c
        else
c         2D
          do ix=2,nx1-1
            x(ix,2    ,1,ie) = c*x(ix,1  ,1,ie) + x(ix,2    ,1,ie)
            x(ix,ny1-1,1,ie) = c*x(ix,ny1,1,ie) + x(ix,ny1-1,1,ie)
          enddo
          do iy=2,ny1-1
            x(2    ,iy,1,ie) = c*x(1  ,iy,1,ie) + x(2    ,iy,1,ie)
            x(nx1-1,iy,1,ie) = c*x(nx1,iy,1,ie) + x(nx1-1,iy,1,ie)
          enddo
        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine dface_ext(x)
c     Extend interior to face of element
c
      include 'SIZE'
      include 'INPUT'
      real x(nx1,ny1,nz1,1)
c
      do ie=1,nelt
c
         if (if3d) then
c
          do iz=2,nz1-1
          do ix=2,nx1-1
            x(ix,1  ,iz,ie) = x(ix,2    ,iz,ie)
            x(ix,ny1,iz,ie) = x(ix,ny1-1,iz,ie)
          enddo
          enddo
c
          do iz=2,nz1-1
          do iy=2,ny1-1
            x(1  ,iy,iz,ie) = x(2    ,iy,iz,ie)
            x(nx1,iy,iz,ie) = x(nx1-1,iy,iz,ie)
          enddo
          enddo
c
          do iy=2,ny1-1
          do ix=2,nx1-1
            x(ix,iy,1  ,ie) = x(ix,iy,2    ,ie)
            x(ix,iy,nz1,ie) = x(ix,iy,nz1-1,ie)
          enddo
          enddo
c
        else
c
          do ix=2,nx1-1
            x(ix,1  ,1,ie) = x(ix,2    ,1,ie)
            x(ix,ny1,1,ie) = x(ix,ny1-1,1,ie)
          enddo
          do iy=2,ny1-1
            x(1  ,iy,1,ie) = x(2    ,iy,1,ie)
            x(nx1,iy,1,ie) = x(nx1-1,iy,1,ie)
          enddo
c
        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine dface_add1si(x,c)
c     Scale interior and add to face of element
c
      include 'SIZE'
      include 'INPUT'
      real x(nx1,ny1,nz1,1)
c
      do ie=1,nelt
c
        if (if3d) then
c
          do iz=2,nz1-1
          do ix=2,nx1-1
            x(ix,1  ,iz,ie) = x(ix,1  ,iz,ie) + c*x(ix,2    ,iz,ie)
            x(ix,ny1,iz,ie) = x(ix,ny1,iz,ie) + c*x(ix,ny1-1,iz,ie)
          enddo
          enddo
c
          do iz=2,nz1-1
          do iy=2,ny1-1
            x(1  ,iy,iz,ie) = x(1  ,iy,iz,ie) + c*x(2    ,iy,iz,ie)
            x(nx1,iy,iz,ie) = x(nx1,iy,iz,ie) + c*x(nx1-1,iy,iz,ie)
          enddo
          enddo
c
          do iy=2,ny1-1
          do ix=2,nx1-1
            x(ix,iy,1  ,ie) = x(ix,iy,1  ,ie) + c*x(ix,iy,2    ,ie)
            x(ix,iy,nz1,ie) = x(ix,iy,nz1,ie) + c*x(ix,iy,nz1-1,ie)
          enddo
          enddo
c
        else
c
c         2D
c
          do ix=2,nx1-1
            x(ix,1  ,1,ie) = x(ix,1  ,1,ie) + c*x(ix,2    ,1,ie)
            x(ix,ny1,1,ie) = x(ix,ny1,1,ie) + c*x(ix,ny1-1,1,ie)
          enddo
          do iy=2,ny1-1
            x(1  ,iy,1,ie) = x(1  ,iy,1,ie) + c*x(2    ,iy,1,ie)
            x(nx1,iy,1,ie) = x(nx1,iy,1,ie) + c*x(nx1-1,iy,1,ie)
          enddo
c
        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine init_weight_op

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      parameter(levb=lelv+lbelv)
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      common /weightop/ w(lx2,lz2,2,3,levb)
      real l
      real w
      integer e,e0,eb

      e0  = 0
      if (ifield.gt.1) e0 = nelt

      n=nx2+1
      if (if3d) then
         do e=1,nelt
            call rzero(l(1,1,1,e),nx1*ny1*nz1)
            do k=2,n
               do j=2,n
                  l(2,j,k,e)=1
                  l(n,j,k,e)=1
               enddo
            enddo
            do k=2,n
               do i=2,n
                  l(i,2,k,e)=1
                  l(i,n,k,e)=1
               enddo
            enddo
            do j=2,n
               do i=2,n
                  l(i,j,2,e)=1
                  l(i,j,n,e)=1
               enddo
            enddo
         enddo
      else
         do e=1,nelt
            call rzero(l(1,1,1,e),nx1*ny1*nz1)
            do j=2,n
               l(2,j,1,e)=1
               l(n,j,1,e)=1
            enddo
            do i=2,n
               l(i,2,1,e)=1
               l(i,n,1,e)=1
            enddo
         enddo
      endif

      call dface_ext(l)
      call dssum(l,nx1,ny1,nz1)
      call dface_add1si(l,-1.)
      call s_face_to_int(l,1.)
c     l now holds the count matrix C on the outer pressure nodes
      if (if3d) then
         do e=1,nelt
            eb = e0+e
            do k=1,nz2
               do j=1,ny2
                  w(j,k,1,1,eb)=1.0/l(2,j+1,k+1,e)
                  w(j,k,2,1,eb)=1.0/l(n,j+1,k+1,e)
               enddo
            enddo
            do k=1,nz2
               do i=1,nx2
                  w(i,k,1,2,eb)=1.0/l(i+1,2,k+1,e)
                  w(i,k,2,2,eb)=1.0/l(i+1,n,k+1,e)
               enddo
            enddo
            do j=1,ny2
               do i=1,nx2
                  w(i,j,1,3,eb)=1.0/l(i+1,j+1,2,e)
                  w(i,j,2,3,eb)=1.0/l(i+1,j+1,n,e)
               enddo
            enddo
         enddo
      else
         do e=1,nelt
            eb = e0+e
            do j=1,ny2
               w(j,1,1,1,eb)=1.0/l(2,j+1,1,e)
               w(j,1,2,1,eb)=1.0/l(n,j+1,1,e)
            enddo
            do i=1,nx2
               w(i,1,1,2,eb)=1.0/l(i+1,2,1,e)
               w(i,1,2,2,eb)=1.0/l(i+1,n,1,e)
            enddo
         enddo
      endif
      end
c-----------------------------------------------------------------------
      subroutine do_weight_op(x)
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      parameter(levb=lelv+lbelv)
      common /weightop/ w(lx2,lz2,2,3,levb)
      real w

      real x(0:nx1-1,0:ny1-1,0:nz1-1,1)
      integer e,e0,eb

      e0  = 0
      if (ifield.gt.1) e0 = nelt

      if (if3d) then
         do e=1,nelt
            eb = e0 + e
            do k=1,nz2
               do j=1,ny2
                  x(  1,j,k,e)=w(j,k,1,1,eb)*x(  1,j,k,e)
                  x(nx2,j,k,e)=w(j,k,2,1,eb)*x(nx2,j,k,e)
               enddo
            enddo
            do k=1,nz2
               do i=2,nx2-1
                  x(i,  1,k,e)=w(i,k,1,2,eb)*x(i,  1,k,e)
                  x(i,ny2,k,e)=w(i,k,2,2,eb)*x(i,ny2,k,e)
               enddo
            enddo
            do j=2,ny2-1
               do i=2,nx2-1
                  x(i,j,  1,e)=w(i,j,1,3,eb)*x(i,j,  1,e)
                  x(i,j,nz2,e)=w(i,j,2,3,eb)*x(i,j,nz2,e)
               enddo
            enddo
         enddo
      else
         do e=1,nelt
            eb = e0 + e
            do j=1,ny2
               x(  1,j,0,e)=w(j,1,1,1,eb)*x(  1,j,0,e)
               x(nx2,j,0,e)=w(j,1,2,1,eb)*x(nx2,j,0,e)
            enddo
            do i=2,nx2-1
               x(i,  1,0,e)=w(i,1,1,2,eb)*x(i,  1,0,e)
               x(i,ny2,0,e)=w(i,1,2,2,eb)*x(i,ny2,0,e)
            enddo
         enddo
      endif
      end
c-----------------------------------------------------------------------
      subroutine calcz(d,e,n,dmax,dmin,z,ierr)
c
c     Num. Rec. 2nd Ed., p.473
c
c     Note:  d(1:n) is the diagonal of the sym. tridiagonal matrix
c            e(1:n) is the upper diagonal of the tridiagonal matrix,
c                   WITH e(n) ARBITRARY (a slight shift from Num.Rec.)
c
c            z(n:n) is the packed array of eigenvectors
c
      real  d(n),e(n),z(n,n)
      real  smalln,small
c
      call ident(z,n)
      one = 1.
c
c     Find smallest number  (pff mod to Num. Rec. 2nd Ed., p.473)
c
      small = 0.5
      do i = 1,100
         smalln = small * small
         if (smalln .eq.0) then
            do j=1,1000
               smalln = 0.5*small
               if (smalln .eq.0) goto 10
               small = smalln
            enddo
         endif
         small = smalln
      enddo
   10 continue
      small = 10.*small
      small = max(small,1e-15)

c     write(6,*) 'this is small:',small
c
      do 15 l=1,n
         iter = 0
c
    1    do 12 m=l,n-1
            dd = abs( d(m) ) + abs( d(m+1) )
            de = e(m) + dd
            df = abs(dd - de)
c           write(6,112) iter,m,'de:',dd,de,df,small
            if ( df .le. small ) goto 2
   12    continue
  112    format(i3,i9,1x,a3,1p4e16.8)
         m = n
c
    2    if ( m .ne. l ) then
            if ( iter .eq. 600 ) then
               write (6,*) 'too many iterations in calc'
c              n10 = min(n,10)
c              do i=1,n
c                 write(6,9) d(i),(z(i,j),j=1,n10)
c              enddo
c   9          format(1pe12.4,' e:',1p10e12.4)
c              call exitt(1)
               ierr=1
               return
            endif
c
            iter = iter + 1
            g = ( d(l+1) - d(l) ) / ( 2.0 * e(l) )
            r = pythag(g,one)
            g = d(m) - d(l) + e(l)/(g+sign(r,g))
            s = 1.0
            c = 1.0
            p = 0.0
c
            do 14 i = m-1,l,-1
               f = s * e(i)
               b = c * e(i)
               r = pythag(f,g)
               e(i+1) = r
               if ( abs(r) .le. small ) then
                  d(i+1) = d(i+1) - p
                  e(m)   = 0.
                  goto 1
               endif
               s = f/r
               c = g/r
               g = d(i+1) - p
               r = ( d(i)-g )*s + 2.*c*b
               p = s*r
               d(i+1) = g + p
               g = c*r - b
c      ...     find eigenvectors ... (new, 11/19/94, pff, p.363. Num.Rec.I.)
               do 13 k=1,n
                  f = z(k,i+1)
                  z(k,i+1) = s*z(k,i)+c*f
                  z(k,i  ) = c*z(k,i)-s*f
   13          continue
c      ...     end of eigenvector section ...
   14       continue
c
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0
            goto 1
         endif
c
   15 continue
c
c     write (8,8) (d(j),j=1,n)
c   8 format('eig:',8f10.4)
c
      dmin = d(1)
      dmax = d(1)
      do 40 i = 1 , n
        dmin = min( d(i) , dmin )
        dmax = max( d(i) , dmax )
   40 continue
c
c     Output eigenvectors
c
c     n10 = min(n,10)
c     do i=1,n
c        write(6,9) d(i),(z(i,j),j=1,n10)
c     enddo
c   9 format(1pe12.4,' e:',1p10e12.4)
c
      ierr=0
      return
      end
c-----------------------------------------------------------------------
      function pythag(a,b)
c
c     Compute sqrt(a*a + b*b) w/o under- or over-flow.
c
      absa=abs(a)
      absb=abs(b)
      if (absa.gt.absb) then
         pythag = absa*sqrt(1. + (absb/absa)**2 )
      else
         if (absb.eq.0.) then
            pythag = 0.
         else
            pythag = absb*sqrt(1. + (absa/absb)**2 )
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine swap(b,ind,n,temp)
      real B(1),TEMP(1)
      integer IND(1)
C***
C***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
C***  INTO ITEM(I), WHERE JJ=IND(I).
C***
      DO 20 I=1,N
         JJ=IND(I)
         TEMP(I)=B(JJ)
   20 CONTINUE
      DO 30 I=1,N
   30 B(I)=TEMP(I)
      RETURN
      END
c=======================================================================
      subroutine gradl_rst_t(u,ur,us,ut,md,if3d)  ! GLL grad-transpose
c
c     Thus routine originally from fsi file: u5.usr (May 2010)
c
c
      include 'SIZE'
      include 'DXYZ'

      real    u(1),ur(1),us(1),ut(1)
      logical if3d

c     dgradl holds GLL-based derivative / interpolation operators

      parameter (ldg=lxd**3,lwkd=2*ldg)
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      m0 = md-1
      call get_dgll_ptr (ip,md,md)
      if (if3d) then
         call local_grad3_t(u,ur,us,ut,m0,1,d(ip),dt(ip),wkd)
      else
         call local_grad2_t(u,ur,us   ,m0,1,d(ip),dt(ip),wkd)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gradl_rst(ur,us,ut,u,md,if3d)  ! GLL-based gradient
c
      include 'SIZE'
      include 'DXYZ'

      real    ur(1),us(1),ut(1),u(1)
      logical if3d

c     dgradl holds GLL-based derivative / interpolation operators

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      m0 = md-1
      call get_dgll_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,u,m0,1,d(ip),dt(ip))
      else
         call local_grad2(ur,us   ,u,m0,1,d(ip),dt(ip))
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t(u,ur,us,ut,N,e,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,e,D,Dt
      real u (0:N,0:N,0:N,1)
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

      call mxm(Dt,m1,ur,m1,u(0,0,0,e),m2)

      do k=0,N
         call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u(0,0,0,e),w,m3)

      call mxm(ut,m2,D ,m1,w,m1)
      call add2(u(0,0,0,e),w,m3)

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad2_t(u,ur,us,N,e,D,Dt,w)
c     Output: ur,us         Input:u,N,e,D,Dt
      real u (0:N,0:N,1)
      real ur(0:N,0:N),us(0:N,0:N),ut(0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1

      call mxm(Dt,m1,ur,m1,u(0,0,e),m1)
      call mxm(us,m1,D ,m1,w         ,m1)
      call add2(u(0,0,e),w,m2)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgll_ptr (ip,mx,md)
c
c     Get pointer to GLL-GLL interpolation dgl() for pair (mx,md)
c
      include 'SIZE'

c     dgradl holds GLL-based derivative / interpolation operators

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgradl/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

c     Pointers into GLL-based derivative / interpolation operators

      parameter (ld=2*lxd)
      common /jgradl/ pd    (0:ld*ld)
     $              , pdg   (0:ld*ld)
     $              , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl

      ij = md + ld*(mx-1)
      ip = pdg (ij)

      if (ip.eq.0) then

         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)

        if (nid.eq.985) write(6,*) nstore,ldg,ij,md,mx,ip,' NSTOR'
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c
         call gen_dgll(d (ip),dt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_dgll(dgl,dgt,mp,np,w)
c
c     Generate derivative from np GL points onto mp GL points
c
c        dgl  = derivative matrix, mapping from velocity nodes to pressure
c        dgt  = transpose of derivative matrix
c        w    = work array of size (3*np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
c
      real dgl(mp,np),dgt(np*mp),w(1)
c
c
      iz = 1
      id = iz + np
c
      call zwgll (w(iz),dgt,np)  ! GL points
      call zwgll (w(id),dgt,mp)  ! GL points
c
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
         enddo
      enddo
c
      call transpose(dgt,np,dgl,mp)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_h1_crs(h1,h2)

      include 'SIZE'
      include 'GEOM'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      include 'ZPER'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      real h1(1),h2(1)

      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex

      integer gs_handle
      integer null_space,e

      character*3 cb
      common /scrxxt/ cmlt(lcr,lelv),mask(lcr,lelv)
      real mask

      common /scrch/   a(2*lx1*ly1*lz1*lelv)
      common /scrxxti/ ia(lcr,lcr,lelv), ja(lcr,lcr,lelv)

      common /scrmgx/ w1(lx1*ly1*lz1*lelv),w2(lx1*ly1*lz1*lelv)

      integer*8 ngv

      t0 = dnekclock()

      nxc     = 2

      if (nid.eq.0) write(6,*) 'setup h1 coarse grid, nxc=', nxc

      ncr     = nxc**ndim
      nxyz_c  = ncr
c
c     Set SEM_to_GLOB
c
      if (nid.eq.0) write(6,*) 'hsmg_setup_dssum, ifgtp=',ifgtp
      if (ifgtp) then
          if (nid.eq.0) write(6,*) 'use regular mesh for multigrid'
          call gen_gtp_vertex    (vertex, ncrnr)
c         call exitt(1)
      else
          call get_vert
      endif

      call set_vert2(se_to_gcrs,ngv,nxc,nelv,vertex,.true.)

c     Set mask
      z=0
      ntot=nelv*nxyz_c
      nzc=1
      if (if3d) nzc=nxc
      call rone(mask,ntot)
      call rone(cmlt,ntot)
      nfaces=2*ndim
c     ifield=1			!c? avo: set in set_overlap through 'TSTEP'?

      if (ifield.eq.1) then
         do ie=1,nelv
         do iface=1,nfaces
            cb=cbc(iface,ie,ifield)
            if (cb.eq.'O  '  .or.  cb.eq.'ON '  .or.  cb.eq.'MM '  .or.
     $          cb.eq.'mm '  .or.  cb.eq.'ms '  .or.  cb.eq.'MS '  .or.
     $          cb.eq.'PEC'                                      )
     $          call facev(mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
         enddo
         enddo
      endif

c     Set global index of dirichlet nodes to zero; xxt will ignore them

      call gs_setup(gs_handle,se_to_gcrs,ntot,nekcomm,np)
      call gs_op   (gs_handle,mask,1,2,0)  !  "*"
      call gs_op   (gs_handle,cmlt,1,1,0)  !  "+"
      call gs_free (gs_handle)
      call set_jl_crs_mask(ntot,mask,se_to_gcrs)

      call invcol1(cmlt,ntot)

c     Setup local SEM-based Neumann operators (for now, just full...)

c     NOTE: a(),h1,...,w2() must all be large enough
      n = nx1*ny1*nz1*nelv
      call get_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)

      call set_mat_ij(ia,ja,ncr,nelv)

      null_space=0
c     if (ifield.eq.1) then
c        if (ifvcor)  null_space=1
c     elseif (ifield.eq.ifldmhd) then
c        if (ifbcor)  null_space=1
c     endif

      nz=ncr*ncr*nelv
      call crs_setup(xxth(ifield),nekcomm,np, ntot,se_to_gcrs,
     $               nz,ia,ja,a, null_space)
c     call crs_stats(xxth(ifield))

      t0 = dnekclock()-t0
      if (nid.eq.0) then
         write(6,*) 'done :: setup h1 coarse grid ',t0, ' sec'
         write(6,*) ' '
      endif

      return
      end
c
c-----------------------------------------------------------------------
      subroutine set_jl_crs_mask(n, mask, se_to_gcrs)
      real mask(1)
      integer*8 se_to_gcrs(1)
      do i=1,n
         if(mask(i).lt.0.1) se_to_gcrs(i)=0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine set_mat_ij(ia,ja,n,ne)
      integer n,ne
      integer ia(n,n,ne), ja(n,n,ne)
c
      integer i,j,ie
      do ie=1,ne
      do j=1,n
      do i=1,n
         ia(i,j,ie)=(ie-1)*n+i-1
         ja(i,j,ie)=(ie-1)*n+j-1
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine get_local_crs_galerkin(a,ncl,nxc,h1,h2,w1,w2)

c     This routine generates Nelv submatrices of order ncl using
c     Galerkin projection

      include 'SIZE'

      real    a(ncl,ncl,1),h1(1),h2(1)
      real    w1(nx1*ny1*nz1,nelv),w2(nx1*ny1*nz1,nelv)

      parameter (lcrd=lx1**ldim)
      common /ctmp1z/ b(lcrd,8)

      integer e

      do j=1,ncl
         call gen_crs_basis(b(1,j),j) ! bi- or tri-linear interpolant
      enddo

      isd  = 1
      imsh = 1

      nxyz = nx1*ny1*nz1
      do j = 1,ncl
         do e = 1,nelv
            call copy(w1(1,e),b(1,j),nxyz)
         enddo

         call axhelm3(w2,w1,h1,h2,imsh,isd)        ! A^e * bj

         do e = 1,nelv
         do i = 1,ncl
            a(i,j,e) = vlsc2(b(1,i),w2(1,e),nxyz)  ! bi^T * A^e * bj
         enddo
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_crs_basis(b,j) ! bi- tri-linear

      include 'SIZE'
      real b(nx1,ny1,nz1)

      real z0(lx1),z1(lx1)
      real zr(lx1),zs(lx1),zt(lx1)

      integer p,q,r

      call zwgll(zr,zs,nx1)

      do i=1,nx1
         z0(i) = .5*(1-zr(i))  ! 1-->0
         z1(i) = .5*(1+zr(i))  ! 0-->1
      enddo

      call copy(zr,z0,nx1)
      call copy(zs,z0,nx1)
      call copy(zt,z0,nx1)

      if (mod(j,2).eq.0)                        call copy(zr,z1,nx1)
      if (j.eq.3.or.j.eq.4.or.j.eq.7.or.j.eq.8) call copy(zs,z1,nx1)
      if (j.gt.4)                               call copy(zt,z1,nx1)

      if (ndim.eq.3) then
         do r=1,nx1
         do q=1,nx1
         do p=1,nx1
            b(p,q,r) = zr(p)*zs(q)*zt(r)
         enddo
         enddo
         enddo
      else
         do q=1,nx1
         do p=1,nx1
            b(p,q,1) = zr(p)*zs(q)
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_crs_basis2(b,j) ! bi- tri-quadratic

      include 'SIZE'
      real b(nx1,ny1,nz1)

      real z0(lx1),z1(lx1),z2(lx1)
      real zr(lx1),zs(lx1),zt(lx1)

      integer p,q,r

      call zwgll(zr,zs,nx1)

      do i=1,nx1
         z0(i) = .5*(zr(i)-1)*zr(i)  ! 1-->0   ! Lagrangian, ordered
         z1(i) = 4.*(1+zr(i))*(1-zr(i))        ! lexicographically
         z2(i) = .5*(zr(i)+1)*zr(i)  ! 0-->1   !
      enddo

      call copy(zr,z0,nx1)
      call copy(zs,z0,nx1)
      call copy(zt,z0,nx1)

      if (mod(j,2).eq.0)                        call copy(zr,z1,nx1)
      if (j.eq.3.or.j.eq.4.or.j.eq.7.or.j.eq.8) call copy(zs,z1,nx1)
      if (j.gt.4)                               call copy(zt,z1,nx1)

      if (ndim.eq.3) then
         do r=1,nx1
         do q=1,nx1
         do p=1,nx1
            b(p,q,r) = zr(p)*zs(q)*zt(r)
         enddo
         enddo
         enddo
      else
         do q=1,nx1
         do p=1,nx1
            b(p,q,1) = zr(p)*zs(q)
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cki(name3,x,n,ii)
      include 'SIZE'
      character*3 name3
      real x(n)

      return

      xx = 0.
      do i=1,n
         xx=xx+x(i)*x(i)
      enddo

      xx   =  glsum(xx,1)
      ntot = iglsum(n ,1)

      xx = xx/ntot
      if (xx.gt.0) xx = sqrt(xx)

      xmn = glmin (x,n)
      xmx = glmax (x,n)
      xma = glamax(x,n)

      if (nid.eq.0) write(6,1) name3,ii,n,xmn,xx,xmx,xma
    1 format('cki ',a3,2i8,1p4e16.8)

      if (ii.eq.0) write(6,*) 'Quit in cki.'
c     nz1=1/(nx1-ny1)
c     if (ii.eq.0) stop

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_mg(h1,h2)
      include 'SIZE'
      include 'TOTAL'
      real h1(1),h2(1)

      call ifill (nelfld,nelv,ldimt)

      call load_semhat_weighted    !   Fills the SEMHAT arrays
      call swap_lengths
      call hsmg_setup(h1,h2)
      call h1mg_setup(h1,h2)

      call set_up_h1_crs(h1,h2)

      return
      end
c-----------------------------------------------------------------------
