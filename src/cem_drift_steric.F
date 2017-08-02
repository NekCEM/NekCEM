c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_steric
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'
      include 'STERIC'

      integer  i,j,i0,isd,imsh,maxit
      real     wk1(lpts),wk2(lpts), tol,glamax
      real     diff_n_dt(lpts), diff_p_dt(lpts),diff_e_dt(lpts)
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/
      
      real wa (lx1*ly1*lz1*lelt)
      real wb (lx1*ly1*lz1*lelt)
      
      real wa_n (lx1*ly1*lz1*lelt)
      real wa_p (lx1*ly1*lz1*lelt)
      real wb_n (lx1*ly1*lz1*lelt)
      real wb_p (lx1*ly1*lz1*lelt)

      npts= nx1*ny1*nz1*nelt 
                             
      if (icalled.eq.0) then

        call cmult2 (h1c(1,1),diff_n,dt,npts) ! h1n= h1n * (diff_n*dt)
        call rone (h2c(1,1),npts)             ! h2n= 1.0

        call cmult2 (h1c(1,2),diff_p,dt,npts) ! h1p= h1p * (diff_p*dt)
        call rone (h2c(1,2),npts)             ! h2p= 1.0

        call cmult2 (h1cc(1,1,1),h1c(1,1),1.0+ste_ac(1,1),npts) !steric control term
        call cmult2 (h1cc(1,1,2),h1c(1,1),    ste_ac(1,2),npts)
        call cmult2 (h1cc(1,2,1),h1c(1,2),    ste_ac(2,1),npts)
        call cmult2 (h1cc(1,2,2),h1c(1,2),1.0+ste_ac(2,2),npts)

        call copy (h2cc(1,1,1),h2c(1,1),npts)
        call rzero(h2cc(1,1,2),npts)
        call rzero(h2cc(1,2,1),npts)
        call copy (h2cc(1,2,2),h2c(1,2),npts)

c        icalled = 1! If not Dirichlet, need uptate each loop 

       endif

c...  set cggo parameters

      isd = 1
      imsh = 1
      maxit = 5000
      tol = param(22)

c...  get rhs for cN and cP

      call cem_drift_rhs_bdf1_steric

C     boundary surface integral
      call bcneusc0(resN,mu_n,cn,npts,2,1)
      call bcneusc0(resP,mu_p,cp,npts,3,1)

      call cem_drift_grad(gradc(1,1,1),gradc(1,2,1),gradc(1,3,1),cN)
      call cem_drift_grad(gradc(1,1,2),gradc(1,2,2),gradc(1,3,2),cP)
      call bcneusc_steric(resN,mu_n,npts,2,1)
      call bcneusc_steric(resP,mu_p,npts,3,1)

C     Robin, add to LHS
      call bcneusc(wa_n,npts,2,-1)
      call col2 (wa_n,h1cc(1,1,1),npts)
      call add2 (h2cc(1,1,1),wa_n,npts)

      call bcneusc(wa_p,npts,3,-1)
      call col2 (wa_p,h1cc(1,2,2),npts)
      call add2 (h2cc(1,2,2),wa_p,npts)

C     Neumann + Robin, add to RHS
      call bcneusc(wa_n,npts,2,1)
      call col2 (wa_n,h1cc(1,1,1),npts)
      call add2 (resN,wa_n,npts)

      call bcneusc(wa_p,npts,3,1)
      call col2 (wa_p,h1cc(1,2,2),npts)
      call add2 (resP,wa_p,npts)

C     Dirichlet seperate
      call bcdirsc (cn,npts,2)
      call bcdirsc (cp,npts,3)
      
      call axhelm_steric(wa_n,wa_p,cn,cp,imsh,isd)

      call sub3 (wb_n,resN,wa_n,npts)
      call sub3 (wb_p,resP,wa_p,npts)

C     dssum on RHS
      call dssum (wb_n,nx1,ny1,nz1) ! Lan
      call col2  (wb_n,dmask(1,2),npts)

      call dssum (wb_p,nx1,ny1,nz1)
      call col2  (wb_p,dmask(1,3),npts)

C     Solve Linear system
      if (ste_ac(1,2)*ste_ac(2,1).eq.0) then 
        call cggo2(wa_n,wb_n,h1cc(1,1,1),h2cc(1,1,1),mult,dmask(1,2)
     $       ,imsh,isd,npts,tol,maxit)
        call cggo2(wa_p,wb_p,h1cc(1,2,2),h2cc(1,2,2),mult,dmask(1,3)
     $       ,imsh,isd,npts,tol,maxit)
      else
        call drift_hmh_gmres_steric(wa_n,wa_p,wb_n,wb_p
     $       ,mult,isd,imsh,npts,tol)
      endif

C     Add Dirichlet back
      call add2(cn,wa_n,npts)
      call add2(cp,wa_p,npts)
    
      return
      end
      
c---------------------------------------------------------------------
      subroutine cem_drift_rhs_bdf1_steric
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'ZPER'
      include 'STERIC'
      real g11_h1(lpts),g12_h1(lpts),g21_h1(lpts),g22_h1(lpts) ! for steric
      real cdc_n(lpts), cdc_p(lpts),h1o(lpts),h2z(lpts), tmp(lpts)

      call rzero(h2z,npts)
      call rone (h1o,npts)

      call cem_drift_rhs_bdf1

C     Adding STERIC terms, Lan
      call cmult2(g11_h1,cn,ste_g(1,1),npts) ! ci * gij
      call cmult2(g12_h1,cn,ste_g(1,2),npts)
      call cmult2(g21_h1,cp,ste_g(2,1),npts)
      call cmult2(g22_h1,cp,ste_g(2,2),npts)

      call add2s2(g11_h1,h1o,-1.0*ste_ac(1,1),npts) ! c1 * g11 - ste_an
      call add2s2(g12_h1,h1o,-1.0*ste_ac(1,2),npts) ! c1 * g12 - ste_as
      call add2s2(g21_h1,h1o,-1.0*ste_ac(2,1),npts) ! c2 * g21 - ste_as
      call add2s2(g22_h1,h1o,-1.0*ste_ac(2,2),npts) ! c2 * g22 - ste_ap
      
      call col2(g11_h1,mu_n,npts) ! mu_i ci g_ij 
      call col2(g12_h1,mu_n,npts)
      call col2(g21_h1,mu_p,npts)
      call col2(g22_h1,mu_p,npts)

      call axhelm3(cdc_n,cn,g11_h1,h2z,1,1)
      call axhelm3(tmp  ,cp,g12_h1,h2z,1,1)
      call add2s2 (cdc_n,tmp,1.0,npts)

      call axhelm3(cdc_p,cn,g21_h1,h2z,1,1)
      call axhelm3(tmp  ,cp,g22_h1,h2z,1,1)
      call add2s2 (cdc_p,tmp,1.0,npts)

      call add2s2(resN,cdc_n,-1.0*dt,npts) 
      call add2s2(resP,cdc_p,-1.0*dt,npts)

      return
      end
      
C-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_steric(un,up,fn,fp
     $           ,wt,isd,imsh,n,tol)
c-----------------------------------------------------------------------
c     For Steric term, involving cross term linear system 
c      Hii = B + dt(1 + ste_aii)*A_Di
c      Hij = dt*ste_aij*A_Di
c    
c      H = (Hij), f = [f1;f2;...fp], u = [u1;u2;...;up], p = # of species
c      
c     solving H*u = f
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      include 'DRIFT'
      include 'STERIC'
C
      integer  n,outer,isd,imsh
      integer  iter,iconv
      integer  i,j,k,m
      real     divex,etime_p,tolpss
      
      real     un(n),up(n),wn(n),wp(n),rn(n),rp(n),xp(n),xn(n)
      real     fn(n),fp(n)
      real     wt(n)
      real     tol,alpha,l,temp,rnorm
      real     glsc3
      real*8   etime1,dnekclock
C

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
      call rzero(xp,n)
      call rzero(xn,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.5000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (rn,fn,n)                  ! r = res
            call copy  (rp,fp,n)                  ! r = res
         else
            !update residual
            call copy   (rn,fn,n)                  ! r = res
            call copy   (rp,fp,n)                  ! r = res
            call axhelm_steric(wn,wp,xn,xp,imsh,isd)! w = A x

            call dssum  (wn,nx1,ny1,nz1) !Lan dssum
            call col2   (wn,dmask(1,2),n)
            call dssum  (wp,nx1,ny1,nz1)
            call col2   (wp,dmask(1,3),n)

            call add2s2 (rn,wn,-1.,n)              ! r = r - w
            call add2s2 (rp,wp,-1.,n)              ! r = r - w
         endif

         gamma(1) = glsc3(rn,rn,wt,n) + glsc3(rp,rp,wt,n)  ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1     include 'SIZE'
         
         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(vn(1,1),rn,temp,n)            !  v  = r / gamma
         call cmult2(vp(1,1),rp,temp,n)            !  v  = r / gamma
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1
            call axhelm_steric(wn,wp,vn(1,j),vp(1,j),imsh,isd)! w = A v 
            call dssum  (wn,nx1,ny1,nz1)
            call col2   (wn,dmask(1,2),n)
            call dssum  (wp,nx1,ny1,nz1)
            call col2   (wp,dmask(1,3),n)

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)= glsc3(wn,vn(1,i),wt,n)      ! h    = (w,v )
     $               + glsc3(wp,vp(1,i),wt,n)      ! i,j       i

               call add2s2(wn,vn(1,i),-h(i,j),n)   ! w = w - h    v
               call add2s2(wp,vp(1,i),-h(i,j),n)   !         i,j  i
            enddo                              

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                             
            alpha = sqrt( glsc3(wn,wn,wt,n)          !            ______
     $                  + glsc3(wp,wp,wt,n) )        ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))!*norm_fac     
c           ratio = rnorm/div0
            if ((nid.eq.0).and.(istep.le.5))
     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' Divergence')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(vn(1,j+1),wn,temp,n)   ! v    = w / alpha
            call cmult2(vp(1,j+1),wp,temp,n)   ! v    = w / alpha
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
            call add2s2(xn,vn(1,i),c(i),n)     ! x = x + c  z
            call add2s2(xp,vp(1,i),c(i),n)     ! x = x + c  z
         enddo                                 !          i  i
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

       call copy(un,xn,n)
       call copy(up,xp,n)

c      call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nid.eq.0) then
          if ((mod(istep,iocomm).eq.0).or.(istep.le.10))
     $    write(6,9999) istep,iter,tolpss,etime_p,etime1
      endif
 9999 format(' ',' ',i9,' gmres   : iteration#',i5,1p3e12.4)

      return
      end
c----------------------------------------------------------------------
      subroutine axhelm_steric(an,ap,xn,xp,imsh,isd)
c----------------------------------------------------------------------
c     Costum Axhelm for steric, coupling many things
c     Use h1cc(1,i,j), h2cc(1,i,j) for  h1*A + h2*B
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'DRIFT'
      include 'STERIC'
      
      real           AN    (LX1,LY1,LZ1,LELT)
     $ ,             AP    (LX1,LY1,LZ1,LELT)
     $ ,             xN    (LX1,LY1,LZ1,LELT)
     $ ,             xP    (LX1,LY1,LZ1,LELT)
     $ ,             ttp   (LX1,LY1,LZ1,LELT) ! =tmp    
      real           term1,term2,dnekclock
      integer isd,imsh,n

      n = nx1*ny1*nz1*nelv ! npts?

      call axhelm3(an ,xn,h1cc(1,1,1),h2cc(1,1,1),imsh,isd) ! (B+A)*cn
      call axhelm3(ttp,xp,h1cc(1,1,2),h2cc(1,1,2),imsh,isd) ! A*cp

      call add2 (an,ttp,n)

      call axhelm3(ttp,xn,h1cc(1,2,1),h2cc(1,2,1),imsh,isd) ! A*cn 
      call axhelm3(ap ,xp,h1cc(1,2,2),h2cc(1,2,2),imsh,isd) ! (B+A)*cp

      call add2 (ap,ttp,n)

      return
      end

c----------------------------------------------------------------------
      subroutine bcneusc_steric(wa,wb,n,ifld,itype)
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'
      include 'STERIC'
      integer n,ifld,itype ! ifld = 1,2
      real wa(n),wb(n) ! wb=mu_n
      real tmp
      integer i,j,k,l,ifld_c
      real tmp_c(n,2)!nion

      call copy(tmp_c(1,1),cN,n) !todo wait for restructure, Lan
      call copy(tmp_c(1,2),cP,n)
      
      ifld_c = ifld - 1 ! local index of cn cp
      
      if (ifld.eq.2) then !cN

        if (ncemface_rob(ifld).gt.0) then
          do i = 1,ncemface_rob(ifld)
            j = cemface_rob(i,ifld)
            k = cemface(j)
c            do l=ifld_c,ifld_c  ! run through cN cP
              tmp = gradc(k,1,1)*unxm(j)
     $            + gradc(k,2,1)*unym(j)+gradc(k,3,1)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k)*ste_g(1,1)-ste_ac(1,1))
     $              *tmp*aream(j)
              tmp = gradc(k,1,2)*unxm(j)
     $            + gradc(k,2,2)*unym(j)+gradc(k,3,2)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k)*ste_g(1,2))
     $              *tmp*aream(j)
c            enddo
          enddo
        endif

        if (ncemface_nmn(ifld).gt.0) then
          do i = 1,ncemface_nmn(ifld)
            j = cemface_nmn(i,ifld)
            k = cemface(j)
c            do l=ifld_c,ifld_c  ! run through cN cP
              tmp = gradc(k,1,1)*unxm(j)
     $            + gradc(k,2,1)*unym(j)+gradc(k,3,1)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k)*ste_g(1,1)-ste_ac(1,1))
     $              *tmp*aream(j)
              tmp = gradc(k,1,2)*unxm(j)
     $            + gradc(k,2,2)*unym(j)+gradc(k,3,2)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k)*ste_g(1,2))
     $              *tmp*aream(j)
c            enddo
          enddo
        endif

      elseif (ifld.eq.3) then

        if (ncemface_rob(ifld).gt.0) then
          do i = 1,ncemface_rob(ifld)
            j = cemface_rob(i,ifld)
            k = cemface(j)
c            do l=ifld_c,ifld_c  ! run through cN cP
              tmp = gradc(k,1,2)*unxm(j)
     $            + gradc(k,2,2)*unym(j)+gradc(k,3,2)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cP(k)*ste_g(2,2)-ste_ac(2,2))
     $              *tmp*aream(j)
              tmp = gradc(k,1,1)*unxm(j)
     $            + gradc(k,2,1)*unym(j)+gradc(k,3,1)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cP(k)*ste_g(2,1))
     $              *tmp*aream(j)
c            enddo
          enddo
        endif

        if (ncemface_nmn(ifld).gt.0) then
          do i = 1,ncemface_nmn(ifld)
            j = cemface_nmn(i,ifld)
            k = cemface(j)
c            do l=ifld_c,ifld_c  ! run through cN cP
              tmp = gradc(k,1,2)*unxm(j)
     $            + gradc(k,2,2)*unym(j)+gradc(k,3,2)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cP(k)*ste_g(2,2)-ste_ac(2,2))
     $              *tmp*aream(j)
              tmp = gradc(k,1,1)*unxm(j)
     $            + gradc(k,2,1)*unym(j)+gradc(k,3,1)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cP(k)*ste_g(2,1))
     $              *tmp*aream(j)
c            enddo
          enddo
        endif

      endif

      return
      end

c----------------------------------------------------------------------
