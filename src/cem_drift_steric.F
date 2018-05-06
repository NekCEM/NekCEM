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

      integer  isd,imsh,maxit
      integer  i,ic,ic2
      real     wk1(lpts),wk2(lpts), tol,glamax
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/
      real wa (lx1*ly1*lz1*lelt)
      real wb (lx1*ly1*lz1*lelt)
      
      real wa_C (lx1*ly1*lz1*lelt,LDIMT-1)
      real wb_C (lx1*ly1*lz1*lelt,LDIMT-1)
      real tmp  (lx1*ly1*lz1*lelt)

      npts= lx1*ly1*lz1*nelt 
                             
      if (icalled.eq.0) then
        do ic = 1,lcdim
          call cmult2 (h1c(1,ic),diff_n(1,ic),dt,npts) ! h1n= h1n * (diff_n*dt)
          call rone (h2c(1,ic),npts)             ! h2n= 1.0
        enddo

        if (ifstevec) then

          do ic = 1,lcdim
            do ic2 = 1,lcdim  
              if (ic.eq.ic2) then
                call cadd2(tmp,s_veca(1,ic,ic2),1.0,npts) ! tmp = s_veca + 1.0
                call col3 (h1s(1,ic,ic2),h1c(1,ic),tmp,npts)
                call copy (h2s(1,ic,ic2),h2c(1,ic),npts)
              else
                call col3 (h1s(1,ic,ic2),h1c(1,ic)
     $               ,s_veca(1,ic,ic2),npts)
                call rzero(h2s(1,ic,ic2),npts)
              endif
            enddo
          enddo
c          icalled = 1! If not Dirichlet, need uptate each loop 

        else

          do ic = 1,lcdim
            do ic2 = 1,lcdim  
              if (ic.eq.ic2) then
                call cmult2 (h1s(1,ic,ic2),h1c(1,ic)
     $               ,1.0+s_a(ic,ic2),npts)
                call copy   (h2s(1,ic,ic2),h2c(1,ic),npts)
              else
                call cmult2 (h1s(1,ic,ic2),h1c(1,ic),s_a(ic,ic2),npts)
                call rzero  (h2s(1,ic,ic2),npts)
              endif
            enddo
          enddo
c          icalled = 1! If not Dirichlet, need uptate each loop 
        endif

      endif

c...  set cggo parameters

      isd = 1
      imsh = 1
      maxit = 5000
      tol = param(22)

c...  get rhs for cN and cP

      call cem_drift_rhs_bdf1_steric ! note: output change to resC and ionC

C     boundary surface integral
      do ic=1,lcdim
         call bcneusc0(rescn(1,ic),mu_n(1,ic),cN(1,ic),npts,ic+1,1)
         call cem_drift_grad(gradc(1,1,ic),gradc(1,2,ic),gradc(1,3,ic)
     $        ,cN(1,ic))
      enddo

      if (ifstevec) then
        do ic=1,lcdim
          call bcneusc_stevec(rescn(1,ic),mu_n(1,ic),lcdim,npts,ic+1,1)
        enddo
      else
        do ic=1,lcdim ! involves coupled gradc
          call bcneusc_steric(rescn(1,ic),mu_n(1,ic),lcdim,npts,ic+1,1)
        enddo
      endif

C     Robin, add to LHS
      do ic=1,lcdim
         call bcneusc(wa_C(1,ic),npts,ic+1,-1)
         call col2 (wa_C(1,ic),h1s(1,ic,ic),npts)
         call add2 (h2s(1,ic,ic),wa_C(1,ic),npts)
      enddo

C     Neumann + Robin, add to RHS
      do ic=1,lcdim
         call bcneusc(wa_C(1,ic),npts,ic+1,1)
         call col2 (wa_C(1,ic),h1s(1,ic,ic),npts)
         call add2 (rescn(1,ic),wa_C(1,ic),npts)
      enddo

C     Dirichlet seperate
      do ic=1,lcdim
         call bcdirsc (cN(1,ic),npts,ic+1)
      enddo
      call axhelm_steric(wa_C,cN,h1s,h2s,lcdim,imsh,isd)
      do ic=1,lcdim
         call sub3(wb_C(1,ic),rescn(1,ic),wa_C(1,ic),npts)
      enddo

C     dssum on RHS
      do ic=1,lcdim
         call dssum (wb_C(1,ic),nx1,ny1,nz1) ! Lan
         call col2  (wb_C(1,ic),dmask(1,ic+1),npts)
      enddo

C     Solve Linear system
      if (ifstevec) then
        call drift_hmh_gmres_steric(wa_C,wb_C,h1s,h2s,lcdim
     $       ,mult,dmask(1,2),isd,imsh,npts,tol)
      else
        if (s_a(1,2)*s_a(2,1).eq.0) then ! Todo: need a flag to tell if diag.
           do ic=1,lcdim
              call cggo2(wa_C(1,ic),wb_C(1,ic),h1s(1,ic,ic),h2s(1,ic,ic)
     $             ,mult,dmask(1,ic+1),imsh,isd,npts,tol,maxit)
           enddo
        else
          call drift_hmh_gmres_steric(wa_C,wb_C,h1s,h2s,lcdim
     $         ,mult,dmask(1,2),isd,imsh,npts,tol)
        endif
      endif

C     Add Dirichlet back
      do ic=1,lcdim
         call add2(cn(1,ic),wa_C(1,ic),npts)
      enddo

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
      real h1o(lpts),h2z(lpts),tmp(lpts)
      real g_h1(lpts,LDIMT-1,LDIMT-1) 
      integer i,j,ic

      call rone (h1o,npts)
      call rzero(h2z,npts)

      do ic = 1,lcdim
        call cem_drift_rhs_bdf1(ic)
      enddo

C     Adding STERIC terms, Lan

      if (ifstevec) then

        do i=1,lcdim
          do j=1,lcdim
            ! g_h1(i,j) = c_i * ste_g(i,j) - ste_sc(i,j)
            call col3   (g_h1(1,i,j),cN(1,i),s_vecg(1,i,j),npts)
            call subcol3(g_h1(1,i,j),h1o,s_veca(1,i,j),npts)
            call col2   (g_h1(1,i,j),mu_n(1,i),npts)
          enddo
        enddo

      else
        do i=1,lcdim
          do j=1,lcdim
            ! g_h1(i,j) = c_i * ste_g(i,j) - ste_sc(i,j)
            call cmult2(g_h1(1,i,j),cN(1,i),s_g(i,j),npts)
            call add2s2(g_h1(1,i,j),h1o,-1.0*s_a(i,j),npts)
            call col2  (g_h1(1,i,j),mu_n(1,i),npts)
          enddo
        enddo
      endif

      ! res_i = res_i - dt * mu_i sum_j{ [c_i ste_g(i,j) - ste_ac(i,j)] * dC_j }
      do i=1,lcdim
         do j=1,lcdim
            call rzero(tmp,npts)
            call axhelm3(tmp,cN(1,j),g_h1(1,i,j),h2z,1,1)
            call add2s2 (rescn(1,i),tmp,-1.0*dt,npts)
         enddo
      enddo

      return
      end
C-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_steric(uC,fC,helm1,helm2,nion
     $           ,wt,mask,isd,imsh,n,tol)
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
      include 'DRIFT'
      include 'MGMRES'
C
      integer  nion,n,outer,isd,imsh
      integer  iter,iconv
      integer  i,j,k,m
      integer  ic ! safe index fo nion
      real     divex,etime_p,tolpss
      
      real     uC(lpts,LDIMT-1)
     $       , fC(lpts,LDIMT-1)
     $       , helm1(lpts,LDIMT-1,LDIMT-1) !h1cc
     $       , helm2(lpts,LDIMT-1,LDIMT-1) !h2cc
     $       , mask(lpts,LDIMT-1)
     $       , wt(lpts)
      real     tol,alpha,l,temp,rnorm
      real     glsc3
      real*8   etime1,dnekclock
C

      if (nid.eq.0) write(6,*) 'start: hmh_gmres_steric2'

      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      call rzero(xxC,lpts*(LDIMT-1))
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.5000)
         outer = outer+1
         if(iter.eq.0) then
            do ic=1,nion
            call copy  (rC(1,ic),fC(1,ic),n)      ! r = res
            enddo
         else
            !update residual
            do ic=1,nion
            call copy  (rC(1,ic),fC(1,ic),n)      ! r = res
            enddo
            call axhelm_steric(wC,xxC,helm1,helm2,nion,imsh,isd)! w = A x

            do ic=1,nion
            call dssum  (wC(1,ic),nx1,ny1,nz1)
            call col2   (wC(1,ic),mask(1,ic),n)
            enddo
            do ic=1,nion
            call add2s2 (rC(1,ic),wC(1,ic),-1.,n)   ! r = r - w
            enddo
         endif
         gamma(1) = 0.0
         do ic=1,nion
           gamma(1) = gamma(1)+glsc3(rC(1,ic),rC(1,ic),wt,n)! gamma  = (r,r)
         enddo
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         write(6,*) 'initial residule',gamma(1)
                                                   !      1     include 'SIZE'
         
         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         do ic=1,nion
         call cmult2(vC(1,ic,1),rC(1,ic),temp,n)            !  v  = r / gamma
         enddo
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1
            call axhelm_steric(wC,vC(1,1,j),helm1,helm2,nion,imsh,isd)! w = A v 
            do ic=1,nion
            call dssum  (wC(1,ic),nx1,ny1,nz1)
            call col2   (wC(1,ic),mask(1,ic),n)
            enddo

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)=0.0
               do ic=1,nion
               h(i,j) = h(i,j) + glsc3(wC(1,ic),vC(1,ic,i),wt,n) ! h    = (w,v )
               enddo                                             ! i,j       i
               do ic=1,nion
               call add2s2(wC(1,ic),vC(1,ic,i),-h(i,j),n)   ! w = w - h    v
               enddo                                        !         i,j  i
            enddo                              

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
            alpha=0.0
            do ic=1,nion                            
            alpha = alpha + glsc3(wC(1,ic),wC(1,ic),wt,n)  !            ______
            enddo                                          ! alpha =  \/ (w,w)
            alpha=sqrt(alpha)

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
            do ic=1,nion
            call cmult2(vC(1,ic,j+1),wC(1,ic),temp,n)   ! v    = w / alpha
            enddo                                       !  j+1
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
            do ic=1,nion
            call add2s2(xxC(1,ic),vC(1,ic,i),c(i),n)  ! x = x + c  z
            enddo                                      !  i  i
         enddo
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue
       
       do ic=1,nion
       call copy(uC(1,ic),xxC(1,ic),n)
       enddo

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
      subroutine axhelm_steric(AU,U,helm1,helm2,nion,imsh,isd)
c----------------------------------------------------------------------
c     Costum Axhelm for steric, coupling with nion
c     Use h1cc(1,i,j), h2cc(1,i,j) for  h1*A + h2*B
      implicit none
      include 'SIZE'
      include 'STERIC' 
      include 'DRIFT'
      integer nion,n      
      real           AU    (lpts,LDIMT-1)
     $ ,             U     (lpts,LDIMT-1)
     $ ,             tmpAU (lpts) ! =tmp    
     $ ,             helm1 (lpts,LDIMT-1,LDIMT-1) ! h1cc    
     $ ,             helm2 (lpts,LDIMT-1,LDIMT-1) ! h2cc   
      integer isd,imsh
      integer i,j
      n = nx1*ny1*nz1*nelv

      call rzero(AU,lpts*nion)      

C     Todo: if diagonal, switch to this loop to speed up
c      do i=1,nion
c         call axhelm3(tmpAU,U(1,i),helm1(1,i,i),helm2(1,i,i),imsh,isd)
c         call add2   (AU(1,i),tmpAU,n)
c      enddo
      do i=1,nion
        do j=1,nion
          call axhelm3(tmpAU,U(1,j),helm1(1,i,j),helm2(1,i,j),imsh,isd)
          call add2   (AU(1,i),tmpAU,n)
        enddo
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine bcneusc_steric(wa,wb,nion,n,ifld,itype)
c----------------------------------------------------------------------
C     This subroutine deal with the surface integral of STERIC term.
C     Since everything is coupling, to avoid too complicated, 
C     we only do single species and put a loop outside the routine
C     input ionC, gradC
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'
      include 'STERIC'
      integer nion,n,ifld,itype 
      real wa(lpts),wb(lpts) ! wb=mu_n
      real tmp
      integer i,j,k,l,ifld_c
      
      ifld_c = ifld - 1 ! local ion index

      if (ncemface_rob(ifld).gt.0) then
        do i = 1,ncemface_rob(ifld)
          j = cemface_rob(i,ifld)
          k = cemface(j)
          do l=1,nion
            if (l.eq.ifld_c) then ! if diag
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k,ifld_c)*s_g(ifld_c,l)-s_a(ifld_c,l))
     $              *tmp*aream(j)
            else
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k,ifld_c)*s_g(ifld_c,l))
     $              *tmp*aream(j)
            endif
          enddo
        enddo
      endif

      if (ncemface_nmn(ifld).gt.0) then
        do i = 1,ncemface_nmn(ifld)
          j = cemface_nmn(i,ifld)
          k = cemface(j)
          do l=1,nion
            if (l.eq.ifld_c) then ! if diag
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k,ifld_c)*s_g(ifld_c,l)-s_a(ifld_c,l))
     $              *tmp*aream(j)
            else
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k,ifld_c)*s_g(ifld_c,l))
     $              *tmp*aream(j)
            endif
          enddo
        enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine bcneusc_stevec(wa,wb,nion,n,ifld,itype)
c----------------------------------------------------------------------
C     This subroutine deal with the surface integral of STERIC term.
C     Using vector s_vecg instead of s_g
C     Since everything is coupling, to avoid too complicated, 
C     we only do single species and put a loop outside the routine
C     input ionC, gradC
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'
      include 'STERIC'
      integer nion,n,ifld,itype 
      real wa(lpts),wb(lpts) ! wb=mu_n
      real tmp
      integer i,j,k,l,ifld_c
      
      ifld_c = ifld - 1 ! local ion index

      if (ncemface_rob(ifld).gt.0) then
        do i = 1,ncemface_rob(ifld)
          j = cemface_rob(i,ifld)
          k = cemface(j)
          do l=1,nion
            if (l.eq.ifld_c) then ! if diag
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *( cN(k,ifld_c)*s_vecg(k,ifld_c,l)
     $              -  s_veca(k,ifld_c,l) )
     $              *tmp*aream(j)
            else
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k,ifld_c)*s_vecg(k,ifld_c,l))
     $              *tmp*aream(j)
            endif
          enddo
        enddo
      endif

      if (ncemface_nmn(ifld).gt.0) then
        do i = 1,ncemface_nmn(ifld)
          j = cemface_nmn(i,ifld)
          k = cemface(j)
          do l=1,nion
            if (l.eq.ifld_c) then ! if diag
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *( cN(k,ifld_c)*s_vecg(k,ifld_c,l)
     $              -  s_veca(k,ifld_c,l) )
     $              *tmp*aream(j)
            else
              tmp = gradc(k,1,l)*unxm(j)
     $            + gradc(k,2,l)*unym(j)+gradc(k,3,l)*unzm(j)
              wa(k) = wa(k) + dt*wb(k)
     $              *(cN(k,ifld_c)*s_vecg(k,ifld_c,l))
     $              *tmp*aream(j)
            endif
          enddo
        enddo
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine chk_diff(s3,ary1,ary2) ! Todo: tempory, wait to remove
c----------------------------------------------------------------------
C     find l2 difference, used in debug, Lan
      include 'SIZE'
      character*3 s3
      parameter (lt=lx1*ly1*lz1*lelt)
      real ary1(lt),ary2(lt)
      real err, aerr(lt) ! working array
      real glsc2

      n=nx1*ny1*nz1*nelt

      call sub3(aerr,ary1,ary2,n)            
      err = glsc2(aerr,aerr,n)

      if (nid.eq.0) write(6,*)'checkl2err: ',s3,err

      return
      end
c----------------------------------------------------------------------
