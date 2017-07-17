c-----------------------------------------------------------------------
c     Inside cem_drift_newton.F we use Jacobi-free Newton Krylov method
c     (JFNK) to solve the steady state of PNP equation
c     By Ping Hsuan 06/27/2017
c-----------------------------------------------------------------------
c     This subroutine aim to initial variables for JFNK
c     Those variables are defined in NEWTON
      subroutine cem_drift_newton_init
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'      
      include 'POISSON'
      include 'NEWTON'

      call rzero(sp_k  ,npts)
      call rzero(sn_k  ,npts)
      call rzero(gp_k  ,npts)
      call rzero(gn_k  ,npts)
      call rzero(cpi   ,npts)
      call rzero(cni   ,npts)
      call rzero(cpo   ,npts)
      call rzero(cno   ,npts)
      call rzero(cp_k  ,npts)
      call rzero(cn_k  ,npts)
      call rzero(fp    ,npts)
      call rzero(fn    ,npts)

      return
      end
c---------------------------------------------------------------------
c     This subroutine is the main core of JFNK method which involves
c     using variables cN and cP in DRIFT
      subroutine cem_drift_newton
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'      
      include 'POISSON'
      include 'NEWTON'
      include 'CTIMER'
      include 'RTIMER'

      real alpha,tt
      real new_tol,tol
      real glmax,glsc3 
      real rnorm,rnorm1,rnorm2,rinorm1,rinorm2,rinorm,ratio
      real fnorm,fnorm1,fnorm2

      integer i,j,k,n
      integer maxit
      real glsum
      real f_pre,f_now

      alpha = 1            ! relaxation parameter \alpha \in (0,1]             
c     nsteps = 200         ! steps for pseudo-transient continuation
      maxit = 20           ! iteration for Newton method
      jaceps = 1e-5        ! perturbation parameter for Jacobi-free
      epsinv = 1./jaceps   ! eps for approximation for jacobian
      tol = 1e-2           ! tolerance for GMRES
      new_tol = 1e-5       ! tolerance for Newton method

      dt_newton = param(1) ! pseudo-transient time step
                           ! it can grow as istep increase 
      dtinv = 1./dt
      dt_newton_inv = 1./dt_newton

      cpu_t = 0.0
      cpu_chk = 0.0

      call cem_drift_newton_init

      call copy (cp_k,cP,npts) ! get the initial condition 
      call copy (cn_k,cN,npts) ! get the initial condition

      call copy(cpi,cp_k,npts)
      call copy(cni,cn_k,npts)

      call compute_gk      ! compute rhs g_k for GMRES

      fnorm1 = glsc3(fp,fp,mult,npts)  
      fnorm2 = glsc3(fn,fn,mult,npts)
      fnorm1 = sqrt(fnorm1)
      fnorm2 = sqrt(fnorm2)
      fnorm  = max(fnorm1,fnorm2)

      f_now = fnorm

      call chsign (gp_k,npts)
      call chsign (gn_k,npts)

      call compute_energy  ! check the energy of the system in the initial
                           ! stage

c     start pseudo-time step
      do istep = 1,nsteps 

         cpu_dtime = dclock()

         call rzero(sn_k,npts)
         call rzero(sp_k,npts)

c        SER for tau^{istep}
         if ((istep .eq. 1)) then
            dt_newton = param(1)
c            dt = dt_newton
         else 
            dt_newton = dt_newton * f_pre/f_now
c            dt = dt_newton
         end if

         dtinv = 1./dt
         dt_newton_inv = 1./dt_newton
         
c        start Newton iteration

         do newton_iter = 1, maxit

            call drift_hmh_gmres_newton   ! solve sn_k
     $      (sn_k,gn_k,cn_k,cp_k,fn,mult,npts,tol,2)
            
            call drift_hmh_gmres_newton   ! solve sp_k
     $      (sp_k,gp_k,cp_k,cn_k,fp,mult,npts,tol,1)
            
            call add2s2(cp_k,sp_k,alpha,npts)   ! cp_k = cp_k + alpha * sp_k
            call add2s2(cn_k,sn_k,alpha,npts)   ! cn_k = cn_k + alpha * sn_k

            call compute_gk   ! compute rhs g_k for checking tolerance
                              ! and also store for next Newton iteration
            call chsign (gp_k,npts)
            call chsign (gn_k,npts)

            if (newton_iter.eq.1) then
               rinorm1 = glsc3(gp_k,gp_k,mult,npts)
               rinorm2 = glsc3(gn_k,gn_k,mult,npts)
               rinorm1 = sqrt(rinorm1)
               rinorm2 = sqrt(rinorm2)
               rinorm  = max(rinorm1,rinorm2)
            else 
               rnorm1 = glsc3(gp_k,gp_k,mult,npts)
               rnorm2 = glsc3(gn_k,gn_k,mult,npts)
               rnorm1 = sqrt(rnorm1)
               rnorm2 = sqrt(rnorm2)
               rnorm  = max(rnorm1,rnorm2)
            endif
            ratio= rnorm/rinorm

            if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) then
               write(6,90) istep,newton_iter,ratio,rnorm1,rnorm2,
     $                     dt_newton,dt
            endif 

           if (ratio.lt.new_tol) goto 900
c          if (rnorm.lt.new_tol) goto 900

         enddo
 90   format('newton iter',2i6,1p5e12.4)   
900   continue
      
      call copy(cpi,cp_k,npts)   ! update cp_i for compute gp 
      call copy(cni,cn_k,npts)   ! update cn_i for compute gn
                                 ! for next pseudo-time step
      call copy(cP,cp_k,npts) ! just for cem_out
      call copy(cN,cn_k,npts) ! just for cem_out

c     compute the norm of f(c_p)^{istep-1} and f(c_n)^{istep-1}
      fnorm1 = glsc3(fp,fp,mult,npts)  
      fnorm2 = glsc3(fn,fn,mult,npts)
      fnorm1 = sqrt(fnorm1)
      fnorm2 = sqrt(fnorm2)
      fnorm  = max(fnorm1,fnorm2)
      f_pre = f_now  ! set up old fnorm to f_pre
      f_now = fnorm  ! set up new fnorm to f_now

c     compute the CPU_time
      cpu_dtime = dclock()-cpu_dtime
      cpu_t = cpu_t+cpu_dtime
      cpu_t_step = cpu_t/istep
      cpu_p_t = glsum(cpu_t_step /npts,1)/np
      
      time = time + dt
  
      call compute_energy

      call userchk

      call cem_out
      
      enddo
      
      call cem_end

      return

      end
c-----------------------------------------------------------------------
c     This routine computes gk for each Newton iteration
c
c          g^{n}_k = u^{n}_k - \delta t f(u^{n}_k) - u^{n}_0
c
c     Input    cpi,cni denote the initial of each Newton iteration 
c     Input    cp_k,cn_k
c     Output   gp_k,gn_k
      subroutine compute_gk
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEWTON'

      call rzero(fp,npts)
      call rzero(fn,npts)

      call compute_f ! compute f(u^{n}_k)
     
      call rzero(gp_k,npts)
      call rzero(gn_k,npts)

      call add2s2(gp_k,cpi,-1.0,npts)
      call add2s2(gn_k,cni,-1.0,npts)

      call add2s2(gp_k,fp,-1.0*dt_newton,npts)
      call add2s2(gn_k,fn,-1.0*dt_newton,npts)

      call add2s2(gp_k,cp_k,1.0,npts)
      call add2s2(gn_k,cn_k,1.0,npts)
      
      call cmult(gp_k,dt_newton_inv,npts)   ! case with divided by dt_newton
      call cmult(gn_k,dt_newton_inv,npts)   ! case with divided by dt_newton

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON' 
      call rzero(fp,npts)
      call rzero(fn,npts)

      call copy(cP,cp_k,npts) 
      call copy(cN,cn_k,npts)

      call cem_drift_source
      call cem_drift_poisson
      call cem_drift_sem_bdf1_newton

      call copy(cpo,cP,npts)
      call copy(cno,cN,npts)

      call sub3(fp,cpo,cp_k,npts)
      call sub3(fn,cno,cn_k,npts)

      call cmult(fp,dtinv,npts)
      call cmult(fn,dtinv,npts)

      return
      end
c-----------------------------------------------------------------------
c     This routine computes Jp where J is the Jacobian matrix and p 
c     is a vector. Note that we do not construct the Jacobian matrix 
c     exactly but use the following formula to compute Jp
c   
c        J_k s_k = s_k - (dt /eps) *  ( f(u_k + eps*s_k) - f(u_k) )
c   
c     where f(u_k) has been store in each Newton iteration
c           f(u_k + eps*s_k) has to be computed in each GMRES iteration
      subroutine Jacobi_free(Jp,u0,u1,f,p,iflag)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(npts), p(npts)
      real     u0(npts), u1(npts)
      real     tmp(npts), tmp1(npts)
      real     fo(npts), f(npts)
      integer  iflag
      real     eps
      real     pnorm,unorm,glsc3
     
      
      pnorm = glsc3(p,p,mult,npts)
      unorm = glsc3(u0,u0,mult,npts)
      pnorm = sqrt(pnorm)
      unorm = sqrt(unorm)

      eps = (1+unorm)*1e-14   ! formula for varing eps
      eps = sqrt(eps)         ! for using the case using 
                              ! "ub" filled with previous
                              !  step info first and 
                              !  Dirichlet info on the boundary
      eps = eps/pnorm         ! this eps will blow up because of 
                              ! pnorm = 0

      epsinv = 1./eps         

      if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) then
          write(6,10) istep,newton_iter,eps,epsinv
      endif

  10  format (2i8,' eps=',1p2e17.7)
      call rzero  (Jp,npts)
      call copy   (tmp,u0,npts)
      call add2s2 (tmp,p,eps,npts) !jaceps -> eps
      
      if (iflag .eq. 1) then
      
      call copy   (cP,tmp,npts)
      call copy   (cN,u1 ,npts)
      
      else if (iflag .eq. 2) then
      
      call copy   (cP,u1,npts)
      call copy   (cN,tmp,npts)

      endif
      
      call cem_drift_source
      call cem_drift_poisson
      call cem_drift_sem_bdf1_newton

      if (iflag .eq. 1) then

      call sub3   (fo,cP,tmp,npts)        ! fo = cP - tmp
      call cmult  (fo,dtinv,npts)         ! fo = fo / dt
      call sub3   (tmp1,fo,f,npts)        ! tmp1 = fo - f
      call cmult  (tmp1,epsinv,npts)      ! tmp1 = (fo-f)/eps
      call cmult  (tmp1,dt_newton,npts)   ! tmp1 = (fo-f)*dt_newton/eps
      call sub3   (Jp,p,tmp1,npts)        ! Jp = p - tmp1
      call cmult  (Jp,dt_newton_inv,npts) ! Jp = p/dt_newton - (fo-f)/eps     
                                          ! case with divided by dt_newton

      else if (iflag .eq. 2) then

      call sub3   (fo,cN,tmp,npts)        ! fo = cN - tmp
      call cmult  (fo,dtinv,npts)         ! fo = fo / dt
      call sub3   (tmp1,fo,f,npts)        ! tmp1 = fo - f
      call cmult  (tmp1,epsinv,npts)      ! tmp1 = (fo-f)/eps
      call cmult  (tmp1,dt_newton,npts)   ! tmp1 = (fo-f)*dt_newton/eps
      call sub3   (Jp,p,tmp1,npts)        ! Jp = p - tmp1
      call cmult  (Jp,dt_newton_inv,npts) ! Jp = p/dt_newton - (fo-f)/eps
                                          ! case with divided by dt_newton
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton
     $           (phi,res,u0,u1,f,wt,n,tol,iflag)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      include 'NEWTON'

      integer  n,outer, iflag
      real     phi(n),res(n),wt(n)
      real     tol,alpha,l,temp
      real     eps, u0(n), u1(n), f(n)
      real*8   etime1,dnekclock

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
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
            call Jacobi_free (w,u0,u1,f,x,iflag)! w = A x
            call add2s2 (r,w,-1.,n) ! r = r - w
         endif

         gamma(1) = glsc3(r,r,wt,n)                ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }

         tolps = 0.1*gamma(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call Jacobi_free (w,u0,u1,f,v(1,j),iflag) ! w = A v

            !modified Gram-Schmidt
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
            alpha = sqrt(glsc3(w,w,wt,n))     ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))

c           if ((nid.eq.0).and.(istep.le.2))
c    $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

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
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,newton_iter,iter,tolpss
      endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_newton
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'

      integer  i,j,i0,isd,imsh,maxit
      real     wk1(lpts),wk2(lpts), tol,glamax
      real     diff_n_dt(lpts), diff_p_dt(lpts),diff_e_dt(lpts)
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/

      real wa (lx1*ly1*lz1*lelt)
      real wb (lx1*ly1*lz1*lelt)

      npts= nx1*ny1*nz1*nelt 
                             
      if (icalled.eq.0) then

        call cmult2 (h1n,diff_n,dt,npts) ! h1n= h1n * (diff_n*dt)
        call rone (h2n,npts)             ! h2n= 1.0

        call cmult2 (h1p,diff_p,dt,npts) ! h1p= h1p * (diff_p*dt)
        call rone (h2p,npts)             ! h2p= 1.0

        icalled = 1

       endif

c     set cggo parameters

      isd = 1
      imsh = 1
      maxit = 5000
      tol = param(22)

c     get rhs for cN and cP

      call cem_drift_rhs_bdf1
      call bcneusc0 (resN,mu_n,cn,npts,2,1)
      call bcneusc0 (resP,mu_p,cp,npts,3,1)

      if (param(2).ne.0) then ! this is trial case

c     solve cn
      call bcneusc (wa,npts,2,-1)
      call col2 (wa,h1n,npts)
      call add2 (h2n,wa,npts)
      call bcdirsc (cn,npts,2)
      call axhelm3 (wa,cn,h1n,h2n,1,1)
      call sub3 (wb,resN,wa,npts)
      call bcneusc (wa,npts,2,1)
      call col2 (wa,h1n,npts)
      call add2 (wb,wa,npts)
 
      call dssum (wb,nx1,ny1,nz1)
      call col2 (wb,dmask(1,2),npts)
      call cggo2(wa,wb,h1n,h2n,mult,dmask(1,2),imsh,isd,npts,tol,maxit)
      call add2 (cn,wa,npts)

c     solve cp: fixme for newton

      call bcneusc (wa,npts,3,-1)
      call col2 (wa,h1p,npts)
      call add2 (h2p,wa,npts)
      call bcdirsc (cp,npts,3)
      call axhelm3 (wa,cp,h1p,h2p,1,1)
      call sub3 (wb,resP,wa,npts)
      call bcneusc (wa,npts,3,1)
      call col2 (wa,h1p,npts)
      call add2 (wb,wa,npts)

      call dssum (wb,nx1,ny1,nz1)
      call col2 (wb,dmask(1,3),npts)
      call cggo2(wa,wb,h1p,h2p,mult,dmask(1,3),imsh,isd,npts,tol,maxit)
      call add2(cp,wa,npts)

      elseif (param(2).eq.0) then ! this is default option

c     solve cn: below support dirichlet inhomogeneous

      call axhelm3 (wa,cn_b,h1n,h2n,1,1)
      call sub3    (wb,resN,wa,npts)

      call dssum (wb,nx1,ny1,nz1)
      call col2 (wb,dmask(1,2),npts)
      call cggo2 (wa,wb,h1n,h2n,mult,dmask(1,2),imsh,isd,npts,tol,maxit)
      call add3 (cn,cn_b,wa,npts)

c     solve cp: below support dirichlet inhomogeneous

      call axhelm3 (wa,cp_b,h1p,h2p,1,1)
      call sub3 (wb,resP,wa,npts)

      call dssum (wb,nx1,ny1,nz1)
      call col2 (wb,dmask(1,3),npts)
      call cggo2 (wa,wb,h1p,h2p,mult,dmask(1,3),imsh,isd,npts,tol,maxit)
      call add3 (cp,cp_b,wa,npts)
     
      endif

      return
      end
c---------------------------------------------------------------------
