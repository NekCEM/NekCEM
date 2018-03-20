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

      integer ic

      do ic = 1,lcdim
        call rzero(sn_k(1,ic),npts)
        call rzero(gn_k(1,ic),npts)
        call rzero(cni (1,ic),npts)
        call rzero(cno (1,ic),npts)
        call rzero(cn_k(1,ic),npts)
        call rzero(fn  (1,ic),npts)
      enddo

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
      real tolNT,tolGMRES
      real glmax,glsc3 
      real rnorm,rnorms(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)
      real temp(lpts1)

      integer i,j,k,n,ic,isep
      integer maxit
      real glsum
      real f_pre,f_now
      real max_dtNT,min_tol

      isep = 1 ! 0=couple, 1=separate, Lan

      alpha = 1            ! relaxation parameter \alpha \in (0,1]             
c     nsteps = 200         ! steps for pseudo-transient continuation
      maxit = 20           ! iteration for Newton method
      jaceps = 1e-5        ! perturbation parameter for Jacobi-free
      epsinv = 1./jaceps   ! eps for approximation for jacobian
      tolGMRES = 1e-10           ! tolerance for GMRES
      tolNT = 1e-5       ! tolerance for Newton method

      dtNT = param(1) ! pseudo-transient time step
                           ! it can grow as istep increase 
      max_dtNT  = 1E5   ! Max. of dtNT, avoid NaN
      min_tol   = 1E-12 ! Min. of tol,  avoid NaN

      dtinv = 1./dt
      dtNTinv = 1./dtNT

      cpu_t = 0.0
      cpu_chk = 0.0

      call cem_drift_newton_init
      
      do ic = 1,lcdim
        call copy (cn_k(1,ic),cN(1,ic),npts) ! get the initial condition
        call copy (cni(1,ic),cn_k(1,ic),npts)
      enddo

      call compute_f_couple(fn,cn_k,npts)
      call compute_gk(gn_k,fn,cn_k,cni,npts) ! compute rhs g_k for GMRES
      do ic = 1,lcdim
        call chsign (gn_k(1,ic),npts)
      enddo

      do ic = 1,lcdim
        fnorms(ic) = glsc3(fn(1,ic),fn(1,ic),mult,npts)
        fnorms(ic) = sqrt(fnorms(ic))
      enddo
      fnorm = glmax(fnorms,lcdim)
      f_now = fnorm

C      call compute_energy  ! check the energy of the system in the initial
                           ! stage
c     start pseudo-time step
      do istep = 1,nsteps 

        cpu_dtime = dclock()

c       SER for tau^{istep}
        if ((istep .eq. 1)) then
          dtNT = param(1)
c          dt = dt_newton
        else
          if (f_now.lt.min_tol) then
            dtNT = max_dtNT
          else 
            dtNT = dtNT * f_pre/f_now
c            dt = dt_newton
          endif
        endif
        if (dtNT.gt.max_dtNT) then
          dtNT = max_dtNT           
        endif
        dtinv = 1./dt
        dtNTinv = 1./dtNT
         
c       start Newton iteration

        do iterNT=1,maxit
          if (ifsteric) then
c            call drift_hmh_gmres_NTsteric(sn_k,gn_k,cn_k,fn
c     $           ,mult,lcdim,npts,tolGMRES)
          else
            do ic = 1,lcdim
              call drift_hmh_gmres_newton(sn_k(1,ic),gn_k(1,ic)
     $             ,cn_k,fn(1,ic),mult,lcdim,npts,tolGMRES,ic,isep)
            enddo
          endif
          do ic = 1,lcdim
            call add2s2(cn_k(1,ic),sn_k(1,ic),alpha,npts) ! cp_k=cp_k+alpha*sp_k
          enddo

          call compute_f_couple(fn,cn_k,npts)
          call compute_gk(gn_k,fn,cn_k,cni,npts)  ! checking tol. + next NT iter
          do ic = 1,lcdim
            call chsign(gn_k(1,ic),npts)
          enddo

          do ic = 1,lcdim
            rnorms(ic) = glsc3(gn_k(1,ic),gn_k(1,ic),mult,npts)
            rnorms(ic) = sqrt(rnorms(ic))
          enddo
          rnorm = glmax(rnorms,lcdim)
          if (iterNT.eq.1) rinorm=rnorm
          ratio=rnorm/rinorm 

          if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) write(6,90)
     $       istep,iterNT,ratio,rnorm,rinorm,dtNT,dt

          if (ratio.lt.tolNT) goto 900
c          if (rnorm.lt.new_tol) goto 900

        enddo
 90     format('newton iter',2i6,1p5e12.4)   
900     continue

        do ic = 1,lcdim
          call copy(cni(1,ic),cn_k(1,ic),npts)   ! update cn_i for compute gn
        enddo                           ! for next pseudo-time step
        do ic = 1,lcdim
          call copy(cN(1,ic),cn_k(1,ic),npts) ! just for cem_out
        enddo

c       compute the norm of f(c_p)^{istep-1} and f(c_n)^{istep-1}
        do ic = 1,lcdim
          fnorms(ic) = glsc3(fn(1,ic),fn(1,ic),mult,npts)  
          fnorms(ic) = sqrt(fnorms(ic))
        enddo
        fnorm  = glmax(fnorms,lcdim)
        f_pre = f_now  ! set up old fnorm to f_pre
        f_now = fnorm  ! set up new fnorm to f_now

c       compute the CPU_time
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
      subroutine compute_gk(gon,fon,ckn,c0n,n)
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' !lcdim
      include 'NEWTON'
      integer n,ic
      real gon(lpts1,lcdim)
     $   , ckn(lpts1,lcdim)
     $   , c0n(lpts1,lcdim)
     $   , fon(lpts1,lcdim)

      do ic = 1,lcdim
        call add3s3(gon(1,ic),ckn(1,ic),fon(1,ic),c0n(1,ic)
     $       ,1.0,-1.0*dtNT,-1.0,n)
        call cmult(gon(1,ic),dtNTinv,n)   ! case with divided by dt_newton
      enddo

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f_sep(fon,cin,iflag,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON' 

      integer n,ic,iflag
      real fon(lpts1) ! f output
      real cin(lpts1,lcdim) ! c input
      real cn_bak(lpts1,lcdim) ! c output (sem_bdf)

      
      call cem_drift_sem_bdf1_newton(cin,n,iflag)

      call sub3(fon,cN(1,iflag),cin(1,iflag),n)
      call cmult(fon,dtinv,n)

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f_couple(fon,cin,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON' 
      integer n,ic
      real fon(lpts1,lcdim) ! f output
      real cin(lpts1,lcdim) ! c input
      real con(lpts1,lcdim) ! c output (sem_bdf)
      real cn_bak(lpts1,lcdim) ! c output (sem_bdf)

c      do ic = 1,lcdim
c        call copy(cN(1,ic),cin(1,ic),n)
c      enddo
      call cem_drift_sem_bdf1_newton(cin,n,0)

      do ic = 1,lcdim
        call sub3(fon(1,ic),cN(1,ic),cin(1,ic),n)
        call cmult(fon(1,ic),dtinv,n)
      enddo

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
      subroutine JacobiMatVec(Jp,p,uc,fi,n,iflag,isep)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(lpts1),p(lpts1)
      real     uc(lpts1,lcdim)      ! variable 
      real     fi(lpts1)            ! parameters
      real     uep(lpts1,lcdim)     ! uc + eps p
      real     foi(lpts1),fon(lpts1,lcdim)
      integer  n,iflag,ic,isep
      real     eps,pnorm,unorm,glsc3 

      pnorm = glsc3(p,p,mult,n)
      unorm = glsc3(uc(1,iflag),uc(1,iflag),mult,n)
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

c     if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) then
c         write(6,10) istep,newton_iter,eps,epsinv
c     endif

c 10  format (2i8,' eps=',1p2e17.7)
      do ic = 1,lcdim
        call copy(uep(1,ic),uc(1,ic),n)
      enddo
      call add3s2 (uep(1,iflag),uc(1,iflag),p,1.0,eps,n) ! uep = u + eps*p
      if (isep.eq.0) then
        call compute_f_couple(fon,uep,n)
      else
        do ic = 1,lcdim
          call compute_f_sep(fon(1,ic),uep,iflag,n)
        enddo
      endif
      call sub3   (foi,fon(1,iflag),fi,n)

      call cmult  (foi,epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
      call sub3   (Jp,p,foi,n)
      call cmult  (Jp,dtNTinv,n)     ! Jp = p/dt_newton - (fo-f)/eps
                                             ! case with divided by dt_newton
      return
      end
c-----------------------------------------------------------------------
      subroutine JacobiMatVec_steric(Jp,p,uc,fi,n)! w = A v 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(lpts1,lcdim),p(lpts1,lcdim)
      real     uc(lpts1,lcdim)      ! variable 
      real     fi(lpts1,lcdim)      ! parameters
      real     uep(lpts1,lcdim)     ! uc + eps p
      real     foi(lpts1,lcdim),fon(lpts1,lcdim)
      integer  n,nion,iflag,ic,isep
      real     pnorms(lcdim),unorms(lcdim)
      real     epss(lcdim),epssinv(lcdim)
      real     glsc3

      n = lx1*ly1*lz1*nelt
      do ic = 1,lcdim
        pnorms(ic) = glsc3(p(1,ic),p(1,ic),mult,n)
        unorms(ic) = glsc3(uc(1,ic),uc(1,ic),mult,n)
        pnorms(ic) = sqrt(pnorms(ic))
        unorms(ic) = sqrt(unorms(ic))
        epss(ic) = (1+unorms(ic))*1e-14   ! formula for varing eps
        epss(ic) = sqrt(epss(ic))         ! for using the case using 
                                          ! "ub" filled with previous
                                    !  step info first and 
                                    !  Dirichlet info on the boundary
        epss(ic) = epss(ic)/pnorms(ic)! this eps will blow up because of 
                                    ! pnorm = 0
        epssinv(ic) = 1./epss(ic)
      enddo

c     if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) then
c         write(6,10) istep,newton_iter,eps,epsinv
c     endif

c 10  format (2i8,' eps=',1p2e17.7)
      do ic = 1,lcdim
        call add3s2 (uep(1,ic),uc(1,ic),p(1,ic),1.0,epss(ic),n) ! uep = u + eps*p
      enddo

      call compute_f_couple(fon,uep,n)

      do ic = 1,lcdim
        call sub3   (foi(1,ic),fon(1,ic),fi(1,ic),n)
        call cmult  (foi(1,ic),epssinv(ic)*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
        call sub3   (Jp(1,ic),p(1,ic),foi(1,ic),n)
        call cmult  (Jp(1,ic),dtNTinv,n)     ! Jp = p/dt_newton - (fo-f)/eps
      enddo                                             ! case with divided by dt_newton
      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton
     $           (phi,res,uc,f,wt,nion,n,tol,iflag,isep)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'DRIFT'
      include 'NEWTON'
      include 'NGMRES'



      integer  n,noin,outer,iflag,isep
      real     phi(lpts1),res(lpts1),wt(lpts1)
      real     tol,alpha,l,temp
      real     eps,uc(lpts1,nion),f(lpts1)
      real*8   etime1,dnekclock

      integer  nion,i,j,k,iconv,iter,m
      real     rnorm,tolpss
      real     glsc3

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      
      call rzero(x_NT,n)
      call rzero(h_NT,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r_NT,res,n)                  ! r = res
         else
            !update residual
            call copy   (r_NT,res,n)                  ! r = res
            call JacobiMatVec(w_NT,x_NT,uc,f,n,iflag,isep)! w = A x
            call add2s2 (r_NT,w_NT,-1.,n) ! r = r - w
         endif

         gamma_NT(1) = glsc3(r_NT,r_NT,wt,n)                ! gamma  = (r,r)
         gamma_NT(1) = sqrt(gamma_NT(1))                 ! gamma  = sqrt{ (r,r) }

         tolps = 0.1*gamma_NT(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable

         !check for lucky convergence
         rnorm = 0.
         if(gamma_NT(1) .eq. 0.) goto 9000
         temp = 1./gamma_NT(1)
         call cmult2(v_NT(1,1),r_NT,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call JacobiMatVec(w_NT,v_NT(1,j),uc,f,n,iflag,isep) ! w = A v

            !modified Gram-Schmidt
            do i=1,j
               h_NT(i,j)=glsc3(w_NT,v_NT(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i
               call add2s2(w_NT,v_NT(1,i),-h_NT(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_NT(i,j)
               h_NT(i  ,j)=  c_NT(i)*temp + s_NT(i)*h_NT(i+1,j)
               h_NT(i+1,j)= -s_NT(i)*temp + c_NT(i)*h_NT(i+1,j)
            enddo
                                              !            ______
            alpha = sqrt(glsc3(w_NT,w_NT,wt,n))     ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h_NT(j,j)*h_NT(j,j)+alpha*alpha)
            temp = 1./l
            c_NT(j) = h_NT(j,j) * temp
            s_NT(j) = alpha  * temp
            h_NT(j,j) = l
            gamma_NT(j+1) = -s_NT(j) * gamma_NT(j)
            gamma_NT(j)   =  c_NT(j) * gamma_NT(j)

            rnorm = abs(gamma_NT(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_NT(1,j+1),w_NT,temp,n)   ! v    = w / alpha
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
            temp = gamma_NT(k)
            do i=j,k+1,-1
               temp = temp - h_NT(k,i)*c_NT(i)
            enddo
            c_NT(k) = temp/h_NT(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_NT,v_NT(1,i),c_NT(i),n)     ! x = x + c  z
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x_NT,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,iterNT,iter,tolpss
      endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end

cC-----------------------------------------------------------------------
c      subroutine drift_hmh_gmres_NTsteric(uC,fC,uref,fi,wt,nion,n,tol)
cc-----------------------------------------------------------------------
cc     For Steric term, involving cross term linear system 
cc      Hii = B + dt(1 + ste_aii)*A_Di
cc      Hij = dt*ste_aij*A_Di
cc    
cc      H = (Hij), f = [f1;f2;...fp], u = [u1;u2;...;up], p = # of species
cc      
cc     solving H*u = f
c      implicit none
c      include 'SIZE'
c      include 'TOTAL'
c      include 'FDMH1'
c      include 'DRIFT'
c      include 'MGMRES'
cC
c      integer  nion,n,outer
c      integer  iter,iconv
c      integer  i,j,k,m
c      integer  ic ! safe index fo nion
c      real     divex,etime_p,tolpss
c
c      real     uC(lpts1,lcdim)
c     $       , fC(lpts1,lcdim)
c     $       , mask(lpts1,lcdim)
c     $       , wt(lpts1)
c     $       , fi(lpts1,lcdim)
c     $       , uref(lpts1,lcdim)
c      real     tol,alpha,l,temp,rnorm
c      real     glsc3
c      real*8   etime1,dnekclock
cC
c
c      if (nid.eq.0) write(6,*) 'start: hmh_gmres_newton_steric'
c
c      n = nx1*ny1*nz1*nelv
c
c      etime1 = dnekclock()
c      etime_p = 0.
c      divex = 0.
c      iter  = 0
c      m     = lgmres
c
c      tolps = tol
c      tolpss= tolps
c      iconv = 0
c      call rzero(xxC,lpts*(LDIMT-1))
c      call rzero(h,m*m)
c
c      outer = 0
c      do while (iconv.eq.0.and.iter.lt.5000)
c         outer = outer+1
c         if(iter.eq.0) then
c            do ic=1,nion
c            call copy  (rC(1,ic),fC(1,ic),n)      ! r = res
c            enddo
c         else
c            !update residual
c            do ic=1,nion
c            call copy  (rC(1,ic),fC(1,ic),n)      ! r = res
c            enddo
c            call JacobiMatVec_steric(wC,xxC,uref,fi,n)! w = A x
c
c            do ic=1,nion
c            call dssum  (wC(1,ic),nx1,ny1,nz1)
c            call col2   (wC(1,ic),mask(1,ic),n)
c            enddo
c            do ic=1,nion
c            call add2s2 (rC(1,ic),wC(1,ic),-1.,n)   ! r = r - w
c            enddo
c         endif
c         gamma(1) = 0.0
c         do ic=1,nion
c           gamma(1) = gamma(1)+glsc3(rC(1,ic),rC(1,ic),wt,n)! gamma  = (r,r)
c         enddo
c         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
c         write(6,*) 'initial residule',gamma(1)
c                                                   !      1     include 'SIZE'
c
c         !check for lucky convergence
c         rnorm = 0.
c         if(gamma(1) .eq. 0.) goto 9000
c         temp = 1./gamma(1)
c         do ic=1,nion
c         call cmult2(vC(1,ic,1),rC(1,ic),temp,n)            !  v  = r / gamma
c         enddo
c         !write(6,*) 'start form m-th krylov subspace'
c         do j=1,m
c            iter = iter+1
c            call JacobiMatVec_steric(wC,vC(1,1,j),uref,fi,n)! w = A v 
c            do ic=1,nion
c            call dssum  (wC(1,ic),nx1,ny1,nz1)
c            call col2   (wC(1,ic),mask(1,ic),n)
c            enddo
c
cc           !modified Gram-Schmidt
c            do i=1,j
c               h(i,j)=0.0
c               do ic=1,nion
c               h(i,j) = h(i,j) + glsc3(wC(1,ic),vC(1,ic,i),wt,n) ! h    = (w,v )
c               enddo                                             ! i,j       i
c               do ic=1,nion
c               call add2s2(wC(1,ic),vC(1,ic,i),-h(i,j),n)   ! w = w - h    v
c               enddo                                        !         i,j  i
c            enddo
c
c            !apply Givens rotations to new column
c            do i=1,j-1
c               temp = h(i,j)
c               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
c               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
c            enddo
c            alpha=0.0
c            do ic=1,nion
c            alpha = alpha + glsc3(wC(1,ic),wC(1,ic),wt,n)  !            ______
c            enddo                                          ! alpha =  \/ (w,w)
c            alpha=sqrt(alpha)
c
c            if(alpha.eq.0.) goto 900  !converged
c            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
c            temp = 1./l
c            c(j) = h(j,j) * temp
c            s(j) = alpha  * temp
c            h(j,j) = l
c            gamma(j+1) = -s(j) * gamma(j)
c            gamma(j)   =  c(j) * gamma(j)
c
c            rnorm = abs(gamma(j+1))!*norm_fac     
cc           ratio = rnorm/div0
c            if ((nid.eq.0).and.(istep.le.5))
c     $           write (6,66) iter,tolpss,rnorm,istep
c   66       format(i5,1p2e12.5,i8,' Divergence')
c
c            if (rnorm .lt. tolps) goto 900 !converged
c            if (j.eq.m) goto 1000 !not converged, restart
c
c            temp = 1./alpha
c            do ic=1,nion
c            call cmult2(vC(1,ic,j+1),wC(1,ic),temp,n)   ! v    = w / alpha
c            enddo                                       !  j+1
c         enddo
cc        write(6,*) 'end of forming m-th krylov subspace'
c  900    iconv = 1
c 1000    continue
c
cc        back substitution
cc             -1
cc        c = H   gamma
cc        write(6,*) 'start solving least squre problem'
c         do k=j,1,-1
c            temp = gamma(k)
c            do i=j,k+1,-1
c               temp = temp - h(k,i)*c(i)
c            enddo
c            c(k) = temp/h(k,k)
c         enddo
c         !sum up Arnoldi vectors
c         do i=1,j
c            do ic=1,nion
c            call add2s2(xxC(1,ic),vC(1,ic,i),c(i),n)  ! x = x + c  z
c            enddo                                      !  i  i
c         enddo
cc        write(6,*) 'end of solving least squre problem'
c      enddo
c 9000 continue
c
c       do ic=1,nion
c       call copy(uC(1,ic),xxC(1,ic),n)
c       enddo
c
cc      call ortho   (res) ! Orthogonalize wrt null space, if present
c
c      etime1 = dnekclock()-etime1
c      if (nid.eq.0) then
c          if ((mod(istep,iocomm).eq.0).or.(istep.le.10))
c     $    write(6,9999) istep,iter,tolpss,etime_p,etime1
c      endif
c 9999 format(' ',' ',i9,' gmres   : iteration#',i5,1p3e12.4)
c
c      return
c      end
c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_newton(cinput,n,iflag)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
c      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'

      integer iflag,n,ic
      real cinput(lpts1,lcdim)

      do ic = 1,lcdim
       call copy(cN(1,ic),cinput(1,ic),n)
      enddo
      if (iflag.eq.0) then ! coupling every thing
        call cem_drift_op_bdf
      else ! seperate
        call cem_drift_source
        call cem_drift_poisson
        call cem_drift_sem_bdf1_init

        call cem_drift_rhs_bdf1(iflag)
        call cem_drift_lhs_bdf1(iflag)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_energy !FIXME: restr later, Lan
c-----------------------------------------------------------------------
c      implicit none
c      include 'SIZE'
c      include 'TOTAL'
c      include 'DRIFT'
c      include 'POISSON'
c
c      real energy_density(npts),energy
c      integer i,j,k
c      real dummy1(npts),dummy2(npts),dummy3(npts)
c      real dummy4(npts),dummy5(npts)
c      logical :: lexist
c      real glsum
c
c      call rzero(energy_density,npts)
c      call rzero(dummy1,npts)
c      call rzero(dummy2,npts)
c      call rzero(dummy3,npts)
c
c      call copy(dummy1,cN,npts)
c      call copy(dummy2,cP,npts)
c      call copy(dummy3,cP,npts)
c
c      call add2s2(dummy3,cN,-1.0,npts)
c      call add2s2(dummy3,ni, 1.0,npts)
c
c      call col2(dummy3,potent,npts)
c      call col2(dummy1,LOG(cN),npts)
c      call col2(dummy2,LOG(cP),npts)
c
c      call add2s2(energy_density,dummy1,1.0,npts)
c      call add2s2(energy_density,dummy2,1.0,npts)
c      call add2s2(energy_density,dummy3,0.5,npts)
c
c      call col2(energy_density,bmn,npts)
c      energy =  glsum(energy_density,npts)
c
c      if ((nid.eq.0).and.(mod(istep,iocomm).eq.0))
c     $    write(6,10) energy,time
c 10   format('     energy: ',1p1e18.7, ' at time=',1p1e12.4)
      return
      end
c-----------------------------------------------------------------------


