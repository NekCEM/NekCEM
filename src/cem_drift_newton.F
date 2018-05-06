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

c     Set ParamNT
      if (nid.eq.0) write(*,*)'ifparamNT',ifparamNT
      if (ifparamNT) then
        ! set parameters inside uservp and turn ifparamNT on
      else  ! load default parameters
        ifsep    = 1       ! 0=couple, 1=separate, Lan

        alpha    = 1.0      ! relaxation parameter \alpha \in (0,1]
c       nsteps   = 200      ! steps for pseudo-transient continuation
        maxit    = 20       ! iteration for Newton method
        jaceps   = 1e-5     ! perturbation parameter for Jacobi-free
        tolGMRES = 1e-10    ! tolerance for GMRES
        tolNT    = 1e-5     ! tolerance for Newton method

        dtNT     = param(1) ! pseudo-transient time step
                            !   it can grow as istep increase 
        f_pre    = 0.0      ! zero to set dtNT (or SER)
        max_dtNT = 1E10     ! Max. of dtNT, avoid NaN
        min_tol  = 1E-12    ! Min. of tol,  avoid NaN
      endif

c     Clear array
      do ic = 1,lcdim
        call rzero(sn_k(1,ic),npts)
        call rzero(gn_k(1,ic),npts)
        call rzero(cni (1,ic),npts)
        call rzero(cno (1,ic),npts)
        call rzero(cn_k(1,ic),npts)
        call rzero(fn  (1,ic),npts)
      enddo
      call izero(BDFcnt,3)  ! count the BDF1 usage

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

      real tt
      real glmax,glsum,glsc3 
      real rnorm,rnorms(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)
      real temp(lpts1)

      integer i,j,k,n,ic
      integer irestart ! print restart file for newton
      integer ifirst,ilast

      irestart = param(83)

c     parameter settings
      call cem_drift_newton_init

      epsinv = 1./jaceps
      dtinv = 1./dt
      dtNTinv = 1./dtNT


c     Adjust restart time step
      if (.not.usesteps) then
         if (ifrestart) then
            nsteps = int((fintim-restart_time)/dt)
            if (nid.eq.0) then
               write(6,*) 'fintim',fintim
               write(6,*) 'restart_time',restart_time
               write(6,*) 'fintim-restart_time',fintim-restart_time
            endif
         else
            nsteps = int((fintim)/dt)
            if (nid.eq.0) write(6,*) 'fintim',fintim
         endif
         if (nid.eq.0) write(6,*) 'adjusted nsteps',nsteps
      endif

      ifirst = irststep+1
      ilast = nsteps+irststep

      if (nid.eq.0) then
         write(6,*) '============================'
         write(6,*) '========  Begin Run ========'
         write(6,*) '============================'
         write(6,3) ifirst-1,ilast
      endif
    3 format(' istep =', i10,'  to',i10)

c     Initial time step info 
      istep = ifirst-1
      call userchk
      call cem_out

      if (nsteps.eq.0) call cem_end


c     compute initial vec and norm
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

C      call compute_energy  ! check the energy of the system in the initial stage


c     Start timer
      cpu_t = 0.0
      cpu_chk = 0.0


c     Start pseudo-time step
      istep = ifirst
      do while (.true.) 
        cpu_dtime = dclock()

c       SER for tau^{istep}
        if (abs(f_pre).lt.min_tol) then ! if not, use f_pre to do SER
          ! Do Nothing, set dtNT in init.
c            dtNT = param(1)
c            dt = dt_newton
          f_pre = 2*min_tol
          if (nid.eq.0) write(*,*)'f_pre',f_pre
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
         
c       Start Newton iteration
        do iterNT=1,maxit
          if (ifsteric) then
            call drift_hmh_gmres_NTsteric(sn_k,gn_k,cn_k,fn
     $           ,mult,lcdim,npts,tolGMRES)
          else
            do ic = 1,lcdim
              call drift_hmh_gmres_newton(sn_k(1,ic),gn_k(1,ic)
     $             ,cn_k,fn(1,ic),mult,lcdim,npts,tolGMRES,ic,ifsep)
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
        enddo                                    ! for next pseudo-time step
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
c        call compute_energy ! check energy


c       compute the CPU_time
        cpu_dtime = dclock()-cpu_dtime
        cpu_t = cpu_t+cpu_dtime
        cpu_t_step = cpu_t/istep
        cpu_p_t = glsum(cpu_t_step /npts,1)/np
        
        time = time + dt
    
        cpu_chk = dclock()
        call userchk
        cpu_chk = dclock()-cpu_chk

        call cem_out

        if (irestart.gt.0) then 
          if (mod(istep,irestart).eq.0) then
            if (nid.eq.0) then
              write(6,*)'Writing newton.restart! at istep=',istep
              call write_paramNT  ! write newton.restart

              write(6,*)'newton BDFcnt: #Poisson solver = ',BDFcnt(1)
              do ic = 1,lcdim
                write(6,991)ic,BDFcnt(ic+1)
              enddo
            endif
          endif
        endif 

        istep = istep+1
        if (lastep.eq.1) goto 901 !exit
        if (usesteps) then
          if (istep.gt.ilast) goto 901 !exit
        else
          if (time.ge.fintim) goto 901 !exit
        endif
      enddo

  901 continue
  991 format(' newton BDFcnt: #cN',I2,' solver = ',i10)

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
      real cin (lpts1,lcdim) ! c input
      real cout(lpts1,lcdim) ! c output
      real cn_bak(lpts1,lcdim) ! c output (sem_bdf)

      
      call cem_drift_sem_bdf1_newton(cin,cout,n,iflag)

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
      real cin (lpts1,lcdim) ! c input
      real cout(lpts1,lcdim) ! c output
      real con (lpts1,lcdim) ! c output (sem_bdf)
      real cn_bak(lpts1,lcdim) ! c output (sem_bdf)

c      do ic = 1,lcdim
c        call copy(cN(1,ic),cin(1,ic),n)
c      enddo
      call cem_drift_sem_bdf1_newton(cin,cout,n,0)

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
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'

      integer  n,noin,outer,iflag,isep
      real     phi(lpts1),res(lpts1),wt(lpts1)
      real     tol,alpha,l,temp
      real     eps,uc(lpts1,nion),f(lpts1)
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
            call JacobiMatVec(w,x,uc,f,n,iflag,isep)! w = A x
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

            call JacobiMatVec(w,v(1,j),uc,f,n,iflag,isep) ! w = A v

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

c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
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

C-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_NTsteric(uC,fC,uref,fi,wt,nion,n,tol)
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
      integer  nion,n,outer
      integer  iter,iconv
      integer  i,j,k,m
      integer  ic ! safe index fo nion
      real     divex,etime_p,tolpss

      real     uC(lpts1,lcdim)
     $       , fC(lpts1,lcdim)
     $       , mask(lpts1,lcdim)
     $       , wt(lpts1)
     $       , fi(lpts1,lcdim)
     $       , uref(lpts1,lcdim)
      real     tol,alpha,l,temp,rnorm
      real     glsc3
      real*8   etime1,dnekclock
C

      if (nid.eq.0) write(6,*) 'start: hmh_gmres_newton_steric'

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
            call JacobiMatVec_steric(wC,xxC,uref,fi,n)! w = A x

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
            call JacobiMatVec_steric(wC,vC(1,1,j),uref,fi,n)! w = A v 
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
c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_newton(cinput,coutput,n,iflag)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
c      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'
      include 'NEWTON' !BDFcnt

      integer iflag,n,ic
      real cinput (lpts1,lcdim)
      real coutput(lpts1,lcdim)

      do ic = 1,lcdim
       call copy(cN(1,ic),cinput(1,ic),n)
      enddo
      if (iflag.eq.0) then ! coupling every thing
        call cem_drift_op_bdf
        do ic = 1,lcdim+1
          BDFcnt(ic) = BDFcnt(ic) + 1
        enddo
      else ! seperate
        call cem_drift_source
        call cem_drift_poisson
        call cem_drift_sem_bdf1_init
        BDFcnt(1) = BDFcnt(1) + 1

        call cem_drift_rhs_bdf1(iflag)
        call cem_drift_lhs_bdf1(iflag)
        BDFcnt(iflag+1) = BDFcnt(iflag+1) + 1
      endif
      do ic = 1,lcdim
       call copy(coutput(1,ic),cN(1,ic),n)
      enddo

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
      subroutine read_paramNT
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      real paramNT(20)
      integer nn,i

C     Read paramNT

      open(unit=11,file='newton.restart')
      nn = 0
      do i=1,20
      read(11,*,end=99)paramNT(i)
      nn=nn+1
      enddo
   99 continue

C     Set paramNT

      ! 1~7 Pseudo-timestep
      ifsep    = int(paramNT(3))   ! (int) 0=couple, 1=separate, Lan
      dtNT     = paramNT(4)        ! pseudo-transient time step
      max_dtNT = paramNT(5)        ! Max. of dtNT, avoid NaN
      min_tol  = paramNT(6)        ! Min. of tol,  avoid NaN
      f_pre    = paramNT(7)        ! zero to set dtNT (or SER)

      ! 8~12 Newton iter
      alpha    = paramNT(8)        ! relaxation parameter \alpha \in (0,1]
      maxit    = int(paramNT(9))   ! (int) iteration for Newton method
      jaceps   = paramNT(10)       ! perturbation parameter for Jacobi-free
      tolNT    = paramNT(11)       ! tolerance for Newton method

      ! 13~20 else
      tolGMRES = paramNT(13)       ! tolerance for GMRES

      close(11)

      return
      end
c-----------------------------------------------------------------------
      subroutine write_paramNT
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'

      integer i
      real paramNT(20)

      open(unit=11,file='newton.restart')

      paramNT(1)  = 0
      paramNT(2)  = 0
      paramNT(3)  = ifsep
      paramNT(4)  = dtNT
      paramNT(5)  = max_dtNT
      paramNT(6)  = min_tol
      paramNT(7)  = f_pre
      paramNT(8)  = alpha
      paramNT(9)  = maxit
      paramNT(10) = jaceps
      paramNT(11) = tolNT
      paramNT(12) = 0
      paramNT(13) = tolGMRES
      paramNT(14) = 0
      paramNT(15) = 0
      paramNT(16) = 0
      paramNT(17) = 0
      paramNT(18) = 0
      paramNT(19) = 0
      paramNT(20) = 0

      write(11,3)int(paramNT(1)),'1: '
      write(11,3)int(paramNT(2)),'2: '
      write(11,3)int(paramNT(3)),
     $           '3: ifsep,(int) 0 = couple, 1 = seperate'
      write(11,7)paramNT(4),'4: dtNT, pseudo-timestep tau'
      write(11,7)paramNT(5),'5: max_dtNT, Max. of dtNT, avoid NaN'
      write(11,7)paramNT(6),'6: min_tol, Min. of tol,  avoid NaN'
      write(11,7)paramNT(7),'7: f_pre, f norm from previous step'
      write(11,7)paramNT(8),
     $           '8: alpha, relaxation parameter \alpha \in (0,1]'
      write(11,3)int(paramNT(9)),
     $           '9: maxit, (int) iteration for Newton method'
      write(11,7)paramNT(10),
     $           '10: jaceps, perturbation parameter for Jacobi-free'
      write(11,7)paramNT(11),'11: tolNT, tolerance for Newton method'
      write(11,7)paramNT(12),'12: '
      write(11,7)paramNT(13),'13: tolGMRES, tolerance for GMRES'
      write(11,7)paramNT(14),'14: '
      write(11,7)paramNT(15),'15: '
      write(11,7)paramNT(16),'16: '
      write(11,7)paramNT(17),'17: '
      write(11,7)paramNT(18),'18: '
      write(11,7)paramNT(19),'19: '
      write(11,7)paramNT(20),'20: '

    3 format(I4,'              ',A)
    7 format(E13.4'     ',A)

      close(11)

      return
      end
c-----------------------------------------------------------------------
