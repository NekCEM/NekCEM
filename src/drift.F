!> \file drift.F 
!!
!! \brief Schrodinger solver main subroutines          
!!
!! A more detailed description of the file


!> \brief Brief description of the subroutine
!!
!! A more detailed descripton goes here. This is an example
!! of a detailed description of a method
!!
!! @param param1 a description of the first parameter
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
!! @see cem_op()
!! @see READATW()
!! @see RDMESH()
!! @see RDMESHW()
c--------------------------------------------------------------------- 
      subroutine cem_drift_temporary                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'  
      include 'POISSON'  
      include 'ZPER'   
      integer  i, npts3, nxzfl3

      npts3 = 3*npts
      nxzfl3= 3*nxzfl

c...  these to be deleted all when IFHEX/IFTET are done

      if (nid.eq.0) write(6,*) 'set variables: cem_drift_temporary'
          call rzero (xmn ,npts)
          call rzero (ymn ,npts)
          call rzero (zmn ,npts)
          call rzero (bmn ,npts)
          call rzero (rxmn,npts)
          call rzero (rymn,npts)
          call rzero (rzmn,npts)
          call rzero (sxmn,npts)
          call rzero (symn,npts)
          call rzero (szmn,npts)
          call rzero (txmn,npts)
          call rzero (tymn,npts)
          call rzero (tzmn,npts)
          call rzero (jacm,npts)

          call rzero (cN   ,npts)
          call rzero (cP   ,npts)
          call rzero (cE   ,npts)
          call rzero (scN  ,npts)
          call rzero (scP  ,npts)
          call rzero (scE  ,npts)
          call rzero (cNQ  ,npts3)
          call rzero (cPQ  ,npts3)
          call rzero (cEQ  ,npts3)
          call rzero (scNQ ,npts3)
          call rzero (scPQ ,npts3)
          call rzero (scEQ ,npts3)
          call rzero (cNJ  ,npts3)
          call rzero (cPJ  ,npts3)
          call rzero (cEJ  ,npts3)
          call rzero (DcN  ,npts)
          call rzero (DcP  ,npts)
          call rzero (DcE  ,npts)
          call rzero (rD   ,npts)
          call rzero (rR   ,npts)
          call rzero (eG   ,npts)
          call rzero (eR   ,npts)
          call rzero (resN ,npts)
          call rzero (resP ,npts)
          call rzero (resE ,npts)
          call rzero (resNQ,npts3)
          call rzero (resPQ,npts3)
          call rzero (resEQ,npts3)
          call rzero (resNJ,npts3)
          call rzero (resPJ,npts3)
          call rzero (resEJ,npts3)

          call rzero (charge,npts)
          call rzero (potent,npts)
          call rzero (rhs   ,npts)
          call rzero (spotent,npts)
          call rzero (epotent,npts)

          call rzero (errN ,npts)
          call rzero (errP ,npts)
          call rzero (errE ,npts)
          call rzero (kcN  ,npts)
          call rzero (kcP  ,npts)
          call rzero (kcE  ,npts)
          call rzero (qbm1 ,npts)

          call rzero (fcN  ,nxzfl)
          call rzero (fcP  ,nxzfl)
          call rzero (fcE  ,nxzfl)
          call rzero (fcNJ ,nxzfl3)
          call rzero (fcPJ ,nxzfl3)
          call rzero (fcEJ ,nxzfl3)

          call rzero (w3mn ,nxyz )
          call rzero (unxm ,nxzfl)
          call rzero (unym ,nxzfl)
          call rzero (unzm ,nxzfl)
          call rzero (aream,nxzfl)

          do i=1,npts
             xmn (i)= xm1 (i,1,1,1)
             ymn (i)= ym1 (i,1,1,1)
             zmn (i)= zm1 (i,1,1,1)
             bmn (i)= bm1 (i,1,1,1)
             rxmn(i)= rxm1(i,1,1,1)
             rymn(i)= rym1(i,1,1,1)
             rzmn(i)= rzm1(i,1,1,1)
             sxmn(i)= sxm1(i,1,1,1)
             symn(i)= sym1(i,1,1,1)
             szmn(i)= szm1(i,1,1,1)
             txmn(i)= txm1(i,1,1,1)
             tymn(i)= tym1(i,1,1,1)
             tzmn(i)= tzm1(i,1,1,1)
             jacm(i)= jacm1(i,1,1,1)
          enddo

          do i=1,nxyz
             w3mn(i)= w3m1 (i,1,1)
          enddo

          do i=1,nxzfl
             unxm (i)= unx (i,1,1,1)
             unym (i)= uny (i,1,1,1)
             unzm (i)= unz (i,1,1,1)
             aream(i)= area(i,1,1,1)
          enddo
      if (nid.eq.0) write(6,*) 'done: temporary variables for geom'

c...  get global numbering on face and set material properties           
        call cem_set_fc_ptr      !global numbering index on face
        call cem_drift_uvp       !set parameters                
        call cem_set_neuman_fc_ptr !fixme by ying

c...  inverse mass matrix including material constants 
        if (ifrk.and.ifse) then
           call copy   (qbm1,bmn,npts)
           call dssum  (qbm1,nx1,ny1,nz1)
           call invers2(qbm1,qbm1,npts)
        else
           call invers2(qbm1,bmn,npts)
        endif
        if (nid.eq.0) write(6,*) 'done: inverse mass matrix'

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_uvp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE' 
      include 'EMWAVE' 
      include 'DRIFT'
      include 'POISSON'
      integer  i,j,k,l,ie,ieg
      real     tmp
c...  FIXME: give correct numbers for parameters  
c...  3/13/14: all these to be defined in user file.
      if(ifsol) then
      call rone(mu_n,npts)
      call rone(mu_p,npts)
      call rone(mu_e,npts)
      call rone(ni,npts)
      K_beta=1.0
      temperature=1.0
      tau_n = 1
      tau_p = 1
      call rone(charge,npts)
      call rone(d_permit,npts)
      call rone(d_permea,npts)
      else
      call rone(mu_n,npts)
      call rone(mu_p,npts)
      call rone(mu_e,npts)
      call rone(ni,npts)
      call cmult(ni, 1.E+26,npts)
      scale_x =1E-9
      scale_mu=1.e-9
      K_beta= 1.38044E-23
      temperature= 300.0
      e_charge=1.602e-19
      epsilon_0=8.854e-12
      constant_0 =e_charge/temperature/K_beta
      r_coulomb=constant_0*e_charge/epsilon_0/4.0/pi
      scale_time=e_charge*scale_x**2/(scale_mu*K_beta*temperature)
      scale_concent=epsilon_0*k_beta*temperature/e_charge**2/scale_x**2
      scale_potent =k_beta*temperature/e_charge
      scale_elect  =scale_potent/scale_x
      tau_n = 1.0E-4/scale_time
      tau_p = 1.0E-4/scale_time
      tau_cE= 1.0E-6/scale_time
      mu_n0=7.7e-9/scale_mu
      mu_p0=5.1e-9/scale_mu
      mu_e0=3.86e-9/scale_mu
      call cmult(mu_n,mu_n0,npts)
      call cmult(mu_p,mu_p0,npts)
      call cmult(mu_e,mu_e0,npts)
      call cmult(ni, 1./scale_concent,npts)
      write(6,*) 'scale_x',scale_x
      write(6,*) 'scale_mu',scale_mu
      write(6,*) 'scale_time',scale_time
      write(6,*) 'scale_concent',scale_concent
      write(6,*) 'scale_potent',scale_potent
      write(6,*) 'scale_elect',scale_elect
      write(6,*) 'tau_n',tau_n
      write(6,*) 'tau_p',tau_p
      write(6,*) 'tau_e',tau_ce
      write(6,*) 'e_charge',e_charge
      write(6,*) 'k_beta',k_beta
      write(6,*) 'constant_0',constant_0
      write(6,*) 'r_coulomb',r_coulomb
      write(6,*) 'mu_n0',mu_n0
      write(6,*) 'mu_p0',mu_p0
      call rone(charge,npts)
c      call cmult(charge,1.602e-19,npts)
      do ie=1,nelt
         if (if_in_region(ie).eq.1) then !donner
            do i=1,nxyz
            j=(ie-1)*nxyz+i
            d_permit(j)= 6.5
            d_permea(j)= 1.0
            enddo
         elseif (if_in_region(ie).eq.-1) then !accepter
            do i=1,nxyz
            j=(ie-1)*nxyz+i
            d_permit(j)= 3.9
            d_permea(j)= 1.0
            enddo
         else
             if (nid.eq.0) write(6,*) if_in_region(ie),
     $        'if_in_region not defined in usr'
         endif
      enddo
      endif
      call copy(diff_n,mu_n,npts)
      call copy(diff_p,mu_p,npts)
      call copy(diff_e,mu_e,npts)
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_op_rk
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      include 'RK5'
      integer  ii


      if (ifse) then
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson
         do ii=1,5
            call rk_c (ii)
            call cem_drift_sem_op
            call rk_drift_ab(ii)
            call col2  (cN,dmask(1,2),npts)
            call col2  (cP,dmask(1,3),npts)
            if (ifexct) call col2(cE,dmask(1,4),npts)
            if (ifpec.and.ifinhom) then
               call add2s2(cN,cN_b,1.0,npts)
               call add2s2(cP,cP_b,1.0,npts)
               if (ifexct) call add2s2(cE,cE_b,1.0,npts)
            endif
         enddo
c     call cem_drift_poisson1 ! variable function for h1 with h2=0
      else
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson
         do ii=1,5
            call rk_c (ii)
            call cem_drift_dg_op
            call rk_drift_ab(ii)
         enddo
      endif

      if (iffilter) call drift_filter(0.01)

      return
      end

c--------------------------------------------------------------------- 
c.... 1st/2nd order semi-implicit BDF scheme
      subroutine cem_drift_op_bdf
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      if     (ifbdf1) then      ! param(17)=-1: BDF1
         if(ifsol) then
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson ! solve potent and electric feild cEE
         call cem_drift_sem_bdf1
         else
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson ! solve potent and electric feild cEE
         call cem_drift_mobility_update
c         call cem_drift_current_density    
         call cem_drift_force  
         call cem_drift_sem_bdf1
         endif
         if (iffilter) call drift_filter(0.01)
         return

      elseif (ifbdf2) then      ! param(17)=-2: BDF2

       if (istep.eq.1) then
         call cem_drift_sem_bdf2_init
         return
       else      
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson ! solve potent and electric feild cEE
         call cem_drift_sem_bdf2
         return
       endif
      endif
 
 
      if (param(17).eq.-3) then
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson
         call cem_drift_sem_bdf2_124
         return
      endif
      if (param(17).eq.-4) then
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson_bdf2_24
         call cem_drift_sem_bdf2_124
         return
      endif
      if (param(17).eq.-5) then
         call cem_drift_poisson_bdf2_24
         call cem_drift_sem_bdf2_124
         return
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_dg_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      common /scrns/ srflxs(12*2*ldim*lx1*lz1*lelt)        
      real           srflxs
      call cem_drift_dg
      call cem_drift_restrict_to_face
      call cem_drift_flux (srflxs)
      call cem_drift_add_flux_to_res (srflxs)
      call cem_drift_source_dg   
      call cem_drift_invqmass_dg
      
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      integer i
      call cem_drift_sem1
      call cem_drift_sem_source1   
      call col2(cN,dmask(1,2),npts)
      call col2(cP,dmask(1,3),npts)
      call col2(cE,dmask(1,4),npts)
      call cem_drift_invqmass_sem
      call col2(cN,dmask(1,2),npts)
      call col2(cP,dmask(1,3),npts)
      call col2(cE,dmask(1,4),npts)
c      call cem_drift_poisson1 ! return potent and electric feild cEE by constant coeff poisson equation
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_sem1               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      real    wk0(lpts4),wk1(lpts4),wk2(lpts4),wk3(lpts4)
      real    wkxyz(lpts4*3),wkx(lpts4),wky(lpts4),wkz(lpts4)
      integer  i,isd,ifld,imsh,maxit
      real    xx,yy,aa,bb,cc,tt
c... 2d case for now
      isd   = 1
      imsh  = 1
c... pre-computation for DtN stuffs
      do ifld=1,nfield
      call rone     (dmask(1,ifld),npts)      ! for grids                         
      call dd_setdmaskp(dmask(1,ifld),npts,ifld)
      enddo

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)

      do i=1,npts
         h1n(i)=-diff_n(i)/charge(i)
         h1p(i)=-diff_p(i)/charge(i)
      enddo
      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call axhelm3(resN,cN,h1n,h2n,imsh,isd)!compute AU=-h1\delta u
      call axhelm3(resP,cP,h1p,h2p,imsh,isd)!compute AU=-h1\delta u
      if (ifdealias) then
      call drift_intp_u_grad_v(wk0,potent,cN)
      call add2s2 (resN,wk0,+1.0,npts)
      call drift_intp_u_grad_v(wk0,potent,cP)
      call add2s2 (resP,wk0,-1.0,npts)
      else
      call axhelm3(wk0,potent,cN,h2n,imsh,isd)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0,+1.0,npts)
      call axhelm3(wk0,potent,cP,h2p,imsh,isd)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0,npts)
      endif

      call dssum  (resN,nx1,ny1,nz1)
      call dssum  (resP,nx1,ny1,nz1)
      if (ifexct) then
      do i=1,npts
         h1e(i)=-diff_e(i)/charge(i)
      enddo
      call rzero(h2e,npts)
      call axhelm3(resE,cE,h1e,h2e,imsh,isd)!compute AU=-h1\delta u
      call dssum  (resE,nx1,ny1,nz1)
      endif

      return
      end
c--------------------------------------------------------------------- 
      subroutine elliptic2d(au,u,a,b,c,imsh,isd)              
C------------------------------------------------------------------
C
C     Compute the 2nd elliptic  matrix-vector product,
C     AU = -\div( a \nabla u)  + c \cdot \grad u + b u , for NEL elements.
C     a, b: scaler, c: vector
C------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,k,imsh,isd
      real tmp1(lpts),wrkdx(lpts),wrkdy(lpts),wrkdz(lpts)
      real a(1),b(1),c(1),u(1),au(1)
      call copy(tmp1,u,npts)
      call axhelm3(au,tmp1,a,b,imsh,isd)!compute AU=-a\delta u+ b u
      call drift_grad(wrkdx,wrkdy,wrkdz,tmp1)  
      k=npts
      if (if3d) then
      do i=1,npts
         tmp1(i)=wrkdx(i)*c(i)+wrkdy(i)*c(i+k)+wrkdz(i)*c(i+2*k)
      enddo
      else
      do i=1,npts
         tmp1(i)=wrkdx(i)*c(i)+wrkdy(i)*c(i+k)
      enddo
      endif
      call col2(tmp1,bmn,npts) !tmp1=tmp1*bmn
      call add2s2 (au,tmp1,1.0,npts) !Au=Au+c \cdot \grad u
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_dg                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i
      real xx,yy,tt,cc,glamax!fixme delte later 
      
      call drift_wght_div  (resN,wNJ(1,1),wNJ(1,2),wNJ(1,3),
     $                      cNJ(1,1),cNJ(1,2),cNJ(1,3))
      call drift_wght_div  (resP,wPJ(1,1),wPJ(1,2),wPJ(1,3),
     $                      cPJ(1,1),cPJ(1,2),cPJ(1,3))
      call drift_wght_div  (resE,wEJ(1,1),wEJ(1,2),wEJ(1,3),
     $                      cEJ(1,1),cEJ(1,2),cEJ(1,3))
      
      call drift_wght_grad (resNQ(1,1),resNQ(1,2),resNQ(1,3),cN)
      call drift_wght_grad (resPQ(1,1),resPQ(1,2),resPQ(1,3),cP)
      call drift_wght_grad (resEQ(1,1),resEQ(1,2),resEQ(1,3),cE)
      call drift_wght_grad (resEN(1,1),resEN(1,2),resEN(1,3),potent)
      if (if3d) then
      do i=1,npts
         resNQ(i,1)=diff_n(i)/charge(i)*resNQ(i,1)
     $              -mu_n(i)*cN(i)*resEN(i,1)
         resNQ(i,2)=diff_n(i)/charge(i)*resNQ(i,2)
     $              -mu_n(i)*cN(i)*resEN(i,2)
         resNQ(i,3)=diff_n(i)/charge(i)*resNQ(i,3)
     $              -mu_n(i)*cN(i)*resEN(i,3)
         resPQ(i,1)=diff_p(i)/charge(i)*resPQ(i,1)
     $              +mu_p(i)*cP(i)*resEN(i,1)
         resPQ(i,2)=diff_p(i)/charge(i)*resPQ(i,2)
     $              +mu_p(i)*cP(i)*resEN(i,2)
         resPQ(i,3)=diff_p(i)/charge(i)*resPQ(i,3)
     $              +mu_p(i)*cP(i)*resEN(i,3)
         resEQ(i,1)=diff_e(i)/charge(i)*resEQ(i,1)
         resEQ(i,2)=diff_e(i)/charge(i)*resEQ(i,2)
         resEQ(i,3)=diff_e(i)/charge(i)*resEQ(i,3)
      enddo
      else
      do i=1,npts
         resNQ(i,1)=diff_n(i)/charge(i)*resNQ(i,1)
     $              -mu_n(i)*cN(i)*resEN(i,1)
         resNQ(i,2)=diff_n(i)/charge(i)*resNQ(i,2)
     $              -mu_n(i)*cN(i)*resEN(i,2)
         resPQ(i,1)=diff_p(i)/charge(i)*resPQ(i,1)
     $              +mu_p(i)*cP(i)*resEN(i,1)
         resPQ(i,2)=diff_p(i)/charge(i)*resPQ(i,2)
     $              +mu_p(i)*cP(i)*resEN(i,2)
         resEQ(i,1)=diff_e(i)/charge(i)*resEQ(i,1)
         resEQ(i,2)=diff_e(i)/charge(i)*resEQ(i,2)
      enddo
      endif

c...  FIXME
c      call cem_drift_force 

      return
      end

c-----------------------------------------------------------------------
      subroutine drift_grad(w1,w2,w3,u1)
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
      subroutine drift_wght_grad(w1,w2,w3,u1)
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
      subroutine drift_wght_div (w0,w1,w2,w3,u1,u2,u3)
      implicit none
c     div:   w = div(u)/jacm
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer  k,i,j,e,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real    ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3
      real    u1rw,u1sw,u1tw 
      real    u2rw,u2sw,u2tw 
      real    u3rw,u3sw,u3tw 

      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(ur1,us1,ut1,u1,nn,e,dxm1,dxtm1)
         call local_grad3(ur2,us2,ut2,u2,nn,e,dxm1,dxtm1)
         call local_grad3(ur3,us3,ut3,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            u1rw = ur1(i)*w3mn(i)
            u1sw = us1(i)*w3mn(i)
            u1tw = ut1(i)*w3mn(i)
            u2rw = ur2(i)*w3mn(i)
            u2sw = us2(i)*w3mn(i)
            u2tw = ut2(i)*w3mn(i)
            u3rw = ur3(i)*w3mn(i)
            u3sw = us3(i)*w3mn(i)
            u3tw = ut3(i)*w3mn(i)

            w1(k) = (  u1rw*rxmn(k)
     $               + u1sw*sxmn(k)
     $               + u1tw*txmn(k))

            w2(k) = (  u2rw*rymn(k)
     $               + u2sw*symn(k)
     $               + u2tw*tymn(k))

            w3(k) = (  u3rw*rzmn(k)
     $               + u3sw*szmn(k)
     $               + u3tw*tzmn(k))

            w0(k) = w1(k)+w2(k)+w3(k)
         enddo

         enddo

      else

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
         call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            u1rw = ur1(i)*w3mn(i)
            u1sw = us1(i)*w3mn(i)
            u2rw = ur2(i)*w3mn(i)
            u2sw = us2(i)*w3mn(i)

            w1(k) = (  u1rw*rxmn(k)
     $               + u1sw*sxmn(k))

            w2(k) = (  u2rw*rymn(k)
     $               + u2sw*symn(k))
         
            w0(k) = w1(k)+w2(k)       
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_restrict_to_face
      implicit none
c     Restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      integer  i,j

      do j=1,ncemface
          i= cemface(j)
           fcNJ(j,1)= cNJ(i,1)
           fcNJ(j,2)= cNJ(i,2)
           fcNJ(j,3)= cNJ(i,3)
           fcPJ(j,1)= cPJ(i,1)
           fcPJ(j,2)= cPJ(i,2)
           fcPJ(j,3)= cPJ(i,3)
           fcEJ(j,1)= cEJ(i,1)
           fcEJ(j,2)= cEJ(i,2)
           fcEJ(j,3)= cEJ(i,3)
           fcN(j)   = cN(i)*diff_n(i)/charge(i)
           fcP(j)   = cP(i)*diff_p(i)/charge(i)
           fcE(j)   = cE(i)*diff_e(i)/charge(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_add_flux_to_res(srflx)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'

      real     srflx(1), a
      integer  i,j,k

      k   = nxzfl            

      if (if3d) then
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resN (i)   = resN (i)   + a*srflx( 0*k+j)
         resP (i)   = resP (i)   + a*srflx( 1*k+j)
         resE (i)   = resE (i)   + a*srflx( 2*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx( 3*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx( 4*k+j)
         resNQ(i,3) = resNQ(i,3) + a*srflx( 5*k+j)
         resPQ(i,1) = resPQ(i,1) + a*srflx( 6*k+j)
         resPQ(i,2) = resPQ(i,2) + a*srflx( 7*k+j)
         resPQ(i,3) = resPQ(i,3) + a*srflx( 8*k+j)
         resEQ(i,1) = resEQ(i,1) + a*srflx( 9*k+j)
         resEQ(i,2) = resEQ(i,2) + a*srflx(10*k+j)
         resEQ(i,3) = resEQ(i,3) + a*srflx(11*k+j)
         enddo
      else
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resN (i)   = resN (i)   + a*srflx( 0*k+j)
         resP (i)   = resP (i)   + a*srflx( 1*k+j)
         resE (i)   = resE (i)   + a*srflx( 2*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx( 3*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx( 4*k+j)
         resPQ(i,1) = resPQ(i,1) + a*srflx( 5*k+j)
         resPQ(i,2) = resPQ(i,2) + a*srflx( 6*k+j)
         resEQ(i,1) = resEQ(i,1) + a*srflx( 7*k+j)
         resEQ(i,2) = resEQ(i,2) + a*srflx( 8*k+j)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source_dg    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'RK5'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      real    dumm1,dummy1(lpts),dummy2(lpts),dummy3(lpts)
      integer i
      real tt,xx,yy,aa,bb,cc !delete later fixme
      if (.not.ifsrc) return
     
      call usersrc (dumm1,wk1,wk2,wk3,dummy1,dummy2,dummy3)
c      call col2    (wk1,bmn,npts) !fixme for dg should including the boundary info or not.. 
      call add2s2  (resN,wk1,1.,npts) 
c      call col2    (wk2,bmn,npts) 
      call add2s2  (resP,wk2,1.,npts) 
c      call col2    (wk3,bmn,npts) 
      call add2s2  (resE,wk3,1.,npts) 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_source1  
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'
      include 'RK5'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      real    dummy1(lpts)
      real    dummy2(lpts)
      real    dummy3(lpts)
      integer i,dumm1
      if (.not.ifsrc) return

      call usersrc (dumm1,wk1,wk2,wk3,dummy1,dummy2,dummy3)
      call dssum   (wk1,nx1,ny1,nz1)
      call dssum   (wk2,nx1,ny1,nz1)
      call add2s2  (resN,wk1,1.,npts) 
      call add2s2  (resP,wk2,1.,npts)
      if(ifexct) then
      call dssum   (wk3,nx1,ny1,nz1)
      call add2s2  (resE,wk3,1.,npts)
      endif 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_invqmass_sem
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'
      call col2(resN,qbm1,npts) !resN=resN*qbm1
      call col2(resP,qbm1,npts) !resP=resP*qbm1
      if (ifexct) call col2(resE,qbm1,npts) !resE=resE*qbm1
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_invqmass_dg
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      call col2(resN      ,qbm1,npts)       !resN =resN*qbm1
      call col2(resNQ(1,1),qbm1,npts) !resNQ=resNQ*qbm1
      call col2(resNQ(1,2),qbm1,npts) !array for copy current
      call col2(resNQ(1,3),qbm1,npts) ! not used for RK45
      call col2(resP      ,qbm1,npts) !resP =resP*qbm1
      call col2(resPQ(1,1),qbm1,npts) !resPQ=resPQ*qbm1
      call col2(resPQ(1,2),qbm1,npts) 
      call col2(resPQ(1,3),qbm1,npts) 
      call col2(resE      ,qbm1,npts) !resE =resE*qbm1
      call col2(resEQ(1,1),qbm1,npts) !resEQ=resEQ*qbm1
      call col2(resEQ(1,2),qbm1,npts) 
      call col2(resEQ(1,3),qbm1,npts) 
      call copy(cNJ,resNQ,npts*3)
      call copy(cPJ,resPQ,npts*3)
      call copy(cEJ,resEQ,npts*3)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_flux (srflx)  
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'DRIFT'
      include 'RK5'   
     
      real*8  dclock,t0
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer i,j,k,e,f,i0
      real    unx0,uny0,unz0
      real    tmp,tmpi1,tmpi2,tmpi3,tmpr1,tmpr2,tmpr3
      real    tmpur,tmpui,grad_r

      integer icalld
      save    icalld
      data    icalld/0/
      character CB*3

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0
      !call usersol(RKtime,sQr,sQi,sQi,sUr,sUi,sUi) !FIXME

      k = nxzfl

      if (if3d) then
        do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
           unz0=unzm(i)
           srflx( 0*k+i)=-(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)+unz0*fcNJ(i,3))
           srflx( 1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)+unz0*fcPJ(i,3))
           srflx( 2*k+i)=-(unx0*fcEJ(i,1)+uny0*fcEJ(i,2)+unz0*fcEJ(i,3))
           srflx( 3*k+i)=-unx0*fcN(i) 
           srflx( 4*k+i)=-uny0*fcN(i) 
           srflx( 5*k+i)=-unz0*fcN(i) 
           srflx( 6*k+i)=-unx0*fcP(i) 
           srflx( 7*k+i)=-uny0*fcP(i) 
           srflx( 8*k+i)=-unz0*fcP(i) 
           srflx( 9*k+i)=-unx0*fcE(i) 
           srflx(10*k+i)=-uny0*fcE(i) 
           srflx(11*k+i)=-unz0*fcE(i) 

         enddo

        else
         do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
C.........FIXME figure out the right sign +/-
           srflx(0*k+i)=-(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)) 
c           srflx(0*k+i)=srflx(0*k+i)+fcN(i)*(unx0+uny0) 
           srflx(1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)) 
           srflx(2*k+i)=-(unx0*fcEJ(i,1)+uny0*fcEJ(i,2)) 
           srflx(3*k+i)=-unx0*fcN(i)  ! NxfHx
           srflx(4*k+i)=-uny0*fcN(i)  ! NxfHx
           srflx(5*k+i)=-unx0*fcP(i)  ! NxfHx
           srflx(6*k+i)=-uny0*fcP(i)  ! NxfHx
           srflx(7*k+i)=-unx0*fcE(i)  ! NxfHx
           srflx(8*k+i)=-uny0*fcE(i)  ! NxfHx

         enddo
        endif

        t0=dclock()

        if (if3d) then
           call gs_op_fields(gsh_face,srflx,nxzfl,12,1,1,0)
        else
           call gs_op_fields(gsh_face,srflx,nxzfl,9,1,1,0)
        endif
        call measure_comm(t0)
        
        if (IFPEC) call cem_drift_flux_pec(srflx)
         
        k = nxzfl

        if (if3d) then
         do i=1,nxzfl
         srflx( 0*k+i)= 0.5*srflx( 0*k+i)
         srflx( 1*k+i)= 0.5*srflx( 1*k+i) 
         srflx( 2*k+i)= 0.5*srflx( 2*k+i) 
         srflx( 3*k+i)= 0.5*srflx( 3*k+i) 
         srflx( 4*k+i)= 0.5*srflx( 4*k+i) 
         srflx( 5*k+i)= 0.5*srflx( 5*k+i) 
         srflx( 6*k+i)= 0.5*srflx( 6*k+i) 
         srflx( 7*k+i)= 0.5*srflx( 7*k+i) 
         srflx( 8*k+i)= 0.5*srflx( 8*k+i) 
         srflx( 9*k+i)= 0.5*srflx( 9*k+i) 
         srflx(10*k+i)= 0.5*srflx(10*k+i) 
         srflx(11*k+i)= 0.5*srflx(11*k+i) 
         enddo
        else

         do i=1,nxzfl
         srflx( 0*k+i)= 0.5*srflx( 0*k+i)
         srflx( 1*k+i)= 0.5*srflx( 1*k+i) 
         srflx( 2*k+i)= 0.5*srflx( 2*k+i) 
         srflx( 3*k+i)= 0.5*srflx( 3*k+i) 
         srflx( 4*k+i)= 0.5*srflx( 4*k+i) 
         srflx( 5*k+i)= 0.5*srflx( 5*k+i) 
         srflx( 6*k+i)= 0.5*srflx( 6*k+i) 
         srflx( 7*k+i)= 0.5*srflx( 7*k+i) 
         srflx( 8*k+i)= 0.5*srflx( 8*k+i) 
         enddo
        endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_flux_pec(srflx)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'DRIFT'  
      include 'PARALLEL'
      include 'NEKUSE'
      include 'RK5'    
 
      real     srflx(1)
      integer  nxyzf,nface,nxz,k,i0,i,j,e,f,ef          
      real     unx0,uny0,unz0
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      k    = nxzfl
       call usersol(RKtime,scn,scp,sce,scn,scp,sce)

      if (if3d) then
         do e= 1,nelt
         do f= 1,nfaces
         do j= 1,nxzf

         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,2)
         if (CB.eq.'PEC') then

             unx0=unxm(i)
             uny0=unym(i)
             unz0=unzm(i)
             srflx( 0*k+i)= 0.d0
             srflx( 1*k+i)= 0.d0
             srflx( 2*k+i)= 0.d0
             srflx( 3*k+i)= 2.d0*(srflx( 3*k+i)+unx0*scN(i0))
             srflx( 4*k+i)= 2.d0*(srflx( 4*k+i)+uny0*scN(i0))
             srflx( 5*k+i)= 2.d0*(srflx( 5*k+i)+unz0*scN(i0))
             srflx( 6*k+i)= 2.d0*(srflx( 6*k+i)+unx0*scP(i0))
             srflx( 7*k+i)= 2.d0*(srflx( 7*k+i)+uny0*scP(i0))
             srflx( 8*k+i)= 2.d0*(srflx( 8*k+i)+unz0*scP(i0))
             srflx( 9*k+i)= 2.d0*(srflx( 9*k+i)+unx0*scE(i0))
             srflx(10*k+i)= 2.d0*(srflx(10*k+i)+uny0*scE(i0))
             srflx(11*k+i)= 2.d0*(srflx(11*k+i)+unz0*scE(i0))

         endif
         enddo
         enddo
         enddo
      else
         do e= 1,nelt
         do f= 1,nfaces
         do j= 1,nxzf

         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,2)
         if (CB.eq.'PEC') then

             unx0=unxm(i)
             uny0=unym(i)
             srflx(0*k+i)= 0.d0
c             srflx(0*k+i)= 2.d0*(srflx( 0*k+i)-(unx0+uny0)*scN(i0))
             srflx(1*k+i)= 0.d0
             srflx(2*k+i)= 0.d0
             srflx(3*k+i)= 2.d0*(srflx(3*k+i)+unx0*scN(i0))
             srflx(4*k+i)= 2.d0*(srflx(4*k+i)+uny0*scN(i0))
             srflx(5*k+i)= 2.d0*(srflx(5*k+i)+unx0*scP(i0))
             srflx(6*k+i)= 2.d0*(srflx(6*k+i)+uny0*scP(i0))
             srflx(7*k+i)= 2.d0*(srflx(7*k+i)+unx0*scE(i0))
             srflx(8*k+i)= 2.d0*(srflx(8*k+i)+uny0*scE(i0))

         endif
           
         enddo
         enddo
         enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine rk_drift_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'
      include 'DRIFT'

      integer  ii   
      real     ca,cb
c     kHX=rk4a(i)*kHX

      ca = rk4a(ii)
      cb = rk4b(ii)
      call rk4_upd(cN,kcN,resN,cb,ca,dt,npts)
      call rk4_upd(cP,kcP,resP,cb,ca,dt,npts)
      if (ifexct) call rk4_upd(cE,kcE,resE,cb,ca,dt,npts)

      return
      end
c----------------------------------------------------------------------
      subroutine cem_drift_op_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'

      if (nid.eq.0) write(6,*) 'compute eigenvalues'

c     call cem_drift_eig        !call drift_eig_arpack
c     call cem_end

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_force  
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real glamax,xx,yy,tt,tmp,tmp0
      real rG(lpts4),rSRH(lpts4),rrad(lpts4),raug(lpts4)
      call rzero(rD,npts)
      call rzero(rR,npts)
      call rzero(eG,npts)
      call rzero(eR,npts)
c1.....for exciton
      if (ifexct) then
         !Onsager dissociation
         call cem_drift_dissociation(rD)
         !generation for exciton
         call cem_drift_photogeneration(eG)
         !relaxation for exciton
         do i=1,npts
         eR(i)=cE(i)/tau_cE
         enddo
      endif
c2.....for electron and hole
c2.1.... electrolyte Dissociation (depend on electrical field and exciton) 
c2.2.....carrier recombination (depend on electron and hole )
c      call rzero(rSRH,npts)
c      call rzero(rrad,npts)
c      call rzero(raug,npts)
c      call cem_drift_recombination_SRH(rSRH)
c      call cem_drift_recombination_rad(rrad)
c      call cem_drift_recombination_aug(raug)
c      call add2s2(rR,rSRH,1.0,npts)
c      call add2s2(rR,rrad,1.0,npts)
c      call add2s2(rR,raug,1.0,npts)
       tmp0=e_charge/epsilon_0*scale_time
     $      *scale_mu*scale_concent**2/scale_concent !tmp0=1.0
       tmp0=1.0
       do i=1,npts
          rR(i)=tmp0*(mu_n(i)+mu_p(i))*cN(i)*cP(i)/d_permit(i)
       enddo
c2.2.....carrier generation (depend on electron and hole )
c      call rzero(rG,npts)
c      call cem_drift_generation(rG)


c......adding to the source term
      call col2(eG,bmn,npts)
      call col2(rD,bmn,npts)
      call col2(rR,bmn,npts)
      call col2(eR,bmn,npts)
      call add2s2(resN,rD, 1.0,npts)
      call add2s2(resN,rR, -1.0,npts)
      call add2s2(resP,rD, 1.0,npts)
      call add2s2(resP,rR, -1.0,npts)
      if (ifexct) then
      call add2s2(resE,eG, 1.0,npts)
      call add2s2(resE,rD, -1.0,npts)
      call add2s2(resE,rR, 0.25,npts)
      call add2s2(resE,eR, -1.0,npts)
      endif
c      write(6,*) 'xxx1',istep,glamax(EG,npts)
c      write(6,*) 'xxx1',istep,glamax(rD,npts)
c      write(6,*) 'xxx1',istep,glamax(rR,npts)
c      write(6,*) 'xxx1',istep,glamax(Er,npts)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_recombination_SRH(out)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     glamax,out(lpts4),wk1(lpts4),wk2(lpts4)
      call copy(out,cN,npts)
      call col2(out,cP,npts)
      call copy(wk1,ni,npts)
      call col2(wk1,wk1,npts)
      call add2s2(out,wk1,-1.0,npts) !out=cN*cP-ni**2

      call copy(wk1,cN,npts)
      call add2s2(wk1,ni,1.0,npts)
      call cmult(wk1,tau_p,npts)
      call copy(wk2,cP,npts)
      call add2s2(wk2,ni,1.0,npts)
      call cmult(wk2,tau_n,npts)
      call add2s2(wk1,wk2,1.0,npts)
      call invcol1(wk1,npts)         !wk1=1.0/(tau_p*(cP+ni)+tau_n*(cN+ni))
      call col2(out,wk1,npts)
      return
      end
C----------------------------------------------------------------------
      subroutine cem_drift_recombination_rad(out)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     out(lpts4)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_recombination_aug(out)
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     out(lpts4)
      return
      end
C----------------------------------------------------------------------
      real function cem_drift_schottky_d_bc(i0)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i0,i,j,ie 
      real     tmp0,tmp1,tmp2
      real     boltz_c
      tmp0=constant_0*r_coulomb*scale_elect
      if (abs(cEE(i0,1)).le.1.e-8) then
      cem_drift_schottky_d_bc=ni(i0)
      else
      tmp1=tmp0*abs(cEE(i0,1))!f
      tmp2=1./tmp1+1./sqrt(tmp1)-1./tmp1*sqrt((1.+2.0*sqrt(tmp1)))!\phi
      cem_drift_schottky_d_bc=4.*tmp2**2*ni(i0)*exp(sqrt(tmp1))
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_generation(out)
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     out(lpts4)
      real     tmp1,tmp2,tmp3,tmp4,tmp5,alpha_n,alpha_p
      tmp1= 2.25*10
      tmp2=-3.21*10**6*scale_x/scale_potent
      tmp3= 3.8
      tmp4=-1.75*10**6*scale_x/scale_potent
      do i=1,npts
         if (abs(cee(i,1)) .lt. 1.e-8) then
            out(i)=0
         else
         alpha_n=tmp1*exp(tmp2/abs(cee(i,1)))
         alpha_p=tmp3*exp(tmp4/abs(cee(i,1)))
         tmp5   =alpha_n*abs(cNJ(i,1))+alpha_p*abs(cPJ(i,1))
         out(i) =tmp5*scale_x*1.e+6
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_photogeneration(out)
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     out(lpts4),glmin,xmin,xx
      real     tmp0,tmp1,gamma0,alpha0
      alpha0=2.e7
      gamma0=4.31e21
      tmp0=-alpha0*scale_x
      tmp1=alpha0*gamma0*scale_time/scale_concent
      xmin=glmin(xm1,npts)
      do i=1,npts
         xx=(xm1(i,1,1,1)-xmin)
         out(i)=tmp1*exp(tmp0*xx)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_dissociation(out)
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     out(lpts4)
      real     gamma_r
      real     tmp0,tmp1,tmp2,tmp3,zr,zi,z1r,z1i
      complex  z,z1
      !out=D(cee,ce)
      tmp0=e_charge/epsilon_0*scale_mu
      tmp0=tmp0*3.0/(4.0*pi)/scale_x**3
      gamma_r=tmp0
      do i=1,npts
         tmp0=gamma_r*(mu_n(i)+mu_p(i))/d_permit(i)
         tmp2=tmp0*cE(i)*exp(-0.5*constant_0)*scale_concent
         tmp3=-0.5*constant_0*r_coulomb
     $        *abs(cEE(i,1))*scale_potent/scale_x !-b(|E|)
         call zsqrt(2.0*tmp3,0,zr,zi)!sqrt(-2b(|E|))
         z  =2.0*cmplx(zr,zi) !2sqrt(-2b(|E|))
         call CBESSJ(z,1, z1) !z1=J_1(z)
         z1r=real(z1)
         z1i=imag(z1)
         tmp1=z1i/zi
         out(i)=tmp2*tmp1*scale_time/scale_concent
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_mobility_update
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real     nd,tmp0,tmp1,tmp2,tmp3,tmp4,glamax
      tmp0=5.E-4*sqrt(scale_elect)
      do i=1,npts
         tmp1=tmp0*sqrt(abs(cee(i,1)))
         tmp2=exp(tmp1)
      !for electrons
         mu_n(i)=mu_n0*tmp2
      !for holes
         mu_p(i)=mu_p0*tmp2
      enddo
c      write(6,*) 'yyy',istep,glamax(cee(1,1),npts)
c      write(6,*) 'yyy',istep,mu_n(1)
c      write(6,*) 'yyy',istep,mu_p(1)
      return
      end
c-----------------------------------------------------------------------
c...  Solve electric potential (potent): direct or iterative methods
c... _required params: param(22,23,116,117,118)
      subroutine cem_drift_poisson
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT' 
      include 'POISSON'
      include 'ZPER'   
      include 'BCS'   
  
      integer icalled
      save    icalled
      data    icalled/0/
      real    wk1(lpts),wk2(lpts)
      integer i,isd,ifld,imsh,maxit    
      real    tol
      real*8  dclock
      real btime2,btime1,glamax
      real  xx,yy,zz
      if (icalled.eq.0) then

c...  setup variable for poisson solve: [ h1*A + h2*B ] phi = B*rhs
          do ifld=1,nfield
          call rone    (dmask(1,ifld),npts)          
          call dd_setdmaskp(dmask(1,ifld),npts,ifld)     
          enddo
          call rone    (mult,npts)
          call dssum   (mult,nx1,ny1,nz1)
          call invcol1 (mult,npts)
          call rone    (h1,npts) 
          if (.not.ifsol) call col2(h1,d_permit,npts)
          call rzero   (h2,npts)

c... setup variables for fast diagonalization: [A - kwave2*B] phi= B*rhs
          if (IFGFDMDD) call gfdm_init(nx1,ny1,nz1,.false.,0.0)
          icalled = 1
      endif
      if(ifpec.and.ifinhom) then
      call axhelm3(wk1,potent_b,h1,h2,1,1) 
      call add2s2(rhs,wk1,-1.0,npts)
      endif 
c...  direct methods   : param(23)<0
      if (IFDIRECT) then       

        if   (IFGFDMDD)    then   ! param(23)= -1
        call gfdm_pres_solv(potent,rhs,wk1,wk2,.false.,0.0)
        endif

c...  iterative methods: param(23)>=0 (0: default, positive#: precondition) 
      else          

        isd   = 1
        imsh  = 1
        maxit = 2000
        tol   = param(22) 

        if     (IFCG)    then   ! param(22)= 1
        call dssum (rhs,nx1,ny1,nz1)
        call col2  (rhs,dmask(1,1),npts)
        call cggo2 (potent,rhs,h1,h2,mult,dmask(1,1)
     $             ,imsh,isd,npts,tol,maxit)

        elseif (IFgmres) then   ! param(22)= 0
        call dssum (rhs,nx1,ny1,nz1)
        call col2  (rhs,dmask(1,1),npts)
        call drift_hmh_gmres
     $             (potent,rhs,h1,h2,mult,dmask(1,1),isd,imsh,npts,tol)
        endif

      endif

      if (ifpec .and.ifinhom) call add2s2(potent,potent_b,1.0,npts)
c...  compute electric fields from electric potential: cEE=\grad (potent) 
      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent)          
      call chsign(cEE(1,1),npts)
      call chsign(cEE(1,2),npts)
      call chsign(cEE(1,3),npts)
      return
      end
c----------------------------------------------------------------------
      subroutine cem_drift_grad(w1,w2,w3,u1)
c----------------------------------------------------------------------
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
             w3(k) =0.0
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      real     wk1(lpts),wk2(lpts),wk3(lpts)
      integer  idummy
      real     glamax
      if (.not.ifsrc) return  ! param(6)=1 to turn on

      call usersrc (idummy,resN,resP,resE,rhs,wk2,wk3)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs2
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
      real wkn(lpts),wkp(lpts),wke(lpts)
      real wk1(lpts),wk2(lpts),wk3(lpts)
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      tt=time+1.0
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call rone (h1n,npts)
      call rzero(h2n,npts)
      if (ifdealias) then
      call drift_intp_u_grad_v(wkn,cN,h1n)
      call drift_intp_u_grad_v(wkp,cP,h1n)
      else
      call axhelm2(wkn,cN,h1n,h2n,1,1)       ! w = A x
      call axhelm2(wkp,cP,h1n,h2n,1,1)       ! w = A x
      endif
      
      do i=1, npts
c below define rhs for potential, which should be given in source FIXME
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(pi*xx)**2*sin(tt)**2   ! u
         cc=-2.0*pi**2*sin(tt)**2*cos(2*pi*xx)
         ee=cc-charge(i)*aa
         rhs(i)=cP(i)-cN(i)
     $          +dt*(wk2(i)-wk1(i))
         rhs(i)=charge(i)*rhs(i)
         rhs(i)=rhs(i)+ee
      enddo       
      call col2(rhs,bm1,npts)
      call add2s2(rhs,wkn, 1.0*dt,npts)
      call add2s2(rhs,wkp,-1.0*dt,npts)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source_semi    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      integer i
      real tt,xx,yy,aa,bb,cc,a1,a2,dd !delete later fixme
c      if (.not.ifsrc) return
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call copy    (resN,cN,npts)
      call add2s2  (resN,wk1,1.*dt,npts) 
      return
      end

c=======================================================================
      subroutine cem_drift_sem_bdf1              
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'  
      integer  i,isd,imsh,maxit
      real     wk1(lpts),wk2(lpts), tol,glamax
      real     diff_n_dt(lpts), diff_p_dt(lpts),diff_e_dt(lpts)
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/

      if (icalled.eq.0) then

c...    elliptic: 2nd + constant+ 1st
        do i=1,npts
        diff_n_dt(i)= diff_n(i)*dt
        diff_p_dt(i)= diff_p(i)*dt
        diff_e_dt(i)= diff_e(i)*dt
        enddo
        call rone    (h2n,npts)           ! h2n= 1.0
        call rone    (h2p,npts)           ! h2p= 1.0
        call rone    (h2e,npts)           ! h2e= 1.0
        call invers2 (h1n,charge,npts)    ! h1n= 1./charge
        call col2   (h1n,diff_n_dt,npts) ! h1n= h1n * (diff_n*dt)
        call invers2 (h1p,charge,npts)    ! h1p= 1./charge
        call col2   (h1p,diff_p_dt,npts) ! h1p= h1p * (diff_p*dt)
        call invers2 (h1e,charge,npts)    ! h1p= 1./charge
        call col2   (h1e,diff_e_dt,npts) ! h1p= h1p * (diff_p*dt)

        icalled = 1

       endif 

c...  direct methods: param(23)<0
      if (IFDIRECT)    then 
         
         if (IFGFDMDD)   then 

         call cem_drift_rhs_bdf1       

         do i=1,npts
         diff_n_dt(i)= diff_n(i)*dt
         diff_p_dt(i)= diff_p(i)*dt
         dtbdi_n(i)  = 1.0/diff_n_dt(i)
         dtbdi_p(i)  = 1.0/diff_p_dt(i)
         enddo
  
         call col2 (resN,dtbdi_n,npts)  ! scaling rhs: rhsN = rhsN/(diff_n_dt)
         call col2 (resP,dtbdi_p,npts)  ! scaling rhs: rhsP = rhsP/(diff_p_dt)

         call gfdm_pres_solv(cN,resN,wk1,wk2,.false.,-dtbdi_n(1)) !FIXME kwave2= 1/diff_n_dt 
         call gfdm_pres_solv(cP,resP,wk1,wk2,.false.,-dtbdi_p(1)) !FIXME kwave2= 1/diff_p_dt
         if (ifexct) then
         diff_e_dt(i)= diff_e(i)*dt
         dtbdi_e(i)  = 1.0/diff_e_dt(i)
         call col2 (resE,dtbdi_e,npts)  ! scaling rhs: rhsP = rhsP/(diff_p_dt)
         call gfdm_pres_solv(cE,resE,wk1,wk2,.false.,-dtbdi_e(1)) !FIXME kwave2= 1/diff_p_dt
         endif
         endif

c...  iterative methods: param(23)>=0 (0: default, positive#: precondition) 
      else 

         isd   = 1 
         imsh  = 1
         tol   = param(22)
         maxit = 2000

         if     (IFCG)    then 

         call cem_drift_rhs_bdf1   

         call dssum (resN,nx1,ny1,nz1)
         call dssum (resP,nx1,ny1,nz1)
         call col2  (resN,dmask(1,2),npts)  
         call col2  (resP,dmask(1,3),npts)  
         call cggo2 (cN,resN,h1n,h2n,mult,dmask(1,2)
     $              ,imsh,isd,npts,tol,maxit)
         call cggo2 (cP,resP,h1p,h2p,mult,dmask(1,3)
     $              ,imsh,isd,npts,tol,maxit)
         if (IFEXCT) then
         call dssum (resE,nx1,ny1,nz1)
         call col2  (resE,dmask(1,4),npts)  
         call cggo2 (cE,resE,h1e,h2e,mult,dmask(1,4)
     $              ,imsh,isd,npts,tol,maxit)
         endif
         elseif (IFGMRES) then

         call cem_drift_rhs_bdf1   

         call dssum (resN,nx1,ny1,nz1)
         call dssum (resP,nx1,ny1,nz1)
         call col2  (resN,dmask(1,2),npts)  
         call col2  (resP,dmask(1,3),npts)  

         call drift_hmh_gmres
     $              (cN,resN,h1n,h2n,mult,dmask(1,2),isd,imsh,npts,tol)
         call drift_hmh_gmres
     $              (cP,resP,h1p,h2p,mult,dmask(1,3),isd,imsh,npts,tol)
         if (ifexct) then
         call dssum (resE,nx1,ny1,nz1)
         call col2  (resE,dmask(1,4),npts)  
         call drift_hmh_gmres
     $              (cE,resE,h1e,h2e,mult,dmask(1,4),isd,imsh,npts,tol)
         endif
         endif
      endif  
      if (ifpec .and.ifinhom) then
         call add2s2(cN,cn_b,1.0,npts)
         call add2s2(cP,cp_b,1.0,npts)
         if (ifexct) call add2s2(cE,ce_b,1.0,npts)
      endif
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_current_density    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      real    wk1(npts,3)
      call cem_drift_grad(cNJ(1,1),cNJ(1,2),cNJ(1,3),CN)
      call copy(wk1,cEE(1,1),npts)
      call col2(wk1,cN,npts)
      call col2(wk1,mu_n,npts)
      call col2(cNJ(1,1),mu_n,npts)
      call add2s2(cNJ(1,1),wk1,-1.0,npts)

      call cem_drift_grad(cPJ(1,1),cPJ(1,2),cPJ(1,3),CP)
      call col2(cPJ(1,1),mu_p,npts)
      call col2(wk1,cP,npts)
      call col2(wk1,mu_p,npts)
      call col2(cPJ(1,1),mu_p,npts)
      call add2s2(cPJ(1,1),wk1, 1.0,npts)
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_rhs_bdf1    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'
      real    wk_cN (lpts), wk_cP (lpts),wk_cE(lpts)
      real    wk_h1n(lpts), wk_h1p(lpts),glamax
      integer i 
      call copy   (wk_cN,cN,npts)
      call copy   (wk_cP,cP,npts)
      call col2   (wk_cN,bm1,npts)  
      call col2   (wk_cP,bm1,npts)  
      call add2s1 (resN,wk_cN,dt,npts) !resN=cN+dt*(resN) with weights
      call add2s1 (resP,wk_cP,dt,npts) !resP=cP+dt*(resP) with weights

c...  resN=resN + wk_cN: (wk_cN= -cn\grad\Phi,\grad v), h1= mu_n*cN, h2=0
c...  resP=resP - wk_cP: (wk_cP=  cp\grad\Phi,\grad v), h1= mu_p*cP, h2=0
c...  compute the nonlinear terms: cN*(grad(phi),cP*(grad(phi)
c...  and add it to the right-hand-side  
      call copy(wk_h1n,cN,npts)
      call col2(wk_h1n,mu_n,npts) ! wk_h1n= mu_n * cN  (=h1)      
      call copy(wk_h1p,cP,npts) 
      call col2(wk_h1p,mu_p,npts) ! wk_h1p= mu_p * cP  (=h1)
      call axhelm3(wk_cN,potent,wk_h1n,h2,1,1) 
      call axhelm3(wk_cP,potent,wk_h1p,h2,1,1)

      call add2s2 (resN,wk_cN, 1.0*dt,npts) 
      call add2s2 (resP,wk_cP,-1.0*dt,npts)
      
      if(ifpec.and.ifinhom) then
      call axhelm3(wk_cN,cN_b,h1n,h2n,1,1) 
      call add2s2(resN,wk_cN,-1.0,npts)
      call axhelm3(wk_cP,cP_b,h1p,h2p,1,1) 
      call add2s2(resP,wk_cP,-1.0,npts)
      if(.not.ifsol) then
      call copy(wk_h1n,cN_b,npts)
      call col2(wk_h1n,mu_n,npts) ! wk_h1n= mu_n * cN  (=h1)      
      call copy(wk_h1p,cP_b,npts) 
      call col2(wk_h1p,mu_p,npts) ! wk_h1p= mu_p * cP  (=h1)
      call axhelm3(wk_cN,potent,wk_h1n,h2,1,1) 
      call axhelm3(wk_cP,potent,wk_h1p,h2,1,1)
      call add2s2 (resN,wk_cN,-1.0*dt,npts) 
      call add2s2 (resP,wk_cP,+1.0*dt,npts)
      endif
      endif 
      if (ifexct) then
      call copy   (wk_cE,cE,npts)
      call col2   (wk_cE,bm1,npts)  
      call add2s1 (resE,wk_cE,dt,npts) !resE=cE+dt*(resE) with weights
      if(ifpec.and.ifinhom) then
      call axhelm3(wk_cE,cE_b,h1e,h2e,1,1) 
      call add2s2(resE,wk_cE,-1.0,npts)
      endif 
      endif 
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_sem_bdf2               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      integer  i,isd,imsh,maxit
      real     wk1(lpts),wk2(lpts),dummy1(lpts)
      real     tol,diff_n_dt(lpts), diff_p_dt(lpts),diff_e_dt(lpts)
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/

      if (icalled.eq.0) then

c...    elliptic: 2nd + constant+ 1st
        do i=1,npts
        diff_n_dt(i)= diff_n(i)*dt
        diff_p_dt(i)= diff_p(i)*dt
        enddo
        call rone    (h2n,npts)
        call rone    (h2p,npts)
        call cmult   (h2n,1.5,npts)       ! h2n= 1.5
        call cmult   (h2p,1.5,npts)       ! h2p= 1.5
        call invers2 (h1n,charge,npts)    ! h1n= 1./charge
        call col2    (h1n,diff_n_dt,npts) ! h1n= 1./charge*(diff_n*dt)
        call invers2 (h1p,charge,npts)    ! h1p= 1./charge
        call col2    (h1p,diff_p_dt,npts) ! h1p= 1./charge*(diff_p*dt)
        if (ifexct) then
        do i=1,npts
        diff_e_dt(i)= diff_e(i)*dt
        enddo
        call rone    (h2e,npts)
        call cmult   (h2e,1.5,npts)       ! h2n= 1.5
        call invers2 (h1e,charge,npts)    ! h1n= 1./charge
        call col2    (h1e,diff_e_dt,npts) ! h1n= 1./charge*(diff_n*dt)
        endif
        icalled = 1

      endif

      if (istep.eq.1) then
        if(nid.eq.1) write(6,*) 'initial first two step first'
        call exitt
        return 
      endif     

c...  direct methods: param(23)<0
      if (IFDIRECT) then 

         if  (IFGFDMDD)     then 
         call cem_drift_rhs_bdf2
         do i=1,npts
         diff_n_dt(i)= diff_n(i)*dt
         diff_p_dt(i)= diff_p(i)*dt
         dtbdi_n(i)  = 1.0/diff_n_dt(i)
         dtbdi_p(i)  = 1.0/diff_p_dt(i)
         enddo
         call col2 (resN,dtbdi_n,npts)  ! scaling rhs: rhsN = rhsN/(diff_n_dt)
         call col2 (resP,dtbdi_p,npts)  ! scaling rhs: rhsP = rhsP/(diff_p_dt)
         do i=1,npts
         dtbdi_n(i)  = 1.5/diff_n_dt(i)
         dtbdi_p(i)  = 1.5/diff_p_dt(i)
         enddo
         call copy(potent1,potent,npts)  ! save potent1 = potent (istep-1)
         call copy(cN1,cN,npts)          ! save cN1= cN(istep-1)
         call copy(cP1,cP,npts)          ! save cP1= cP(istep-1) 
         
         call gfdm_pres_solv(cN,resN,wk1,wk2,.false.,-dtbdi_n(1)) !fixme kwave2= 1/diff_n_dt 
         call gfdm_pres_solv(cP,resP,wk1,wk2,.false.,-dtbdi_p(1)) !fixme kwave2= 1/diff_p_dt
         if (ifexct) then
         do i=1,npts
         diff_e_dt(i)= diff_e(i)*dt
         dtbdi_e(i)  = 1.0/diff_e_dt(i)
         enddo
         call col2 (resE,dtbdi_e,npts)  ! scaling rhs: rhsE = rhsE/(diff_e_dt)
         do i=1,npts
         dtbdi_e(i)  = 1.5/diff_e_dt(i)
         enddo
         call copy(cE1,cE,npts)          ! save cE1= cE(istep-1) 
         call gfdm_pres_solv(cE,resE,wk1,wk2,.false.,-dtbdi_e(1)) !fixme kwave2= 1/diff_e_dt
         endif
         endif

c...  iterative methods: param(23)>=0 (0: default, positive#: precondition) 
      else

         isd   = 1
         imsh  = 1
         tol   = param(22)
         maxit = 2000

         if     (IFCG)    then 

         call cem_drift_rhs_bdf2   

         call dssum (resN,nx1,ny1,nz1)
         call col2  (resN,dmask(1,2),npts)  
         call dssum (resP,nx1,ny1,nz1)
         call col2  (resP,dmask(1,3),npts)  

         call copy(potent1,potent,npts)  ! save potent1 = potent (istep-1)
         call copy(cN1,cN,npts)          ! save cN1= cN(istep-1)
         call copy(cP1,cP,npts)          ! save cP1= cP(istep-1) 
         
         call cggo2 (cN,resN,h1n,h2n,mult,dmask(1,2)
     $              ,imsh,isd,npts,tol,maxit)
         call cggo2 (cP,resP,h1p,h2p,mult,dmask(1,3)
     $              ,imsh,isd,npts,tol,maxit)
         if (ifexct) then
         call dssum (resE,nx1,ny1,nz1)
         call col2  (resE,dmask(1,4),npts)  
         call copy(cE1,cE,npts)          ! save cP1= cP(istep-1) 
         call cggo2 (cE,resE,h1e,h2e,mult,dmask(1,4)
     $              ,imsh,isd,npts,tol,maxit)
         endif
         elseif (IFGMRES) then

         call cem_drift_rhs_bdf2   

         call dssum (resN,nx1,ny1,nz1)
         call col2  (resN,dmask(1,2),npts)  
         call dssum (resP,nx1,ny1,nz1)
         call col2  (resP,dmask(1,3),npts)  

         call copy(potent1,potent,npts)  ! save potent1 = potent (istep-1)
         call copy(cN1,cN,npts)          ! save cN1= cN(istep-1)
         call copy(cP1,cP,npts)          ! save cP1= cP(istep-1) 
         
         call drift_hmh_gmres
     $              (cN,resN,h1n,h2n,mult,dmask(1,2),isd,imsh,npts,tol)
         call drift_hmh_gmres
     $              (cP,resP,h1p,h2p,mult,dmask(1,3),isd,imsh,npts,tol)
         if (ifexct) then
         call dssum (resE,nx1,ny1,nz1)
         call col2  (resE,dmask(1,4),npts)  
         call copy(cE1,cE,npts)          ! save cP1= cP(istep-1) 
         call drift_hmh_gmres
     $              (cE,resE,h1e,h2e,mult,dmask(1,4),isd,imsh,npts,tol)
         endif
         endif

      endif
      if (ifpec .and.ifinhom) then
         call add2s2(cN,cn_b,1.0,npts)
         call add2s2(cP,cp_b,1.0,npts)
         if (ifexct) call add2s2(cE,ce_b,1.0,npts)
      endif

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_rhs_bdf2    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'

      real wk_pot (lpts)
      real wk_cN (lpts), wk_cP (lpts), wk_cE (lpts)
      real wk_cN1(lpts), wk_cP1(lpts)        
      real wk_h1n(lpts), wk_h1p(lpts)                 


      call copy   (wk_cN,cN,npts)
      call cmult  (wk_cN,2.0,npts)
      call add2s2 (wk_cN,cN1,-0.5,npts)
      call col2   (wk_cN,bm1,npts)      !wk_cN=2.0*cN-0.5*cN1 with weights information
      call copy   (wk_cP,cP,npts)
      call cmult  (wk_cP,2.0,npts)
      call add2s2 (wk_cP,cP1,-0.5,npts) !wk_cP=2.0*cP-0.5*cP1 with weights information
      call col2   (wk_cP,bm1,npts)  

      
      call cmult  (resN,dt ,npts)       !resN = resN*dt 
      call add2s2 (resN,wk_cN,1.0,npts) !resN = resN + 2.0*cN-0.5*cN1 

      call cmult  (resP,dt ,npts)       !resP = resP*dt 
      call add2s2 (resP,wk_cP,1.0,npts) !resP = resP + 2.0*cP-0.5*cP1 

c...  nonlinear terms
      call copy   (wk_cN,cN ,npts)
      call cmult  (wk_cN,2.0,npts)      !wk_cN= 2*cN
      call add2s2 (wk_cN,cN1,-1.0,npts) !wk_cN= 2*cN-cN1

      call copy   (wk_cP,cP ,npts)                   
      call cmult  (wk_cP,2.0,npts)      !wk_cP= 2*cP      
      call add2s2 (wk_cP,cP1,-1.0,npts) !wk_cP= 2*cP-cP1

      call copy   (wk_pot,potent,npts)                   
      call cmult  (wk_pot,2.0    ,npts)        !wk_pot= 2*potent  
      call add2s2 (wk_pot,potent1,-1.0,npts)   !wk_pot= 2*potent-potent1
    
      call copy   (wk_h1n,wk_cN,npts)
      call col2   (wk_h1n,mu_n,npts)     ! wk_h1n= mu_n * cN  (=h1)      
      call copy   (wk_h1p,wk_cP,npts) 
      call col2   (wk_h1p,mu_p,npts)     ! wk_h1p= mu_p * cP  (=h1)
      call axhelm3(wk_cN1,wk_pot,wk_h1n,h2,1,1)!wk_cN1=(cn\grad\Phi,\grad v)
      call axhelm3(wk_cP1,wk_pot,wk_h1p,h2,1,1)!wk_cP1=(cp\grad\Phi,\grad v)

c      call add2s2 (resN,wk_cN1, dt,npts)
c      call add2s2 (resP,wk_cP1,-dt,npts)
     
c...  nonlinear terms
      call axhelm3(wk_cN,potent,cN,h2,1,1)!wk_cN1=(cn\grad\Phi,\grad v)
      call axhelm3(wk_cP,potent,cP,h2,1,1)!wk_cN1=(cn\grad\Phi,\grad v)
      call axhelm3(wk_cN1,potent1,cN1,h2,1,1)!wk_cN1=(cn\grad\Phi,\grad v)
      call axhelm3(wk_cP1,potent1,cP1,h2,1,1)!wk_cN1=(cn\grad\Phi,\grad v)
      
      call add2s2 (resN,wk_cN, 2.0*dt,npts)
      call add2s2 (resP,wk_cP,-2.0*dt,npts)

      call add2s2 (resN,wk_cN1,-1.0*dt,npts)
      call add2s2 (resP,wk_cP1,+1.0*dt,npts)
     
      if(ifpec.and.ifinhom) then
      call axhelm3(wk_cN,cN_b,h1n,h2n,1,1) 
      call add2s2(resN,wk_cN,-1.0,npts)
      call axhelm3(wk_cP,cP_b,h1p,h2p,1,1) 
      call add2s2(resP,wk_cP,-1.0,npts)
      endif 
      if (ifexct) then
      call copy   (wk_cE,cE,npts)
      call cmult  (wk_cE,2.0,npts)
      call add2s2 (wk_cE,cE1,-0.5,npts) !wk_cE=2.0*cE-0.5*cE1 with weights information
      call col2   (wk_cE,bm1,npts)  
      call cmult  (resE,dt ,npts)       !resE = resE*dt 
      call add2s2 (resE,wk_cE,1.0,npts) !resE = resE + 2.0*cE-0.5*cE1 
      if(ifpec.and.ifinhom) then
      call axhelm3(wk_cE,cE_b,h1e,h2e,1,1) 
      call add2s2(resE,wk_cE,-1.0,npts)
      endif
      endif 
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_div(wk0,wkx,wky,wkz)   
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      real    wkx(1),wky(1),wkz(1),wk0(1)

      if (if3d) then
      call drift_grad(wk1,wk2,wk3,wkx)
      call copy(wk0,wk1,npts)
      call drift_grad(wk1,wk2,wk3,wky)
      call add2s2(wk0,wk2,1.0,npts)
      call drift_grad(wk1,wk2,wk3,wky)
      call add2s2(wk0,wk3,1.0,npts)
      else
      call drift_grad(wk1,wk2,wk3,wkx)
      call copy(wk0,wk1,npts)
      call drift_grad(wk1,wk2,wk3,wky)
      call add2s2(wk0,wk2,1.0,npts)
      endif
      return
      end
c--------------------------------------------------------------------- 
      subroutine drift_filter(wght)
c--------------------------------------------------------------------- 
c
c     filter vx,vy,vz, and p by simple interpolation
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
c
c
c     These are the dimensions that we interpolate onto for v and p:
      parameter(lxv=lx1-1)
      parameter(lxp=lx2-1)
c
      real intdv(lx1,lx1)
      real intuv(lx1,lx1)
      real intdp(lx1,lx1)
      real intup(lx1,lx1)
      real intv(lx1,lx1)
      real intp(lx1,lx1)
c
      save intdv
      save intuv
      save intdp
      save intup
      save intv
      save intp

      common /ctmp0/ intw,intt
     $             , wk1,wk2
     $             , zgmv,wgtv,zgmp,wgtp,tmax(100),omax(103)

      real intw(lx1,lx1)
      real intt(lx1,lx1)
      real wk1  (lx1,lx1,lx1,lelt)
      real wk2  (lx1,lx1,lx1)
      real zgmv(lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)


      character*18 sfmt

      integer icalld
      save    icalld
      data    icalld /0/

      imax = nid
      imax = iglmax(imax,1)
      jmax = iglmax(imax,1)
      if (icalld.eq.0) then
         icalld = 1
         ncut = param(85) ! 1 or 2   
         call build_new_filter(intv,zgm1,nx1,ncut,wght,nid)
      endif

      call filterq(cN,intv,nx1,nz1,wk1,wk2,intt,if3d,exmx)
      call filterq(cP,intv,nx1,nz1,wk1,wk2,intt,if3d,eymx)
      call filterq(cE,intv,nx1,nz1,wk1,wk2,intt,if3d,ezmx)


c     write(6,1) istep,time,exmx,eymx,ezmx,hxmx,hymx,hzmx
c   1 format(i8,1p7e10.3,' filt')


      return
      end
c-----------------------------------------------------------------------
      subroutine drift_wght_div_d(w0,w1,w2,w3,u1,u2,u3) 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE' 
      include 'TOTAL' 
c     include 'EMWAVE' 
      include 'DEALIAS' 
c     ddddd
c     Combined weighted div:   w0 = B*div(u1,u2,u3) on a finer mesh
c     with dealiasing
c
      real     w0(1)                  
      real     u1(1),u2(1),u3(1)
      real     w1(1),w2(1),w3(1)
                                    
      integer    ldd,md
      parameter (ldd = lxd**ldim  )
      parameter ( md = lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)
 
      real         ud1r,ud1s,ud1t
     $           , ud2r,ud2s,ud2t
     $           , ud3r,ud3s,ud3t

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw

      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd0(ldd),wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd0,wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)                     
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)
      real            dxmd,dxtmd,w3md,wgld,zmd

      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy, nn, mm, mx
      save    icalld 
      data    icalld /0/
 
      mx    = nx1
      nn    = nx1-1   
      mm    = md -1
      nxyd  = md**ndim
      nxy   = mx**ndim

      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

          if (if3d) then
             l = 0
             do k=1,md
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)*wgld(k)
             enddo
             enddo
             enddo
          else
             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)              
             enddo
             enddo
          endif
                                           
          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)

      endif

      do e=1,nelt

         j  = nxyd *(e-1) 
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u1(j1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,u2(j1),mx,md,if3d, 0)
         call intp_rstd (ju3,u3(j1),mx,md,if3d, 0)

         if (if3d) then

            call local_grad3 (ud1r,ud1s,ud1t,ju1,mm,1,dxmd,dxtmd)
            call local_grad3 (ud2r,ud2s,ud2t,ju2,mm,1,dxmd,dxtmd)
            call local_grad3 (ud3r,ud3s,ud3t,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j
 
            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u1tw = ud1t(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)
            u2tw = ud2t(i)*w3md(i)
            u3rw = ud3r(i)*w3md(i)
            u3sw = ud3s(i)*w3md(i)
            u3tw = ud3t(i)*w3md(i)

            wd1(i) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k)
     $               + u1tw*txmd(k))

            wd2(i) = ( u2rw*rymd(k)
     $               + u2sw*symd(k)
     $               + u2tw*tymd(k))

            wd3(i) = ( u3rw*rzmd(k)
     $               + u3sw*szmd(k)
     $               + u3tw*tzmd(k))

            wd0(i) = wd1(i)+wd2(i)+wd3(i)

            enddo

        else

            call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
            call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
            call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j

            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)

            wd1(i) = (  u1rw*rxmd(k)
     $               + u1sw*sxmd(k))

            wd2(i) = (  u2rw*rymd(k)
     $               + u2sw*symd(k))

            wd0(i) = wd1(i)+wd2(i)      

           enddo

         endif

         call intp_rstd (w0(j1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_div_grad_d(w0,w1,w2,w3,u1,u2,u3,u4) 
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE' 
      include 'TOTAL' 
c     include 'EMWAVE' 
      include 'DEALIAS' 
c     ddddd
c     Combined weighted div(u4*(u1,u2,u3)^t):   w0 = B*div(u4*(u1,u2,u3)^t) on a finer mesh
c     with dealiasing
c
      real     w0(1)                  
      real     u1(1),u2(1),u3(1),u4(1)
      real     w1(1),w2(1),w3(1)
                                    
      integer    ldd,md
      parameter (ldd = lxd**ldim  )
      parameter ( md = lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)
 
      real         ud1r,ud1s,ud1t
     $           , ud2r,ud2s,ud2t
     $           , ud3r,ud3s,ud3t

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      real    ju4(ldd)
      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd0(ldd),wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd0,wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)                     
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)
      real            dxmd,dxtmd,w3md,wgld,zmd

      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy, nn, mm, mx
      save    icalld 
      data    icalld /0/
 
      mx    = nx1
      nn    = nx1-1   
      mm    = md -1
      nxyd  = md**ndim
      nxy   = mx**ndim

      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

          if (if3d) then
             l = 0
             do k=1,md
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)*wgld(k)
             enddo
             enddo
             enddo
          else
             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)              
             enddo
             enddo
          endif
                                           
          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)

      endif

      do e=1,nelt

         j  = nxyd *(e-1) 
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u1(j1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,u2(j1),mx,md,if3d, 0)
         call intp_rstd (ju3,u3(j1),mx,md,if3d, 0)
         call intp_rstd (ju4,u4(j1),mx,md,if3d, 0)
         call col2(ju1,ju4,nxyd)
         call col2(ju2,ju4,nxyd)
         call col2(ju3,ju4,nxyd)

         if (if3d) then

            call local_grad3 (ud1r,ud1s,ud1t,ju1,mm,1,dxmd,dxtmd)
            call local_grad3 (ud2r,ud2s,ud2t,ju2,mm,1,dxmd,dxtmd)
            call local_grad3 (ud3r,ud3s,ud3t,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j
 
            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u1tw = ud1t(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)
            u2tw = ud2t(i)*w3md(i)
            u3rw = ud3r(i)*w3md(i)
            u3sw = ud3s(i)*w3md(i)
            u3tw = ud3t(i)*w3md(i)

            wd1(i) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k)
     $               + u1tw*txmd(k))

            wd2(i) = ( u2rw*rymd(k)
     $               + u2sw*symd(k)
     $               + u2tw*tymd(k))

            wd3(i) = ( u3rw*rzmd(k)
     $               + u3sw*szmd(k)
     $               + u3tw*tzmd(k))

            wd0(i) = wd1(i)+wd2(i)+wd3(i)

            enddo

        else

            call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
            call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
            call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j

            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)

            wd1(i) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k))

            wd2(i) = ( u2rw*rymd(k)
     $               + u2sw*symd(k))

            wd0(i) = wd1(i)+wd2(i)      

           enddo

         endif

         call intp_rstd (w0(j1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_intp_u_grad_v(w0,u,v) 
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE' 
      include 'TOTAL' 
c     include 'EMWAVE' 
      include 'DEALIAS' 
c     ddddd
c     Combined weighted div(u4*(u1,u2,u3)^t):   w0 = B*div(u4*(u1,u2,u3)^t) on a finer mesh
c     with dealiasing
c
      real     w0(lx1,ly1,lz1,1)                  
      real     u(lx1,ly1,lz1,1)                  
      real     v(lx1,ly1,lz1,1)                  
                                    
      integer    ldd,md
      parameter (ldd = lxd**ldim  )
      parameter ( md = lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)
 
      real         ud1r,ud1s,ud1t
     $           , ud2r,ud2s,ud2t
     $           , ud3r,ud3s,ud3t

      common /t1t2/dudrd(lxd,lyd,lzd)
     $            ,dudsd(lxd,lyd,lzd)
     $            ,tmp1d(lxd,lyd,lzd)
     $            ,tmp2d(lxd,lyd,lzd)
 
      real         dudrd,dudsd,tmp1d,tmp2d
      real         tm1(lxd,lyd,lzd)
      real         tm2(lxd,lyd,lzd)
      
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      real    ju4(ldd)
      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd0(ldd),wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd0,wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)                     
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)
      real            dxmd,dxtmd,w3md,wgld,zmd
      common /ddtmp6/ g1m1d(lxd**ldim*lelt)
     $               ,g2m1d(lxd**ldim*lelt)
     $               ,g4m1d(lxd**ldim*lelt)
     $               ,jacmd(lxd**ldim*lelt)
     $               ,wjd(lxd**ldim*lelt)
      real            g1m1d,g2m1d,g4m1d,jacmd,wjd
      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy, nn, mm, mx,nptsd
      save    icalld 
      data    icalld /0/
      EQUIVALENCE    (DUDRd,TM1),(DUDSd,TM2)
      mx    = nx1
      nn    = nx1-1   
      mm    = md -1
      nxyd  = md**ndim
      nxy   = mx**ndim
      nptsd = nxyd*nelt
      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)              
             enddo
             enddo
                                           
          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (jacmd(jd),jacm1(j1,1,1,1),mx,md,if3d, 0)
          CALL INVERS2 (WJd(jd),JACMd(jd),nxyd)
         CALL VDOT2(G1M1d(jd),RXMd(jd),RYMd(jd),RXMd(jd),RYMd(jd),nxyd)
         CALL VDOT2(G2M1d(jd),sXMd(jd),sYMd(jd),sXMd(jd),sYMd(jd),nxyd)
         CALL VDOT2(G4M1d(jd),RXMd(jd),RYMd(jd),SXMd(jd),SYMd(jd),nxyd)
          CALL COL2  (G1M1d(jd),WJd(jd),nxyd)
          CALL COL2  (G2M1d(jd),WJd(jd),nxyd)
          CALL COL2  (G4M1d(jd),WJd(jd),nxyd)
          CALL COL2 (G1M1d(jd),W3Md,nxyd)
          CALL COL2 (G2M1d(jd),W3Md,nxyd)
          CALL COL2 (G4M1d(jd),W3Md,nxyd)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)
      endif
         
      
      call rzero(w0,npts)
      do e=1,nelt
         call rzero(wd0,nxyd)
         j  = nxyd *(e-1) + 1 
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u(j1,1,1,1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,v(j1,1,1,1),mx,md,if3d, 0)
         
         CALL MXM  (DXMd,md,ju1,md,DUDRd,md)
         CALL MXM  (ju1,md,DxTMd,md,DUDSd,md)
         CALL COL3 (TMP1d,DUDRd,G1M1d(j),nxyd)
         CALL COL3 (TMP2d,DUDSd,G2M1d(j),nxyd)

         CALL ADDCOL3 (TMP1d,DUDSd,G4M1d(j),NXYd)
         CALL ADDCOL3 (TMP2d,DUDRd,G4M1d(j),NXYd)
         CALL COL2 (TMP1d,ju2,NXYd)
         CALL COL2 (TMP2d,ju2,NXYd)
         CALL MXM  (DXTMd,md,TMP1d,md,TM1,md)
         CALL MXM  (TMP2d,md,DxMd,md,TM2,md)
         CALL ADD2 (wd0,TM1,NXYd)
         CALL ADD2 (wd0,TM2,NXYd)


         call intp_rstd (w0(j1,1,1,1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c-----------------------------------------------------------------------
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
            call axhelm3(w,x,h1,h2,imsh,isd)       ! w = A x
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

            call axhelm3(w,v(1,j),h1,h2,imsh,isd)  ! w = A v
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
            if (nid.eq.0) 
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
      subroutine cem_drift_sem_bdf2_124               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,isd,ifld,imsh,maxit
      real tol,tt
      integer  icalled
      save     icalled
      data     icalled/0/

      if (icalled.eq.0) then
         do ifld=1,nfield
         call rone    (dmask(1,ifld),npts )      ! for grids                         
         call dd_setdmaskp(dmask(1,ifld),npts,ifld)
         enddo
         call rone    (mult,npts)
         call dssum   (mult,nx1,ny1,nz1)
         call invcol1 (mult,npts)
         icalled=1
      endif
c.... I just keep this temporarily
      if (istep.eq.1) then
        tt=time+dt
        call usersol(tt, scn, scp, sce, scn, scp, sce)
        call copy(cN1,cN,npts)
        call copy(cP1,cP,npts)
        call copy(potent1,potent,npts)
        call copy(cN,scN,npts)
        call copy(cP,scP,npts)
        return 
      endif     
 
c---- rhs=rhs*dt+ u^n
      call cem_drift_source_bdf2_124   
      call dssum   (resN,nx1,ny1,nz1)
      call col2    (resN,dmask(1,2),npts)  
      call dssum   (resP,nx1,ny1,nz1)
      call col2    (resP,dmask(1,3),npts)  

c...  elliptic: 2nd + constant+ 1st
      do i=1,npts
         h1n(i)=diff_n(i)/charge(i)*dt
         h2n(i)=1.5 ! mu_n*\delta\phi*dt+1.0
      enddo
      do i=1,npts
         h1p(i)=diff_p(i)/charge(i)*dt
         h2p(i)=1.5 !- mu_p*\delta\phi
      enddo

      call copy(cN1,cN,npts)
      call copy(cP1,cP,npts)
      call copy(potent1,potent,npts)

      isd   = 1
      imsh  = 1
      tol   = param(22)
      maxit=2000

      if (IFCG) then
         call cggo2 (cN,resN,h1n,h2n,mult,dmask(1,2)
     $              ,imsh,isd,npts,tol,maxit)
         call cggo2 (cP,resP,h1p,h2p,mult,dmask(1,3)
     $              ,imsh,isd,npts,tol,maxit)
      endif
      if (IFGMRES) then
         call drift_hmh_gmres(cN,resN,h1n,h2n
     $                   ,mult,dmask(1,2),isd,imsh,npts,tol)
         call drift_hmh_gmres(cP,resP,h1p,h2p
     $                   ,mult,dmask(1,3),isd,imsh,npts,tol)
         if (nid.eq.0) write(6,*) 'done: elliptic_gmres3'
      endif
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_drift_source_bdf2_124    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk0(lpts*3),wk1(lpts),wk2(lpts),wk3(lpts)
      integer i,npts3
      real tt,xx,yy,aa,bb,cc,a1,a2,dd !delete later fixme

      if (param(17).eq.-3) then
      call rzero(resN,npts)
      call add2s2(resN,cN,2.0,npts)
      call add2s2(resN,cN1,-0.5,npts)
      call rzero(resP,npts)
      call add2s2(resP,cP,2.0,npts)
      call add2s2(resP,cP1,-0.5,npts)
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call add2s2 (resN,wk1,1.0*dt,npts)
      call add2s2 (resP,wk2,1.0*dt,npts)
      call col2(resN,bm1,npts)
      call col2(resP,bm1,npts)

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call rzero  (wk1 ,npts)
      call add2s2 (wk1 ,cN , 2.0,npts)
      call add2s2 (wk1 ,cN1,-1.0,npts)
      call rzero  (wk2 ,npts)
      call add2s2 (wk2 ,cP , 2.0,npts)
      call add2s2 (wk2 ,cP1,-1.0,npts)
      call rzero  (wk3 ,npts)
      call add2s2 (wk3 ,potent  , 1.0,npts)
      call axhelm3(wk0 ,wk3,wk1,h2n,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0, 1.0*dt,npts)
      call axhelm3(wk0 ,wk3,wk2,h2p,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0*dt,npts)
      elseif (param(17).eq.-4) then
      call rzero(resN,npts)
      call add2s2(resN,cN,2.0,npts)
      call add2s2(resN,cN1,-0.5,npts)
      call rzero(resP,npts)
      call add2s2(resP,cP,2.0,npts)
      call add2s2(resP,cP1,-0.5,npts)
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call add2s2 (resN,wk1,1.0*dt,npts)
      call add2s2 (resP,wk2,1.0*dt,npts)
      call col2(resN,bm1,npts)
      call col2(resP,bm1,npts)

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call rzero  (wk1 ,npts)
      call add2s2 (wk1 ,cN , 2.0,npts)
      call add2s2 (wk1 ,cN1,-1.0,npts)
      call rzero  (wk2 ,npts)
      call add2s2 (wk2 ,cP , 2.0,npts)
      call add2s2 (wk2 ,cP1,-1.0,npts)
      call rzero  (wk3 ,npts)
      call add2s2 (wk3 ,potent  , 1.0,npts)
      call axhelm3(wk0 ,wk3,wk1,h2n,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0, 1.0*dt,npts)
      call axhelm3(wk0 ,wk3,wk2,h2p,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0*dt,npts)
      elseif (param(17).eq.-5) then
      call rzero(resN,npts)
      call add2s2(resN,cN,2.0,npts)
      call add2s2(resN,cN1,-0.5,npts)
      call rzero(resP,npts)
      call add2s2(resP,cP,2.0,npts)
      call add2s2(resP,cP1,-0.5,npts)
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call add2s2 (resN,wk1,1.0*dt,npts)
      call add2s2 (resP,wk2,1.0*dt,npts)
      call col2(resN,bm1,npts)
      call col2(resP,bm1,npts)

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call rzero  (wk1 ,npts)
      call add2s2 (wk1 ,cN , 2.0,npts)
      call add2s2 (wk1 ,cN1,-1.0,npts)
      call rzero  (wk2 ,npts)
      call add2s2 (wk2 ,cP , 2.0,npts)
      call add2s2 (wk2 ,cP1,-1.0,npts)
      call rzero  (wk3 ,npts)
      call add2s2 (wk3 ,potent  , 2.0,npts)
      call add2s2 (wk3 ,potent1 ,-1.0,npts)
      call axhelm3(wk0 ,wk3,wk1,h2n,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0, 1.0*dt,npts)
      call axhelm3(wk0 ,wk3,wk2,h2p,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0*dt,npts)
      else
      call exitt
      endif
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_poisson_bdf2_24
c----------------------------------------------------------------------
c...  implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT' 
      include 'POISSON'
      integer i
      real    tt !fixme delete finally
      integer icalled
      save    icalled
      data    icalled/0/

c     if (nid.eq.0) write(6,*) '\n call cem_drift_poissonq'
      npts =  nx1*ny1*nz1*nelt

      if (icalled.eq.0) then
          do ifld=1,nfield
          call rone    (dmask(1,ifld),npts )      ! for grids                         
          call dd_setdmaskp(dmask(1,ifld),npts,ifld)
          enddo
          call rone    (mult,npts)
          call dssum   (mult,nx1,ny1,nz1)
          call invcol1 (mult,npts)

          call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
          call rzero   (h2,npts)
          icalled = 1
      endif

      call cem_drift_poisson_rhs24  
      call col2    (rhs,bm1,npts)            ! B*rhs
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,dmask(1,1),npts)

      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = param(22) 

      if (IFCG) then
        call cggo2 (potent,rhs,h1,h2,mult,dmask(1,1)
     $             ,imsh,isd,npts,tol,maxit)
      endif
      if (IFgmres) then
         call drift_hmh_gmres(potent,rhs,h1,h2
     $                   ,mult,dmask(1,1),isd,imsh,npts,tol)
      endif

      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent) !cEE=\grad potent 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs24
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      tt=time+dt+1.0
      if (istep.eq. 0) tt=time+1.0
      do i=1,npts
         rhs(i)=charge(i)*(cP(i)-cN(i))*2.0
         rhs(i)=rhs(i)-charge(i)*(cP1(i)-cN1(i))
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(pi*xx)**2*sin(tt)**2   ! u
         cc=-2.0*pi**2*sin(tt)**2*cos(2*pi*xx)
         ee=cc-charge(i)*aa
         rhs(i)=rhs(i)+ee
      enddo
     
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson1
c----------------------------------------------------------------------
c...  implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT' 
      include 'POISSON'
      include 'ZPER'   
      integer i
      real    tt !fixme delete finally
      integer icalled
      save    icalled
      data    icalled/0/
      real    wk0(lpts),wk1(lpts),wk2(lpts)

c     if (nid.eq.0) write(6,*) '\n call cem_drift_poissonq'
      npts =  nx1*ny1*nz1*nelt

      if (icalled.eq.0) then
         call rone    (dmask(1,1),npts*nfield )      ! for grids                         
         call dd_setdmaskp(dmask(1,1),npts,1)
         call dd_setdmaskp(dmask(1,2),npts,2)
         call dd_setdmaskp(dmask(1,3),npts,3)
         call dd_setdmaskp(dmask(1,4),npts,4)

          call rone    (mult,npts)
          call dssum   (mult,nx1,ny1,nz1)
          call invcol1 (mult,npts)

          call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
          call rzero   (h2,npts)

          if (IFGFDMDD) call gfdm_init(nx1,ny1,nz1,.false.,0.0)      

          icalled = 1
      endif

c      dt1=rktime1-rktime
      dt1=dt
      call cem_drift_poisson_rhs1_rk  
      call usersrc (57,wk0,wk1,wk2,wk0,wk1,wk2)
      call add2s2(rhs,wk0,-1.0*dt1,npts)
      call add2s2(rhs,wk1, 1.0*dt1,npts)
      call col2    (rhs,bm1,npts)            ! B*rhs
      call axhelm3(wk0,cN,h1,h2,1,1)
      call add2s2(rhs,wk0, 1.0*dt1,npts)
      call axhelm3(wk0,cP,h1,h2,1,1)
      call add2s2(rhs,wk0,-1.0*dt1,npts)
      
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,dmask(1,1),npts)
      do i=1,npts
         h1(i)=1.0+charge(i)*dt1*(mu_n(i)*cN(i)+mu_p(i)*cP(i))
      enddo
      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = param(22) 

      if (IFPRECON) then
         call gfdm_pres_solv(potent,rhs,wk1,wk2,.false.,0.0) 
      endif
      
      if (IFCG) then
        call cggo2 (potent,rhs,h1,h2,mult,dmask(1,1)
     $             ,imsh,isd,npts,tol,maxit)
      endif
      if (IFgmres) then
         call drift_hmh_gmres(potent,rhs,h1,h2
     $                   ,mult,dmask(1,1),isd,imsh,npts,tol)
      endif
   
      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent) !cEE=\grad potent 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs1_rk
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      if (ifrk) then
         tt=rktime1+1.0
      endif
      if (ifbdf) then
         tt=time+dt+1.0
         if (istep.eq. 0) tt=time+1.0
      endif
      do i=1,npts
         rhs(i)=charge(i)*(cP(i)-cN(i))
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(pi*xx)**2*sin(tt)**2   ! u
         cc=-2.0*pi**2*sin(tt)**2*cos(2*pi*xx)
         ee=cc-charge(i)*aa
         rhs(i)=rhs(i)+ee
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine fast_poisson_test(mode)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      include 'ZPER'  ! lelx,lely,lelz in nek5000

      common /rwork/ ff(lpts),uu(lpts),wwk(lpts),wwk1(lpts),wwk2(lpts)
      real    ff,uu,wwk,wwk1,wwk2,wmax,glamax
      real    kwave2
      real    pi2,xx,yy
      integer mode,i

      nelx = param(116)+0.001
      nely = param(117)+0.001
      nelz = param(118)+0.001

      if (nelx.gt.lelx) call exitti('lelx too small:$',lelx)
      if (nely.gt.lely) call exitti('lely too small:$',lely)
      if (nelz.gt.lelz) call exitti('lelz too small:$',lelz)


      pi2 = pi/2.
      kwave2 =0.0 !100 !0.0 ! pi*pi

      call gfdm_init(nx1,ny1,nz1,.false.,0.0)

      if (mode.eq.1) then
        do i=1,npts
         xx = xm1(i,1,1,1) !*i
         yy = ym1(i,1,1,1) !*i
         ff(i) = pi*pi2*cos(pi2*xx)*cos(pi2*yy)
     $         - kwave2*cos(pi2*xx)*cos(pi2*yy)
         uu(i) = cos(pi2*xx)*cos(pi2*yy)
         enddo
      else
        do i=1,npts
         xx = xm1(i,1,1,1) !*i
         yy = ym1(i,1,1,1) !*i
         ff(i) = pi*pi2*sin(pi2*xx)*sin(pi2*yy)
     $         - kwave2*sin(pi2*xx)*sin(pi2*yy)
         uu(i) = sin(pi2*xx)*sin(pi2*yy)
         enddo
      endif

      call col2(ff,bm1,npts)

      call gfdm_pres_solv(wwk,ff,wwk1,wwk2,.false.,0.0) ! (A - kwave2*B)z = r

      do i=1,npts
         write(10,10) i,uu(i),wwk(i)
      enddo
  10  format(i8,2e25.15)

      call sub3(wwk1,uu,wwk,npts)

      wmax=glamax(wwk1,npts)

      if (nid.eq.0) write(6,*) 'gdfm:: maxerr',wmax

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_bdf2_init 
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'
      integer kstep,i,j
      real    dt1,glamax
      real    dummy1(lpts)
      
      if (ifsol) then 
         call cem_drift_source           ! param(6)= 1   (to turn on source)
         call cem_drift_poisson          ! solve potent and electric feild cEE
         call copy (potent1,potent,npts) ! save potent1= potent(istep=0)
         call usersol(time+dt, scn, scp, sce, dummy1, dummy1, dummy1)
         call copy (cN1,cN,npts)         ! save cN1= cN (istep=0)
         call copy (cP1,cP,npts)         ! save cP1= cP (istep=0)
         call copy (cN,scN,npts)         ! get  cN at istep=1 
         call copy (cP,scP,npts)         ! get  cP at istep=1
         if (ifexct) then
         call copy (cE1,cE,npts)         ! save cP1= cP (istep=0)
         call copy (cE,scE,npts)         ! get  cN at istep=1 
         endif
       else
         dt1=dt
         dt =dt**2
         kstep=dt1/dt
      
         call copy (cN1,cN,npts)         ! save cN1= cN (istep=0)
         call copy (cP1,cP,npts)         ! save cP1= cP (istep=0)
         if (ifexct) then
         call copy (cE1,cE,npts)         ! save cP1= cP (istep=0)
         endif
      
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson ! solve potent and electric feild cEE
         call cem_drift_sem_bdf1
         time=time+dt
         call copy (potent1,potent,npts) ! save potent1= potent(istep=0)

         do i=2,kstep
         call cem_drift_source  ! param(6)= 1   (to turn on source)
         call cem_drift_poisson ! solve potent and electric feild cEE
         call cem_drift_sem_bdf1
         time=time+dt
         enddo

         time=time-dt1
         dt  =dt1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_neuman_fc_ptr
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BCS'
      integer  ifld,i,i0,i1,j,l1,l2,l3,e,f
      character CB*3

      if (ifflow) then
          i0= 1
          i1= 1
      endif
      if (ifheat) then
          i0= 1
          i1= nfield
      endif

      do ifld= i0,i1
         l1=0
         l2=0
         l3=0
      do e= 1,nelt
      do f= 1,nfaces
      do j= 1,nxzf
         CB= CBC(f,e,ifld+1)
         i =(e-1)*nfaces*nxzf+(f-1)*nxzf+j !global face numbering
         if  (CB.eq.'N  ') then
            l1= l1+1
            cemface_nmn(l1,ifld)= i
         elseif  (CB.eq.'PEC') then
            l2= l2+1
            cemface_pec(l2,ifld)= i
         elseif  (CB.eq.'R  ') then
            l3= l3+1
            cemface_r  (l3,ifld)= i
         endif
      enddo
      enddo
      enddo

      ncemface_nmn(ifld)  = l1
      ncemface_pec(ifld)  = l2
      ncemface_r  (ifld)  = l3
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bd_nmn(bdin,bdout,ifld)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BCS'
     
      real    bdin(1),bdout(1)
      integer ifld,i,i0,j
      
      do i=1,ncemface_nmn(ifld)
         j=cemface_nmn(i,ifld)
         i0=cemface(j)
         bdout(i0)=bdout(i0)+bdin(i)*aream(j) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bd_pec(bdin,bdout,ifld)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BCS'
      real    bdin(1),bdout(1)
      integer ifld,i,i0,j        
      
      do i=1,ncemface_pec(ifld)
         j =cemface_pec(i,ifld)
         i0=cemface(j)
         bdout(i0)=bdout(i0)+bdin(i)*aream(j) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bd_r  (bdin,bdout,ifld)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'BCS'
      real    bdin(1),bdout(1)
      integer ifld,i,i0,j        
      
      do i=1,ncemface_r (ifld) 
         j =cemface_r  (i,ifld)
         i0=cemface(j)
         bdout(i0)=bdout(i0)+bdin(i)*aream(j) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_scaling(flag)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      integer flag
      integer i,j,i0
      if (flag.eq.1) then
         call cmult(cN,scale_concent,npts)      
         call cmult(cP,scale_concent,npts)      
         call cmult(potent,scale_potent,npts)      
         call cmult(cEE,scale_potent,npts*3)      
      else
         
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine dd_setdmaskp(dmask,n,ifld)
c-----------------------------------------------------------------------
C     set dmask for poisson problem
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer   e,f,i,j,n,ifld
      real      dmask(1)
      character CB*3
      do e=1,nelt
      do f=1,nfaces

         CB =CBC(f,e,ifld+1)
c        Homogeneous Dirichlet boundary conditions
         if (CB.eq.'PEC') call facev (dmask(1),e,f,0.0,nx1,ny1,nz1)
         if (CB.eq.'PML') call facev (dmask(1),e,f,0.0,nx1,ny1,nz1) ! fischer/mmin: 1/23/07

      enddo
      enddo

      call dsop(dmask(1),'MUL',nx1,ny1,nz1)

c     if (nid.eq.0) write(6,*) ' done: setdmaskp'

      return
      end
c-----------------------------------------------------------------------
      real*8 Function Fact(K)
c-----------------------------------------------------------------------
      Integer i
      Real*8  f
      F=1.d0
      do i=2, k 
         f=f*dfloat(i)
      end do
      Fact=f
      return
      End
c-----------------------------------------------------------------------
!*******************************************
!*           FUNCTION  GAMMA(X)            *
!* --------------------------------------- *
!* Returns the value of Gamma(x) in double *
!* precision as EXP(LN(GAMMA(X))) for X>0. *
!*******************************************
      real*8 Function Gamma(xx)
c-----------------------------------------------------------------------
      parameter(ONE=1.d0,FPF=5.5d0,HALF=0.5d0)
      real*8 xx
      real*8 cof(6)
      real*8 stp,x,tmp,ser
      integer j
      cof(1)=76.18009173d0
      cof(2)=-86.50532033d0
      cof(3)=24.01409822d0
      cof(4)=-1.231739516d0
      cof(5)=0.120858003d-2
      cof(6)=-0.536382d-5
      stp=2.50662827465d0
  
      x=xx-ONE
      tmp=x+FPF
      tmp=(x+HALF)*LOG(tmp)-tmp
      ser=ONE
      do j=1, 6
         x=x+ONE
         ser=ser+cof(j)/x
      end do
      Gamma = EXP(tmp+LOG(stp*ser))
      return
      End
c-----------------------------------------------------------------------
      Subroutine CBESSJ(z, nu, z1)
!---------------------------------------------------
!                       inf.     (-z^2/4)^k
!   Jnu(z) = (z/2)^nu x Sum  ------------------
!                       k=0  k! x Gamma(nu+k+1)
!  (nu must be >= 0).
!---------------------------------------------------
      Parameter(MAXK=20,ZERO=0.d0)
      Complex z,z1
      Integer k
      Complex sum,tmp
      Real*8 Fact, Gamma
      sum = CMPLX(ZERO,ZERO)
      do k=0, MAXK
      !calculate (-z**2/4)**k
        tmp = (-z*z/4.d0)**k
      !divide by k!
        tmp = tmp / Fact(k)
      !divide by Gamma(nu+k+1)
        tmp = tmp / Gamma(dfloat(nu+k+1))
      !actualize sum
        sum = sum + tmp
      end do
      !calculate (z/2)**nu
      tmp = (z/2)**nu
      !multiply (z/2)**nu by sum
      z1 = tmp*sum
      return
      End

!end of file cbessj.f90      
      

      
      
      

