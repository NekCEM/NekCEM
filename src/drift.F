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

c...  inverse mass matrix including material constants 
        call invers2(qbm1,bmn,npts)
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
      include 'DRIFT'
      include 'POISSON'
      integer  i,j,k,l,ie,ieg

c...  FIXME: give correct numbers for parameters  
c...  3/13/14: all these to be defined in user file.

      mu_n =  1.0   
      mu_p =  1.0
      mu_e =  1.0
      K_beta= 1.0
      temperature= 1.0
      tau_cE= 1.0

      diff_n= K_beta*temperature*mu_n  ! charge to be divided later
      diff_p= K_beta*temperature*mu_p 
      diff_e= K_beta*temperature*mu_e 

      do i=1,npts
         d_permit(i)= 1.0
         d_permea(i)= 1.0
         charge  (i)= 1.0
      enddo
 
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
      include 'RK5'
      integer  ii
      integer i
      real xx,yy,tt,aa,bb,cc,glamax !fixme delete


      if (ifse) then
         do ii=1,5
            call rk_c (ii)
            call cem_drift_sem_op
            call rk_drift_ab(ii)
         enddo
      else
         do ii=1,5
            call rk_c (ii)
            call cem_drift_dg_op
            call rk_drift_ab(ii)
         enddo
      endif
      call drift_filter(0.05)

c-----compute the error at each time step
      tt=time+dt
      call usersol(tt, scn, scp, sce, scn, scp, sce)
      call copy(errn,cn,npts)
      call add2s2(errn,scn,-1.0,npts) !errn=cn-scn
      cc= glamax(Errn(1),Npts)
      if (nid.eq.0) write(6,*) 'ppp-rk45n',istep,cc
      call copy(errp,cp,npts)
      call add2s2(errp,scp,-1.0,npts) !errn=cn-scn
      cc= glamax(Errp(1),Npts)
      if (nid.eq.0) write(6,*) 'ppp-rk45p',istep,cc

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_drift_op_bdf
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      integer  ii
      integer i,ie,j
      real xx,yy,tt,aa,bb,cc,glamax !fixme delete


c.... 1st order semi-implicit scheme

c      call cem_drift_sem_semi
      call cem_drift_sem_semi_mix

c.... 2nd order semi-implicit scheme
c      call cem_drift_sem_bdf2
      !call cem_fsi_sem_semi !FIXME locally defined, delete later


c...  compute the error at each time step
c      tt=time+dt
c      call usersol(tt, scn, scp, sce, scn, scp, sce)
c      call copy(errn,cn,npts)
c      call add2s2(errn,scn,-1.0,npts) !errn=cn-scn
c      call copy(errn,scn,npts)
c      call add2s2(errn,potent,-1.0,npts) !errn=cn-scn
c      call copy(errn,cp,npts)
c      call add2s2(errn,cn,-2.0,npts) !errn=cn-scn
c      cc= glamax(Errn(1),Npts)
c      write(6,*) 'ppp-bdf1',istep,cc
c      write(6,*) 'ppp-bdf2',istep,cc

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
      call cem_drift_poisson ! return potent and electric field cEE
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

      call cem_drift_sem
      call cem_drift_sem_source   
      call cem_drift_invqmass_sem
c      call cem_drift_poisson ! return potent and electric feild cEE
      call cem_drift_potential_sem ! return potent and electric feild cEE

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_sem               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      real    wk0(lpts4)
      integer  i,isd,imsh
c... 2d case for now
      isd   = 1
      imsh  = 1
c... pre-computation for DtN stuffs
      call rone    (mask,npts )      ! for grids                         
      call setmaskp(mask,npts)
      if (nid.eq.0) write(6,*) ' done: setmaskp'

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) ' done: mult'
      call cem_drift_div(wk0,cEE(1,1),cEE(1,2),cEE(1,3))
      do i=1,npts
         h1n(i)=-diff_n/charge(i)
         h2n(i)=-mu_n*wk0(i)
         h3n(i)=-mu_n*cEE(i,1)
         h3n(i+npts)=-mu_n*cEE(i,2)
         h3n(i+2*npts)=0.0!-mu_n*cEE(i,3)
      enddo
      do i=1,npts
         h1p(i)=-diff_p/charge(i)
         h2p(i)=mu_p*wk0(i)
         h3p(i)=mu_p*cEE(i,1)
         h3p(i+npts)=mu_p*cEE(i,2)
         h3p(i+2*npts)=0.0!-mu_p*cEE(i,3)
      enddo
      call elliptic2d(resN,cN,h1n,h2n,h3n,imsh,isd)
      call dssum  (resN,nx1,ny1,nz1)
      call col2   (resN,mask,npts)
      call elliptic2d(resP,cP,h1p,h2p,h3p,imsh,isd)
      call dssum  (resP,nx1,ny1,nz1)
      call col2   (resP,mask,npts)
c      call elliptic2d(resE,cE,h1e,h2e,h3e,imsh,isd)
c      call dssum  (resE,nx1,ny1,nz1)
c      call col2   (resE,mask,npts)

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_sem_semi               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,isd,imsh,maxit
      real tol
c... 2d case for now
      isd   = 1
      imsh  = 1
c... pre-computation for DtN stuffs
      call rone    (mask,npts )      ! for grids                         
      call setmaskp(mask,npts)
      if (nid.eq.0) write(6,*) ' done: setmaskp'

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) ' done: mult'
 
c     h1=D_n, h3=-\grad D_n+\mu_n\grad\phi, h2=\mu_n\delta \phi+1.0/dt
c     if \phi=0, D_n constant, then h2=h3=0
c...  elliptic: 2nd + constant+ 1st
      do i=1,npts
         h1n(i)=diff_n/charge(i)*dt
         h2n(i)=1.0 ! mu_n*\delta\phi*dt+1.0
         h3n(i)=1.0*dt !-dx(D_n)+\mu_n*cEN(i,1)
         h3n(i+npts)=1.0*dt !-dy(D_n)+\mu_n*cEN(i,2)
         h3n(i+2*npts)=1.0*dt!-dz(D_n)+\mu_n*cEN(i,3)
      enddo
c      call rzero(h3n,3*npts) !turn off the first order term
c      do i=1,npts
c         h1p(i)=-diff_p/charge(i)
c         h2p(i)=0.0 !- mu_p*\delta\phi
c         h3p(i)=0.0 ! dx(D_p)-\mu_p*cEp(i,1)
c         h3p(i+npts)=0.0 ! dy(D_p)-\mu_p*cEp(i,2)
c         h3p(i+2*npts)=0.0 ! dz(D_p)-\mu_p*cEp(i,3)
c      enddo
c      do i=1,npts
c         h1e(i)=-diff_e/charge(i)
c         h2e(i)=0.0 !- mu_e*\delta\phi
c         h3e(i)=0.0 ! dx(D_e)-\mu_e*cEe(i,1)
c         h3e(i+npts)=0.0 ! dy(D_e)-\mu_e*cEe(i,2)
c         h3e(i+2*npts)=0.0 ! dz(D_e)-\mu_n*cEe(i,3)
c      enddo
c----- rhs=rhs*dt+ u^n
      call cem_drift_source_semi   
      call col2    (resN,bm1,npts)            ! B*rhs
      call dssum   (resN,nx1,ny1,nz1)
      call col2    (resN,mask,npts)  
c      call col2    (resP,bm1,npts)            ! B*rhs
c      call dssum   (resP,nx1,ny1,nz1)
c      call col2    (resP,mask,npts)  
c      call col2    (resE,bm1,npts)            ! B*rhs
c      call dssum   (resE,nx1,ny1,nz1)
c      call col2    (resE,mask,npts)  

      tol   = 1.e-16
      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES
      if (IFGMRES) then
c         call rzero     (cN,npts)
         call elliptic_gmres3(cN,resN,h1n,h2n,h3n
     $                   ,mult,mask,isd,imsh,npts,tol)
c         call rzero     (cP,npts)
c         call elliptic_gmres3(cP,resP,h1n,h2n,h3n
c     $                   ,mult,mask,isd,imsh,npts,tol)
c         call rzero     (cE,npts)
c         call elliptic_gmres3(cE,resE,h1n,h2n,h3n
c     $                   ,mult,mask,isd,imsh,npts,tol)
         if (nid.eq.0) write(6,*) 'done: elliptic_gmres3'
      else
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
         call exitt
      endif
c      if (IFCG) then
c          maxit = 2000
c          call cggo2 (cN,resN,h1n,h2n,mult,mask,imsh,isd,npts,maxit)
c      else
c         if (nid.eq.0) write(6,*) 'need to turn on the CG solver'
c         call exitt
c      endif
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
      call axhelm3(au,tmp1,a,b,imesh,isd)!compute AU=-a\delta u+ b u
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
c      call drift_wght_div  (resE,wEJ(1,1),wEJ(1,2),wEJ(1,3),
c     $                      cEJ(1,1),cEJ(1,2),cEJ(1,3))
      
      call drift_wght_grad (resNQ(1,1),resNQ(1,2),resNQ(1,3),cN)
      call drift_wght_grad (resPQ(1,1),resPQ(1,2),resPQ(1,3),cP)
c      call drift_wght_grad (resEQ(1,1),resEQ(1,2),resEQ(1,3),cE)
      call drift_wght_grad (resEN(1,1),resEN(1,2),resEN(1,3),potent)

      do i=1,npts
         resNQ(i,1)=diff_n/charge(i)*resNQ(i,1)-mu_n*cN(i)*resEN(i,1)
         resNQ(i,2)=diff_n/charge(i)*resNQ(i,2)-mu_n*cN(i)*resEN(i,2)
         resNQ(i,3)=diff_n/charge(i)*resNQ(i,3)-mu_n*cN(i)*resEN(i,3)
         resPQ(i,1)=diff_p/charge(i)*resPQ(i,1)+mu_p*cP(i)*resEN(i,1)
         resPQ(i,2)=diff_p/charge(i)*resPQ(i,2)+mu_p*cP(i)*resEN(i,2)
         resPQ(i,3)=diff_p/charge(i)*resPQ(i,3)+mu_p*cP(i)*resEN(i,3)
c         resEQ(i,1)=diff_e/charge(i)*resEQ(i,1)
c         resEQ(i,2)=diff_e/charge(i)*resEQ(i,2)
c         resEQ(i,3)=diff_e/charge(i)*resEQ(i,3)
      enddo


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
c           fcPJ(j,1)= cPJ(i,1)
c           fcPJ(j,2)= cPJ(i,2)
c           fcPJ(j,3)= cPJ(i,3)
c           fcEJ(j,1)= cEJ(i,1)
c           fcEJ(j,2)= cEJ(i,2)
c           fcEJ(j,3)= cEJ(i,3)
           fcN(j)   = cN(i)*diff_n/charge(i)
c           fcP(j)   = cP(i)*diff_p/charge(i)
c           fcE(j)   = cE(i)*diff_e/charge(i)
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
      integer i
      real tt,xx,yy,aa,bb,cc !delete later fixme
c      if (.not.ifsrc) return
     
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call col2    (wk1,bmn,npts) 
      call add2s2  (resN,wk1,1.,npts) 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_source  
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
      integer i
      real tt,xx,yy,aa,bb,cc !delete later fixme
c      if (.not.ifsrc) return

      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call col2    (wk1,bmn,npts) 
      call dssum   (wk1,nx1,ny1,nz1)
      call col2    (wk1,mask,npts)
      call add2s2  (resN,wk1,1.,npts) 
      call col2    (wk2,bmn,npts) 
      call dssum   (wk2,nx1,ny1,nz1)
      call col2    (wk2,mask,npts)
      call add2s2  (resP,wk2,1.,npts) 
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
      call col2(resN,qbm1,npts) !resN=resN*qbm1
      call col2(resP,qbm1,npts) !resP=resP*qbm1
      call col2(resE,qbm1,npts) !resE=resE*qbm1
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
           srflx(0*k+i)=srflx(0*k+i)+fcN(i)*(unx0+uny0) 
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
      call rk4_upd(cE,kcE,resE,cb,ca,dt,npts)

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
      subroutine cem_drift_force_dg  
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'EMWAVE'
      include 'RK5'
      integer  i,j,ie 
      real xx,yy,tt,tmp
c... FIXME  should be multiplied with bm1 for weighted results                   
      do i=1,npts
         rD(i) = 0.0 
         rR(i) = 0.0 !charge*(mu_n+mu_p)*cN(i)*cP(i)/d_permit(i)
         eG(i) = 0.0 
         eR(i) = 0.0 !cE(i)/tau_cE
      enddo
      do i=1,npts
         resN(i)=resN(i)*qbm1(i)
         resP(i)=resP(i)*qbm1(i)
         resE(i)=resE(i)*qbm1(i)
      enddo
      tt=RKtime
      do ie=1,nelt
         if      (if_in_region(ie).eq. 1) then
             do i=1,nxyz
             j=(ie-1)*nxyz+i
             xx=xm1(i,1,1,ie)
             yy=ym1(i,1,1,ie)

             !resN(j)=0
             tmp=sin(pi*xx)*sin(2*pi*yy)
             resN(j)=resN(j)+tmp*(cos(tt)+5.0*pi**2*sin(tt))
c             tmp=sin(pi*xx)*sin(2*pi*yy)
c             resP(j)=resP(j)+tmp*(cos(tt)+5.0*pi**2*sin(tt))
c             tmp=sin(pi*xx)*sin(pi*yy)
c             resE(j)=resE(j)+tmp*(cos(tt)+2.0*pi**2*sin(tt))
             
             resN(j)=resN(j)+rD(j)-rR(j)
c             resP(j)=resP(j)+rD(j)-rR(j)
c             resE(j)=resE(j)+eG(j)+1.0/4.0*rR(j)-eR(j)-rD(j)
             enddo
         elseif (if_in_region(ie).eq.-1) then
             do i=1,nxyz
             j=(ie-1)*nxyz+i
             xx=xm1(i,1,1,ie)
             yy=ym1(i,1,1,ie)

             tmp=sin(pi*xx)*sin(2*pi*yy)
             resN(j)=resN(j)+tmp*(cos(tt)+5.0*pi**2*sin(tt))
             !resP(j)=0.0
c             tmp=sin(pi*xx)*sin(2*pi*yy)
c             resP(j)=resP(j)+tmp*(cos(tt)+5.0*pi**2*sin(tt))
c             tmp=sin(pi*xx)*sin(pi*yy)
c             resE(j)=resE(j)+tmp*(cos(tt)+2.0*pi**2*sin(tt))
c             
c             resN(j)=resN(j)+rD(j)-rR(j)
c             resP(j)=resP(j)+rD(j)-rR(j)
c             resE(j)=resE(j)+eG(j)+1.0/4.0*rR(j)-eR(j)-rD(j)
             enddo
         endif
       enddo
      do i=1,npts
         resN(i)=resN(i)*bmn(i)
c         resP(i)=resP(i)*bmn(i)
c         resE(i)=resE(i)*bmn(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson
c----------------------------------------------------------------------
c...  implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT' 
      include 'POISSON'
      integer i
      real    tt !fixme delete finally
      if (nid.eq.0) write(6,*) '\n call cem_drift_poissonq'
      npts =  nx1*ny1*nz1*nelt
      call rone    (mask,npts)      ! mask is one everywhere                 
      call setmaskp(mask,npts)      ! mask is zero for PEC, PML 
      if (nid.eq.0) write(6,*) 'done: setmaskp'

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) 'done: dssum mult in poissonq'

      call cem_drift_poisson_rhs  
      call col2    (rhs,bm1,npts)            ! B*rhs
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,npts)
      if (nid.eq.0) write(6,*) 'done:: dssum rsh in poissonq'

      call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call rzero   (h2,npts)
      if (nid.eq.0) write(6,*) 'done:: set h1 and h2 in poissonq'


      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = 1.e-18 ! 1.e-10


      if (nid.eq.0) write(6,*) 'start IFCG/IFGMRES',IFCG,IFGMRES

c      if (IFCG) then
         call cggo2 (potent,rhs,h1,h2,mult,mask,imsh,isd,npts,maxit)
      bmax=glmax(potent,npts)
      bmin=glmin(potent,npts)
      if (nid.eq.0) write(6,*) 'phi_max/phi_min',bmax,bmin
c      endif

      if (nid.eq.0) write(6,*) 'done:: IFCG/IFGMRES',IFCG,IFGMRES
      if (nid.eq.0) write(6,*) 'done:: cem_drift_poissonq'
   
c      tt=time+dt
c      call usersol(tt, scn, scp, sce, scn, scp, sce)
c      call copy(potent,scn,npts)
      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent) !cEE=\grad potent 
      return
      end
c-----------------------------------------------------------------------
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
             w3(k)=w2(k)
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
      real wkn(lpts),wkp(lpts),wke(lpts)
      real wk1(lpts),wk2(lpts),wk3(lpts)
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      tt=time+dt
      if (istep .eq. 1) then
         call copy(rhs,cP,npts)
         call add2s2(rhs,cN,-1.0,npts)
         do i=1,npts
         rhs(i)=charge(i)*rhs(i)
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(2*pi*xx)*sin(2*pi*yy)
         aa=aa**2*sin(tt+1.0)**2   ! u
         cc=-8.0*pi**2*sin(tt+1.0)**2
     &      *(cos(4*pi*xx)*sin(2*pi*yy)**2+cos(4*pi*yy)*sin(2*pi*xx)**2)
         ee=cc-charge(i)*aa
         rhs(i)=rhs(i)+ee
         enddo
         return
      endif
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call cem_drift_div(wkN,cNJ(1,1),cNJ(1,2),cNJ(1,3))
      call cem_drift_div(wkP,cPJ(1,1),cPJ(1,2),cPJ(1,3))
      call cem_drift_div(wkE,cEJ(1,1),cEJ(1,2),cEJ(1,3))
      do i=1, npts
         rhs(i)=cP(i)-cN(i)
     $          +dt*(wkp(i)-wkn(i)+wk2(i)-wk1(i))
         rhs(i)=charge(i)*rhs(i)
c below define rhs for potential, which should be given in source FIXME
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(2*pi*xx)*sin(2*pi*yy)
         aa=aa**2*sin(tt+1.0)**2   ! u
         cc=-8.0*pi**2*sin(tt+1.0)**2
     &      *(cos(4*pi*xx)*sin(2*pi*yy)**2+cos(4*pi*yy)*sin(2*pi*xx)**2)
         ee=cc-charge(i)*aa
         rhs(i)=rhs(i)+ee
      enddo       
      return
      end
c-----------------------------------------------------------------------
      subroutine elliptic_gmres3(
     $           phi,res,h1,h2,h3,wt,mask,isd,imsh,n,tol)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,outer,isd,imsh
      real     phi(n),res(n),h1(n),h2(n),h3(3*n),wt(n),mask(n)
      real     tol,alpha,l,temp
      real*8   etime1,dnekclock

      if (nid.eq.0) write(6,*) 'start: elliptic_gmres3'

      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
c      call rzero(x,n)
      call copy(x,phi,n)
      iter  = 1
      call rzero(h,m*m)
      outer = 0
      do while (iconv.eq.0.and.iter.lt.1000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r,res,n)                  ! r = res
         else
            !update residual
            call copy   (r,res,n)                  ! r = res
            !call axhelm2(w,x,h1,h2,imsh,isd)       ! w = A x
            call elliptic2d(w,x,h1,h2,h3,imsh,isd)
            call dssum  (w,nx1,ny1,nz1)
            call col2   (w,mask,n)
            call add2s2 (r,w,-1.,n)                ! r = r - w
         endif
         
         gamma(1) = glsc3(r,r,wt,n)                ! gamma  = (r,r) 
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
c         write(6,*) 'initial residule',gamma(1)
                                                   !      1

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
c         write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            !call axhelm2(w,v(1,j),h1,h2,imsh,isd)  ! w = A v
            call elliptic2d(w,v(1,j),h1,h2,h3,imsh,isd)
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
        !write(6,*) 'end of forming m-th krylov subspace'
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

c--------------------------------------------------------------------- 
      subroutine cem_drift_sem_bdf2               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,isd,imsh,maxit
      real tol

      
      if(istep.eq.1) then
        call copy(cNO,cN,npts)
        call cem_drift_sem_semi   
        return 
      endif     
c... 2d case for now
      isd   = 1
      imsh  = 1
c... pre-computation for DtN stuffs
      call rone    (mask,npts )      ! for grids                         
      call setmaskp(mask,npts)
      if (nid.eq.0) write(6,*) ' done: setmaskp'

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) ' done: mult'
 
c     h1=D_n, h3=-\grad D_n+\mu_n\grad\phi, h2=\mu_n\delta \phi+1.0/dt
c     if \phi=0, D_n constant, then h2=h3=0
c...  elliptic: 2nd + constant+ 1st
      do i=1,npts
         h1n(i)=diff_n/charge(i)*dt
         h2n(i)=1.5 ! mu_n*\delta\phi*dt+1.0
         h3n(i)=1.0*dt !-dx(D_n)+\mu_n*cEN(i,1)
         h3n(i+npts)=1.0*dt !-dy(D_n)+\mu_n*cEN(i,2)
         h3n(i+2*npts)=1.0*dt!-dz(D_n)+\mu_n*cEN(i,3)
      enddo
c      call rzero(h3n,3*npts) !turn off the first order term
c      do i=1,npts
c         h1p(i)=-diff_p/charge(i)
c         h2p(i)=0.0 !- mu_p*\delta\phi
c         h3p(i)=0.0 ! dx(D_p)-\mu_p*cEp(i,1)
c         h3p(i+npts)=0.0 ! dy(D_p)-\mu_p*cEp(i,2)
c         h3p(i+2*npts)=0.0 ! dz(D_p)-\mu_p*cEp(i,3)
c      enddo
c      do i=1,npts
c         h1e(i)=-diff_e/charge(i)
c         h2e(i)=0.0 !- mu_e*\delta\phi
c         h3e(i)=0.0 ! dx(D_e)-\mu_e*cEe(i,1)
c         h3e(i+npts)=0.0 ! dy(D_e)-\mu_e*cEe(i,2)
c         h3e(i+2*npts)=0.0 ! dz(D_e)-\mu_n*cEe(i,3)
c      enddo
c----- rhs=rhs*dt+ u^n
      call cem_drift_source_bdf2   
      call add2s2  (resN,cNO,-0.5,npts) 
      call col2    (resN,bm1,npts)            ! B*rhs
      call dssum   (resN,nx1,ny1,nz1)
      call col2    (resN,mask,npts)  
c      call col2    (resP,bm1,npts)            ! B*rhs
c      call dssum   (resP,nx1,ny1,nz1)
c      call col2    (resP,mask,npts)  
c      call col2    (resE,bm1,npts)            ! B*rhs
c      call dssum   (resE,nx1,ny1,nz1)
c      call col2    (resE,mask,npts)  

      call copy(cNO,cN,npts)
      tol   = 1.e-16
      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES
      if (IFGMRES) then
c         call rzero     (cN,npts)
         call elliptic_gmres3(cN,resN,h1n,h2n,h3n
     $                   ,mult,mask,isd,imsh,npts,tol)
c         call rzero     (cP,npts)
c         call elliptic_gmres3(cP,resP,h1n,h2n,h3n
c     $                   ,mult,mask,isd,imsh,npts,tol)
c         call rzero     (cE,npts)
c         call elliptic_gmres3(cE,resE,h1n,h2n,h3n
c     $                   ,mult,mask,isd,imsh,npts,tol)
         if (nid.eq.0) write(6,*) 'done: elliptic_gmres3'
      else
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
         call exitt
      endif
c      if (IFCG) then
c          maxit = 2000
c          call cggo2 (cN,resN,h1n,h2n,mult,mask,imsh,isd,npts,maxit)
c      else
c         if (nid.eq.0) write(6,*) 'need to turn on the CG solver'
c         call exitt
c      endif
      return
      end
c=======================================================================
      subroutine cem_drift_source_bdf2    
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
      call rzero   (resN,npts)
      call add2s2  (resN,wk1,1.*dt,npts) 
      call add2s2  (resN,cN,2.0,npts) 
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_drift_sem_semi_mix               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,isd,imsh,maxit
      real wk1(lpts),wk2(lpts),wk3(lpts),wk0(lpts)
      real tol
      real tt,glamax !fixme delete finally


c... 2d case for now
      isd   = 1
      imsh  = 1
c... pre-computation for DtN stuffs
      call rone    (mask,npts )      ! for grids                         
      call setmaskp(mask,npts)
      if (nid.eq.0) write(6,*) ' done: setmaskp'

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) ' done: mult'
c...  elliptic: 2nd + constant+ 1st
      call cem_drift_div(wk0,cEE(1,1),cEE(1,2),cEE(1,3))
      do i=1,npts
         h1n(i)=diff_n/charge(i)*dt
         h2n(i)=1.0+mu_n*wk0(i)*dt
         h3n(i)=mu_n*cEE(i,1)*dt
         h3n(i+npts)=mu_n*cEE(i,2)*dt
         h3n(i+2*npts)=0.0!mu_n*cEE(i,3)*dt
      enddo
      call rzero(h3n,3*npts) !turn off the first order term
      do i=1,npts
         h1p(i)=diff_p/charge(i)*dt
         h2p(i)=1.0-mu_p*wk0(i)*dt
         h3p(i)=-mu_p*cEE(i,1)*dt
         h3p(i+npts)=-mu_p*cEE(i,2)*dt
         h3p(i+2*npts)=0.0!-mu_p*cEE(i,3)*dt
      enddo
      call rzero(h3p,3*npts) !turn off the first order term
c      do i=1,npts
c         h1e(i)=-diff_e/charge(i)
c         h2e(i)=0.0 !- mu_e*\delta\phi
c         h3e(i)=0.0 ! dx(D_e)-\mu_e*cEe(i,1)
c         h3e(i+npts)=0.0 ! dy(D_e)-\mu_e*cEe(i,2)
c         h3e(i+2*npts)=0.0 ! dz(D_e)-\mu_n*cEe(i,3)
c      enddo
c----- rhs=rhs*dt+ u^n
      call cem_drift_source_semi_mix   
      call col2    (resN,bm1,npts)            ! B*rhs
      call dssum   (resN,nx1,ny1,nz1)
      call col2    (resN,mask,npts)  
      call col2    (resP,bm1,npts)            ! B*rhs
      call dssum   (resP,nx1,ny1,nz1)
      call col2    (resP,mask,npts)  
c      call col2    (resE,bm1,npts)            ! B*rhs
c      call dssum   (resE,nx1,ny1,nz1)
c      call col2    (resE,mask,npts)  

      call copy(cNO,cN,npts)
      tol   = 1.e-16
      if (nid.eq.0) write(6,*) 'IFCG/IFGMRES',IFCG,IFGMRES
      if (IFGMRES) then
         call elliptic_gmres3(cN,resN,h1n,h2n,h3n
     $                   ,mult,mask,isd,imsh,npts,tol)
         call elliptic_gmres3(cP,resP,h1p,h2p,h3p
     $                   ,mult,mask,isd,imsh,npts,tol)
c         call elliptic_gmres3(cE,resE,h1n,h2n,h3n
c     $                   ,mult,mask,isd,imsh,npts,tol)
         if (nid.eq.0) write(6,*) 'done: elliptic_gmres3'
      else
         if (nid.eq.0) write(6,*) 'need to turn on the gmres solver'
         call exitt
      endif
      call cem_drift_potential_sem ! return potent and electric feild cEE
      call drift_grad(cNJ(1,1),cNJ(1,2),cNJ(1,3),cN)
      call drift_grad(cPJ(1,1),cPJ(1,2),cPJ(1,3),cP)
c      call drift_grad(cEJ(1,1),cEJ(1,2),cEJ(1,3),cE)

      do i=1,npts
      cNJ(i,1)=diff_n/charge(i)*cNJ(i,1)-mu_n*cN(i)*cEE(i,1)
      cNJ(i,2)=diff_n/charge(i)*cNJ(i,2)-mu_n*cN(i)*cEE(i,2)
      cNJ(i,3)=diff_n/charge(i)*cNJ(i,3)-mu_n*cN(i)*cEE(i,3)
      cPJ(i,1)=diff_p/charge(i)*cPJ(i,1)+mu_p*cP(i)*cEE(i,1)
      cPJ(i,2)=diff_p/charge(i)*cPJ(i,2)+mu_p*cP(i)*cEE(i,2)
      cPJ(i,3)=diff_p/charge(i)*cPJ(i,3)+mu_p*cP(i)*cEE(i,3)
      enddo 
      tt=time+dt
      call usersol(tt, scn, scp, sce, scn, scp, sce)
      call copy(errn,cn,npts)
      call add2s2(errn,scn,-1.0,npts)
      write(6,*) 'xxxx',istep,glamax(errn,npts)
      call copy(errp,cp,npts)
      call add2s2(errp,scp,-1.0,npts)
      write(6,*) 'xxxx',istep,glamax(errp,npts)
c      call exitt
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_source_semi_mix    
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
c      if (.not.ifsrc) return
      npts3=npts*3
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)

      call drift_grad(wk0(1),wk0(npts+1),wk0(2*npts+1),cN)
      call col2    (wk0,cEE,npts3)            ! B*rhs
      call copy    (resN,cN,npts)
      call add2s2  (resN,wk0,-1.*dt,npts) 
      call add2s2  (resN,wk1,1.*dt,npts) 
      
      call drift_grad(wk0(1),wk0(npts+1),wk0(2*npts+1),cP)
      call col2    (wk0,cEE,npts3)            ! B*rhs
      call copy    (resP,cP,npts)
      call add2s2  (resP,wk0,1.*dt,npts) 
      call add2s2  (resP,wk2,1.*dt,npts) 

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
      subroutine cem_drift_potential_sem               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,isd,imsh,maxit
      real wk1(lpts),wk2(lpts),wk3(lpts),wk0(lpts)
      real wkx(lpts),wky(lpts),wkz(lpts)
      real tol
      real tt,glamax !fixme delete finally
    
c... 2d case for now
      isd   = 1
      imsh  = 1
c... pre-computation for DtN stuffs
      call rone    (mask,npts )      ! for grids                         
      call setmaskp(mask,npts)
      if (nid.eq.0) write(6,*) ' done: setmaskp'

      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)
      if (nid.eq.0) write(6,*) ' done: mult'
 
c...  elliptic: 2nd + 1st
      if (istep .eq. 1) then
      call rone(h1,npts)
      else
      do i=1,npts
         h1(i)=1.0+dt*(mu_n*cN(i)+mu_p*cP(i))*charge(i)
      enddo
      endif
      call rzero(h2,npts)
      call cem_drift_poisson_rhs
      call col2    (rhs,bm1,npts)            ! B*rhs
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,npts)  
      
      tol   = 1.e-16
      maxit = 2000
      call cggo2 (potent,rhs,h1,h2,mult,mask,imsh,isd,npts,maxit)
      call drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent)
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
      implicit none
      include 'SIZE' 
      include 'TOTAL' 
c     include 'EMWAVE' 
      include 'DEALIAS' 
c     ddddd
c     Combined weighted divl:   w0 = B*div(u1,u2,u3) on a finer mesh
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

            wd1(k) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k)
     $               + u1tw*txmd(k))

            wd2(k) = (  u2rw*rymd(k)
     $               + u2sw*symd(k)
     $               + u2tw*tymd(k))

            wd3(k) = (  u3rw*rzmd(k)
     $               + u3sw*szmd(k)
     $               + u3tw*tzmd(k))

            wd0(k) = wd1(k)+wd2(k)+wd3(k)

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

            wd1(k) = (  u1rw*rxmd(k)
     $               + u1sw*sxmd(k))

            wd2(k) = (  u2rw*rymd(k)
     $               + u2sw*symd(k))

            wd0(k) = wd1(k)+wd2(k)      

           enddo

         endif

         call intp_rstd (w0(j1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
