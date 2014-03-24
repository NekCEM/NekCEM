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
      integer i,ie,j
      real xx,yy,tt,aa,bb,cc,glamax !fixme delete
      real myshx(lx1,ly1,lz1,lelt)
      real myshy(lx1,ly1,lz1,lelt)
      real myshz(lx1,ly1,lz1,lelt)
      real mysex(lx1,ly1,lz1,lelt)
      real mysey(lx1,ly1,lz1,lelt)
      real mysez(lx1,ly1,lz1,lelt)

      do ii=1,5
         call rk_c (ii)
c         call cem_drift_op
c         write(6,*) 'xxx--using DG'
         call cem_drift_sem_op
         write(6,*) 'xxx--using SEM'
         call rk_drift_ab(ii)
      enddo
         tt=rktime
         call usersol(tt, scn, scp, sce, scn, scp, sce)
         do i=1,npts
         errn(i)=cn(i)-scn(i)
         enddo
         cc= glamax(Errn(1),Npts)
         write(6,*) 'ppp',istep,cc
c      call exitt
      
      return
      end

c-----------------------------------------------------------------------
      subroutine cem_drift_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      common /scrns/ srflxs(12*2*ldim*lx1*lz1*lelt)        
      real           srflxs
      real xx,yy,tt,cc,glamax !fixme delet after testing
      integer i  !fixme delet after testing
c      call cem_drift_poisson ! return potent and electric field cEE
      !FIXME DONE ON TESTING POISSON SOLVER, start testing on DD solver
      call cem_drift
      call cem_drift_restrict_to_face
      call cem_drift_flux (srflxs)
      call cem_drift_add_flux_to_res (srflxs)
      call cem_drift_source   
      call cem_drift_invqmass
      
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

c      call cem_drift_poisson ! return potent and electric feild cEE
      call cem_drift_sem
      call cem_drift_source   
      call cem_drift_invqmass_sem

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
      integer  i,isd,imsh
      real h1n(lpts4),h2n(lpts4),h3n(lpts4*3)
      real h1p(lpts4),h2p(lpts4),h3p(lpts4*3)
      real h1e(lpts4),h2e(lpts4),h3e(lpts4*3)
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
 
      call rone    (h1n,npts)        ! for  [ h1*u'' + h3*u'+h2*u ] 
      call chsign  (h1n,npts)        ! negative for helmholtz, positive fo     r poisson
      call rone    (h2n,npts)       
      call rone    (h3n,npts*3)      
c     h1=-D_n, h3=-\grad D_n+\mu_n\grad\phi, h2=\mu_n\delta \phi
c     if \phi=0, D_n constant, then h2=h3=0
c...  elliptic: 2nd + constant+ 1st
      call rzero   (h2n,npts)
      call rzero   (h3n,npts*3)
      do i=1,npts
         h1n(i)=-diff_n/charge(i)
         h2n(i)=0.0 ! mu_n*\delta\phi
         h3n(i)=0.0 !-dx(D_n)+\mu_n*cEN(i,1)
         h3n(i+npts)=0.0 !-dy(D_n)+\mu_n*cEN(i,2)
         h3n(i+2*npts)=0.0 !-dz(D_n)+\mu_n*cEN(i,3)
      enddo
      call elliptic2d(resN,cN,h1n,h2n,h3n,imsh,isd)
c      call elliptic2d(resP,cP,h1p,h2p,h3p,imsh,isd)
c      call elliptic2d(resE,cE,h1e,h2e,h3e,imsh,isd)

      return
      end
c--------------------------------------------------------------------- 
      subroutine elliptic2d(au,u,a,b,c,imsh,isd)              
c--------------------------------------------------------------------- 
C------------------------------------------------------------------
C
C     Compute the 2nd elliptic  matrix-vector product,
C     AU = a \delta u  + c \cdot \grad u + b u , for NEL elements.
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
      call axhelm2(au,u,a,b,imesh,isd)
      call drift_grad(wrkdx,wrkdy,wrkdz,u)  
      k=npts
      if (if3d) then
      do i=1,npts
         tmp1(i)=wrkdx(i)*c(i)+wrkdy(i)*c(i+k)+wrkdz(i)*c(i+2*k)
         tmp1(i)=tmp1(i)*bmn(i)
      enddo
      else
      do i=1,npts
         tmp1(i)=wrkdx(i)*c(i)+wrkdy(i)*c(i+k)
         tmp1(i)=tmp1(i)*bmn(i)
      enddo
      endif
      call add2s2 (au,tmp1,1.0,npts) 
      call dssum  (au,nx1,ny1,nz1)
      call col2   (au,mask,npts)
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift                
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
c      call drift_wght_div  (resP,wPJ(1,1),wPJ(1,2),wPJ(1,3),
c     $                      cPJ(1,1),cPJ(1,2),cPJ(1,3))
c      call drift_wght_div  (resE,wEJ(1,1),wEJ(1,2),wEJ(1,3),
c     $                      cEJ(1,1),cEJ(1,2),cEJ(1,3))
      
      call drift_wght_grad (resNQ(1,1),resNQ(1,2),resNQ(1,3),cN)
c      call drift_wght_grad (resPQ(1,1),resPQ(1,2),resPQ(1,3),cP)
c      call drift_wght_grad (resEQ(1,1),resEQ(1,2),resEQ(1,3),cE)
c      call drift_wght_grad (resEN(1,1),resEN(1,2),resEN(1,3),potent)

      do i=1,npts
         resNQ(i,1)=diff_n/charge(i)*resNQ(i,1)!-mu_n*cN(i)*resEN(i,1)
         resNQ(i,2)=diff_n/charge(i)*resNQ(i,2)!-mu_n*cN(i)*resEN(i,2)
         resNQ(i,3)=diff_n/charge(i)*resNQ(i,3)!-mu_n*cN(i)*resEN(i,3)
c         resPQ(i,1)=diff_p/charge(i)*resPQ(i,1)!+mu_p*cP(i)*resEN(i,1)
c         resPQ(i,2)=diff_p/charge(i)*resPQ(i,2)!+mu_p*cP(i)*resEN(i,2)
c         resPQ(i,3)=diff_p/charge(i)*resPQ(i,3)!+mu_p*cP(i)*resEN(i,3)
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
c         resP (i)   = resP (i)   + a*srflx( 1*k+j)
c         resE (i)   = resE (i)   + a*srflx( 2*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx( 3*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx( 4*k+j)
         resNQ(i,3) = resNQ(i,3) + a*srflx( 5*k+j)
c         resPQ(i,1) = resPQ(i,1) + a*srflx( 6*k+j)
c         resPQ(i,2) = resPQ(i,2) + a*srflx( 7*k+j)
c         resPQ(i,3) = resPQ(i,3) + a*srflx( 8*k+j)
c         resEQ(i,1) = resEQ(i,1) + a*srflx( 9*k+j)
c         resEQ(i,2) = resEQ(i,2) + a*srflx(10*k+j)
c         resEQ(i,3) = resEQ(i,3) + a*srflx(11*k+j)
         enddo
      else
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resN (i)   = resN (i)   + a*srflx( 0*k+j)
c         resP (i)   = resP (i)   + a*srflx( 1*k+j)
c         resE (i)   = resE (i)   + a*srflx( 2*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx( 3*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx( 4*k+j)
c         resPQ(i,1) = resPQ(i,1) + a*srflx( 5*k+j)
c         resPQ(i,2) = resPQ(i,2) + a*srflx( 6*k+j)
c         resEQ(i,1) = resEQ(i,1) + a*srflx( 7*k+j)
c         resEQ(i,2) = resEQ(i,2) + a*srflx( 8*k+j)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source    
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
      integer i
      real tt,xx,yy,aa,bb,cc !delete later fixme
      if (.not.ifsrc) return

c...  FIXME ying/misun 10/09/2012
      tt=RKtime+1.0
      do i=1,npts
         wnx(i)=bmn(i) !cN(i)
         wny(i)=bmn(i) !cP(i)
         wnz(i)=bmn(i)       
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(2*pi*xx)*sin(2*pi*yy)
         bb=cos(tt)*aa
         cc=8.0*pi**2*aa*sin(tt) 
         resn(i)=resn(i)+(bb+cc)*bmn(i)
c         rese(i)=rese(i)+(bb+cc)*bmn(i)
         aa=2.0*sin(2*pi*xx)*sin(2*pi*yy)
         bb=cos(tt)*aa
         cc=8.0*pi**2*aa*sin(tt) 
         resp(i)=resp(i)+(bb+cc)*bmn(i)
      enddo
c      call usersrc (57,resN,resP,wnz,wnx,wny,wnz)
c      call exitt                                                           
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
      integer i
      do i=1,npts
         resN (i)  =resN(i)   *qbm1(i)
c         resP (i)  =resP(i)   *qbm1(i)
c         resE (i)  =resE(i)   *qbm1(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_invqmass
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      integer i

      do i=1,npts
         resN (i)  =resN(i)   *qbm1(i)
c         resP (i)  =resP(i)   *qbm1(i)
c         resE (i)  =resE(i)   *qbm1(i)
         resNQ(i,1)=resNQ(i,1)*qbm1(i) ! arrays for copy to currents
         resNQ(i,2)=resNQ(i,2)*qbm1(i) ! not used for RK45      
         resNQ(i,3)=resNQ(i,3)*qbm1(i)
c         resPQ(i,1)=resPQ(i,1)*qbm1(i)
c         resPQ(i,2)=resPQ(i,2)*qbm1(i)
c         resPQ(i,3)=resPQ(i,3)*qbm1(i)
c         resEQ(i,1)=resEQ(i,1)*qbm1(i)
c         resEQ(i,2)=resEQ(i,2)*qbm1(i)
c         resEQ(i,3)=resEQ(i,3)*qbm1(i)
         cNJ(i,1)  =resNQ(i,1)           
         cNJ(i,2)  =resNQ(i,2)
         cNJ(i,3)  =resNQ(i,3)
c         cPJ(i,1)  =resPQ(i,1)
c         cPJ(i,2)  =resPQ(i,2)
c         cPJ(i,3)  =resPQ(i,3)            
c         cEJ(i,1)  =resEQ(i,1)
c         cEJ(i,2)  =resEQ(i,2)
c         cEJ(i,3)  =resEQ(i,3)            
      enddo
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
c           srflx( 1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)+unz0*fcPJ(i,3))
c           srflx( 2*k+i)=-(unx0*fcEJ(i,1)+uny0*fcEJ(i,2)+unz0*fcEJ(i,3))
           srflx( 3*k+i)=-unx0*fcN(i) 
           srflx( 4*k+i)=-uny0*fcN(i) 
           srflx( 5*k+i)=-unz0*fcN(i) 
c           srflx( 6*k+i)=-unx0*fcP(i) 
c           srflx( 7*k+i)=-uny0*fcP(i) 
c           srflx( 8*k+i)=-unz0*fcP(i) 
c           srflx( 9*k+i)=-unx0*fcE(i) 
c           srflx(10*k+i)=-uny0*fcE(i) 
c           srflx(11*k+i)=-unz0*fcE(i) 

         enddo

        else
         do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
C.........FIXME figure out the right sign +/-
           srflx(0*k+i)=-(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)) 
c           srflx(1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)) 
c           srflx(2*k+i)=-(unx0*fcEJ(i,1)+uny0*fcEJ(i,2)) 
           srflx(3*k+i)=-unx0*fcN(i)  ! NxfHx
           srflx(4*k+i)=-uny0*fcN(i)  ! NxfHx
c           srflx(5*k+i)=-unx0*fcP(i)  ! NxfHx
c           srflx(6*k+i)=-uny0*fcP(i)  ! NxfHx
c           srflx(7*k+i)=-unx0*fcE(i)  ! NxfHx
c           srflx(8*k+i)=-uny0*fcE(i)  ! NxfHx

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
             srflx( 3*k+i)= 2.d0*(srflx( 3*k+i)+unx0*fcN(i0))
             srflx( 4*k+i)= 2.d0*(srflx( 4*k+i)+uny0*fcN(i0))
             srflx( 5*k+i)= 2.d0*(srflx( 5*k+i)+unz0*fcN(i0))
c             srflx( 6*k+i)= 2.d0*(srflx( 6*k+i)+unx0*fcP(i0))
c             srflx( 7*k+i)= 2.d0*(srflx( 7*k+i)+uny0*fcP(i0))
c             srflx( 8*k+i)= 2.d0*(srflx( 8*k+i)+unz0*fcP(i0))
c             srflx( 9*k+i)= 2.d0*(srflx( 9*k+i)+unx0*fcE(i0))
c             srflx(10*k+i)= 2.d0*(srflx(10*k+i)+uny0*fcE(i0))
c             srflx(11*k+i)= 2.d0*(srflx(11*k+i)+unz0*fcE(i0))

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
             srflx(1*k+i)= 0.d0
             srflx(2*k+i)= 0.d0
             srflx(3*k+i)= 2.d0*(srflx(3*k+i)+unx0*fcN(i0))
             srflx(4*k+i)= 2.d0*(srflx(4*k+i)+uny0*fcN(i0))
c             srflx(5*k+i)= 2.d0*(srflx(5*k+i)+unx0*fcP(i0))
c             srflx(6*k+i)= 2.d0*(srflx(6*k+i)+uny0*fcP(i0))
c             srflx(7*k+i)= 2.d0*(srflx(7*k+i)+unx0*fcE(i0))
c             srflx(8*k+i)= 2.d0*(srflx(8*k+i)+uny0*fcE(i0))

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
c      call rk4_upd(cP,kcP,resP,cb,ca,dt,npts)
c      call rk4_upd(cE,kcE,resE,cb,ca,dt,npts)

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
      real xx,yy 

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
      call invcol2 (rhs,d_permit,npts)   ! B*rhs
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

      if (IFCG) then
         call cggo2 (potent,rhs,h1,h2,mult,mask,imsh,isd,npts,maxit)
      bmax=glmax(potent,npts)
      bmin=glmin(potent,npts)
      if (nid.eq.0) write(6,*) 'phi_max/phi_min',bmax,bmin
      endif
      if (IFGMRES) then
         call  hmh_gmres2(potent,rhs,h1,h2,mult,mask,isd,imsh,npts,tol)
      endif

      if (nid.eq.0) write(6,*) 'done:: IFCG/IFGMRES',IFCG,IFGMRES
      if (nid.eq.0) write(6,*) 'done:: cem_drift_poissonq'

      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent) !cEE=\grad potent 
      !FIXME TEST CASE DELTE FINALLY
c      do i=1,npts
c         xx=xm1(i,1,1,1)
c         yy=ym1(i,1,1,1)
c         spotent(i)=sin(pi*xx)*sin(pi*yy)
c         write(6,*) i,spotent(i)-potent(i)
c         spotent(i)=pi*cos(pi*xx)*sin(pi*yy)
c         write(6,*) i,spotent(i)-cEE(i,1)
c         spotent(i)=pi*sin(pi*xx)*cos(pi*yy)
c         write(6,*) i,spotent(i)-cEE(i,2)
c      enddo
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
      real xx,yy
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      do i=1, npts
         rhs(i)=charge(i)*(cP(i)-cN(i))
c.... FIXME for testing, delete finally
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         charge(i)=2.0*pi**2
         rhs(i)=charge(i)*(sin(pi*xx)*sin(pi*yy))
      enddo       
      return
      end
c-----------------------------------------------------------------------
