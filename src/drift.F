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

      if (nid.eq.0) write(6,*) 'set temporary variables for geom'
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
          call rzero (scN  ,npts)
          call rzero (scP  ,npts)
          call rzero (cNQ  ,npts3)
          call rzero (cPQ  ,npts3)
          call rzero (scNQ ,npts3)
          call rzero (scPQ ,npts3)
          call rzero (cNJ  ,npts3)
          call rzero (cPJ  ,npts3)
          call rzero (DcN  ,npts)
          call rzero (DcP  ,npts)
          call rzero (RG   ,npts)
          call rzero (RU   ,npts)
          call rzero (resN ,npts)
          call rzero (resP ,npts)
          call rzero (resNQ,npts3)
          call rzero (resPQ,npts3)
          call rzero (resNJ,npts3)
          call rzero (resPJ,npts3)
          call rzero (errN ,npts)
          call rzero (errP ,npts)
          call rzero (kcN  ,npts)
          call rzero (kcP  ,npts)
          call rzero (qbm1N,npts)
          call rzero (qbm1P,npts)

          call rzero (fcN  ,nxzfl)
          call rzero (fcP  ,nxzfl)
          call rzero (fcNJ ,nxzfl3)
          call rzero (fcPJ ,nxzfl3)

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
        call cem_set_fc_ptr       !global numbering index on face
        call cem_drift_uvp       !returns planck(i)             

c...  inverse mass matrix including material constants 
        call invers2(qbm1N,bmn,npts)
        call invers2(qbm1P,bmn,npts)
        if (nid.eq.0) write(6,*) 'done: inverse mass matrix'

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_uvp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      integer  i,ie
      real     e_mobility_max,e_mobility_min
      real     h_mobility_max,h_mobility_min
      real     glmax,glmin

c...  nekuvp returns planck(i), using planck_c defined in uservp

      do ie=1,nelt
          call nekuvp(ie)  
      enddo 

      e_mobility_max = glmax(e_mobility,npts) 
      e_mobility_min = glmin(e_mobility,npts) 
      h_mobility_max = glmax(h_mobility,npts) 
      h_mobility_min = glmin(h_mobility,npts) 

      if (nid.eq.0) write(6,*) 'max e_mobility=',e_mobility_max
      if (nid.eq.0) write(6,*) 'min e_mobility=',e_mobility_min
      if (nid.eq.0) write(6,*) 'max h_mobility=',h_mobility_max
      if (nid.eq.0) write(6,*) 'min h_mobility=',h_mobility_min

      if ((e_mobility_max.le.0).or.(e_mobility_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define mu_n in uservp(.usr)'
          call exitt
      endif
      if ((h_mobility_max.le.0).or.(h_mobility_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define mu_p in uservp(.usr)'
          call exitt
      endif
                
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
      integer  ii

      do ii=1,5
         call rk_c (ii)
         call cem_drift_op
         call rk_drift_ab(ii)
      enddo
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
      common /scrns/ srflxs(8*2*ldim*lx1*lz1*lelt)        
      real           srflxs

      call cem_drift
      call cem_drift_restrict_to_face
      call cem_drift_flux (srflxs)
      call cem_drift_add_flux_to_res (srflxs)
      call cem_drift_source   
      call cem_drift_invqmass

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'

      call drift_wght_grad (resNQ(1,1),resNQ(1,2),resNQ(1,3),cN)
      call drift_wght_grad (resPQ(1,1),resPQ(1,2),resPQ(1,3),cP)
      !call current_density :cNJ=-charge*e_mobility*cN*grad(potent) +charge*DcN*resNQ
      !call current_density :cPJ=-charge*h_mobility*cJ*grad(potent) +charge*DcP*resPQ
      call drift_wght_div  (resN,wNJ(1,1),wNJ(1,2),wNJ(1,3),
     $                           cNJ(1,1),cNJ(1,2),cNJ(1,3))
      call drift_wght_div  (resP,wPJ(1,1),wPJ(1,2),wPJ(1,3),
     $                           cPJ(1,1),cPJ(1,2),cPJ(1,3))

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
      integer  i,j

      do j=1,ncemface
          i= cemface(j)
           fcNJ(j,1)= cNJ(i,1)
           fcNJ(j,2)= cNJ(i,2)
           fcNJ(j,3)= cNJ(i,3)
           fcPJ(j,1)= cPJ(i,1)
           fcPJ(j,2)= cPJ(i,2)
           fcPJ(j,3)= cPJ(i,3)
           fcN(j)   = cN(i)
           fcP(j)   = cP(i)
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
         resN (i)   = resN (i)   + a*srflx(0*k+j)
         resP (i)   = resP (i)   + a*srflx(1*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx(2*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx(3*k+j)
         resNQ(i,3) = resNQ(i,3) + a*srflx(4*k+j)
         resPQ(i,1) = resPQ(i,1) + a*srflx(5*k+j)
         resPQ(i,2) = resPQ(i,2) + a*srflx(6*k+j)
         resPQ(i,3) = resPQ(i,3) + a*srflx(7*k+j)
         enddo
      else
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resN (i)   = resN (i)   + a*srflx(0*k+j)
         resP (i)   = resP (i)   + a*srflx(1*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx(2*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx(3*k+j)
         resPQ(i,1) = resPQ(i,1) + a*srflx(4*k+j)
         resPQ(i,2) = resPQ(i,2) + a*srflx(5*k+j)
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

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      integer i

      if (.not.ifsrc) return

c...  FIXME ying/misun 10/09/2012
      do i=1,npts
         wnx(i)=bmn(i) !cN(i)
         wny(i)=bmn(i) !cP(i)
         wnz(i)=bmn(i)       
      enddo
      call usersrc (57,resN,resP,wnz,wnx,wny,wnz)
                                                                 
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
         resN (i)  =resN(i)   *qbm1N(i)
         resP (i)  =resP(i)   *qbm1P(i)
         resNJ(i,1)=resNJ(i,1)*qbm1N(i) ! arrays for copy to currents
         resNJ(i,2)=resNJ(i,2)*qbm1N(i) ! not used for RK45      
         resNJ(i,3)=resNJ(i,3)*qbm1N(i)
         resPJ(i,1)=resPJ(i,1)*qbm1P(i)
         resPJ(i,2)=resPJ(i,2)*qbm1P(i)
         resPJ(i,3)=resPJ(i,3)*qbm1P(i)
         cNJ(i,1)   =resNJ(i,1)           
         cNJ(i,2)   =resNJ(i,2)
         cNJ(i,3)   =resNJ(i,3)
         cPJ(i,1)   =resPJ(i,1)
         cPJ(i,2)   =resPJ(i,2)
         cPJ(i,3)   =resPJ(i,3)            
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
     
      real*8  dclock,t0
      real    srflx(1)
      integer i,j,k,e,f,i0
      real    unx0,uny0,unz0
      real    tmp,tmpi1,tmpi2,tmpi3,tmpr1,tmpr2,tmpr3
      real    tmpur,tmpui,grad_r
      real    fjump_condi ,fjump_condr 
      real    fjump_condi1,fjump_condr1 
      real    fjump_condi2,fjump_condr2 
      real    fjump_condi3,fjump_condr3 

      integer icalld
      save    icalld
      data    icalld/0/
      character CB*3

      !call usersol(RKtime,sQr,sQi,sQi,sUr,sUi,sUi) !FIXME

      k = nxzfl

      if (if3d) then
        do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
           unz0=unzm(i)
           srflx(0*k+i)=-(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)+unz0*fcNJ(i,3))
           srflx(1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)+unz0*fcPJ(i,3))
           srflx(2*k+i)= unx0*fcN(i) 
           srflx(3*k+i)= uny0*fcN(i) 
           srflx(4*k+i)= unz0*fcN(i) 
           srflx(5*k+i)= unx0*fcP(i) 
           srflx(6*k+i)= uny0*fcP(i) 
           srflx(7*k+i)= unz0*fcP(i) 

         enddo

        else
         do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)

           srflx(0*k+i)= -(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)) 
           srflx(1*k+i)= +(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)) 
           srflx(2*k+i)= - unx0*fcN(i)  ! NxfHx
           srflx(3*k+i)= - uny0*fcN(i)  ! NxfHx
           srflx(4*k+i)= - unx0*fcP(i)  ! NxfHx
           srflx(5*k+i)= - uny0*fcP(i)  ! NxfHx

         enddo
        endif

        t0=dclock()

        if (if3d) then
           call gs_op_fields(gsh_face,srflx,nxzfl,8,1,1,0)
        else
           call gs_op_fields(gsh_face,srflx,nxzfl,6,1,1,0)
        endif
        call measure_comm(t0)
        
        k = nxzfl

        if (if3d) then
         do i=1,nxzfl
         srflx(0*k+i)= 0.5*srflx(0*k+i)
         srflx(1*k+i)= 0.5*srflx(1*k+i) 
         srflx(2*k+i)= 0.5*srflx(2*k+i) 
         srflx(3*k+i)= 0.5*srflx(3*k+i) 
         srflx(4*k+i)= 0.5*srflx(4*k+i) 
         srflx(5*k+i)= 0.5*srflx(5*k+i) 
         srflx(6*k+i)= 0.5*srflx(6*k+i) 
         srflx(7*k+i)= 0.5*srflx(7*k+i) 
         enddo
        else
         do i=1,nxzfl
         srflx(0*k+i)= 0.5*srflx(0*k+i)
         srflx(1*k+i)= 0.5*srflx(1*k+i) 
         srflx(2*k+i)= 0.5*srflx(2*k+i) 
         srflx(3*k+i)= 0.5*srflx(3*k+i) 
         srflx(4*k+i)= 0.5*srflx(4*k+i) 
         srflx(5*k+i)= 0.5*srflx(5*k+i) 
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
      include 'PARALLEL'
      include 'NEKUSE'
      include 'SCHROD'
 
      real     srflx(1)
      integer  nxyzf,nface,nxz,k,i0,i,j,e,f,ef          
      real     unx0,uny0,unz0
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      k    = nxzfl
      call usersol(RKtime,sQr,sQi,sQi,sUr,sUi,sUi)  !FIXME

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
             srflx(0*k+i)= 0.d0
             srflx(1*k+i)= 0.d0
             srflx(2*k+i)= 2.d0*(srflx(2*k+i)+unx0*sUr(i0))
             srflx(3*k+i)= 2.d0*(srflx(3*k+i)+uny0*sUr(i0))
             srflx(4*k+i)= 2.d0*(srflx(4*k+i)+unz0*sUr(i0))
             srflx(5*k+i)= 2.d0*(srflx(5*k+i)+unx0*sUi(i0))
             srflx(6*k+i)= 2.d0*(srflx(6*k+i)+uny0*sUi(i0))
             srflx(7*k+i)= 2.d0*(srflx(7*k+i)+unz0*sUi(i0))

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
             srflx(2*k+i)= 2.d0*(srflx(2*k+i)+unx0*sUr(i0))
             srflx(3*k+i)= 2.d0*(srflx(3*k+i)+uny0*sUr(i0))
             srflx(4*k+i)= 2.d0*(srflx(4*k+i)+unx0*sUi(i0))
             srflx(5*k+i)= 2.d0*(srflx(5*k+i)+uny0*sUi(i0))

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

      call cem_drift_eig        !call drift_eig_arpack
      call cem_end

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_drift_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      integer       lpts_eig2
      parameter    (lpts_eig2=lpts_eig*lpts_eig)

      common /EIGC/ res(  lpts_eig2)
     $            , ev (  lpts_eig2)
     $            , er (  lpts_eig)
     $            , ei (  lpts_eig)
     $            , rew(4*lpts_eig)
      real          res,ev,er,ei,rew
      integer       k1,k,i,j,jj,nn,nn2,npts2,npts3,npts4,info

      common /scrns/  srflxs(8*2*ldim*lx1*lz1*lelt)
      real            srflxs

      nfield = 2
      npts2= npts*2
      npts3= npts*3
      npts4= npts*4

      if (nfield*npts.gt.lpts_eig) then
         write(6,*) 'lpts_eig',lpts_eig, nfield*npts 
         write(6,*) 'insufficient memory in cem_drift_eig, abort'     
         call exitt
      endif

      if (nid.eq.0) write(6,*) 'start cem_drift_eig'
      nn  = npts*nfield
      nn2 = nn*nn
      call rzero(res,nn2)

      k1 = 0
      do jj=1,nfield
      do j =1,npts

         k1 =  k1+1
         k  = (k1-1)*nn+1

         call rzero(Ur   ,npts )
         call rzero(resUi,npts )
         call rzero(resQr,npts3)

         call rzero(Ui   ,npts )
         call rzero(resUr,npts )
         call rzero(resQi,npts3)

         if (jj.eq.1) Ur(j)= 1.  ! jth col of identity matrix
         if (jj.eq.2) Ui(j)= 1.

         !if (jj.eq.3) Qr(j,1)= 1.
         !if (jj.eq.4) Qr(j,2)= 1.
         !if (jj.eq.5) Qr(j,3)= 1.
         !if (jj.eq.2) Qi(j,1)= 1.
         !if (jj.eq.3) Qi(j,2)= 1.
         !if (jj.eq.4) Qi(j,3)= 1.

         call cem_drift
         call cem_drift_restrict_to_face
         call cem_drift_flux (srflxs)
         call cem_drift_add_flux_to_res (srflxs)
         
         do i=1,npts

            if (abs(resUr(i))  .lt.1.e-13) resUr(i)  =0.
            if (abs(resUi(i))  .lt.1.e-13) resUi(i)  =0.
            if (abs(resQr(i,1)).lt.1.e-13) resQr(i,1)=0.
            if (abs(resQr(i,2)).lt.1.e-13) resQr(i,2)=0.
            if (abs(resQr(i,3)).lt.1.e-13) resQr(i,3)=0.
            if (abs(resQi(i,1)).lt.1.e-13) resQi(i,1)=0.
            if (abs(resQi(i,2)).lt.1.e-13) resQi(i,2)=0.
            if (abs(resQi(i,3)).lt.1.e-13) resQi(i,3)=0.

         enddo

         call copy(res(k+0*npts),resUr,npts)
         call copy(res(k+1*npts),resUi,npts)

        !call copy(res(k+2*npts),resQr(1,1),npts)
        !call copy(res(k+3*npts),resQr(1,2),npts)
        !call copy(res(k+4*npts),resQr(1,3),npts)
        !call copy(res(k+1*npts),resQi(1,1),npts)
        !call copy(res(k+2*npts),resQi(1,2),npts)
        !call copy(res(k+3*npts),resQi(1,3),npts)

      enddo
      enddo

      if (nid.eq.0) write(6,*) 'obtain spatial operator'
      if (nelv.eq.1) call outmat(res,nx1*3,nx1*3,'res   ',1)

      open(50,file='eig_data.dat')
      do i=1,nn2
         write(50,*) res(i)
      enddo
      close(50)

      write(6,*) '                                      '
      write(6,*) 'IFEIG T :: computing eigenvalues      '
      write(6,*) '                                      '

      call dgeev('N','V',nn,res,nn,er,ei,ev,nn,ev,nn,rew,4*nn,info)

      open(60,file='eigenvalue.dat')
      do i=1,nn
        write(60,*)  er(i),ei(i)
      enddo
      close(60)

    1 format(1p1e19.12)

      if (nid.eq.0) write(6,*) 'eigenvalue computation completed!!!'
      call exitt


      return
      end

c-----------------------------------------------------------------------

