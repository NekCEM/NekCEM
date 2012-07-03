c--------------------------------------------------------------------- 
      subroutine cem_schrod_op_rk
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'
      include 'SCHROD'

      do ii=1,5
         call rk_c (ii)
         call cem_schrod_op
         call rk_schrod_ab(ii)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_schrod_op
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'SCHROD'
      common /scrns/ srflxs(8*2*ldim*lx1*lz1*lelt)        

      call cem_schrod
      call cem_schrod_restrict_to_face
      call cem_schrod_flux (srflxs)
      call cem_schrod_add_flux_to_res (srflxs)
      call cem_schrod_invqmass

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_schrod                
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'
      include 'SCHROD'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)

      COMMON /SCRLB1/ WX1  (LX1*LY1*LZ1*LELT)
     $,               WY1  (LX1*LY1*LZ1*LELT)
     $,               WZ1  (LX1*LY1*LZ1*LELT)
      COMMON /SCRLB2/ WX2  (LX1*LY1*LZ1*LELT)
     $,               WY2  (LX1*LY1*LZ1*LELT)
     $,               WZ2  (LX1*LY1*LZ1*LELT)

      n = nx1*ny1*nz1*nelt

      call schrod_wght_grad (resQr(1,1),resQr(1,2),resQr(1,3),Ur)
      call schrod_wght_grad (resQi(1,1),resQi(1,2),resQi(1,3),Ui)
      call schrod_wght_div  (resUi,wx1,wy1,wz1,qr(1,1),qr(1,2),qr(1,3))
      call schrod_wght_div  (resUr,wx2,wy2,wz2,qi(1,1),qi(1,2),qi(1,3))
      call chsign           (resUi,n)

      return
      end

c-----------------------------------------------------------------------
      subroutine schrod_grad(w1,w2,w3,u1)
c     Non-weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
c
      integer mode, imode
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1

      if (if3d) then

         do e=1,nelt
         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (    u1r(i)*rxm1(k,1,1,1)
     $                 + u1s(i)*sxm1(k,1,1,1)
     $                 + u1t(i)*txm1(k,1,1,1))/jacm1(i,1,1,e)

            w2(k) = (    u1r(i)*rym1(k,1,1,1)
     $                 + u1s(i)*sym1(k,1,1,1)
     $                 + u1t(i)*tym1(k,1,1,1))/jacm1(i,1,1,e)

            w3(k) = (    u1r(i)*rzm1(k,1,1,1)
     $                 + u1s(i)*szm1(k,1,1,1)
     $                 + u1t(i)*tzm1(k,1,1,1))/jacm1(i,1,1,e)

          enddo

          enddo

       else

c        ...2d grad 
          do e= 1,nelt

           j= nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j

             w1(k) = ( u1r(i)*rxm1(k,1,1,1)
     $               + u1s(i)*sxm1(k,1,1,1))/jacm1(i,1,1,e)

             w2(k) =(  u1r(i)*rym1(k,1,1,1)
     $               + u1s(i)*sym1(k,1,1,1))/jacm1(i,1,1,e)
             w3(k)=w2(k)
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine schrod_wght_grad(w1,w2,w3,u1)
c     Weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
c
      integer mode, imode
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1

      if (if3d) then

         do e=1,nelt
         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
c...        rym1,sym1,tym1,...jac are multiplied. w3m1 is purely weight.
c...        area(i,1,1,1) is multiplied with jacobina weight.
            u1rw = u1r(i)*w3m1(i,1,1)
            u1sw = u1s(i)*w3m1(i,1,1)
            u1tw = u1t(i)*w3m1(i,1,1)
            w1(k) =      u1rw*rxm1(k,1,1,1)
     $                 + u1sw*sxm1(k,1,1,1)
     $                 + u1tw*txm1(k,1,1,1)

            w2(k) =      u1rw*rym1(k,1,1,1)
     $                 + u1sw*sym1(k,1,1,1)
     $                 + u1tw*tym1(k,1,1,1)

            w3(k) =      u1rw*rzm1(k,1,1,1)
     $                 + u1sw*szm1(k,1,1,1)
     $                 + u1tw*tzm1(k,1,1,1)

          enddo

          enddo

       else

c        ...2d grad 
          do e= 1,nelt

           j= nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j


             u1rw = u1r(i)*w3m1(i,1,1)
             u1sw = u1s(i)*w3m1(i,1,1)

c             u1rw = u1r(i) !*w3m1(i,1,1)
c             u1sw = u1s(i) !*w3m1(i,1,1)

             w1(k) = ( u1rw*rxm1(k,1,1,1)
     $               + u1sw*sxm1(k,1,1,1))!/jacm1(i,1,1,e)

             w2(k) =( u1rw*rym1(k,1,1,1)
     $               + u1sw*sym1(k,1,1,1))!/jacm1(i,1,1,e)
             w3(k)=w2(k)
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine schrod_wght_div (w0,w1,w2,w3,u1,u2,u3)
c     div:   w = div(u)/jacm1
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer  e

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)

      nxyz  = nx1*ny1*nz1
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

            u1rw = ur1(i)*w3m1(i,1,1)
            u1sw = us1(i)*w3m1(i,1,1)
            u1tw = ut1(i)*w3m1(i,1,1)
            u2rw = ur2(i)*w3m1(i,1,1)
            u2sw = us2(i)*w3m1(i,1,1)
            u2tw = ut2(i)*w3m1(i,1,1)
            u3rw = ur3(i)*w3m1(i,1,1)
            u3sw = us3(i)*w3m1(i,1,1)
            u3tw = ut3(i)*w3m1(i,1,1)

            w1(k) = (  u1rw*rxm1(k,1,1,1)
     $               + u1sw*sxm1(k,1,1,1)
     $               + u1tw*txm1(k,1,1,1))

            w2(k) = (  u2rw*rym1(k,1,1,1)
     $               + u2sw*sym1(k,1,1,1)
     $               + u2tw*tym1(k,1,1,1))

            w3(k) = (  u3rw*rzm1(k,1,1,1)
     $               + u3sw*szm1(k,1,1,1)
     $               + u3tw*tzm1(k,1,1,1))

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
            u1rw = ur1(i)*w3m1(i,1,1)
            u1sw = us1(i)*w3m1(i,1,1)
            u2rw = ur2(i)*w3m1(i,1,1)
            u2sw = us2(i)*w3m1(i,1,1)

            w1(k) = (  u1rw*rxm1(k,1,1,1)
     $               + u1sw*sxm1(k,1,1,1))

            w2(k) = (  u2rw*rym1(k,1,1,1)
     $               + u2sw*sym1(k,1,1,1))
         
            w0(k) = w1(k)+w2(k)       
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_schrod_restrict_to_face
c     Restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) call cem_set_fc_ptr
      icalld = 1

       do j=1,ncemface
          i= cemface(j)
           fQr(j,1) = Qr(i,1)
           fQr(j,2) = Qr(i,2)
           fQr(j,3) = Qr(i,3)
           fQi(j,1) = Qi(i,1)
           fQi(j,2) = Qi(i,2)
           fQi(j,3) = Qi(i,3)
           fUr(j)   = Ur(i)
           fUi(j)   = Ui(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_schrod_add_flux_to_res(srflx)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'SCHROD'

      real srflx(1)

      k= nx1*nz1*2*ndim*nelt

      n = nx1*ny1*nz1*nelt

      if (if3d) then
         do j = 1,ncemface
         i = cemface(j)
         a = area(j,1,1,1)
         resUr(i)   = resUr(i)   + a*srflx(0*k+j)
         resUi(i)   = resUi(i)   + a*srflx(1*k+j)
         resQr(i,1) = resQr(i,1) + a*srflx(2*k+j)
         resQr(i,2) = resQr(i,2) + a*srflx(3*k+j)
         resQr(i,3) = resQr(i,3) + a*srflx(4*k+j)
         resQi(i,1) = resQi(i,1) + a*srflx(5*k+j)
         resQi(i,2) = resQi(i,2) + a*srflx(6*k+j)
         resQi(i,3) = resQi(i,3) + a*srflx(7*k+j)
         enddo
      else
         do j = 1,ncemface
         i = cemface(j)
         a = area(j,1,1,1)
         resUr(i)   = resUr(i)   + a*srflx(0*k+j)
         resUi(i)   = resUi(i)   + a*srflx(1*k+j)
         resQr(i,1) = resQr(i,1) + a*srflx(2*k+j)
         resQr(i,2) = resQr(i,2) + a*srflx(3*k+j)
         resQi(i,1) = resQi(i,1) + a*srflx(4*k+j)
         resQi(i,2) = resQi(i,2) + a*srflx(5*k+j)
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_schrod_invqmass
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'SCHROD'

      integer icalld
      save    icalld
      data    icalld  /0/
      real*8  dclock

      integer n

      n = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call col3   (sqbm1,planck,bm1,n)
          call invcol1(sqbm1,n)
          call invers2(sqbm2,bm1,n)
          icalld = 1
      endif

      do i=1,n
         resUr(i)  =resUr(i)  *sqbm1(i)
         resUi(i)  =resUi(i)  *sqbm1(i)
         resQr(i,1)=resQr(i,1)*sqbm2(i)
         resQr(i,2)=resQr(i,2)*sqbm2(i)
         resQr(i,3)=resQr(i,3)*sqbm2(i)
         resQi(i,1)=resQi(i,1)*sqbm2(i)
         resQi(i,2)=resQi(i,2)*sqbm2(i)
         resQi(i,3)=resQi(i,3)*sqbm2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_schrod_flux (srflx)  

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'SCHROD'
     
      real*8  dclock,t0
      real    srflx(1)
      integer nxyzf,nxz,nface
      integer i,j,k,e,f

      integer icalld
      save    icalld
      data    icalld/0/
      character CB*3

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1
      nface = 2*ndim

      if (icalld.eq.0) then
          call planck_restrict_to_face !obtain fw
          call copy (srflxc,fw,nxyzf)
          call gs_op_fields(gsh_face,srflxc,nxyzf,1,1,1,0)
          do i=1,nxyzf
             fw(i)=srflxc(i)-2*fw(i)
          enddo

          do e= 1,nelt
          do f= 1,nface
          do j= 1,nxz

          i = (e-1)*nxz*nface + nxz*(f-1) + j  ! face numbering

          CB = CBC(f,e,2)
          if (CB.eq.'E  ')then
          else 
             fw(i)=0.0
          endif 
          enddo
          enddo
          enddo

          icalld=1
      endif

      call usersol(RKtime,sQr,sQi,sQi,sUr,sUi,sUi) !FIXME

        k = nxyzf
        if (if3d) then
         do i=1,nxyzf                
           unx0=unx(i,1,1,1)
           uny0=uny(i,1,1,1)
           unz0=unz(i,1,1,1)
           srflx(0*k+i)= -(unx0*fQi(i,1)+uny0*fQi(i,2)+unz0*fQi(i,3)) 
           srflx(1*k+i)=  (unx0*fQr(i,1)+uny0*fQr(i,2)+unz0*fQr(i,3)) 
           srflx(2*k+i)= - unx0*fUr(i) 
           srflx(3*k+i)= - uny0*fUr(i) 
           srflx(4*k+i)= - unz0*fUr(i) 
           srflx(5*k+i)= - unx0*fUi(i) 
           srflx(6*k+i)= - uny0*fUi(i) 
           srflx(7*k+i)= - unz0*fUi(i) 

c**********this is for 1d discontinous case
c           tmp1=fw(i)*0.5*(RKtime+7.0)
c           tmp1=unx(i,1,1,1)*(cos(tmp1)-1.0)
c           tmp2=sin(fw(i)*0.5*(RKtime+7.0))
c           tmp2=unx(i,1,1,1)*tmp2
c           srflx(0*k+i)=srflx(0*k+i)
c     $                  -1.0*(fQi(i,1,1,1)*tmp1+fQr(i,1,1,1)*tmp2)
c           srflx(1*k+i)=srflx(1*k+i)
c     $                  +1.0*(fQr(i,1,1,1)*tmp1-fQi(i,1,1,1)*tmp2)
c           srflx(2*k+i)=srflx(2*k+i)
c     $                  -1.0*(fUr(i,1,1)  *tmp1-fUi(i,1,1)  *tmp2)
c           srflx(5*k+i)=srflx(5*k+i)
c     $                  -1.0*(fUi(i,1,1)  *tmp1+fUr(i,1,1)  *tmp2)
c**********this is for 3d cylinder discontinous case

           i0 = cemface(i)  ! i0= global numbering per core
           tmp= fw(i)  

           if (abs(tmp).gt.0)  then ! if nonzero case               

           tmpi1= sQi(i0,1)
           tmpi2= sQi(i0,2)
           tmpi3= sQi(i0,3)
           tmpr1= sQr(i0,1)
           tmpr2= sQr(i0,2)
           tmpr3= sQr(i0,3)
           tmpUr= sUr(i0)
           tmpUi= sUi(i0)

           grad_r = frx(i)**2+fry(i)**2 ! tmp3=1 if circle
           grad_r = sqrt(grad_r)

           fjump_condi  = grad_r*(unx0*tmpi1+uny0*tmpi2+unz0*tmpi3)
           fjump_condr  = grad_r*(unx0*tmpr1+uny0*tmpr2+unz0*tmpr3)
           fjump_condr1 = unx0*tmpUr
           fjump_condr2 = uny0*tmpUr
           fjump_condr3 = unz0*tmpUr
           fjump_condi1 = unx0*tmpUi
           fjump_condi2 = uny0*tmpUi
           fjump_condi3 = unz0*tmpUi

           srflx(0*k+i)= srflx(0*k+i) + fjump_condi
           srflx(1*k+i)= srflx(1*k+i) - fjump_condr
           srflx(2*k+i)= srflx(2*k+i) + fjump_condr1
           srflx(3*k+i)= srflx(3*k+i) + fjump_condr2
           srflx(4*k+i)= srflx(4*k+i) + fjump_condr3
           srflx(5*k+i)= srflx(5*k+i) + fjump_condi1
           srflx(6*k+i)= srflx(6*k+i) + fjump_condi2
           srflx(7*k+i)= srflx(7*k+i) + fjump_condi3

           endif
         enddo

        else
         do i=1,nxyzf                
           unx0=unx(i,1,1,1)
           uny0=uny(i,1,1,1)

           srflx(0*k+i)= -(unx0*fQi(i,1)+uny0*fQi(i,2)) 
           srflx(1*k+i)= +(unx0*fQr(i,1)+uny0*fQr(i,2)) 
           srflx(2*k+i)= - unx0*fUr(i)  ! NxfHx
           srflx(3*k+i)= - uny0*fUr(i)  ! NxfHx
           srflx(4*k+i)= - unx0*fUi(i)  ! NxfHx
           srflx(5*k+i)= - uny0*fUi(i)  ! NxfHx

c          tmp1=fw(i)*0.5*(RKtime+7.0)
c          tmp1=unx(i,1,1,1)*(cos(tmp1)-1.0)
c          tmp2=sin(fw(i)*0.5*(RKtime+7.0))
c          tmp2=unx(i,1,1,1)*tmp2
 
           i0 = cemface(i)  ! i0= global numbering per core
           tmp= fw(i)

           if (abs(tmp).gt.0)  then ! if nonzero case               


           tmpi1= sQi(i0,1)
           tmpi2= sQi(i0,2)
           tmpr1= sQr(i0,1)
           tmpr2= sQr(i0,2)
           tmpUr= sUr(i0)
           tmpUi= sUi(i0)

           grad_r = frx(i)**2+fry(i)**2 ! tmp3=1 if circle
           grad_r = sqrt(grad_r)

           fjump_condi  = grad_r*(unx0*tmpi1+uny0*tmpi2)
           fjump_condr  = grad_r*(unx0*tmpr1+uny0*tmpr2)
           fjump_condr1 = unx0*tmpUr
           fjump_condr2 = uny0*tmpUr
           fjump_condi1 = unx0*tmpUi
           fjump_condi2 = uny0*tmpUi

           srflx(0*k+i)= srflx(0*k+i) + fjump_condi
           srflx(1*k+i)= srflx(1*k+i) - fjump_condr
           srflx(2*k+i)= srflx(2*k+i) + fjump_condr1
           srflx(3*k+i)= srflx(3*k+i) + fjump_condr2
           srflx(4*k+i)= srflx(4*k+i) + fjump_condi1
           srflx(5*k+i)= srflx(5*k+i) + fjump_condi2

           endif
         enddo
        endif

        t0=dclock()

        if (if3d) then
           call gs_op_fields(gsh_face,srflx,nxyzf,8,1,1,0)
        else
           call gs_op_fields(gsh_face,srflx,nxyzf,6,1,1,0)
        endif
        call measure_comm(t0)
        
        if (IFPEC) call cem_schrod_flux_pec(srflx)

        k = nxyzf

        if (if3d) then
         do i=1,nxyzf
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
         do i=1,nxyzf
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
      subroutine cem_schrod_flux_pec(srflx)                              

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
      integer  nxyzf,nface,nxz,k,i,e,f,ef          
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      nxyzf= nx1*nz1*2*ndim*nelt
      nface= 2*ndim
      nxz  = nx1*nz1
      k    = nxyzf

      call usersol(RKtime,sQr,sQi,sQi,sUr,sUi,sUi)  !FIXME

      if (if3d) then
         do e= 1,nelt
         do f= 1,nface
         do j= 1,nxz

         i = (e-1)*nxz*nface + nxz*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,2)
         if (CB.eq.'PEC') then

             unx0=unx(i,1,1,1)
             uny0=uny(i,1,1,1)
             unz0=unz(i,1,1,1)
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
         do f= 1,nface
         do j= 1,nxz

         i = (e-1)*nxz*nface + nxz*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,2)
         if (CB.eq.'PEC') then

             unx0=unx(i,1,1,1)
             uny0=uny(i,1,1,1)
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
      subroutine planck_restrict_to_face
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      do j = 1,ncemface
         i = cemface(j)
         fW(j)= Planck(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine rk_schrod_ab(ii)
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'
      include 'SCHROD'

      integer istep

      n=nx1*ny1*nz1*nelt

c     kHX=rk4a(i)*kHX

      ca = rk4a(ii)
      cb = rk4b(ii)
      call rk4_upd(Ur,kUr,resUr,cb,ca,dt,n)
      call rk4_upd(Ui,kUi,resUi,cb,ca,dt,n)

      return
      end

