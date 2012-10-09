!> \file schrod.F 
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
      subroutine cem_schrod_temporary                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD' 
      integer  i, npts3

      npts3= 3*npts

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

          call rzero (sQr,npts3)
          call rzero (sQi,npts3)
          call rzero (sUr,npts)
          call rzero (sUi,npts)
          call rzero (Qr,npts3)
          call rzero (Qi,npts3)
          call rzero (Ur,npts)
          call rzero (Ui,npts)
          call rzero (resQr,npts3)
          call rzero (resQi,npts3)
          call rzero (resUr,npts)
          call rzero (resUi,npts)
          call rzero (errQr,npts3)
          call rzero (errQi,npts3)
          call rzero (errUr,npts)
          call rzero (errUi,npts)
          call rzero (kUr,npts)
          call rzero (kUi,npts)

          call rzero (planck,npts)
          call rzero (sqbm1 ,npts)
          call rzero (sqbm2 ,npts)

          call rzero (w3mn,nxyz)
          call rzero (unxm,nxzfl)
          call rzero (unym,nxzfl)
          call rzero (unzm,nxzfl)
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
        call cem_schrod_uvp       !returns planck(i)             
        call cem_schrod_materials !other material related setting

c...  inverse mass matrix including material constants 
        call col3   (sqbm1,planck,bmn,npts)
        call invcol1(sqbm1,npts)
        call invers2(sqbm2,bmn,npts)
        if (nid.eq.0) write(6,*) 'done: inverse mass matrix'

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_schrod_uvp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      integer  ie
      real     planck_max,planck_min
      real     glmax,glmin

c...  nekuvp returns planck(i), using planck_c defined in uservp

      do ie=1,nelt
          call nekuvp(ie)  
      enddo
      planck_max=glmax(planck,npts)
      planck_min=glmin(planck,npts)

      if (nid.eq.0) write(6,*) 'maxmum  planck constant=',planck_max
      if (nid.eq.0) write(6,*) 'minimum planck constant=',planck_min
      if ((planck_max.le.0).or.(planck_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define plank_c in uservp (.usr)'
          call exitt
      endif
                
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_schrod_materials
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      real     xx,yy,zz,rr 
      integer  i,j

c...  compute radius and save the array  
      do i = 1,npts
         xx= xmn(i)
         yy= ymn(i)
         zz= zmn(i)
         rr= xx**2+yy**2
         rad(i) = sqrt(rr)
      enddo

c...  gradient of radius on the surfaces    
      call schrod_grad (wn(1,1),wn(1,2),wn(1,3),rad)

      do j = 1,ncemface
         i = cemface(j)
         frx(j) = wn(i,1)
         fry(j) = wn(i,2)
         frz(j) = wn(i,3)
      enddo

      if (nid.eq.0) write(6,*) 'done: (frx,fry,frz)=grad(rr)'

      ! add stuffs when needed

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_schrod_op_rk
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      integer  ii

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
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      common /scrns/ srflxs(8*2*ldim*lx1*lz1*lelt)        
      real           srflxs

      call cem_schrod
      call cem_schrod_restrict_to_face
      call cem_schrod_flux (srflxs)
      call cem_schrod_add_flux_to_res (srflxs)
      call cem_schrod_source   
      call cem_schrod_invqmass

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_schrod                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      call schrod_wght_grad (resQr(1,1),resQr(1,2),resQr(1,3),Ur)
      call schrod_wght_grad (resQi(1,1),resQi(1,2),resQi(1,3),Ui)
      call schrod_wght_div  (resUi,wn(1,1),wn(1,2),wn(1,3),
     $                             qr(1,1),qr(1,2),qr(1,3))
      call schrod_wght_div  (resUr,wn(1,1),wn(1,2),wn(1,3),
     $                             qi(1,1),qi(1,2),qi(1,3))
      call chsign           (resUi,npts)

      return
      end

c-----------------------------------------------------------------------
      subroutine schrod_grad(w1,w2,w3,u1)
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
      subroutine schrod_wght_grad(w1,w2,w3,u1)
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
      subroutine schrod_wght_div (w0,w1,w2,w3,u1,u2,u3)
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
      subroutine cem_schrod_restrict_to_face
      implicit none
c     Restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      integer  i,j

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
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'SCHROD'

      real     srflx(1), a
      integer  i,j,k

      k   = nxzfl            

      if (if3d) then
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
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
         a = aream(j)
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
      subroutine cem_schrod_source    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'SCHROD'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      integer i

      if (.not.ifsrc) return

      do i=1,npts
         wnx(i)=bmn(i)*Ur(i)
         wny(i)=bmn(i)*Ui(i)
         wnz(i)=bmn(i)*Ui(i)
      enddo
      call usersrc (57,resUr,resUi,resUi,wnx,wny,wnz)
                                                                 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_schrod_invqmass
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'SCHROD'
      integer i

      do i=1,npts
         resUr(i)  =resUr(i)  *sqbm1(i)
         resUi(i)  =resUi(i)  *sqbm1(i)
         resQr(i,1)=resQr(i,1)*sqbm2(i)
         resQr(i,2)=resQr(i,2)*sqbm2(i)
         resQr(i,3)=resQr(i,3)*sqbm2(i)
         resQi(i,1)=resQi(i,1)*sqbm2(i)
         resQi(i,2)=resQi(i,2)*sqbm2(i)
         resQi(i,3)=resQi(i,3)*sqbm2(i)
         Qr(i,1)   =resQr(i,1)           
         Qr(i,2)   =resQr(i,2)
         Qr(i,3)   =resQr(i,3)
         Qi(i,1)   =resQi(i,1)
         Qi(i,2)   =resQi(i,2)
         Qi(i,3)   =resQi(i,3)            
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_schrod_flux (srflx)  
      implicit none
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

      if (icalld.eq.0) then
          call planck_restrict_to_face !obtain fw
          call copy (srflxc,fw,nxzfl)
          call gs_op_fields(gsh_face,srflxc,nxzfl,1,1,1,0)
          do i=1,nxzfl 
             fw(i)=srflxc(i)-2*fw(i)
          enddo

          do e= 1,nelt
          do f= 1,nfaces
          do j= 1,nxzf

          i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering

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

        k = nxzfl
        if (if3d) then
         do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
           unz0=unzm(i)
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

           if (abs(tmp).gt.0.01)  then ! if nonzero case               

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
           grad_r=1.0
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
         do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)

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
           call gs_op_fields(gsh_face,srflx,nxzfl,8,1,1,0)
        else
           call gs_op_fields(gsh_face,srflx,nxzfl,6,1,1,0)
        endif
        call measure_comm(t0)
        
        if (IFPEC) call cem_schrod_flux_pec(srflx)

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
      subroutine cem_schrod_flux_pec(srflx)
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
      subroutine planck_restrict_to_face
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      integer  i,j

      do j = 1,ncemface
         i = cemface(j)
         fW(j)= Planck(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine rk_schrod_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'
      include 'SCHROD'

      integer  ii   
      real     ca,cb
c     kHX=rk4a(i)*kHX

      ca = rk4a(ii)
      cb = rk4b(ii)
      call rk4_upd(Ur,kUr,resUr,cb,ca,dt,npts)
      call rk4_upd(Ui,kUi,resUi,cb,ca,dt,npts)

      return
      end
c----------------------------------------------------------------------
      subroutine cem_schrod_op_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      if (nid.eq.0) write(6,*) 'compute eigenvalues'

      call cem_schrod_eig        !call schrod_eig_arpack
      call cem_end

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_schrod_eig
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
         write(6,*) 'insufficient memory in cem_schrod_eig, abort'     
         call exitt
      endif

      if (nid.eq.0) write(6,*) 'start cem_schrod_eig'
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

         call cem_schrod
         call cem_schrod_restrict_to_face
         call cem_schrod_flux (srflxs)
         call cem_schrod_add_flux_to_res (srflxs)
         
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

