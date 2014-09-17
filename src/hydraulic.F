c-----------------------------------------------------------------------
      subroutine cem_hydraulic_op_rk
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  ii

      real*8  dclock
      real    stime0,stime5

      stime5=dclock()
      do ii=1,5

         stime0=dclock()
         call rk_c (ii)
         call cem_hydraulic_op
         call rk_hydraulic_ab(ii)
         stime0=dclock()-stime0

      enddo
      stime5=dclock()-stime5

#ifdef TIMER
      if  (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if  (nid.eq.0) then
          write(6,*) '----------------------------------::'
          write(6,*) 'cem_hydraulic_op_rk            (sec)::',stime0
          write(6,*) 'cem_hydraulic_op_rk per istep  (sec)::',stime5
          write(6,*) '----------------------------------::'
      endif
      endif
#endif 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_hydraulic_op 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'    
      integer  i
      real     bgamma
      real     bmpn1,bmpn2,bmpn3
      real     betaqn1,betaqn2,betaqn3
      real     omegen1,omegen2,omegen3
      real     xx,yy,zz
      real     tmp1,tmp2,tmp3    
      real     hydro_beta2, hydro_gamma2, hydro_omega2

      integer  icalld
      save     icalld
      data     icalld/0/

      common /scrnsh/ srflxh(4*lxzfl10)
      real            srflxh

      
      call hydraulic_div(QN1(1,1),JN1(1,1),JN1(1,2),JN1(1,3))
      call hydraulic_wght_grad
     $     (resPN1(1,1),resPN1(1,2),resPN1(1,3),QN1(1,1))
      call hydraulic_wght_div(resQN1(1,1),JN1(1,1),JN1(1,2),JN1(1,3))
      call cem_hydraulic_restrict_to_face  !fHN,fEN,fQN1,fJN1
      call cem_hydraulic_flux(srflxh)
      call cem_hydraulic_add_flux_to_res(srflxh)
      call cem_hydraulic_model                    
      call cem_hydraulic_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_hydraulic_model    
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'    
      integer  i
      real     bgamma
      real     bmpn1,bmpn2,bmpn3
      real     betaqn1,betaqn2,betaqn3
      real     omegen1,omegen2,omegen3
      real     xx,yy,zz
      real     tmp1,tmp2,tmp3    
      real     hydro_beta2, hydro_gamma2, hydro_omega2

      integer  icalld
      save     icalld
      data     icalld/0/

      common /scrnsh/ srflxh(4*lxzfl10)
      real            srflxh

      if (icalld.eq.0) then
          if (nid.eq.0) write(6,*) 'hydro_beta =', hydro_beta
          if (nid.eq.0) write(6,*) 'hydro_gamma=', hydro_gamma
          if (nid.eq.0) write(6,*) 'hydro_omega=', hydro_omega
          icalld = 1
      endif
      
      hydro_beta2  = hydro_beta**2

      do i= 1,npts

         xx      = xmn(i)
         yy      = ymn(i)
         zz      = zmn(i)
         bgamma  = bmn(i) * hydro_gamma
         hydro_omega2 = permittivity(i)*hydro_omega**2

         bmPN1   = bgamma * PN1(i,1)
         bmPN2   = bgamma * PN1(i,2)
         bmPN3   = bgamma * PN1(i,3)
         betaQN1 = hydro_beta2  * resPN1(i,1) 
         betaQN2 = hydro_beta2  * resPN1(i,2) 
         betaQN3 = hydro_beta2  * resPN1(i,3) 
         omegEN1 = hydro_omega2 * resEN (i,1)
         omegEN2 = hydro_omega2 * resEN (i,2)
         omegEN3 = hydro_omega2 * resEN (i,3)
     
         resPN1(i,1) = -bmPN1 + betaQN1 + omegEN1
         resPN1(i,2) = -bmPN2 + betaQN2 + omegEN2
         resPN1(i,3) = -bmPN3 + betaQN3 + omegEN3

         resJN1(i,1) =  PN1(i,1)*bmn(i)
         resJN1(i,2) =  PN1(i,2)*bmn(i)
         resJN1(i,3) =  PN1(i,3)*bmn(i)

         resPN1(i,1) =   betaQN1 
         resPN1(i,2) =   betaQN2
         resPN1(i,3) =   betaQN3

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_hydraulic_flux (srflx)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

      real*8  dclock,t0
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k,e,f

      C0= 1.0 !central flux

c...  fluxes for QN     
      k = nxzfl
      do i= 1,nxzfl
         j= cemface(i)
         srflx(0*k+i)= -unx(i,1,1,1)*fJN1(i,1)-uny(i,1,1,1)*fJN1(i,2) !resQN=resQN+div(J)
         srflx(1*k+i)= -unx(i,1,1,1)*fPN1(i,1)  ! resPN=resPN+grad(Q)
         srflx(2*k+i)= -uny(i,1,1,1)*fPN1(i,2)  ! resPN=resPN+grad(Q)
      enddo

      call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)
      call cem_hydraulic_flux_pec(srflx)

      do i=1,nxzfl
         srflx(0*k+i)= C0 *srflx(0*k+i)
         srflx(1*k+i)= 0.0*srflx(1*k+i)
         srflx(2*k+i)= 0.0*srflx(2*k+i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_hydraulic_restrict_to_face
c-----------------------------------------------------------------------
c     restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

c...  sacrify efficiency for 2D for simplicity
      do j=1, ncemface
           i= cemface(j)
           fJN1(j,1) = JN1(i,1)
           fJN1(j,2) = JN1(i,2)
           fQN1(j,1) = QN1(i,1)
      enddo
     
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_hydraulic_add_flux_to_res(srflx)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  i,j,k
      real     srflx(1),a

      k = nxzfl

      do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resQN1(i,1) = resQN1(i,1) + a*srflx(0*k+j)
         resPN1(i,1) = resPN1(i,1) + a*srflx(1*k+j)
         resPN1(i,2) = resPN1(i,2) + a*srflx(2*k+j)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_hydraulic_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'
      include 'PML'
      real     ca,cb
      integer  ii
      real*8 dclock
      real   stime0

      stime0=dclock()

c     kHX=rk4a(i)*kHX
      ca = rk4a(ii)
      cb = rk4b(ii)

      if (if3d) then
          call rk4_upd(JN1(1,1),kJN1(1,1),resJN1(1,1),cb,ca,dt,npts)
          call rk4_upd(JN1(1,2),kJN1(1,2),resJN1(1,2),cb,ca,dt,npts)
          call rk4_upd(JN1(1,3),kJN1(1,3),resJN1(1,3),cb,ca,dt,npts)
          call rk4_upd(PN1(1,1),kPN1(1,1),resPN1(1,1),cb,ca,dt,npts)
          call rk4_upd(PN1(1,2),kPN1(1,2),resPN1(1,2),cb,ca,dt,npts)
          call rk4_upd(PN1(1,3),kPN1(1,3),resPN1(1,3),cb,ca,dt,npts)
      else
          call rk4_upd(JN1(1,1),kJN1(1,1),resJN1(1,1),cb,ca,dt,npts)
          call rk4_upd(JN1(1,2),kJN1(1,2),resJN1(1,2),cb,ca,dt,npts)
          call rk4_upd(PN1(1,1),kPN1(1,1),resPN1(1,1),cb,ca,dt,npts)
          call rk4_upd(PN1(1,2),kPN1(1,2),resPN1(1,2),cb,ca,dt,npts)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hydraulic_wght_grad(w1,w2,w3,u1)
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


             w1(k) = ( u1rw*rxmn(k)
     $               + u1sw*sxmn(k))

             w2(k) =(  u1rw*rymn(k)
     $               + u1sw*symn(k))
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hydraulic_wght_div (w0,u1,u2,u3)
      implicit none
c     div:   w = div(u)/jacm
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),u1(1),u2(1),u3(1)
      real     w1,w2,w3
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

            w1   = (  u1rw*rxmn(k)
     $              + u1sw*sxmn(k)
     $              + u1tw*txmn(k))

            w2   = (  u2rw*rymn(k)
     $              + u2sw*symn(k)
     $              + u2tw*tymn(k))

            w3  = (   u3rw*rzmn(k)
     $              + u3sw*szmn(k)
     $              + u3tw*tzmn(k))

            w0(k) = w1 + w2 + w3
         enddo

         enddo

      else

         do e = 1,nelt

            j = nxyz*(e-1)

            call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
            call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

         do i = 1,nxyz

            k = i+j
            u1rw = ur1(i)*w3mn(i)
            u1sw = us1(i)*w3mn(i)
            u2rw = ur2(i)*w3mn(i)
            u2sw = us2(i)*w3mn(i)

            w1   = ( u1rw*rxmn(k)
     $           +   u1sw*sxmn(k) )

            w2   = ( u2rw*rymn(k)
     $           +   u2sw*symn(k) )
         
            w0(k) =  w1 + w2 
         enddo

         enddo
         

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hydraulic_div (w0,u1,u2,u3)
      implicit none
c     div:   w = div(u)/jacm
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),u1(1),u2(1),u3(1)
      real     w1,w2,w3
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

            u1rw = ur1(i)/jacm(k)
            u1sw = us1(i)/jacm(k)
            u1tw = ut1(i)/jacm(k)
            u2rw = ur2(i)/jacm(k)
            u2sw = us2(i)/jacm(k)
            u2tw = ut2(i)/jacm(k)
            u3rw = ur3(i)/jacm(k)
            u3sw = us3(i)/jacm(k)
            u3tw = ut3(i)/jacm(k)

            w1   = (  u1rw*rxmn(k)
     $              + u1sw*sxmn(k)
     $              + u1tw*txmn(k))

            w2   = (  u2rw*rymn(k)
     $              + u2sw*symn(k)
     $              + u2tw*tymn(k))

            w3  = (   u3rw*rzmn(k)
     $              + u3sw*szmn(k)
     $              + u3tw*tzmn(k))

            w0(k) = w1 + w2 + w3
         enddo

         enddo

      else

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
         call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            u1rw = ur1(i)/jacm(k)
            u1sw = us1(i)/jacm(k)
            u2rw = ur2(i)/jacm(k)
            u2sw = us2(i)/jacm(k)

            w1 = (  u1rw *rxmn(k)
     $            + u1sw *sxmn(k) )
            w2 = (  u2rw *rymn(k)
     $            + u2sw *symn(k) )
         
            w0(k) =  w1 + w2 
         enddo

         enddo

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine cem_hydraulic_flux_pec(srflx)                          
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
 
      real      srflx(1)
      integer   k,i,j,l,e,f,ef          
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      common  /BDRY1/ cempec(6*lx1*lz1*2*ldim*lelt)
      common  /BDRY2/ ncempec
      integer  cempec,ncempec

      k = nxzfl
      l = 0

      if (icalld.eq.0) then

         do e= 1,nelt            
         do f= 1,nfaces
         do j= 1,nxzf
            CB= CBC(f,e,1)
            i =(e-1)*nfaces*nxzf+(f-1)*nxzf+j
            if ((CB.eq.'PEC').or.(CB.eq.'PML')) then
              l= l+1
              cempec(l)= i
              ncempec  = l
            endif
         enddo
         enddo
         enddo
         icalld =1

      endif

      if (ncempec.eq.0) return
 
      k = nxzfl
      do j = 1,ncempec         
         i = cempec(j)
         srflx(0*k+i)=0.0*srflx(0*k+i)
         srflx(1*k+i)=0.0*srflx(1*k+i)
c        srflx(2*k+i)=0.0*srflx(2*k+i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_hydraulic_invqmass
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer  i

      do i=1,npts
         resJN1(i,1) = resJN1(i,1)*bminv(i)
         resJN1(i,2) = resJN1(i,2)*bminv(i)
         resJN1(i,3) = resJN1(i,3)*bminv(i)
         resPN1(i,1) = resPN1(i,1)*bminv(i)
         resPN1(i,2) = resPN1(i,2)*bminv(i)
         resPN1(i,3) = resPN1(i,3)*bminv(i)
         resQN1(i,1) = resQN1(i,1)*bminv(i)
         QN1   (i,1) = resQN1(i,1)
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine cem_hydraulic_op_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'

      if (nid.eq.0) write(6,*) 'compute eigenvalues: hydraulic'

      call cem_hydraulic_eig        
      call cem_end

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_hydraulic_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer       lpts_eig2
      parameter    (lpts_eig2=lpts_eig*lpts_eig)

      common /EIGC/ res( lpts_eig2)
     $            , ev ( lpts_eig2)
     $            , er ( lpts_eig )
     $            , ei ( lpts_eig )
     $            , rew(4*lpts_eig)
      real          res,ev,er,ei,rew
      integer       k1,k,i,ii,j,jj,nn,nn2,npts2,npts3,npts4,info

      real    glmax,glmin
      real    eigmax_r,eigmin_r,eigmax_i,eigmin_i
      real    tol,xx,yy

      common /scrnsh/ srflxh(4*lxzfl10)
      real            srflxh

      nfield = 4
      npts2= npts*2
      npts3= npts*3
      npts4= npts*4

      if (nfield*npts.gt.lpts_eig) then
         write(6,*) 'lpts_eig',lpts_eig, nfield*npts 
         write(6,*) 'insufficient memory in cem_hydraulic_eig, abort' 
         call exitt
      endif

      if (nid.eq.0) write(6,*) 'start cem_hydraulic_eig'
      nn  = npts*nfield
      nn2 = nn*nn
      call rzero(res,nn2)

      k1 = 0
      do jj = 1,nfield
      do j  = 1,npts

         k1 =  k1+1
         k  = (k1-1)*nn+1

         call rzero(PN1   ,npts3)
         call rzero(JN1   ,npts3)
         call rzero(QN1   ,npts3)
         call rzero(resPN1,npts3)
         call rzero(resJN1,npts3)
         call rzero(resQN1,npts3)

         if (jj.eq.1) JN1(j,1)= 1.  ! jth col of identity matrix
         if (jj.eq.2) JN1(j,2)= 1.
         if (jj.eq.3) PN1(j,1)= 1.
         if (jj.eq.4) PN1(j,2)= 1.

         call hydraulic_wght_div(resQN1(1,1),JN1(1,1),JN1(1,2),JN1(1,3))

         call cem_hydraulic_restrict_to_face         ! fQN,fJN
         call cem_hydraulic_flux(srflxh)             ! resPN=fQN, resQN=fJN
         call cem_hydraulic_add_flux_to_res(srflxh)  ! resPN,resQN
         call cem_hydraulic_invqmass                 ! added
         call copy (QN1(1,1),resQN1,npts)                

c        call rzero(resQN1,npts3)
         call rzero(resPN1,npts3)
         call rzero(resJN1,npts3)

         call hydraulic_wght_grad
     $        (resPN1(1,1),resPN1(1,2),resPN1(1,3),QN1(1,1))
C        call hydraulic_wght_div(resQN1(1,1),JN1(1,1),JN1(1,2),JN1(1,3))
         call cem_hydraulic_restrict_to_face  !fHN,fEN,fQN1,fJN1
         call cem_hydraulic_flux(srflxh)
         call cem_hydraulic_add_flux_to_res(srflxh)
         call cem_hydraulic_model  
         call cem_hydraulic_invqmass ! added
         
         do ii = 1,npts
 
            tol= 1.e-13
            if (abs(resPN1(ii,1)).lt.tol) resPN1(ii,1)=0.
            if (abs(resPN1(ii,2)).lt.tol) resPN1(ii,2)=0.
            if (abs(resJN1(ii,1)).lt.tol) resJN1(ii,1)=0.
            if (abs(resJN1(ii,2)).lt.tol) resJN1(ii,2)=0.
            if (abs(resQN1(ii,1)).lt.tol) resQN1(ii,1)=0.
 
            write(6,22) jj,j,ii,resPN1(ii,1),resPN1(ii,2),
     $                          resJN1(ii,1),resJN1(ii,2),resQN1(ii,1)  
         enddo
   22   format('tt=',3i5,5e13.5)

         call copy(res(k+0*npts),resPN1(1,1),npts)
         call copy(res(k+1*npts),resPN1(1,2),npts)
         call copy(res(k+2*npts),resJN1(1,1),npts)
         call copy(res(k+3*npts),resJN1(1,2),npts)

      enddo
      enddo

      if (nid.eq.0) write(6,*) 'obtain spatial operator'
      if (nelv.eq.1) call outmat(res,nx1*3,nx1*3,'res   ',1)

      open(50,file='spatial_matrix.dat')
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

      eigmax_r= glmax(er,nn)
      eigmin_r= glmin(er,nn)
      eigmax_i= glmax(ei,nn)
      eigmin_i= glmin(ei,nn)

      write(6,*) 'eigenvalue: max(real)=',eigmax_r
      write(6,*) 'eigenvalue: min(real)=',eigmin_r
      write(6,*) 'eigenvalue: max(imag)=',eigmax_i
      write(6,*) 'eigenvalue: min(imag)=',eigmin_i

    1 format(1p1e19.12)

      if (nid.eq.0) write(6,*) 'eigenvalue computation completed!!!'
      call exitt


      return
      end

c-----------------------------------------------------------------------
