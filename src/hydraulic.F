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

c ... 5-stage 4th-order RK
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
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
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
      integer  i
      real     gamm,bbeta,omegap2,bmm,bmv,gammbm,omegap
      real     bmpn1,bmpn2,bmpn3
      real     betaqn1,betaqn2,betaqn3
      real     omegen1,omegen2,omegen3
      real     xx,yy,zz

      common /scrnsh/ srflxh(4*lxzfl10)
      real            srflxh

      integer icalld
      save    icalld
      data    icalld/0/

      if (icalld.eq.0) then

          do i=1,npts
             xx=xm1(i,1,1,1)
             yy=ym1(i,1,1,1)
             zz=zm1(i,1,1,1)
             JN1(i,1)= sin(xx-time)
             JN1(i,2)= cos(yy-time)
             JN1(i,3)= 0.0        
             PN1(i,1)= -cos(xx-time)
             PN1(i,2)= sin(yy-time)
             PN1(i,3)= 0.0           
             QN1(i,1)= cos(xx-time)-sin(yy-time)
             resJN1(i,1)= 0.0        
             resJN1(i,2)= 0.0        
             resJN1(i,3)= 0.0        
             resPN1(i,1)= 0.0        
             resPN1(i,2)= 0.0        
             resPN1(i,3)= 0.0        
             resQN1(i,1)= 0.0        
             resQN1(i,2)= 0.0        
             resQN1(i,3)= 0.0        
          enddo 
          icalld=1

      endif

      call hydraulic_wght_grad
     $  (resQN1(1,1),resQN1(1,2),resQN1(1,3),QN1(1,1),npts)
      call hydraulic_wght_div (resWN1,JN1(1,1),JN1(1,2),JN1(1,3),npts)
      call cem_hydraulic_restrict_to_face     !fHN,fEN,fQN1,fJN1
      call cem_hydraulic_flux(srflxh)
      call cem_hydraulic_add_flux_to_res(srflxh)

      gamm    = 1.0 ! gamma1
      bbeta   = 1.0 ! beta**2 
      omegap2 = 0.0 ! omegap**2

      do i= 1,npts

         bmm     = bm1(i,1,1,1)
         bmv     = 1.0 /bmm
         gammbm  = gamm*bmm
         omegap  = omegap2*permittivity(i)

         bmPN1   = gammbm * PN1(i,1)
         bmPN2   = gammbm * PN1(i,2)
         bmPN3   = gammbm * PN1(i,3)
         betaQN1 = bbeta  * resQN1(i,1) 
         betaQN2 = bbeta  * resQN1(i,2) 
         betaQN3 = bbeta  * resQN1(i,3) 
         omegEN1 = omegap * resEN (i,1)
         omegEN2 = omegap * resEN (i,2)
         omegEN3 = omegap * resEN (i,3)
     
         resPN1(i,1) = -bmPN1 + betaQN1 + omegEN1
         resPN1(i,2) = -bmPN2 + betaQN2 + omegEN2
         resPN1(i,3) = -bmPN3 + betaQN3 + omegEN3

         resJN1(i,1) =  bmPN1
         resJN1(i,2) =  bmPN2                      
         resJN1(i,3) =  bmPN3                 
         QN1   (i,1) = resWN1(i,1)*bmv  ! grad(J)
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

      C0= 0.5 !central flux

c...  fluxes for H and E

c...  fluxes for QN     
      k = nxzfl
      if (if3d) then
         do i= 1,nxzfl
            j= cemface(i)
            srflx(0*k+i)= -unxm(i)*fQN1(i,1)  ! fQN1 
            srflx(1*k+i)= -unym(i)*fQN1(i,1)  ! fQN1 
            srflx(2*k+i)= -unzm(i)*fQN1(i,1)  ! fQN1 
            srflx(3*k+i)= -unxm(i)*fJN1(i,1)  ! fJN1 
     $                    -unym(i)*fJN1(i,2)  ! fJN1 
     $                    -unzm(i)*fJN1(i,3)  ! fJN1 
         enddo
         call gs_op_fields(gsh_face,srflx,nxzfl,4,1,1,0)
         do i=1,nxzfl
            srflx(0*k+i)= 0.5*srflx(0*k+i)
            srflx(1*k+i)= 0.5*srflx(1*k+i)
            srflx(2*k+i)= 0.5*srflx(2*k+i)
            srflx(3*k+i)= 0.5*srflx(3*k+i)
         enddo
      else
         do i= 1,nxzfl
            j= cemface(i)
            srflx(0*k+i)= -unxm(i)*fQN1(i,1)  ! fQN1  
            srflx(1*k+i)= -unym(i)*fQN1(i,1)  ! fQN1  
            srflx(2*k+i)= -unxm(i)*fJN1(i,1)  ! fJN1 
     $                    -unym(i)*fJN1(i,2)  ! fJN1 
         enddo
         call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)
         do i=1,nxzfl
            srflx(0*k+i)= 0.5*srflx(0*k+i)
            srflx(1*k+i)= 0.5*srflx(1*k+i)
            srflx(2*k+i)= 0.5*srflx(2*k+i)
         enddo
      endif

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
           fHN (j,1) = HN (i,1)
           fHN (j,2) = HN (i,2)
           fHN (j,3) = HN (i,3)
           fEN (j,1) = EN (i,1)
           fEN (j,2) = EN (i,2)
           fEN (j,3) = EN (i,3)
           fJN1(j,1) = JN1(i,1)
           fJN1(j,2) = JN1(i,2)
           fJN1(j,3) = JN1(i,3)
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

      if     (imode.eq.3) then !IF3D
        do j = 1,ncemface
           i = cemface(j)
           a = aream(j)
c          resHN (i,1) = resHN (i,1) + a*srflx(0*k+j)
c          resHN (i,2) = resHN (i,2) + a*srflx(1*k+j)
c          resHN (i,3) = resHN (i,3) + a*srflx(2*k+j)
c          resEN (i,1) = resEN (i,1) + a*srflx(3*k+j)
c          resEN (i,2) = resEN (i,2) + a*srflx(4*k+j)
c          resEN (i,3) = resEN (i,3) + a*srflx(5*k+j)
           resQN1(i,1) = resQN1(i,1) + a*srflx(0*k+j)
           resQN1(i,2) = resQN1(i,2) + a*srflx(1*k+j)
           resQN1(i,3) = resQN1(i,3) + a*srflx(2*k+j)
           resWN1(i,1) = resWN1(i,1) + a*srflx(3*k+j)
        enddo
      elseif (imode.eq.2) then !IFTM
        do j=  1,ncemface
           i = cemface(j)
           a = aream(j)
c          resHN (i,1) = resHN (i,1) + a*srflx(0*k+j)
c          resHN (i,2) = resHN (i,2) + a*srflx(1*k+j)
c          resEN (i,3) = resEN (i,3) + a*srflx(2*k+j)
           resQN1(i,1) = resQN1(i,1) + a*srflx(0*k+j)
           resQN1(i,2) = resQN1(i,2) + a*srflx(1*k+j)
           resWN1(i,1) = resWN1(i,1) + a*srflx(2*k+j)
        enddo
      elseif (imode.eq.1) then !IFTE
        do j = 1,ncemface
           i = cemface(j)
           a = aream(j)
c          resEN (i,1) = resEN (i,1) + a*srflx(0*k+j)
c          resEN (i,2) = resEN (i,2) + a*srflx(1*k+j)
c          resHN (i,3) = resHN (i,3) + a*srflx(2*k+j)
           resQN1(i,1) = resQN1(i,1) + a*srflx(0*k+j)
           resQN1(i,2) = resQN1(i,2) + a*srflx(1*k+j)
           resWN1(i,1) = resWN1(i,1) + a*srflx(2*k+j)
        enddo
      endif

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

      call rk4_upd(JN1(1,1),kJN1(1,1),resJN1(1,1),cb,ca,dt,npts)
      call rk4_upd(JN1(1,2),kJN1(1,2),resJN1(1,2),cb,ca,dt,npts)
      call rk4_upd(JN1(1,3),kJN1(1,3),resJN1(1,3),cb,ca,dt,npts)
      call rk4_upd(PN1(1,1),kPN1(1,1),resPN1(1,1),cb,ca,dt,npts)
      call rk4_upd(PN1(1,2),kPN1(1,2),resPN1(1,2),cb,ca,dt,npts)
      call rk4_upd(PN1(1,3),kPN1(1,3),resPN1(1,3),cb,ca,dt,npts)

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

            w3  = (  u3rw*rzmn(k)
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
            u1rw = ur1(i)*w3mn(i)
            u1sw = us1(i)*w3mn(i)
            u2rw = ur2(i)*w3mn(i)
            u2sw = us2(i)*w3mn(i)

            w1   = (  u1rw*rxmn(k)
     $              + u1sw*sxmn(k))

            w2   = (  u2rw*rymn(k)
     $              + u2sw*symn(k))
         
            w0(k) = w1 + w2 
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
