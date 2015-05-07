!> \file cem_dg.F 
!!
!! \brief CEM main subroutines          
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
      subroutine cem (igeom)           
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  igeom

      if     (IFSCHROD ) then

          if (ifrk ) call cem_schrod_op_rk
          if (ifeig) call cem_schrod_op_eig

      elseif (IFDRIFT  ) then

      !    if (ifrk ) call cem_drift_op_rk
          if (ifrk ) call cem_drift_dg_op_rk
          if (ifbdf) call cem_drift_op_bdf
          if (ifeig) call cem_drift_op_eig

      elseif (IFHYDRO  ) then

          if (ifrk ) call cem_hydraulic_op_rk
          if (ifeig) call cem_hydraulic_op_eig

      else !  Maxwell Solver
          if (ifrk ) call cem_maxwell_op_rk
          if (ifeig) call cem_maxwell_op_eig
          if (ifexp) call cem_maxwell_op_exp

      endif

      return
      end
C----------------------------------------------------------------------
      subroutine cem_drude 
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer i
   
c...  simpler with explicit writing than calling mat1 routines many times
      do i=1,npts
        !cdrude1 = DRUDEA(i)!*bm1(i,1,1,1)
        !cdrude2 = DRUDEB(i)!*bm1(i,1,1,1)
 
        !resJN(i,1)= cdrude1*JN(i,1)+cdrude2*EN(i,1)
        !resJN(i,2)= cdrude1*JN(i,2)+cdrude2*EN(i,2)
        !resJN(i,3)= cdrude1*JN(i,3)+cdrude2*EN(i,3)
 
        !resJX(i,1,1,1,1)= resJx(i,1,1,1,1)+incEx(i,1,1,1)
        !resJY(i,1,1,1,1)= resJy(i,1,1,1,1)+incEy(i,1,1,1)
        !resJZ(i,1,1,1,1)= resJz(i,1,1,1,1)+incEz(i,1,1,1)

        !resEN(i,1) = resEN(i,1)-JN(i,1)*bmn(i)
        !resEN(i,2) = resEN(i,2)-JN(i,2)*bmn(i)
        !resEN(i,3) = resEN(i,3)-JN(i,3)*bmn(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_lorentz
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer i

      do i=1,npts

c        clorenA1 = LORENA(i,1)!*bm1(i,1,1,1)
c        clorenB1 = LORENB(i,1)!*bm1(i,1,1,1)
c        clorenC1 = LORENC(i,1)!*bm1(i,1,1,1)
c
c        clorenA1_Jx= clorenA1*JN(i,1)                
c        clorenA1_Jy= clorenA1*JN(i,2)                
c        clorenA1_Jz= clorenA1*JN(i,3)                
c
c        clorenB1_Px= clorenB1*PN(i,1)                
c        clorenB1_Py= clorenB1*PN(i,2)                
c        clorenB1_Pz= clorenB1*PN(i,3)                
c
c        clorenC1_Ex= clorenC1*EN(i)                
c        clorenC1_Ey= clorenC1*EN(i)                
c        clorenC1_Ez= clorenC1*EN(i)                
c
c        resJN(i,1)= clorenA1_Jx+clorenB1_Px+clorenC1_Ex 
c        resJN(i,2)= c1orenA1_Jy+clorenB1_Py+clorenC1_Ey
c        resJN(i,3)= clorenA1_Jz+clorenB1_Pz+clorenC1_Ez
c
c        clorenA2 = LORENA(i,2)!*bm1(i,1,1,1)
c        clorenB2 = LORENB(i,2)!*bm1(i,1,1,1)
c        clorenC2 = LORENC(i,2)!*bm1(i,1,1,1)
c
c        clorenA2_Jx= clorenA2*JN(i,2)                
c        clorenA2_Jy= clorenA2*JN(i,2)                
c        clorenA2_Jz= clorenA2*JN(i,2)                
c
c        clorenB2_Px= clorenB2*PN(i,2)                
c        clorenB2_Py= clorenB2*PN(i,2)                
c        clorenB2_Pz= clorenB2*PN(i,2)                
c
c        clorenC2_Ex= clorenC2*EN(i)                
c        clorenC2_Ey= clorenC2*EN(i)                
c        clorenC2_Ez= clorenC2*EN(i)                
c
c        resJx(i,2)= clorenA2_Jx+clorenB2_Px+clorenC2_Ex 
c        resJy(i,2)= c1orenA2_Jy+clorenB2_Py+clorenC2_Ey
c        resJz(i,2)= clorenA2_Jz+clorenB2_Pz+clorenC2_Ez
c
c        clorenA3 = LORENA(i,3)!*bm1(i,1,1,1)
c        clorenB3 = LORENB(i,3)!*bm1(i,1,1,1)
c        clorenC3 = LORENC(i,3)!*bm1(i,1,1,1)

c        clorenA3_Jx= clorenA3*JN(i,3)                
c        clorenA3_Jy= clorenA3*JN(i,3)                
c        clorenA3_Jz= clorenA3*JN(i,3)                
c
c        clorenB3_Px= clorenB3*PN(i,3)                
c        clorenB3_Py= clorenB3*PN(i,3)                
c        clorenB3_Pz= clorenB3*PN(i,3)                
c
c        clorenC3_Ex= clorenC3*EN(i)                
c        clorenC3_Ey= clorenC3*EN(i)                
c        clorenC3_Ez= clorenC3*EN(i)                
c
c        resJN(i,3)= clorenA3_Jx+clorenB3_Px+clorenC3_Ex 
c        resJN(i,3)= c1orenA3_Jy+clorenB3_Py+clorenC3_Ey
c        resJN(i,3)= clorenA3_Jz+clorenB3_Pz+clorenC3_Ez
c
c        cJx=(JN(i,1)+JN(i,2)+JN(i,3))*bmn(i)
c        cJy=(JN(i,1)+JN(i,2)+JN(i,3))*bmn(i)
c        cJz=(JN(i,1)+JN(i,2)+JN(i,3))*bmn(i)
c
c        resEN(i)= resEN(i)-cJx                 
c        resEN(i)= resEN(i)-cJy                 
c        resEN(i)= resEN(i)-cJz            
c
c        resPN(i,1)= JN(i,1)                     
c        resPN(i,1)= JN(i,1)                     
c        resPN(i,1)= JN(i,1)                
c        resPN(i,2)= JN(i,2)                     
c        resPN(i,2)= JN(i,2)                     
c        resPN(i,2)= JN(i,2)                
c        resPN(i,3)= JN(i,3)                     
c        resPN(i,3)= JN(i,3)                     
c        resPN(i,3)= JN(i,3)                

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_c(i)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      integer  i

      RKtime=time+dt*rk4c(i)


      return
      end
c-----------------------------------------------------------------------
      subroutine rk4_upd(h,kh,resh,cb,ca,dt,n)
c-----------------------------------------------------------------------
      implicit none
      real     h(1),kh(1),resh(1)
      real     cb,ca,dt
      integer  n,nresid,iand,i

#ifdef _OPENACC
!$ACC DATA PRESENT(kh,resh,h)
!$ACC PARALLEL LOOP GANG VECTOR
      do i=1,n          
         kh(i  ) = ca*kh(i  ) + dt*resh(i  )
         h (i  ) =     h(i  ) + cb*kh  (i  )
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
#else
c     below routine calculates 4 consecutive ones at one time
      nresid = iand(n,3)

      do i=1,n-nresid,4

         kh(i  ) = ca*kh(i  ) + dt*resh(i  )
         kh(i+1) = ca*kh(i+1) + dt*resh(i+1)
         kh(i+2) = ca*kh(i+2) + dt*resh(i+2)
         kh(i+3) = ca*kh(i+3) + dt*resh(i+3)

         h (i  ) =     h(i  ) + cb*kh  (i  )
         h (i+1) =     h(i+1) + cb*kh  (i+1)
         h (i+2) =     h(i+2) + cb*kh  (i+2)
         h (i+3) =     h(i+3) + cb*kh  (i+3)

      enddo

      if (nresid.eq.0) then

         return

      elseif (nresid.eq.1) then
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n  ) =     h(n  ) + cb*kh  (n  )
      elseif (nresid.eq.2) then
         kh(n-1) = ca*kh(n-1) + dt*resh(n-1)
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n-1) =     h(n-1) + cb*kh  (n-1)
         h (n  ) =     h(n  ) + cb*kh  (n  )
      else
         kh(n-2) = ca*kh(n-2) + dt*resh(n-2)
         kh(n-1) = ca*kh(n-1) + dt*resh(n-1)
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n-2) =     h(n-2) + cb*kh  (n-2)
         h (n-1) =     h(n-1) + cb*kh  (n-1)
         h (n  ) =     h(n  ) + cb*kh  (n  )
      endif
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine rk_storage
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      if     (IFRK45) then

      rk4a(1) =   0.0
      rk4a(2) =  -567301805773.0/1357537059087.0
      rk4a(3) =  -2404267990393.0/2016746695238.0
      rk4a(4) =  -3550918686646.0/2091501179385.0
      rk4a(5) =  -1275806237668.0/842570457699.0

      rk4b(1) =   1432997174477.0/9575080441755.0
      rk4b(2) =   5161836677717.0/13612068292357.0
      rk4b(3) =   1720146321549.0/2090206949498.0
      rk4b(4) =   3134564353537.0/4481467310338.0
      rk4b(5) =   2277821191437.0/14882151754819.0

      rk4c(1) =   0.0
      rk4c(2) =   1432997174477.0/9575080441755.0
      rk4c(3) =   2526269341429.0/6820363962896.0
      rk4c(4) =   2006345519317.0/3224310063776.0
      rk4c(5) =   2802321613138.0/2924317926251.0
      rk4c(6) =   1.

      elseif (IFRK22) then

      rk4a(1) =   0.0 
      rk4a(2) =  -1.0                         

      rk4b(1) =   1.0 
      rk4b(2) =   1.0/2.0                       

c     rk4a(1) =   0.0 
c     rk4a(2) =  -11/15 !756391.0/934407.0
c     rk4a(3) =  -5/3   !36441873.0/15625003.0
c     rk4a(4) =  -1     !1953125.0/1085297.0
c     !rk4a(5) =   0.0
c
c     rk4b(1) =   1/3  !8.0/141.0 
c     rk4b(2) =   5/6  !6627.0/2000.0
c     rk4b(3) =   3/5  !609375.0/1085297.0
c     rk4b(4) =   1/4  !198961.0/526383.0
c     !rk4b(5) =   0.0
c
c     rk4c(1) =   0.0
c     rk4c(2) =   1.0/3.0
c     rk4c(3) =   5.0/9.0
c     rk4c(4) =   8.0/9.0
c     !rk4c(5) =   1.0
c
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_inc                            
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      !call userinc                                           

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_source                         
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      if (.not.ifsrc) return

      call usersrc (57,reshn(1,1),reshn(1,2),reshn(1,3),
     $                 resen(1,1),resen(1,2),resen(1,3))

      return
      end
c-----------------------------------------------------------------------
      subroutine rectangular_to_cylindrical(vr,vt,vz1,vx,vy,vz,theta)   
      implicit none
      real       vr,vt,vz1,vx,vy,vz,theta   

      vr= vx*cos(theta)+vy*sin(theta)
      vt=-vx*sin(theta)+vy*cos(theta)
      vz1= vz

      return
      end

c-----------------------------------------------------------------------
      subroutine cylindrical_to_rectangular(vx,vy,vz1,vr,vt,vz,theta)   
      implicit none
      real       vr,vt,vz1,vx,vy,vz,theta   

      vx= vr*cos(theta)-vt*sin(theta)
      vy= vr*sin(theta)+vt*cos(theta)
      vz1= vz

      return
      end
c-----------------------------------------------------------------------
      subroutine energy_field (u,v1,v2,v3,n)                 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer  i,n
      real     u(n),v1(n),v2(n),v3(n)
      real     tmp,const

      const = 1.0

      do i=1,n 
         tmp = v1(i)**2+v2(i)**2+v3(i)**2
         tmp = const*tmp
         u(i)= sqrt (tmp)
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine outfldrn (x,txt,kk)
      implicit none
      include 'SIZE'
      real        x(lx1,ly1,lz1,lelt)
      character*7 txt
      integer     kk,ie,iz,i,j,i1,i0 
C
      do ie=1,nelv,4
         do iz=1,nz1,1
            if (iz.eq.1) write(6,106) txt,iz,ie,kk
            if (iz.gt.1) write(6,107) 
            i1 = ie+1
            do j=ny1,1,-1
               write(6,105) ((x(i,j,iz,ie+i0),i=1,nx1),i0=0,3)
            enddo
         enddo
      enddo
C
  107 FORMAT(' ')
  105 FORMAT(4(4f6.0,5x))
c 105 FORMAT(4f6.0,20x,4f6.0)
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',2I4,/,
     $           5X,'       X            ')
C
      return
      end

c-----------------------------------------------------------------------
      subroutine outmat (u,m,n,name6,ie)
      implicit none
      include 'SIZE'
      integer     i,j,m,n,ie,n12     
      real        u(m,n)
      character*6 name6
c
c     Print out copies of a global matrix
c
      write(6,*)
      write(6,*) ie,' matrix: ',name6,m,n
      n12 = min(n,12)
      do i=1,m
         write(6,6) ie,name6,(u(i,j),j=1,n12)
      enddo
    6 format(i3,1x,a6,20e15.7)
      write(6,*)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_fc_ptr
c
c     Set up pointer to restrict u to faces ! NOTE: compact
c
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
c
      integer k,l,i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2

#ifdef GPU        
#ifdef _OPENACC
      integer p
  
      COMMON /cemfce_acc/ cemface2(2*ldim*lx1*lz1*lelt)
     $     ,              cemface_ptr(0:2*ldim*lx1*lz1*lelt)
     $     ,              cemface_ind(2*ldim*lx1*lz1*lelt)
      integer cemface2, cemface_ptr, cemface_ind
#endif
#endif

      if (nid.eq.0) then
         write(6,*) 'start: cem_set_fc_ptr, volume index j=cemface(i)'
      endif

      k = 0
      l = 0
      do e=1,nelt
      do f=1,nfaces

         ef     = eface(f)
         js1    = skpdat(1,f)
         jf1    = skpdat(2,f)
         jskip1 = skpdat(3,f)
         js2    = skpdat(4,f)
         jf2    = skpdat(5,f)
         jskip2 = skpdat(6,f)

         i = 0
         do j2=js2,jf2,jskip2
         do j1=js1,jf1,jskip1

            i = i+1
            l = l+1
            k = i+nxzf*(ef-1)+nxzf*nfaces*(e-1)           ! face   numbering
            cemface(k) = j1+nx1*(j2-1)+nxyz*(e-1) ! global numbering
            ncemface   = l
            !write(6,*) 'fc',ncemface,e,k,f,ef!,srfout(k)!glo_num3d(cemface(k))
            ! glo_num(cemface(k))  get global numbering
            ! cemface gets the Ed's numbering at faces following pff order

         enddo
         enddo

      enddo
      enddo

#ifdef GPU       
#ifdef _OPENACC
      call icopy(cemface_ind,cemface,ncemface)
      call isort(cemface_ind,cemface2,ncemface)

C Used in GPU calculation
      p = 1
      cemface_ptr(1) = cemface_ind(1)
      do i=2,ncemface
         if (cemface_ind(i) .ne. cemface_ind(i-1)) then
            p = p+1
            cemface_ptr(p) = i
         endif
      enddo
      cemface_ptr(p+1) = ncemface+1
      cemface_ptr(0) = p
#endif
#endif

      if (nid.eq.0) then
         write(6,*) 'done: cem_set_fc_ptr, ncemface= ',ncemface
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine cem_alternate_flux_ptr
c
c     Set up pointer to restrict u to faces ! NOTE: compact
c
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
c
      integer k,l,i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2
      integer k_opp, ef_opp

      if (nid.eq.0) then
          write(6,*) 'start: cem_alternatve_flux_ptr-> index setting'
      endif

      k = 0
      l = 0
      do e = 1,nelt
      do f = 1,nfaces,2
         ef    = eface(f  )
         ef_opp= eface(f+1) ! opposite
         do i=1,nxzf    
            k     = i+nxzf*(ef-1)    +nxzf*nfaces*(e-1)   ! face   numbering
            k_opp = i+nxzf*(ef_opp-1)+nxzf*nfaces*(e-1)   ! face   numbering
            cemface_alt(k_opp)= cemface(k    )
            cemface_alt(k    )= cemface(k_opp)
c           write(6,*) 'k/k_opp',k_opp,cemface_alt(k_opp) 
c    $                          ,k,cemface_alt(k)
         enddo

      enddo
      enddo

c     do i=1,ncemface
c        write(6,*) 'cemface/cemface_alt',i,cemface(i),cemface_alt(i)
c     enddo

      if (nid.eq.0) then
         write(6,*) 'done: cem_alternate_flux_ptr'
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine measure_time(name8,ireset,t0,time_sum,iprint)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'RTIMER'   
      real*8      dclock
      real        t0, dtime, time_init, time_sum
      character*8 name8
      integer     ireset,iprint,ifreq
      integer icalld
      save    icalld
      data    icalld/0/
      
      dtime   = dclock()-t0
      time_sum= time_sum+dtime
      t0      = dtime

      ifreq   = 5
      icalld  = icalld+1

#ifdef TIMER
      if (iprint.ne.0) then
      if (istep.le.10.or.mod(istep,iocomm).eq.0)   then
      if ((nid.eq.0).and.(mod(icalld,ifreq).eq.0)) then
          if (ireset.eq.0) then
              write(6,1000) istep,name8,t0    
          else
              write(6,1001) istep,name8,t0,time_sum        
          endif
      endif
      endif
      endif
#endif

 1000 format(i8,1x,a8,1x,2e15.7)
 1001 format(i8,1x,a8,1x,3e15.7)

      return
      end
c-----------------------------------------------------------------------
      subroutine measure_comm(t0)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'RTIMER'  
      integer icalld
      save    icalld
      data    icalld/0/
      real*8  dclock,t0,dtime,dtime5
      integer ifreq

      if (icalld.eq.0) comm_t=0.0

      dtime = dclock()-t0
      comm_t= comm_t+dtime
      dtime5= comm_t/istep
      icalld= icalld+1

      ifreq=5             
      if (ifexp) ifreq=larnol

#ifdef TIMER
      if (istep.le.10.or.mod(istep,iocomm).eq.0)   then
      if ((nid.eq.0).and.(mod(icalld,ifreq).eq.0)) then
       write(6,*) 'cpu: gs_op_fields per call   (sec)::',dtime  
       write(6,*) 'cpu: gs_op_fields per step   (sec)::',dtime5  
       write(6,*) 'cpu: gs_op_fields total      (sec)::',comm_t  
      endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine measure_comm_acc(t0)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'RTIMER'  
      integer icalld
      save    icalld
      data    icalld/0/
      real*8  dclock,t0,dtime,dtime5
      integer ifreq

      if (icalld.eq.0) comm_t_acc=0.0
      if (icalld.eq.0) comm_t_gpu=0.0

      dtime     = dclock()-t0
#ifdef GPU      
#ifdef _OPENACC
      comm_t_acc= comm_t_acc+dtime
      dtime5    = comm_t_acc/istep
#else
      comm_t_gpu= comm_t_gpu+dtime
      dtime5    = comm_t_gpu/istep
#endif
#endif
      icalld    = icalld+1

      ifreq=5             
      if (ifexp) ifreq=larnol

#ifdef TIMER
      if (istep.le.10.or.mod(istep,iocomm).eq.0)   then
      if ((nid.eq.0).and.(mod(icalld,ifreq).eq.0)) then
       write(6,*) 'acc: gs_op_fields_acc per call (sec)::',dtime  
       write(6,*) 'acc: gs_op_fields_acc per step (sec)::',dtime5  
       write(6,*) 'acc: gs_op_fields_acc total    (sec)::',comm_t_acc  
      endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)
c-----------------------------------------------------------------------
C
C     -  Spectral interpolation from A to B via tensor products
C     -  scratch arrays: w(na*na*nb + nb*nb*na)
C
C     5/3/00  -- this routine replaces specmp in navier1.f, which
c                has a potential memory problem
C
C
      logical if3d
c
      real b(nb,nb,nb),a(na,na,na)
      real w(ldw)
      real ba(nb,na), ab(na,nb)
c
      ltest = na*nb
      if (if3d) ltest = na*na*nb + nb*na*na
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine specmpnface(b,nb,a,na,ba,ab,if3d,w,ldw)
            
c     -  Spectral interpolation from A to B via tensor products
c     -  scratch arrays: w(na*na*nb + nb*nb*na)
c
c     5/3/00  -- this routine replaces specmp in navier1.f, which        
c                has a potential memory problem                           
          
      logical if3d                                                     
c          
      real b(nb,nb),a(na,na)
      real w(ldw)
c          
      ltest = na*nb
           
      if (if3d)      ltest = na*na*nb + nb*na*na
          
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c         
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine intp_rstd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      integer mx, md, idir
      logical if3d
c
      parameter       (ldg=lxd**3,lwkd=4*lx1)
c      parameter       (ldg=1,lwkd=1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real            jgl,jgt,wkd
c
      parameter       (ld=2*lxd)
c      parameter       (ld=1    )
      common /ctmp0/  w(ld**ldim,2)

      integer   i
      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
 
      ldw = 2*(ld**ldim)
  
      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
         call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else
         call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif

      return
      end
             
c-----------------------------------------------------------------------
      subroutine intp_rsd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      integer mx, md, idir
      logical if3d
c
      parameter       (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/   d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                 jgt(ldg), wkd(lwkd)
      real             jgl,jgt,wkd
 
      parameter       (ld=2*lxd)
      common /ctmp0/   w(ld**ldim,2)

      integer   i
      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
 
      ldw = 2*(ld**ldim)
  
      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
          call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else 
          call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine gen_int(jgl,jgt,mp,np,w)
 
c     Generate interpolation from np GLL points to mp GL points
c
c        jgl  = interpolation matrix, mapping from velocity nodes to pressure
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c

      real jgl(mp,np),jgt(np*mp),w(1)
c
      iz = 1
      id = iz + np
c
      call zwgll (w(iz),jgt,np)
      call zwgll (w(id),jgt,mp)
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo
c
      call transpose(jgt,np,jgl,mp)
c
      return
      end
          
c-----------------------------------------------------------------------
      subroutine gen_dgl(dgl,dgt,mp,np,w)
c
c     Generate derivative from np GL points onto mp GL points
c
c        dgl  = interpolation matrix, mapping from velocity nodes to pressure
c        dgt  = transpose of interpolation matrix
c        w    = work array of size (3*np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
c
      real dgl(mp,np),dgt(np*mp),w(1)
c
c
      iz = 1
      id = iz + np  
c
      call zwgll  (w(iz),dgt,np)  ! GL points
      call zwgll  (w(id),dgt,mp)  ! GL points
c
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c
      n  = np -1 
      do i=1,mp   
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np  
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
         enddo
      enddo
c
      call transpose(dgt,np,dgl,mp)

c
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
      include 'SIZE'            ! need nid
      character*5  avar5,lvar5
      character*10 sub_name10
c
      if (n.gt.m) then
         write(6,*)             
         write(6,*) nid,n,m
c        write(6,*) nid,n,m,avar5,lvar5,sub_name10
c        write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i9,2(1x,a5),1x,a10)
         call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md)
c
c     Get pointer to jgl() for interpolation pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip, mx, md
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c
      if (ip.eq.0) then
c
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c
         nwrkd = mx + md

         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
c
         call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgl_ptr (ip,mx,md)
c
c     Get pointer to GL-GL interpolation dgl() for pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip, mx, md
c
      parameter (ld=2*lxd)
c      parameter (ld=1    )
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pdg (ij)
c
      if (ip.eq.0) then
c
         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c
         call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)

      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ifoption(ifif)        
      implicit none
      include 'SIZE'
      include 'EMWAVE'
      include 'PARALLEL'

      integer ifip,iglmax,imaxg
      logical ifif
 
      ifip =0
      if (ifif) ifip =1

      imaxg = iglmax(ifip,1)
      if (imaxg.eq.1) ifif =.true.
   
      return
      end
c-----------------------------------------------------------------------
      subroutine fcsum2(xsum,asum,x,e,f)
c
c     Compute the weighted sum of X over face f of element e
c
c     x is an (NX,NY,NZ) data structure
c     f  is in the preprocessor notation 
c
c     xsum is sum (X*area)
c     asum is sum (area)

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'
      real x(lx1,ly1,lz1,1)
      integer e,f,fd

      asum = 0.
      xsum = 0.

c     Set up counters ;  fd is the dssum notation.
      call dsset(nx1,ny1,nz1)
      fd     = eface1(f)
      js1    = skpdat(1,fd)
      jf1    = skpdat(2,fd)
      jskip1 = skpdat(3,fd)
      js2    = skpdat(4,fd)
      jf2    = skpdat(5,fd)
      jskip2 = skpdat(6,fd)

      i = 0
      do j2=js2,jf2,jskip2
      do j1=js1,jf1,jskip1
         i = i+1
         xsum = xsum+area(i,1,f,e)*x(j1,j2,1,e)
         asum = asum+area(i,1,f,e)
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      real function geom_xyradius()
      implicit none   
      ! return the maximum radius in the x-y plane 
      ! of the nodes in the mesh.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n,i
      real     radius_squared,radhere_squared
      real     gop_work
      
      n   = 8*nelt
      radius_squared = 0

      do i=1,n
         radhere_squared = xc(i,1)*xc(i,1) + yc(i,1)*yc(i,1)
         if (radhere_squared.ge.radius_squared) then
             radius_squared = radhere_squared
         endif
      enddo

      call gop(radius_squared, gop_work, 'M  ', 1)
      geom_xyradius = sqrt(radius_squared)

      return
      end
c ---------------------------------------------------------
      real function final_geom_xyradius()
      implicit none
      ! return the maximum radius in the x-y plane 
      ! of the nodes in the mesh.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer n,i
      real    radius_squared,radhere_squared
      real    gop_work
      
      n = nx1*ny1*nz1*nelt
      radius_squared = 0

      do i=1,n
         radhere_squared = xm1(i,1,1,1)*xm1(i,1,1,1) 
     $                   + ym1(i,1,1,1)*ym1(i,1,1,1)
         if (radhere_squared.ge.radius_squared) then
             radius_squared = radhere_squared
         endif
      enddo

      call gop(radius_squared, gop_work, 'M  ', 1)
      final_geom_xyradius = sqrt(radius_squared)

      return
      end
c-----------------------------------------------------------------------
      subroutine gh_face_extend(x,zg,n,gh_type,e,v)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      real     x(1),zg(1),e(1),v(1)
      integer  gh_type
      integer  n

      if (ndim.eq.2) then
         call gh_face_extend_2d(x,zg,n,gh_type,e,v)
      else
         call gh_face_extend_3d(x,zg,n,gh_type,e,v)
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine gh_face_extend_2d(x,zg,n,gh_type,e,v)
c-----------------------------------------------------------------------
      implicit none
c
c     Extend 2D faces into interior via gordon hall
c
c     gh_type:  1 - vertex only
c               2 - vertex and faces
c
c
      integer n,gh_type,ntot,i,j,ii,jj
      real    x(n,n)
      real    zg(n)
      real    e(n,n)
      real    v(n,n)
      real    si,sj,hi,hj
c
c     Build vertex interpolant
c
      ntot=n*n
      call rzero(v,ntot)
      do jj=1,n,n-1
      do ii=1,n,n-1
         do j=1,n
         do i=1,n
            si     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            sj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            v(i,j) = v(i,j) + si*sj*x(ii,jj)
         enddo
         enddo
      enddo
      enddo
      if (gh_type.eq.1) then
         call copy(x,v,ntot)
         return
      endif


c     Extend 4 edges
      call rzero(e,ntot)
c
c     x-edges
c
      do jj=1,n,n-1
         do j=1,n
         do i=1,n
            hj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            e(i,j) = e(i,j) + hj*(x(i,jj)-v(i,jj))
         enddo
         enddo
      enddo
c
c     y-edges
c
      do ii=1,n,n-1
         do j=1,n
         do i=1,n
            hi     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            e(i,j) = e(i,j) + hi*(x(ii,j)-v(ii,j))
         enddo
         enddo
      enddo

      call add3(x,e,v,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine gh_face_extend_3d(x,zg,n,gh_type,e,v)
c-----------------------------------------------------------------------
c
c     Extend faces into interior via gordon hall
c
c     gh_type:  1 - vertex only
c               2 - vertex and edges
c               3 - vertex, edges, and faces
c
c
      real x(n,n,n)
      real zg(n)
      real e(n,n,n)
      real v(n,n,n)
      integer gh_type
c
c     Build vertex interpolant
c
      ntot=n*n*n
      call rzero(v,ntot)
      do kk=1,n,n-1
      do jj=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            si       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            sj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            sk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            v(i,j,k) = v(i,j,k) + si*sj*sk*x(ii,jj,kk)
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
      if (gh_type.eq.1) then
         call copy(x,v,ntot)
         return
      endif
c
c
c     Extend 12 edges
      call rzero(e,ntot)
c
c     x-edges
c
      do kk=1,n,n-1
      do jj=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            e(i,j,k) = e(i,j,k) + hj*hk*(x(i,jj,kk)-v(i,jj,kk))
         enddo
         enddo
         enddo
      enddo
      enddo
c
c     y-edges
c
      do kk=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            e(i,j,k) = e(i,j,k) + hi*hk*(x(ii,j,kk)-v(ii,j,kk))
         enddo
         enddo
         enddo
      enddo
      enddo
c
c     z-edges
c
      do jj=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            e(i,j,k) = e(i,j,k) + hi*hj*(x(ii,jj,k)-v(ii,jj,k))
         enddo
         enddo
         enddo
      enddo
      enddo
c
      call add2(e,v,ntot)
c
      if (gh_type.eq.2) then
         call copy(x,e,ntot)
         return
      endif
c
c     Extend faces
c
      call rzero(v,ntot)
c
c     x-edges
c
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            v(i,j,k) = v(i,j,k) + hi*(x(ii,j,k)-e(ii,j,k))
         enddo
         enddo
         enddo
      enddo
c
c     y-edges
c
      do jj=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            v(i,j,k) = v(i,j,k) + hj*(x(i,jj,k)-e(i,jj,k))
         enddo
         enddo
         enddo
      enddo
c
c     z-edges
c
      do kk=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            v(i,j,k) = v(i,j,k) + hk*(x(i,j,kk)-e(i,j,kk))
         enddo
         enddo
         enddo
      enddo
c
      call add2(v,e,ntot)
      call copy(x,v,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine map_to_crs(a,na,b,nb,if3d,w,ldw)
c-----------------------------------------------------------------------
      implicit none
c
c     Input:   b
c     Output:  a
c
      real     a(1),b(1),w(1)
      logical  if3d
      integer  na,nb,ldw
c
      integer   lx
      parameter(lx=40)
      real      za(lx),zb(lx)
c
      real    iba(lx*lx),ibat(lx*lx)
      save    iba,ibat
c
      integer nao,nbo
      save    nao,nbo
      data    nao,nbo  / -9, -9/
c
      if (na.gt.lx.or.nb.gt.lx) then
         write(6,*)'ERROR: increase lx in map_to_crs to max:',na,nb
         call exitt
      endif
c
      if (na.ne.nao  .or.   nb.ne.nbo) then
         nao = na
         nbo = nb
         call zwgll(za,w,na)
         call zwgll(zb,w,nb)
         call igllm(iba,ibat,zb,za,nb,na,nb,na)
      endif
c
      call specmpn(a,na,b,nb,iba,ibat,if3d,w,ldw)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_blowup(val)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real     val

      if (val.gt.1e20) then
        if (nid.eq.0) then
        write(6,*) 'solutions BLOW UP > 1e20 => simulation terminated!'
        endif
        call exitt
      endif

      return
      end
c-----------------------------------------------------------------------

#ifdef GPU         
#ifdef _OPENACC

c-----------------------------------------------------------------------
      subroutine rk4_upd_acc(h,kh,resh,cb,ca,dt,n)
c-----------------------------------------------------------------------
      implicit none
      real     h(n),kh(n),resh(n)
      real     cb,ca,dt,ktmp
      integer  n,nresid,iand,i
 
!$ACC DATA PRESENT(kh,resh,h)
!$ACC PARALLEL LOOP GANG VECTOR
      do i=1,n          
         kh(i  ) = ca*kh(i  ) + dt*resh(i  )
         h (i  ) =     h(i  ) + cb*kh  (i  )
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA

      return
      end

#endif 
#endif 
