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

      real     h(n),kh(n),resh(n)
      real     cb,ca,dt,ktmp
      integer  n,nresid,iand,i

#ifdef _OPENACC
!$ACC DATA PRESENT(kh,resh,h)
!$ACC PARALLEL LOOP GANG WORKER VECTOR
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
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      if (ifrk45) then
         rk4a(1) = 0.0
         rk4a(2) = -567301805773.0/1357537059087.0
         rk4a(3) = -2404267990393.0/2016746695238.0
         rk4a(4) = -3550918686646.0/2091501179385.0
         rk4a(5) = -1275806237668.0/842570457699.0

         rk4b(1) = 1432997174477.0/9575080441755.0
         rk4b(2) = 5161836677717.0/13612068292357.0
         rk4b(3) = 1720146321549.0/2090206949498.0
         rk4b(4) = 3134564353537.0/4481467310338.0
         rk4b(5) = 2277821191437.0/14882151754819.0

         rk4c(1) = 0.0
         rk4c(2) = 1432997174477.0/9575080441755.0
         rk4c(3) = 2526269341429.0/6820363962896.0
         rk4c(4) = 2006345519317.0/3224310063776.0
         rk4c(5) = 2802321613138.0/2924317926251.0
         rk4c(6) = 1.0
      elseif (ifrk22) then
         rk4a(1) = 0.0
         rk4a(2) = -1.0

         rk4b(1) = 1.0
         rk4b(2) = 1.0/2.0
      endif

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
c-----------------------------------------------------------------------
c     Set up pointer to restrict fields to faces.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      common /cemfce_acc/
     $     cemface2(2*ldim*lx1*lz1*lelt),
     $     cemface_ptr(0:2*ldim*lx1*lz1*lelt),
     $     cemface_ind(2*ldim*lx1*lz1*lelt)
      integer cemface2,cemface_ptr,cemface_ind

      integer k,l,i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2,p

      if (nid.eq.0) then
         write(6,*) 'start: cem_set_fc_ptr, volume index j=cemface(i)'
      endif

      k = 0
      l = 0
      do e = 1,nelt
         do f = 1,nfaces
            ef = eface(f)
            js1 = skpdat(1,f)
            jf1 = skpdat(2,f)
            jskip1 = skpdat(3,f)
            js2 = skpdat(4,f)
            jf2 = skpdat(5,f)
            jskip2 = skpdat(6,f)

            i = 0
            do j2 = js2,jf2,jskip2
               do j1 = js1,jf1,jskip1
                  i = i+1
                  l = l+1
c     Face numbering.
                  k = i+nxzf*(ef-1)+nxzf*nfaces*(e-1)
c     Global number; `cemface` gets Ed's numbering at faces following
c     pff order.
                  cemface(k) = j1+nx1*(j2-1)+nxyz*(e-1)
               enddo
            enddo
            ncemface = l
         enddo
      enddo

#ifdef _OPENACC
      call icopy(cemface_ind,cemface,ncemface)
      call isort(cemface_ind,cemface2,ncemface)

      p = 1
      cemface_ptr(1) = cemface_ind(1)
      do i = 2,ncemface
         if (cemface_ind(i).ne.cemface_ind(i-1)) then
            p = p+1
            cemface_ptr(p) = i
         endif
      enddo
      cemface_ptr(p+1) = ncemface+1
      cemface_ptr(0) = p
#endif

      if (nid.eq.0) then
         write(6,*) 'done: cem_set_fc_ptr, ncemface= ',ncemface
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

#ifdef _OPENACC
      comm_t_acc= comm_t_acc+dtime
      dtime5    = comm_t_acc/istep
#else
      comm_t_gpu= comm_t_gpu+dtime
      dtime5    = comm_t_gpu/istep
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
         call exitt(1)
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
         call exitt(1)
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
         call exitt(1)
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
         call exitt(1)
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
        call exitt(1)
      endif

      return
      end
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
      subroutine cem_grad (w1,w2,w3,u1)
c-----------------------------------------------------------------------
c     non-weighted grad:   (w1,w2,w3) = grad(u1)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w1(lpts),w2(lpts),w3(lpts),u1(lpts)

      real u1rw,u1sw,u1tw
      real u2rw,u2sw,u2tw
      real u3rw,u3sw,u3tw

      real rxmn_k,sxmn_k,txmn_k
      real rymn_k,symn_k,tymn_k
      real rzmn_k,szmn_k,tzmn_k
      integer e,i,j,k,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real u1r,u1s,u1t

      nn = nx1-1

      if (if3d) then
         do e = 1,nelt
            j = nxyz*(e-1)

            call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

            do i = 1,nxyz
               k = i+j

               u1rw = u1r(i)
               u1sw = u1s(i)
               u1tw = u1t(i)

               rxmn_k = rxmn(k)/jacm(k)
               sxmn_k = sxmn(k)/jacm(k)
               txmn_k = txmn(k)/jacm(k)
               rymn_k = rymn(k)/jacm(k)
               symn_k = symn(k)/jacm(k)
               tymn_k = tymn(k)/jacm(k)
               rzmn_k = rzmn(k)/jacm(k)
               szmn_k = szmn(k)/jacm(k)
               tzmn_k = tzmn(k)/jacm(k)

               w1(k) = u1rw*rxmn_k+u1sw*sxmn_k+u1tw*txmn_k
               w2(k) = u1rw*rymn_k+u1sw*symn_k+u1tw*tymn_k
               w3(k) = u1rw*rzmn_k+u1sw*szmn_k+u1tw*tzmn_k
            enddo
         enddo
      else
         do e = 1,nelt
            j = nxyz*(e-1)

            call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

            do i = 1,nxyz
               k = i+j

               u1rw = u1r(i)
               u1sw = u1s(i)

               rxmn_k= rxmn(k)/jacm(k)
               sxmn_k= sxmn(k)/jacm(k)

               rymn_k= rymn(k)/jacm(k)
               symn_k= symn(k)/jacm(k)

               w1(k) = u1rw*rxmn_k+u1sw*sxmn_k
               w2(k) = u1rw*rymn_k+u1sw*symn_k
               w3(k) = 0
            enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_div(w0,w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
c     div: w = div(u)/jacm
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3

      nn = nx1-1
c
      if (if3d) then
         do e = 1,nelt
            j = nxyz*(e-1)

            call local_grad3(ur1,us1,ut1,u1,nn,e,dxm1,dxtm1)
            call local_grad3(ur2,us2,ut2,u2,nn,e,dxm1,dxtm1)
            call local_grad3(ur3,us3,ut3,u3,nn,e,dxm1,dxtm1)

            do i=1,nxyz
               k = i+j
               w1(k) = (ur1(i)*rxmn(k)
     $                + us1(i)*sxmn(k)
     $                + ut1(i)*txmn(k))/jacm(k)

               w2(k) = (ur2(i)*rymn(k)
     $                + us2(i)*symn(k)
     $                + ut2(i)*tymn(k))/jacm(k)

               w3(k) = (ur3(i)*rzmn(k)
     $                + us3(i)*szmn(k)
     $                + ut3(i)*tzmn(k))/jacm(k)

               w0(k) = w1(k)+w2(k)+w3(k)
            enddo
         enddo
      else
         do e = 1,nelt

            j = nxyz*(e-1)

            call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
            call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

            do i = 1,nxyz

               k = i+j
               w1(k) = (ur1(i)*rxmn(k)
     $                + us1(i)*sxmn(k))/jacm(k)
               w2(k) = (ur2(i)*rymn(k)
     $                + us2(i)*symn(k))/jacm(k)

               w0(k) = w1(k)+w2(k)
            enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_error(u,exact,error,n,l2,linf)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'MASS'

      integer n
      real u(n),exact(n),error(n)
      real l2,linf

      real glsc3,glamax

      call sub3(error,exact,u,n)

      l2 = glsc3(error,bm1,error,n)/volvm1
      if (l2.gt.0.0) l2 = sqrt(l2)

      linf = glamax(error,n)

      return
      end
c-----------------------------------------------------------------------
