c---------------------------------------------------------------------
c...  This file includes the main routines for solving Maxwell eqs.
c...  Globally defined variable: npts, npts3, nxyz, nxzf,nxzfl, nfaces
c---------------------------------------------------------------------

c...  preproessor compile options: CPU/OMP/CUDA/OpenACC

c-----------------------------------------------------------------------
c...  this routine masure time for single or multiple communication call
c...  and automatically choose the faster one as an option
      subroutine cem_communication_cost_check
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      real*8   dclock
      real     t0,t1,t2
      real     srflx_tmp1(lxzfl*6) ! single vector array
      real     srflx_tmp2(lxzfl,6) ! 6 dimensions
      integer  nxzfl6,maxcomm,iglmax

      nxzfl6= nxzfl*6
      call rone (srflx_tmp1,nxzfl6)
      call rone (srflx_tmp2,nxzfl6)

      t0 = dclock()
      call gs_op_fields(gsh_face,srflx_tmp1,nxzfl,6,1,1,0)
      t1 = dclock()-t0

      t0= dclock()
      call gs_op_fields(gsh_face,srflx_tmp2(1,1),nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,srflx_tmp2(1,2),nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,srflx_tmp2(1,3),nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,srflx_tmp2(1,4),nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,srflx_tmp2(1,5),nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,srflx_tmp2(1,6),nxzfl,1,1,1,0)
      t2 = dclock()-t0

      if (nid.eq.0) then
         write(6,*) 'nid=0: gs_op_fields: single vector',t1,'sec'
         write(6,*) 'nid=0: gs_op_fields: multi  vector',t2,'sec'
      endif

      icomm_option = 0 ! default
      if (t1.gt.t2) icomm_option = 1
      maxcomm = iglmax(icomm_option,1)
      call bcast(maxcomm,isize)
      icomm_option = maxcomm

      if (nid.eq.0) then
         if (icomm_option.eq.0)
     $        write(6,*) 'gs_op_fields: chose single vector'
         if (icomm_option.eq.1)
     $        write(6,*) 'gs_op_fields: chose multi vector '
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine cem_maxwell_init
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'
      integer  i,npts3,nxzfl3

      npts3 = 3*npts
      nxzfl3 = 3*nxzfl

c     These to be deleted all when IFHEX/IFTET are done

      if (nid.eq.0) write(6,*) 'start: initialize geom and fields'
c     Coordinates, geometry, field variables
      call rzero (xmn,npts)
      call rzero (ymn,npts)
      call rzero (zmn,npts)
      call rzero (bmn,npts)
      call rzero (bminv,npts)
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
      call rzero (shn,npts3)
      call rzero (sen,npts3)
      call rzero (hn,npts3)
      call rzero (en,npts3)
      call rzero (khn,npts3)
      call rzero (ken,npts3)
      call rzero (errhn,npts3)
      call rzero (erren,npts3)
      call rzero (reshn,npts3)
      call rzero (resen,npts3)

c     Geometry related
      call rzero (w3mn ,nxyz)
      call rzero (unxm ,nxzfl)
      call rzero (unym ,nxzfl)
      call rzero (unzm ,nxzfl)
      call rzero (aream,nxzfl)

c     Material properties
      call rzero (permittivity,npts)
      call rzero (permeability,npts)
      call rzero (impede,npts)
      call rzero (conduc,npts)
      call rzero (faceif,npts)

c     pml related arrays
      call rzero (pmlsigma,lx1*ly1*lz1*lelt*3)
      call rzero (pmldn,lx1*ly1*lz1*lelt*3)
      call rzero (pmlbn,lx1*ly1*lz1*lelt*3)

      call rzero (respmldn,lx1*ly1*lz1*lelt*3)
      call rzero (respmlbn,lx1*ly1*lz1*lelt*3)

      call rzero (respmlen,lx1*ly1*lz1*lelt*3)
      call rzero (respmlhn,lx1*ly1*lz1*lelt*3)
      call rzero (kpmlbn,lx1*ly1*lz1*lelt*3)
      call rzero (kpmldn,lx1*ly1*lz1*lelt*3)
      call izero (pmlptr,lelt)

      do i=1,npts
         xmn (i) = xm1  (i,1,1,1)
         ymn (i) = ym1  (i,1,1,1)
         zmn (i) = zm1  (i,1,1,1)
         bmn (i) = bm1  (i,1,1,1)
         rxmn(i) = rxm1 (i,1,1,1)
         rymn(i) = rym1 (i,1,1,1)
         rzmn(i) = rzm1 (i,1,1,1)
         sxmn(i) = sxm1 (i,1,1,1)
         symn(i) = sym1 (i,1,1,1)
         szmn(i) = szm1 (i,1,1,1)
         txmn(i) = txm1 (i,1,1,1)
         tymn(i) = tym1 (i,1,1,1)
         tzmn(i) = tzm1 (i,1,1,1)
         jacm(i) = jacm1(i,1,1,1)
         bminv(i)= 1.0/bmn(i)
      enddo

      do i=1,nxyz
         w3mn(i) = w3m1 (i,1,1)
      enddo

      do i=1,nxzfl
         unxm (i)= unx (i,1,1,1)
         unym (i)= uny (i,1,1,1)
         unzm (i)= unz (i,1,1,1)
         aream(i)= area(i,1,1,1)
      enddo
      if (nid.eq.0) write(6,*) 'done: initialization geom/fields'

      call cem_set_fc_ptr       ! global numbering index on face
      if (ifpml) then
         call pml_errchk
         call pml_fill_faceary(faceary,pmlthick)
         call pml_extent_and_tags(pmlinner,pmlouter,pmltag,faceary)
      endif
      call cem_maxwell_uvp      ! returns material properties
      call cem_maxwell_materials ! other material related settings
      if (ifpml) then
         call pml_calc_sigma(pmlinner,pmlouter,pmltag,pmlorder,
     $                       pmlreferr)
      endif

c     Communication routines
      call cem_communication_cost_check

c     Inverse mass matrix including material constants
      call col3    (ebm1,permittivity,bmn,npts)
      call col3    (hbm1,permeability,bmn,npts)
      call invcol1 (ebm1,npts)
      call invcol1 (hbm1,npts)
      call invers2 (jacmi,jacm1,npts)
      if (nid.eq.0) write(6,*) 'done: inverse mass matrix'

      return
      end
c---------------------------------------------------------------------
      subroutine cem_maxwell_init_fields
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      if (ifrestart) then
         call restart_swap
         call userchk
         call cem_out
      else
         call userini(time,HN(1,1),HN(1,2),HN(1,3),
     $                     EN(1,1),EN(1,2),EN(1,3))
         call userchk
         call cem_out
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine cem_maxwell_uvp
c---------------------------------------------------------------------
c     Get material parameters from `uservp` and check whether they are
c     sensible.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer ie,npts3
      real glmax,glmin
      real permit_max,permit_min
      real permea_max,permea_min

      npts3 = 3*npts
c     All the arguments are just dummy arguments. This is a temporary
c     measure for Maxwell until we move to a system where `uservp` is
c     called only once and looping is handled in the .usr file. See
c
c     https://github.com/NekCEM/NekCEM/issues/12
c
      call uservp(1,1,1,1)

c     Check permittivity/permeability setting
      permit_max = glmax(permittivity,npts)
      permit_min = glmin(permittivity,npts)
      permea_max = glmax(permeability,npts)
      permea_min = glmin(permeability,npts)

      if (nid.eq.0) write(6,*) 'permit_max=',permit_max
      if (nid.eq.0) write(6,*) 'permit_min=',permit_min
      if (nid.eq.0) write(6,*) 'permea_max=',permea_max
      if (nid.eq.0) write(6,*) 'permea_min=',permea_min

      if ((permit_max.le.0).or.(permit_min.le.0)) then
         if (nid.eq.0) write(6,*) 'define permit in uservp (.usr)'
         call exitt(1)
      endif
      if ((permea_max.le.0).or.(permea_min.le.0)) then
         if (nid.eq.0) write(6,*) 'define permea in uservp (.usr)'
         call exitt(1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_materials
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer e,f,i,j
      character CB*3
      real const

      if (nid.eq.0) write(6,*) 'start: cem_maxwell_materials'

c     Compute impedence/conductance.
      do i = 1,npts
         impede(i) = sqrt(permeability(i)/permittivity(i))
         conduc(i) = sqrt(permittivity(i)/permeability(i))
      enddo

c     Restrict impedence/conductance to face.
      do j = 1,ncemface
         i = cemface(j)
         zimpede(j) = impede(i)
         yconduc(j) = conduc(i)
      enddo

      call copy(Z_0,zimpede,nxzfl)
      call copy(Y_0,yconduc,nxzfl)

c     Addition of the face values: Z_0 = Z_0^+ + Z_0^- and
c     Y_0 = Y_0^+ + Y_0^-. Note that gsh_face follows Ed's face
c     numbering.
      call gs_op_fields(gsh_face,Z_0,nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,Y_0,nxzfl,1,1,1,0)

c     Special case for PEC boundary: Z_0 = 2*Z_0^-, Y_0 = 2*Y_0^-
      do e = 1,nelt
         do f = 1,nfaces
            do j = 1,nxzf
               i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j ! face numbering
               CB = CBC(f,e,1)
               if (CB.eq.'PEC') then
                  Y_0(i) = 2.0*Y_0(i)
                  Y_1(i) = 2.0*Y_1(i)
                  Z_0(i) = 2.0*Z_0(i)
                  Z_1(i) = 2.0*Z_1(i)
               endif
            enddo
         enddo
      enddo

c     Z_1 = Z_0 - zimpede = Z_0^+
      call sub3 (Z_1,Z_0,zimpede,nxzfl)
      call sub3 (Y_1,Y_0,yconduc,nxzfl)

c     Z_0 = 0.5*Z_0
      const = 0.5
      call cmult (Z_0,const,nxzfl)
      call cmult (Y_0,const,nxzfl)

      if (nid.eq.0) write(6,*) 'done: cem_maxwell_materials'

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_rk
c-----------------------------------------------------------------------
c     Time step the PDE using a 5-stage 4th-order RK method.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      do rkstep = 1,5
         call rk_c(rkstep)
         call cem_maxwell_op
         call rk_maxwell_ab(rkstep)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_dummy
c-----------------------------------------------------------------------
c...  curl evaluation -> flux
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

      common /nonblock/ gs_aggress_nonblock
      integer           gs_aggress_nonblock

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx
      integer icalld
      save    icalld
      data    icalld/0/


c.... this pairs with nonblocking_4

      if (param(88).eq.1) then

         gs_aggress_nonblock= 1
         if (nid.eq.0) write(6,*) 'gs_aggress_nonblock:ON, param(88)=1'

         if (icalld.eq.0) then
            call cem_maxwell_restrict_to_face
            call gs_dummy_flux3d(srflx)
            call gs_op_fields_irecv(gsh_face,srflx,nxzfl,6,1,1,0)
            icalld=1
         else
            call gs_op_fields_isend(gsh_face,srflx,nxzfl,6,1,1,0)
            call gs_op_fields_wait (gsh_face,srflx,nxzfl,6,1,1,0)
         endif

      else

         gs_aggress_nonblock= 0
         if (nid.eq.0) write(6,*) 'gs_aggress_nonblock:OFF, param(88)=0'

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_timing1
c-----------------------------------------------------------------------
c...  curl evaluation -> flux
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7

      call cem_maxwell
      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux            (srflx)
      call cem_maxwell_add_flux_to_res (srflx)
      call cem_maxwell_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_timing2
c-----------------------------------------------------------------------
c...  flux -> curl evaluation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7

      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux            (srflx)
      call cem_maxwell
      call cem_maxwell_add_flux_to_res (srflx)
      call cem_maxwell_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_timing3
c-----------------------------------------------------------------------
c...  flux -> curl evaluation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7

      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux_nonblocking(srflx)
      call cem_maxwell_add_flux_to_res (srflx)
      call cem_maxwell_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_timing4
c-----------------------------------------------------------------------
c...  flux -> curl evaluation
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7

      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux_nonblocking_4(srflx)
      call cem_maxwell_add_flux_to_res (srflx)
      call cem_maxwell_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real srflx

      call cem_maxwell
      call cem_maxwell_restrict_to_face
      call userinc(rktime,fhn(1,1),fhn(1,2),fhn(1,3),
     $                    fen(1,1),fen(1,2),fen(1,3))
      call cem_maxwell_flux(srflx)
      call cem_maxwell_add_flux_to_res(srflx)
      if (ifpml) call pml_step
      call usersrc(rktime,reshn(1,1),reshn(1,2),reshn(1,3),
     $                    resen(1,1),resen(1,2),resen(1,3))
      call cem_maxwell_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer   i

c...  sacrify 2D optimization for simplicity::
      if (ifdealias) then

         if     (imode.eq.3) then !IF3D

            call maxwell_wght_dcurl
     $        (resEN(1,1),resEN(1,2),resEN(1,3),HN(1,1),HN(1,2),HN(1,3))
            call maxwell_wght_dcurl
     $        (resHN(1,1),resHN(1,2),resHN(1,3),EN(1,1),EN(1,2),EN(1,3))
            call chsign (resHN(1,1),npts)
            call chsign (resHN(1,2),npts)
            call chsign (resHN(1,3),npts)

         elseif (imode.eq.2) then !IFTM

            call maxwell_wght_dcurl
     $        (resHN(1,1),resHN(1,2),resEN(1,3),HN(1,1),HN(1,2),EN(1,3))
            call chsign (resHN(1,1),npts)
            call chsign (resHN(1,2),npts)

         elseif (imode.eq.1) then !IFTE
            call maxwell_wght_dcurl
     $        (resEN(1,1),resEN(1,2),resHN(1,3),EN(1,1),EN(1,2),HN(1,3))
            call chsign (resHN(1,3),npts)

         endif

      else

#ifdef _OPENACC

         if     (imode.eq.3) then !IF3D

            call maxwell_wght_curl_acc
     $        (resEN(1,1),resEN(1,2),resEN(1,3),HN(1,1),HN(1,2),HN(1,3))
            call maxwell_wght_curl_acc
     $        (resHN(1,1),resHN(1,2),resHN(1,3),EN(1,1),EN(1,2),EN(1,3))
            call chsign (resHN(1,1),npts)
            call chsign (resHN(1,2),npts)
            call chsign (resHN(1,3),npts)

         elseif (imode.eq.2) then !IFTM

            call maxwell_wght_curl_acc
     $        (resHN(1,1),resHN(1,2),resEN(1,3),HN(1,1),HN(1,2),EN(1,3))
            call chsign (resHN(1,1),npts)
            call chsign (resHN(1,2),npts)

         elseif (imode.eq.1) then !IFTE
            call maxwell_wght_curl_acc
     $        (resEN(1,1),resEN(1,2),resHN(1,3),EN(1,1),EN(1,2),HN(1,3))
            call chsign (resHN(1,3),npts)

         endif
#else

         if     (imode.eq.3) then !IF3D

            call maxwell_wght_curl
     $        (resEN(1,1),resEN(1,2),resEN(1,3),HN(1,1),HN(1,2),HN(1,3))
            call maxwell_wght_curl
     $        (resHN(1,1),resHN(1,2),resHN(1,3),EN(1,1),EN(1,2),EN(1,3))
            call chsign (resHN(1,1),npts)
            call chsign (resHN(1,2),npts)
            call chsign (resHN(1,3),npts)

         elseif (imode.eq.2) then !IFTM

            call maxwell_wght_curl
     $        (resHN(1,1),resHN(1,2),resEN(1,3),HN(1,1),HN(1,2),EN(1,3))
            call chsign (resHN(1,1),npts)
            call chsign (resHN(1,2),npts)

         elseif (imode.eq.1) then !IFTE
            call maxwell_wght_curl
     $        (resEN(1,1),resEN(1,2),resHN(1,3),EN(1,1),EN(1,2),HN(1,3))
            call chsign (resHN(1,3),npts)

         endif

#endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_restrict_to_face
c-----------------------------------------------------------------------
c     restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

!$ACC DATA PRESENT (fHN,fEN,HN,EN,cemface)
      if (imode.eq.3) then
!$ACC PARALLEL LOOP INDEPENDENT
         do j = 1,ncemface
            i = cemface(j)
            fHN(j,1) = HN(i,1)
            fHN(j,2) = HN(i,2)
            fHN(j,3) = HN(i,3)
            fEN(j,1) = EN(i,1)
            fEN(j,2) = EN(i,2)
            fEN(j,3) = EN(i,3)
         enddo
!$ACC END PARALLEL LOOP
      elseif (imode.eq.2) then
!$ACC PARALLEL LOOP INDEPENDENT
         do j = 1, ncemface
            i = cemface(j)
            fHN(j,1) = HN(i,1)
            fHN(j,2) = HN(i,2)
            fHN(j,3) = 0.0
            fEN(j,1) = 0.0
            fEN(j,2) = 0.0
            fEN(j,3) = EN(i,3)
         enddo
!$ACC END PARALLEL LOOP
      elseif (imode.eq.1) then
!$ACC PARALLEL LOOP INDEPENDENT
         do j = 1, ncemface
            i = cemface(j)
            fHN(j,1) = 0.0
            fHN(j,2) = 0.0
            fHN(j,3) = HN(i,3)
            fEN(j,1) = EN(i,1)
            fEN(j,2) = EN(i,2)
            fEN(j,3) = 0.0
         enddo
!$ACC END PARALLEL LOOP
      endif
!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_add_flux_to_res(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      common /cemfce_acc/
     $     cemface2(2*ldim*lx1*lz1*lelt),
     $     cemface_ptr(0:2*ldim*lx1*lz1*lelt),
     $     cemface_ind(2*ldim*lx1*lz1*lelt)
      integer cemface2,cemface_ptr,cemface_ind

      real srflx(6*2*ldim*lx1*lz1*lelt)

      integer i,j,k
#ifdef _OPENACC
      integer l,p,nptr
#endif
      real a

      k = nxzfl
#ifdef _OPENACC
      nptr = cemface_ptr(0)
!$ACC DATA PRESENT(resHN,resEN,aream,srflx,cemface2,cemface_ptr,
!$ACC&             cemface_ind)
      if (imode.eq.3) then      ! 3D
!$ACC PARALLEL LOOP
         do p = 1,nptr
!$ACC LOOP SEQ
            do l = cemface_ptr(p),cemface_ptr(p+1)-1
               i = cemface_ind(l)
               j = cemface2(l)
               a = aream(j)
               resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
               resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
               resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
               resEN(i,1) = resEN(i,1) + a*srflx(3*k+j)
               resEN(i,2) = resEN(i,2) + a*srflx(4*k+j)
               resEN(i,3) = resEN(i,3) + a*srflx(5*k+j)
            enddo
         enddo
      elseif (imode.eq.2) then  ! TM
!$ACC PARALLEL LOOP
         do p = 1,nptr
!$ACC LOOP SEQ
            do l = cemface_ptr(p),cemface_ptr(p+1)-1
               i = cemface_ind(l)
               j = cemface2(l)
               a = aream(j)
               resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
               resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
               resEN(i,3) = resEN(i,3) + a*srflx(2*k+j)
            enddo
         enddo
      elseif (imode.eq.1) then  ! TE
!$ACC PARALLEL LOOP
         do p = 1,nptr
!$ACC LOOP SEQ
            do l = cemface_ptr(p),cemface_ptr(p+1)-1
               i = cemface_ind(l)
               j = cemface2(l)
               a = aream(j)
               resEN(i,1) = resEN(i,1) + a*srflx(0*k+j)
               resEN(i,2) = resEN(i,2) + a*srflx(1*k+j)
               resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
            enddo
         enddo
      endif
!$ACC END DATA
#else
      if (imode.eq.3) then  ! 3D
         do j = 1,ncemface
            i = cemface(j)
            a = aream(j)
            resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
            resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
            resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
            resEN(i,1) = resEN(i,1) + a*srflx(3*k+j)
            resEN(i,2) = resEN(i,2) + a*srflx(4*k+j)
            resEN(i,3) = resEN(i,3) + a*srflx(5*k+j)
         enddo
      elseif (imode.eq.2) then  ! TM
         do j = 1,ncemface
            i = cemface(j)
            a = aream(j)
            resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
            resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
            resEN(i,3) = resEN(i,3) + a*srflx(2*k+j)
         enddo
      elseif (imode.eq.1) then  ! TE
         do j = 1,ncemface
            i = cemface(j)
            a = aream(j)
            resEN(i,1) = resEN(i,1) + a*srflx(0*k+j)
            resEN(i,2) = resEN(i,2) + a*srflx(1*k+j)
            resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
         enddo
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux(srflx)
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      real srflx(6*2*ldim*lx1*lz1*lelt)

      if (if3d) then
         call cem_maxwell_flux3d(srflx)
      else
         call cem_maxwell_flux2d(srflx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux_nonblocking(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      real     srflx(1)

      if (if3d) then
         call cem_maxwell_flux3d_nonblocking (srflx)
      else
         call cem_maxwell_flux2d_nonblocking (srflx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux_nonblocking_4(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      real     srflx(1)

      if (if3d) then
         call cem_maxwell_flux3d_nonblocking_4 (srflx)
      else
         call cem_maxwell_flux2d_nonblocking (srflx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux2d(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'RK5'

      real dclock,t0
      real srflx(6*2*ldim*lx1*lz1*lelt)
      real C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k,e,f

      nxzfl = nx1*nz1*nelt*2*ndim

      if (ifcentral) C0 = 0.0
      if (ifupwind) C0 = 1.0

!$ACC  DATA PRESENT(srflx,cemface,unxm,unym,fEN,fHN)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)

      if (imode.eq.2) then


c     TM case
         k = nxzfl
!$ACC PARALLEL LOOP
         do i = 1,nxzfl
c     -(n x E)_x
            srflx(0*k+i) = -unym(i)*fEN(i,3)
c     -(n x E)_y
            srflx(1*k+i) = unxm(i)*fEN(i,3)
c     -(n x H)_z
            srflx(2*k+i) = -unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
         enddo
!$ACC END PARALLEL LOOP

         call userfsrc(rktime,srflx(3*k+1),srflx(4*k+1),srflx(2*k+1),
     $                        srflx(0*k+1),srflx(1*k+1),srflx(5*k+1))

!$ACC UPDATE HOST(srflx)
         t0 = dclock()
         call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)
         call measure_comm(t0)

         if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)
!$ACC UPDATE DEVICE(srflx)

         k = nxzfl
!$ACC PARALLEL LOOP
         do i = 1,nxzfl
            Y0 = Y_0(i)
            Y1 = Y_1(i)
            Z0 = Z_0(i)
            Z1 = Z_1(i)

            fu1 = unym(i)*srflx(2*k+i)
            fu2 = -unxm(i)*srflx(2*k+i)
            fw3 = unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

            srflx(0*k+i) = 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)
            srflx(1*k+i) = 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)
            srflx(2*k+i) = 0.5/Z0*(Z1*srflx(2*k+i)-C0*fw3)
         enddo
!ACC END PARALLEL LOOP


      elseif (imode.eq.1) then


c     TE case
         k = nxzfl
!$ACC PARALLEL LOOP
         do i = 1,nxzfl
c     -(n x H)_x
            srflx(0*k+i) = -unym(i)*fHN(i,3)
c     -(n x H)_y
            srflx(1*k+i) = unxm(i)*fHN(i,3)
c     -(n x E)_z
            srflx(2*k+i) = -unxm(i)*fEN(i,2)+unym(i)*fEN(i,1)
         enddo
!$ACC END PARALLEL LOOP

         call userfsrc(rktime,srflx(0*k+1),srflx(1*k+1),srflx(3*k+1),
     $                        srflx(4*k+1),srflx(5*k+1),srflx(2*k+1))


!$ACC UPDATE HOST(srflx)
         t0 = dclock()
         call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)
         call measure_comm(t0)

         if (ifpml.or.ifpec) call cem_maxwell_flux_pec(srflx)
!$ACC UPDATE DEVICE(srflx)

         k = nxzfl
!$ACC PARALLEL LOOP
         do i = 1,nxzfl
            Y0 = Y_0(i)
            Y1 = Y_1(i)
            Z0 = Z_0(i)
            Z1 = Z_1(i)
c     (n x (n x [[E]]))_x
            fw1 =  unym(i)*srflx(2*k+i)
c     (n x (n x [[E]]))_y
            fw2 = -unxm(i)*srflx(2*k+i)
c     n x (n x [[H]])
            fu3 =  unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)
            srflx(0*k+i) = 0.5/Z0*(Z1*srflx(0*k+i)-C0*fw1)
            srflx(1*k+i) = 0.5/Z0*(Z1*srflx(1*k+i)-C0*fw2)
            srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i)-C0*fu3)
         enddo
!$ACC END PARALLEL LOOP
      endif

!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux3d(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'RK5'

      real dclock,t0
      real srflx(6*2*ldim*lx1*lz1*lelt)
      real C0,y0,y1,z0,z1,y02,z02,c02y,c02z,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k

      if (ifcentral) C0 = 0.0
      if (ifupwind) C0 = 1.0

!$ACC  DATA PRESENT(srflx,fEN,fHN,unxm,unym,unzm)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)
      k = nxzfl
!$ACC PARALLEL LOOP
      do i = 1,nxzfl
c     - n x E
         srflx(0*k+i) = -unym(i)*fEN(i,3)+unzm(i)*fEN(i,2)
         srflx(1*k+i) = -unzm(i)*fEN(i,1)+unxm(i)*fEN(i,3)
         srflx(2*k+i) = -unxm(i)*fEN(i,2)+unym(i)*fEN(i,1)
c     - n x H
         srflx(3*k+i) = -unym(i)*fHN(i,3)+unzm(i)*fHN(i,2)
         srflx(4*k+i) = -unzm(i)*fHN(i,1)+unxm(i)*fHN(i,3)
         srflx(5*k+i) = -unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
      enddo
!$ACC END PARALLEL LOOP

      call userfsrc(rktime,srflx(3*k+1),srflx(4*k+1),srflx(5*k+1),
     $                     srflx(0*k+1),srflx(1*k+1),srflx(2*k+1))

      t0 = dclock()
      call gs_op_fields(gsh_face,srflx,nxzfl,6,1,1,0)
      call measure_comm(t0)

c     The array srflx now contains n x [E] and n x [H]. (un^+ = - un^-.)
      if (ifpec.or.ifpml) call cem_maxwell_flux_pec(srflx)

      k = nxzfl
!$ACC PARALLEL LOOP
      do i = 1,nxzfl
         Y0 = Y_0(i)
         Y1 = Y_1(i)
         Z0 = Z_0(i)
         Z1 = Z_1(i)

         Y02 = -0.5/Y0*Y1
         Z02 = 0.5/Z0*Z1
         C02Y = 0.5/Y0*C0
         C02Z = 0.5/Z0*C0

c     fu = n x (n x [H])
         fu1 =unym(i)*srflx(5*k+i)-unzm(i)*srflx(4*k+i)
         fu2 =unzm(i)*srflx(3*k+i)-unxm(i)*srflx(5*k+i)
         fu3 =unxm(i)*srflx(4*k+i)-unym(i)*srflx(3*k+i)

c     fw = n x (n x [E])
         fw1 =unym(i)*srflx(2*k+i)-unzm(i)*srflx(1*k+i)
         fw2 =unzm(i)*srflx(0*k+i)-unxm(i)*srflx(2*k+i)
         fw3 =unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

         srflx(0*k+i) = Y02*srflx(0*k+i)-C02Y*fu1
         srflx(1*k+i) = Y02*srflx(1*k+i)-C02Y*fu2
         srflx(2*k+i) = Y02*srflx(2*k+i)-C02Y*fu3
         srflx(3*k+i) = Z02*srflx(3*k+i)-C02Z*fw1
         srflx(4*k+i) = Z02*srflx(4*k+i)-C02Z*fw2
         srflx(5*k+i) = Z02*srflx(5*k+i)-C02Z*fw3
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux2d_nonblocking (srflx)
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


      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0

!$ACC  DATA PRESENT(srflx,cemface,unxm,unym,fEN,fHN)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)
      if     (imode.eq.2) then !IFTM

         k = nxzfl
!$ACC PARALLEL LOOP
         do i=1,nxzfl
            j=cemface(i)
            srflx(0*k+i)=-unym(i)*fEN(i,3) ! NxfHx
            srflx(1*k+i)= unxm(i)*fEN(i,3) ! NxfHy
            srflx(2*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1) ! NxfEz
         enddo
!$ACC END PARALLEL LOOP

c Multi vector gs_op_fields option removed in svn r2015
         t0=dclock()
         call gs_op_fields_irecv(gsh_face,srflx,nxzfl,3,1,1,0)
         call gs_op_fields_isend(gsh_face,srflx,nxzfl,3,1,1,0)
         call cem_maxwell
         call gs_op_fields_wait (gsh_face,srflx,nxzfl,3,1,1,0)

         call measure_comm(t0)

         if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)

         k = nxzfl
!$ACC PARALLEL LOOP
         do i=1,nxzfl

            Y0 = Y_0(i)
            Y1 = Y_1(i)
            Z0 = Z_0(i)
            Z1 = Z_1(i)

            fu1 = unym(i)*srflx(2*k+i)
            fu2 =-unxm(i)*srflx(2*k+i)
            fw3 = unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

            srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)
            srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)
            srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)
         enddo
!ACC END PARALLEL LOOP

      elseif (imode.eq.1) then  !IFTE

         k = nxzfl
!$ACC PARALLEL LOOP
         do i=1,nxzfl
            srflx(0*k+i)= -unym(i)*fHN(i,3) ! NxfHx
            srflx(1*k+i)=  unxm(i)*fHN(i,3) ! NxfHy
            srflx(2*k+i)= -unxm(i)*fEN(i,2)+unym(i)*fEN(i,1) ! NxfEz
         enddo
!$ACC END PARALLEL LOOP

         t0=dclock()
         call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)

         call measure_comm(t0)

         if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)

         k = nxzfl
!$ACC PARALLEL LOOP
         do i=1,nxzfl

            Y0 = Y_0(i)
            Y1 = Y_1(i)
            Z0 = Z_0(i)
            Z1 = Z_1(i)

            fw1 = unym(i)*srflx(2*k+i)
            fw2 =-unxm(i)*srflx(2*k+i)
            fu3 = unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

            srflx(0*k+i) = 0.5/Z0*( Z1*srflx(0*k+i) - C0*fw1)
            srflx(1*k+i) = 0.5/Z0*( Z1*srflx(1*k+i) - C0*fw2)
            srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3)

         enddo
!$ACC END PARALLEL LOOP
      endif
!$ACC END DATA
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux3d_nonblocking(srflx)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

c      n x dEx - c* n x n x dHx
c
      real*8  dclock,t0
      real    srflx(1)
      real    C0,y0,y1,z0,z1,y02,z02,c02y,c02z,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k
      real    srflx_tmp(lxzfl1,6)

      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex

!$ACC  DATA PRESENT(srflx,fEN,fHN,unxm,unym,unzm)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)
      k = nxzfl
!$ACC PARALLEL LOOP
      do i= 1,nxzfl
         ! - n x E
         srflx(0*k+i)=-unym(i)*fEN(i,3)+unzm(i)*fEN(i,2)
         srflx(1*k+i)=-unzm(i)*fEN(i,1)+unxm(i)*fEN(i,3)
         srflx(2*k+i)=-unxm(i)*fEN(i,2)+unym(i)*fEN(i,1)
         ! - n x H
         srflx(3*k+i)=-unym(i)*fHN(i,3)+unzm(i)*fHN(i,2)
         srflx(4*k+i)=-unzm(i)*fHN(i,1)+unxm(i)*fHN(i,3)
         srflx(5*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
      enddo
!$ACC END PARALLEL LOOP

      call gs_op_fields_irecv(gsh_face,srflx,nxzfl,6,1,1,0)
      call gs_op_fields_isend(gsh_face,srflx,nxzfl,6,1,1,0)
      call cem_maxwell
      call gs_op_fields_wait (gsh_face,srflx,nxzfl,6,1,1,0)

      !srflx now contains n x [E] and n x [H],with [F] = F^+ - F^-;(observe un^+ = - un^-)
      if (IFPEC.or.IFPML)  call cem_maxwell_flux_pec(srflx)

      k = nxzfl
!$ACC PARALLEL LOOP
      do i=1,nxzfl
         Y0   = Y_0(i)
         Y1   = Y_1(i)
         Z0   = Z_0(i)
         Z1   = Z_1(i)
         Y02  =-0.5/Y0*Y1
         Z02  = 0.5/Z0*Z1
         C02Y = 0.5/Y0*C0
         C02Z = 0.5/Z0*C0

c     fu = n x n x [H]
         fu1 =unym(i)*srflx(5*k+i)-unzm(i)*srflx(4*k+i)
         fu2 =unzm(i)*srflx(3*k+i)-unxm(i)*srflx(5*k+i)
         fu3 =unxm(i)*srflx(4*k+i)-unym(i)*srflx(3*k+i)

c     fw = n x n x [E]
         fw1 =unym(i)*srflx(2*k+i)-unzm(i)*srflx(1*k+i)
         fw2 =unzm(i)*srflx(0*k+i)-unxm(i)*srflx(2*k+i)
         fw3 =unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

         srflx(0*k+i) = Y02*srflx(0*k+i) - C02Y*fu1
         srflx(1*k+i) = Y02*srflx(1*k+i) - C02Y*fu2
         srflx(2*k+i) = Y02*srflx(2*k+i) - C02Y*fu3
         srflx(3*k+i) = Z02*srflx(3*k+i) - C02Z*fw1
         srflx(4*k+i) = Z02*srflx(4*k+i) - C02Z*fw2
         srflx(5*k+i) = Z02*srflx(5*k+i) - C02Z*fw3

      enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_dummy_flux3d(srflx)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

      common /nonblock/ gs_aggress_nonblock
      integer           gs_aggress_nonblock

      integer icalld
      save    icalld
      data    icalld/0/

c      n x dEx - c* n x n x dHx
c
      real*8  dclock,t0
      real    srflx(1)
      real    C0,y0,y1,z0,z1,y02,z02,c02y,c02z,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k


      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex

!$ACC  DATA PRESENT(srflx,fEN,fHN,unxm,unym,unzm)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)
      k = nxzfl
!$ACC PARALLEL LOOP
      do i= 1,nxzfl
         ! - n x E
         srflx(0*k+i)=-unym(i)*fEN(i,3)+unzm(i)*fEN(i,2)
         srflx(1*k+i)=-unzm(i)*fEN(i,1)+unxm(i)*fEN(i,3)
         srflx(2*k+i)=-unxm(i)*fEN(i,2)+unym(i)*fEN(i,1)
         ! - n x H
         srflx(3*k+i)=-unym(i)*fHN(i,3)+unzm(i)*fHN(i,2)
         srflx(4*k+i)=-unzm(i)*fHN(i,1)+unxm(i)*fHN(i,3)
         srflx(5*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux3d_nonblocking_4(srflx)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

      common /nonblock/ gs_aggress_nonblock
      integer           gs_aggress_nonblock

      integer icalld
      save    icalld
      data    icalld/0/

c      n x dEx - c* n x n x dHx
c
      real*8  dclock,t0
      real    srflx(1)
      real    C0,y0,y1,z0,z1,y02,z02,c02y,c02z,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k


      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex

!$ACC  DATA PRESENT(srflx,fEN,fHN,unxm,unym,unzm)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)
      k = nxzfl
!$ACC PARALLEL LOOP
      do i= 1,nxzfl
         ! - n x E
         srflx(0*k+i)=-unym(i)*fEN(i,3)+unzm(i)*fEN(i,2)
         srflx(1*k+i)=-unzm(i)*fEN(i,1)+unxm(i)*fEN(i,3)
         srflx(2*k+i)=-unxm(i)*fEN(i,2)+unym(i)*fEN(i,1)
         ! - n x H
         srflx(3*k+i)=-unym(i)*fHN(i,3)+unzm(i)*fHN(i,2)
         srflx(4*k+i)=-unzm(i)*fHN(i,1)+unxm(i)*fHN(i,3)
         srflx(5*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
      enddo
!$ACC END PARALLEL LOOP

      call gs_op_fields_isend(gsh_face,srflx,nxzfl,6,1,1,0)
      call cem_maxwell
      call gs_op_fields_wait (gsh_face,srflx,nxzfl,6,1,1,0)
      call gs_op_fields_irecv(gsh_face,srflx,nxzfl,6,1,1,0)


      !srflx now contains n x [E] and n x [H],with [F] = F^+ - F^-;(observe un^+ = - un^-)
      if (IFPEC.or.IFPML)  call cem_maxwell_flux_pec(srflx)

      k = nxzfl
!$ACC PARALLEL LOOP
      do i=1,nxzfl
         Y0   = Y_0(i)
         Y1   = Y_1(i)
         Z0   = Z_0(i)
         Z1   = Z_1(i)
         Y02  =-0.5/Y0*Y1
         Z02  = 0.5/Z0*Z1
         C02Y = 0.5/Y0*C0
         C02Z = 0.5/Z0*C0

c     fu = n x n x [H]
         fu1 =unym(i)*srflx(5*k+i)-unzm(i)*srflx(4*k+i)
         fu2 =unzm(i)*srflx(3*k+i)-unxm(i)*srflx(5*k+i)
         fu3 =unxm(i)*srflx(4*k+i)-unym(i)*srflx(3*k+i)

c     fw = n x n x [E]
         fw1 =unym(i)*srflx(2*k+i)-unzm(i)*srflx(1*k+i)
         fw2 =unzm(i)*srflx(0*k+i)-unxm(i)*srflx(2*k+i)
         fw3 =unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

         srflx(0*k+i) = Y02*srflx(0*k+i) - C02Y*fu1
         srflx(1*k+i) = Y02*srflx(1*k+i) - C02Y*fu2
         srflx(2*k+i) = Y02*srflx(2*k+i) - C02Y*fu3
         srflx(3*k+i) = Z02*srflx(3*k+i) - C02Z*fw1
         srflx(4*k+i) = Z02*srflx(4*k+i) - C02Z*fw2
         srflx(5*k+i) = Z02*srflx(5*k+i) - C02Z*fw3

      enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA


      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux_pec(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      common /bdry1/ cempec(6*lx1*lz1*2*ldim*lelt),ncempec
      integer cempec,ncempec

      real srflx(1)
      integer k,i,j,l,e,f,ef,ii
      character CB*3

      integer icalld
      save icalld
      data icalld /0/

      k = nx1*nz1*2*ndim*nelt

      if (icalld.eq.0) then
         l = 0
         do e= 1,nelt
            do f= 1,nfaces
               do j= 1,nx1*nz1
                  CB= CBC(f,e,1)
                  i =(e-1)*2*ndim*nx1*nz1+(f-1)*nx1*nz1+j
                  if ((CB.eq.'PEC').or.(CB.eq.'PML')) then
                     l= l+1
                     cempec(l)= i
                     ncempec  = l
                  endif
               enddo
            enddo
         enddo
         icalld = 1
      endif


c     Here srflx is (n x -E, n x -H) since we're on the boundary and
c     fgs_gop left us untouched. We want mirror conditions, i.e.
c
c     n x E+ = -n x E-
c     n x H+ =  n x H-
c
c     (cf. Hesthaven, Warburton, HO-/SM on Unst. Grids, I. Time-Domain
c     Solution of Maxwell's Eqns., p. 14). We have
c
c     n x [E] = n x (E+ - E-) = n x -2E-
c     n x [H] = n x (H+ - H-) = 0


      if (ncempec.ge.1) then
         k = nx1*nz1*2*ndim*nelt
         if (imode.eq.3) then !IF3D
            do j = 1,ncempec
               i = cempec(j)
               srflx(0*k+i) = 2.0*srflx(0*k+i)
               srflx(1*k+i) = 2.0*srflx(1*k+i)
               srflx(2*k+i) = 2.0*srflx(2*k+i)
               srflx(3*k+i) = 0
               srflx(4*k+i) = 0
               srflx(5*k+i) = 0
            enddo
         elseif (imode.eq.2) then ! IFTM
            do j = 1,ncempec
               i = cempec(j)
               srflx(0*k+i) = 2.0*srflx(0*k+i)
               srflx(1*k+i) = 2.0*srflx(1*k+i)
               srflx(2*k+i) = 0
            enddo
         elseif (imode.eq.1) then ! IFTE
            do j = 1,ncempec
               i = cempec(j)
               srflx(0*k+i) = 0
               srflx(1*k+i) = 0
               srflx(2*k+i) = 2.0*srflx(2*k+i)
            enddo
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine maxwell_wght_curl(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real    w1(lpts),w2(lpts),w3(lpts),u1(lpts),u2(lpts),u3(lpts)

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw

      real    rxmn_k, sxmn_k, txmn_k
      real    rymn_k, symn_k, tymn_k
      real    rzmn_k, szmn_k, tzmn_k
      integer e,i,j,k, nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)

      real u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t

      nn    = nx1-1

      if (if3d) then

         do e=1,nelt

            j = nxyz*(e-1)

            call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
            call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
            call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

            do i=1,nxyz

               k = i+j

               u1rw = u1r(i) * w3mn(i)
               u1sw = u1s(i) * w3mn(i)
               u1tw = u1t(i) * w3mn(i)
               u2rw = u2r(i) * w3mn(i)
               u2sw = u2s(i) * w3mn(i)
               u2tw = u2t(i) * w3mn(i)
               u3rw = u3r(i) * w3mn(i)
               u3sw = u3s(i) * w3mn(i)
               u3tw = u3t(i) * w3mn(i)

               rxmn_k= rxmn(k)
               sxmn_k= sxmn(k)
               txmn_k= txmn(k)
               rymn_k= rymn(k)
               symn_k= symn(k)
               tymn_k= tymn(k)
               rzmn_k= rzmn(k)
               szmn_k= szmn(k)
               tzmn_k= tzmn(k)

               w1(k)= u3rw * rymn_k + u3sw * symn_k + u3tw * tymn_k
     $              - u2rw * rzmn_k - u2sw * szmn_k - u2tw * tzmn_k
               w2(k)= u1rw * rzmn_k + u1sw * szmn_k + u1tw * tzmn_k
     $              - u3rw * rxmn_k - u3sw * sxmn_k - u3tw * txmn_k
               w3(k)= u2rw * rxmn_k + u2sw * sxmn_k + u2tw * txmn_k
     $              - u1rw * rymn_k - u1sw * symn_k - u1tw * tymn_k

            enddo

         enddo
      else

         do e= 1,nelt

            j= nxyz*(e-1)

            call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
            call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1)
            call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

            do i=1,nxyz

               k = i+j

c...         sgn1/sgn2 is assigned "-" or "+" when calling; sgn1/sgn2=1/-1 if TM,
c...         sgn1/sgn2 is assigned "-" or "+" when calling; sgn1/sgn2=-1/1 if TE,

               u1rw = u1r(i) * w3mn(i)
               u1sw = u1s(i) * w3mn(i)
               u2rw = u2r(i) * w3mn(i)
               u2sw = u2s(i) * w3mn(i)

               u3rw = u3r(i) * w3mn(i)
               u3sw = u3s(i) * w3mn(i)

               rxmn_k= rxmn(k)
               sxmn_k= sxmn(k)
               rymn_k= rymn(k)
               symn_k= symn(k)

               w1(k) = ( u3rw * rymn_k + u3sw * symn_k )
               w2(k) =-( u3rw * rxmn_k + u3sw * sxmn_k )
               w3(k) = ( u2rw * rxmn_k + u2sw * sxmn_k
     $                 - u1rw * rymn_k - u1sw * symn_k )

            enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine maxwell_wght_dcurl(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
c     Combined weighted:   (w1,w2,w3) = B*curl(u1,u2,u3) on a finer mesh
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

      real    rxmn_k, sxmn_k, txmn_k
      real    rymn_k, symn_k, tymn_k
      real    rzmn_k, szmn_k, tzmn_k

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
c     This will need to be dynamic if we have p-refinement
         do e=1,nelt

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

         call intp_rstd (ju1,u1(j1),mx,md,if3d, 0) ! 0 ==> forward
         call intp_rstd (ju2,u2(j1),mx,md,if3d, 0)
         call intp_rstd (ju3,u3(j1),mx,md,if3d, 0)

         if (if3d) then

            call local_grad3 (ud1r,ud1s,ud1t,ju1,mm,1,dxmd,dxtmd)
            call local_grad3 (ud2r,ud2s,ud2t,ju2,mm,1,dxmd,dxtmd)
            call local_grad3 (ud3r,ud3s,ud3t,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

               k = i+j

               u1rw = ud1r(i) * w3md(i)
               u1sw = ud1s(i) * w3md(i)
               u1tw = ud1t(i) * w3md(i)
               u2rw = ud2r(i) * w3md(i)
               u2sw = ud2s(i) * w3md(i)
               u2tw = ud2t(i) * w3md(i)
               u3rw = ud3r(i) * w3md(i)
               u3sw = ud3s(i) * w3md(i)
               u3tw = ud3t(i) * w3md(i)

               rxmn_k= rxmd(k)
               sxmn_k= sxmd(k)
               txmn_k= txmd(k)
               rymn_k= rymd(k)
               symn_k= symd(k)
               tymn_k= tymd(k)
               rzmn_k= rzmd(k)
               szmn_k= szmd(k)
               tzmn_k= tzmd(k)

               wd1(k)= u3rw * rymn_k + u3sw * symn_k + u3tw * tymn_k
     $               - u2rw * rzmn_k - u2sw * szmn_k - u2tw * tzmn_k
               wd2(k)= u1rw * rzmn_k + u1sw * szmn_k + u1tw * tzmn_k
     $               - u3rw * rxmn_k - u3sw * sxmn_k - u3tw * txmn_k
               wd3(k)= u2rw * rxmn_k + u2sw * sxmn_k + u2tw * txmn_k
     $               - u1rw * rymn_k - u1sw * symn_k - u1tw * tymn_k

            enddo


         else

            call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
            call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
            call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

               k = i+j

               u1rw = ud1r(i) * w3md(i)
               u1sw = ud1s(i) * w3md(i)
               u2rw = ud2r(i) * w3md(i)
               u2sw = ud2s(i) * w3md(i)
               u3rw = ud3r(i) * w3md(i)
               u3sw = ud3s(i) * w3md(i)

               rxmn_k= rxmd(k)
               sxmn_k= sxmd(k)
               rymn_k= rymd(k)
               symn_k= symd(k)

               wd1(k) = ( u3rw * rymn_k + u3sw * symn_k )
               wd2(k) =-( u3rw * rxmn_k + u3sw * sxmn_k )
               wd3(k) = ( u2rw * rxmn_k + u2sw * sxmn_k
     $                  - u1rw * rymn_k - u1sw * symn_k )

            enddo

         endif

         call intp_rstd (w1(j1),wd1,mx,md,if3d,1) ! 1 --> backward
         call intp_rstd (w2(j1),wd2,mx,md,if3d,1) ! 1 --> backward
         call intp_rstd (w3(j1),wd3,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine maxwell_curl(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
c     non-weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real    w1(lpts),w2(lpts),w3(lpts),u1(lpts),u2(lpts),u3(lpts)

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw

      real    rxmn_k, sxmn_k, txmn_k
      real    rymn_k, symn_k, tymn_k
      real    rzmn_k, szmn_k, tzmn_k
      integer e,i,j,k,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)

      real u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t

      nn    = nx1-1

      if (if3d) then
         do e = 1,nelt
            j = nxyz*(e-1)

            call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
            call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
            call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

            do i = 1,nxyz
               k = i+j

               u1rw = u1r(i)
               u1sw = u1s(i)
               u1tw = u1t(i)
               u2rw = u2r(i)
               u2sw = u2s(i)
               u2tw = u2t(i)
               u3rw = u3r(i)
               u3sw = u3s(i)
               u3tw = u3t(i)

               rxmn_k = rxmn(k) / jacm(k)
               sxmn_k = sxmn(k) / jacm(k)
               txmn_k = txmn(k) / jacm(k)
               rymn_k = rymn(k) / jacm(k)
               symn_k = symn(k) / jacm(k)
               tymn_k = tymn(k) / jacm(k)
               rzmn_k = rzmn(k) / jacm(k)
               szmn_k = szmn(k) / jacm(k)
               tzmn_k = tzmn(k) / jacm(k)

               w1(k)= u3rw * rymn_k + u3sw * symn_k + u3tw * tymn_k
     $              - u2rw * rzmn_k - u2sw * szmn_k - u2tw * tzmn_k
               w2(k)= u1rw * rzmn_k + u1sw * szmn_k + u1tw * tzmn_k
     $              - u3rw * rxmn_k - u3sw * sxmn_k - u3tw * txmn_k
               w3(k)= u2rw * rxmn_k + u2sw * sxmn_k + u2tw * txmn_k
     $              - u1rw * rymn_k - u1sw * symn_k - u1tw * tymn_k
            enddo
         enddo
      else
         do e = 1,nelt
            j = nxyz*(e-1)

            call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
            call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1)
            call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

            do i = 1,nxyz
               k = i+j

               u1rw = u1r(i)
               u1sw = u1s(i)
               u2rw = u2r(i)
               u2sw = u2s(i)
               u3rw = u3r(i)
               u3sw = u3s(i)

               rxmn_k= rxmn(k) / jacm(k)
               sxmn_k= sxmn(k) / jacm(k)
               rymn_k= rymn(k) / jacm(k)
               symn_k= symn(k) / jacm(k)

               w1(k) = ( u3rw * rymn_k + u3sw * symn_k )
               w2(k) =-( u3rw * rxmn_k + u3sw * sxmn_k )
               w3(k) = ( u2rw * rxmn_k + u2sw * sxmn_k
     $                 - u1rw * rymn_k - u1sw * symn_k )

            enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_invqmass
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'
      integer  i

#ifdef _OPENACC
!$ACC DATA PRESENT(resHN,resEN,hbm1,ebm1)
      if     (imode.eq.3) then !IF3D
!$ACC PARALLEL LOOP
         do i=1,npts
            resHN(i,1)=resHN(i,1)*hbm1(i)
            resHN(i,2)=resHN(i,2)*hbm1(i)
            resHN(i,3)=resHN(i,3)*hbm1(i)
            resEN(i,1)=resEN(i,1)*ebm1(i)
            resEN(i,2)=resEN(i,2)*ebm1(i)
            resEN(i,3)=resEN(i,3)*ebm1(i)
         enddo
!$ACC END PARALLEL LOOP
      elseif (imode.eq.2) then  !IFTM
!$ACC PARALLEL LOOP
         do i=1,npts
            resHN(i,1)=reshN(i,1)*hbm1(i)
            resHN(i,2)=reshN(i,2)*hbm1(i)
            resEN(i,3)=reseN(i,3)*ebm1(i)
         enddo
!$ACC END PARALLEL LOOP
      elseif (imode.eq.1) then  !IFTE
!$ACC PARALLEL LOOP
         do i=1,npts
            resEN(i,1)=resEN(i,1)*ebm1(i)
            resEN(i,2)=resEN(i,2)*ebm1(i)
            resHN(i,3)=resHN(i,3)*hbm1(i)
         enddo
!$ACC END PARALLEL LOOP
      endif

!$ACC END DATA


#else

      if     (imode.eq.3) then !IF3D
         do i=1,npts
            resHN(i,1)=resHN(i,1)*hbm1(i)
            resHN(i,2)=resHN(i,2)*hbm1(i)
            resHN(i,3)=resHN(i,3)*hbm1(i)
            resEN(i,1)=resEN(i,1)*ebm1(i)
            resEN(i,2)=resEN(i,2)*ebm1(i)
            resEN(i,3)=resEN(i,3)*ebm1(i)
         enddo
      elseif (imode.eq.2) then  !IFTM
         do i=1,npts
            resHN(i,1)=reshN(i,1)*hbm1(i)
            resHN(i,2)=reshN(i,2)*hbm1(i)
            resEN(i,3)=reseN(i,3)*ebm1(i)
         enddo
      elseif (imode.eq.1) then  !IFTE
         do i=1,npts
            resEN(i,1)=resEN(i,1)*ebm1(i)
            resEN(i,2)=resEN(i,2)*ebm1(i)
            resHN(i,3)=resHN(i,3)*hbm1(i)
         enddo
      endif

#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_maxwell_ab(ii)
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

      if (imode.eq.3) then
         call rk4_upd(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
         call rk4_upd(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
         call rk4_upd(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)
         call rk4_upd(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
         call rk4_upd(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)
         call rk4_upd(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)
      elseif (imode.eq.2) then
         call rk4_upd(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
         call rk4_upd(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
         call rk4_upd(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)
      elseif (imode.eq.1) then
         call rk4_upd(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
         call rk4_upd(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)
         call rk4_upd(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)
      endif

      if (ifpml) then
         call rk4_upd(pmlbn(1,1),kpmlbn(1,1),respmlbn(1,1),cb,ca,dt
     $        ,npts)
         call rk4_upd(pmlbn(1,2),kpmlbn(1,2),respmlbn(1,2),cb,ca,dt
     $        ,npts)
         call rk4_upd(pmlbn(1,3),kpmlbn(1,3),respmlbn(1,3),cb,ca,dt
     $        ,npts)
         call rk4_upd(pmldn(1,1),kpmldn(1,1),respmldn(1,1),cb,ca,dt
     $        ,npts)
         call rk4_upd(pmldn(1,2),kpmldn(1,2),respmldn(1,2),cb,ca,dt
     $        ,npts)
         call rk4_upd(pmldn(1,3),kpmldn(1,3),respmldn(1,3),cb,ca,dt
     $        ,npts)
      endif

      stime0=dclock()-stime0

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
         if (nid.eq.0) then
            write(6,*) 'rk_maxwell_ab                (sec)::',stime0
         endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      if (nid.eq.0) write(6,*) 'compute eigenvalues'

      call cem_maxwell_eig        !call schrod_eig_arpack
      call cem_end

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer       lpts_eig2
      parameter    (lpts_eig2=lpts_eig*lpts_eig)

      common /EIGC/ res(  lpts_eig2)
     $            , ev (  lpts_eig2)
     $            , er (  lpts_eig)
     $            , ei (  lpts_eig)
     $            , rew(4*lpts_eig)
      real          res,ev,er,ei,rew
      integer       info,i,j,jj,k,k1,nn,nn2,npts2,npts3,npts4,npts6


      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx

      nfield=3*(ndim-1)

      npts2= npts*2
      npts3= npts*3
      npts4= npts*4
      npts6= npts*6

      if (nfield*npts.gt.lpts_eig) then
         write(6,*) 'lpts_eig',lpts_eig, nfield*npts
         write(6,*) 'insufficient memory in cem_maxwell_eig, abort'
         call exitt(1)
      endif

      if (nid.eq.0) write(6,*) 'start cem_maxwell_eig',imode
      nn  = npts*nfield
      nn2 = nn*nn
      call rzero(res,nn2)

c     if3d
      if (imode.eq.3) then

         k1 = 0
         do jj= 1,nfield
            do j = 1,npts

               k1 =  k1+1
               k  = (k1-1)*nn+1

               call rzero(HN,npts3)
               call rzero(EN,npts3)

               if (jj.eq.1) HN(j,1)= 1. ! jth col of identity matrix
               if (jj.eq.2) HN(j,2)= 1. ! jth col of identity matrix
               if (jj.eq.3) HN(j,3)= 1. ! jth col of identity matrix
               if (jj.eq.4) EN(j,1)= 1. ! jth col of identity matrix
               if (jj.eq.5) EN(j,2)= 1. ! jth col of identity matrix
               if (jj.eq.6) EN(j,3)= 1. ! jth col of identity matrix

               call cem_maxwell
               call cem_maxwell_restrict_to_face
               call cem_maxwell_flux (srflx)
               call cem_maxwell_add_flux_to_res (srflx)
               !call cem_maxwell_invqmass

               call copy(res(k+0*npts),resHN(1,1),npts)
               call copy(res(k+1*npts),resHN(1,2),npts)
               call copy(res(k+2*npts),resHN(1,3),npts)
               call copy(res(k+3*npts),resEN(1,1),npts)
               call copy(res(k+4*npts),resEN(1,2),npts)
               call copy(res(k+5*npts),resEN(1,3),npts)

            enddo
         enddo

      elseif (imode.eq.2) then

         k1 = 0
         do jj = 1,nfield
            do j  = 1,npts

               k1 =  k1+1
               k  = (k1-1)*nn+1

               call rzero(HN,npts3)
               call rzero(EN,npts3)

               if (jj.eq.1) HN(j,1)= 1. ! jth col of identity matrix
               if (jj.eq.2) HN(j,2)= 1. ! jth col of identity matrix
               if (jj.eq.3) EN(j,3)= 1. ! jth col of identity matrix

               call cem_maxwell
               call cem_maxwell_restrict_to_face
               call cem_maxwell_flux (srflx)
               call cem_maxwell_add_flux_to_res (srflx)
               !call cem_maxwell_invqmass

               call copy(res(k+0*npts),resHN(1,1),npts)
               call copy(res(k+1*npts),resHN(1,2),npts)
               call copy(res(k+2*npts),resEN(1,3),npts)

            enddo
         enddo

      elseif (imode.eq.1) then

         k1 = 0
         do jj = 1,nfield
            do j  = 1,npts

               k1 =  k1+1
               k  = (k1-1)*nn+1

               call rzero(EN(1,1),npts)
               call rzero(EN(1,2),npts)
               call rzero(HN(1,3),npts)

               if (jj.eq.1) EN(j,1)= 1. ! jth col of identity matrix
               if (jj.eq.2) EN(j,2)= 1. ! jth col of identity matrix
               if (jj.eq.3) HN(j,3)= 1. ! jth col of identity matrix

               call cem_maxwell
               call cem_maxwell_restrict_to_face
               call cem_maxwell_flux (srflx)
               call cem_maxwell_add_flux_to_res (srflx)
               !call cem_maxwell_invqmass

               call copy(res(k+0*npts),resEN(1,1),npts)
               call copy(res(k+1*npts),resEN(1,2),npts)
               call copy(res(k+2*npts),resHN(1,3),npts)

            enddo
         enddo

      endif

      if (nid.eq.0) write(6,*) 'obtain spatial operator'

      do i=1,nn2
         if (abs(res(i)).lt.1.e-13) res(i)=0.0
      enddo

      open(50,file='spatial_matrix.dat')

      do i=1,nn2
         write(50,*) res(i)
      enddo
      close(50)


      if (nid.eq.0) then
         write(6,*) '                                      '
         write(6,*) 'IFEIG T :: computing eigenvalues      '
         write(6,*) '                                      '
      endif

      call dgeev('N','V',nn,res,nn,er,ei,ev,nn,ev,nn,rew,4*nn,info)

      open(60,file='eigenvalue.dat')
      do i=1,nn
         write(60,*)  er(i),ei(i)
      enddo
      close(60)

    1 format(1p1e19.12)

      if (nid.eq.0) write(6,*) 'eigenvalue computation completed!!!'
      call exitt(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_exp
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      integer  icalld
      save     icalld
      data     icalld/0/
      real     glsc2,beta,tmp
      integer  i,j,k,n6,n,m,info

      n6 = npts*(3+3*(ndim-2))
      n  = npts
      m  = larnol

      if (icalld.eq.0) then
         call rzero(ru ,n6)
         call rzero(rv0,n6*(m+1))
         call rzero(rv ,n6*m)
         call rzero(rvt,n6*m)
         call rzero(rh0,(m+1)*m)
         call rzero(rh ,m*m)
         call rzero(rw ,m)
         call rzero(wsp,2*m*(m+2))
         icalld = 1
      endif

      if     (imode.eq.3) then
         call copy(ru(0*npts+1),hn(1,1),npts)
         call copy(ru(1*npts+1),hn(1,2),npts)
         call copy(ru(2*npts+1),hn(1,3),npts)
         call copy(ru(3*npts+1),en(1,1),npts)
         call copy(ru(4*npts+1),en(1,2),npts)
         call copy(ru(5*npts+1),en(1,3),npts)
      elseif (imode.eq.2) then
         call copy(ru(0*npts+1),hn(1,1),npts)
         call copy(ru(1*npts+1),hn(1,2),npts)
         call copy(ru(2*npts+1),en(1,3),npts)
      elseif (imode.eq.1) then
         call copy(ru(0*npts+1),en(1,1),npts)
         call copy(ru(1*npts+1),en(1,2),npts)
         call copy(ru(2*npts+1),hn(1,3),npts)
      endif

      beta= glsc2(ru,ru,n6)
      beta= sqrt (beta)

      call arnoldi_vec (rv0,rh0,n6,m,ru)
      call arnoldi_mat (rvt,rv,rh,rv0,rh0,m,n6)

      call copy_r2z    (zh,rh,m*m)
      call zgeev  ('N','V',m,zh,m,eval,vri,m,vr,m,awork,2*m,rwork,info)
      call copy_z2z    (vri,vr,m*m)

      call zgetrf      (m,m,vri,m,ipiv,info) !vri= LU factorization
      call zgetri      (m,vri,m,ipiv,awork,m,info)!vri= X^(-1)

      call z_exp       (eval,dt,m)      !eval=exp(dt*eval)
      call z_rowscale  (vri,eval,m)     !vri =eval*X^(-1);
      call c_mxm       (vr,m,vri,m,zh,m)         !zh =X*vri

      call copy_z2r    (rw,zh(1,1),m)                 !zh(i,1) (i=1:m)
      call mxm         (rv,n6,rw,m,ru,1)
      call cmult       (ru,beta,n6)

      if     (imode.eq.3) then
         call copy(hn(1,1),ru(0*npts+1),npts)
         call copy(hn(1,2),ru(1*npts+1),npts)
         call copy(hn(1,3),ru(2*npts+1),npts)
         call copy(en(1,1),ru(3*npts+1),npts)
         call copy(en(1,2),ru(4*npts+1),npts)
         call copy(en(1,3),ru(5*npts+1),npts)
      elseif (imode.eq.2) then
         call copy(hn(1,1),ru(0*npts+1),npts)
         call copy(hn(1,2),ru(1*npts+1),npts)
         call copy(en(1,3),ru(2*npts+1),npts)
      elseif (imode.eq.1) then
         call copy(en(1,1),ru(0*npts+1),npts)
         call copy(en(1,2),ru(1*npts+1),npts)
         call copy(hn(1,3),ru(2*npts+1),npts)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine arnoldi_mat (vt,v,h,v0,h0,m,n)
      implicit none
      integer  m,n,i,j,k
      real     h0(m+1,m),h(m,m)
      real     v0(n,m+1)
      real     v (n,m)
      real     vt(m,n)

      do j=1,m
         do i=1,n
            v (i,j)= v0(i,j)
            vt(j,i)= v0(i,j)
         enddo
      enddo

      do i=1,m
         do j=1,m
            h(i,j) = h0(i,j)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine arnoldi_vec(v,h,n,m,b)
      implicit none
      integer  i,j,m,n
      real     v(n,m+1),b(n)
      real     h(m+1,m)
      real     eps, s1, si, glsc2, cnt

      eps = 1.e-13

      s1 = glsc2(b,b,n)
      si = 1./sqrt(s1) ! should add a check and a comment here

      call cmult2(v(1,1),b,si,n)

      do j=1,m

         call amult (v(1,j+1),v(1,j),n)

         do i=1,j
            h(i,j) = glsc2 (v(1,j+1),v(1,i),n)
            cnt    =       -h(i,j  )
            call add2s2    (v(1,j+1),v(1,i),cnt,n)
         enddo

         s1 = glsc2(v(1,j+1),v(1,j+1),n)
         s1 = sqrt(s1)
         h(j+1,j) = s1

         if (abs(h(j+1,j)).le. eps) then
            write(6,*) 'h', h(j+1,j)
            call exitt(1)
            return              !  should add a comment here
         else
            si = 1.0/s1
            call cmult(v(1,j+1),si,n)
         endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine amult(w,v,n6)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      integer  n6
      real     w(n6),v(n6)
      common  /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real            srflx


      if        (imode.eq.3) then
         call copy(hn(1,1), v(1+0*npts), npts)
         call copy(hn(1,2), v(1+1*npts), npts)
         call copy(hn(1,3), v(1+2*npts), npts)
         call copy(en(1,1), v(1+3*npts), npts)
         call copy(en(1,2), v(1+4*npts), npts)
         call copy(en(1,3), v(1+5*npts), npts)
      elseif    (imode.eq.2) then
         call copy(hn(1,1), v(1+0*npts), npts)
         call copy(hn(1,2), v(1+1*npts), npts)
         call copy(en(1,3), v(1+2*npts), npts)
      elseif    (imode.eq.1) then
         call copy(en(1,1), v(1+0*npts), npts)
         call copy(en(1,2), v(1+1*npts), npts)
         call copy(hn(1,3), v(1+2*npts), npts)
      endif

      call cem_maxwell
      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux (srflx)
      call cem_maxwell_add_flux_to_res (srflx)
      call usersrc(rktime,reshn,resen)
      call cem_maxwell_invqmass

      if        (imode.eq.3) then
         call copy(w(1+0*npts), reshn(1,1), npts)
         call copy(w(1+1*npts), reshn(1,2), npts)
         call copy(w(1+2*npts), reshn(1,3), npts)
         call copy(w(1+3*npts), resen(1,1), npts)
         call copy(w(1+4*npts), resen(1,2), npts)
         call copy(w(1+5*npts), resen(1,3), npts)
      elseif    (imode.eq.2) then
         call copy(w(1+0*npts), reshn(1,1), npts)
         call copy(w(1+1*npts), reshn(1,2), npts)
         call copy(w(1+2*npts), resen(1,3), npts)
      elseif    (imode.eq.1) then
         call copy(w(1+0*npts), resen(1,1), npts)
         call copy(w(1+1*npts), resen(1,2), npts)
         call copy(w(1+2*npts), reshn(1,3), npts)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_r2z(z,r,n)
      implicit none
      integer    n,i
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         z(i)=r(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2r(r,z,n)
      implicit none
      integer    n,i
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         r(i)=real(z(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2z(z2,z1,n)
      implicit none
      integer    n,i
      complex*16 z1(n),z2(n)

      do i=1,n
         z2(i)=z1(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c_mxm(a,n1,b,n2,c,n3)
c-----------------------------------------------------------------------
c     Complex matrix multiplication
      implicit none
      include 'SIZE'
      include 'PARALLEL'

      integer n1,n2,n3,i

      complex a(1),b(1),c(1)
      real c_one(2),c_zero(2)

      save c_one,c_zero
      data c_one  / 1. , 0. /
      data c_zero / 0. , 0. /

      call zgemm( 'N','N',n1,n3,n2,c_one,a,n1,b,n2,c_zero,c,n1)

      return
      end
c-----------------------------------------------------------------------
      subroutine z_rowscale(a,s,n)
c-----------------------------------------------------------------------
      implicit none
      integer    n,i,j
      complex*16 a(n,n),s(n)

      do i=1,n
         do j=1,n
            a(i,j) = s(i)*a(i,j)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_exp(f,dt,n)
c-----------------------------------------------------------------------
      implicit none
      integer    n,i
      real*8     dt
      complex*16 f(n),a

      do i=1,n
         a    = dt*f(i)
         f(i) = exp(a)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_rowscale2(a,b,c,dt,n)
c-----------------------------------------------------------------------
      implicit none
      integer    n,i
      complex*16 a(n),b(n),c(n)
      real       dt

      do i=1,n
         a(i)=b(i)*exp(c(i)*dt)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wght_div (w0,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     div:     w0 = div(u)=du1/dx+du2/dy+du3/dz; u=(u1,u2,u3)
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),u1(1),u2(1),u3(1)
      integer  e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real  ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3

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

               w0(k) = (ur1(i)*rxmn(k)+ us1(i)*sxmn(k)+ ut1(i)*txmn(k))
     $              + (ur2(i)*rymn(k)+ us2(i)*symn(k)+ ut2(i)*tymn(k))
     $              + (ur3(i)*rzmn(k)+ us3(i)*szmn(k)+ ut3(i)*tzmn(k))
            enddo

         enddo

      else

         do e=1,nelt

            j = nxyz*(e-1)

            call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
            call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

            do i=1,nxyz

               k = i+j

               w0(k) = (  ur1(i)*rxmn(k) + us1(i)*sxmn(k) ) ! du/dx
     $              + (  ur2(i)*rymn(k) + us2(i)*symn(k) ) ! du/dy
            enddo

         enddo

      endif

      return
      end

#ifdef _OPENACC
c-----------------------------------------------------------------------
      subroutine curl_acc(u1r, u1s, u1t, u2r, u2s, u2t, u3r, u3s, u3t,
     $                    rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn,
     $                    w1,  w2,  w3,  w3mn,if3d)
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'
      include 'WZ'
      logical if3d
      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)
      real w1 (lx1,ly1,lz1,lelt)
      real w2 (lx1,ly1,lz1,lelt)
      real w3 (lx1,ly1,lz1,lelt)

      real w3mn(lx1,ly1,lz1)
      real rxmn(lx1,ly1,lz1,lelt)
      real sxmn(lx1,ly1,lz1,lelt)
      real txmn(lx1,ly1,lz1,lelt)
      real rymn(lx1,ly1,lz1,lelt)
      real symn(lx1,ly1,lz1,lelt)
      real tymn(lx1,ly1,lz1,lelt)
      real rzmn(lx1,ly1,lz1,lelt)
      real szmn(lx1,ly1,lz1,lelt)
      real tzmn(lx1,ly1,lz1,lelt)

      integer i,j,k,l,e


!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,w1,w2,w3)
!$ACC&     PRESENT(w3mn,rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn)

      if (if3d) then
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
         do e = 1,nelt
            do k = 1,nz1
               do j = 1,ny1
                  do i = 1,nx1
                     w1(i,j,k,e)= (u3r(i,j,k,e)*rymn(i,j,k,e)
     $                    + u3s(i,j,k,e)*symn(i,j,k,e)
     $                    + u3t(i,j,k,e)*tymn(i,j,k,e)
     $                    - u2r(i,j,k,e)*rzmn(i,j,k,e)
     $                    - u2s(i,j,k,e)*szmn(i,j,k,e)
     $                    - u2t(i,j,k,e)*tzmn(i,j,k,e))*w3mn(i,j,k)

                     w2(i,j,k,e)= (u1r(i,j,k,e)*rzmn(i,j,k,e)
     $                    + u1s(i,j,k,e)*szmn(i,j,k,e)
     $                    + u1t(i,j,k,e)*tzmn(i,j,k,e)
     $                    - u3r(i,j,k,e)*rxmn(i,j,k,e)
     $                    - u3s(i,j,k,e)*sxmn(i,j,k,e)
     $                    - u3t(i,j,k,e)*txmn(i,j,k,e))*w3mn(i,j,k)

                     w3(i,j,k,e)= (u2r(i,j,k,e)*rxmn(i,j,k,e)
     $                    + u2s(i,j,k,e)*sxmn(i,j,k,e)
     $                    + u2t(i,j,k,e)*txmn(i,j,k,e)
     $                    - u1r(i,j,k,e)*rymn(i,j,k,e)
     $                    - u1s(i,j,k,e)*symn(i,j,k,e)
     $                    - u1t(i,j,k,e)*tymn(i,j,k,e))*w3mn(i,j,k)
                  enddo
               enddo
            enddo
         enddo
!$ACC END PARALLEL LOOP
      else
!$ACC PARALLEL LOOP COLLAPSE(3) GANG WORKER VECTOR
         do e = 1,nelt
            do j = 1,ny1
               do i = 1,nx1
                  w1(i,j,1,e)= (u3r(i,j,1,e)*rymn(i,j,1,e)
     $                 + u3s(i,j,1,e)*symn(i,j,1,e))*w3mn(i,j,1)

                  w2(i,j,1,e)= (
     $                 - u3r(i,j,1,e)*rxmn(i,j,1,e)
     $                 - u3s(i,j,1,e)*sxmn(i,j,1,e))*w3mn(i,j,1)

                  w3(i,j,1,e)= (u2r(i,j,1,e)*rxmn(i,j,1,e)
     $                 + u2s(i,j,1,e)*sxmn(i,j,1,e)
     $                 - u1r(i,j,1,e)*rymn(i,j,1,e)
     $                 - u1s(i,j,1,e)*symn(i,j,1,e))*w3mn(i,j,1)
               enddo
            enddo
         enddo
!$ACC END PARALLEL LOOP
      endif


!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine maxwell_wght_curl_acc(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      include 'SIZE'
      include 'TOTAL'

      real    w1(lpts),w2(lpts),w3(lpts),u1(lpts),u2(lpts),u3(lpts)

      common /ctmp1/ u1r(lpts1),u1s(lpts1),u1t(lpts1)
     $             , u2r(lpts1),u2s(lpts1),u2t(lpts1)
     $             , u3r(lpts1),u3s(lpts1),u3t(lpts1)
      real           u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t


!$ACC DATA PRESENT(w1,w2,w3,u1,u2,u3)
!$ACC&     PRESENT(w3mn,rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn)
!$ACC&     PRESENT(dxm1,dxtm1)
!$ACC&     PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)

c     call test_global_curl_grad3_acc
      if (if3d) then
         call global_curl_grad3_acc
     $        (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,u1,u2,u3,dxm1)
      else
         call global_curl_grad2_acc
     $        (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,u1,u2,u3,dxm1)
      endif

c     call test_curl_acc(u1r, u1s, u1t, u2r, u2s, u2t, u3r, u3s, u3t,
      call curl_acc (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,
     $               rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn,
     $               w1,w2,w3,w3mn,if3d)
!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine global_curl_grad3_acc
     $           (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,u1,u2,u3,d)
c-----------------------------------------------------------------------
      include 'SIZE'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)

      real u1 (lx1,ly1,lz1,lelt)
      real u2 (lx1,ly1,lz1,lelt)
      real u3 (lx1,ly1,lz1,lelt)

      real d(lx1,ly1)

      real tmpr1,tmps1,tmpt1
      real tmpr2,tmps2,tmpt2
      real tmpr3,tmps3,tmpt3
      integer i,j,k,l,e

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)
!$ACC&    PRESENT(u1,u2,u3,d)

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&     private(tmpr1,tmpr2,tmpr3,tmps1,tmps2,tmps3,
!$ACC&             tmpt1,tmpt2,tmpt3)
!dir$ NOBLOCKING
      do e = 1,nelt
         do k = 1,nz1
            do j = 1,ny1
               do i = 1,nx1
                  tmpr1 = 0.0
                  tmpr2 = 0.0
                  tmpr3 = 0.0
                  tmps1 = 0.0
                  tmps2 = 0.0
                  tmps3 = 0.0
                  tmpt1 = 0.0
                  tmpt2 = 0.0
                  tmpt3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpr1=tmpr1+d(i,l)*u1(l,j,k,e)

                     tmpr2=tmpr2+d(i,l)*u2(l,j,k,e)
                     tmpr3=tmpr3+d(i,l)*u3(l,j,k,e)

                     tmps1=tmps1+d(j,l)*u1(i,l,k,e)

                     tmps2=tmps2+d(j,l)*u2(i,l,k,e)
                     tmps3=tmps3+d(j,l)*u3(i,l,k,e)

                     tmpt1=tmpt1+d(k,l)*u1(i,j,l,e)

                     tmpt2=tmpt2+d(k,l)*u2(i,j,l,e)
                     tmpt3=tmpt3+d(k,l)*u3(i,j,l,e)
                  enddo
                  u1r(i,j,k,e) = tmpr1
                  u2r(i,j,k,e) = tmpr2
                  u3r(i,j,k,e) = tmpr3
                  u1s(i,j,k,e) = tmps1
                  u2s(i,j,k,e) = tmps2
                  u3s(i,j,k,e) = tmps3
                  u1t(i,j,k,e) = tmpt1
                  u2t(i,j,k,e) = tmpt2
                  u3t(i,j,k,e) = tmpt3
               enddo
            enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine global_curl_grad2_acc
     $           (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,u1,u2,u3,d)
c-----------------------------------------------------------------------
      include 'SIZE'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)

      real u1 (lx1,ly1,lz1,lelt)
      real u2 (lx1,ly1,lz1,lelt)
      real u3 (lx1,ly1,lz1,lelt)

      real d(lx1,ly1)

      real tmpr1,tmps1,tmpt1
      real tmpr2,tmps2,tmpt2
      real tmpr3,tmps3,tmpt3
      integer i,j,k,l,e

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)
!$ACC&    PRESENT(u1,u2,u3,d)

!$ACC PARALLEL LOOP COLLAPSE(3) GANG WORKER VECTOR
!$ACC&     private(tmpr1,tmpr2,tmpr3,tmps1,tmps2,tmps3)
!dir$ NOBLOCKING
      do e = 1,nelt
         do j = 1,ny1
            do i = 1,nx1
               tmpr1 = 0.0
               tmpr2 = 0.0
               tmpr3 = 0.0
               tmps1 = 0.0
               tmps2 = 0.0
               tmps3 = 0.0
!$ACC LOOP SEQ
               do l=1,nx1
                  tmpr1=tmpr1+d(i,l)*u1(l,j,1,e)

                  tmpr2=tmpr2+d(i,l)*u2(l,j,1,e)
                  tmpr3=tmpr3+d(i,l)*u3(l,j,1,e)

                  tmps1=tmps1+d(j,l)*u1(i,l,1,e)

                  tmps2=tmps2+d(j,l)*u2(i,l,1,e)
                  tmps3=tmps3+d(j,l)*u3(i,l,1,e)
               enddo
               u1r(i,j,1,e) = tmpr1
               u2r(i,j,1,e) = tmpr2
               u3r(i,j,1,e) = tmpr3
               u1s(i,j,1,e) = tmps1
               u2s(i,j,1,e) = tmps2
               u3s(i,j,1,e) = tmps3
               u1t(i,j,1,e) = 0
               u2t(i,j,1,e) = 0
               u3t(i,j,1,e) = 0
            enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA

      return
      end
c----------------------------------------------------------------------
#endif
c-----------------------------------------------------------------------
      subroutine cem_3d_graphene_current
     $     (fjn,kfjn,resfjn,params,gindex,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      real fjn(lxzfl,3,6),kfjn(lxzfl,3,6),resfjn(lxzfl,3,6)
      real params(lxzfl,12)
      integer gindex(lxzfl),n

      integer i,j
      real nH_x,nH_y,nH_z,ndotE,nEn_x,nEn_y,nEn_z
      real tmp1,tmp2,tmp3,f1,f2,f3,Yfac,cpfac,jnfac
      real a_d,b_d,b_cp1,a_211,a_221,b_11,b_21,b_cp2,a_212,a_222,b_12
     $     ,b_22
      real ca,cb

      do i = 1,n
         j = gindex(i)

         a_d = params(j,1)
         b_d = params(j,2)
         b_cp1 = params(j,3)
         a_211 = params(j,4)
         a_221 = params(j,5)
         b_11 = params(j,6)
         b_21 = params(j,7)
         b_cp2 = params(j,8)
         a_212 = params(j,9)
         a_222 = params(j,10)
         b_12 = params(j,11)
         b_22 = params(j,12)
c     Three components of -n x H
         nH_x = -unym(j)*fHN(j,3)+unzm(j)*fHN(j,2)
         nH_y = unxm(j)*fHN(j,3)-unzm(j)*fHN(j,1)
         nH_z = -unxm(j)*fHN(j,2)+unym(j)*fHN(j,1)
c     Three components of n x (E x n)
         ndotE = unxm(j)*fEN(j,1)+unym(j)*fEN(j,2)+unzm(j)*fEN(j,3)
         nEn_x = fEN(j,1)-unxm(j)*ndotE
         nEn_y = fEN(j,2)-unym(j)*ndotE
         nEn_z = fEN(j,3)-unzm(j)*ndotE

         Yfac = 0.5/Y_0(j)
         tmp1 = Yfac*(nH_x+yconduc(j)*nEn_x)
         tmp2 = Yfac*(nH_y+yconduc(j)*nEn_y)
         tmp3 = Yfac*(nH_z+yconduc(j)*nEn_z)
         cpfac = b_cp1+b_cp2
         jnfac = 1.0-cpfac*Yfac
         fjn(j,1,1) = (fjn(j,1,2)+fjn(j,1,3)+fjn(j,1,5)-cpfac*tmp1)
     $        /jnfac
         fjn(j,2,1) = (fjn(j,2,2)+fjn(j,2,3)+fjn(j,2,5)-cpfac*tmp2)
     $        /jnfac
         fjn(j,3,1) = (fjn(j,3,2)+fjn(j,3,3)+fjn(j,3,5)-cpfac*tmp3)
     $        /jnfac
c     Forcing term
         f1 = tmp1-Yfac*fjn(j,1,1)
         f2 = tmp2-Yfac*fjn(j,2,1)
         f3 = tmp3-Yfac*fjn(j,3,1)
c     Drude term
         resfjn(j,1,2) = -a_d*fjn(j,1,2)+b_d*f1
         resfjn(j,2,2) = -a_d*fjn(j,2,2)+b_d*f2
         resfjn(j,3,2) = -a_d*fjn(j,3,2)+b_d*f3
c     First critical point term
         resfjn(j,1,3) = fjn(j,1,4)+b_11*f1
         resfjn(j,2,3) = fjn(j,2,4)+b_11*f2
         resfjn(j,3,3) = fjn(j,3,4)+b_11*f3
         resfjn(j,1,4) = -a_211*fjn(j,1,3)-a_221*fjn(j,1,4)+b_21*f1
         resfjn(j,2,4) = -a_211*fjn(j,2,3)-a_221*fjn(j,2,4)+b_21*f2
         resfjn(j,3,4) = -a_211*fjn(j,3,3)-a_221*fjn(j,3,4)+b_21*f3
c     Second critical point term
         resfjn(j,1,5) = fjn(j,1,6)+b_12*f1
         resfjn(j,2,5) = fjn(j,2,6)+b_12*f2
         resfjn(j,3,5) = fjn(j,3,6)+b_12*f3
         resfjn(j,1,6) = -a_212*fjn(j,1,5)-a_222*fjn(j,1,6)+b_22*f1
         resfjn(j,2,6) = -a_212*fjn(j,2,5)-a_222*fjn(j,2,6)+b_22*f2
         resfjn(j,3,6) = -a_212*fjn(j,3,5)-a_222*fjn(j,3,6)+b_22*f3
      enddo

      ca = rk4a(rkstep)
      cb = rk4b(rkstep)
      call rk4_upd(fjn(1,1,2),kfjn(1,1,2),resfjn(1,1,2),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,2),kfjn(1,2,2),resfjn(1,2,2),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,2),kfjn(1,3,2),resfjn(1,3,2),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,3),kfjn(1,1,3),resfjn(1,1,3),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,3),kfjn(1,2,3),resfjn(1,2,3),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,3),kfjn(1,3,3),resfjn(1,3,3),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,4),kfjn(1,1,4),resfjn(1,1,4),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,4),kfjn(1,2,4),resfjn(1,2,4),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,4),kfjn(1,3,4),resfjn(1,3,4),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,5),kfjn(1,1,5),resfjn(1,1,5),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,5),kfjn(1,2,5),resfjn(1,2,5),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,5),kfjn(1,3,5),resfjn(1,3,5),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,6),kfjn(1,1,6),resfjn(1,1,6),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,6),kfjn(1,2,6),resfjn(1,2,6),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,6),kfjn(1,3,6),resfjn(1,3,6),cb,ca,dt,nxzfl)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_te_graphene_current
     $     (fjn,kfjn,resfjn,params,gindex,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      real fjn(lxzfl,3,6),kfjn(lxzfl,3,6),resfjn(lxzfl,3,6)
      real params(lxzfl,12)
      integer gindex(lxzfl),n

      integer i,j
      real nH_x,nH_y,nEn_x,nEn_y
      real f1,f2,tmp1,tmp2,Yfac,cpfac,jnfac
      real a_d,b_d,b_cp1,a_211,a_221,b_11,b_21,b_cp2,a_212,a_222,b_12
     $     ,b_22
      real ca,cb

      do i = 1,n
c     Global face number for the i-th graphene face.
         j = gindex(i)

         a_d = params(j,1)
         b_d = params(j,2)
         b_cp1 = params(j,3)
         a_211 = params(j,4)
         a_221 = params(j,5)
         b_11 = params(j,6)
         b_21 = params(j,7)
         b_cp2 = params(j,8)
         a_212 = params(j,9)
         a_222 = params(j,10)
         b_12 = params(j,11)
         b_22 = params(j,12)
c     -(n x H)_x
         nH_x = -unym(j)*fHN(j,3)
c     -(n x H)_y
         nH_y = unxm(j)*fHN(j,3)
c     (n x (E x n))_x
         nEn_x = unym(j)**2*fEN(j,1)-unxm(j)*unym(j)*fEN(j,2)
c     (n x (E x n))_y
         nEn_y = unxm(j)**2*fEN(j,2)-unxm(j)*unym(j)*fEN(j,1)
         Yfac = 0.5/Y_0(j)
         tmp1 = Yfac*(nH_x+yconduc(j)*nEn_x)
         tmp2 = Yfac*(nH_y+yconduc(j)*nEn_y)
         cpfac = b_cp1+b_cp2
         jnfac = 1.0-cpfac*Yfac
         fjn(j,1,1) = (fjn(j,1,2)+fjn(j,1,3)+fjn(j,1,5)-cpfac*tmp1)
     $        /jnfac
         fjn(j,2,1) = (fjn(j,2,2)+fjn(j,2,3)+fjn(j,2,5)-cpfac*tmp2)
     $        /jnfac
c     Forcing term
         f1 = tmp1-Yfac*fjn(j,1,1)
         f2 = tmp2-Yfac*fjn(j,2,1)
c     Drude term
         resfjn(j,1,2) = -a_d*fjn(j,1,2)+b_d*f1
         resfjn(j,2,2) = -a_d*fjn(j,2,2)+b_d*f2
c     First critical point term
         resfjn(j,1,3) = fjn(j,1,4)+b_11*f1
         resfjn(j,2,3) = fjn(j,2,4)+b_11*f2
         resfjn(j,1,4) = -a_211*fjn(j,1,3)-a_221*fjn(j,1,4)+b_21*f1
         resfjn(j,2,4) = -a_211*fjn(j,2,3)-a_221*fjn(j,2,4)+b_21*f2
c     Second critical point term
         resfjn(j,1,5) = fjn(j,1,6)+b_12*f1
         resfjn(j,2,5) = fjn(j,2,6)+b_12*f2
         resfjn(j,1,6) = -a_212*fjn(j,1,5)-a_222*fjn(j,1,6)+b_22*f1
         resfjn(j,2,6) = -a_212*fjn(j,2,5)-a_222*fjn(j,2,6)+b_22*f2
      enddo

      ca = rk4a(rkstep)
      cb = rk4b(rkstep)
      call rk4_upd(fjn(1,1,2),kfjn(1,1,2),resfjn(1,1,2),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,2),kfjn(1,2,2),resfjn(1,2,2),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,3),kfjn(1,1,3),resfjn(1,1,3),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,3),kfjn(1,2,3),resfjn(1,2,3),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,4),kfjn(1,1,4),resfjn(1,1,4),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,4),kfjn(1,2,4),resfjn(1,2,4),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,5),kfjn(1,1,5),resfjn(1,1,5),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,5),kfjn(1,2,5),resfjn(1,2,5),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,1,6),kfjn(1,1,6),resfjn(1,1,6),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,2,6),kfjn(1,2,6),resfjn(1,2,6),cb,ca,dt,nxzfl)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_tm_graphene_current
     $     (fjn,kfjn,resfjn,params,gindex,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      real fjn(lxzfl,3,6),kfjn(lxzfl,3,6),resfjn(lxzfl,3,6)
      real params(lxzfl,12)
      integer gindex(lxzfl),n

      integer i,j
      real nH_z
      real Yfac,tmp,f,cpfac,jnfac
      real a_d,b_d,b_cp1,a_211,a_221,b_11,b_21,b_cp2,a_212,a_222,b_12
     $     ,b_22
      real ca,cb

      do i = 1,n
         j = gindex(i)

         a_d = params(j,1)
         b_d = params(j,2)
         b_cp1 = params(j,3)
         a_211 = params(j,4)
         a_221 = params(j,5)
         b_11 = params(j,6)
         b_21 = params(j,7)
         b_cp2 = params(j,8)
         a_212 = params(j,9)
         a_222 = params(j,10)
         b_12 = params(j,11)
         b_22 = params(j,12)
c     -(n x H)_z
         nH_z = -unxm(j)*fHN(j,2)+unym(j)*fHN(j,1)
c     Note that n x (E x n) = E in the TM case.
         Yfac = 0.5/Y_0(j)
         tmp = Yfac*(nH_z+yconduc(j)*fEN(j,3))
         cpfac = b_cp1 + b_cp2
         jnfac = 1.0-cpfac*Yfac
         fjn(j,3,1) = (fjn(j,3,2)+fjn(j,3,3)+fjn(j,3,5)-cpfac*tmp)
     $        /jnfac
c     Forcing term
         f = tmp-Yfac*fjn(j,3,1)
c     Drude term
         resfjn(j,3,2) = -a_d*fjn(j,3,2)+b_d*f
c     First critical point term
         resfjn(j,3,3) = fjn(j,3,4)+b_11*f
         resfjn(j,3,4) = -a_211*fjn(j,3,3)-a_221*fjn(j,3,4)+b_21*f
c     Second critical point term
         resfjn(j,3,5) = fjn(j,3,6)+b_12*f
         resfjn(j,3,6) = -a_212*fjn(j,3,5)-a_222*fjn(j,3,6)+b_22*f
      enddo

      ca = rk4a(rkstep)
      cb = rk4b(rkstep)
      call rk4_upd(fjn(1,3,2),kfjn(1,3,2),resfjn(1,3,2),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,3),kfjn(1,3,3),resfjn(1,3,3),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,4),kfjn(1,3,4),resfjn(1,3,4),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,5),kfjn(1,3,5),resfjn(1,3,5),cb,ca,dt,nxzfl)
      call rk4_upd(fjn(1,3,6),kfjn(1,3,6),resfjn(1,3,6),cb,ca,dt,nxzfl)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_drude(jn,kjn,resjn,params,dindex,n)
c-----------------------------------------------------------------------
c     Drude model. If
c
c     a = params(-,1)
c     b = params(-,2)
c
c     then the dielectric function is
c
c     eps(omega) = eps_oo-b/(omega*(omega+i*a))
c
c     where eps_oo is the value stored in `permittivity`.
c
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      real jn(lpts,3),kjn(lpts,3),resjn(lpts,3)
      real params(lpts,2)
      integer dindex(lpts)
      integer n

      integer i,j
      real a,b
      real ca,cb

      do i = 1,n
         j = dindex(i)
         a = params(j,1)
         b = params(j,2)

         resen(j,1) = resen(j,1)-jn(j,1)*bmn(j)
         resen(j,2) = resen(j,2)-jn(j,2)*bmn(j)
         resen(j,3) = resen(j,3)-jn(j,3)*bmn(j)

         resjn(j,1) = -a*jn(j,1)+b*en(j,1)
         resjn(j,2) = -a*jn(j,2)+b*en(j,2)
         resjn(j,3) = -a*jn(j,3)+b*en(j,3)
      enddo

      ca = rk4a(rkstep)
      cb = rk4b(rkstep)
      call rk4_upd(jn(1,1),kjn(1,1),resjn(1,1),cb,ca,dt,npts)
      call rk4_upd(jn(1,2),kjn(1,2),resjn(1,2),cb,ca,dt,npts)
      call rk4_upd(jn(1,3),kjn(1,3),resjn(1,3),cb,ca,dt,npts)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_lorentz(jn,kjn,resjn,params,lindex,n)
c-----------------------------------------------------------------------
c     Lorentz model. If
c
c     a = lorenzparams(-,1)
c     b = lorenzparams(-,2)
c     c = lorenzparams(-,3).
c
c     This corresponds to a dielectric function of
c
c     eps(omega) = eps_oo+c/(b-i*a*omega-omega**2)
c
c     where eps_oo is the value stored in `permittivity`.
c
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      real jn(lpts,3,2),kjn(lpts,3,2),resjn(lpts,3,2)
      real params(lpts,3)
      integer lindex(lpts)
      integer n

      integer i,j
      real a,b,c
      real ca,cb

      do i = 1,n
         j = lindex(i)

         a = params(j,1)
         b = params(j,2)
         c = params(j,3)

         resen(j,1) = resen(j,1)-jn(j,1,1)*bmn(j)
         resen(j,2) = resen(j,2)-jn(j,2,1)*bmn(j)
         resen(j,3) = resen(j,3)-jn(j,3,1)*bmn(j)

         resjn(j,1,1) = -a*jn(j,1,1)-b*jn(j,1,2)+c*en(j,1)
         resjn(j,2,1) = -a*jn(j,2,1)-b*jn(j,2,2)+c*en(j,2)
         resjn(j,3,1) = -a*jn(j,3,1)-b*jn(j,3,2)+c*en(j,3)

         resjn(j,1,2) = jn(j,1,1)
         resjn(j,2,2) = jn(j,2,1)
         resjn(j,3,2) = jn(j,3,1)
      enddo

      ca = rk4a(rkstep)
      cb = rk4b(rkstep)
      call rk4_upd(jn(1,1,1),kjn(1,1,1),resjn(1,1,1),cb,ca,dt,npts)
      call rk4_upd(jn(1,2,1),kjn(1,2,1),resjn(1,2,1),cb,ca,dt,npts)
      call rk4_upd(jn(1,3,1),kjn(1,3,1),resjn(1,3,1),cb,ca,dt,npts)
      call rk4_upd(jn(1,1,2),kjn(1,1,2),resjn(1,1,2),cb,ca,dt,npts)
      call rk4_upd(jn(1,2,2),kjn(1,2,2),resjn(1,2,2),cb,ca,dt,npts)
      call rk4_upd(jn(1,3,2),kjn(1,3,2),resjn(1,3,2),cb,ca,dt,npts)

      return
      end
c-----------------------------------------------------------------------
