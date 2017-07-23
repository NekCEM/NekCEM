c----------------------------------------------------------------------
c     NekCEM: Spectral Element Discontinuous Galerkin code for
c     computational electromagnetics developed by Argonne National
c     Laboratory.
c----------------------------------------------------------------------
      program nekcem
c----------------------------------------------------------------------
      implicit none

      call cem_init
      call cem_solve
      call cem_end

      end
c----------------------------------------------------------------------
      subroutine cem_init
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'RTIMER'

      real dclock
      integer igeom

      call iniproc
      call opcount(1)

      call walltime(total_wtime1)
      call walltime(setup_wtime1)
      setup_time1 = dclock()
      total_time1 = dclock()

      call initdim
      call initdat
      call files
      call readat

      call setvar               ! read parameters from rea file
      call set_logics           ! read parameters from rea file

      call setup_topo

      call usrdat               ! modify element vertices
      call setlog
      igeom = 2
      call gengeom(igeom)
      call usrdat2              ! modify mesh coordinates
      call geom_reset(1)
      call set_unr              ! normal and tangent vectors
      call vrdsmsh              ! validation mesh
      call bcmask               ! boundary mapping
      call load_semhat_weighted ! semHAT

      call time00
      call opcount(2)
      call dofcnt

      if (ifrk) call rk_storage

      call walltime(setup_wtime2)
      setup_time2 = dclock()

      return
      end
c----------------------------------------------------------------------
      subroutine walltime(seconds)
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RTIMER'

      real seconds

      real MPI_Wtime

#ifdef MPI
      seconds = MPI_Wtime()
      if (nid.eq.0) write(6,*) 'MPI walltime: ',seconds
#endif

      return
      end
c----------------------------------------------------------------------
      subroutine cputime(seconds,icalld)
c----------------------------------------------------------------------
      implicit none

      integer icalld
      real dclock,seconds

      if (icalld.eq.0) then
         seconds = dclock()
      else
         seconds = dclock()-seconds
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine cem_solve
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RTIMER'

      real dclock,pretime,endtime
      real dxmin

      call walltime(solve_wtime1)
      solve_time1 = dclock()

c     If requested, compute `dt` from the CFL number
      call get_dxmin(dxmin)
      call set_dt(dxmin)

c     Initial conditions for non-restart/restart cases
      istep = 0
      irststep = 0              ! restart step from previous run
      time = 0.0

#if defined SCHROD
      call cem_schrod_init
      call time_advancing_pde
#elif defined DRIFT
      call cem_drift_init
      call time_advancing_pde
#elif defined MAXWELL
      call cem_maxwell_init
      call acc_copy_all_in
      call cem_maxwell_init_fields
      call time_advancing_pde
      call acc_copy_all_out
#else
#error "Invalid application"
#endif

      call walltime(solve_wtime2)
      solve_time2 = dclock()

      return
      end
c----------------------------------------------------------------------
      subroutine preprocessor_options
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'

      if (nid.eq.0) then
         write(6,*) ' Preprocessor options:  '
#ifdef NO_MPI
         write(6,*) '  MPI off with -DNO_MPI'
#else
         write(6,*) '  MPI on: -DNO_MPI not defined'
#endif
#ifdef GPU
         write(6,*) '  GPU on with -DGPU'
#ifdef _OPENACC
         write(6,*) '  GPU on with OpenACC'
#else
         write(6,*) '  OpenACC off: _OPENACC not defined'
#endif
#else
         write(6,*) '  GPU off at compile time'
#endif
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine cem_end
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RTIMER'

      real dclock,t0
      real pretime,endtime
      real glsum,glmax,glmin

      istep = istep-1

c     cputime
      total_time2 = dclock()
      setup_time = setup_time2-setup_time1
      solve_time = solve_time2-solve_time1
      total_time = total_time2-total_time1

      if (nid.eq.0) then
         write(6,66) setup_time2,setup_time1,setup_time
         write(6,67) solve_time2,solve_time1,solve_time
         write(6,68) total_time2,total_time1,total_time
      endif
 66   format(' setup: time2/time1/time2-time1::',1p3e11.4,' sec ')
 67   format(' solve: time2/time1/time2-time1::',1p3e11.4,' sec ')
 68   format(' total: time2/time1/time2-time1::',1p3e11.4,' sec ')

      if (istep.gt.0) cpu_t_step = cpu_t /istep
      if (istep.gt.0) comm_t_step = comm_t/istep
      if (dumpno.gt.0) cpu_io_step = cpu_t_io/dumpno
      if (dumpno.eq.0) cpu_io_step = 0.0

      setup_tmax = glmax(setup_time,1)
      solve_tmax = glmax(solve_time,1)
      total_tmax = glmax(total_time,1)

      setup_tmin = glmin(setup_time,1)
      solve_tmin = glmin(solve_time,1)
      total_tmin = glmin(total_time,1)

      setup_time = glsum(setup_time,1)/np
      solve_time = glsum(solve_time,1)/np
      total_time = glsum(total_time,1)/np

      cpu_p_t = glsum(cpu_t_step /npts,1)/np
      cpu_io_step_pts = glsum(cpu_io_step/npts,1)/np

c     Average total computation time per timestep.
      cpu_t_step = glsum(cpu_t_step ,1)/np
c     Average total communication time per timestep.
      comm_t_step = glsum(comm_t_step,1)/np
c     Average total io time per iostep.
      cpu_io_step = glsum(cpu_io_step,1)/np

c     Average total computation time over all timesteps.
      cpu_t = glsum(cpu_t,1)/np
c     Average total communication time over all timesteps.
      comm_t = glsum(comm_t,1)/np
c     Average total io time over all timesteps.
      cpu_t_io = glsum(cpu_t_io,1)/np
c     Average total io time over all timesteps.
      cpu_chk = glsum(cpu_chk,1)/np
c     Average total io time over all timesteps.
      swaptime = glsum(swaptime,1)/np

      if (nid.eq.0) write(6,*)
     $ '==============================================================='

      if (nid.eq.0) write(6,5) cpu_t
      if (nid.eq.0) write(6,6) comm_t
      if (nid.eq.0) write(6,7) cpu_t_io
      if (nid.eq.0) write(6,8) cpu_chk

      if (nid.eq.0) write(6,9) cpu_t_step
      if (nid.eq.0) write(6,10) comm_t_step
      if (nid.eq.0) write(6,11) cpu_io_step

      if (nid.eq.0) write(6,12) cpu_p_t
      if (nid.eq.0) write(6,13) cpu_io_step_pts
      if (nid.eq.0) write(6,14) swaptime

      if (nsteps.ne.0) then
         if (nid.eq.0) write(6,15) cpu_io_step_pts/cpu_p_t
         if (nid.eq.0) write(6,16) comm_t/cpu_t*100.0
      endif

      if (nid.eq.0) write(6,*)
     $ '==============================================================='

    5 format('    total computation         ::',  1pe11.4,' sec ')
    6 format('    total communication       ::',  1pe11.4,' sec ')
    7 format('    total io time             ::',  1pe11.4,' sec ')
    8 format('    total userchk time        ::',  1pe11.4,' sec ')
    9 format('    total computation   /step ::',  1pe11.4,' sec ')
   10 format('    total communication /step ::',  1pe11.4,' sec ')
   11 format('    total io time       /step ::',  1pe11.4,' sec ')
   12 format('    computation time/step/pts ::',  1pe11.4,' sec ')
   13 format('    io          time/step/pts ::',  1pe11.4,' sec ')
   14 format('    io          swaptime      ::',  1pe11.4,' sec ')
   15 format('    io/comp     time/step/pts)::',  1pe11.4,' %   ')
   16 format('    communication/computation ::',  1pe11.4,' %   ')

      call exitt(0)

      return
      end
c-----------------------------------------------------------------------
      subroutine setvar
c----------------------------------------------------------------------
c     Initialize variables
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WZ'
      include 'TSTEP'

      integer maxhis

c     Enforce splitting/Uzawa according to the way the code was compiled
      if (lx1.eq.lx2) ifsplit = .true.
      if (lx1.ne.lx2) ifsplit = .false.

      fintim = param(10)
      nsteps = param(11)
      dt = abs(param(12))
      iocomm = param(13)
      timeio = param(14)
      iostep = param(15)
      maxhis = param(16)
      lastep = 0
      nbdinp = param(27)

c     If positive, `fintim` overrides nsteps.
      if (fintim.gt.0.0) then
         usesteps = .false.
      else
         usesteps = .true.
      endif

c     Print interval defaults to 1
      if (iocomm.eq.0) iocomm = nsteps+1

c     Max history dumps defaults to 1000
      if (maxhis.le.0.or.maxhis.gt.1000) then
         maxhis = 1000
         param(16) = maxhis
      endif

      nbd = 0

c     Useful constant
      pi = 4.0*atan(1.0)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_dt(dxmin)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real dxmin

      real cfl

#if defined DRIFT
      if (param(12).lt.0) then
         dt = abs(param(12))
         cfl = dt/dxmin**2
      elseif (param(12).gt.0) then
         cfl = param(12)
         dt = cfl*dxmin**2
      else
         if (nid.eq.0) write(6,*) 'set param(12) with a nonzero value'
         call exitt(1)
      endif
#else
      if (param(12).lt.0) then
         dt = abs(param(12))
         cfl = dt/dxmin
      elseif (param(12).gt.0) then
         cfl = param(12)
         dt = cfl*dxmin
      else
         if (nid.eq.0) write(6,*) 'set param(12) with a nonzero value'
         call exitt(1)
      endif
#endif

      if (nid.eq.0) then
         write(6,*) ' dxmin = ',dxmin
         write(6,*) 'dt = ',dt, ', cfl = ',cfl
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine acc_copy_all_in
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

#ifdef _OPENACC
      common /ctmp1/   u1r(lpts1),u1s(lpts1),u1t(lpts1)
     $               , u2r(lpts1),u2s(lpts1),u2t(lpts1)
     $               , u3r(lpts1),u3s(lpts1),u3t(lpts1)
      real             u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t

      common /gsh_tmp/ ug(lpts*12)
      real             ug

      common /scrns/   srflx(6*2*ldim*lx1*lz1*lelt)
      real             srflx

      common /gsh_acc/ gsh_face_acc,ids_lgl1(-1:lpts*2)
     $               , ids_ptr(-1:lpts*2)
      integer gsh_face_acc,ids_lgl1,ids_ptr

      common /cemfce_acc/
     $                 cemface2(2*ldim*lx1*lz1*lelt)
     $               , cemface_ptr(0:2*ldim*lx1*lz1*lelt)
     $               , cemface_ind(2*ldim*lx1*lz1*lelt)
      integer          cemface2,cemface_ptr,cemface_ind

      common /bdry1/ cempec(6*lx1*lz1*2*ldim*lelt),ncempec
      integer cempec,ncempec

!$ACC ENTER DATA COPYIN (xm1,ym1,zm1)
!$ACC ENTER DATA COPYIN (dxm1,dxtm1,w3mn,w3m1)
!$ACC ENTER DATA COPYIN (jacm1,jacm)
!$ACC ENTER DATA COPYIN (unx,uny,unz)
!$ACC ENTER DATA COPYIN (rxm1,sxm1,txm1)
!$ACC ENTER DATA COPYIN (rym1,sym1,tym1)
!$ACC ENTER DATA COPYIN (rzm1,szm1,tzm1)
!$ACC ENTER DATA COPYIN (unxm,unym,unzm)
!$ACC ENTER DATA COPYIN (rxmn,sxmn,txmn)
!$ACC ENTER DATA COPYIN (rymn,symn,tymn)
!$ACC ENTER DATA COPYIN (rzmn,szmn,tzmn)
!$ACC ENTER DATA COPYIN (bmn,bminv,bm1,hbm1,ebm1)
!$ACC ENTER DATA COPYIN (Y_0,Y_1,Z_0,Z_1,yconduc)
!$ACC ENTER DATA COPYIN (permittivity,permeability)
!$ACC ENTER DATA COPYIN (pmlptr,pmlsigma)
!$ACC ENTER DATA COPYIN (respmlhn,respmlen)
!$ACC ENTER DATA COPYIN (respmlbn,respmldn)
!$ACC ENTER DATA COPYIN (kpmlbn,kpmldn)
!$ACC ENTER DATA COPYIN (pmlbn,pmldn)

!$ACC ENTER DATA COPYIN (cemface,cempec,aream,area)
!$ACC ENTER DATA COPYIN (cemface2,cemface_ptr,cemface_ind)

!$ACC ENTER DATA COPYIN (resHN,resEN,kHN,kEN)
!$ACC ENTER DATA COPYIN (ids_lgl1,ids_ptr,ug)

!$ACC ENTER DATA CREATE (u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t)
!$ACC ENTER DATA CREATE (fHN,fEN,srflx)
!$ACC ENTER DATA COPYIN (HN,EN,SHN,SEN)

c!$ACC ENTER DATA COPYIN (dym1,dzm1,dymt1,dztm1,bmn,binvm1)
c!$ACC ENTER DATA COPYIN (potent_b,cN_b,cP_b,cE_b)
c!$ACC ENTER DATA COPYIN (resN,resP,resE,resNQ,resPQ,resEQ)
c!$ACC ENTER DATA COPYIN (resNJ,resPJ,resEJ)
c!$ACC ENTER DATA COPYIN (mu_p,mu_n,mu_e,diff_n,diff_p,diff_e)
c!$ACC ENTER DATA COPYIN (charge,ni)
c!$ACC ENTER DATA COPYIN (d_permit,d_permea)
c!$ACC ENTER DATA COPYIN (cemface_nmn,cemface_pec,cemface_r)
c!$ACC ENTER DATA COPYIN (cN,cP,cE,cEE)
c!$ACC ENTER DATA COPYIN (cN1,cP1,cE1,scn,scp,sce)
c!$ACC ENTER DATA COPYIN (potent,potent1,rhs)

c!$ACC ENTER DATA COPYIN (approx_o)
c!$ACC ENTER DATA COPYIN (cNJ,cPJ)
c!$ACC ENTER DATA COPYIN (G1M1,G2M1)
c!$ACC ENTER DATA COPYIN (G3M1,G4M1)
c!$ACC ENTER DATA COPYIN (G5M1,G6M1)
#endif

      return
      end
c----------------------------------------------------------------------
      subroutine acc_copy_all_out()
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

#ifdef _OPENACC
      common /ctmp1/   u1r(lpts1),u1s(lpts1),u1t(lpts1)
     $               , u2r(lpts1),u2s(lpts1),u2t(lpts1)
     $               , u3r(lpts1),u3s(lpts1),u3t(lpts1)
      real             u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t

      common /gsh_tmp/ ug(lpts*12)
      real             ug

      common /scrns/   srflx(6*2*ldim*lx1*lz1*lelt)
      real             srflx

      common /gsh_acc/ gsh_face_acc,ids_lgl1(-1:lpts*2)
     $               , ids_ptr(-1:lpts*2)
      integer gsh_face_acc,ids_lgl1,ids_ptr

      COMMON /cemfce_acc/
     $                 cemface2(2*ldim*lx1*lz1*lelt)
     $               , cemface_ptr(0:2*ldim*lx1*lz1*lelt)
     $               , cemface_ind(2*ldim*lx1*lz1*lelt)
      integer          cemface2,cemface_ptr,cemface_ind

      common /bdry1/ cempec(6*lx1*lz1*2*ldim*lelt),ncempec
      integer cempec,ncempec

!$ACC EXIT DATA DELETE (xm1,ym1,zm1)
!$ACC EXIT DATA DELETE (dxm1,dxtm1,w3mn,w3m1)
!$ACC EXIT DATA DELETE (jacm1,jacm)
!$ACC EXIT DATA DELETE (unx,uny,unz)
!$ACC EXIT DATA DELETE (rxm1,sxm1,txm1)
!$ACC EXIT DATA DELETE (rym1,sym1,tym1)
!$ACC EXIT DATA DELETE (rzm1,szm1,tzm1)
!$ACC EXIT DATA DELETE (unxm,unym,unzm)
!$ACC EXIT DATA DELETE (rxmn,sxmn,txmn)
!$ACC EXIT DATA DELETE (rymn,symn,tymn)
!$ACC EXIT DATA DELETE (rzmn,szmn,tzmn)
!$ACC EXIT DATA DELETE (bmn,bminv,bm1,hbm1,ebm1)
!$ACC EXIT DATA DELETE (Y_0,Y_1,Z_0,Z_1,yconduc)
!$ACC EXIT DATA DELETE (permittivity,permeability)

!$ACC EXIT DATA DELETE (pmlptr,pmlsigma)
!$ACC EXIT DATA DELETE (respmlhn,respmlen)
!$ACC EXIT DATA DELETE (respmlbn,respmldn)
!$ACC EXIT DATA DELETE (respmlbn,respmldn)
!$ACC EXIT DATA DELETE (kpmlbn,kpmldn)
!$ACC EXIT DATA DELETE (pmlbn,pmldn)

!$ACC EXIT DATA DELETE (cemface,cempec,aream,area)
!$ACC EXIT DATA DELETE (cemface2,cemface_ptr,cemface_ind)
!$ACC EXIT DATA DELETE (resHN,resEN,kHN,kEN)
!$ACC EXIT DATA DELETE (ids_lgl1,ids_ptr,ug)

!$ACC EXIT DATA DELETE (u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t)
!$ACC EXIT DATA DELETE (fHN,fEN,srflx)
!$ACC EXIT DATA COPYOUT (HN,EN,SHN,SEN)

c!$ACC EXIT DATA DELETE (DYM1,DZM1,DYTM1,DZTM1,bmn,binvm1)
c!$ACC EXIT DATA COPYOUT (potent_b,cN_b,cP_b,cE_b)
c!$ACC EXIT DATA DELETE (resN,resP,resE,resNQ,resPQ,resEQ)
c!$ACC EXIT DATA DELETE (resNJ,resPJ,resEJ)
c!$ACC EXIT DATA DELETE (mu_p,mu_n,mu_e,diff_n,diff_p,diff_e)
c!$ACC EXIT DATA DELETE (charge,ni)
c!$ACC EXIT DATA DELETE (d_permit,d_permea)
c!$ACC EXIT DATA DELETE (cemface_nmn,cemface_pec,cemface_r)
c!$ACC EXIT DATA COPYOUT (cN,cP,cE,cEE)
c!$ACC EXIT DATA COPYOUT (cN1,cP1,cE1,scn,scp,sce)
c!$ACC EXIT DATA COPYOUT (potent,potent1,rhs)
c
c!$ACC EXIT DATA DELETE (approx_o)
c!$ACC EXIT DATA COPYOUT (cNJ,cPJ)
c!$ACC EXIT DATA DELETE (G1M1,G2M1)
c!$ACC EXIT DATA DELETE (G3M1,G4M1)
c!$ACC EXIT DATA DELETE (G5M1,G6M1)
#endif

      return
      end
c----------------------------------------------------------------------
      subroutine time_advancing_pde
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'RTIMER'
      real glsum
      integer ifirst,ilast

      if (.not.usesteps) then
         nsteps = int(fintim/dt)
      endif
      ifirst = irststep+1
      ilast = nsteps+irststep

      if (nid.eq.0) then
         write(6,*) '============================'
         write(6,*) '========  Begin Run ========'
         write(6,*) '============================'
         write(6,3) ifirst-1,ilast
      endif
    3 format(' istep =', i6,'   to',i10)

      call userchk
      call cem_out
      if (nsteps.eq.0) call cem_end

c     Start time stepping timers
      cpu_t = 0.0
      acc_t = 0.0
      cpu_chk = 0.0
      gpu_t = 0.0

      istep = ifirst
      do while (.true.)
         cpu_dtime = dclock()

#if defined SCHROD
         if (ifrk) call cem_schrod_op_rk
         if (ifeig) call cem_schrod_op_eig
#elif defined DRIFT
         if (ifbdf) call cem_drift_op_bdf
         if (ifeig) call cem_drift_op_eig
#elif defined MAXWELL
         if (ifrk) call cem_maxwell_op_rk
         if (ifeig) call cem_maxwell_op_eig
         if (ifexp) call cem_maxwell_op_exp
#else
#error "Invalid application"
#endif
         cpu_dtime = dclock()-cpu_dtime
         cpu_t = cpu_t+cpu_dtime
         cpu_t_step = cpu_t/istep
         cpu_p_t = glsum(cpu_t_step /npts,1)/np

         time = time+dt

         cpu_chk = dclock()
         call userchk
         cpu_chk = dclock()-cpu_chk

         call cem_out

         istep = istep+1
         if (lastep.eq.1) exit
         if (usesteps) then
            if (istep.gt.nsteps) exit
         else
            if (time.ge.fintim) exit
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine initdim
c-----------------------------------------------------------------------
c     Transfer array dimensions to common blocks.
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'ZPER'
      include 'DEALIAS'
      include 'EXPONE'

      nx1 = lx1
      ny1 = ly1
      nz1 = lz1

      nx2 = lx2
      ny2 = ly2
      nz2 = lz2

      nx3 = lx3
      ny3 = ly3
      nz3 = lz3

      nelt = lelt
      nelv = lelv
      ndim = ldim

      marnol = larnol           ! For Arnoldi

      npert = lpert

      nxd = lxd
      nyd = lyd
      nzd = lzd

      nelx = lelx
      nely = lely
      nelz = lelz

      if (ldim.eq.3) then
         nxyz = nx1*ny1*nz1
         nxzf = nx1*nz1
         nfaces = ndim*2
      else
         nxyz = nx1*ny1
         nxzf = nx1*nz1
         nfaces = ndim*2
      endif
      npts = nxyz*nelt
      nxzfl = nxzf*nfaces*nelt

      return
      end
c--------------------------------------------------------------------
      subroutine initdat
c--------------------------------------------------------------------
c     Initialize and set default values.
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /doit/ ifdoit
      logical ifdoit

      integer lcbc,nel8

c     Set default logicals
      ifflow = .false.
      ifmvbd = .false.
      ifheat = .false.
      ifsplit = .false.
      ifdoit = .false.
      iftet = .false.
      ifhex = .false.
      ifdg = .false.

      call preprocessor_options

      if (mesh.eq.0) ifhex = .true.
      if (mesh.eq.1) iftet = .true.

      if (nid.eq.0) write(6,*) '  Mesh options: '
      if (nid.eq.0) write(6,*) '  ifhex = ',ifhex,' (hexahedral)'
      if (nid.eq.0) write(6,*) '  iftet = ',iftet,' (tetrahedral)'

      ifdg = .true.

c     Turn off (on) diagnostics for communication
      ifgprnt = .false.

      call rzero(param,200)

c     The initialization of `cbc` is done in `readat`
      lcbc = 3*6*lelt*(ldimt1+1)
      call blank(cbc,lcbc)

      call blank(ccurve,8*lelt)
      nel8 = 8*lelt
      call rzero(xc,nel8)
      call rzero(yc,nel8)
      call rzero(zc,nel8)

      return
      end
c------------------------------------------------------------------------
      subroutine gengeom(igeom)
c------------------------------------------------------------------------
c     Generate geometry data.
      implicit none
      include 'SIZE'
      include 'GEOM'

      integer igeom

      COMMON /scruz/ xm3(lx3,ly3,lz3,lelt)
     $ ,             ym3(lx3,ly3,lz3,lelt)
     $ ,             zm3(lx3,ly3,lz3,lelt)
      real xm3,ym3,zm3

      call genwz
      call gencoor(xm3,ym3,zm3)
      call geom1(xm3,ym3,zm3)

      call volume

      call setdef

      return
      end
c----------------------------------------------------------------------
      subroutine files
c----------------------------------------------------------------------
c     Get the name of the `.rea` file and the `.re2` file.
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer narg,iargc

      narg = iargc()
      if (narg.gt.1) then
         if (nid.eq.0) then
            write(*,*) 'ERROR: files: too many command line arguments'
         endif
         call exitt(1)
      elseif (narg.eq.0) then
         if (nid.eq.0) then
            write(*,*) 'ERROR: files: pass read file as first argument'
         endif
         call exitt(1)
      endif
      call getarg(1,session)

      call blank(reafle,132)
      call blank(re2fle,132)
      reafle = trim(session)//'.rea'
      re2fle = trim(session)//'.re2'

      if (nid.eq.0) then
         write(*,*) ' Beginning session: ',trim(reafle)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine err_chk(ierr,string)
C------------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      character(*) string(132)
      integer ierr

      character ostring(132)
      character(10) s10
      integer iglsum,indx1,k,len

      ierr = iglsum(ierr,1)
      if (ierr.eq.0) return

      len = indx1(string,'$',1)
      call blank(ostring,132)
      write(s10,11) ierr
 11   format(1x,' ierr = ',i3)

      call chcopy(ostring,string,len-1)
      call chcopy(ostring(len),s10,10)

      if (nid.eq.0) write(6,1) (ostring(k),k=1,len+10)
    1 format('ERROR: ',132a1)

      call exitt(ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine time00
C------------------------------------------------------------------------
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'CTIMER'
C
      nmxmf=0
      nmxms=0
      ndsum=0
      nsett=0
      ncdtp=0
      npres=0
      nmltd=0
      ngsum=0
      nprep=0
      ndsnd=0
      ndadd=0
      nhmhz=0
      naxhm=0
      ngop =0
      nusbc=0
      ncopy=0
      ninvc=0
      ninv3=0
      nsolv=0
      nslvb=0
      nddsl=0
      ncrsl=0
      ndott=0
      nbsol=0
c
      tmxmf=0.0
      tmxms=0.0
      tdsum=0.0
      tdsmn=9.9e9
      tdsmx=0.0
      tsett=0.0
      tcdtp=0.0
      tpres=0.0
      tmltd=0.0
      tgsum=0.0
      tgsmn=9.9e9
      tgsmx=0.0
      tprep=0.0
      tdsnd=0.0
      tdadd=0.0
      thmhz=0.0
      taxhm=0.0
      tgop =0.0
      tusbc=0.0
      tcopy=0.0
      tinvc=0.0
      tinv3=0.0
      tsolv=0.0
      tslvb=0.0
      tddsl=0.0
      tcrsl=0.0
      tdott=0.0
      tbsol=0.0
      tbso2=0.0
      etims0= dclock()
C
      return
      end
C------------------------------------------------------------------------
      subroutine opcount(ICALL)
C------------------------------------------------------------------------
      include 'SIZE'
      include 'OPCTR'
      character*6 sname(maxrts)
      integer     ind  (maxrts)
      integer     idum (maxrts)
C
      if (icall.eq.1) then
         nrout=0
      endif
      if (icall.eq.1.or.icall.eq.2) then
         dcount = 0.0
         do i=1,maxrts
            ncall(i) = 0
            dct(i)   = 0.0
         enddo
      endif
      if (icall.eq.3) then
C
C        Sort and print out diagnostics
C
         write(6,*) nid,' opcount',dcount
         dhc = dcount
         call gop(dhc,dwork,'+  ',1)
         if (nid.eq.0) then
            write(6,*) nid,' TOTAL OPCOUNT',dhc
         endif
C
         CALL DRCOPY(rct,dct,nrout)
         CALL SORT(rct,ind,nrout)
         CALL CHSWAPR(rname,6,ind,nrout,sname)
         call iswap(ncall,ind,nrout,idum)
C
         if (nid.eq.0) then
            do i=1,nrout
               write(6,201) nid,rname(i),rct(i),ncall(i)
            enddo
  201       format(2x,' opnode',i4,2x,a6,g18.7,i12)
         endif
      endif
      return
      end
C
C------------------------------------------------------------------------
      subroutine dofcnt
C------------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      COMMON /CTMP0/ DUMMY0(LCTMP0)
      COMMON /CTMP1/ DUMMY1(LCTMP1)
      COMMON /SCRNS/ WORK(LCTMP1)
C
      ntot1=nx1*ny1*nz1*nelv
      ntot2=nx2*ny2*nz2*nelv
C
      call col3 (work,vmult,v1mask,ntot1)
      vpts = glsum(work,ntot1) + .1
      nvtot=vpts
      work(1)=ntot2
      ppts = glsum(work,1) + .1
      nptot=ppts
C
      work(1)=0.0
      do i=1,ntot1
         if (vmult(i,1,1,1).lt.0.5) work(1)=work(1)+vmult(i,1,1,1)
      enddo
      epts = glsum(work,1) + .1
      netot=epts
      if (nid.eq.0) write(6,*) 'dofs:',nvtot,nptot,netot
      return
      end
c-----------------------------------------------------------------------
      subroutine geom_reset(icall)
C------------------------------------------------------------------------
C
C     Generate geometry data
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      include 'WZ'
      include 'MASS'

c
      COMMON /scruz/ XM3 (LX1,LY1,LZ1,LELT)
     $ ,             YM3 (LX1,LY1,LZ1,LELT)
     $ ,             ZM3 (LX1,LY1,LZ1,LELT)
C
c
      integer icall

      ntot = nx1*ny1*nz1*nelv
c
      if (lx3.eq.lx1) then
         call copy(xm3,xm1,ntot)
         call copy(ym3,ym1,ntot)
         call copy(zm3,zm1,ntot)
      else
         call map13_all(xm3,xm1)
         call map13_all(ym3,ym1)
         if (ldim.eq.3) call map13_all(zm3,zm1)
      endif
c
      call geom1 (xm3,ym3,zm3)
      call geom2
c     call updmsys (1)
      call volume
      call setinvm
      call setdef
c     call sfastax
c
      return
      end
c-----------------------------------------------------------------------
      subroutine map13_all(x3,x1)
C------------------------------------------------------------------------
c
      include 'SIZE'
      include 'TOTAL'
c
      real x3(lx3,ly3,lz3,lelt)
      real x1(lx1,ly1,lz1,lelt)
c
      integer e
c
      do e=1,nelv
         call map13 (x3(1,1,1,e),x1(1,1,1,e), e)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
