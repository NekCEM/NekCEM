c#####################################################################
c     NEKCEM: Spectral Element/Spectral Element Discontinuous Galerkin
c     (SEDG) Code for Computational Electromagnetics developed by
c     Argonne National Laboratory
c#####################################################################

c ... NEKCEM: Main Program

      PROGRAM NEKTON

      call cem_presetup
      call cem_solve
      call cem_end

      end

c----------------------------------------------------------------------
c ... presetup: initialization
      subroutine cem_presetup
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'RTIMER'
      real*8   dclock

      call proc_init2

      call walltime(total_wtime1)
      call walltime(setup_wtime1)
      setup_time1= dclock()
      total_time1= dclock()

      call data_init2
      call read_init2
      call para_init2

      call geom_init2
      call op_counts

      call walltime(setup_wtime2)
      setup_time2= dclock()

      return
      end
c----------------------------------------------------------------------
      subroutine walltime(seconds)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RTIMER'
      real*8   omp_get_wtime,seconds
      integer  OMP_GET_THREAD_NUM
      integer  OMP_GET_NUM_THREADS
      real*8   MPI_Wtime

#ifdef MPI
      seconds = MPI_Wtime()
      if (nid.eq.0) write(6,*) 'MPI walltime:',seconds
#endif

#ifdef OMP
c$OMP PARALLEL PRIVATE(TID)
      TID    = OMP_GET_THREAD_NUM()
      TIDNUM = OMP_GET_NUM_THREADS()
      seconds= omp_get_wtime()
      if (nid.eq.0) then
         write(6,*) 'OMP walltime: nid=',nid,seconds
         if (tid.eq.0) write(6,*) 'OMP walltime: tid= 0',seconds
         if (tid.eq.1) write(6,*) 'OMP walltime: tid= 1',seconds
         if (tid.eq.2) write(6,*) 'OMP walltime: tid= 2',seconds
         if (tid.eq.3) write(6,*) 'OMP walltime: tid= 3',seconds
         if (tid.eq.4) write(6,*) 'OMP walltime: tid= 4',seconds
         if (tid.eq.5) write(6,*) 'OMP walltime: tid= 5',seconds
         if (tid.eq.6) write(6,*) 'OMP walltime: tid= 6',seconds
         if (tid.eq.7) write(6,*) 'OMP walltime: tid= 7',seconds
      endif
c$OMP END PARALLEL
#endif

      return
      end
c----------------------------------------------------------------------
      subroutine cputime(seconds,icalld)
      implicit none
      integer icalld
      real*8  dclock, seconds

      if (icalld.eq.0) then
         seconds = dclock()
      else
         seconds = dclock() - seconds
      endif

      return
      end
c----------------------------------------------------------------------
c ... main solver
      subroutine cem_solve
      implicit none
c----------------------------------------------------------------------
      include 'SIZE'
      include 'RTIMER'
      real*8   dclock,pretime,endtime

      call preprocessor_options(1)
      call walltime(solve_wtime1)
      solve_time1= dclock()

      !call cpu_time(totalu0) !totalu = secnds(0.0)

c ... main solver: spatial and temporal integration
      call cem_solver

      call walltime(solve_wtime2)
      solve_time2 = dclock()
      call preprocessor_options(2)

      !call cpu_time(totalu1)           !totalu = secnds(totalu)
      totalu = totalu1-totalu0

c     if (nid.eq.0)  write(6,*)  'cputime:',totalu0,totalu1,totalu

      return
      end
c----------------------------------------------------------------------
      subroutine preprocessor_options(ioption)
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      integer  ioption

      if (ioption.eq.0) then
         if (nid.eq.0)  write(6,*) ' Preprocessor options::  '
#ifdef NO_MPI
         if (nid.eq.0) write(6,*) '  MPI off with -DNO_MPI'
#else
         if (nid.eq.0) write(6,*) '  MPI on: -DNO_MPI not defined'
#endif
#ifdef GPU
         if (nid.eq.0) write(6,*) '  GPU on with -DGPU'
#ifdef _OPENACC
         if (nid.eq.0) write(6,*) '  GPU on with OpenACC'
#else
         if (nid.eq.0) write(6,*) '  OpenACC off: _OPENACC not defined'
#endif
#else
         if (nid.eq.0) write(6,*) '  GPU off at compile time'
#endif
      elseif (ioption.eq.1) then
#ifdef HPM
         if (nid.eq.0)  write(6,*) '  hpm on with -DHPM'
         call hpm_setup         ! check hpm libray package on system
#else
         if (nid.eq.0)  write(6,*) '  hpm off: -DHPM not defined'
#endif
      elseif (ioption.eq.2) then
#ifdef HPM
         call hpm_release
         if (nid.eq.0)  write(6,*) 'done: hpm_release'
#endif
      else

c...  additional

      endif

      return
      end
C----------------------------------------------------------------------
c ... finalize
      subroutine cem_end
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RTIMER'
      real*8   dclock,t0

      real     pretime,endtime
      real     setup_tmax,solve_tmax,total_tmax
      real     setup_tmin,solve_tmin,total_tmin
      real     cpu_t_step,comm_t_step,cpu_p_t,cpu_io_step_pts
      real     gpu_t_step,comm_t_gpu_step
      real     acc_t_step,comm_t_acc_step
      real     glsum,glmax,glmin

      istep= istep-1

c...  cputime
      total_time2 = dclock()
      setup_time = setup_time2  - setup_time1
      solve_time = solve_time2  - solve_time1
      total_time = total_time2  - total_time1

      if (nid.eq.0) then
         write(6,66) setup_time2,setup_time1,setup_time
         write(6,67) solve_time2,solve_time1,solve_time
         write(6,68) total_time2,total_time1,total_time
      endif

 66   format(' setup: time2/time1/time2-time1::',1p3e11.4,' sec ')
 67   format(' solve: time2/time1/time2-time1::',1p3e11.4,' sec ')
 68   format(' total: time2/time1/time2-time1::',1p3e11.4,' sec ')

      if (istep .gt.0) cpu_t_step = cpu_t /istep
      if (istep .gt.0) comm_t_step= comm_t/istep
      if (dumpno.gt.0) cpu_io_step= cpu_t_io/dumpno
      if (dumpno.eq.0) cpu_io_step= 0.0

      setup_tmax = glmax(setup_time,1)
      solve_tmax = glmax(solve_time,1)
      total_tmax = glmax(total_time,1)

      setup_tmin = glmin(setup_time,1)
      solve_tmin = glmin(solve_time,1)
      total_tmin = glmin(total_time,1)

      setup_time = glsum(setup_time,1)/np
      solve_time = glsum(solve_time,1)/np
      total_time = glsum(total_time,1)/np

      cpu_p_t         = glsum(cpu_t_step /npts,1)/np
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

      if (nid.eq.0) write(6,5)  cpu_t
      if (nid.eq.0) write(6,6)  comm_t
      if (nid.eq.0) write(6,7)  cpu_t_io
      if (nid.eq.0) write(6,8)  cpu_chk

      if (nid.eq.0) write(6,9)  cpu_t_step
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

    5 format('    total computation         ::',  1pe11.4,' sec ') !purely computation over timesteps
    6 format('    total communication       ::',  1pe11.4,' sec ') !communication in flux over timesteps 
    7 format('    total io time             ::',  1pe11.4,' sec ') !io time over timesteps
    8 format('    total userchk time        ::',  1pe11.4,' sec ') !io time over timesteps
    9 format('    total computation/step    ::',  1pe11.4,' sec ') !computation per step                
   10 format('    total communication /step ::',  1pe11.4,' sec ') !communication per step
   11 format('    total io time       /step ::',  1pe11.4,' sec ') !io time per step          

   12 format('    computation time/step/pts ::',  1pe11.4,' sec ') !computation per step per point 
   13 format('    io          time/step/pts ::',  1pe11.4,' sec ') !io time per step per point
   14 format('    io          swaptime      ::',  1pe11.4,' sec ') !communication vs. computation
   15 format('    io/comp     time/step/pts)::',  1pe11.4,' %   ') !communication time per step
   16 format('    communication/computation ::',  1pe11.4,' %   ') !communication vs. computation

      call exitt(0)

      return
      end
c-----------------------------------------------------------------------
c     processor initialization
      subroutine proc_init2
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      include 'OPCTR'
      include 'CTIMER'

      call iniproc2
      call opcount(1)

      return
      end
c-----------------------------------------------------------------------
c     data nitialization
      subroutine data_init2
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      real*8   dclock,t0
      if (nid.eq.0) write(6,*) 'call data_init2'

      t0= dclock()

      call initdim
      call initdat
      call files2  !read files

      if (nid.eq.0) write(6,*) 'done: data_init2::', dclock()-t0,' sec'

      return
      end
c-----------------------------------------------------------------------
c     read rea and map files
      subroutine read_init2
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      real*8   dclock,t0

      if (nid.eq.0) write(6,*) 'call read_init2'

      t0= dclock()

      if (IFHEX) call readat2
      if (nid.eq.0) write(6,*) 'done: readat2 ::',dclock()-t0,' sec'

      return
      end
c-----------------------------------------------------------------------
c     set parameters and echo
      subroutine para_init2
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'EMWAVE'
      INCLUDE 'ZPER'

      if (IFODE) return

      if (nid.eq.0) write(6,*) 'call para_init2'
      call setvar
      if (nid.eq.0) write(6,*) 'done: para_init2'

      return
      end
c-----------------------------------------------------------------------
c     geometry initialization
      subroutine geom_init2
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      integer  igeom
      integer  e

      if (IFODE) return

      igeom =  2
      call setup_topo      !glonum,dsset(skpdat),setedge,nelt,etc

      if     (IFHEX) then
         if (if3d) then
            nxyz  = nx1*ny1*nz1
            nxzf  = nx1*nz1
            nfaces= ndim*2
         else
            nxyz  = nx1*ny1
            nxzf  = nx1*nz1
            nfaces= ndim*2
         endif
      elseif (IFTET) then
         if (if3d) then
            nxyz  = nx1*(nx1+1)*(nx1+2)/6
            nxzf  = nx1*(nx1+1)/2
            nfaces= ndim+1
         else
            nxyz  = nx1*(nx1+1)/2
            nxzf  = nx1
            nfaces= ndim+1
         endif
      endif

      npts = nxyz*nelt
      nxzfl= nxzf*nfaces*nelt

      if (nid.eq.0) write(6,10) npts,nelt,nxyz,nxzfl,nfield
  10  format ('npts/nelt/nxyz/nxzfl/nfield',4i9,i5)

      if (nid.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if (nid.eq.0) write(6,'(A,/)') ' done :: usrdat'

c     if (nid.eq.0) write(6,*) 'call set_logics'
c     call set_logics      !option logical value setting
c     if (nid.eq.0) write(6,'(A,/)') ' done :: set_logics'

      if (nid.eq.0) write(6,*) 'call setlog'
      call setlog          !field initialization
      if (nid.eq.0) write(6,'(A,/)') ' done :: setlog'

      if (nid.eq.0) write(6,*) 'call gengeom'
      call gengeom (igeom) !temorary blocking for hmholtz 11/06/06 misun
      if(nid.eq.0) write(6,'(A,/)') ' done :: gengeom'

      if (nid.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat2'

      if (nid.eq.0) write(6,*) 'call geom_reset'
      call geom_reset(1)   !ic.f:  not using this for XM3
      if(nid.eq.0) write(6,'(A,/)') ' done :: geom_reset'

      if (nid.eq.0) write(6,*) 'call set_unr   '
      call set_unr         !normal, tangent vectors
      if(nid.eq.0) write(6,'(A,/)') ' done :: set_unr'

      call vrdsmsh         !validation mesh
      call bcmask          !bdry.f - boundary mapping

      if(nid.eq.0) write(6,*) 'call load_semhat_weighted'
      call load_semhat_weighted !semHAT
      if(nid.eq.0) write(6,'(A,/)') ' done :: load_semhat_weighted'

      !call setics         !validation mesh

      return
      end
c----------------------------------------------------------------------
c     countings
      subroutine op_counts
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'

      call time00
      call opcount(2)
      call dofcnt

      return
      end
c-----------------------------------------------------------------------
c     dt  setup: default
      subroutine get_dtsize(dxmin)
      implicit none
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      real     dxmin,cfl

#if defined DRIFT
      if (param(12).lt.0) then
         dt = abs(param(12))
         cfl = dt/dxmin**2
      elseif (param(12).gt.0) then
         cfl = param(12)
         dt = cfl*dxmin**2
      else
         if (nid.eq.0) write(6,*)
     $        'set param(12) with a nonzero value'
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
         if (nid.eq.0) write(6,*)
     $        'set param(12) with a nonzero value'
         call exitt(1)
      endif
#endif

      if (nid.eq.0) then
         write(6,*) ' dxmin=',dxmin
         write(6,*) 'dt=',dt, ', CFL=',cfl
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_solver
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'OPCTR'
      include 'CTIMER'
      include 'RTIMER'

      integer lxr,lyr,lzr,lxyzr,lpsc9,lzdt,lctmp1n
      parameter (lxr = 1)
      parameter (lyr = 1)
      parameter (lzr = 1)
      parameter (lxyzr = lxr*lyr*lzr)
      parameter (lpsc9 = 1)
      parameter (lzdt = 1)
      parameter (lctmp1n = lxyzr*lpsc9)

      common zdt(lzdt)
      real zdt,dxmin

c     Set parameters
      call set_logics
      call errchk_logics

c     Initialize variables, materials, and applications
#if defined SCHROD
      call cem_schrod_temporary
#elif defined DRIFT
      call cem_drift_temporary
#elif defined MAXWELL
      call cem_maxwell_temporary
#else
#error "Invalid application"
#endif

c     setting fintime, dt, cfl
      call get_dxmin (dxmin)
      call get_dtsize(dxmin)

c     initial conditions for non-restart/restart cases
 100  istep   = 0 !
      irststep= 0 ! restart step from previous run
      time    = 0.0
      call set_initial_fields

c     timestep options and constants
      call errchk_timestep_options
      call rk_storage

c     begin time iteration
#if defined SCHROD
      call time_advancing_pde
#elif defined DRIFT
      call time_advancing_pde
#elif defined MAXWELL
      call gs_dummy
      call time_advancing_pde
      call gs_dummy
#else
#error "Invalid application"
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine errchk_logics
c---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      imode=0

      if (if3d) then
         imode=3
         return
      endif

      if (iftm) imode=2
      if (ifte) imode=1

      if (imode.eq.0) then
         if (nid.eq.0) write(6,*) 'param(4): TE or TM not defiend'
         call exitt(1)
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine errchk_timestep_options
c---------------------------------------------------------------------
c     Error check to see if timestepping method is defined.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      if (ifeig.or.ifexp.or.ifrk.or.ifbdf) then
         if (nid.eq.0) write(6,*) 'param(17): defined as ', param(17)
      else
         if (nid.eq.0) write(6,*) 'param(17): not defined', param(17)
         call exitt(1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine time_advancing_pde
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'RTIMER'
      INCLUDE 'DRIFT'
      INCLUDE 'POISSON'
      INCLUDE 'BCS'
      integer  ifirst, ilast, igeom
      real     fintime

#ifdef GPU
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
      integer          cemface2, cemface_ptr, cemface_ind

!$ACC  DATA COPYIN(cemface,aream,dxm1,dxtm1,wxm1,unxm,unym,unzm)
!$ACC&      COPYIN(w3mn,rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn)
!$ACC&      COPYIN(bminv,bm1,hbm1,ebm1,Y_0,Y_1,Z_0,Z_1,XM1,YM1,ZM1)
!$ACC&      COPYIN(cemface2,cemface_ptr,cemface_ind)
!$ACC&      COPYIN(ids_lgl1,ids_ptr)
!$ACC&      CREATE(u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t)
!$ACC&      COPYIN(resHN,resEN,HN,EN,kHN,kEN,fHN,fEN,srflx)
!$ACC&      CREATE(ErrHN,ErrEN,SHN,SEN)
!$ACC&      CREATE(ug)

!$ACC&      COPYIN(DYM1,DZM1,DYTM1,DZTM1,bmn,binvm1)
!$ACC&      COPYIN(resN,resP,resE,resNQ,resPQ,resEQ)
!$ACC&      COPYIN(resNJ,resPJ,resEJ)
!$ACC&      COPYIN(mu_p,mu_n,mu_e,charge,ni,diff_n,diff_p,diff_e)
!$ACC&      COPYIN(d_permit,d_permea)
!$ACC&      COPYIN(cemface_nmn,cemface_pec,cemface_r)
!$ACC&      COPYIN(cN,cP,cE,cEE)
!$ACC&      COPYIN(cN1,cP1,cE1,scn,scp,sce,potent,potent1,rhs)
!$ACC&      CREATE(potent_b,cN_b,cP_b,cE_b)
!$ACC&      COPYIN(approx_o)
!$ACC&      COPYIN(cNJ,cPJ)
!$ACC&      COPYIN(G1M1,G2M1)
!$ACC&      COPYIN(G3M1,G4M1)
!$ACC&      COPYIN(G5M1,G6M1)


!$ACC UPDATE DEVICE(HN,EN)
#endif
#endif

c     counting and printing
      call print_steps(ifirst,ilast,fintime)

c     start time stepping timers
      cpu_t = 0.0
      acc_t = 0.0
      cpu_chk = 0.0
      gpu_t = 0.0


      do istep = ifirst,ilast
         cpu_dtime = dclock()
#if defined SCHROD
         if (ifrk) call cem_schrod_op_rk
         if (ifeig) call cem_schrod_op_eig
#elif defined DRIFT
         if (ifrk) call cem_drift_dg_op_rk
         if (ifbdf) call cem_drift_op_bdf
         if (ifeig) call cem_drift_op_eig
#elif defined MAXWELL
         if (ifrk) call cem_maxwell_op_rk
         if (ifeig) call cem_maxwell_op_eig
         if (ifexp) call cem_maxwell_op_exp
#else
#error "Invalid application"
#endif
         time = time+dt
         cpu_dtime = dclock()-cpu_dtime
         cpu_t = cpu_t+cpu_dtime

         cpu_chk = dclock()
#ifdef GPU
#ifdef _OPENACC
         call userchk_acc
#endif
#else
         call userchk
#endif

         cpu_chk = dclock()-cpu_chk
#ifndef  NO_IO
         call cem_out
#endif
         if (lastep.eq.1)     goto 1001
         if (time.ge.fintime) goto 1001
      enddo
 1000 continue
 1001 continue


#ifdef GPU
#ifdef _OPENACC
!$ACC END DATA
#endif
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine print_steps(ifirst,ilast,fintime)
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      integer  ifirst,ilast
      real     fintime

c ... irststep, nsteps: global variables
      fintime= param(10)
      nsteps = param(11)
      iocomm = param(13)
      iostep = param(15)

      if (nid.eq.0) then
         write(6,*) 'initial time=',time,', final time=',fintime
         write(6,*) 'total nsteps=',nsteps
         write(6,*) 'iocomm=',iocomm,', iostep=',iostep
      endif

      ifirst= irststep + 1 !irstep: global varale.
      ilast = nsteps + irststep

      if (nsteps.eq.0) call cem_end      !param(11)

      if (nid.eq.0) then
         write(6,*) '. '
         write(6,*) '============================'
         write(6,*) '========  BEGIN RUN ========'
         write(6,*) '============================'
         write(6,3) ifirst,ilast
      endif
    3 format(' istep=', i6,'    to',i10)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_initial_fields
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'EMWAVE'
      INCLUDE 'SCHROD'
      INCLUDE 'DRIFT'
      INCLUDE 'POISSON'
      INCLUDE 'RTIMER'

      if (ifrestart) then
c ...    read/update initial fields from previous step: restarting case
#ifndef  NO_IO
         call restart_swap
#endif
         if (nid.eq.0) write(6,*) 'restart: done, time=',time
         call userchk
         if (nid.eq.0) write(6,*) 'userchk: done, ifrestart=',ifrestart
      else

         if (ifsol) then
#if defined SCHROD
            call usersol(time,sQr,sQi,sQi,sUr,sUi,sUi)
#elif defined DRIFT
            call usersol(time,scNQ,scPQ,scPQ,scN,scP,scP)
#elif defined MAXWELL
            call usersol(time,sHN(1,1),sHN(1,2),sHN(1,3),
     $                        sEN(1,1),sEN(1,2),sEN(1,3))
#else
#error "Invalid application"
#endif
            if (nid.eq.0) write(6,*) 'usersol: done, IFSOL=',ifsol
         endif

#if defined SCHROD
         call userini(time,Qr,Qi,Qi,Ur,Ui,Ui)
         call schrod_grad (Qr(1,1),Qr(1,2),Qr(1,3),Ur)
         call schrod_grad (Qi(1,1),Qi(1,2),Qi(1,3),Ui)
#elif defined DRIFT
         call userini(time,cN,cP,cE,cN,cP,cE)
#elif defined MAXWELL
         call userini(time,HN(1,1),HN(1,2),HN(1,3),
     $                     EN(1,1),EN(1,2),EN(1,3))
#else
#error "Invalid application"
#endif
         if (nid.eq.0) write(6,*) 'userini: done, IFSOL=',ifsol
         call userchk

         if (nid.eq.0) write(6,*) 'userchk: done, ifrestart=',ifrestart
#ifndef  NO_IO
         call cem_out
#endif

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine initdim
C-------------------------------------------------------------------
C
C     Transfer array dimensions to common
C
C-------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'ZPER'
      INCLUDE 'DEALIAS'
      INCLUDE 'EXPONE'
C
      NX1=LX1
      NY1=LY1
      NZ1=LZ1
C
      NX2=LX2
      NY2=LY2
      NZ2=LZ2
C
      NX3=LX3
      NY3=LY3
      NZ3=LZ3
C
      NELT=LELT
      NELV=LELV
      NDIM=LDIM

      MARNOL=LARNOL  ! arnoldi

      NPERT = LPERT

      NXD= LXD
      NYD= LYD
      NZD= LZD

      NELX=LELX
      NELY=LELY
      NELZ=LELZ

      RETURN
      END
C
      subroutine initdat
C--------------------------------------------------------------------
C
C     Initialize and set default values.
C
C--------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      COMMON /DOIT/ IFDOIT
      LOGICAL       IFDOIT
      integer lcbc, nel8
C
C     Set default logicals
C
      IFFLOW   = .FALSE.
      IFMVBD   = .FALSE.
      IFHEAT   = .FALSE.
      IFSPLIT  = .FALSE.
      IFDOIT   = .FALSE.
      IFTET    = .FALSE.
      IFHEX    = .FALSE.
      IFDG     = .FALSE.
      IFNEDELEC= .FALSE.
      IFODE    = .FALSE.
      IFPDE    = .FALSE.
      IFCEM    = .TRUE. !default

      if (nid.eq.0) write(6,*) 'IFCEM= ',IFCEM

      call preprocessor_options(0)

      if ((level.eq.1).and.(lnumqd.eq.1).and.(lnumsp.eq.1)) then
         IFODE= .false.
         IFPDE= .true.
      else
         IFODE= .true.
         IFPDE= .false.
      endif

      if (nid.eq.0) write(6,*) '  PDE/ODE Options:: '
      if (nid.eq.0) write(6,*) '  IFPDE= ',IFPDE
      if (nid.eq.0) write(6,*) '  IFODE= ',IFODE,' (quantum solver)'

      if (mesh.eq.0) IFHEX = .true.
      if (mesh.eq.1) IFTET = .true.

      if (nid.eq.0) write(6,*) '  Mesh options:: '
      if (nid.eq.0) write(6,*) '  IFHEX= ',IFHEX,' (hexahedral )'
      if (nid.eq.0) write(6,*) '  IFTET= ',IFTET,' (tetrahedral)'

      IFDG = .true. !default
      if (nedelec.eq.1) IFNEDELEC = .true.
      if (nid.eq.0) write(6,*) '  Available basis:: '
      if (nid.eq.0) write(6,*) '  IFDG     = ',IFDG     ,' (SEDG   )'
      if (nid.eq.0) write(6,*) '  IFNEDELEC= ',IFNEDELEC,' (Nedelec)'

C     Turn off (on) diagnostics for communication
C
      IFGPRNT= .FALSE.
C
      CALL RZERO (PARAM,200)
C
C     The initialization of CBC is done in READAT
C
      LCBC = 3*6*LELT*(LDIMT1+1)
      CALL BLANK(CBC,LCBC)
C
      CALL BLANK(CCURVE ,8*LELT)
      NEL8 = 8*LELT
      CALL RZERO(XC,NEL8)
      CALL RZERO(YC,NEL8)
      CALL RZERO(ZC,NEL8)
C
      RETURN
      END
C------------------------------------------------------------------------
      subroutine setvar
C------------------------------------------------------------------------
C
C     Initialize variables
C
C------------------------------------------------------------------------
      implicit none
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'WZ'
      INCLUDE 'TSTEP'
      INCLUDE 'EMWAVE'
      integer  nfldt,mfield,maxhis,nabmsh,iadv,ifld1
      real     tlfac,one

C
C     Enforce splitting/Uzawa according to the way the code was compiled
C
      IF (LX1.EQ.LX2) IFSPLIT=.TRUE.
      IF (LX1.NE.LX2) IFSPLIT=.FALSE.


C     IF (IFHEAT) THEN
C        NFIELD = 2 + nfield
C        NFLDT  = 1 + nfield
C     ELSE
C        NFIELD = 1 + nfield
C        NFLDT  = 1 + nfiled
C     ENDIF


c     DO 100 IFIELD=MFIELD,NFIELD
c        IF (IFTMSH(IFIELD)) THEN
c            NELFLD(IFIELD) = NELT
c        ELSE
c            NELFLD(IFIELD) = NELV
c        ENDIF
c100  CONTINUE

      NMXH   = 1000
      NMXP   = 2000
      NMXE   = 1000
      NMXNL  = 100
C
      !PARAM(86) = 0
C
      BETAG  = PARAM(3)
      GTHETA = PARAM(4)
      DT     = abs(PARAM(12))
      DTINIT = DT
      FINTIM = PARAM(10)
      NSTEPS = PARAM(11)
      IOCOMM = PARAM(13)
      TIMEIO = PARAM(14)
      IOSTEP = PARAM(15)
      MAXHIS = PARAM(16)
      LASTEP = 0
      TOLPDF = abs(PARAM(21))
      TOLHDF = abs(PARAM(22))
      TOLREL = abs(PARAM(24))
      TOLABS = abs(PARAM(25))
      CTARG  = PARAM(26)
      NBDINP = PARAM(27)
      NABMSH = PARAM(28)


C     Check accuracy requested.
C
      IF (TOLREL.LE.0.) TOLREL = 0.01
C
C     Relaxed pressure iteration; maximum decrease in the residual.
C
      PRELAX = 0.1*TOLREL
      IF (.NOT.IFTRAN .AND. .NOT.IFNAV) PRELAX = 1.E-5
C
C     Tolerance for nonlinear iteration
C
      TOLNL  = 1.E-4
C
C     Fintim overrides nsteps
C
C     IF (FINTIM.NE.0.) NSTEPS = 200000  ! FIXME
c     IF (.NOT.IFTRAN ) NSTEPS = 1
C
C     Print interval defaults to 1
C
      IF (IOCOMM.EQ.0)  IOCOMM = nsteps+1
C
C     Max history dumps defaults to 1000
C
      IF (MAXHIS.LE.0 .OR. MAXHIS.GT.1000) THEN
         MAXHIS    = 1000
         PARAM(16) = ( MAXHIS )
      ENDIF
C
C     Set logical for natural convection
C
      IFNATC = .FALSE.
      IF (BETAG.GT.0.) IFNATC=.TRUE.
C
C     Set default for mesh integration scheme
C
      IF (NABMSH.LE.0 .OR. NABMSH.GT.3) THEN
         NABMSH    = NBDINP
         PARAM(28) = (NABMSH)
      ENDIF
C
C     Set default for mixing length factor
C
      TLFAC = 0.14
      IF (PARAM(49) .LE. 0.0) PARAM(49) = TLFAC
C
C     Courant number only applicable if convection in ANY field.
C
      IADV  = 0
      IFLD1 = 1
      IF (.NOT.IFFLOW) IFLD1 = 2
c     DO 200 IFIELD=IFLD1,NFIELD
c        IF (IFADVC(IFIELD)) IADV = 1
c200  CONTINUE
C
C     If characteristics, need number of sub-timesteps (DT/DS).
C     Current sub-timeintegration scheme: RK4.
C     If not characteristics, i.e. standard semi-implicit scheme,
C     check user-defined Courant number.
C
C     IF (IADV.EQ.1) CALL SETCHAR
C
C     Initialize order of time-stepping scheme (BD)
C     Initialize time step array.
C
      NBD    = 0
      CALL RZERO (DTLAG,10)
C
C     Useful constants
C
      one = 1.
      PI  = 4.*ATAN(one)

      RETURN
      END
C----------------------------------------------------------------------
      subroutine gengeom (igeom)
C------------------------------------------------------------------------
      implicit none
c
c     Generate geometry data
c

      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'GEOM'
      include 'WZ'
      include 'MASS'
      integer  igeom

      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
      real           xm3,ym3,zm3

      call genwz
      if (nid.eq.0) write(6,*) 'genwz'
      call gencoor (XM3,YM3,ZM3) !generate xm1
      if (nid.eq.0) write(6,*) 'gencoor'
      call geom1   (XM3,YM3,ZM3) !jacm1,rxm1,sxm1,rzm1,..bm1,..,area,..
      if (nid.eq.0) write(6,*) 'geom1'

c     call geom2    !xm2,ym2,....
c     CALL UPDMSYS (1)
      call volume
      if (nid.eq.0) write(6,*) 'volumn'
c     call setinvm !cem doesn't use this inverse variables

      call setdef
      if (nid.eq.0) write(6,*) 'setdef'

      if (nid.eq.0.and.istep.le.1) then
         write(6,*) 'done :: generate geomerty data'
         write(6,*) ' '
      endif

      RETURN
      END
c----------------------------------------------------------------------
      subroutine files2
C------------------------------------------------------------------------
c
c     Defines machine specific input and output file names.
c
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
C
      CHARACTER*132 NAME
      CHARACTER*1   SESS1(132),PATH1(132),NAM1(132)
      EQUIVALENCE  (SESSION,SESS1)
      EQUIVALENCE  (PATH,PATH1)
      EQUIVALENCE  (NAME,NAM1)
      CHARACTER*1  DMP(4),FLD(4),REA(4),HIS(4),SCH(4) ,ORE(4), NRE(4)
      CHARACTER*1  RE2(4)
      CHARACTER*4  DMP4  ,FLD4  ,REA4  ,HIS4  ,SCH4   ,ORE4  , NRE4
      CHARACTER*4  RE24
      EQUIVALENCE (DMP,DMP4), (FLD,FLD4), (REA,REA4), (HIS,HIS4)
     $          , (SCH,SCH4), (ORE,ORE4), (NRE,NRE4)
     $          , (RE2,RE24)
      DATA DMP4,FLD4,REA4 /'.dmp','.fld','.rea'/
      DATA HIS4,SCH4      /'.his','.sch'/
      DATA ORE4,NRE4      /'.ore','.nre'/
      DATA RE24           /'.re2'       /
      CHARACTER*78  STRING
C
C     Find out the session name:
C
      CALL BLANK(SESSION,132)
      CALL BLANK(PATH   ,132)

      ierr = 0

c...  FIXME: not sure yet if want to read .rea in one node. misun 7/13/10
      IF(NID.EQ.0) THEN

         OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD')
         READ(8,10) SESSION
         READ(8,10) PATH
         CLOSE(UNIT=8)
 10      FORMAT(A132)

         GOTO 23
 24      ierr = 1

      ENDIF
  23  call err_chk(ierr,' Cannot open SESSION.NAME!$')

      len = ltrunc(path,132)
      if(indx1(path1(len),'/',1).lt.1) then
         call chcopy(path1(len+1),'/',1)
      endif

      call bcast(SESSION,132*CSIZE)
      call bcast(PATH,132*CSIZE)

      CALL BLANK(REAFLE,132)
      CALL BLANK(RE2FLE,132)
      CALL BLANK(FLDFLE,132)
      CALL BLANK(HISFLE,132)
      CALL BLANK(SCHFLE,132)
      CALL BLANK(DMPFLE,132)
      CALL BLANK(OREFLE,132)
      CALL BLANK(NREFLE,132)
      CALL BLANK(NAME  ,132)
C
C     Construct file names containing full path to host:
C
      LS=LTRUNC(SESSION,132)
      LPP=LTRUNC(PATH,132)
      LSP=LS+LPP
c
      call chcopy(nam1(    1),path1,lpp)
      call chcopy(nam1(lpp+1),sess1,ls )
      l1 = lpp+ls+1
      ln = lpp+ls+4
c
c
c .rea file
      call chcopy(nam1  (l1),rea , 4)
      call chcopy(reafle    ,nam1,ln)
c      write(6,*) 'reafile:',reafle
c
c .re2 file
      call chcopy(nam1  (l1),re2 , 4)
      call chcopy(re2fle    ,nam1,ln)
c
c .fld file
      call chcopy(nam1  (l1),fld , 4)
      call chcopy(fldfle    ,nam1,ln)
c
c .his file
      call chcopy(nam1  (l1),his , 4)
      call chcopy(hisfle    ,nam1,ln)
c
c .sch file
      call chcopy(nam1  (l1),sch , 4)
      call chcopy(schfle    ,nam1,ln)
c
c
c .dmp file
      call chcopy(nam1  (l1),dmp , 4)
      call chcopy(dmpfle    ,nam1,ln)
c
c .ore file
      call chcopy(nam1  (l1),ore , 4)
      call chcopy(orefle    ,nam1,ln)
c
c .nre file
      call chcopy(nam1  (l1),nre , 4)
      call chcopy(nrefle    ,nam1,ln)
c
C     Write the name of the .rea file to the logfile.
C
      IF (NID.EQ.0) THEN
         CALL CHCOPY(STRING,REAFLE,78)
         WRITE(6,1000) STRING
         WRITE(6,1001)
 1000    FORMAT(//,2X,'Beginning session:',/,2X,A78)
 1001    FORMAT(/,' ')
      ENDIF
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine err_chk(ierr,string)
C------------------------------------------------------------------------
      character*1 string(132)
      character*1 ostring(132)
      character*10 s10
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      ierr = iglsum(ierr,1)
      if(ierr.eq.0) return

      len = indx1(string,'$',1)
      call blank(ostring,132)
      write(s10,11) ierr
   11 format(1x,' ierr=',i3)

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
c         call gdsum(dhc,1,dwork)
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
c     INCLUDE 'SOLN'
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
         if (if3d) call map13_all(zm3,zm1)
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
