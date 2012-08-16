c#####################################################################
c     NEKCEM: Spectral Element Discontinuous Galerkin (SEDG) Code for 
c     Computational Electromagnetics developed by Argonne National Lab
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
      include 'SIZE'
      include 'INPUT'
      real*8   dclock,pretime,endtime                
      common  /tottime/ setup_time,solve_time,total_time
      real              setup_time,solve_time,total_time

      common  /wtimes/  setup_wtime,solve_wtime,total_wtime
      real              setup_wtime,solve_wtime,total_wtime

 
      call walltime(total_wtime,0)
      call walltime(setup_wtime,0)

      pretime = dclock()

      call data_init2 
      call read_init2 
      call para_init2 

c     Map BCs !fill CBC arrays
      if (ifmoab) then
#ifdef MOAB
       if (nid.eq.0) write(6,*) 'call   nekMOAB_bcs'
       call nekMOAB_bcs
       if (nid.eq.0) write(6,*) 'done:: nekMOAB_bcs'
#endif
      endif

      call geom_init2
      call op_counts

      endtime = dclock()-pretime

      setup_time= endtime                 
      total_time= pretime 

      call walltime(setup_wtime,1)

      if (nid.eq.0)  write(6,3)  endtime            
    3 format(' total presetup::',1pe11.4,' sec ')

      return
      end
c----------------------------------------------------------------------
      subroutine walltime(seconds,icalld)
      real    seconds
      integer icalld                   
    
#ifdef MPI
      if (icalld.eq.0) then
          seconds = mpi_wtime ( )
      else
          seconds = mpi_wtime ( ) - seconds;
      endif    
      if (nid.eq.0) 'MPI walltime:',icalld,seconds 
#endif    
#ifdef OMP
      if (icalld.eq.0) then
          seconds = omp_get_wtime ( )
      else
          seconds = omp_get_wtime ( ) - seconds;
      endif    
      if (nid.eq.0) 'OMP walltime:',icalld,seconds 
#endif    
      
      return
      end
c----------------------------------------------------------------------
c ... main solver                      
      subroutine cem_solve                  
c----------------------------------------------------------------------
      include 'SIZE'
      real*8   dclock,pretime,endtime                
      common  /tottime/ setup_time,solve_time,total_time
      real              setup_time,solve_time,total_time

      common  /wtimes/  setup_wtime,solve_wtime,total_wtime
      real              setup_wtime,solve_wtime,total_wtime

      call walltime(solve_wtime,0)
      call preprocessor_options(1)
      
c ... main solver: spatial and temporal integration     
      pretime = dclock()

      call cem_solver

      endtime = dclock()-pretime
      solve_time = endtime                 

      call preprocessor_options(2)
      call walltime(solve_wtime,1)

      return
      end
c----------------------------------------------------------------------
      subroutine preprocessor_options(ioption)
c----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      integer  ioption

      if (nid.eq.0)  write(6,*) '\n Preprocessor options::  '
      if     (ioption.eq.0) then

#ifdef MOAB       
          IFMOAB=.true.
          if (nid.eq.0) write(6,*) '  MOAB on with -DMOAB'
#else
          IFMOAB=.false.
          if (nid.eq.0) write(6,*) '  MOAB off at compile time'
#endif


#ifdef NO_MPI     
          if (nid.eq.0) write(6,*) '  MPI off with -DNO_MPI'
#else
          if (nid.eq.0) write(6,*) '  MPI on: -DNO_MPI not defined'
#endif


#ifdef GPU        
          if (nid.eq.0) write(6,*) '  GPU on with -DGPU'
#else
          if (nid.eq.0) write(6,*) '  GPU off at compile time'
#endif


#ifdef OMP        
          if (nid.eq.0) write(6,*) '  OpenMP on with -DOMP'
          call omp_setup
#else
          if (nid.eq.0) write(6,*) '  OpenMP off: -DOMP not defined'
#endif


      elseif (ioption.eq.1) then

#ifdef HPM        
          if (nid.eq.0)  write(6,*) '  hpm on with -DHPM'
          call hpm_setup !check hpm libray package on system
#else
          if (nid.eq.0)  write(6,*) '  hpm off: -DHPM not defined'
#endif
          if (nid.eq.0)  write(6,*) '\n'


      elseif (ioption.eq.2) then

#ifdef HPM       
          call hpm_release    
          if (nid.eq.0)  write(6,*) 'done: hpm_release'   
          if (nid.eq.0)  write(6,*) '\n'
#endif

      else

c...  additional

      endif

      return
      end
C----------------------------------------------------------------------
      subroutine omp_setup
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      integer icalld
      save    icalld
      data    icalld/0/
c     integer NTHREADS,TID,TIDNUM,TIDMAX: moved into PARALLEL
      integer OMP_GET_THREAD_NUM
      integer OMP_GET_NUM_THREADS
      integer OMP_GET_MAX_THREADS
      integer OMP_IN_PARALLEL

#ifdef OMP

c$OMP PARALLEL PRIVATE(TID)
        TID    = OMP_GET_THREAD_NUM()
        TIDNUM = OMP_GET_NUM_THREADS()
        TIDMAX = OMP_GET_MAX_THREADS()
        NTHREADS= TIDNUM

        if (nid.eq.0) then 
        if (TID .eq. 0) then ! master thread
        write(*,*) "\n"
        write(*,*) "Nthreads           ::", tidnum
        write(*,*) "tid, tidnum, tidmax::", tid,tidnum,tidmax
        write(*,*) "OMP_GET_THREAD_NUM ::", OMP_GET_THREAD_NUM()
        write(*,*) "OMP_GET_NUM_THREADS::", OMP_GET_NUM_THREADS()
        write(*,*) "OMP_GET_MAX_THREADS::", OMP_GET_MAX_THREADS()
        write(*,*) "OMP_IN_PARALLEL   :: ", OMP_IN_PARALLEL()
        endif
        endif

c$OMP END PARALLEL

#endif

      return
      end
c----------------------------------------------------------------------
c ... finalize                         
      subroutine cem_end                            
c----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      real             pretime,endtime                
      common /tottime/ setup_time,solve_time,total_time
      real             setup_time,solve_time,total_time
      common /wtimes/  setup_wtime,solve_wtime,total_wtime
      real             setup_wtime,solve_wtime,total_wtime
      common /ccpu  /  cpu_t,cpu_dtime,cpu_chk                  
      common /commt /  comm_t                        
      common /iotime/  cpu_t_io, cpu_io_step, swaptime
      real             cpu_t,cpu_dtime,cpu_chk,comm_t
      real             cpu_t_io, cpu_io_step, swaptime
      common /iocnt/   dumpno
      integer          dumpno
      real*8           dclock                                

      istep= istep-1    
      call print_options(-1)

      total_time= dclock()-total_time
      call walltime(total_wtime,1)

      if (istep .gt.0) cpu_t_step = cpu_t /istep 
      if (istep .gt.0) comm_t_step= comm_t/istep  
      if (dumpno.gt.0) cpu_io_step= cpu_t_io/dumpno 
      if (dumpno.eq.0) cpu_io_step=0.0 
      
      npts=nx1*ny1*nz1*nelt

      setup_time     = glsum(setup_time,1)/np
      solve_time     = glsum(solve_time,1)/np
      total_time     = glsum(total_time,1)/np

      cpu_p_t        = glsum(cpu_t_step /npts,1)/np !avg total computation time per step per grid
      cpu_io_step_pts= glsum(cpu_io_step/npts,1)/np !avg total io time per step per grid

      cpu_t_step = glsum(cpu_t_step ,1)/np !avg total computation time per timestep
      comm_t_step= glsum(comm_t_step,1)/np !avg total communication time per timestep    
      cpu_io_step= glsum(cpu_io_step,1)/np !avg total io time per iostep

      cpu_t    = glsum(cpu_t   ,1)/np !avg total computation time over timesteps 
      comm_t   = glsum(comm_t  ,1)/np !avg total communication time over timesteps  
      cpu_t_io = glsum(cpu_t_io,1)/np !avg total io time over timesteps
      cpu_chk  = glsum(cpu_chk ,1)/np !avg total io time over timesteps
      swaptime = glsum(swaptime,1)/np !avg total io time over timesteps

      if (nid.eq.0)  write(6,*) '====================================='
      if (nid.eq.0)  write(6,*) '   Job successfully completed !!!!   '
      if (nid.eq.0)  write(6,*) '====================================='
      if (nid.eq.0)  write(6,*) 'time averaged over np=',np

      if (nid.eq.0)  write(6,33) setup_wtime   
      if (nid.eq.0)  write(6,44) solve_wtime   
      if (nid.eq.0)  write(6,55) total_wtime   

      if (nid.eq.0)  write(6,3)  setup_time   
      if (nid.eq.0)  write(6,4)  solve_time   
      if (nid.eq.0)  write(6,5)  total_time   

      if (nid.eq.0)  write(6,6)  cpu_t                         
      if (nid.eq.0)  write(6,7)  comm_t                           
      if (nid.eq.0)  write(6,8)  cpu_t_io                         
      if (nid.eq.0)  write(6,17) cpu_chk                          

      if (nid.eq.0)  write(6,9)  cpu_t_step                          
      if (nid.eq.0)  write(6,10) comm_t_step                           
      if (nid.eq.0)  write(6,11) cpu_io_step                      

      if (nid.eq.0)  write(6,12) cpu_p_t                          
      if (nid.eq.0)  write(6,13) cpu_io_step_pts                      
      if (nid.eq.0)  write(6,14) swaptime                              

      if (nsteps.ne.0) then
      if (nid.eq.0)  write(6,15) cpu_io_step_pts/cpu_p_t              
      if (nid.eq.0)  write(6,16) comm_t/cpu_t*100.0                    
      endif

      if (nid.eq.0)  write(6,*) '====================================='
      if (nid.eq.0)  write(6,*) '\n'    

#ifdef OMP
      setup_time     = glsum(setup_time,1)/(np*nthreads)
      solve_time     = glsum(solve_time,1)/(np*nthreads)
      total_time     = glsum(total_time,1)/(np*nthreads)

      cpu_p_t        = glsum(cpu_t_step /npts,1)/(np*nthreads)
      cpu_io_step_pts= glsum(cpu_io_step/npts,1)/(np*nthreads)

      cpu_t_step = glsum(cpu_t_step ,1)/(np*nthreads)
      comm_t_step= glsum(comm_t_step,1)/(np*nthreads)
      cpu_io_step= glsum(cpu_io_step,1)/(np*nthreads)

      cpu_t    = glsum(cpu_t   ,1)/(np*nthreads)
      comm_t   = glsum(comm_t  ,1)/(np*nthreads)
      cpu_t_io = glsum(cpu_t_io,1)/(np*nthreads)
      cpu_chk  = glsum(cpu_chk ,1)/(np*nthreads)
      swaptime = glsum(swaptime,1)/(np*nthreads)

      if (nid.eq.0)  write(6,*) '====================================='
      if (nid.eq.0)  write(6,*) 'time averaged over (np*nthreads)     '
      if (nid.eq.0)  write(6,*) 'np*nthreads=',np*nthreads
      if (nid.eq.0)  write(6,*) 'np,nthreads=',np,',',nthreads
      if (nid.eq.0)  write(6,*) '====================================='
      if (nid.eq.0)  write(6,33) setup_wtime   
      if (nid.eq.0)  write(6,44) solve_wtime   
      if (nid.eq.0)  write(6,55) total_wtime   

      if (nid.eq.0)  write(6,3)  setup_time, setup_wtime
      if (nid.eq.0)  write(6,4)  solve_time, solve_wtime
      if (nid.eq.0)  write(6,5)  total_time, total_wtime

      if (nid.eq.0)  write(6,6)  cpu_t
      if (nid.eq.0)  write(6,7)  comm_t
      if (nid.eq.0)  write(6,8)  cpu_t_io
      if (nid.eq.0)  write(6,17) cpu_chk

      if (nid.eq.0)  write(6,9)  cpu_t_step
      if (nid.eq.0)  write(6,10) comm_t_step
      if (nid.eq.0)  write(6,11) cpu_io_step

      if (nid.eq.0)  write(6,12) cpu_p_t
      if (nid.eq.0)  write(6,13) cpu_io_step_pts
      if (nid.eq.0)  write(6,14) swaptime

      if (nsteps.ne.0) then
      if (nid.eq.0)  write(6,15) cpu_io_step_pts/cpu_p_t
      if (nid.eq.0)  write(6,16) comm_t/cpu_t*100.0
      if (nid.eq.0)  write(6,*) '====================================='
      endif
#endif

   33 format(' setup walltime      ::',  1p2e11.4,' sec ')
   44 format(' solve walltime      ::',  1p2e11.4,' sec ') !computation+io+overhead
   55 format(' total walltime      ::',  1p2e11.4,' sec ') !presolve+solver
    3 format(' setup cputime       ::',  1p2e11.4,' sec ')
    4 format(' solve cputime (w/io)::',  1p2e11.4,' sec ') !computation+io+overhead
    5 format(' total cputime       ::',  1p2e11.4,' sec ') !presolve+solver
    6 format(' total comp          ::',  1pe11.4,' sec ') !purely computation over timesteps
    7 format(' total communication ::',  1pe11.4,' sec ') !communication in flux over timesteps 
    8 format(' total io time       ::',  1pe11.4,' sec ') !io time over timesteps
   17 format(' total chk time      ::',  1pe11.4,' sec ') !io time over timesteps
    9 format(' computation   /step ::',  1pe11.4,' sec ') !computation per step                
   10 format(' communication /step ::',  1pe11.4,' sec ') !communication per step
   11 format(' io time       /step ::',  1pe11.4,' sec ') !io time per step          

   12 format(' comp time /step/pts ::',  1pe11.4,' sec ') !computation per step per point 
   13 format(' io   time /step/pts ::',  1pe11.4,' sec ') !io time per step per point
   14 format(' io swaptime         ::',  1pe11.4,' sec ') !communication vs. computation
   15 format(' io/comp  (/step/pts)::',  1pe11.4,' %   ') !communication time per step
   16 format(' communication/comp  ::',  1pe11.4,' %   ') !communication vs. computation

c     call printout     
      call exitt     !mpi_finalize is needed for darshan                                   

      return
      end
c-----------------------------------------------------------------------
c     processor initialization
      subroutine proc_init2
c----------------------------------------------------------------------
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
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

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
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

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
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'EMWAVE'
      INCLUDE 'ZPER'

      if (nid.eq.0) write(6,*) 'call para_init2'
      call setvar
      !call echopar2
      if (nid.eq.0) write(6,*) 'done: para_init2'

      return
      end
c-----------------------------------------------------------------------
c     geometry initialization
      subroutine geom_init2
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
c     INCLUDE 'WAKEPOT'
c     INCLUDE 'MOVWIN'

      igeom = 1
      call setup_topo    ! glonum,dsset(skpdat),setedge,etc

      if(nid.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat'

      call set_options   ! option logical value setting
      call setlog        ! Field initialization
      call gengeom(igeom)! temorary blocking for hmholtz 11/06/06 misun 

      if(nid.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat2'

      if(nid.eq.0) write(6,*) 'call geom_reset'
      call geom_reset(1) ! ic.f:  not using this for XM3
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat2'

      if(nid.eq.0) write(6,*) 'call set_unr   '
      call set_unr       ! normal, tangent vectors
      if(nid.eq.0) write(6,'(A,/)') ' done :: set_unr'

      !call vrdsmsh      ! validation mesh
      !call bcmask       ! bdry.f - boundary mapping

      if(nid.eq.0) write(6,*) 'call load_semhat_weighted'
      call load_semhat_weighted !semHAT
      if(nid.eq.0) write(6,'(A,/)') ' done :: load_semhat_weighted'

      !call setics       ! validation mesh

      return
      end
c----------------------------------------------------------------------
c     countings
      subroutine op_counts
c----------------------------------------------------------------------
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
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      real     dxmin
 
      if     (param(12).lt.0) then
          dt =abs(param(12))
          cfl=dt/dxmin
      elseif (param(12).gt.0) then
          cfl=param(12)
          dt =cfl*dxmin
      else
          if (nid.eq.0) write(6,*) 'set param(12) with a nonzero value'
          call exitt
      endif

      if (nid.eq.0) then
         write(6,*) '\n dxmin=',dxmin
         write(6,*) 'dt=',dt, ', CFL=',cfl
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_solver     
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
c     INCLUDE 'WAKEPOT'
c     INCLUDE 'MOVWIN'

      integer    lzdt
      parameter (LXR=1)
      parameter (LYR=1)
      parameter (LZR=1)
      parameter (LXYZR=LXR*LYR*LZR)
      parameter (LPSC9=1)
      parameter (LZDT=1)
      parameter (lctmp1n = LXYZR*LPSC9)

      common        zdt(LZDT)
      common /ccpu  /  cpu_t,cpu_dtime,cpu_chk

c...  setting logical options read from param
      call set_logics
      call errchk_logics

c...  initialization variables, materials, applications                        
c ... FIXME: update once the presetup for tet is done.                                              
      if (ifschrod) then
          call cem_schrod_temporary
      else
          call cem_maxwell_temporary
      endif

      call cem_apps_init ! presetup for each application         

c...  setting fintime, dt, cfl               
      call get_dxmin (dxmin)
      call get_dtsize(dxmin)

c ... initial conditions for non-restart/restart cases                                   
 100  istep   = 0 !   
      irststep= 0 ! restart step from previous run
      time    = 0.0
      call set_initial_fields

c ... timestep options and constants                                                           
      call errchk_timestep_options
      call rk_storage

c ... begin time iteration                                                                   
      call time_advancing

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_apps_init                
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
c     INCLUDE 'WAKEPOT'
c     INCLUDE 'MOVWIN'

      nmscale=1.0 
      if (IFNM)   nmscale=1.e-9
      if (ifpml ) call pml_setup
C     if (ifpml ) call cem_outer_pmlsrf !to be improved
C     if (ifbeam) call beam_setup

      return
      end
c---------------------------------------------------------------------
      subroutine errchk_logics
c---------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'

      imode=0

      if (if3d) then
          imode=3
          return
      endif 

      if (iftm) imode=2 
      if (ifte) imode=1 

      if (imode.eq.0) then
          if (nid.eq.0) write(6,*) 'param(4): TE or TM not defiend'
          call exitt
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine errchk_timestep_options
c---------------------------------------------------------------------
c ... error check for timestepping defined
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'

      if (ifrk22) ifrk=.true.
      if (ifrk33) ifrk=.true.
      if (ifrk44) ifrk=.true.
      if (ifrk45) ifrk=.true.
      if (ifsim1) ifsim=.true.
      if (ifsim2) ifsim=.true.
      if (ifsim3) ifsim=.true.
      if (ifsim4) ifsim=.true.
      if (ifsim5) ifsim=.true.


      if (ifeig.or.ifexp.or.ifrk.or.ifsim) then
         if (nid.eq.0) write(6,*) 'param(17): defined as ', param(17)
      else
         if (nid.eq.0) write(6,*) 'param(17): not defined', param(17)
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine time_advancing
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'SCHROD'
      common  /ccpu/ cpu_t,cpu_dtime,cpu_chk                 
      integer  ifirst, ilast
      real     fintime      

c ... counting and printing                               
      call print_steps(ifirst,ilast,fintime)

c ... start time stepping              
      cpu_t = 0.0
      cpu_chk= 0.0
      do 1000 istep= ifirst,ilast

         !if (ifmovw) call movwin_setup 

         cpu_dtime= dclock()

         call cem(igeom)

         time = time+dt

         cpu_dtime= dclock()-cpu_dtime
         cpu_t= cpu_t+cpu_dtime             

         cpu_chk= dclock()
         call userchk
         cpu_chk= dclock()-cpu_chk

         call cem_out

         if (lastep.eq.1)     goto 1001
         if (time.ge.fintime) goto 1001

         call time_remaining(cpu_t)

 1000 continue
 1001 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine print_steps(ifirst,ilast,fintime)                
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
c     INCLUDE 'WAKEPOT'
c     INCLUDE 'MOVWIN'
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
      subroutine time_remaining(cpu_t)
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      real cpu_t

      rtime= (cpu_t/istep)*(nsteps-istep)/60.0
      if ((nid.eq.0).and.(mod1(istep,iostep).eq.iostep)) then
               write(6,3) nsteps, rtime
      endif
    3 format(/,' max timestep=',i9,
     $        '; estimated remaining time (min)=',f10.2)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_initial_fields 
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'SCHROD' 
c     INCLUDE 'WAKEPOT'
c     INCLUDE 'MOVWIN'

      if (ifrestart) then
c ...    read/update initial fields from previous step: restarting case
         call restart_swap
         if (nid.eq.0) write(6,*) 'restart: done, time=',time 
         call userchk
         if (nid.eq.0) write(6,*) 'userchk: done, ifrestart=',ifrestart
      else

         if (ifsol) then
           if (IFSCHROD) then
             call usersol(time,sQr,sQi,sQi,sUr,sUi,sUi)
           else
             call usersol(time,sHN(1,1),sHN(1,2),sHN(1,3),
     $                         sEN(1,1),sEN(1,2),sEN(1,3))
           endif
           if (nid.eq.0) write(6,*) 'usersol: done, IFSOL=',ifsol
         endif

         if (IFSCHROD) then
             call userini(time,Qr,Qi,Qi,Ur,Ui,Ui)
             call schrod_grad (Qr(1,1),Qr(1,2),Qr(1,3),Ur)
             call schrod_grad (Qi(1,1),Qi(1,2),Qi(1,3),Ui)
         else
             call userini(time, HN(1,1), HN(1,2), HN(1,3),
     $                          EN(1,1), EN(1,2), EN(1,3))
         endif
         if (nid.eq.0) write(6,*) 'userini: done, IFSOL=',ifsol
    
         call userchk
         if (nid.eq.0) write(6,*) 'userchk: done, ifrestart=',ifrestart
         call cem_out

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine adjust_dt
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
c     INCLUDE 'WAKEPOT'
c     INCLUDE 'MOVWIN' 


c     if (param(14).gt.0) call dt_adjusted_nano


c     if (param(50).ne.0) then
c        n=nx1*ny1*nz1*nelt
c        dist=cspeed*dt
c        if     (param(50).eq.1) then
c          xmax=glmax(xm1,n)
c          xmin=glmin(xm1,n)
c          !nsteps=(xmax-xmin)/dist
c        elseif (param(50).eq.2) then
c          ymax=glmax(ym1,n)
c          ymin=glmin(ym1,n)
c          !nsteps=(ymax-ymin)/dist
c        elseif (param(50).eq.3) then
c          zmax=glmax(zm1,n)
c          zmin=glmin(zm1,n)
c          !nsteps=10*(zmax-zmin)/dist
c        else
c          if (nid.eq.0) write(6,*) 'param(50)=0,1,2?',param(50)
c          call exitt
c        endif
c        if (nid.eq.0) write(6,*) 'param(50)=',param(50)
c        if (nid.eq.0) write(6,*) 'updated nsteps for FT',nsteps
c     endif
c

c     if (ifwake) then
c        mstep = param(75) ! preset for wake
c        if (nid.eq.0) write(6,*) ' ifwake--true' 
c        if (nid.eq.0) write(6,*) 'after mstep=',mstep
c        call cem_wakepotential! preset for wake
c        dt = dzt/mstep        ! reset dt for wake
c
c        if (nid.eq.0) write(6,*)
c        if (nid.eq.0) write(6,*) 'check and adjust dt >>>>>>'
c        if (nid.eq.0) write(6,*) 'c*dt=',cspeed*dt, dxmin*cfl
c        if (cspeed*dt.gt.dxmin*cfl) then
c           if (nid.eq.0) then
c              write(6,*) 'small luniform or mstep',luniform,dt,dxmin
c              !call exitt
c           endif
c        endif
c        if (nid.eq.0) write(6,*) ' ifwake--true--end' 
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine chkmovwin
c-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
c     INCLUDE 'MOVWIN'

c     ifmovw=.false.
c     if (lmov.eq.1) then 
c         ifmovw=.true.
c         if (nid.eq.0) write(6,*) 'IFMOVW=',ifmovw
c     endif

      return
      end

c-----------------------------------------------------------------------
      subroutine printout
C-------------------------------------------------------------------
C
C     Transfer array dimensions to common
C
C-------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'

      if (nid.eq.0) then 

       if (if3d) then
         if (IFCENTRAL) write(6,335)
         if (IFCENTRAL) write(6,331)
         if (IFCENTRAL) write(6,335)
         if (IFUPWIND)  write(6,335)
         if (IFUPWIND)  write(6,332)
         if (IFUPWIND)  write(6,335)
         if ((.not.IFCENTRAL).and.(.not.IFUPWIND)) then
               write(6,336) 
               call exitt
         endif
       else
         if (IFCENTRAL) write(6,335)
         if (IFCENTRAL) write(6,333)
         if (IFCENTRAL) write(6,335)
         if (IFUPWIND)  write(6,335)
         if (IFUPWIND)  write(6,334)
         if (IFUPWIND)  write(6,335)
         if ((.not.IFCENTRAL).and.(.not.IFUPWIND)) then               
               write(6,336) 
               call exitt
         endif
       endif

      endif

 229  format ('   BEGIN with New DT    ',e15.7) 
 330  format ('         BEGIN          ') 
 331  format ('3D CENTRAL FLUX used !!!')
 332  format ('3D UPWIND  FLUX used !!!')
 333  format ('2D CENTRAL FLUX used !!!')
 334  format ('2D UPWIND  FLUX used !!!')
 335  format ('------------------------') 
 336  format ('ERR: NOT defined IFCENTRAL, IFUPWIND') 

c     call opcount(3)
c     call hmt_opcount_end()
c     call TIMEOUT
 
c----------------------------------------------------------------------
c     Time stepping loop end
c----------------------------------------------------------------------
 
c     CALL PREPOST (.true.,'   ')
c
      IF (NSTEPS.EQ.0) THEN
c             CALL PREPOST (.true.,'   ')
              CALL cem_out                          
      ENDIF

c     IF (instep.EQ.0) THEN
c        LASTEP=1
c        t0 = dclock()
c        CALL PREPOST (.true.,'   ')
c        t2 = t2 + (dclock()-t0)
c        NSTEPS=0
c     ENDIF
 

c     if (nid.eq.0) then
c        write(6,*) 'prepost time ::',t2,' seconds'
c     endif
 
c     call COMMENT
c     call DIAGNOS
c     call xxt_stats(xxth) 

      return
      end
c-----------------------------------------------------------------------
      subroutine initdim
C-------------------------------------------------------------------
C
C     Transfer array dimensions to common
C
C-------------------------------------------------------------------
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
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      COMMON /DOIT/ IFDOIT
      LOGICAL       IFDOIT
C
C     Set default logicals
C
      IFFLOW   = .FALSE.
      IFMVBD   = .FALSE.
      IFHEAT   = .FALSE.
      IFSPLIT  = .FALSE.
      IFDOIT   = .FALSE.
      IFMOAB   = .FALSE.
      IFTET    = .FALSE.
      IFHEX    = .FALSE.
      IFDG     = .FALSE.
      IFNEDELEC= .FALSE.
      IFCEM    = .TRUE. !default 
      if (nid.eq.0) write(6,*) 'IFCEM= ',IFCEM      

      call preprocessor_options(0)    !IFMOAB

      if (mesh.eq.0) IFHEX = .true.      
      if (mesh.eq.1) IFTET = .true.      

      if (nid.eq.0) write(6,*) '\n Mesh options:: '
      if (nid.eq.0) write(6,*) '  IFHEX= ',IFHEX,' (hexahedral )'
      if (nid.eq.0) write(6,*) '  IFTET= ',IFTET,' (tetrahedral)'      

      IFDG = .true. !default
      if (nedelec.eq.1) IFNEDELEC = .true.      
      if (nid.eq.0) write(6,*) '\n Available basis:: '
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
C
      subroutine comment
C---------------------------------------------------------------------
C
C     No need to comment !!
C
C---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      LOGICAL  IFCOUR
      SAVE     IFCOUR
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      REAL*8 ETIME0,ETIME1,ETIME2
      SAVE   ETIME0,ETIME1,ETIME2
      DATA   ETIME0,ETIME1,ETIME2 /0.0, 0.0, 0.0/
      REAL*8 DCLOCK
C
C     Only node zero makes comments.
      IF (NID.NE.0) RETURN
C
C
      IF (ETIME0.EQ.0.0) ETIME0=DCLOCK()
      ETIME1=ETIME2
      ETIME2=DCLOCK()
C
      IF (ISTEP.EQ.0) THEN
         IFCOUR  = .FALSE.
         DO 10 IFIELD=1,NFIELD
            IF (IFADVC(IFIELD)) IFCOUR = .TRUE.
 10      CONTINUE
         IF (IFWCNO) IFCOUR = .TRUE.
         WRITE (6,*) ' '
         WRITE (6,*) 'Initialization successfully completed'
         IF (TIME.NE.0.0) WRITE (6,*) 'Initial time is:',TIME
         WRITE (6,*) ' '
         WRITE (6,*) 'START OF SIMULATION'
         WRITE (6,*) ' '
      ELSEIF (ISTEP.GT.0 .AND. LASTEP.EQ.0 .AND. IFTRAN) THEN
         ETIME=ETIME2-ETIME1
         TTIME=ETIME2-ETIME0
         IF (     IFCOUR) 
     $      WRITE (6,100) ISTEP,TIME,DT,COURNO,TTIME,ETIME
         IF (.NOT.IFCOUR) WRITE (6,101) ISTEP,TIME,DT
      ELSEIF (LASTEP.EQ.1) THEN
         WRITE (6,*) ' '
         WRITE (6,*) 'Simulation successfully completed'
      ENDIF
 100  FORMAT('Step',I6,', t=',E14.7,', DT=',E14.7
     $,', C=',F7.3,2E11.4)
 101  FORMAT('Step',I6,', time=',E12.5,', DT=',E11.3)
C      call flush_io()
      RETURN
      END
C
      subroutine exit2
C     This is here because calling Sun-4 Fortran's EXIT causes a core dump (!)
      call exitt
      END
C
C------------------------------------------------------------------------
      subroutine setvar
C------------------------------------------------------------------------
C
C     Initialize variables
C
C------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      INCLUDE 'EMWAVE'
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
C
C     Enforce splitting/Uzawa according to the way the code was compiled
C
      IF (LX1.EQ.LX2) IFSPLIT=.TRUE.
      IF (LX1.NE.LX2) IFSPLIT=.FALSE.
C
C     Geometry on Mesh 3 or 1?
C
      IFGMSH3 = .TRUE.
      IF ( IFSTRS )           IFGMSH3 = .FALSE.
      IF (.NOT.IFFLOW)        IFGMSH3 = .FALSE.
C
      NFIELD = 1
      IF (IFHEAT) THEN
         NFIELD = 2 + NPSCAL
         NFLDT  = 1 + NPSCAL
      ENDIF

      MFIELD = 1
      IF (IFMVBD) MFIELD = 0
C
      DO 100 IFIELD=MFIELD,NFIELD
         IF (IFTMSH(IFIELD)) THEN
             NELFLD(IFIELD) = NELT
         ELSE
             NELFLD(IFIELD) = NELV
         ENDIF
 100  CONTINUE
C
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
      IF (.NOT.IFTRAN ) NSTEPS = 1
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
      DO 200 IFIELD=IFLD1,NFIELD
         IF (IFADVC(IFIELD)) IADV = 1
 200  CONTINUE
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
C
C------------------------------------------------------------------------
      subroutine echopar2
C------------------------------------------------------------------------
C
C     Echo the nonzero parameters from the readfile to the logfile
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'EMWAVE'
      INCLUDE 'ZPER'

      CHARACTER*132 STRING
      CHARACTER*1  STRING1(132)
      EQUIVALENCE (STRING,STRING1)
C
      IF (nid.ne.0) RETURN
C
      OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')
      REWIND(UNIT=9)
C
C
      READ(9,*,ERR=400)
      READ(9,*,ERR=400) VNEKTON
      NKTONV=VNEKTON
      VNEKMIN=2.5
      IF(VNEKTON.LT.VNEKMIN)THEN
         PRINT*,' Error: This NEKTON Solver Requires a .rea file'
         PRINT*,' from prenek version ',VNEKMIN,' or higher'
         PRINT*,' Please run the session through the preprocessor'
         PRINT*,' to bring the .rea file up to date.'
         call exitt
      ENDIF
      READ(9,*,ERR=400) NDIM
c     error check
      IF(NDIM.NE.LDIM)THEN
         WRITE(6,10) LDIM,NDIM
   10       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $              /,2X,'       for spatial dimension equal to',I2,'.'
     $              /,2X,'       The data file has dimension',I2,'.')
         CALL exitt
      ENDIF
C
      CALL BLANK(STRING,132)
      CALL CHCOPY(STRING,REAFLE,132)
      Ls=LTRUNC(STRING,132)
      READ(9,*,ERR=400) NPARAM
      WRITE(6,82) NPARAM,(STRING1(j),j=1,Ls)
C
      DO 20 I=1,NPARAM
         CALL BLANK(STRING,132)
         READ(9,80,ERR=400) STRING
         Ls=LTRUNC(STRING,132)
         !IF (PARAM(i).ne.0.0) 
         WRITE(6,81) I,(STRING1(j),j=1,Ls)
      ! write(6,*) 'do pppram--',param(86)
   20 CONTINUE
   80 FORMAT(A132)
   81 FORMAT(I4,3X,132A1)
   82 FORMAT(I4,3X,'Parameters from file:',132A1)
      CLOSE (UNIT=9)
      write(6,*) ' '

c     if (param(62).gt.0) then
c        if(nid.eq.0) write(6,*)
c    &      'enable byte swap for output'
c        call set_bytesw_write(1)
c     endif
c
      return
C
C     Error handling:
C
  400 CONTINUE
      WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA'
     $    ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
      CALL exitt
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA'
     $    ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
      CALL exitt
C
      RETURN
      END

C------------------------------------------------------------------------
      subroutine echopar
C------------------------------------------------------------------------
C
C     Echo the nonzero parameters from the readfile to the logfile
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      CHARACTER*80 STRING
      CHARACTER*1  STRING1(80)
      EQUIVALENCE (STRING,STRING1)
C
      IF (nid.ne.0) RETURN
C
      OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')
      REWIND(UNIT=9)
C
C
      READ(9,*,ERR=400)
      READ(9,*,ERR=400) VNEKTON
      NKTONV=VNEKTON
      VNEKMIN=2.5
      IF(VNEKTON.LT.VNEKMIN)THEN
         PRINT*,' Error: This NEKTON Solver Requires a .rea file'
         PRINT*,' from prenek version ',VNEKMIN,' or higher'
         PRINT*,' Please run the session through the preprocessor'
         PRINT*,' to bring the .rea file up to date.'
         call exitt
      ENDIF
      READ(9,*,ERR=400) NDIM
c     error check
      IF(NDIM.NE.LDIM)THEN
         WRITE(6,10) LDIM,NDIM
   10       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $              /,2X,'       for spatial dimension equal to',I2,'.'
     $              /,2X,'       The data file has dimension',I2,'.')
         CALL exitt
      ENDIF
C
      CALL BLANK(STRING,80)
      CALL CHCOPY(STRING,REAFLE,80)
      Ls=LTRUNC(STRING,80)
      READ(9,*,ERR=400) NPARAM
      WRITE(6,82) NPARAM,(STRING1(j),j=1,Ls)
C
      DO 20 I=1,NPARAM
         CALL BLANK(STRING,80)
         READ(9,80,ERR=400) STRING
         Ls=LTRUNC(STRING,80)
         IF (PARAM(i).ne.0.0) WRITE(6,81) I,(STRING1(j),j=1,Ls)
   20 CONTINUE
   80 FORMAT(A80) 
   81 FORMAT(I4,3X,80A1)
   82 FORMAT(I4,3X,'Parameters from file:',80A1)
      CLOSE (UNIT=9)
      RETURN
C
C     Error handling:
C
  400 CONTINUE
      WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA'
     $    ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
      CALL exitt
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA'
     $    ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
      CALL exitt
C
      RETURN
      END
C
C----------------------------------------------------------------------
      subroutine gengeom (igeom)
C------------------------------------------------------------------------
c
c     Generate geometry data
c

      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'GEOM'
      include 'WZ'
      include 'MASS'
 
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
 
      call genwz
      call gencoor (XM3,YM3,ZM3) !generate xm1
      call geom1   (XM3,YM3,ZM3) !jacm1,rxm1,sxm1,rzm1,..bm1,..,area,..

c     call geom2    !xm2,ym2,....
c     CALL UPDMSYS (1)
      call volume
c     call setinvm !cem doesn't use this inverse variables

      call setdef

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
   10 FORMAT(A132)

      GOTO 23
  24  ierr = 1

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

      call exitt

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
C
      subroutine timeout
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'CTIMER'
      dimension vdsum(lp),vgsum(lp),vbsol(lp),vusbc(lp)
      dimension vdsnd(lp),vdadd(lp),work(lp),vgop (lp)
      dimension vdsmn(lp),vdsmx(lp),vgsmn(lp),vgsmx(lp)
C
      real min_dsum, max_dsum, avg_dsum
      real min_gop,  max_gop,  avg_gop
      real min_crsl, max_crsl, avg_crsl
c
      real dhc, dwork
C
      tstop=dclock()
      ttotal=tstop-etimes
      tttstp=tstop-etims0
c
      if (nid.eq.0) then
         write(6,44) taaaa,tbbbb,tcccc,tdddd,teeee
   44    format(1p5e12.4,' HSMG time')
      endif
c
c      call gdhigh(tusbc,1,wwork)
c
      min_dsum = tdsum
      call gop(min_dsum,wwork,'m  ',1)
c      call gdlow (min_dsum,1,wwork)
      max_dsum = tdsum
      call gop(max_dsum,wwork,'M  ',1)
c      call gdhigh(max_dsum,1,wwork)
      avg_dsum = tdsum
      call gop(avg_dsum,wwork,'+  ',1)
c      call gdsum (avg_dsum,1,wwork)
      avg_dsum = avg_dsum/np
c
      min_gop = tgop
      call gop(min_gop,wwork,'m  ',1)
c      call gdlow (min_gop,1,wwork)
      max_gop = tgop
      call gop(max_gop,wwork,'M  ',1)
c      call gdhigh(max_gop,1,wwork)
      avg_gop = tgop
      call gop(avg_gop,wwork,'+  ',1)
c      call gdsum (avg_gop,1,wwork)
      avg_gop = avg_gop/np
c
      min_crsl = tcrsl
      call gop(min_crsl,wwork,'m  ',1)
c      call gdlow (min_crsl,1,wwork)
      max_crsl = tcrsl
      call gop(max_crsl,wwork,'M  ',1)
c      call gdhigh(max_crsl,1,wwork)
      avg_crsl = tcrsl
      call gop(avg_crsl,wwork,'+  ',1)
c      call gdsum (avg_crsl,1,wwork)
      avg_crsl = avg_crsl/np
c
      tttstp = tttstp + 1e-7
      if (nid.eq.0) then
         write(6,*) 'total time',ttotal,tttstp
         ttotal=tttstp
         pcopy=tcopy/tttstp
         write(6,*) 'copy time',ncopy,tcopy,pcopy
         pmxmf=tmxmf/tttstp
         write(6,*) 'mxmf time',nmxmf,tmxmf,pmxmf
         pmxms=tmxms/tttstp
         write(6,*) 'mxms time',nmxms,tmxms,pmxms
         pinv3=tinv3/tttstp
         write(6,*) 'inv3 time',ninv3,tinv3,pinv3
         pinvc=tinvc/tttstp
         write(6,*) 'invc time',ninvc,tinvc,pinvc
         pmltd=tmltd/tttstp
         write(6,*) 'mltd time',nmltd,tmltd,pmltd
         pcdtp=tcdtp/tttstp
         write(6,*) 'cdtp time',ncdtp,tcdtp,pcdtp
         ppres=tpres/tttstp
         write(6,*) 'pres time',npres,tpres,ppres
         phmhz=thmhz/tttstp
         write(6,*) 'hmhz time',nhmhz,thmhz,phmhz
         pusbc=tusbc/tttstp
         write(6,*) 'usbc time',nusbc,tusbc,pusbc
         paxhm=taxhm/tttstp
         write(6,*) 'axhm time',naxhm,taxhm,paxhm
c
         pgop =tgop /tttstp
         write(6,*) 'gop  time',ngop ,tgop ,pgop 
         write(6,*) 'gop  min ',min_gop 
         write(6,*) 'gop  max ',max_gop 
         write(6,*) 'gop  avg ',avg_gop 
c
         pdsum=tdsum/tttstp
         write(6,*) 'dsum time',ndsum,tdsum,pdsum
         write(6,*) 'dsum min ',min_dsum
         write(6,*) 'dsum max ',max_dsum
         write(6,*) 'dsum avg ',avg_dsum
c
         pgsum=tgsum/tttstp
         write(6,*) 'gsum time',ngsum,tgsum,pgsum
         pdsnd=tdsnd/tttstp
         write(6,*) 'dsnd time',ndsnd,tdsnd,pdsnd
         pdadd=tdadd/tttstp
         write(6,*) 'dadd time',ndadd,tdadd,pdadd
         pdsmx=tdsmx/tttstp
         write(6,*) 'dsmx time',ndsmx,tdsmx,pdsmx
         pdsmn=tdsmn/tttstp
         write(6,*) 'dsmn time',ndsmn,tdsmn,pdsmn
         pgsmx=tgsmx/tttstp
         write(6,*) 'gsmx time',ngsmx,tgsmx,pgsmx
         pgsmn=tgsmn/tttstp
         write(6,*) 'gsmn time',ngsmn,tgsmn,pgsmn
         pslvb=tslvb/tttstp
         write(6,*) 'slvb time',nslvb,tslvb,pslvb
         pddsl=tddsl/tttstp
         write(6,*) 'ddsl time',nddsl,tddsl,pddsl
c
         pcrsl=tcrsl/tttstp
         write(6,*) 'crsl time',ncrsl,tcrsl,pcrsl
         write(6,*) 'crsl min ',min_crsl
         write(6,*) 'crsl max ',max_crsl
         write(6,*) 'crsl avg ',avg_crsl
c
         psolv=tsolv/tttstp
         write(6,*) 'solv time',nsolv,tsolv,psolv
         psett=tsett/tttstp
         write(6,*) 'sett time',nsett,tsett,psett
         pprep=tprep/tttstp
         write(6,*) 'prep time',nprep,tprep,pprep
         pbsol=tbsol/tttstp
         write(6,*) 'bsol time',nbsol,tbsol,pbsol
         pbso2=tbso2/tttstp
         write(6,*) 'bso2 time',nbso2,tbso2,pbso2
      endif
      if (np.gt.0) then
c        call rzero(vbsol,np)
c        vbsol(node)=tbsol
c        call gop(vbsol,work,'+  ',np)
c
         call rzero(vusbc,np)
         vusbc(node)=tusbc
         call gop(vusbc,work,'+  ',np)
         write(6,*) nid,' nusbc',nusbc,pusbc
c
         call rzero(vdsnd,np)
         vdsnd(node)=tdsnd
         call gop(vdsnd,work,'+  ',np)
c
         call rzero(vdadd,np)
         vdadd(node)=tdadd
         call gop(vdadd,work,'+  ',np)
c
         call rzero(vgsum,np)
         vgsum(node)=tgsum
         call gop(vgsum,work,'+  ',np)
c
         call rzero(vdsum,np)
         vdsum(node)=tdsum
         call gop(vdsum,work,'+  ',np)
c
         call rzero(vgop ,np)
         vgop (node)=tgop 
         call gop(vgop ,work,'+  ',np)
c
         call rzero(vdsmx,np)
         vdsmx(node)=tdsmx
         call gop(vdsmx,work,'+  ',np)
c
         call rzero(vdsmn,np)
         vdsmn(node)=tdsmn
         call gop(vdsmn,work,'+  ',np)
c
         call rzero(vgsmx,np)
         vgsmx(node)=tgsmx
         call gop(vgsmx,work,'+  ',np)
c
         call rzero(vgsmn,np)
         vgsmn(node)=tgsmn
         call gop(vgsmn,work,'+  ',np)
c
         if (nid.eq.0) then
c
            write(6,202) np,nelgv,tttstp
            write(6,203) ndsum,ndsum,nbsol
            write(6,*) 'qqq ip tdsum tdsnd tdadd tgsum tgop tusbs'
            do 100 ip=1,np
               write(6,204) ip,vdsum(ip),vdsnd(ip),vdadd(ip)
     $                         ,vgsum(ip),vgop(ip),vusbc(ip)
  100       continue
            write(6,*) 'qqq ip dsavg tdsmn tdsmx tgsmn tgsmx'
            do 200 ip=1,np
               dsavg = vdsum(ip)/(ndsum)
               write(6,204) ip,dsavg,vdsmn(ip),vdsmx(ip)
     $                        ,vgsmn(ip),vgsmx(ip)
  200       continue
  202       format('qqq  np,nel,tttstp:',2i8,f12.5)
  203       format('qqq  num procs',/,' dot,dsum,bsol:',3i8)
  204       format('qqq', i4,6f12.5)
C
C
            write(6,*) 'qqq ip tdsum tdsnd tdadd tgsum tgop'
            do 110 ip=1,np
               rdsum=vdsum(ip)/tttstp
               rdsnd=vdsnd(ip)/tttstp
               rdadd=vdadd(ip)/tttstp
               rgsum=vgsum(ip)/tttstp
               rgop =vgop (ip)/tttstp
               rusbc=vusbc(ip)/tttstp
               write(6,204) ip,rdsum,rdsnd,rdadd,rgsum,rgop,rusbc
  110       continue
         endif
      endif
C
      return
      end


C------------------------------------------------------------------------
      subroutine timeout_OLD
C------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'CTIMER'
      dimension vdsum(lp),vgsum(lp),vbsol(lp),vusbc(lp)
      dimension vdsnd(lp),vdadd(lp),work(lp),vgop (lp)
      dimension vdsmn(lp),vdsmx(lp),vgsmn(lp),vgsmx(lp)
C
      tstop=dclock()
      ttotal=tstop-etimes
      tttstp=tstop-etims0
c
      if (nid.eq.0) then
         write(6,*) 'total time',ttotal,tttstp
         ttotal=tttstp
         pcopy=tcopy/tttstp
         write(6,*) 'copy time',ncopy,tcopy,pcopy
         pmxmf=tmxmf/tttstp
         write(6,*) 'mxmf time',nmxmf,tmxmf,pmxmf
         pmxms=tmxms/tttstp
         write(6,*) 'mxms time',nmxms,tmxms,pmxms
         pinv3=tinv3/tttstp
         write(6,*) 'inv3 time',ninv3,tinv3,pinv3
         pinvc=tinvc/tttstp
         write(6,*) 'invc time',ninvc,tinvc,pinvc
         pmltd=tmltd/tttstp
         write(6,*) 'mltd time',nmltd,tmltd,pmltd
         pcdtp=tcdtp/tttstp
         write(6,*) 'cdtp time',ncdtp,tcdtp,pcdtp
         ppres=tpres/tttstp
         write(6,*) 'pres time',npres,tpres,ppres
         phmhz=thmhz/tttstp
         write(6,*) 'hmhz time',nhmhz,thmhz,phmhz
         pusbc=tusbc/tttstp
         write(6,*) 'usbc time',nusbc,tusbc,pusbc
         paxhm=taxhm/tttstp
         write(6,*) 'axhm time',naxhm,taxhm,paxhm
         pgop =tgop /tttstp
         write(6,*) 'gop  time',ngop ,tgop ,pgop 
         pdsum=tdsum/tttstp
         write(6,*) 'dsum time',ndsum,tdsum,pdsum
         pgsum=tgsum/tttstp
         write(6,*) 'gsum time',ngsum,tgsum,pgsum
         pdsnd=tdsnd/tttstp
         write(6,*) 'dsnd time',ndsnd,tdsnd,pdsnd
         pdadd=tdadd/tttstp
         write(6,*) 'dadd time',ndadd,tdadd,pdadd
         pdsmx=tdsmx/tttstp
         write(6,*) 'dsmx time',ndsmx,tdsmx,pdsmx
         pdsmn=tdsmn/tttstp
         write(6,*) 'dsmn time',ndsmn,tdsmn,pdsmn
         pgsmx=tgsmx/tttstp
         write(6,*) 'gsmx time',ngsmx,tgsmx,pgsmx
         pgsmn=tgsmn/tttstp
         write(6,*) 'gsmn time',ngsmn,tgsmn,pgsmn
         pslvb=tslvb/tttstp
         write(6,*) 'slvb time',nslvb,tslvb,pslvb
         pddsl=tddsl/tttstp
         write(6,*) 'ddsl time',nddsl,tddsl,pddsl
         pcrsl=tcrsl/tttstp
         write(6,*) 'crsl time',ncrsl,tcrsl,pcrsl
         psolv=tsolv/tttstp
         write(6,*) 'solv time',nsolv,tsolv,psolv
         psett=tsett/tttstp
         write(6,*) 'sett time',nsett,tsett,psett
         pprep=tprep/tttstp
         write(6,*) 'prep time',nprep,tprep,pprep
         pbsol=tbsol/tttstp
         write(6,*) 'bsol time',nbsol,tbsol,pbsol
         pbso2=tbso2/tttstp
         write(6,*) 'bso2 time',nbso2,tbso2,pbso2
      endif
      if (np.gt.0) then
c        call rzero(vbsol,np)
c        vbsol(node)=tbsol
c        call gop(vbsol,work,'+  ',np)
c
         call rzero(vusbc,np)
         vusbc(node)=tusbc
         call gop(vusbc,work,'+  ',np)
         write(6,*) nid,' nusbc',nusbc,pusbc
c
         call rzero(vdsnd,np)
         vdsnd(node)=tdsnd
         call gop(vdsnd,work,'+  ',np)
c
         call rzero(vdadd,np)
         vdadd(node)=tdadd
         call gop(vdadd,work,'+  ',np)
c
         call rzero(vgsum,np)
         vgsum(node)=tgsum
         call gop(vgsum,work,'+  ',np)
c
         call rzero(vdsum,np)
         vdsum(node)=tdsum
         call gop(vdsum,work,'+  ',np)
c
         call rzero(vgop ,np)
         vgop (node)=tgop 
         call gop(vgop ,work,'+  ',np)
c
         call rzero(vdsmx,np)
         vdsmx(node)=tdsmx
         call gop(vdsmx,work,'+  ',np)
c
         call rzero(vdsmn,np)
         vdsmn(node)=tdsmn
         call gop(vdsmn,work,'+  ',np)
c
         call rzero(vgsmx,np)
         vgsmx(node)=tgsmx
         call gop(vgsmx,work,'+  ',np)
c
         call rzero(vgsmn,np)
         vgsmn(node)=tgsmn
         call gop(vgsmn,work,'+  ',np)
c
         if (nid.eq.0) then
c
            write(6,202) np,nelgv,tttstp
            write(6,203) ndsum,ndsum,nbsol
            write(6,*) 'qqq ip tdsum tdsnd tdadd tgsum tgop tusbs'
            do 100 ip=1,np
               write(6,204) ip,vdsum(ip),vdsnd(ip),vdadd(ip)
     $                         ,vgsum(ip),vgop(ip),vusbc(ip)
  100       continue
            write(6,*) 'qqq ip dsavg tdsmn tdsmx tgsmn tgsmx'
            do 200 ip=1,np
               dsavg = vdsum(ip)/(ndsum)
               write(6,204) ip,dsavg,vdsmn(ip),vdsmx(ip)
     $                        ,vgsmn(ip),vgsmx(ip)
  200       continue
  202       format('qqq  np,nel,tttstp:',2i8,f12.5)
  203       format('qqq  num procs',/,' dot,dsum,bsol:',3i8)
  204       format('qqq', i4,6f12.5)
C
C
            write(6,*) 'qqq ip tdsum tdsnd tdadd tgsum tgop'
            do 110 ip=1,np
               rdsum=vdsum(ip)/tttstp
               rdsnd=vdsnd(ip)/tttstp
               rdadd=vdadd(ip)/tttstp
               rgsum=vgsum(ip)/tttstp
               rgop =vgop (ip)/tttstp
               rusbc=vusbc(ip)/tttstp
               write(6,204) ip,rdsum,rdsnd,rdadd,rgsum,rgop,rusbc
  110       continue
         endif
      endif
C
      return
      end
C------------------------------------------------------------------------
      subroutine opcount(ICALL)
C------------------------------------------------------------------------
C
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
         do 100 i=1,maxrts
            ncall(i) = 0
            dct(i)   = 0.0
  100    continue
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
            do 200 i=1,nrout
               write(6,201) nid,rname(i),rct(i),ncall(i)
  200       continue
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
      do 10 i=1,ntot1
         if (vmult(i,1,1,1).lt.0.5) work(1)=work(1)+vmult(i,1,1,1)
   10 continue
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
