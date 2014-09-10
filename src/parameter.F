
c-----------------------------------------------------------------------
      subroutine set_logics
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      include 'NEKUSE'

      IFDG     =.true.    ! spectral element DG method: default 
      IFSE     =.false.   ! spectral element method   :  
      IFMAXWELL=.false.   ! maxwell solver
      IFSCHROD =.false.   ! schrodinger sholver
      IFDRIFT  =.false.   ! drift-diffusion solver
      IFEXCT   =.false.   ! drift-diffusion solver
      IFVARDT  =.false.   ! variable timestep size   
      IFFILTER =.false.   ! filter
      IFDIELEC =.false.   ! dielectric material
      IFIC     =.false.   ! initial condition
      IFRESTART=.false.   ! restart output file
      IFDTCTL  =.false.   ! dt control
      IFSCAT   =.false.   ! scattered field mode
      IFSFTF   =.false.   ! scattered-total field mode
      IFDRUDE  =.false.   ! drude model
      IFLOREN  =.false.   ! lorentz model
      IFHYDRO  =.false.   ! hydrodynamic model
      IFBEAM   =.false.   ! beam simulation
      IFWAKE   =.false.   ! wakefield calculation
      IFNAPY   =.false.   ! napoly contour for wake field calculation
      IFWAVEGUIDE=.false. ! waveguide problems
      IFPBG    =.false.   ! photonic band gap application
      IFWAKET  =.false.   ! wake field calculation 
      IFVTK    =.false.   ! vtk format output files

      IFDEALIAS=.false.   ! dealiasing
      IFUPWIND =.false.   ! upwind flux for DG 
      IFCENTRAL=.false.   ! central flux for DG

      IFSOL    =.false.   ! the cases with exact solutions 
      IFSRC    =.false.   ! source term
      IFPOISSON=.false.   ! poisson solver
      IFGMRES  =.false.   ! GMRES iterative solver
      IFCG     =.false.   ! CG iterative solver        
      IFSEMG   =.false.   ! GMRES multigird            
      IFFDM    =.false.   ! diagonal preconditioner             
      IFPRECON =.false.   ! fast diagonalization method for preconditioning
      IFGFDMDD =.false.   ! fast diagonalization method
      IFDIRECT =.false.

      IFRK45   =.false.   ! Runge-Kutta (RK) 5-state 4-th order
      IFRK44   =.false.   ! RK 4-stage 4th order
      IFRK33   =.false.   ! RK 3-stage 3rd order
      IFRK22   =.false.   ! RK 2-stage 2nd order
      IFRK     =.false.   ! All Runge-Kutta methods
      IFBDF    =.false.   ! All BDF methods
      IFBDF1   =.false.   ! BDF1 time stepping
      IFBDF2   =.false.   ! BDF2 time stepping
      IFSIM5   =.false.   ! 5th-order simplitic time stepping
      IFSIM4   =.false.   ! 4th-order simplitic time stepping
      IFSIM3   =.false.   ! 3rd-order simplitic time stepping
      IFSIM2   =.false.   ! 2nd-order simplitic time stepping
      IFSIM1   =.false.   ! 1st-order simplitic time stepping
      IFSIM    =.false.   ! all simplicit time stepping
      IFEXP    =.false.   ! exponential time stepping
      IFEIG    =.false.   ! eigenvalue solver

      IFTE     =.false.   ! transeverse electric mode
      IFTM     =.false.   ! transverse magnetic mode
      IFNM     =.false.   ! nanoscale 
      IFINHOM  =.false.   ! inhomogeneous boundary 

      call ifoption (IFRESTART)
      call ifoption (IFPP )

c     if (param( 14).lt.0)  IFPOISSON=.true. : use it for nanoperiodicity with (-)
      if (param(  4).eq.1)  IFTE     =.true.
      if (param(  4).eq.2)  IFTM     =.true.
      if (param(  5).eq.1)  IFDRUDE  =.true.
      if (param(  5).eq.2)  IFLOREN  =.true.
      if (param(  5).eq.3)  IFHYDRO  =.true.
      if (param(  5).eq.4)  IFBEAM   =.true.
      if (param(  5).eq.5 ) IFBEAM   =.true. !to be updated misun 7/9/2012

      if (param( 5).lt.10) IFMAXWELL=.true. !to be updated misun 7/9/2012
      if ((param(5).ge.10).and.(param(5).lt.20)) IFSCHROD =.true. 
      if ((param(5).ge.20).and.(param(5).lt.40)) IFDRIFT  =.true.
      if ((param(5).eq.21).or. (param(5).eq.31)) IFEXCT   =.true.
      if ((param(5).ge.30).and.(param(5).lt.40)) IFSE     =.true.

      if (param(  6).eq. 1) IFSRC    =.true.
      if (param(  7).eq. 1) IFINHOM  =.true.

c     if (param( 14).eq.0)  IFPOISSON=.true.
c     if (param( 14).eq.1)  IFSCAT   =.true.
c     if (param( 14).eq.3)  IFSFTF   =.true.

      if (param( 16).eq. 1) IFSOL    =.true.
      if (param( 17).eq.10) IFEIG    =.true.
      if (param( 17).eq. 0) IFRK45   =.true.
      if (param( 17).eq. 1) IFEXP    =.true.
      if (param( 17).eq.45) IFRK45   =.true.
      if (param( 17).eq.44) IFRK44   =.true.
      if (param( 17).eq.33) IFRK33   =.true.
      if (param( 17).eq.22) IFRK22   =.true.
      if (param( 17).eq.555)IFSIM5   =.true.
      if (param( 17).eq.444)IFSIM4   =.true.
      if (param( 17).eq.333)IFSIM3   =.true.
      if (param( 17).eq.222)IFSIM2   =.true.
      if (param( 17).eq.111)IFSIM1   =.true.
      if (param( 17).eq.-1) IFBDF1   =.true.
      if (param( 17).eq.-2) IFBDF2   =.true.
      if (param( 17).lt. 0) IFBDF    =.true.

      if (param( 18).eq. 1) IFFILTER =.true.
      if (param( 18).eq.-1) IFDEALIAS=.true.
      if (param( 19).eq. 0) IFUPWIND =.true.
      if (param( 19).eq. 1) IFCENTRAL=.true.
      if (param( 84).ne. 0) IFRESTART=.true.
      if (param(103).eq. 1) IFNM     =.true.

      if (param( 21).eq. 0) IFGMRES  =.true.
      if (param( 21).eq. 1) IFCG     =.true.
      if (param( 23).eq. 1) IFFDM    =.true. 
      if (param( 24).eq. 1) IFSEMG   =.true.
      if (param( 25).eq. 1) IFGFDMDD =.true.

      if (ifrk22.or.ifrk33.or.ifrk44.or.ifrk45) ifrk=.true.
      if (ifsim1.or.ifsim2.or.ifsim3.or.ifsim4.or.ifsim5) ifsim=.true.

      if (nid.eq.0) then

      write(6,*) 'IFMAXWELL   ',IFMAXWELL
      write(6,*) 'IFSCHROD    ',IFSCHROD
      write(6,*) 'IFDRIFT     ',IFDRIFT  
      write(6,*) 'IFDRUDE     ',IFDRUDE
      write(6,*) 'IFLOREN     ',IFLOREN
      write(6,*) 'IFBEAM      ',IFBEAM
      write(6,*) 'IFWAKE      ',IFWAKE
      write(6,*) 'IFMOVW      ',IFMOVW
      write(6,*) 'IFFILTER    ',IFFILTER
      write(6,*) 'IFNAPY      ',IFNAPY
      write(6,*) 'IFWAKET     ',IFWAKET
      write(6,*) 'IFVARDT     ',IFVARDT
      write(6,*) 'IFPBG       ',IFPBG
      write(6,*) 'IFDIELEC    ',IFDIELEC
      write(6,*) 'IFPOISSON   ',IFPOISSON
      write(6,*) 'IFGMRES     ',IFGMRES   
      write(6,*) 'IFCG        ',IFCG      
      write(6,*) 'IFSEMG      ',IFSEMG    
      write(6,*) 'IFFDM       ',IFFDM     
      write(6,*) 'IFIC        ',IFIC
      write(6,*) 'IFSOL       ',IFSOL
      write(6,*) 'IFSRC       ',IFSRC
      write(6,*) 'IFSCAT      ',IFSCAT
      write(6,*) 'IFSFTF      ',IFSFTF
      write(6,*) 'IFCENTRAL   ',IFCENTRAL
      write(6,*) 'IFUPWIND    ',IFUPWIND
      write(6,*) 'IFDEALIAS   ',IFDEALIAS
      write(6,*) 'IFTM        ',IFTM
      write(6,*) 'IFTE        ',IFTE
      write(6,*) 'IFEIG       ',IFEIG
      write(6,*) 'IFEXP       ',IFEXP
      write(6,*) 'IFRK45      ',IFRK45
      write(6,*) 'IFRK44      ',IFRK44
      write(6,*) 'IFRK33      ',IFRK33
      write(6,*) 'IFRK22      ',IFRK22
      write(6,*) 'IFSIM5      ',IFSIM5 
      write(6,*) 'IFSIM4      ',IFSIM4 
      write(6,*) 'IFSIM3      ',IFSIM3 
      write(6,*) 'IFSIM2      ',IFSIM2 
      write(6,*) 'IFSIM1      ',IFSIM1 
      write(6,*) 'IFDTCTL     ',IFDTCTL
      write(6,*) 'IFPEC       ',IFPEC
      write(6,*) 'IFINHOM     ',IFINHOM
      write(6,*) 'IFPML       ',IFPML
      write(6,*) 'IFDTN       ',IFDTN
      write(6,*) 'IFNM        ',IFNM
      write(6,*) 'IFRESTART   ',IFRESTART

      if (IFEXP) write(6,*) 'Krylov_Dim =',larnol

      endif

      if (IFNM) then

          cspeed  = param(1)
          permit_0= param(2)
          permea_0= param(3)

          xyzscale= 1.0
          if (IFDRUDE) xyzscale=1.e+9 ! for visual tool
          if (IFBEAM ) xyzscale=1.0
          if (IFWAKE ) xyzscale=1.0

          nmscale = 1.0
          if (IFDRUDE) nmscale =1.e-9 ! for visual tool
          if (IFBEAM ) nmscale =1.0
          if (IFWAKE ) nmscale =1.0

      else
          cspeed  = 1.0
          permit_0= 1.0
          permea_0= 1.0

          xyzscale= 1.0
          nmscale = 1.0
      endif

      if (nid.eq.0) then

      write(6,*)
      write(6,*) 'cspeed   =',cspeed
      write(6,*) 'permit_0 =',permit_0
      write(6,*) 'permea_0 =',permea_0
      write(6,*) 'nmscale  =',nmscale
      write(6,*)

      endif

      if (if3d) then
        hxactive = .true.
        hyactive = .true.
        hzactive = .true.
        exactive = .true.
        eyactive = .true.
        ezactive = .true.

        fldidx(1) = 1
        fldidx(2) = 2
        fldidx(3) = 3
        fldidx(4) = 4
        fldidx(5) = 5
        fldidx(6) = 6
        fldcnt = 6
      elseif (iftm) then
        hxactive = .true.
        hyactive = .true.
        hzactive = .false.
        exactive = .false.
        eyactive = .false.
        ezactive = .true.

        fldidx(1) = 1 !hx
        fldidx(2) = 2
        fldidx(3) = 0 !hz
        fldidx(4) = 0 !ex
        fldidx(5) = 0
        fldidx(6) = 3 !ez
        fldcnt = 3
      elseif (ifte) then
        hxactive = .false.
        hyactive = .false.
        hzactive = .true.
        exactive = .true.
        eyactive = .true.
        ezactive = .false.

        fldidx(1) = 0 !hx
        fldidx(2) = 0
        fldidx(3) = 3 !hz
        fldidx(4) = 1 !ex
        fldidx(5) = 2
        fldidx(6) = 0 !ez
        fldcnt = 3
      endif

      return
      end  
