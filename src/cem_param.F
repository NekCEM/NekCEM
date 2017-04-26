c-----------------------------------------------------------------------
      subroutine set_logics
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      IFDG     =.true.    ! spectral element DG method: default
      IFSE     =.false.   ! spectral element method   :
      IFEXCT   =.false.   ! drift-diffusion solver
      IFVARDT  =.false.   ! variable timestep size
      IFFILTER =.false.   ! filter
      IFIC     =.false.   ! initial condition
      IFRESTART=.false.   ! restart output file
      IFSCAT   =.false.   ! scattered field mode
      IFSFTF   =.false.   ! scattered-total field mode
      IFGRAPH  =.false.   ! graphene
      IFNAPY   =.false.   ! napoly contour for wake field calculation
      IFPBG    =.false.   ! photonic band gap application
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
      IFINHOM  =.false.   ! inhomogeneous boundary

      call ifoption (IFRESTART)
      call ifoption (IFPP )

c     if (param( 14).lt.0)  IFPOISSON=.true. : use it for nanoperiodicity with (-)
      if (param(  4).eq.1)  IFTE     =.true.
      if (param(  4).eq.2)  IFTM     =.true.
      if (param( 60).eq.1)  IFGRAPH  =.true.

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

      if (abs(param( 21)).eq. 1) IFGMRES  =.true.
      if (abs(param( 21)).eq. 2) IFCG     =.true.
      if (abs(param( 21)).eq. 3) IFSEMG   =.true.
      if (param( 23).eq. 1) IFFDM    =.true.
      if (param( 25).eq. 1) IFGFDMDD =.true.

      if (ifrk22.or.ifrk33.or.ifrk44.or.ifrk45) ifrk=.true.
      if (ifsim1.or.ifsim2.or.ifsim3.or.ifsim4.or.ifsim5) ifsim=.true.

      if (nid.eq.0) then

      write(6,*) 'IFGRAPH     ',IFGRAPH
      write(6,*) 'IFMOVW      ',IFMOVW
      write(6,*) 'IFFILTER    ',IFFILTER
      write(6,*) 'IFNAPY      ',IFNAPY
      write(6,*) 'IFVARDT     ',IFVARDT
      write(6,*) 'IFPBG       ',IFPBG
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
      write(6,*) 'IFPEC       ',IFPEC
      write(6,*) 'IFROB       ',IFROB
      write(6,*) 'IFNEU       ',IFNEU
      write(6,*) 'IFINHOM     ',IFINHOM
      write(6,*) 'IFPML       ',IFPML
      write(6,*) 'IFCUR       ',IFCUR
      write(6,*) 'IFDTN       ',IFDTN
      write(6,*) 'IFRESTART   ',IFRESTART

      if (IFEXP) write(6,*) 'Krylov_Dim =',larnol

      endif

      xyzscale = 1.0

      return
      end
c-----------------------------------------------------------------------
