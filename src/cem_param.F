c-----------------------------------------------------------------------
      subroutine set_logics
c-----------------------------------------------------------------------
c     Read from the rea file and do some checks that the parameters are
c     sane.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      include 'PML'

      ifdg = .true.             ! spectral element DG method: default
      ifse = .false.            ! spectral element method
      ifexct = .false.          ! drift-diffusion solver
      ifsteric = .false.        ! drift-diffusion w/ steric term
      iffilter = .false.        ! filter
      ific = .false.            ! initial condition
      ifrestart = .false.       ! restart output file

      ifdealias = .false.       ! dealiasing
      ifupwind = .false.        ! upwind flux for DG
      ifcentral = .false.       ! central flux for DG

      ifsol = .false.           ! the cases with exact solutions
      ifsrc = .false.           ! source term
      ifpoisson = .false.       ! poisson solver
      ifgmres = .false.         ! GMRES iterative solver
      ifcg = .false.            ! CG iterative solver
      ifsemg = .false.          ! GMRES multigird
      iffdm = .false.           ! diagonal preconditioner
      ifprecon = .false.        ! fast diagonalization method for preconditioning
      ifgfdmdd = .false.        ! fast diagonalization method

      ifrk45 = .false.          ! Runge-Kutta (RK) 5-state 4-th order
      ifrk22 = .false.          ! RK 2-stage 2nd order
      ifbdf = .false.           ! All BDF methods
      ifbdf1 = .false.          ! BDF1 time stepping
      ifbdf1r = .false.         ! BDF1 time stepping with Richardson extrapolation
      ifbdf2 = .false.          ! BDF2 time stepping
      ifexp = .false.           ! exponential time stepping
      ifeig = .false.           ! eigenvalue solver
c TODO: default ifte = true
      ifte = .true.            ! transeverse electric mode
      iftm = .false.            ! transverse magnetic mode
      ifinhom = .false.         ! inhomogeneous boundary

      call ifoption(ifrestart)
      call ifoption(ifpp)

      if (param(4).eq.1) ifte = .true.
      if (param(4).eq.2) iftm = .true.

      if ((param(5).eq.21).or.(param(5).eq.31)) ifexct = .true.
      if ((param(5).ge.30).and.(param(5).lt.40)) ifse = .true.

      if (param(6).eq.1) ifsrc = .true.
      if (param(7).eq.1) ifinhom = .true.

      if (param(16).eq.1) ifsol = .true.
      if (param(17).eq.10) ifeig = .true.
      if (param(17).eq.0) ifrk45 = .true.
      if (param(17).eq.1) ifexp = .true.
      if (param(17).eq.45) ifrk45 = .true.
      if (param(17).eq.22) ifrk22 = .true.
      if (param(17).eq.-1) ifbdf1 = .true.
      if (param(17).eq.-11) ifbdf1r = .true.
      if (param(17).eq.-2) ifbdf2 = .true.
      if (param(17).lt.0) ifbdf = .true.

      if (param(18).eq.1) iffilter = .true.
      if (param(18).eq.-1) ifdealias = .true.
      if (param(19).eq.0) ifupwind = .true.
      if (param(19).eq.1) ifcentral = .true.
      if (param(84).ne.0) ifrestart = .true.

      if (abs(param(21)).eq.1) ifgmres = .true.
      if (abs(param(21)).eq.2) ifcg = .true.
      if (abs(param(21)).eq.3) ifsemg = .true.
      if (param(23).eq.1) iffdm = .true.
      if (param(25).eq.1) ifgfdmdd = .true.

c     Set PML parameters for the Maxwell solver
      pmlthick = int(param(77))
      pmlorder = param(78)
      pmlreferr = param(79)

      if (ifrk22.or.ifrk45) then
         ifrk= .true.
      else
         ifrk = .false.
      endif

      if (nid_io.eq.0) then
         write(6,*) 'iffilter    ',iffilter
         write(6,*) 'ifpoisson   ',ifpoisson
         write(6,*) 'ifgmres     ',ifgmres
         write(6,*) 'ifcg        ',ifcg
         write(6,*) 'ifsemg      ',ifsemg
         write(6,*) 'iffdm       ',iffdm
         write(6,*) 'ific        ',ific
         write(6,*) 'ifsol       ',ifsol
         write(6,*) 'ifsrc       ',ifsrc
         write(6,*) 'ifcentral   ',ifcentral
         write(6,*) 'ifupwind    ',ifupwind
         write(6,*) 'ifdealias   ',ifdealias
         write(6,*) 'iftm        ',iftm
         write(6,*) 'ifte        ',ifte
         write(6,*) 'ifeig       ',ifeig
         write(6,*) 'ifexp       ',ifexp
         write(6,*) 'ifrk45      ',ifrk45
         write(6,*) 'ifrk22      ',ifrk22
         write(6,*) 'ifpec       ',ifpec
         write(6,*) 'ifrob       ',ifrob
         write(6,*) 'ifneu       ',ifneu
         write(6,*) 'ifinhom     ',ifinhom
         write(6,*) 'ifpml       ',ifpml
         write(6,*) 'ifdtn       ',ifdtn
         write(6,*) 'ifrestart   ',ifrestart

         if (ifexp) write(6,*) 'Krylov_Dim =',larnol
      endif

c     Do we have a valid imode?
      imode = 0
      if (ldim.eq.3) then
         imode = 3
         return
      endif
      if (iftm) imode = 2
      if (ifte) imode = 1
      if (imode.eq.0) then
         if (nid_io.eq.0) write(6,*) 'param(4): TE or TM not defined'
         call exitt(1)
      endif
c     Do we have a valid timestepping method?
      if (ifeig.or.ifexp.or.ifrk.or.ifbdf) then
         if (nid_io.eq.0) write(6,*) 'param(17): defined as ', param(17)
      else
         if (nid_io.eq.0) write(6,*) 'param(17): not defined', param(17)
         call exitt(1)
      endif

      return
      end
c-----------------------------------------------------------------------
