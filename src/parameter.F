
c-----------------------------------------------------------------------
      subroutine set_logics
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      include 'NEKUSE'

      IFSE     =.false.
      IFMAXWELL=.false.
      IFSCHROD =.false.
      IFDRIFT  =.false.
      IFVARDT  =.false.
      IFFILTER =.false.
      IFDIELEC =.false.
      IFIC     =.false.
      IFRESTART=.false.
      IFDTCTL  =.false.
      IFSCAT   =.false.
      IFSFTF   =.false.
      IFDRUDE  =.false.
      IFLOREN  =.false.
      IFHYDRO  =.false.
      IFBEAM   =.false.
      IFWAKE   =.false.
      IFNAPY   =.false.
      IFWAVEGUIDE=.false.
      IFPBG    =.false.
      IFWAKET  =.false.
      IFVTK    =.false.

      IFRESTART=.false.
      IFFILTER =.false.
      IFDEALIAS=.false.
      IFUPWIND =.false.
      IFCENTRAL=.false.

      IFSOL    =.false.
      IFSRC    =.false.
      IFPOISSON=.false.
      IFGMRES  =.false.
      IFCG     =.false.
      IFGFDM   =.false.

      IFRK45   =.false.
      IFRK44   =.false.
      IFRK33   =.false.
      IFRK22   =.false.
      IFRK     =.false.
      IFBDF    =.false.
      IFBDF1   =.false.
      IFBDF2   =.false.
      IFSIM5   =.false.
      IFSIM4   =.false.
      IFSIM3   =.false.
      IFSIM2   =.false.
      IFSIM1   =.false.
      IFSIM    =.false.
      IFEXP    =.false.
      IFEIG    =.false.

      IFTE     =.false.
      IFTM     =.false.
      IFNM     =.false.

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
      if ((param(5).ge.30).and.(param(5).lt.40)) IFSE     =.true.

c     if (param( 14).eq.0)  IFPOISSON=.true.
c     if (param( 14).eq.1)  IFSCAT   =.true.
c     if (param( 14).eq.3)  IFSFTF   =.true.

      if (param( 16).eq. 1) IFSOL    =.true.
      if (param( 16).ne. 1) IFSRC    =.true.

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
      if (IFBDF1.or.IFBDF2) IFBDF    =.true.

      if (param( 18).eq. 1) IFFILTER =.true.
      if (param( 18).eq.-1) IFDEALIAS=.true.
      if (param( 19).eq. 0) IFUPWIND =.true.
      if (param( 19).eq. 1) IFCENTRAL=.true.
      if (param( 21).eq. 0) IFGMRES  =.true.
      if (param( 21).eq. 1) IFCG     =.true.
      if (param( 84).ne. 0) IFRESTART=.true.
      if (param(103).eq. 1) IFNM     =.true.
c     if (param(116).ne. 0) IFGFDM   =.true.

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
