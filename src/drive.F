!> \file drive.F 
!! 
!! \brief main driver file for NEKCEM 
!! 
!! This file contains the application entry point for NEKCEM 


!> \brief A brief description goes here.
!! 
!! A more detailed descripton goes here. This is an example 
!!
!! @param param1 a description of the first parameter 
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
!! @see gengeom() 

#ifdef NO_C_MAIN
      PROGRAM NEKTON
#else
      SUBROUTINE DRIVE_NEKTON
#endif

C
C     CEM Simulation
C

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'WAKEPOT'
      INCLUDE 'MOVWIN' 
C
C     Declare scratch arrays
C
c...  start: to be kept temporarily
      INTEGER    lzdt  ! temp
      PARAMETER (LXR=1)
      PARAMETER (LYR=1)
      PARAMETER (LZR=1)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LPSC9=1)
      PARAMETER (LZDT=1) !temp
      parameter (lctmp1n = LXYZR*LPSC9)
      common  zdt(LZDT)      !temp

c     PARAMETER (LXR=LX1+6)
c     PARAMETER (LYR=LY1+6)
c     PARAMETER (LZR=LZ1+6)
c     PARAMETER (LXYZR=LXR*LYR*LZR)
c     PARAMETER (LPSC9=LDIMT+9)
c     PARAMETER (LZDT=100000) !temp

      COMMON /CTMP1/ DUMMY1(LCTMP1n)
      COMMON /CTMP0/ DUMMY0(LCTMP0)
      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)
      real*4  dummy1
c...  end: to be deleted                
 
      common /ccpu/ cpu_t,cpu_t0
      REAL    e, oe
      integer WDS
      integer IGEOM
      integer NTOTT, NSRF 
      real*8  t0,t2,pi

      integer icalld0
      save    icalld0
      data    icalld0/0/
 
c ... Processor initialization 
      call iniproc
      call comm_init()
      if (nid.eq.0) write(6,*) 'Number of Processors ::',np
      call opcount(1)

      time0  = dclock()
      etimes = dclock()

c ... Data initialization
      t0     = dclock()
      call initdim      ! nx=lx1, ny=ly1, nz=lz1 
      call initdat      ! ifoptions  
      call files        ! read files
      if (nid.eq.0) write(6,*) 'initialization',dclock()-t0,' sec'

c ... Read data file                                          
      t0 = dclock()
      call readat 
      if (nid.eq.0) write(6,*) 'readat time ::',dclock()-t0,' sec'

c ... Set parameters                                         
      call setvar
      call echopar

c ... Check for negative IOSTEP to turn on dumping feature
      instep  = 1
      istep   = 0
      fintime = param(10)

      if (nsteps.eq.0 .and. fintime.eq.0.) instep=0 ! not using
 
c ... Geometry initialization
      igeom   = 2
      call connect       ! glonum,dsset(skpdat),setedge,etc
c ... call GENWZ         ! block: redundant call in connect.
      call usrdat        ! usr file
      call set_options   ! option logical value setting
      call SETLOG        ! Field initialization
      call gengeom(igeom)! temorary blocking for hmholtz 11/06/06 misun 
      call usrdat2       ! deformation
      call geom_reset(1) ! ic.f:  not using this for XM3
      call set_unr       ! normal, tangent vectors
      call vrdsmsh       ! validation mesh
      call bcmask        ! bdry.f - boundary mapping
      call load_semhat_weighted !semHAT


c ... CEM START 
      ntot = nx1*ny1*nz1*nelv
      NSRF = 2*ndim*nx1*nz1*nelv  ! to check diff between elmts

c ... restart for moving window    
 100  istep= 0
      time = 0.0

      call cem_init      ! initialize
      call cem_uvp       ! assign material property
   
c ... setup for some options
      if (ifpml ) call pml_setup                                     
      if (ifpml ) call cem_outer_pmlsrf
      if (ifbeam) call beam_setup 

c ... initial conditions                                   
      call userini(time,hx,hy,hz,ex,ey,ez) 

c ... initialize
      if (icalld0.eq.0) then
          dt= param(12)          ! defalut
          call get_dxmin (dxmin) ! get minimum mesh size
          call adjust_dt         !           
          icalld0= 1
      endif

        
      if (ifsol) call usersol(time,hx,hy,hz,ex,ey,ez) ! if analytic solutions                
      if (ific)  call setics                          ! restart:only for aps mesh

c ... print and output initial fields 
      call userchk
      call cem_out

c ... countings and rk-steup
      call time00
      call opcount(2)
      call dofcnt
      call rk_storage 

c ... auto nsteps=0
      !if (param(14).eq.0) nsteps= param(11) 
      !if (param(14).eq.1) nsteps= nsteps     
      !if (nid.eq.0)  write(6,*) 'nsteps=',nsteps
 
c ... Start time stepping              
      cpu_t=0.0          

      DO 1000 istep=1,nsteps 

         !if (ifmovw) call movwin_setup 

         cpu_t0=dclock()
         call cem (igeom)
         cpu_t= cpu_t+dclock()-cpu_t0

         time=time+dt

         call cem_out
         call userchk

         if (lastep.eq.1)     goto 1001
         if (time.ge.fintime) goto 1001

         rtime= (cpu_t/istep)*(nsteps-istep)/60.0            
         if ((nid.eq.0).and.(mod1(istep,iostep).eq.iostep)) 
     $        write(6,1111) nsteps, rtime  

 1000 CONTINUE
 1001 CONTINUE

 1111 format(/,' max timestep=',i9,
     $        '; estimated remaining time (min)=',f10.2)

c ... End time stepping
      call userchk                                        

c ... Print                           
      call printout          

      END
 
c-----------------------------------------------------------------------
      subroutine adjust_dt
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'WAKEPOT'
      INCLUDE 'MOVWIN' 

      call get_dxmin (dxmin) 

      if (nid.eq.0) write(6,*) 'dxmin =',dxmin
      if (nid.eq.0) write(6,*) 'cspeed=',cspeed
      if (nid.eq.0) write(6,*) 'before mstep=',mstep
      if (nid.eq.0) write(6,*) 'default dt=',dt

      cfl = 1.0
      if (param(21).gt.0)  cfl= param(21)

      dt = cfl*dxmin/cspeed
      if (nid.eq.0) write(6,*) 'cfl=',cfl
      if (nid.eq.0) write(6,*) 'dt=cfl*dxmin/cspeed=',dt

      if (ifvardt)  then
           dt = cfl*dxmin/cspeed ! auto dt
           fldmax= param(18) ! auto dt
           idtchg= 1
      endif

      if (param(14).gt.0) call dt_adjusted_nano
      if (param(21).eq.0) then
          nsteps=param(11)
          dt= param(12)
      endif

      if (param(50).ne.0) then
         n=nx1*ny1*nz1*nelt
         dist=cspeed*dt
         if     (param(50).eq.1) then
           xmax=glmax(xm1,n)
           xmin=glmin(xm1,n)
           !nsteps=(xmax-xmin)/dist
         elseif (param(50).eq.2) then
           ymax=glmax(ym1,n)
           ymin=glmin(ym1,n)
           !nsteps=(ymax-ymin)/dist
         elseif (param(50).eq.3) then
           zmax=glmax(zm1,n)
           zmin=glmin(zm1,n)
           !nsteps=10*(zmax-zmin)/dist
         else
           if (nid.eq.0) write(6,*) 'param(50)=0,1,2?',param(50)
           call exitt
         endif
         if (nid.eq.0) write(6,*) 'param(50)=',param(50)
         if (nid.eq.0) write(6,*) 'updated nsteps for FT',nsteps
      endif

      if (nid.eq.0) then
          write(6,*) 'used DT=',dt
          write(6,*) '             '                         
          write(6,*) '-------------'
          write(6,*) '    BEGIN    '
          write(6,*) '-------------'
      endif

      if (ifwake) then
         mstep = param(75) ! preset for wake
         if (nid.eq.0) write(6,*) 'after mstep=',mstep
         call cem_wakepotential! preset for wake
         dt = dzt/mstep        ! reset dt for wake

         if (nid.eq.0) write(6,*)
         if (nid.eq.0) write(6,*) 'check and adjust dt >>>>>>'
         if (nid.eq.0) write(6,*) 'c*dt=',cspeed*dt, dxmin*cfl
         if (cspeed*dt.gt.dxmin*cfl) then
            if (nid.eq.0) then
               write(6,*) 'small luniform or mstep',luniform,dt,dxmin
               !call exitt
            endif
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine chkmovwin
c-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'MOVWIN'

      ifmovw=.false.
      if (lmov.eq.1) then 
          ifmovw=.true.
          if (nid.eq.0) write(6,*) 'IFMOVW=',ifmovw
      endif

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
c     call exitt

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

      NXD = LXD 
      NYD = LYD 
      NZD = LZD 

      MARNOL=LARNOL  ! arnoldi 

      NPERT = LPERT

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
      IFFLOW = .FALSE.
      IFMVBD = .FALSE.
      IFHEAT = .FALSE.
      IFSPLIT= .FALSE.
      IFDOIT = .FALSE.
      IFDG   = .FALSE.
      IFCEM  = .FALSE.
C


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
      PARAM(86) = 0
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
      subroutine echopar
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
      subroutine gengeom (igeom)
C----------------------------------------------------------------------
C
C     Generate geometry data
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'GEOM'
      include 'WZ'
      include 'MASS'
C
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C
      ntot=85750                    
      CALL GENCOOR (XM3,YM3,ZM3)
      CALL GEOM1   (XM3,YM3,ZM3)
      CALL GEOM2
C     CALL UPDMSYS (1)
      CALL VOLUME
      CALL SETINVM
C
      RETURN
      END
C
      subroutine files
C----------------------------------------------------------------------
C
C     Defines machine specific input and output file names.
C
C----------------------------------------------------------------------
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
      OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD')
      READ(8,10) SESSION
      READ(8,10) PATH
      CLOSE(UNIT=8)
   10 FORMAT(A132)


C
c     call hmt_fix_path(132,PATH1)
C
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
      call hmt_set_file_names(132,reafle)
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
C
      subroutine time00
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

      subroutine timeout_OLD
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
      subroutine opcount(ICALL)
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
      subroutine dofcnt
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
      if (nid.eq.0) write(6,*) ' dofs:',nvtot,nptot,netot
      return
      end
c-----------------------------------------------------------------------
      subroutine geom_reset(icall)
C
C     Generate geometry data
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
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
      CALL GEOM1 (XM3,YM3,ZM3)
      CALL GEOM2
c     CALL UPDMSYS (1)
      CALL VOLUME
      CALL SETINVM
      CALL SETDEF
c     CALL SFASTAX
c
      return
      end
c-----------------------------------------------------------------------
      subroutine map13_all(x3,x1)
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
