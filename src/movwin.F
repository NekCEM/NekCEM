!> \file movwin.F 
!! 
!! \brief Moving window  test ----                        
!! 
!! Moving window implementation                 

!> \brief Brief description of the subroutine 
!! 
!! A more detailed descripton goes here. This is an example 
!! of a detailed description of a method
!!
!! @param param1 a description of the first parameter 
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
!! @see READATW() 
!! @see RDMESH() 
!! @see RDMESHW() 
c---------------------------------------------------------------------
      subroutine movwin_setup             
c---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'WAKEPOT'
      INCLUDE 'MOVWIN'

      n = nx1*ny1*nz1*nelt
      mode=nuniform*param(75)

      if (mod1(istep,mode).eq.mode) then 
         if (nid.eq.0) write(6,*) 'SETUP: MOVING WINDOW FRAME',ifmovw,
     $        mwcnt
      else
         return
      endif

      wetime00=dclock()
      call mwsave_pre 
      call setmovw      ! or call setmovw1 ???
      call copy (hx(1,1,1,1),whx(1,1,1,1),n)
      call copy (hy(1,1,1,1),why(1,1,1,1),n)
      call copy (hz(1,1,1,1),whz(1,1,1,1),n)
      call copy (ex(1,1,1,1),wex(1,1,1,1),n)
      call copy (ey(1,1,1,1),wey(1,1,1,1),n)
      call copy (ez(1,1,1,1),wez(1,1,1,1),n)
      wetime11=dclock()
      if (nid.eq.0) write(6,*) 'MOVING WINDOW SETUP TIME',
     $                          wetime11-wetime00

      return
      end
c---------------------------------------------------------------------
      subroutine READATW            
c---------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      INCLUDE 'ZPER'
      include 'PARALLEL'
      include 'EMWAVE'
      include 'MOVWIN'

      open (unit=9,FILE=REAFLE,STATUS='OLD')
 
      wetime1=dclock()
      call RDPARAM   
      wetime2=dclock()
      !write(6,*) 'RDparam=',wetime2-wetime1
      call RDMESHW
      wetime3=dclock()
      !write(6,*) 'RDmesh =',wetime3-wetime2
      call RDCURVW
      wetime4=dclock()
      !write(6,*) 'RDcurv =',wetime4-wetime3
      call RDBDRYW
      wetime5=dclock()
      !write(6,*) 'RDbdry =',wetime5-wetime4
      call RDICDF  
      wetime6=dclock()
      !write(6,*) 'RDicdf =',wetime6-wetime5
      call RDMATP 
      wetime7=dclock()
      !write(6,*) 'RDmat  =',wetime7-wetime6
      call RDHIST 
      wetime8=dclock()
      !write(6,*) 'RDhist =',wetime8-wetime7
      call RDOUT  
      wetime9=dclock()
      !write(6,*) 'RDout  =',wetime9-wetime8
      call RDOBJ  
      wetime10=dclock()
      !write(6,*) 'RDobj  =',wetime10-wetime9

      close(unit=9)

      return
      end 


!> \brief Read moving window mesh 
!! 
!! Read the moving window mesh. Skips elements not
!! in the current window.
!!
!! @param param1 a description of the first parameter 
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
c---------------------------------------------------------------------
      subroutine RDMESHW
c---------------------------------------------------------------------
c     read mesh for moving window
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'SCRCT'
      include 'TSTEP'
      include 'ZPER'
      include 'EMWAVE'
      include 'MOVWIN'
      include 'WAKEPOT'

      integer  icalld 
      save     icalld 
      data     icalld/0/ 

      READ(9,*,ERR=400,END=400) XFAC,YFAC,XZERO,YZERO
      READ(9,*,ERR=500)
      READ(9,*,ERR=500,END=500) NELGT,NDIM,NELGV

c ... basic setup to read mesh,bdry,curve
      msize = param(86)
      nwin2d= lelg2d
      mw1   = nwin2d* icalld+1
      mw2   = nwin2d*(icalld+msize)
      icalld= icalld+1
      if (nid.eq.0) write(6,*) 'number MOVWIN calld=',icalld
      mwcnt = icalld
      
      nelgt_mw= nelgt
      call izero(gllel_mw,nelgt_mw)

      iel=0
      do ieg= 1,nelgt_mw
         gllel_mw(ieg)=0
         if ((ieg.ge.mw1).and.(ieg.le.mw2))  then
            iel  = iel+1
            gllel_mw(ieg)=iel
            nelgt=iel    ! assign nelgt newly
            nelgv=iel    ! assign nelgt newly
         endif
      enddo
      
c ... read map                  
      ! gllnid(ieg)=0   ! temp 
      ! gllnid has to be assigned with xxt_elm_to_procw 
      ! when it reads the map data: 1st column    

      CALL MAPELPR() 

c ... read moving window data
      nsides= ndim*2

      do ieg= 1,nelgt_mw                ! global elt num
 
           if (GLLNID_MW(IEG).EQ.NID) THEN
             
            iel= gllel_mw(ieg)                                         
            iel= gllel(iel)                                           

            !if (nid.eq.0) write(6,*) 'gll_mw::',nid,iel,ieg            
            READ(9,30,ERR=500,END=500) IGROUP(IEL)
   30       FORMAT(43X,I5)

            IF     (NDIM.EQ.2)THEN
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=1,4)
               call rzero                (zc(1 ,iel)     ,4)
            ELSE IF(NDIM.EQ.3)THEN
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (ZC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=5,8)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=5,8)
               READ(9,*,ERR=500,END=500) (ZC(IC,IEL),IC=5,8)
               !write(6,*) 'mw::',nid,ieg,iel,(ZC(IC,IEL),IC=1,1)
            ENDIF


         else ! Skip over data for element NOT on this window          

            READ(9,41,ERR=500,END=500) ADUM
            IF     (NDIM.EQ.2)THEN
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
            ELSE IF(NDIM.EQ.3)THEN
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
            ENDIF

         endif
      enddo

   41 format(A1)
 
C     End of mesh read.
 
      RETURN
 
C     Error handling:
 
  400 CONTINUE
      WRITE(6,401)
  401 FORMAT(2X,'ERROR READING MOVING WINDOW, CHECK READ FILE'
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt
C
  500 CONTINUE
      WRITE(6,501) IEG
  501 FORMAT(2X,'ERROR READING MOVING MESH DATA NEAR ELEMENT',I7
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt
C
      RETURN


      return
      end

c---------------------------------------------------------------------
      subroutine RDCURVW
C     Read curve side data
      INCLUDE    'SIZE'
      INCLUDE    'INPUT'
      INCLUDE    'PARALLEL'
      INCLUDE    'MOVWIN'
      CHARACTER*1 ANS

      READ(9,*)
      READ(9,*)NCURVE
      CALL RZERO(CURVE ,48*LELT)
      CALL BLANK(CCURVE, 8*LELT)

      IF (NCURVE.GT.0) THEN

         DO 50 ICURVE=1,NCURVE

            IF     (NELGT_MW.LT.1000) THEN
               READ(9,60,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSEIF (NELGT_MW.LT.1000000) THEN
               READ(9,61,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSEIF (NELGT_MW.LT.10000000) THEN
               READ(9,62,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ENDIF

   60       FORMAT(I3,I3,5G14.6,1X,A1)
   61       FORMAT(I2,I6,5G14.6,1X,A1)
   62       FORMAT(I1,I7,5G14.6,1X,A1)

           if (GLLNID_MW(IEG).EQ.NID) THEN
             
                IEL= GLLEL_MW(IEG)
                IEL= GLLEL(IEL)   ! (IEL)=GLLEL_MW(IEG)? to be checked.         
                CURVE (1,IEDG,IEL)=R1
                CURVE (2,IEDG,IEL)=R2
                CURVE (3,IEDG,IEL)=R3
                CURVE (4,IEDG,IEL)=R4
                CURVE (5,IEDG,IEL)=R5
                CCURVE(  IEDG,IEL)=ANS
               ! if (nid.eq.0) write(6,*) 'CURVE::',IEG,IEL,R1
            ENDIF

   50    CONTINUE

      ENDIF

      RETURN

C     Error handling:
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING MOVING WINDOW CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
      RETURN

      end

c---------------------------------------------------------------------
      subroutine RDBDRYW
c---------------------------------------------------------------------
C     .Read Boundary Conditions (and connectivity data)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      INCLUDE 'EMWAVE'
      INCLUDE 'MOVWIN'

      CHARACTER    CBC1*1,CBC3*3,CHTEMP*1,CHTMP3*3
      EQUIVALENCE (CHTEMP,CHTMP3)
      character*80 string
 
C     Set up TEMPORARY value for NFIELD - NFLDT
      NFLDT  = 1 
      IF (IFHEAT) NFLDT=2+NPSCAL
      NBCS   = NFLDT
      IBCS   = 2
      IF (IFFLOW) IBCS = 1
      NSIDES = 2*NDIM

C     Read boundary conditions for all fields
      LCBC=18*LELT*(LDIMT1 + 1)
      LRBC=30*LELT*(LDIMT1 + 1)

      CALL RZERO(BC ,LRBC)
      CALL BLANK(CBC,LCBC)
 
      READ(9,*,ERR=500,END=500)  !   ***** BOUNDARY CONDITIONS *****
      ibcnew =1 

      write(6,*) 'MOVWIN-------BDRY'                   
      DO 100 IFIELD=ibcnew,NBCS  !   DO 100 IFIELD=IBCS,NBCS
        NEL=NELGT_MW
        read(9,81) string        !   ***** FLUID   BOUNDARY CONDITIONS *****
        call capit(string,80)

        IF (indx1(string,'NO ',3).eq.0) then ! we have acitve bc info
         IF(VNEKTON .LE. 2.52) NBCREA = 3
         IF(VNEKTON .GE. 2.55) NBCREA = 5
 
         DO 80 IEG=1,NEL         
         DO 80 ISIDE=1,NSIDES

            IF (GLLNID_MW(IEG).EQ.NID) THEN
               IEG0=GLLEL_MW(IEG)
               IEL=GLLEL(IEG0)

               IF     (NELGT_MW.LT.1000) THEN
                  READ(9,50,ERR=500,END=500)
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
                  READ(6,50,ERR=500,END=500)
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
               ELSEIF (NELGT_MW.LT.100000) THEN
                  READ(9,51,ERR=500,END=500)
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
               ELSEIF (NELGT_MW.LT.10000000) THEN
                  READ(9,52,ERR=500,END=500)
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
               ENDIF

   50          FORMAT(A1,A3,2I3,5G14.7)
   51          FORMAT(A1,A3,I5,I1,5G14.7)
   52          FORMAT(A1,A3,I6,5G14.7)

               IF (CHTEMP.NE.' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP 

               !check for fortran function as denoted by lower case bc's:
               CBC1=CBC(ISIDE,IEL,IFIELD)
               CBC3=CBC(ISIDE,IEL,IFIELD)
               ICBC1=ICHAR(CBC1)
               IF (ICBC1.GE.97.AND.ICBC1.LE.122) THEN
                  IF(CBC3(3:3).NE.'i')NLINES=BC(1,ISIDE,IEL,IFIELD)
                  IF(CBC3(3:3).EQ.'i')NLINES=BC(4,ISIDE,IEL,IFIELD)
                  DO 60 I=1,NLINES
   60             READ(9,*,ERR=500,END=500)
               ENDIF
            ELSE
               READ(9,*,ERR=500,END=500)   cbc1  ! dummy read                          
            ENDIF

   80    CONTINUE
        ENDIF

   81   FORMAT(a80)
  100 CONTINUE
      
      write(6,*) 'MOVWIN-------BDRY end'                        
C     END OF BC READ

      !DO IEG= 1,NELGT_MW         
         !IEL= GLLEL_MW(IEG)

      DO IEG= 1,NELGT         

         IEL= GLLEL(IEG)

      DO ISIDE=1,NSIDES

         if     ((IEG.LE.nwin2d).AND.(ISIDE.EQ.5)) THEN

             CBC(ISIDE,IEL,2)= 'NA '             
             DO II=1,NBCREA
                BC (II,ISIDE,IEL,2)=0.0                 
             ENDDO
             write(6,111) nelgt,iel,iside,CBC(ISIDE,IEL,2)

         elseif ((IEG.GT.NELGT-nwin2d).AND.(ISIDE.EQ.6)) THEN

             CBC(ISIDE,IEL,2)= 'NA ' ! temporary             
             DO II=1,NBCREA
                BC (II,ISIDE,IEL,2)=0.0                  
             ENDDO
             write(6,112) nelgt,iel,iside,CBC(ISIDE,IEL,2)

         endif    

      ENDDO
      ENDDO
  111 FORMAT('BC face 5::',3I5,'::  ',A3)
  112 FORMAT('BC face 6::',3I5,'::  ',A3)

C     Check for dummy line:  "NO THERMAL B.C.'S"
      IF (NFLDT.EQ.1) READ(9,*,ERR=500,END=500)
      RETURN
 
C     Error handling:
  500 CONTINUE
      WRITE(6,501) IFIELD,IEG
  501 FORMAT(2X,'ERROR READING MOVING MESH BDRY CONDITIONS FOR FIELD',
     $       I4,I6,/,2X,'ABORTING IN ROUTINE RDBDRY.')
      call exitt

      RETURN
      end

c---------------------------------------------------------------------
      subroutine mwsave_pre 
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'MOVWIN'

      integer  icalld 
      save     icalld 
      data     icalld/0/
      integer  e,e1,e2
   
      n   = nx1*ny1*nz1*nelt
      nxyz= nx1*ny1*nz1              

      nwin2d= lelg2d
      if (nid.eq.0) write(6,*) 'nwin2d=',nwin2d,nelgt,nelt

      do ieg = nwin2d+1,nelgt               
         k = nwin2d
         e1 = gllel(ieg)
         e2 = gllel(ieg-k)
         do i= 1,nxyz
            whx(i,1,1,e2 )= hx(i,1,1,e1) 
            why(i,1,1,e2 )= hy(i,1,1,e1) 
            whz(i,1,1,e2 )= hz(i,1,1,e1) 
            wex(i,1,1,e2 )= ex(i,1,1,e1) 
            wey(i,1,1,e2 )= ey(i,1,1,e1) 
            wez(i,1,1,e2 )= ez(i,1,1,e1) 
           !wez(i,1,1,e-k)= ez(i,1,1,e) 
         enddo
      enddo

      do ieg  = nelgt-nwin2d+1,nelgt
            e = gllel(ieg)
         do i = 1,nxyz
            whx(i,1,1,e)= 0           
            why(i,1,1,e)= 0           
            whz(i,1,1,e)= 0           
            wex(i,1,1,e)= 0           
            wey(i,1,1,e)= 0           
            wez(i,1,1,e)= 0           
         enddo
      enddo

      if (mwcnt.eq.4) then 
      call cem_out
      if (nid.eq.0) write(6,*) 'mwcnt=',mwcnt
      call exitt
      endif        
    
      return
      end
c---------------------------------------------------------------------
      subroutine mwsave_aft 
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'MOVWIN'

      integer  icalld 
      save     icalld 
      data     icalld/0/
      integer  e
   
      n   = nx1*ny1*nz1*nelt
      nxyz= nx1*ny1*nz1              

      nwin2d= lelg2d
      write(6,*) 'nwin2d=',nwin2d,nelt

      filtw=0.5
      call copy (hx,whx,n)
      call copy (hy,why,n)
      call copy (hz,whz,n)
      call copy (ex,wex,n)
      call copy (ey,wey,n)
      call copy (ez,wez,n)

      filtw=0.5
c     call cem_filter(filtw)

      return 
      end   
c---------------------------------------------------------------------
      subroutine setmovw1
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'WAKEPOT'
C
C     Declare scratch arrays
C
      integer    lzdt  ! temp

      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LPSC9=LDIMT+9)
      PARAMETER (LZDT=100000) !temp
      parameter (lctmp1n = LXYZR*LPSC9)
      COMMON /CTMP1/ DUMMY1(LCTMP1n)
      common /ccpu/ cpu_t,cpu_t0
      real*4  dummy1

      common  zdt(LZDT)      !temp
C
      COMMON /CTMP0/ DUMMY0(LCTMP0)
      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)

      REAL    e, oe
      integer WDS
      integer IGEOM
      integer NTOTT, NSRF
      REAL*8  t0,t2,pi
      REAL*8  HXer, HYer,EZer
      REAL*8  EXer, EYer,HZer
      REAL*8  HXgp, HYgp,EZgp
      REAL*8  EXgp, EYgp,HZgp
      real    dt0
      save    dt0

      integer icalld0
      save    icalld0
      data    icalld0/0/

      
c ... Data initialization

      call INITDAT      ! ifoptions  
      call FILES        ! read files

      call READATW

      call SETVAR
      call ECHOPAR

      call CONNECT

      igeom = 2
      call gengeom (igeom)!! temorary blocking for hmholtz 11/06/06 misun 

      call geom_reset(1)  ! not using this for XM3
      call set_options

      call usrdat2

      call SET_UNR   
      call VRDSMSH
      call SETLOG
      call BCMASK
      call load_semhat_weighted

      call RK_STORAGE

      ntot = nx1*ny1*nz1*nelv
      nsrf = 2*ndim*nx1*nz1*nelv  ! to check diff between elmts

      call cem_uvp

      if (ifpml ) call pml_setup
      if (ifbeam) call beam_setup

      call get_dxmin (dxmin)
c
c     dt = param(12)    ! defalut
c     twopi = 2.0*pi
c     cfl = 1.0/(sqrt(twopi))
c
      if (ifwake)   mstep = param(75)     ! preset for wake
      if (ifwake)   call cem_wakepotential! preset for wake
      if (ifwake)   dt = dzt/mstep        ! reset dt for wake
      if (nid.eq.0) write(6,*) 'default DT=',dt
 
      if (param(14).eq.0) nsteps= param(11)
      if (param(14).eq.1) nsteps= nsteps


      return
      end
c---------------------------------------------------------------------
      subroutine setmovw
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'WAKEPOT'
C
C     Declare scratch arrays
C
      integer    lzdt  ! temp

      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LPSC9=LDIMT+9)
      PARAMETER (LZDT=100000) !temp
      parameter (lctmp1n = LXYZR*LPSC9)
      COMMON /CTMP1/ DUMMY1(LCTMP1n)
      common /ccpu/ cpu_t,cpu_t0
      real*4  dummy1

      common  zdt(LZDT)      !temp
C
      COMMON /CTMP0/ DUMMY0(LCTMP0)
      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)

      REAL    e, oe
      integer WDS
      integer IGEOM
      integer NTOTT, NSRF
      REAL*8  t0,t2,pi
      REAL*8  HXer, HYer,EZer
      REAL*8  EXer, EYer,HZer
      REAL*8  HXgp, HYgp,EZgp
      REAL*8  EXgp, EYgp,HZgp
      real    dt0
      save    dt0

      integer icalld0
      save    icalld0
      data    icalld0/0/

      
c ... Processor initialization 

c     call INIPROC
c     call comm_init()
c
c     if (nid.eq.0) then
c         write(6,*) 'Number of Processors ::',np
c     endif


c ... Data initialization

      call INITDAT      ! ifoptions  
      call FILES        ! read files

      call READATW

      if (nid.eq.0) then
          write(6,*) 'readat time ::',dclock()-t0,' seconds'
      endif

c     call SETVAR

      call ECHOPAR

c ... Check for negative IOSTEP to turn on dumping feature

      instep  = 1
      fintime = param(10)

      if (nid.eq.0) write(6,*) '22  mwcnt=',mwcnt
      if (nsteps.eq.0 .and. fintime.eq.0.) instep=0

c ... Geometry initialization

      IGEOM   = 2

      if (nid.eq.0) write(6,*) 'if3d: ',if3d,nelgv,nelt,nelv,nx1

      if (nid.eq.0) write(6,*) '33 mwcnt=',mwcnt
c ... connect1.f-initds(face ordering),dsset(skpdat),setedge,genxyzl,geom1,geodat
      call CONNECT

c ... call GENWZ ! block due to redundant call;been called in connect.
      call usrdat

c ... preset for arrays: beam and wake
c     call arrays_setup   ! need for usrdat2 ???          

      if (nid.eq.0) write(6,*) '44 mwcnt=',mwcnt
c ... coef.f - lagmass(BM1LAG),geom1,setinvm(BINVM1/BINVTM1),setdef(deform)
      call gengeom (igeom)!! temorary blocking for hmholtz 11/06/06 misun 

c ... option checking  
      call set_options

c ... deformation
      call usrdat2

c ... ic.f
      call geom_reset(1)  ! not using this for XM3

c ... coef.f - normal, tangent vectors
      call SET_UNR   ! check rxm1, rym1, sxm1,

      if (nid.eq.0) write(6,*) nid,' call vrdsms',instep
      call VRDSMSH

c ... Field initialization

      if (nid.eq.0) write(6,*) nid,' call setlog',instep
      call SETLOG

c ... bdry.f - boundary mapping
      call BCMASK

c ... semHAT
      call load_semhat_weighted

c ... Need eigenvalues to set tolerances in presolve (SETICS)

      if (nid.eq.0) write(6,*) 'this is ifflow:',ifflow,nsteps

c ... IF (FINTIM.NE.0.0.OR.NSTEPS.NE.0) CALL GENEIG (IGEOM)

c ... cem start 
      ntot = nx1*ny1*nz1*nelv
      nsrf = 2*ndim*nx1*nz1*nelv  ! to check diff between elmts

c ... restart    

c ... initialize
c     call cem_init

c ... assign material property
      call cem_uvp

c ... initialize option setup
      !ifpml=.true.  ! TEMP
      if (ifpml ) call pml_setup

      !ifpml=.false.  ! TEMP
      if (ifbeam) call beam_setup

      !user param= zsigma,rsigma
      !auto:zposition,zsigma,rsigma,zcut,wakerad1,wakerad2                   

c ... initial conditions                                   
c     call userini(time,hx,hy,hz,ex,ey,ez)

c ... initialize
      if (icalld0.eq.0) then

c ... get minimum mesh size
      call get_dxmin (dxmin)

      if (nid.eq.0)  write(6,*) 'dxmin=',dxmin
      if (nid.eq.0)  write(6,*) 'cspeed=',cspeed
 
c ... get dt                       
      dt = param(12)    ! defalut
      twopi = 2.0*pi
      cfl = 1.0/(sqrt(twopi))

      if (nid.eq.0)  write(6,*) 'dt=',dt
      if (nid.eq.0)  write(6,*) 'before mstep=',mstep

      if (ifvardt)  dt = dxmin/cspeed     ! auto dt
      if (ifvardt)  fldmax= param(18)     ! auto dt
      if (ifvardt)  idtchg= 1

      if (ifwake)   mstep = param(75)     ! preset for wake
      if (nid.eq.0)  write(6,*) 'after mstep=',mstep

      if (ifwake)   call cem_wakepotential! preset for wake
      if (ifwake)   dt = dzt/mstep        ! reset dt for wake
      if (ifwake)   then
         if (nid.eq.0) write(6,*)
         if (nid.eq.0) write(6,*) 'check and adjust dt >>>>>>'
         if (nid.eq.0) write(6,*) 'c*dt=',cspeed*dt, dxmin*cfl
         if (cspeed*dt.gt.dxmin*cfl) then
            if (nid.eq.0) then
               write(6,*) 'small luniform or mstep',luniform,dt,dxmin
               call exitt
            endif
         endif
      endif
      if (nid.eq.0) write(6,*) 'default DT=',dt
      if (nid.eq.0) write(6,*)
      icalld0=1

      endif

c ... update initial if analytic solutions                
c     if (ifsol) call usersol(time,hx,hy,hz,ex,ey,ez)

c ... restart: right now using only for aps mesh
      if (ific)  call setics

c ... print and output initial fields 
      if (.not.ifvardt) call userchk

c ... reduce dt by a half dt                                
      if ((ifvardt).and.(idtchg.gt.1)) dt = dt/2.0
      if ((ifvardt).and.(idtchg.gt.1).and.(nid.eq.0))  write(6,*)
      if ((ifvardt).and.(idtchg.gt.1).and.(nid.eq.0))
     $     write(6,*) ' new dt=',dt,idtchg

c ... Start time stepping

      !call TIME00
      call opcount(2)
      call dofcnt
      call RK_STORAGE

c ... write steatement showing start

      if  (nid.eq.0)  write(6,*) '-------------'
      if  (nid.eq.0)  write(6,*) '    BEGIN    '
      if  (nid.eq.0)  write(6,*) '-------------'

c ... time info                        
      !cpu_t0=dclock()

c ... auto nsteps=0
      if (param(14).eq.0) nsteps= param(11)
      if (param(14).eq.1) nsteps= nsteps

C     write(6,*) 'SETUP: MOVING WINDOW FRAME',ifmovw

      return
      end
c---------------------------------------------------------------------

