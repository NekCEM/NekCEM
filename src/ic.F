c-----------------------------------------------------------------------
      subroutine setics
C-----------------------------------------------------------------------
C
C     Set initial conditions.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'IXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'MASS'
      INCLUDE 'MVGEOM'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
      INCLUDE 'EMWAVE'
c
      logical  iffort(  ldimt1,0:lpert)
     $       , ifrest(0:ldimt1,0:lpert)
     $       , ifprsl(  ldimt1,0:lpert)
c
      LOGICAL  IFANYP
      LOGICAL  IFPERT  ! FIXME temporary misun 6/4/07

      COMMON /RDUMP/ NTDUMP
      COMMON /CTMP1/ WORK(LX1,LY1,LZ1,LELV)
     $ ,             TA1 (LX2,LY1,LZ1)
     $ ,             TA2 (LX2,LY2,LZ1)
      common /VPTS1/ vxp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              ,vyp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              ,vzp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              , tp(lpx1*lpy1*lpz1*lpelt,ldimt,lpert)


      ifpert=.false. ! FIXME temporary   : misun 6/4/07

C     Initialize all fields:
C
      NTOTP=NX2*NY2*NZ2*NELV
      NTOTV=NX1*NY1*NZ1*NELV
      NTOTT=NX1*NY1*NZ1*NELT
c
      CALL RZERO(VX,NTOTV)
      CALL RZERO(VY,NTOTV)
      CALL RZERO(VZ,NTOTV)

      DO 10 IFLD=1,LDIMT
         CALL RZERO(T(1,1,1,1,IFLD),NTOTT)
   10 CONTINUE
c
      jp = 0
C
      irst = param(46)
      if (irst.gt.0) return
c
c
c     If moving geometry then add a perturbation to the
c     mesh coordinates (see Subroutine INIGEOM)
c
C
C     Find out what type of i.c. is requested
C     Current options: 
C
C     (1) - User specified fortran function (default is zero i.c.)
C     (2) - Restart from file(s)
C     (3) - Activate pre-solver => steady diffusion / steady Stokes
C
C     If option (2) is requested, also return with the name of the
C     restart file(s) together with the associated dump number
C
      call slogic (iffort,ifrest,ifprsl,nfiles)
C
C     Set up proper initial values for turbulence model arrays
C
C
C      ***** TEMPERATURE AND PASSIVE SCALARS ******
C
C     Check if any pre-solv necessary for temperature/passive scalars
C
      IFANYP = .FALSE.
      DO 100 IFIELD=2,NFIELD
         IF (IFPRSL(IFIELD,jp)) THEN
            IF (NID.EQ.0) WRITE(6,101) IFIELD
            IFANYP = .TRUE.
         ENDIF
  100 CONTINUE
  101 FORMAT(2X,'Using PRESOLVE option for field',I2,'.')
C
C
C     If any pre-solv, do pre-solv for all temperatur/passive scalar fields
C
C
C     Fortran function initial conditions for temp/pass. scalars.
C
      MAXFLD = NFIELD
      IF (IFMODEL.AND.IFKEPS) MAXFLD = NFIELD-2
c
      jp = 0
      do 200 ifield=2,maxfld
         if (nid.eq.0) write(6,*) 'nekuic::',ifield,iffort(ifield,jp)
         if (iffort(ifield,jp)) call nekuic
 200  continue
c
      if (ifpert) then
         ! call exitt
         ifield=2
         do jp=1,npert
         if (nid.eq.0) write(6,*) 'nekuicP',ifield,jp,iffort(ifield,jp)
            if (iffort(ifield,jp)) call nekuic
         enddo
      endif

      jp = 0
      
C     Restart files
C
 
      call restart(nfiles)
C
C
C      ***** VELOCITY ******
C
C     (If restarting for V, we're done,
C     ...else, do pre-solv for fluid if requested.)
C
      IFIELD = 1
C
C
C     Fortran function initial conditions for velocity.
C
      if (nid.eq.0) write(6,*) 'call nekuic for vel:',iffort(1,jp),jp
      ifield = 1
      if (iffort(ifield,jp)) call nekuic
c
      if (ifpert) then
          call exitt
         ifield=1
         do jp=1,npert
            if (iffort(ifield,jp)) call nekuic
            if (nid.eq.0) write(6,*) 'ic vel pert:',iffort(1,jp),jp
         enddo
      endif
      jp = 0
c
      ntotv = nx1*ny1*nz1*nelv
      vxmax = glamax(vx,ntotv)
      vymax = glamax(vy,ntotv)
C
C     Fortran function initial conditions for turbulence k-e model
C
      if (ifmodel .and. ifkeps) then
         mfldt = nfield - 1
         do 300 ifield=mfldt,nfield
            if (iffort(ifield,jp)) call nekuic
 300     continue
      endif
C
C     Initial mesh velocities
C
C
C     Compute additional initial values for turbulence model arrays
C     based on I.C.
C
C
C     If convection-diffusion of a passive scalar with a fixed velocity field,
C     make sure to fill up lagged arrays since this will not be done in
C     the time-stepping procedure (no flow calculation) (01/18/91 -EMR).
C
      IF (.NOT.IFFLOW.AND.IFHEAT) THEN
         ITEST=0
         DO 400 IFIELD=2,NFIELD
            IF (IFADVC(IFIELD)) ITEST=1
 400     CONTINUE
         IF (ITEST.EQ.1) THEN
            NBDMAX = 3
            NBDSAV = NBDINP
            NBDINP = NBDMAX
            DO 500 I=1,NBDMAX
c              CALL LAGVEL
 500        CONTINUE
            NBDINP = NBDSAV
         ENDIF
      ENDIF
C     
C     Ensure that all processors have the same time as node 0.
C
      IF (NID.NE.0) TIME=0.0
      TIME=GLSUM(TIME,1)
      NTDUMP=0
      IF (TIMEIO.NE.0.0) NTDUMP = INT( TIME/TIMEIO )
C
C     Ensure that initial field is continuous!
C
      NXYZ1=NX1*NY1*NZ1
      NTOTT=NELT*NXYZ1
      NTOTV=NELV*NXYZ1
      NTOTG=NELGV*NXYZ1
c
C     first.. a test...
      ifield = 2
      if (ifflow) ifield = 1
      call rone(work,ntotv)
      ifield = 1
c                  i've no idea why this next stmt is here... pff 6/28/99
      flag_gs_init = 1
      CALL DSSUM(work,NX1,NY1,NZ1)
c     flag_gs_init = 0
      CALL COL2(work,VMULT,NTOTV)
      rtot  = glsum(work,ntotv)
      rtotv = ntotg
      rdif  = (rtot-rtotv)/rtotv
      if (nid.eq.0) write(6,51) ntotg,rtot,rdif
      if (nid.eq.0) write(6,52) rtotv,rtot,rdif
   51 format('dssum test1:',i12,3g13.5)
   52 format('dssum testf:',3g13.5)
c
      ntot2 = nx2*ny2*nz2*nelv
      vxmax = glamax(vx,ntotv)
      vymax = glamax(vy,ntotv)
      vzmax = glamax(vz,ntotv)
      prmax = glamax(pr,ntot2)
      ttmax = glamax(t ,ntott)
      if (nid.eq.0) write(6,6) vxmax,vymax,vzmax,prmax,ttmax
    6 format('uvwpt max:',5g13.5)
c
      small=1.0E-15
      if (vxmax.eq.0) call perturb(vx,'v',small)
      if (vymax.eq.0) call perturb(vy,'v',small)
      if (vzmax.eq.0) call perturb(vz,'v',small)
      if (prmax.eq.0) call perturb(pr,'p',small)
      if (ttmax.eq.0) call perturb(t ,'t',small)
c
      if (nid.eq.0) write(6,*) 'npscal',npscal
      do i=1,npscal
         psmax = glamax(t(1,1,1,1,1+i),ntott)
         if (psmax.eq.0) call perturb(t(1,1,1,1,1+i),'t',small)
      enddo
c
      if (ifflow) then
         ifield = 1
         call dssum(vx,nx1,ny1,nz1)
         call dssum(vy,nx1,ny1,nz1)
         call dssum(vz,nx1,ny1,nz1)
         call col2 (vx,vmult,ntotv)
         call col2 (vy,vmult,ntotv)
         call col2 (vz,vmult,ntotv)
      endif
c
      if (ifheat) then
         ifield = 2
         call dssum(t ,nx1,ny1,nz1)
         call col2 (t ,tmult,ntott)
         do i=1,npscal
            call dssum(t(1,1,1,1,1+i),nx1,ny1,nz1)
            call col2 (t(1,1,1,1,1+i),tmult,ntott)
         enddo
      endif
c
      if (ifpert) then
         do jp=1,npert
            ifield = 1
            call dssum(vxp(1,jp),nx1,ny1,nz1)
            call dssum(vyp(1,jp),nx1,ny1,nz1)
            call dssum(vzp(1,jp),nx1,ny1,nz1)
            call col2 (vxp(1,jp),vmult,ntotv)
            call col2 (vyp(1,jp),vmult,ntotv)
            call col2 (vzp(1,jp),vmult,ntotv)
            ifield = 2
            call dssum(tp(1,1,jp),nx1,ny1,nz1)
            call col2 (tp(1,1,jp),tmult,ntotv)
c           note... must be updated for addl pass. scal's. pff 4/26/04
            vxmax = glamax(vxp(1,jp),ntotv)
            vymax = glamax(vyp(1,jp),ntotv)
            if (nid.eq.0) write(6,111) jp,vxmax,vymax
  111       format(i5,1p2e12.4,' max pert vel')
         enddo
      endif
      jp = 0
c
c
      vxmax = glamax(vx,ntotv)
      vymax = glamax(vy,ntotv)
      vzmax = glamax(vz,ntotv)
      prmax = glamax(pr,ntot2)
      ttmax = glamax(t ,ntott)
      if (nid.eq.0) write(6,16) vxmax,vymax,vzmax,prmax,ttmax
   16 format('uvwpt max2',5g13.5)
c
      xxmin = glmin(xm1,ntotv)
      yymin = glmin(ym1,ntotv)
      zzmin = glmin(zm1,ntotv)
      xxmax = glmax(xm1,ntotv)
      yymax = glmax(ym1,ntotv)
      zzmax = glmax(zm1,ntotv)
      if (nid.eq.0) write(6,7) xxmin,yymin,zzmin,xxmax,yymax,zzmax
    7 format('xyz minmx:',6g13.5)
c
      if (ifrest(0,jp)) then
c        mesh has been read in.  recompute geometric factors
         if (nid.eq.0) write(6,*) 'New geometry -- remapping.'
c        call dsavg(xm1) ! doesn't work for periodic !
c        call dsavg(ym1) ! doesn't work for periodic !
c        call dsavg(zm1) ! doesn't work for periodic !
         call geom_reset(1)
      endif
c
      return
      end
C            
c-----------------------------------------------------------------------
      subroutine slogic (iffort,ifrest,ifprsl,nfiles)
C---------------------------------------------------------------------
C
C     Set up logicals for initial conditions.
C
C---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
c
      logical  iffort(  ldimt1,0:lpert)
     $       , ifrest(0:ldimt1,0:lpert)
     $       , ifprsl(  ldimt1,0:lpert)
c
      logical ifgetx,ifgetz,ifgetu,ifgetw,ifgetp,ifgett,ifgtps(ldimt1)
     $       ,ifgtim
      character*80 line,fname,cdum
      character*2  s2
      character*1  line1(80)
      equivalence (line1,line)
C
C     Default is user specified fortran function (=0 if not specified)
C
      nfldt = nfield

      do jp=0,npert
         ifrest(0,jp) = .false.
         do ifld=1,nfldt
            iffort(ifld,jp) = .true.
            ifrest(ifld,jp) = .false.
            ifprsl(ifld,jp) = .false.
         enddo
      enddo

      jp = 0
      nfiles=0
C
C     Check for Presolve options     
C
      DO 1000 ILINE=1,15 
         LINE=INITC(ILINE)
         CALL CAPIT(LINE,80)
         IF (INDX1(LINE,'PRESOLV',7).NE.0) THEN
C           found a presolve request
            CALL BLANK(INITC(ILINE),80)
            CALL LJUST(LINE)
            CALL CSPLIT(CDUM,LINE,' ',1)
C
            IF (LTRUNC(LINE,80).EQ.0) THEN
               IF (NID.EQ.0) WRITE(6,700)
  700          FORMAT(/,2X,'Presolve options: ALL')
C              default - all fields are presolved.
               DO 800 IFIELD=1,nfldt
                  ifprsl(ifield,jp) = .true.
                  iffort(ifield,jp) = .false.
  800          CONTINUE
            ELSE
C           check line for arguments
C
               LL=LTRUNC(LINE,80)
               IF (NID.EQ.0) WRITE(6,810) (LINE1(L),L=1,LL)
  810          FORMAT(/,2X,'Presolve options: ',80A1)
C
               IF (INDX2(LINE,'U',1).NE.0) THEN
                  ifprsl(1,jp) = .true.
                  iffort(1,jp) = .false.
               ENDIF
C
               IF (INDX2(LINE,'T',1).NE.0) THEN
                  ifprsl(2,jp) = .true.
                  iffort(2,jp) = .false.
               ENDIF
C
               DO 900 IFIELD=3,NPSCAL+2
                  IP=IFIELD-2
                  WRITE(S2,901) IP
                  IF (INDX2(LINE,S2,2).NE.0) THEN
                     ifprsl(ifield,jp) = .true.
                     iffort(ifield,jp) = .false.
                  ENDIF
  900          CONTINUE
  901          FORMAT('P',I1)
            ENDIF
         ENDIF
 1000    CONTINUE
C
C     Check for restart options
C
      jp = 0
      iline=1 
      DO 2000 ILINE=1,15
c        if (ifpert)  jp=iline-1  ! FIXME misun 6/4/07        
         LINE=INITC(ILINE)
         IF (LTRUNC(LINE,80).NE.0) THEN
C           found a filename
            NFILES=NFILES+1
            INITC(NFILES)=LINE
C
            IF (NID.EQ.0.AND.NFILES.EQ.1) WRITE(6,1010)
 1010       FORMAT(/,2X,'Restart options:')
            IF (NID.EQ.0) WRITE(6,'(A80)') LINE
C
C           Parse restart options 
            call sioflag(ifgetx,ifgetz,ifgetu,ifgetw,ifgetp
     $                  ,ifgett,ifgtps,ifgtim,ndumps,fname,line)
            IF (IFGETX) THEN
               IFREST(0,jp) = .TRUE.
            ENDIF
            IF (IFGETU) THEN
               iffort(1,jp) = .false.
               ifprsl(1,jp) = .false.
               ifrest(1,jp) = .true.
            ENDIF
            IF (IFGETT) THEN
               iffort(2,jp) = .false.
               ifprsl(2,jp) = .false.
               ifrest(2,jp) = .true.
            ENDIF
            DO 1900 IFIELD=3,nfldt
c              write(6,*) 'ifgetps:',(ifgtps(k),k=1,nfield)
               IF (IFGTPS(IFIELD-2)) THEN
                  iffort(ifield,jp) = .false.
                  ifprsl(ifield,jp) = .false.
                  ifrest(ifield,jp) = .true.
               ENDIF
 1900       CONTINUE
         ENDIF
 2000 CONTINUE
C
      return
      end

c-----------------------------------------------------------------------
      subroutine restart(nfiles)
C----------------------------------------------------------------------
C
C     (1) Open restart file(s)
C     (2) Check previous spatial discretization 
C     (3) Map (K1,N1) => (K2,N2) if necessary
C
C     nfiles > 1 has several implications:
C
C     i.   For std. run, data is taken from last file in list, unless
C          explicitly specified in argument list of filename
C
C     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
C          subsequent files are for B-field or perturbation fields
C
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LXYZT=LX1*LY1*LZ1*LELT)
      PARAMETER (LPSC9=LDIMT+9)

      common /poiss/  phi (ltot),rhs (ltot)
     $             ,  h1  (ltot),h2  (ltot)
     $             ,  mask(ltot),mult(ltot)
     $             ,  binv(ltot)
      real    mask,mult

c
      COMMON /SCRNS/ SDUMP(LXYZT,7)
      integer mesg(40)
c
c
      COMMON /SCRUZ/ SDMP2(LXYZT,4)
      COMMON /SCRMG/ SDMP3(LXYZT,4)
c
C     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 9.
      COMMON /CTMP1/ TDUMP(LXYZR,LPSC9)
      real*4         tdump
      character*4    adump(lxyzr,lpsc9)
      equivalence   (tdump,adump)
c
c
c     cdump comes in via PARALLEL (->TOTAL)
c     COMMON /CBLELG/ CDUMP(LELG)
c
      CHARACTER*20 FRMAT
      CHARACTER*2  EXCODER(10)
      CHARACTER*80 FNAME
      CHARACTER*1  FNAME1(80)
      EQUIVALENCE (FNAME1,FNAME)
C
      INTEGER      HNAMI (20)
      CHARACTER*80 HNAME
      CHARACTER*1  HNAME1(80)
      EQUIVALENCE (HNAME,HNAME1)
      EQUIVALENCE (HNAME,HNAMI )

      logical  ifpert !FIXME misun 6/4/07

C
C     Local logical flags to determine whether to copy data or not.
C
      logical ifgetx,ifgetz,ifgetu,ifgetw,ifgetp,ifgett,ifgtps(ldimt1)
     $       ,ifgtim,ifok,iffmat
      integer iposx,iposz,iposu,iposw,iposp,ipost,ipsps(ldimt1)
C
      logical ifbytsw, if_byte_swap_test
      real*4   bytetest

      common /VPTS1/ vxp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              ,vyp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              ,vzp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              ,prp(lpx1*lpy1*lpz1*lpelv,lpert)
     $              , tp(lpx1*lpy1*lpz1*lpelt,ldimt,lpert)
      common /VPTS2/  pm(lbx2,lby2,lbz2,lbelv)
c
      ifpert = .false.
      ifok=.false.
      iffmat=.true.
      frmat='NORMAL'
      if (param(32).eq.1.0) iffmat=.false.
      IF (PARAM(32).EQ.1.0) FRMAT='UNFORMATTED'
      IF (PARAM(32).EQ.2.0) FRMAT='UNIVERSAL'
      if (param(32).GE.3.0) iffmat=.false.
      IF (PARAM(32).ge.3.0) FRMAT='UNFORMAT2'
      ifbytsw = .false.

      if (param(32).eq.6.0) then
         do ifile=1,nfiles
            call  mlti_file_input(ifile)
         enddo
         return
      endif
c
      DO 6000 IFILE=1,NFILES
C
C        Loop over the requested restart files
C
         call sioflag(ifgetx,ifgetz,ifgetu,ifgetw,ifgetp,ifgett,
     $                ifgtps,ifgtim,ndumps,fname,initc(ifile))
C
C
C       Standard Case:   Requested inputs are in the file.
C
        if (nid.eq.0) then
          if (iffmat) then
            open (unit=91,file=fname,status='old',err=500)
          elseif (param(32).eq.1) then
            open (unit=91,file=fname
     $           ,form='unformatted',status='old',err=500)
          elseif (param(32).eq.3) then
c
            len= ltrunc(fname,79)
            call izero (hnami,20)
            call chcopy(hname1,fname,len)
c           test for presence of .fld file
            OPEN (UNIT=92,FILE=HNAME
     $           ,FORM='FORMATTED',STATUS='OLD',ERR=500)
            close(UNIT=92)
            call byte_open(hname)
c
c           Open header file
c
            ihd= indx1 (fname,'.fld',4)
            call chcopy(hname1(ihd),'.fhd',4)
            open (unit=92,file=hname
     $           ,form='formatted',status='old',err=500)
          elseif (param(32).eq.4) then
c
            len= ltrunc(fname,79)
            call izero (hnami,20)
            call chcopy(hname1,fname,len)
c           test for presence of file
            open (unit=91,file=hname
     $           ,form='unformatted',status='old',err=500)
            close(UNIT=91)
            call byte_open(hname)
c
          endif 
          ifok = .true.
        endif
  500   continue

        call lbcast(ifok)

        if (.not.ifok) goto 5000
         
C
c        IF (NDUMPS.EQ.0) NDUMPS=10000
c
         ndumps = 1
C
C        Only NODE 0 reads from the disk.
C
         DO 1000 IDUMP=1,NDUMPS
C
            IF (NID.EQ.0) THEN
               if (iffmat) then

                 if (nelgt.lt.10000) then
                    read(91,91,err=1500,end=1500)
     $              neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10)
   91               format(4i4,1x,g13.4,i5,1x,10a2)
                 else
                    read(91,92,err=1500,end=1500)
     $              neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10)
   92               format(i10,3i4,1P1e18.9,i9,1x,15a2)
                 endif

               elseif (param(32).eq.3) then
                 READ(92,'(4I4,1X,G13.4,I5,1X,10A2)',ERR=1500,END=1500)
     $           NELTR,NXR,NYR,NZR,RSTIME,ISTEPR,(EXCODER(I),I=1,10)
                 call byte_read(bytetest,1)
                 ifbytsw = if_byte_swap_test(bytetest)
               elseif (param(32).eq.4) then
                 call byte_read(hnami,20)
                 if (nelgt.lt.10000) then
                    read(hname,'(4i4,1x,g13.4,i5,1x,10a2)'
     $              ,err=1500,end=1500)
     $               neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10)
                 else
                    read(hname,'(i10,3i4,1P1e18.9,i9,1x,15a2)'
     $              ,err=1500,end=1500)
     $               neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10)
                 endif
                 call byte_read(bytetest,1)
                 ifbytsw = if_byte_swap_test(bytetest)
               ELSE
                 READ(91,ERR=1500,END=1500)
     $           NELTR,NXR,NYR,NZR,RSTIME,ISTEPR,(EXCODER(I),I=1,10)
               ENDIF
               mesg(1) = neltr
               mesg(2) = nxr
               mesg(3) = nyr
               mesg(4) = nzr
               write(6,333) nelgt,neltr,nxr,nyr,nzr,' NELR',param(32)
  333          format(2i9,3i4,a5,1p1e12.4)
               call chcopy(mesg(5),excoder,20)
               len  = 14*isize
            endif
c
            IF (IDUMP.EQ.1) THEN
               len  = 14*isize
               call bcast(mesg,len)
               neltr = mesg(1)
               nxr   = mesg(2)
               nyr   = mesg(3)
               nzr   = mesg(4)
               call   chcopy(excoder,mesg(5),20)
c              write(6,*) nid,' excoder: ',(excoder(k),k=1,10)
c
               call lbcast(ifbytsw)
C
C              Find out the position of the data via the header info.
C
C              Assumptions about position of data:
C
C              .If X is found then Y (Z) are implied
C              .If U is found then V (W) are implied
C              .Column data is space or , delimited.
C
               IF (NELTR.NE.NELGT.AND.NID.EQ.0) THEN
                  PRINT*,'ERROR: RESTART FILE MUST HAVE SAME'
                  PRINT*,'NUMBER OF ELEMENTS AS NEW SESSION'
                  PRINT*,' OLD ',NELTR,' NEW ',NELGT
                  call exitt               
               ENDIF
C
C              Bounds checking on mapped data.
               IF (NXR.GT.LXR) THEN
                  WRITE(6,20) NXR,NX1
   20             FORMAT(//,2X,
     $            'WARNING:  Attempt to map from',I3,
     $            ' to N=',I3,'.',/,2X,
     $            'NEKTON currently supports mapping from N+2 or less.'
     $            ,/,2X,'Increase N or LXR in IC.FOR.')
                  CALL EMERXIT
               ENDIF
C
C
C              Figure out position of data in file "IFILE"
C
               NOUTS=0
               IPOSX=0
               IPOSY=0
               IPOSZ=0
               IPOSU=0
               IPOSV=0
               IPOSW=0
               IPOSP=0
               IPOST=0
               DO 40 I=1,NPSCAL
                  IPSPS(I)=0
   40          CONTINUE
C
               DO 50 I=1,10
                  IF (EXCODER(I).EQ.'X') THEN
                     NOUTS=NOUTS + 1
                     IPOSX=NOUTS
                     NOUTS=NOUTS+1
                     IPOSY=NOUTS
                     IF (IF3D) THEN
                        NOUTS=NOUTS + 1
                        IPOSZ=NOUTS
                     ENDIF
                  ENDIF
                  IF (EXCODER(I).EQ.'U') THEN
                     NOUTS=NOUTS + 1
                     IPOSU=NOUTS
                     NOUTS=NOUTS+1
                     IPOSV=NOUTS
                     IF (IF3D) THEN
                        NOUTS=NOUTS + 1
                        IPOSW=NOUTS
                     ENDIF
                  ENDIF
                  IF (EXCODER(I).EQ.'P') THEN
                     NOUTS=NOUTS + 1
                     IPOSP=NOUTS
                  ENDIF
                  IF (EXCODER(I).EQ.'T') THEN
                     NOUTS=NOUTS + 1
                     IPOST=NOUTS
                  ENDIF
                  IF (EXCODER(I).EQ.'1') THEN
                     NOUTS=NOUTS + 1
                     IPSPS(1)=NOUTS
                  ENDIF
                  IF (EXCODER(I).EQ.'2') THEN
                     NOUTS=NOUTS + 1
                     IPSPS(2)=NOUTS
                  ENDIF
                  IF (EXCODER(I).EQ.'3') THEN
                     NOUTS=NOUTS + 1
                     IPSPS(3)=NOUTS
                  ENDIF
                  IF (EXCODER(I).EQ.'4') THEN
                     NOUTS=NOUTS + 1
                     IPSPS(4)=NOUTS
                  ENDIF
   50          CONTINUE
C
               LNAME=LTRUNC(FNAME,80)
               IF (NID.EQ.0) WRITE(6,61) (FNAME1(I),I=1,LNAME)
               IF (NID.EQ.0) WRITE(6,62) 
     $             IPOSU,IPOSV,IPOSW,IPOSP,IPOST,NOUTS
   61          FORMAT(/,2X,'Restarting from file ',80A1)
   62          FORMAT(2X,'Columns for restart data, U,V,W,P,T,N: ',6I4)
C
C              Make sure the requested data is present in this file....
               IF (IPOSX.EQ.0) IFGETX=.FALSE.
               IF (IPOSY.EQ.0) IFGETX=.FALSE.
               IF (IPOSZ.EQ.0) IFGETZ=.FALSE.
               IF (IPOSU.EQ.0) IFGETU=.FALSE.
               IF (IPOSV.EQ.0) IFGETU=.FALSE.
               IF (IPOSW.EQ.0) IFGETW=.FALSE.
               IF (IPOSP.EQ.0) IFGETP=.FALSE.
               IF (IPOST.EQ.0) IFGETT=.FALSE.
               DO 65 I=2,NPSCAL
                  IF (IPSPS(I).EQ.0) IFGTPS(I)=.FALSE.
   65          CONTINUE
C
C              End of restart file header evaluation.
C
            ENDIF
C
C           Read the error estimators
C
            IF(NID.EQ.0)THEN
               if (iffmat) then
                  READ(91,'(6G11.4)',END=1500)(CDUMP(I),I=1,NELTR)
               ELSEIF (FRMAT.eq.'UNFORMATTED') then
                  READ(91,END=1500)(CDUMP(I),I=1,NELTR)
c              ELSEIF (FRMAT.eq.'UNFORMAT2') then
c                 No header files! :)
               ENDIF
            ENDIF
C
C           Read the current dump, double buffer so that we can
C           fit the data on a distributed memory machine,
C           and so we won't have to read the restart file twice
C           in case of an incomplete data file.
C
            NXYZR = NXR*NYR*NZR
C
C           Read the data
C
            DO 200 IEG=1,NELGT
               ifok = .false.
               IF (NID.EQ.0) THEN
                 IF (MOD(IEG,100).EQ.1) WRITE(6,*) 'Reading',IEG
                 IF (FRMAT.EQ.'NORMAL') THEN
                    READ(91,*,ERR=1500,END=1500)
     $              ((tdump(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
                 ELSEIF (FRMAT.EQ.'UNFORMATTED') THEN
                    READ(91,ERR=1500,END=1500)
     $              ((tdump(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
                 ELSEIF (FRMAT.EQ.'UNFORMAT2') THEN
                    do ii=1,nouts
                       call byte_read(tdump(1,II),nxyzr)
                    enddo
                 ELSEIF (FRMAT.EQ.'UNIVERSAL') THEN
c                   fixed for real*8 conversion... pff 5/21/96
                    READ(91,201,ERR=1500,END=1500)
     $              ((adump(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
  201               FORMAT(20A4)
                 ENDIF
                 IFOK=.TRUE.
               ENDIF
C
C              Notify other processors that we've read the data OK.
C
               call lbcast(ifok)
               IF (.NOT.IFOK) GOTO 1600
C
C              MAPDMP maps data from NXR to NX1
C              (and sends data to the appropriate processor.)
C
C              The buffer SDUMP is used so that if an incomplete dump
C              file is found (e.g. due to UNIX io buffering!), then
C              the previous read data stored in VX,VY,.., is not corrupted.
C
               CALL MAPDMP
     $         (SDUMP(1,1),TDUMP(1,IPOSX),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDUMP(1,2),TDUMP(1,IPOSY),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDUMP(1,3),TDUMP(1,IPOSZ),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDUMP(1,4),TDUMP(1,IPOSU),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDUMP(1,5),TDUMP(1,IPOSV),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDUMP(1,6),TDUMP(1,IPOSW),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDUMP(1,7),TDUMP(1,IPOSP),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
               CALL MAPDMP
     $         (SDMP2(1,1),TDUMP(1,IPOST),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
C              passive scalars
               NP3=MIN(NPSCAL,3)
               DO 100 IPS=1,NP3
                  IF (IFGTPS(IPS)) CALL MAPDMP(SDMP2(1,IPS+1)
     $               ,TDUMP(1,IPSPS(IPS)),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
  100          CONTINUE               
               NP4=NPSCAL-3
               NP4=MIN(NP4,4)
               DO 110 IPS=1,NP4
                  IP4=IPS+3
                  IF (IFGTPS(IP4)) CALL MAPDMP(SDMP3(1,IPS)
     $               ,TDUMP(1,IPSPS(IP4)),IEG,NXR,NYR,NZR,FRMAT,ifbytsw)
  110          CONTINUE               
  200       CONTINUE               
C
C           Successfully read a complete field, store it.
C
            NXYZ2=NX2*NY2*NZ2
            NXYZ1=NX1*NY1*NZ1
            NTOTT=NELT*NXYZ1
            NTOTV=NELV*NXYZ1
c
c           if (ifpert.and.ifile.ge.2) then  ! FIXME this is original line:6/4/07 misun

             if (ifpert) then
               j=ifile-1  ! pointer to perturbation field
               if (ifgetu) call copy(vxp(1,j),sdump(1,4),ntotv)
               if (ifgetu) call copy(vyp(1,j),sdump(1,5),ntotv)
               if (ifgetw) call copy(vzp(1,j),sdump(1,6),ntotv)
               if (ifgetp) then
                  if (nid.eq.0) write(6,*) 'getting restart pressure'
                  if (ifsplit) then
                     call copy(prp(1,j),sdump(1,7),ntotv)
                  else
                     do ie=1,nelv
                        ie1 = (ie-1)*nxyz1+1
                        ie2 = (ie-1)*nxyz2+1
                        call map12 (prp(ie2,j),sdump(ie1,7),ie)
                     enddo
                  endif
               endif
               if (ifgett) call copy(tp(1,1,j),sdmp2(1,1),ntott)
C              passive scalars
               np3=min(3,npscal)
               do ips=1,np3
                  if (ifgtps(ips))
     $            call copy(tp(1,ips+1,j),sdmp2(1,ips+1),ntott)
               enddo
               np4=npscal-3
               np4=min(np4,4)
               do ips=1,np4
                  ip4=ips+3
                  if (ifgtps(ip4))
     $            call copy(tp(1,ip4+1,j),sdmp3(1,ips),ntott)
               enddo
c
            else  ! Std. Case
               call copy(xm1,sdump(1,1),ntott)
               call copy(ym1,sdump(1,2),ntott)
               call copy(zm1,sdump(1,3),ntott)
               call copy(vx ,sdump(1,4),ntotv)
               call copy(vy ,sdump(1,5),ntotv)
               call copy(vz ,sdump(1,6),ntotv)
               call copy(ex,vx,ntotv)
               call copy(ey,vy,ntotv)
               call copy(ez,vz,ntotv)
               if (ifgetp) then
                  if (nid.eq.0) write(6,*) 'getting restart pressure'
                  if (ifsplit) then
                     call copy(pr,sdump(1,7),ntotv)
                  else
                     do iel=1,nelv
                        iiel = (iel-1)*nxyz1+1
                        call map12 (pr(1,1,1,iel),sdump(iiel,7),iel)
                     enddo
                  endif
               endif

               call copy(t,sdmp2(1,1),ntott)
               call copy(phi,t,ntott)
C              passive scalars
               np3=min(3,npscal)
               do ips=1,np3
                  if (ifgtps(ips))
     $            call copy(t(1,1,1,1,ips+1),sdmp2(1,ips+1),ntott)
               enddo
               np4=npscal-3
               np4=min(np4,4)
               do ips=1,np4
                  ip4=ips+3
                  if (ifgtps(ip4))
     $            call copy(t(1,1,1,1,ip4+1),sdmp3(1,ips),ntott)
               enddo
c
               if (ifgtim) time=rstime
c
            endif
 1000    CONTINUE
         GOTO 1600
C
C           Else, end of file found - notify other processors.
 1500       CONTINUE
            IFOK=.FALSE.
            call lbcast(ifok)
C
 1600    CONTINUE
C
         IF (IDUMP.EQ.1.AND.NID.EQ.0) THEN
            WRITE(6,1700) FNAME
            WRITE(6,1701) IEG,IXYZ
            WRITE(6,1702) 
     $            ((TDUMP(JXYZ,II),II=1,NOUTS),JXYZ=IXYZ-1,IXYZ)
 1700       FORMAT(5X,'WARNING:  No data read in for file ',A80)
 1701       FORMAT(5X,'Failed on  element',I4,',  point',I5,'.')
 1702       FORMAT(5X,'Last read dump:',/,5G15.7)
            write(6,*) nid,'call exitt 1702a',idump
            call exitt
         ELSEIF (IDUMP.EQ.1) THEN
            write(6,*) nid,'call exitt 1702b',idump
            call exitt
         ELSE
            IDUMP=IDUMP-1
            IF (NID.EQ.0) WRITE(6,1800) IDUMP
 1800       FORMAT(2X,'Successfully read data from dump number',I3,'.')
         ENDIF
         if (frmat.eq.'UNFORMAT2') then
            if (nid.eq.0.and.param(32).eq.3.) close(unit=92)
            if (nid.eq.0) call byte_close()
         else
            if (nid.eq.0) close(unit=91)
         endif
         GOTO 6000
C
C        Can't open file...
 5000    CONTINUE
         IF (NID.EQ.0) WRITE(6,5001) FNAME 
 5001    FORMAT(2X,'   *******   ERROR   *******    '
     $       ,/,2X,'   *******   ERROR   *******    '
     $       ,/,2X,'   Could not open restart file:'
     $       ,/,A80
     $      ,//,2X,'Quitting in routine RESTART.')
         CLOSE(UNIT=91)
         if (nid.eq.0.and.param(32).eq.3.) close(unit=92)
         call exitt
 5002    CONTINUE
         IF (NID.EQ.0) WRITE(6,5001) HNAME 
         call exitt
C
C
C     End of IFILE loop
 6000 CONTINUE
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine sioflag(ifgetx,ifgetz,ifgetu,ifgetw,ifgetp,ifgett,
     $                   ifgtps,ifgtim,ndumps,fname,rsopts)
C
C     Set IO flags according to Restart Options File, RSOPTS
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
C
      CHARACTER*80 RSOPTS,FNAME
      CHARACTER*2  S2
      LOGICAL IFGETX,IFGETZ,IFGETU,IFGETW,IFGETP,IFGETT,IFGTPS(LDIMT1)
     $       ,IFGTIM
      LOGICAL IFGTRL
C
C     Scratch variables..
      LOGICAL IFDEFT,IFANYC
      CHARACTER*80 RSOPT     ,LINE
      CHARACTER*1  RSOPT1(80),LINE1(80)
      EQUIVALENCE (RSOPT1,RSOPT)
      EQUIVALENCE (LINE1,LINE)
C
C     Parse filename
C
C        CSPLIT splits S1 into two parts, delimited by S2.  
C        S1 returns with 2nd part of S1.  CSPLIT returns 1st part.
C
      RSOPT=RSOPTS
      CALL LJUST(RSOPT)
      CALL CSPLIT(FNAME,RSOPT,' ',1)
C     check fname for user supplied extension.
      IF (INDX1(FNAME,'.',1).EQ.0) THEN
         LEN=LTRUNC(FNAME,80)
         LEN1=LEN+1
         LEN4=LEN+4
         FNAME(LEN1:LEN4)='.fld'
      ENDIF
C
C     Parse restart options
C
C        set default flags
C
      IFGETX=.FALSE.
      IFGETZ=.FALSE.
      IFGETU=.FALSE.
      IFGETW=.FALSE.
      IFGETP=.FALSE.
      IFGETT=.FALSE.
      DO 100 I=1,LDIMT1
         IFGTPS(I)=.FALSE.
  100 CONTINUE
      IFGTIM=.TRUE.
      NDUMPS=0
C
C     Check for default case - just a filename given, no i/o options specified
C
      IFDEFT=.TRUE.
C
C     Parse file for i/o options and/or dump number
C
      CALL CAPIT(RSOPT,80)
      IF (LTRUNC(RSOPT,80).NE.0) THEN
C
C        Check for explicit specification of restart TIME.
C
         ITO=INDX2(RSOPT,'TIME',4)
         IFGTIM=.TRUE.
         IF (ITO.NE.0) THEN
C           user has specified the time explicitly.
            IT1=INDX2(RSOPT,'=',1)
            IT8=80-IT1
            CALL BLANK(LINE,80)
            CALL CHCOPY(LINE,RSOPT1(IT1),IT8)
            IF (IFGTRL(TTIME,LINE)) THEN
               IFGTIM=.FALSE.
               TIME=TTIME
            ENDIF
C           remove the user specified time from the RS options line.
            ITA=80-ITO+1
            CALL BLANK(RSOPT1(ITO),ITA)
            CALL LJUST(LINE)
            IT1=INDX1(LINE,' ',1)
            ITB=80-IT1+1
            CALL CHCOPY(RSOPT1(ITO),LINE1(IT1),ITB)
         ENDIF
C
C        Parse field specifications.
C
         IXO=INDX2(RSOPT,'X',1)
         IF (IXO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETX=.TRUE.
            IF (IF3D) IFGETZ=.TRUE.
         ENDIF
C
         IVO=INDX2(RSOPT,'U',1)
         IF (IVO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETU=.TRUE.
            IF (IF3D) IFGETW=.TRUE.
         ENDIF
C
         IPO=INDX2(RSOPT,'P',1)
         IF (IPO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETP=.TRUE.
         ENDIF
C
         ITO=INDX2(RSOPT,'T',1)
         IF (ITO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETT=.TRUE.
         ENDIF
C
         DO 300 I=2,NPSCAL
            WRITE (S2,301) I
            IPO=INDX2(RSOPT,S2,2)
            IF (IPO.NE.0) THEN
               IFDEFT=.FALSE.
               IFGTPS(I)=.TRUE.
            ENDIF
  300    CONTINUE
  301    FORMAT('P',I1)
C
C        Get number of dumps from remainder of user supplied line.
C
         IF (IFGTRL(TDUMPS,RSOPT)) NDUMPS=INT(TDUMPS)
      ENDIF
C
C     If no fields were explicitly specified, assume getting all fields. 
C
      IF (IFDEFT) THEN
         IF (IFXYO) THEN
            IFGETX=.TRUE.
            IF (IF3D) IFGETZ=.TRUE.
         ENDIF
         IFANYC=.FALSE.
         DO 400 I=1,NFIELD
            IF (IFADVC(I)) IFANYC=.TRUE.
  400    CONTINUE
         IF (IFFLOW.OR.IFANYC) THEN
            IFGETU=.TRUE.
            IF (IF3D) IFGETW=.TRUE.
         ENDIF
         IF (IFFLOW) IFGETP=.TRUE.
         IF (IFHEAT) IFGETT=.TRUE.
         DO 410 I=1,NPSCAL
            IFGTPS(I)=.TRUE.
  410    CONTINUE
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine mapdmp(sdump,tdump,ieg,nxr,nyr,nzr,frmat,if_byte_sw)
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      CHARACTER*20 FRMAT
C
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
C
      REAL   SDUMP(LXYZ1,LELT)
      REAL*4 TDUMP(LXYZR)
c
      logical if_byte_sw
c
      NXYZ=NX1*NY1*NZ1
      NXYR=NXR*NYR*NZR
C
C     Serial processor code:
C
      IF (NP.EQ.1) THEN
C
         IF (FRMAT.EQ.'UNIVERSAL') CALL VRNVERT     (TDUMP,NXYR)
c        IF (if_byte_sw)           call vbyte_swap  (TDUMP,NXYR)
         IF (if_byte_sw)           CALL byte_reverse(TDUMP,NXYR)
         IF (NXR.EQ.NX1.AND.NYR.EQ.NY1.AND.NZR.EQ.NZ1) THEN
            CALL COPY4r(SDUMP(1,IEG),TDUMP,NXYZ)
         ELSE
C           do the map    (assumes that NX=NY=NZ, or NX=NY, NZ=1)
            CALL MAPAB4r(SDUMP(1,IEG),TDUMP,NXR,1)
         ENDIF
C
       ELSE
C
C     Parallel code - send data to appropriate processor and map.
C
         JNID=GLLNID(IEG)
         MTYPE=3333+IEG
         LEN=4*NXYR
         LE1=4
         IF (NID.EQ.0.AND.JNID.NE.0) THEN
c           hand-shake
            CALL CSEND(MTYPE,TDUMP,LE1,JNID,NULLPID)
            CALL CRECV(MTYPE,dummy,LE1)
            CALL CSEND(MTYPE,TDUMP,LEN,JNID,NULLPID)
         ELSEIF (NID.NE.0.AND.JNID.EQ.NID) THEN
C           Receive data from node 0
            CALL CRECV(MTYPE,dummy,LE1)
            CALL CSEND(MTYPE,TDUMP,LE1,0,NULLPID)
            CALL CRECV(MTYPE,TDUMP,LEN)
         ENDIF
C
C        If the data is targeted for this processor, then map 
C        to appropriate element.
C
         IF (JNID.EQ.NID) THEN
            IE=GLLEL(IEG)
            IF (FRMAT.EQ.'UNIVERSAL') CALL VRNVERT     (TDUMP,NXYR)
c           IF (if_byte_sw)           call vbyte_swap  (TDUMP,NXYR)
            IF (if_byte_sw)           CALL byte_reverse(TDUMP,NXYR)
            IF (NXR.EQ.NX1.AND.NYR.EQ.NY1.AND.NZR.EQ.NZ1) THEN
               CALL COPY4r(SDUMP(1,IE),TDUMP,NXYZ)
            ELSE
               CALL MAPAB4r(SDUMP(1,IE),TDUMP,NXR,1)
            ENDIF
         ENDIF
C
C        End of parallel distribution/map routine.
C
      ENDIF
      return
      END
c-----------------------------------------------------------------------
      subroutine mapab(x,y,nxr,nel)
C---------------------------------------------------------------
C
C     Interpolate Y(NXR,NYR,NZR,NEL) to X(NX1,NY1,NZ1,NEL)
C     (assumes that NXR=NYR=NZR, or NXR=NYR, NZR=1)
C---------------------------------------------------------------
C
      INCLUDE 'SIZE'
      INCLUDE 'IXYZ'
      INCLUDE 'WZ'
C
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      DIMENSION X(NX1,NY1,NZ1,NEL)
      DIMENSION Y(NXR,NXR,NXR,NEL)
C
      COMMON /CTMP0/ XA(LXYZR)      ,XB(LX1,LY1,LZR)
     $              ,IRES(LXR,LXR)  ,ITRES(LXR,LXR)
     $              ,ZGMR(20)       ,WGTR(20)
      REAL XA,XB,IRES,ITRES,ZGMR,WGTR
C
      INTEGER NOLD,ICALLD
      SAVE    NOLD,ICALLD
      DATA    NOLD,ICALLD /0,0/
C
      NZR = NXR
      IF(NZ1.EQ.1) NZR=1
      NYZR = NXR*NZR
      NXY1 = NX1*NY1
C
      IF (NXR.NE.NOLD) THEN
         NOLD=NXR
         CALL ZWGLL   (ZGMR,WGTR,NXR)
         CALL IGLLM   (IRES,ITRES,ZGMR,ZGM1,NXR,NX1,NXR,NX1)      
         IF (NID.EQ.0) WRITE(6,10) NXR,NX1
   10       FORMAT(2X,'Mapping restart data from Nold=',I2
     $               ,' to Nnew=',I2,'.')
      ENDIF
C
      DO 1000 IE=1,NEL
         CALL MXM (IRES,NX1,Y(1,1,1,IE),NXR,XA,NYZR)
         DO 100 IZ=1,NZR
            IZOFF = 1 + (IZ-1)*NX1*NXR
            CALL MXM (XA(IZOFF),NX1,ITRES,NXR,XB(1,1,IZ),NY1)
  100    CONTINUE
         IF (NDIM.EQ.3) THEN
            CALL MXM (XB,NXY1,ITRES,NZR,X(1,1,1,IE),NZ1)
         ELSE
            CALL COPY(X(1,1,1,IE),XB,NXY1)
         ENDIF
 1000 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine mapab4R(x,y,nxr,nel)
C---------------------------------------------------------------
C
C     Interpolate Y(NXR,NYR,NZR,NEL) to X(NX1,NY1,NZ1,NEL)
C     (assumes that NXR=NYR=NZR, or NXR=NYR, NZR=1)
c
c     Input:  real*4,  Output:  default precision
c
C---------------------------------------------------------------
C
      INCLUDE 'SIZE'
      INCLUDE 'IXYZ'
      INCLUDE 'WZ'
C
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      REAL*4 X(NX1,NY1,NZ1,NEL)
      REAL   Y(NXR,NXR,NXR,NEL)
C
      COMMON /CTMP0/ XA(LXYZR)      ,XB(LX1,LY1,LZR) ,XC(LXYZR)      
     $              ,IRES(LXR,LXR)  ,ITRES(LXR,LXR)
     $              ,ZGMR(20)       ,WGTR(20)
      REAL XA,XB,IRES,ITRES,ZGMR,WGTR
C
      INTEGER NOLD,ICALLD
      SAVE    NOLD,ICALLD
      DATA    NOLD,ICALLD /0,0/
C
      NZR = NXR
      IF(NZ1.EQ.1) NZR=1
      NYZR = NXR*NZR
      NXY1 = NX1*NY1
      nxyzr = nxr*nxr*nzr
C
      IF (NXR.NE.NOLD) THEN
         NOLD=NXR
         CALL ZWGLL   (ZGMR,WGTR,NXR)
         CALL IGLLM   (IRES,ITRES,ZGMR,ZGM1,NXR,NX1,NXR,NX1)      
         IF (NID.EQ.0) WRITE(6,10) NXR,NX1
   10       FORMAT(2X,'Mapping restart data from Nold=',I2
     $               ,' to Nnew=',I2,'.')
      ENDIF
C
      DO 1000 IE=1,NEL
         call copy4r(xc,y(1,1,1,ie),nxyzr)
         CALL MXM (IRES,NX1,xc,NXR,XA,NYZR)
         DO 100 IZ=1,NZR
            IZOFF = 1 + (IZ-1)*NX1*NXR
            CALL MXM (XA(IZOFF),NX1,ITRES,NXR,XB(1,1,IZ),NY1)
  100    CONTINUE
         IF (NDIM.EQ.3) THEN
            CALL MXM (XB,NXY1,ITRES,NZR,X(1,1,1,IE),NZ1)
         ELSE
            CALL COPY(X(1,1,1,IE),XB,NXY1)
         ENDIF
 1000 CONTINUE
C
      return
      END
C
c-----------------------------------------------------------------------
      subroutine nekuic
C------------------------------------------------------------------
C
C     User specified fortran function (=0 if not specified)
C
C------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'TURBO'
      INCLUDE 'PARALLEL'
      INCLUDE 'NEKUSE'
C
      NEL   = NELFLD(IFIELD)
C
      IF (IFMODEL .AND. IFKEPS .AND. IFIELD.EQ.IFLDK) THEN
C
      DO 100 IEL=1,NEL
         IEG = LGLEL(IEL,NODE)
         DO 100 K=1,NZ1
         DO 100 J=1,NY1
         DO 100 I=1,NX1
            CALL NEKASGN (I,J,K,IEL)
            CALL USERIC  (I,J,K,IEG)
            T(I,J,K,IEL,IFIELD-1) = TURBK
 100  CONTINUE
C
      ELSEIF (IFMODEL .AND. IFKEPS .AND. IFIELD.EQ.IFLDE) THEN
C
      DO 200 IEL=1,NEL
         IEG = LGLEL(IEL,NODE)
         DO 200 K=1,NZ1
         DO 200 J=1,NY1
         DO 200 I=1,NX1
            CALL NEKASGN (I,J,K,IEL)
            CALL USERIC  (I,J,K,IEG)
            T(I,J,K,IEL,IFIELD-1) = TURBE
 200  CONTINUE
C
      ELSE
C
      DO 300 IEL=1,NEL
         IEG = LGLEL(IEL,NODE)
         DO 300 K=1,NZ1
         DO 300 J=1,NY1
         DO 300 I=1,NX1
           CALL NEKASGN (I,J,K,IEL)
           CALL USERIC  (I,J,K,IEG)
           if (jp.eq.0) then
             IF (IFIELD.EQ.1) THEN
               VX(I,J,K,IEL) = UX
               VY(I,J,K,IEL) = UY
               VZ(I,J,K,IEL) = UZ
c              BX(I,J,K,IEL) = UX
c              BY(I,J,K,IEL) = UY
c              BZ(I,J,K,IEL) = UZ
             ELSE
               T(I,J,K,IEL,IFIELD-1) = TEMP
             ENDIF
           else
             ijke = i+nx1*((j-1)+ny1*((k-1) + nz1*(iel-1)))
             IF (IFIELD.EQ.1) THEN
c              VXP(IJKE,jp) = UX
c              VYP(IJKE,jp) = UY
c              VZP(IJKE,jp) = UZ
             ELSE
c              TP(IJKE,IFIELD-1,jp) = TEMP
             ENDIF
           endif

 300  CONTINUE
C
      ENDIF
c
      return
      END
c-----------------------------------------------------------------------
      LOGICAL FUNCTION IFGTRL(VALUE,LINE)
C
C     Read VALUE from LINE and set IFGTRL to .TRUE. if successful,
C                                  IFGTRL to .FALSE. otherwise.
C
C     This complicated function is necessary thanks to the Ardent,
C     which won't allow free formatted reads (*) from internal strings!
C
      CHARACTER*80 LINE
      CHARACTER*80 WORK
      CHARACTER*8  FMAT
C
C     Note that the format Fn.0 is appropriate for fields of type:
C          34   34.0  34.0e+00
C     The only difficulty would be with '34' but since we identify
C     the field width exactly, there is no problem.
C
      IFGTRL=.FALSE.
      VALUE=0.0
C
      WORK=LINE
      CALL LJUST(WORK)
      IFLDW=INDX1(WORK,' ',1)-1
C
      IF (IFLDW.GT.0) THEN
         WRITE(FMAT,10) IFLDW
   10    FORMAT('(F',I3.3,'.0)')
         READ(WORK,FMAT,ERR=100,END=100) TVAL
         VALUE=TVAL
         IFGTRL=.TRUE.
         return
      ENDIF
C
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine perturb(tt,cc,eps)
      include 'SIZE'
      include 'TOTAL'
c
      real tt(1)
      character*1 cc
c
      if (cc.eq.'v') then
         ifielt = ifield
         ifield = 1
         n = nx1*ny1*nz1*nelv
         call vcospf(tt,bm1,n)
         call cmult(tt,eps,n)
         call dssum(tt,nx1,ny1,nz1)
         ifield = ifielt
      else if (cc.eq.'t') then
         ifielt = ifield
         ifield = 2
         n = nx1*ny1*nz1*nelt
         call vcospf(tt,bm1,n)
         call cmult(tt,eps,n)
         call dssum(tt,nx1,ny1,nz1)
         ifield = ifielt
      else if (nx1.eq.nx2) then
         ifielt = ifield
         ifield = 1
         n = nx1*ny1*nz1*nelv
         call vcospf(tt,bm1,n)
         call cmult(tt,eps,n)
         call dssum(tt,nx1,ny1,nz1)
         ifield = ifielt
      else
         n = nx2*ny2*nz2*nelv
         call vcospf(tt,bm2,n)
         call cmult(tt,eps,n)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine vcospf(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = cos(1000.*y(i))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine vbyte_swap(x,n)
      character*1 x(0:3,1),tmp0,tmp1
      character*1 in (0:3), out(0:3)
      real*4      in4     , out4
      equivalence (in ,in4 )
      equivalence (out,out4)
c
      do i=1,n
         do j=0,3
            in (j) = x(j,i)
         enddo
         tmp0   = x(0,i)
         tmp1   = x(1,i)
         x(0,i) = x(3,i)
         x(1,i) = x(2,i)
         x(2,i) = tmp1
         x(3,i) = tmp0
         do j=0,3
            out(j) = x(j,i)
         enddo
         write(6,*) 'swap:',i,in4,out4
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)
      include 'SIZE'
c
      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern
c
      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.
c
      test2 = bytetest
      call byte_reverse(test2,1)
      if (nid.eq.0) 
     $   write(6,*) 'Byte swap:',if_byte_swap_test,bytetest,test2
      return
      end
c-----------------------------------------------------------------------
      subroutine dsavg(u)
c
c
      include 'SIZE'
      include 'TOTAL'
      real u(lx1,ly1,lz1,lelt)
c
c
c     Take direct stiffness avg of u
c
c
      ifieldo = ifield
      if (ifflow) then
         ifield = 1
         ntot = nx1*ny1*nz1*nelv
         call dssum(u,nx1,ny1,nz1)
         call col2 (u,vmult,ntot)
      else
         ifield = 2
         ntot = nx1*ny1*nz1*nelt
         call dssum(u,nx1,ny1,nz1)
         call col2 (u,tmult,ntot)
      endif
      ifield = ifieldo
c
      return
      end
c-----------------------------------------------------------------------

c     Parallel file i/o reader  (param(66) = 6)
c
c     6/14/2006  pff
c
c     File format:

c     132 byte header  (containing nel_block, etc.)
c     4-byte swap test
c     nel_block 4-byte integers, indicating elements contained in file
c     nxyz . nel_block vec data
c     nxyz . nel_block vec data
c     nxyz . nel_block scalar data


c-----------------------------------------------------------------------
      subroutine mfi_gets(u,wk,lwk)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,lelv)
      real*4 wk(lwk)
      integer e,eg,msg_id(lelt)

      common /scruz/ w2(4*lx1*ly1*lz1*lelv)

      call gsync() ! clear outstanding message queues.

      num_recv  = nxr*nyr*nzr*nelv
      num_avail = lwk*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_gets a')

      num_read  = nxr*nyr*nzr
      num_avail = 4*lx1*ly1*lz1*lelv*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_gets b')

      len   = nxr*nyr*nzr*wdsizr
      nxyzr = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzr = 2*nxyzr
      l = 1
      do e=1,nelv
         eg = lglel(e,node)
         msg_id(e) = irecv(eg,wk(l),len)
         l = l+nxyzr
      enddo

      if (nid.eq.pid0) then ! open file fid0
c        call byte_open_indx(mfi_fname,indx)
         do e=1,nelr
            jnid = gllnid(er(e))
            call byte_read(w2,nxyzr)
            call csend(er(e),w2,len,jnid,0)  ! blocking send
         enddo
c        call byte_close_indx(indx)
      endif

c     We don't currently support byte_swap for 64-bit reals
      if (if_byte_sw.and.wdsizr.eq.8) write(6,*) 'ERROR mfi_gets a'
      if (if_byte_sw.and.wdsizr.eq.8) call exitt

      nxyzr = nxr*nyr*nzr
      nxyzw = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzw = 2*nxyzw

      l = 1
      do e=1,nelv
         call msgwait(msg_id(e))
         if (if_byte_sw) call byte_reverse(wk(l),nxyzr)
         if (nxr.eq.nx1.and.nyr.eq.ny1.and.nzr.eq.nz1) then
            if (wdsizr.eq.4) then         ! COPY
               call copy4r(u(1,e),wk(l),nxyzr)
            else
               call copy  (u(1,e),wk(l),nxyzr)
            endif
         else                             ! INTERPOLATE
            if (wdsizr.eq.4) then
               call mapab4r(u(1,e),wk(l),nxr,1)
            else
               call mapab  (u(1,e),wk(l),nxr,1)
            endif
         endif
         l = l+nxyzr
      enddo

      call gsync() ! clear outstanding message queues.

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_getv(u,v,w,wk,lwk)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,lelv),v(lx1*ly1*lz1,lelv),w(lx1*ly1*lz1,lelv)
      real*4 wk(lwk)
      integer e,eg,msg_id(lelt)

      common /scruz/ w2(4*lx1*ly1*lz1*lelv)

      call gsync() ! clear outstanding message queues.

      num_recv  = ndim*nxr*nyr*nzr*nelv
      num_avail = lwk*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv a')

      num_read  = ndim*nxr*nyr*nzr
      num_avail = 4*lx1*ly1*lz1*lelv*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv b')

      len   = ndim*nxr*nyr*nzr*wdsizr
      nxyzr = ndim*nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzr = 2*nxyzr

      l = 1
      do e=1,nelv
         eg = lglel(e,node)
         msg_id(e) = irecv(eg,wk(l),len)
         l = l+nxyzr
      enddo

      if (nid.eq.pid0) then ! open file fid0
c        call byte_open_indx(mfi_fname,indx)
         do e=1,nelr
            jnid = gllnid(er(e))
            call byte_read(w2,nxyzr)
            call csend(er(e),w2,len,jnid,0)  ! blocking send
         enddo
c        call byte_close_indx(indx)
      endif

      if (if_byte_sw.and.wdsizr.eq.8) write(6,*) 'ERROR mfi_getv a'
      if (if_byte_sw.and.wdsizr.eq.8) call exitt

      nxyzr = nxr*nyr*nzr
      nxyzv = ndim*nxr*nyr*nzr
      nxyzw = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzw = 2*nxyzw

      l = 1
      do e=1,nelv

         call msgwait(msg_id(e))

         if (if_byte_sw) call byte_reverse(wk(l),nxyzv)

         if (nxr.eq.nx1.and.nyr.eq.ny1.and.nzr.eq.nz1) then
            if (wdsizr.eq.4) then         ! COPY
               call copy4r(u(1,e),wk(l        ),nxyzr)
               call copy4r(v(1,e),wk(l+  nxyzw),nxyzr)
               if (if3d) 
     $         call copy4r(w(1,e),wk(l+2*nxyzw),nxyzr)
            else
               call copy  (u(1,e),wk(l        ),nxyzr)
               call copy  (v(1,e),wk(l+  nxyzw),nxyzr)
               if (if3d) 
     $         call copy  (w(1,e),wk(l+2*nxyzw),nxyzr)
            endif
         else                             ! INTERPOLATE
            if (wdsizr.eq.4) then
               call mapab4r(u(1,e),wk(l        ),nxr,1)
               call mapab4r(v(1,e),wk(l+  nxyzw),nxr,1)
               if (if3d) 
     $         call mapab4r(w(1,e),wk(l+2*nxyzw),nxr,1)
            else
               call mapab  (u(1,e),wk(l        ),nxr,1)
               call mapab  (v(1,e),wk(l+  nxyzw),nxr,1)
               if (if3d) 
     $         call mapab  (w(1,e),wk(l+2*nxyzw),nxr,1)
            endif
         endif
         l = l+ndim*nxyzw
      enddo

      call gsync() ! clear outstanding message queues.

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_hdr(hdr)
      include 'SIZE'

      character*132 hdr

      if (indx1(hdr,'#std',4).eq.1) then
          call parse_std_hdr(hdr)
      else
         if (nid.eq.0) write(6,80) hdr
         if (nid.eq.0) write(6,80) 'NONSTD HDR, parse_hdr, abort.'
   80    format(a80)
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr(hdr)
      include 'SIZE'
      include 'INPUT'
      include 'RESTART'

      character*132 hdr

c                4  7  10  13   23    33    53    62     68     74
      read(hdr,1) wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         , ifiler,nfiler
     $         , (rdcode1(k),k=1,20)                   ! 74+20=94
    1 format(4x,i2,3i3,2i10,e20.13,i9,2i6,20a1)

c     Assign read conditions, according to rdcode
c          NOTE:  This will be extended to general case in future.
c                 For now, what you see in file is what you get.

      ifgetx = .false.
      ifgetu = .false.
      ifgetp = .false.
      ifgett = .false.
      do k=1,npscal
         ifgtps(k) = .false.
      enddo

      ifgtim = .true.  ! this is the default

      if (rdcode1(1).eq.'X') ifgetx = .true.
      if (rdcode1(2).eq.'U') ifgetu = .true.
      if (rdcode1(3).eq.'P') ifgetp = .true.
      if (rdcode1(4).eq.'T') ifgett = .true.
      do k=1,npscal
         if (rdcode1(4+k).ne.' ') ifgtps(k) = .true.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_get_hdr0(hdr,fname)
c
      include 'SIZE'
      include 'PARALLEL'
      character*132 hdr
      character*80  fname

      if (nid.eq.0) then
         call mbyte_open(fname,0) ! open  blah000.fldnn
         call blank     (hdr,132)
         call byte_read (hdr, 33) ! 4 x 33 = 132
         call byte_close()
      endif
      call bcast(hdr,132)

      return
      end
c-----------------------------------------------------------------------
      subroutine mlti_file_input(ifile)
c
c     (1) Open restart file(s)
c     (2) Check previous spatial discretization 
c     (3) Map (K1,N1) => (K2,N2) if necessary
c
c     nfiles > 1 has several implications:
c
c     i.   For std. run, data is taken from last file in list, unless
c          explicitly specified in argument list of filename
c
c     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
c          subsequent files are for B-field or perturbation fields
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      character*132 hdr

      parameter (lwk = 7*lx1*ly1*lz1*lelv)
      common /scrns/ wk(lwk)

      if (nid.eq.0) write(6,*) 'IFILE:',ifile
      if (nid.eq.0) write(6,*) initc(ifile)
      call mfi_get_hdr0(hdr,initc(ifile))
      call parse_hdr   (hdr)

      call set_pid(initc(ifile))! determine reader nodes; open files

      if (ifgetx) call mfi_getv(xm1,ym1,zm1,wk,lwk)
      if (ifgetu) then
         call mfi_getv(vx,vy,vz,wk,lwk)
      endif
      if (ifgett) call mfi_gets(t,wk,lwk)
      do k=1,npscal
         if (ifgtps(k)) call mfi_gets(t(1,1,1,1,k+1),wk,lwk)
      enddo

      if (ifgtim) time = timer

      if (nid.eq.pid0) call byte_close()

      return
      end
c-----------------------------------------------------------------------
      subroutine mbyte_open(hname,fid) ! open  blah000.fldnn
      include 'SIZE'
      include 'TSTEP'
c
      integer fid
      character*80 hname

      character*6  six,fmt,s6
      save         six
      data         six / '!!!!!!' /
      !data         six / '??????' /

      character*80 fname
      character*1  fname1(80)
      equivalence (fname1,fname)

      integer      iname(20)
      equivalence (iname,fname)

      call  izero (iname,20)
      len = ltrunc(hname,80)
      call chcopy (fname,hname,len)

      do ipass=1,2      ! 2nd pass, in case 1 file/directory
         do k=6,1,-1
            i1 = indx1(fname,six,k)
            if (i1.ne.0) then
               write(fmt,1) k,k
    1          format('(i',i1,'.',i1,')')
               write(s6,fmt) fid
               call chcopy(fname1(i1),s6,k)
               goto 10
            endif
         enddo
   10    continue
      enddo
      
      write(6,6) nid,istep,len,(fname1(k),k=1,len)
    6 format(2i8,i3,' OPEN: ',80a1)

      call byte_open(fname)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_pid(hname)  ! determine which nodes are readers
      character*80 hname

      include 'SIZE'
      include 'PARALLEL'
      include 'RESTART'

      integer stride
      character*132 hdr
      logical if_byte_swap_test
      real*4 bytetest

      stride = np / nfiler
c     write(6,*) nid,' STRIDE:',stride,np,nfiler

      if (stride.lt.1) then
         write(6,*) nfiler,np,'  TOO MANY FILES, set_pid abort'
         call exitt
      endif

      if (mod(nid,stride).eq.0) then
         pid0 = nid
         pid1 = nid + stride
         fid0 = nid / stride

         call mbyte_open(hname,fid0) ! open  blah000.fldnn

         call blank     (hdr,132)
         call byte_read (hdr, 33) ! 4 x 33 = 132
         call parse_hdr (hdr)     ! Re-parse hdr, for nelr info

         call byte_read (bytetest,1)
         if_byte_sw = if_byte_swap_test(bytetest)

         call byte_read (er,nelr) ! Get list of elements to be read

      endif

      call lbcast(if_byte_sw)   ! broadcast byte swap from node 0

      return
      end
c-----------------------------------------------------------------------
