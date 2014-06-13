      subroutine prepost(ifdoin,prefix)
C-------------------------------------------------------------------
C
C     Store results for later postprocessing
C
C------------------------------------------------------------------
C
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'GEOM'
      INCLUDE 'CTIMER'
      INCLUDE 'IXYZ'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'PARALLEL'
      INCLUDE 'SOLN'
C
C     Work arrays and temporary arrays
C
      COMMON /SCRNS/ VXAX   (LX1,LY1,LELV)
     $             , VYAX   (LX1,LY1,LELV)
     $             , PRAX   (LX2,LY2,LELV)
     $             , PM1    (LX1,LY1,LZ1,LELV)
     $             , YAX    (LX1,LY1,LELT)
     $             , TAX    (LX1,LY1,LELT,LDIMT)
c
C     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 3.
      PARAMETER (LXYZ=LX1*LY1*LZ1) 
      PARAMETER (LPSC9=LDIMT+9)
      COMMON /CTMP1/ TDUMP(LXYZ,LPSC9)
      real*4         tdump
      CHARACTER*4    ADUMP(LXYZ,LPSC9)
      EQUIVALENCE   (ADUMP,TDUMP)
      real*4         test_pattern
c
      character*3    prefix
c
      character*1 fhdfle1(80)
      character*80 fhdfle
      equivalence (fhdfle,fhdfle1)
      character*1  fldfile2(120)
      integer      fldfilei( 60)
      equivalence (fldfilei,fldfile2)
c
      real           tdmp(4)
      EQUIVALENCE   (TDUMP,tdmp)
c
      COMMON /DOIT/ IFDOIT
      LOGICAL       IFDOIT
      LOGICAL       IFDOIN
C     
      REAL AXISM1 (LX1,LY1)
      REAL AXISM2 (LX2,LY2)
      REAL PA     (LX1,LY2,LZ2)
      REAL PB     (LX1,LY1,LZ2)
      real*8  AXM18,AXM28,PA8,PB8
      EQUIVALENCE (AXM18,AXISM1)
      EQUIVALENCE (AXM28,AXISM2)
      EQUIVALENCE (PA8,PA)
      EQUIVALENCE (PB8,PB)
C
      REAL DMPTIM(1000),HDUMP(25)
      REAL XPART(10),YPART(10),ZPART(10)
      CHARACTER*2 EXCODE(15)
      CHARACTER*10 FRMAT
      INTEGER NOPEN(99)
      SAVE    NOPEN
      DATA    NOPEN /99*0/
c     INTEGER ICALLD
c     SAVE    ICALLD
c     DATA    ICALLD /0/
      COMMON /RDUMP/ NTDUMP
      DATA NDUMPS / 0 /
c
      integer maxstep
      save    maxstep
      data    maxstep /99999999/
C
C     On first call, open field file
      etime1=dclock()
c
c     Hardcode param(65)=1 so that we *always* get .fld01, .fld02, etc.
c     pff 7/21/98
c
      param(65)=1
c
      IF (NID.EQ.0.AND.ICALLD.EQ.0) THEN
         IF (PARAM(65).EQ.0.0) THEN
           IF (PARAM(66).EQ.1.0) THEN
             OPEN(UNIT=24,FILE=FLDFLE,FORM='UNFORMATTED',STATUS='NEW')
           ELSEIF (PARAM(66).ge.3.0) THEN
              call  izero    (fldfilei,33)
              len = ltrunc   (fldfle,131)
              call chcopy    (fldfile2,fldfle,len)
              call byte_open (fldfile2)
              if (param(66).eq.3.0) then
c                open header file
                 ihd = indx1(fldfle,'.fld',4)
                 call chcopy(fhdfle1,fldfle,80)
                 call chcopy(fhdfle1(ihd),'.fhd',4)
                 OPEN(UNIT=27,FILE=FHDFLE,FORM='FORMATTED'
     $               ,STATUS='UNKNOWN')
              else
c                write header as character string
                 call blank(fhdfle,80)
              endif
c
              if (param(66).eq.5.0) then
                 write(6,*) 'param(66)=5 not supported, reset to 4'
                 param(66)=4.
              endif
c
           ELSE
             OPEN(UNIT=24,FILE=FLDFLE,FORM='FORMATTED',STATUS='UNKNOWN')
           ENDIF
         ENDIF
c        write(6,*) 'hisfile:',hisfle
         OPEN(UNIT=25,FILE=HISFLE,FORM='FORMATTED',STATUS='NEW')
         OPEN(UNIT=26,FILE=SCHFLE,FORM='FORMATTED',STATUS='NEW')
      ENDIF
C
      DENSITY = AVTRAN(1)
C
      IF (IFAXIS) THEN
         NTOTM1 = NX1*NY1*NELT
         CALL COPY (YAX,YM1,NTOTM1)
         DO 5 IEL=1,NELT
            IF (IFRZER(IEL)) THEN
               CALL MXM  (YM1(1,1,1,IEL),NX1,IATJL1,NY1,AXISM1,NY1)
               CALL COPY (YM1(1,1,1,IEL),AXISM1,NX1*NY1)
            ENDIF
    5    CONTINUE
         IF (IFFLOW) THEN
            NTOTM1 = NX1*NY1*NELV
            NTOTM2 = NX2*NY2*NELV
            CALL COPY (VXAX,VX,NTOTM1)
            CALL COPY (VYAX,VY,NTOTM1)
            CALL COPY (PRAX,PR,NTOTM2)
            DO 10 IEL=1,NELV
               IF (IFRZER(IEL)) THEN
                  CALL MXM  (VX(1,1,1,IEL),NX1,IATJL1,NY1,AXISM1,NY1)
                  CALL COPY (VX(1,1,1,IEL),AXISM1,NX1*NY1)
                  CALL MXM  (VY(1,1,1,IEL),NX1,IATJL1,NY1,AXISM1,NY1)
                  CALL COPY (VY(1,1,1,IEL),AXISM1,NX1*NY1)
                  CALL MXM  (PR(1,1,1,IEL),NX2,IATJL2,NY2,AXISM2,NY2)
                  CALL COPY (PR(1,1,1,IEL),AXISM2,NX2*NY2)
               ENDIF
 10         CONTINUE
         ENDIF
         IF (IFHEAT) THEN
            NTOTM1 = NX1*NY1*NELT
            DO 15 IFLDT=1,NPSCAL+1
               CALL COPY (TAX(1,1,1,IFLDT),T(1,1,1,1,IFLDT),NTOTM1)
 15         CONTINUE
            DO 30 IEL=1,NELT
               IF (IFRZER(IEL)) THEN
                 DO 25 IFLDT=1,NPSCAL+1
                   CALL MXM  (T(1,1,1,IEL,IFLDT),NX1,IATJL1,NY1,
     $                                               AXISM1,NY1)
                   CALL COPY (T(1,1,1,IEL,IFLDT),AXISM1,NX1*NY1)
 25              CONTINUE
               ENDIF
 30         CONTINUE
         ENDIF
      ENDIF
C
C     Map the pressure onto the velocity mesh
C
      NGLOB1 = NX1*NY1*NZ1*NELV
      NYZ2   = NY2*NZ2
      NXY1   = NX1*NY1
      NXYZ   = NX1*NY1*NZ1
C
      IF (IFSPLIT) THEN
         CALL COPY(PM1,PR,NGLOB1)
      ELSE
         DO 1000 IEL=1,NELV
            CALL MXM (IXM21,NX1,PR(1,1,1,IEL),NX2,PA(1,1,1),NYZ2)        
            DO 100 IZ=1,NZ2
               CALL MXM (PA(1,1,IZ),NX1,IYTM21,NY2,PB(1,1,IZ),NY1)
 100        CONTINUE
            CALL MXM (PB(1,1,1),NXY1,IZTM21,NZ2,PM1(1,1,1,IEL),NZ1)
 1000    CONTINUE
      ENDIF
C
C      IF(IFFLOW)THEN
C        Average the pressure on elemental boundaries
C        NOT IN THE NEW VERSION !!! (ER)
C
C         IFIELD=1
C         CALL DSSUM (PM1,NX1,NY1,NZ1)
C         CALL COL2 (PM1(1,1,1,1),VMULT(1,1,1,1),NGLOB1)
C      ENDIF
C
C     Store the Gauss-Legendre Lobatto collocation points, 
C     and the corresponding velocities and pressure.
C   
c
c     Changed, 8/5/99   pff
c
c     IF(NSTEPS.EQ.0) NSTEPS=1
c     IF(FINTIM .EQ. 0.0 .AND. ISTEP .EQ. NSTEPS) LASTEP=1
      IF(ISTEP .GE. NSTEPS) LASTEP=1
c
c
      TIMDUMP=0
      IF(TIMEIO.NE.0.0)THEN
         IF(TIME .GE. (NTDUMP + 1) * TIMEIO) THEN
            TIMDUMP=1.
            NTDUMP=NTDUMP+1
         ENDIF
      ENDIF
      IIIDMP=0
      IF(ISTEP.GT.0 .AND. IOSTEP.GT.0)THEN
         IF(MOD(ISTEP,IOSTEP) .EQ. 0) ifdoit=.true.
      ENDIF
C
C--------------------------------------------------------
C     Check for io request in file 'ioinfo'
C
      if (nid.eq.0) then
         open(unit=87,file='ioinfo',status='old',err=88)
         read(87,*,end=87,err=87) idummy
         if (iiidmp.eq.0) iiidmp=idummy
C
         if (idummy.ne.0) then
C           overwrite last i/o request!
            rewind(87)
            write(87,86)
   86       format(' 0')
         endif
C
   87    continue
         close(unit=87)
   88    continue
	 if (iiidmp.ne.0) write(6,*) 'Output:',iiidmp
      endif
C
      tdmp(1)=iiidmp
      call gop(tdmp,tdmp(3),'+  ',1)
      iiidmp= tdmp(1)
      if (iiidmp.lt.0) maxstep=abs(iiidmp)
      if (istep.ge.maxstep.or.iiidmp.eq.-2) lastep=1
      if (iiidmp.eq.-2) return
      if (iiidmp.lt.0) iiidmp = 0
C--------------------------------------------------------
C
C
C
      IF (IFDOIN) IFDOIT=.TRUE.
      IF (IIIDMP.NE.0.OR.LASTEP.EQ.1.OR.TIMDUMP.EQ.1.) IFDOIT=.TRUE.
      IF (IFDOIT) THEN
C
C        We will output a dump.
C
         iprefix = i_find_prefix(prefix,99)
         IF (PARAM(65).NE.0.0.AND.NID.EQ.0) THEN
C          Open new file for each dump on /cfs
           NOPEN(iprefix)=NOPEN(iprefix)+1
           CALL FILE2(NOPEN(iprefix),prefix)
           IF (PARAM(66).EQ.1.0) THEN
             OPEN(UNIT=24,FILE=FLDFLE,FORM='UNFORMATTED',STATUS='NEW')
           ELSEIF (PARAM(66).ge.3.0) THEN
              call  izero    (fldfilei,33)
              len = ltrunc   (fldfle,131)
              call chcopy    (fldfile2,fldfle,len)
              call byte_open (fldfile2)
c
              if (param(66).eq.3.0) then
c                open header file
                 ihd = indx1(fldfle,'.fld',4)
                 call chcopy(fhdfle1,fldfle,80)
                 call chcopy(fhdfle1(ihd),'.fhd',4)
                 OPEN(UNIT=27,FILE=FHDFLE,FORM='FORMATTED'
     $               ,STATUS='UNKNOWN')
              else
c                write header as character string
                 call blank(fhdfle,80)
              endif
c
              if (param(66).eq.5.0) then
                 write(6,*) 'param(66)=5 not supported, reset to 4'
                 param(66)=4.
              endif
           ELSE
             OPEN(UNIT=24,FILE=FLDFLE,FORM='FORMATTED',STATUS='UNKNOWN')
           ENDIF
         ENDIF
C        Figure out what goes in EXCODE
         NDUMPS=NDUMPS+1
         DMPTIM(NDUMPS)=TIME
         CALL BLANK(EXCODE,30)
         IF(IFXYO)THEN
            EXCODE(1)='X'
            EXCODE(2)='Y'
            IF(IF3D)EXCODE(3)='Z'
         ENDIF
         IF(IFVO)THEN
            EXCODE(4)='U'
         ENDIF
         IF(IFPO)EXCODE(5)='P'
         IF(IFTO)EXCODE(6)='T'
         IF(NPSCAL.GT.0)THEN
            DO 1905 IIP=1,NPSCAL
               IF(IFPSCO(IIP))WRITE(EXCODE(IIP+6),'(I1)')IIP
1905        CONTINUE
         ENDIF
c
C        Dump
c
         IF (NID.EQ.0) call dump_header(excode)
c
         call get_id(id)
c
         do ieg=1,nelgt
c
            jnid = gllnid(ieg)
            ie   = gllel (ieg)
c
            if (nid.eq.0) then
               if (jnid.eq.0) then
                  call fill_tmp(tdump,id,ie)
               else
                  mtype=2000+ieg
                  len=4*id*nxyz
                  dum1=0.
                  call csend(mtype,dum1,wdsize,jnid,nullpid)
                  call crecv(mtype,tdump,len)
               endif
               call out_tmp(id)
            elseif (nid.eq.jnid) then
               call fill_tmp(tdump,id,ie)
               dum1=0.
c
               mtype=2000+ieg
               len=4*id*nxyz
               call crecv(mtype,dum1,wdsize)
               call csend(mtype,tdump,len,node0,nullpid)
            endif
         enddo
c
         if (nid.eq.0) call close_fld
c
      endif
C
C
C
      IF (LASTEP.EQ.1) THEN
         MAXHIS = PARAM(16)
         IF (NHIS.GT.0) CALL SKIPHIS (ISTEP,MAXHIS)
         IF (NID.EQ.0) THEN
            CLOSE(UNIT=25)
            CLOSE(UNIT=26)
c
            IF (param(65).eq.0.0) then
               IF (PARAM(66).ge.3) then
                  call byte_close()
                  IF (PARAM(66).eq.3) CLOSE(UNIT=27)
               else
                  CLOSE(UNIT=24)
               endif
            endif
c
         ENDIF
      ENDIF
C
      ICALLD = ICALLD+1
      nprep=icalld
      tprep=tprep+dclock()-etime1
      IFDOIT=.FALSE.

      return
      END
c=======================================================================
      subroutine skiphis (kstep,maxhis)
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
C
      DIMENSION PARTDMP(3),HISDMP(7),QINTDMP(7)
      CHARACTER CHDMP*10
C
      integer icalld
      save    icalld
      data    icalld /0/
C
      IF (NID.NE.0) return
      IF (ICALLD.NE.0) return
      ICALLD=1
      REWIND 25
      REWIND 26
C
      NSKIP  = INT( (KSTEP - 1)/MAXHIS ) + 1
C
      DO 100 IS=1,KSTEP
      DO 100 IH=1,NHIS
C
         IF (HCODE(10,IH).EQ.'P') THEN
            READ (26,'(4G14.6,A10)',ERR=500,END=500)
     $       TIME,(PARTDMP(I),I=1,3),CHDMP
            GOTO 500
         ENDIF
         IF (HCODE(10,IH).EQ.'H') THEN
            ISKH=0
            DO 200 I=1,6
               IF (HCODE(I,IH).NE.' ') ISKH=ISKH + 1
  200       CONTINUE
            READ (26,'(5G14.6)',ERR=500,END=500)
     $      TIME,(HISDMP(I),I=1,ISKH)
            GOTO 500
         ENDIF
         IF (HCODE(10,IH).EQ.'I') THEN
            ISKQ=0
            DO 300 I=1,3
               IF (HCODE(I,IH).NE.' ') ISKQ=ISKQ + 1
  300       CONTINUE
            DO 310 I=5,6
               IF (HCODE(I,IH).NE.' ') ISKQ=ISKQ + 1
  310       CONTINUE
            READ (26,'(5G14.6)',ERR=500,END=500)
     $      TIME,(QINTDMP(I),I=1,ISKQ)
         ENDIF
C
  500    IF ( MOD(IS,NSKIP).NE.0 ) GOTO 100
C
         IF (HCODE(10,IH).EQ.'P') THEN
            WRITE (25,'(4G14.6,A10)') TIME,(PARTDMP(I),I=1,3),CHDMP
            GOTO 100
         ENDIF
         IF (HCODE(10,IH).EQ.'H') THEN
            WRITE (25,'(1p5e16.8)')TIME,(HISDMP(I),I=1,ISKH)
            GOTO 100
         ENDIF
         IF (HCODE(10,IH).EQ.'I') THEN
            WRITE (25,'(5G14.6)')TIME,(QINTDMP(I),I=1,ISKQ)
         ENDIF
C
  100 CONTINUE
C
      return
      END
c=======================================================================
      subroutine file2(nopen,PREFIX)
C----------------------------------------------------------------------
C
C     Defines machine specific input and output file names.
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'PARALLEL'
C
      CHARACTER*132 NAME
      CHARACTER*1   SESS1(132),PATH1(132),NAM1(132)
      EQUIVALENCE  (SESSION,SESS1)
      EQUIVALENCE  (PATH,PATH1)
      EQUIVALENCE  (NAME,NAM1)
      CHARACTER*1  DMP(4),FLD(4),REA(4),HIS(4),SCH(4) ,ORE(4), NRE(4)
      CHARACTER*4  DMP4  ,FLD4  ,REA4  ,HIS4  ,SCH4   ,ORE4  , NRE4
      EQUIVALENCE (DMP,DMP4), (FLD,FLD4), (REA,REA4), (HIS,HIS4)
     $          , (SCH,SCH4), (ORE,ORE4), (NRE,NRE4)
      CHARACTER*1  NUMRL(0:9)
      DATA DMP4,FLD4,REA4 /'.dmp','.fld','.rea'/
      DATA HIS4,SCH4      /'.his','.sch'/
      DATA ORE4,NRE4      /'.ore','.nre'/
      DATA NUMRL          /'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER*78  STRING
c
      character*1    prefix(3)
C
      call blank(name  ,132)
      call blank(fldfle,132)
C
C     Construct file names containing full path to host:
C
      LS=LTRUNC(SESSION,132)
      LPP=LTRUNC(PATH,132)
      LSP=LS+LPP
      l = 0
      DO 100 I=1,LPP
         l = l+1
         NAM1(l)=PATH1(I)
  100 CONTINUE
c
      if (prefix(1).ne.' '.and.prefix(2).ne.' '.and.
     $     prefix(3).ne.' ') then
         do i=1,3
            l = l+1
            NAM1(l)=prefix(i)
         enddo
      endif
C
      DO 200 I=1,LS
         l = l+1
         NAM1(l)=SESS1(I)
  200 CONTINUE
C
C .fld file
      DO 300 I=1,4
         l = l+1
         NAM1(l)=FLD(I)
  300 CONTINUE
      if (nopen.lt.100) then
C        less than 100 dumps....
         ITEN=NOPEN/10
         l = l+1
         NAM1(l)=NUMRL(ITEN)
         IONE=MOD(NOPEN,10)
         l = l+1
         NAM1(l)=NUMRL(IONE)
      elseif (nopen.lt.1000) then
C        less than 1000 dumps....
         IHUN=NOPEN/100
         l = l+1
         NAM1(l)=NUMRL(IHUN)
         ITEN=MOD(NOPEN,100)/10
         l = l+1
         NAM1(l)=NUMRL(ITEN)
         IONE=MOD(NOPEN,10)
         l = l+1
         NAM1(l)=NUMRL(IONE)
      elseif (nopen.lt.10000) then
C        less than 10000 dumps....
         ITHO=NOPEN/1000
         l = l+1
         NAM1(l)=NUMRL(ITHO)
         IHUN=MOD(NOPEN,1000)/100
         l = l+1
         NAM1(l)=NUMRL(IHUN)
         ITEN=MOD(NOPEN,100)/10
         l = l+1
         NAM1(l)=NUMRL(ITEN)
         IONE=MOD(NOPEN,10)
         l = l+1
         NAM1(l)=NUMRL(IONE)
      endif
      FLDFLE=NAME
C
C     Write the name of the .fld file to the logfile.
C
      IF (NID.EQ.0) THEN
         CALL CHCOPY(STRING,FLDFLE,78)
         WRITE(6,1000) istep,time,STRING
 1000    FORMAT('step ',i0,', t=',1pe12.4,', file: ',A)
      ENDIF
C
      return
      END
c=======================================================================
      subroutine vcnvert44(a,n)
c     fixed for real*8 input  .... 5/21/96
      REAL*4 A(1)
      REAL*4 B
      CHARACTER*4 CB
      EQUIVALENCE (CB,B)
C
      DO 100 I=1,N
         CALL CONVERT(CB,A(i))
         A(I)=B
  100 CONTINUE
      return
      END
c=======================================================================
      subroutine vcnvert(a,n)
c     fixed for real*8 input  .... 5/21/96
      REAL   A(1)
      REAL*4 B,a4
      CHARACTER*4 CB
      EQUIVALENCE (CB,B)
C
      DO 100 I=1,N
         A4  = A(i)
         CALL CONVERT(CB,A4)
         A(I)=B
  100 CONTINUE
      return
      END
c=======================================================================
      subroutine rzero4(a,n)
      real*4 A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END
c=======================================================================
      subroutine copyX4(a,b,n)
      REAL*4 A(1)
      REAL   B(1)
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c=======================================================================
      subroutine copy4r(a,b,n)
      real   a(1)
      real*4 b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end
c=======================================================================
      function i_find_prefix(prefix,imax)
c
      character*3 prefix
      character*3 prefixes(99)
      save        prefixes
      data        prefixes /99*'...'/
c
      integer nprefix
      save    nprefix
      data    nprefix /0/
c
c     Scan existing list of prefixes for a match to "prefix"
c
      do i=1,nprefix
         if (prefix.eq.prefixes(i)) then
            i_find_prefix = i
            return
         endif
      enddo
c
c     If we're here, we didn't find a match.. bump list and return
c
      nprefix                = nprefix + 1
      prefixes(nprefix)      = prefix
      i_find_prefix          = nprefix
c
c     Array bounds check on prefix list
c
      if (nprefix.gt.99.or.nprefix.gt.imax) then
         write(6,*) 'Hey! nprefix too big! ABORT in i_find_prefix'
     $      ,nprefix,imax
         call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dump_header(excode)
c
      include 'SIZE'
      include 'TOTAL'
c
      CHARACTER*2 EXCODE(15)
c
      real*4         test_pattern
c
      character*1 fhdfle1(80)
      character*80 fhdfle
      equivalence (fhdfle,fhdfle1)
c
      ikstep = istep
      do ik=1,10
         if (ikstep.gt.99999) ikstep = ikstep/10
      enddo
c
c     WRITE(6,'(4I4,1pe14.7,I5,1X,15A2,1X,A12)')
c    $  NELGT,NX1,NY1,NZ1,TIME,ikstep,(EXCODE(I),I=1,15),
c    $  'NELT,NX,NY,N'
c
      IF (PARAM(66).EQ.1.0) THEN
C       unformatted i/o
        WRITE(24)
     $  NELGT,NX1,NY1,NZ1,TIME,ISTEP,(EXCODE(I),I=1,15)
      ELSEIF (PARAM(66).EQ.3.0) THEN
C       formatted i/o to header file
        WRITE(27,'(4I4,1pe14.7,I5,1X,15A2,1X,A12)')
     $  NELGT,NX1,NY1,NZ1,TIME,ikstep,(EXCODE(I),I=1,15),
     $  'NELT,NX,NY,N'
      ELSEIF (PARAM(66).eq.4.0) THEN
C       formatted i/o to header file
        WRITE(fhdfle,'(4I4,1pe14.7,I5,1X,15A2,1X,A12)')
     $  NELGT,NX1,NY1,NZ1,TIME,ikstep,(EXCODE(I),I=1,15),
     $  ' 4 NELT,NX,NY,N'
        call byte_write(fhdfle,20)
      ELSEIF (PARAM(66).eq.5.0) THEN
C       formatted i/o to header file
        WRITE(fhdfle,'(4I4,1pe14.7,I5,1X,15A2,1X,A12)')
     $  NELGT,NX1,NY1,NZ1,TIME,ikstep,(EXCODE(I),I=1,15),
     $  ' 8 NELT,NX,NY,N'
        call byte_write(fhdfle,20)
      ELSE
C       formatted i/o
        WRITE(24,'(4I4,1pe14.7,I5,1X,15A2,1X,A12)')
     $  NELGT,NX1,NY1,NZ1,TIME,ikstep,(EXCODE(I),I=1,15),
     $  'NELT,NX,NY,N'
      ENDIF
C     cdrror is a dummy cerror value for now.
      CDRROR=0.0
      IF (PARAM(66).EQ.1.0) THEN
C       unformatted i/o
        WRITE(24)(CDRROR,I=1,NELGT)
      ELSEIF (PARAM(66).eq.3. .or. param(66).eq.4.0) then
C       write byte-ordering test pattern to byte file...
        test_pattern = 6.54321
        call byte_write(test_pattern,1)
      ELSEIF (PARAM(66).eq.5.) then
        test_pattern8 = 6.54321
        call byte_write(test_pattern8,2)
      ELSE
C       formatted i/o
        WRITE(24,'(6G11.4)')(CDRROR,I=1,NELGT)
      ENDIF
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fill_tmp(tdump,id,ie)
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
c
      COMMON /SCRNS/ VXAX   (LX1,LY1,LELV)
     $             , VYAX   (LX1,LY1,LELV)
     $             , PRAX   (LX2,LY2,LELV)
     $             , PM1    (LX1,LY1,LZ1,LELV)
     $             , YAX    (LX1,LY1,LELT)
     $             , TAX    (LX1,LY1,LELT,LDIMT)
C
C     Fill work array
C
      PARAMETER (LXYZ=LX1*LY1*LZ1)
      PARAMETER (LPSC9=LDIMT+9)
      real*4 tdump(lxyz,lpsc9)
C
      nxyz = nx1*ny1*nz1
c
      ID=0
      IF(IFXYO)THEN
         ID=ID+1
         CALL COPYx4(TDUMP(1,ID),XM1(1,1,1,IE),NXYZ)
         ID=ID+1
         CALL COPYx4(TDUMP(1,ID),YM1(1,1,1,IE),NXYZ)
         IF(IF3D) THEN
            ID=ID+1
            CALL COPYx4(TDUMP(1,ID),ZM1(1,1,1,IE),NXYZ)
         ENDIF
      ENDIF
c
      IF(IFVO)THEN
         IF (IE.LE.NELV) THEN
            ID=ID+1
            CALL COPYx4(TDUMP(1,ID),VX(1,1,1,IE),NXYZ)
            ID=ID+1
            CALL COPYx4(TDUMP(1,ID),VY(1,1,1,IE),NXYZ)
            IF(IF3D)THEN
               ID=ID+1
               CALL COPYx4(TDUMP(1,ID),VZ(1,1,1,IE),NXYZ)
            ENDIF
         ELSE
            ID=ID+1
            CALL RZERO4(TDUMP(1,ID),NXYZ)
            ID=ID+1
            CALL RZERO4(TDUMP(1,ID),NXYZ)
            IF(IF3D)THEN
               ID=ID+1
               CALL RZERO4(TDUMP(1,ID),NXYZ)
            ENDIF
         ENDIF
      ENDIF
      IF(IFPO)THEN
         IF (IE.LE.NELV) THEN
            ID=ID+1
            CALL COPYx4(TDUMP(1,ID),PM1(1,1,1,IE),NXYZ)
         ELSE
            ID=ID+1
            CALL RZERO4(TDUMP(1,ID),NXYZ)
         ENDIF
      ENDIF
      IF(IFTO)THEN
         ID=ID+1
         CALL COPYx4(TDUMP(1,ID),T(1,1,1,IE,1),NXYZ)
      ENDIF
C     PASSIVE SCALARS
      IF(NPSCAL.GT.0)THEN
         DO IIP=1,NPSCAL
            IF(IFPSCO(IIP))THEN
               ID=ID+1
               CALL COPYx4(TDUMP(1,ID),T(1,1,1,IE,IIP+1),NXYZ)
           ENDIF             
         ENDDO
      ENDIF             
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_id(id)
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
C
C     Count up amount of data to be shipped
C
      ID=0
      IF(IFXYO)THEN
         ID=ID+1
         ID=ID+1
         IF(IF3D) THEN
            ID=ID+1
         ENDIF
      ENDIF
c
      IF(IFVO)THEN
         ID=ID+1
         ID=ID+1
         IF(IF3D)THEN
            ID=ID+1
         ENDIF
      ENDIF
      IF(IFPO) ID=ID+1
      IF(IFTO) ID=ID+1
C     PASSIVE SCALARS
      IF(NPSCAL.GT.0)THEN
         DO IIP=1,NPSCAL
            IF(IFPSCO(IIP))THEN
               ID=ID+1
           ENDIF             
         ENDDO
      ENDIF             
c
      return
      end
c-----------------------------------------------------------------------
      subroutine close_fld
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
c
      IF (PARAM(65).NE.0.0.AND.NID.EQ.0) then
         IF (PARAM(66).eq.3) then
            CLOSE(UNIT=27)
            call byte_close()
         elseif (param(66).eq.4 .or. param(66).eq.5) then
            call byte_close()
         else
            CLOSE(UNIT=24)
         endif
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_tmp(id)
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
c
      PARAMETER (LXYZ=LX1*LY1*LZ1)
      PARAMETER (LPSC9=LDIMT+9)
c
      COMMON /CTMP1/ TDUMP(LXYZ,LPSC9)
      real*4         tdump
      CHARACTER*4    ADUMP(LXYZ,LPSC9)
      EQUIVALENCE   (ADUMP,TDUMP)
c
      CHARACTER*11 FRMAT
c
      nxyz = nx1*ny1*nz1
c
      call blank(frmat,11)
      if (id.le.9) then
         WRITE(FRMAT,1801) ID
 1801    FORMAT('(1p',I1,'e14.6)')
      else
         WRITE(FRMAT,1802) ID
 1802    FORMAT('(1p',I2,'e14.6)')
      endif
c
      IF (PARAM(66).EQ.0.0) THEN
c
C       formatted i/o
        WRITE(24,FRMAT)
     $      ((TDUMP(I,II),II=1,ID),I=1,NXYZ)
c
      ELSEIF (PARAM(66).EQ.1.0) THEN
c
C        unformatted i/o
         WRITE(24)((TDUMP(I,II),II=1,ID),I=1,NXYZ)
c
      ELSEIF (PARAM(66).EQ.3.0.or.param(66).eq.4.) THEN
c
         do ii=1,id
            call byte_write(tdump(1,ii),nxyz)
         enddo
c
      ELSEIF (PARAM(66).EQ.2.0) THEN
c
C        universal format i/o
         NWORDS=ID*NXYZ
         CALL VCNVERT44(TDUMP,NWORDS)
         WRITE(24,'(20A4)')
     $       ((ADUMP(I,II),II=1,ID),I=1,NXYZ)
c
      ENDIF
c
      return
      end
c-----------------------------------------------------------------------
      subroutine prepost_map(isave) ! isave=0-->fwd, isave=1-->bkwd

c     Store results for later postprocessing

      include 'SIZE'
      include 'TOTAL'
C
C     Work arrays and temporary arrays
C
      common /scruz/ vxax   (lx1,ly1,lelv)
     $             , vyax   (lx1,ly1,lelv)
     $             , prax   (lx2,ly2,lelv)
     $             , yax    (lx1,ly1,lelt)
      common /scrmg/ tax    (lx1,ly1,lelt,ldimt)
      common /scrcg/ pm1    (lx1,ly1,lz1,lelv)
C
c
      common /prepst/ pa(lx1,ly2,lz2),pb(lx1,ly1,lz2)
      integer e

      if (isave.eq.0) then ! map to GLL grid

         if (ifaxis) then
            ntotm1 = nx1*ny1*nelt
            call copy (yax,ym1,ntotm1)
            do 5 e=1,nelt
               if (ifrzer(e)) then
                  call mxm  (ym1(1,1,1,e),nx1,iatjl1,ny1,pb,ny1)
                  call copy (ym1(1,1,1,e),pb,nx1*ny1)
               endif
    5       continue
            if (ifflow) then
               ntotm1 = nx1*ny1*nelv
               ntotm2 = nx2*ny2*nelv
               call copy (vxax,vx,ntotm1)
               call copy (vyax,vy,ntotm1)
               call copy (prax,pr,ntotm2)
               do 10 e=1,nelv
                  if (ifrzer(e)) then
                     call mxm  (vx(1,1,1,e),nx1,iatjl1,ny1,pb,ny1)
                     call copy (vx(1,1,1,e),pb,nx1*ny1)
                     call mxm  (vy(1,1,1,e),nx1,iatjl1,ny1,pb,ny1)
                     call copy (vy(1,1,1,e),pb,nx1*ny1)
                     call mxm  (pr(1,1,1,e),nx2,iatjl2,ny2,pb,ny2)
                     call copy (pr(1,1,1,e),pb,nx2*ny2)
                  endif
 10            continue
            endif
            if (ifheat) then
               ntotm1 = nx1*ny1*nelt
               do 15 ifldt=1,npscal+1
                  call copy (tax(1,1,1,ifldt),t(1,1,1,1,ifldt),ntotm1)
 15            continue
               do 30 e=1,nelt
                  if (ifrzer(e)) then
                    do 25 ifldt=1,npscal+1
                      call mxm  (t(1,1,1,e,ifldt),nx1,iatjl1,ny1,
     $                                                  pb,ny1)
                      call copy (t(1,1,1,e,ifldt),pb,nx1*ny1)
 25                 continue
                  endif
 30            continue
            endif
         endif
C        Map the pressure onto the velocity mesh
C
         ntot1 = nx1*ny1*nz1*nelv
         nyz2  = ny2*nz2
         nxy1  = nx1*ny1
         nxyz  = nx1*ny1*nz1
         nxyz2 = nx2*ny2*nz2
C
         if (ifsplit) then
            call copy(pm1,pr,ntot1)
c        elseif (if_full_pres) then
c           call rzero(pm1,ntot1)
c           do e=1,nelv
c              call copy(pm1(1,1,1,e),pr(1,1,1,e),nxyz2)
c           enddo
         else
            do 1000 e=1,nelv
               call mxm (ixm21,nx1,pr(1,1,1,e),nx2,pa(1,1,1),nyz2)        
               do 100 iz=1,nz2
                  call mxm (pa(1,1,iz),nx1,iytm21,ny2,pb(1,1,iz),ny1)
  100          continue
               call mxm (pb(1,1,1),nxy1,iztm21,nz2,pm1(1,1,1,e),nz1)
 1000       continue
         endif

      else       ! map back

         if (ifaxis) then
            ntot1 = nx1*ny1*nelt
            call copy (ym1,yax,ntot1)
            if (ifflow) then
               ntot1 = nx1*ny1*nelv
               ntot2 = nx2*ny2*nelv
               call copy (vx,vxax,ntot1)
               call copy (vy,vyax,ntot1)
               call copy (pr,prax,ntot2)
            endif
            if (ifheat) then
               ntot1 = nx1*ny1*nelt
               do 3000 ifldt=1,npscal+1
                  call copy (t(1,1,1,1,ifldt),tax(1,1,1,ifldt),ntot1)
 3000          continue
            endif
         endif

      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine outpost(v1,v2,v3,vp,vt,name3)

      include 'SIZE'
      include 'INPUT'

      real v1(1),v2(1),v3(1),vp(1),vt(1)
      character*3 name3


      itmp=0
      if (ifto) itmp=1
      call outpost2(v1,v2,v3,vp,vt,itmp,name3)

      return
      end

      subroutine outpost2(v1,v2,v3,vp,vt,nfldt,name3)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'

      parameter(ltot1=lx1*ly1*lz1*lelt)
      parameter(ltot2=lx2*ly2*lz2*lelv)
      common /outtmp/  w1(ltot1),w2(ltot1),w3(ltot1),wp(ltot2)
     &                ,wt(ltot1,ldimt)
c
      real v1(1),v2(1),v3(1),vp(1),vt(ltot1,1)
      character*3 name3
      logical if_save(ldimt)
c
      ntot1  = nx1*ny1*nz1*nelv
      ntot1t = nx1*ny1*nz1*nelt
      ntot2  = nx2*ny2*nz2*nelv

      if(nfldt.gt.ldimt) then
        write(6,*) 'ABORT: outpost data too large (nfldt>ldimt)!'
        call exitt
      endif

c store solution
      call copy(w1,vx,ntot1)
      call copy(w2,vy,ntot1)
      call copy(w3,vz,ntot1)
      call copy(wp,pr,ntot2)
      do i = 1,nfldt
         call copy(wt(1,i),t(1,1,1,1,i),ntot1t)
      enddo

c swap with data to dump
      call copy(vx,v1,ntot1)
      call copy(vy,v2,ntot1)
      call copy(vz,v3,ntot1)
      call copy(pr,vp,ntot2)
      do i = 1,nfldt
         call copy(t(1,1,1,1,i),vt(1,i),ntot1t)
      enddo

c dump data
      if_save(1) = ifto
      ifto = .false.
      if(nfldt.gt.0) ifto = .true.
      do i = 1,ldimt-1
         if_save(i+1) = ifpsco(i)
         ifpsco(i) = .false.
         if(i+1.le.nfldt) ifpsco(i) = .true.
      enddo

      call prepost(.true.,name3)

      ifto = if_save(1)
      do i = 1,ldimt-1
         ifpsco(i) = if_save(i+1)
      enddo

c restore solution data
      call copy(vx,w1,ntot1)
      call copy(vy,w2,ntot1)
      call copy(vz,w3,ntot1)
      call copy(pr,wp,ntot2)
      do i = 1,nfldt
         call copy(t(1,1,1,1,i),wt(1,i),ntot1t)
      enddo

      return
      end

