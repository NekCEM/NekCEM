      SUBROUTINE READAT
C
C     Read in data supplied by preprocessor and
C     (eventually) echo check.
C
C     NOTE:  This routine has been broken up into several submodules
C            for ease of comprehension and in anticipation of future
C            restructuring of the NEKTON input.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'ZPER'
C
      OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')
C
C     Read Parameters
C
      CALL RDPARAM
C
C     Read Mesh Data and Group ID
C
      CALL RDMESH(ifgtp)

      if (ifgtp) then
c
         call genbox
c
      else
c
C        Read curved side data
C
         CALL RDCURVE
C
C        Read Boundary Conditions (and connectivity data)
C
         CALL RDBDRY
C
      endif
C
C     Read Initial Conditions / Drive Force
C
      CALL RDICDF
C
C     Read materials property data
C
      CALL RDMATP
C
C     Read history data
C
      CALL RDHIST
C
C     Read output specs
C
      CALL RDOUT
C
C     Read objects
C
      CALL RDOBJ
C
C     End of input data, close read file.
C
      CLOSE(UNIT=9)
      if (.not.IFFMTIN) CLOSE(UNIT=8)
C
      RETURN
      END

      SUBROUTINE RDPARAM
C
C     .Read in parameters supplied by preprocessor and
C      (eventually) echo check.
C
C     .Broadcast run parameters to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'ZPER'
      INCLUDE 'EMWAVE'
      INCLUDE 'NEKUSE'
C
      READ(9,*,ERR=400)    !
      READ(9,*,ERR=400) VNEKTON
      NKTONV=(VNEKTON+0.5)
      IF (IFCEM) THEN
      VNEKMIN=0.1
      ELSE
      VNEKMIN=2.5
      IF(VNEKTON.LT.VNEKMIN)THEN
         PRINT*,' Error: This NEKTON Solver Requires a .rea file'
         PRINT*,' from prenek version ',VNEKMIN,' or higher'
         PRINT*,' Please run the session through the preprocessor'
         PRINT*,' to bring the .rea file up to date.'
         call exitt
      ENDIF
      ENDIF
      READ(9,*,ERR=400) NDIM
c     error check
      IF(NDIM.NE.LDIM)THEN
         WRITE(6,10) LDIM,NDIM
   10       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $              /,2X,'       for spatial dimension equal to',I2,'.'
     $              /,2X,'       The data file has dimension',I2,'.')
         call exitt
      ENDIF

      IF (NDIM.EQ.3) IF3D=.TRUE.
      IF (NDIM.NE.3) IF3D=.FALSE.

      READ(9,*,ERR=400) NPARAM
      DO 20 I=1,NPARAM
         READ(9,*,ERR=400)PARAM(I)
   20 CONTINUE

      NPSCAL=INT(PARAM(23))
      NPSCL1=NPSCAL+1
      NPSCL2=NPSCAL+2

      IF (NPSCL1.GT.LDIMT) THEN
         WRITE(6,21) LDIMT,NPSCL1
   21    FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $           /,2X,'       for',I3,' passive scalars.  This run'
     $           /,2X,'       requires that LDIMT be set to',I3,'.')
         call exitt
      ENDIF
    
      READ(9,*,ERR=400) NSKIP

      DO 25 I=1,NSKIP
            READ(9,*,ERR=500)
   25 CONTINUE

C     Read logical equation type descriptors....

      READ(9,*,ERR=500) NLOGIC
      READ(9,*,ERR=500) IFFLOW  
      READ(9,*,ERR=500) IFHEAT
      READ(9,*,ERR=500) IFTRAN 

      IFSRC    =.FALSE.  
      IFCENTRAL=.FALSE.  
      IFUPWIND =.FALSE.
      IFTM     =.FALSE.  
      IFTE     =.FALSE.
      IFDEALIAS=.FALSE.
      IFRK4    =.FALSE.  
      IFEXP    =.FALSE.
      IFEIG    =.FALSE.
      IFNM     =.FALSE.

      IF (NLOGIC.GE.4 ) READ(9,*,ERR=500) IFSRC  
      IF (NLOGIC.GE.5 ) READ(9,*,ERR=500) IFCENTRAL
      IF (NLOGIC.GE.6 ) READ(9,*,ERR=500) IFUPWIND 
      IF (NLOGIC.GE.7 ) READ(9,*,ERR=500) IFTM   
      IF (NLOGIC.GE.8 ) READ(9,*,ERR=500) IFTE     
      IF (NLOGIC.GE.9 ) READ(9,*,ERR=500) IFDEALIAS
      IF (NLOGIC.GE.10) READ(9,*,ERR=500) IFRK4 
      IF (NLOGIC.GE.11) READ(9,*,ERR=500) IFEXP  
      IF (NLOGIC.GE.12) READ(9,*,ERR=500) IFEIG 
      IF (NLOGIC.GE.13) READ(9,*,ERR=500) IFNM   

      MAXLOG=15
      IF(NLOGIC.GT.MAXLOG)THEN
         DO 30 IL=1,NLOGIC-MAXLOG
            READ(9,*,ERR=500)
   30    CONTINUE
      ENDIF
 
      nelx   = abs(param(116))   ! check for global tensor-product structure
      nely   = abs(param(117))
      nelz   = abs(param(118))
      n_o    = 0
 
      if (n_o.eq.0) then
         ifzper=.false.
         ifgfdm=.false.
         if (nelz.gt.0) ifzper=.true.
         if (nelx.gt.0) ifgfdm=.true.
         if (nelx.gt.0) ifzper=.false.
      endif
 
C     End of PARAMETER read for all processors.
C
      RETURN
  
C     Error handling:

  400 CONTINUE
      WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDPARAM.')
      call exitt
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDPARAM.')
      call exitt
C
      RETURN
      END


      SUBROUTINE RDMESH(ifgtp)

C     .Read number of elements
C
C     .Construct sequential element-processor partition according
C      to number of elements and processors
C
C     .Selectively read mesh (defined by element vertices, and group numbers)
C      on each processor
 
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      CHARACTER*1 ADUM
C
      REAL    DUM(4)
c
      logical ifgtp  ! if-global-tensor-product
C
      IFFMTIN = .TRUE.
C
C     Read Elemental Mesh data
C
      READ(9,*,ERR=400,END=400) XFAC,YFAC,XZERO,YZERO
      READ(9,*,ERR=500)
      READ(9,*,ERR=500,END=500) NELGT,NDIM,NELGV
C
      ifgtp = .false.
      if (ndim.lt.0) then
          ifgtp = .true.     ! domain is a global tensor product
          return
      endif
C
      IF (NELGT.GT.0) THEN
C
C-----------------------------------------------------------------
C Read Formatted Mesh data
C-----------------------------------------------------------------
C          
      NELTMX=NP*LELT
      NELVMX=NP*LELV

      IF (NELGT.GT.NELTMX.OR.NELGV.GT.NELVMX) THEN

      if (nid.eq.0) write (6,*) 'help:',lp,np,nelvmx,nelgv,neltmx,nelgt
      if (nid.eq.0) write (6,*) 'help:',lelt,lelv
         IF (NP.EQ.1) WRITE(6,11) LELT,NELGV,NELGT
         IF (NP.GT.1) WRITE(6,12) LELT,NELGV,NELGT
   11       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $             ,/,2X,'       for',I4,' elements.'
     $             ,/,2X,'       The data file has dimensions',2I4,'.')
   12       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $             ,/,2X,'       for',I4,' elements per processor.'
     $             ,/,2X,'       The data file has dimensions',2I4,'.'
     $             ,/,6X,'   Rerun with more processors or recompile.')
         call exitt

      ENDIF
 
C     Construct a tempory element to processor mapping so that we
C     can assimilate all the read data here for a parallel implementation.
C     (This part will need to be slightly modified when we include the 
C     variable K restart in order to read and map the restart mesh.)
 
      CALL MAPELPR()          ! parallel distribution
 
      NSIDES=NDIM*2
      DO IEG=1,NELGT          ! global elt num
         IF (GLLNID(IEG).EQ.NID) THEN
            IEL=GLLEL(IEG)    ! get local elt num
            READ(9,30,ERR=500,END=500) IGROUP(IEL)
   30       FORMAT(43X,I5)
C           Read Corner data
            IF     (NDIM.EQ.2)THEN
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=1,4)
                              call rzero (zc(1 ,iel)     ,4)
            ELSE IF(NDIM.EQ.3)THEN
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (ZC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=5,8)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=5,8)
               READ(9,*,ERR=500,END=500) (ZC(IC,IEL),IC=5,8)
            ENDIF
         ELSE
C           Skip over this data for element NOT on this processor
            READ(9,41,ERR=500,END=500) ADUM
C           Read Corner data
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
         ENDIF
      ENDDO   
   41 FORMAT(A1)
C
C     End of mesh read.
C
      RETURN
C
C     Error handling:
C
  400 CONTINUE
      WRITE(6,401) 
  401 FORMAT(2X,'ERROR READING SCALE FACTORS, CHECK READ FILE'
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt
C
  500 CONTINUE
      WRITE(6,501) IEG
  501 FORMAT(2X,'ERROR READING MESH DATA NEAR ELEMENT',I7
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt
C
      RETURN
C
C-----------------------------------------------------------------
C Else ....  Unformatted read
C-----------------------------------------------------------------
C
      ELSE
C
      IFFMTIN = .FALSE.
      IF (nid.eq.0) write(6,*) 'opening ',re2fle
      OPEN(UNIT=8,FILE=RE2FLE,FORM='UNFORMATTED',STATUS='OLD')
C
      NELGT = -NELGT
      NELTMX=NP*LELT
      NELVMX=NP*LELV
      IF (NELGT.GT.NELTMX.OR.NELGV.GT.NELVMX) THEN
      if (nid.eq.0) write (6,*) 'help:',lp,np,nelvmx,nelgv,neltmx,nelgt
      if (nid.eq.0) write (6,*) 'help:',lelt,lelv
         IF (NP.EQ.1) WRITE(6,11) LELT,NELGV,NELGT
         IF (NP.GT.1) WRITE(6,12) LELT,NELGV,NELGT
         call exitt
      ENDIF
C
C     Construct a tempory element to processor mapping so that we
C     can assimilate all the read data here for a parallel implementation.
C     (This part will need to be slightly modified when we include the 
C      variable K restart in order to read and map the restart mesh.)
C
      CALL MAPELPR()
C
      NSIDES=NDIM*2
      DO 1040 IEG=1,NELGT
         IF (GLLNID(IEG).EQ.NID) THEN
            IEL=GLLEL(IEG)
            READ(8,ERR=1500,END=1500) IGROUP(IEL)
C           Read Corner data
            IF(NDIM.EQ.2)THEN
               READ(8,ERR=1500,END=1500) (XC(IC,IEL),IC=1,4)
               READ(8,ERR=1500,END=1500) (YC(IC,IEL),IC=1,4)
            ELSE IF(NDIM.EQ.3)THEN
               READ(8,ERR=1500,END=1500) (XC(IC,IEL),IC=1,4)
               READ(8,ERR=1500,END=1500) (YC(IC,IEL),IC=1,4)
               READ(8,ERR=1500,END=1500) (ZC(IC,IEL),IC=1,4)
               READ(8,ERR=1500,END=1500) (XC(IC,IEL),IC=5,8)
               READ(8,ERR=1500,END=1500) (YC(IC,IEL),IC=5,8)
               READ(8,ERR=1500,END=1500) (ZC(IC,IEL),IC=5,8)
            ENDIF
         ELSE
C           Skip over this data for element NOT on this processor
            READ(8,ERR=1500,END=1500) IDUM
C           Read Corner data
            IF(NDIM.EQ.2)THEN
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
            ELSE IF(NDIM.EQ.3)THEN
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
               READ(8,ERR=1500,END=1500) (DUM(j),j=1,4)
            ENDIF
         ENDIF
 1040 CONTINUE
C
C     End of unformatted mesh read.
      ENDIF
C
      RETURN
C
 1500 CONTINUE
      WRITE(6,1501) IEG
 1501 FORMAT(2X,'ERROR READING unformatted MESH DATA NEAR ELEMENT',I7
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt
C
      RETURN
      END

      SUBROUTINE RDCURVE
C
C     .Read curve side data
C
C     .Disperse curve side data to all processors according 
C      to sequential partition scheme
C
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      CHARACTER*1 ANS
C
c     write(6,*) IFFMTIN=.true. ; current setting
C
      IF (IFFMTIN) THEN   ! FMTIN=ForMaTtedINput
C
C     Read formatted curve side data 
C
      READ(9,*)
      READ(9,*)NCURVE
      CALL RZERO(CURVE ,48*LELT)
      CALL BLANK(CCURVE, 8*LELT)
      IF (NCURVE.GT.0) THEN
         DO 50 ICURVE=1,NCURVE
            IF     (NELGT.LT.1000) THEN
               READ(9,60,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSEIF (NELGT.LT.1000000) THEN
               READ(9,61,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSEIF (NELGT.LT.10000000) THEN
               READ(9,62,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ENDIF
   60       FORMAT(I3,I3,5G14.6,1X,A1)
   61       FORMAT(I2,I6,5G14.6,1X,A1)
   62       FORMAT(I1,I7,5G14.6,1X,A1)
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               CURVE (1,IEDG,IEL)=R1
               CURVE (2,IEDG,IEL)=R2
               CURVE (3,IEDG,IEL)=R3
               CURVE (4,IEDG,IEL)=R4
               CURVE (5,IEDG,IEL)=R5
               CCURVE(  IEDG,IEL)=ANS
            ENDIF
   50    CONTINUE
      ENDIF
      RETURN
C
C     Error handling:
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
      RETURN
C
      ELSE
C
C     Read unformatted curve side data 
C
      READ(8) NCURVE
      CALL RZERO(CURVE ,48*LELT)
      CALL BLANK(CCURVE, 8*LELT)
      IF (NCURVE.GT.0) THEN
         DO 1050 ICURVE=1,NCURVE
            READ(8,ERR=1500,END=1500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               CURVE (1,IEDG,IEL)=R1
               CURVE (2,IEDG,IEL)=R2
               CURVE (3,IEDG,IEL)=R3
               CURVE (4,IEDG,IEL)=R4
               CURVE (5,IEDG,IEL)=R5
               CCURVE(  IEDG,IEL)=ANS
            ENDIF
 1050    CONTINUE
      ENDIF
      RETURN
C
C     Error handling:
C
 1500 CONTINUE
      WRITE(6,1501)
 1501 FORMAT(2X,'ERROR READING unformatted CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
C
      RETURN
      ENDIF
      END

      SUBROUTINE RDBDRY
C
C     .Read Boundary Conditions (and connectivity data)
C
C     .Disperse boundary condition data to all processors 
C      according to sequential partition scheme
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      CHARACTER    CBC1*1,CBC3*3,CHTEMP*1,CHTMP3*3
      EQUIVALENCE (CHTEMP,CHTMP3)
      character*80 string
C
C     Set up TEMPORARY value for NFIELD - NFLDT
C
      NFLDT  = 1                
      IF (IFHEAT) NFLDT=2+NPSCAL
      NBCS   = NFLDT
      IBCS   = 2
      IF (IFFLOW) IBCS = 1
      NSIDES = 2*NDIM
C
C     Read boundary conditions for all fields
C
      LCBC=18*LELT*(LDIMT1 + 1)
      LRBC=30*LELT*(LDIMT1 + 1)

      CALL RZERO(BC ,LRBC)
      CALL BLANK(CBC,LCBC)
C
C-----------------------------------------------------------------
C  Formatted Reads
C-----------------------------------------------------------------
C
      IF (IFFMTIN) THEN
C
      READ(9,*,ERR=500,END=500)  !   ***** BOUNDARY CONDITIONS *****
      ibcnew = 1

      DO 100 IFIELD=ibcnew,NBCS  !     DO 100 IFIELD=IBCS,NBCS

        NEL=NELGT
C       Fluid and/or thermal
        read(9,81) string        !  ***** FLUID   BOUNDARY CONDITIONS *****
        call capit(string,80)

        if (indx1(string,'NO ',3).eq.0) then ! we have acitve bc info
 
         IF(VNEKTON .LE. 2.52) NBCREA = 3
         IF(VNEKTON .GE. 2.55) NBCREA = 5
 
         DO 80 IEG=1,NEL
         DO 80 ISIDE=1,NSIDES
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               IF (NELGT.LT.1000) THEN
                  READ(9,50,ERR=500,END=500)    
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
c  write(6,*) 'cbc',CBC(ISIDE,IEL,IFIELD),ID1,ID2,ISIDE,IEL,IFIELD !<--read 
c                 write(6,50)
c    $            CHTEMP,
c    $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
c    $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
   50             FORMAT(A1,A3,2I3,5G14.7)
               ELSEIF (NELGT.LT.100000) THEN
                  READ(9,51,ERR=500,END=500)    
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
   51             FORMAT(A1,A3,I5,I1,5G14.7)
               ELSEIF (NELGT.LT.10000000) THEN
                  READ(9,52,ERR=500,END=500)    
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
   52             FORMAT(A1,A3,I6,5G14.7)
               ENDIF
C              Mesh B.C.'s in 1st column of 1st field
               IF (CHTEMP.NE.' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP
C              check for fortran function as denoted by lower case bc's:
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
               READ(9,*,ERR=500,END=500)   cbc1  ! dummy read, pff 4/28/05
c              IEL=1
c              READ(9,65,ERR=500,END=500)    
c    $         CBCS(ISIDE,IEL),ID1,ID2,(BCS(II,ISIDE,IEL),II=1,NBCREA)
c  65          FORMAT(1X,A3,2I3,5G14.7)
c              check for fortran function as denoted by lower case bc's:
c              CBC1=CBCS(ISIDE,IEL)
c              CBC3=CBCS(ISIDE,IEL)
c              ICBC1=ICHAR(CBC1)
c              IF (ICBC1.GE.97.AND.ICBC1.LE.122) THEN
c                 IF(CBC3(3:3).NE.'i')NLINES=BCS(1,ISIDE,IEL)
c                 IF(CBC3(3:3).EQ.'i')NLINES=BCS(4,ISIDE,IEL)
c                 DO 70 I=1,NLINES
c  70             READ(9,*,ERR=500,END=500)
c              ENDIF
            ENDIF
   80    CONTINUE
        endif
   81   format(a80)
  100 CONTINUE
C
C     END OF BC READ
C
C     Check for dummy line:  "NO THERMAL B.C.'S"
      IF (NFLDT.EQ.1) READ(9,*,ERR=500,END=500)
C
      RETURN
C
C     Error handling:
C
  500 CONTINUE
      WRITE(6,501) IFIELD,IEG
  501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I6
     $    ,/,2X,'ABORTING IN ROUTINE RDBDRY.')
      call exitt
      RETURN
C
C
      ELSE
C
C-----------------------------------------------------------------
C  UNformatted Reads
C-----------------------------------------------------------------
C
c     READ(8,ERR=500,END=500)

      DO 1100 IFIELD=IBCS,NBCS
         NEL=NELGT
C        Fluid and/or thermal
         NBCREA = 5
C
         DO 1080 IEG=1,NEL
         DO 1080 ISIDE=1,NSIDES
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               READ(8,ERR=1500,END=1500)    
     $         CHTMP3,
     $         CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $         (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
C
C              Mesh B.C.'s in 1st column of 1st field
               IF (CHTEMP.NE.' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP
C              check for fortran function as denoted by lower case bc's:
            ELSE
               IEL=1
               READ(8,ERR=1500,END=1500) CHTMP3,
     $         CBCS(ISIDE,IEL),ID1,ID2,(BCS(II,ISIDE,IEL),II=1,NBCREA)
C              check for fortran function as denoted by lower case bc's:
            ENDIF
 1080    CONTINUE
 1100 CONTINUE
C
C     END OF BC READ
C
      RETURN
C
C     Error handling:
C
 1500 CONTINUE
      WRITE(6,1501) IFIELD,IEG
 1501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I6
     $    ,/,2X,'(unformatted) ABORTING IN ROUTINE RDBDRY.')
      call exitt
      ENDIF
C
      RETURN
      END

      SUBROUTINE RDICDF
C
C     .Read Initial Conditions / Drive Force
C
C     .Broadcast ICFILE to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      CHARACTER*80 LINE
      LOGICAL      IFGTIL
C
C     Read Initial Conditions/Restart Files
      CALL BLANK(INITC,1200)
      READ(9,80,ERR=200,END=200) LINE
      IF (INDX1(LINE,'RESTART',7).NE.0) THEN
         IF (.NOT.IFGTIL(NSKIP,LINE)) GOTO 200
c        READ(LINE,*,ERR=200,END=200) NSKIP
         DO 50 I=1,NSKIP
            READ(9,80,ERR=200,END=200) INITC(I)
   50    CONTINUE
         READ(9,80,ERR=200,END=200) LINE
      ENDIF
   80 FORMAT(A80)
      IF (.NOT.IFGTIL(NSKIP,LINE)) GOTO 200
c     READ(LINE,*,ERR=200,END=200)NSKIP
      DO 100 I=1,NSKIP
         READ(9,80,ERR=200,END=200) LINE
  100 CONTINUE
C     Read drive force data
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NSKIP
      DO 110 I=1,NSKIP
        READ(9,80,ERR=200,END=200) LINE
  110 CONTINUE
      RETURN
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING INITIAL CONDITION/DRIVE FORCE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDICDF.')
      call exitt
C
      RETURN
      END
      SUBROUTINE RDMATP
C
C     .Read materials property data
C
C     .Disperse material properties to all processors according 
C      to sequential partition scheme
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      CHARACTER*80 LINE
C
      CALL IZERO(MATYPE,16*LDIMT1)
      CALL RZERO(CPGRP ,48*LDIMT1)
C
C     Read material property data
C
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NSKIP
      READ(9,*,ERR=200,END=200) NPACKS
      DO 100 IIG=1,NPACKS
         IFVPS=.TRUE.
         READ(9,*)IGRP,IFLD,ITYPE
         MATYPE(IGRP,IFLD)=ITYPE
         DO 100 IPROP=1,3
            IF(ITYPE.EQ.1) READ(9,* ) CPGRP(IGRP,IFLD,IPROP)
            IF(ITYPE.EQ.2) READ(9,80) LINE
   80 FORMAT(A80)
  100 CONTINUE
C
      RETURN
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING MATERIAL PROPERTIES DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDMATP.')
      call exitt
C
      RETURN
      END
      SUBROUTINE RDHIST
C
C     .Read history data
C
C     .Broadcast to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
C
      ! FIXME
      ! lhis never defined, set to zero. ak, 2006-07-31
      lhis = 0

      CALL BLANK (HCODE ,11*lhis)
      CALL IZERO (LOCHIS, 4*lhis)
C
C     Read history data
      READ (9,*)
      READ (9,*,ERR=200,END=200) NHIS
      if (nhis.gt.lhis) then
         write(6,*) nid,' Too many history pts. RESET LHIS.',nhis,lhis
         call exitt
      endif
c
C     HCODE(10) IS WHETHER IT IS HISTORY, STREAKLINE, PARTICLE, ETC.
      if (nhis.gt.0) then
         do i=1,nhis
            if (nelgt.lt.100000) then
               read(9,130,err=200,end=200)
     $         (hcode(ii,i),ii=1,11),(lochis(i2,i),i2=1,4)
  130          format(1x,11a1,1x,4i5)
            else
               read(9,131,err=200,end=200)
     $         (hcode(ii,i),ii=1,11),(lochis(i2,i),i2=1,4)
  131          format(1x,11a1,1x,3i5,i10)
            endif
c
c           threshold lochis locations to allow easy specification of "NX,NY,NZ"
c           pff 1/7/97
c
            if (hcode(10,i).eq.'H') then
               lochis(1,i) = min(lochis(1,i),nx1)
               lochis(2,i) = min(lochis(2,i),ny1)
               lochis(3,i) = min(lochis(3,i),nz1)
c
c              if lochis_k = -1, set it to nxk/2   pff 8/21/03
c
               if (lochis(1,i).eq.-1) lochis(1,i) = (nx1+1)/2
               if (lochis(2,i).eq.-1) lochis(2,i) = (ny1+1)/2
               if (lochis(3,i).eq.-1) lochis(3,i) = (nz1+1)/2
            endif
         enddo
      endif
C
      return
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING HISTORY DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDHIST.')
      call exitt
C
      RETURN
      END
      SUBROUTINE RDOUT
C
C     .Read output specs
C
C     .Broadcast to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
C
C     Read output specs
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NOUTS
      READ(9,*,ERR=200,END=200) IFXYO
      READ(9,*,ERR=200,END=200) IFVO
      READ(9,*,ERR=200,END=200) IFPO
      READ(9,*,ERR=200,END=200) IFTO
      ! FIXME: Our files only seem to have five of these. Hope this
      ! doesn't matter. Gfortran code (correctly) bails out here.
      ! ak, 8/10/06
      ! READ(9,*,ERR=200,END=200) IFBO   !  IFTGO
      ! READ(9,*,ERR=200,END=200) IPSCO
      !IF (IPSCO.GT.0) THEN
      !   DO I=1,IPSCO
      !      READ(9,*,ERR=200,END=200) IFPSCO(I)
      !   ENDDO
      !ENDIF
C
      RETURN
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING OUTPUT SPECIFICATION DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDOUT.')
      call exitt
C
      RETURN
      END

      SUBROUTINE RDOBJ
C
C     .Read objects
C
C     .Broadcast to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
C
C     Default if no data is read No Objects
C
      NOBJ=0

      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NOBJ
C
      IF(NOBJ.GT.MAXOBJ)THEN
         write(6,*) nid,'ERROR, too many objects:',nobj,maxobj
         call exitt
      ENDIF
C
      DO 10 IOBJ = 1,NOBJ
         READ(9,*,ERR=200,END=200) NMEMBER(IOBJ)
         IF(NMEMBER(IOBJ).GT.MAXMBR)THEN
            PRINT*,'ERROR: Too many members in object ',IOBJ
            call exitt
         ENDIF
         DO 5 MEMBER=1,NMEMBER(IOBJ)
            READ(9,*,ERR=200,END=200) OBJECT(IOBJ,MEMBER,1),
     $                                OBJECT(IOBJ,MEMBER,2)
    5    CONTINUE
   10 CONTINUE
      if (nid.eq.0) write(6,*) nobj,' objects found,'
     $                             ,(nmember(k),k=1,nobj)
C
      RETURN
C
C     Error handling:  For old versions, default to no objects
C
  200 CONTINUE
      NOBJ=0
C
      RETURN
      END


      SUBROUTINE VRDSMSH
C=====================================================================
C     Verify that mesh and dssum are properly defined by performing
C        a direct stiffness operation on the X,Y and Z coordinates.
C     Note that periodic faces are not checked here.
C=====================================================================
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT),TB(LX1,LY1,LZ1,LELT)
     $           ,QMASK(LX1,LY1,LZ1,LELT),tmp(2)
      CHARACTER*3 CB
C
      IERR      = 0
      EPS       = 1.0e-04
      IFIELD    = 1
      IF (IFHEAT) IFIELD = 2
      NXYZ1     = NX1*NY1*NZ1
      NTOT      = NX1*NY1*NZ1*NELT
      NFACES    = 2*NDIM
C
C     First check - use 1/Multiplicity
C
      IF (IFHEAT) THEN 
         CALL COPY(TA,TMULT,NTOT)
      ELSE
         CALL COPY(TA,VMULT,NTOT)
      ENDIF

c     write(6,1) 
c    $(nid,'tab4',lglel(ie,node),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
c   1 format(i3,a4,i3,16f5.2)
c
      CALL DSSUM(TA,NX1,NY1,NZ1)
c
c     write(6,1) 
c    $(nid,'taaf',lglel(ie,node),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
c
      CALL RONE (TB,NTOT)
      CALL SUB2 (TB,TA,NTOT)

      DO 1000 IE=1,NELT
      IEG=LGLEL(IE,NODE)
      DO 1000 IZ=1,NZ1
      DO 1000 IY=1,NY1
      DO 1000 IX=1,NX1
         IF (ABS(TB(IX,IY,IZ,IE)).GT.EPS ) THEN
            WRITE(6,1005) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TA(IX,IY,IZ,IE),eps
c           WRITE(7,1005) IX,IY,IZ,IEG
c    $      ,XM1(IX,IY,IZ,IE),TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE)
c    $      ,QMASK(IX,IY,IZ,IE)
 1005       FORMAT(2X,'WARNING: DSSUM problem at:',/
     $            ,2X,'I,J,K,IE:',4I5,/
     $            ,2X,'Near X =',3G16.8,', d:',2G16.8)
            IERR=4
         ENDIF
 1000 CONTINUE
C
C     Set up QMASK quickly to annihilate checks on periodic bc's
C
      CALL RONE(QMASK,NTOT)
      DO 100 IEL=1,NELT
      DO 100 IFACE=1,NFACES
         CB =CBC(IFACE,IEL,IFIELD)
         IF (CB.EQ.'P  ') 
     $         CALL FACEV(QMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
  100 CONTINUE
      CALL DSOP(QMASK,'MUL',NX1,NY1,NZ1)
C
C     X-component
C
      CALL COPY(TA,XM1,NTOT)
      CALL COPY(TB,XM1,NTOT)
      CALL DSOP(TA,'MIN',NX1,NY1,NZ1)
      CALL DSOP(TB,'MAX',NX1,NY1,NZ1)
      CALL SUB2(TA,XM1,NTOT)
      CALL SUB2(TB,XM1,NTOT)
      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)
      DO 1100 IE=1,NELT
         XSCMAX = VLMAX(XM1(1,1,1,IE),NXYZ1)
         XSCMIN = VLMIN(XM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(XSCMAX-XSCMIN)
         SCAL2=ABS(XSCMAX)
         SCAL3=ABS(XSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         XSCALE = 1./SCAL1
         IEG=LGLEL(IE,NODE)
         DO 1100 IZ=1,NZ1
         DO 1100 IY=1,NY1
         DO 1100 IX=1,NX1
         IF (ABS(TA(IX,IY,IZ,IE)*XSCALE).GT.EPS .OR.
     $       ABS(TB(IX,IY,IZ,IE)*XSCALE).GT.EPS ) THEN
            WRITE(6,1105) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),xscale
 1105       FORMAT(1X,'WARNING: Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',4I5,/
     $            ,1X,'Near X =',3G16.8,', d:',3G16.8)
            IERR=1
         ENDIF
 1100 CONTINUE
C
C     Y-component
C
      CALL COPY(TA,YM1,NTOT)
      CALL COPY(TB,YM1,NTOT)
      CALL DSOP(TA,'MIN',NX1,NY1,NZ1)
      CALL DSOP(TB,'MAX',NX1,NY1,NZ1)
      CALL SUB2(TA,YM1,NTOT)
      CALL SUB2(TB,YM1,NTOT)
      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)
      DO 1200 IE=1,NELT
         YSCMAX = VLMAX(YM1(1,1,1,IE),NXYZ1)
         YSCMIN = VLMIN(YM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(YSCMAX-YSCMIN)
         SCAL2=ABS(YSCMAX)
         SCAL3=ABS(YSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         YSCALE = 1./SCAL1
         IEG=LGLEL(IE,NODE)
         DO 1200 IZ=1,NZ1
         DO 1200 IY=1,NY1
         DO 1200 IX=1,NX1
         IF (ABS(TA(IX,IY,IZ,IE)*YSCALE).GT.EPS .OR.
     $       ABS(TB(IX,IY,IZ,IE)*YSCALE).GT.EPS ) THEN
            WRITE(6,1205) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),yscale
 1205       FORMAT(1X,'WARNING: Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',4I5,/
     $            ,1X,'Near Y =',3G16.8,', d:',3G16.8)
            IERR=2
         ENDIF
 1200 CONTINUE
C
C     Z-component
C
      IF (IF3D) THEN
       CALL COPY(TA,ZM1,NTOT)
       CALL COPY(TB,ZM1,NTOT)
       CALL DSOP(TA,'MIN',NX1,NY1,NZ1)
       CALL DSOP(TB,'MAX',NX1,NY1,NZ1)
       CALL SUB2(TA,ZM1,NTOT)
       CALL SUB2(TB,ZM1,NTOT)
       CALL COL2(TA,QMASK,NTOT)
       CALL COL2(TB,QMASK,NTOT)
       DO 1300 IE=1,NELT
          ZSCMAX = VLMAX(ZM1(1,1,1,IE),NXYZ1)
          ZSCMIN = VLMIN(ZM1(1,1,1,IE),NXYZ1)
          SCAL1=ABS(ZSCMAX-ZSCMIN)
          SCAL2=ABS(ZSCMAX)
          SCAL3=ABS(ZSCMIN)
          SCAL1=MAX(SCAL1,SCAL2)
          SCAL1=MAX(SCAL1,SCAL3)
          ZSCALE = 1./SCAL1
          IEG=LGLEL(IE,NODE)
          DO 1300 IZ=1,NZ1
          DO 1300 IY=1,NY1
          DO 1300 IX=1,NX1
          IF (ABS(TA(IX,IY,IZ,IE)*ZSCALE).GT.EPS .OR.
     $        ABS(TB(IX,IY,IZ,IE)*ZSCALE).GT.EPS ) THEN
           WRITE(6,1305) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),zscale
 1305       FORMAT(1X,'WARNING: Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',4I5,/
     $            ,1X,'Near Z =',3G16.8,', d:',3G16.8)
            IERR=3
          ENDIF
 1300  CONTINUE
      ENDIF
C
      IF (IERR.gt.0) THEN
         WRITE(6,1400) 
 1400    FORMAT
     $   ('  Mesh consistency check failed.  EXITING in VRDSMSH.')
      ENDIF
C
      tmp(1)=ierr
      CALL GOP(tmp,tmp(2),'M  ',1)
      IF (tmp(1).ge.4.0) THEN
         call exitt
      ENDIF
C
      IF (NID.EQ.0) WRITE(6,1500) 
 1500 FORMAT('  Mesh consistency check, OK.')
      RETURN
      END
      SUBROUTINE ROTAT2(XYZ,ANGLE,NPTS)
C
C     Rotate NPTS through ANGLE (in two directions IF3D).
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      DIMENSION XYZ(3,1)
      COMMON /CTMP0/ RMTRX(3,3),RX(3,3),RZ(3,3),XYZN(3,10)
C
      SINA=SIN(ANGLE)
      COSA=COS(ANGLE)
      CALL RZERO(RX,9)
      CALL RZERO(RZ,9)
      RX(1,1)=COSA
      RX(2,2)=COSA
      RX(1,2)=SINA
      RX(2,1)=-SINA
      RX(3,3)=1.0
      IF (IF3D) THEN
         RZ(1,1)=COSA
         RZ(3,3)=COSA
         RZ(1,3)=SINA
         RZ(3,1)=-SINA
         RZ(2,2)=1.0
      ELSE
         RZ(1,1)=1.0
         RZ(2,2)=1.0
         RZ(3,3)=1.0
      ENDIF
      CALL MXM(RX,3,RZ,3,RMTRX,3)
C
C     Strip mine mxms in chunks of 10:
      DO 100 I=1,NPTS-10,10
         CALL MXM(RMTRX,3,XYZ(1,I),3,XYZN,10)
         CALL COPY(XYZ(1,I),XYZN,30)
  100 CONTINUE
      N10=MOD1(NPTS,10)
      I=NPTS-N10+1
      CALL RZERO(XYZN,30)
      IF (N10.GT.0) THEN
         CALL MXM(RMTRX,3,XYZ(1,I),3,XYZN,N10)
         CALL COPY(XYZ(1,I),XYZN,3*N10)
      ENDIF
C
      RETURN
      END

      SUBROUTINE SCALE(XYZL,NL)
C
C     Rescale XYZL such that the mean value of IXX=IYY=IZZ for each element.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      DIMENSION XYZL(3,8,LELT)
      COMMON /CTMP0/ VO(LELT),XYZI(3,LELT),CG(3,LELT)
     $              ,TI(6),WORK(6)
C
C     Compute volumes -
C
      CALL VOLUME2(VO,XYZL,NL)
      VTOT=GLSUM (VO,NL)
C
C     Compute (weighted) average inertia for each element.
C
      NCRNR=2**NDIM
      CALL RZERO(TI,6)
      DO 100 IL=1,NL
         VO0 = VO(IL)/VTOT
         CALL INRTIA(XYZI(1,IL),CG(1,IL),XYZL(1,1,IL),NCRNR,1)
         TI(1)=TI(1)+XYZI(1,IL)*VO0
         TI(2)=TI(2)+XYZI(2,IL)*VO0
         TI(3)=TI(3)+XYZI(3,IL)*VO0
         TI(4)=TI(4)+CG(1,IL)  *VO0
         TI(5)=TI(5)+CG(2,IL)  *VO0
         TI(6)=TI(6)+CG(3,IL)  *VO0
  100 CONTINUE
      CALL GOP(TI,WORK,'+  ',6)
      XI  =SQRT(TI(1))
      YI  =SQRT(TI(2))
      ZI  =1.0
      IF (IF3D) ZI=SQRT(TI(3))
C
C     Rescale ( & shift to a nearly mean zero )
C
      DO 200 IL=1,NL
      DO 200 IC=1,NCRNR
         XYZL(1,IC,IL)=(XYZL(1,IC,IL)-TI(4))/XI
         XYZL(2,IC,IL)=(XYZL(2,IC,IL)-TI(5))/YI
         XYZL(3,IC,IL)=(XYZL(3,IC,IL)-TI(6))/ZI
  200 CONTINUE
C
      RETURN
      END
      SUBROUTINE INRTIA(XYZI,CG,XYZL,N,ITYPE)
C
C     Compute cg and inertia for a collection of unit point masses.
C     This is a global (multiprocessor) operation, only IF itype=2.
C
      DIMENSION XYZI(3),CG(3),XYZL(3,1)
      DIMENSION TI(4),WORK(4)
C
      TI(1)=0.0
      TI(2)=0.0
      TI(3)=0.0
      TI(4)=N
      DO 100 I=1,N
         TI(1)=TI(1)+XYZL(1,I)
         TI(2)=TI(2)+XYZL(2,I)
         TI(3)=TI(3)+XYZL(3,I)
  100 CONTINUE
      IF (ITYPE.EQ.2) CALL GOP(TI,WORK,'+  ',4)
      IF (TI(4).EQ.0.0) TI(4)=1.0
      CG(1)=TI(1)/TI(4)
      CG(2)=TI(2)/TI(4)
      CG(3)=TI(3)/TI(4)
C
      TI(1)=0.0
      TI(2)=0.0
      TI(3)=0.0
      DO 200 I=1,N
         TI(1)=TI(1)+( XYZL(1,I)-CG(1) )**2
         TI(2)=TI(2)+( XYZL(2,I)-CG(2) )**2
         TI(3)=TI(3)+( XYZL(3,I)-CG(3) )**2
  200 CONTINUE
      IF (ITYPE.EQ.2) CALL GOP(TI,WORK,'+  ',3)
      TI(1)=TI(1)/TI(4)
      TI(2)=TI(2)/TI(4)
      TI(3)=TI(3)/TI(4)
      IF (ITYPE.EQ.2) THEN
C        std. def'n of inertia.
         XYZI(1)=TI(2)+TI(3)
         XYZI(2)=TI(3)+TI(1)
         XYZI(3)=TI(1)+TI(2)
      ELSE
         XYZI(1)=TI(1)
         XYZI(2)=TI(2)
         XYZI(3)=TI(3)
      ENDIF
C
      RETURN
      END
      SUBROUTINE VOLUME2(VOL,XYZ,N)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      DIMENSION XYZ(3,2,2,2,1)
      DIMENSION VOL(1)
C
      DO 1000 IE=1,N
         VOL(IE)=0.0
         IF (IF3D) THEN        
           DO 20 K=1,2
           DO 20 J=1,2
           DO 20 I=1,2
              VOL1 = (XYZ(1,2,J,K,IE)-XYZ(1,1,J,K,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,I,J,2,IE)-XYZ(3,I,J,1,IE))
              VOL2 = (XYZ(1,2,J,K,IE)-XYZ(1,1,J,K,IE))
     $             * (XYZ(2,I,J,2,IE)-XYZ(2,I,J,1,IE))
     $             * (XYZ(3,I,2,K,IE)-XYZ(3,I,1,K,IE))
              VOL3 = (XYZ(1,I,2,K,IE)-XYZ(1,I,1,K,IE))
     $             * (XYZ(2,2,J,K,IE)-XYZ(2,1,J,K,IE))
     $             * (XYZ(3,I,J,2,IE)-XYZ(3,I,J,1,IE))
              VOL4 = (XYZ(1,I,J,2,IE)-XYZ(1,I,J,1,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,I,2,K,IE)-XYZ(3,I,1,K,IE))
              VOL5 = (XYZ(1,I,2,K,IE)-XYZ(1,I,1,K,IE))
     $             * (XYZ(2,I,J,2,IE)-XYZ(2,I,J,1,IE))
     $             * (XYZ(3,2,J,K,IE)-XYZ(3,1,J,K,IE))
              VOL6 = (XYZ(1,I,J,2,IE)-XYZ(1,I,J,1,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,2,J,K,IE)-XYZ(3,1,J,K,IE))
              VOL(IE) = VOL(IE)+VOL1+VOL2+VOL3+VOL4+VOL5+VOL6
   20      CONTINUE
           VOL(IE)=VOL(IE)/8.0
         ELSE
C     2-D:
            DO 40 J=1,2
            DO 40 I=1,2
              VOL1 = (XYZ(1,2,J,1,IE)-XYZ(1,1,J,1,IE))
     $             * (XYZ(2,I,2,1,IE)-XYZ(2,I,1,1,IE))
              VOL3 = (XYZ(1,I,2,1,IE)-XYZ(1,I,1,1,IE))
     $             * (XYZ(2,2,J,1,IE)-XYZ(2,1,J,1,IE))
              VOL(IE)=VOL(IE)+VOL1+VOL3
   40      CONTINUE
           VOL(IE)=VOL(IE)/4.0
         ENDIF
         VOL(IE)=ABS(VOL(IE))
 1000 CONTINUE
C
      RETURN
      END
      SUBROUTINE FINDCG(CG,XYZ,N)
C
C     Compute cg for N elements.
C
      INCLUDE 'SIZE'
      DIMENSION CG(3,1),XYZ(3,8,1)
C
      NCRNR=2**NDIM
      CALL RZERO(CG,3*N)
      DO 100 I =1,N
      DO 100 IC=1,NCRNR
         CG(1,I)=CG(1,I)+XYZ(1,IC,I)
         CG(2,I)=CG(2,I)+XYZ(2,IC,I)
         CG(3,I)=CG(3,I)+XYZ(3,IC,I)
  100 CONTINUE
      TMP=1.0/(NCRNR)
      CALL CMULT(CG,TMP,3*N)
      RETURN
      END
      SUBROUTINE DIVIDE(LIST1,LIST2,NL1,NL2,IFOK,LIST,NL,XYZI,CG,wgt)
C
C     Divide the elements associated with this subdomain according to
C     the direction having the smallest moment of inertia (the "long"
C     direction).
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
C
      DIMENSION LIST(LELT),LIST1(LELT),LIST2(LELT)
      DIMENSION XYZI(3),CG(3,LELT),wgt(1)
      COMMON /CTMP0/ XCG(LELT),YCG(LELT),ZCG(LELT)
      REAL IXX,IYY,IZZ
      INTEGER WORK(2),WRK2(2)
      LOGICAL IFOK
C
C     Choose "long" direction:
C
      IXX=XYZI(1)
      IYY=XYZI(2)
      IZZ=XYZI(3)
      IF (IF3D) THEN
         IF (IXX.LE.IYY.AND.IXX.LE.IZZ) THEN
            DO 104 IE=1,NL
               XCG(IE)=CG(1,IE)
               YCG(IE)=CG(2,IE)
               ZCG(IE)=CG(3,IE)
  104       CONTINUE
         ELSEIF (IYY.LE.IXX.AND.IYY.LE.IZZ) THEN
            DO 106 IE=1,NL
               XCG(IE)=CG(2,IE)
               YCG(IE)=CG(3,IE)
               ZCG(IE)=CG(1,IE)
  106       CONTINUE
         ELSEIF (IZZ.LE.IXX.AND.IZZ.LE.IYY) THEN
            DO 108 IE=1,NL
               XCG(IE)=CG(3,IE)
               YCG(IE)=CG(1,IE)
               ZCG(IE)=CG(2,IE)
  108       CONTINUE
         ENDIF
      ELSE
C     2-D:
         IF (IXX.LE.IYY) THEN
            DO 114 IE=1,NL
               XCG(IE)=CG(1,IE)
               YCG(IE)=CG(2,IE)
  114       CONTINUE
         ELSE
            DO 116 IE=1,NL
               XCG(IE)=CG(2,IE)
               YCG(IE)=CG(1,IE)
  116       CONTINUE
         ENDIF
      ENDIF
      call col2(xcg,wgt,nl)
      call col2(ycg,wgt,nl)
      call col2(zcg,wgt,nl)
C
C     Find median value of CG to determine dividing point:
C
      XM=FMDIAN(XCG,NL,IFOK)
      YM=FMDIAN(YCG,NL,IFOK)
      ZM=0.0
      IF (IF3D) ZM=FMDIAN(ZCG,NL,IFOK)
C
C     Diagnostics
C
      IF (.NOT.IFOK) THEN
         WRITE(6,130) NID,NL,XM,YM,ZM
         DO 120 IL=1,NL
            WRITE(6,135) NID,IL,XCG(IL),YCG(IL),ZCG(IL)
  120    CONTINUE
  130    FORMAT(I3,'DIVIDE: NL,XM,YM,ZM',I3,3F12.5)
  135    FORMAT(I3,'DIVIDE: NID,IL,XC,YC,ZCG',I4,3F12.5)
      ENDIF
C
C=============================================================
C     Divide LIST into LIST1 (XCG < XM) and LIST2 (XCG>XM).
C=============================================================
C
      NL1=0
      NL2=0
      DO 200 IE=1,NL
         IF (XCG(IE).LT.XM) THEN
            NL1=NL1+1
            LIST1(NL1)=LIST(IE)
         ENDIF
         IF (XCG(IE).GT.XM) THEN
            NL2=NL2+1
            LIST2(NL2)=LIST(IE)
         ENDIF
         IF (XCG(IE).EQ.XM) THEN
C
C           We have to look at the other directions to arrive at
C           a unique subdivision algortithm.
C
C
C           More Diagnostics
C
            IF (.NOT.IFOK) WRITE(6,201) NID,IE,XCG(IE),XM
  201    FORMAT(I3,'DIVIDE: IE,XCG,XM:',I4,3F12.5)
C
            IF (YCG(IE).LT.YM) THEN
               NL1=NL1+1
               LIST1(NL1)=LIST(IE)
            ENDIF
            IF (YCG(IE).GT.YM) THEN
               NL2=NL2+1
               LIST2(NL2)=LIST(IE)
            ENDIF
            IF (YCG(IE).EQ.YM) THEN
C              look at 3rd direction.  
               IF (IF3D .AND. ZCG(IE).LT.ZM) THEN
                  NL1=NL1+1
                  LIST1(NL1)=LIST(IE)
               ELSE IF (IF3D .AND. ZCG(IE).GT.ZM) THEN
                  NL2=NL2+1
                  LIST2(NL2)=LIST(IE)
               ELSE 
C                 for 2- or 3-D intdeterminate case:
                  NL1=NL1+1
                  LIST1(NL1)=LIST(IE)
               ENDIF
            ENDIF
C
         ENDIF
  200 CONTINUE
C
C     Check for an even distribution (i.e. - not different by
C     more than 1):
C
      IFOK=.TRUE.  
      WORK(1)=NL1 
      WORK(2)=NL2
      CALL IGOP(WORK,WRK2,'+  ',2)
      IF (ABS(WORK(1)-WORK(2)).GT.1) IFOK=.FALSE.
C
      RETURN
      END
      SUBROUTINE GRAY(IPG,IPG1,ITMP,N)
      DIMENSION IPG(1),IPG1(1),ITMP(1)
      INTEGER CD
C
      CD=LOG2(N)
      NPR = 1
      IPG(1) = 1
      DO 100 I=1,CD
         N1  = NPR+1
         N2  = NPR*2
         CALL ICOPY(IPG(N1),IPG(1),NPR)
         CALL IADD (IPG(N1),NPR,NPR)
         CALL IFLIP(IPG(N1),NPR)
         NPR = 2*NPR
  100 CONTINUE
      N1=NPR+1
c     IPG(0)=IPG(NPR)
c     IPG(N1)=IPG(1)
C
C     Compute the inverse of the gray code mapping.  In other words,
C     given the ith entry on the gray code list, IPG(i) returns the
C     appropriate Node number (=nid+1).  Given a node number, NODE,
C     IPG1(NODE) returns the location in the gray code (i.e. ring) map.
C
      CALL ICOPY(ITMP,IPG ,NPR)
      CALL ISORT(ITMP,IPG1,NPR)
C
      RETURN
      END
c-----------------------------------------------------------------------
