!> \file bdry.F   
!! 
!! \brief brief description of the file 
!! 
!! boundary                                     


!> \brief Brief description of the subroutine 
!! 
!! A more detailed descripton goes here. This is an example 
!! of a detailed description of a method
!!
!! @param param1 a description of the first parameter 
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
c----------------------------------------------------------------------
      subroutine SETLOG           
c----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'TURBO'
      INCLUDE 'EMWAVE'
                  
      COMMON  /CPRINT/ IFPRINT
      COMMON  /nekcb/ cb
      CHARACTER CB*3
      LOGICAL  IFALGN,IFNORX,IFNORY,IFNORZ,IFPRINT
                
      NFACE  = 2*NDIM
      NMXV   = NFACE*NELV
      NMXT   = NFACE*NELT
               
      IFPRINT= .TRUE.
      IFVCOR = .TRUE.
      IFGEOM = .FALSE.
      IFINTQ = .FALSE.
      IFSURT = .FALSE.
      IFWCNO = .FALSE.
      IFSWALL= .FALSE.
              
      DO 10 IFIELD=1,NFIELD
            IFNONL(IFIELD) = .FALSE.
 10   CONTINUE
C
C     CALL LFALSE (IFEPPM,NMXV)
c     CALL LFALSE (IFQINP,NMXV)
C
C     IF (IFMODEL) CALL SETSHL
C
      IF (IFMVBD) THEN
         IFGEOM = .TRUE.
         IF   ( IFFLOW .AND. .NOT.IFNAV  )  IFWCNO = .TRUE.
      ENDIF

      IF (IFCEM) THEN
C
      IFIELD = 2
        DO IEL=1,NELV                
          DO IFC=1,NFACE
            CB=CBC(IFC,IEL,IFIELD)
            IF      (CB.EQ.'PEC' .OR. CB.EQ.'pec') THEN
                    IFPEC = .TRUE.
            ENDIF
            IF  (CB.EQ.'PMC' .OR. CB.EQ.'pmc') THEN
                    IFPMC = .TRUE.
            ENDIF
            IF  (CB.EQ.'PML' .OR. CB.EQ.'pml') THEN
                    IFPML = .TRUE.
            ENDIF
            IF  (CB.EQ.'DTN' .OR. CB.EQ.'dtn') THEN
                    IFDTN = .TRUE.
            ENDIF
            IF  (CB.EQ.'N  ' .OR. CB.EQ.'n  ') THEN
                    IFNEU = .TRUE.
            ENDIF

          ENDDO
        ENDDO
C
      ENDIF
     
      IF (IFFLOW) THEN
         IFIELD = 1
         DO 100 IEL=1,NELV
         DO 100 IFC=1,NFACE
            CB = CBC(IFC,IEL,IFIELD)
C           CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
C           IF ( .NOT.IFSTRS ) CALL CHKCBC  (CB,IEL,IFC,IFALGN)
            IF  (CB.EQ.'O  ' .OR. CB.EQ.'o  ' .OR.
     $           CB.EQ.'ON ' .OR. CB.EQ.'on ' .OR.
     $           CB.EQ.'S  ' .OR. CB.EQ.'s  ' .OR.
     $           CB.EQ.'SL ' .OR. CB.EQ.'sl ' .OR.
     $           CB.EQ.'MS ' .OR. CB.EQ.'ms ')  THEN
                                              IFVCOR          = .FALSE.
                                              IFEPPM(IFC,IEL) = .TRUE.
            ENDIF
            IF  (CB.EQ.'VL ' .OR. CB.EQ.'vl ' .OR.
     $           CB.EQ.'WSL' .OR. CB.EQ.'wsl' .OR.
     $           CB.EQ.'SL ' .OR. CB.EQ.'sl ' .OR.
     $           CB.EQ.'SHL' .OR. CB.EQ.'shl' .OR.
     $           CB.EQ.'MS ' .OR. CB.EQ.'ms ' .OR.
     $           CB.EQ.'O  ' .OR. CB.EQ.'o  ' .OR.
     $           CB.EQ.'ON ' .OR. CB.EQ.'on ')  THEN
                                              IFQINP(IFC,IEL) = .TRUE.
            ENDIF
            IF  (CB.EQ.'MS ' .OR. CB.EQ.'ms ' .OR.
     $           CB.EQ.'MSI' .OR. CB.EQ.'msi' ) THEN
                                              IFSURT          = .TRUE.
            ENDIF
            IF  (CB.EQ.'WS ' .OR. CB.EQ.'ws ' .OR.
     $           CB.EQ.'WSL' .OR. CB.EQ.'wsl') THEN
                                              IFSWALL         = .TRUE.
                                              IFCWUZ          = .TRUE.
            ENDIF
  100    CONTINUE
      ENDIF
C
      IF (IFHEAT) THEN
C
         DO 250 IFIELD=2,NFIELD
         DO 250 IEL=1,NELFLD(IFIELD)
         DO 250 IFC=1,NFACE
            CB=CBC(IFC,IEL,IFIELD)
            IF  (CB.EQ.'P  ' .OR. CB.EQ.'P  ') THEN
                                              IFNONL(IFIELD)  = .TRUE.
            ENDIF
  250    CONTINUE
C
      ENDIF
C
      IF (NHIS.GT.0) THEN
         IQ = 0
         DO 300 IH=1,NHIS
            IF ( HCODE(10,IH) .EQ. 'I' ) THEN
               IFINTQ = .TRUE.
               IOBJ   = LOCHIS(1,IH)
               IQ     = IQ + 1
               IF (IOBJ.GT.NOBJ .OR. IOBJ.LT.0)  THEN
                  WRITE (6,*) 
     $            'ERROR : Undefined Object for integral',IQ
                  STOP
               ENDIF
            ENDIF
  300    CONTINUE
      ENDIF
C
C     Establish global consistency of LOGICALS amongst all processors.
C
      CALL GLLOG(IFVCOR , .FALSE.)
      CALL GLLOG(IFSURT , .TRUE. )
      CALL GLLOG(IFSWALL, .TRUE. )
      CALL GLLOG(IFCWUZ , .TRUE. )
      CALL GLLOG(IFWCNO , .TRUE. )
      CALL GLLOG(IFDTN  , .TRUE. )
      CALL GLLOG(IFNEU  , .TRUE. )
      CALL GLLOG(IFPEC  , .TRUE. )
      CALL GLLOG(IFPML  , .TRUE. )
      CALL GLLOG(IFPMC  , .TRUE. )

      DO 400 IFIELD=2,NFIELD
         CALL GLLOG(IFNONL(IFIELD),.TRUE.)
  400 CONTINUE
C
      IF (NID.EQ.0) THEN
         WRITE (6,*) 'IFSTRS =',IFSTRS
         WRITE (6,*) 'IFCHAR =',IFCHAR
         WRITE (6,*) 'IFAXIS =',IFAXIS
         WRITE (6,*) 'IFMVBD =',IFMVBD
         WRITE (6,*) 'IFMELT =',IFMELT
         WRITE (6,*) 'IFMODEL=',IFMODEL
         WRITE (6,*) 'IFKEPS =',IFKEPS
         WRITE (6,*) '  '
         WRITE (6,*) 'IFVCOR =',IFVCOR
         WRITE (6,*) 'IFINTQ =',IFINTQ
         WRITE (6,*) 'IFCWUZ =',IFCWUZ
         WRITE (6,*) 'IFSWALL=',IFSWALL
         WRITE (6,*) 'IFGEOM =',IFGEOM
         WRITE (6,*) 'IFSURT =',IFSURT
         WRITE (6,*) 'IFWCNO =',IFWCNO
         WRITE (6,*) 'IFCEM  =',IFCEM     
         DO 500 IFIELD=1,NFIELD
            WRITE (6,*) '  '
            WRITE (6,*) 'IFADVC for field',IFIELD,'   = ',IFADVC(IFIELD)
            WRITE (6,*) 'IFNONL for field',IFIELD,'   = ',IFNONL(IFIELD)
 500     CONTINUE
      ENDIF
C
      RETURN
      END
C
