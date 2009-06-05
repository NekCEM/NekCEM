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
      CALL LFALSE (IFEPPM,NMXV)
      CALL LFALSE (IFQINP,NMXV)
C
      IF (IFMODEL) CALL SETSHL
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

          ENDDO
        ENDDO
C
      ENDIF
     
      IF (IFFLOW) THEN
         IFIELD = 1
         DO 100 IEL=1,NELV
         DO 100 IFC=1,NFACE
            CB = CBC(IFC,IEL,IFIELD)
            CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
            IF ( .NOT.IFSTRS ) CALL CHKCBC  (CB,IEL,IFC,IFALGN)
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
      SUBROUTINE SETRZER
C-------------------------------------------------------------------
C
C     Check for axisymmetric case.
C     Are some of the elements close to the axis?
C
C-------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
C
C     Single or double precision???
C
      DELTA = 1.E-9
      X     = 1.+DELTA
      Y     = 1.
      DIFF  = ABS(X-Y)
      IF (DIFF.EQ.0.) EPS = 1.E-7
      IF (DIFF.GT.0.) EPS = 1.E-14
C
      DO 100 IEL=1,NELT
         IFRZER(IEL) = .FALSE.
         IF (IFAXIS) THEN
            NVERT = 0
            DO 10 IC=1,4
               IF(ABS(YC(IC,IEL)).LT.EPS) NVERT = NVERT+1
 10         CONTINUE
         ENDIF
         IEDGE = 1
         IF ((NVERT.EQ.2).AND.(CCURVE(IEDGE,IEL).EQ.' '))
     $       IFRZER(IEL) = .TRUE.
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
C
C      Check direction of normal of an element face for
C      alignment with the X, Y, or Z axis.
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
C
      LOGICAL IFALGN,IFNORX,IFNORY,IFNORZ
C
      SUMX    = 0.0
      SUMY    = 0.0
      SUMZ    = 0.0
      TOLNOR  = 1.0e-3
      IFALGN  = .FALSE.
      IFNORX  = .FALSE.
      IFNORY  = .FALSE.
      IFNORZ  = .FALSE.
C
      IF (NDIM.EQ.2) THEN
C
         NCPF = NX1
         DO 100 IX=1,NX1
            SUMX = SUMX + ABS( ABS(UNX(IX,1,IFC,IEL)) - 1.0 )
            SUMY = SUMY + ABS( ABS(UNY(IX,1,IFC,IEL)) - 1.0 )
  100    CONTINUE
         SUMX = SUMX / NCPF
         SUMY = SUMY / NCPF
         IF ( SUMX.LT.TOLNOR ) THEN
            IFNORX  = .TRUE.
            IFALGN = .TRUE.
         ENDIF
         IF ( SUMY.LT.TOLNOR ) THEN
            IFNORY  = .TRUE.
            IFALGN = .TRUE.
         ENDIF
C
      ELSE
C
         NCPF = NX1*NX1
         DO 200 IX=1,NX1
         DO 200 IY=1,NY1
            SUMX = SUMX + ABS( ABS(UNX(IX,IY,IFC,IEL)) - 1.0 )
            SUMY = SUMY + ABS( ABS(UNY(IX,IY,IFC,IEL)) - 1.0 )
            SUMZ = SUMZ + ABS( ABS(UNZ(IX,IY,IFC,IEL)) - 1.0 )
  200    CONTINUE
         SUMX = SUMX / NCPF
         SUMY = SUMY / NCPF
         SUMZ = SUMZ / NCPF
         IF ( SUMX.LT.TOLNOR ) THEN
            IFNORX  = .TRUE.
            IFALGN = .TRUE.
         ENDIF
         IF ( SUMY.LT.TOLNOR ) THEN
            IFNORY  = .TRUE.
            IFALGN = .TRUE.
         ENDIF
         IF ( SUMZ.LT.TOLNOR ) THEN
            IFNORZ  = .TRUE.
            IFALGN = .TRUE.
         ENDIF
C
      ENDIF
C
      RETURN
      END


      SUBROUTINE CHKAXCB
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      CHARACTER CB*3
C
      IFLD  = 1
      NFACE = 2*NDIM
C
      DO 100 IEL=1,NELV
      DO 100 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFLD)
         IF  (CB.EQ.'A  ' .AND. IFC.NE.1)  GOTO 9000
  100 CONTINUE
C
      RETURN
C
 9000 WRITE (6,*) ' Element face on the axis of symmetry must be FACE 1'
      WRITE (6,*) ' Element',IEL,'   face',IFC,'  is on the axis.'
      STOP
C
      END
      SUBROUTINE CHKCBC (CB,IEL,IFC,IFALGN)
C
C     Check for illegal boundary conditions
C
      CHARACTER CB*3
      LOGICAL IFALGN
C
C     Laplacian formulation only
C
      IF  (CB.EQ.'SH ' .OR.  CB.EQ.'sh ' .OR.
     $     CB.EQ.'SHL' .OR.  CB.EQ.'shl' .OR.
     $     CB.EQ.'S  ' .OR.  CB.EQ.'s  ' .OR.
     $     CB.EQ.'SL ' .OR.  CB.EQ.'sl ' .OR.
     $     CB.EQ.'MS ' .OR.  CB.EQ.'ms ' .OR.
     $     CB.EQ.'MSI' .OR.  CB.EQ.'msi'    )                GOTO 9001
      IF ( .NOT.IFALGN .AND.
     $    (CB.EQ.'ON ' .OR.  CB.EQ.'on ' .OR. CB.EQ.'SYM') ) GOTO 9010
      RETURN
C
 9001 WRITE (6,*) ' Illegal traction boundary conditions detected for'
      GOTO 9999
 9010 WRITE (6,*) ' Mixed B.C. on a side nonaligned with either the X,Y,
     $ or Z axis detected for'
 9999 WRITE (6,*) ' Element',IEL,'   side',IFC,'.'
      WRITE (6,*) ' Selected option only allowed for STRESS FORMULATION'
      WRITE (6,*) ' Execution terminates'
      STOP
      END

      SUBROUTINE BCMASK
C
C     Zero out masks corresponding to Dirichlet boundary points.
C
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'MVGEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
C
      common  /nekcb/ cb
      CHARACTER CB*3
      LOGICAL IFALGN,IFNORX,IFNORY,IFNORZ
C
      NFACES=2*NDIM
      NXYZ  =NX1*NY1*NZ1
C
C     Masks for moving mesh
C
c     IF (IFMVBD) THEN
c        IFIELD = 0
c        CALL STSMASK (W1MASK,W2MASK,W3MASK)
c     ENDIF
C
C     Masks for flow variables
C
      IF (IFFLOW) THEN
         IFIELD = 1
         NEL    = NELFLD(IFIELD)
         NTOT   = NXYZ*NEL
C
C        Pressure mask
C
         CALL RONE(PMASK,NTOT)
         DO 50 IEL=1,NELV
         DO 50 IFACE=1,NFACES
            CB=CBC(IFACE,IEL,IFIELD)
            IF (CB.EQ.'O  ' .OR. CB.EQ.'ON ')
     $         CALL FACEV(PMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
   50    CONTINUE
C
C        Zero out mask at Neumann-Dirichlet interfaces
C
         CALL DSOP(PMASK,'MUL',NX1,NY1,NZ1)
C
C        Velocity masks
C
         IF (IFSTRS) THEN
c          CALL STSMASK (V1MASK,V2MASK,V3MASK)
         ELSE
C
           CALL RONE(V1MASK,NTOT)
           CALL RONE(V2MASK,NTOT)
           CALL RONE(V3MASK,NTOT)
C
           DO 100 IEL=1,NELV
           DO 100 IFACE=1,NFACES
              CB =CBC(IFACE,IEL,IFIELD)
              CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFACE,IEL)
C
C            All-Dirichlet boundary conditions
C
           IF (CB.EQ.'v  ' .OR. CB.EQ.'V  ' .OR. CB.EQ.'vl ' .OR.
     $       CB.EQ.'VL ' .OR. CB.EQ.'W  ') THEN
             CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             GOTO 100
         ENDIF
          
C
C        Mixed-Dirichlet-Neumann boundary conditions
C
         IF (CB.EQ.'SYM') THEN
             IF ( .NOT.IFALGN .OR. IFNORX )
     $            CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             IF ( IFNORY )
     $            CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             IF ( IFNORZ )
     $            CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             GOTO 100
         ENDIF
         IF (CB.EQ.'ON ') THEN
             IF ( IFNORY .OR. IFNORZ )
     $            CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             IF ( .NOT.IFALGN .OR. IFNORX .OR. IFNORZ )
     $            CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             IF ( .NOT.IFALGN .OR. IFNORX .OR. IFNORY )
     $            CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             GOTO 100
         ENDIF
         IF (CB.EQ.'A  ') THEN
             CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
         ENDIF
  100    CONTINUE
C
         CALL DSOP(V1MASK,'MUL',NX1,NY1,NZ1)
         CALL DSOP(V2MASK,'MUL',NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSOP(V3MASK,'MUL',NX1,NY1,NZ1)
C
       ENDIF
C
      ENDIF
C
C     Masks for passive scalars +
C     k and e if k-e turbulence modem:
C     k = nfield-1
C     e = nfield
C
      IF (IFHEAT) THEN
C
         DO 1200 IFIELD=2,NFIELD
            IPSCAL = IFIELD-1
            NEL    = NELFLD(IFIELD)
            NTOT   = NXYZ*NEL
            CALL RONE (TMASK(1,1,1,1,IPSCAL),NTOT)
         DO 1100 IEL=1,NEL
         DO 1100 IFACE=1,NFACES
            CB =CBC(IFACE,IEL,IFIELD)
C
C           Assign mask values.
C
            IF  (CB.EQ.'T  ' .OR. CB.EQ.'t  ' .OR. 
     $           CB.EQ.'MCI' .OR. CB.EQ.'MLI' .OR.
     $           CB.EQ.'KD ' .OR. CB.EQ.'kd ' .OR.
     $           CB.EQ.'ED ' .OR. CB.EQ.'ed ' .OR.
     $           CB.EQ.'KW ' .OR. CB.EQ.'KWS' .OR. CB.EQ.'EWS')
     $           CALL FACEV (TMASK(1,1,1,1,IPSCAL),
     $                       IEL,IFACE,0.0,NX1,NY1,NZ1)
 1100       CONTINUE
         CALL DSOP (TMASK(1,1,1,1,IPSCAL),'MUL',NX1,NY1,NZ1)
 1200    CONTINUE
C
      ENDIF
C
      RETURN
      END


      SUBROUTINE BCPEC  (V1,V2,V3)
C
C     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3).
C     Use IFIELD as a guide to which boundary conditions are to be applied.
C
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
      INCLUDE 'CTIMER'
      INCLUDE 'EMWAVE'

      COMMON /SCRUZ/ TMP1(LX1,LY1,LZ1,LELV)
     $             , TMP2(LX1,LY1,LZ1,LELV)
     $             , TMP3(LX1,LY1,LZ1,LELV)
      COMMON /SCRMG/ TMQ1(LX1,LY1,LZ1,LELV)
     $             , TMQ2(LX1,LY1,LZ1,LELV)
     $             , TMQ3(LX1,LY1,LZ1,LELV)
C
      DIMENSION V1(NX1,NY1,NZ1,LELV),V2(NX1,NY1,NZ1,LELV)
     $         ,V3(NX1,NY1,NZ1,LELV)
c
      common  /nekcb/ cb
      character cb*3
c
      logical ifonbc
c
      ifonbc = .false.
c
      if (icalld.eq.0) then
         tusbc=0.0
         nusbc=0
         icalld=icalld+1
      endif
      nusbc=nusbc+1
      etime1=dclock()
C
C
      NFACES=2*NDIM
      NXYZ  =NX1*NY1*NZ1
      NEL   =NELFLD(IFIELD)
      NTOT  =NXYZ*NEL
C
      CALL RZERO(TMP1,NTOT)
      CALL RZERO(TMP2,NTOT)

      IF (IF3D) CALL RZERO(TMP3,NTOT)
C
C     Velocity boundary conditions
C
      DO 2100 ISWEEP=1,2
         DO 2000 IE=1,NEL
         DO 2000 IFACE=1,NFACES
            CB  = CBC(IFACE,IE,IFIELD)
            BC1 = BC(1,IFACE,IE,IFIELD)
            BC2 = BC(2,IFACE,IE,IFIELD)
            BC3 = BC(3,IFACE,IE,IFIELD)

            IF (CB.EQ.'P  ') THEN   ! 2/24/06  mmin

               CALL FACEV (TMP1,IE,IFACE,BC1,NX1,NY1,NZ1)
               CALL FACEV (TMP2,IE,IFACE,BC2,NX1,NY1,NZ1)

               IF (IF3D) CALL FACEV (TMP3,IE,IFACE,BC3,NX1,NY1,NZ1)
               IF ( IFQINP(IFACE,IE) )
     $         CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                       TMP3(1,1,1,IE),IE,IFACE)
            ENDIF

            IF (CB.EQ.'v  ' .OR. CB.EQ.'vl ' .OR. 
     $          CB.EQ.'ws ' .OR. CB.EQ.'wsl' .OR.
     $          CB.EQ.'mv ' .OR. CB.EQ.'mvn') THEN
                CALL FACEIV (CB,TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                       TMP3(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
                IF ( IFQINP(IFACE,IE) )
     $          CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                        TMP3(1,1,1,IE),IE,IFACE)
            ENDIF

            IF (CB.EQ.'ON ' .OR. CB.EQ.'on ') then   ! 5/21/01 pff
                ifonbc =.true.
                CALL FACEIV ('v  ',TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                       TMP3(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
            ENDIF

 2000    CONTINUE

         DO 2010 IE=1,NEL
         DO 2010 IFACE=1,NFACES
            IF (CBC(IFACE,IE,IFIELD).EQ.'W  ') THEN
               CALL FACEV (TMP1,IE,IFACE,0.0,NX1,NY1,NZ1)
               CALL FACEV (TMP2,IE,IFACE,0.0,NX1,NY1,NZ1)
               IF (IF3D) CALL FACEV (TMP3,IE,IFACE,0.0,NX1,NY1,NZ1)
            ENDIF
 2010    CONTINUE
 2100 CONTINUE

      RETURN
      END


      SUBROUTINE BCDIRSC(S)
C
C     Apply Dirichlet boundary conditions to surface of scalar, S.
C     Use IFIELD as a guide to which boundary conditions are to be applied.
C
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TOPOL'
C
      DIMENSION S(LX1,LY1,LZ1,LELT)
      COMMON /SCRSF/ TMP(LX1,LY1,LZ1,LELT)
     $             , TMA(LX1,LY1,LZ1,LELT)
     $             , SMU(LX1,LY1,LZ1,LELT)
      common  /nekcb/ cb
      CHARACTER CB*3
C
      IFLD   = 1
      NFACES = 2*NDIM
      NXYZ   = NX1*NY1*NZ1
      NEL    = NELFLD(IFIELD)
      NTOT   = NXYZ*NEL
      NFLDT  = NFIELD - 1
C
      CALL RZERO(TMP,NTOT)
C
C     Temperature boundary condition
C
      DO 2100 ISWEEP=1,2
C
c        IF (IFMODEL .AND. IFKEPS .AND. IFIELD.GE.NFLDT)
c    $       CALL TURBWBC (TMP,TMA,SMU)
C
         DO 2010 IE=1,NEL
         DO 2010 IFACE=1,NFACES
            CB=CBC(IFACE,IE,IFIELD)
            BC1=BC(1,IFACE,IE,IFIELD)
            BC2=BC(2,IFACE,IE,IFIELD)
            BC3=BC(3,IFACE,IE,IFIELD)
            BC4=BC(4,IFACE,IE,IFIELD)
            BCK=BC(4,IFACE,IE,IFLD)
            BCE=BC(5,IFACE,IE,IFLD)
            IF (CB.EQ.'T  ') CALL FACEV (TMP,IE,IFACE,BC1,NX1,NY1,NZ1)
            IF (CB.EQ.'MCI') CALL FACEV (TMP,IE,IFACE,BC4,NX1,NY1,NZ1)
            IF (CB.EQ.'MLI') CALL FACEV (TMP,IE,IFACE,BC4,NX1,NY1,NZ1)
            IF (CB.EQ.'KD ') CALL FACEV (TMP,IE,IFACE,BCK,NX1,NY1,NZ1)
            IF (CB.EQ.'ED ') CALL FACEV (TMP,IE,IFACE,BCE,NX1,NY1,NZ1)
            IF (CB.EQ.'t  ' .OR. CB.EQ.'kd ' .OR. CB.EQ.'ed ') 
     $         CALL FACEIS (CB,TMP(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
 2010    CONTINUE
C
C        Take care of Neumann-Dirichlet shared edges...
C
         IF (ISWEEP.EQ.1) CALL DSOP(TMP,'MXA',NX1,NY1,NZ1)
         IF (ISWEEP.EQ.2) CALL DSOP(TMP,'MNA',NX1,NY1,NZ1)
 2100 CONTINUE
C
C     Copy temporary array to temperature array.
C
      CALL COL2(S,TMASK(1,1,1,1,IFIELD-1),NTOT)
      CALL ADD2(S,TMP,NTOT)
      RETURN
      END

C
      SUBROUTINE BCNEUSC(S,ITYPE)
C
C     Apply Neumann boundary conditions to surface of scalar, S.
C     Use IFIELD as a guide to which boundary conditions are to be applied.
C
C     If ITYPE = 1, then S is returned as the rhs contribution to the 
C                   volumetric flux.
C
C     If ITYPE =-1, then S is returned as the lhs contribution to the 
C                   diagonal of A.
C
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'NEKUSE'
C
      DIMENSION S(LX1,LY1,LZ1,LELT)
      common  /nekcb/ cb
      CHARACTER CB*3
C
      NFACES=2*NDIM
      NXYZ  =NX1*NY1*NZ1
      NEL   =NELFLD(IFIELD)
      NTOT  =NXYZ*NEL
      CALL RZERO(S,NTOT)
C
      IF (ITYPE.EQ.-1) THEN
C
C        Compute diagonal contributions to accomodate Robin boundary conditions
C
         DO 1000 IE=1,NEL
         DO 1000 IFACE=1,NFACES
            IEG=LGLEL(IE,NODE)
            CB =CBC(IFACE,IE,IFIELD)
            IF (CB.EQ.'C  ' .OR. CB.EQ.'c  ' .OR.
     $          CB.EQ.'R  ' .OR. CB.EQ.'r  ') THEN
C
               IF (CB.EQ.'C  ') HC   = BC(2,IFACE,IE,IFIELD)
               IF (CB.EQ.'R  ') THEN
                                TINF = BC(1,IFACE,IE,IFIELD)
                                HRAD = BC(2,IFACE,IE,IFIELD)
               ENDIF
               IA=0
C
C IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
C
               CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
               DO 100 IZ=KZ1,KZ2
               DO 100 IY=KY1,KY2
               DO 100 IX=KX1,KX2
                  IA = IA + 1
                  TS = T(IX,IY,IZ,IE,IFIELD-1)
                  IF (CB.EQ.'c  ' .OR. CB.EQ.'r  ') THEN
                     CALL NEKASGN (IX,IY,IZ,IE)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                  ENDIF
                  IF (CB.EQ.'r  ' .OR. CB.EQ.'R  ') 
     $               HC = HRAD * (TINF**2 + TS**2) * (TINF + TS)
                  S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE) +
     $               HC*AREA(IA,1,IFACE,IE)/BM1(IX,IY,IZ,IE)
  100          CONTINUE
            ENDIF
 1000    CONTINUE
      ENDIF
      IF (ITYPE.EQ.1) THEN
C
C        Add passive scalar fluxes to rhs
C
         DO 2000 IE=1,NEL
         DO 2000 IFACE=1,NFACES
            IEG=LGLEL(IE,NODE)
            CB =CBC(IFACE,IE,IFIELD)
            IF (CB.EQ.'F  ' .OR. CB.EQ.'f  ' .OR. 
     $          CB.EQ.'C  ' .OR. CB.EQ.'c  ' .OR. 
     $          CB.EQ.'R  ' .OR. CB.EQ.'r  ' ) THEN
C
                IF (CB.EQ.'F  ') FLUX=BC(1,IFACE,IE,IFIELD)
                IF (CB.EQ.'C  ') FLUX=BC(1,IFACE,IE,IFIELD)
     $                               *BC(2,IFACE,IE,IFIELD)
                IF (CB.EQ.'R  ') THEN
                                 TINF=BC(1,IFACE,IE,IFIELD)
                                 HRAD=BC(2,IFACE,IE,IFIELD)
                ENDIF
C
C              Add local weighted flux values to rhs, S.
C
C IA is areal counter, assumes advancing fastest index first. (IX...IY...IZ)
               IA=0
               CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFACE)
               DO 200 IZ=KZ1,KZ2
               DO 200 IY=KY1,KY2
               DO 200 IX=KX1,KX2
                  IA = IA + 1
                  TS = T(IX,IY,IZ,IE,IFIELD-1)
                  IF (CB.EQ.'f  ') THEN
                     CALL NEKASGN (IX,IY,IZ,IE)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                  ENDIF
                  IF (CB.EQ.'c  ') THEN
                     CALL NEKASGN (IX,IY,IZ,IE)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                     FLUX = TINF*HC
                  ENDIF
                  IF (CB.EQ.'r  ') THEN
                     CALL NEKASGN (IX,IY,IZ,IE)
                     CALL USERBC  (IX,IY,IZ,IFACE,IEG)
                  ENDIF
                  IF (CB.EQ.'R  ' .OR. CB.EQ.'r  ') 
     $               FLUX = HRAD*(TINF**2 + TS**2)*(TINF + TS) * TINF
C
C                 Add computed fluxes to boundary surfaces:
C
                  S(IX,IY,IZ,IE) = S(IX,IY,IZ,IE)
     $                           + FLUX*AREA(IA,1,IFACE,IE)
  200          CONTINUE
            ENDIF
 2000    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE FACEIS (CB,S,IEL,IFACE,NX,NY,NZ)
C
C     Assign inflow boundary conditions to face(IE,IFACE)
C     for scalar S.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'NEKUSE'
C
      DIMENSION S(LX1,LY1,LZ1)
      CHARACTER CB*3
c
      common  /nekcb/ cb3
      character*3 cb3
      cb3 = cb
c
C
C     Passive scalar term
C
      IEG = LGLEL(IEL,NODE)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
C
      IF (CB.EQ.'t  ') THEN
C
         DO 100 IZ=KZ1,KZ2
         DO 100 IY=KY1,KY2
         DO 100 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            S(IX,IY,IZ) = TEMP
  100    CONTINUE
         RETURN
C
      ELSEIF (CB.EQ.'ms ' .OR. CB.EQ.'msi') THEN
C
         DO 200 IZ=KZ1,KZ2
         DO 200 IY=KY1,KY2
         DO 200 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            S(IX,IY,IZ) = SIGMA
  200    CONTINUE
C
      ELSEIF (CB.EQ.'kd ') THEN
C
         DO 300 IZ=KZ1,KZ2
         DO 300 IY=KY1,KY2
         DO 300 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            S(IX,IY,IZ) = TURBK
  300    CONTINUE
C
      ELSEIF (CB.EQ.'ed ') THEN
C
         DO 400 IZ=KZ1,KZ2
         DO 400 IY=KY1,KY2
         DO 400 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            S(IX,IY,IZ) = TURBE
  400    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE FACEIV (CB,V1,V2,V3,IEL,IFACE,NX,NY,NZ)
C
C     Assign fortran function boundary conditions to 
C     face IFACE of element IEL for vector (V1,V2,V3).
C
      INCLUDE 'SIZE'
      INCLUDE 'NEKUSE'
      INCLUDE 'PARALLEL'
C
      DIMENSION V1(NX,NY,NZ),V2(NX,NY,NZ),V3(NX,NY,NZ)
      CHARACTER CB*3
c
      common  /nekcb/ cb3
      character*3 cb3
      cb3 = cb
c
      IEG = LGLEL(IEL,NODE)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
C
      IF (CB.EQ.'v  ' .OR. CB.EQ.'ws ' .OR. CB.EQ.'mv '.OR. 
     $    CB.EQ.'mvn') THEN
C
         DO 100 IZ=KZ1,KZ2
         DO 100 IY=KY1,KY2
         DO 100 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = UX
            V2(IX,IY,IZ) = UY
            V3(IX,IY,IZ) = UZ
  100    CONTINUE
         RETURN
C
      ELSEIF (CB.EQ.'vl ' .OR. CB.EQ.'wsl') THEN
C
         DO 120 IZ=KZ1,KZ2
         DO 120 IY=KY1,KY2
         DO 120 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = UN
            V2(IX,IY,IZ) = U1
            V3(IX,IY,IZ) = U2
  120    CONTINUE
         RETURN
C
      ELSEIF (CB.EQ.'s  ' .OR. CB.EQ.'sh ') THEN
C
         DO 200 IZ=KZ1,KZ2
         DO 200 IY=KY1,KY2
         DO 200 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = TRX
            V2(IX,IY,IZ) = TRY
            V3(IX,IY,IZ) = TRZ
  200    CONTINUE
         RETURN
C
      ELSEIF (CB.EQ.'sl ' .OR. CB.EQ.'shl') THEN
C
         DO 220 IZ=KZ1,KZ2
         DO 220 IY=KY1,KY2
         DO 220 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = TRN
            V2(IX,IY,IZ) = TR1
            V3(IX,IY,IZ) = TR2
  220    CONTINUE
C
      ELSEIF (CB.EQ.'ms ') THEN
C
         DO 240 IZ=KZ1,KZ2
         DO 240 IY=KY1,KY2
         DO 240 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = -PA
            V2(IX,IY,IZ) = TR1
            V3(IX,IY,IZ) = TR2
  240    CONTINUE
C
      ELSEIF (CB.EQ.'on ' .OR. CB.EQ.'o  ') THEN
C
         DO 270 IZ=KZ1,KZ2
         DO 270 IY=KY1,KY2
         DO 270 IX=KX1,KX2
            CALL NEKASGN (IX,IY,IZ,IEL)
            CALL USERBC  (IX,IY,IZ,IFACE,IEG)
            V1(IX,IY,IZ) = -PA
            V2(IX,IY,IZ) = 0.0
            V3(IX,IY,IZ) = 0.0
  270    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE NEKASGN (IX,IY,IZ,IEL)
C
C     Assign NEKTON variables for definition (by user) of
C     boundary conditions at collocation point (IX,IY,IZ)
C     of element IEL.
C
C       X             X-coordinate
C       Y             Y-coordinate
C       Z             Z-coordinate
C       UX            X-velocity
C       UY            Y-velocity
C       UZ            Z-velocity
C       TEMP          Temperature
C       PS1           Passive scalar No. 1
C       PS2           Passive scalar No. 2
C        .             .
C        .             .
C       PS9           Passive scalar No. 9
C       SI2           Strainrate invariant II
C       SI3           Strainrate invariant III
C
C     Variables to be defined by user for imposition of
C     boundary conditions :
C
C       SH1           Shear component No. 1
C       SH2           Shear component No. 2
C       TRX           X-traction
C       TRY           Y-traction
C       TRZ           Z-traction
C       SIGMA         Surface-tension coefficient
C       FLUX          Flux
C       HC            Convection heat transfer coefficient
C       HRAD          Radiation  heat transfer coefficient
C       TINF          Temperature at infinity
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'INPUT'
      INCLUDE 'NEKUSE'
c
      common  /nekcb/ cb
      CHARACTER CB*3
C
      COMMON /SCREV / SII (LX1,LY1,LZ1,LELT)
     $              , SIII(LX1,LY1,LZ1,LELT)
C
        X     = XM1(IX,IY,IZ,IEL)
        Y     = YM1(IX,IY,IZ,IEL)
        Z     = ZM1(IX,IY,IZ,IEL)
        R     = X**2+Y**2
        IF (R.GT.0.0) R=SQRT(R)
        IF (X.NE.0.0 .OR. Y.NE.0.0) THETA = ATAN2(Y,X)
C
        UX    = VX(IX,IY,IZ,IEL)
        UY    = VY(IX,IY,IZ,IEL)
        UZ    = VZ(IX,IY,IZ,IEL)
        TEMP  = T(IX,IY,IZ,IEL,1)
        DO 100 IPS=1,NPSCAL
           PS(IPS) = T(IX,IY,IZ,IEL,IPS+1)
 100    CONTINUE
        SI2   = SII (IX,IY,IZ,IEL)
        SI3   = SIII(IX,IY,IZ,IEL)
c
        cbu   = cb
C
      RETURN
      END
      SUBROUTINE BCNEUTR
C
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      COMMON /SCRSF/ TRX(LX1,LY1,LZ1)
     $             , TRY(LX1,LY1,LZ1)
     $             , TRZ(LX1,LY1,LZ1)
      COMMON /CTMP0/ STC(LX1,LY1,LZ1)
      REAL SIGST(LX1,LY1)
C
      LOGICAL IFALGN,IFNORX,IFNORY,IFNORZ
      common  /nekcb/ cb
      CHARACTER CB*3
C
      IFLD  = 1
      NFACE = 2*NDIM
      NXY1  = NX1*NY1
      NXYZ1 = NX1*NY1*NZ1
C
      DO 100 IEL=1,NELV
      DO 100 IFC=1,NFACE
C
         CB  = CBC (IFC,IEL,IFLD)
         BC1 = BC(1,IFC,IEL,IFLD)
         BC2 = BC(2,IFC,IEL,IFLD)
         BC3 = BC(3,IFC,IEL,IFLD)
         BC4 = BC(4,IFC,IEL,IFLD)
         CALL RZERO3 (TRX,TRY,TRZ,NXYZ1)
C
C        Prescribed tractions and shear tractions
C
         IF (CB.EQ.'S  ' .OR. CB.EQ.'SL ' .OR.
     $       CB.EQ.'SH ' .OR. CB.EQ.'SHL' ) THEN
             CALL TRCON (TRX,TRY,TRZ,BC1,BC2,BC3,IEL,IFC)
             IF (IFQINP(IFC,IEL)) CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
             GOTO 120
         ENDIF
         IF (CB.EQ.'s  ' .OR. CB.EQ.'sl ' .OR.
     $       CB.EQ.'sh ' .OR. CB.EQ.'shl' ) THEN
             CALL FACEIV (CB,TRX,TRY,TRZ,IEL,IFC,NX1,NY1,NZ1)
             CALL FACCVS (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
             IF (IFQINP(IFC,IEL)) CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
             GOTO 120
         ENDIF
C
C        Prescribed outflow ambient pressure
C
         IF (CB.EQ.'ON ' .OR. CB.EQ.'O  ') THEN
             BCN = -BC1
             BC2 =  0.0
             BC3 =  0.0
             CALL TRCON   (TRX,TRY,TRZ,BCN,BC2,BC3,IEL,IFC)
             CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
             GOTO 120
         ENDIF
         IF (CB.EQ.'on ' .OR. CB.EQ.'o  ') THEN
             CALL FACEIV  (CB,TRX,TRY,TRZ,IEL,IFC,NX1,NY1,NZ1)
             CALL FACCVS  (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
             CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
             GOTO 120
         ENDIF
C
C     Surface-tension
C
         IF (CB.EQ.'MS ' .OR. CB.EQ.'MSI' .OR.
     $       CB.EQ.'ms ' .OR. CB.EQ.'msi') THEN
             IF (CB.EQ.'MS ') THEN
                BCN = -BC1
                CALL TRCON   (TRX,TRY,TRZ,BCN,BC2,BC3,IEL,IFC)
                CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
             ENDIF
             IF (CB.EQ.'ms ') THEN
                CALL FACEIV  (CB,TRX,TRY,TRZ,IEL,IFC,NX1,NY1,NZ1)
                CALL FACCVS  (TRX,TRY,TRZ,AREA(1,1,IFC,IEL),IFC)
                CALL GLOBROT (TRX,TRY,TRZ,IEL,IFC)
             ENDIF
             IF (CB(1:1).EQ.'M') THEN
                CALL CFILL  (SIGST,BC4,NXY1)
             ELSE
                CALL FACEIS (CB,STC,IEL,IFC,NX1,NY1,NZ1)
                CALL FACEXS (SIGST,STC,IFC,0)
             ENDIF
             IF (IFAXIS) THEN
                CALL TRSTAX (TRX,TRY,SIGST,IEL,IFC)
             ELSEIF (NDIM.EQ.2) THEN
                CALL TRST2D (TRX,TRY,SIGST,IEL,IFC)
             ELSE
                CALL TRST3D (TRX,TRY,TRZ,SIGST,IEL,IFC)
             ENDIF
         ENDIF
C
c 120    CALL ADD2 (BFX(1,1,1,IEL),TRX,NXYZ1)
c        CALL ADD2 (BFY(1,1,1,IEL),TRY,NXYZ1)
c        IF (NDIM.EQ.3) CALL ADD2 (BFZ(1,1,1,IEL),TRZ,NXYZ1)
C
  120 CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRCON (TRX,TRY,TRZ,TR1,TR2,TR3,IEL,IFC)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TOPOL'
C
      DIMENSION TRX(LX1,LY1,LZ1)
     $        , TRY(LX1,LY1,LZ1)
     $        , TRZ(LX1,LY1,LZ1)
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFC)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
      I = 0
C
      IF (NDIM.EQ.2) THEN
         DO 100 J2=JS2,JF2,JSKIP2
         DO 100 J1=JS1,JF1,JSKIP1
            I = I + 1
            TRX(J1,J2,1) = TR1*AREA(I,1,IFC,IEL)
            TRY(J1,J2,1) = TR2*AREA(I,1,IFC,IEL)
  100    CONTINUE
      ELSE
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            I = I + 1
            TRX(J1,J2,1) = TR1*AREA(I,1,IFC,IEL)
            TRY(J1,J2,1) = TR2*AREA(I,1,IFC,IEL)
            TRZ(J1,J2,1) = TR3*AREA(I,1,IFC,IEL)
  200    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE TRST2D (TRX,TRY,SIGST,IEL,IFC)
C
C     Compute taction due to surface tension (2D)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'DXYZ'
      INCLUDE 'TOPOL'
      INCLUDE 'WZ'
      COMMON /CTMP1/ A1X(LX1),A1Y(LX1),STX(LX1),STY(LX1)
C
      DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),SIGST(LX1,1)
      DIMENSION CANG(2),SANG(2)
      DIMENSION IXN(2),IYN(2),IAN(2)
C
      DO 100 IX=1,NX1
         AA = SIGST(IX,1) * WXM1(IX)
         STX(IX) = T1X(IX,1,IFC,IEL) * AA
         STY(IX) = T1Y(IX,1,IFC,IEL) * AA
  100 CONTINUE
C 
      IF (IFC.EQ.3 .OR. IFC.EQ.4) THEN
         CALL CHSIGN (STX,NX1)
         CALL CHSIGN (STY,NX1)
      ENDIF
C
      IF (IFC.EQ.1 .OR. IFC.EQ.3) THEN
         CALL MXM (DXTM1,NX1,STX,NX1,A1X,1)
         CALL MXM (DXTM1,NX1,STY,NX1,A1Y,1)
      ELSE
         CALL MXM (DYTM1,NY1,STX,NY1,A1X,1)
         CALL MXM (DYTM1,NY1,STY,NY1,A1Y,1)
      ENDIF
C
      CALL DSSET (NX1,NY1,NZ1)
      IFACE  = EFACE1(IFC)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
      I = 0
C
      DO 200 J2=JS2,JF2,JSKIP2
      DO 200 J1=JS1,JF1,JSKIP1
         I = I + 1
         TRX(J1,J2,1) = TRX(J1,J2,1) - A1X(I)
         TRY(J1,J2,1) = TRY(J1,J2,1) - A1Y(I)
  200 CONTINUE
C
C     Contact angle corrections
C
      CALL CTANG2D (CANG,SANG,IXN,IYN,IAN,IFC,IEL)
      DO 500 I=1,2
         IX = IXN(I)
         IY = IYN(I)
         IA = IAN(I)
         TRX(IX,IY,1)=TRX(IX,IY,1) + SIGST(IA,1)*CANG(I)
         TRY(IX,IY,1)=TRY(IX,IY,1) + SIGST(IA,1)*SANG(I)
  500 CONTINUE
C
      RETURN
      END
      SUBROUTINE TRSTAX (TRX,TRY,SIGST,IEL,IFC)
C
C     Compute taction due to surface tension (axisymmetric)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'DXYZ'
      INCLUDE 'TOPOL'
      INCLUDE 'WZ'
      COMMON /CTMP1/ A1X(LX1),A1Y(LX1),A2X(LX1),A2Y(LX1)
     $             , STX(LX1),STY(LX1),XJM1(LX1)
      COMMON /CTMP0/ XFM1(LX1),YFM1(LX1),T1XF(LX1),T1YF(LX1)
C
      DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),SIGST(LX1,LY1)
      DIMENSION CANG(2),SANG(2)
      DIMENSION IXN(2),IYN(2),IAN(2)
      LOGICAL IFGLJ
C
      IFGLJ = .FALSE.
      IF ( IFRZER(IEL) .AND. (IFC.EQ.2 .OR. IFC.EQ.4) ) IFGLJ = .TRUE.
      CALL FACEC2 (XFM1,YFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL),IFC)
C
      IF (IFGLJ) THEN
         CALL MXM (DAM1,NY1,XFM1,NY1,T1XF,1)
         CALL MXM (DAM1,NY1,YFM1,NY1,T1YF,1)
         YS0 = T1YF(1)
      ELSE
         CALL MXM (DXM1,NX1,XFM1,NX1,T1XF,1)
         CALL MXM (DXM1,NX1,YFM1,NX1,T1YF,1)
      ENDIF
C
      DO 10 IX=1,NX1
         XJM1(IX)=SQRT( T1XF(IX)**2 + T1YF(IX)**2 )
         T1XF(IX)=T1XF(IX) / XJM1(IX)
         T1YF(IX)=T1YF(IX) / XJM1(IX)
   10 CONTINUE
C
      IF ( IFGLJ ) THEN
         CALL MXM (DAM1,1,T1XF,NY1,T1XS0,1)
         CALL MXM (DAM1,1,UNY(1,1,IFC,IEL),NY1,UNYS0,1)
         DDX    = WAM1(1)*SIGST(1,1)*T1XS0*YS0
         DDY    = WAM1(1)*SIGST(1,1)*T1YF(1)*YS0*2.0
         A2X(1) = WAM1(1)*SIGST(1,1)*XJM1(1)*UNX(1,1,IFC,IEL)*UNYS0
         A2Y(1) = 0.0
         STX(1) = 0.0
         STY(1) = 0.0
         DO 100 IY=2,NY1
            AA = WAM1(IY) * SIGST(IY,1) / (1.0 + ZAM1(IY))
            STX(IY) = T1XF(IY) * AA
            STY(IY) = T1YF(IY) * AA
            AA = AA * XJM1(IY) * UNY(IY,1,IFC,IEL)
            A2X(IY) = UNX(IY,1,IFC,IEL) * AA
            A2Y(IY) = UNY(IY,1,IFC,IEL) * AA
  100    CONTINUE
      ELSE
         DO 200 IX=1,NX1
            AA = SIGST(IX,1) * WXM1(IX)
            STX(IX) = T1XF(IX) * AA
            STY(IX) = T1YF(IX) * AA
            AA = AA * XJM1(IX) * UNY(IX,1,IFC,IEL)
            A2X(IX) = UNX(IX,1,IFC,IEL) * AA
            A2Y(IX) = UNY(IX,1,IFC,IEL) * AA
  200    CONTINUE
      ENDIF
C
      IF (IFGLJ) THEN
         DO 220 IY=1,NY1
            YSIY = T1YF(IY)*XJM1(IY)
            DTX1 = 0.0
            DTY1 = DATM1(IY,1)*DDY
            DTX2 = YSIY*STX(IY)
            DTY2 = YSIY*STY(IY)
            DTY3 = 0.0
            DO 240 J=2,NY1
               DTYS = DATM1(IY,J)*YFM1(J)
               DTX1 = DTX1 + DTYS*STX(J)
               DTY3 = DTY3 + DTYS*STY(J)
  240       CONTINUE
            A1X(IY) = DTX1 + DTX2
            A1Y(IY) = DTY1 + DTY2 + DTY3
  220    CONTINUE
            A1X(1)  = A1X(1) + DDX
      ELSE
         CALL MXM  (DXTM1,NX1,STX,NX1,A1X,1)
         CALL MXM  (DXTM1,NX1,STY,NX1,A1Y,1)
         CALL COL2 (A1X,YFM1,NX1)
         CALL COL2 (A1Y,YFM1,NX1)
      ENDIF
C
      CALL DSSET (NX1,NY1,NZ1)
      IFACE  = EFACE1(IFC)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
      I = 0
C
      DO 300 J2=JS2,JF2,JSKIP2
      DO 300 J1=JS1,JF1,JSKIP1
         I  = I + 1
         TRX(J1,J2,1) = TRX(J1,J2,1) - A2X(I) - A1X(I)
         TRY(J1,J2,1) = TRY(J1,J2,1) - A2Y(I) - A1Y(I)
  300 CONTINUE
C
C     Contact angle corrections
C
      CALL CTANG2D (CANG,SANG,IXN,IYN,IAN,IFC,IEL)
      DO 500 I=1,2
         IX = IXN(I)
         IY = IYN(I)
         IA = IAN(I)
         AA = SIGST(IA,1)*YM1(IX,IY,1,IEL)
         TRX(IX,IY,1)=TRX(IX,IY,1) + AA*CANG(I)
         TRY(IX,IY,1)=TRY(IX,IY,1) + AA*SANG(I)
  500 CONTINUE
C
      RETURN
      END
      SUBROUTINE CTANG2D (CANG,SANG,IXN,IYN,IAN,IFC,IEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'INPUT'
C
      DIMENSION CANG(2),SANG(2)
      DIMENSION IXN(2),IYN(2),IAN(2),ISN(2),NEBPT(4,2)
      CHARACTER CBN*3
C
      DATA NEBPT /4,1,2,3, 2,3,4,1/
      IFLD = 1
      EPS  = 1.e-6
C
      DO 100 I=1,2
         IFCN    = NEBPT(IFC,I)
         CBN     = CBC(IFCN,IEL,IFLD)
         IXN(I)  = 1
         IYN(I)  = 1
         IAN(I)  = 1
         ISN(I)  = 1
         CANG(I) = 0.0
         SANG(I) = 0.0
         IF (CBN.EQ.'E  '    .OR. CBN.EQ.'P  '    .OR. 
     $       CBN(1:1).EQ.'M' .OR. CBN(1:1).EQ.'m') GOTO 100
         NC = IFC
         IF (I.EQ.2) NC=IFCN
         IF (NC  .EQ.2 .OR. NC  .EQ.3) IXN(I) = NX1
         IF (NC  .EQ.3 .OR. NC  .EQ.4) IYN(I) = NY1
         IF (IFC .EQ.2 .OR. IFC .EQ.3) ISN(I) = NX1
         IF (IFCN.EQ.2 .OR. IFCN.EQ.3) IAN(I) = NX1
         IX = IXN(I)
         IY = IYN(I)
         IA = IAN(I)
         IS = ISN(I)
         IF (CBN(1:1).EQ.'V'   .OR. CBN(1:1).EQ.'v'   .OR. 
     $       CBN     .EQ.'S  ' .OR. CBN     .EQ.'s  ' .OR. 
     $       CBN     .EQ.'SL ' .OR. CBN     .EQ.'sl ' .OR. 
     $       CBN(1:1).EQ.'O'   .OR. CBN(1:1).EQ.'o' ) THEN
             UX=VX(IX,IY,1,IEL)
             UY=VY(IX,IY,1,IEL)
             UM=UX**2 + UY**2
             IF (UM.GT.EPS) THEN
                 UNLX=UNX(IS,1,IFCN,IEL)
                 UNLY=UNY(IS,1,IFCN,IEL)
                 UM=SQRT(UM)
                 DOT =UX*UNLX + UY*UNLY
                 IF (DOT.LT.0.0) UM=-UM
                 CANG(I)=UX/UM
                 SANG(I)=UY/UM
                 GOTO 100
             ENDIF
         ENDIF
         CANG(I)=UNX(IS,1,IFCN,IEL)
         SANG(I)=UNY(IS,1,IFCN,IEL)
  100 CONTINUE
C
      RETURN      
      END
      SUBROUTINE TRST3D (TRX,TRY,TRZ,SIGST,IEL,IFC)
C
C     Compute taction due to surface tension (3D)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'WZ'
      COMMON /CTMP0/  XFM1(LX1,LY1),YFM1(LX1,LY1),ZFM1(LX1,LY1)
      COMMON /CTMP1/  DRM1(LX1,LX1),DRTM1(LX1,LY1)
     $             ,  DSM1(LX1,LX1),DSTM1(LX1,LY1)
     $             ,  WGS(LX1,LY1)
      COMMON /SCRMG/  XRM1(LX1,LY1),YRM1(LX1,LY1),ZRM1(LX1,LY1)
     $             ,  XSM1(LX1,LY1),YSM1(LX1,LY1),ZSM1(LX1,LY1)
      COMMON /SCRUZ/  S1X(LX1,LY1),S1Y(LX1,LY1),S1Z(LX1,LY1)
     $             ,  S2X(LX1,LY1),S2Y(LX1,LY1),S2Z(LX1,LY1)
      COMMON /SCRNS/  G1X(LX1,LY1),G1Y(LX1,LY1),G1Z(LX1,LY1)
     $             ,  G2X(LX1,LY1),G2Y(LX1,LY1),G2Z(LX1,LY1)
     $             ,  GBS(LX1,LY1),GB1L(LX1,LY1),GB2L(LX1,LY1)
C
      DIMENSION TRX(LX1,LY1,LZ1),TRY(LX1,LY1,LZ1),TRZ(LX1,LY1,LZ1)
      DIMENSION SIGST(LX1,LY1)
C
      NXY1 = NX1*NY1
C
      CALL RZERO3 (S1X,S1Y,S1Z,NXY1)
      CALL RZERO3 (S2X,S2Y,S2Z,NXY1)
      CALL FACEXV (XFM1,YFM1,ZFM1,XM1(1,1,1,IEL),YM1(1,1,1,IEL),
     $             ZM1(1,1,1,IEL),IFC,0)
      CALL SETDRS (DRM1,DRTM1,DSM1,DSTM1,IFC)
C
      CALL MXM (DRM1,NX1, XFM1,NX1,XRM1,NY1)
      CALL MXM (DRM1,NX1, YFM1,NX1,YRM1,NY1)
      CALL MXM (DRM1,NX1, ZFM1,NX1,ZRM1,NY1)
      CALL MXM (XFM1,NX1,DSTM1,NY1,XSM1,NY1)
      CALL MXM (YFM1,NX1,DSTM1,NY1,YSM1,NY1)
      CALL MXM (ZFM1,NX1,DSTM1,NY1,ZSM1,NY1)
C
      DO 100 IX=1,NX1
      DO 100 IY=1,NY1
         GB1X=XRM1(IX,IY)
         GB1Y=YRM1(IX,IY)
         GB1Z=ZRM1(IX,IY)
         GB2X=XSM1(IX,IY)
         GB2Y=YSM1(IX,IY)
         GB2Z=ZSM1(IX,IY)
         GB11=GB1X*GB1X + GB1Y*GB1Y + GB1Z*GB1Z
         GB12=GB1X*GB2X + GB1Y*GB2Y + GB1Z*GB2Z
         GB22=GB2X*GB2X + GB2Y*GB2Y + GB2Z*GB2Z
         GDET=GB11*GB22 - GB12*GB12
         IF (GDET .LT. 1.E-20) GO TO 9001
         GT11= GB22/GDET
         GT12=-GB12/GDET
         GT22= GB11/GDET
         GB1L(IX,IY)=SQRT(GB11)
         GB2L(IX,IY)=SQRT(GB22)
         GBS (IX,IY)=SQRT(GDET)
         WGS (IX,IY)=WXM1(IX)*WYM1(IY)*SIGST(IX,IY)
         BB = GBS(IX,IY) * WGS(IX,IY)
         G1X(IX,IY) = BB * ( GT11*GB1X + GT12*GB2X )
         G1Y(IX,IY) = BB * ( GT11*GB1Y + GT12*GB2Y )
         G1Z(IX,IY) = BB * ( GT11*GB1Z + GT12*GB2Z )
         G2X(IX,IY) = BB * ( GT12*GB1X + GT22*GB2X )
         G2Y(IX,IY) = BB * ( GT12*GB1Y + GT22*GB2Y )
         G2Z(IX,IY) = BB * ( GT12*GB1Z + GT22*GB2Z )
  100    CONTINUE
C
      CALL MXM (DRTM1,NX1,G1X,NX1,S1X,NY1)
      CALL MXM (DRTM1,NX1,G1Y,NX1,S1Y,NY1)
      CALL MXM (DRTM1,NX1,G1Z,NX1,S1Z,NY1)
C
      CALL MXM (G2X,NX1,DSM1,NY1,S2X,NY1)
      CALL MXM (G2Y,NX1,DSM1,NY1,S2Y,NY1)
      CALL MXM (G2Z,NX1,DSM1,NY1,S2Z,NY1)
C
      CALL ADD2 (S1X,S2X,NXY1)
      CALL ADD2 (S1Y,S2Y,NXY1)
      CALL ADD2 (S1Z,S2Z,NXY1)
C
C     Contact angle option on hold
C
C      ICONTAC=INT(BC2)
C      IF (ICONTAC.NE.0) THEN
C         IX=1
C         IY=1
C         IF (ICONTAC.GE.3) IY=NY1
C         IF (ICONTAC.EQ.2 .OR. ICONTAC.EQ.3) IX=NX1
C         ANG = BC3 * PI / 180.00
C         RR  = YM1(IX,IY,IZ,IEL)
C         TRX(IX,IY,IZ)=TRX(IX,IY,IZ) + RR*SIGST*COS( ANG )
C         TRY(IX,IY,IZ)=TRY(IX,IY,IZ) + RR*SIGST*SIN( ANG )
C      ENDIF
C
      CALL FACSUB2 (TRX,TRY,TRZ,S1X,S1Y,S1Z,IFC)
C
      RETURN
C
 9001 WRITE ( 6,*) 'Zero area for Element=',IEL,'    Face=',IFC
      STOP
C
      END
      SUBROUTINE SETDRS (DRM1,DRTM1,DSM1,DSTM1,IFC)
C
      INCLUDE 'SIZE'
      INCLUDE 'DXYZ'
C
      DIMENSION DRM1(LX1,LX1),DRTM1(LX1,LX1)
     $        , DSM1(LY1,LY1),DSTM1(LY1,LY1)
C
      NXY1=NX1*NY1
C
      IF (IFC.EQ.5 .OR. IFC.EQ.6) THEN
         CALL COPY (DRM1 ,DXM1 ,NXY1)
         CALL COPY (DSM1 ,DYM1 ,NXY1)
         CALL COPY (DRTM1,DXTM1,NXY1)
         CALL COPY (DSTM1,DYTM1,NXY1)
      ELSEIF (IFC.EQ.2 .OR. IFC.EQ.4) THEN
         CALL COPY (DRM1 ,DYM1 ,NXY1)
         CALL COPY (DSM1 ,DZM1 ,NXY1)
         CALL COPY (DRTM1,DYTM1,NXY1)
         CALL COPY (DSTM1,DZTM1 ,NXY1)
      ELSE     
         CALL COPY (DRM1 ,DZM1 ,NXY1)
         CALL COPY (DSM1 ,DXM1 ,NXY1)
         CALL COPY (DRTM1,DZTM1,NXY1)
         CALL COPY (DSTM1,DXTM1,NXY1)
      ENDIF
C
      RETURN
      END

      SUBROUTINE GLOBROT (R1,R2,R3,IEL,IFC)
C
C     Rotate vector components R1,R2,R3 at face IFC 
C     of element IEL from local to global system.
C
C     R1, R2, R3 have the (NX,NY,NZ) data structure
C     IFACE1 is in the preprocessor notation 
C     IFACE  is the dssum notation.
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TOPOL'
C
      DIMENSION R1(LX1,LY1,LZ1)
     $        , R2(LX1,LY1,LZ1)
     $        , R3(LX1,LY1,LZ1)
C
      CALL DSSET (NX1,NY1,NZ1)
      IFACE  = EFACE1(IFC)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
      I = 0
C
      IF (NDIM.EQ.2) THEN
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            I = I+1
            RNORL = R1(J1,J2,1)
            RTAN1 = R2(J1,J2,1)
            R1(J1,J2,1) = RNORL*UNX(I,1,IFC,IEL) +
     $                    RTAN1*T1X(I,1,IFC,IEL)
            R2(J1,J2,1) = RNORL*UNY(I,1,IFC,IEL) +
     $                    RTAN1*T1Y(I,1,IFC,IEL)
  200    CONTINUE
      ELSE
         DO 300 J2=JS2,JF2,JSKIP2
         DO 300 J1=JS1,JF1,JSKIP1
            I = I+1
            RNORL = R1(J1,J2,1)    
            RTAN1 = R2(J1,J2,1)    
            RTAN2 = R3(J1,J2,1)    
            R1(J1,J2,1) = RNORL*UNX(I,1,IFC,IEL) +
     $                    RTAN1*T1X(I,1,IFC,IEL) +
     $                    RTAN2*T2X(I,1,IFC,IEL)
            R2(J1,J2,1) = RNORL*UNY(I,1,IFC,IEL) +
     $                    RTAN1*T1Y(I,1,IFC,IEL) +
     $                    RTAN2*T2Y(I,1,IFC,IEL)
            R3(J1,J2,1) = RNORL*UNZ(I,1,IFC,IEL) +
     $                    RTAN1*T1Z(I,1,IFC,IEL) +
     $                    RTAN2*T2Z(I,1,IFC,IEL)
  300       CONTINUE
         ENDIF
C
      RETURN
      END

      SUBROUTINE FACEC2 (A1,A2,B1,B2,IFC)
C
C     2-D Geometry only
C     Extract A1,A2 from B1,B2 on surface IFC.
C
C     A1, A2 have the (NX1,  1,NFACE) data structure
C     B1, B2 have the (NX1,NY1,    1) data structure
C
      INCLUDE 'SIZE'
C
      DIMENSION A1(LX1),A2(LX1),B1(LX1,LY1),B2(LX1,LY1)
C
      IX=1
      IY=1
      IF (IFC.EQ.1 .OR. IFC.EQ.3) THEN
         IF (IFC.EQ.3) IY = NY1
         DO 10 IX=1,NX1
            A1(IX)=B1(IX,IY)
            A2(IX)=B2(IX,IY)
   10    CONTINUE
      ELSE
         IF (IFC.EQ.2) IX = NX1
         DO 20 IY=1,NY1
            A1(IY)=B1(IX,IY)
            A2(IY)=B2(IX,IY)
   20   CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE LFALSE (IFA,N)
      LOGICAL IFA(1)
      DO 100 I=1,N
      IFA(I)=.FALSE.
  100 CONTINUE
      RETURN
      END
      SUBROUTINE SETSHL
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      COMMON /SCRMG/ V1(LX1,LY1,LZ1,LELV)
     $             , V2(LX1,LY1,LZ1,LELV)
     $             , V3(LX1,LY1,LZ1,LELV)
     $             , VV(LX1,LY1,LZ1,LELV)
C
      common  /nekcb/ cb
      CHARACTER CB*3
C
      IFIELD = 1
      NFACE  = 2*NDIM
      NTOT1  = NX1*NY1*NZ1*NELV
      DELTA  = 1.E-9
      X      = 1.+DELTA
      Y      = 1.
      DIFF   = ABS(X-Y)
      IF (DIFF.EQ.0.) EPSA = 1.E-06
      IF (DIFF.GT.0.) EPSA = 1.E-13
C
      CALL RZERO3  (V1,V2,V3,NTOT1)
      CALL BCTWALL (V1,V2,V3)
      CALL OPDOT   (VV,V1,V2,V3,V1,V2,V3,NTOT1)
      VDOT  = GLMAX(VV,NTOT1)
      VMAX  = SQRT(VDOT)
      IF (VMAX .LT. EPSA) VMAX = -EPSA
C
      DO 100 IEL=1,NELV
      DO 100 IFC=1,NFACE
         CB=CBC(IFC,IEL,IFIELD)
         IF (CB.NE.'V  ' .AND. CB.NE.'v  '  .AND. CB.NE.'VL ' .AND. 
     $       CB.NE.'vl ') GOTO 100
             IF (VMAX .GT. 0.0) THEN
                 CALL CHKZVN (VMAX,IEL,IFC,IVNORL)
                 IF (IVNORL.EQ.1) GOTO 100
             ENDIF
             IF (CB.EQ.'V  ') CBC(IFC,IEL,IFIELD)='WS '
             IF (CB.EQ.'VL ') CBC(IFC,IEL,IFIELD)='WSL'
             IF (CB.EQ.'v  ') CBC(IFC,IEL,IFIELD)='ws '
             IF (CB.EQ.'vl ') CBC(IFC,IEL,IFIELD)='wsl'
 100  CONTINUE
C
      RETURN
      END
      SUBROUTINE CHKZVN (VMAX,IEL,IFC,IVNORL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      COMMON /SCRMG/ V1(LX1,LY1,LZ1,LELV)
     $             , V2(LX1,LY1,LZ1,LELV)
     $             , V3(LX1,LY1,LZ1,LELV)
     $             , VV(LX1,LY1,LZ1,LELV)
C
      NXZ1  = NX1*NZ1
      TOLV  = 0.01*VMAX
C
      VNOR1 = FACDOT(V1(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
      VNOR2 = FACDOT(V2(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
      VNOR  = VNOR1 + VNOR2
      IF (NDIM.EQ.3) THEN
          VNOR3 = FACDOT(V3(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
          VNOR  = VNOR + VNOR3
      ENDIF
      VNOR = ABS(VNOR) / NXZ1
C
      IVNORL = 1
      IF (VNOR .LT. TOLV) IVNORL = 0
C
      RETURN
      END
      SUBROUTINE BCTWALL (TMP1,TMP2,TMP3)
C
C     Apply Dirichlet boundary conditions to surface of vector (V1,V2,V3)
C     (No antimask operation is applied).
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
C
      DIMENSION TMP1(NX1,NY1,NZ1,1)
     $        , TMP2(NX1,NY1,NZ1,1)
     $        , TMP3(NX1,NY1,NZ1,1)
      common  /nekcb/ cb
      CHARACTER CB*3
C
      NFACE = 2*NDIM
      NTOT1 = NX1*NY1*NZ1*NELV
C
      CALL RZERO (TMP1,NTOT1)
      CALL RZERO (TMP2,NTOT1)
      IF (IF3D) CALL RZERO (TMP3,NTOT1)
C
      DO 2000 IEL=1,NELV
      DO 2000 IFC=1,NFACE
         CB  = CBC (IFC,IEL,IFIELD)
         BC1 = BC(1,IFC,IEL,IFIELD)
         BC2 = BC(2,IFC,IEL,IFIELD)
         BC3 = BC(3,IFC,IEL,IFIELD)
         IF (CB.EQ.'V  ' .OR. CB.EQ.'VL '  .OR.
     $       CB.EQ.'WS ' .OR. CB.EQ.'WSL') THEN
             CALL FACEV (TMP1,IEL,IFC,BC1,NX1,NY1,NZ1)
             CALL FACEV (TMP2,IEL,IFC,BC2,NX1,NY1,NZ1)
             IF (NDIM.EQ.3) CALL FACEV (TMP3,IEL,IFC,BC3,NX1,NY1,NZ1)
             IF (CB.EQ.'VL ' .OR. CB.EQ.'WSL')
     $       CALL GLOBROT (TMP1(1,1,1,IEL),TMP2(1,1,1,IEL),
     $                     TMP3(1,1,1,IEL),IEL,IFC)
         ENDIF
         IF (CB.EQ.'v  ' .OR. CB.EQ.'vl ' .OR. 
     $       CB.EQ.'ws ' .OR. CB.EQ.'wsl' .OR.
     $       CB.EQ.'mv ' .OR. CB.EQ.'mvn') THEN
             CALL FACEIV (CB,TMP1(1,1,1,IEL),TMP2(1,1,1,IEL),
     $                    TMP3(1,1,1,IEL),IEL,IFC,NX1,NY1,NZ1)
             IF (CB.EQ.'vl ' .OR. CB.EQ.'wsl')
     $       CALL GLOBROT (TMP1(1,1,1,IEL),TMP2(1,1,1,IEL),
     $                     TMP3(1,1,1,IEL),IEL,IFC)
         ENDIF
 2000 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ANTIMSK1(X,XMASK,N)
C------------------------------------------------------------------
C
C     Return only Dirichlet boundary values of X
C
C-------------------------------------------------------------------
      REAL  X(1),XMASK(1)
      include 'OPCTR'
C
      DO 100 I=1,N
         X(I) = X(I)*(1.-XMASK(I))
 100  CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine output(u,n)
      integer  n
      real     u(n)

      do i=1,n
         write(6,*) i,u(i)
      enddo
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine output2(u,v,n)
      integer  n
      real     u(n)
      real     v(n)

      do i=1,n
         write(6,*) i,u(i),v(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine einit_poissonq(eix,eiy,eiz)
c
c     solve -nabla^2 phi = q , q from userq()
c     set ex,ey,ez = -grad phi
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'POISSON'

      real    eix(ltot),eiy(ltot),eiz(ltot)

      n =  nx1*ny1*nz1*nelv

      call rone    (mask,n)           ! mask is one everywhere                 
      call setmaskp(mask)             ! mask is zero for PEC, PML 

      call rone    (mult,n)
      call dssum   (mult,nx1,ny1,nz1)

      call setqvol (rhs)              ! rhs
      !do i=1,n
      !    write(6,*) 'm-',i,mult(i),mask(i)           
      !enddo
      !call chsign  (rhs,n)
      call col2    (rhs,bm1,n)        ! B*rhs
      !call exitt

      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,mask,n)

      call rone    (h1,n)             ! Solve  [ h1*A + h2*B ]  phi = B*rhs
      call rzero   (h2,n)

      call invcol1 (mult,n)
      !do i=1,n
      !   write(6,*) 'h1 h2',h1(i),h2(i),mult(i)
      !enddo
      !call exitt

      call copy    (binv,bm1,n)
      call dssum   (binv,nx1,ny1,nz1)
      call invcol1 (binv,n)

      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = 1.e-18 ! 1.e-10
             
c     call cggo 
c    $    (phi,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binv,'init')

      call cggo2 (phi,rhs,mult,mask,r,w,p,z,imsh,isd,n,maxit) 
      if (nid.eq.0) write(6,*) 'POISSON Solver Done!'
     
      !c1 = -1.0
      call cem_grad(eix,eiy,eiz,phi,c1) 

      !call outpost (eix,eiy,eiz,p,phi,'   ')

      return
      end
c------------------------------------------------------------------
c     subroutine axhelm (au,u,helm1,helm2,imesh,isd)
c
c     Compute the (Helmholtz) matrix-vector product,
c     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
c
c------------------------------------------------------------------
c     include 'SIZE'
c     include 'WZ'
c     include 'DXYZ'
c     include 'GEOM'
c     include 'MASS'
c     include 'INPUT'
c     include 'PARALLEL'
c     include 'CTIMER'
c
c     COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
c     COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
c     LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
c
c     REAL           AU    (LX1,LY1,LZ1,1)
c
c     call chsign  (phi,n)
c     call gradm1  (eix,eiy,eiz,phi) 
c     call outpost (eix,eiy,eiz,p,phi,'   ')
c
c     return
c     end
c------------------------------------------------------------------
      subroutine axhelm (au,u,helm1,helm2,imesh,isd)
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'CTIMER'
C
      COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)
      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
      COMMON /CTMP0/ TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)
      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)

      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)
C
      IF(IMESH.EQ.1) NEL=NELV
      IF(IMESH.EQ.2) NEL=NELT

      NXY=NX1*NY1
      NYZ=NY1*NZ1
      NXZ=NX1*NZ1
      NXYZ=NX1*NY1*NZ1
      NTOT=NXYZ*NEL
C
      if (icalld.eq.0) taxhm=0.0
      icalld=icalld+1
      naxhm=icalld
      etime1=dclock()
C
      IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
      CALL RZERO (AU,NTOT)
C
      DO 100 IEL=1,NEL
        ie=iel
C
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IEL) )
C
        IF (NDIM.EQ.2) THEN
C
C       2-d case ...............
C
           IF (IFFAST(IEL)) THEN
C
C          Fast 2-d mode: constant properties and undeformed element
C
           H1 = HELM1(1,1,1,IEL)
           CALL MXM   (WDDX,NX1,U(1,1,1,IEL),NX1,TM1,NYZ)
           CALL MXM   (U(1,1,1,IEL),NX1,WDDYT,NY1,TM2,NY1)
           CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ)
           CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)
           CALL ADD3  (AU(1,1,1,IEL),TM1,TM2,NXYZ)
           CALL CMULT (AU(1,1,1,IEL),H1,NXYZ)
C
           ELSE
C
C          General case, speed-up for undeformed elements
C
           CALL MXM  (DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           CALL MXM  (U(1,1,1,IEL),NX1,DYTM1,NY1,DUDS,NY1)
           CALL COL3 (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3 (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)

           IF (IFDFRM(IEL)) THEN
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
           ENDIF
           CALL COL2 (TMP1,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,HELM1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           CALL MXM  (TMP2,NX1,DYM1,NY1,TM2,NY1)
           CALL ADD2 (AU(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM2,NXYZ)
C
           ENDIF
C
        ELSE
C
C       3-d case ...............
C
           IF (IFFAST(IEL)) THEN
C
C          Fast 3-d mode: constant properties and undeformed element
C
           H1 = HELM1(1,1,1,IEL)
           CALL MXM   (WDDX,NX1,U(1,1,1,IEL),NX1,TM1,NYZ)
           DO 5 IZ=1,NZ1
           CALL MXM   (U(1,1,IZ,IEL),NX1,WDDYT,NY1,TM2(1,1,IZ),NY1)
 5         CONTINUE
           CALL MXM   (U(1,1,1,IEL),NXY,WDDZT,NZ1,TM3,NZ1)
           CALL COL2  (TM1,G4M1(1,1,1,IEL),NXYZ)
           CALL COL2  (TM2,G5M1(1,1,1,IEL),NXYZ)
           CALL COL2  (TM3,G6M1(1,1,1,IEL),NXYZ)
           CALL ADD3  (AU(1,1,1,IEL),TM1,TM2,NXYZ)
           CALL ADD2  (AU(1,1,1,IEL),TM3,NXYZ)
           CALL CMULT (AU(1,1,1,IEL),H1,NXYZ)
C
           ELSE
C
C          General case, speed-up for undeformed elements
C
           CALL MXM(DXM1,NX1,U(1,1,1,IEL),NX1,DUDR,NYZ)
           DO 10 IZ=1,NZ1
              CALL MXM(U(1,1,IZ,IEL),NX1,DYTM1,NY1,DUDS(1,1,IZ),NY1)
   10      CONTINUE
           CALL MXM     (U(1,1,1,IEL),NXY,DZTM1,NZ1,DUDT,NZ1)
           CALL COL3    (TMP1,DUDR,G1M1(1,1,1,IEL),NXYZ)
           CALL COL3    (TMP2,DUDS,G2M1(1,1,1,IEL),NXYZ)
           CALL COL3    (TMP3,DUDT,G3M1(1,1,1,IEL),NXYZ)
           IF (IFDFRM(IEL)) THEN
              CALL ADDCOL3 (TMP1,DUDS,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP1,DUDT,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDR,G4M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP2,DUDT,G6M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDR,G5M1(1,1,1,IEL),NXYZ)
              CALL ADDCOL3 (TMP3,DUDS,G6M1(1,1,1,IEL),NXYZ)
           ENDIF
           CALL COL2 (TMP1,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP2,HELM1(1,1,1,IEL),NXYZ)
           CALL COL2 (TMP3,HELM1(1,1,1,IEL),NXYZ)
           CALL MXM  (DXTM1,NX1,TMP1,NX1,TM1,NYZ)
           DO 20 IZ=1,NZ1
              CALL MXM(TMP2(1,1,IZ),NX1,DYM1,NY1,TM2(1,1,IZ),NY1)
   20      CONTINUE
           CALL MXM  (TMP3,NXY,DZM1,NZ1,TM3,NZ1)
           CALL ADD2 (AU(1,1,1,IEL),TM1,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM2,NXYZ)
           CALL ADD2 (AU(1,1,1,IEL),TM3,NXYZ)
C
           ENDIF
C
        ENDIF
C
 100  CONTINUE
C
      IF (IFH2) CALL ADDCOL4 (AU,HELM2,BM1,U,NTOT)
C
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
      IF (IFAXIS.AND.(ISD.EQ.2)) THEN
         DO 200 IEL=1,NEL
C
            IF (IFRZER(IEL)) THEN
               CALL MXM(U  (1,1,1,IEL),NX1,DATM1,NY1,DUAX,1)
               CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
            ENDIF
C
            DO 190 J=1,NY1
            DO 190 I=1,NX1
               IF (YM1(I,J,1,IEL).NE.0.) THEN
                  TERM1 = BM1(I,J,1,IEL)*U(I,J,1,IEL)/YM1(I,J,1,IEL)**2
                  IF (IFRZER(IEL)) THEN
                     TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J)*DUAX(I)
     $                       *JACM1(I,1,1,IEL)/YSM1(I)
                  ELSE
                     TERM2 = 0.
                  ENDIF
                  AU(I,J,1,IEL) = AU(I,J,1,IEL)
     $                          + HELM1(I,J,1,IEL)*(TERM1+TERM2)
               ENDIF
  190       CONTINUE
  200    CONTINUE
      ENDIF
C
      taxhm=taxhm+(dclock()-etime1)
      return
      end
C
c=======================================================================
      subroutine setfast (helm1,helm2,imesh)
C-------------------------------------------------------------------
C
C     Set logicals for fast evaluation of A*x
C
C-------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      REAL HELM1(NX1,NY1,NZ1,1), HELM2(NX1,NY1,NZ1,1)
C
      IF (IMESH.EQ.1) NEL=NELV
      IF (IMESH.EQ.2) NEL=NELT
      NXYZ = NX1*NY1*NZ1
      NTOT = NXYZ*NEL
C
      DELTA = 1.E-9
      X    = 1.+DELTA
      Y    = 1.
      DIFF = ABS(X-Y)
      IF (DIFF.EQ.0.0) EPSM = 1.E-6
      IF (DIFF.GT.0.0) EPSM = 1.E-13
C
      DO 100 ie=1,NEL
         IFFAST(ie) = .FALSE.
         IF (IFDFRM(ie).OR.IFAXIS .OR. IFMODEL ) THEN
            IFFAST(ie) = .FALSE.
         ELSE
           H1MIN  = VLMIN(HELM1(1,1,1,ie),NXYZ)
           H1MAX  = VLMAX(HELM1(1,1,1,ie),NXYZ)
           TESTH1 = ABS((H1MAX-H1MIN)/(H1MAX+H1MIN))
           IF (TESTH1.LT.EPSM) IFFAST(ie) = .TRUE.
         ENDIF
 100  CONTINUE
c
      IFH2   = .FALSE.
      TESTH2 =  VLAMAX(HELM2,NTOT)
      IF (TESTH2.GT.0.) IFH2 = .TRUE.
      return
      end
C
c=======================================================================
      subroutine sfastax
C----------------------------------------------------------------------
C
C     For undeformed elements, set up appropriate elemental matrices
C     and geometric factors for fast evaluation of Ax.
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      COMMON /FASTAX/ WDDX(LX1,LY1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      LOGICAL IFIRST
      SAVE    IFIRST
      DATA    IFIRST /.TRUE./
C
      NXX=NX1*NX1
      IF (IFIRST) THEN
         CALL RZERO(WDDX,NXX)
         DO 100 I=1,NX1
         DO 100 J=1,NX1
         DO 100 IP=1,NX1
            WDDX(I,J) = WDDX(I,J) + WXM1(IP)*DXM1(IP,I)*DXM1(IP,J)
  100    CONTINUE
         NYY=NY1*NY1
         CALL RZERO(WDDYT,NYY)
         DO 200 I=1,NY1
         DO 200 J=1,NY1
         DO 200 IP=1,NY1
            WDDYT(J,I) = WDDYT(J,I) + WYM1(IP)*DYM1(IP,I)*DYM1(IP,J)
  200    CONTINUE
         NZZ=NZ1*NZ1
         CALL RZERO(WDDZT,NZZ)
         DO 300 I=1,NZ1
         DO 300 J=1,NZ1
         DO 300 IP=1,NZ1
            WDDZT(J,I) = WDDZT(J,I) + WZM1(IP)*DZM1(IP,I)*DZM1(IP,J)
  300    CONTINUE
         IFIRST=.FALSE.
      ENDIF
C
      IF (NDIM.EQ.3) THEN
         DO 1001 IE=1,NELT
            IF (.NOT.IFDFRM(IE)) THEN
               DO 1000 IZ=1,NZ1
               DO 1000 IY=1,NY1
               DO 1000 IX=1,NX1
                  G4M1(IX,IY,IZ,IE)=G1M1(IX,IY,IZ,IE)/WXM1(IX)
                  G5M1(IX,IY,IZ,IE)=G2M1(IX,IY,IZ,IE)/WYM1(IY)
                  G6M1(IX,IY,IZ,IE)=G3M1(IX,IY,IZ,IE)/WZM1(IZ)
 1000          CONTINUE
            ENDIF
 1001    CONTINUE
      ELSE
         DO 2001 IE=1,NELT
            IF (.NOT.IFDFRM(IE)) THEN
               DO 2000 IY=1,NY1
               DO 2000 IX=1,NX1
                  G4M1(IX,IY,1,IE)=G1M1(IX,IY,1,IE)/WXM1(IX)
                  G5M1(IX,IY,1,IE)=G2M1(IX,IY,1,IE)/WYM1(IY)
 2000          CONTINUE
            ENDIF
 2001    CONTINUE
      ENDIF
      return
      end
C
c=======================================================================
      subroutine setprec (dpcm1,helm1,helm2,imesh,isd)
C-------------------------------------------------------------------
C
C     Generate diagonal preconditioner for the Helmholtz operator.
C
C-------------------------------------------------------------------
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
      REAL            DPCM1 (LX1,LY1,LZ1,1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      REAL            HELM1(NX1,NY1,NZ1,1), HELM2(NX1,NY1,NZ1,1)
      REAL YSM1(LY1)
C
      IF(IMESH.EQ.1) NEL=NELV
      IF(IMESH.EQ.2) NEL=NELT

      NTOT = NEL*NX1*NY1*NZ1
C
C     The following 4 lines provide a convenient debugging option
       IF (PARAM( 104 ).NE.0) THEN
          CALL RONE(DPCM1,NTOT)
          return
       ENDIF
C
      CALL RZERO(DPCM1,NTOT)
      DO 1000 IE=1,NEL
C
        IF (IFAXIS) CALL SETAXDY ( IFRZER(IE) )
C
        DO 320 IQ=1,NX1
        DO 320 IZ=1,NZ1
        DO 320 IY=1,NY1
        DO 320 IX=1,NX1
           DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + 
     $                          G1M1(IQ,IY,IZ,IE) * DXTM1(IX,IQ)**2
           !write(6,*) 'g1m1',G1M1(IQ,IY,IZ,IE),DPCM1(IX,IY,IZ,IE)
  320      CONTINUE

        DO 340 IQ=1,NY1
        DO 340 IZ=1,NZ1
        DO 340 IY=1,NY1
        DO 340 IX=1,NX1
           DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + 
     $                          G2M1(IX,IQ,IZ,IE) * DYTM1(IY,IQ)**2
           !write(6,*) 'g2m1',G2M1(IQ,IY,IZ,IE),DPCM1(IX,IY,IZ,IE)
  340      CONTINUE

        IF (NDIM.EQ.3) THEN
           DO 360 IQ=1,NZ1
           DO 360 IZ=1,NZ1
           DO 360 IY=1,NY1
           DO 360 IX=1,NX1
              DPCM1(IX,IY,IZ,IE) = DPCM1(IX,IY,IZ,IE) + 
     $                             G3M1(IX,IY,IQ,IE) * DZTM1(IZ,IQ)**2
  360      CONTINUE
C
C       Add cross terms if element is deformed.
C
        IF (IFDFRM(IE)) THEN
           DO 600 IY=1,NY1
           DO 600 IZ=1,NZ1
           DPCM1(1,IY,IZ,IE) = DPCM1(1,IY,IZ,IE)
     $            + G4M1(1,IY,IZ,IE) * DXTM1(1,1)*DYTM1(IY,IY)
     $            + G5M1(1,IY,IZ,IE) * DXTM1(1,1)*DZTM1(IZ,IZ)
           DPCM1(NX1,IY,IZ,IE) = DPCM1(NX1,IY,IZ,IE)
     $            + G4M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DYTM1(IY,IY)
     $            + G5M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DZTM1(IZ,IZ)
  600      CONTINUE
           DO 700 IX=1,NX1
           DO 700 IZ=1,NZ1
             DPCM1(IX,1,IZ,IE) = DPCM1(IX,1,IZ,IE)
     $            + G4M1(IX,1,IZ,IE) * DYTM1(1,1)*DXTM1(IX,IX)
     $            + G6M1(IX,1,IZ,IE) * DYTM1(1,1)*DZTM1(IZ,IZ)
             DPCM1(IX,NY1,IZ,IE) = DPCM1(IX,NY1,IZ,IE)
     $            + G4M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DXTM1(IX,IX)
     $            + G6M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DZTM1(IZ,IZ)
  700      CONTINUE
           DO 800 IX=1,NX1
           DO 800 IY=1,NY1
             DPCM1(IX,IY,1,IE) = DPCM1(IX,IY,1,IE)
     $            + G5M1(IX,IY,1,IE) * DZTM1(1,1)*DXTM1(IX,IX)
     $            + G6M1(IX,IY,1,IE) * DZTM1(1,1)*DYTM1(IY,IY)
             DPCM1(IX,IY,NZ1,IE) = DPCM1(IX,IY,NZ1,IE)
     $            + G5M1(IX,IY,NZ1,IE) * DZTM1(NZ1,NZ1)*DXTM1(IX,IX)
     $            + G6M1(IX,IY,NZ1,IE) * DZTM1(NZ1,NZ1)*DYTM1(IY,IY)
  800      CONTINUE
        ENDIF
      ELSE
C
       IF (IFDFRM(IE)) THEN
           IZ=1
           DO 602 IY=1,NY1
             DPCM1(1,IY,IZ,IE) = DPCM1(1,IY,IZ,IE)
     $            + G4M1(1,IY,IZ,IE) * DXTM1(1,1)*DYTM1(IY,IY)
             DPCM1(NX1,IY,IZ,IE) = DPCM1(NX1,IY,IZ,IE)
     $            + G4M1(NX1,IY,IZ,IE) * DXTM1(NX1,NX1)*DYTM1(IY,IY)
             !write(6,*) 'grs',G4M1(NX1,IY,IZ,IE),DPCM1(NX1,IY,IZ,IE)
  602      CONTINUE

           DO 702 IX=1,NX1
           DO 702 IZ=1,NZ1
             DPCM1(IX,1,IZ,IE) = DPCM1(IX,1,IZ,IE)
     $            + G4M1(IX,1,IZ,IE) * DYTM1(1,1)*DXTM1(IX,IX)
             DPCM1(IX,NY1,IZ,IE) = DPCM1(IX,NY1,IZ,IE)
     $            + G4M1(IX,NY1,IZ,IE) * DYTM1(NY1,NY1)*DXTM1(IX,IX)
             !write(6,*) 'grs',G4M1(IX,NY1,IZ,IE),DPCM1(IX,NY1,IZ,IE)
  702      CONTINUE
         ENDIF
      ENDIF
 1000 CONTINUE
C
      CALL COL2    (DPCM1,HELM1,NTOT)
      CALL ADDCOL3 (DPCM1,HELM2,BM1,NTOT)
C
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
      IF (IFAXIS.AND.(ISD.EQ.2)) THEN
         DO 1200 IEL=1,NEL
C
            IF (IFRZER(IEL)) THEN
               CALL MXM(YM1(1,1,1,IEL),NX1,DATM1,NY1,YSM1,1)
            ENDIF
C
            DO 1190 J=1,NY1
            DO 1190 I=1,NX1
               IF (YM1(I,J,1,IEL).NE.0.) THEN
                  TERM1 = BM1(I,J,1,IEL)/YM1(I,J,1,IEL)**2
                  IF (IFRZER(IEL)) THEN
                     TERM2 =  WXM1(I)*WAM1(1)*DAM1(1,J)
     $                       *JACM1(I,1,1,IEL)/YSM1(I)
                  ELSE
                     TERM2 = 0.
                  ENDIF
                  DPCM1(I,J,1,IEL) = DPCM1(I,J,1,IEL)
     $                             + HELM1(I,J,1,IEL)*(TERM1+TERM2)
               ENDIF
 1190       CONTINUE
 1200    CONTINUE
      ENDIF
C
      call dssum   (dpcm1,nx1,ny1,nz1)
      call invcol1 (dpcm1,ntot)

      !do i=1,ntot
      !   write(6,*) i,dpcm1(i,1,1,1) !,mask(i)
      !enddo
      !call exitt
C
      return
      end
C
c=======================================================================
      subroutine chktcg1 (tol,res,h1,h2,mask,mult,imesh,isd)
C-------------------------------------------------------------------
C
C     Check that the tolerances are not too small for the CG-solver.
C     Important when calling the CG-solver (Gauss-Lobatto mesh) with
C     zero Neumann b.c.
C
C-------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'EIGEN'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      COMMON /CTMP0/ W1   (LX1,LY1,LZ1,LELT)
     $ ,             W2   (LX1,LY1,LZ1,LELT)
      REAL RES  (LX1,LY1,LZ1,1)
      REAL H1   (LX1,LY1,LZ1,1)
      REAL H2   (LX1,LY1,LZ1,1)
      REAL MULT (LX1,LY1,LZ1,1)
      REAL MASK (LX1,LY1,LZ1,1)
C
      IF (EIGAA.NE.0.) THEN
         ACONDNO = EIGGA/EIGAA
      ELSE
         ACONDNO = 10.
      ENDIF
C
C     Single or double precision???
C
      DELTA = 1.E-9
      X     = 1.+DELTA
      Y     = 1.
      DIFF  = ABS(X-Y)
      IF (DIFF.EQ.0.) EPS = 1.E-6
      IF (DIFF.GT.0.) EPS = 1.E-13
C
      IF (IMESH.EQ.1) THEN
          NL  = NELV
          VOL = VOLVM1
      ELSEIF (IMESH.EQ.2) THEN
          NL  = NELT
          VOL = VOLTM1
      ENDIF

      NTOT1 = NX1*NY1*NZ1*NL
      CALL COPY (W1,RES,NTOT1)
C
      IF (IMESH.EQ.1) THEN
         CALL COL3 (W2,BINVM1,W1,NTOT1)
         RINIT  = SQRT(GLSC3 (W2,W1,MULT,NTOT1)/VOLVM1)
      ELSE
         CALL COL3 (W2,BINTM1,W1,NTOT1)
         RINIT  = SQRT(GLSC3 (W2,W1,MULT,NTOT1)/VOLTM1)
      ENDIF
      RMIN   = EPS*RINIT
      IF (TOL.LT.RMIN) THEN
         IF (NID.EQ.0.AND.IFPRINT)
     $   WRITE (6,*) 'New CG1-tolerance (RINIT*epsm) = ',RMIN,TOL
         TOL = RMIN
      ENDIF
C
      CALL RONE (W1,NTOT1)
      BCNEU1 = GLSC3(W1,MASK,MULT,NTOT1)
      BCNEU2 = GLSC3(W1,W1  ,MULT,NTOT1)
      BCTEST = ABS(BCNEU1-BCNEU2)
C
      CALL AXHELM (W2,W1,H1,H2,IMESH,ISD)
      CALL COL2   (W2,W2,NTOT1)
      CALL COL2   (W2,BM1,NTOT1)
      BCROB  = SQRT(GLSUM(W2,NTOT1)/VOL)
C
      IF ((BCTEST .LT. .1).AND.(BCROB.LT.(EPS*ACONDNO))) THEN
C         OTR = GLSC3 (W1,RES,MULT,NTOT1)
         TOLMIN = RINIT*EPS*10.
         IF (TOL .LT. TOLMIN) THEN
             TOL = TOLMIN
             IF (NID.EQ.0.AND.IFPRINT)
     $       WRITE(6,*) 'New CG1-tolerance (Neumann) = ',TOLMIN
         ENDIF
      ENDIF
C
      return
      end

c=======================================================================
c ggggg
      subroutine cggo2(x,f,mult,mask,r,w,p,z,imsh,isd,n,niter)
      include 'SIZE'

c     Solve Ax=f where A is SPD and is invoked by the routine ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided routine ax(w,z,n) returns  w := Az,  
c
c     User-provided routine solveM(z,r,n) ) returns  z := M^-1 r,  
c
c     User-provided array wght() is used to scale inner-products
c     of the form (p,z) = p'*wght*z
c
      !common /srcsf/ grs(3,lx1,ly1,lelt)
      !common /srccg/ ur(lx1*ly1),us(lx1*ly1)
      !$             , wr(lx1*ly1),ws(lx1*ly1)

      integer n, niter
      real    eps, rhsnorm
      real    x(n),f(n),mult(n),mask(n),r(n),w(n),p(n),z(n)

      logical ifmcor


c     set machine tolerances

      one = 1.
      eps = 1.e-20
      rhsnorm = sqrt(glsc2(f,f,n))
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7
      eps= eps*rhsnorm

      !tube-9: just remain above: no additional setting for eps

      if (nid.eq.0) write(6,5) eps,n,niter
   5  format(' eps=',e15.7,', n=',i6, ', iter=',i6)

      pap = 0.0

      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)

      ifmcor= .false.
      skmin =  glmin(mask,n)

      if (skmin.gt.0) ifmcor = .true.

      if (ifmcor) then
          smean = -1./glsum(mult,n)
          rmean = smean*glsc2(r,mult,n)
          call cadd(r,rmean,n)
      endif

      rnorm= sqrt(glsc3(r,mult,r,n)) ! inner product

      iter = 0
      if (nid.eq.0) write(6,*) 'iter;rnorm',iter,rnorm,ifmcor
      if (nid.eq.0) write(6,*) '.'

      miter = niter

      rtz1=1.0
      do 1000 iter=1,miter

         call solveM(z,r,n)         ! Invert preconditioner here

         if (ifmcor) then
             rmean = smean*glsc2(z,mult,n)
             call cadd(z,rmean,n)
         endif

         rtz2 = rtz1
         rtz1 = glsc3(r,mult,z,n)   ! parallel inner product
         beta = rtz1/rtz2

         if (iter.eq.1) beta=0.0

         call add2s1 (p,z,beta,n)    ! p=p*beta+z
         call axhelm (w,p,h1,h2,imsh,isd)
         call dssum  (w,nx1,ny1,nz1)
         call col2   (w,mask,n)

         pap = glsc3 (w,mult,p,n)

         alpha=  rtz1/pap
         alphm= -alpha

         call add2s2(x,p,alpha,n)
         call add2s2(r,w,alphm,n)

         rtr = glsc2(r,r,n)

         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)

         if (nid.eq.0)  write(6,6) iter,rnorm,alpha,beta,pap
c        if (iter.le.10.or.mod(iter,10).eq.0) write(6,6) iter,rnorm
    6    format(' cg:',i8,1p4e12.4)

         if (rtr.le.rlim2) goto 1001

 1000 continue
 1001 continue

      return
      end
c=======================================================================
      subroutine cggo(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,binv,name)
             ! (phi,rhs,h1,h2,mask,mult,imsh,tol,maxit,isd,binv,'init')
C-------------------------------------------------------------------------
C
C     Solve the Helmholtz equation, H*U = RHS,
C     using preconditioned conjugate gradient iteration.
C     Preconditioner: diag(H).
C
C------------------------------------------------------------------------
      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
c
      COMMON  /CPRINT/ IFPRINT, IFHZPC
      LOGICAL          IFPRINT, IFHZPC
C
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      logical ifmcor
C
      real x(1),f(1),h1(1),h2(1),mask(1),mult(1),binv(1)
      parameter        (lg=lx1*ly1*lz1*lelt)
      COMMON /SCRCG/ d (lg) , scalar(2)
      common /SCRMG/ r (lg) , w (lg) , p (lg) , z (lg)
c
      parameter (maxcg=900)
      common /tdarray/ diagt(maxcg),upper(maxcg)
      common /iterhm/ niterhm
      character*4 name
c
c **  zero out stuff for Lanczos eigenvalue estimator

      call rzero(diagt,maxcg)
      call rzero(upper,maxcg)
C
C     Initialization
C
      NXYZ   = NX1*NY1*NZ1
      NEL    = NELV
      VOL    = VOLVM1

      IF (IMSH.EQ.2) NEL=NELT
      IF (IMSH.EQ.2) VOL=VOLTM1
      n  = NEL*NXYZ
c
      tol=tin
      if (param(22).ne.0) tol=abs(param(22))
      niter = min(maxit,maxcg)
      niter = 2000
C
C     Speed-up for undeformed elements and constant properties.
C     Construct diagonal preconditioner.
      !if (.not.ifsolv) then
      !   call setfast(h1,h2,imesh)
      !   ifsolv = .true.
      !endif
C
C     Set up diag preconditioner.
C
      !call setprec(D,h1,h2,imsh,isd)
      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)
c
c     Check for non-trivial null-space
c
      ifmcor= .false.
      h2max = glmax(h2  ,n)
      skmin = glmin(mask,n)
      if (skmin.gt.0.and.h2max.eq.0) ifmcor = .true.
C
      if (ifmcor) then
         smean = -1./glsum(mult,n)
         rmean = smean*glsc2(r,mult,n)
         call cadd(r,rmean,n)
         !write(6,*) 'cggo 0'
         !call exitt
      endif
C
      krylov = 0
      rtz1=1.0
      niterhm = 0

      do 1000 iter=1,niter

         call col3(z,r,d,n)    ! Jacobi Preconditioner

         if (ifmcor) then
            rmean = smean*glsc2(z,mult,n)
            call cadd(z,rmean,n)
            ! write(6,*) 'cggo 1'
            !call exitt
         endif

         rtz2=rtz1
         scalar(1)=vlsc3 (z,r,mult,n)
         scalar(2)=vlsc32(r,mult,binv,n)
         beta=rtz1/rtz2
         call gop(scalar,w,'+  ',2)
         rtz1=scalar(1)
         rbn2=sqrt(scalar(2)/vol)

         if (iter.eq.1) rbn0 = rbn2
         if (param(22).lt.0) tol=abs(param(22))*rbn0
         
         tol=1.e-14  ! FIXME temporary    02/18/07 misun

         !IF (IFPRINT.AND.NID.EQ.0.AND.PARAM(74).NE.0) THEN
         !   WRITE(6,3002) istep,ITER,name,ifmcor,rbn2,TOL,h1(1),h2(1)
         !ENDIF
 
c        Always take at least one iteration   (for projection) pff 11/23/98
         IF (rbn2.LE.TOL.and.(iter.gt.1 .or. istep.le.5)) THEN
            NITER = ITER-1
            IF(NID.EQ.0.AND.((.NOT.IFHZPC).OR.IFPRINT))
     $      WRITE(6,3000) istep,name,NITER,rbn2,rbn0,TOL
            GO TO 9999
         ENDIF
c
         beta = rtz1/rtz2

         if (iter.eq.1) beta=0.0

         call add2s1 (p,z,beta,n)
         call axhelm (w,p,h1,h2,imsh,isd)
         call dssum  (w,nx1,ny1,nz1)
         call col2   (w,mask,n)

         rho0 = rho
         rho  = glsc3(w,p,mult,n)
         alpha= rtz1/rho
         alphm=-alpha

         call add2s2(x,p ,alpha,n)
         call add2s2(r,w ,alphm,n)

c        Generate tridiagonal matrix for Lanczos scheme
         if (iter.eq.1) then
            krylov = krylov+1
            diagt(iter) = rho/rtz1
         elseif (iter.le.maxcg) then
            krylov = krylov+1
            diagt(iter)    = (beta**2 * rho0 + rho ) / rtz1
            upper(iter-1)  = -beta * rho0 / sqrt(rtz2 * rtz1)
         endif
 1000 continue
      niter = iter-1
c
      if (nid.eq.0) write (6,3001) istep,niter,name,rbn2,rbn0,tol
 3000 format(i9,4x,'Hmholtz ',a4,': ',I6,1p6E13.4)
 3001 format(2i6,' **ERROR**: Failed in HMHOLTZ: ',a4,1p6E13.4)
 3002 format(i3,i6,' Helmholtz ',a4,1x,l4,':',1p6E13.4)
 9999 continue

      niterhm = niter
      ifsolv = .false.
c
c
c     Call eigenvalue routine for Lanczos scheme:
c          two work arrays are req'd if you want to save "diag & upper"
c
c     if (iter.ge.3) then
c        niter = iter-1
c        call calc (diagt,upper,w,z,krylov,dmax,dmin)
c        cond = dmax/dmin
c        if (nid.eq.0) write(6,6) istep,cond,dmin,dmax,' lambda'
c     endif
c   6 format(i9,1p3e12.4,4x,a7)
c
c     if (n.gt.0) write(6,*) 'quit in cggo'
c     if (n.gt.0) call exitt
c     call exitt
      return
      end
c=======================================================================
      function vlsc3(r,b,m,n)
      real r(1),b(1),m(1)
      s = 0.
      do i=1,n
         s = s + r(i)*b(i)*m(i)
      enddo
      vlsc3 = s
      return
      end
c=======================================================================
      function vlsc32(r,b,m,n)
      real r(1),b(1),m(1)
      s = 0.
      do i=1,n
         s = s + b(i)*m(i)*r(i)*r(i)
      enddo
      vlsc32 = s
      return
      end
c=======================================================================
      subroutine setmaskp(mask)
C
C     set mask for poisson problem
C
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      real mask(lx1,ly1,lz1,lelv)

      common    /nekcb/ cb
      character cb*3
      integer   e,f

      ! specially for cbc it works for ifield =2, why ??
      ! misun 10/30/06

      if (IFCEM)  ifield=2
      nfaces = 2*ndim

      do e=1,nelv
      do f=1,nfaces
         cb =cbc(f,e,ifield)
c        Homogeneous Dirichlet boundary conditions
         if (cb.eq.'PEC') call facev (mask,e,f,0.0,nx1,ny1,nz1)
         if (cb.eq.'PML') call facev (mask,e,f,0.0,nx1,ny1,nz1) ! fischer/mmin: 1/23/07
       
      enddo
      enddo
      call dsop(mask,'MUL',nx1,ny1,nz1)

      return
      end
c=======================================================================
      subroutine setqvol(bql)
C-----------------------------------------------------------------------
C
C     Set user specified volumetric forcing function (e.g. heat source).
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
c
      real bql(lx1*ly1*lz1,lelt,1)
C
      NEL   = NELFLD(IFIELD)
      NXYZ1 = NX1*NY1*NZ1
      NTOT1 = NXYZ1*NEL
C
      DO 100 IEL=1,NEL
         CALL NEKUQ (bql,IEL)
  100 CONTINUE
C   
      RETURN
      END
C
C------------------------------------------------------------------
      subroutine nekuq (bql,iel)
C
C     Generate user-specified volumetric source term (temp./p.s.)
C
C------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'PARALLEL'
      include 'TSTEP'
      include 'NEKUSE'
c
      integer iel
      real bql(lx1,ly1,lz1,lelt,1)
c
      IELG = LGLEL(IEL,NODE)
      DO 10 K=1,NZ1
      DO 10 J=1,NY1
      DO 10 I=1,NX1
         CALL NEKASGN (I,J,K,IEL)
         CALL USERQ   (I,J,K,IELG)
         BQL(I,J,K,IEL,1) = QVOL
 10   CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine makeuq
C
C     Fill up user defined forcing function and collocate will the
C     mass matrix on the Gauss-Lobatto mesh.
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      ntot = nx1*ny1*nz1*nelv
c
      call rzero   ( bq ,    ntot)
      call setqvol ( bq          )
      call col2    ( bq ,bm1,ntot)

      return
      end
c-----------------------------------------------------------------------
