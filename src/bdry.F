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
        DO 250 IFIELD=1,NFIELD
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
  250   CONTINUE
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
c             CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFACE,IEL)
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
c              IF ( IFQINP(IFACE,IE) )
c    $         CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE),
c    $                       TMP3(1,1,1,IE),IE,IFACE)
            ENDIF

            IF (CB.EQ.'v  ' .OR. CB.EQ.'vl ' .OR. 
     $          CB.EQ.'ws ' .OR. CB.EQ.'wsl' .OR.
     $          CB.EQ.'mv ' .OR. CB.EQ.'mvn') THEN
                CALL FACEIV (CB,TMP1(1,1,1,IE),TMP2(1,1,1,IE),
     $                       TMP3(1,1,1,IE),IE,IFACE,NX1,NY1,NZ1)
c               IF ( IFQINP(IFACE,IE) )
c    $          CALL GLOBROT (TMP1(1,1,1,IE),TMP2(1,1,1,IE),
c    $                        TMP3(1,1,1,IE),IE,IFACE)
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
            IEG=LGLEL(IE)
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
            IEG=LGLEL(IE)
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
      IEG = LGLEL(IEL)
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
