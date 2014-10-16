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
              
      DO 10 IFLD=1,NFIELD
            IFNONL(IFLD) = .FALSE.
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

      ifield0= 1 
      NFLD   = nfield


      if (nid.eq.0) then
          write(6,*) 'ifflow/ifheat/nfield',ifflow,ifheat,nfield
          write(6,*) 'ifield0/nfld',ifield0,nfld   
      endif

      IF (IFCEM) THEN

        DO 250 IFLD=ifield0,NFLD  
        DO IEL=1,NELV                
          DO IFC=1,NFACE
            CB=CBC(IFC,IEL,IFLD)
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

      DO 400 IFLD=ifield0,NFLD   
         CALL GLLOG(IFNONL(IFLD),.TRUE.)
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
         DO 500 IFLD=ifield0,NFLD
            WRITE (6,*) '  '
            WRITE (6,*) 'IFADVC for field',IFLD,'   = ',IFADVC(IFLD)
            WRITE (6,*) 'IFNONL for field',IFLD,'   = ',IFNONL(IFLD)
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
c        IFLD = 0
c        CALL STSMASK (W1MASK,W2MASK,W3MASK)
c     ENDIF
C
C     Masks for flow variables
C
      IF (IFFLOW) THEN
         IFLD = 1
         NEL    = NELFLD(IFLD)
         NTOT   = NXYZ*NEL
C
C        Pressure mask
C
         CALL RONE(PMASK,NTOT)
         DO 50 IEL=1,NELV
         DO 50 IFACE=1,NFACES
            CB=CBC(IFACE,IEL,IFLD)
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
C        IF (IFSTRS) THEN
c          CALL STSMASK (V1MASK,V2MASK,V3MASK)
C        ELSE
C
           CALL RONE(V1MASK,NTOT)
           CALL RONE(V2MASK,NTOT)
           CALL RONE(V3MASK,NTOT)
C
           DO 100 IEL=1,NELV
           DO 100 IFACE=1,NFACES
              CB =CBC(IFACE,IEL,IFLD)
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
c            IF ( .NOT.IFALGN .OR. IFNORX )
c    $            CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
c            IF ( IFNORY )
c    $            CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
c            IF ( IFNORZ )
c    $            CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             GOTO 100
         ENDIF
         IF (CB.EQ.'ON ') THEN
c            IF ( IFNORY .OR. IFNORZ )
c    $            CALL FACEV (V1MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
c            IF ( .NOT.IFALGN .OR. IFNORX .OR. IFNORZ )
c    $            CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
c            IF ( .NOT.IFALGN .OR. IFNORX .OR. IFNORY )
c    $            CALL FACEV (V3MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
             GOTO 100
         ENDIF
         IF (CB.EQ.'A  ') THEN
c            CALL FACEV (V2MASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
         ENDIF
  100    CONTINUE
C
         CALL DSOP(V1MASK,'MUL',NX1,NY1,NZ1)
         CALL DSOP(V2MASK,'MUL',NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSOP(V3MASK,'MUL',NX1,NY1,NZ1)
C
       ENDIF
C
C     ENDIF
C
C     Masks for passive scalars +
C     k and e if k-e turbulence modem:
C     k = nfield-1
C     e = nfield
C
C
      RETURN
      END


C
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
