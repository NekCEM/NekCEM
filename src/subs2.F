      SUBROUTINE OPDOT (DP,A1,A2,A3,B1,B2,B3,N)
C
      INCLUDE 'SIZE'
C
      DIMENSION DP(LX1,LY1,LZ1,1)
     $        , A1(LX1,LY1,LZ1,1)
     $        , A2(LX1,LY1,LZ1,1)
     $        , A3(LX1,LY1,LZ1,1)
     $        , B1(LX1,LY1,LZ1,1)
     $        , B2(LX1,LY1,LZ1,1)
     $        , B3(LX1,LY1,LZ1,1)
C
      IF (NDIM.EQ.2) THEN
         CALL VDOT2 (DP,A1,A2,B1,B2,N)
      ELSE
         CALL VDOT3 (DP,A1,A2,A3,B1,B2,B3,N)
      ENDIF
C
      RETURN
      END
      SUBROUTINE OPADDS (A1,A2,A3,B1,B2,B3,CONST,N,ISC)
C
      INCLUDE 'SIZE'
C
      DIMENSION A1(LX1,LY1,LZ1,1)
     $        , A2(LX1,LY1,LZ1,1)
     $        , A3(LX1,LY1,LZ1,1)
     $        , B1(LX1,LY1,LZ1,1)
     $        , B2(LX1,LY1,LZ1,1)
     $        , B3(LX1,LY1,LZ1,1)
C
      IF (ISC.EQ.1) THEN
         CALL ADD2S1 (A1,B1,CONST,N)
         CALL ADD2S1 (A2,B2,CONST,N)
         IF (NDIM.EQ.3) CALL ADD2S1 (A3,B3,CONST,N)
      ELSEIF (ISC.EQ.2) THEN
         CALL ADD2S2 (A1,B1,CONST,N)
         CALL ADD2S2 (A2,B2,CONST,N)
         IF (NDIM.EQ.3) CALL ADD2S2 (A3,B3,CONST,N)
      ENDIF
C
      RETURN
      END
      SUBROUTINE FACEXS (A,B,IFACE1,IOP)
C
C     IOP = 0
C     Extract scalar A from B on face IFACE1.
C
C     IOP = 1
C     Extract scalar B from A on face IFACE1.
C
C     A has the (NX,NY,NFACE) data structure
C     B has the (NX,NY,NZ)    data structure
C     IFACE1 is in the preprocessor notation 
C     IFACE  is the dssum notation.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
      DIMENSION A(LX1,LY1),B(LX1,LY1,LZ1)
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      IF (IOP.EQ.0) THEN
         DO 100 J2=JS2,JF2,JSKIP2
         DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A(I,1) = B(J1,J2,1)
  100    CONTINUE
      ELSE
         DO 150 J2=JS2,JF2,JSKIP2
         DO 150 J1=JS1,JF1,JSKIP1
            I = I+1
            B(J1,J2,1) = A(I,1)
  150    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE FACEXV (A1,A2,A3,B1,B2,B3,IFACE1,IOP)
C
C     IOP = 0
C     Extract vector (A1,A2,A3) from (B1,B2,B3) on face IFACE1.
C
C     IOP = 1
C     Extract vector (B1,B2,B3) from (A1,A2,A3) on face IFACE1.
C
C     A1, A2, A3 have the (NX,NY,NFACE) data structure
C     B1, B2, B3 have the (NX,NY,NZ)    data structure
C     IFACE1 is in the preprocessor notation 
C     IFACE  is the dssum notation.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
      DIMENSION A1(LX1,LY1),A2(LX1,LY1),A3(LX1,LY1),
     $          B1(LX1,LY1,LZ1),B2(LX1,LY1,LZ1),B3(LX1,LY1,LZ1)
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
      I = 0
C
      IF (IOP.EQ.0) THEN
         DO 100 J2=JS2,JF2,JSKIP2
         DO 100 J1=JS1,JF1,JSKIP1
            I = I+1
            A1(I,1) = B1(J1,J2,1)
            A2(I,1) = B2(J1,J2,1)
            A3(I,1) = B3(J1,J2,1)
  100    CONTINUE
      ELSE
         DO 150 J2=JS2,JF2,JSKIP2
         DO 150 J1=JS1,JF1,JSKIP1
            I = I+1
            B1(J1,J2,1) = A1(I,1)
            B2(J1,J2,1) = A2(I,1)
            B3(J1,J2,1) = A3(I,1)
  150    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE FACSUB2 (A1,A2,A3,B1,B2,B3,IFACE1)
C
C     Subtract B1,B2,B3 from A1,A2,A3 on surface IFACE1 of element IE.
C
C     A1, A2, A3 have the (NX,NY,NZ)    data structure
C     B1, B2, B3 have the (NX,NY,NFACE) data structure
C     IFACE1 is in the preprocessor notation 
C     IFACE  is the dssum notation.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
      DIMENSION A1(LX1,LY1,LZ1),A2(LX1,LY1,LZ1),A3(LX1,LY1,LZ1),
     $          B1(LX1,LY1),B2(LX1,LY1),B3(LX1,LY1)
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A1(J1,J2,1) = A1(J1,J2,1) - B1(I,1)
         A2(J1,J2,1) = A2(J1,J2,1) - B2(I,1)
         A3(J1,J2,1) = A3(J1,J2,1) - B3(I,1)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CMULT2 (A,B,CONST,N)
      DIMENSION A(1),B(1)
      DO 100 I=1,N
         A(I)=B(I)*CONST
 100  CONTINUE
      RETURN
      END
      SUBROUTINE ADD3S (A,B,C,CONST,N)
      DIMENSION A(1),B(1),C(1)
      DO 100 I=1,N
        A(I)=B(I)+CONST*C(I)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE EMERXIT
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      INCLUDE 'PARALLEL'
C
C     Try to hang in there on the first few time steps (pff 8/92)
      IF (IFTRAN.AND.ISTEP.LT.9) RETURN
C
      LASTEP = 1
      CALL PREPOST(.true.,'   ')
C
      IF (NP.EQ.1) THEN
         WRITE (6,*) '       '
         WRITE (6,*) 
     $   ' Emergency exit at step =',ISTEP,'   time =',TIME
         WRITE (6,*) 
     $   ' Latest solution and data are dumped for post-processing.'
         WRITE (6,*) ' *** STOP ***'
      ELSE
         WRITE (6,*) '       '
         WRITE (6,*) NID,
     $   ' Emergency exit at step =',ISTEP,'   time =',TIME
         WRITE (6,*) 
     $   ' Latest solution and data are dumped for post-processing.'
         WRITE (6,*) ' *** STOP ***'
      ENDIF
C
      call exitt
      RETURN
      END
      SUBROUTINE FACCVS (A1,A2,A3,B,IFACE1)
C
C     Collocate scalar B with vector A, components A1,A2,A3,
C     on the surface IFACE1 of an element.
C
C         A1,A2,A3 have the (NX,NY,NZ) data structure
C         B has the (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A1(LX1,LY1,LZ1),A2(LX1,LY1,LZ1),A3(LX1,LY1,LZ1),
     $          B(LX1,LY1)
C
C     Set up counters
C
      CALL DSSET(NX1,NY1,NZ1)
      IFACE  = EFACE1(IFACE1)
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
            I = I+1
            A1(J1,J2,1) = A1(J1,J2,1)*B(I,1)
            A2(J1,J2,1) = A2(J1,J2,1)*B(I,1)
  100    CONTINUE
      ELSE
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            I = I+1
            A1(J1,J2,1) = A1(J1,J2,1)*B(I,1)
            A2(J1,J2,1) = A2(J1,J2,1)*B(I,1)
            A3(J1,J2,1) = A3(J1,J2,1)*B(I,1)
  200    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE OUTF1 (X,TXT,IEL,IFC)
      INCLUDE 'SIZE'
      DIMENSION X(LX1,LZ1,6,LELT)
      CHARACTER*10 TXT
C
         NFACE = 2*NDIM
         NZI   = NZ1
         NZJ   =  1
         NZINC = -1
         NXI   =  1
         NXJ   = NX1
         NXINC =  1
C
         WRITE(21,106) TXT,IFC,NFACE
         DO 100 J=NZI,NZJ,NZINC
         WRITE(21,105) (X(I,J,IFC,IEL),I=NXI,NXJ,NXINC)
  100    CONTINUE
C
  105 FORMAT(5E15.6)
  106 FORMAT(///,5X,'     ^              ',/,
     $           5X,'   S |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',I2,/,
     $           5X,'       R            ',/)
C
      RETURN
      END
      SUBROUTINE OUTM1 (X,TXT,NP,IEL,IP)
      INCLUDE 'SIZE'
      DIMENSION X(LX1,LY1,LZ1,LELT)
      CHARACTER*10 TXT
C
         NYI   = NY1
         NYJ   =  1
         NYINC = -1
         NXI   =  1
         NXJ   = NX1
         NXINC =  1
C
         WRITE(6,106) TXT,IP,NP
         DO 100 J=NYI,NYJ,NYINC
         WRITE(6,105) (X(I,J,IP,IEL),I=NXI,NXJ,NXINC)
  100    CONTINUE
C
c 105 FORMAT(1p8e10.3)
  105 FORMAT(8f10.3)
  106 FORMAT(///,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',I2,/,
     $           5X,'       X            ',/)
C
      RETURN
      END
      SUBROUTINE OUTM2 (X,TXT,NP,IEL,IP)
      INCLUDE 'SIZE'
      DIMENSION X(LX2,LY2,LZ2,LELV)
      CHARACTER*10 TXT
C
         NYI   = NY2
         NYJ   =  1
         NYINC = -1
         NXI   =  1
         NXJ   = NX2
         NXINC =  1
C
         WRITE(21,106) TXT,IP,NP
         DO 100 J=NYI,NYJ,NYINC
         WRITE(21,105) (X(I,J,IP,IEL),I=NXI,NXJ,NXINC)
  100    CONTINUE
C
  105 FORMAT(5E15.6)
  106 FORMAT(///,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',I2,/,
     $           5X,'       X            ',/)
C
      RETURN
      END
      SUBROUTINE SKIPCNR (NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'INPUT'
      COMMON /INDXFC/ MCRFC(4,6)
     $              , MFCCR(3,8)
     $              , MEGCR(3,8)
     $              , MFCEG(2,12)
     $              , MCREG(2,12)
     $              , MCRRST(3,8)
     $              , MIDRST(3,12)
     $              , MCRIND(8)
     $              , MEDIND(2,4)
     $              , NTEFC(2,12)
     $              , NTCRF(2,3)     
C
      NFACE = 2*NDIM
      NCRFC = NFACE - 2
      NMXCR = 8*NEL
      CALL LFALSE (IFNSKP,NMXCR)
C
      DO 100 IEL=1,NEL
      DO 100 IFC=1,NFACE
         IF (CDOF(IFC,IEL).EQ.'1') THEN
             ICR=MCRFC(1,IFC)
             IFNSKP(ICR,IEL)=.TRUE.
         ELSEIF (CDOF(IFC,IEL).EQ.'2') THEN
             ICR=MCRFC(2,IFC)
             IFNSKP(ICR,IEL)=.TRUE.
         ELSEIF (CDOF(IFC,IEL).EQ.'3') THEN
             ICR=MCRFC(3,IFC)
             IFNSKP(ICR,IEL)=.TRUE.
         ELSEIF (CDOF(IFC,IEL).EQ.'4') THEN
             ICR=MCRFC(4,IFC)
             IFNSKP(ICR,IEL)=.TRUE.
         ENDIF
         IF (CDOF(IFC,IEL).EQ.'*') THEN
         DO 160 ICRFC=1,NCRFC
            ICR=MCRFC(ICRFC,IFC)
            IFNSKP(ICR,IEL)=.TRUE.
  160    CONTINUE
         ENDIF
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE SETMASK (C1MASK,C2MASK,C3MASK,HVMASK,NEL) 
C 
      INCLUDE 'SIZE'
C
      DIMENSION HVMASK (LX1,LY1,LZ1,1)
     $        , C1MASK(LX1,LY1,LZ1,1)
     $        , C2MASK(LX1,LY1,LZ1,1)
     $        , C3MASK(LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
      CALL RZERO3  (C1MASK,C2MASK,C3MASK,NTOT1)
C
      DO 100 IEL=1,NEL
      DO 100 IZ=1,NZ1
      DO 100 IY=1,NY1
      DO 100 IX=1,NX1
         HMV=ABS( HVMASK(IX,IY,IZ,IEL) )
         IF (HMV .GT. 2.9) THEN
             C1MASK(IX,IY,IZ,IEL) = 1.0
         ENDIF
         IF ((HMV.GT.1.9 .AND. HMV.LT.2.1) .OR. HMV.GT.4.9) THEN
             C2MASK(IX,IY,IZ,IEL) = 1.0
         ENDIF
  100 CONTINUE
C
      IF (NDIM.EQ.3) CALL COPY (C3MASK,C2MASK,NTOT1)
C
      RETURN
      END
      SUBROUTINE SETCDAT 
C
      INCLUDE 'SIZE'
      COMMON /INDXFC/ MCRFC(4,6)
     $              , MFCCR(3,8)
     $              , MEGCR(3,8)
     $              , MFCEG(2,12)
     $              , MCREG(2,12)
     $              , MCRRST(3,8)
     $              , MIDRST(3,12)
     $              , MCRIND(8)
     $              , MEDIND(2,4)
     $              , NTEFC(2,12)
     $              , NTCRF(2,3)     
C
      NMID = (NX1 +1)/2
C
C     Corners on faces
C
      MCRFC(1,1) = 1
      MCRFC(2,1) = 2
      MCRFC(3,1) = 6
      MCRFC(4,1) = 5
      MCRFC(1,2) = 2
      MCRFC(2,2) = 3
      MCRFC(3,2) = 7
      MCRFC(4,2) = 6
      MCRFC(1,3) = 3
      MCRFC(2,3) = 4
      MCRFC(3,3) = 8
      MCRFC(4,3) = 7
      MCRFC(1,4) = 4
      MCRFC(2,4) = 1
      MCRFC(3,4) = 5
      MCRFC(4,4) = 8
      MCRFC(1,5) = 4
      MCRFC(2,5) = 3
      MCRFC(3,5) = 2
      MCRFC(4,5) = 1
      MCRFC(1,6) = 5
      MCRFC(2,6) = 6
      MCRFC(3,6) = 7
      MCRFC(4,6) = 8
C
C     Faces at corners
C
      MFCCR(1,1) = 4
      MFCCR(2,1) = 1
      MFCCR(3,1) = 5
      MFCCR(1,2) = 1
      MFCCR(2,2) = 2
      MFCCR(3,2) = 5
      MFCCR(1,3) = 2
      MFCCR(2,3) = 3
      MFCCR(3,3) = 5
      MFCCR(1,4) = 3
      MFCCR(2,4) = 4
      MFCCR(3,4) = 5
      MFCCR(1,5) = 4
      MFCCR(2,5) = 1
      MFCCR(3,5) = 6
      MFCCR(1,6) = 1
      MFCCR(2,6) = 2
      MFCCR(3,6) = 6
      MFCCR(1,7) = 2
      MFCCR(2,7) = 3
      MFCCR(3,7) = 6
      MFCCR(1,8) = 3
      MFCCR(2,8) = 4
      MFCCR(3,8) = 6
C
C     Edges at corners
C
      MEGCR(1,1) = 4
      MEGCR(2,1) = 1
      MEGCR(3,1) = 9
      MEGCR(1,2) = 1
      MEGCR(2,2) = 2
      MEGCR(3,2) = 10
      MEGCR(1,3) = 2
      MEGCR(2,3) = 3
      MEGCR(3,3) = 11
      MEGCR(1,4) = 3
      MEGCR(2,4) = 4
      MEGCR(3,4) = 12
      MEGCR(1,5) = 8
      MEGCR(2,5) = 5
      MEGCR(3,5) = 9
      MEGCR(1,6) = 5
      MEGCR(2,6) = 6
      MEGCR(3,6) = 10
      MEGCR(1,7) = 6
      MEGCR(2,7) = 7
      MEGCR(3,7) = 11
      MEGCR(1,8) = 7
      MEGCR(2,8) = 8
      MEGCR(3,8) = 12
C
C     Faces on edges
C
      MFCEG(1,1)  = 1
      MFCEG(2,1)  = 5
      MFCEG(1,2)  = 2
      MFCEG(2,2)  = 5
      MFCEG(1,3)  = 3
      MFCEG(2,3)  = 5
      MFCEG(1,4)  = 4
      MFCEG(2,4)  = 5
      MFCEG(1,5)  = 1
      MFCEG(2,5)  = 6
      MFCEG(1,6)  = 2
      MFCEG(2,6)  = 6
      MFCEG(1,7)  = 3
      MFCEG(2,7)  = 6
      MFCEG(1,8)  = 4
      MFCEG(2,8)  = 6
      MFCEG(1,9)  = 4
      MFCEG(2,9)  = 1
      MFCEG(1,10) = 1
      MFCEG(2,10) = 2
      MFCEG(1,11) = 2
      MFCEG(2,11) = 3
      MFCEG(1,12) = 3
      MFCEG(2,12) = 4
C
C     Corners at edges
C
      MCREG(1,1)  = 1
      MCREG(2,1)  = 2
      MCREG(1,2)  = 2
      MCREG(2,2)  = 3
      MCREG(1,3)  = 4
      MCREG(2,3)  = 3
      MCREG(1,4)  = 1
      MCREG(2,4)  = 4
      MCREG(1,5)  = 5
      MCREG(2,5)  = 6
      MCREG(1,6)  = 6
      MCREG(2,6)  = 7
      MCREG(1,7)  = 8
      MCREG(2,7)  = 7
      MCREG(1,8)  = 5
      MCREG(2,8)  = 8
      MCREG(1,9)  = 1
      MCREG(2,9)  = 5
      MCREG(1,10) = 2
      MCREG(2,10) = 6
      MCREG(1,11) = 3
      MCREG(2,11) = 7
      MCREG(1,12) = 4
      MCREG(2,12) = 8
C
C     Corner indices (Vol array)
C
      MCRRST(1,1) = 1
      MCRRST(2,1) = 1
      MCRRST(3,1) = 1
      MCRRST(1,2) = NX1
      MCRRST(2,2) = 1
      MCRRST(3,2) = 1
      MCRRST(1,3) = NX1
      MCRRST(2,3) = NX1
      MCRRST(3,3) = 1
      MCRRST(1,4) = 1
      MCRRST(2,4) = NX1
      MCRRST(3,4) = 1
      MCRRST(1,5) = 1
      MCRRST(2,5) = 1
      MCRRST(3,5) = NX1
      MCRRST(1,6) = NX1
      MCRRST(2,6) = 1
      MCRRST(3,6) = NX1
      MCRRST(1,7) = NX1
      MCRRST(2,7) = NX1
      MCRRST(3,7) = NX1
      MCRRST(1,8) = 1
      MCRRST(2,8) = NX1
      MCRRST(3,8) = NX1
C
C     Mid-edge indcies (Vol array)
C 
      MIDRST(1,1)  = NMID
      MIDRST(1,2)  = NX1
      MIDRST(1,3)  = NMID
      MIDRST(1,4)  = 1
      MIDRST(1,5)  = NMID
      MIDRST(1,6)  = NX1
      MIDRST(1,7)  = NMID
      MIDRST(1,8)  = 1
      MIDRST(1,9)  = 1
      MIDRST(1,10) = NX1
      MIDRST(1,11) = NX1
      MIDRST(1,12) = 1
      MIDRST(2,1)  = 1
      MIDRST(2,2)  = NMID
      MIDRST(2,3)  = NX1
      MIDRST(2,4)  = NMID
      MIDRST(2,5)  = 1
      MIDRST(2,6)  = NMID
      MIDRST(2,7)  = NX1
      MIDRST(2,8)  = NMID
      MIDRST(2,9)  = 1
      MIDRST(2,10) = 1
      MIDRST(2,11) = NX1
      MIDRST(2,12) = NX1
      MIDRST(3,1)  = 1
      MIDRST(3,2)  = 1
      MIDRST(3,3)  = 1
      MIDRST(3,4)  = 1
      MIDRST(3,5)  = NX1
      MIDRST(3,6)  = NX1
      MIDRST(3,7)  = NX1
      MIDRST(3,8)  = NX1
      MIDRST(3,9)  = NMID
      MIDRST(3,10) = NMID
      MIDRST(3,11) = NMID
      MIDRST(3,12) = NMID
C
C     1-D corners indices (Vol array)
C
      MCRIND(1) = 1
      MCRIND(2) = NX1
      MCRIND(3) = NX1**2
      MCRIND(7) = NX1**3
      MCRIND(4) = MCRIND(3) - NX1 + 1
      MCRIND(5) = MCRIND(7) - MCRIND(3) + 1
      MCRIND(6) = MCRIND(5) + NX1 - 1
      MCRIND(8) = MCRIND(7) - NX1 + 1
C
C     1-D  edge indices (Face array)
C
      MEDIND(1,1) = 1
      MEDIND(2,1) = NX1
      MEDIND(1,2) = NX1**2 - NX1 + 1
      MEDIND(2,2) = NX1**2
      MEDIND(1,3) = 1
      MEDIND(2,3) = MEDIND(1,2)
      MEDIND(1,4) = NX1
      MEDIND(2,4) = NX1**2
C
C     1-D edge index type (Face array)
C
      NTEFC(1,1)  = 1
      NTEFC(2,1)  = 1
      NTEFC(1,2)  = 1
      NTEFC(2,2)  = 4
      NTEFC(1,3)  = 1
      NTEFC(2,3)  = 2
      NTEFC(1,4)  = 1
      NTEFC(2,4)  = 3
      NTEFC(1,5)  = 2
      NTEFC(2,5)  = 1
      NTEFC(1,6)  = 2
      NTEFC(2,6)  = 4
      NTEFC(1,7)  = 2
      NTEFC(2,7)  = 2
      NTEFC(1,8)  = 2
      NTEFC(2,8)  = 3
      NTEFC(1,9)  = 3
      NTEFC(2,9)  = 3
      NTEFC(1,10) = 4
      NTEFC(2,10) = 3
      NTEFC(1,11) = 4
      NTEFC(2,11) = 4
      NTEFC(1,12) = 3
      NTEFC(2,12) = 4
C
C     Corner index address on face in MCRRST
C
      NTCRF(1,1) = 1
      NTCRF(2,1) = 3
      NTCRF(1,2) = 2
      NTCRF(2,2) = 3
      NTCRF(1,3) = 1
      NTCRF(2,3) = 2
C
      RETURN
      END
      SUBROUTINE EDGINDF (LF1,LF2,LFSKIP,ISD,IFCN)
C
      INCLUDE 'SIZE'
      COMMON /INDXFC/ MCRFC(4,6)
     $              , MFCCR(3,8)
     $              , MEGCR(3,8)
     $              , MFCEG(2,12)
     $              , MCREG(2,12)
     $              , MCRRST(3,8)
     $              , MIDRST(3,12)
     $              , MCRIND(8)
     $              , MEDIND(2,4)
     $              , NTEFC(2,12)
     $              , NTCRF(2,3)     
C
      ITYP = NTEFC(IFCN,ISD)
C
      LF1 = MEDIND(1,ITYP)
      LF2 = MEDIND(2,ITYP)
C
      LFSKIP = 1
      IF (ITYP .GE. 3) LFSKIP = NX1
C
      RETURN
      END
      SUBROUTINE EDGINDV (LV1,LV2,LVSKIP,ISD)
C
      INCLUDE 'SIZE'
      COMMON /INDXFC/ MCRFC(4,6)
     $              , MFCCR(3,8)
     $              , MEGCR(3,8)
     $              , MFCEG(2,12)
     $              , MCREG(2,12)
     $              , MCRRST(3,8)
     $              , MIDRST(3,12)
     $              , MCRIND(8)
     $              , MEDIND(2,4)
     $              , NTEFC(2,12)
     $              , NTCRF(2,3)     
C
      IODD = ISD - ISD/2*2
      ICR1 = MCREG(1,ISD)
      ICR2 = MCREG(2,ISD)
C
      LV1  = MCRIND(ICR1)
      LV2  = MCRIND(ICR2)
C
      IF (ISD .GE. 9) THEN
         LVSKIP = NX1**2
      ELSE
         IF (IODD.EQ.0) THEN
            LVSKIP = NX1
         ELSE
            LVSKIP = 1
         ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE AMASK (VB1,VB2,VB3,V1,V2,V3,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      COMMON /SCRSF/ A1MASK(LX1,LY1,LZ1,LELT)
     $             , A2MASK(LX1,LY1,LZ1,LELT)
     $             , A3MASK(LX1,LY1,LZ1,LELT)
      COMMON /CTMP0/ WA(LX1,LY1,LZ1,LELT)
C
      DIMENSION VB1(LX1,LY1,LZ1,1)
     $        , VB2(LX1,LY1,LZ1,1)
     $        , VB3(LX1,LY1,LZ1,1)
     $        , V1(LX1,LY1,LZ1,1)
     $        , V2(LX1,LY1,LZ1,1)
     $        , V3(LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
      CALL RONE (WA,NTOT1)
      CALL COPY (VB1,V1,NTOT1)
      CALL COPY (VB2,V2,NTOT1)
      IF (NDIM.EQ.3) CALL COPY (VB3,V3,NTOT1)
C
      IF (IFIELD.EQ.1) THEN
         CALL SUB3  (A1MASK,WA,V1MASK,NTOT1)
         CALL SUB3  (A2MASK,WA,V2MASK,NTOT1)
         IF (NDIM.EQ.3) CALL SUB3 (A3MASK,WA,V3MASK,NTOT1)
      ELSE
         CALL SUB3  (A1MASK,WA,W1MASK,NTOT1)
         CALL SUB3  (A2MASK,WA,W2MASK,NTOT1)
         IF (NDIM.EQ.3) CALL SUB3 (A3MASK,WA,W3MASK,NTOT1)
      ENDIF
C
      CALL QMASK (VB1,VB2,VB3,A1MASK,A2MASK,A3MASK,NEL)
C
      RETURN
      END
      SUBROUTINE RMASK (R1,R2,R3,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'MVGEOM'
C
      DIMENSION R1  (LX1,LY1,LZ1,1)
     $        , R2  (LX1,LY1,LZ1,1)
     $        , R3  (LX1,LY1,LZ1,1)
C
      IF (IFIELD.EQ.1) THEN
         CALL QMASK (R1,R2,R3,V1MASK,V2MASK,V3MASK,NEL)
      ELSE
         CALL QMASK (R1,R2,R3,W1MASK,W2MASK,W3MASK,NEL)
      ENDIF
C
      RETURN
      END
      SUBROUTINE QMASK (R1,R2,R3,R1MASK,R2MASK,R3MASK,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      COMMON /CTMP1/ S1(LX1,LY1,LZ1,LELT)
     $             , S2(LX1,LY1,LZ1,LELT)
     $             , S3(LX1,LY1,LZ1,LELT)
C
      DIMENSION R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , R3(LX1,LY1,LZ1,1)
     $        , R1MASK(LX1,LY1,LZ1,1)
     $        , R2MASK(LX1,LY1,LZ1,1)
     $        , R3MASK(LX1,LY1,LZ1,1)
C
      NTOT1 = NX1*NY1*NZ1*NEL
C
C     (0) Collocate Volume Mask
C
      CALL COPY  (S1,R1,NTOT1)
      CALL COPY  (S2,R2,NTOT1)
      CALL COL2  (R1,R1MASK,NTOT1)
      CALL COL2  (R2,R2MASK,NTOT1)
      IF (NDIM.EQ.3) THEN
         CALL COPY (S3,R3,NTOT1)
         CALL COL2 (R3,R3MASK,NTOT1)
      ENDIF
C
C     (1) Face Mask
C
      IF (IFLMSF(IFIELD)) THEN
         IF (NDIM.EQ.2) THEN
            CALL FCMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)
         ELSE
            CALL FCMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
         ENDIF
      ENDIF
C
C     (2) Edge Mask  (3-D only)
C
      IF (NDIM.EQ.3 .AND. IFLMSE(IFIELD)) 
     $   CALL EGMASK (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
C
C     (3) Corner Mask
C
      IF (IFLMSC(IFIELD)) THEN
         IF (NDIM.EQ.2) THEN
            CALL CRMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)
         ELSE
            CALL CRMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
         ENDIF
      ENDIF
C
      RETURN
      END
      SUBROUTINE FCMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      DIMENSION R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , S1(LX1,LY1,LZ1,1)
     $        , S2(LX1,LY1,LZ1,1)
     $        , R1MASK(LX1,LY1,LZ1,1)
     $        , R2MASK(LX1,LY1,LZ1,1)
C
      NFACE = 2*NDIM
C
      DO 100 IEL=1,NEL
      DO 100 IFC=1,NFACE
         IF (.NOT.IFMSFC(IFC,IEL,IFIELD)) GO TO 100
         CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
         DO 120 J2=JS2,JF2,JSKIP2
         DO 120 J1=JS1,JF1,JSKIP1
            RNOR = ( S1(J1,J2,1,IEL)*VNX(J1,J2,1,IEL) +
     $               S2(J1,J2,1,IEL)*VNY(J1,J2,1,IEL) ) *
     $               R1MASK(J1,J2,1,IEL)
            RTN1 = ( S1(J1,J2,1,IEL)*V1X(J1,J2,1,IEL) +
     $               S2(J1,J2,1,IEL)*V1Y(J1,J2,1,IEL) ) *
     $               R2MASK(J1,J2,1,IEL)
            R1(J1,J2,1,IEL) = RNOR*VNX(J1,J2,1,IEL) +
     $                        RTN1*V1X(J1,J2,1,IEL)
            R2(J1,J2,1,IEL) = RNOR*VNY(J1,J2,1,IEL) +
     $                        RTN1*V1Y(J1,J2,1,IEL)
            write(6,*) '=====FCMASK===',R1(J1,J2,1,IEL),R2(J1,J2,1,IEL)
  120       CONTINUE
  100    CONTINUE
C
      RETURN
      END
      SUBROUTINE FCMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      DIMENSION R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , R3(LX1,LY1,LZ1,1)
     $        , S1(LX1,LY1,LZ1,1)
     $        , S2(LX1,LY1,LZ1,1)
     $        , S3(LX1,LY1,LZ1,1)
     $        , R1MASK(LX1,LY1,LZ1,1)
     $        , R2MASK(LX1,LY1,LZ1,1)
     $        , R3MASK(LX1,LY1,LZ1,1)
C
      NFACE = 2*NDIM
C
      DO 100 IEL=1,NEL
      DO 100 IFC=1,NFACE
         IF (.NOT.IFMSFC(IFC,IEL,IFIELD)) GO TO 100
         CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
         DO 120 J2=JS2,JF2,JSKIP2
         DO 120 J1=JS1,JF1,JSKIP1
            RNOR = ( S1(J1,J2,1,IEL)*VNX(J1,J2,1,IEL) +
     $               S2(J1,J2,1,IEL)*VNY(J1,J2,1,IEL) +
     $               S3(J1,J2,1,IEL)*VNZ(J1,J2,1,IEL) ) *
     $               R1MASK(J1,J2,1,IEL)
            RTN1 = ( S1(J1,J2,1,IEL)*V1X(J1,J2,1,IEL) +
     $               S2(J1,J2,1,IEL)*V1Y(J1,J2,1,IEL) +
     $               S3(J1,J2,1,IEL)*V1Z(J1,J2,1,IEL) ) *
     $               R2MASK(J1,J2,1,IEL)
            RTN2 = ( S1(J1,J2,1,IEL)*V2X(J1,J2,1,IEL) +
     $               S2(J1,J2,1,IEL)*V2Y(J1,J2,1,IEL) +
     $               S3(J1,J2,1,IEL)*V2Z(J1,J2,1,IEL) ) *
     $               R3MASK(J1,J2,1,IEL)
            R1(J1,J2,1,IEL) = RNOR*VNX(J1,J2,1,IEL) +
     $                        RTN1*V1X(J1,J2,1,IEL) +
     $                        RTN2*V2X(J1,J2,1,IEL)
            R2(J1,J2,1,IEL) = RNOR*VNY(J1,J2,1,IEL) +
     $                        RTN1*V1Y(J1,J2,1,IEL) +
     $                        RTN2*V2Y(J1,J2,1,IEL)
            R3(J1,J2,1,IEL) = RNOR*VNZ(J1,J2,1,IEL) +
     $                        RTN1*V1Z(J1,J2,1,IEL) +
     $                        RTN2*V2Z(J1,J2,1,IEL) 
            write(6,*) '=====FCMASK3===',R1(J1,J2,1,IEL),R2(J1,J2,1,IEL)
  120       CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE EGMASK (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      DIMENSION R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , R3(LX1,LY1,LZ1,1)
     $        , S1(LX1,LY1,LZ1,1)
     $        , S2(LX1,LY1,LZ1,1)
     $        , S3(LX1,LY1,LZ1,1)
     $        , R1MASK(LX1,LY1,LZ1,1)
     $        , R2MASK(LX1,LY1,LZ1,1)
     $        , R3MASK(LX1,LY1,LZ1,1)
C
      NEDGE = 12
C
      DO 100 IEL=1,NEL
      DO 100 ISD=1,NEDGE
         IF (.NOT.IFMSEG(ISD,IEL,IFIELD)) GOTO 100
         CALL EDGINDV (LV1,LV2,LVSKIP,ISD)
         DO 120 LV=LV1,LV2,LVSKIP
            RNOR = ( S1(LV,1,1,IEL)*VNX(LV,1,1,IEL) +
     $               S2(LV,1,1,IEL)*VNY(LV,1,1,IEL) +
     $               S3(LV,1,1,IEL)*VNZ(LV,1,1,IEL) ) *
     $               R1MASK(LV,1,1,IEL)
            RTN1 = ( S1(LV,1,1,IEL)*V1X(LV,1,1,IEL) +
     $               S2(LV,1,1,IEL)*V1Y(LV,1,1,IEL) +
     $               S3(LV,1,1,IEL)*V1Z(LV,1,1,IEL) ) *
     $               R2MASK(LV,1,1,IEL)
            RTN2 = ( S1(LV,1,1,IEL)*V2X(LV,1,1,IEL) +
     $               S2(LV,1,1,IEL)*V2Y(LV,1,1,IEL) +
     $               S3(LV,1,1,IEL)*V2Z(LV,1,1,IEL) ) *
     $               R3MASK(LV,1,1,IEL)
            R1(LV,1,1,IEL) = RNOR*VNX(LV,1,1,IEL) +
     $                       RTN1*V1X(LV,1,1,IEL) +
     $                       RTN2*V2X(LV,1,1,IEL)
            R2(LV,1,1,IEL) = RNOR*VNY(LV,1,1,IEL) +
     $                       RTN1*V1Y(LV,1,1,IEL) +
     $                       RTN2*V2Y(LV,1,1,IEL)
            R3(LV,1,1,IEL) = RNOR*VNZ(LV,1,1,IEL) +
     $                       RTN1*V1Z(LV,1,1,IEL) +
     $                       RTN2*V2Z(LV,1,1,IEL) 
            write(6,*) '=====EGMASK===',R1(LV,1,1,IEL) , R2(LV,1,1,IEL)
  120    CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CRMSK2 (R1,R2,S1,S2,R1MASK,R2MASK,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      COMMON /INDXFC/ MCRFC(4,6)
     $              , MFCCR(3,8)
     $              , MEGCR(3,8)
     $              , MFCEG(2,12)
     $              , MCREG(2,12)
     $              , MCRRST(3,8)
     $              , MIDRST(3,12)
     $              , MCRIND(8)
     $              , MEDIND(2,4)
     $              , NTEFC(2,12)
     $              , NTCRF(2,3)     
      DIMENSION R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , S1(LX1,LY1,LZ1,1)
     $        , S2(LX1,LY1,LZ1,1)
     $        , R1MASK(LX1,LY1,LZ1,1)
     $        , R2MASK(LX1,LY1,LZ1,1)
C
      NCRNR = 2**NDIM
C
      DO 100 IEL=1,NEL
      DO 100 ICR=1,NCRNR
         IF (.NOT.IFMSCR(ICR,IEL,IFIELD)) GO TO 100
         IX = MCRRST(1,ICR)
         IY = MCRRST(2,ICR)
         IZ = 1
         RNOR = ( S1(IX,IY,IZ,IEL)*VNX(IX,IY,IZ,IEL) +
     $            S2(IX,IY,IZ,IEL)*VNY(IX,IY,IZ,IEL) ) *
     $            R1MASK(IX,IY,IZ,IEL)
         RTN1 = ( S1(IX,IY,IZ,IEL)*V1X(IX,IY,IZ,IEL) +
     $            S2(IX,IY,IZ,IEL)*V1Y(IX,IY,IZ,IEL) ) *
     $            R2MASK(IX,IY,IZ,IEL)
         R1(IX,IY,IZ,IEL) = RNOR*VNX(IX,IY,IZ,IEL) +
     $                      RTN1*V1X(IX,IY,IZ,IEL)
         R2(IX,IY,IZ,IEL) = RNOR*VNY(IX,IY,IZ,IEL) +
     $                      RTN1*V1Y(IX,IY,IZ,IEL)
            write(6,*) '==CRMASK2===', R1(IX,IY,IZ,IEL),R2(IX,IY,IZ,IEL)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE CRMSK3 (R1,R2,R3,S1,S2,S3,R1MASK,R2MASK,R3MASK,NEL)
C
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TSTEP'
      COMMON /INDXFC/ MCRFC(4,6)
     $              , MFCCR(3,8)
     $              , MEGCR(3,8)
     $              , MFCEG(2,12)
     $              , MCREG(2,12)
     $              , MCRRST(3,8)
     $              , MIDRST(3,12)
     $              , MCRIND(8)
     $              , MEDIND(2,4)
     $              , NTEFC(2,12)
     $              , NTCRF(2,3)     
      DIMENSION R1(LX1,LY1,LZ1,1)
     $        , R2(LX1,LY1,LZ1,1)
     $        , R3(LX1,LY1,LZ1,1)
     $        , S1(LX1,LY1,LZ1,1)
     $        , S2(LX1,LY1,LZ1,1)
     $        , S3(LX1,LY1,LZ1,1)
     $        , R1MASK(LX1,LY1,LZ1,1)
     $        , R2MASK(LX1,LY1,LZ1,1)
     $        , R3MASK(LX1,LY1,LZ1,1)
C
      NCRNR = 2**NDIM
C
      DO 100 IEL=1,NEL
      DO 100 ICR=1,NCRNR
         IF (.NOT.IFMSCR(ICR,IEL,IFIELD)) GO TO 100
         IX = MCRRST(1,ICR)
         IY = MCRRST(2,ICR)
         IZ = MCRRST(3,ICR)
         RNOR = ( S1(IX,IY,IZ,IEL)*VNX(IX,IY,IZ,IEL) +
     $            S2(IX,IY,IZ,IEL)*VNY(IX,IY,IZ,IEL) +
     $            S3(IX,IY,IZ,IEL)*VNZ(IX,IY,IZ,IEL) ) *
     $            R1MASK(IX,IY,IZ,IEL)
         RTN1 = ( S1(IX,IY,IZ,IEL)*V1X(IX,IY,IZ,IEL) +
     $            S2(IX,IY,IZ,IEL)*V1Y(IX,IY,IZ,IEL) +
     $            S3(IX,IY,IZ,IEL)*V1Z(IX,IY,IZ,IEL) ) *
     $            R2MASK(IX,IY,IZ,IEL)
         RTN2 = ( S1(IX,IY,IZ,IEL)*V2X(IX,IY,IZ,IEL) +
     $            S2(IX,IY,IZ,IEL)*V2Y(IX,IY,IZ,IEL) +
     $            S3(IX,IY,IZ,IEL)*V2Z(IX,IY,IZ,IEL) ) *
     $            R3MASK(IX,IY,IZ,IEL)
         R1(IX,IY,IZ,IEL) = RNOR*VNX(IX,IY,IZ,IEL) +
     $                      RTN1*V1X(IX,IY,IZ,IEL) +
     $                      RTN2*V2X(IX,IY,IZ,IEL)
         R2(IX,IY,IZ,IEL) = RNOR*VNY(IX,IY,IZ,IEL) +
     $                      RTN1*V1Y(IX,IY,IZ,IEL) +
     $                      RTN2*V2Y(IX,IY,IZ,IEL)
         R3(IX,IY,IZ,IEL) = RNOR*VNZ(IX,IY,IZ,IEL) +
     $                      RTN1*V1Z(IX,IY,IZ,IEL) +
     $                      RTN2*V2Z(IX,IY,IZ,IEL) 
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SETAXDY (IFAXDY)
C
      INCLUDE 'SIZE'
      INCLUDE 'DXYZ'
C
      LOGICAL IFAXDY
C
      IF (IFAXDY) THEN
         CALL COPY (DYM1 ,DAM1 ,NY1*NY1)
         CALL COPY (DYTM1,DATM1,NY1*NY1)
      ELSE
         CALL COPY (DYM1 ,DCM1 ,NY1*NY1)
         CALL COPY (DYTM1,DCTM1,NY1*NY1)
      ENDIF
C
      RETURN
      END
      SUBROUTINE SETAXW1 (IFAXWG)
C
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
C
      LOGICAL IFAXWG
C
      IF (IFAXWG) THEN
         CALL COPY (W3M1,W2AM1,NX1*NY1)
      ELSE
         CALL COPY (W3M1,W2CM1,NX1*NY1)
      ENDIF
C
      RETURN
      END
      SUBROUTINE SETAXW2 (IFAXWG)
C
      INCLUDE 'SIZE'
      INCLUDE 'WZ'
C
      LOGICAL IFAXWG
C
      IF (IFAXWG) THEN
         CALL COPY (W3M2,W2AM2,NX2*NY2)
      ELSE
         CALL COPY (W3M2,W2CM2,NX2*NY2)
      ENDIF
C
      RETURN
      END
