      SUBROUTINE FCAVER(XAVER,A,IEL,IFACE1)
C------------------------------------------------------------------------
C
C     Compute the average of A over the face IFACE1 in element IEL.
C
C         A is a (NX,NY,NZ) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'GEOM'
      INCLUDE 'TOPOL'
      REAL A(LX1,LY1,LZ1,1)
C
      FCAREA = 0.
      XAVER  = 0.
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
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         FCAREA = FCAREA+AREA(I,1,IFACE1,IEL)
         XAVER  = XAVER +AREA(I,1,IFACE1,IEL)*A(J1,J2,1,IEL)
  100 CONTINUE
C
      XAVER = XAVER/FCAREA
      RETURN
      END
      SUBROUTINE FACCL2(A,B,IFACE1)
C
C     Collocate B with A on the surface IFACE1 of element IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         5 Jan 1989 15:12:22      PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1)
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
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = A(J1,J2,1)*B(I,1)
  100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE FACCL3(A,B,C,IFACE1)
C
C     Collocate B with A on the surface IFACE1 of element IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         5 Jan 1989 15:12:22      PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1)
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
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = B(J1,J2,1)*C(I,1)
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE FADDCL3(A,B,C,IFACE1)
C
C     Collocate B with C and add to A on the surface IFACE1 of element IE.
C
C         A is a (NX,NY,NZ) data structure
C         B is a (NX,NY,NZ) data structure
C         C is a (NX,NY,IFACE) data structure
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         29 Jan 1990 18:00 PST   PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1,LZ1),C(LX1,LY1)
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
C
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         A(J1,J2,1) = A(J1,J2,1) + B(J1,J2,1)*C(I,1)
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      LOGICAL FUNCTION IFGTIL(IVALUE,LINE)
C
C     Read IVALUE from LINE and set IFGTIL to .TRUE. if successful,
C                                   IFGTIL to .FALSE. otherwise.
C
C     This complicated function is necessary thanks to the Ardent,
C     which won't allow free formatted reads (*) from internal strings!
C
      CHARACTER*80 LINE
      CHARACTER*80 WORK
      CHARACTER*8  FMAT
C
      IFGTIL=.FALSE.
      IVALUE=0
C
      WORK=LINE
      CALL LJUST(WORK)
      IFLDW=INDX1(WORK,' ',1)-1
C
      IF (IFLDW.GT.0) THEN
         WRITE(FMAT,10) IFLDW
   10    FORMAT('(F',I3.3,'.0)')
         READ(WORK,FMAT,ERR=100,END=100) TVAL
         IVALUE=INT(TVAL)
         IFGTIL=.TRUE.
         RETURN
      ENDIF
C
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE VRNVERT(A,N)
      REAL*4 A(1)
      real*4 a4
      real*4 revert
c     DO 100 I=1,N
c        a4   = a(i)
c        a4   = REVERT(a4)
c        a(i) = a4
c 100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION FACDOT(A,B,IFACE1)
C
C     Take the dot product of A and B on the surface IFACE1 of element IE.
C
C         IFACE1 is in the preprocessor notation 
C         IFACE  is the dssum notation.
C         5 Jan 1989 15:12:22      PFF
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
      DIMENSION A(LX1,LY1,LZ1),B(LX1,LY1)
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
C
      SUM=0.0
      I = 0
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         I = I+1
         SUM = SUM + A(J1,J2,1)*B(I,1)
  100 CONTINUE
C
      FACDOT = SUM
C
      RETURN
      END
C
      SUBROUTINE NEKUVP (IEL)
C------------------------------------------------------------------
C
C     Generate user-specified material properties
C
C------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'PARALLEL'
      INCLUDE 'NEKUSE'
      INCLUDE 'EMWAVE'

      IELG = LGLEL(IEL)

      DO 10 K=1,NZ1
      DO 10 J=1,NY1
      DO 10 I=1,NX1

        !CALL NEKASGN (I,J,K,IEL )
         CALL USERVP  (I,J,K,IELG)

         EPSILO (I,J,K,IEL) = permit
         MU     (I,J,K,IEL) = permea
         DRUDEA (I,J,K,IEL) = drude_alpha
         DRUDEB (I,J,K,IEL) = drude_beta
c        LORENA(I,J,K,IEL,1) = loren_gamma1
c        LORENA(I,J,K,IEL,2) = loren_gamma2
c        LORENA(I,J,K,IEL,3) = loren_gamma3
c        LORENB(I,J,K,IEL,1) = loren_omega1
c        LORENB(I,J,K,IEL,2) = loren_omega2
c        LORENB(I,J,K,IEL,3) = loren_omega3
c        LORENC(I,J,K,IEL,1) = loren_eps1
c        LORENC(I,J,K,IEL,2) = loren_eps2
c        LORENC(I,J,K,IEL,3) = loren_eps3

 10   CONTINUE
      RETURN
      END
C
