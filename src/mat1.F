!> \file mat1.F 
!! 
!! \brief matrix utility routines
!! 
!! This file contains various matrix utility routines, such as copy 

      SUBROUTINE BLANK(A,N)
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
C
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE

      RETURN
      END
C
      SUBROUTINE VSQ (A,N)
      DIMENSION  A(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vsq   '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I = 1, N
 100     A(I) = A(I)**2
      RETURN
      END
C
      SUBROUTINE VSQRT(A,N)
      DIMENSION  A(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vsqrt '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I = 1, N
 100     A(I) = SQRT(A(I))
      RETURN
      END
C
      SUBROUTINE ADD2S1(A,B,C1,N)
      DIMENSION A(1),B(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s1'
      endif
      isbcnt = 2*N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE ADD2S2(A,B,C1,N)
      DIMENSION A(1),B(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s2'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE ADD3S2(A,B,C,C1,C2,N)
      DIMENSION A(1),B(1),C(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add3s2'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE ADD3S3(A,B,C,D,C1,C2,C3,N)
      DIMENSION A(1),B(1),C(1),D(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add3s2'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)+C3*D(I)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE ADD4(A,B,C,D,N)
      REAL A(1),B(1),C(1),D(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add4  '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=B(I)+C(I)+D(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE INVERS2(A,B,N)
      REAL A(1),B(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'inver2'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=1./B(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE INVCOL1(A,N)
      REAL A(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl1'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=1./A(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE INVCOL2(A,B,N)
C
      REAL A(1),B(1)
      include 'CTIMER'
      include 'OPCTR'
C
      if (icalld.eq.0) tinvc=0.0
      icalld=icalld+1
      ninvc=icalld
      etime1=dclock()
C
C
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl2'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)/B(I)
 100  CONTINUE
      tinvc=tinvc+(dclock()-etime1)
      RETURN
      END
C
      SUBROUTINE INVCOL3(A,B,C,N)
      REAL A(1),B(1),C(1)
C
      include 'OPCTR'
      include 'CTIMER'
      if (icalld.eq.0) tinv3=0.0
      icalld=icalld+1
      ninv3=icalld
      etime1=dclock()
C
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl3'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=B(I)/C(I)
 100  CONTINUE
      tinv3=tinv3+(dclock()-etime1)
      RETURN
      END
C
      SUBROUTINE COL4(A,B,C,D,N)
      REAL A(1),B(1),C(1),D(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col4  '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=B(I)*C(I)*D(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE ADDCOL3(A,B,C,N)
      REAL A(1),B(1),C(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE ADDCOL4(A,B,C,D,N)
      REAL A(1),B(1),C(1),D(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl4'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE ASCOL5 (A,B,C,D,E,N)
      REAL A(1),B(1),C(1),D(1),E(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ascol5'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I) = B(I)*C(I)-D(I)*E(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE SUB2(A,B,N)
      REAL A(1),B(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub2  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE SUB3(A,B,C,N)
      REAL A(1),B(1),C(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub3  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=B(I)-C(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE SUBCOL3(A,B,C,N)
      REAL A(1),B(1),C(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'subcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE SUBCOL4(A,B,C,D,N)
      REAL A(1),B(1),C(1),D(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'subcl4'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)*D(I)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE CZERO(A,N)
      COMPLEX    A(1)
      DO 100 I = 1, N
 100     A(I ) = (0,0)
      RETURN
      END

      SUBROUTINE RZERO(A,N)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      RETURN
      END
C
      SUBROUTINE IZERO(A,N)
      INTEGER A(1)
C
      DO 100 I = 1, N
 100     A(I ) = 0
      RETURN
      END

      SUBROUTINE IZERO8(A,N)
      INTEGER*8 A(1)
C
      DO 100 I = 1, N
 100     A(I ) = 0
      RETURN
      END
C
      SUBROUTINE IONE(A,N)
      INTEGER   A(1)
      DO 100 I = 1, N
 100     A(I ) = 1
      RETURN
      END
C
      SUBROUTINE RONE(A,N)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 1.0
      RETURN
      END
C
C> Function comment blah 
C> Another comment line
      SUBROUTINE CFILL(A,B,N)
C> input parameter 
      DIMENSION  A(1)
C
      DO 100 I = 1, N
 100     A(I) = B
      RETURN
      END
C
      SUBROUTINE IFILL(IA,IB,N)
      DIMENSION IA(1)
C
      DO 100 I = 1, N
 100     IA(I) = IB
      RETURN
      END
C

!> \brief Copy N REAL elements from B to A 
!! 
!! this routine copies some stuff 
!!
!! @param B the source for the copy
!! @param A the destination for the copy 
!! @todo make this more general than just working with REAL
      SUBROUTINE COPY(A,B,N)
      REAL A(1), B(1)
      include 'CTIMER'
C
      if (icalld.eq.0) tcopy=0.0
      icalld=icalld+1
      ncopy=icalld
      etime1=dclock()
C
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      tcopy=tcopy+(dclock()-etime1)
      RETURN
      END
C
      SUBROUTINE CHCOPY(A,B,N)
      CHARACTER*1 A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      RETURN
      END
C
      SUBROUTINE ICOPY(A,B,N)
      INTEGER A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      RETURN
      END

      SUBROUTINE ICOPY8(A,B,N)
      INTEGER*8 A(1), B(1)
C
      DO 100 I = 1, N
 100     A(I) = B(I)
      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine copyi4(a,b,n)
      integer a(1)
      real    b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine opcopy (a1,a2,a3,b1,b2,b3)
      include 'SIZE'
      REAL A1(1),A2(1),A3(1),B1(1),B2(1),B3(1)
      NTOT1=NX1*NY1*NZ1*NELV
      CALL COPY(A1,B1,NTOT1)
      CALL COPY(A2,B2,NTOT1)
      IF(NDIM.EQ.3)CALL COPY(A3,B3,NTOT1)
      return
      END
c-----------------------------------------------------------------------

      SUBROUTINE CHSIGN(A,N)
      REAL A(1)
C
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE CMULT(A,CONST,N)
      REAL A(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cmult '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)*CONST
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE CADD(A,CONST,N)
      REAL A(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cadd  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      RETURN
      END
      SUBROUTINE IADD(I1,ISCAL,N)
      DIMENSION I1(1)
C
      DO 10 I=1,N
         I1(I)=I1(I)+ISCAL
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CADD2(A,B,CONST,N)
      REAL A(1),B(1)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cadd2 '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=B(I)+CONST
 100  CONTINUE
      RETURN
      END
C
      REAL FUNCTION VLMIN(VEC,N)
      REAL VEC(1)
      TMIN = 99.0E20
C
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      VLMIN = TMIN
      RETURN
      END
C
      integer function ivlmin(vec,n)
      integer vec(1),tmin
      if (n.eq.0) then
         ivlmin=0
         return
      endif
      tmin = 8888888
      do i=1,n
         tmin = min(tmin,vec(i))
      enddo
      ivlmin = tmin
      return
      end
C
      integer function ivlmax(vec,n)
      integer vec(1),tmax
      if (n.eq.0) then
         ivlmax=0
         return
      endif
      TMAX =-8888888
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      Ivlmax = tmax
      return
      end
C
      REAL FUNCTION VLMAX(VEC,N)
      REAL VEC(1)
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      VLMAX = TMAX
      RETURN
      END
C
      REAL FUNCTION VLAMAX(VEC,N)
      REAL VEC(1)
      TAMAX = 0.0
C
      DO 100 I=1,N
         TAMAX = MAX(TAMAX,ABS(VEC(I)))
 100  CONTINUE
      VLAMAX = TAMAX
      RETURN
      END
C
      REAL FUNCTION VLSUM(VEC,N)
      REAL VEC(1)
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vlsum '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      SUM = 0.
C
      DO 100 I=1,N
         SUM=SUM+VEC(I)
 100  CONTINUE
      VLSUM = SUM
      RETURN
      END
C
      SUBROUTINE VCROSS (U1,U2,U3,V1,V2,V3,W1,W2,W3,N)
C
C     Compute a Cartesian vector cross product.
C
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
      DIMENSION W1(1),W2(1),W3(1)
C
C
      DO 100 I=1,N
         U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
         U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
         U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE VDOT2 (DOT,U1,U2,V1,V2,N)
C
C     Compute a Cartesian vector dot product. 2-d version
C
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1)
      DIMENSION V1(1),V2(1)
C
C
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) 
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE VDOT3 (DOT,U1,U2,U3,V1,V2,V3,N)
C
C     Compute a Cartesian vector dot product. 3-d version
C
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
C
C
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE ADDTNSR(S,H1,H2,H3,NX,NY,NZ)
C
C     Map and add to S a tensor product form of the three functions H1,H2,H3.
C     This is a single element routine used for deforming geometry.
C
      DIMENSION H1(1),H2(1),H3(1)
      DIMENSION S(NX,NY,NZ)
C
      DO 200 IZ=1,NZ
      DO 200 IY=1,NY
         HH = H2(IY)*H3(IZ)
         DO 100 IX=1,NX
            S(IX,IY,IZ)=S(IX,IY,IZ)+HH*H1(IX)
  100    CONTINUE
  200 CONTINUE
      RETURN
      END
      FUNCTION LTRUNC(STRING,L)
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
      DATA BLNK/' '/
C
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      RETURN
      END
      FUNCTION MOD1(I,N)
C
C     Yields MOD(I,N) with the exception that if I=K*N, result is N.
C
      MOD1=0
      IF (I.EQ.0) THEN
         RETURN
      ENDIF
      IF (N.EQ.0) THEN
         WRITE(6,*) 
     $  'WARNING:  Attempt to take MOD(I,0) in FUNCTION MOD1.'
         RETURN
      ENDIF
      II = I+N-1
      MOD1 = MOD(II,N)+1
      RETURN
      END
      INTEGER FUNCTION LOG2(K)
      RK=(K)
      RLOG=LOG10(RK)
      RLOG2=LOG10(2.0)
      RLOG=RLOG/RLOG2+0.5
      LOG2=INT(RLOG)
      RETURN
      END
      SUBROUTINE IFLIP(I1,N)
      DIMENSION I1(1)
      N1=N+1
      N2=N/2
      DO 10 I=1,N2
         ILAST=N1-I
         ITMP=I1(ILAST)
         I1(ILAST)=I1(I)
         I1(I)=ITMP
   10 CONTINUE
      RETURN
      END
      SUBROUTINE ISWAP(B,IND,N,TEMP)
      INTEGER B(1),IND(1),TEMP(1)
C***
C***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
C***  INTO ITEM(I), WHERE JJ=IND(I).
C***
      DO 20 I=1,N
         JJ=IND(I)
         TEMP(I)=B(JJ)
   20 CONTINUE
      DO 30 I=1,N
   30 B(I)=TEMP(I)
      RETURN
      END
C
C----------------------------------------------------------------------------
C
C     Vector reduction routines which require communication 
C     on a parallel machine. These routines must be substituted with
C     appropriate routines which take into account the specific architecture.
C
C----------------------------------------------------------------------------
C
      FUNCTION GLSUM (X,N)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      GLSUM = TMP(1)
      RETURN
      END

      REAL FUNCTION GLAMAX(A,N)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMAX = 0.0
      DO 100 I=1,N
         TMAX = MAX(TMAX,ABS(A(I)))
 100  CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLAMAX=ABS(TMP(1))
      RETURN
      END
c-----------------------------------------------------------------------
      function iglmin(a,n)
      integer a(1),tmin
      integer tmp(1),work(1)
      tmin=  7777777
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call igop(tmp,work,'m  ',1)
      iglmin=tmp(1)
      if (iglmin .eq. 7777777) iglmin=0
      return
      end
c-----------------------------------------------------------------------
      function iglmax(a,n)
      integer a(1),tmax
      integer tmp(1),work(1)
      tmax= -999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call igop(tmp,work,'M  ',1)
      iglmax=tmp(1)
      if (iglmax .eq. -999999) iglmax=0
      return
      end
c-----------------------------------------------------------------------
      function iglsum(a,n)
      integer a(1),tsum
      integer tmp(1),work(1)
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call igop(tmp,work,'+  ',1)
      iglsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      FUNCTION GLMAX(A,N)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLMAX=TMP(1)
      RETURN
      END
C
      FUNCTION GLMIN(A,N)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL GOP(TMP,WORK,'m  ',1)
      GLMIN = TMP(1)
      RETURN
      END
C
      SUBROUTINE GLLOG(LA,LB)
C
C     If ANY LA=LB, then ALL LA=LB.
C
      LOGICAL LA,LB
      DIMENSION TMP(1),WORK(1)
C
      TMP(1)=1.0
      IF (LB) THEN
         IF (LA) TMP(1)=0.0
      ELSE
         IF (.NOT.LA) TMP(1)=0.0
      ENDIF
      CALL GOP(TMP,WORK,'*  ',1)
      IF (TMP(1).EQ.0.0) LA=LB
      RETURN
      END
C
      FUNCTION GLSC3(A,B,MULT,N)
C
C     Perform inner-product in double precision
C
      REAL A(1),B(1),MULT(1)
      REAL DT,TMP(1),WORK(1)
      REAL Ta,Tb,Tm
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc3 '
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DT = 0.0
      TMP(1) = 0.0

      DO 10 I=1,N
         Ta = A(I)
         Tb = B(I)
         Tm = MULT(I)
         DT = DT+Ta*Tb*Tm
c        TMP(1) = TMP(1)+ A(I)*B(I)*MULT(I)
 10   CONTINUE

      TMP(1)=DT
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC3 = TMP(1)

      RETURN
      END
C
      FUNCTION GLSC2(X,Y,N)
C
C     Perform inner-product in double precision
C
      DIMENSION X(1), Y(1)
      DIMENSION TMP(1),WORK(1)
      real*8 ds,dx,dy
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
c     tmp(1)=sdot(n,x,1,y,1)
      ds = 0.0
      do 10 i=1,n
         dx=x(i)
         dy=y(i)
         ds=ds+dx*dy
   10 continue
      tmp(1)=ds
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC2 = TMP(1)
      RETURN
      END

      FUNCTION FMDIAN(A,N,IFOK)
C     find the Median of the (global) set A
      INCLUDE 'SIZE'
      DIMENSION A(1)
      DIMENSION WORK1(5),WORK2(5)
      DIMENSION GUES(100)
      LOGICAL IFOK
C
      AMP  =1.5
      AFAC =1.5
      GMIN =GLMIN(A,N)
      GMAX =GLMAX(A,N)
      GMIN0=GLMIN(A,N)
      GMAX0=GLMAX(A,N)
      GUESS=(GMAX+GMIN)/2.0
      EPS  =(GMAX-GMIN)
      IF (EPS.EQ.0.0) THEN
         FMDIAN=GMAX
         RETURN
      ENDIF
      WORK1(1)=N
      CALL GOP(WORK1,WORK2,'+  ',1)
      NTOT=WORK1(1)
      N2 = (NTOT+1)/2
      IF (.NOT.IFOK) THEN
        WRITE(6,8) NID,N,(A(I),I=1,N)
        WRITE(6,9) NID,NTOT,N2,N,GMIN,GMAX
    8   FORMAT(I5,'N,A:',I5,10(6F10.5,/)) 
    9   FORMAT(I5,'mnx:',3I6,2F10.5)
      ENDIF
C
C     This is the trial loop
C
      ITRY=-1
   10 CONTINUE
      ITRY=ITRY+1
      II=ITRY+1
      IF (II.LE.100) GUES(II)=GUESS
C     error check for infinite loop
      IF (ITRY.GT.2*NTOT) GOTO 9000
      CALL RZERO(WORK1,5)
      NLT=0
      NGT=0
      CLT=GMIN0
      CGT=GMAX0
      DO 100 I=1,N
         AA=A(I)
         IF (AA.NE.GUESS) THEN
            IF (AA.LT.GUESS) THEN
               NLT=NLT+1
C              CLT - closest value to GUESS Less Than GUESS
               IF (AA.GT.CLT) CLT=AA
            ENDIF
            IF (AA.GT.GUESS) THEN
               NGT=NGT+1
C              CGT - closest value to GUESS Greater Than GUESS
               IF (AA.LT.CGT) CGT=AA
            ENDIF
            DUM=1./(EPS+ABS(AA-GUESS))
            WORK1(1)=WORK1(1)+DUM
            WORK1(2)=WORK1(2)+DUM*AA
         ELSE
C           detected values equaling the guess.
            WORK1(5)=WORK1(5)+1.0
         ENDIF
  100 CONTINUE
C     Invoke vector reduction across processors:
      WORK2(1)=CLT
      CLT=GLMAX(WORK2,1)
      WORK2(1)=CGT
      CGT=GLMIN(WORK2,1)
      WORK1(3)=NLT
      WORK1(4)=NGT
      CALL GOP(WORK1,WORK2,'+  ',5)
      NLT=WORK1(3)
      NGT=WORK1(4)
      IF (.NOT.IFOK) THEN
         WRITE(6,101) NID,GUESS,CLT,CGT
         WRITE(6,102) NID,(WORK1(I),I=1,5)
  101    FORMAT(I5,'Glg:',3F12.5)
  102    FORMAT(I5,'WORK1:',5F12.5)
      ENDIF
C
C     Done?
C
      IF (NLT.GT.N2.OR.NGT.GT.N2) THEN
C        we're not done.....
         IF (NGT.GT.NLT) THEN
C           guess is too low
            GMIN=CGT
            G2=CGT+MAX(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2.GT.GMAX) G2=0.5*(GUESS+GMAX)
            EPS=AFAC*ABS(G2-GUESS)
C           see that we move at least as far as the next closest value.
            GUESS=MAX(G2,CGT)
            GOTO 10
         ELSE IF (NLT.GT.NGT) THEN
C           guess is too high
            GMAX=CLT
            G2=CLT+MIN(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2.LT.GMIN) G2=0.5*(GUESS+GMIN)
            EPS=AFAC*ABS(G2-GUESS)
C           see that we move at least as far as the next closest value.
            GUESS=MIN(G2,CLT)
            GOTO 10
         ENDIF
      ELSE
C
C        we're done....
         IF (WORK1(5).NE.0) THEN
C           the median is (usually) one of the values 
            FMDIAN=GUESS
            IF (WORK1(5).EQ.1.0) THEN
               IF (MOD(NTOT,2).EQ.0) THEN
                  IF (NGT.GT.NLT) THEN
                     FMDIAN=0.5*(GUESS+CGT)
                  ELSE
                     FMDIAN=0.5*(GUESS+CLT)
                  ENDIF
               ELSE
                  IF (NGT.EQ.NLT) THEN
                     FMDIAN=GUESS
                  ELSE IF(NGT.GT.NLT) THEN
                     FMDIAN=CGT
                  ELSE
                     FMDIAN=CLT
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF (MOD(NTOT,2).EQ.0) THEN
               IF (NGT.EQ.NLT) THEN
                  FMDIAN=0.5*(CLT+CGT)
               ELSE IF(NGT.GT.NLT) THEN
                  FMDIAN=0.5*(GUESS+CGT)
               ELSE
                  FMDIAN=0.5*(GUESS+CLT)
               ENDIF
            ELSE
               IF (NGT.EQ.NLT) THEN
                  FMDIAN=GUESS
               ELSE IF(NGT.GT.NLT) THEN
                  FMDIAN=CGT
               ELSE
                  FMDIAN=CLT
               ENDIF
           ENDIF
         ENDIF
C
      ENDIF
       IF (.NOT.IFOK) WRITE(6,*) NID,'FMDIAN2',FMDIAN,(A(I),I=1,N)
      RETURN
C
C     Error handling
C
 9000 CONTINUE
      WRITE(6,11) NTOT,GMIN0,GMAX0,GUESS
   11 FORMAT('ABORTING IN FMDIAN: N,AMIN,AMAX:',I6,3G14.6)
      DO 13 I1=1,N,5
        IN=I1+5 
        IN=MIN(IN,N)
        WRITE(6,12) NID,(A(I),I=I1,IN)
   12   FORMAT(I4,' FMA:',5G14.6)
   13 CONTINUE
      DO 15 I1=1,ITRY,5
        IN=I1+5
        IN=MIN(IN,ITRY)
        WRITE(6,14) NID,(GUES(I),I=I1,IN)
   14   FORMAT(I4,' FMG:',5G14.6)
   15 CONTINUE
      STOP
      END
C
C========================================================================
C     Double precision matrix and vector routines
C========================================================================
C
      SUBROUTINE DCADD(A,CONST,N)
      real*8 A(1),CONST
C
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE DSUB2(A,B,N)
      real*8 A(1), B(1)
C
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE DADD2(A,B,N)
      real*8 A(1), B(1)
C
      DO 100 I=1,N
         A(I)=A(I)+B(I)
 100  CONTINUE
      RETURN
      END
      SUBROUTINE CHSWAPR(B,l,IND,N,TEMP)
      INTEGER IND(1)
      CHARACTER*6 B(1),TEMP(1)
C***
C***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
C***  INTO ITEM(I), WHERE JJ=IND(I).
C***
      DO 20 I=1,N
         JJ=IND(I)
         TEMP(I)=B(JJ)
   20 CONTINUE
      DO 30 I=1,N
   30 B(I)=TEMP(I)
      RETURN
      END
      SUBROUTINE DRCOPY(R,D,n)
      real*8    d(1)
      dimension r(1)
      do 10 i=1,n
         r(i)=d(i)
   10 continue
      return
      end
      subroutine sorts(xout,xin,work,n)
      real xout(1),xin(1),work(1)
      call copy(xout,xin,n)
      call sort(xout,work,n)
      return
      end
C
      FUNCTION IVLSUM(A,N)
      INTEGER A(1)
      INTEGER TSUM
      if (n.eq.0) then
         ivlsum = 0
         return
      endif
      TSUM=A(1)
      DO 100 I=2,N
         TSUM=TSUM+A(I)
  100 CONTINUE
      IVLSUM=TSUM
      RETURN
      END
C
      SUBROUTINE ICADD(A,C,N)
      INTEGER A(1),C
      DO 100 I = 1, N
 100     A(I) = A(I) + C
      RETURN
      END
      subroutine isort(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(1),ind(1)
      integer aa
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end

      subroutine sort(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real a(1),aa
      integer ind(1)
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
C
c-----------------------------------------------------------------------
      REAL FUNCTION VLSC2(X,Y,N)
      REAL X(1),Y(1)
      include 'SIZE'
      include 'OPCTR'
      include 'PARALLEL'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      s = 0.
      do i=1,n
         s = s + x(i)*y(i)
      enddo
      vlsc2=s
      return
      end
c-----------------------------------------------------------------------
      real function vlsc21(x,y,n)
      real x(1),y(1)
      include 'SIZE'
      include 'OPCTR'
      include 'PARALLEL'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC21'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      s = 0.
      do i=1,n
         s = s + x(i)*x(i)*y(i)
      enddo
      vlsc21=s
      return
      end
c-----------------------------------------------------------------------
      subroutine iswap_ip(x,p,n)
      integer x(1),xstart
      integer p(1)
c
c     In-place permutation: x' = x(p)
c
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! iswap_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine iswapt_ip(x,p,n)
      integer x(1),t1,t2
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine swap_ip(x,p,n)
      real    x(1),xstart
      integer p(1)
c
c     In-place permutation: x' = x(p)
c
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! swap_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine swapt_ip(x,p,n)
      real    x(1),t1,t2
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! swapt_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ident(a,n)
      real  a(n,n)
      call rzero(a,n*n)
      do i=1,n
         a(i,i) = 1.0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
c
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine gaujordf(a,m,n,indr,indc,ipiv,ierr,rmult)
C
C     Gauss-Jordan matrix inversion with full pivoting
c
c     Num. Rec. p. 30, 2nd Ed., Fortran
c
C
C     a     is an m x n matrix
C     rmult is a  work array of dimension m
C
c
      real a(m,n),rmult(m)
      integer indr(m),indc(n),ipiv(n)
c
      ierr = 0
      eps = 1.e-9
c
      call izero(ipiv,m)
c
      do k=1,m
         amx=0.
c
c        Pivot search
         do i=1,m
            if (ipiv(i).ne.1) then
               do j=1,m
                  if (ipiv(j).eq.0) then
                    if (abs(a(i,j)).ge.amx) then
                       amx = abs(a(i,j))
                       ir  = i
                       jc  = j
                    endif
                 elseif (ipiv(j).gt.1) then
                    ierr = -ipiv(j)
                    return
                 endif
              enddo
           endif
        enddo
        ipiv(jc) = ipiv(jc) + 1
c
c       Swap rows
        if (ir.ne.jc) then
           do j=1,n
              tmp     = a(ir,j)
              a(ir,j) = a(jc,j)
              a(jc,j) = tmp
           enddo
        endif
        indr(k)=ir
        indc(k)=jc
c       write(6 ,*) k,' Piv:',jc,a(jc,jc)
c       write(28,*) k,' Piv:',jc,a(jc,jc)
        if (abs(a(jc,jc)).lt.eps) then
           write(6 ,*) 'small Gauss Jordan Piv:',jc,a(jc,jc)
           write(28,*) 'small Gauss Jordan Piv:',jc,a(jc,jc)
           ierr = jc
           return
        endif
        piv = 1./a(jc,jc)
        a(jc,jc)=1.
        do j=1,n
           a(jc,j) = a(jc,j)*piv
        enddo
c
        do j=1,n
           work    = a(jc,j)
           a(jc,j) = a(1 ,j)
           a(1 ,j) = work
        enddo
        do i=2,m
           rmult(i) = a(i,jc)
           a(i,jc)  = 0.
        enddo
c
        do j=1,n
        do i=2,m
           a(i,j) = a(i,j) - rmult(i)*a(1,j)
        enddo
        enddo
c
        do j=1,n
           work    = a(jc,j)
           a(jc,j) = a(1 ,j)
           a(1 ,j) = work
        enddo

      enddo
 
c     unscramble matrix
      do j=m,1,-1
         if (indr(j).ne.indc(j)) then
            do i=1,m
               tmp=a(i,indr(j))
               a(i,indr(j))=a(i,indc(j))
               a(i,indc(j))=tmp
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE RZERO3 (A,B,C,N)
      DIMENSION A(1),B(1),C(1)
      DO 100 I=1,N
         A(I)=0.0
         B(I)=0.0
         C(I)=0.0
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE INVCHK2 (A,B,N)
      DIMENSION A(1),B(1)
      DO 100 I=1,N
         IF (B(I) .GT. 1.e-20) THEN
             A(I)=A(I)/B(I)
         ELSE
             A(I)=0.0
         ENDIF
 100  CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
      CALL DSSET (NX1,NY1,NZ1)
      IFACE  = EFACE1(IFC)
      JS1    = SKPDAT(1,IFACE)
      JF1    = SKPDAT(2,IFACE)
      JSKIP1 = SKPDAT(3,IFACE)
      JS2    = SKPDAT(4,IFACE)
      JF2    = SKPDAT(5,IFACE)
      JSKIP2 = SKPDAT(6,IFACE)
C
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX132(S1,S2,L2)
      CHARACTER*132 S1,S2
C
      N1=80-L2+1
      INDX132=0
      IF (N1.LT.1) RETURN
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX132=I
            RETURN
         ENDIF
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      integer function indx2(s1,l1,s2,l2)
      character*132 s1,s2

      n1=l1-l2+1
      indx2=0
      if (n1.lt.1) return

      do i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx2=i
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*132 S1,S2
C
      N1=132-L2+1
      INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine csplit(s0,s1,s2,l0)
      CHARACTER*132 S0,S1,S2
C     split string S1 into two parts, delimited by S2.
C
      I2=INDX_CUT(S1,S2,L0)
      IF (I2.EQ.0) return
C
      I1=I2-1
      CALL BLANK(S0,132)
      S0(1:I1)=S1(1:I1)
      CALL LSHFT(S1,I2)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine lshft(string,ipt)
C     shift string from IPT to the left
C     INPUT : "abcde......    test    "
C     OUTPUT: "e......    test        "     if ipt.eq.5
      CHARACTER*1 STRING(132)
C
      DO 20 J=1,133-IPT
         IJ=IPT+J-1
         STRING(J)=STRING(IJ)
   20 CONTINUE
      DO 30 J=134-IPT,132
         STRING(J)=' '
   30 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ljust(string)
C     left justify string
      CHARACTER*1 STRING(132)
C
      IF (STRING(1).NE.' ') return
C
      DO 100 I=2,132
C
         IF (STRING(I).NE.' ') THEN
            DO 20 J=1,133-I
               IJ=I+J-1
               STRING(J)=STRING(IJ)
   20       CONTINUE
            DO 30 J=134-I,132
               STRING(J)=' '
   30       CONTINUE
            return
         ENDIF
C
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE UNITVEC (X,Y,Z,N)
      DIMENSION X(1),Y(1),Z(1)
      DO 100 I=1,N
      XLNGTH = SQRT( X(I)**2 + Y(I)**2 + Z(I)**2 )
      IF (XLNGTH.NE.0.0) THEN
         X(I) = X(I)/XLNGTH
         Y(I) = Y(I)/XLNGTH
         Z(I) = Z(I)/XLNGTH
      ENDIF
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CAPIT(LETTRS,N)
C     Capitalizes string of length n
      CHARACTER LETTRS(N)
C
      DO 5 I=1,N
         INT=ICHAR(LETTRS(I))
         IF(INT.GE.97 .AND. INT.LE.122) THEN
            INT=INT-32
            LETTRS(I)=CHAR(INT)
         ENDIF
5     CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine irank_vec(ind,nn,a,m,n,key,nkey,aa)
c
c     Compute rank of each unique entry a(1,i) 
c
c     Output:   ind(i)  i=1,...,n    (global) rank of entry a(*,i)
c               nn  = max(rank)
c               a(j,i) is destroyed
c
c     Input:    a(j,i) j=1,...,m;  i=1,...,n  
c               m      :   leading dim. of v  (ldv must be .ge. m)
c               key    :   sort key
c               nkey   :   
c
c     Although not mandatory, this ranking procedure is probably
c     most effectively employed when the keys are pre-sorted. Thus,
c     the option is provided to sort vi() prior to the ranking.
c
c
      integer ind(n),a(m,n)
      integer key(nkey),aa(m)
      logical iftuple_ianeb,a_ne_b
c
      if (m.eq.1) then
c
         write(6,*) 
     $        'WARNING: For single key, not clear that rank is unique!'
         call irank(a,ind,n)
         return
      endif
c
c
      nk = min(nkey,m)
      call ituple_sort(a,m,n,key,nk,ind,aa)
c
c     Find unique a's
c
      nn=1
c
      call icopy(aa,a,m)
      a(1,1) = nn
      a(2,1)=ind(1)
c
      do i=2,n
         a_ne_b = iftuple_ianeb(aa,a(1,i),key,nk)
         if (a_ne_b) then
            call icopy(aa,a(1,i),m)
            nn = nn+1
         endif
         a(1,i) = nn
         a(2,i) = ind(i)
      enddo
c
c     Set ind() to rank
c
      do i=1,n
         iold=a(2,i)
         ind(iold) = a(1,i)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine ituple_sort(a,lda,n,key,nkey,ind,aa)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_ialtb
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
c           aa  = a  (l)
            call icopy(aa,a(1,l),lda)
            ii  = ind(l)
         else
c           aa =   a(ir)
            call icopy(aa,a(1,ir),lda)
            ii = ind(ir)
c           a(ir) =   a( 1)
            call icopy(a(1,ir),a(1,1),lda)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
c              a(1) = aa
               call icopy(a(1,1),aa,lda)
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if (iftuple_ialtb(a(1,j),a(1,j+1),key,nkey)) j=j+1
            endif
            if (iftuple_ialtb(aa,a(1,j),key,nkey)) then
c              a(i) = a(j)
               call icopy(a(1,i),a(1,j),lda)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
c        a(i) = aa
         call icopy(a(1,i),aa,lda)
         ind(i) = ii
      GOTO 100
      end
c
c-----------------------------------------------------------------------
c
      subroutine tuple_sort(a,lda,n,key,nkey,ind,aa)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_altb
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
c           aa  = a  (l)
            call copy(aa,a(1,l),lda)
            ii  = ind(l)
         else
c           aa =   a(ir)
            call copy(aa,a(1,ir),lda)
            ii = ind(ir)
c           a(ir) =   a( 1)
            call copy(a(1,ir),a(1,1),lda)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
c              a(1) = aa
               call copy(a(1,1),aa,lda)
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
c              if ( a(j).lt.a(j+1) ) j=j+1
               if (iftuple_altb(a(1,j),a(1,j+1),key,nkey)) j=j+1
            endif
c           if (aa.lt.a(j)) then
            if (iftuple_altb(aa,a(1,j),key,nkey)) then
c              a(i) = a(j)
               call copy(a(1,i),a(1,j),lda)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
c        a(i) = aa
         call copy(a(1,i),aa,lda)
         ind(i) = ii
      GOTO 100
      end
c
c-----------------------------------------------------------------------
c
      logical function iftuple_ialtb(a,b,key,nkey)
      integer a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).lt.b(k)) then
            iftuple_ialtb = .true.
            return
         elseif (a(k).gt.b(k)) then
            iftuple_ialtb = .false.
            return
         endif
      enddo
      iftuple_ialtb = .false.
      return
      end
c
c-----------------------------------------------------------------------
c
      logical function iftuple_altb(a,b,key,nkey)
      real a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).lt.b(k)) then
            iftuple_altb = .true.
            return
         elseif (a(k).gt.b(k)) then
            iftuple_altb = .false.
            return
         endif
      enddo
      iftuple_altb = .false.
      return
      end
c
c-----------------------------------------------------------------------
c
      logical function iftuple_ianeb(a,b,key,nkey)
      integer a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).ne.b(k)) then
            iftuple_ianeb = .true.
            return
         endif
      enddo
      iftuple_ianeb = .false.
      return
      end
c
c-----------------------------------------------------------------------
      subroutine irank(A,IND,N)
C
C     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
C
      integer A(1),IND(1)
C
      if (n.le.1) return
      DO 10 J=1,N
         IND(j)=j
   10 continue
C
      if (n.eq.1) return
      L=n/2+1
      ir=n
  100 continue
         IF (l.gt.1) THEN
            l=l-1
            indx=ind(l)
            q=a(indx)
         ELSE
            indx=ind(ir)
            q=a(indx)
            ind(ir)=ind(1)
            ir=ir-1
            if (ir.eq.1) then
               ind(1)=indx
               return
            endif
         ENDIF
         i=l
         j=l+l
  200    continue
         IF (J.le.IR) THEN
            IF (J.lt.IR) THEN
               IF ( A(IND(j)).lt.A(IND(j+1)) ) j=j+1
            ENDIF
            IF (q.lt.A(IND(j))) THEN
               IND(I)=IND(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            ENDIF
         GOTO 200
         ENDIF
         IND(I)=INDX
      GOTO 100
      END
c
c-----------------------------------------------------------------------
c
      subroutine iranku(r,input,n,w,ind)
c
c     Return the rank of each input value, and the maximum rank.
c
c     OUTPUT:    r(k) = rank of each entry,  k=1,..,n
c                maxr = max( r )
c                w(i) = sorted & compressed list of input values
c
      integer r(1),input(1),ind(1),w(1)
c
      call icopy(r,input,n)
      call isort(r,ind,n)
c
      maxr  = 1
      rlast = r(1) 
      do i=1,n
c        Bump rank only when r_i changes
         if (r(i).ne.rlast) then
            rlast = r(i)
            maxr  = maxr + 1
         endif
         r(i) = maxr
      enddo
      call iunswap(r,ind,n,w)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine iunswap(b,ind,n,temp)
      integer b(1),ind(1),temp(1)
c
c     sort associated elements by putting item(jj)
c     into item(i), where jj=ind(i).
c
      do 20 i=1,n
         jj=ind(i)
         temp(jj)=b(i)
   20 continue
      do 30 i=1,n
   30 b(i)=temp(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine irank_vec_tally(ind,nn,a,m,n,key,nkey,key2,aa)
c
c     Compute rank of each unique entry a(1,i) 
c
c     Output:   ind(i)  i=1,...,n    (global) rank of entry a(*,i)
c               nn  = max(rank)
c               a(j,i) is destroyed
c               a(1,i) tally of preceding structure values
c
c     Input:    a(j,i) j=1,...,m;  i=1,...,n  
c               m      :   leading dim. of v  (ldv must be .ge. m)
c               key    :   sort key
c               nkey   :   
c
c     Although not mandatory, this ranking procedure is probably
c     most effectively employed when the keys are pre-sorted. Thus,
c     the option is provided to sort vi() prior to the ranking.
c
c
      integer ind(n),a(m,n)
      integer key(nkey),key2(0:3),aa(m)
      logical iftuple_ianeb,a_ne_b
c
c
      nk = min(nkey,m)
      call ituple_sort(a,m,n,key,nk,ind,aa)
c     do i=1,n
c        write(6,*) i,' sort:',(a(k,i),k=1,3)
c     enddo
c
c
c     Find unique a's
c
      call icopy(aa,a,m)
      nn=1
      mm=0
c
      a(1,1) = nn
      a(2,1)=ind(1)
      a(3,1)=mm
c
      do i=2,n
         a_ne_b = iftuple_ianeb(aa,a(1,i),key,nk)
         if (a_ne_b) then              ! new structure
            ms = aa(3)                 ! structure type
            if (aa(2).eq.0) ms = aa(2) ! structure type
            mm = mm+key2(ms)           ! n dofs
            call icopy(aa,a(1,i),m)
            nn = nn+1
         endif
         a(1,i) = nn
         a(2,i) = ind(i)
         a(3,i) = mm
      enddo
      ms = aa(3)
      if (aa(2).eq.0) ms = aa(2) ! structure type
      nn = mm+key2(ms)
c
c     Set ind() to rank
c
      do i=1,n
         iold=a(2,i)
         ind(iold) = a(1,i)
      enddo
c
c     Set a1() to number of preceding dofs
c
      do i=1,n
         iold=a(2,i)
         a(1,iold) = a(3,i)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine col2s2(x,y,s,n)
      real x(n),y(n)
c
      do i=1,n
         x(i)=s*x(i)*y(i)
      enddo
c
      return
      end
c-----------------------------------------------------------------------


#ifdef _OPENACC
c-----------------------------------------------------------------------
      SUBROUTINE CHSIGN_ACC(A,N)
      REAL A(N)
C
!$ACC DATA PRESENT(A(1:N))
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA

      RETURN
      END
C
c-----------------------------------------------------------------------
      SUBROUTINE RZERO_ACC(A,N)
      REAL A(N)
C
!$ACC DATA PRESENT(A(1:N))
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I) = 0.
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA

      RETURN
      END
C
c-----------------------------------------------------------------------
      SUBROUTINE SUB3_ACC(A,B,C,N)
      REAL A(N),B(N),C(N)
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub3  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
!$ACC DATA PRESENT(A(1:N),B(1:N),C(1:N))
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I)=B(I)-C(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA

      RETURN
      END

C-----------------------------------------------------------------------
      FUNCTION GLMAX_ACC(A,N)
      REAL A(N)
      DIMENSION TMP(1),WORK(1)
      TMAX=-99.0e20
!$ACC DATA PRESENT(A(1:N))
!$ACC PARALLEL LOOP REDUCTION(MAX:TMAX)
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA

      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLMAX_ACC=TMP(1)
      RETURN
      END

C-----------------------------------------------------------------------
      REAL FUNCTION GLAMAX_ACC(A,N)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMAX = 0.0
!$ACC DATA PRESENT(A(1:N))
!$ACC PARALLEL LOOP REDUCTION(MAX:TMAX)
      DO 100 I=1,N
         TMAX = MAX(TMAX,ABS(A(I)))
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA

      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLAMAX_ACC=ABS(TMP(1))
      RETURN
      END

C-----------------------------------------------------------------------
      FUNCTION GLSC3_ACC(A,B,MULT,N)
C
C     Perform inner-product in double precision
C
      REAL A(N),B(N),MULT(N)
      REAL DT,TMP,WORK(1)
      REAL Ta,Tb,Tm
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc3 '
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
C      DT = 0.0
      TMP = 0.0

!$ACC DATA PRESENT(A(1:N))
!$ACC PARALLEL LOOP REDUCTION(+:TMP)
      DO 10 I=1,N
C         Ta = A(I)
C         Tb = B(I)
C         Tm = MULT(I)
C         DT = DT+Ta*Tb*Tm
        TMP = TMP + A(I)*B(I)*MULT(I)
 10   CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA

C      TMP(1)=DT
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC3_ACC = TMP

      RETURN
      END

#endif
