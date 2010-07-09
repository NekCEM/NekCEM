      SUBROUTINE COL2(A,B,N)
      DIMENSION A(1),B(1)
      include 'OPCTR'
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
c     IF (N.GT.36) THEN
c        N12=N-MOD(N,12)
c        CALL COL2A(A,B,N12)
c        DO 100 I=N12+1,N
c           A(I)=A(I)*B(I)
c 100    CONTINUE
c     ELSE
         DO 200 I=1,N
            A(I)=A(I)*B(I)
  200    CONTINUE
c     ENDIF
      RETURN
      END
      SUBROUTINE COL3(A,B,C,N)
      DIMENSION A(1),B(1),C(1)
      include 'OPCTR'
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
c     IF (N.GT.36) THEN
c        N12=N-MOD(N,12)
c        CALL COL3A(A,B,C,N12)
c        DO 100 I=N12*1,N
c           A(I)=B(I)*C(I)
c 100    CONTINUE
c     ELSE
         DO 200 I=1,N
            A(I)=B(I)*C(I)
  200    CONTINUE
c     ENDIF
      RETURN
      END
      SUBROUTINE ADD2(A,B,N)
      DIMENSION A(1),B(1)
      include 'OPCTR'
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
      DO 200 I=1,N
         A(I)=A(I)+B(I)
  200 CONTINUE
      RETURN
      END
      SUBROUTINE ADD3(A,B,C,N)
      DIMENSION A(1),B(1),C(1)
      include 'OPCTR'
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
c     IF (N.GT.36) THEN
c        N12=N-MOD(N,12)
c        CALL ADD3A(A,B,C,N12)
c        DO 100 I=N12+1,N
c           A(I)=B(I)+C(I)
c 100    CONTINUE
c     ELSE
         DO 200 I=1,N
            A(I)=B(I)+C(I)
  200    CONTINUE
c     ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine fcol3(x,y,z,n)
      real x(1),y(1),z(1)
      do i=1,n
         x(i) = y(i)*z(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine faddcol4(x,y,z,q,n)
      real x(1),y(1),z(1),q(1)
      do i=1,n
         x(i) = x(i) + y(i)*z(i)*q(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine fadd2(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = x(i)+y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      integer*8 function i8glmax(a,n)
      integer*8 a(1),tmax
      integer*8 tmp(1),work(1)
      tmax= -999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      !call i8gop(tmp,work,'M  ',1) temporary

      call igop(tmp,work,'M  ',1)
      i8glmax=tmp(1)
      if (i8glmax .eq. -999999) i8glmax=0
      return
      end

