c-----------------------------------------------------------------------
      subroutine setics
C-----------------------------------------------------------------------
C
C     Set initial conditions.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'DEALIAS'
      INCLUDE 'INPUT'
      INCLUDE 'IXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'MASS'
      INCLUDE 'MVGEOM'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
c
      logical  iffort(  ldimt1,0:lpert)
     $       , ifrest(0:ldimt1,0:lpert)
     $       , ifprsl(  ldimt1,0:lpert)
c
      LOGICAL  IFANYP
      common /rdump/ ntdump
      common /inelr/ nelrr
      common /ctmp1/ work(lx1,ly1,lz1,lelv)
     $ ,             ta1 (lx2,ly1,lz1)
     $ ,             ta2 (lx2,ly2,lz1)
      integer*8 ntotg,nn

      real psmax(LDIMT)

      if (nid.eq.0) then
          write(6,*) 'set initial conditions'
      endif
C
C     Initialize all fields:
C
      nxyz2=nx2*ny2*nz2
      ntot2=nxyz2*nelv
      nxyz1=nx1*ny1*nz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1
c
c
      CALL RZERO(VX,NTOTT)
      CALL RZERO(VY,NTOTT)
      CALL RZERO(VZ,NTOTT)
      CALL RZERO(PR,nxyz2*nelt)
      DO 10 IFLD=1,LDIMT
         CALL RZERO(T(1,1,1,1,IFLD),NTOTT)
   10 CONTINUE

      jp = 0                  ! set counter for perturbation analysis

      irst = param(46)        ! for lee's restart (rarely used)
c     if (irst.gt.0)  call setup_convect(2)


c     If moving geometry then add a perturbation to the
c     mesh coordinates (see Subroutine INIGEOM)

c     IF (IFMVBD) CALL PTBGEOM      

C
C     Find out what type of i.c. is requested
C     Current options: 
C
C     (1) - User specified fortran function (default is zero i.c.)
C     (2) - Restart from file(s)
C     (3) - Activate pre-solver => steady diffusion / steady Stokes
C
C     If option (2) is requested, also return with the name of the
C     restart file(s) together with the associated dump number
C
c     call slogic (iffort,ifrest,ifprsl,nfiles)
C
C     Set up proper initial values for turbulence model arrays
C

c     IF (IFMODEL) CALL PRETMIC

C
C      ***** TEMPERATURE AND PASSIVE SCALARS ******
C
C     Check if any pre-solv necessary for temperature/passive scalars
C
      IFANYP = .FALSE.
      DO 100 IFLD=2,NFIELD
         IF (IFPRSL(IFLD,jp)) THEN
            IF (NID.EQ.0) WRITE(6,101) IFLD
            IFANYP = .TRUE.
         ENDIF
  100 CONTINUE
  101 FORMAT(2X,'Using PRESOLVE option for field',I2,'.')
C
C
C     If any pre-solv, do pre-solv for all temperatur/passive scalar fields
C
c     IF (IFANYP) CALL PRSOLVT

C
C     Fortran function initial conditions for temp/pass. scalars.
C
c     call nekgsync()
c     call restart(nfiles)
c     call nekgsync()
C
C
C      ***** VELOCITY ******
C
C     (If restarting for V, we're done,
C     ...else, do pre-solv for fluid if requested.)
C
c     IFIELD = 1
c     IF (IFPRSL(IFLD,jp)) CALL PRSOLVV
C
C
C     Fortran function initial conditions for velocity.
C
c     ifield = 1
c     if (iffort(ifield,jp)) then
c        if (nid.eq.0) write(6,*) 'call nekuic for vel  '
c        call nekuic
c     endif
c
c     if (ifpert) then
c        ifield=1
c        do jp=1,npert
c           if (iffort(ifield,jp)) call nekuic
c           if (nid.eq.0) write(6,*) 'ic vel pert:',iffort(1,jp),jp
c        enddo
c     endif
c     jp = 0
c
c     ntotv = nx1*ny1*nz1*nelv
c
c     Fortran function initial conditions for turbulence k-e model
c
c     if (ifmodel .and. ifkeps) then
c        mfldt = nfield - 1
c        do 300 ifield=mfldt,nfield
c           if (iffort(ifield,jp)) call nekuic
c300     continue
c     endif

C....          
C
C     Initial mesh velocities
C
c     IF (IFMVBD) CALL OPCOPY (WX,WY,WZ,VX,VY,VZ)
C      IF (IFMVBD.AND..NOT.IFREST(0,jp)) CALL MESHV (2)

C
C....      
C     Compute additional initial values for turbulence model arrays
C     based on I.C.
C
C     IF (IFMODEL) CALL POSTMIC
C....      

C
C     If convection-diffusion of a passive scalar with a fixed velocity field,
C     make sure to fill up lagged arrays since this will not be done in
C     the time-stepping procedure (no flow calculation) (01/18/91 -EMR).
C
c     IF (.NOT.IFFLOW.AND.IFHEAT) THEN
c        ITEST=0
c        DO 400 IFLD=2,NFIELD
c           IF (IFADVC(IFLD)) ITEST=1
c400     CONTINUE
c        IF (ITEST.EQ.1) THEN
c           NBDMAX = 3
c           NBDSAV = NBDINP
c           NBDINP = NBDMAX
c           DO 500 I=1,NBDMAX
c...           CALL LAGVEL
c500        CONTINUE
c           NBDINP = NBDSAV
c        ENDIF
c     ENDIF
c     
c     Ensure that all processors have the same time as node 0.
c
c     IF (NID.NE.0) TIME=0.0
c     TIME=GLSUM(TIME,1)
c     NTDUMP=0
c     IF (TIMEIO.NE.0.0) NTDUMP = INT( TIME/TIMEIO )
c
c     Ensure that initial field is continuous!
c
      nxyz1=nx1*ny1*nz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1
      nn = nxyz1
      ntotg=nelgv*nn
c
      ifield = 2
      if (ifflow) ifield = 1
      call rone(work,ntotv)
      ifield = 1
      CALL DSSUM(work,NX1,NY1,NZ1)
      CALL COL2(work,VMULT,NTOTV)
      rdif = glsum(work,ntotv)
      rtotg = ntotg
      rdif = (rdif-rtotg)/rtotg
      if (abs(rdif).gt.1e-14) then
         if (nid.eq.0) write(*,*) 'ERROR: dssum test has failed!',rdif
         call exitt
      endif

      vxmax = glamax(vx,ntotv)
      vymax = glamax(vy,ntotv)
      vzmax = glamax(vz,ntotv)
      prmax = glamax(pr,ntot2)

      ntot = nxyz1*nelfld(2)
      ttmax = glamax(t ,ntot)

      do i=1,NPSCAL
         ntot = nx1*ny1*nz1*nelfld(i+2)
         psmax(i) = glamax(T(1,1,1,1,i+1),ntot)
      enddo

c 
c     small=1.0E-20
c     ifldsave = ifield
c     if (vxmax.eq.0.0) call perturb(vx,1,small)
c     if (vymax.eq.0.0) call perturb(vy,1,small)
c     if (vzmax.eq.0.0) call perturb(vz,1,small)
c     if (prmax.eq.0.0.and.ifsplit) call perturb(pr,1,small)
c     if (ttmax.eq.0.0) call perturb(t ,2,small)
c
C print min values
      xxmax = glmin(xm1,ntott)
      yymax = glmin(ym1,ntott)
      zzmax = glmin(zm1,ntott)

      vxmax = glmin(vx,ntotv)
      vymax = glmin(vy,ntotv)
      vzmax = glmin(vz,ntotv)
      prmax = glmin(pr,ntot2)

      ntot = nxyz1*nelfld(2)
      ttmax = glmin(t ,ntott)

      do i=1,LDIMT-1
         ntot = nxyz1*nelfld(i+2)
         psmax(i) = glmin(T(1,1,1,1,i+1),ntot)
      enddo

      if (nid.eq.0) then
         write(6,19) xxmax,yymax,zzmax
   19    format(' xyz min  ',5g13.5)
      endif
      if (nid.eq.0) then
         write(6,20) vxmax,vymax,vzmax,prmax,ttmax
   20    format(' uvwpt min',5g13.5)
      endif
      if (LDIMT-1.gt.0) then
         if (nid.eq.0) write(6,21) (psmax(i),i=1,LDIMT-1)
   21    format(' PS min   ',50g13.5)
      endif

c print max values
      xxmax = glmax(xm1,ntott)
      yymax = glmax(ym1,ntott)
      zzmax = glmax(zm1,ntott)

      vxmax = glmax(vx,ntotv)
      vymax = glmax(vy,ntotv)
      vzmax = glmax(vz,ntotv)
      prmax = glmax(pr,ntot2)

      ntot = nxyz1*nelfld(2)
      ttmax = glmax(t ,ntott)

      if (nid.eq.0) then
         write(6,16) xxmax,yymax,zzmax
   16    format(' xyz max  ',5g13.5)
      endif

      if (nid.eq.0) then
         write(6,17) vxmax,vymax,vzmax,prmax,ttmax
   17    format(' uvwpt max',5g13.5)
      endif

c     call setup_convect(2)

      if(nid.eq.0) then
        write(6,*) 'done :: set initial conditions'
        write(6,*) ' '
      endif

      return
      end
C            
c-----------------------------------------------------------------------
      function i1_from_char(s1)
      character*1 s1

      character*10 n10
      save         n10
      data         n10 / '0123456789' /

      i1_from_char = indx2(n10,10,s1,1)-1

      return
      end
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX_CUT(S1,S2,L2)
C
C     INDX_CUT is returned with the location of S2 in S1 (0 if not found)
C     S1     is returned with 1st occurance of S2 removed.
C
      CHARACTER*1 S1(132),S2(132)
C
      I1=INDX1(S1,S2,L2)
C
      IF (I1.NE.0) THEN
C
         N1=132-L2
         DO 100 I=I1,N1
            I2=I+L2
C           remove the 1st occurance of S2 from S1.
            S1(I)=S1(I2)
  100    CONTINUE
         N2=N1+1
         DO 200 I=N2,132
            S1(I)=' '
  200    CONTINUE
      ENDIF
C
      INDX_CUT=I1
      return
      END
c-----------------------------------------------------------------------
      subroutine chknorm (ifzero)
C----------------------------------------------------------------------
C
C     Check if trivial user specified initial conditions
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      LOGICAL IFZERO
C
      IFZERO = .TRUE.
C
      IF (IFFLOW) THEN
c        IFIELD = 1
         IMESH  = 1
C        CALL UNORM
         IF (VNRML8.GT.0.) IFZERO = .FALSE.
      ENDIF
      IF (IFHEAT) THEN
         DO 100 IFLD=2,NFIELD
            IMESH = 1
            IF (IFTMSH(IFLD)) IMESH = 2
C           CALL UNORM
            IF (TNRML8(IFLD).GT.0.) IFZERO = .FALSE.
 100     CONTINUE
      ENDIF
c
      return
      END
C
c-----------------------------------------------------------------------
      subroutine prsolvt
C----------------------------------------------------------------------
C
C     Use steady state solution as initial condition 
C     for temperatur/passive scalar
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      LOGICAL  IFSAV1,IFSAV2(LDIMT1)
C
      IF (NID.EQ.0) WRITE(6,*) ' '
      IF (NID.EQ.0) WRITE(6,*) 'Conduction pre-solver activated'
C
C     Set logical IFTRAN to false (steady state)
C     Save logicals for convection
C     Turn convection off
C
      NFIELD = 1
      IFSAV1 = IFTRAN
      IFTRAN = .FALSE.
      DO 100 IFLD=2,NFIELD
         IFSAV2(IFLD) = IFADVC(IFLD)
         IFADVC(IFLD) = .FALSE.
 100  CONTINUE
C

c...  CALL SETPROP
c...  CALL SETSOLV

C
      IF(NID.EQ.0)WRITE(6,*)'Steady conduction/passive scalar problem'
C
      DO 200 IGEOM=1,2
C        CALL HEAT (IGEOM)
 200  CONTINUE

C
C     Set IFTRAN to true again
C     Turn convection on again
C
      IFTRAN = IFSAV1
      DO 300 IFLD=2,NFIELD
         IFADVC(IFLD) = IFSAV2(IFLD)
 300  CONTINUE
C
      return
      END
C
c-----------------------------------------------------------------------
      LOGICAL FUNCTION IFGTRL(VALUE,LINE)
C
C     Read VALUE from LINE and set IFGTRL to .TRUE. if successful,
C                                  IFGTRL to .FALSE. otherwise.
C
C     This complicated function is necessary thanks to the Ardent,
C     which won't allow free formatted reads (*) from internal strings!
C
      CHARACTER*132 LINE
      CHARACTER*132 WORK
      CHARACTER*8  FMAT
C
C     Note that the format Fn.0 is appropriate for fields of type:
C          34   34.0  34.0e+00
C     The only difficulty would be with '34' but since we identify
C     the field width exactly, there is no problem.
C
      IFGTRL=.FALSE.
      VALUE=0.0
C
      WORK=LINE
      CALL LJUST(WORK)
      IFLDW=INDX1(WORK,' ',1)-1
C
      IF (IFLDW.GT.0) THEN
         WRITE(FMAT,10) IFLDW
   10    FORMAT('(F',I3.3,'.0)')
         READ(WORK,FMAT,ERR=100,END=100) TVAL
         VALUE=TVAL
         IFGTRL=.TRUE.
         return
      ENDIF
C
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine perturb(tt,ifld,eps)
      include 'SIZE'
      include 'TOTAL'
c
      real tt(1)
      integer ifld

      ifield = ifld

      n = nx1*ny1*nz1*nelfld(ifield)
      call vcospf(tt,bm1,n)
      call cmult(tt,eps,n)
      call dssum(tt,nx1,ny1,nz1)

      return
      end
c-----------------------------------------------------------------------
      subroutine vcospf(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = cos(1000.*y(i))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine vbyte_swap(x,n)
      character*1 x(0:3,1),tmp0,tmp1
      character*1 in (0:3), out(0:3)
      real*4      in4     , out4
      equivalence (in ,in4 )
      equivalence (out,out4)
c
      do i=1,n
         do j=0,3
            in (j) = x(j,i)
         enddo
         tmp0   = x(0,i)
         tmp1   = x(1,i)
         x(0,i) = x(3,i)
         x(1,i) = x(2,i)
         x(2,i) = tmp1
         x(3,i) = tmp0
         do j=0,3
            out(j) = x(j,i)
         enddo
         write(6,*) 'swap:',i,in4,out4
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)
      include 'SIZE'
c
      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern
c
      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.
c
      test2 = bytetest
      call byte_reverse(test2,1,ierr)
      if (nid.eq.0) 
     $   write(6,*) 'byte swap:',if_byte_swap_test,bytetest,test2
      return
      end
c-----------------------------------------------------------------------
