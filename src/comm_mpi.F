      subroutine iniproc2
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      logical flag

c...  this is a dummy subroutine: to be deleted misun 07/10

      return
      end
c-----------------------------------------------------------------------
      subroutine create_comm(icomm)
      include 'mpif.h'

c     mpi_comm_world: communicator (handle)
c     itmp: subset of the group of  mpi_comm_world (handle)
c     icomm: new communicator (handle)
      call mpi_comm_group (mpi_comm_world,itmp,ierr)
      call mpi_comm_create (mpi_comm_world,itmp,icomm,ierr)
      call mpi_group_free (itmp,ierr)
c     write(6,*) 'nekcomm:',nekcomm

      return
      end
c-----------------------------------------------------------------------
      subroutine iniproc
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      include 'mpif.h'
C
      COMMON /CTMP0/ IPG(LP),IPG1(LP),ITMP(LP)
      INTEGER CD1
C
C
C     Initialize MPI
C
      call mpi_init      (ierr)
      call init_nek_comm (nid,np,wdsize)
      flag_gs_init = 0
c
      if (nid.eq.0) call printHeader

      CD  = LOG2(NP)
      PID = 0
      NODE= NID+1
C
C     Check cube dimension 
C
      IF (NP.GT.LP) THEN
         WRITE(6,*) 
     $   'ERROR: Code compiled for a max of',LP,' processors.'
         WRITE(6,*) 
     $   'Recompile with LP =',NP,' or run with fewer processors.'
         WRITE(6,*) 
     $   'Aborting in routine INIPROC.'
         call exitt
      ENDIF
C
c     These flags added for native 64-bit machines  (pff 10/1/98)
c
      ifdblas = .false.
      if (wdsize.eq.8) ifdblas = .true.
c
c     This is to determine the message length of integers
c
      isize = 4
c
c     On the Cray:   ifdblas = .false., since sblas = 64 bit
c                    isize  is 8
c
c     ifdblas = .false.
c     isize   = 8
c
c
      MANAGER=0
      ALLNODES=-1
      NULLPID=0
      NODE0=0
C
C
C     Initialize ring pass data
C
      CALL GRAY(IPG,IPG1,ITMP,NP)
      IGNODE=IPG1(NODE)-1
      DO 100 IP=1,NP
         IGNODE=IGNODE+1
         IF (IGNODE.GT.NP) IGNODE=IGNODE-NP
         IPRING(IP)=IPG(IGNODE)-1
  100 CONTINUE
      LFTNBR=IPRING(NP)
      RGTNBR=IPRING(2)
      if (nid.eq.0) then
         write(6,*) 'Number of processors:',np
         WRITE(6,*) 'REAL    wdsize      :',WDSIZE
         WRITE(6,*) 'INTEGER wdsize      :',ISIZE
      endif
C
C     All done
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine init_nek_comm2
      include 'mpif.h'
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
C
      call create_comm(nekcomm) ! set up nekton specific communicator
c
      nid_  = mynode()
      np_   = numnodes()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine init_nek_comm(nido,npo,wdsize)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer wdsize
C
      call nek_comm            ! set up nekton specific communicator
c
      nid  = mynode()
      np   = numnodes()
      nido = nid
      npo  = np
c
      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) wdsize=8
      nekreal = mpi_real
      if (wdsize.eq.8) nekreal = mpi_double_precision
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gop( x, w, op, n)
c
c     Global vector commutative operation using spanning tree.
c
      INCLUDE 'CTIMER'
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
      real x(n), w(n)
      character*3 op
c
      if (icalld.eq.0) then
        tgop =0.0d0
        ngop =0
        icalld=1
      endif
      ngop = ngop + 1
      etime1=dclock()
c
      if (op.eq.'+  ') then
c        call mpi_allreduce_(x,w,n,nekreal,mpi_sum ,nekcomm,ierr)
         call mpi_allreduce (x,w,n,nekreal,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
c        call mpi_allreduce_(x,w,n,nekreal,mpi_max ,nekcomm,ierr)
         call mpi_allreduce (x,w,n,nekreal,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
c        call mpi_allreduce_(x,w,n,nekreal,mpi_min ,nekcomm,ierr)
         call mpi_allreduce (x,w,n,nekreal,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
c        call mpi_allreduce_(x,w,n,nekreal,mpi_prod,nekcomm,ierr)
         call mpi_allreduce (x,w,n,nekreal,mpi_prod,nekcomm,ierr)
      else
         write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
         call exitt
      endif
c
      call copy(x,w,n)
c
      tgop =tgop +(dclock()-etime1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine igop( x, w, op, n)
c
c     Global vector commutative operation using spanning tree.
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
      integer x(n), w(n)
      character*3 op
c
      if (op.eq.'+  ') then
c       call mpi_allreduce_(x,w,n,mpi_integer,mpi_sum ,nekcomm,ierr)
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
c       call mpi_allreduce_(x,w,n,mpi_integer,mpi_max ,nekcomm,ierr)
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
c       call mpi_allreduce_(x,w,n,mpi_integer,mpi_min ,nekcomm,ierr)
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
c       call mpi_allreduce_(x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
        call exitt
      endif
c
      call icopy(x,w,n)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine csend(mtype,buf,len,jnid,jpid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)
c     call mpi_send_(buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)
      call mpi_send (buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)
      return
      END
c-----------------------------------------------------------------------
      subroutine crecv(   mtype, buf ,lenm )
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

c     call mpi_recv_(buf,len,mpi_byte
c    $            ,jnid,mtype,nekcomm,status,ierr)
      call mpi_recv (buf,len,mpi_byte
     $            ,jnid,mtype,nekcomm,status,ierr)
c
      if (len.gt.lenm) 
     $    write(6,*) nid,'long message in mpi_crecv:',len,lenm
c
      return
      END
c-----------------------------------------------------------------------
      integer function numnodes()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c     call mpi_comm_size_(nekcomm, numnodes , ierr)
      call mpi_comm_size (nekcomm, numnodes , ierr)
      return
      end
c-----------------------------------------------------------------------
      integer function mynode()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer myid
c     call mpi_comm_rank_(nekcomm, myid, ierr)
      call mpi_comm_rank (nekcomm, myid, ierr)
      mynode = myid
      return
      end
c-----------------------------------------------------------------------
      real*8 function dclock()
      include 'mpif.h'
c
      real*4 etime,q(2)
      save q
      data q /0.,0./
c
c     dclock=mpi_wtime_()
      dclock=mpi_wtime ()
c     dclock=etime(q)    ! for alpha
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lbcast(ifif)
C
C     Broadcast logical variable to all processors.
C
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'
C
      logical ifif
C
      if (np.eq.1) return
C
      item=0
      if (ifif) item=1
      call bcast(item,isize)
      ifif=.false.
      if (item.eq.1) ifif=.true.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine bcast(buf,len)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)
c     call mpi_bcast_(buf,len,mpi_byte,0,nekcomm,ierr)
      call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine nek_comm
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c     call mpi_comm_group_(mpi_comm_world,nekgroup,ierr)
c     call mpi_comm_create_(mpi_comm_world,nekgroup,nekcomm,ierr)
c     call mpi_group_free_(nekgroup,ierr)
      call mpi_comm_group (mpi_comm_world,nekgroup,ierr)
      call mpi_comm_create (mpi_comm_world,nekgroup,nekcomm,ierr)
       write(6,*) 'nekcomm',nekcomm,mpi_comm_world,nekgroup,ierr
        call exitt
      call mpi_group_free (nekgroup,ierr)
      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
c     call mpi_isend_(x,len,mpi_byte,jnid,msgtag
c    $       ,nekcomm,imsg,ierr)
      call mpi_isend (x,len,mpi_byte,jnid,msgtag
     $       ,nekcomm,imsg,ierr)
      isend = imsg
c     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)
c
      return
      end
c-----------------------------------------------------------------------
      function irecv(msgtag,x,len)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
c     call mpi_irecv_(x,len,mpi_byte,mpi_any_source,msgtag
c    $       ,nekcomm,imsg,ierr)
      call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,nekcomm,imsg,ierr)
      irecv = imsg
c     write(6,*) nid,' irecv:',imsg,msgtag,len
c
      return
      end
c-----------------------------------------------------------------------
      subroutine msgwait(imsg)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
c
c     write(6,*) nid,' msgwait:',imsg
c
c     call mpi_wait_(imsg,status,ierr)
      call mpi_wait (imsg,status,ierr)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gsync   
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer x(1),w(1)
c
      x(1) = 0
c     call mpi_allreduce_(x,w,n,mpi_integer,mpi_sum,nekcomm,ierr)
      call mpi_allreduce (x,w,n,mpi_integer,mpi_sum,nekcomm,ierr)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rring(work,x,n,irg)
      DIMENSION X(N),WORK(N)
C
C     Pass data X around a ring network, and receive corresponding
C     data into WORK.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      include 'mpif.h'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /scrns/ wrk(lx1*ly1*lz1*lelt)
      common /ctmp1/ istat(lx1*lelt)
C
      IF (IRG.EQ.1) THEN
         CALL COPY(WORK,X,N)
      ELSE
         LEN   = WDSIZE*N
         ITYPE = NID
         JTYPE = RGTNBR
         CALL copy (wrk,work,n)
c
c        call mpi_irecv_(work,2*len,mpi_byte,rgtnbr,jtype
         call mpi_irecv (work,2*len,mpi_byte,rgtnbr,jtype
     $        ,nekcomm,msg,ierr)
c
c        call mpi_send_(wrk ,len,mpi_byte,lftnbr,itype
         call mpi_send (wrk ,len,mpi_byte,lftnbr,itype
     $       ,nekcomm,ierr)
c        call mpi_wait_(msg,istat,ierr)
         call mpi_wait (msg,istat,ierr)
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine cring(work,x,n,irg)
C
C     Pass data X around a ring network, and receive corresponding
C     data into WORK.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      include 'mpif.h'
      DIMENSION X(N),WORK(N)
      CHARACTER*1 X,WORK
C
      PARAMETER (LTOT2=LX1*LY1*LZ1*LELT*2)
      PARAMETER (LTOT8=LTOT2*4)
      COMMON /CTMP0/ CW1(LTOT8)
      CHARACTER*1    CW1
      DIMENSION      RW1(LTOT2)
      EQUIVALENCE   (RW1,CW1)
c
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /scrns/ wrk(lx1*ly1*lz1*lelt)
      common /ctmp1/ istat(lx1*lelt)
C
      IF (IRG.EQ.1) THEN
         CALL chcopy(WORK,X,N)
      ELSE
C
C        enough work space?
         IF (N.GT.LTOT8) THEN
            WRITE(6,100) NID,N,LTOT8
  100       FORMAT(2X,I5,
     $     'WARNING: In routine CRING, not enough work space.'
     $     ,/,2X,'Required # of words:',I7,' Supplied:',I7)
            CALL EXITT
         ENDIF
C
c
         CALL chcopy(CW1,WORK,N)
         LEN   = N
         ITYPE = NID
         JTYPE = RGTNBR
c
c        call mpi_irecv_(work,2*len,mpi_byte,rgtnbr,jtype
         call mpi_irecv (work,2*len,mpi_byte,rgtnbr,jtype
     $       ,nekcomm,msg   ,ierr)
c        call mpi_send_(rw1 ,len,mpi_byte,lftnbr,itype
         call mpi_send (rw1 ,len,mpi_byte,lftnbr,itype
     $       ,nekcomm,ierr)
c        call mpi_wait_(msg,istat,ierr)
         call mpi_wait (msg,istat,ierr)
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine iring(work,x,n,irg)
      INTEGER X(N),WORK(N)
C
C     Pass data X around a ring network, and receive corresponding
C     data into WORK.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      include 'mpif.h'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /scrns/ wrk(lx1*ly1*lz1*lelt)
      common /ctmp1/ istat(lx1*lelt)
C
      IF (IRG.EQ.1) THEN
         CALL ICOPY(WORK,X,N)
      ELSE
         LEN   = isize*N
         ITYPE = NID
         JTYPE = RGTNBR
         CALL icopy (wrk,work,n)
c
c        call mpi_irecv_(work,2*len,mpi_byte,rgtnbr,jtype
         call mpi_irecv (work,2*len,mpi_byte,rgtnbr,jtype
     $       ,nekcomm,msg ,ierr)
c        call mpi_send_(wrk ,len,mpi_byte,lftnbr,itype
         call mpi_send (wrk ,len,mpi_byte,lftnbr,itype
     $       ,nekcomm, ierr)
c        call mpi_wait_(msg,istat,ierr)
         call mpi_wait (msg,istat,ierr)
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine close_unit(io)
      close (unit=io)
      return
      end
c-----------------------------------------------------------------------
      subroutine exitt
c     common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'
c
      x = 1.0
      x = glsum(x,1)
c
      write(6,*) nid,' normal exit.'
      call flush_io
c
c     if (nid.eq.0) then              ! For the NERSC IBM
c        do k=16,24
c           call close_unit(k)
c        enddo
c     endif
c
c
c     call mpi_finalize_(ierr)
      call mpi_finalize (ierr)
      call exit
c
      return
      end
c-----------------------------------------------------------------------
      subroutine printHeader

      INCLUDE 'HEADER'

      return
      end
c-----------------------------------------------------------------------
      function igl_running_sum(in)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      integer x,w,r

      x = in  ! running sum
      w = in  ! working buff
      r = 0   ! recv buff

      call mpi_scan(x,r,1,mpi_integer,mpi_sum,nekcomm,ierr)
      igl_running_sum = r

      return
      end

