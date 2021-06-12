c----------------------------------------------------------------------
      subroutine setupcomm()
      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'
      include 'INPUT'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /happycallflag/ icall
      integer nid_global_root(0:nsessmax-1)
      character*132 session_mult(0:nsessmax-1),path_mult(0:nsessmax-1)

      logical ifhigh

      ! Init MPI
      call mpi_initialized(mpi_is_initialized, ierr)

      if (mpi_is_initialized .eq. 0 ) call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,np_global,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,nid_global,ierr)


      ! check upper tag size limit
      call mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,nval,flag,ierr)

      if (nval.lt.(10000+lp)) then
         if(nid_global.eq.0) write(6,*) 'ABORT: MPI_TAG_UB too small!'
         call exitt(1)
      endif

      ! set defaults
      nid         = nid_global
      nekcomm     = MPI_COMM_WORLD
      iglobalcomm = MPI_COMM_WORLD

      ifneknek    = .false.
      nsessions   = 1
      session_id  = 1

      ierr = 0
      nlin = 0

      if (nid.eq.0) then
         write(6,*) 'Reading session file ...'
         open (unit=8,file='SESSION.NAME',status='old',err=24)
 21      read (8,*,END=22)
         nlin = nlin + 1
         goto 21
 22      rewind(8)

         if (nlin.gt.2) read(8,*,err=24) nsessions

         call blank(session_mult(0),132)
         call blank(path_mult(0)   ,132)
         read(8,11,err=24) session_mult(0)
         read(8,11,err=24) path_mult(0)
         if (nsessions.gt.1) read(8,*,err=24) npsess(0)    
         ! fixed with npesss(0): each session using same # cores

 11      format(a132)
         close(8)
         write(6,*) 'Number of sessions   :',nsessions
         goto 23
 24      ierr = 1
      endif
 23   continue

      call err_chk(ierr,' Error while reading SESSION.NAME!$')

      call bcast(nsessions,ISIZE)
      call bcast(session_mult(0),132*CSIZE)
      call bcast(path_mult(0),132*CSIZE)
      call bcast(npsess(0),ISIZE) 

      if (nsessions .gt. 1) ifneknek = .true.

      if (nsessions .gt. nsessmax) then
          call exitti('nsessmax in SIZE too low!$',nsessmax,1)
      endif

      ! single session run
      if (.not.ifneknek) then
         ifneknekc = .false.
         session   = session_mult(0)
         path      = path_mult(0)
         call mpi_comm_dup(mpi_comm_world,iglobalcomm,ierr)
         intracomm = iglobalcomm
         return
      endif

c     Check if specified number of ranks in each session is consistent
c     with the total number of ranks
      npall=0
      do n=0,nsessions-1
         npall=npall+npsess(0) !fixed
      enddo

      if (npall.ne.np_global) then
         call exitti('Number of ranks does not match!$',npall,1)
      endif

c     Assign key for splitting into multiple groups
      nid_global_root_next=0
      do n=0,nsessions-1
         nid_global_root(n)=nid_global_root_next
         nid_global_root_next=nid_global_root(n)+npsess(0) !fixed
         if (nid_global.ge.nid_global_root(n).and.
     $       nid_global.lt.nid_global_root_next) idsess = n
      enddo

      call mpi_comm_split(mpi_comm_world,idsess,nid,intracomm,ierr)

      session_id = idsess+1
      session = session_mult(0)      
      path    = path_mult   (0)     

      return
      end
c----------------------------------------------------------------------
      subroutine setupcomm_neknek()
      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'
      include 'INPUT'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /happycallflag/ icall
      integer nid_global_root(0:nsessmax-1)
      character*132 session_mult(0:nsessmax-1),path_mult(0:nsessmax-1)

      logical ifhigh

      ! Init MPI
      call mpi_initialized(mpi_is_initialized, ierr)

      if (mpi_is_initialized .eq. 0 ) call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD,np_global,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,nid_global,ierr)


      ! check upper tag size limit
      call mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,nval,flag,ierr)

      if (nval.lt.(10000+lp)) then
         if(nid_global.eq.0) write(6,*) 'ABORT: MPI_TAG_UB too small!'
         call exitt(1)
      endif

      ! set defaults
      nid         = nid_global
      nekcomm     = MPI_COMM_WORLD
      iglobalcomm = MPI_COMM_WORLD

      ifneknek    = .false.
      ifneknekc   = .false. ! session are uncoupled
      ifneknekm   = .false. ! not moving
      nsessions   = 1

      ierr = 0
      nlin = 0

      if (nid.eq.0) then
         write(6,*) 'Reading session file ...'
         open (unit=8,file='SESSION.NAME',status='old',err=24)
 21      read (8,*,END=22)
         nlin = nlin + 1
         goto 21
 22      rewind(8)
         if (nlin.gt.2) read(8,*,err=24) nsessions

         if (nsessions.gt.1) read(8,*,err=24) ifneknekc

         do n=0,nsessions-1
            call blank(session_mult(n),132)
            call blank(path_mult(n)   ,132)
            read(8,11,err=24) session_mult(n)
            read(8,11,err=24) path_mult(n)
            if (nsessions.gt.1) read(8,*,err=24)  npsess(n)
c           write(6,*) 'n,npsess(n)',n,npsess(n)
         enddo

 11      format(a132)
         close(8)
         write(6,*) 'Number of sessions:',nsessions
         goto 23
 24      ierr = 1
      endif
 23   continue
      call err_chk(ierr,' Error while reading SESSION.NAME!$')

      call bcast(nsessions,ISIZE)

      do n = 0,nsessions-1
         call bcast(npsess(n),ISIZE)
         call bcast(session_mult(n),132*CSIZE)
         call bcast(path_mult(n),132*CSIZE)
      enddo

      if (nsessions .gt. 1) ifneknek = .true.

      if (nsessions .gt. nsessmax) then
          call exitti('nsessmax in SIZE too low!$',nsessmax,1)
      endif

      ! single session run
      if (.not.ifneknek) then
         ifneknekc = .false.
         session   = session_mult(0)
         path      = path_mult(0)
         call mpi_comm_dup(mpi_comm_world,iglobalcomm,ierr)
         intracomm = iglobalcomm
         return
      endif

c     Check if specified number of ranks in each session is consistent
c     with the total number of ranks
      npall=0
      do n=0,nsessions-1
         npall=npall+npsess(n)
      enddo

      if (npall.ne.np_global) then
         call exitti('Number of ranks does not match!$',npall,1)
      endif

c     Assign key for splitting into multiple groups
      nid_global_root_next=0
      do n=0,nsessions-1
         nid_global_root_next=nid_global_root(n)+npsess(n)
         nid_global_root_next=nid_global_root(n)+npsess(0) !fixed
         if (nid_global.ge.nid_global_root(n).and.
     $       nid_global.lt.nid_global_root_next) idsess = n
      enddo

      call mpi_comm_split(mpi_comm_world,idsess,nid,intracomm,ierr)

      session_id = idsess+1
      session = session_mult(idsess)
      path    = path_mult   (idsess)

c...  nid is global

      return
      end
c----------------------------------------------------------------------
      subroutine iniproc
      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'mpif.h'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      logical flag
      integer*8 isize_mpi, lb

      nid  = mynode()
      nid_ = nid
      np   = numnodes()
      np_  = np

c     nio = -1             ! Default io flag


c...  nid is set in each session
c     if (nid.eq.0) then
c        if (ifneknek) then
c          call set_stdout(' ',idsess)
c        else
c          call set_stdout(' ',-1)
c        endif
c     endif
 
      if (nid_io.eq.0) call printHeader


c.... check 
c     if (nid.eq.0) nio=0  ! Only node 0 writes

      if (np.gt.lp) then 
      if (nid_io.eq.0) then
         write(6,*)
     $   'ERROR: Code compiled for a max of',LP,' processors.'
         write(6,*)
     $   'Recompile with LP =',NP,' or run with fewer processors.'
         write(6,*)
     $   'Aborting in routine INIPROC.'
         call exitt(1)
      endif
      endif

      if (wdsize .eq. 4) then
      if (nid_io.eq.0) then
         write(6,*) 'ABORT: single precision mode not supported!'
         call exitt(1)
      endif
      endif
      nekreal = mpi_real
      if (wdsize.eq.8) nekreal = mpi_double_precision

      call MPI_Type_Get_Extent(MPI_DOUBLE_PRECISION,lb,isize_mpi,ierr)
      if (isize_mpi .ne. wdsize) then
         call exitti('MPI real size does not match$',isize_mpi)
      endif

      call MPI_Type_Get_Extent(MPI_INTEGER,lb,isize_mpi,ierr)
       if (isize_mpi .ne. isize) then
         call exitti('MPI integer size does not match$',isize_mpi)
      endif

      call MPI_Type_Get_Extent(MPI_INTEGER8,lb,isize_mpi,ierr)
       if (isize_mpi .ne. isize8) then
         call exitti('MPI integer8 size does not match$',isize_mpi)
      endif



c
      PID = 0
      NULLPID=0
      NODE0=0
      NODE= NID+1

C     Test timer accuracy
      edif = 0.0
      do i = 1,10
         e1 = dclock()
         e2 = dclock()
         edif = edif + e2-e1
      enddo
      edif = edif/10.

      if (nid_io.eq.0) then
         write(6,*) 'Number of processors:',np
         WRITE(6,*) 'REAL    wdsize      :',WDSIZE
         WRITE(6,*) 'INTEGER wdsize      :',ISIZE
         WRITE(6,'(A,1pE8.2)') ' Timer accuracy      : ',edif
         WRITE(6,*) ' '
      endif

      call crystal_setup(cr_h,nekcomm,np)  ! set cr handle to new instance

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

      return
      end
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
         call exitt(1)
      endif
c
      call copy(x,w,n)
c
      tgop =tgop +(dclock()-etime1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine igop(x,w,op,n)
c-----------------------------------------------------------------------
c     Global vector commutative operation using spanning tree.
      implicit none
      include 'mpif.h'

      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer nid,np,nekcomm,nekgroup,nekreal

      integer x(n),w(n)
      integer n
      character*3 op

      integer ierr

      if (op.eq.'+  ') then
        call mpi_allreduce(x,w,n,mpi_integer,mpi_sum,nekcomm,ierr)
      elseif (op.eq.'M  ') then
        call mpi_allreduce(x,w,n,mpi_integer,mpi_max,nekcomm,ierr)
      elseif (op.eq.'m  ') then
        call mpi_allreduce(x,w,n,mpi_integer,mpi_min,nekcomm,ierr)
      elseif (op.eq.'*  ') then
        call mpi_allreduce(x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
      else
         write(6,*) nid,' OP ',op,' not supported. ABORT in GOP.'
        call exitt(1)
      endif

      call icopy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine i8gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      integer*8 x(n), w(n)
      character*3 op

      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer8,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
      endif

      call i8copy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine csend(mtype,buf,len,jnid,jpid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real   buf(1)
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
      real   buf(1)
      len = lenm
      jnid = mpi_any_source

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
      real   etime,q(2)
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
      SUBROUTINE LBCAST(IFIF)
C
C  Broadcast logical variable to all processors.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'mpif.h'

      LOGICAL IFIF
C
      IF (NP.EQ.1) RETURN
C
      ITEM=0
      IF (IFIF) ITEM=1
      call bcast(item,isize)
      IFIF=.FALSE.
      IF (ITEM.EQ.1) IFIF=.TRUE.
      RETURN
      END

c-----------------------------------------------------------------------
c     subroutine bcast(buf,len)
c-----------------------------------------------------------------------
c     include 'mpif.h'
c     common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c     real*4 buf(1)
c
c     call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)
c
c     return
c     end
c-----------------------------------------------------------------------
      subroutine bcast(buf,len)
c-----------------------------------------------------------------------
      include 'mpif.h'
 
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4   buf(len)
 
      call mpi_bcast(buf,len,mpi_byte,0,nekcomm,ierr)
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
      call mpi_group_free (nekgroup,ierr)
      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      integer x(len)
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
      integer x(len)
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
      subroutine close_unit(io)
      close (unit=io)
      return
      end
c-----------------------------------------------------------------------
      subroutine exitti(stringi,idata,irc)
      character*1 stringi(132)
      character*1 stringo(132)
      character*11 s11
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      call blank(stringo,132)
      call chcopy(stringo,stringi,132)
      len = indx1(stringo,'$',1)
      write(s11,11) idata
   11 format(1x,i10)
      call chcopy(stringo(len),s11,11)

      if (nid.eq.0) write(6,1) (stringo(k),k=1,len+10)
    1 format('EXIT: ',132a1)

      call exitt(irc)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt(rc)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'mpif.h'

      integer rc

      integer ierr,tmp(1),irc
      integer iglmax
      real tstop,dnekclock
      real papi_flops,papi_mflops,dtmp1,dtmp2,dgp

c     Make the return code the maximum of the return codes
      tmp(1) = rc
      irc = iglmax(tmp,1)

c     Communicate unhappiness to the other session
      call nekgsync

      tstop = dnekclock()
      ttotal = tstop-etimes
      nxyz = nx1*ny1*nz1

      if (nid_io.eq.0) then
         dtmp1 = 0
         dtmp2 = 0
         if(istep.gt.0) then
            dgp = nvtot
            dgp = max(dgp,1.0)
            dtmp1 = np*ttime/(dgp*max(istep,1))
            dtmp2 = ttime/max(istep,1)
         endif
         write(6,*) ' '
         write(6,'(A)') 'call exitt: dying ...'
         write(6,*) ' '
      endif

      call flush_io
      call mpi_finalize(ierr)
      call exit(irc)

      return
      end
c-----------------------------------------------------------------------
      subroutine printHeader

      INCLUDE 'HEADER'

      return
      end
c-----------------------------------------------------------------------
      integer function igl_running_sum(in)
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
c-----------------------------------------------------------------------
      subroutine nekgsync()

      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_barrier(nekcomm,ierr)

      return
      end

c-----------------------------------------------------------------------
      real*8 function dnekclock()
      include 'mpif.h'
c
      dnekclock=mpi_wtime()
c
      return
      end
c-----------------------------------------------------------------------
      real*8 function dnekclock_sync()
      include 'mpif.h'
c
      call nekgsync()
      dnekclock_sync=mpi_wtime()
c
      return
      end

c-----------------------------------------------------------------------
      subroutine nek_flops(flops,mflops)
      real*4 rtime,ptime,mflops
      integer*8 flops

c     call getflops_papi(flops,mflops)

      return
      end

c-----------------------------------------------------------------------
c      subroutine getflops_papi(flpops,mflops)
c#ifdef PAPI
c      include 'f77papi.h'
c      real*4 rtime,ptime,mflops
c      integer*8 flpops
c
cc     call papif_flops(rtime,ptime,flpops,mflops,ierr)
c      if(ierr.ne.0) then
c        flpops = -1
c        mflops = -1
c      endif
c#else
c      flpops = -1
c      mflops = -0
c#endif
c
c      return
c      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine crecv2(mtype,buf,lenm,jnid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm

      call mpi_recv (buf,len,mpi_byte
     $              ,jnid,mtype,nekcomm,status,ierr)
c
      if (len.gt.lenm) then
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine crecv3(mtype,buf,len,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call mpi_recv (buf,len,mpi_byte
     $            ,jnid,mtype,nekcomm,status,ierr)
      call mpi_get_count (status,mpi_byte,len,ierr)
c
      if (len.gt.lenm) then
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
