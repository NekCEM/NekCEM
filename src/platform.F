c-----------------------------------------------------------------------
      subroutine platform_timer(ivb) ! mxm, ping-pong, and all_reduce timer

      include 'SIZE'
      include 'TOTAL'


      call mxm_test_all(nid,ivb)  ! measure mxm times
c     call exitti('done mxm_test_all$',ivb)

      call comm_test(ivb)         ! measure message-passing and all-reduce times

      return
      end
c-----------------------------------------------------------------------
      subroutine comm_test(ivb) ! measure message-passing and all-reduce times
                                ! ivb = 0 --> minimal verbosity
                                ! ivb = 1 --> fully verbose
                                ! ivb = 2 --> smaller sample set(shorter)

      include 'SIZE'
      include 'PARALLEL'

      call gop_test(ivb)   ! added, Jan. 8, 2008

      log_np=log2(np)
      np2 = 2**log_np
      if (np2.eq.np) call gp2_test(ivb)   ! added, Jan. 8, 2008

      io = 6
      n512 = min(512,np-1)

c     do nodeb=1,n512
c        call pingpongo(alphas,betas,0,nodeb,.0005,io,ivb)
c        if (nid.eq.0) write(6,2) nodeb,np,alphas,betas
c   2    format(2i10,1p2e15.7,' alpha betao')
c     enddo

      do kk=0,0
      do nodeb=1,n512
         call pingpong (alphas,betas,0,nodeb,.0005,io,ivb,kk)
         if (nid.eq.0) write(6,1) nodeb,np,alphas,betas,kk
    1    format(2i10,1p2e15.7,' alpha beta',i1)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pingpong(alphas,betas,nodea,nodeb,dt,io,ivb,kk)

      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd = 3*lt)
      common /scrns/ x(mwd),y(mwd)

      include 'mpif.h'
      integer status(mpi_status_size)

      character*10 fname

      if (nid.eq.nodea) then
         write(fname,3) np,nodeb
    3    format('t',i4.4,'.',i4.4)
         if (io.ne.6) open (unit=io,file=fname)
      endif

      call nekgsync
      call get_msg_vol(msg_vol,dt,nodea,nodeb) ! Est. msg vol for dt s

      nwds = 0
      if (nid.eq.nodea.and.ivb.gt.0) write(io,*)

      betas = 0  ! Reported inverse bandwidth
      count = 0

      do itest = 1,500

         nloop = msg_vol/(nwds+2)
         nloop = min(nloop,1000)
         nloop = max(nloop,4)

         len   = 8*nwds
     
         if (kk.eq.0)
     $      call ping_loop (t1,t0,len,nloop,nodea,nodeb,nid,x,y,x,y)
         if (kk.eq.1)
     $      call ping_loop1(t1,t0,len,nloop,nodea,nodeb,nid,x,y,x,y)
         if (kk.eq.2)
     $      call ping_loop2(t1,t0,len,nloop,nodea,nodeb,nid,x,y,x,y)

         if (nid.eq.nodea) then
            tmsg = (t1-t0)/(2*nloop)   ! 2*nloop--> Double Buffer
            tmsg = tmsg / 2.           ! one-way cost = 1/2 round-trip
            tpwd = tmsg                ! time-per-word
            if (nwds.gt.0) tpwd = tmsg/nwds
            if (ivb.gt.0) write(io,1) nodeb,np,nloop,nwds,tmsg,tpwd,kk
    1       format(3i6,i12,1p2e16.8,' pgn',i1)

            if (nwds.eq.1) then
               alphas = tmsg
            elseif (nwds.gt.10000) then   ! "average" beta
               betas = (betas*count + tpwd)/(count+1)
               count = count + 1
            endif
         endif

         if (ivb.eq.2) then
            nwds = (nwds+1)*1.25
         else
            nwds = (nwds+1)*1.016
         endif
         if (nwds.gt.mwd) then
c        if (nwds.gt.1024) then
            if (nid.eq.nodea.and.io.ne.6) close(unit=io)
            call nekgsync
            return
         endif

      enddo

      if (nid.eq.nodea.and.io.ne.6) close(unit=io)
      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine pingpongo(alphas,betas,nodea,nodeb,dt,io,ivb)

      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd = 3*lt)
      common /scrns/ x(mwd),y(mwd)

      include 'mpif.h'
      integer status(mpi_status_size)

      character*10 fname

      if (nid.eq.nodea) then
         write(fname,3) np,nodeb
    3    format('t',i4.4,'.',i4.4)
         if (io.ne.6) open (unit=io,file=fname)
      endif

      call nekgsync
      call get_msg_vol(msg_vol,dt,nodea,nodeb) ! Est. msg vol for dt s

      nwds = 0
      if (nid.eq.nodea.and.ivb.gt.0) write(io,*)

      betas = 0  ! Reported inverse bandwidth
      count = 0

      do itest = 1,500
         call nekgsync
         nloop = msg_vol/(nwds+2)
         nloop = min(nloop,1000)
         nloop = max(nloop,1)

         len   = 8*nwds
         jnid = mpi_any_source

         if (nid.eq.nodea) then

            msg  = irecv(itest,y,1)
            call csend(itest,x,1,nodeb,0)   ! Initiate send, to synch.
            call msgwait(msg)

            t0 = mpi_wtime ()
            do i=1,nloop
               call mpi_irecv(y,len,mpi_byte,mpi_any_source,i
     $                        ,nekcomm,msg,ierr)
               call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)
               call mpi_wait (msg,status,ierr)
            enddo
            t1 = mpi_wtime ()
            tmsg = (t1-t0)/nloop
            tmsg = tmsg / 2.       ! Round-trip message time = twice one-way
            tpwd = tmsg
            if (nwds.gt.0) tpwd = tmsg/nwds
            if (ivb.gt.0) write(io,1) nodeb,np,nloop,nwds,tmsg,tpwd
    1       format(3i6,i12,1p2e16.8,' pgo')

            if (nwds.eq.1) then
               alphas = tmsg
            elseif (nwds.gt.10000) then
               betas = (betas*count + tpwd)/(count+1)
               count = count + 1
            endif

         elseif (nid.eq.nodeb) then

            call crecv(itest,y,1)           ! Initiate send, to synch.
            call csend(itest,x,1,nodea,0)

            t0 = dnekclock()
            do i=1,nloop
               call mpi_recv (y,len,mpi_byte
     $               ,jnid,i,nekcomm,status,ierr)
               call mpi_send (x,len,mpi_byte,nodea,i,nekcomm,ierr)
            enddo
            t1 = dnekclock()
            tmsg = (t1-t0)/nloop

         endif

         nwds = (nwds+1)*1.016
         if (nwds.gt.mwd) then
            if (nid.eq.nodea.and.io.ne.6) close(unit=io)
            call nekgsync
            return
         endif

      enddo

      if (nid.eq.nodea.and.io.ne.6) close(unit=io)
      call nekgsync

      return
      end
c-----------------------------------------------------------------------
      subroutine get_msg_vol(msg_vol,dt,nodea,nodeb)
      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrns/ x(3*lt),y(3*lt)
!
!     Est. msg vol for dt s
!
      msg_vol = 1000

      nwds  = min(1000,lt)
      nloop = 50
 
      tmsg = 0.
      call gop(tmsg,t1,'+  ',1)

      len = 8*nwds
      if (nid.eq.nodea) then

         msg  = irecv(1,y,1)
         call csend(1,x,1,nodeb,0)   ! Initiate send, to synch.
         call msgwait(msg)

         t0 = dnekclock()
         do i=1,nloop
            msg  = irecv(i,y,len)
            call csend(i,x,len,nodeb,0)
            call msgwait(msg)
         enddo
         t1   = dnekclock()
         tmsg = (t1-t0)/nloop
         tpwd = tmsg/nwds

      elseif (nid.eq.nodeb) then

         call crecv(1,y,1)           ! Initiate send, to synch.
         call csend(1,x,1,nodea,0)

         t0 = dnekclock()
         do i=1,nloop
            call crecv(i,y,len)
            call csend(i,x,len,nodea,0)
         enddo
         t1   = dnekclock()
         tmsg = (t1-t0)/nloop
         tmsg = 0.

      endif

      call gop(tmsg,t1,'+  ',1)
      msg_vol = nwds*(dt/tmsg)
c     if (nid.eq.nodea) write(6,*) nid,msg_vol,nwds,dt,tmsg,' msgvol'

      return
      end
c-----------------------------------------------------------------------
      subroutine gop_test(ivb)
      include 'SIZE'
      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      include 'mpif.h'
      integer status(mpi_status_size)

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd = 3*lt)
      common /scrns/ x(mwd),y(mwd)
      common /scruz/ times(2,500)
      common /scrcg/ nwd(500)

      nwds  = 1
      mtest = 0
      do itest = 1,500
         nwds = (nwds+1)*1.016
         if (nwds.gt.mwd) goto 100
         mtest = mtest+1
         nwd(mtest) = nwds
      enddo
  100 continue

      nwds = 1
      do itest = mtest,1,-1

         tiny = 1.e-27
         call cfill(x,tiny,mwd)
         nwds = nwd(itest)
         call nekgsync

         t0 = mpi_wtime ()
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         call gop(x,y,'+  ',nwds)
         t1 = mpi_wtime ()

         tmsg = (t1-t0)/6 ! six calls
         tpwd = tmsg
         if (nwds.gt.0) tpwd = tmsg/nwds
         times(1,itest) = tmsg
         times(2,itest) = tpwd

      enddo
  101 continue


      if (nid.eq.0) then
         nwds = 1
         do itest=1,500
            if (ivb.gt.0.or.itest.eq.1) 
     $         write(6,1) np,nwds,(times(k,itest),k=1,2)
    1       format(i12,i12,1p2e16.8,' gop')
            nwds = (nwds+1)*1.016
            if (nwds.gt.mwd) goto 102
         enddo
  102    continue
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gp2_test(ivb)

      include 'SIZE'
      include 'mpif.h'

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)

      parameter  (lt=lx1*ly1*lz1*lelt)
      parameter (mwd = 3*lt)
      common /scrns/ x(mwd),y(mwd)
      common /scruz/ times(2,500)

      call rzero(x,mwd)

      nwds = 1
      do itest = 1,500
         call gp2(x,y,'+  ',1,nid,np)

         t0 = mpi_wtime ()
         call gp2(x,y,'+  ',nwds,nid,np)
         call gp2(x,y,'+  ',nwds,nid,np)
         call gp2(x,y,'+  ',nwds,nid,np)
         call gp2(x,y,'+  ',nwds,nid,np)
         t1 = mpi_wtime ()

         tmsg = (t1-t0)/4 ! four calls
         tpwd = tmsg
         if (nwds.gt.0) tpwd = tmsg/nwds
         times(1,itest) = tmsg
         times(2,itest) = tpwd

         nwds = (nwds+1)*1.016
         if (nwds.gt.mwd) goto 101
      enddo
  101 continue


      if (nid.eq.0) then
         nwds = 1
         do itest=1,500
            if (ivb.gt.0.or.itest.eq.1) 
     $         write(6,1) np,nwds,(times(k,itest),k=1,2)
    1       format(i12,i12,1p2e16.8,' gp2')
            nwds = (nwds+1)*1.016
            if (nwds.gt.mwd) goto 102
         enddo
  102    continue
      endif

      return
      end
c-----------------------------------------------------------------------
      integer function xor(m,n)
c
c  If NOT running on a parallel processor, it is sufficient to
c  have this routine return a value of XOR=1.
c
c  Pick one of the following:
c
c  UNIX 4.2, f77:
       XOR = OR(M,N)-AND(M,N)
c
c  Intel FTN286:
c     XOR = M.NEQV.N
c
c  Ryan-McFarland Fortran
C      XOR = IEOR(M,N)
c
c     XOR = 0
c     IF(M.EQ.1 .OR.  N.EQ.1) XOR=1
c     IF(M.EQ.0 .AND. N.EQ.0) XOR=0
c     IF(M.EQ.1 .AND. N.EQ.1) XOR=0
c     IF(M.GT.1 .OR.N.GT.1 .OR.M.LT.0.OR.N.LT.0) THEN
c        PRINT*,'ERROR IN XOR'
c        STOP
c     ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine gp2( x, w, op, n, nid, np)
c
c     Global vector commutative operation using spanning tree.
c
c     Std. fan-in/fan-out

      real x(n), w(n)
      character*3 op

      integer bit, bytes, cnt, diff, spsize, i, 
     *   parent, troot, xor, root, lnp, log2
      logical ifgot

      integer type
      save    type
      data    type  /998/

      type  = type+100
      if (type.gt.9992) type=type-998
      typer = type-1
      bytes = 8*n

      root    = 0
      troot   = max0((nid/np)*np, root)
      diff    = xor(nid,troot)
      nullpid = 0

c     Accumulate contributions from children, if any
      level2=1
    5 continue
         level=level2
         level2=level+level
         if (mod(nid,level2).ne.0) goto 20
            call crecv(type,w,bytes)
            if (op.eq.'+  ') then
               do i=1,n
                  x(i) = x(i) + w(i)
               enddo
            elseif (op.eq.'*  ') then
               do i=1,n
                  x(i) = x(i) * w(i)
               enddo
            elseif (op.eq.'M  ') then
               do i=1,n
                  x(i) = max(x(i),w(i))
               enddo
            elseif (op.eq.'m  ') then
               do i=1,n
                  x(i) = min(x(i),w(i))
               enddo
            endif
         if (level2.lt.np) goto 5

c     Pass result back to parent
   20 parent = nid-level
      if (nid .ne. 0) call csend(type,x,bytes,parent,nullpid)

c     Await final answer from node 0 via log_2 fan out
      level=np/2
      ifgot=.false.
      if (nid.eq.root) ifgot=.true.

      lnp = log2(np)
      do i=1,lnp
        if (ifgot) then
           jnid=nid+level
           call csend(typer,x,bytes,jnid,nullpid)
        elseif (mod(nid,level).eq.0) then
           call crecv(typer,x,bytes)
           ifgot=.true.
        endif
        level=level/2
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ping_loop1(t1,t0,len,nloop,nodea,nodeb,nid,x,y)

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real x(1),y(1)

      include 'mpif.h'
      integer status(mpi_status_size)

      i=0
      if (nid.eq.nodea) then
         call nekgsync
         call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr)    ! 1b
         call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
c        call mpi_rsend(x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
         call msgwait(msg)                                          ! 1b

         t0 = mpi_wtime ()
         do i=1,nloop
            call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr) ! 2b
            call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
c           call mpi_rsend(x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
            call mpi_wait (msg,status,ierr)                         ! 2b
         enddo
         t1 = mpi_wtime ()

      elseif (nid.eq.nodeb) then

         call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr)    ! 1a
         call nekgsync
         call mpi_wait (msg,status,ierr)                            ! 1a

         j=i
         do i=1,nloop
            call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr) ! 2a
c           call mpi_rsend(x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_wait (msg,status,ierr)                         ! 2a
            j=i
         enddo
c        call mpi_rsend(x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb
         call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb

      else
         call nekgsync
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ping_loop2(t1,t0,len,nloop,nodea,nodeb,nid,x,y)

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real x(1),y(1)

      include 'mpif.h'
      integer status(mpi_status_size)

      i=0
      if (nid.eq.nodea) then
         call nekgsync
         call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr)    ! 1b
         call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)        ! 1a
         call msgwait(msg)                                          ! 1b

         t0 = mpi_wtime ()
         do i=1,nloop
            call mpi_send (x,len,mpi_byte,nodeb,i,nekcomm,ierr)     ! 2a
            call mpi_irecv(y,len,mpi_byte,nodeb,i,nekcomm,msg,ierr) ! 2b
            call mpi_wait (msg,status,ierr)                         ! 2b
         enddo
         t1 = mpi_wtime ()

      elseif (nid.eq.nodeb) then

         call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr)    ! 1a
         call nekgsync
         call mpi_wait (msg,status,ierr)                            ! 1a

         j=i
         do i=1,nloop
            call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)     ! 1b
            call mpi_irecv(y,len,mpi_byte,nodea,i,nekcomm,msg,ierr) ! 2a
            call mpi_wait (msg,status,ierr)                         ! 2a
            j=i
         enddo
         call mpi_send (x,len,mpi_byte,nodea,j,nekcomm,ierr)        ! nb

      else
         call nekgsync
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ping_loop(t1,t0,len,nloop,nodea,nodeb,nid,x1,y1,x2,y2)
c     Double Buffer : does 2*nloop timings

      common /nekmpi/ mid,np,nekcomm,nekgroup,nekreal

      real x1(1),y1(1),x2(1),y2(1)

      include 'mpif.h'
      integer status(mpi_status_size)

      itag=1
      if (nid.eq.nodea) then
         call mpi_irecv(y1,len,mpi_byte,nodeb,itag,nekcomm,msg1,ierr)   ! 1b 
         call nekgsync


         t0 = mpi_wtime ()
         do i=1,nloop
            call mpi_send (x1,len,mpi_byte,nodeb,itag,nekcomm,ierr)     ! 1a 
            call mpi_irecv(y2,len,mpi_byte,nodeb,itag,nekcomm,msg2,ierr)! 2b 
            call mpi_wait (msg1,status,ierr)                            ! 1b
            call mpi_send (x2,len,mpi_byte,nodeb,itag,nekcomm,ierr)     ! 2a 
            call mpi_irecv(y1,len,mpi_byte,nodeb,itag,nekcomm,msg1,ierr)! 3b 
            call mpi_wait (msg2,status,ierr)                            ! 2b
         enddo
         t1 = mpi_wtime ()
         call mpi_send (x1,len,mpi_byte,nodeb,itag,nekcomm,ierr)        ! nb
         call mpi_wait (msg1,status,ierr)                              ! nb

      elseif (nid.eq.nodeb) then

         call mpi_irecv(y1,len,mpi_byte,nodea,itag,nekcomm,msg1,ierr)   ! nb 
         call nekgsync


         do i=1,nloop
            call mpi_wait (msg1,status,ierr)                            ! 1a
            call mpi_send (x1,len,mpi_byte,nodea,itag,nekcomm,ierr)     ! 1b
            call mpi_irecv(y2,len,mpi_byte,nodea,itag,nekcomm,msg2,ierr)! 2a
            call mpi_wait (msg2,status,ierr)                            ! 2a 
            call mpi_send (x2,len,mpi_byte,nodea,itag,nekcomm,ierr)     ! 2b
            call mpi_irecv(y1,len,mpi_byte,nodea,itag,nekcomm,msg1,ierr)! 3a
         enddo
         call mpi_wait (msg1,status,ierr)                            ! 2a 
         call mpi_send (x1,len,mpi_byte,nodea,itag,nekcomm,ierr)        ! nb

      else
         call nekgsync
      endif

      return
      end
c-----------------------------------------------------------------------
      integer*8 function i8gl_running_sum(in)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      integer*8 x,r

      x = in  ! running sum
      r = 0   ! recv buff

      call mpi_scan(x,r,1,mpi_integer8,mpi_sum,nekcomm,ierr)
      i8gl_running_sum = r

      return
      end
c-----------------------------------------------------------------------
