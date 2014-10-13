c-----------------------------------------------------------------------
      subroutine outmesh
      include 'SIZE'
      include 'TOTAL'
      integer e,eg

      common /cmesh/ xt(2**ldim,ldim)

      len = wdsize*ndim*(2**ndim)

      if (nid.eq.0) open(unit=29,file='rea.new')

      do eg=1,nelgt
         mtype = eg
         call nekgsync()          !  belt
         jnid = gllnid(eg)
         e    = gllel (eg)
         if (jnid.eq.0 .and. nid.eq.0) then
            call get_el(xt,xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e))
            call out_el(xt,eg)
         elseif (nid.eq.0) then
            call crecv(mtype,xt,len)
            call out_el(xt,eg)
         elseif (jnid.eq.nid) then
            call get_el(xt,xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e))
            call csend(mtype,xt,len,0,0)
         endif
         call nekgsync()          !  suspenders
      enddo

      if (nid.eq.0) close(29)
      call nekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine out_el(xt,e)
      include 'SIZE'
      include 'TOTAL'

      real xt(2**ldim,ldim)
      integer e

      integer ed(8)
      save    ed
      data    ed  / 1,2,4,3 , 5,6,8,7 /

      write(29,1) e
      write(29,2) ((xt(ed(k),j),k=1,4),j=1,ndim)
      write(29,2) ((xt(ed(k),j),k=5,8),j=1,ndim)

    1 format(12x,'ELEMENT',i6,' [    1 ]    GROUP     0')
    2 format(1p4e18.10)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_el(xt,x,y,z)
      include 'SIZE'
      include 'TOTAL'

      real xt(2**ldim,ldim)
      real x(nx1,ny1,nz1),y(nx1,ny1,nz1),z(nx1,ny1,nz1)

      l = 0
      do k=1,nz1,nz1-1
      do j=1,ny1,ny1-1
      do i=1,nx1,nx1-1
         l = l+1
         xt(l,1) = x(i,j,k)
         xt(l,2) = y(i,j,k)
         xt(l,3) = z(i,j,k)
      enddo
      enddo
      enddo

      return
      end
