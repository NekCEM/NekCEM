c-----------------------------------------------------------------------
      subroutine setup_topo
c
c     Parallel compatible routine to find 
c     connectivity of element structure.
c
c     On Processor 0:
c
c     .Verify right-handedness of elements.
c     .Verify element-to-element reciprocity of BC's
c     .Verify correlation between E-E BC's and physical coincidence
c     .Set rotations
c     .Determine multiplicity
c     .Set up direct stiffness summation arrays.
c
c     All Processors:
c
c     .Disperse/Receive BC and MULT temporary data read from
c     preprocessor.
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'NONCON'
      include 'ZPER'
      include 'SCRCT'
c
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
c
      common /c_is1/ glo_num(1*lx1*ly1*lz1*lelv)
      integer*8      glo_num
      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer        vertex

      if (nid.eq.0) write(6,*) 'setup mesh topology'
c
c     Initialize key arrays for Direct Stiffness SUM.
c
      NXL=3  
      NYL=3   
      NZL=1+2*(NDIM-2)

      call initds
      call dsset (nx1,ny1,nz1)
      call setedge
c
c=================================================
c     Establish (global) domain topology
c=================================================
c
c     .Generate topologically correct mesh data.
c     .Set up element centers, face centers, etc. 
c     .Check  right handedness of elements.
c     .Check  element boundary conditions.
c     .Establish Element-Element rotations
c     .Construct the element to processor map and

      call genxyzl
      call setside
      call verify

      ncrnr = 2**ndim

      if (nelgv.eq.nelgt) then

         if (ifgtp) then
            call gen_gtp_vertex    (vertex, ncrnr)
         else
            call get_vert
         endif
         call setup_ds
     $    (gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
           gsh_fld(2)=gsh_fld(1)

      else
c
c       For conjugate heat transfer, it is assumed that fluid
c       elements are listed both globally and locally with lower
c       element numbers than the solid elements.
c
c       We currently assume that there is at least one fluid elem.
c       per processor.

        call get_vert
c       call outmati(vertex,4,nelv,'vrtx V')
        call setup_ds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)

c       call get_vert  (vertex, ncrnr, nelgt, '.mp2')  !  LATER !
c       call outmati(vertex,4,nelt,'vrtx T')

        call setup_ds(gsh_fld(2),nx1,ny1,nz1,nelt,nelgt,vertex,glo_num)

c       check if there is a least one fluid element on each processor
        do iel = 1,nelt
           ieg = lglel(iel)
           if (ieg.le.nelgv) goto 101
        enddo

        if(nid.eq.0) write(6,*)
     &    'ERROR: each domain must contain at least one element!'
        call exitt
 101  continue

      endif

#ifdef _OPENACC
      call setup_dgds2_acc
     $        (gsh_fld(3),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
#else
      call setup_dgds2
     $        (gsh_fld(3),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
#endif

      if(nid.eq.0) then
        write(6,*) 'done :: setup mesh topology'
        write(6,*) '.'
      endif


      ntotv = nx1*ny1*nz1*nelv
      ntott = nx1*ny1*nz1*nelt
c
c     if (ifflow) then
c        ifield = 1
         call rone    (vmult,ntotv)
         call dssum   (vmult,nx1,ny1,nz1)
         vmltmax = glmax(vmult,ntotv)
         ivmltmax= vmltmax
         if (nid.eq.0) write(6,*) ivmltmax,' max multiplicity'

         call invcol1 (vmult,ntotv)
c     endif


c     if (ifheat) then
c        ifield = 2
         call rone    (tmult,ntott)
         call dssum   (tmult,nx1,ny1,nz1)
         call invcol1 (tmult,ntott)
c     endif


      return
      end
c-----------------------------------------------------------------------
      subroutine setup_ds(gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'

      integer gs_handle
      integer vertex(1)
      integer*8 glo_num(1),ngv

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      t0 = dclock()

c     Global-to-local mapping for gs
      call set_vert2(glo_num,ngv,nx,nel,vertex,.false.)

c     Initialize gather-scatter code
      if     (IFHEX) then
        if (if3d) then
          nxyz  = nx*ny*nz
          nxzf  = nx*nz
          nfaces= ndim*2
        else
          nxyz  = nx*ny
          nxzf  = nx*nz
          nfaces= ndim*2
        endif
      elseif (IFTET) then
        if (if3d) then
          nxyz  = nx*(nx+1)*(nx+2)/6
          nxzf  = nx*(nx+1)/2
          nfaces= ndim+1
        else
          nxyz  = nx*(nx+1)/2
          nxzf  = nx
          nfaces= ndim+1
        endif
      endif

      npts = nxyz*nelt
      nxzfl= nxzf*nfaces*nelt

      call gs_setup(gs_handle,glo_num,npts,nekcomm,np)
      !write(6,*) 'gs_handle',gs_handle

c     call gs_chkr(glo_num)

      t1 = dclock() - t0
      if (nid.eq.0) then
         write(6,1) t1,gs_handle,nx,ngv,melg
    1    format('   setup_ds time',1pe11.4,' seconds ',2i3,2i12)
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine set_vert2(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     Given global array, vertex, pointing to hex vertices, set up
c     a new array of global pointers for an nx^ndim set of elements.
c
      include 'SIZE'
      include 'INPUT'
c
      integer*8 glo_num(1),ngv
      integer vertex(1),nx
      logical ifcenter

      if (if3d) then
         call setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)
      else
         call setvert2d(glo_num,ngv,nx,nel,vertex,ifcenter)
      endif

      if(nid.eq.0) write(6,*) '  call usrsetvert'
      call usrsetvert(glo_num,nel,nx,nx,nx)
      if(nid.eq.0) write(6,'(A,/)') '   done :: usrsetvert'

      return
      end
c-----------------------------------------------------------------------
      subroutine gbtuple_rank(tuple,m,n,nmax,cr_h,nid,np,ind)
c
c     Return a unique rank for each matched tuple set. Global.  Balanced.
c
c     tuple is destroyed.
c
c     By "balanced" we mean that none of the tuple entries is likely to
c     be much more uniquely populated than any other, so that any of
c     the tuples can serve as an initial (parallel) sort key
c
c     First two slots in tuple(:,i) assumed empty
c
      integer ind(nmax),tuple(m,nmax),cr_h

      parameter (mmax=40)
      integer key(mmax),wtuple(mmax)

      if (m.gt.mmax) then
         write(6,*) nid,m,mmax,' gbtuple_rank fail'
         call exitt
      endif

      do i=1,n
         tuple(1,i) = mod(tuple(3,i),np) ! destination processor
         tuple(2,i) = i                  ! return location
      enddo

      ni= n
      ky= 1  !Assumes crystal_new already called
      
C     do i=1,n
C        write(6,*) 'tup:',
C    $   nid,tuple(1,i),tuple(2,i),tuple(3,i),tuple(4,1)
C     enddo
      
      call crystal_ituple_transfer(cr_h, tuple,m,ni,nmax, ky)

      nimx = iglmax(ni,1)
      if (ni.gt.nmax)   write(6,*) ni,nmax,n,'cr_xfer problem, A'
      if (nimx.gt.nmax) call exitt

      nkey = m-2
      do k=1,nkey
         key(k) = k+2
      enddo

      call irank_vecn(ind,nu,tuple,m,ni,key,nkey,wtuple)! tuple re-ordered,
                                                        ! but contents same

      nu_tot   = igl_running_sum(nu) ! running sum over P processors
      nu_prior = nu_tot - nu

      do i=1,ni
         tuple(3,i) = ind(i) + nu_prior  ! global ranking
      enddo

      call crystal_ituple_transfer(cr_h, tuple,m,ni,nmax, ky)

      nk = 1  ! restore to original order, local rank: 2; global: 3
      ky = 2
      call ituple_sort(tuple,m,n,ky,nk,ind,wtuple)


      return
      end
c-----------------------------------------------------------------------
      subroutine setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     setup unique ids for dssum  
c     note:
c     total number of unique vertices, edges and faces has to be smaller 
c     than 2**31 (integer-4 limit).
c     if nelgt < 2**31/12 we're ok for sure (independent of N)! 
c
      include 'SIZE'
      include 'CTIMER'
      include 'PARALLEL'
      include 'TOPOL'
      include 'GEOM'

      integer*8 glo_num(1),ngv
      integer   vertex(0:1,0:1,0:1,1),nx
      logical   ifcenter

      integer  edge(0:1,0:1,0:1,3,lelt),enum(12,lelt),fnum(6,lelt)
      common  /scrmg/ edge,enum,fnum

      integer  etuple(4,5*12*lelt),ftuple(5,6,5*lelt)
      integer  ind(5*12*lelt)
      common  /scrns/ ind,etuple
      equivalence (etuple,ftuple)

      integer  gvf(4),facet(4),aa(3),key(3),e
      logical  ifij
      
      integer*8 igv,ig0
      integer*8 ngvv,ngve,ngvs,ngvi,ngvm
      integer*8 n_on_edge,n_on_face,n_in_interior
      integer*8 i8glmax
c
      ny   = nx
      nz   = nx
      nxyz = nx*ny*nz
c
      key(1)=1
      key(2)=2
      key(3)=3
c
c     Assign hypercube ordering of vertices
c     -------------------------------------
c
c     Count number of unique vertices
      nlv  = 2**ndim
      ngvv = iglmax(vertex,nlv*nel)
c
      do e=1,nel
         do k=0,1
         do j=0,1
         do i=0,1
c           Local to global node number (vertex)
            il  = 1 + (nx-1)*i + nx*(nx-1)*j + nx*nx*(nx-1)*k
            ile = il + nx*ny*nz*(e-1)
            glo_num(ile)   = vertex(i,j,k,e)
         enddo
         enddo
         enddo
      enddo
      ngv  = ngvv
c
      if (nx.eq.2) return
c
c     Assign global vertex numbers to SEM nodes on each edge
c     ------------------------------------------------------
c
c     Assign edge labels by bounding vertices.  
      do e=1,nel
         do k=0,1
         do j=0,1
         do i=0,1
            edge(i,j,k,1,e) = vertex(i,j,k,e)  ! r-edge
            edge(j,i,k,2,e) = vertex(i,j,k,e)  ! s-edge
            edge(k,i,j,3,e) = vertex(i,j,k,e)  ! t-edge
         enddo
         enddo
         enddo
      enddo
c
c     Sort edges by bounding vertices.
      do i=0,12*nel-1
         if (edge(0,i,0,1,1).gt.edge(1,i,0,1,1)) then
            kswap = edge(0,i,0,1,1)
            edge(0,i,0,1,1) = edge(1,i,0,1,1)
            edge(1,i,0,1,1) = kswap
         endif
         etuple(3,i+1) = edge(0,i,0,1,1)
         etuple(4,i+1) = edge(1,i,0,1,1)
      enddo
c
c     Assign a number (rank) to each unique edge
      m    = 4
      n    = 12*nel
      nmax = 12*lelt*5  ! 5x for crystal router factor of safety

      call gbtuple_rank(etuple,m,n,nmax,cr_h,nid,np,ind)
      do i=1,12*nel
         enum(i,1) = etuple(3,i)
      enddo
      n_unique_edges = iglmax(enum,12*nel)
c
      n_on_edge = nx-2
      ngve      = n_unique_edges*n_on_edge
      do e=1,nel
         iedg_loc = 0
c
c        Edges 1-4
         do k=0,1
         do j=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = nx*(nx-1)*j + nx*nx*(nx-1)*k
            i0e = i0 + nxyz*(e-1)
            if (glo_num(i0e+1).lt.glo_num(i0e+nx)) then
               do i=2,nx-1                                   ! std forward case
                  glo_num(i0e+i) = igv + i-1
               enddo
            else
               do i=2,nx-1                                   ! backward case
                  glo_num(i0e+i) = igv + 1 + n_on_edge-(i-1)
               enddo
            endif
            iedg_loc = iedg_loc + 1
         enddo
         enddo
c
c        Edges 5-8
         do k=0,1
         do i=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = 1+(nx-1)*i + nx*nx*(nx-1)*k
            i0e = i0 + nxyz*(e-1)
            if (glo_num(i0e).lt.glo_num(i0e+nx*(nx-1))) then
               do j=2,nx-1                                   ! std forward case
                  glo_num(i0e+(j-1)*nx) = igv + j-1
               enddo
            else
               do j=2,nx-1                                   ! backward case
                  glo_num(i0e+(j-1)*nx) = igv + 1 + n_on_edge-(j-1)
               enddo
            endif
            iedg_loc = iedg_loc + 1
         enddo
         enddo
c
c        Edges 9-12
         do j=0,1
         do i=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = 1 + (nx-1)*i + nx*(nx-1)*j
            i0e = i0 + nxyz*(e-1)
            if (glo_num(i0e).lt.glo_num(i0e+nx*nx*(nx-1))) then
               do k=2,nx-1                                   ! std forward case
                  glo_num(i0e+(k-1)*nx*nx) = igv + k-1
               enddo
            else
               do k=2,nx-1                                   ! backward case
                  glo_num(i0e+(k-1)*nx*nx) = igv + 1 + n_on_edge-(k-1)
               enddo
            endif
            iedg_loc = iedg_loc + 1
         enddo
         enddo
      enddo
      ngv   = ngv + ngve
c
c     Asign global node numbers on the interior of each face
c     ------------------------------------------------------ 
c
c     Assign faces by 3-tuples 
c
c     (The following variables all take the symmetric 
c     notation of IFACE as arguments:)
c
c     ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
c                         as depicted below, e.g. ICFACE(i,2)=2,4,6,8.
c
c                        3+-----+4    ^ Y
c                        /  2  /|     |
c     Edge 1 extends    /     / |     |
c       from vertex   7+-----+8 +2    +----> X
c       1 to 2.        |  4  | /     /
c                      |     |/     /
c                     5+-----+6    Z
c                         3
c
      nfaces=ndim*2
      ncrnr =2**(ndim-1)
      do e=1,nel
         do ifac=1,nfaces
            do icrn=1,ncrnr
               i                  = icface(icrn,ifac)-1
               facet(icrn)        = vertex(i,0,0,e)
            enddo
            call isort(facet,ind,ncrnr)
            call icopy(ftuple(3,ifac,e),facet,ncrnr-1)
         enddo
      enddo

c     Assign a number (rank) to each unique face
      m    = 5
      n    = 6*nel
      nmax = 6*lelt*5  ! 5x for crystal router factor of safety
      if (nid.eq.0) write(6,*) nid,nel,nmax,n,' Fgbtuple',nx
      call gbtuple_rank(ftuple,m,n,nmax,cr_h,nid,np,ind)
      do i=1,6*nel
         fnum(i,1) = ftuple(3,i,1)
      enddo
      n_unique_faces = iglmax(fnum,6*nel)
c
      call dsset (nx,ny,nz)
      do e=1,nel
       do iface=1,nfaces
         i0 = skpdat(1,iface)
         i1 = skpdat(2,iface)
         is = skpdat(3,iface)
         j0 = skpdat(4,iface)
         j1 = skpdat(5,iface)
         js = skpdat(6,iface)
c
c        On each face, count from minimum global vertex number,
c        towards smallest adjacent vertex number.  e.g., suppose
c        the face is defined by the following global vertex numbers:
c
c
c                    11+--------+81
c                      |c      d|
c                      |        |
c                      |        |
c                      |a      b|
c                    15+--------+62
c                          
c        We would count from c-->a, then towards d.
c
         gvf(1) = glo_num(i0+nx*(j0-1)+nxyz*(e-1))
         gvf(2) = glo_num(i1+nx*(j0-1)+nxyz*(e-1))
         gvf(3) = glo_num(i0+nx*(j1-1)+nxyz*(e-1))
         gvf(4) = glo_num(i1+nx*(j1-1)+nxyz*(e-1))
c
         call irank(gvf,ind,4)
c
c        ind(1) tells which element of gvf() is smallest.
c
         ifij = .false.
         if (ind(1).eq.1) then
            idir =  1
            jdir =  1
            if (gvf(2).lt.gvf(3)) ifij = .true.
         elseif (ind(1).eq.2) then
            idir = -1
            jdir =  1
            if (gvf(1).lt.gvf(4)) ifij = .true.
         elseif (ind(1).eq.3) then
            idir =  1
            jdir = -1
            if (gvf(4).lt.gvf(1)) ifij = .true.
         elseif (ind(1).eq.4) then
            idir = -1
            jdir = -1
            if (gvf(3).lt.gvf(2)) ifij = .true.
         endif
c
         if (idir.lt.0) then
            it=i0
            i0=i1
            i1=it
            is=-is
         endif
c
         if (jdir.lt.0) then
            jt=j0
            j0=j1
            j1=jt
            js=-js
         endif
c
         nxx = nx*nx
         n_on_face = (nx-2)*(ny-2)
         ngvs  = n_unique_faces*n_on_face
         ig0 = ngv + n_on_face*(fnum(iface,e)-1)
         if (ifij) then
            k=0
            l=0
            do j=j0,j1,js
            do i=i0,i1,is
               k=k+1
c              this is a serious kludge to stay on the face interior
               if (k.gt.nx.and.k.lt.nxx-nx .and.
     $            mod(k,nx).ne.1.and.mod(k,nx).ne.0) then
c                 interior
                  l = l+1
                  glo_num(i+nx*(j-1)+nxyz*(e-1)) = l + ig0
               endif
            enddo
            enddo
         else
            k=0
            l=0
            do i=i0,i1,is
            do j=j0,j1,js
               k=k+1
c              this is a serious kludge to stay on the face interior
               if (k.gt.nx.and.k.lt.nxx-nx .and.
     $            mod(k,nx).ne.1.and.mod(k,nx).ne.0) then
c                 interior
                  l = l+1
                  glo_num(i+nx*(j-1)+nxyz*(e-1)) = l + ig0
               endif
            enddo
            enddo
         endif
       enddo
      enddo
      ngv   = ngv + ngvs
c
c     Finally,  number interiors (only ifcenter=.true.)
c     -------------------------------------------------
c
      n_in_interior = (nx-2)*(ny-2)*(nz-2)
      ngvi = n_in_interior*nelgt
      if (ifcenter) then
         do e=1,nel
            ig0 = ngv + n_in_interior*(lglel(e)-1)
            l = 0
            do k=2,nz-1
            do j=2,ny-1
            do i=2,nx-1
               l = l+1
               glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(e-1)) = ig0+l
            enddo
            enddo
            enddo
         enddo
         ngv = ngv + ngvi
      else
         do e=1,nel
            l = 0
            do k=2,nz-1
            do j=2,ny-1
            do i=2,nx-1
               l = l+1
               glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(e-1)) = 0
            enddo
            enddo
            enddo
         enddo
      endif
c
c     Quick check on maximum #dofs:
      m    = nxyz*nelt
      ngvm = i8glmax(glo_num,m)
      ngvv = ngvv + ngve + ngvs  ! number of unique ids w/o interior 
      ngvi = ngvi + ngvv         ! total number of unique ids 
      if (nid.eq.0) write(6,1) nx,ngvv,ngvi,ngv,ngvm
    1 format('   setvert3d:',i4,4i12)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine setvert2d(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     setup unique ids for dssum  
c
      include 'SIZE'
      include 'CTIMER'
      include 'PARALLEL'
      include 'TOPOL'
      include 'GEOM'

      integer*8 glo_num(1),ngv
      integer   vertex(0:1,0:1,1),nx
      logical   ifcenter

      integer  edge(0:1,0:1,2,lelt),enum(4,lelt)
      common  /scrmg/ edge,enum

      integer  etuple(4,4*lelt*5)
      integer  ind(4*lelt*5)
      common  /scrns/ ind,etuple

      integer  gvf(4),aa(3),key(3),e,eg
      logical  ifij

      integer*8 igv,ig0
      integer*8 ngvv,ngve,ngvs,ngvi,ngvm
      integer*8 n_on_edge,n_on_face,n_in_interior
      integer*8 i8glmax
c
c     memory check...
c
      ny   = nx
      nz   = 1
      nxyz = nx*ny*nz
c
      key(1)=1
      key(2)=2
      key(3)=3
c
c     Count number of unique vertices
      nlv  = 2**ndim
      ngvv = iglmax(vertex,nlv*nel)
      ngv  = ngvv
c
c     Assign hypercube ordering of vertices.
      do e=1,nel
         do j=0,1
         do i=0,1
c           Local to global node number (vertex)
            il  = 1 + (nx-1)*i + nx*(nx-1)*j
            ile = il + nx*ny*(e-1)
            glo_num(ile)   = vertex(i,j,e)
         enddo
         enddo
      enddo
      if (nx.eq.2) return
c
c     Assign edge labels by bounding vertices.  
      do e=1,nel
         do j=0,1
         do i=0,1
            edge(i,j,1,e) = vertex(i,j,e)  ! r-edge
            edge(j,i,2,e) = vertex(i,j,e)  ! s-edge
       ! write(6,*) 'edge:0 ',e,edge(i,j,1,e),edge(j,i,2,e)
         enddo
         enddo
      enddo
      !call exitt

c     Sort edges by bounding vertices.
      do i=0,4*nel-1
         if (edge(0,i,1,1).gt.edge(1,i,1,1)) then
            kswap = edge(0,i,1,1)
            edge(0,i,1,1) = edge(1,i,1,1)
            edge(1,i,1,1) = kswap
         endif
         etuple(3,i+1) = edge(0,i,1,1)
         etuple(4,i+1) = edge(1,i,1,1)
      enddo

c     Assign a number (rank) to each unique edge
      m    = 4
      n    = 4*nel
      nmax = 4*lelt*5  ! 5x for crystal router factor of safety

      call gbtuple_rank(etuple,m,n,nmax,cr_h,nid,np,ind)
      do i=1,4*nel
         enum(i,1) = etuple(3,i)
      enddo
      n_unique_edges = iglmax(enum,4*nel)

c= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c     Assign global vertex numbers to SEM nodes on each edge
      n_on_edge = nx-2
      do e=1,nel

         iedg_loc = 0

c        Edges 1-2
         do j=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = nx*(nx-1)*j
            i0e = i0 + nxyz*(e-1)
            if (glo_num(i0e+1).lt.glo_num(i0e+nx)) then
               do i=2,nx-1                                   ! std forward case
                  glo_num(i0e+i) = igv + i-1
               enddo
            else
               do i=2,nx-1                                   ! backward case
                  glo_num(i0e+i) = igv + 1 + n_on_edge-(i-1)
               enddo
            endif
            iedg_loc = iedg_loc + 1
         enddo
c
c        Edges 3-4
         do i=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,e)-1)
            i0  = 1+(nx-1)*i
            i0e = i0 + nxyz*(e-1)
            if (glo_num(i0e).lt.glo_num(i0e+nx*(nx-1))) then
               do j=2,nx-1                                   ! std forward case
                  glo_num(i0e+(j-1)*nx) = igv + j-1
               enddo
            else
               do j=2,nx-1                                   ! backward case
                  glo_num(i0e+(j-1)*nx) = igv + 1 + n_on_edge-(j-1)
               enddo
            endif
            iedg_loc = iedg_loc + 1
         enddo
      enddo
 
      ngve = n_unique_edges*n_on_edge
      ngv  = ngv + ngve    
c
c     Finally,  number interiors  
c
      n_in_interior = (nx-2)*(ny-2)
      ngvi          = n_in_interior*nelgt
      if (ifcenter) then
         do e=1,nel
            ig0 = ngv + n_in_interior*(lglel(e)-1)
            l = 0
            do j=2,ny-1
            do i=2,nx-1
               l = l+1
               glo_num(i+nx*(j-1)+nxyz*(e-1)) = ig0+l
            enddo
            enddo
         enddo
         ngv = ngv + ngvi
      else
         do e=1,nel
            l = 0
            do j=2,ny-1
            do i=2,nx-1
               l = l+1
               glo_num(i+nx*(j-1)+nxyz*(e-1)) = 0
            enddo
            enddo
         enddo
      endif

c     do i=1,nel*nx*ny
c        write(6,*) 'glo--',i,glo_num(i)
c     enddo

c     Quick check on maximum #dofs:
      m    = nxyz*nelt
      ngvm = i8glmax(glo_num,m)
      ngvv = ngvv + ngve         ! number of unique ids w/o interior 
      ngvi = ngvi + ngvv         ! total number of unique ids 
      if (nid.eq.0) write(6,1) nx,ngvv,ngvi,ngv,ngvm
    1 format('   setvert2d:',i4,4i12)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine irank_vecn(ind,nn,a,m,n,key,nkey,aa)
c
c     Compute rank of each unique entry a(1,i) 
c
c     Output:   ind(i)  i=1,...,n    (global) rank of entry a(*,i)
c               nn  = max(rank)
c               a(j,i) is permuted
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

      nk = min(nkey,m)
      call ituple_sort(a,m,n,key,nk,ind,aa)

c     Find unique a's
      call icopy(aa,a,m)
      nn     = 1
      ind(1) = nn
c
      do i=2,n
         a_ne_b = iftuple_ianeb(aa,a(1,i),key,nk)
         if (a_ne_b) then
            call icopy(aa,a(1,i),m)
            nn = nn+1
         endif
         ind(i) = nn ! set ind() to rank
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer vertex

      integer e,eg

      integer icalld
      save    icalld
      data    icalld  /0/

      if (icalld.gt.0) return
      icalld = 1

      if (ifhex) ncrnr = 2**ndim
      if (iftet) ncrnr = ndim + 1
      if (nid.eq.0) write(6,*) 'get_vert, ncrnr=',ncrnr 

      call get_vert_map(vertex, ncrnr, nelgt, '.map',ifgfdm)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_map(vertex, nlv, nel, suffix, ifgfdm)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      logical ifgfdm
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer vertex(nlv,1)
      character*4 suffix

      parameter(mdw=2+2**ldim)
      parameter(ndw=7*lx1*ly1*lz1*lelv/mdw)
      common /scrns/ wk(mdw,ndw)   ! room for long ints, if desired
      integer wk,e,eg,eg0,eg1

      character*132 mapfle
      character*1   mapfle1(132)
      equivalence  (mapfle,mapfle1)

      iok = 0
      if (nid.eq.0) then
         lfname = ltrunc(reafle,132) - 4
         call blank (mapfle,132)
         call chcopy(mapfle,reafle,lfname)
         call chcopy(mapfle1(lfname+1),suffix,4)
         open(unit=80,file=mapfle,status='old',err=99)
         read(80,*,err=99) neli,nnzi
         iok = 1
      endif
   99 continue
      iok = iglmax(iok,1)
      if (iok.eq.0) goto 999     ! Mapfile not found

      if (nid.eq.0) then
         neli = iglmax(neli,1)   ! communicate to all procs
      else
         neli = 0
         neli = iglmax(neli,1)   ! communicate neli to all procs
      endif

      npass = 1 + (neli/ndw)

      if (npass.gt.np) then
         if (nid.eq.0) write(6,*) npass,np,neli,ndw,'Error get_vert_map'
         call exitt
      endif 

      len = 4*mdw*ndw
      if (nid.gt.0.and.nid.lt.npass) msg_id=irecv(nid,wk,len)
      call nekgsync

      if (nid.eq.0) then
         eg0 = 0
         do ipass=1,npass
            eg1 = min(eg0+ndw,neli)
            m   = 0
            do eg=eg0+1,eg1
               m = m+1
               read(80,*,end=998) (wk(k,m),k=2,mdw)
               if (.not.ifgfdm) gllnid(eg) = wk(2,m)  !proc map,  must still be divided
               wk(1,m) = eg                           !eg=global node id
            !  write(6,*) (wk(k,m),k=1,mdw)
            enddo
            if (ipass.lt.npass) call csend(ipass,wk,len,ipass,0) !send to ipass
            eg0 = eg1
         enddo
         close(80)
         ntuple = m
      elseif (nid.lt.npass) then
         call msgwait(msg_id)
         ntuple = ndw
      else
         ntuple = 0
      endif

c     Distribute and assign partitions

      if (nid.eq.0) write(6,*) 'ifgfdm',ifgfdm

      if (.not.ifgfdm) then 
      lng = isize*neli
      call bcast(gllnid,lng)
      call assign_gllnid(gllnid,gllel,nelgt,nelgv,np) ! gllel is used as scratch
      endif

c      if(nid.eq.0) then
c         write(99,*) (gllnid(i),i=1,nelgt)
c         write( 6,*) (gllnid(i),i=1,nelgt)
C      endif

      nelt=0 !     Count number of elements on this processor
      nelv=0
      do eg=1,neli
         if (gllnid(eg).eq.nid) then
            if (eg.le.nelgv) nelv=nelv+1
            if (eg.le.nelgt) nelt=nelt+1
         endif
      enddo
      if (np.le.64) write(6,33) nid,nelv,nelt,nelgv,nelgt        
  33  format(5i12,' NELV')                                       

c     now: crystal route vertex by processor id
      do i=1,ntuple
         eg=wk(1,i)
         wk(2,i)=gllnid(eg)        ! processor id for element eg
      !  write(6,*) 'ntuple',nid,ntuple,i,eg,gllnid(eg)
      enddo

      key = 2 ! processor id is in wk(2,:)
      call crystal_ituple_transfer(cr_h,wk,mdw,ntuple,ndw,key)
 
      if (.not.ifgfdm) then            ! no sorting for gfdm?
         key = 1  ! Sort tuple list by eg := wk(1,:)
         nkey= 1
         call crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)
      endif

      iflag = 0
      if (ntuple.ne.nelt) then
         write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FAIL'
         iflag=1
      else
         nv = 2**ndim
         do e=1,nelt
            call icopy(vertex(1,e),wk(3,e),nv)
         enddo
      endif

      iflag = iglmax(iflag,1)
      if (iflag.gt.0) then
         do mid=0,np-1
            call nekgsync
            if (mid.eq.nid)
     $      write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FB'
            call nekgsync
         enddo
         call nekgsync
         call exitt
      endif

      return

  999 continue
      if (nid.eq.0) write(6,*) 'ABORT: Could not find map file ',mapfle
      call exitt

  998 continue
      if (nid.eq.0) write(6,*)ipass,npass,eg0,eg1,mdw,m,eg,'get vX fail'
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine initds
C
C          -- Direct Stiffness Initialization Routine --
C
C     Set up required data for packing data on faces of spectral cubes.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
C     Nominal ordering for direct stiffness summation of faces
C
      J=0
      DO 5 IDIM=1,NDIM
      DO 5 IFACE=1,2
         J=J+1
         NOMLIS(IFACE,IDIM)=J
         !write(6,*) 'nom',IFACE,NDIM,NOMLIS(IFACE,IDIM)
    5 CONTINUE
C
C     Assign Ed's numbering scheme to PF's scheme.
C
      EFACE(1)=4
      EFACE(2)=2
      EFACE(3)=1
      EFACE(4)=3
      EFACE(5)=5
      EFACE(6)=6
C
C     Assign inverse of Ed's numbering scheme to PF's scheme.
C
      EFACE1(1)=3
      EFACE1(2)=2
      EFACE1(3)=4
      EFACE1(4)=1
      EFACE1(5)=5
      EFACE1(6)=6
C
C     Assign group designation to each face to determine ordering of indices.
C
      GROUP(1)=0
      GROUP(2)=1
      GROUP(3)=1
      GROUP(4)=0
      GROUP(5)=0
      GROUP(6)=1
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine setedge
C
C     .Initialize EDGE arrays for face and edge specific tasks.
C
C     .NOTE: Sevaral arrays in common are initialized via 
C            BLOCKDATA EDGEC
C
C     Computed arrays: 
C
C     IEDGE  -  Minimal list of wire frame nodes.
C               Used to search for all physical coincidences.
C
C
C     IEDGEF - .Ordered list of wire frame nodes 
C               associated with faces 1 through 6.
C              .Each of 4 sides of square frame stored 
C               individually so that rotations are 
C               readily handled.
C              .Two types of node orderings stored -
C               (0) is clockwise marching
C               (1) is counter-clockwise marching 
C                   for image face.
C                               
C
C     IFACE         - indicates the face number.  Two notations
C                     are currently in use:
C
C                     i) Preprocessor notation:
C                
C                                       +--------+     ^ S
C                                      /        /|     |
C                                     /    3   / |     |
C                               4--> /        /  |     |
C                                   +--------+ 2 +     +----> R
C                                   |        |  /     /
C                                   |    6   | /     /
C                                   |        |/     /
C                                   +--------+     T
C                                       1
C
C                    ii) Symmetric notation:
C                     
C                                       +--------+     ^ S
C                                      /        /|     |
C                                     /    4   / |     |
C                               1--> /        /  |     |
C                                   +--------+ 2 +     +----> R
C                                   |        |  /     /
C                                   |    6   | /     /
C                                   |        |/     /
C                                   +--------+     T
C                                       3
C
C     EFACE(IFACE)  -   Given face number IFACE in symmetric notation, 
C                       returns preprocessor notation face number.
C
C     EFACE1(IFACE) -   Given face number IFACE in preprocessor notation, 
C                       returns symmetric notation face number.
C
C  The following variables all take the symmetric 
C  notation of IFACE as arguments:
C
C     ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
C                         as depicted below, e.g. ICFACE(i,2)=2,4,6,8.
C
C                        3+-----+4    ^ Y
C                        /  2  /|     |
C     Edge 1 extends    /     / |     |
C       from vertex   7+-----+8 +2    +----> X
C       1 to 2.        |  4  | /     /
C                      |     |/     /
C                     5+-----+6    Z
C                         3
C
C     IEDGFC(i,IFACE) -   Gives the 4 edges which border the face IFACE
C                         Edge numbering is as follows:
C                            Edge = 1,2,3,4     run in +r direction
C                            Edge = 5,6,7,8     run in +s direction
C                            Edge = 9,10,11,12  run in +t direction
C
C                         Ordering of each edge is such that a monotonically
C                         increasing sequence of vertices is associated with
C                         the start point of a corresponding set of 
C                         monotonically increasing edge numbers, e.g.,
C
C     ICEDG(i,IEDGE)  -   Gives 3 variables for determining the stride along
C                         a given edge, IEDGE;  i=1 gives the starting vertex
C                                               i=2 gives the stopping vertex
C                                               i=3 gives the stride size.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
      COMMON /CTMP0/ ITMP(3,3,3)
      INTEGER ORDER
C
      NXL=3
      NYL=3
      NZL=1+2*(NDIM-2)
      NXY   =NXL*NYL
      NXYZ  =NXL*NYL*NZL
      NFACES=2*NDIM
C
C----------------------------------------------------------------------
C     Set up edge arrays (temporary - required only for defining DS)
C----------------------------------------------------------------------
C
C     Fill corners - 1 through 8.
C
      I3D=1
      IF (NDIM.EQ.2) I3D=0
C
      I=0
      DO 10 I3=0,I3D
         IZ=1+(NZL-1)*I3
         DO 10 I2=0,1
            IY=1+(NYL-1)*I2
            DO 10 I1=0,1
               IX=1+(NXL-1)*I1
               I=I+1
               IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
               !write(6,*) 'IEDGE',I,IEDGE(I)
   10 CONTINUE

C     Fill X-direction edges.
C
      DO 20 I3=0,I3D
         IZ=1+(NZL-1)*I3
         DO 20 I2=0,1
            IY=1+(NYL-1)*I2
            DO 20 IX=2,NXL-1
               I=I+1
               IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
               !write(6,*) 'x IEDGE',I,IEDGE(I)
   20 CONTINUE
C
C     Fill Y-direction edges.
C
      DO 30 I3=0,I3D
         IZ=1+(NZL-1)*I3
         DO 30 I1=0,1
            IX=1+(NXL-1)*I1
            DO 30 IY=2,NYL-1
               I=I+1
               IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
               !write(6,*) 'y IEDGE',I,IEDGE(I)
   30 CONTINUE
C
C     Fill Z-direction edges.
C
      IF (NDIM.EQ.3) THEN
         DO 40 I2=0,1
            IY=1+(NYL-1)*I2
            DO 40 I1=0,1
               IX=1+(NXL-1)*I1
               DO 40 IZ=2,NZL-1
                  I=I+1
                  IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
                  !write(6,*) 'z IEDGE',I,IEDGE(I)
   40    CONTINUE
      ENDIF
C
      CALL IZERO(INVEDG,27)

      DO 44 II=1,20
         IX=IEDGE(II)
         INVEDG(IX)=II
         !write(6,*) IX,INVEDG(IX) : not used 
   44 CONTINUE
C
C
C     GENERAL FACE, GENERAL ROTATION EDGE NUMBERS.
C
      IF (NDIM.EQ.3) THEN
C
C        Pack 3-D edge numbering:
C   
C        Fill temporary array with local index numbers:
C
         DO 50 IX=1,NXYZ
            ITMP(IX,1,1)=IX
   50    CONTINUE
C
C        Two sets are required, the base cube and the image cube 
C        which is being summed with it.
C
         DO 1000 IMAGE=0,1
C
C        Pack edges for each face, no rotation.
C
         DO 500 IFACE=1,NFACES
            JS1    = SKPDAT(1,IFACE)
            JF1    = SKPDAT(2,IFACE)
            JSKIP1 = SKPDAT(3,IFACE)
            JS2    = SKPDAT(4,IFACE)
            JF2    = SKPDAT(5,IFACE)
            JSKIP2 = SKPDAT(6,IFACE)
C
C           Choose proper indexing order according to face type and image.
C
            ORDER = (-1)**(GROUP(IFACE)+IMAGE)

            IF (ORDER.EQ.1) THEN
C
C              Forward ordering:
C
C            +-------------+    ^ v1
C            | --------->| |    |
C            | ^    2    | |    +-->
C            | |         | |      v2
C            | |1       3| |
C            | |    4    V |
C            | |<--------- |
C            F-------------I     F is fiducial node.
C
C                                I is location of fiducial node for
C                                     image face.
C
C           Load edge 1:
C
            J=0
            J2=JS2
            DO 100 J1=JS1,JF1-JSKIP1,JSKIP1
               J=J+1
               IEDGEF(J,1,IFACE,IMAGE)=ITMP(J1,J2,1)
  100       CONTINUE
C
C           Load edge 2:
C
            J=0
            J1=JF1
            DO 200 J2=JS2,JF2-JSKIP2,JSKIP2
               J=J+1
               IEDGEF(J,2,IFACE,IMAGE)=ITMP(J1,J2,1)
  200       CONTINUE
C
C
C           Load edge 3:
C
            J=0
            J2=JF2
            DO 300 J1=JF1,JS1+JSKIP1,-JSKIP1
               J=J+1
               IEDGEF(J,3,IFACE,IMAGE)=ITMP(J1,J2,1)
  300       CONTINUE
C
C           Load edge 4:
C
            J=0
            J1=JS1
            DO 400 J2=JF2,JS2+JSKIP2,-JSKIP2
               J=J+1
               IEDGEF(J,4,IFACE,IMAGE)=ITMP(J1,J2,1)
  400       CONTINUE
C
            ELSE
C
C           Reverse ordering:
C
C            +-------------+
C            | |<--------- |       ^ v2
C            | |    2    ^ |       |
C            | |         | |    <--+
C            | |3       1| |     v1
C            | V    4    | |
C            | --------->| |
C            I-------------F     F is fiducial node.
C
C                                I is location of fiducial node for
C                                     image face.
C
C           Load edge 1:
C   
            J=0
            J1=JS1
            DO 105 J2=JS2,JF2-JSKIP2,JSKIP2
               J=J+1
               IEDGEF(J,1,IFACE,IMAGE)=ITMP(J1,J2,1)
  105       CONTINUE
C
C           Load edge 2:
C
            J=0
            J2=JF2
            DO 205 J1=JS1,JF1-JSKIP1,JSKIP1
               J=J+1
               IEDGEF(J,2,IFACE,IMAGE)=ITMP(J1,J2,1)
  205       CONTINUE
C
C           Load edge 3:
C
            J=0
            J1=JF1
            DO 305 J2=JF2,JS2+JSKIP2,-JSKIP2
                J=J+1
               IEDGEF(J,3,IFACE,IMAGE)=ITMP(J1,J2,1)
  305       CONTINUE
C
C           Load edge 4:
C
            J=0
            J2=JS2
            DO 405 J1=JF1,JS1+JSKIP1,-JSKIP1
               J=J+1
              IEDGEF(J,4,IFACE,IMAGE)=ITMP(J1,J2,1)
  405       CONTINUE
            ENDIF
C
  500    CONTINUE
 1000    CONTINUE
      ELSE
C
C        Load edge information for 2-D case
C
         IEDGEF(1,1,1,0) = NXY - NXL + 1
         IEDGEF(1,2,1,0) = 1
         IEDGEF(1,1,2,0) = NXL
         IEDGEF(1,2,2,0) = NXY
         IEDGEF(1,1,3,0) = 1
         IEDGEF(1,2,3,0) = NXL
         IEDGEF(1,1,4,0) = NXY
         IEDGEF(1,2,4,0) = NXY - NXL + 1
C
         IEDGEF(1,1,1,1) = 1
         IEDGEF(1,2,1,1) = NXY - NXL + 1
         IEDGEF(1,1,2,1) = NXY
         IEDGEF(1,2,2,1) = NXL
         IEDGEF(1,1,3,1) = NXL
         IEDGEF(1,2,3,1) = 1
         IEDGEF(1,1,4,1) = NXY - NXL + 1
         IEDGEF(1,2,4,1) = NXY
c        write(6,*) 'IEDGEF = '
c    $,  IEDGEF(1,1,1,0),IEDGEF(1,2,1,0),IEDGEF(1,1,2,0)
c    $,  IEDGEF(1,2,2,0),IEDGEF(1,1,3,0),IEDGEF(1,2,3,0)
c    $,  IEDGEF(1,1,4,0),IEDGEF(1,2,4,0)
c        write(6,*) 'IEDGEF = '
c    $,  IEDGEF(1,1,1,1),IEDGEF(1,2,1,1),IEDGEF(1,1,2,1)
c    $,  IEDGEF(1,2,2,1),IEDGEF(1,1,3,1),IEDGEF(1,2,3,1)
c    $,  IEDGEF(1,1,4,1),IEDGEF(1,2,4,1)

      ENDIF
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine dsset(nx,ny,nz)
C
C     Set up arrays IXCN,ESKIP,SKPDAT,NEDG,NOFFST for new NX,NY,NZ
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TOPOL'

      INTEGER NXO,NYO,NZO
      SAVE    NXO,NYO,NZO
      DATA    NXO,NYO,NZO /3*0/
C
C     Check if element surface counters are already set from last call...
C
      IF (NXO.EQ.NX.AND.NYO.EQ.NY.AND.NZO.EQ.NZ) RETURN
C
C     else, proceed....
C
      NXO = NX
      NYO = NY
      NZO = NZ
C
C     Establish corner to elemental node number mappings
C
      DUMMY2 = 0
      IC=0
      DO 10 ICZ=0,1
      DO 10 ICY=0,1
      DO 10 ICX=0,1
C       Supress vectorization to 
        IF(ICX.EQ.0)DUMMY=0        
        IF(ICX.EQ.1)DUMMY=1
        DUMMY2=DUMMY2+DUMMY

        IC=IC+1
        IXCN(IC)= 1 + (NX-1)*ICX + NX*(NY-1)*ICY + NX*NY*(NZ-1)*ICZ
   10   CONTINUE
 
C     Assign indices for direct stiffness summation of arbitrary faces.
C
C
C     Y-Z Planes (Faces 1 and 2)
C
      SKPDAT(1,1)=1
      SKPDAT(2,1)=NX*(NY-1)+1
      SKPDAT(3,1)=NX
      SKPDAT(4,1)=1
      SKPDAT(5,1)=NY*(NZ-1)+1
      SKPDAT(6,1)=NY
      
C
      SKPDAT(1,2)=1             + (NX-1)
      SKPDAT(2,2)=NX*(NY-1)+1   + (NX-1)
      SKPDAT(3,2)=NX
      SKPDAT(4,2)=1
      SKPDAT(5,2)=NY*(NZ-1)+1
      SKPDAT(6,2)=NY
C
C     X-Z Planes (Faces 3 and 4)
C
      SKPDAT(1,3)=1
      SKPDAT(2,3)=NX
      SKPDAT(3,3)=1
      SKPDAT(4,3)=1
      SKPDAT(5,3)=NY*(NZ-1)+1
      SKPDAT(6,3)=NY
C
      SKPDAT(1,4)=1           + NX*(NY-1)
      SKPDAT(2,4)=NX          + NX*(NY-1)
      SKPDAT(3,4)=1
      SKPDAT(4,4)=1
      SKPDAT(5,4)=NY*(NZ-1)+1
      SKPDAT(6,4)=NY
C
C     X-Y Planes (Faces 5 and 6)
C
      SKPDAT(1,5)=1
      SKPDAT(2,5)=NX
      SKPDAT(3,5)=1
      SKPDAT(4,5)=1
      SKPDAT(5,5)=NY
      SKPDAT(6,5)=1
C
      SKPDAT(1,6)=1           + NX*NY*(NZ-1)
      SKPDAT(2,6)=NX          + NX*NY*(NZ-1)
      SKPDAT(3,6)=1
      SKPDAT(4,6)=1
      SKPDAT(5,6)=NY
      SKPDAT(6,6)=1
       
C
C     Set up skip indices for each of the 12 edges 
C
C         Note that NXY = NX*NY even for 2-D since
C         this branch does not apply to the 2D case anyway.
C
C     ESKIP(*,1) = start location
C     ESKIP(*,2) = end 
C     ESKIP(*,3) = stride
C

      NXY=NX*NY
      ESKIP( 1,1) = IXCN(1) + 1
      ESKIP( 1,2) = IXCN(2) - 1
      ESKIP( 1,3) = 1
      ESKIP( 2,1) = IXCN(3) + 1
      ESKIP( 2,2) = IXCN(4) - 1
      ESKIP( 2,3) = 1
      ESKIP( 3,1) = IXCN(5) + 1
      ESKIP( 3,2) = IXCN(6) - 1
      ESKIP( 3,3) = 1
      ESKIP( 4,1) = IXCN(7) + 1
      ESKIP( 4,2) = IXCN(8) - 1
      ESKIP( 4,3) = 1
      ESKIP( 5,1) = IXCN(1) + NX
      ESKIP( 5,2) = IXCN(3) - NX
      ESKIP( 5,3) = NX
      ESKIP( 6,1) = IXCN(2) + NX
      ESKIP( 6,2) = IXCN(4) - NX
      ESKIP( 6,3) = NX
      ESKIP( 7,1) = IXCN(5) + NX
      ESKIP( 7,2) = IXCN(7) - NX
      ESKIP( 7,3) = NX
      ESKIP( 8,1) = IXCN(6) + NX
      ESKIP( 8,2) = IXCN(8) - NX
      ESKIP( 8,3) = NX
      ESKIP( 9,1) = IXCN(1) + NXY
      ESKIP( 9,2) = IXCN(5) - NXY
      ESKIP( 9,3) = NXY
      ESKIP(10,1) = IXCN(2) + NXY
      ESKIP(10,2) = IXCN(6) - NXY
      ESKIP(10,3) = NXY
      ESKIP(11,1) = IXCN(3) + NXY
      ESKIP(11,2) = IXCN(7) - NXY
      ESKIP(11,3) = NXY
      ESKIP(12,1) = IXCN(4) + NXY
      ESKIP(12,2) = IXCN(8) - NXY
      ESKIP(12,3) = NXY
C
C     Load reverse direction edge arrays for reverse mappings...
C
      DO 20 IED=1,12
         IEDM=-IED
         ESKIP(IEDM,1) =  ESKIP(IED,2)
         ESKIP(IEDM,2) =  ESKIP(IED,1)
         ESKIP(IEDM,3) = -ESKIP(IED,3)
   20 CONTINUE
C
C     Compute offset for global edge vector given current element
C     dimensions.....
C
C     NGSPED(ITE,ICMP) = number of global (ie, distinct) special edges
C                        of type ITE (1,2, or 3)  for field ICMP.
C
C                        ITE = 1 implies an "X" edge
C                        ITE = 2 implies an "Y" edge
C                        ITE = 3 implies an "Z" edge
C
C     Set up number of nodes along each of the 3 types of edges
C     (endpoints excluded).
C
      NEDG(1)=NX-2
      NEDG(2)=NY-2
      NEDG(3)=NZ-2
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine genxyzl
C
C     Generate xyz coordinates 
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SCRCT'
      COMMON /CTMP0/ XCB(2,2,2),YCB(2,2,2),ZCB(2,2,2),H(3,3,2),INDX(8)
C
      NXL=3   
      NYL=3   
      NZL=1+2*(NDIM-2)
      NTOT3=NXL*NYL*NZL*NELT
C
C   Preprocessor Corner notation:      Symmetric Corner notation:
C                
C           4+-----+3    ^ s                    3+-----+4    ^ s        
C           /     /|     |                      /     /|     |          
C          /     / |     |                     /     / |     |          
C        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r   
C         |     | /     /                     |     | /     /           
C         |     |/     /                      |     |/     /            
C        5+-----+6    t                      5+-----+6    t             
C
      DO 10 IX=1,NXL
         H(IX,1,1)=0.5*FLOAT(3-IX)
         H(IX,1,2)=0.5*FLOAT(IX-1)
   10 CONTINUE
      DO 20 IY=1,NYL
         H(IY,2,1)=0.5*FLOAT(3-IY)
         H(IY,2,2)=0.5*FLOAT(IY-1)
   20 CONTINUE
      DO 30 IZ=1,NZL
         H(IZ,3,1)=0.5*FLOAT(3-IZ)
         H(IZ,3,2)=0.5*FLOAT(IZ-1)
   30 CONTINUE
C
      INDX(1)=1
      INDX(2)=2
      INDX(3)=4
      INDX(4)=3
      INDX(5)=5
      INDX(6)=6
      INDX(7)=8
      INDX(8)=7
C
      CALL RZERO(XML,NTOT3)
      CALL RZERO(YML,NTOT3)
      CALL RZERO(ZML,NTOT3)
      CALL RZERO(XCB,8)
      CALL RZERO(YCB,8)
      CALL RZERO(ZCB,8)
C
      DO 5000 IE=1,NELT
C
         NDIM2 = 2**NDIM
         DO 50 IX=1,NDIM2
            I=INDX(IX)
            XCB(IX,1,1)=XC(I,IE)
            YCB(IX,1,1)=YC(I,IE)
            ZCB(IX,1,1)=ZC(I,IE)
   50    CONTINUE
C
C        Map R-S-T space into physical X-Y-Z space.
C
         DO 100 IZT=1,ndim-1
         DO 100 IYT=1,2
         DO 100 IXT=1,2
C
         DO 100 IZ=1,NZL
         DO 100 IY=1,NYL
         DO 100 IX=1,NXL
            XML(IX,IY,IZ,IE)=XML(IX,IY,IZ,IE)+
     $        H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*XCB(IXT,IYT,IZT)
            YML(IX,IY,IZ,IE)=YML(IX,IY,IZ,IE)+
     $        H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*YCB(IXT,IYT,IZT)
            ZML(IX,IY,IZ,IE)=ZML(IX,IY,IZ,IE)+
     $        H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*ZCB(IXT,IYT,IZT)
  100    CONTINUE
C
 5000 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine verify
C
C     .Verify right-handedness of elements.
C     .Verify element-to-element reciprocity of BC's
C     .Verify correlation between E-E BC's and physical coincidence
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'INPUT'
      INCLUDE 'SCRCT'
C
      CALL VERRHE

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine setside
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TOPOL'
      INCLUDE 'SCRCT'
C
C     SIDE(i,IFACE,IE) -  Physical (xyz) location of element side midpoint.
C                         i=1,2,3 gives x,y,z value, respectively.
C                         i=4  gives average dimension of face for setting
C                              tolerances.
C
      INDX(1)=1
      INDX(2)=2
      INDX(3)=4
      INDX(4)=3
      INDX(5)=5
      INDX(6)=6
      INDX(7)=8
      INDX(8)=7
C
C     Flip vertex array structure
C
      CALL RZERO(XYZ,24*NELT)

      call rzero(xyz,24*nelt)
      if (if3d) then
         do ie=1,nelt
         do j=1,8
            ivtx = indx(j)
            xyz(1,ivtx,ie) = xc(j,ie)
            xyz(2,ivtx,ie) = yc(j,ie)
            xyz(3,ivtx,ie) = zc(j,ie)
c           write(6,1) ie,j,ivtx,xc(j,ie),yc(j,ie),zc(j,ie),' xcz'
c           write(6,1) ie,j,ivtx,(xyz(k,ivtx,ie),k=1,3),' vtx'
c   1       format(3i5,1p3e12.4,a4)
         enddo
         enddo
      else
         do ie=1,nelt
         do j=1,4
            ivtx = indx(j)
            xyz(1,ivtx,ie) = xc(j,ie)
            xyz(2,ivtx,ie) = yc(j,ie)
            xyz(3,ivtx,ie) = 0.0
         enddo
         enddo
      endif

C
C     Compute location of center and "diameter" of each element side.
C
      NFACES=NDIM*2
      NCRNR =2**(NDIM-1)
      CALL RZERO(SIDE,24*NELT)
      DO 500 ICRN=1,NCRNR
      DO 500 IFAC=1,NFACES
         IVTX = ICFACE(ICRN,IFAC)
         ICR1 = NCRNR+(ICRN-1)
         ICR1 = MOD1(ICR1,NCRNR)
         IVT1 = ICFACE(ICR1,IFAC)
         DO 400 IE=1,NELT
            DO 300 IDIM=1,NDIM
               SIDE(IDIM,IFAC,IE)=SIDE(IDIM,IFAC,IE)+XYZ(IDIM,IVTX,IE)
               SIDE(   4,IFAC,IE)=SIDE(   4,IFAC,IE)+
     $                     ( XYZ(IDIM,IVTX,IE)-XYZ(IDIM,IVT1,IE) )**2
  300       CONTINUE
            SIDE(4,IFAC,IE)=SQRT( SIDE(4,IFAC,IE) )
  400    CONTINUE
  500 CONTINUE
      AVWGHT=1.0/FLOAT(NCRNR)
      CALL CMULT(SIDE,AVWGHT,24*NELT)
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine verrhe
C     8 Mar 1989 21:58:26   PFF
C     Verify right-handedness of given elements. 
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      INCLUDE 'TOPOL'
      LOGICAL IFYES,IFCSTT
C
      IFCSTT=.TRUE.

      IF (.NOT.IF3D) THEN
      DO 1000 IE=1,NELT
C
C        CRSS2D(A,B,O) = (A-O) X (B-O)
C
         C1=CRSS2D(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,1,IE))
         C2=CRSS2D(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,2,IE))
         C3=CRSS2D(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,3,IE))
         C4=CRSS2D(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,4,IE))
C
         IF (C1.LE.0.0.OR.C2.LE.0.0.OR.
     $       C3.LE.0.0.OR.C4.LE.0.0 ) THEN
C
            IEG=LGLEL(IE)
            WRITE(6,800) IEG,C1,C2,C3,C4
  800       FORMAT(/,2X,'WARNING:  Detected non-right-handed element.',
     $      /,2X,'Number',I4,'  C1-4:',4E12.4)
            IFCSTT=.FALSE.
         ENDIF
 1000 CONTINUE
C
C     Else 3-D:
C
      ELSE      

      DO 2000 IE=1,NELT
C
C        VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
C
         V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
         V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
         V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
         V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
         V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
         V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
         V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
         V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
C
         IF (V1.LE.0.0.OR.V2.LE.0.0.OR.
     $       V3.LE.0.0.OR.V4.LE.0.0.OR.
     $       V5.LE.0.0.OR.V6.LE.0.0.OR.
     $       V7.LE.0.0.OR.V8.LE.0.0    ) THEN
C
            IEG=LGLEL(IE)
            WRITE(6,1800) IEG,V1,V2,V3,V4,V5,V6,V7,V8
 1800       FORMAT(/,2X,'WARNING:  Detected non-right-handed element.',
     $      /,2X,'Number',I4,'  V1-8:',4E12.4
     $      /,2X,'      ',4X,'       ',4E12.4)
            IFCSTT=.FALSE.
         ENDIF
 2000 CONTINUE
      ENDIF
C
 9000 CONTINUE
C
C     Print out results from right-handed check
C
      IF (.NOT.IFCSTT) WRITE(6,2001)
C
C     Check consistency accross all processors.
C
      CALL GLLOG(IFCSTT,.FALSE.)
      !if (nid.eq.0) write(6,*) 'EE------t4'

      IF (.NOT.IFCSTT) THEN
         IF (NID.EQ.0) WRITE(6,2003) NELGT
         call exitt
      ELSE
         IF (NID.EQ.0) WRITE(6,2002) NELGT
      ENDIF
C
 2001 FORMAT(//,'  Elemental geometry not right-handed, ABORTING'
     $      ,' in routine VERRHE.')
 2002 FORMAT('   Right-handed check complete for',I4,' elements. OK.')
 2003 FORMAT('   Right-handed check failed for',I4,' elements.'
     $      ,'  Exiting in routine VERRHE.')

      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION VOLUM0(P1,P2,P3,P0)
C
C                           3
C     Given four points in R , (P1,P2,P3,P0), VOLUM0 returns
C     the volume enclosed by the parallelagram defined by the
C     vectors { (P1-P0),(P2-P0),(P3-P0) }.  This routine has
C     the nice feature that if the 3 vectors so defined are
C     not right-handed then the volume returned is negative.
C
      REAL P1(3),P2(3),P3(3),P0(3)
C
         U1=P1(1)-P0(1)
         U2=P1(2)-P0(2)
         U3=P1(3)-P0(3)
C
         V1=P2(1)-P0(1)
         V2=P2(2)-P0(2)
         V3=P2(3)-P0(3)
C
         W1=P3(1)-P0(1)
         W2=P3(2)-P0(2)
         W3=P3(3)-P0(3)
C
         CROSS1 = U2*V3-U3*V2
         CROSS2 = U3*V1-U1*V3
         CROSS3 = U1*V2-U2*V1
C
         VOLUM0  = W1*CROSS1 + W2*CROSS2 + W3*CROSS3
         
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION CRSS2D(XY1,XY2,XY0)
      REAL XY1(2),XY2(2),XY0(2)
C
         V1X=XY1(1)-XY0(1)
         V2X=XY2(1)-XY0(1)
         V1Y=XY1(2)-XY0(2)
         V2Y=XY2(2)-XY0(2)
         CRSS2D = V1X*V2Y - V1Y*V2X
C
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
       KX1=1
       KY1=1
       KZ1=1
       KX2=NX1
       KY2=NY1
       KZ2=NZ1    ! Ed's order 
       IF (IFACE.EQ.1) KY2=1
       IF (IFACE.EQ.2) KX1=NX1
       IF (IFACE.EQ.3) KY1=NY1
       IF (IFACE.EQ.4) KX2=1
       IF (IFACE.EQ.5) KZ2=1
       IF (IFACE.EQ.6) KZ1=NZ1
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine facev(a,ie,iface,val,nx,ny,nz)
C
C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      INCLUDE 'SIZE'
      DIMENSION A(NX,NY,NZ,LELT)

      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
         A(IX,IY,IZ,IE)=VAL
  100 CONTINUE

      RETURN
      END

c-----------------------------------------------------------------------
      subroutine ifacev(a,ie,iface,val,nx,ny,nz)
C
C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      include 'SIZE'
      integer a(nx,ny,nz,lelt),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz,ie)=val
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine facec(a,b,ie,iface,nx,ny,nz,nel)
C
C     Copy the face (IFACE) of B to A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      DIMENSION A(NX,NY,NZ,NEL)
      DIMENSION B(NX,NY,NZ,NEL)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
         A(IX,IY,IZ,IE)=B(IX,IY,IZ,IE)
  100 CONTINUE

      RETURN
      END

c-----------------------------------------------------------------------
      subroutine outfldio (x,txt)
      INCLUDE 'SIZE'
      integer x(lx1,ly1,lz1,lelt)
      character*10 txt
C
      do ie=1,nelv,2
         do iz=1,nz1,1
            if (iz.eq.1) write(6,106) txt,iz,ie
            if (iz.gt.1) write(6,107) 
            i1 = ie+1
            do j=ny1,1,-1
               write(6,105) (x(i,j,iz,ie),i=1,nx1)
     $                    , (x(i,j,iz,i1),i=1,nx1)
            enddo
         enddo
      enddo
C
  107 FORMAT(' ')
  105 FORMAT(4i6,20x,4i6)
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',I2,/,
     $           5X,'       X            ')
C
      return
      end

c-----------------------------------------------------------------------
      subroutine outfldi (x,txt)
      INCLUDE 'SIZE'
      integer x(lx1,ly1,lz1,lelt)
      character*10 txt
c
      character*6 s(20,20)
c
      if (lx1.ne.4 .or. nelv.gt.3) return
c
      ! FIXME
      ! ie defined below, set to zero here. ak, 2006-07-31
      ie = 0

      write(6,106) txt,ie,ie
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','elem. = ',I2,'/',I2,/,
     $           5X,'       X            ')
C
C
      call blank(s,6*20*20)
      do ie=1,3
         if (ie.eq.1) then
            jstart  = 1
            istart  = 1
            istride = 3
         else
            jstart  = 2 + lx1
            istart  = 1
            istride = 1
            if (ie.eq.2) istart = 7
         endif
c
         i=istart
         do iy=ny1,1,-1
            j=jstart
            do ix=1,nx1
               write(s(i,j),6) x(ix,iy,1,ie)
               j=j+1
            enddo
            i=i+istride
         enddo
    6    format(20i6)
      enddo
c
      do i=1,10
         write(6,7) (s(i,l),l=1,j-1)
      enddo
    7 format(20a6)

      write(6,*)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldr (x,txt)
      INCLUDE 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt
c
      character*6 s(20,20)
c
      if (lx1.ne.4 .or. nelv.gt.3) return

      ! FIXME
      ! ie defined below, set to zero here. ak, 2006-07-31
      ie = 0

      write(6,106) txt,ie,ie
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','elem. = ',I2,'/',I2,/,
     $           5X,'       X            ')
      call blank(s,6*20*20)
C
C
      do ie=1,2
         if (ie.eq.1) then
            jstart  = 1
            istart  = 1
            istride = 3
         else
            jstart  = 2 + lx1
            istart  = 1
            istride = 1
            if (ie.eq.2) istart = 7
         endif
c
         i=istart
         do iy=ny1,1,-1
            j=jstart
            do ix=1,nx1
               write(s(i,j),6) x(ix,iy,1,ie)
               j=j+1
            enddo
            i=i+istride
         enddo
    6    format(f6.2)
      enddo
c
      do i=1,10
         write(6,7) (s(i,l),l=1,j-1)
      enddo
    7 format(20a6)
c
      write(6,*)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine checkit(idum)
      write(6,*) 'continue?'
      read (5,*) idum

      return
      end

c-----------------------------------------------------------------------
      subroutine outfldface (x,txt,ichk)

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      real x(nx1*nz1,2*ndim,lelt)
      character*10 txt
c
      integer idum,e,ef
      save idum
      data idum /3/

      integer f
   
      rnel = nelv
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelv-ne+1

      write(6,117) 

      do ie = ne1,1,-ne
         l  = ie-1
         write(6,117)
         do f=1, 4
            ef=eface(f)
            if (nx1.eq.5) write(6,115) ((x(i,ef,e+l),i=1,nx1),e=1,ne)
            if (nx1.eq.6) write(6,116) ((x(i,ef,e+l),i=1,nx1),e=1,ne)
         enddo
      enddo

  115 FORMAT(3(5f10.7,10x))
  116 FORMAT(4(6f10.7,10x))
  117 FORMAT(' ')

      return
      end

c-----------------------------------------------------------------------
      subroutine setup_dgds2
     $           (gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      integer*8 glo_num(1),ngv
      integer   gs_handle
      integer   vertex(1)
      real*8    dclock,t0,t00

      common /c_is2/ glo_num1((lx1+2)*(ly1+2)*(lz1+2)*lelt)
      integer*8 glo_num1

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      t0 = dclock()

c     Global-to-local mapping for gs
      mx = nx+2

      if (nid.eq.0) write(6,*) 'setup_dgds2:: init'

      t00= dclock()           
      call set_vert2(glo_num1,ngv,mx,nel,vertex,.false.) !2nd call
      if (nid.eq.0)  write(6,2)  dclock()-t00
    2 format(' set_vert time::',1pe11.4,' seconds ')

      ntot = (nx1+2)*(ny1+2)*nelt
      if (if3d) ntot = (nx1+2)*(ny1+2)*(nz1+2)*nelt
c
      mz0 = 1
      mz1 = 1

      if (if3d) mz0 = 0
      if (if3d) mz1 = nz1+1
 
c     call icopy8          (glo_num1,glo_num,ntot)
      call iface_vert_int8 (glo_num,glo_num1,mz0,mz1,nelt) 

      nface= 2*ndim
      ntot = nx1*nz1*nface*nelt !total grids on faces

      t00= dclock()
      call gs_setup(gsh_face,glo_num,ntot,nekcomm,np)
      if (nid.eq.0)  write(6,3)  dclock()-t00
    3 format('   gs_setup time::',1pe11.4,' seconds ')

      if (nid.eq.0) then
         write(6,5) dclock()-t0,gsh_face,nx,ngv,melg
    5    format(' setup_dsdg2 time::',1pe11.4,' seconds ',2i3,2i12)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine iface_vert2d(fa,va,jz0,jz1,nel)
      include 'SIZE'
      integer  e,f,jz0,jz1,nel
      integer  fa(nx1*nz1,2*ndim,nel)
      integer  va(1:nx1,1:ny1,jz0:jz1,nel)
c
      ntot = nx1*nz1*2*ndim*nel
      call izero(fa,ntot)
c
      mx1 = nx1 +2
      my1 = ny1 +2
      mz1 = nz1 +2

      if (ndim.eq.2) mz1=1
c
      nface= 2*ndim
      do e = 1,nel
      do f = 1,nface

         call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
 
         if     (f.eq.1) then ! EB notation
            ky1=ky1 !-1
            ky2=ky1
         elseif (f.eq.2) then
            kx1=kx1 !+1
            kx2=kx1
         elseif (f.eq.3) then
            ky1=ky1 !+1
            ky2=ky1
         elseif (f.eq.4) then
            kx1=kx1 !-1
            kx2=kx1
         elseif (f.eq.5) then
            kz1=kz1 !-1
            kz2=kz1
         elseif (f.eq.6) then
            kz1=kz1 !+1
            kz2=kz1
         endif
 
         i = 0
         do iz= kz1,kz2
         do iy= ky1,ky2
         do ix= kx1,kx2
            i = i+1
            fa(i,f,e)= va(ix,iy,iz,e)
            !write(6,*) 'fa',i,f,e,fa(i,f,e)
         enddo
         enddo
         enddo

      enddo
      enddo
c                                                                                                                    
      return
      end
c-----------------------------------------------------------------------
      subroutine iface_vert_int8(fa,va,jz0,jz1,nel)
      include 'SIZE'
      integer*8 fa(nx1*nz1,2*ndim,nel),va(0:nx1+1,0:ny1+1,jz0:jz1,nel)
c     integer fa(1,1,1),va(1)
      integer e,f
c
      ntot = nx1*nz1*2*ndim*nel
      call izero8(fa,ntot)
c
      mx1 = nx1+2
      my1 = ny1+2
      mz1 = nz1+2

      if (ndim.eq.2) mz1=1
c
      nface = 2*ndim
      do e=1,nel
      do f=1,nface
         call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
c
         if     (f.eq.1) then ! EB notation
            ky1=ky1-1
            ky2=ky1
         elseif (f.eq.2) then
            kx1=kx1+1
            kx2=kx1
         elseif (f.eq.3) then
            ky1=ky1+1
            ky2=ky1
         elseif (f.eq.4) then
            kx1=kx1-1
            kx2=kx1
         elseif (f.eq.5) then
            kz1=kz1-1
            kz2=kz1
         elseif (f.eq.6) then
            kz1=kz1+1
            kz2=kz1
         endif
c
         i = 0
         do iz=kz1,kz2
         do iy=ky1,ky2
         do ix=kx1,kx2
            i = i+1
            fa(i,f,e)=va(ix,iy,iz,e)
           !write(6,*) 'va',ix,iy,iz,fa(i,f,e)
         enddo
         enddo
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine iface_vert2(fa,va,jz0,jz1,nel)
      include 'SIZE'
      integer fa(nx1*nz1,2*ndim,nel),va(0:nx1+1,0:ny1+1,jz0:jz1,nel)
      integer e,f
c
      ntot = nx1*nz1*2*ndim*nel
      call izero(fa,ntot)
c
      mx1 = nx1+2
      my1 = ny1+2
      mz1 = nz1+2

      if (ndim.eq.2) mz1=1
c
      nface = 2*ndim
      do e=1,nel
      do f=1,nface
         call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
c
         if     (f.eq.1) then ! EB notation
            ky1=ky1-1
            ky2=ky1
         elseif (f.eq.2) then
            kx1=kx1+1
            kx2=kx1
         elseif (f.eq.3) then
            ky1=ky1+1
            ky2=ky1
         elseif (f.eq.4) then
            kx1=kx1-1
            kx2=kx1
         elseif (f.eq.5) then
            kz1=kz1-1
            kz2=kz1
         elseif (f.eq.6) then
            kz1=kz1+1
            kz2=kz1
         endif
c
         i = 0
         do iz=kz1,kz2
         do iy=ky1,ky2
         do ix=kx1,kx2
            i = i+1
            fa(i,f,e)=va(ix,iy,iz,e)
C           write(6,*) 'va',ix,iy,iz,fa(i,f,e)
         enddo
         enddo
         enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine f77_get_vert (vertex,nlv,nel,suffix)
      include 'SIZE'
      include 'INPUT'
      integer vertex(nlv,nel)
      character*4 suffix
c
      character*132 mapfle
      character*1   mapfle1(132)
      equivalence  (mapfle,mapfle1)
c
      lfname = ltrunc(reafle,132) - 4
      call chcopy(mapfle,reafle,lfname)
      call chcopy(mapfle1(lfname+1),suffix,4)
c
      write(6,*) 'lfname:',lfname
      write(6,*) 'suffix:',suffix
      write(6,*) 'mapfle:',mapfle
      write(6,*) 'reafle:',reafle
      open(unit=80,file=mapfle,status='old',err=999)
c
      read(80,*) neli,nnzi
c
c     Get vertices from rsb routine... hypercubic ordering
c
      do ie=1,neli
         read(80,*) iproc,(vertex(k,ie),k=1,nlv)
      enddo
c
      close(unit=80)
      if (neli.eq.nel) return
c
   99 continue
      write(6,*) 'ABORT in f77_get_vert',neli,nel
      call exitt
c
  999 continue
      write(6,*) 'ABORT in f77_get_vert',mapfle
      call exitt
      return
      end
c-----------------------------------------------------------------------
      subroutine bindec(bin_in)
      integer bin_in,d,b,b2
c
      keep  = bin_in
      d  = bin_in
      b2 = 1
      b  = 0
      do l=1,12
         b  = b + b2*mod(d,10)
         d  = d/10
         b2 = b2*2
         if (d.eq.0) goto 1
      enddo
    1 continue
      bin_in = b
      return
      end
c
c-----------------------------------------------------------------------
      subroutine ifacev_redef(a,ie,iface,val,nx,ny,nz)
C
C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      include 'SIZE'
      integer a(nx,ny,nz,lelt),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz,ie)=val
  100 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine assign_gllnid(gllnid,iunsort,nelgt,nelgv,np)
c
      integer gllnid(1),iunsort(1),nelgt,np 
      integer e,eg


      log2p = log2(np)
      np2   = 2**log2p
      if (np2.eq.np.and.nelgv.eq.nelgt) then   ! std power of 2 case

         npstar = ivlmax(gllnid,nelgt)+1
         nnpstr = npstar/np
         do eg=1,nelgt
            gllnid(eg) = gllnid(eg)/nnpstr
         enddo

         return

      elseif (np2.eq.np) then   ! std power of 2 case, conjugate heat xfer

c        Assign fluid elements
         npstar = max(np,ivlmax(gllnid,nelgv)+1)
         nnpstr = npstar/np
         do eg=1,nelgv
            gllnid(eg) = gllnid(eg)/nnpstr
         enddo

c        Assign solid elements
         nelgs  = nelgt-nelgv  ! number of solid elements
         npstar = max(np,ivlmax(gllnid(nelgv+1),nelgs)+1)
         nnpstr = npstar/np
         do eg=nelgv+1,nelgt
            gllnid(eg) = gllnid(eg)/nnpstr
         enddo

         return

      elseif (nelgv.ne.nelgt) then
c        call exitti
c    $       ('Conjugate heat transfer requires P=power of 2.$',np)
      endif


c  Below is the code for P a non-power of two:

c  Split the sorted gllnid array (read from .map file) 
c  into np contiguous partitions. 

c  To load balance the partitions in case of mod(nelgt,np)>0 
c  add 1 contiguous entry out of the sorted list to NODE_i 
c  where i = np-mod(nelgt,np) ... np


      nel   = nelgt/np       ! number of elements per processor
      nmod  = mod(nelgt,np)  ! bounded between 1 ... np-1
      npp   = np - nmod      ! how many paritions of size nel 
 
      ! sort gllnid  
      call isort(gllnid,iunsort,nelgt)

      ! setup partitions of size nel 
      k   = 0
      do ip = 0,npp-1
         do e = 1,nel  
            k = k + 1 
            gllnid(k) = ip
         enddo
      enddo
      ! setup partitions of size nel+1
      if(nmod.gt.0) then 
        do ip = npp,np-1
           do e = 1,nel+1  
              k = k + 1 
              gllnid(k) = ip
           enddo
        enddo 
      endif

      ! unddo sorting to restore initial ordering by
      ! global element number
      call iswapt_ip(gllnid,iunsort,nelgt)

      return
      end
c-----------------------------------------------------------------------
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end

c-----------------------------------------------------------------------


C=======================================================================
#ifdef GPU

#ifdef _OPENACC

c-----------------------------------------------------------------------
      subroutine isort8(a,ind,n) 
c     Sort routine for a = int*8, ind=int.
c     Uses heap sort (p 231 Num. Rec., 1st Ed.)

      integer*8 a(1),aa
      integer ind(1)

      do j=1,n
         ind(j)=j
      enddo

      if (n.le.1) return
      L=n/2+1
      ir=n
 100   continue
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
 200         continue
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
         goto 200
         endif
           a(i) = aa
         ind(i) = ii
      goto 100
      end

c-----------------------------------------------------------------------
      subroutine setup_dgds2_acc
     $           (gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      integer*8 glo_num(1),ngv
      integer   gs_handle
      integer   vertex(1)
      real*8    dclock,t0,t00

      common /c_is2/ glo_num1((lx1+2)*(ly1+2)*(lz1+2)*lelt)
      integer*8 glo_num1

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

C     For GPU
      common /gsh_acc/ gsh_face_acc,ids_lgl1(-1:lpts*2)
     $     ,           ids_ptr(-1:lpts*2),ids_lgl2(-1:lpts*2)
      integer gsh_face_acc,ids_lgl1,ids_ptr,ids_lgl2

      common /dum/ n_nonlocal
      integer      n_nonlocal
      
      integer*8 wk(lpts)
      real      u(lpts), v(lpts)
      integer   idx(lpts),ndssum,nglobl,pcount

      t0 = dclock()

c     Global-to-local mapping for gs
      mx = nx+2

      if (nid.eq.0) write(6,*) 'setup_dgds2:: init'

      t00= dclock()           
      call set_vert2(glo_num1,ngv,mx,nel,vertex,.false.) !2nd call
      if (nid.eq.0)  write(6,2)  dclock()-t00
    2 format(' set_vert time::',1pe11.4,' seconds ')

      ntot = (nx1+2)*(ny1+2)*nelt
      if (if3d) ntot = (nx1+2)*(ny1+2)*(nz1+2)*nelt
c
      mz0 = 1
      mz1 = 1

      if (if3d) mz0 = 0
      if (if3d) mz1 = nz1+1
 
c     call icopy8          (glo_num1,glo_num,ntot)
      call iface_vert_int8 (glo_num,glo_num1,mz0,mz1,nelt) 

      nface= 2*ndim
      ntot = nx1*nz1*nface*nelt !total grids on faces

      t00= dclock()
c      call gs_setup(gsh_face,glo_num,ntot,nekcomm,np)

      if (nid.eq.0)  write(6,3)  dclock()-t00
    3 format('   gs_setup time::',1pe11.4,' seconds ')

      
C     added for gs_op_acc
      do ipass = 1,2
         call icopy8(wk,glo_num,ntot)
         call isort8(wk,idx,ntot)
         
         ig = 0
         ip = 0
         glo_num_last = 0
         
         do i = 1,ntot
            if (wk(i).ne.0) then
               ip = ip + 1
               il = idx(i)
               if (ig.eq.0) then
                  ig=ig+1
                  ids_lgl1(ip) = il
                  ids_lgl2(ip) = ig
               elseif( wk(i) .eq. glo_num_last) then
                  ids_lgl1(ip) = il
                  ids_lgl2(ip) = ig
               else
                  ig = ig+1
                  ids_lgl1(ip) = il
                  ids_lgl2(ip) = ig
               endif
               glo_num_last = wk(i)
               wk(i) = abs(wk(i))
            endif
         enddo

         ids_lgl1(0) = ip
         ids_lgl2(0) = ig
         ndssum      = ip
         nglobl      = ig

         if (ipass.eq.1) then ! eliminate singletons
            call gs_setup(gsh_face,glo_num,ntot,nekcomm,mp) 
            call rone(u,ntot)
            call gs_op(gsh_face,u,1,1,0) ! 1 ===> +
            call rone(v,ntot)
            call ldssum(v,wk)

            n_nonlocal = 0
            do i = 1,ndssum   ! Identifiy singletons and nonlocals
               il = ids_lgl1(i)
               ig = ids_lgl2(i)
               if (u(il) .lt.1.1) then        ! Singleton
                  glo_num(il) = 0
               elseif (u(il).ne.v(il)) then   ! Nonlocal
                  n_nonlocal = n_nonlocal+1
                  glo_num(il) = -glo_num(il)
               endif
            enddo
         else ! ipass = 2
            igl = 0
            nnl = 0
            do i = 1,n_nonlocal ! Identify nonlocal in ug()
               ig = ids_lgl2(i)
               if (ig.ne.igl) nnl = nnl+1
               wk(nnl) = wk(i)
               igl = ig
            enddo
            n_nonlocal = nnl
            call gs_setup(gsh_face_acc,wk,n_nonlocal,nekcomm,mp)
         endif
      enddo

      ids_lgl1(-1) = gsh_face
      ids_lgl2(-1) = n_nonlocal

C Used in GPU calculation
      pcount = 1
      ids_ptr(-1) = n_nonlocal
      ids_ptr(0)  = ids_lgl2(0)
      ids_ptr(1)  = 1
      do i=2,ndssum
         if (ids_lgl2(i) .ne. ids_lgl2(i-1)) then
            pcount = pcount+1
            ids_ptr(pcount) = i
         endif
      enddo
      ids_ptr(pcount+1) = ndssum + 1
c
c     Check  
c      write(*,*) pcount, ids_lgl1(0), ids_lgl2(0), n_nonlocal
c      if (pcount .ne. ids_lgl2(0)) stop

      if (nid.eq.0) then
         write(6,5) dclock()-t0,gsh_face,nx,ngv,melg
    5    format(' setup_dsdg2 time::',1pe11.4,' seconds ',2i3,2i12)
      endif
c
      return
      end

#endif 

#endif

