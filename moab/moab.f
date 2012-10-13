c     Nek-MOAB interface notes
c     ------------------------
c Files, repo:
c Most of the Nek-MOAB interface is in the trunk/nek/3rdparty/moab.f source file, which is
c called into from nek's rdparam, readat, and usrchk subroutines.  From 37k ft, nek's datastructure
c stores fields on 3d (nxnxn) arrays on elements, corresponding to a (n-1)th order Guass-Lobatto
c quadrature.  Data is duplicated in those arrays between neighboring elements.
c 
c Global-local ids, etc:
c - in Nek, vertices not really identified as entities, just by coordinates; fields and
c   vertex coordinates are stored in 3d arrays, in lexicographical order, on elements
c - corner vertex positions also held, in separate arrays, in fe order
c - each proc stores global to local processor map (gllnid) and element number (gllel) arrays
c   for element to processor/local id location
c - 2 possible materials (fluid, solid), with fluid elements before solid elements locally
c - side boundary conditions stored in 3d array (elem, side#=6, field), where field is temp or vel
c
c Flow:
c - rdparam (connect2.f): read ifmoab, ifcoup, ifvcoup
c - readat (connect2.f): 
c   . read h5mfle, block/sideset info
c   . bcast h5mfle, block/sideset info
c   . call moab_dat
c     - call nekMOAB_load:
c       . instantiate iMesh, iMeshP, get root set
c       . load file
c       . get tag handles for material, neumann sets
c     - call nekMOAB_get_elems:
c       . foreach block (fluid + solid), init elem iter, elem counts, starts
c       . check against hard-set sizes from SIZE
c       . foreach block (fluid + solid), call nekMOAB_gllnid:
c         - set gllnid(e), global id to proc index
c     - call chk_nel(connect2.f): check local sizes against global hard-set sizes from SIZE
c     - call nekMOAB_create_tags: create results fields tags, set SEM_DIMS on root set
c     - call mapelpr(map2.f): 
c       . check #procs against nelgt, exit if too few elems
c       . call set_proc_map(map2.f)-get_map(map2.f)-get_vert(navier8.f)-nekMOAB_loadConn:
c         - get connectivity in vertex(8,elem) in lex order in terms of global ids
c       . compute, distribute gllel (global to local elem map)
c     - call moab_geometry(xm1,ym1,zm1):
c       . call nekMOAB_loadCoord: for each block (fluid+solid):
c         - get connectivity pointer
c         - get vertex coordinates in blocked-storage array (xxx... yyy... zzz...)
c       . call permute_and_map for each coord array: permute coords, then compute spectral pts
c         using high-order (hex27) vertices
c     - call xml2xc: set corner vertex positions (xc, yc, zc), flip corners from lex to fe ordering
c
c     - call ifill: fill bc array moabbc with -1s
c     - call nekMOAB_BC(moabbc):
c       . for each sideset:
c         - get sideset id, field #
c         - call nekMOAB_intBC(moabbc->bcdata):
c           . get all faces in set recursively; foreach face:
c             - get connectivity of face, hexes adj to face; foreach hex:
c               . get hex connectivity, side # of face wrt hex
c               . call nekMOAB_getElNo: get local hex #, through global ids, gllel
c               . set bcdata(side_no, elem_no, field) to sideset id
c         - end (nekMOAB_intBC)
c     - end (nekMOAB_BC)
c   . end (moab_dat)
c - setlog (bdry.f): print values of ifmoab, ifcoup, ifvcoup
c (solve)
c - usrchk (zero.usr, xxx.usr):
c   . call usr_moab_output (zero.usr, xxx.usr):
c     - call nekMOAB_export_vars (moab.f): write fields to MOAB tags
c     - if an io step, call iMesh_save
c 
#ifdef PTRSIZE8
#define IMESH_HANDLE integer*8
#else
#define IMESH_HANDLE integer*4
#endif

#define IMESH_NULL 0
#define IMESH_ASSERT \
      if (ierr.ne.0) call imesh_err(ierr,imeshh,'moab.f ',__LINE__)

#define IMESH_NULLSTRIP(s) s(:index(s, char(0))-1)

#define MYLOC LOC

c-----------------------------------------------------------------------
      subroutine nekMOAB_init(comm_f, imesh_instance, partn_handle, 
     $     fileset_handle, ierr)

      implicit none
#include "NEKMOAB"
#include "mpif.h"
      IBASE_HANDLE_T imesh_instance, comm_c
      iBase_EntitySetHandle fileset_handle, partn_handle
      integer comm_f, ierr, comm_sz, comm_rk

      if (imesh_instance .eq. 0) then
c     !Initialize imesh and load file
         imeshh = IMESH_NULL
         call iMesh_newMesh(" ", imeshh, ierr)
         IMESH_ASSERT
         iCreatedImesh = 1
      else
         imeshh = imesh_instance
         iCreatedImesh = 0
      endif

#ifdef MPI
      call MPI_Comm_size(comm_f, comm_sz, ierr)
      if (ierr .ne. MPI_SUCCESS) return

      if (comm_sz .gt. 1 .and. partn_handle .eq. 0) then
         call moab_comm_f2c(comm_f, comm_c)
         call iMeshP_createPartitionAll(%VAL(imeshh), 
     1        %VAL(comm_c), hPartn, ierr)
         IMESH_ASSERT
         iCreatedPartn = 1
      else
         hPartn = partn_handle
         iCreatedPartn = 0
      endif
#endif

      if (fileset_handle .ne. 0) fileset = fileset_handle

      return
      end

c-----------------------------------------------------------------------
      subroutine nekMOAB_import
      implicit none
#include "NEKMOAB"
      include 'PARALLEL'
      include 'GEOM'
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer nekcomm, nekgroup, nekreal, nid_, np_

      integer i,ierr,e

      if (imeshh .eq. 0 .or. hPartn .eq. 0 .or.
     $     fileset .eq. 0) then
         call nekMOAB_init(nekcomm, imeshh, hPartn, fileset, 
     $        ierr)
         if (ierr .ne. iBase_SUCCESS) return
      endif
         
      if (fileset .eq. 0) then
         call nekMOAB_load   ! read mesh using MOAB
      endif

#ifdef MPI
      partsSize = 0
      rpParts = IMESH_NULL
      call iMeshP_getLocalParts(%VAL(imeshh), %VAL(hPartn), 
     1     rpParts, partsSize, partsSize, ierr)
#endif

      call nekMOAB_create_tags             ! allocate MOAB tags to reflect Nek variables

      call nekMOAB_get_elems              ! read material sets and establish mapping
      call chk_nel
      call mapelpr2                        ! create gllel mapping 

      if (ifhex) call moab_geometry(xm1,ym1,zm1)     ! fill xm1,ym1,zm1
      if (iftet) call moab_geometry_tetstris (xm1,ym1,zm1)

      call xml2xc                         ! fill xc,yc,zc

c      do e = 1, nelt
c         do i=1,4
c         write(6,*) e,xc(i,e),yc(i,e),zc(i,e)
c         enddo
c      enddo
c      call exitt
 
      call nekMOAB_BC             ! read MOAB BCs 

c     call nekMOAB_compute_diagnostics

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_get_instance(imesh_instance, partn, 
     $     fileset_handle)

      implicit none
#include "NEKMOAB"
#include "mpif.h"
      IBASE_HANDLE_T imesh_instance, partn
      iBase_EntitySetHandle fileset_handle

      imesh_instance = imeshh
      fileset_handle = fileset
      partn = hPartn

      return
      end

c-----------------------------------------------------------------------
      subroutine nekMOAB_create_tags
c
      implicit none
#include "NEKMOAB"

      integer ierr, ntot
      integer semdim(3)
      iBase_TagHandle tagh

      ntot = nx1*ny1*nz1

c tags used for initializing model, should already be there
      globalIdTag = 0
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "GLOBAL_ID",       !/*in*/ const char* tag_name,
     $     globalIdTag,       !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT
      
      matsetTag = 0
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "MATERIAL_SET", !/*in*/ const char* tag_name,
     $     matSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT

      neusetTag = 0
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "NEUMANN_SET", !/*in*/ const char* tag_name,
     $     neuSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT

c create a tag to store SEM dimensions, and set it on the file set
      semdim(1) = nx1
      semdim(2) = ny1
      semdim(3) = nz1
      call nekMOAB_create_find_tag("SEM_DIMS", 
     $     " moab:TAG_STORAGE_TYPE=SPARSE", 3, iBase_INTEGER, tagh,
     $     .true., ierr)
      IMESH_ASSERT
      call iMesh_setEntSetData(%VAL(imeshh), %VAL(fileset), %VAL(tagh), 
     1     semdim, 12, ierr)
      IMESH_ASSERT

c tags for results variables
      call nekMOAB_create_find_tag("SEM_X", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, xm1Tag, .true., ierr)
      IMESH_ASSERT
      call nekMOAB_create_find_tag("SEM_Y", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, ym1Tag, .true., ierr)
      IMESH_ASSERT
      call nekMOAB_create_find_tag("SEM_Z", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, zm1Tag, .true., ierr)
      IMESH_ASSERT

      call nekMOAB_create_find_tag("VX", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, vxTag, .true., ierr)
      IMESH_ASSERT
      call nekMOAB_create_find_tag("VY", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, vyTag, .true., ierr)
      IMESH_ASSERT
      call nekMOAB_create_find_tag("VZ", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, vzTag, .true., ierr)
      IMESH_ASSERT

      call nekMOAB_create_find_tag("TEMP", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $     ntot, iBase_DOUBLE, tTag, .true., ierr)
      IMESH_ASSERT

      if (nx2.eq.nx1 .and. ny2.eq.ny1 .and. nz2.eq.nz1) then
         call nekMOAB_create_find_tag("PRESS", 
     $        " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     $        ntot, iBase_DOUBLE, pTag, .true., ierr)
         IMESH_ASSERT
      endif

c     initialize these tags to zero since their use depends on user input
      vtTag = 0
      vpTag = 0
      vdTag = 0
      vpowTag = 0
      dTag = 0
      powTag = 0
      call nekMOAB_create_find_tag("VTEMP", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", 
     $     1, iBase_DOUBLE, vtTag, .true., ierr)
      IMESH_ASSERT

      if (nx2.eq.nx1 .and. ny2.eq.ny1 .and. nz2.eq.nz1) then
         call nekMOAB_create_find_tag("VPRESS", 
     $       " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", 
     $        1, iBase_DOUBLE, vpTag, .true., ierr)
         IMESH_ASSERT
      endif

c may or may not have these tags, depending on coupler state
      call nekMOAB_create_find_tag("VDENSITY", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", 
     $     1, iBase_DOUBLE, vdTag, .true., ierr)
      IMESH_ASSERT

      call nekMOAB_create_find_tag("VPOWER", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", 
     $     1, iBase_DOUBLE, vpowTag, .true., ierr)
      IMESH_ASSERT

c may or may not have these tags, depending on coupler state
      call nekMOAB_create_find_tag("DENSITY", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", 
     $     ntot, iBase_DOUBLE, dTag, .true., ierr)
      IMESH_ASSERT

      call nekMOAB_create_find_tag("POWER", 
     $     " moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=-1.0", 
     $     ntot, iBase_DOUBLE, powTag, .true., ierr)
      IMESH_ASSERT

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_create_find_tag(tag_name, tag_opts, tag_size, 
     $     tag_datatype, tag_handle, create_if_missing, ierr)
c attempt to get a tag handle of a specified name; if not found, create one; if found
c verify tag characteristics
      implicit none
#include "NEKMOAB"
      include 'mpif.h'

      character*(*) tag_name, tag_opts
      integer tag_size, tag_datatype, ierr, dum_int
      logical create_if_missing
      iBase_TagHandle tag_handle

      call iMesh_getTagHandle(%VAL(imeshh), tag_name, tag_handle, 
     $     ierr)
      if (iBase_SUCCESS .ne. ierr .and. create_if_missing) then
c need to create it
         call iMesh_createTagWithOptions(%VAL(imeshh), tag_name,
     $        tag_opts, %VAL(tag_size), %VAL(tag_datatype), 
     $        tag_handle, ierr)
         IMESH_ASSERT
      else if (create_if_missing) then
c verify characteristics: size, datatype
         call iMesh_getTagSizeValues(%VAL(imeshh), %VAL(tag_handle), 
     $        dum_int, ierr)
         if (ierr .ne. iBase_SUCCESS .or. dum_int .ne. tag_size) then
            ierr = iBase_INVALID_TAG_HANDLE
            return
         endif
         call iMesh_getTagType(%VAL(imeshh), %VAL(tag_handle), 
     $        dum_int, ierr)
        if (ierr .ne. iBase_SUCCESS .or. dum_int .ne. tag_datatype) then
           ierr = iBase_INVALID_TAG_HANDLE
           return
        endif
      endif

      return
      end
c
c-----------------------------------------------------------------------
      subroutine nekMOAB_load
c
c     Load "filename" into imesh/moab, store imesh handle 
c     in /nekmoab/ common block
c
      implicit none
#include "NEKMOAB"
      include 'mpif.h'
c two forms of load options, depending whether we\'re running serial or parallel      
      character*(*) parLoadOpt, serLoadOpt
      parameter(parLoadOpt=" moab:PARALLEL=READ_PART   moab:PARTITION=PA
     $RALLEL_PARTITION moab:PARALLEL_RESOLVE_SHARED_ENTS moab:PARTITION_
     $DISTRIBUTE moab:CPUTIME")
      parameter(serLoadOpt = " ")
c      parameter(saveOut = "meshdata.txt")
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal
      integer nekcomm, nekgroup, nekreal, nid_, np

      integer ierr
      IBASE_HANDLE_T ccomm

c     create a file set to load into
      call iMesh_createEntSet(%VAL(imeshh), %VAL(0), fileset, ierr)
      IMESH_ASSERT

#ifdef MPI
      if (np .gt. 1) then
         call iMeshP_loadAll(%VAL(imeshh), %VAL(hPartn),%VAL(fileset),
     $        H5MFLE, parLoadOpt, ierr)
         IMESH_ASSERT

c      if  (nid .eq. 0) then         
c         call iMeshP_saveAll(%VAL(imeshh), %VAL(hPartn),%VAL(fileset),
c     $        saveOut," ", ierr,12,1)
c         IMESH_ASSERT
c      endif

      else
#endif
         call iMesh_load(%VAL(imeshh), %VAL(fileset),
     $        H5MFLE, serLoadOpt, ierr)
         IMESH_ASSERT
#ifdef MPI
      endif
#endif
c      write(6,*) 'hPartn :', hPartn 
c      if (nid .eq. 0) then
c         call iMeshP_saveAll(%VAL(imeshh), %VAL(hPartn),%VAL(fileset),
c     $        saveOut," ", ierr)
c      endif
      return
      end  
c-----------------------------------------------------------------------
      subroutine nekMOAB_get_elems()
c
c get fluid/solid elements and establish global id mapping
c
      implicit none
#include "NEKMOAB"
      include 'PARALLEL'
      include 'SCRCT'
      include 'ZPER'

      include 'mpif.h'

      integer iglsum, i, ierr
      iBase_EntitySetHandle dumsets(numsts)
      IBASE_HANDLE_T valptr, setsptr
      integer dumval, dumalloc, dumsize, dumnum, ns, ilast, atend, 
     $        count, tmpcount, ietmp1(numsts), ietmp2(numsts), npass, 
     $        ipass, m, k, iwork
      common /ctmp0/ iwork(lelt)
      common /nekmpi/  nid_,np_,nekcomm,nekgroup,nekreal
      integer          nid_,np_,nekcomm,nekgroup,nekreal

      integer iBase_DIM,iMesh_STRUC

c assign mesh structure and dimension 
      if (if3d) then
         if (IFHEX) iMesh_STRUC= iMesh_HEXAHEDRON   !=7 
         if (IFTET) iMesh_STRUC= iMesh_TETRAHEDRON  !=6 
         iBase_DIM= ndim 
      else
         if (IFHEX) iMesh_STRUC= iMesh_QUADRILATERAL!=4   
         if (IFTET) iMesh_STRUC= iMesh_TRIANGLE     !=3 
         iBase_DIM= ndim 
      endif
      if (nid.eq.0) 
     $  write(6,*) 'call nekMOAB_get_elems: iMesh_Struc=',iMesh_STRUC
        write(6,*) 'call nekMOAB_get_elems: iBase_DIM  =',iBase_DIM

c get fluid, other material sets, and count elements in them
      valptr  = loc(dumval)
      setsptr = loc(dumsets(1))
      dumalloc= numsts
      ilast= 0
      do i = 1, numflu+numoth
         dumval = matids(i)
         dumsize= numsts
c get the set by matset number
         call iMesh_getEntSetsByTagsRec(%VAL(imeshh), %VAL(fileset),
     $        matsetTag, valptr, %VAL(1), %VAL(1), 
     $        setsptr, dumsize, dumsize, ierr)
         if (dumsize .gt. 1) then
            call exitti('More than one material set with id ', dumval)
         endif
         IMESH_ASSERT
c get the number of hexes
         if (dumsize .gt. 0) then
            call iMesh_getNumOfTopoRec(%VAL(imeshh), %VAL(dumsets(1)), 
     $           %VAL(iMesh_STRUC), %VAL(1), dumnum, ierr)
            IMESH_ASSERT
            matsets(i) = dumsets(1)
            iestart(i) = ilast + 1
            ilast = ilast + dumnum
            iecount(i) = dumnum
c     get an iterator for this set, used later
            call iMesh_initEntArrIterRec(%VAL(imeshh), %VAL(dumsets(1)),
     $           %VAL(iBase_DIM), %VAL(iMesh_STRUC),
     $           %VAL(dumnum), %VAL(0), %VAL(1), ieiter(i), ierr)
         else
            matsets(i) = 0
            iestart(i) = ilast + 1
            iecount(i) = 0
         endif

c set total number if nec
         if (i .eq. numflu) then
            nelv = ilast
         endif
c this is if, not elseif, to handle numoth=0
         if (i .eq. numflu+numoth) then
            nelt = ilast
         endif
      enddo

c set remaining values to default values
      do i = numflu+numoth+1, numsts
         iecount(i)= -1
         iestart(i)= -1
         ieiter(i) = 0
         matsets(i)= 0
      enddo

c check local size
      if(nelt .le. 0 .or. nelv .gt. lelv) then
         print *, 'ABORT: nelv is invalid in nekmoab_proc_map'
         print *, 'nelv, lelv = ', nelv, lelv
         call exitt
      endif

c reduce to get global numbers of fluid, other elements, and check size
      nelgv = iglsum(nelv,1)
      nelgt = iglsum(nelt,1)
      if (NELGT .gt. LELG) then
         print *, 'ABORT; increase lelg ',nelgv,lelg
         call exitt
      endif

c do global scan to allow computation of global element ids
      ietmp1(1) = 0
      do i = 1, numflu
         ietmp1(1) = ietmp1(1) + iecount(i)
      enddo
      ietmp1(2) = 0
      do i = numflu+1, numflu+numoth
         ietmp1(2) = ietmp1(2) + iecount(i)
      enddo
      call mpi_scan(ietmp1,ietmp2,numflu+numoth,MPI_INTEGER,MPI_SUM,
     $     nekcomm, ierr)
      if (ierr .ne. MPI_SUCCESS) ierr = iBase_FAILURE
      IMESH_ASSERT
c     set returned nums to exclusive start - 1
      ietmp2(1) = ietmp2(1) - ietmp1(1)
      ietmp2(2) = ietmp2(2) - ietmp1(2)
c     set gids for local fluid, other elems
      call izero(lglel, nelt)
      call izero(gllel, nelgt)
      call izero(gllnid,nelgt)

      do i = 1, nelv
         lglel(i) = ietmp2(1)+i
         gllel(lglel(i)) = i
         gllnid(lglel(i)) = nid_
      enddo
      do i = nelv+1, nelt
         lglel(i) = nelgv+ietmp2(2)-nelv+i
         gllel(lglel(i)) = i
         gllnid(lglel(i))= nid_
      enddo
c     now communicate to other procs (taken from set_proc_map in map2.f)
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllel(k),iwork,'+  ',m)
         k = k+m
      enddo

      call igop(GLLNID, iwork, 'M  ', NELGT)

c     set the global id tag on elements

      do i = 1, numflu+numoth
         atend = 0
         count = 0
         if (iecount(i) .ne. 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           ierr)
            IMESH_ASSERT
         else
            atend = 1
         endif

         do while (atend .eq. 0)
c use the same iterator for all variables, since the elems are the same
            call nekMOAB_set_int_tag(ieiter(i), globalIdTag, 1, 
     $           tmpcount, lglel(iestart(i)+count))
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           %VAL(tmpcount), atend, ierr)
            IMESH_ASSERT
            count = count + tmpcount
         enddo
      enddo

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_gllnid(matset, iter, count)
      implicit none
c
c initialize global ids for hexes in this set
      iBase_EntitySetHandle matset
      iBase_EntityArrIterator iter
      integer itmp, i, j, ierr, gid, count, atend
      IBASE_HANDLE_T tag_ptr
      pointer(tag_ptr, gid(1))

#include "NEKMOAB"
      include 'PARALLEL'
c 
      call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(iter), ierr)
      IMESH_ASSERT
      i = 0
      atend = 0
      do while (atend .eq. 0)
         call iMesh_tagIterate(%VAL(imeshh), %VAL(globalIdTag), 
     $        %VAL(iter), tag_ptr, itmp, ierr)
         IMESH_ASSERT

c set the global ids to this proc
         do j = 1, itmp
            if (gid(j) .gt. nelgt) 
     $           call exitti('Global id greater than NELGT', gid(i))
            gllnid(gid(j)) = nid
         enddo
c step the iterator
         call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(iter), %VAL(itmp), 
     $        atend, ierr)
         IMESH_ASSERT
         i = i + itmp
      enddo

c assert we got the right number of elements
      if (i .lt. count) call exitti(
     $     'Wrong number of entities in region iterator', i)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_loadConn(vertex, nelgt, ncrnr)
c
c     Fill the vertex array with connectivity data from imesh
c
c     vertex(ncrnr, nelt): int array, global id of vertices in element nelt
c     nelgt: int, global number of elements
c     ncrnr: int, number of corner vertices per element (should be 8)
c
      implicit none
#include "NEKMOAB"
      integer ncrnr,nelgt
      integer vertex(ncrnr, *), i

c get corner vertex gids
      integer e_in_set, eid, j, k, nv, ierr, e_in_chunk, v_per_e,e
      integer gids (27)  ! FIXME 7/27/2012:  27 will hold general case?
      iBase_EntityArrIterator iter
      IBASE_HANDLE_T connect_ptr
      iBase_EntityHandle connect
      pointer (connect_ptr, connect(0:1))
c      integer   indv  (4,lelt)
      integer   indv  (4)
      data    indv / 1, 2, 3, 4/              

      integer l2c(8),l2c_tet(4)
      save    l2c,l2c_tet
      data    l2c     / 1, 2, 4, 3, 5, 6, 8, 7 / !hexahedral
      data    l2c_tet / 1, 2, 3, 4/              !tetrahedron

c     assign variable for each mesh structures
      if (if3d) then
          if (ifhex) nv = 8
          if (iftet) nv = 4
      else
          if (ifhex) nv = 4
          if (iftet) nv = 3
      endif
      write(6,*) 'nv is: ',nv
c     assign corner verticies from moab
      do i = 1, numflu+numoth
         e_in_set = 0
         eid = iestart(i)
         do while (e_in_set .lt. iecount(i))
            if (e_in_set .eq. 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           ierr)
            IMESH_ASSERT
            endif

c     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), 
     $           connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

c     for each element
            do j = 0, e_in_chunk-1
c     get vertex gids for this e
               call iMesh_getIntArrData(%VAL(imeshh), !iMesh_Instance instance,
     $              connect(j*v_per_e), %VAL(nv), %VAL(globalIdTag), 
     $              loc(gids), nv, nv, ierr)
               IMESH_ASSERT
c     permute into vertex array
               if (ifhex) then
                  do k=1, nv
                  vertex(k, eid) = gids(l2c(k))
                  enddo
               else
                  do k=1, nv
                  vertex(k, eid) = gids(l2c_tet(k))
c                  write(6,*) 'vertex:',eid,vertex(k,eid) 
                  enddo
               endif
               eid = eid + 1
            enddo

            e_in_set = e_in_set + e_in_chunk
         enddo
      enddo

c ------TESTING----------
      do e=1,nelgt
      write(6,*) 'b4moab', vertex(1,e),vertex(2,e),vertex(3,e),
     $                     vertex(4,e)
         call isort(vertex(1,e),indv,4)
      write(6,*) 'afmoab', vertex(1,e),vertex(2,e),vertex(3,e),
     $                     vertex(4,e)
      write(6,*) 'indv', indv(1),indv(2),indv(3),
     $                   indv(4)
      enddo
c----TESTING-------------
      return
      end 
c-----------------------------------------------------------------------------
      subroutine nekMOAB_loadCoord(xmc, ymc, zmc, lc)
c     
c     stuff the xyz coords of the 27 verts of each local element -- 
c     shared vertex coords are stored redundantly
c     
      implicit none
#include "NEKMOAB"
      integer lc    
      real    xmc(lc,*), ymc(lc,*), zmc(lc,*)
      IBASE_HANDLE_T connect_ptr
      iBase_EntityHandle connect
      pointer(connect_ptr, connect(0:1))
      integer i, j, k, ierr, e_in_chunk, e_in_set, v_per_e, e_tot

      e_tot = 0
      do i = 1, numflu+numoth
         e_in_set = 0
         do while (e_in_set .lt. iecount(i))
            if (e_in_set .eq. 0) then
               call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $              ierr)
               IMESH_ASSERT
            endif

c     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), 
     $           connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

c     for each element
            do j = 0, e_in_chunk-1
c     get vertex gids for this e
               do k = 1, lc !lc=TWENTYSEVEN,TEN,NINE,SIX
                  call iMesh_getVtxCoord(%VAL(imeshh),
     $                 %VAL(connect(j*v_per_e+k-1)), 
     $                 xmc(k,e_tot+j+1), ymc(k,e_tot+j+1), 
     $                 zmc(k,e_tot+j+1), ierr)
             IMESH_ASSERT
               enddo
            enddo
            e_tot = e_tot + e_in_chunk
            e_in_set = e_in_set + e_in_chunk
         enddo
      enddo

      return
      end 
c-----------------------------------------------------------------------------
      subroutine nekMOAB_BC()
c
c     Mark the "side" (hex/face) with the index of the corresponding set in ibcsts array
c

      implicit none
#include "NEKMOAB"
      common /mbc/ moabbc(6,lelt,ldimt1) 
      integer moabbc

      IBASE_HANDLE_T hentSet(*)
      pointer (rpentSet, hentSet)
      integer entSetSize

      integer ierr, i, j, set_idx, set_ids(numsts)

      call ifill(moabbc, -1, 6*lelt*ldimt1)

      !Sidesets in cubit come in as entity sets with the NEUMANN_SET -- see sample file
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "NEUMANN_SET", neuSetTag, ierr)
      IMESH_ASSERT

      rpentSet = IMESH_NULL
      entSetSize      = 0
      call iMesh_getEntSetsByTagsRec(%VAL(imeshh),
     $     %VAL(fileset), neuSetTag, %VAL(IMESH_NULL), %VAL(1), %VAL(0),
     $     rpentSet, entSetSize, entSetSize, ierr)
      IMESH_ASSERT

      !print *, 'H3', entSetSize
c     get the set ids
      if (entSetSize .gt. numsts) then
c     too many sets, need to bail
         ierr = iBase_FAILURE
         IMESH_ASSERT
      endif
      do i = 1, entSetSize
         call iMesh_getEntSetIntData(%VAL(imeshh),
     $     %VAL(hentSet(i)), %VAL(neuSetTag), set_ids(i), ierr)
         IMESH_ASSERT
      enddo

c     loop over local neusets, will be <= numsts
      do i=1, entSetSize
c     find the right set
         set_idx = -1
         do j = 1, numsts
            if (set_ids(i) .eq. ibcsts(j)) then
               set_idx = j
               if (set_idx .ne. -1) then
                  call nekMOAB_intBC(moabbc, hentSet(i), set_ids(i), 
     $                 bcf(set_idx))
               endif
            endif
         enddo
      enddo

      call free(rpentSet)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_intBC(bcdata, setHandle, setId, field)
c
c     Internal function, don\'t call directly
c

      implicit none
#include "NEKMOAB"

      integer bcdata(6,lelt,ldimt1), field

      iBase_EntitySetHandle setHandle

      IBASE_HANDLE_T hset
      integer setId !coming from cubit

      IBASE_HANDLE_T faces(*)
      pointer (rpfaces, faces)
      integer facesSize

      IBASE_HANDLE_T ahex(*)
      pointer (rpahex, ahex)
      integer ahexSize

      IBASE_HANDLE_T fvtx(*)
      pointer (rpfvtx, fvtx)
      integer fvtxSize

      IBASE_HANDLE_T hvtx(*)
      pointer (rphvtx, hvtx)
      integer hvtxSize

      integer ierr, i, j, elno

      integer side_no, side_offset, side_sense
      integer hexMBCNType
      integer num_sides
      integer elnos(2)
      integer numids

      integer iMesh_FACE_STRUC,iBase_DIM,iMesh_FACE_VERTNUM
      integer iMesh_FACE_SHAPE,iMesh_STRUC

c assign mesh structure and dimension 
      if (if3d) then

         if (ifhex) iMesh_STRUC   = iMesh_HEXAHEDRON   
         if (ifhex) iMesh_FACE_STRUC  = iMesh_QUADRILATERAL
         if (ifhex) iMesh_FACE_VERTNUM = 4

         if (iftet) iMesh_STRUC   = iMesh_TETRAHEDRON   
         if (iftet) iMesh_FACE_STRUC= iMesh_TRIANGLE      
         if (iftet) iMesh_FACE_VERTNUM = 3

         iMesh_FACE_SHAPE= iBase_FACE
         iBase_DIM= ndim
      else
         if (ifhex) iMesh_STRUC   = iMesh_QUADRILATERAL
         if (ifhex) iMesh_FACE_STRUC= iMesh_LINE_SEGMENT 
         if (ifhex) iMesh_FACE_VERTNUM =2 

         if (iftet) iMesh_STRUC   = iMesh_TRIANGLE     
         if (iftet) iMesh_FACE_STRUC= iMesh_LINE_SEGMENT  
         if (iftet) iMesh_FACE_VERTNUM = 2

         iMesh_FACE_SHAPE= iBase_EDGE
         iBase_DIM= ndim
      endif
      if (nid.eq.0)
     $  write(6,*) 'call nekMOAB_get_elems: iMesh_Struc=',iMesh_STRUC


      numids = 0

      call iMesh_MBCNType(%VAL(iMesh_STRUC), hexMBCNType)

      !what faces are in this set?
      facesSize = 0
      rpfaces = IMESH_NULL
      call iMesh_getEntitiesRec(%VAL(imeshh), 
     $  %VAL(setHandle),%VAL(iMesh_FACE_SHAPE),%VAL(iMesh_FACE_STRUC),
     $  %VAL(1), rpfaces, facesSize, facesSize, ierr)
      IMESH_ASSERT

      num_sides = 0

      do i=1, facesSize
         !get vertices defining the face
         fvtxSize = 0
         rpfvtx = IMESH_NULL
         call iMesh_getEntAdj(%VAL(imeshh), %VAL(faces(i)), 
     $        %VAL(iBase_VERTEX), rpfvtx, fvtxSize, fvtxSize, ierr)
         IMESH_ASSERT         

         !get hexes adjacent to the face (1 or 2) 
         ! -- hopefully only returns hexes on local proc, but untested
         ahexSize = 0
         rpahex = IMESH_NULL                  
         call iMesh_getEntAdj(%VAL(imeshh), %VAL(faces(i)), 
     $        %VAL(iBase_DIM), rpahex, ahexSize, ahexSize, ierr)
         IMESH_ASSERT

         do j=1, ahexSize                                
            !get verts adjacent to the hex
            hvtxSize = 0
            rphvtx = IMESH_NULL                  
            call iMesh_getEntAdj(%VAL(imeshh), 
     $           %VAL(ahex(j)), %VAL(iBase_VERTEX), 
     $           rphvtx, hvtxSize, hvtxSize, ierr)
            IMESH_ASSERT

            !Get the side number
            call MBCN_SideNumberUlong(%VAL(rphvtx),
     $           %VAL(hexMBCNType), %VAL(rpfvtx), 
     $           %VAL(iMesh_FACE_VERTNUM), %VAL(iBase_DIM-1),
     $           side_no, side_sense, side_offset) 
           if ((side_no .lt. 0) .or. (side_no .gt. 5)) ierr = 1
           IMESH_ASSERT

           side_no = side_no + 1 !moab side number is zero based
           num_sides = num_sides + 1

           !call nekMOAB_getGlobElNo(ahex(j), elno)
           !print *, 'BLAH', elno, side_no, setId

           call nekMOAB_getElNo(ahex(j), elno)
           if (ahexSize .eq. 2) elnos(j) = elno

           if (bcdata(side_no, elno, field) .ne. -1) 
     $          print *, 'Warning: resetting BC, bcno, elno, sideno = ',
     $            setId, elno, side_no 
           bcdata(side_no, elno, field) = setId
           numids = numids + 1

           call free(rphvtx)
         enddo

         if (ahexSize .eq. 2) 
     $        print *, 'Warning: face shared by 2 hexes: ', elnos(1), 
     $        elnos(2)

         call free(rpahex)
         call free(rpfvtx)

      enddo

      call free(rpfaces)

      if (nid.eq.0) write(6,*) 'Setid, numids = ', nid,setId, numids

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_getElNo(handle, elno)
c
c     Given the imesh handle of an element (hex), 
c     return its local nek element number in elno
c     Cannot be used until the GLLEL array has been set (in map2.f)
c
      implicit none
#include "NEKMOAB"
      include 'PARALLEL'

      IBASE_HANDLE_T handle
      integer elno

      !first get global id
      call nekMOAB_getGlobElNo(handle, elno)
      !then convert to local id
      elno = GLLEL(elno)

      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_getGlobElNo(handle, elno)

      implicit none
#include "NEKMOAB"
      include 'PARALLEL'

      IBASE_HANDLE_T handle
      integer elno
      
      integer ierr

      call iMesh_getIntData(%VAL(imeshh), %VAL(handle), 
     $     %VAL(globalIdTag), elno, ierr)
      IMESH_ASSERT
      
      if (GLLNID(elno) .ne. NID) then
         print *, 'Wrong proc for element; gid, proc, gllnid = ',
     $        elno, nid, gllnid(elno)
         call exitt
      endif
                          
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_bcs() ! fill the nek cbc arrays
      implicit none
c
#include "NEKMOAB"      
      common /mbc/ moabbc(6,lelt,ldimt1)
      integer moabbc

      integer e,f,l, j
      character*3 cbi

      integer ibcs(3), i, lcbc, nface
      data ibcs / 0, 0, 0 /

c------------------------------------------------------
c-----TESTING in the case of a non-uniform mesh-------- 
c      open(unit=13,file='bcdata0.txt',status='unknown')
c      open(unit=14,file='bcdata1.txt',status='unknown')
c      open(unit=15,file='bcdata2.txt',status='unknown')
c      open(unit=16,file='bcdata3.txt',status='unknown')
c------------------------------------------------------

      lcbc=18*lelt*(ldimt1 + 1)
      call blank(cbc,lcbc)

      nface = 2*ndim
      do l=2,2      ! For nek5: l=1,nfield
         do e=1,nelt
            do f=1,nface
               if (moabbc(f,e,l) .ne. -1) then
c     moabbc stores local set index, retrieve the character bc type
                  do j = 1, numsts
                     if (moabbc(f,e,l) .eq. ibcsts(j)
     $                    .and. bcf(j) .eq. l) then
                        cbc(f, e, l) = bctyps(j)
c----------------------TESTING---------------------------
c                        if (nid .eq. 0) then
c                        write(13,*) 'BC data:', nid,f,e,l,
c     $                              cbc(f,e,l)
c                        elseif (nid .eq. 1) then
c                        write(14,*) 'BC data:', nid,f,e,l,
c     $                              cbc(f,e,l)
c                       elseif (nid .eq. 2) then
c                        write(15,*) 'BC data:', nid,f,e,l,
c     $                              cbc(f,e,l)
c                        elseif (nid .eq. 3) then
c                        write(16,*) 'BC data:', nid,f,e,l,
c     $                              cbc(f,e,l)
c                        endif
c--------------------END TESTING---------------------------
                     endif
                  enddo
               else
                  cbc(f, e, l) = 'E  '
               endif
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine moab_geometry (xmlo,ymlo,zmlo)

      implicit none
#include "NEKMOAB"      

      real      xmlo(lx1*ly1*lz1,*)
     $        , ymlo(lx1*ly1*lz1,*)
     $        , zmlo(lx1*ly1*lz1,*)
      integer   e, nmoab,i

      integer    lc, nface
      parameter (lc= 3**ldim)
      common /tcrmg/ xcm(lc,lelt), ycm(lc,lelt), zcm(lc,lelt)
      real           xcm,ycm,zcm  
       
c     get coords from moab 
      call nekMOAB_loadCoord(xcm, ycm, zcm, lc) 

c     if (ifhex) nface= 2*ndim
c     if (iftet) nface= ndim+1
c     call outmat(xmc,mc,nface,'xmcdat',nelt)
c     call outmat(ymc,mc,nface,'ymcdat',nelt)
c     call outmat(zmc,mc,nface,'zmcdat',nelt)

c     nmoab = 3     ! not used, we replaced nmoab by mc

c     Interpolate for each element
      do e=1,nelt   
         call permute_and_map(xmlo(1,e),xcm(1,e),nx1,nmoab,e)
         call permute_and_map(ymlo(1,e),ycm(1,e),ny1,nmoab,e)
         if (if3d) call permute_and_map(zmlo(1,e),zcm(1,e),nz1,nmoab,e)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine xml2xc
      
      implicit none
#include "NEKMOAB"
      include 'GEOM'
      integer i, j, k, l,zc3,zc4
      real    dtmp

      integer e
      if (ifhex) then
      
      if (if3d) then
         do e=1,nelt
            l=0
         do k=1,nz1,nz1-1
         do j=1,ny1,ny1-1
         do i=1,nx1,nx1-1
            l = l+1
            xc(l,e) = xm1(i,j,k,e)
            yc(l,e) = ym1(i,j,k,e)
            zc(l,e) = zm1(i,j,k,e)
         enddo
         enddo
         enddo
         enddo

      else  ! 2D
         do e=1,nelt
            l=0
         do j=1,ny1,ny1-1
         do i=1,nx1,nx1-1
            l = l+1
            xc(l,e) = xm1(i,j,1,e)
            yc(l,e) = ym1(i,j,1,e)
c           if (nid.eq.0) write(6,*) e,l,' : ',xc(l,e),yc(l,e)

         enddo
         enddo
         l = 0
         enddo
      endif

      do e=1,nelt ! flip corners back to pre-proc notation
         dtmp = xc(3,e)
         xc(3,e) = xc(4,e)
         xc(4,e) = dtmp
         dtmp = yc(3,e)
         yc(3,e) = yc(4,e)
         yc(4,e) = dtmp
         dtmp = zc(3,e)
         zc(3,e) = zc(4,e)
         zc(4,e) = dtmp

         if(if3d) then
           dtmp = xc(7,e)
           xc(7,e) = xc(8,e)
           xc(8,e) = dtmp
           dtmp = yc(7,e)
           yc(7,e) = yc(8,e)
           yc(8,e) = dtmp
           dtmp = zc(7,e)
           zc(7,e) = zc(8,e)
           zc(8,e) = dtmp
         endif
      enddo
      endif


      if (iftet) then
      zc3 = (nx1)*(nx1+1)/2      
      zc4 = (nx1)*(nx1+1)*(nx1+2)/6
      if (if3d) then !needs work for tets!

         do e=1,nelt
            write(6,*) 'test for tets!!!!!!!'

            xc(1,e) = xm1(1,1,1,e)
            yc(1,e) = ym1(1,1,1,e)
            zc(1,e) = zm1(1,1,1,e)

            xc(2,e) = xm1(nx1,1,1,e)
            yc(2,e) = ym1(nx1,1,1,e)
            zc(2,e) = zm1(nx1,1,1,e)

            xc(3,e) = xm1(zc3,1,1,e)
            yc(3,e) = ym1(zc3,1,1,e)
            zc(3,e) = zm1(zc3,1,1,e)
                        
            xc(4,e) = xm1(zc4,1,1,e)
            yc(4,e) = ym1(zc4,1,1,e)
            zc(4,e) = zm1(zc4,1,1,e)

         enddo

      else  ! 2D
         do e=1,nelt

            xc(1,e) = xm1(1,1,1,e)
            xc(2,e) = xm1(nx1,1,1,e)
            xc(3,e) = xm1(zc3,1,1,e) 
            
            yc(1,e) = ym1(1,1,1,e)
            yc(2,e) = ym1(nx1,1,1,e)
            yc(3,e) = ym1(zc3,1,1,e)
            

c            write(6,*) e,l,' : ',xc(3,e),yc(3,e)

         enddo
      endif

      endif
   
      return
      end
c-----------------------------------------------------------------------
      subroutine permute_and_map(x,xcm,nx,nmoab,e)
c
      implicit none
#include "NEKMOAB"
      integer nx, nmoab, e, ncoord

      real x(1),xcm(0:1)

      common /ctmp0/ z3(3),zpt(lx1)
     $             , xt(3,3,3)             
     $             , wk(3*lx1*ly1*lz1)
     $             , xw(3*lx1*ly1*lz1)
      real z3, zpt, xt, wk, xw

      real interp(lx1*3),interpt(lx1*3),zl(lx1*3)
      save interp,interpt

      integer moabmap (27) 
      integer moabmap1(27),moabmap2(9),moabmap3(10),moabmap4(6)
      save    moabmap1,moabmap2,moabmap3,moabmap4
      data    moabmap1
     $                 /  0,  8,  1, 11, 24,  9,  3, 10,  2  
     $                 , 12, 20, 13, 23, 26, 21, 15, 22, 14 
     $                 ,  4, 16,  5, 19, 25, 17,  7, 18,  6/ !hexahedron 27
      data    moabmap2 /  0, 4, 1, 7, 8, 5, 3, 6, 2 /        !quadrilateral 9
      data    moabmap3 /  0, 4, 1, 6, 5, 2, 7, 8, 9, 3 /     !tetrahedron 10
      data    moabmap4 /  0, 3, 1, 5, 4, 2 /                 !triangle 6

      integer nmlast,nxlast
      save    nmlast,nxlast
      data    nmlast,nxlast / 0,0 /

      real gh_edge (3,3,3),gh_vtx(3,3,3),zgh(3)
      save zgh
      data zgh / -1., 0., 1. /
      integer gh_type, i, j, ldw

      if (nx.ne.nxlast) then ! define interp. op
          nxlast = nx
          call zwgll (zl,interp,nx)    !interp is dummy here
          call igllm (interp,interpt,zgh,zl,3,nx,3,nx)
      endif

      if (ifhex) ncoord=3**ndim
      if (iftet) ncoord=(ndim+1)*(ndim+2)/2

      if (if3d) then
          if (ifhex) call copy(moabmap,moabmap1,ncoord)
          if (iftet) call copy(moabmap,moabmap3,ncoord)
      else
          if (ifhex) call copy(moabmap,moabmap2,ncoord)
          if (iftet) call copy(moabmap,moabmap4,ncoord)
      endif

      do i=1,ncoord      ! currently support only 3x3x3 in moab
         j = moabmap(i)
         xt(i,1,1)= xcm(j)
      enddo

      if (ifhex) then

          gh_type = 2 ! edge extension
          call gh_face_extend(xt,zgh,3,gh_type,gh_edge,gh_vtx)

c         Interpolate from 3x3x3 to (nx1 x ny1 x nz1) SEM mesh
          ldw = 3*lx1*ly1*lz1
          call map_to_crs(x,nx1,xt,3,if3d,wk,ldw)

      else

        !saumil's temporary....

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_export_vars
      implicit none
#include "NEKMOAB"
      include 'GEOM'
      include 'SOLN'
      integer i, j, ierr, ntot, tmpcount, count, atend
      real    vtrans,tag_ptr(1)

      ntot = nx1*ny1*nz1
      do i = 1, numflu+numoth
         atend = 0
         count = 0
         if (iecount(i) .ne. 0) then
            call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           ierr)
            IMESH_ASSERT
         else
            atend = 1
         endif

         do while (atend .eq. 0)
c     use the same iterator for all variables, since the elems are the same
            if (xm1Tag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), xm1Tag, ntot, tmpcount, 
     $           xm1)
            if (ym1Tag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), ym1Tag, ntot, tmpcount, 
     $           ym1)
            if (zm1Tag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), zm1Tag, ntot, tmpcount, 
     $           zm1)

            if (vxTag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), vxTag, ntot, tmpcount, 
     $           vx)
            if (vyTag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), vyTag, ntot, tmpcount, 
     $           vy)
            if (vzTag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), vzTag, ntot, tmpcount, 
     $           vz)

            if (tTag .ne. 0)
     $           call nekMOAB_set_tag(ieiter(i), tTag, ntot, tmpcount, 
     $           t)

            if (tTag .ne. 0 .and. nx2.eq.nx1 .and. ny2.eq.ny1 .and. 
     $           nz2.eq.nz1) 
     $           call nekMOAB_set_tag(ieiter(i), pTag, ntot, tmpcount, 
     $           pr)

            if (dTag .ne. 0) 
     $           call nekMOAB_set_tag(ieiter(i), dTag, ntot, tmpcount, 
     $           vtrans)

c     vertex-based temperature
            if (vtTag .ne. 0)
     $           call nekMOAB_set_vertex_tag(ieiter(i), vtTag, tmpcount, 
     $           t)
c     vertex-based density
            if (vdTag .ne. 0)
     $           call nekMOAB_set_vertex_tag(ieiter(i), vdTag, tmpcount, 
     $           vtrans)
c     vertex-based pressure, but only if its there
            if (vpTag .ne. 0 .and.
     $           nx2.eq.nx1 .and. ny2.eq.ny1 .and. nz2.eq.nz1) 
     $           call nekMOAB_set_vertex_tag(ieiter(i), vpTag, 
     $           tmpcount, pr)

c     step the iterator
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           %VAL(tmpcount), atend, ierr)
            IMESH_ASSERT

            count = count + tmpcount
         enddo

c     double-check the total number of elements in this set
         if (count .ne. iecount(i)) then
            call exitti('Wrong no of elems iterating over matset ', 
     $           matids(i))
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_import_vars(tagh, is_v, field)
      implicit none
#include "NEKMOAB"

      real field(lx1,ly1,lz1,lelt)

      include 'GEOM'
      include 'SOLN'
      integer i, j, ierr, ntot, tmpct, count, atend, is_v
      iBase_TagHandle tagh
      real tag_ptr(1)

      ntot = nx1*ny1*nz1

      if (tagh .eq. 0) return

      do i = 1, numflu+numoth
         call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), ierr)
         IMESH_ASSERT

         atend = 0
         count = 0
         do while (atend .eq. 0)
c use the same iterator for all variables, since the elems are the same
            if (is_v .eq. 1) then
               call nekMOAB_get_vertex_tag(ieiter(i), tagh, tmpct, 
     $              field, ierr)
            else
               call nekMOAB_get_tag(ieiter(i), tagh, ntot, tmpct, 
     $              field, ierr)
            endif

c     step the iterator
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           %VAL(tmpct), atend, ierr)
            IMESH_ASSERT

            count = count + tmpct
         enddo

c double-check the total number of elements in this set
         if (count .ne. iecount(i)) then
            call exitti('Wrong no of elems iterating over matset ', 
     $           matids(i))
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_set_tag(iter, tagh, size, count, vals)
      implicit none

#include "NEKMOAB"      
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      integer ierr, i, ivals, size, count
      real vals(*), tag_vals
      pointer(tag_ptr, tag_vals(1))

      call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), 
     $     %VAL(iter), tag_ptr, count, ierr)
      IMESH_ASSERT

c set the tag vals
      ivals = size * count
      do i = 1, ivals
         tag_vals(i) = vals(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_set_int_tag(iter, tagh, size, count, vals)
      implicit none

#include "NEKMOAB"      
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      integer ierr, i, ivals, size, count
      integer vals(*), tag_vals
      pointer(tag_ptr, tag_vals(1))

      call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), 
     $     %VAL(iter), tag_ptr, count, ierr)
      IMESH_ASSERT

c set the tag vals
      ivals = size * count
      do i = 1, ivals
         tag_vals(i) = vals(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_set_vertex_tag(iter, tagh, count, vals)
      implicit none

#include "NEKMOAB"      
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      integer ierr, i, ivals, size, count, v_per_e, ntot
      real vals(*), tag_vals(27)
      iBase_EntityHandle connect
      pointer (connect_ptr, connect(1))

      call iMesh_connectIterate(%VAL(imeshh), %VAL(iter), 
     $     connect_ptr, v_per_e, count, ierr)
      IMESH_ASSERT

c only works if nx, ny, nz are equal, and if vpere is 27
      if (nx1 .ne. ny1 .or. nx1 .ne. nz1 .or. v_per_e .ne. 27) then
         ierr = iBase_FAILURE
         IMESH_ASSERT
      endif

      ntot = nx1 * ny1 * nz1

c set the tag vals
      ivals = 1
      do i = 1, count
c        transfer spectral variable to vertex variable
c        corners only
#define INDEX(i,j,k) nx1*ny1*k + nx1*j + i
         tag_vals(1) = vals(1+INDEX(0,0,0))
         tag_vals(2) = vals(1+INDEX(nx1-1,0,0))
         tag_vals(3) = vals(1+INDEX(nx1-1,ny1-1,0))
         tag_vals(4) = vals(1+INDEX(0,ny1-1,0))
         tag_vals(5) = vals(1+INDEX(0,0,nz1-1))
         tag_vals(6) = vals(1+INDEX(nx1-1,0,nz1-1))
         tag_vals(7) = vals(1+INDEX(nx1-1,ny1-1,nz1-1))
         tag_vals(8) = vals(1+INDEX(0,ny1-1,nz1-1))
#undef INDEX
         call iMesh_setDblArrData(%VAL(imeshh), 
     $        connect(ivals), %VAL(8), %VAL(tagh), tag_vals(1), 
     $        %VAL(8), ierr)
         IMESH_ASSERT
         ivals = ivals + v_per_e
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_get_tag(iter, tagh, size, count, vals, ierr)
      implicit none

#include "NEKMOAB"      
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      integer ierr, i, ivals, size, count
      real vals(*), tag_vals
      pointer(tag_ptr, tag_vals(1))

      call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), 
     $     %VAL(iter), tag_ptr, count, ierr)
c don't assert here, just return
      if (iBase_SUCCESS .ne. ierr) return

c set the tag vals
      ivals = size * count
      do i = 1, ivals
         vals(i) = tag_vals(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_get_vertex_tag(iter, tagh, count, vals, ierr)
      implicit none

#include "NEKMOAB"      
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      integer ierr, i, j, ivals, size, count, v_per_e, ntot
      real vals(*), tag_vals(27)
      iBase_EntityHandle connect
      pointer (connect_ptr, connect(1))
      real avg

      call iMesh_connectIterate(%VAL(imeshh), %VAL(iter), 
     $     connect_ptr, v_per_e, count, ierr)
c don't assert here, just return
      if (iBase_SUCCESS .ne. ierr) return

c only works if nx1, ny1, nz1 are equal, and if v_per_e is 27
      if (nx1 .ne. ny1 .or. nx1 .ne. nz1 .or. v_per_e .ne. 27) then
         ierr = iBase_FAILURE
         return
      endif

      ntot = nx1 * ny1 * nz1

c set the tag vals
      ivals = 1
      do i = 1, count
c        transfer spectral variable fromvertex variable
c        corners only
         call iMesh_getDblArrData(%VAL(imeshh), 
     $        connect(ivals), %VAL(8), %VAL(tagh), tag_vals(1), 
     $        %VAL(8), ierr)
c don't assert here, just return
         if (iBase_SUCCESS .ne. ierr) return
         ivals = ivals + v_per_e

#define INDEX(i,j,k) lx1*ly1*k + lx1*j + i
         avg = 0.0d0
         do j = 1, 8
            avg = avg + tag_vals(j)
         enddo
         vals(1+INDEX(0,0,0)) = avg
         vals(1+INDEX(nx1-1,0,0)) = avg
         vals(1+INDEX(nx1-1,ny1-1,0)) = avg
         vals(1+INDEX(0,ny1-1,0)) = avg
         vals(1+INDEX(0,0,nz1-1)) = avg
         vals(1+INDEX(nx1-1,0,nz1-1)) = avg
         vals(1+INDEX(nx1-1,ny1-1,nz1-1)) = avg
         vals(1+INDEX(0,ny1-1,nz1-1)) = avg
#undef INDEX
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_compute_diagnostics()
      implicit none

#include "NEKMOAB"
      include 'GEOM'

      integer i, j, k, l

      integer e, intv
      real avg(3)

      intv = nelt / 10
      e = 1
      do while (e .lt. nelt)
         avg(1) = 0.0
         avg(2) = 0.0
         avg(3) = 0.0
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  avg(1) = avg(1) + xm1(i,j,k,e)
                  avg(2) = avg(2) + ym1(i,j,k,e)
                  if (if3d)
     $                 avg(3) = avg(3) + zm1(i,j,k,e)
               enddo
            enddo
         enddo
         avg(1) = avg(1) / (nx1*ny1*nz1)
         avg(2) = avg(2) / (nx1*ny1*nz1)
         avg(3) = avg(3) / (nx1*ny1*nz1)
         print *, "Average for e is ", e, avg(1), avg(2), avg(3)

         e = e + intv
      enddo

      return
      end

c-----------------------------------------------------------------------
      block data nekMOABdata
#include "NEKMOAB"
      data imeshh/0/, hPartn/0/, fileset/0/, 
     $     rpParts/0/, rpHexes/0/, rpxm1/0/, rpym1/0/, rpzm1/0/, 
     $     rpvx/0/, rpvy/0/, rpvz/0/, rpt/0/, rpp/0/,
     $     globalIdTag/0/, matsetTag/0/, neusetTag/0/,
     $     matsets/numsts*0/, ieiter/numsts*0/, 
     $     xm1Tag/0/, ym1Tag/0/, zm1Tag/0/, vxTag/0/, vyTag/0/, 
     $     vzTag/0/, tTag/0/, 
     $     pTag/0/, dTag/0/, powTag/0/, vtTag/0/, vpTag/0/, vdTag/0/, 
     $     vpowTag/0/, senseTag/0/, 
     $     iCreatedImesh/0/, iCreatedPartn/0/, iCreatedFileset/0/, 
     $     iestart/numsts*0/, iecount/numsts*0/, 
     $     partsSize/0/, hexesSize/0/
      end
c-----------------------------------------------------------------------
      subroutine moab_geometry_tetstris (xmlo,ymlo,zmlo)

      implicit none
#include "NEKMOAB"      
      include 'GEOM'
      include 'PARALLEL'
      real      xmlo(lx1*ly1*lz1,*)
     $        , ymlo(lx1*ly1*lz1,*)
     $        , zmlo(lx1*ly1*lz1,*)
      integer   e,nmoab,i,j,k
      integer   ncrn

      common /ivrtx/ vtx ((ldim+1),lelt)
      integer   vtx
      integer   lc2
      parameter(lc2=(ldim+1))!*(ldim+2)/2) ! lc2=10 (3d), lc2=6 (2d)

      integer indv(4)
      data indv /1,2,3,4/

      common /tcrmg2/ xcm(lc2,lelt), ycm(lc2,lelt), zcm(lc2,lelt)
      common /tcrmg3/ txcm(lc2,lelt), tycm(lc2,lelt), tzcm(lc2,lelt)
      real            xcm, ycm, zcm
      real            txcm,tycm,tzcm

c      open(unit=11,file='mesh1.txt',status='unknown')
      ncrn = 3
      nxyz = nx1*(nx1+1)/2
      if (if3d) nxyz = nx1*(nx1+1)*(nx1+2)/6
      if (if3d) ncrn = 4 

      call nekMOAB_loadCoord(xcm,ycm,zcm,lc2)   !Get coords from moab
      call nekMOAB_loadConn (vtx, nelgt, ncrn)

c-------ARRANGE for new permutation-------
         do e = 1, nelt      
         call isort(vtx(1,e),indv,ncrn)
         do i = 1, lc2
            txcm(i,e) = xcm(indv(i),e)
c            xcm(i,e) = txcm(i,e)
            tycm(i,e) = ycm(indv(i),e)
c            ycm(i,e) = tycm(i,e)
            tzcm(i,e) = zcm(indv(i),e)     
c            zcm(i,e) = tzcm(i,e)
         enddo
         enddo

       if (if3d) then

          do e=1,nelt   

c         call readandmesh3d_fekete(xmlo(1,e),txcm(1,e),nx1-1,nxyz,e)
c         call readandmesh3d_fekete(ymlo(1,e),tycm(1,e),nx1-1,nxyz,e)
c         call readandmesh3d_fekete(zmlo(1,e),tzcm(1,e),nx1-1,nxyz,e)

          enddo

       else

          do e=1,nelt
  
c         call readandmesh2d_fekete(xmlo(1,e),txcm(1,e),nx1-1,nxyz,e)
c         call readandmesh2d_fekete(ymlo(1,e),tycm(1,e),nx1-1,nxyz,e)
     
          enddo

       endif


c       do e = 1,nelt
c       do i = 1,nxyz 
c          write(11,23) e,i,xmlo(i,e),ymlo(i,e),zmlo(i,e)
c23        format(i3,' ',i3,' ',f7.3,f7.3,f7.3)
c       enddo 
c       enddo

      return
      end
c--------------------------------------------------------------------------
      subroutine nekMOAB_loadCoord2(xmc, ymc, zmc, lc)
c     
c     stuff the xyz coords of the 27 verts of each local element -- 
c     shared vertex coords are stored redundantly
c     
      implicit none
#include "NEKMOAB"
      integer lc    
      real    xmc(lc,*), ymc(lc,*), zmc(lc,*)
      IBASE_HANDLE_T connect_ptr
      iBase_EntityHandle connect
      pointer(connect_ptr, connect(0:1))
      integer i, j, k, ierr, e_in_chunk, e_in_set, v_per_e, e_tot

      e_tot = 0
      do i = 1, numflu+numoth
         e_in_set = 0
         do while (e_in_set .lt. iecount(i))
            if (e_in_set .eq. 0) then
               call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $              ierr)
               IMESH_ASSERT
            endif

c     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), 
     $           connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

c     for each element
            do j = 0, e_in_chunk-1
c     get vertex gids for this e
               do k = 1, lc !lc=TWENTYSEVEN,TEN,NINE,SIX
                  call iMesh_getVtxCoord(%VAL(imeshh),
     $                 %VAL(connect(j*v_per_e+k-1)), 
     $                 xmc(k,e_tot+j+1), ymc(k,e_tot+j+1), 
     $                 zmc(k,e_tot+j+1), ierr)
             IMESH_ASSERT
               enddo
            enddo
            e_tot = e_tot + e_in_chunk
            e_in_set = e_in_set + e_in_chunk
         enddo
      enddo

      return
      end
c-------------------------------------------------------------------------
