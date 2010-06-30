C-------------------------------------------------------------------tttt
      subroutine connect
C
C     Parallel compatible routine to find 
C     connectivity of element structure.
C
C     On Processor 0:
C
C     .Verify right-handedness of elements.
C     .Verify element-to-element reciprocity of BC's
C     .Verify correlation between E-E BC's and physical coincidence
C     .Set rotations
C     .Determine multiplicity
C     .Set up direct stiffness summation arrays.
C
C      All Processors:
C
C     .Disperse/Receive BC and MULT temporary data read from preprocessor.
C
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'NONCON'
      INCLUDE 'SCRCT'
      INCLUDE 'EMWAVE'
      INCLUDE 'MOVWIN'
c
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C
      common /scrmg/ glo_num(4*lx1*ly1*lz1*lelv)
      common /scrns/ vertex (7*lx1*ly1*lz1*lelv)
      integer glo_num,vertex
C
C     Initialize key arrays for Direct Stiffness SUM.
C
      NXL=NX1
      NYL=NY1
      NZL=1+2*(NDIM-2)

      CALL INITDS              ! Ed & PF ordering for edges
      CALL DSSET (NX1,NY1,NZ1) ! SKPDAT corner ordering
      CALL SETEDGE             ! Edge ordering
 
C=================================================
C     Establish (global) domain topology
C=================================================
C
C     .Generate topologically correct mesh data.
C     .Set up element centers, face centers, etc. 
C     .Check  right handedness of elements.
C     .Check  element boundary conditions.
C     .Establish Element-Element rotations
C     .Construct the element to processor map and
C
      CALL GENXYZL    ! Generate physical mesh
      CALL SETSIDE
      !CALL VERIFY

C     After verifying the mesh, it's necessary to map the 
C     elements to processors and set the rotations
C
c     CALL MAPELPR(2)
c
c     IF (NP.GT.1) THEN
c        CALL READAT
c        CALL SETVAR
c        CALL GENXYZL
c        CALL SETSIDE
c     ENDIF
c
c     CALL SETCDOF 
c     IF (IFAXIS            ) CALL SETRZER
c     IF (IFMVBD            ) CALL CBCMESH
c     IF (IFMODEL.AND.IFKEPS) CALL CBCTURB
c     CALL CHKAXCB
c
      CALL GENWZ
      CALL GENCOOR (XM3,YM3,ZM3)
      CALL GEOM1 (XM3,YM3,ZM3) ! Geo-mapping : CXRM1,RX,W,JAC...(coef.f)
      CALL GEOM2

C========================================================================
C     Set up multiplicity and direct stiffness arrays for each IFIELD
C========================================================================
C
      MFIELD=2

      IF (IFFLOW) MFIELD=1
      IF (IFMVBD) MFIELD=0
c
c     must read the map file before calling fxxt_ivertex_map
c
      ncrnr = 2**ndim
c
c     First call, for masks, boundaries, etc.
c
      if (nelgv.eq.nelgt) then
         if (ifgtp) then
           call gen_gtp_vertex    (vertex, ncrnr)
         else
           if (ifmovw) then
              if (np.eq.1) then
              call fxxt_ivertex_mapw1(vertex, nelgt_mw, ncrnr)
              else
              call fxxt_ivertex_mapw (vertex, nelgt_mw, ncrnr)
              !do i = 1,(mw2-mw1+1)*ncrnr 
              !   write(6,*) 'vertex--', ncrnr,i,vertex(i)
              !enddo
              endif
           else
              !vertex is assigned from map file
              call fxxt_ivertex_map  (vertex, nelgt, ncrnr)
           endif
         endif
         ! pass vertex and get glo_num
         call setupds1(gsh_fld(1),nx1,ny1,nz1,nelv,vertex,glo_num)
         call setupds2(gsh_fld(1),nx1,ny1,nz1,nelv,vertex,glo_num)
         gsh_fld(2)=gsh_fld(1)
      else
         write(6,*) 'if',ifflow,ifheat
         call exitt
         call f77_get_vert  (vertex, ncrnr, nelgv, '.map')
         call setupds1(gsh_fld(1),nx1,ny1,nz1,nelv,vertex,glo_num)
         call setupds2(gsh_fld(1),nx1,ny1,nz1,nelv,vertex,glo_num)
c
         call f77_get_vert  (vertex, ncrnr, nelgt, '.mp2')
         call setupds1(gsh_fld(2),nx1,ny1,nz1,nelt,vertex,glo_num)
         call setupds2(gsh_fld(2),nx1,ny1,nz1,nelt,vertex,glo_num)
      endif
c
      call setup_dgds(vertex)
c
      ntotv = nx1*ny1*nz1*nelv
      ntott = nx1*ny1*nz1*nelt
c
      if (ifflow) then
         ifield = 1
         call rone    (vmult,ntotv)
         call dssum   (vmult,nx1,ny1,nz1)
         call invcol1 (vmult,ntotv)
      endif

      
      if (ifheat) then
         ifield = 2
         call rone    (tmult,ntott)
         call dssum   (tmult,nx1,ny1,nz1)
         call invcol1 (tmult,ntott)
      endif
  
      if (ifmvbd)  call copy (wmult,vmult,ntott)

      do ifield=3,nfield                  ! Additional pass. scalrs.
         if (nelg(ifield).eq.nelgv) then
            gsh_fld(ifield) = gsh_fld(1)
            call copy (tmult(1,1,1,1,ifield-1),vmult,ntotv)
         else
            gsh_fld(ifield) = gsh_fld(2)
            call copy (tmult(1,1,1,1,ifield-1),tmult,ntott)
         endif
      enddo
c
      RETURN
      END
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
        !write(6,*) 'DUM ----',DUMMY2,DUMMY     

        IC=IC+1
        IXCN(IC)= 1 + (NX-1)*ICX + NX*(NY-1)*ICY + NX*NY*(NZ-1)*ICZ
        !write(6,*) 'IC',IC,IXCN(IC)
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
       
c     write(6,*) 'SKP 1:',SKPDAT(1,1),SKPDAT(2,1),SKPDAT(3,1)
c    $                   ,SKPDAT(4,1),SKPDAT(5,1),SKPDAT(6,1)
c     write(6,*) 'SKP 2:',SKPDAT(1,2),SKPDAT(2,2),SKPDAT(3,2)
c    $                   ,SKPDAT(4,2),SKPDAT(5,2),SKPDAT(6,2)
c     write(6,*) 'SKP 3:',SKPDAT(1,3),SKPDAT(2,3),SKPDAT(3,3)
c    $                   ,SKPDAT(4,3),SKPDAT(5,3),SKPDAT(6,3)
c     write(6,*) 'SKP 4:',SKPDAT(1,4),SKPDAT(2,4),SKPDAT(3,4)
c    $                   ,SKPDAT(4,4),SKPDAT(5,4),SKPDAT(6,4)
c     write(6,*) 'SKP 5:',SKPDAT(1,5),SKPDAT(2,5),SKPDAT(3,5)
c    $                   ,SKPDAT(4,5),SKPDAT(5,5),SKPDAT(6,5)
c     write(6,*) 'SKP 6:',SKPDAT(1,6),SKPDAT(2,6),SKPDAT(3,6)
c    $                   ,SKPDAT(4,6),SKPDAT(5,6),SKPDAT(6,6)
 
C     Set up skip indices for each of the 12 edges 
C
C         Note that NXY = NX*NY even for 2-D since
C         this branch does not apply to the 2D case anyway.
C
C     ESKIP(*,1) = start location
C     ESKIP(*,2) = end 
C     ESKIP(*,3) = stride
C
      ! ESKIP not used !!!
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
      NXL=NX1
      NYL=NY1
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
      INCLUDE 'SCRCT'
C
      CALL VERRHE
      if (np.eq.1) then
c        Note this routine *can* be called in parallel.  It is simply
c        turned off to reduce start up time.  pff 2/7/98
         CALL VERBCS
      elseif (nid.eq.0) then
         write(6,*) 'NOTE:  No more E-E BC verification in parallel'
      endif
C
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

      IF (IF3D) THEN
         DO 100 J=1,8
         IVTX = INDX(J)
         DO 100 IE=1,NELT
            XYZ(1,IVTX,IE) = XC(J,IE)
            XYZ(2,IVTX,IE) = YC(J,IE)
            XYZ(3,IVTX,IE) = ZC(J,IE)
  100    CONTINUE
      ELSE
         DO 200 J=1,4
         IVTX = INDX(J)
         DO 200 IE=1,NELT
            XYZ(1,IVTX,IE) = XC(J,IE)
            XYZ(2,IVTX,IE) = YC(J,IE)
            XYZ(3,IVTX,IE) = 0.0
  200    CONTINUE
      ENDIF
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
C
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
            IEG=LGLEL(IE,NODE)
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
            IEG=LGLEL(IE,NODE)
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
C
      IF (.NOT.IFCSTT) THEN
         IF (NID.EQ.0) WRITE(6,2003) NELGT
         call exitt
      ELSE
         IF (NID.EQ.0) WRITE(6,2002) NELGT
      ENDIF
C
 2001 FORMAT(//,'  Elemental geometry not right-handed, ABORTING'
     $      ,' in routine VERRHE.')
 2002 FORMAT('  Right-handed check complete for',I4,' elements. OK.')
 2003 FORMAT('  Right-handed check failed for',I4,' elements.'
     $      ,'  Exiting in routine VERRHE.')
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine verbcs
C     8 Mar 1989 21:58:26   PFF
C     Verify reciprocity of boundary conditions.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      INCLUDE 'TOPOL'
      INCLUDE 'TSTEP'
C
      LOGICAL IFCLOS,IFRECP
      CHARACTER*3 CB1,CBJ,CBCHAR(8,2)
      DATA CBCHAR /'E  ','P  ','MSI','MPI','msi','mpi','MLI','MCI' 
     $            ,'E  ','P  ','MPI','MSI','mpi','msi','MCI','MLI'/
      IFCLOS=.TRUE.
      IFRECP=.TRUE.
C
C     In order to check E-E connectivity, we need to know the location
C     of the faces.  Find midpoint of each face or edge by averaging.
C
C     Check Element-Element reciprocity
C
C
C     PARALLEL COMPATIBLE VERSION: 
C       .Set consistency to .FALSE. for each element-face.
C       .Pass each data set around to other processors, checking off
C        each side that is indeed consistent.
C       .Scan list at completion to verify consistency of all faces.
C
      MFIELD=2
c     write(6,*) 'IFFLOW=',IFFLOW
c      IF (IFFLOW) MFIELD=1
      DO 6000 IFIELD=MFIELD,NFIELD
C
         NFACES=2*NDIM
         NEL=NELFLD(IFIELD)
         NELM=NEL
         DO 10 I=1,NP
            NELM = MAX(NELM,NELP(IFIELD,I))
   10    CONTINUE
         DO 100 IFACE=1,NFACES
         DO 100 IE =1,NEL
         IFCNST(IFACE,IE)=.TRUE.
         CB1 = CBC(IFACE,IE,IFIELD)
         IF (CB1.EQ.'P  ' .OR. CB1.EQ.'E  ' .OR. 
     $       CB1.EQ.'MSI' .OR. CB1.EQ.'msi' .OR.
     $       CB1.EQ.'MPI' .OR. CB1.EQ.'mpi' .OR.
     $       CB1.EQ.'MLI' .OR. CB1.EQ.'MCI') 
     $       IFCNST(IFACE,IE)=.FALSE.
         write(6,*) 'CBC==',IFACE,IE,IFIELD,CB1
  100    CONTINUE
         !call exitt
C
         DO 2000 IRG=1,NP
            NXTNID=IPRING(IRG)
            NXNODE=NXTNID+1
            NEL  = NELFLD(IFIELD)
            NELJ = NELP(IFIELD,NXNODE)
C
C           The following data needs to be passed around the ring
C           in order to verify E-E consistency:
C
            CALL CRING(CBCS,CBC(1,1,IFIELD),18*NELM,IRG)
            CALL RRING(BCS,BC(1,1,1,IFIELD),30*NELM,IRG)
            CALL RRING(SIDES,SIDE,          24*NELM,IRG) 
C
C           See how many IFCNST's we can set to .TRUE. on this pass.
C
            DO 1000 ICHR=1,8
            DO 1000 IE=1,NEL
            DO 1000 IFACE=1,NFACES
               IEG = LGLEL(IE,NODE)
               CB1 = CBC(IFACE,IE,IFIELD)
               IF (.NOT.IFCNST(IFACE,IE).AND.
     $            CB1.EQ.CBCHAR(ICHR,1)) THEN
                  JEG    = INT( BC(1,IFACE,IE,IFIELD) )
                  JFACE  = INT( BC(2,IFACE,IE,IFIELD) )
                  JNID   = GLLNID(JEG)
C
C                 Is the opposite face/element a member of the newly arrived
C                 set of data??
C
                  IF (JNID.EQ.NXTNID) THEN
                     JE   = GLLEL(JEG)
                     CBJ  = CBCS(JFACE,JE)
                     JBC1 = INT( BCS(1,JFACE,JE) )
                     JBC2 = INT( BCS(2,JFACE,JE) )
C
C                    Check for reciprocity of E-E bc.s
C
                     IFRECP=.TRUE.
                     IF (CBJ.NE.CBCHAR(ICHR,2) .OR. JBC1.NE.IEG .OR.
     $                  JBC2.NE.IFACE) THEN
                        IFRECP=.FALSE.
                        IF (NP.GT.1) THEN
                           WRITE(6,510) NID,IEG,IFACE,IEG,IFACE
     $                    ,CB1,JEG,JFACE,JEG,JFACE,CBJ,JBC1,JBC2
                        ELSE
                           WRITE(6,511) IEG,IFACE,IEG,IFACE
     $                    ,CB1,JEG,JFACE,JEG,JFACE,CBJ,JBC1,JBC2
                        ENDIF
  510                   FORMAT(1X,I5,' WARNING: Detected inconsistency'
     $                  ,' for element',I4,', side',I2,':'
     $                  ,/,1X,' IE  SIDE',5X,'CBC',5X,'BC1',5X,'BC2'
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3   )
  511                   FORMAT(1X,' WARNING: Detected inconsistency'
     $                  ,' for element',I4,', side',I2,':'
     $                  ,/,1X,' IE  SIDE',5X,'CBC',5X,'BC1',5X,'BC2'
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3   )
                     ENDIF
C
C                    Check for closeness of sides (w/ semi-check for periodic)
C
                     IFC = EFACE1(IFACE)
                     JFC = EFACE1(JFACE)
                     IF (IF3D) THEN
                        DIST = 
     $                     SQRT( (SIDE(1,IFC,IE)-SIDES(1,JFC,JE))**2 +
     $                           (SIDE(2,IFC,IE)-SIDES(2,JFC,JE))**2 +
     $                           (SIDE(3,IFC,IE)-SIDES(3,JFC,JE))**2 )
                     ELSE
                        DIST = 
     $                     SQRT( (SIDE(1,IFC,IE)-SIDES(1,JFC,JE))**2 +
     $                           (SIDE(2,IFC,IE)-SIDES(2,JFC,JE))**2 )
                     ENDIF
                     EPS =.0005*MIN(SIDE(4,IFC,IE),SIDES(4,JFC,JE))
C
c                    WRITE(7,512) IE,JE,IFC,JFC,IFACE,JFACE,DIST,EPS,
c    $               SIDE(1,IFC,IE),SIDE(1,IFC,IE),SIDE(1,IFC,IE),
c    $               SIDES(1,JFC,JE),SIDES(1,JFC,JE),SIDES(1,JFC,JE)
  512                FORMAT(6I4,2E11.4,/,4X,3E14.4,/,4X,3E14.4)
C
                     IFCLOS=.TRUE.
                     IF (ICHR.NE.2 .AND. DIST.GT.EPS) THEN
                        IFCLOS=.FALSE.
                        IF (NP.GT.1) THEN
                           WRITE(6,610) NID,IEG,IFACE,DIST,EPS
     $                                     ,JEG,JFACE
  610                      FORMAT(2X,'WARNING:  Sides too far '
     $                    ,'apart for internal boundary condition.'
     $                    ,/,2X,' NID   EL  FACE  DIST      TOL.'
     $                    ,/,2X,  I4,   I6,  I5,  G10.3,   G10.3
     $                    ,/,2X,  4X,   I6,  I5                 )
                        ELSE
c                          WRITE(7,611)     IEG,IFACE,DIST,EPS
c    $                                     ,JEG,JFACE
                           WRITE(6,611)     IEG,IFACE,DIST,EPS
     $                                     ,JEG,JFACE
  611                      FORMAT(2X,'WARNING:  Sides too far '
     $                    ,'apart for internal boundary condition.'
     $                    ,/,2X,'   EL  FACE  DIST      TOL.'
     $                    ,/,2X,    I6,  I5,  G10.3,   G10.3
     $                    ,/,2X,    I6,  I5                 )
                        ENDIF
                     ENDIF
C                    consistent?
                     IF (IFRECP.AND.IFCLOS) IFCNST(IFACE,IE)=.TRUE.
C                    End of element-element match for this pass of the RING.
                  ENDIF
C                 End of consistency check for this element.
               ENDIF
 1000       CONTINUE
 2000    CONTINUE
C
C        Ring check is completed, now verify that all BC's for this field
C        are consistent:
C
         DO 4000 IFACE=1,NFACES
         DO 4000 IE =1,NEL
         IF (.NOT.IFCNST(IFACE,IE)) THEN
            IEG=LGLEL(IE,NODE)
            WRITE(6,4001) NID,IFIELD,IEG,IFACE
            call exitt
            RETURN
         ENDIF
 4000    CONTINUE
C     End of fields loop
 6000 CONTINUE
C
C     Formats.
C
 4001 FORMAT(1X,I5,1X
     $   ,'Boundary condition consistency check failed for '
     $   ,'field',I2,', element',I5,', face',I2,'.'
     $   ,/,2X,'ABORTING in routine VERBCS.')
 4002 FORMAT('  Boundary condition consistency check, OK.')
C
C     Must be OK:
C
      IF (NID.EQ.0) WRITE(6,4002)
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
      subroutine setupds(gs_handle,nx,ny,nz,nel,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      integer gs_handle,fgs_init
      integer vertex(1),glo_num(1)
c
      call set_vert(glo_num,ngv,nx,nel,vertex,.false.)
c     call outfldi (glo_num,'glo_num   ')

c     Initialize gather-scatter code
 
      nvec      = 3             ! number of vectors in gs for CEM
      if (if3d) nvec = 6
c
      t0 = dclock()
                            
      call fgs_init_vec_sz(nvec)

      ntot      = nx*ny*nz*nel
      gs_handle = fgs_init(glo_num,ntot,np)

      t1 = dclock()
      et = t1-t0
c
      if (nid.eq.0) then
         write(6,1) et,t0,t1,gs_handle
    1    format('fgs_init time ::',1pe11.4,' seconds ',2e11.4,i10)
      endif
c
      return
      end

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
      subroutine setupds1(gs_handle,nx,ny,nz,nel,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      integer gs_handle,fgs_init
      integer vertex(1),glo_num(1)

      ! write(6,*) 'nx, ny, nz,nel,glo_num::',nx,ny,nz,nel              
      call set_vert(glo_num,ngv,nx1,nel,vertex,.false.)  ! 1st call
c
      return
      end
c-----------------------------------------------------------------------
      subroutine setupds2(gs_handle,nx,ny,nz,nel,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      integer gs_handle,fgs_init
      integer vertex(1),glo_num(1)

c     Initialize gather-scatter code
 
      nvec      = 3        ! number of vectors in gs for CEM
      if (if3d) nvec = 6
c
      t0 = dclock()
      call fgs_init_vec_sz(nvec)

      ntot = nx*ny*nz*nel
      gs_handle = fgs_init(glo_num,ntot,np)

      !do i=1,ntot           
      !    write(6,*) 'ggg',i,glo_num(i)
      !enddo
      !call exitt

      t1 = dclock()
      et = t1-t0
c
      if (nid.eq.0) then
         write(6,1) et,t0,t1,gs_handle
    1    format('fgs_init time ::',1pe11.4,' seconds ',2e11.4,i10)
      endif
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
      subroutine outfld (x,txt,nx1,ny1,nz1,nel)

      integer e,nx1,ny1,nz1,nel
      real x(nx1,ny1,nz1,nel)
      character*10 txt

c
      mtot = nx1*ny1*nz1*nel 

      do ie = ne1,1,-ne
         l  = ie-1
         do k = 1,1
            do j=ny1,1,-1

              if (nx1.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.8) write(6,108) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
            enddo
         enddo
      enddo
c
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(4f9.5,2x))
  104 FORMAT(3(5f7.3,2x))
  105 FORMAT(3(5f9.5,10x))
c pff 105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(4(6f9.5,5x))
c 106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
   
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldro (x,txt,ichk)

      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(nx1,ny1,nz1,lelt)
      character*10 txt
c
      integer idum,e
      save idum
      data idum /3/

      if (idum.lt.0) return
c
      mtot = nx1*ny1*nz1*nelv

      if (nx1.gt.7.or.nelv.gt.16) return

      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelv
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelv-ne+1

      do ie = ne1,1,-ne
         l  = ie-1
         do k = 1,1
            if (ie.eq.ne1) write(6,116) txt,k,ie,xmin,xmax,istep,time
            write(6,117)    
            do j=ny1,1,-1
              if (nx1.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
              if (nx1.eq.8) write(6,108) ((x(i,j,k,e+l),i=1,nx1),e=1,ne)
            enddo
         enddo
      enddo
c
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(4f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(3(5f9.5,10x))
c pff 105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(4(6f9.5,5x))
c 106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^  ',/,
     $    5X,'   Y |         ',/,
     $    5X,'     |         ',A10,/,
     $    5X,'     +---->    ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X       ','Step  =',I9,f15.5)
  117 FORMAT(' ')
 
c     if (ichk.eq.1.and.idum.gt.0) call checkit(idum)

      return
      end

c-----------------------------------------------------------------------
      subroutine outfldrp (x,txt,ichk)
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(nx2,ny2,nz2,lelt)
      character*10 txt
c
      integer idum,e
      save idum
      data idum /3/
      if (idum.lt.0) return
c
      mtot = nx2*ny2*nz2*nelv
      if (nx2.gt.7.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelv
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelv-ne+1
      do ie=ne1,1,-ne
         l=ie-1
         do k=1,1
            if (ie.eq.ne1) write(6,116) txt,k,ie,xmin,xmax,istep,time
            write(6,117) 
            do j=ny2,1,-1
              if (nx2.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.8) write(6,108) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
            enddo
         enddo
      enddo
c
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c

c     if (ichk.eq.1.and.idum.gt.0) call checkit(idum)

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_dgds(vertex)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
c
      integer vertex(1)
c
      common /scruz/ glo_num(4*lx1*ly1*lz1*lelt)
      integer glo_num
      integer fgs_init
c
      mx = nx1+2

      call set_vert (glo_num,ngv,mx,nelt,vertex,.false.) ! 2nd call

      ! glo_num has extended index with mx : global info

      ntot = (nx1+2)*(ny1+2)*nelt
      if (if3d) ntot = (nx1+2)*(ny1+2)*(nz1+2)*nelt
c
      mz0 = 1
      mz1 = 1

      if (if3d) mz0 = 0
      if (if3d) mz1 = nz1+1
 
c     Use the below call for development only for x-y 1 elt periodic.
      if (param(26).eq.1) then
          call setvert_periodic(glo_num,ngv,mx,nelt,vertex,.false.) 
          if (np.gt.1) then
           if (nid.eq.0) write(6,*) 'setvert_p:np should be 1, np=',np
           call exitt
          endif
c     !!! ATTENTION: 
c     !!! don't run it in parallel <-- this is only serial runs
c     !!! since the range of the glo_num are conflicting. 11/03/08 misun
      endif

      call icopy       (vertex,glo_num,ntot)
      call iface_vert2 (glo_num,vertex,mz0,mz1,nelt) 
 
c     Initialize gather-scatter code
 
      t0  = dclock()
      ntot= nx1*nz1*2*ndim*nelt

      if (nid.eq.0)  write(6,*) 'fgs_ini:',ntot,np
c
      nvec = 3
      if (if3d) nvec = 6

      call fgs_init_vec_sz(nvec)
      gsh_face  = fgs_init(glo_num,ntot,np)  

      t1 = dclock()
      et = t1-t0
c
      if (nid.eq.0) then
         write(6,1) et,t0,t1,gsh_face
    1    format('fgs_init face ::',1pe11.4,' seconds ',2e11.4,i10)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine iface_vert2d(fa,va,jz0,jz1,nel)
C
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
      subroutine iface_vert2(fa,va,jz0,jz1,nel)
C
      include 'SIZE'
      integer fa(nx1*nz1,2*ndim,nel),va(0:nx1+1,0:ny1+1,jz0:jz1,nel)
      integer e,f
c
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
            !write(6,*) 'va',ix,iy,iz,fa(i,f,e)
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
c
      subroutine set_vert(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c
c     Given global array, vertex, pointing to hex vertices, set up
c     a new array of global pointers for an nx^ndim set of elements.
c
c     No communication required, because vertex has all the data.
c
c     The output can go straight into gs_init:
c
c         call        gs_init_vec_sz(ndim)
c         gs_handle = gs_init(glo_num,n,NP)
c
c     where n := nx^ndim * nel.
c
      include 'SIZE'
      include 'INPUT'
c
      integer glo_num(1),vertex(1),ngv,nx
      logical ifcenter
c
      if (if3d) then
         call setvert3d (glo_num,ngv,nx,nel,vertex,ifcenter)
      else
         if (nid.eq.0) write(6,*) 'SETVERT2D A:',nel
         call setvert2d (glo_num,ngv,nx,nel,vertex,ifcenter)
         if (nid.eq.0) write(6,*) 'SETVERT2D B:',nel
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setvert2d(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     set up gsexch interfaces for direct stiffness summation.  pff 2/3/98
c     hmt revisited 12/10/01
c
c
      include 'SIZE'
      include 'CTIMER'
      include 'PARALLEL'
      include 'TOPOL'
      include 'GEOM'
c
      integer glo_num(1),vertex(0:1,0:1,nelgt),ngv,nx,nel
      logical ifcenter
c
c     4 vertices, 4 edges, 1 center
c
      integer  s(3,9,lelg),snum(9*lelg)
      common  /vptsol/ s,snum
c
      integer gvf(6),facet(6),aa(6),key(3),key2(0:3)
c
c     do i=1,nel
c        write(6,1) (vertex(k,0,i),k=0,3),i,nel
c  1     format(4i9,2i5,' VERTEX')
c     enddo
c
c     memory check...
c
      ny=nx
      if (nid.eq.0) write(6,*) 'in setvert2dn',nx,ny
c
      isz_have= 27*lx1*ly1*lz1*lelv
      isz_need= 16*nelgt
      if (isz_need.gt.isz_have) then
         if (nid.eq.0) write(6,*)
c    $      'Abort in setvert2dn.  Increase storage,',isz_need,isz_have
c        call exitt
      endif
c
      key(1)=1
      key(2)=2
      key(3)=3
c
c     build big structure list, vertices, edges, faces
c
      call izero(s,27*nelgv)
      do ie=1,nelgv
         l=0
         do j=0,1           !VERTICES
         do i=0,1
            l=l+1
            s(1,l,ie) = vertex(i,j,ie) 
            s(2,l,ie) = vertex(i,j,ie)
            s(3,l,ie) = 1   ! vertex
         enddo
         enddo
c
         do j=0,1           !EDGES
            l=l+1
            do i=0,1                                    ! x-edge
               i0 = min(vertex(0,j,ie),vertex(1,j,ie)) 
               i1 = max(vertex(0,j,ie),vertex(1,j,ie))
               s(i+1,l,ie) = i*i1 + (1-i)*i0            ! store w/ min leading
            enddo
            s(3,l,ie) = 2   ! vertex
         enddo
         do j=0,1
            l=l+1
            do i=0,1                                    ! y-edge
               i0 = min(vertex(j,0,ie),vertex(j,1,ie)) 
               i1 = max(vertex(j,0,ie),vertex(j,1,ie))
               s(i+1,l,ie) = i*i1 + (1-i)*i0            ! store w/ min leading
            enddo
            s(3,l,ie) = 2   ! vertex
         enddo
c
         imin = vertex(0,0,ie)
         do j=0,1           !CENTER
         do i=0,1
            if (vertex(i,j,ie).le.imin) then
               i0 = 1-i
               j0 = 1-j
               iopp = vertex(i0,j0,ie)
               imin = vertex(i ,j ,ie)
            endif
         enddo
         enddo
         l=l+1
         s(1,l,ie) = imin     ! store w/ min leading
c        s(2,l,ie) = iopp     ! store w/ min leading
c        s(3,l,ie) = 0        ! center
         s(2,l,ie) = 0        ! store w/ min leading
         s(3,l,ie) = ie       ! center
c
      enddo
c
c     compute running tally of number of dof's 
c
      key2(0) = (nx-2)**2
      key2(1) = 1
      key2(2) = nx-2
c
      nstruct = 9*nelgt
      call irank_vec_tally(snum,n_tally,s,3,nstruct,key,3,key2,aa)
c
c     Assign vertices.  Assume hypercube ordering.
c
      do ie =1,nel
         ieg=lglel(ie,node)
         l =0
         do j=0,1
         do i=0,1
            l=l+1
            il  = 1 + (nx-1)*i + nx*(nx-1)*j
            ile = il + nx*nx*(ie-1)
            glo_num(ile)   = s(1,l,ieg)+1
            !write(6,*) 'setvert2d: glo-',il,ile,glo_num(ile)               
            !set global numbring for vertices
         enddo
         enddo
      enddo
c
      ! FIXME
      ! iedg_loc never set, set to zero here. ak, 2006-07-31
      ! set global numbring for edges         

      iedg_loc = 0

      n_on_edge = nx-2
      do ie=1,nel
         ieg = lglel(ie,node)
         l=2**ndim
c
c        Edges 1-2
         do j=0,1                                            ! x-edge
            l   = l+1
            igv = s(1,l,ieg)
            i0  = nx*(nx-1)*j
            i0e = i0 + nx*nx*(ie-1)
            if (glo_num(i0e+1).gt.glo_num(i0e+nx)) then
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
            l   = l+1
            igv = s(1,l,ieg)
            i0  = 1+(nx-1)*i
            i0e = i0 + nx*nx*(ie-1)
            if (glo_num(i0e).gt.glo_num(i0e+nx*(nx-1))) then
               do j=2,nx-1                                   ! std forward case
                  glo_num(i0e+(j-1)*nx) = igv + j-1
               enddo
            else
               do j=2,nx-1                                   ! backward case
                  glo_num(i0e+(j-1)*nx) = igv + 1 + n_on_edge-(j-1)
               enddo
            endif
         enddo
      enddo
c
      if(ifcenter) then 
         do ie=1,nel           ! Centers
            ieg = lglel(ie,node)
            igv = s(1,9,ieg)
            k   = 0
            do j=2,nx-1
               do i=2,nx-1
                  ii = i+(j-1)*nx + nx*nx*(ie-1)
                  k  = k+1
                  glo_num(ii) = igv + k
               enddo
            enddo
         enddo
      else
         do ie=1,nel           ! Centers
            ieg = lglel(ie,node)
            igv = s(1,9,ieg)
            k   = 0
            do j=2,nx-1
               do i=2,nx-1
                  ii = i+(j-1)*nx + nx*nx*(ie-1)
                  k  = k+1
                  glo_num(ii) = 0
               enddo
            enddo
         enddo
      endif
c
      ngv = n_tally
c
c     Quick check on maximum #dofs:
      m   = nx*nx*nelt
      ngvm= iglmax(glo_num,m)

c.... glo_num defined at vertices and edges

      if (nid.eq.0) write(6,*) 'NUM Coarse verts:',ngvm,ngv,n_tally
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setvert3d(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     set up gsexch interfaces for direct stiffness summation.  pff 2/3/98
c     hmt revisited 12/10/01
c
c
      include 'SIZE'
      include 'CTIMER'
      include 'PARALLEL'
      include 'TOPOL'
      include 'GEOM'
c
      integer  glo_num(1),vertex(0:1,0:1,0:1,nelgt),ngv,nx,nel
      logical  ifcenter
c
      integer  edge(0:1,0:1,0:1,3,lelg),enum(12,lelg)
      integer  face(3,6,lelg),fnum(6,lelg)
      common  /vptsol/ edge,enum
      equivalence  (edge,face),(enum,fnum)
c
      integer gvf(4),ind(4),facet(4),aa(3),key(3)
      logical ifij
 
c     memory check...
c
c     do i=1,nel     
c        write(6,11) (vertex(k,0,0,i),k=0,7),i
c 11     format(9i8,' vertex')
c     enddo
c     call exitt
 
      ny   = nx
      nz   = nx
      nxyz = nx*ny*nz

      if (nid.eq.0) write(6,*) 'in setvert3d',nx,ny,nz

      isz_have= 27*lx1*ly1*lz1*lelv
      isz_need= 36*nelgt

      if (isz_need.gt.isz_have) then
c     if (nid.eq.0) write(6,*)
c     $'Abort in setvert3d.  Increase storage,',isz_need,isz_have
c     call exitt
      endif
 
      key(1)=1
      key(2)=2
      key(3)=3
c
c     Count number of unique vertices
      nlv = 2**ndim
      ngv = ivlmax(vertex,nlv*nelgt)
c
c     Assign hypercube ordering of vertices.
      do ie=1,nel   
         ieg = lglel(ie,node)
         do k=0,1
         do j=0,1
         do i=0,1
c           Local to global node number (vertex)
            il  = 1 + (nx-1)*i + nx*(nx-1)*j + nx*nx*(nx-1)*k
            ile = il + nx*ny*nz*(ie-1)
            glo_num(ile) = vertex(i,j,k,ieg)
           !write(6,*) 'glo_num:: ',nelgt,nel,ieg,ile,vertex(i,j,k,ieg)
         enddo
         enddo
         enddo
      enddo

      if (nx.eq.2) return  ! nx1=2
c
c     Assign edge labels by bounding vertices.  
      do ie=1,nelgt
         do k=0,1
         do j=0,1
         do i=0,1
            edge(i,j,k,1,ie) = vertex(i,j,k,ie)  ! r-edge
            edge(j,i,k,2,ie) = vertex(i,j,k,ie)  ! s-edge
            edge(k,i,j,3,ie) = vertex(i,j,k,ie)  ! t-edge
         enddo
         enddo
         enddo
      enddo
c
c     Sort edges by bounding vertices.
      do i=0,12*nelgt-1
         if (edge(0,i,0,1,1).gt.edge(1,i,0,1,1)) then
            kswap = edge(0,i,0,1,1)
            edge(0,i,0,1,1) = edge(1,i,0,1,1)
            edge(1,i,0,1,1) = kswap
         endif
      enddo

c     Assign a number (rank) to each unique edge

      call irank_vec(enum,n_unique_edges,edge,2,12*nelgt,key,2,aa)
c
c= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c     Assign global vertex numbers to SEM nodes on each edge
      n_on_edge = nx-2
      do ie=1,nel
c
         ieg = lglel(ie,node)
         iedg_loc = 0
c
c        Edges 1-4
         do k=0,1
         do j=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,ieg)-1)
            i0  = nx*(nx-1)*j + nx*nx*(nx-1)*k
            i0e = i0 + nxyz*(ie-1)
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
            igv = ngv + n_on_edge*(enum(iedg_loc+1,ieg)-1)
            i0  = 1+(nx-1)*i + nx*nx*(nx-1)*k
            i0e = i0 + nxyz*(ie-1)
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
            igv = ngv + n_on_edge*(enum(iedg_loc+1,ieg)-1)
            i0  = 1 + (nx-1)*i + nx*(nx-1)*j
            i0e = i0 + nxyz*(ie-1)
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
c
c     Currently assigned number of vertices
      ngvv  = ngv
      ngv   = ngv + n_unique_edges*n_on_edge
c
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
c
      nfaces=ndim*2
      ncrnr =2**(ndim-1)
      do ie=1,nelgt
         do ifac=1,nfaces
            do icrn=1,ncrnr
               i                  = icface(icrn,ifac)-1
               facet(icrn)        = vertex(i,0,0,ie)
            enddo
            call isort(facet,gvf,ncrnr)
            call icopy(face(1,ifac,ie),facet,ncrnr-1)
         enddo
      enddo
c
c     Assign a number (rank) to each unique face
      call irank_vec(fnum,n_unique_faces,face,3,6*nelgt,key,3,aa)
c
c     Now assign global node numbers on the interior of each face
c
      call dsset (nx,ny,nz)
      do ie=1,nel
       ieg= lglel(ie,node)
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
         gvf(1) = glo_num(i0+nx*(j0-1)+nxyz*(ie-1))
         gvf(2) = glo_num(i1+nx*(j0-1)+nxyz*(ie-1))
         gvf(3) = glo_num(i0+nx*(j1-1)+nxyz*(ie-1))
         gvf(4) = glo_num(i1+nx*(j1-1)+nxyz*(ie-1))
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
         ig0 = ngv + n_on_face*(fnum(iface,ieg)-1)
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
                  glo_num(i+nx*(j-1)+nxyz*(ie-1)) = l + ig0
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
                  glo_num(i+nx*(j-1)+nxyz*(ie-1)) = l + ig0
               endif
            enddo
            enddo
         endif
       enddo
      enddo
c
c     Finally,  number interiors  
c     ngvs := number of global vertices on surface of subdomains
c
      ngve  = ngv
      ngv   = ngv + n_unique_faces*n_on_face
      ngvs  = ngv
c
      n_in_interior = (nx-2)*(ny-2)*(nz-2)
      if (ifcenter) then
         do ie=1,nel
            ig0 = ngv + n_in_interior*(lglel(ie,node)-1)
            l = 0
            do k=2,nz-1
               do j=2,ny-1
                  do i=2,nx-1
                     l = l+1
                     glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(ie-1)) = ig0+l
                  enddo
               enddo
            enddo
         enddo
      else
         do ie=1,nel
            ig0 = ngv + n_in_interior*(lglel(ie,node)-1)
            l = 0
            do k=2,nz-1
               do j=2,ny-1
                  do i=2,nx-1
                     l = l+1
                     glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(ie-1)) = 0
                  enddo
               enddo
            enddo
         enddo
      endif
c
      ngv = ngv + n_in_interior*nelgt
c
c     Quick check on maximum #dofs:

      m    = nxyz*nelt
      ngvm = iglmax(glo_num,m)

c     do i=1,m
c        write(6,*) 'glo_nnn',i,glo_num(i)
c     enddo

      if (nid.eq.0) write(6,1) nx,ngvv,ngve,ngvs,ngv,ngvm
    1 format('setupds3d:',9i9)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine setvert3dw(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     set up gsexch interfaces for direct stiffness summation.  pff 2/3/98
c     hmt revisited 12/10/01
c
c
      include 'SIZE'
      include 'CTIMER'
      include 'PARALLEL'
      include 'TOPOL'
      include 'GEOM'
c
      integer  glo_num(1),vertex(0:1,0:1,0:1,nelgt),ngv,nx,nel
      logical  ifcenter
c
      integer  edge(0:1,0:1,0:1,3,lelg),enum(12,lelg)
      integer  face(3,6,lelg),fnum(6,lelg)
      common  /vptsol/ edge,enum
      equivalence  (edge,face),(enum,fnum)
c
      integer gvf(4),ind(4),facet(4),aa(3),key(3)
      logical ifij
 
c     memory check...
c
      do i=1,nel   
         write(6,11) (vertex(k,0,0,i),k=0,7),i
  11     format(9i8,' vertex')
      enddo
 
      ny   = nx
      nz   = nx
      nxyz = nx*ny*nz

      if (nid.eq.0) write(6,*) 'in setvert3d',nx,ny,nz

      isz_have= 27*lx1*ly1*lz1*lelv
      isz_need= 36*nelgt

      if (isz_need.gt.isz_have) then
c     if (nid.eq.0) write(6,*)
c     $'Abort in setvert3d.  Increase storage,',isz_need,isz_have
c     call exitt
      endif
 
      key(1)=1
      key(2)=2
      key(3)=3
c
c     Count number of unique vertices
      nlv = 2**ndim
      ngv = ivlmax(vertex,nlv*nelgt)
c
c     Assign hypercube ordering of vertices.
      do ie=1,nelgt 
         ieg = lglel(ie,node)
         do k=0,1
         do j=0,1
         do i=0,1
c           Local to global node number (vertex)
            il  = 1 + (nx-1)*i + nx*(nx-1)*j + nx*nx*(nx-1)*k
            ile = il + nx*ny*nz*(ie-1)
            glo_num(ile) = vertex(i,j,k,ieg)
            write(6,*) 'glo_num:: ',ieg,ile,vertex(i,j,k,ieg)
         enddo
         enddo
         enddo
      enddo


      if (nx.eq.2) return  ! nx1=2
c
c     Assign edge labels by bounding vertices.  
      do ie=1,nelgt
         do k=0,1
         do j=0,1
         do i=0,1
            edge(i,j,k,1,ie) = vertex(i,j,k,ie)  ! r-edge
            edge(j,i,k,2,ie) = vertex(i,j,k,ie)  ! s-edge
            edge(k,i,j,3,ie) = vertex(i,j,k,ie)  ! t-edge
         enddo
         enddo
         enddo
      enddo
c
c     Sort edges by bounding vertices.
      do i=0,12*nelgt-1
         if (edge(0,i,0,1,1).gt.edge(1,i,0,1,1)) then
            kswap = edge(0,i,0,1,1)
            edge(0,i,0,1,1) = edge(1,i,0,1,1)
            edge(1,i,0,1,1) = kswap
         endif
      enddo

c     Assign a number (rank) to each unique edge

      call irank_vec(enum,n_unique_edges,edge,2,12*nelgt,key,2,aa)
c
c= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c     Assign global vertex numbers to SEM nodes on each edge
      n_on_edge = nx-2
      do ie=1,nel
c
         ieg = lglel(ie,node)
         iedg_loc = 0
c
c        Edges 1-4
         do k=0,1
         do j=0,1
            igv = ngv + n_on_edge*(enum(iedg_loc+1,ieg)-1)
            i0  = nx*(nx-1)*j + nx*nx*(nx-1)*k
            i0e = i0 + nxyz*(ie-1)
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
            igv = ngv + n_on_edge*(enum(iedg_loc+1,ieg)-1)
            i0  = 1+(nx-1)*i + nx*nx*(nx-1)*k
            i0e = i0 + nxyz*(ie-1)
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
            igv = ngv + n_on_edge*(enum(iedg_loc+1,ieg)-1)
            i0  = 1 + (nx-1)*i + nx*(nx-1)*j
            i0e = i0 + nxyz*(ie-1)
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
c
c     Currently assigned number of vertices
      ngvv  = ngv
      ngv   = ngv + n_unique_edges*n_on_edge
c
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
c
      nfaces=ndim*2
      ncrnr =2**(ndim-1)
      do ie=1,nelgt
         do ifac=1,nfaces
            do icrn=1,ncrnr
               i                  = icface(icrn,ifac)-1
               facet(icrn)        = vertex(i,0,0,ie)
            enddo
            call isort(facet,gvf,ncrnr)
            call icopy(face(1,ifac,ie),facet,ncrnr-1)
         enddo
      enddo
c
c     Assign a number (rank) to each unique face
      call irank_vec(fnum,n_unique_faces,face,3,6*nelgt,key,3,aa)
c
c     Now assign global node numbers on the interior of each face
c
      call dsset (nx,ny,nz)
      do ie=1,nel
       ieg= lglel(ie,node)
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
         gvf(1) = glo_num(i0+nx*(j0-1)+nxyz*(ie-1))
         gvf(2) = glo_num(i1+nx*(j0-1)+nxyz*(ie-1))
         gvf(3) = glo_num(i0+nx*(j1-1)+nxyz*(ie-1))
         gvf(4) = glo_num(i1+nx*(j1-1)+nxyz*(ie-1))
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
         ig0 = ngv + n_on_face*(fnum(iface,ieg)-1)
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
                  glo_num(i+nx*(j-1)+nxyz*(ie-1)) = l + ig0
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
                  glo_num(i+nx*(j-1)+nxyz*(ie-1)) = l + ig0
               endif
            enddo
            enddo
         endif
       enddo
      enddo
c
c     Finally,  number interiors  
c     ngvs := number of global vertices on surface of subdomains
c
      ngve  = ngv
      ngv   = ngv + n_unique_faces*n_on_face
      ngvs  = ngv
c
      n_in_interior = (nx-2)*(ny-2)*(nz-2)
      if (ifcenter) then
         do ie=1,nel
            ig0 = ngv + n_in_interior*(lglel(ie,node)-1)
            l = 0
            do k=2,nz-1
               do j=2,ny-1
                  do i=2,nx-1
                     l = l+1
                     glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(ie-1)) = ig0+l
                  enddo
               enddo
            enddo
         enddo
      else
         do ie=1,nel
            ig0 = ngv + n_in_interior*(lglel(ie,node)-1)
            l = 0
            do k=2,nz-1
               do j=2,ny-1
                  do i=2,nx-1
                     l = l+1
                     glo_num(i+nx*(j-1)+nx*ny*(k-1)+nxyz*(ie-1)) = 0
                  enddo
               enddo
            enddo
         enddo
      endif
c
      ngv = ngv + n_in_interior*nelgt
c
c     Quick check on maximum #dofs:

      m    = nxyz*nelt
      ngvm = iglmax(glo_num,m)
      if (nid.eq.0) write(6,1) nx,ngvv,ngve,ngvs,ngv,ngvm
    1 format('setupds3d:',9i9)
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine setvert_periodic(glo_num,ngv,nx,nel,vertex,ifcenter)
c
c     set up gsexch interfaces for direct stiffness summation.  pff 2/3/98
c     hmt revisited 12/10/01
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
c     include 'CTIMER'
c     include 'PARALLEL'
c     include 'TOPOL'
c     include 'GEOM'
c
      integer glo_num(1),vertex(0:1,0:1,nelgt),ngv,nx,nel
      logical ifcenter
c
c     4 vertices, 4 edges, 1 center
c
      integer  s(3,9,lelg),snum(9*lelg)
      common  /vptsol/ s,snum
c
      ny=nx

      if (nid.eq.0) write(6,*) 'in setvert2dn',nx,ny

      if (nel.eq.1) then

      if (if3d)   then
 
       do k=1,nx     
       do j=1,nx
       do i=1,nx

          ii = i+(j-1)*nx+(k-1)*nx*nx
          i0 = i+1+(k-1)*nx*nx                      
          j0 = i+(k-1)*nx*nx                      
          k0 = i+(j-1)*nx                         

          if (k.ge.1) then
           if (glo_num(ii).eq.0) then
               glo_num(ii) = glo_num(ii)                           
           else
               glo_num(ii) = glo_num(ii)+ (j-1)*nx*nx
     $                                  + (k-1)*nx*nx*nx    
               if (j.eq.nx) then
                   glo_num(ii)=glo_num(j0)
               endif
           endif
          endif

          if (k.eq.nx) then
              glo_num(ii)=glo_num(k0)
          endif
          !write(6,*)  'x-y-z elt=1 periodic 3d',ii,i,j,k,glo_num(ii)                   

       enddo
       enddo
       enddo
       !call exitt

      else

       do j=1,nx
       do i=1,nx

          i0 = 1+(j-1)*nx
          j0 = i
          ii = i+(j-1)*nx

          glo_num(ii) = glo_num(ii)

          if (i.eq.nx) glo_num(ii) = glo_num(i0)

          if  (j.ge.2) then

          if  (j.eq.nx) then
               glo_num(ii) = glo_num(j0)
          else

            if ((i.eq.1 ).or.(i.eq.nx))  then
               glo_num(ii) = nx*nx*nx+j
               glo_num(ii) = glo_num(ii)
            endif

          endif
          endif
          !write(6,*)  ' x-y elt=1 periodic 2d',ii,i,j,k,glo_num(ii)                   


       enddo
       enddo

      endif

      else
 
       IFPP=.true.
       if (IFPP) then

       do ie=1,nelt   
       do k=1,nx     
       do j=1,nx
       do i=1,nx

          ii = i+(j-1)*nx+(k-1)*nx*nx+(ie-1)*nx*nx*nx
          i0 = i+1+(k-1)*nx*nx                      
          j0 = i+(k-1)*nx*nx                      
          k0 = i+(j-1)*nx+(nx-1)*nx*nx+(ie-2)*nx*nx*nx             

          if (glo_num(ii).eq.0) then
              glo_num(ii) = glo_num(ii)                           
          else
              glo_num(ii) = glo_num(ii)+ (j-1)*nx*nx
     $               +(k-1)*nx*nx*nx +(ie-1)*nx*nx*nx*nx
              if ((k.eq.1).and.(ie.ge.2)) then
                   glo_num(ii)=glo_num(k0)
              endif
          endif
          !write(6,*) 'x-y elt=1 periodic: 3d ',ii,ie,i,j,k,glo_num(ii)            

       enddo
       enddo
       enddo
       enddo

       endif

       call setvert_periodic2(glo_num,nx)

       !do ie=1,nelt   
       !do k=1,nx     
       !do j=1,nx
       !do i=1,nx
       !   ii = i+(j-1)*nx+(k-1)*nx*nx+(ie-1)*nx*nx*nx
       !   write(6,*) 'x-y',ii,ie,i,j,k,glo_num(ii)            
       !enddo
       !enddo
       !enddo
       !enddo
       !call exitt

      endif

      return
      end 

c-----------------------------------------------------------------------
      subroutine setvert_periodic2(glo_num,nx)
c
c     set up gsexch interfaces for:                      (10/18/08)
c                       . DG + periodic in x-y 
c     
c
      include 'SIZE'
      include 'TOTAL'

      integer glo_num(nx,nx,nx,nelv)
      integer e

      if (if3d) then
         do e=1,nelt
            do k=1,nx
            do i=1,nx
               glo_num(i,nx,k,e) = glo_num(i,1,k,e)
               glo_num(nx,i,k,e) = glo_num(1,i,k,e)
            enddo
            enddo
         enddo
      else
         do e=1,nelt
            do i=1,nx
               glo_num(i,nx,e,1) = glo_num(i,1,e,1)
               glo_num(nx,i,e,1) = glo_num(1,i,e,1)
            enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
