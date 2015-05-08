      subroutine GLO_GRAD(D,U,U1,U2,U3)
      implicit none
      include 'SIZE'
      integer i,j,k,l,e
      real D (lx1,lx1)
      real U (lx1,ly1,lz1,lelt)
      real U1(lx1,ly1,lz1,lelt)
      real U2(lx1,ly1,lz1,lelt)
      real U3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (D(NX1,NX1))
!$ACC&     PRESENT (U(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (U1(NX1,NY1,NZ1,NELT),U2(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (U3(NX1,NY1,NZ1,NELT))
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    private(tmpu1,tmpu2,tmpu3)
!dir$ NOBLOCKING
      do e=1,nelt
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  tmpu1 = 0.0
                  tmpu2 = 0.0
                  tmpu3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpu1 = tmpu1 + D(i,l)*U(l,j,k,e)
                     tmpu2 = tmpu2 + D(j,l)*U(i,l,k,e)
                     tmpu3 = tmpu3 + D(k,l)*U(i,j,l,e)
                  enddo
                  U1(i,j,k,e) = tmpu1
                  U2(i,j,k,e) = tmpu2
                  U3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END LOOP
!$ACC END DATA
      return
      end
c---------------------------------------------------------------
      subroutine GLO_2GRAD(D,U1,U2,U3,V1,V2,V3)
      implicit none
      include 'SIZE'
      integer i,j,k,l,e
      real D (lx1,ly1)
      real U1(lx1,ly1,lz1,lelt)
      real U2(lx1,ly1,lz1,lelt)
      real U3(lx1,ly1,lz1,lelt)
      real V1(lx1,ly1,lz1,lelt)
      real V2(lx1,ly1,lz1,lelt)
      real V3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (D (NX1,NX1))
!$ACC&     PRESENT (U1(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (U2(NX1,NY1,NZ1,NELT),U3(NX1,NY1,NZ1,NELT))
!$ACC$     PRESENT (V1(NX1,NY1,NZ1,NELT),V2(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (V3(NX1,NY1,NZ1,NELT))
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    private(tmpu1,tmpu2,tmpu3)
!dir$ NOBLOCKING
      do e=1,nelt
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  tmpu1 = 0.0
                  tmpu2 = 0.0
                  tmpu3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpu1 = tmpu1 + D(i,l)*U1(l,j,k,e)
                     tmpu2 = tmpu2 + D(j,l)*U2(i,l,k,e)
                     tmpu3 = tmpu3 + D(k,l)*U3(i,j,l,e)
                  enddo
                  V1(i,j,k,e) = tmpu1
                  V2(i,j,k,e) = tmpu2
                  V3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END LOOP
!$ACC END DATA
      return
      end
