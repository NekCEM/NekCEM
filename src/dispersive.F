c-----------------------------------------------------------------------
c     Compute currents for dispersive materials. Current only the Drude
c     and Lorentz models are supported.
c-----------------------------------------------------------------------
      subroutine cem_dispersive_setup
c-----------------------------------------------------------------------
c     Count the number of dispersive elements and build an array
c     pointing to their global numbers.
c
c     To determine which elements are dispersive, it uses the
c     user-provided subroutine `uservp`, which should fill the array
c     `vdispmarker` with markers. The key is
c
c     0 -> non-dispersive material
c     1 -> Drude material
c     2 -> Lorentz material
c
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i,j,k

      npts = nx1*ny1*nz1*nelt

      i = 0
      j = 0
      do k = 1,npts
         if (vdispmarker(k).eq.0) then
            cycle
         else if (vdispmarker(k).eq.1) then
            cemdrude(i) = k
            i = i+1
         else if (vdispmarker(k).eq.2) then
            cemlorentz(j) = k
            j = j+1
         else
            write(6,*) "Invalid volume marker in cem_dispersive_setup"
            call exitt
         endif
      enddo
      ncemdrude = i
      ncemlorentz = j

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_drude
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i,j
      real bmv
      real ca1,cb1

      do i = 1,ncemdrude
         j = cemdrude(i)
         bmv = bmn(j)
         ca1 = drudea(j)
         cb1 = drudeb(j)

         resEN(j,1) = resEN(j,1)-JN1(j,1)*bmv
         resEN(j,2) = resEN(j,2)-JN1(j,2)*bmv
         resEN(j,3) = resEN(j,3)-JN1(j,3)*bmv

         resJN1(j,1) = -ca1*JN1(j,1)+cb1*EN(j,1)
         resJN1(j,2) = -ca1*JN1(j,2)+cb1*EN(j,2)
         resJN1(j,3) = -ca1*JN1(j,3)+cb1*EN(j,3)

         resPN1(j,1) = JN1(j,1)
         resPN1(j,2) = JN1(j,2)
         resPN1(j,3) = JN1(j,3)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_lorentz
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i
      real bmv
      real ca1,cb1,cc1
      real ca2,cb2,cc2
      real ca3,cb3,cc3

      do i = 1,npts
         bmv = bmn(i)
         ca1 = lorena(i,1)
         cb1 = lorenb(i,1)
         cc1 = lorenc(i,1)
         ca2 = lorena(i,2)
         cb2 = lorenb(i,2)
         cc2 = lorenc(i,2)
         ca3 = lorena(i,3)
         cb3 = lorenb(i,3)
         cc3 = lorenc(i,3)

         resEN(i,1) = resEN(i,1)-(JN1(i,1)+JN2(i,1)+JN3(i,1))*bmv
         resEN(i,2) = resEN(i,2)-(JN1(i,1)+JN2(i,2)+JN3(i,2))*bmv
         resEN(i,3) = resEN(i,3)-(JN1(i,1)+JN2(i,3)+JN3(i,3))*bmv

         resJN1(i,1) = ca1*JN1(i,1)+cb1*PN1(i,1)+cc1*EN(i,1)
         resJN1(i,2) = ca1*JN1(i,2)+cb1*PN1(i,2)+cc1*EN(i,2)
         resJN1(i,3) = ca1*JN1(i,3)+cb1*PN1(i,3)+cc1*EN(i,3)
         resJN2(i,1) = ca2*JN2(i,1)+cb2*PN2(i,1)+cc2*EN(i,1)
         resJN2(i,2) = ca2*JN2(i,2)+cb2*PN2(i,2)+cc2*EN(i,2)
         resJN2(i,3) = ca2*JN2(i,3)+cb2*PN2(i,3)+cc2*EN(i,3)
         resJN3(i,1) = ca3*JN3(i,1)+cb3*PN3(i,1)+cc3*EN(i,1)
         resJN3(i,2) = ca3*JN3(i,2)+cb3*PN3(i,2)+cc3*EN(i,2)
         resJN3(i,3) = ca3*JN3(i,3)+cb3*PN3(i,3)+cc3*EN(i,3)

         resPN1(i,1) = JN1(i,1)
         resPN1(i,2) = JN1(i,2)
         resPN1(i,3) = JN1(i,3)
         resPN2(i,1) = JN2(i,1)
         resPN2(i,2) = JN2(i,2)
         resPN2(i,3) = JN2(i,3)
         resPN3(i,1) = JN3(i,1)
         resPN3(i,2) = JN3(i,2)
         resPN3(i,3) = JN3(i,3)
      enddo

      return
      end
c-----------------------------------------------------------------------
