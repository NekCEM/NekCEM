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
c     Drude model. If
c
c     a = drudea(-)
c     b = drudeb(-)
c
c     then the dielectric function is
c
c     eps(omega) = eps_oo-b/(omega*(omega+i*a))
c
c     where eps_oo is the value stored in `permittivity`.
c
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
c     Lorentz model. This allows for up to three Lorentz poles; the
c     parameters for the i-th pole are stored in
c
c     a = lorena(-,i)
c     b = lorenb(-,i)
c     c = lorenc(-,i).
c
c     This corresponds to a dielectric function of
c
c     eps(omega) = eps_oo+c/(b-i*a*omega-omega**2)
c
c     where eps_oo is the value stored in `permittivity`.
c
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i,j
      real bmv
      real ca1,cb1,cc1
      real ca2,cb2,cc2
      real ca3,cb3,cc3

      do i = 1,ncemlorentz
         j = cemlorentz(i)

         bmv = bmn(j)
         ca1 = lorena(j,1)
         cb1 = lorenb(j,1)
         cc1 = lorenc(j,1)
         ca2 = lorena(j,2)
         cb2 = lorenb(j,2)
         cc2 = lorenc(j,2)
         ca3 = lorena(j,3)
         cb3 = lorenb(j,3)
         cc3 = lorenc(j,3)

         resEN(j,1) = resEN(j,1)-(JN1(j,1)+JN2(j,1)+JN3(j,1))*bmv
         resEN(j,2) = resEN(j,2)-(JN1(j,2)+JN2(j,2)+JN3(j,2))*bmv
         resEN(j,3) = resEN(j,3)-(JN1(j,3)+JN2(j,3)+JN3(j,3))*bmv

         resJN1(j,1) = -ca1*JN1(j,1)-cb1*PN1(j,1)+cc1*EN(j,1)
         resJN1(j,2) = -ca1*JN1(j,2)-cb1*PN1(j,2)+cc1*EN(j,2)
         resJN1(j,3) = -ca1*JN1(j,3)-cb1*PN1(j,3)+cc1*EN(j,3)
         resJN2(j,1) = -ca2*JN2(j,1)-cb2*PN2(j,1)+cc2*EN(j,1)
         resJN2(j,2) = -ca2*JN2(j,2)-cb2*PN2(j,2)+cc2*EN(j,2)
         resJN2(j,3) = -ca2*JN2(j,3)-cb2*PN2(j,3)+cc2*EN(j,3)
         resJN3(j,1) = -ca3*JN3(j,1)-cb3*PN3(j,1)+cc3*EN(j,1)
         resJN3(j,2) = -ca3*JN3(j,2)-cb3*PN3(j,2)+cc3*EN(j,2)
         resJN3(j,3) = -ca3*JN3(j,3)-cb3*PN3(j,3)+cc3*EN(j,3)

         resPN1(j,1) = JN1(j,1)
         resPN1(j,2) = JN1(j,2)
         resPN1(j,3) = JN1(j,3)
         resPN2(j,1) = JN2(j,1)
         resPN2(j,2) = JN2(j,2)
         resPN2(j,3) = JN2(j,3)
         resPN3(j,1) = JN3(j,1)
         resPN3(j,2) = JN3(j,2)
         resPN3(j,3) = JN3(j,3)
      enddo

      return
      end
c-----------------------------------------------------------------------
