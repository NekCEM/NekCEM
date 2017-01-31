c-----------------------------------------------------------------------
c     Compute surface currents for 2D materials.
c-----------------------------------------------------------------------
      subroutine cem_current_bc_ptr
c-----------------------------------------------------------------------
c     Count the number of conducting faces and build an array pointing
c     to their global face numbers.
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'

      integer i,j,k

      k = 0
      do j = 1,ncemface
         if (facef(j).eq.40.0) then
            k = k+1
            cemcur(k) = j
         endif
      enddo
      ncemcur = k

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_volume_current
c-----------------------------------------------------------------------
c     Set volume currents with the current from the faces. Useful for
c     visualizing the currents in Paraview.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer i,j,k

      do i = 1,ncemcur
c     j is the global face number
         j = cemcur(i)
c     k is the volume global number associated with face j.
         k = cemface(j)

         JN(k,1) = fJN(j,1)
         JN(k,2) = fJN(j,2)
         JN(k,3) = fJN(j,3)
         JN1(k,1) = fJN1(j,1)
         JN1(k,2) = fJN1(j,2)
         JN1(k,3) = fJN1(j,3)
         PN1(k,1) = fPN1(j,1)
         PN1(k,2) = fPN1(j,2)
         PN1(k,3) = fPN1(j,3)
         PN2(k,1) = fPN2(j,1)
         PN2(k,2) = fPN2(j,2)
         PN2(k,3) = fPN2(j,3)
         PN3(k,1) = fPN3(j,1)
         PN3(k,2) = fPN3(j,2)
         PN3(k,3) = fPN3(j,3)
         QN1(k,1) = fQN1(j,1)
         QN1(k,2) = fQN1(j,2)
         QN1(k,3) = fQN1(j,3)
         QN2(k,1) = fQN2(j,1)
         QN2(k,2) = fQN2(j,2)
         QN2(k,3) = fQN2(j,3)
         QN3(k,1) = fQN3(k,1)
         QN3(k,2) = fQN3(k,2)
         QN3(k,3) = fQN3(k,3)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_graphene_current
c-----------------------------------------------------------------------
      implicit none

      call cem_te_graphene_current

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_te_graphene_current
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i,j
      real f,tmp,Yfac
      real Y0,Y1
      real a_d,b_d
      real b_cp1,b_cp2
      real a_211,a_221,b_11,b_21
      real a_212,a_222,b_12,b_22

c     FIXME: figure out Drude coefficients
      a_d = 1.0
      b_d = 1.0

      a_211 = 0.463
      a_221 = 1.091
      a_211 = 1.0
      a_221 = 1.0
      b_11 = -0.0012
      b_11 = 0.0
      b_21 = -0.00026
      b_21 = 1.0

      a_212 = 0.585
      a_222 = 1.238
      b_12 = -0.0029
      b_22 = -0.00049

      b_cp1 = -0.0022
      b_cp2 = -0.0051

      do i = 1,ncemcur
c     Global face number for the i-th graphene face.
         j = cemcur(i)
c     Build up the total current
         Yfac = 0.5/Y_0(j)
         tmp = Yfac*(fHN(j,3)+zimpede(j)*fEN(j,1))
         fJN(j,1) = fJN1(j,1)+fPN1(j,1)+fPN2(j,1)-(b_cp1+b_cp2)*tmp
         fJN(j,1) = fJN(j,1)/(1.0 + (b_cp1 + b_cp2)*Yfac)
c     Forcing term
         f = tmp-Yfac*fJN(j,1)
c     Drude term
         resfJN1(j,1) = -a_d*fJN1(j,1)+b_d*f
c     First critical point term
         resfPN2(j,1) = fQN2(j,1)+b_11*f
	 resfQN2(j,1) = -a_211*fPN2(j,1)-a_221*fQN2(j,1)+b_21*f
c     Second critical point term
         resfPN3(j,1) = fQN3(j,1)+b_12*f
	 resfQN3(j,1) = -a_212*fPN3(j,1)-a_222*fQN3(j,1)+b_22*f
      enddo

      i = 1
      j = cemcur(i)
      write(*,*) 'JJJ',fJN1(j,1),fPN2(j,1),fPN3(j,1)
      
      return
      end
c-----------------------------------------------------------------------
