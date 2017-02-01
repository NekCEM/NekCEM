c-----------------------------------------------------------------------
c     Compute surface currents for 2D materials.
c-----------------------------------------------------------------------
      subroutine cem_surface_current_setup
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
      include '2DMAT'

      integer i,j
      real nH_x,nH_y,nEn_x,nEn_y
      real f1,f2,tmp1,tmp2,Yfac,cpfac,jnfac

      do i = 1,ncemcur
c     Global face number for the i-th graphene face.
         j = cemcur(i)
c     (n x H)_x
         nH_x = unym(j)*fHN(j,3)
c     (n x H)_y
         nH_y = -unxm(j)*fHN(j,3)
c     (n x (E x n))_x
         nEn_x = unym(j)**2*fEN(j,1)-unxm(j)*unym(j)*fEN(j,2)
c     (n x (E x n))_y
         nEn_y = unxm(j)**2*fEN(j,2)-unxm(j)*unym(j)*fEN(j,1)
         Yfac = 0.5/Y_0(j)
         tmp1 = Yfac*(nH_x + yconduc(j)*nEn_x)
         tmp2 = Yfac*(nH_y + yconduc(j)*nEn_y)
         cpfac = b_cp1 + b_cp2
         jnfac = 1.0-cpfac*Yfac
         fJN(j,1) = (fJN1(j,1)+fPN1(j,1)+fPN2(j,1)-cpfac*tmp1)/jnfac
         fJN(j,2) = (fJN1(j,2)+fPN1(j,2)+fPN2(j,2)-cpfac*tmp2)/jnfac
c     Forcing term
         f1 = tmp1-Yfac*fJN(j,1)
         f2 = tmp2-Yfac*fJN(j,2)
c     Drude term
         resfJN1(j,1) = -a_d*fJN1(j,1)+b_d*f1
         resfJN1(j,2) = -a_d*fJN1(j,2)+b_d*f2
c     First critical point term
         resfPN2(j,1) = fQN2(j,1)+b_11*f1
         resfPN2(j,2) = fQN2(j,2)+b_11*f2
	 resfQN2(j,1) = -a_211*fPN2(j,1)-a_221*fQN2(j,1)+b_21*f1
         resfQN2(j,2) = -a_211*fPN2(j,2)-a_221*fQN2(j,2)+b_21*f2
c     Second critical point term
         resfPN3(j,1) = fQN3(j,1)+b_12*f1
         resfPN3(j,2) = fQN3(j,2)+b_12*f2
	 resfQN3(j,1) = -a_212*fPN3(j,1)-a_222*fQN3(j,1)+b_22*f1
         resfQN3(j,2) = -a_212*fPN3(j,2)-a_222*fQN3(j,2)+b_22*f2
      enddo

      return
      end
c-----------------------------------------------------------------------
