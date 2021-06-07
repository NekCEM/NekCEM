c-----------------------------------------------------------------------
c
c     Box geometry with periodic boundary conditions.
c
c-----------------------------------------------------------------------
      subroutine userinc(tt,incfhx,incfhy,incfhz,incfex,incfey,incfez)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'

      real tt
      real incfhx(lxzfl),incfhy(lxzfl),incfhz(lxzfl)
      real incfex(lxzfl),incfey(lxzfl),incfez(lxzfl)

      return
      end
c-----------------------------------------------------------------------
      subroutine userini(tt,hx,hy,hz,ex,ey,ez)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real tt
      real hx(lpts),hy(lpts),hz(lpts)
      real ex(lpts),ey(lpts),ez(lpts)

      call usersol(tt,hx,hy,hz,ex,ey,ez)

      return
      end
c-----------------------------------------------------------------------
      subroutine usersol(tt,solhx,solhy,solhz,solex,soley,solez)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real tt
      real solhx(lpts),solhy(lpts),solhz(lpts)
      real solex(lpts),soley(lpts),solez(lpts)

      real omega,tmph,tmpe
      real xx,yy,zz
      integer i,n

      n = nx1*ny1*nz1*nelt

      omega = sqrt(2.0)

!$ACC DATA PRESENT(solhx,solhy,solhz,solex,soley,solez)
!$ACC&     PRESENT(xm1,ym1,zm1)

      if (iftm) then
         tmph = sin(omega*tt)/omega
         tmpe = cos(omega*tt)

!$ACC PARALLEL LOOP
         do i = 1,n
            xx = xm1(i,1,1,1)
            yy = ym1(i,1,1,1)
            zz = zm1(i,1,1,1)

            solhx(i) = cos(xx)*sin(yy)*tmph
            solhy(i) = -sin(xx)*cos(yy)*tmph
            solhz(i) = 0.0
            solex(i) = 0.0
            soley(i) = 0.0
            solez(i) = cos(xx)*cos(yy)*tmpe
         enddo
!$ACC END PARALLEL LOOP

      elseif (ifte) then
         tmph = cos(omega*tt)
         tmpe = sin(omega*tt)/omega

!$ACC PARALLEL LOOP
         do i = 1,n
            xx = xm1(i,1,1,1)
            yy = ym1(i,1,1,1)
            zz = zm1(i,1,1,1)

            solhx(i) = 0.0
            solhy(i) = 0.0
            solhz(i) = sin(xx)*sin(yy)*tmph
            solex(i) = sin(xx)*cos(yy)*tmpe
            soley(i) = -cos(xx)*sin(yy)*tmpe
            solez(i) = 0.0
         enddo
!$ACC END PARALLEL LOOP

      else
         write(*,*) 'ERROR: usersol: invalid imode'
         call exitt(1)
      endif

!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine usersrc(tt,srchx,srchy,srchz,srcex,srcey,srcez)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'

      real tt
      real srchx(lpts),srchy(lpts),srchz(lpts)
      real srcex(lpts),srcey(lpts),srcez(lpts)

      return
      end
c-----------------------------------------------------------------------
      subroutine userfsrc(tt,srcfhx,srcfhy,srcfhz,srcfex,srcfey,srcfez)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'

      real tt
      real srcfhx(lxzfl),srcfhy(lxzfl),srcfhz(lxzfl)
      real srcfex(lxzfl),srcfey(lxzfl),srcfez(lxzfl)

      return
      end
c-----------------------------------------------------------------------
      subroutine uservp(ix,iy,iz,iel)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

c     These don't do anything! This is a temporary measure until
c
c     https://github.com/NekCEM/NekCEM/issues/12
c
c     is resolved.
      integer ix,iy,iz,iel

      integer i

      do i = 1,npts
         permittivity(i) = 1.0
         permeability(i) = 1.0
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
c-----------------------------------------------------------------------
      implicit none

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer i,n
      real sx,sy,sz
      real glmin,glmax
      real xmin,ymin,zmin
      real xmax,ymax,zmax

      n = nx1*ny1*nz1*nelt

      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)

      sx = 2.0*pi/(xmax-xmin)
      sy = 2.0*pi/(ymax-ymin)

      do i = 1,n
         xm1(i,1,1,1) = sx*(xm1(i,1,1,1)-xmin)
         ym1(i,1,1,1) = sy*(ym1(i,1,1,1)-ymin)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RTIMER'

      integer i
      real l2(6),linf(6)
      real l2tol(6),linftol(6)

      if (ifte) then
         l2tol(1) = 0.0
         l2tol(2) = 0.0
         l2tol(3) = 5e-8
         l2tol(4) = 5e-8
         l2tol(5) = 5e-8
         l2tol(6) = 0.0

         linftol(1) = 0.0
         linftol(2) = 0.0
         linftol(3) = 5e-7
         linftol(4) = 5e-7
         linftol(5) = 5e-7
         linftol(6) = 0.0
      elseif (iftm) then
         l2tol(1) = 5e-8
         l2tol(2) = 5e-8
         l2tol(3) = 0.0
         l2tol(4) = 0.0
         l2tol(5) = 0.0
         l2tol(6) = 5e-8

         linftol(1) = 5e-7
         linftol(2) = 5e-7
         linftol(3) = 0.0
         linftol(4) = 0.0
         linftol(5) = 0.0
         linftol(6) = 5e-7
      else
         write(*,*) 'ERROR: userchk: invalid imode'
         call exitt(1)
      endif

      if (istep.le.10.or.mod(istep,iocomm).eq.0) then
         call usersol
     $     (time,shn(1,1),shn(1,2),shn(1,3),sen(1,1),sen(1,2),sen(1,3))

!$ACC UPDATE HOST(hn,en,shn,sen)
         call cem_error(hn(1,1),shn(1,1),errhn(1,1),npts,l2(1),linf(1))
         call cem_error(hn(1,2),shn(1,2),errhn(1,2),npts,l2(2),linf(2))
         call cem_error(hn(1,3),shn(1,3),errhn(1,3),npts,l2(3),linf(3))
         call cem_error(en(1,1),sen(1,1),erren(1,1),npts,l2(4),linf(4))
         call cem_error(en(1,2),sen(1,2),erren(1,2),npts,l2(5),linf(5))
         call cem_error(en(1,3),sen(1,3),erren(1,3),npts,l2(6),linf(6))

         call userprint(istep,time,dt,l2,linf,cpu_t,cpu_p_t)

         do i = 1,6
            if (l2(i).gt.l2tol(i)) call exitt(1)
            if (linf(i).gt.linftol(i)) call exitt(1)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userprint(istep,tt,dt,l2,linf,t1,t2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'

      integer istep
      real tt,dt,t1,t2
      real l2(6),linf(6)

      integer k

      if (nid.eq.0) then
         write(6,101) istep,nelt,nx1-1,npts,tt,dt,(l2(k),k=1,6),t1,t2
         write(6,102) istep,nelt,nx1-1,npts,tt,dt,(linf(k),k=1,6),t1,t2
      endif

 101  format(/,i10,i6,i4,i9,1p9e10.3,e9.2,' CPU: L2')
 102  format(  i10,i6,i4,i9,1p9e10.3,e9.2,' CPU: Linf')

      return
      end
c-----------------------------------------------------------------------

c automatically added by makenek
      subroutine userq(ix,iy,iz,eg)

      return
      end

c automatically added by makenek
      subroutine useric(ix,iy,iz,eg) 

      return
      end

c automatically added by makenek
      subroutine userbc(ix,iy,iz,iside,eg)

      return
      end

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrdat3 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end
