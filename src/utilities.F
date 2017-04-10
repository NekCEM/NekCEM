c-----------------------------------------------------------------------
c     Various utility functions.
c-----------------------------------------------------------------------
      subroutine cem_error(field,sol,error,l2error,linferror)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real field(lpts1,3),sol(lpts1,3),error(lpts1,3)
      real l2error(3),linferror(3)

      integer ntott
      real errx,erry,errz
      real glsc3,glamax

      ntott = nx1*ny1*nz1*nelt

      call sub3(error(1,1),sol(1,1),field(1,1),ntott)
      call sub3(error(1,2),sol(1,2),field(1,2),ntott)
      call sub3(error(1,3),sol(1,3),field(1,3),ntott)
      
      errx = glsc3(error(1,1),bm1,error(1,1),ntott)/volvm1
      erry = glsc3(error(1,2),bm1,error(1,2),ntott)/volvm1
      errz = glsc3(error(1,3),bm1,error(1,3),ntott)/volvm1

      if (errx.gt.0.0) errx = sqrt(errx)
      if (erry.gt.0.0) erry = sqrt(erry)
      if (errz.gt.0.0) errz = sqrt(errz)

      l2error(1) = errx
      l2error(2) = erry
      l2error(3) = errz

      linferror(1) = glamax(error(1,1),ntott)
      linferror(2) = glamax(error(1,2),ntott)
      linferror(3) = glamax(error(1,3),ntott)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_errorchk(field,sol,error,l2tol,linftol)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real field(lpts1,3),sol(lpts1,3),error(lpts1,3)
      real l2tol(3),linftol(3)

      integer i
      real l2error(3),linferror(3)

      call cem_error(field,sol,error,l2error,linferror)
      do i = 1,3
         if (l2error(i).gt.l2tol(i)) stop 1
         if (linferror(i).gt.linftol(i)) stop 1
      enddo

      return
      end
c-----------------------------------------------------------------------
