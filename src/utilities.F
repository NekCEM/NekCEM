c-----------------------------------------------------------------------
c     Various utility functions.
c-----------------------------------------------------------------------
      subroutine cem_error(u,exact,error,n,l2,linf)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'MASS'

      integer n
      real u(n),exact(n),error(n)
      real l2,linf

      real glsc3,glamax

      call sub3(error,exact,u,n)

      l2 = glsc3(error,bm1,error,n)/volvm1
      if (l2.gt.0.0) l2 = sqrt(l2)

      linf = glamax(error,n)

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
