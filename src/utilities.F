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
