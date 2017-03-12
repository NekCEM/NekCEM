c-----------------------------------------------------------------------
      subroutine face_deriv_2d(ut,u,f,e)
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'DXYZ'
      real     ut(lx1),u(lx1,ly1,lelt)
      real     s    
      integer  e,f

c
c     Face derivative:
c
c     Face numbering in reference element, Omega-hat
c
c            3
c        +---^----+
c        |  s|    |
c       4|   +->r |2
c        |        |
c        +--------+
c            1


      common /ctmp0/ ur(lx1),us(lx1),ul(lx1)

      if     (f.eq.3) then

         call mxm(dxm1 ,nx1,u(  1,ny1,e),nx1,ut,1)
         do i=1,lx1
            s = sqrt(rxm1(i,ny1,1,e)**2+rym1(i,ny1,1,e)**2)
            ut(i)=ut(i)*s
         enddo

      elseif (f.eq.1) then

         call mxm(dxm1 ,nx1,u(  1,  1,e),nx1,ul,1)
         do i=1,lx1
            s = -sqrt(rxm1(i,1,1,e)**2+rym1(i,1,1,e)**2)
            ut(1+lx1-i)=ul(i)*s
         enddo

      elseif (f.eq.2) then

         call mxm(u(nx1,  1,e),nx1,dytm1,nx1,ul,1)
         do i=1,lx1
            s = -sqrt(sxm1(nx1,i,1,e)**2+sym1(nx1,i,1,e)**2)
            ut(1+lx1-i)=ul(i)*s
         enddo

      elseif (f.eq.4) then

         call mxm(u(  1,  1,e),nx1,dytm1,nx1,ul,1)
         do i=1,lx1
            s = sqrt(sxm1(1,i,1,e)**2+sym1(1,i,1,e)**2)
            ut(i)=ul(i)*s
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
