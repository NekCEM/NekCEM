c
c-----------------------------------------------------------------------
c
      subroutine ax_csr(u,v,a,ia,ja,n)
c
c     Compute u = Av, a in csr format, n rows
c
      real a(1),u(1),v(1),t
      integer ia(1),ja(1)
c
      do i=1,n
         j1 = ia(i)
         j2 = ia(i+1)-1
         t = 0.0
         do j=j1,j2
            t = t+a(j)*v(ja(j))
         enddo
         u(i) = t
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine atx_csr(u,v,a,ia,ja,n)
c
c     Compute u = A^T v, a in csr format, n rows
c
      real a(1),u(1),v(1)
      integer ia(1),ja(1)
c
      call rzero(u,n)
c
      do i=1,n
         j1 = ia(i)
         j2 = ia(i+1)-1
         t  = v(i)
         do j=j1,j2
            u(ja(j)) = u(ja(j)) + a(j)*t
         enddo
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
