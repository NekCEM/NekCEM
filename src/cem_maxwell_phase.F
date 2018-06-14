c-----------------------------------------------------------------------
      subroutine fftshift(u,N)
c     mimic MATLAB fftshift
c-----------------------------------------------------------------------
      implicit none
      integer N, N_mid, i
      real u(1)
      real tmp(N)

      N_mid = floor(real(N/2)) + 1

      do i=N_mid,N
        tmp(i-N_mid+1) = u(i)
      enddo
      do i=1,N_mid-1
        tmp(N-N_mid+1+i) = u(i)
      enddo
     
      call copy(u,tmp,N)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine Zfftshift(u,N)
c     mimic MATLAB fftshift
c-----------------------------------------------------------------------
      implicit none
      integer N, N_mid, i
      complex*16 u(1)
      complex*16 tmp(N)

      N_mid = floor(real(N/2)) + 1

      do i=N_mid,N
        tmp(i-N_mid+1) = u(i)
      enddo
      do i=1,N_mid-1
        tmp(N-N_mid+1+i) = u(i)
      enddo
     
      call Z_copy(u,tmp,N)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine Z_conj(u,N)
c-----------------------------------------------------------------------
      implicit none
      integer N
      complex*16 u(1)
      real tmp_r(N),tmp_i(N)

      call Z2R_copy(tmp_r,tmp_i,u,n)
      call chsign(tmp_i,n)
      call R2Z_copy(u,tmp_r,tmp_i,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine tmp_ZCOL2(a,b,N)
c-----------------------------------------------------------------------
      implicit none
      integer N,i
      complex*16 a(1),b(1)

      do i = 1,N
        a(i) = a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tmp_ZCOL3(a,b,c,N)
c-----------------------------------------------------------------------
      implicit none
      integer N,i
      complex*16 a(1),b(1),c(1)

      do i = 1,N
        a(i) = b(i)*c(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tmp_ZMULTR(a,const,N)
c-----------------------------------------------------------------------
      implicit none
      integer N,i
      complex*16 a(1)
      real const

      do i = 1,N
        a(i) = a(i)*const
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tmp_ZpowI3(a,b,const,N)
c-----------------------------------------------------------------------
      implicit none
      integer N,i
      complex*16 a(1),b(1)
      integer const

      do i = 1,N
        a(i) = b(i)**const
      enddo

      return
      end
c-----------------------------------------------------------------------
c     evaluate polynomial function f(z)
      subroutine Zpolyval(fz,z,ifconj,p,order,Np,N)
c     ifconj = -1: eval f(conj(z)) 
c            else: eval f(z) 
c     order  = -1: fz = sum( p(1:Np).*z.^(1:Np) )
c            else: fz = sum( p(1:Np).*z.^(Np:1) )
c
c     p      = [c_Np,...,c_2,c_1,c_0], of length Np+1
c-----------------------------------------------------------------------
      implicit none
      integer Np,N,i,j,j0,order,ifconj
      complex*16 p(Np),z(N),zw(N),fz(N),tmp
      real const

      if (Np.gt.100) then
        write(*,*)'WARNING: Np is large for Zpolyval'
      endif


      if (ifconj.eq.-1) then
        call Z_copy(zw,z,N)
        call Z_conj(zw,N)
      else
        call Z_copy(zw,z,N)
      endif
      call Z_zero(fz,N)


      if (order.eq.-1) then

        do i = 1,N
          tmp = complex(0.0,0.0)
          do j = 1,Np+1
            j0 = j-1  ! j0=0:N
            tmp = p(j)*( zw(i)**(j0) )
            fz(i) = fz(i) + tmp 
          enddo
        enddo

      else

        do i = 1,N
          tmp = complex(0.0,0.0)
          do j = 1,Np+1
            j0 = j-1  ! j0=0:N
            tmp = p(j)*( zw(i)**(Np-j0) )
            fz(i) = fz(i) + tmp 
          enddo
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine arg_vlmax(idx,vmax,v,N)
c     input: v, N
c     output: vmax, idx    such that vmax = v(idx) = max(v)
c-----------------------------------------------------------------------
      implicit none
      real v(1),vmax
      integer idx,N,i

      vmax = -99.9E20

      do i=1,N
        if (v(i).gt.vmax) then
          vmax = v(i)
          idx = i
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
