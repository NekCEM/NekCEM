      subroutine mxmf2_gpu(A,N1,B,N2,C,N3)
c
c     unrolled loop version 
c
      real a(n1,n2),b(n2,n3),c(n1,n3)

#ifdef GPU
      call mxm_gpu(a,n1,b,n2,c,n3)
#else 
      if (n2.le.8) then
         if (n2.eq.1) then
            call mxf1(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxf2(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxf3(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxf4(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxf5(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxf6(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxf7(a,n1,b,n2,c,n3)
         else
            call mxf8(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxf9(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxf10(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxf11(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxf12(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxf13(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxf14(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxf15(a,n1,b,n2,c,n3)
         else
            call mxf16(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.24) then
         if (n2.eq.17) then
            call mxf17(a,n1,b,n2,c,n3)
         elseif (n2.eq.18) then
            call mxf18(a,n1,b,n2,c,n3)
         elseif (n2.eq.19) then
            call mxf19(a,n1,b,n2,c,n3)
         elseif (n2.eq.20) then
            call mxf20(a,n1,b,n2,c,n3)
         elseif (n2.eq.21) then
            call mxf21(a,n1,b,n2,c,n3)
         elseif (n2.eq.22) then
            call mxf22(a,n1,b,n2,c,n3)
         elseif (n2.eq.23) then
            call mxf23(a,n1,b,n2,c,n3)
         elseif (n2.eq.24) then
            call mxf24(a,n1,b,n2,c,n3)
         endif
      else
         call mxm44_0(a,n1,b,n2,c,n3)
      endif
#endif
c      do j=1,n3
c        do i=1,n1
c          if ((c(i,j)-g(i,j)) .ge. 1.00E-13) then
c            write (*,900) 'mxm:',i,j,c(i,j),g(i,j),c(i,j)-g(i,j)
c          endif
c        enddo
c      enddo

      return
      end
c-----------------------------------------------------------------------
