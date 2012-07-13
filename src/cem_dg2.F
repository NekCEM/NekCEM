!> \file cem_dg.F 
!!
!! \brief CEM main subroutines          
!!
!! A more detailed description of the file


!> \brief Brief description of the subroutine
!!
!! A more detailed descripton goes here. This is an example
!! of a detailed description of a method
!!
!! @param param1 a description of the first parameter
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
!! @see cem_op()
!! @see READATW()
!! @see RDMESH()
!! @see RDMESHW()
c--------------------------------------------------------------------- 
      subroutine cem (igeom)           
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      if   (IFSCHROD ) then

          if (ifrk ) call cem_schrod_op_rk
          if (ifeig) call cem_schrod_op_eig

      else
       if (param(5).lt.0) then

          if (ifrk ) call cem_maxwell_op_rk
       else
          if (ifeig) call cem_op_eig  
          if (ifexp) call cem_op_exp    
          if (ifrk ) call cem_op_rk
      
c         if (ifsim) call cem_op_sim
c         if (iffilter) call q_filter(0.05) 
C         if (ifwake  ) call cem_divergence
C         if (ifscat  ) call total_field    
C         if (param(14).gt.0) call cem_energy        
C         if (ifwake  ) call cem_wakepotential
       endif
      endif

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_op_rk                   
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'

c ... 5-Stage/4th-order RK
      if     (ifrk45) then
          do ii=1,5             
             call rk_c (ii)
             call cem_op   
             call rk_ab(ii)
          enddo
c ... Standard 2nd/3rd/4th-order RK steps
      !call  cem_op_rk44(ioption)
      elseif (ifrk22) then
             call cem_op_std_rk(2)
      elseif (ifrk33) then
             call cem_op_std_rk(3)
      elseif (ifrk44) then
             call cem_op_std_rk(4)
      endif

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_op_sim                  
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'

      common  /SIMSTP/ isimstep
      integer          isimstep
      integer  icalld
      save     icalld
      data     icalld/0/

c ... Simplectic steps setting
      if (icalld.eq.0) then
          if (ifsim1) isimstep=1
          if (ifsim2) isimstep=2
          if (ifsim3) isimstep=3
          if (ifsim4) isimstep=4
          if (ifsim5) isimstep=5
          icalld=1
      endif

      call cem_op_std_sim(isimstep)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_op_std_rk(irkstep)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'
      common /STEPS/ resWN (6*lx1*ly1*lz1*lelt,lw)

      npts = nx1*ny1*nz1*nelt
      n6   = 3*npts                
      if (imode.eq.3) n6= 6*npts               

c ... RK Standard 4th order 
      if     (irkstep.eq.4) then

      call rk_save    (reswn(1,5),hx,hy,hz,ex,ey,ez,1.0,npts,n6,imode) !un=reswn
      call cem_op                         !f(tn,un)
      call rk_save    (reswn(1,1),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)    !K1=dt*f(tn,un)
      call rk_substeps(reswn(1,1),reswn(1,5),hx,hy,hz,ex,ey,ez,
     $                 0.5,npts,n6,imode) !  un+0.5*K1
      call cem_op                         !f(un+0.5*K1) 
      call rk_save    (reswn(1,2),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)    !K2=dt*f(un+0.5*K1)
      call rk_substeps(reswn(1,2),reswn(1,5),hx,hy,hz,ex,ey,ez,
     $                 0.5,npts,n6,imode) !  un+0.5*K2
      call cem_op                         !f(un+0.5*K2)
      call rk_save    (reswn(1,3),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)    !K3=dt*f(un+0.5*K2)
      call rk_substeps(reswn(1,3),reswn(1,5),hx,hy,hz,ex,ey,ez,
     $                 1.0,npts,n6,imode) !  un+K3 
      call cem_op     !f(un+K3)
      call rk_save    (reswn(1,4),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)    !K4=dt*f(un+K3) 
 
      call rk_fourth  (reswn(1,1),reswn(1,2),reswn(1,3),reswn(1,4),
     $                 reswn(1,5),hx,hy,hz,ex,ey,ez,npts,n6,imode)
                      !u'= un+(K1+2*K2+2*K3+K4)/6

c ... RK Standard 3rd order 
      elseif (irkstep.eq.3) then

      call rk_save    (reswn(1,5),hx,hy,hz,ex,ey,ez,1.0,npts,n6,imode)
      call cem_op     !K1/h
      call rk_save    (reswn(1,1),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)
      call rk_substeps(reswn(1,1),reswn(1,5),hx,hy,hz,ex,ey,ez,0.5,
     $                 npts,n6,imode)
      call cem_op     !K2/h
      call rk_save    (reswn(1,2),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)
      call rk_save2   (reswn(1,3),reswn(1,1),reswn(1,2),
     $                 -1.0,2.0,npts,n6,imode)
      call rk_substeps(reswn(1,3),reswn(1,5),hx,hy,hz,ex,ey,ez,
     $                 1.0,npts,n6,imode)
      call cem_op     !K3/h
      call rk_save    (reswn(1,4),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)
      call rk_third   (reswn(1,1),reswn(1,2),reswn(1,4),reswn(1,5),
     $                 hx,hy,hz,ex,ey,ez,npts,n6,imode)

c ... RK Standard 2nd order 
      elseif (irkstep.eq.2) then

      call rk_save    (reswn(1,3),hx,hy,hz,ex,ey,ez,1.0,npts,n6,imode)
      call cem_op     !K1/h
      call rk_save    (reswn(1,1),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)
      call rk_substeps(reswn(1,1),reswn(1,3),hx,hy,hz,ex,ey,ez,0.5,
     $                 npts,n6,imode)
      call cem_op     !K2/h
      call rk_save    (reswn(1,2),reshx,reshy,reshz,resex,resey,resez,
     $                 dt,npts,n6,imode)
      call rk_second  (reswn(1,2),reswn(1,3),hx,hy,hz,ex,ey,ez,
     $                 npts,n6,imode) 
                      !u'=un+K2                 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_op_std_sim(isimstep)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'
      common /STEPS/ resWN (6*lx1*ly1*lz1*lelt,lw)

c     real    a(5),b(5)
c
c     if     (isimstep.eq.1) then
c         a(1)= 1.0
c         b(1)= 1.0
c     elseif (isimstep.eq.2) then
c         a(1)= 0.5     
c         a(2)= a(1) 
c         b(1)= 0.0  
c         b(2)= 1.0  
c     elseif (isimstep.eq.3) then
c         a(1)= 2.0/3.0
c         a(2)=-a(1) 
c         a(3)= 1.0  
c         b(1)= 7.0/24.0
c         b(2)= 3.0/4.0
c         b(3)=-1.0/24.0
c     elseif (isimstep.eq.4) then
c         a(1)= 1.0/6.0*(2.0+2.0**(1/3)+2.0**(-1/3))
c         a(2)= 1.0/6.0*(1.0-2.0**(1/3)-2.0**(-1/3))
c         a(3)= a(2)                                   
c         a(4)= a(1)                                  
c         b(1)= 0.0           
c         b(2)= 1.0/(2.0-2**(1/3))
c         b(3)= 1.0/(1.0-2**(2/3))
c         b(4)= b(2)                    
c     endif

      npts = nx1*ny1*nz1*nelt
      n6   = 3*npts                
      if (imode.eq.3) n6= 6*npts               

      call rk_save  (reswn(1,2),hx,hy,hz,ex,ey,ez,1.0,npts,n6,3)
      call cem_op    !K1/h
      call rk_save  (reswn(1,1),reshx,reshy,reshz,resex,resey,resez,
     $               dt,npts,n6,3)
      call rk_substeps(reswn(1,1),reswn(1,2),hx,hy,hz,ex,ey,ez,
     $                 1.0,npts,n6,3)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_op_rk44(ioption)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'

      common /STEPS/ resWN (6*lx1*ly1*lz1*lelt,lw)

      a1=0.5*dt
      b1=0.5*dt
      c1=1.0*dt
      const=1.0/6.0

      npts = nx1*ny1*nz1*nelt
      n6   = 3*npts                

      if (imode.eq.3) n6= 6*npts               

      do i=1,npts
         reswn(i+0*npts,5)=hx(i,1,1,1)
         reswn(i+1*npts,5)=hy(i,1,1,1)
         reswn(i+2*npts,5)=hz(i,1,1,1)
         reswn(i+3*npts,5)=ex(i,1,1,1)
         reswn(i+4*npts,5)=ey(i,1,1,1)
         reswn(i+5*npts,5)=ez(i,1,1,1)
      enddo

      call cem_op !K1/h

      do i=1,npts
         reswn(i+0*npts,1)=dt*reshx(i,1,1,1)
         reswn(i+1*npts,1)=dt*reshy(i,1,1,1)
         reswn(i+2*npts,1)=dt*reshz(i,1,1,1)
         reswn(i+3*npts,1)=dt*resex(i,1,1,1)
         reswn(i+4*npts,1)=dt*resey(i,1,1,1)
         reswn(i+5*npts,1)=dt*resez(i,1,1,1)
 
         hx(i,1,1,1) = reswn(i+0*npts,5)+0.5*reswn(i+0*npts,1)
         hy(i,1,1,1) = reswn(i+1*npts,5)+0.5*reswn(i+1*npts,1)
         hz(i,1,1,1) = reswn(i+2*npts,5)+0.5*reswn(i+2*npts,1)
         ex(i,1,1,1) = reswn(i+3*npts,5)+0.5*reswn(i+3*npts,1)
         ey(i,1,1,1) = reswn(i+4*npts,5)+0.5*reswn(i+4*npts,1)
         ez(i,1,1,1) = reswn(i+5*npts,5)+0.5*reswn(i+5*npts,1)
      enddo

      call cem_op !K2/h

      do i=1,npts
         reswn(i+0*npts,2)=dt*reshx(i,1,1,1)
         reswn(i+1*npts,2)=dt*reshy(i,1,1,1)
         reswn(i+2*npts,2)=dt*reshz(i,1,1,1)
         reswn(i+3*npts,2)=dt*resex(i,1,1,1)
         reswn(i+4*npts,2)=dt*resey(i,1,1,1)
         reswn(i+5*npts,2)=dt*resez(i,1,1,1)

         hx(i,1,1,1) = reswn(i+0*npts,5)+0.5*reswn(i+0*npts,2)
         hy(i,1,1,1) = reswn(i+1*npts,5)+0.5*reswn(i+1*npts,2)
         hz(i,1,1,1) = reswn(i+2*npts,5)+0.5*reswn(i+2*npts,2)
         ex(i,1,1,1) = reswn(i+3*npts,5)+0.5*reswn(i+3*npts,2)
         ey(i,1,1,1) = reswn(i+4*npts,5)+0.5*reswn(i+4*npts,2)
         ez(i,1,1,1) = reswn(i+5*npts,5)+0.5*reswn(i+5*npts,2)
      enddo

      call cem_op !K3/h

      do i=1,npts
         reswn(i+0*npts,3)=dt*reshx(i,1,1,1)
         reswn(i+1*npts,3)=dt*reshy(i,1,1,1)
         reswn(i+2*npts,3)=dt*reshz(i,1,1,1)
         reswn(i+3*npts,3)=dt*resex(i,1,1,1)
         reswn(i+4*npts,3)=dt*resey(i,1,1,1)
         reswn(i+5*npts,3)=dt*resez(i,1,1,1)

         hx(i,1,1,1) = reswn(i+0*npts,5)+reswn(i+0*npts,3)
         hy(i,1,1,1) = reswn(i+1*npts,5)+reswn(i+1*npts,3)
         hz(i,1,1,1) = reswn(i+2*npts,5)+reswn(i+2*npts,3)
         ex(i,1,1,1) = reswn(i+3*npts,5)+reswn(i+3*npts,3)
         ey(i,1,1,1) = reswn(i+4*npts,5)+reswn(i+4*npts,3)
         ez(i,1,1,1) = reswn(i+5*npts,5)+reswn(i+5*npts,3)
      enddo

      call cem_op !K4/h

      do i=1,npts
         reswn(i+0*npts,4)=dt*reshx(i,1,1,1)
         reswn(i+1*npts,4)=dt*reshy(i,1,1,1)
         reswn(i+2*npts,4)=dt*reshz(i,1,1,1)
         reswn(i+3*npts,4)=dt*resex(i,1,1,1)
         reswn(i+4*npts,4)=dt*resey(i,1,1,1)
         reswn(i+5*npts,4)=dt*resez(i,1,1,1)
      enddo

      do i=1,npts
       hx(i,1,1,1) = reswn(i+0*npts,5)
     $  +const*(  reswn(i+0*npts,1)+2*reswn(i+0*npts,2)
     $         +2*reswn(i+0*npts,3)+  reswn(i+0*npts,4))
       hy(i,1,1,1) = reswn(i+1*npts,5)
     $  +const*(  reswn(i+1*npts,1)+2*reswn(i+1*npts,2)
     $         +2*reswn(i+1*npts,3)+  reswn(i+1*npts,4))
       hz(i,1,1,1) = reswn(i+2*npts,5)
     $  +const*(  reswn(i+2*npts,1)+2*reswn(i+2*npts,2)
     $         +2*reswn(i+2*npts,3)+  reswn(i+2*npts,4))
       ex(i,1,1,1) = reswn(i+3*npts,5)
     $  +const*(  reswn(i+3*npts,1)+2*reswn(i+3*npts,2)
     $         +2*reswn(i+3*npts,3)+  reswn(i+3*npts,4))
       ey(i,1,1,1) = reswn(i+4*npts,5)
     $  +const*(  reswn(i+4*npts,1)+2*reswn(i+4*npts,2)
     $         +2*reswn(i+4*npts,3)+  reswn(i+4*npts,4))
       ez(i,1,1,1) = reswn(i+5*npts,5)
     $  +const*(  reswn(i+5*npts,1)+2*reswn(i+5*npts,2)
     $         +2*reswn(i+5*npts,3)+  reswn(i+5*npts,4))

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_save(u,ux1,uy1,uz1,ux2,uy2,uz2,a,n,n6,imode)
c     include 'TOTAL'
      integer  n,n6,imode
      real     a
      real     u(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)

      if      (imode.eq. 3) then    !if3D

       do i=1,n    
         u(i+0*n)=a*ux1(i)
         u(i+1*n)=a*uy1(i)
         u(i+2*n)=a*uz1(i)
         u(i+3*n)=a*ux2(i)
         u(i+4*n)=a*uy2(i)
         u(i+5*n)=a*uz2(i)
       enddo

      elseif  (imode.eq. 2) then    !ifTM

       do i=1,n    
         u(i+0*n)=a*ux1(i)
         u(i+1*n)=a*uy1(i)
         u(i+2*n)=a*uz2(i)
       enddo

      elseif  (imode.eq. 1) then    !ifTM

       do i=1,n    
         u(i+0*n)=a*ux2(i)
         u(i+1*n)=a*uy2(i)
         u(i+2*n)=a*uz1(i)
       enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_save2(u3,u1,u2,a,b,n,n6,imode)
c     include 'TOTAL'
      integer  n,n6,imode
      real     a,b
      real     u1(n6),u2(n6),u3(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)

      if      (imode.eq. 3) then    !if3D

       do i=1,n    
         u3(i+0*n)=a*u1(i+0*n)+b*u2(i+0*n)
         u3(i+1*n)=a*u1(i+1*n)+b*u2(i+1*n)
         u3(i+2*n)=a*u1(i+2*n)+b*u2(i+2*n)
         u3(i+3*n)=a*u1(i+3*n)+b*u2(i+3*n)
         u3(i+4*n)=a*u1(i+4*n)+b*u2(i+4*n)
         u3(i+5*n)=a*u1(i+5*n)+b*u2(i+5*n)
       enddo

      elseif  (imode.eq. 2) then    !ifTM

       do i=1,n    
         u3(i+0*n)=a*u1(i+0*n)+b*u2(i+0*n)
         u3(i+1*n)=a*u1(i+1*n)+b*u2(i+1*n)
         u3(i+2*n)=a*u1(i+2*n)+b*u2(i+2*n)
       enddo

      elseif  (imode.eq. 1) then    !ifTM

       do i=1,n    
         u3(i+0*n)=a*u1(i+0*n)+b*u2(i+0*n)
         u3(i+1*n)=a*u1(i+1*n)+b*u2(i+1*n)
         u3(i+2*n)=a*u1(i+2*n)+b*u2(i+2*n)
       enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_substeps(u1,u2,ux1,uy1,uz1,ux2,uy2,uz2,b,n,n6,imode)
c     include 'TOTAL'
      real     b
      integer  n,n6,imode
      real     u1(n6)
      real     u2(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)

      if      (imode.eq. 3) then    !if3D

       do i=1,n    
         ux1(i) = u2(i+0*n)+b*u1(i+0*n)
         uy1(i) = u2(i+1*n)+b*u1(i+1*n)
         uz1(i) = u2(i+2*n)+b*u1(i+2*n)
         ux2(i) = u2(i+3*n)+b*u1(i+3*n)
         uy2(i) = u2(i+4*n)+b*u1(i+4*n)
         uz2(i) = u2(i+5*n)+b*u1(i+5*n)
       enddo

      elseif  (imode.eq. 2) then    !ifTM

       do i=1,n    
         ux1(i) = u2(i+0*n)+b*u1(i+0*n)
         uy1(i) = u2(i+1*n)+b*u1(i+1*n)
         uz2(i) = u2(i+2*n)+b*u1(i+2*n)
       enddo

      elseif  (imode.eq. 1) then    !ifTE

       do i=1,n    
         ux2(i) = u2(i+0*n)+b*u1(i+0*n)
         uy2(i) = u2(i+1*n)+b*u1(i+1*n)
         uz1(i) = u2(i+2*n)+b*u1(i+2*n)
       enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_second(u1,u2,ux1,uy1,uz1,ux2,uy2,uz2,n,n6,imode)
      integer  n,n6,imode
      real     u1(n6),u2(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)

      if     (imode.eq.3) then    !if3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          i3=i+3*n
          i4=i+4*n
          i5=i+5*n
          ux1(i)=u1(i0)+u2(i0)          
          uy1(i)=u1(i1)+u2(i1)
          uz1(i)=u1(i2)+u2(i2)
          ux2(i)=u1(i3)+u2(i3)
          uy2(i)=u1(i4)+u2(i4)
          uz2(i)=u1(i5)+u2(i5)
        enddo
      elseif (imode.eq.2) then    !if3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          ux1(i)=u1(i0)+u2(i0)
          uy1(i)=u1(i1)+u2(i1)
          uz2(i)=u1(i2)+u2(i2)
        enddo
      elseif (imode.eq.1) then    !if3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          ux2(i)=u1(i0)+u2(i0)
          uy2(i)=u1(i1)+u2(i1)
          uz1(i)=u1(i2)+u2(i2)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine 
     $     rk_third(u1,u2,u3,u5,ux1,uy1,uz1,ux2,uy2,uz2,n,n6,imode)
c     include 'TOTAL'
      integer  n,n6,m,imode
      real     u1(n6)
      real     u2(n6)
      real     u3(n6)
      real     u5(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)
      real     const

      cst=1.0/6.0

      if     (imode.eq.3) then    !if3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          i3=i+3*n
          i4=i+4*n
          i5=i+5*n
          ux1(i)=u5(i0)+cst*(u1(i0)+4*u2(i0)+u3(i0))
          uy1(i)=u5(i1)+cst*(u1(i1)+4*u2(i1)+u3(i1))
          uz1(i)=u5(i2)+cst*(u1(i2)+4*u2(i2)+u3(i2))
          ux2(i)=u5(i3)+cst*(u1(i3)+4*u2(i3)+u3(i3))
          uy2(i)=u5(i4)+cst*(u1(i4)+4*u2(i4)+u3(i4))
          uz2(i)=u5(i5)+cst*(u1(i5)+4*u2(i5)+u3(i5))
        enddo
      elseif (imode.eq.2) then    !if3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          ux1(i)=u5(i0)+cst*(u1(i0)+4*u2(i0)+u3(i0))
          uy1(i)=u5(i1)+cst*(u1(i1)+4*u2(i1)+u3(i1))
          uz2(i)=u5(i2)+cst*(u1(i2)+4*u2(i2)+u3(i2))
        enddo
      elseif (imode.eq.1) then    !if3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          ux2(i)=u5(i0)+cst*(u1(i0)+4*u2(i0)+u3(i0))
          uy2(i)=u5(i1)+cst*(u1(i1)+4*u2(i1)+u3(i1))
          uz1(i)=u5(i2)+cst*(u1(i2)+4*u2(i2)+u3(i2))
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine 
     $     rk_fourth(u1,u2,u3,u4,u5,ux1,uy1,uz1,ux2,uy2,uz2,n,n6,imode)
c     include 'TOTAL'
      integer n,n6,imode
      real    u1(n6),u2(n6),u3(n6)
      real    u4(n6),u5(n6)
      real    ux1(n),uy1(n),uz1(n)
      real    ux2(n),uy2(n),uz2(n)
      real    cst

      cst=1.0/6.0

      if     (imode.eq.3) then    !IF3D
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          i3=i+3*n
          i4=i+4*n
          i5=i+5*n
          ux1(i)=u5(i0)+cst*(u1(i0)+2*u2(i0)+2*u3(i0)+u4(i0))
          uy1(i)=u5(i1)+cst*(u1(i1)+2*u2(i1)+2*u3(i1)+u4(i1))
          uz1(i)=u5(i2)+cst*(u1(i2)+2*u2(i2)+2*u3(i2)+u4(i2))
          ux2(i)=u5(i3)+cst*(u1(i3)+2*u2(i3)+2*u3(i3)+u4(i3))
          uy2(i)=u5(i4)+cst*(u1(i4)+2*u2(i4)+2*u3(i4)+u4(i4))
          uz2(i)=u5(i5)+cst*(u1(i5)+2*u2(i5)+2*u3(i5)+u4(i5))
        enddo
      elseif (imode.eq.2) then    !IFTM
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          ux1(i)=u5(i0)+cst*(u1(i0)+2*u2(i0)+2*u3(i0)+u4(i0))
          uy1(i)=u5(i1)+cst*(u1(i1)+2*u2(i1)+2*u3(i1)+u4(i1))
          uz2(i)=u5(i2)+cst*(u1(i2)+2*u2(i2)+2*u3(i2)+u4(i2))
        enddo
      elseif (imode.eq.1) then    !IFTE
        do i=1,n    
          i0=i+0*n
          i1=i+1*n
          i2=i+2*n
          ux2(i)=u5(i0)+cst*(u1(i0)+2*u2(i0)+2*u3(i0)+u4(i0))
          uy2(i)=u5(i1)+cst*(u1(i1)+2*u2(i1)+2*u3(i1)+u4(i1))
          uz1(i)=u5(i2)+cst*(u1(i2)+2*u2(i2)+2*u3(i2)+u4(i2))
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_op_eig
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      call cem_eig   !call cem_eig_arpack

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_eig
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

#if  !defined(CEM_EIG_MEMORY_SMALL) &&    ! defined(CEM_EIG_MEMORY_BIG)
      parameter  (lelm = 1, lq1= 1  )     ! For IFRK4 runs
      parameter  (mm =lq1*lq1*6*lelm)     ! For IFEIG runs: 2D
      parameter  (mm2=mm*mm         )     ! For IFEIG runs: 2D
#endif

#ifdef CEM_EIG_MEMORY_SMALL
     !parameter  (lelm = 18,lq1=lx1 )     ! For IFEIG runs (make sure lelm > nelv for ev test)
      parameter  (lelm =  1, lq1=  1)     ! For IFEIG runs (make sure lelm > nelv for ev test)
      parameter  (mm =  1)!q1*lq1*6*lelm)     ! For IFEIG runs: 2D
      parameter  (mm2=  1)!mm*mm         )     ! For IFEIG runs: 2D
#endif

#ifdef CEM_EIG_MEMORY_BIG
     !parameter  (lelm = 3,lq1=lx1)       ! For IFEIG runs  3D (make sure lelm > nelv for ev test)
     !parameter  (lelm = 1,lq1=lx1)       ! For IFEIG runs  3D (make sure lelm > nelv for ev test)
      parameter  (lelm = 1, lq1=  1)     ! For IFEIG runs (make sure lelm > nelv for ev test)
      parameter  (mm =lq1*lq1*lq1*6*lelm) ! For IFEIG runs: 3D
      parameter  (mm2=mm*mm*mm)           ! For IFEIG runs: 3D
#endif
                                    
      common /ecem/   res( mm2)
     $            ,   ev ( mm2)
     $            ,   er ( mm )
     $            ,   ei ( mm )
     $            ,   rw1(4*mm)
      common /scrns/  srflx(6*2*ldim*lx1*lz1*lelt)

      real            res, evecres, evalresr, evalresi,rw1        

      integer         info,ii,i,jj,j
      integer         k,k1,n,n3,n6,nn,nn2   
      integer         icalld
      save            icalld
      data            icalld /0/


      n   = nx1*ny1*nz1*nelv
      n3  = nx1*ny1*nz1*nelv*3
      n6  = nx1*ny1*nz1*nelv*6

      if (n.gt.mm) then
         write(6,*) 'insufficient memory in cem_eig, abort',mm,n
         call exitt
      endif
      
      if (if3D)  then
          nn  = n6
      else
          nn  = n3
      endif

      nn2 = nn*nn 

      call rzero(res,nn2) 

      if     (imode.eq.3) then 

         k1 = 0

         do jj=1,6                  
         do j =1,n

         k1 =  k1+1
         k  = (k1-1)*nn+1  

         call rzero(hx,n) 
         call rzero(hy,n) 
         call rzero(hz,n) 
         call rzero(ex,n) 
         call rzero(ey,n) 
         call rzero(ez,n) 

         call rzero(reshx,n) 
         call rzero(reshy,n) 
         call rzero(reshz,n) 
         call rzero(resex,n) 
         call rzero(resey,n) 
         call rzero(resez,n) 

         if (jj.eq.1) hx(j,1,1,1) = 1.  ! jth col of identity matrix
         if (jj.eq.2) hy(j,1,1,1) = 1.
         if (jj.eq.3) hz(j,1,1,1) = 1.
         if (jj.eq.4) ex(j,1,1,1) = 1.                         
         if (jj.eq.5) ey(j,1,1,1) = 1.
         if (jj.eq.6) ez(j,1,1,1) = 1.

         call cem_curl
         call cem_restrict_to_face
         call cem_flux (srflx)
         call cem_add_flux_to_res (srflx)

         do i=1,n

            if (abs(reshx(i,1,1,1)).lt.1.e-13) reshx(i,1,1,1)=0.
            if (abs(reshy(i,1,1,1)).lt.1.e-13) reshy(i,1,1,1)=0.
            if (abs(reshz(i,1,1,1)).lt.1.e-13) reshz(i,1,1,1)=0.
            if (abs(resex(i,1,1,1)).lt.1.e-13) resex(i,1,1,1)=0.
            if (abs(resey(i,1,1,1)).lt.1.e-13) resey(i,1,1,1)=0.
            if (abs(resez(i,1,1,1)).lt.1.e-13) resez(i,1,1,1)=0.

         enddo

         call copy(res(k+0*n),reshx,n)
         call copy(res(k+1*n),reshy,n)
         call copy(res(k+2*n),reshz,n)
         call copy(res(k+3*n),resex,n)
         call copy(res(k+4*n),resey,n)
         call copy(res(k+5*n),resez,n)

         enddo
         enddo

      elseif (imode.eq.2) then

         k1 = 0

         do jj=1,3                  

         do j =1,n

         k1 =  k1+1
         k  = (k1-1)*nn+1  

         call rzero(hx,n) 
         call rzero(hy,n) 
         call rzero(ez,n) 

         call rzero(reshx,n) 
         call rzero(reshy,n) 
         call rzero(resez,n) 

         if (jj.eq.1) hx(j,1,1,1) = 1.  ! jth col of identity matrix
         if (jj.eq.2) hy(j,1,1,1) = 1.
         if (jj.eq.3) ez(j,1,1,1) = 1.

         call cem_curl
         call cem_restrict_to_face
         call cem_flux (srflx)
         call cem_add_flux_to_res (srflx)
            
         do i=1,n
            if (abs(reshx(i,1,1,1)).lt.1.e-13) reshx(i,1,1,1)=0.
            if (abs(reshy(i,1,1,1)).lt.1.e-13) reshy(i,1,1,1)=0.
            if (abs(resez(i,1,1,1)).lt.1.e-13) resez(i,1,1,1)=0.
         enddo
          
         call copy(res(k+0*n),reshx,n)
         call copy(res(k+1*n),reshy,n)
         call copy(res(k+2*n),resez,n)
          
         enddo
         enddo
           
      elseif (imode.eq.1) then

         k1 = 0

         do jj=1,3                  
         do j =1,n

         k1 =  k1+1
         k  = (k1-1)*nn+1  

         call rzero(ex,n) 
         call rzero(ey,n) 
         call rzero(hz,n) 

         call rzero(resex,n) 
         call rzero(resey,n) 
         call rzero(reshz,n) 

         if (jj.eq.1) ex(j,1,1,1) = 1.
         if (jj.eq.2) ey(j,1,1,1) = 1.  ! jth col of identity matrix
         if (jj.eq.3) hz(j,1,1,1) = 1.

         call cem_curl
         call cem_restrict_to_face
         call cem_flux (srflx)
         call cem_add_flux_to_res (srflx)

         do i=1,n
            if (abs(resex(i,1,1,1)).lt.1.e-13) resex(i,1,1,1)=0.
            if (abs(resey(i,1,1,1)).lt.1.e-13) resey(i,1,1,1)=0.
            if (abs(reshz(i,1,1,1)).lt.1.e-13) reshz(i,1,1,1)=0.
         enddo

         call copy (res(k+0*n),resex,n)
         call copy (res(k+1*n),resey,n)
         call copy (res(k+2*n),reshz,n)

         enddo
         enddo

      endif 

      if (nelv.eq.1) call outmat(res,nx1*3,nx1*3,'res   ',1)

      do i=1,nn2
         write(60,*) res(i)                  
      enddo

      write(6,*) '                                      '
      write(6,*) 'IFEIG T :: computing eigenvalues      '
      write(6,*) '                                      '

      call dgeev('N','V',nn,res,nn,er,ei,ev,nn,ev,nn,rw1,4*nn,info)

      do i=1,nn
         write(59,*)  er(i),ei(i)                  
      enddo

    1 format(1p1e19.12)

      if (icalld.eq.0) then
        write(6,*) ' IF3D =',if3d,', TM=',iftm, ', CENTRAL=',
     $       ifcentral ,', PEC=',ifpec,', RK4=',IFRK4, 
     $       ', DEALIAS=',ifdealias,', SCAT=',ifscat,', IFNM =',ifnm
        write(6,*) ' IFSRC=', ifsrc,', TE=', ifte,', UPWIND =', 
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude,
     $       ', LORENTZ=',ifloren

        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=', 
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter 
        icalld= 1     
      endif

      write(6,*) 'eigenvalue computation successfully completed!!!'
      call exitt


      return
      end

c-----------------------------------------------------------------------
      subroutine vec_energy2(a,b,c,d,n)
      real     a(1),b(1),c(1)
      real     d
      integer  n

      do i=1,n
         a(i)= a(i)+d*(b(i)**2+c(i)**2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine vec_energy3(a,b,c,d,e,n)
      real     a(1),b(1),c(1),d(1)
      real     e
      integer  n
      
      do i = 1,n
         a(i)= a(i)+e*(b(i)**2+c(i)**2+d(i)**2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_eig_arpack
C----------------------------------------------------------------------
#if USE_ARPACK
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      ! calculation temporaries
      real       srflx(6*2*ldim*lx1*lz1*lelt)

      ! problem size parameters
      integer    ldofperfld, ldof
      parameter (ldofperfld   = lx1*ly1*lz1*lelv)
      parameter (ldof = ldofperfld*6)

      integer    ndofperfld, ndof

      integer    eigvcount,ritzvcount, lworkl
      parameter (eigvcount = 10)
      parameter (ritzvcount = eigvcount*2+20)
      parameter (lworkl = 3*ritzvcount*ritzvcount + 6*ritzvcount)

      ! spectral shift parameters
      real       shift
      parameter (shift = 1.)

      ! ARPACK dnaupd arguments
      character  abmat*1, awhich*2
      integer    aido, ainfo, aldv, alworkl, an, ancv, anev
      real       atol

      integer    aiparam(11), aipntr(14)
      real       aresid(ldof), av(ldof,ritzvcount)
      real       aworkd(3*ldof), aworkl(lworkl)

      ! ARPACK dneupd extra arguments
      logical   arvec
      character ahowmny*1
      logical   aselect(ritzvcount)
      real      adr(eigvcount+1)
      real      adi(eigvcount+1)
      real      az(ldof, eigvcount+1)
      integer   aldz
      real      asigmar, asigmai
      real      aworkev(3*ritzvcount)

      ! loop variable
      integer   i, loopcnt

      ! start of code

      ndofperfld = nx1*ny1*nz1*nelv

      if (if3d) then
          ndof = 6*ndofperfld
      else
          ndof = 3*ndofperfld
      endif

      write(*,*) 'entering arpack eigenvalue computation'
      write(*,*) ndofperfld, ' dofs per field'
      write(*,*) ndof, ' dofs total'

      ! set parameters for arpack
      aido    =  0    ! initial call to RCI
      abmat   = 'I'   ! standard eigenvalue problem
      an      =  ndof
      awhich  = 'LR'  ! want evalues with Largest Real part
      anev    =  eigvcount ! how many eigenvalues
      atol    =  1e-5 ! permitted relative error
      ainfo   =  0    ! use random initial residual
      ancv    =  ritzvcount ! how many ritz vectors (.ge. arpnev+2)
      aldv    =  ldof ! leading dimension of arpv

      aiparam(1) = 1  ! use exact shifts
      aiparam(3) = 100 * ritzvcount ! maximum # of iterations
      aiparam(4) = 1
      aiparam(7) = 1  ! mode 1: Ax = lambda x

      alworkl = lworkl

      loopcnt = 0
      do while (aido.ne.99)
        loopcnt = loopcnt + 1
        if (mod(loopcnt,100).eq.0) then
          write (*,*) 'arpack iteration', loopcnt
        endif
        call dnaupd
     &  (aido, abmat, an, awhich, anev, atol, aresid, 
     &  ancv, av, aldv, aiparam, aipntr,
     &  aworkd, aworkl, alworkl, ainfo)

        if (ainfo.ne.0) then
          write (*,*) 'arpack dnaupd error number:', ainfo
          write (*,*) 'converged eigenvalues     :', aiparam(5)
          call exitt
        endif

        if ((aido.ne.99).and.((aido.eq.1).or.(aido.eq.-1))) then
          call rzero(reshx, ndofperfld) 
          call rzero(reshy, ndofperfld) 
          call rzero(reshz, ndofperfld) 
          call rzero(resex, ndofperfld) 
          call rzero(resey, ndofperfld) 
          call rzero(resez, ndofperfld) 

          call scatter_fields(hx, hy, hz, ex, ey, ez, aworkd(aipntr(1)))

          call cem_curl
          call cem_restrict_to_face
          call cem_flux (srflx)
          call cem_add_flux_to_res (srflx)

          call gather_fields(aworkd(aipntr(2)), 
     &         reshx, reshy, reshz, resex, resey, resez)
          if (shift.ne.0) then
              do i=0,ndof-1
                 aworkd(aipntr(2)+i) = aworkd(aipntr(2)+i) +
     &                                 shift*aworkd(aipntr(1)+i)
              enddo
          endif
        endif
      enddo

      arvec = .false. ! no eigenvectors
      ahowmny = 'R' ! Ritz vectors, but doesn't matter (see arvec)
      aldz = ldof

      call dneupd
     &(arvec, ahowmny, aselect, adr, adi, az, aldz, asigmar, asigmai, 
     & aworkev, abmat,
     & an, awhich, anev, atol, aresid, ancv, av, aldv, aiparam, aipntr, 
     & aworkd, aworkl, alworkl, ainfo)
      if (ainfo.ne.0) then
        write (*,*) 'arpack dneupd error number:', ainfo
        call exitt
      endif

      do i=1,eigvcount
         write(59,*)  adr(i)-shift, adi(i)                  
      enddo
      write (*,*) 'eigenvalue computation complete.'
      write (*,*) 'eigenvalue with max real part:'
      write (*,*) adr(1)-shift, adi(1)
#endif 
c// USE_ARPACK

      return
      end

cccccc  
c-----------------------------------------------------------------------
      subroutine cem_op
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)

      call cem_curl
      call cem_restrict_to_face
      call cem_flux (srflx)
      call cem_add_flux_to_res (srflx)
      call cem_source
      call cem_apps    
      call cem_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_apps 
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'

      if (ifdrude)  call cem_drude
      if (ifloren)  call cem_lorentz

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_op_test                 
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'

      integer n
      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)

      n= nx1*ny1*nz1*nelt
    
c... FIXME
c     if (ifsrc) then
c        call cem_inc
c        call cem_pmlsrf_to_face ! not used                 
c     endif
c... FIXME

      call cem_curl
      call cem_restrict_to_face
      call cem_flux (srflx)
      call cem_add_flux_to_res (srflx)

c... FIXME
c     if (ifpml) call pml_step
c     if (ifsrc) then 
c        if (ifsftf) then
c            call cem_source ! FIXME check out: nano pml, aps pml     
c        else
c            call cem_source ! FIXME check out: nano pml, aps pml     
c        endif  
c     endif
c
c     if (ifdrude) call cem_drude
c... FIXME
 
      call cem_invqmass

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_op_exp
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      integer icalled
      save    icalled
      data    icalled/0/
      common  /exptmp/ zh3(larnol)
      real    zh3
      complex*16 ttmp

      n6 = nx1*ny1*nz1*nelt*(3+3*(ndim-2))
      n  = nx1*ny1*nz1*nelt
      m  = larnol
      marnol  = larnol

      if (icalld.eq.0) then
          call rzero(ru ,n6)
          call rzero(rv0,n6*(m+1))
          call rzero(rv ,n6*m)
          call rzero(rvt,n6*m)
          call rzero(rh0,(m+1)*m)
          call rzero(rh ,m*m)
          call rzero(rw ,m)
          call rzero(wsp,2*m*(m+2))
          icalld = 1
      endif

      if (if3d) ioption=3 
      if (iftm) ioption=2 
      if (ifte) ioption=1 

      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,ioption)

      call arnoldi_vec (rv0,rh0,n6,m,ru,beta)
      call arnoldi_mat (rvt,rv,rh,rv0,rh0,m,n6)          

      call copy_r2z    (zh,rh,m*m)
      call zgeev  ('N','V',m,zh,m,eval,vri,m,vr,m,awork,2*m,rwork,info)
      call copy_z2z    (vri,vr,m*m)
      
      call zgetrf      (m,m,vri,m,ipiv,info)      !vri= LU factorization                
      call zgetri      (m,vri,m,ipiv,awork,m,info)!vri= X^(-1)

      call z_exp       (eval,dt,m)          !eval=exp(dt*eval)
      call z_rowscale  (vri,eval,m)         !vri =eval*X^(-1);vri=X^(-1)
      call c_mxm       (vr,m,vri,m,zh,m)
      call copy_z2r    (zh3,zh,m*m)

      call mxm         (rv,n6,zh3,m,ru,1)
      call cmult       (ru,beta,n6)      

      if (if3d) ioption=-3 
      if (iftm) ioption=-2 
      if (ifte) ioption=-1 

      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,ioption)

      return
      end
c-----------------------------------------------------------------------
      subroutine flds_to_vec(u,ux1,uy1,uz1,ux2,uy2,uz2,n6,n,ioption)
c     include 'TOTAL'
      integer  n6,n,ioption
      real     u(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)

      if      (ioption.eq. 3) then          !if3D
           call copy (u(1 + 0*n),ux1,n)
           call copy (u(1 + 1*n),uy1,n)
           call copy (u(1 + 2*n),uz1,n)
           call copy (u(1 + 3*n),ux2,n)
           call copy (u(1 + 4*n),uy2,n)
           call copy (u(1 + 5*n),uz2,n)
      elseif  (ioption.eq. 2)   then        !ifTM
           call copy (u(1 + 0*n),ux1,n)
           call copy (u(1 + 1*n),uy1,n)
           call copy (u(1 + 2*n),uz2,n)
      elseif  (ioption.eq. 1)   then        !ifTE
           call copy (u(1 + 0*n),ux2,n)
           call copy (u(1 + 1*n),uy2,n)
           call copy (u(1 + 2*n),uz1,n)
      elseif  (ioption.eq.-3)   then        !if3D
           call copy (ux1,u(1 + 0*n),n)
           call copy (uy1,u(1 + 1*n),n)
           call copy (uz1,u(1 + 2*n),n)
           call copy (ux2,u(1 + 3*n),n)
           call copy (uy2,u(1 + 4*n),n)
           call copy (uz2,u(1 + 5*n),n)
      elseif  (ioption.eq.-2)   then        !ifTM
           call copy (ux1,u(1 + 0*n),n)
           call copy (uy1,u(1 + 1*n),n)
           call copy (uz2,u(1 + 2*n),n)
      elseif  (ioption.eq.-1)   then        !ifTE
           call copy (ux2,u(1 + 0*n),n)
           call copy (uy2,u(1 + 1*n),n)
           call copy (uz1,u(1 + 2*n),n)
      else
         write(6,*) 'ERROR:: define ioption 3,2,1,-3,-2,-1'
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine arnoldi_mat (vt,v,h,v0,h0,m,n)
      integer  m,n, i,j,k
      real     h0(m+1,m),h(m,m)
      real     v0(n,m+1)
      real     v (n,m)
      real     vt(m,n)     

      do j=1,m 
      do i=1,n
         v (i,j)= v0(i,j)
      enddo
      enddo

      do i=1,m
      do j=1,m 
         h(i,j) = h0(i,j)
      enddo
      enddo

      do j=1,m
      do i=1,n
         vt(j,i) = v(i,j)
      enddo
      enddo
     
      return
      end
c-----------------------------------------------------------------------
      subroutine arnoldi_vec(v,h,n,m,b,beta)
      integer  e,ii,jj, i,j,k
      real     h(m+1,m)         
      real     v(n,m+1),b(n)
      real     beta, cnt, eps, s1, si


      eps = 1.e-13

      beta= glsc2(b,b,n)  
      beta= sqrt(beta)    
      si  = 1./beta     ! should add a check and a comment here

      call cmult2(v(1,1),b,si,n)

      do j=1,m 

         call cem_amult (v(1,j+1),v(1,j),n)                              

         do i=1,j
            h(i,j) = glsc2 (v(1,j+1),v(1,i),n)
            cnt    =       -h(i,j  )
            call add2s2    (v(1,j+1),v(1,i),cnt,n)
         enddo

         s1 = glsc2(v(1,j+1),v(1,j+1),n)
         s1 = sqrt(s1)                                
         h(j+1,j) = s1                                      

         if (abs(h(j+1,j)).le. eps) then
             write(6,*) 'h', h(j+1,j)
             call exitt
             return   !  should add a comment here
         else
             si = 1.0/s1 
             call cmult(v(1,j+1),si,n)
         endif 

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_amult(aw,av,n6)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer  n,n6                 
      common  /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real     aw(n6),av(n6)

      n = nx1*ny1*nz1*nelt 

      if        (IF3D) then
            call copy(hx,av(1 + 0*n),n)
            call copy(hy,av(1 + 1*n),n)
            call copy(hz,av(1 + 2*n),n)
            call copy(ex,av(1 + 3*n),n)
            call copy(ey,av(1 + 4*n),n)
            call copy(ez,av(1 + 5*n),n)
      else
         if     (IFTM) then
            call copy(hx,av(1 + 0*n),n)
            call copy(hy,av(1 + 1*n),n)
            call copy(ez,av(1 + 2*n),n)
         elseif (IFTE) then
            call copy(ex,av(1 + 0*n),n)
            call copy(ey,av(1 + 1*n),n)
            call copy(hz,av(1 + 2*n),n)
         else
            write(6,*) 'ERROR::IFTM or IFTE not defined'
            call exitt
         endif
      endif

      call cem_curl
      call cem_restrict_to_face
      call cem_flux (srflx)
      call cem_add_flux_to_res (srflx)
      call cem_invqmass

      if        (IF3D) then
            call copy(aw(1 + 0*n),reshx,n)
            call copy(aw(1 + 1*n),reshy,n)
            call copy(aw(1 + 2*n),reshz,n)
            call copy(aw(1 + 3*n),resex,n)
            call copy(aw(1 + 4*n),resey,n)
            call copy(aw(1 + 5*n),resez,n)
      else
         if     (IFTM) then
            call copy(aw(1 + 0*n),reshx,n)
            call copy(aw(1 + 1*n),reshy,n)
            call copy(aw(1 + 2*n),resez,n)
         elseif (IFTE) then
            call copy(aw(1 + 0*n),resex,n)
            call copy(aw(1 + 1*n),resey,n)
            call copy(aw(1 + 2*n),reshz,n)
         else
            write(6,*) 'ERROR::IFTM or IFTE not defined'
            call exitt
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine amult(w,v,n6)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      common  /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real     w(n6),v(n6)
      integer  n                 

      n = nx1*ny1*nz1*nelt 

      if        (IF3D) then
            call copy(hx,v(1 + 0*n),n)
            call copy(hy,v(1 + 1*n),n)
            call copy(hz,v(1 + 2*n),n)
            call copy(ex,v(1 + 3*n),n)
            call copy(ey,v(1 + 4*n),n)
            call copy(ez,v(1 + 5*n),n)
      else
         if     (IFTM) then
            call copy(hx,v(1 + 0*n),n)
            call copy(hy,v(1 + 1*n),n)
            call copy(ez,v(1 + 2*n),n)
         elseif (IFTE) then
            call copy(ex,v(1 + 0*n),n)
            call copy(ey,v(1 + 1*n),n)
            call copy(hz,v(1 + 2*n),n)
         else
            write(6,*) 'ERROR::IFTM or IFTE not defined'
            call exitt
         endif
      endif

      call cem_curl
      call cem_restrict_to_face
      call cem_flux (srflx)
      call cem_add_flux_to_res (srflx)
      call cem_invqmass

      if        (IF3D) then
            call copy(w(1 + 0*n),reshx,n)
            call copy(w(1 + 1*n),reshy,n)
            call copy(w(1 + 2*n),reshz,n)
            call copy(w(1 + 3*n),resex,n)
            call copy(w(1 + 4*n),resey,n)
            call copy(w(1 + 5*n),resez,n)
      else
         if     (IFTM) then
            call copy(w(1 + 0*n),reshx,n)
            call copy(w(1 + 1*n),reshy,n)
            call copy(w(1 + 2*n),resez,n)
         elseif (IFTE) then
            call copy(w(1 + 0*n),resex,n)
            call copy(w(1 + 1*n),resey,n)
            call copy(w(1 + 2*n),reshz,n)
         else
            write(6,*) 'ERROR::IFTM or IFTE not defined'
            call exitt
         endif
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine copy_r2z(z,r,n)
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         z(i)=r(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2r(r,z,n)
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         r(i)=real(z(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2z(z2,z1,n)
      complex*16 z1(n),z2(n)

      do i=1,n
         z2(i)=z1(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_col2(a,b,n)
      complex*16 a(n),b(n)

      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c_mxm(a,n1,b,n2,c,n3) ! complex mxm
      include 'SIZE'
      include 'PARALLEL'   
c
      real a(1),b(1),c(1)   
      real c_one(2),c_zero(2)

      save c_one   ,c_zero
      data c_one  / 1. , 0. /
      data c_zero / 0. , 0. /

      !c=n1xn3  b=n2xn3  a=n1xn2

       call zgemm( 'N','N',n1,n3,n2,c_one,a,n1,b,n2,c_zero,c,n1) !(ifblas)
      !call cgemm( 'N','N',n1,n3,n2,c_one,a,n1,b,n2,c_zero,c,n1)

      return
      end
c-----------------------------------------------------------------------
      real function anorm(a,n,m)
      real a(n,m)

      anorm =0.0
      do j=1,m
         tmp=0.0
      do i=1,n
         tmp=tmp+abs(a(i,j))
      enddo      
      if (tmp.ge.anorm) anorm=tmp
      enddo      

      return
      end
c-----------------------------------------------------------------------
      subroutine z_rowscale(a,s,n)
      integer    n
      complex*16 a(n,n),s(n)

      do i=1,n
      do j=1,n
         a(i,j) = s(i)*a(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_exp(f,dt,n)
      integer    n
      real*8     dt
      complex*16 f(n),a

      do i=1,n
         a    = dt*f(i)
         f(i) = exp(a)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine z_exp_arnoldi(eh,h,dt,m,ch,eval,evec,eveci,ipiv,w1,w2)

      include 'SIZE'
      include 'PARALLEL'

      real     eh(m,m),h(m,m)       
      real     ch(2,m,m),eval(2*m),evec(2*m*m),eveci(2*m*m)
      real     w1(4*m*m),w2(4*m*m)                  
      real     dt
      integer  ipiv(m) 
      integer  info,i,j
      
      if (wdsize.eq.4) then
          write(6,*) 'real4 not supported in z_exp_arnoldi'
          call exitt
      endif

      do j=1,m
      do i=1,m
         ch(1,i,j) = h(i,j)
         ch(2,i,j) = 0.
      enddo
      enddo

      call zgeev('N','V',m,ch,m,eval,evec,m,evec,m,w1,2*m,w2,info)

      if (info.ne.0) then
         write(6,*) 'exp_arnoldi fail:',info,m
      endif

      call copy       (eveci,evec,2*m*m)            
      call zgetri     (m,eveci,m,ipiv,w2,2*m,info) ! inverse of evec = X^(-1)
      call z_exp      (eval,dt,m )                 ! dt*lambda
      call z_rowscale (eveci,eval,m )              ! eveci=exp(dt*eval)*X^(-1)
      call c_mxm      (evec,m,eveci,m,ch,m)        ! ch=X*exp(dt*lambda)*X^(-1)  

      do j=1,m
      do i=1,m
         eh(i,j) = ch(1,i,j)
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_uvp           
C----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'EMWAVE'

      integer e

      do e=1,nelt
         call nekuvp(e)
      enddo

      call cem_mater1 
      call cem_mater2
      call cem_mater3

      return
      end
C----------------------------------------------------------------------
      subroutine cem_energy         
C----------------------------------------------------------------------
c...  this works OK with param(81)=3: needs to be updated
c...  by completely replacing old subroutine cem_out_fields
c...  09/11/2006 misun

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      
      return
      end

c ...mpi file routines starts here...
C----------------------------------------------------------------------
      subroutine cem_out_fields6(dumpno)
C----------------------------------------------------------------------
c...  this works OK with param(81)=6,-6,5,8:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-07-08 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      integer n,dumpno
      integer maxnelt,numgroups,numfields,maxnumnodes
      save    maxnelt,numgroups,numfields,maxnumnodes

      integer icalld        
      save    icalld          
      data    icalld /0/               

      common  ncell,ncellsize,ncelltype     
      integer ncell,ncellsize,ncelltype

      common /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real            vtkxyz
      common /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer         fpcnt
      common /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      if (icalld.eq.0) then
          maxnelt = iglmax(nelt,1)
          call bcast(maxnelt,isize)
          maxnumnodes = nx1*ny1*nz1*maxnelt ! max # nodes per processor

          if (nid.eq.0) then
              numfields   = param(80)           ! this number has to be exact.                            
              numgroups   = param(82)           ! param(82) = # groups 
          endif
          call bcast(numfields,isize)
          call bcast(numgroups,isize)

          call vtk_dump_header_nonswap()!save nodes and cells 
          icalld=1
      endif

      npts = nx1*ny1*nz1*nelt

c ... VTK FORMAT BINARY

      io_option= 6
      call print_options(io_option)

      call initrbio    (numgroups,numfields,maxnumnodes) !rbio_nekcem.c
      call openfile6   (dumpno,nid)                     
      call writeheader6()
      call writenodes6 (vtkxyz,npts)

      if  (if3d) then
          call write3dcells6(fpcnt,nelt,ncell,npts)
      else
          call write2dcells6(fpcnt,nelt,ncell,npts)
      endif

c...  number of fields are four: 
      if (ifschrod) then
      call writefield6 (11,vtkout1,npts) 
      call writefield6 (12,vtkout2,npts) 
      call writefield6 (13,vtkout3,npts)
      call writefield6 (14,vtkout4,npts)
      else
      call writefield6 (1,vtkout1,npts) 
      call writefield6 (2,vtkout2,npts) 
      call writefield6 (3,vtkout3,npts)
      call writefield6 (4,vtkout4,npts)
      endif

      inumfields=4

      if (inumfields.ne.numfields) then
          if (nid.eq.0) write(6,*) 'ERR: param(80) should be set as 4'
          call exitt
      endif 

      call closefile6  ()

c     call pvtk_xml(dumpno)
c     if(param(81).eq.8) then
c       if(nid.eq.0) call pvtk_nmm(dumpno)
c     elseif(param(81).eq.18) then
c       if(nid.eq.0) call pvtk_nmm(dumpno)
c     elseif(param(81).eq.5) then
c       if(nid.eq.0) call pvtk_nm(dumpno)
c     endif

      return
      end

C----------------------------------------------------------------------
      subroutine cem_restart_out(dumpno)    
C----------------------------------------------------------------------
c...  this works OK with param(81)=4:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-06-28 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      include 'SCHROD'

      integer icalld
      save    icalld            
      data    icalld /0/
      integer n,dumpno,idump

      common  ncell,ncellsize,ncelltype     
      integer ncell,ncellsize,ncelltype
      common  /s_rvtk/ vtkxyz_s(3*lx1*ly1*lz1*lelt)
      real             vtkxyz_s
      common  /s_ivtk/ fpcnt_s (8*lx1*ly1*lz1*lelt)
      integer          fpcnt_s
      common  /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               , vtkout2(3*lx1*ly1*lz1*lelt)
     $               , vtkout3(3*lx1*ly1*lz1*lelt)
     $               , vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      npts = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
         call vtk_dump_header_swap()     !save nodes and cells 
         icalld=1
      endif

      call openfile_restart (dumpno,nid) !vtkbin.c

c...  calling c-routines
      call writeheader4(istep,dumpno,nx1,np,time,dt)
      call writenodes4 (vtkxyz_s,npts)
      if (if3d) then
          call write3dcells4_swap(fpcnt_s,nelt,ncell,npts)
      else
          call write2dcells4_swap(fpcnt_s,nelt,ncell,npts)  
      endif

      if (ifschrod) then

        if (ifdouble) then
          call writefield4_double (11,vtkout1,npts)           
          call writefield4_double (12,vtkout2,npts) 
          call writefield4_double (13,vtkout3,npts)           
          call writefield4_double (14,vtkout4,npts) 
        else
          call writefield4 (11,vtkout1,npts)           
          call writefield4 (12,vtkout2,npts) 
          call writefield4 (13,vtkout3,npts)           
          call writefield4 (14,vtkout4,npts) 
        endif

      else

        if (ifdouble) then
          call writefield4_double (1,vtkout1,npts)           
          call writefield4_double (2,vtkout2,npts) 
          call writefield4_double (3,vtkout3,npts)           
          call writefield4_double (4,vtkout4,npts) 
        else
          call writefield4 (1,vtkout1,npts)           
          call writefield4 (2,vtkout2,npts) 
          call writefield4 (3,vtkout3,npts)           
          call writefield4 (4,vtkout4,npts) 
        endif

      endif

      call closefile_restart () 

c     if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
c       write(6,*) ' restart-file-write completed  ::::'
c     endif

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields4(dumpno)
C----------------------------------------------------------------------
c...  this works OK with param(81)=4:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-06-28 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'SCHROD'

      integer icalld
      save    icalld         
      data    icalld /0/                 
      integer n,dumpno,idump

      common  ncell,ncellsize,ncelltype     
      integer ncell,ncellsize,ncelltype

      common /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real            vtkxyz
      common /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer         fpcnt
      common /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      io_option= 4
      call print_options(io_option)

      npts = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call vtk_dump_header_nonswap() !save nodes and cells 
          icalld=1
      endif

c...  calling c-routines
      call openfile4   (dumpno,nid)   ! vtkbin.c
      call writeheader4(istep,dumpno,nx1,np,time,dt)
      call writenodes4 (vtkxyz,npts)
      if (if3d) then
          call write3dcells4(fpcnt,nelt,ncell,npts)
      else
          call write2dcells4(fpcnt,nelt,ncell,npts)  
      endif

      if (ifschrod) then

      call writefield4 (11,vtkout1,npts) 
      call writefield4 (12,vtkout2,npts) 
      call writefield4 (13,vtkout3,npts)
      call writefield4 (14,vtkout4,npts)

      else
      call writefield4 (1,vtkout1,npts) 
      call writefield4 (2,vtkout2,npts) 
      call writefield4 (3,vtkout3,npts)
      call writefield4 (4,vtkout4,npts)
      endif

      call closefile4  () 

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields3(dumpno)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'SCHROD'

      integer icalld
      save    icalld         
      data    icalld /0/
      integer n,dumpno,idump

      common  ncell,ncellsize,ncelltype
      integer ncell,ncellsize,ncelltype

      common /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real            vtkxyz
      common /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer         fpcnt
      common /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)

      io_option= 3
      call print_options(io_option)

      npts = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call vtk_dump_header_nonswap() !save nodes and cells
          icalld=1
      endif

c...  calling c-routines
      call openfile   (dumpno,nid)     ! vtkbin.c
      call writeheader()
      call writenodes (vtkxyz,npts)
      if (if3d) then
          call write3dcells(fpcnt,nelt,ncell,npts)
      else
          call write2dcells(fpcnt,nelt,ncell,npts)
      endif

      if (ifschrod) then
      call writefield (11,vtkout1,npts) 
      call writefield (12,vtkout2,npts) 
      call writefield (13,vtkout3,npts)
      call writefield (14,vtkout4,npts)
      else
      call writefield (1,vtkout1,npts) 
      call writefield (2,vtkout2,npts) 
      call writefield (3,vtkout3,npts)
      call writefield (4,vtkout4,npts)
      endif

      call closefile()  ! vtkbin.c

      if(nid.eq.0) call pvtk_xml3(dumpno)

      return 
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields2(dumpno)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'SCHROD'
      
      integer    fnum, dumpno
      parameter (fnum=111)
      logical    outex, outey, outez, outhx, outhy, outhz
      integer    n

      character*33 filename
      data         filename /'./vtk/ascii-NN-p000000-t00000.vtk'/
      CHARACTER*1  fnarray(33)
      equivalence  (filename,fnarray)

      integer   icalld
      save      icalld          
      data      icalld /0/             

      character*1  numrl(0:9)
      data numrl /'0','1','2','3','4','5','6','7','8','9'/

      fnarray(29)= numrl(mod(dumpno,10))
      fnarray(28)= numrl(mod(dumpno,100)/10)
      fnarray(27)= numrl(mod(dumpno,1000)/100)
      fnarray(26)= numrl(mod(dumpno,10000)/1000)
      fnarray(25)= numrl(mod(dumpno,100000)/10000)

      fnarray(22)= numrl(mod(nid,10))
      fnarray(21)= numrl(mod(nid,100)/10)
      fnarray(20)= numrl(mod(nid,1000)/100)
      fnarray(19)= numrl(mod(nid,10000)/1000)
      fnarray(18)= numrl(mod(nid,100000)/10000)
      fnarray(17)= numrl(mod(nid,1000000)/100000)

      io_option= 2
      call print_options(io_option)

      n = nx1*ny1*nz1*nelt

c ... VTK FORMAT: ASCII

      open(fnum, file=filename)      
      call vtk_dump_header2(fnum)

      outex = exactive
      outey = eyactive
      outez = ezactive
      outhx = hxactive
      outhy = hyactive
      outhz = hzactive

      if (ifscat) then

        call vtk_dump_field2(fnum,'totH  ',tothx,tothy,tothz,istep,time)
        call vtk_dump_field2(fnum,'totE  ',totex,totey,totez,istep,time)
        call vtk_dump_field2(fnum,'scatH ',hx,hy,hz,istep,time)
        call vtk_dump_field2(fnum,'scatE ',ex,ey,ez,istep,time)
        call vtk_dump_field2(fnum,'incH  ',inchx,inchy,inchz,istep,time)
        call vtk_dump_field2(fnum,'incE  ',incex,incey,incez,istep,time)

      else

        call vtk_dump_field2(fnum,'H     ',hx,hy,hz,istep,time)
        call vtk_dump_field2(fnum,'E     ',ex,ey,ez,istep,time)

      endif

      call vtk_dump_field2
     $              (fnum,'region',faceif,faceif,faceif,istep,time)!51 region
      call vtk_dump_field2
     $              (fnum,'face  ',faceid,faceid,faceid,istep,time)!52 face line
   
      if (ifsol) then

        if (outhx) call sub3(errhx, hx, shx, n)
        if (outhy) call sub3(errhy, hy, shy, n)
        if (outhz) call sub3(errhz, hz, shz, n)
        if (outex) call sub3(errex, ex, sex, n)
        if (outey) call sub3(errey, ey, sey, n)
        if (outez) call sub3(errez, ez, sez, n)

       call vtk_dump_field2(fnum,'errH  ',errhx,errhy,errhz,istep,time)
       call vtk_dump_field2(fnum,'errE  ',errex,errey,errez,istep,time)
       call vtk_dump_field2(fnum,'solH  ',shx,shy,shz,istep,time)
       call vtk_dump_field2(fnum,'solE  ',sex,sey,sez,istep,time)
          
      endif

      close(fnum)

      if(nid.eq.0) call pvtk_xml2(dumpno)   

      return 
      end

C----------------------------------------------------------------------
      subroutine cem_out_fields(dumpno)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'SCHROD'
      
      logical outex, outey, outez, outhx, outhy, outhz
      integer n,dumpno
      integer icalld
      save    icalld
      data    icalld /0/


      ifxyo=.true.
      if (istep.gt.iostep) ifxyo = .false.  ! dump geom for 1st file only
      ifvo  = .true.
      ifpo  = .false.
      ifto  = .false.

      io_option= 0
      call print_options(io_option)

      n = nx1*ny1*nz1*nelt

        if (NSTEPS.eq.0) then
C           call outpost(hx,hy,hz,pr,hx,'hxy')
C           call outpost(ex,ey,ez,pr,ey,'exy')
            return
        endif

        if (ifsol) then

c         call outpost(hx,hy,hz,pr,hx,'hxy')
c         call outpost(ex,ey,ez,pr,ey,'exy')

        else 

          if (ifscat) then

c             call outpost(tothx,tothy,tothz,pr,hx,'h-t')
c             call outpost(totex,totey,totez,pr,ey,'e-t')
c             call outpost(hx,hy,hz,pr,hx,'hxy')
c             call outpost(ex,ey,ez,pr,ey,'exy')
c             call outpost(engyh,engyh,engyh,pr,ey,'enh')
c             call outpost(engye,engye,engye,pr,ex,'ene')

          else

c            call outpost(hx,hy,hz,pr,hx,'hxy')
c            call outpost(ex,ey,ez,pr,ey,'exy')

          endif ! ifscat

        endif ! ifsol

      return 
      end

C----------------------------------------------------------------------
      subroutine print_options(io_option)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      
      integer  io_option

      if ((mod(istep,iostep).eq.0).or.(istep.eq.nsteps)) then
      if (nid.eq.0) then
        write(6,*)
        write(6,34) ifsrc,iftm,ifte,ifcentral,ifupwind,ifpml,ifpec,
     $  ifeig,ifexp,ifrk45,ifrk44,ifrk33,ifrk22,ifsim5,ifsim4,ifsim3,
     $  ifsim2,ifsim1                     

        write(6,35) ifdealias,ifscat,ifsftf,ifnm,ifdielec,iffilter,
     $  ifsol,ifpoisson,ifdrude,ifloren,ifwake,ifnapy,ifbeam,lxd,
     $  marnol

c...   not printint due to check if (io_option.eq.0) for now, it ok
        if (io_option.lt.0) return
        if (io_option.eq.0) write(6,3) istep,io_option
        if (io_option.eq.0) write(6,4) istep,io_option
        if (io_option.eq.0) write(6,5) istep,io_option
        if (io_option.eq.0) write(6,6) istep,io_option
        if (io_option.eq.0) write(6,7) istep,io_option
      endif

      endif

   33 format(a12,i10)
   34 format(' SRC=',L1,', TM=',L1,', TE=',L1,', CENTRAL=',L1,
     $ ', UPWIND=',L1,', PML=',L1,', PEC=',L1,' EIG=',L1,', EXP=',L1,
     $ ', RK45=',L1,', RK4=',L1,', RK3=',L1,', RK2=',L1,'SIM5=',L1,
     $ ', SIM4=',L1,', SIM3=',L1,', SIM2=',L1,', SIM1=',L1)
   35 format(' DEALIAS=',L1,', SCAT=',L1,', SFTF=',L1,', NM=',L1,
     $ ', DIELEC=',L1,', FILTER=',L1,', SOL=',L1,', POISSON=',L1,
     $ ', DRUDE=',L1,', LORENTZ=',L1,', WAKE=',L1,', NAPOLY=',L1,
     $ ', BEAM=',L1,', Dealias_M=',i3,', Krylov_Dim=',i3,'\n')     

   3  format(' ---> dumping postx format, istep=',i8,', iooption=',i3)
   4  format(' ---> dumping ascii  vtk, istep=',i8,', iooption=',i3)
   5  format(' ---> dumping binary vtk, istep=',i8,', iooption=',i3)
   6  format(' ---> dumping binary vtk, istep=',i8,', iooption=',i3)
   7  format(' ---> dumping binary vtk, istep=',i8,', iooption=',i3)

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      common  /outvtk/vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      common /ccpu  / cpu_t,cpu_dtime,cpu_chk
      common /commt / comm_t
      common /iotime/ cpu_t_io, cpu_io_step, swaptime
      real            cpu_t,cpu_dtime,cpu_chk,comm_t 
      real            cpu_t_io, cpu_io_step, swaptime

      common /iocnt/  dumpno
      integer         dumpno
      integer iooption,numfiles,irestart,idouble
      save    iooption,numfiles,irestart,idouble
      integer icomptrace,iotrace
      save    icomptrace,iotrace        
      integer iparam1,iparam2,iparam3,iparam4,iparam5,iparam6

      integer icalld
      save    icalld
      data    icalld/0/
      real*8  dclock

c...  dumpno for restarting case comes from previous step
      if (icalld.eq.0) then
         cpu_t_io= 0.0

         if (nid.eq.0) then
             iparam1=param(81)
             iparam2=param(82)
             iparam3=param(83)      !io-frequency for restart files     
             iparam4=0              !restart dumpno             
             iparam5=param(85)      !write compute trace                  
             iparam6=param(86)      !write io trace                       
             iparam7=param(87)      !write io trace                       
         else
             iparam1=0             
             iparam2=0            
             iparam3=0              !io-frequency for restart files     
             iparam4=0                 
             iparam5=0              !write compute trace                  
             iparam6=0              !write io trace                       
             iparam7=0              !write io trace                       
         endif
         if (ifrestart) then
            if (nid.eq.0) then
             iparam4= param(84)+1  !restart dumpno update
            else
             iparam4= 0                     
            endif
         endif
c...     iglsum is not doing smart bcast'ing, revert back to bcast -- jingfu
         call bcast(iparam1,isize) 
         call bcast(iparam2,isize) 
         call bcast(iparam3,isize) 
         call bcast(iparam4,isize) 
         call bcast(iparam5,isize) 
         call bcast(iparam6,isize) 
         call bcast(iparam7,isize) 

         iooption    =iparam1      
         numfiles    =iparam2        
         irestart    =iparam3                                       
         dumpno      =iparam4
         icomptrace  =iparam5
         iotrace     =iparam6
         idouble     =iparam7

         IFDOUBLE=.true.  
         if (idouble.ne.0) IFDOUBLE=.false. 
         if (nid.eq.0) write(6,*) 'IFDOUBLE=',IFDOUBLE
         if (nid.eq.0) write(6,*) 'param(81)=',iparam1
         if (nid.eq.0) write(6,*) 'param(82)=',iparam2
         if (nid.eq.0) write(6,*) 'param(83)=',irestart            
         if (nid.eq.0) write(6,*) 'dumpno   =',dumpno

         icalld  = 1
         call pass_io_params(iooption,numfiles)
      endif

c...  param(85) decide if compute trace should be generated
      if (icomptrace.gt.0) call writecomputetrace(istep,cpu_dtime,cpu_t)

c...  no output/restart files will be generated
      if (iooption.eq.0) return 

c...  this produces output for (io-option=iparam) with timing info
      if (iostep.gt.0) then

      if ((mod(istep,iostep).eq.0).or.(istep.eq.nsteps)) then
        
         if (ifschrod) then

         call vtk_nonswap_field( Ur, Ui, Ui,vtkout1)
         call vtk_nonswap_field( Qr(1,1),Qr(1,2),Qr(1,3),vtkout2)
         call vtk_nonswap_field( Qi(1,1),Qi(1,2),Qi(1,3),vtkout3)
         call vtk_nonswap_field( sUr,sQr(1,1),sQi(1,1),vtkout4)

         else

         call vtk_nonswap_field( HN(1,1), HN(1,2), HN(1,3),vtkout1)
         call vtk_nonswap_field( EN(1,1), EN(1,2), EN(1,3),vtkout2)
c        call vtk_nonswap_field(faceif,faceif,faceif,vtkout3)
         call vtk_nonswap_field(sHN(1,1),sHN(1,2),sHN(1,3),vtkout3)
         call vtk_nonswap_field(sEN(1,1),sEN(1,2),sEN(1,3),vtkout4)

         endif

         pretime = dclock()
         call starttiming
         call checkpoint_write(iooption,dumpno)
c...     maybe need a thread_join function to ensure the last io step finished
         call endtiming
         endtime = dclock()-pretime
         cpu_t_io= cpu_t_io+endtime

         call printio(iooption,istep)

c...  produce restart files    
         if (irestart.gt.0) then !io-frequency for restart files >= nonzero

          if (mod(istep,irestart*iostep).eq.0) then

#ifdef NO_MPI
           write(6,*) 'non-mpi run: no restart file produced'
#else
           swaptime=dclock() 

           if (ifschrod) then
               call vtk_swap_field( Ur, Ui, Ui,vtkout1)
               call vtk_swap_field( Qr(1,1),Qr(1,2),Qr(1,3),vtkout2)
               call vtk_swap_field( Qi(1,1),Qi(1,2),Qi(1,3),vtkout3)
               call vtk_swap_field( sUr,sQr(1,1),sQi(1,1),vtkout4)
           else
               call vtk_swap_field(HN(1,1),HN(1,2),HN(1,3),vtkout1)
               call vtk_swap_field(EN(1,1),EN(1,2),EN(1,3),vtkout2)
               call vtk_swap_field(sHN(1,1),sHN(1,2),sHN(1,3),vtkout3)
               call vtk_swap_field(sEN(1,1),sEN(1,2),sEN(1,3),vtkout4)
           endif

           swaptime=dclock()-swaptime 
           swaptime=swaptime/4.0
           i99=99

          !For threading, restart frequency should be zero on BG/P: fixme -- misun/jing 3/23/2012 
          !To be tested the same case on XK6: not likely to cause the same situation.                      
           if (nid.eq.0) write(6,3) dumpno,irestart*iostep,swaptime
           call checkpoint_write(i99,dumpno)
#endif

          endif     
         endif     
    3    format(' ---> write restart files at dumpno=',i10,
     $          ', irestart*iostep=',i10,', swaptime=',e11.4,'\n')
         dumpno = dumpno + 1  

      endif
      endif
c
      return
      end
c---------------------------------------------------------------
      subroutine cem_geombound
c---------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'GEOMBOUND'

      n = nx1*ny1*nz1*nelv
      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)
      zmin = glmin(zm1,n)
      zmax = glmax(zm1,n)

      xsize = xmax-xmin
      ysize = ymax-ymin
      zsize = zmax-zmin

      xcenter = 0.5*(xmin+xmax)
      ycenter = 0.5*(ymin+ymax)
      zcenter = 0.5*(zmin+zmax)

      return
      end
c---------------------------------------------------------------
      subroutine cem_init
c---------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'

      n    = nx1*ny1*nz1*nelv
      n6   = nx1*ny1*nz1*nelv*(3+3*(ndim-2))
      nsrf = nx1*nz1*nelv*2*ndim
      m    = marnol

      call rzero(shx,n)
      call rzero(shy,n)
      call rzero(shz,n)
      call rzero(sex,n)
      call rzero(sey,n)
      call rzero(sez,n)

      call rzero(hx,n)
      call rzero(hy,n)
      call rzero(hz,n)
      call rzero(ex,n)
      call rzero(ey,n)
      call rzero(ez,n)

      call rzero(jx,3*n)
      call rzero(jy,3*n)
      call rzero(jz,3*n)
      call rzero(px,3*n)
      call rzero(py,3*n)
      call rzero(pz,3*n)

      call rzero(ww1,NSRF)
      call rzero(ww2,NSRF)
      call rzero(ww3,NSRF)
      call rzero(uu1,NSRF)
      call rzero(uu2,NSRF)
      call rzero(uu3,NSRF)

c     call rzero(Fhx,NSRF)
c     call rzero(Fhy,NSRF)
c     call rzero(Fhz,NSRF)
c     call rzero(Fex,NSRF)
c     call rzero(Fey,NSRF)
c     call rzero(Fez,NSRF)
c     call rzero(Fjx,NSRF)
c     call rzero(Fjy,NSRF)
c     call rzero(Fjz,NSRF)
c     call rzero(fincex,NSRF)
c     call rzero(fincey,NSRF)
c     call rzero(fincez,NSRF)
c     call rzero(finchx,NSRF)
c     call rzero(finchy,NSRF)
c     call rzero(finchz,NSRF)
c     call rzero(srftmp,6*NSRF)
c     call rzero(srfinc,6*NSRF)
c     call rzero(srfin ,6*NSRF)
c     call rzero(srfout,6*NSRF)
c     call rzero(srflxi,6*NSRF)
c     call rzero(srflxo,6*NSRF)

      call rzero(khx,n)
      call rzero(khy,n)
      call rzero(khz,n)
      call rzero(kex,n)
      call rzero(key,n)
      call rzero(kez,n)

      call rzero(kjx,3*n)
      call rzero(kjy,3*n)
      call rzero(kjz,3*n)
      call rzero(kpx,3*n)
      call rzero(kpy,3*n)
      call rzero(kpz,3*n)

      call rzero(kbx,n)
      call rzero(kby,n)
      call rzero(kbz,n)
      call rzero(kdx,n)
      call rzero(kdy,n)
      call rzero(kdz,n)

      call rzero(reshx,n)
      call rzero(reshy,n)
      call rzero(reshz,n)
      call rzero(resex,n)
      call rzero(resey,n)
      call rzero(resez,n)

      call rzero(resjx,3*n)
      call rzero(resjy,3*n)
      call rzero(resjz,3*n)
      call rzero(respx,3*n)
      call rzero(respy,3*n)
      call rzero(respz,3*n)

      call rzero(lorena,3*n)
      call rzero(lorenb,3*n)
      call rzero(lorenc,3*n)

      call rzero(epsilo,n)
      call rzero(mu,    n)
      call rzero(impede,n)
      call rzero(conduc,n)
      call rzero(drudea,n)
      call rzero(drudeb,n)
      call rzero(faceif,n)

      call rzero(errex,n)
      call rzero(errey,n)
      call rzero(errez,n)
      call rzero(errhx,n)
      call rzero(errhy,n)
      call rzero(errhz,n)

      nxyzl=nx1*ny1*nz1*nelt*ndim 
      call rzero(pmlsigma,nxyzl)
      call rzero(pmldx,n)
      call rzero(pmldy,n)
      call rzero(pmldz,n)
      call rzero(pmlbx,n)
      call rzero(pmlby,n)
      call rzero(pmlbz,n)

      call rzero(respmldx,n)                      
      call rzero(respmldy,n)                      
      call rzero(respmldz,n)                      
      call rzero(respmlbx,n)                      
      call rzero(respmlby,n)                      
      call rzero(respmlbz,n)                      

      call rzero(respmlex,n)                      
      call rzero(respmley,n)                      
      call rzero(respmlez,n)                      
      call rzero(respmlhx,n)                      
      call rzero(respmlhy,n)                      
      call rzero(respmlhz,n)                      

      call rzero(kpmlbx,n)
      call rzero(kpmlby,n)
      call rzero(kpmlbz,n)
      call rzero(kpmldx,n)
      call rzero(kpmldy,n)
      call rzero(kpmldz,n)

      call izero(pmlptr,nelt)

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_mater2
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i, n
      real    tmpmu, tmpep, tmpt1, tmpt2

      n   = nx1*ny1*nz1*nelt

      do i=1,n 

         tmpmu = MU      (i,1,1,1)                      
         tmpep = EPSILO  (i,1,1,1)
         tmpt1 = sqrt(tmpmu/tmpep)
         tmpt2 = sqrt(tmpep/tmpmu)

         IMPEDE(i,1,1,1) = tmpt1                
         CONDUC(i,1,1,1) = tmpt2                                       

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_mater3 
c
c     Restrict u to faces
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
 
      integer nxyzf
      integer icalld
      save    icalld
      data    icalld  /0/

      real const
      character CB*3
      integer   e,f,i,j

      nxyzf = nx1*nz1*2*ndim*nelt
      nface = 2*ndim
      nxz   = nx1*nz1

      if (icalld.eq.0) then
         call cem_set_fc_ptr
         icalld = 1
      endif

      nxyzf = nx1*nz1*2*ndim*nelt

      do j=1,ncemface
         i              = cemface(j)
         ZIMPEDE(j,1,1) = IMPEDE(i,1,1,1)
         YCONDUC(j,1,1) = CONDUC(i,1,1,1)
         FACEF  (j,1,1) = FACEIF(i,1,1,1) 
      enddo

      call copy (Z_0,ZIMPEDE,nxyzf)
      call copy (Y_0,YCONDUC,nxyzf)

      call gs_op_fields(gsh_face,Z_0,nxyzf,1,1,1,0)
      call gs_op_fields(gsh_face,Y_0,nxyzf,1,1,1,0)

      do e= 1,nelt            
      do f= 1,nface
      do j= 1,nxz

         i   = (e-1)*nxz*nface + nxz*(f-1) + j  ! face numbering

         CB = CBC(f,e,2)

         if (CB.eq.'PEC') then
             Y_0(i,1,1) = 2.0*Y_0(i,1,1)
             Y_1(i,1,1) = 2.0*Y_1(i,1,1)
             Z_0(i,1,1) = 2.0*Z_0(i,1,1)
             Z_1(i,1,1) = 2.0*Z_1(i,1,1)
         endif

      enddo
      enddo
      enddo

      call sub3 (Z_1,Z_0,ZIMPEDE,NXYZF)
      call sub3 (Y_1,Y_0,YCONDUC,NXYZF)

      const = 0.5
 
      call cmult (Z_0,const,nxyzf)
      call cmult (Y_0,const,nxyzf)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_swap_across_faces(fld)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real fld    (LX1*LZ1,2*LDIM,LELT)
      real fldcpy (LX1*LZ1,2*LDIM,LELT)
      real origfld(LX1*LZ1,2*LDIM,LELT)

      call copy(origfld, fld, ncemface)
      call copy( fldcpy, fld, ncemface)
      call gs_op(gsh_face,fldcpy,1,1,0)

      ! the additional copy was necessary to avoid aliasing fld 
      ! fld = fldcpy - origfld

      call sub3 (fld, fldcpy, origfld, ncemface)

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_verify_comm

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real rndfld  (lx1*lz1,2*ldim,lelt)
      real rndcopy (lx1*lz1,2*ldim,lelt)
      real facexfld(lx1*lz1,2*ldim,lelt)
      real faceyfld(lx1*lz1,2*ldim,lelt)
      real facezfld(lx1*lz1,2*ldim,lelt)
      real fullxfld(lx1,ly1,lz1,lelt)
      real fullyfld(lx1,ly1,lz1,lelt)
      real fullzfld(lx1,ly1,lz1,lelt)
      real samefld (lx1,ly1,lz1,lelt)

      integer n,j,i

      n = nx1*ny1*nz1*nelv

      ! fill field with random garbage
      do j=1,ncemface
         rndfld (j,1,1) = zbqlu01(0.)
         rndcopy(j,1,1) = rndfld(j,1,1)
      enddo

      ! perform communication
      call gs_op(gsh_face,rndfld,1,1,0)
      !call cem_swap_across_faces(rndfld)

      ! compare to original, flag identical values
      call rzero(samefld, n)
      do j=1,ncemface
         i              = cemface(j)
         if (abs(rndfld(j,1,1)-rndcopy(j,1,1)).le.1e-12) then
           samefld(i,1,1,1) = 1
         endif
      enddo

      ! fill fields of coordinates of faces
      do j=1,ncemface
         i              = cemface(j)
         facexfld(j,1,1) = xm1(i,1,1,1)
         faceyfld(j,1,1) = ym1(i,1,1,1)
         facezfld(j,1,1) = zm1(i,1,1,1)
      enddo

      ! swap them
      call cem_swap_across_faces(facexfld)
      call cem_swap_across_faces(faceyfld)
      call cem_swap_across_faces(facezfld)

      call rzero(fullxfld, n)
      call rzero(fullyfld, n)
      call rzero(fullzfld, n)

      ! compute the difference, put in full-space field
      do j=1,ncemface
         i              = cemface(j)
         fullxfld(i,1,1,1) = facexfld(j,1,1) - xm1(i,1,1,1)
         fullyfld(i,1,1,1) = faceyfld(j,1,1) - ym1(i,1,1,1)
         fullzfld(i,1,1,1) = facezfld(j,1,1) - zm1(i,1,1,1)
      enddo

      ! and dump for visualization
      call vtk_dump_geometry
      call vtk_dump_field(fullxfld, 'CMVX  ', 0, 0.)
      call vtk_dump_field(fullyfld, 'CMVY  ', 0, 0.)
      call vtk_dump_field(fullzfld, 'CMVZ  ', 0, 0.)
      call vtk_dump_field(samefld, 'SMFX  ', 0, 0.)
      call vtk_dump_field(samefld, 'SMFY  ', 0, 0.)
      call vtk_dump_field(samefld, 'SMFZ  ', 0, 0.)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_mater1
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real     sum1, sum2, sum3, sum4, sum5, sum6
      real     old1, old2, old3, old4, old5, old6
      real     new1, new2, new3, new4, new5, new6
      integer  n,i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2
      integer icalld                                               
      save    icalld 
      data    icalld /0/

      n   = nx1*ny1*nz1

      do e=1,nelv

         do i=1,n
            FACEIF(i,1,1,e)= if_in_region(e)                     
         enddo
                
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE CEM_CURL2        
c---------------------------------------------------------------

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'

      COMMON /SCRMG/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)

      call curl_3d(resEX,resEY,resEZ,HX,HY,HZ,.false.,w1,w2) ! check done
      call curl_3d(resHX,resHY,resHZ,EX,EY,EZ,.false.,w1,w2) ! check done
 
      n = nx1*ny1*nz1*nelt

      call chsign(resHX,n)
      call chsign(resHY,n)
      call chsign(resHZ,n)

      RETURN
      END

c-----------------------------------------------------------------------
      subroutine curl_3d(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
c
      include 'SIZE'
      include 'TOTAL'
c
      logical ifavg
c
      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
c
      ntot  = nx1*ny1*nz1*nelt

      call dudxyz(work1,u3,rym1,sym1,tym1,jacm1,1,2) ! work1=dw/dy ; work2=dv/dz
      call dudxyz(work2,u2,rzm1,szm1,tzm1,jacm1,1,3)
      call sub3  (w1,work1,work2,ntot)

      call dudxyz(work1,u1,rzm1,szm1,tzm1,jacm1,1,3) ! work1=du/dz ; work2=dw/dx
      call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
      call sub3  (w2,work1,work2,ntot)

      call dudxyz(work1,u2,rxm1,sxm1,txm1,jacm1,1,1) ! work1=dv/dx ; work2=du/dy
      call dudxyz(work2,u1,rym1,sym1,tym1,jacm1,1,2)
      call sub3  (w3,work1,work2,ntot)
 
      return
      end
c-----------------------------------------------------------------------
c     subroutine cem_divergence  ! pointwise                     
c
c     include 'SIZE'
c     include 'TOTAL'
c     include 'EMWAVE'
c
c     n = nx1*ny1*nz1*nelt
c
c     call cem_div(divex,divey,divez,ex,ey,ez)
c     call cem_div(divhx,divhy,divhz,hx,hy,hz)
c
c     do i=1,n
c        divE(i,1,1,1)= divex(i,1,1,1)
c    $                 +divey(i,1,1,1)
c    $                 +divez(i,1,1,1)
c        divH(i,1,1,1)= divhx(i,1,1,1)
c    $                 +divhy(i,1,1,1)
c    $                 +divhz(i,1,1,1)
c        
c        errDivE(i,1,1,1)=abs(divE(i,1,1,1)+chgq(i,1,1,1))
c        errDivH(i,1,1,1)=abs(divH(i,1,1,1))
c     enddo
c
c     if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
c       divEmax =glmax(errDivE,n)
c       divHmax =glmax(errDivH,n)
c       write(6,*) '.'
c       write(6,*) ' divEmax',divEmax,', divHmax',divHmax 
c       write(6,*) '.'
c     endif

c     return
c     end
c-----------------------------------------------------------------------
      subroutine cem_dcurl ! weighted form of curl

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      real*8   dclock
      COMMON /SCRMG/ w1(lx1,ly1,lz1,lelv),w2(lx1,ly1,lz1,lelv)

      n = nx1*ny1*nz1*nelt

      if     (imode.eq.3) then !IF3D

          call wght_dcurl  (resEX,resEY,resEZ,HX,HY,HZ,.false.,w1,w2)
          call wght_dcurl  (resHX,resHY,resHZ,EX,EY,EZ,.false.,w1,w2)
          call chsign (resHX,n)
          call chsign (resHY,n)
          call chsign (resHZ,n)

      elseif (imode.eq.2) then !IFTM

          call wght_dcurl(resHX,resHY,resEZ,HX,HY,EZ,.false.,w1,w2)
          call chsign (resHX,n)
          call chsign (resHY,n)

      elseif (imode.eq.1) then !IFTM
          call wght_dcurl(resEX,resEY,resHZ,EX,EY,HZ,.false.,w1,w2)
          call chsign (resHZ,n)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_curl      !weighted form of curl with proper sign
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real cplus, cminus

      cplus = 1.0
      cminus=-1.0

      if     (imode.eq.3) then !IF3D
          call maxwell_curl (resEX,resEY,resEZ,HX,HY,HZ,cplus ,cplus )
          call maxwell_curl (resHX,resHY,resHZ,EX,EY,EZ,cminus,cminus)

      elseif (imode.eq.2) then !IFTM
          call maxwell_curl (resHX,resHY,resEZ,HX,HY,EZ,cplus ,cminus)

      elseif (imode.eq.1) then !IFTE
          call maxwell_curl (resEX,resEY,resHZ,EX,EY,HZ,cminus,cplus )
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drude 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      n = nx1*ny1*nz1*nelt    

c...  simpler with explicit writing than calling mat1 routines many times
      do i=1,n
         cdrude1 = DRUDEA(i,1,1,1)!*bm1(i,1,1,1)
         cdrude2 = DRUDEB(i,1,1,1)!*bm1(i,1,1,1)

         resJX(i,1,1,1,1)= cdrude1*Jx(i,1,1,1,1)+cdrude2*Ex(i,1,1,1)
         resJY(i,1,1,1,1)= cdrude1*Jy(i,1,1,1,1)+cdrude2*Ey(i,1,1,1)
         resJZ(i,1,1,1,1)= cdrude1*Jz(i,1,1,1,1)+cdrude2*Ez(i,1,1,1)

        !resJX(i,1,1,1,1)= resJx(i,1,1,1,1)+incEx(i,1,1,1)
        !resJY(i,1,1,1,1)= resJy(i,1,1,1,1)+incEy(i,1,1,1)
        !resJZ(i,1,1,1,1)= resJz(i,1,1,1,1)+incEz(i,1,1,1)

         resEX(i,1,1,1) = resEX(i,1,1,1)-Jx(i,1,1,1,1)*bm1(i,1,1,1)
         resEY(i,1,1,1) = resEY(i,1,1,1)-Jy(i,1,1,1,1)*bm1(i,1,1,1)
         resEZ(i,1,1,1) = resEZ(i,1,1,1)-Jz(i,1,1,1,1)*bm1(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_lorentz
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      n = nx1*ny1*nz1*nelt    

      do i=1,n

         clorenA1 = LORENA(i,1,1,1,1)!*bm1(i,1,1,1)
         clorenB1 = LORENB(i,1,1,1,1)!*bm1(i,1,1,1)
         clorenC1 = LORENC(i,1,1,1,1)!*bm1(i,1,1,1)

         clorenA1_Jx= clorenA1*Jx(i,1,1,1,1)                
         clorenA1_Jy= clorenA1*Jy(i,1,1,1,1)                
         clorenA1_Jz= clorenA1*Jz(i,1,1,1,1)                

         clorenB1_Px= clorenB1*Px(i,1,1,1,1)                
         clorenB1_Py= clorenB1*Py(i,1,1,1,1)                
         clorenB1_Pz= clorenB1*Pz(i,1,1,1,1)                

         clorenC1_Ex= clorenC1*Ex(i,1,1,1)                
         clorenC1_Ey= clorenC1*Ey(i,1,1,1)                
         clorenC1_Ez= clorenC1*Ez(i,1,1,1)                

         resJx(i,1,1,1,1)= clorenA1_Jx+clorenB1_Px+clorenC1_Ex 
         resJy(i,1,1,1,1)= c1orenA1_Jy+clorenB1_Py+clorenC1_Ey
         resJz(i,1,1,1,1)= clorenA1_Jz+clorenB1_Pz+clorenC1_Ez

         clorenA2 = LORENA(i,1,1,1,2)!*bm1(i,1,1,1)
         clorenB2 = LORENB(i,1,1,1,2)!*bm1(i,1,1,1)
         clorenC2 = LORENC(i,1,1,1,2)!*bm1(i,1,1,1)

         clorenA2_Jx= clorenA2*Jx(i,1,1,1,2)                
         clorenA2_Jy= clorenA2*Jy(i,1,1,1,2)                
         clorenA2_Jz= clorenA2*Jz(i,1,1,1,2)                

         clorenB2_Px= clorenB2*Px(i,1,1,1,2)                
         clorenB2_Py= clorenB2*Py(i,1,1,1,2)                
         clorenB2_Pz= clorenB2*Pz(i,1,1,1,2)                

         clorenC2_Ex= clorenC2*Ex(i,1,1,1)                
         clorenC2_Ey= clorenC2*Ey(i,1,1,1)                
         clorenC2_Ez= clorenC2*Ez(i,1,1,1)                

         resJx(i,1,1,1,2)= clorenA2_Jx+clorenB2_Px+clorenC2_Ex 
         resJy(i,1,1,1,2)= c1orenA2_Jy+clorenB2_Py+clorenC2_Ey
         resJz(i,1,1,1,2)= clorenA2_Jz+clorenB2_Pz+clorenC2_Ez

         clorenA3 = LORENA(i,1,1,1,3)!*bm1(i,1,1,1)
         clorenB3 = LORENB(i,1,1,1,3)!*bm1(i,1,1,1)
         clorenC3 = LORENC(i,1,1,1,3)!*bm1(i,1,1,1)

         clorenA3_Jx= clorenA3*Jx(i,1,1,1,3)                
         clorenA3_Jy= clorenA3*Jy(i,1,1,1,3)                
         clorenA3_Jz= clorenA3*Jz(i,1,1,1,3)                

         clorenB3_Px= clorenB3*Px(i,1,1,1,3)                
         clorenB3_Py= clorenB3*Py(i,1,1,1,3)                
         clorenB3_Pz= clorenB3*Pz(i,1,1,1,3)                

         clorenC3_Ex= clorenC3*Ex(i,1,1,1)                
         clorenC3_Ey= clorenC3*Ey(i,1,1,1)                
         clorenC3_Ez= clorenC3*Ez(i,1,1,1)                

         resJx(i,1,1,1,3)= clorenA3_Jx+clorenB3_Px+clorenC3_Ex 
         resJy(i,1,1,1,3)= c1orenA3_Jy+clorenB3_Py+clorenC3_Ey
         resJz(i,1,1,1,3)= clorenA3_Jz+clorenB3_Pz+clorenC3_Ez

         cJx=(Jx(i,1,1,1,1)+Jx(i,1,1,1,2)+Jx(i,1,1,1,3))*bm1(i,1,1,1)
         cJy=(Jy(i,1,1,1,1)+Jy(i,1,1,1,2)+Jy(i,1,1,1,3))*bm1(i,1,1,1)
         cJz=(Jz(i,1,1,1,1)+Jz(i,1,1,1,2)+Jz(i,1,1,1,3))*bm1(i,1,1,1)

         resEX(i,1,1,1)= resEX(i,1,1,1)-cJx                 
         resEY(i,1,1,1)= resEY(i,1,1,1)-cJy                 
         resEZ(i,1,1,1)= resEZ(i,1,1,1)-cJz            

         resPX(i,1,1,1,1)= JX(i,1,1,1,1)                     
         resPY(i,1,1,1,1)= JY(i,1,1,1,1)                     
         resPZ(i,1,1,1,1)= JZ(i,1,1,1,1)                
         resPX(i,1,1,1,2)= JX(i,1,1,1,2)                     
         resPY(i,1,1,1,2)= JY(i,1,1,1,2)                     
         resPZ(i,1,1,1,2)= JZ(i,1,1,1,2)                
         resPX(i,1,1,1,3)= JX(i,1,1,1,3)                     
         resPY(i,1,1,1,3)= JY(i,1,1,1,3)                     
         resPZ(i,1,1,1,3)= JZ(i,1,1,1,3)                

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_invqmass
C---------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'

      integer icalld
      save    icalld
      data    icalld  /0/
      real*8  dclock

      integer n 

      n = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call col3(qbm1,epsilo,bm1,n)
          call col3(qbm2,mu    ,bm1,n)
          call invcol1(qbm1,n)
          call invcol1(qbm2,n)
          call invers2(qbm3,bm1,n)
          icalld = 1
      endif

      if     (imode.eq.3) then !IF3D
       
          do i=1,n
             reshx(i,1,1,1)=reshx(i,1,1,1)*qbm2(i,1,1,1) ! call col2 (reshx,qbm2,n)
             reshy(i,1,1,1)=reshy(i,1,1,1)*qbm2(i,1,1,1)
             reshz(i,1,1,1)=reshz(i,1,1,1)*qbm2(i,1,1,1)
             resex(i,1,1,1)=resex(i,1,1,1)*qbm1(i,1,1,1)
             resey(i,1,1,1)=resey(i,1,1,1)*qbm1(i,1,1,1)
             resez(i,1,1,1)=resez(i,1,1,1)*qbm1(i,1,1,1)
          enddo

      elseif (imode.eq.2) then !IFTM

          do i=1,n
             reshx(i,1,1,1)=reshx(i,1,1,1)*qbm2(i,1,1,1) ! call col2 (reshx,qbm2,n)
             reshy(i,1,1,1)=reshy(i,1,1,1)*qbm2(i,1,1,1)
             resez(i,1,1,1)=resez(i,1,1,1)*qbm1(i,1,1,1)
          enddo

      elseif (imode.eq.1) then !IFTE

           do i=1,n
             reshz(i,1,1,1)=reshz(i,1,1,1)*qbm2(i,1,1,1)
             resex(i,1,1,1)=resex(i,1,1,1)*qbm1(i,1,1,1)
             resey(i,1,1,1)=resey(i,1,1,1)*qbm1(i,1,1,1)
           enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_l2err
C---------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n
      real     l2err

      N = NX1*NY1*NZ1*NELT

      l2err =0.0
#if 0
      write(6,*) 'test 11'
      l2err = l2err+ GLSC3(resHX,QBM2,HX,N)
      l2err = l2err+ GLSC3(resHY,QBM2,HY,N)
      l2err = l2err+ GLSC3(resHZ,QBM2,HZ,N)
      l2err = l2err+ GLSC3(resEX,QBM1,EX,N)
      write(6,*) 'test 22'
      l2err = l2err+ GLSC3(resEY,QBM1,EY,N)
      l2err = l2err+ GLSC3(resEZ,QBM1,EZ,N)
      write(6,*) 'test 33'

      l2err = l2err+ GLSC3(HX,BM2,HX,N)
      l2err = l2err+ GLSC3(HY,BM2,HY,N)
      l2err = l2err+ GLSC3(HZ,BM2,HZ,N)
      l2err = l2err+ GLSC3(EX,BM1,EX,N)
      write(6,*) 'test 44'
      l2err = l2err+ GLSC3(EY,BM1,EY,N)
      l2err = l2err+ GLSC3(EZ,BM1,EZ,N)

      l2err = sqrt(l2err)
      energy=l2err
#endif
      return
      end


c-----------------------------------------------------------------------
      subroutine rk_c(i)
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      integer  i

      RKtime=time+dt*rk4c(i)

c     if (IFRK44) ilast=5
c     if (IFRK33) ilast=4

c     if (i.ne.ilast) then
c       nextsubdt = rk4c(i+1)
c     else
c       nextsubdt = rk4c(1)
c     endif

      return
      end


c-----------------------------------------------------------------------
      subroutine rk_ab(ii)
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'

      integer istep

      n=nx1*ny1*nz1*nelt

c     kHX=rk4a(i)*kHX

      ca = rk4a(ii)
      cb = rk4b(ii)

      if       (imode.eq.3)  then !IF3D

         call rk4_upd(hx,khx,reshx,cb,ca,dt,n)
         call rk4_upd(hy,khy,reshy,cb,ca,dt,n)
         call rk4_upd(hz,khz,reshz,cb,ca,dt,n)
         call rk4_upd(ex,kex,resex,cb,ca,dt,n)
         call rk4_upd(ey,key,resey,cb,ca,dt,n)
         call rk4_upd(ez,kez,resez,cb,ca,dt,n)
 
         if (ifdrude) then
         call rk4_upd(jx(1,1,1,1,1),kjx(1,1,1,1,1),resjx(1,1,1,1,1),
     $                cb,ca,dt,n)
         call rk4_upd(jy(1,1,1,1,1),kjy(1,1,1,1,1),resjy(1,1,1,1,1),
     $                cb,ca,dt,n)
         call rk4_upd(jz(1,1,1,1,1),kjz(1,1,1,1,1),resjz(1,1,1,1,1),
     $                cb,ca,dt,n)
         endif

         if (ifloren) then
         call rk4_upd(jx,kjx,resjx,cb,ca,dt,3*n)
         call rk4_upd(jy,kjy,resjy,cb,ca,dt,3*n)
         call rk4_upd(jz,kjz,resjz,cb,ca,dt,3*n)
         call rk4_upd(px,kpx,respx,cb,ca,dt,3*n)
         call rk4_upd(py,kpy,respy,cb,ca,dt,3*n)
         call rk4_upd(pz,kpz,respz,cb,ca,dt,3*n)
         endif

         if (ifpml) then
         call rk4_upd(pmlbx,kpmlbx,respmlbx,cb,ca,dt,n)
         call rk4_upd(pmlby,kpmlby,respmlby,cb,ca,dt,n)
         call rk4_upd(pmlbz,kpmlbz,respmlbz,cb,ca,dt,n)
         call rk4_upd(pmldx,kpmldx,respmldx,cb,ca,dt,n)
         call rk4_upd(pmldy,kpmldy,respmldy,cb,ca,dt,n)
         call rk4_upd(pmldz,kpmldz,respmldz,cb,ca,dt,n)
         endif

      elseif   (imode.eq.2) then !IFTM

         call rk4_upd(hx,khx,reshx,cb,ca,dt,n)
         call rk4_upd(hy,khy,reshy,cb,ca,dt,n)
         call rk4_upd(ez,kez,resez,cb,ca,dt,n)

         if (ifdrude) then
            call rk4_upd(jx(1,1,1,1,1),kjx(1,1,1,1,1),resjx(1,1,1,1,1),
     $                   cb,ca,dt,n)
            call rk4_upd(jy(1,1,1,1,1),kjy(1,1,1,1,1),resjy(1,1,1,1,1),
     $                   cb,ca,dt,n)
         endif

         if (ifloren) then
            call rk4_upd(jx,kjx,resjx,cb,ca,dt,3*n)
            call rk4_upd(jy,kjy,resjy,cb,ca,dt,3*n)
            call rk4_upd(px,kpx,respx,cb,ca,dt,3*n)
            call rk4_upd(py,kpy,respy,cb,ca,dt,3*n)
         endif

         if (ifpml) then
            call rk4_upd(pmlbx,kpmlbx,respmlbx,cb,ca,dt,n)
            call rk4_upd(pmlby,kpmlby,respmlby,cb,ca,dt,n)
            call rk4_upd(pmldz,kpmldz,respmldz,cb,ca,dt,n)
         endif

      elseif   (imode.eq.1) then !IFTE

         call rk4_upd(hz,khz,reshz,cb,ca,dt,n)
         call rk4_upd(ey,key,resey,cb,ca,dt,n)
         call rk4_upd(ex,kex,resex,cb,ca,dt,n)

         if (ifdrude) then
            call rk4_upd(jx(1,1,1,1,1),kjx(1,1,1,1,1),resjx(1,1,1,1,1),
     $                   cb,ca,dt,n)
            call rk4_upd(jy(1,1,1,1,1),kjy(1,1,1,1,1),resjy(1,1,1,1,1),
     $                   cb,ca,dt,n)
         endif

         if (ifloren) then
            call rk4_upd(jx,kjx,resjx,cb,ca,dt,3*n)
            call rk4_upd(jy,kjy,resjy,cb,ca,dt,3*n)
            call rk4_upd(px,kpx,respx,cb,ca,dt,3*n)
            call rk4_upd(py,kpy,respy,cb,ca,dt,3*n)
         endif

         if (ifpml) then
            call rk4_upd(pmlbz,kpmlbz,respmlbz,cb,ca,dt,n)
            call rk4_upd(pmldx,kpmldx,respmldx,cb,ca,dt,n)
            call rk4_upd(pmldy,kpmldy,respmldy,cb,ca,dt,n)
         endif

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine rk4_upd(h,kh,resh,cb,ca,dt,n)
c-----------------------------------------------------------------------
      real h(1),kh(1),resh(1)

c     do i=1,n          
c        kh(i  ) = ca*kh(i  ) + dt*resh(i  )
c        h (i  ) =     h(i  ) + cb*kh  (i  )
c     enddo
c     return


c     below routine calculates 4 consecutive ones at one time
      nresid = iand(n,3)

      do i=1,n-nresid,4

         kh(i  ) = ca*kh(i  ) + dt*resh(i  )
         kh(i+1) = ca*kh(i+1) + dt*resh(i+1)
         kh(i+2) = ca*kh(i+2) + dt*resh(i+2)
         kh(i+3) = ca*kh(i+3) + dt*resh(i+3)

         h (i  ) =     h(i  ) + cb*kh  (i  )
         h (i+1) =     h(i+1) + cb*kh  (i+1)
         h (i+2) =     h(i+2) + cb*kh  (i+2)
         h (i+3) =     h(i+3) + cb*kh  (i+3)

      enddo

      if (nresid.eq.0) then

         return

      elseif (nresid.eq.1) then
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n  ) =     h(n  ) + cb*kh  (n  )
      elseif (nresid.eq.2) then
         kh(n-1) = ca*kh(n-1) + dt*resh(n-1)
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n-1) =     h(n-1) + cb*kh  (n-1)
         h (n  ) =     h(n  ) + cb*kh  (n  )
      else
         kh(n-2) = ca*kh(n-2) + dt*resh(n-2)
         kh(n-1) = ca*kh(n-1) + dt*resh(n-1)
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n-2) =     h(n-2) + cb*kh  (n-2)
         h (n-1) =     h(n-1) + cb*kh  (n-1)
         h (n  ) =     h(n  ) + cb*kh  (n  )
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_storage
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      if     (IFRK45) then

      rk4a(1) =   0.0
      rk4a(2) =  -567301805773.0/1357537059087.0
      rk4a(3) =  -2404267990393.0/2016746695238.0
      rk4a(4) =  -3550918686646.0/2091501179385.0
      rk4a(5) =  -1275806237668.0/842570457699.0

      rk4b(1) =   1432997174477.0/9575080441755.0
      rk4b(2) =   5161836677717.0/13612068292357.0
      rk4b(3) =   1720146321549.0/2090206949498.0
      rk4b(4) =   3134564353537.0/4481467310338.0
      rk4b(5) =   2277821191437.0/14882151754819.0

      rk4c(1) =   0.0
      rk4c(2) =   1432997174477.0/9575080441755.0
      rk4c(3) =   2526269341429.0/6820363962896.0
      rk4c(4) =   2006345519317.0/3224310063776.0
      rk4c(5) =   2802321613138.0/2924317926251.0
      rk4c(6) =   1.

      elseif (IFRK22) then

      rk4a(1) =   0.0 
      rk4a(2) =  -1.0                         

      rk4b(1) =   1.0 
      rk4b(2) =   1.0/2.0                       

c     rk4a(1) =   0.0 
c     rk4a(2) =  -11/15 !756391.0/934407.0
c     rk4a(3) =  -5/3   !36441873.0/15625003.0
c     rk4a(4) =  -1     !1953125.0/1085297.0
c     !rk4a(5) =   0.0
c
c     rk4b(1) =   1/3  !8.0/141.0 
c     rk4b(2) =   5/6  !6627.0/2000.0
c     rk4b(3) =   3/5  !609375.0/1085297.0
c     rk4b(4) =   1/4  !198961.0/526383.0
c     !rk4b(5) =   0.0
c
c     rk4c(1) =   0.0
c     rk4c(2) =   1.0/3.0
c     rk4c(3) =   5.0/9.0
c     rk4c(4) =   8.0/9.0
c     !rk4c(5) =   1.0
c
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine restart_swap
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD' 

      common  /RSTTMP/  xmr1(lx1,ly1,lz1,lelt)
     $                , ymr1(lx1,ly1,lz1,lelt)
     $                , zmr1(lx1,ly1,lz1,lelt)
     $                , exr (lx1,ly1,lz1,lelt)
     $                , eyr (lx1,ly1,lz1,lelt)
     $                , ezr (lx1,ly1,lz1,lelt)
     $                , hxr (lx1,ly1,lz1,lelt)
     $                , hyr (lx1,ly1,lz1,lelt)
     $                , hzr (lx1,ly1,lz1,lelt)
      common  /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real             vtkxyz
      common  /c_ivtk/ fpcnt (8*lx1*lx1*lx1*lelt)
      integer          fpcnt

      !write in float --> read in double: doesn't work
      !write in float --> read in float : recover up to 7-8 digits
      !write in double--> read in double: recover up to 16 digits

      integer n,dumpno,idump,idummy
      integer iooption,numfiles     
      real*8  dclock,t0

      integer e,eg

      integer icalld
      save    icalld
      data    icalld/0/

      npts = nx1*ny1*nz1*nelt !                            
      nxyz = nx1*ny1*nz1      !                               
      nxyz3= nx1*ny1*nz1*3    !

      t0= dclock()

      if (nid.eq.0) then
          iooption    = param(81)
          numfiles    = param(82)
          irstdumpno  = param(84) 
          idouble     = param(87) 
      else
          iooption    = 0           
          numfiles    = 0        
          irstdumpno  = 0         
          idouble     = 0         
      endif
 
      call bcast(iooption  ,isize)
      call bcast(numfiles  ,isize)
      call bcast(irstdumpno,isize)
      call bcast(idouble   ,isize)

      IFDOUBLE=.true.                    ! zero: double
      if (idouble.ne.0) IFDOUBLE=.false. ! nonzero: float
      if (nid.eq.0) write(6,*) 'IFDOUBLE=',IFDOUBLE
      if (nid.eq.0) write(6,*) 'restart format param(87)=',idouble  
      if (nid.eq.0) write(6,*) 'restart param(84)=',irstdumpno

      t0= dclock()
      call pass_io_params  (iooption,numfiles)
      call openfile_restart(irstdumpno, nid)
      call readheader4     (irststep,idump,nn1,nn2,restart_time,dtvtk)

      if (nid.eq.0) then
          time = time + restart_time
      else
          time = 0.0
      endif
      time = glsum(time,1)
      if (nid.eq.0) write(6,*) 'restart at time=time+rst_time',time

      !call printchecks  (idump)
      call readnodes4   (vtkxyz,npts) ! dummy reading but necessary     
      call swap_real_backward(vtkxyz,nxyz3) 
      call save2vectors (xmr1,ymr1,zmr1,vtkxyz,npts)
      
      if (if3d) then
           call read3dcells4(idummy,idummy,idummy,npts)  !no return
      else
           call read2dcells4(idummy,idummy,idummy,npts)  !no return
      endif

      if (ifschrod) then

          if (ifdouble) then
          call readfield4_double  (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (ur,ui,ui,vtkxyz,npts)
          call readfield4_double  (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qr(1,1),qr(1,2),qr(1,3),vtkxyz,npts)
          call readfield4_double  (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qi(1,1),qi(1,2),qi(1,3),vtkxyz,npts)
          else
          call readfield4         (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (ur,ui,ui,vtkxyz,npts)
          call readfield4         (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qr(1,1),qr(1,2),qr(1,3),vtkxyz,npts)
          call readfield4         (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qi(1,1),qi(1,2),qi(1,3),vtkxyz,npts)
          endif

      else

          if (ifdouble) then
          call readfield4_double   (idummy,vtkxyz,npts)   
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (hx,hy,hz,vtkxyz,npts)
          call readfield4_double   (idummy,vtkxyz,npts)   
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (ex,ey,ez,vtkxyz,npts)
          else
          call readfield4          (idummy,vtkxyz,npts)   
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (hx,hy,hz,vtkxyz,npts)
          call readfield4          (idummy,vtkxyz,npts)    
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (ex,ey,ez,vtkxyz,npts)
          endif

      endif
      call closefile_restart

      return
      end
c-----------------------------------------------------------------------
      subroutine save2vectors (ux1,uy1,uz1,vtkxyz,n)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer i,k,n
      real    ux1(lx1,ly1,lz1,lelt)
      real    uy1(lx1,ly1,lz1,lelt)
      real    uz1(lx1,ly1,lz1,lelt)
      real    vtkxyz(3*lx1*ly1*lz1*lelt)

      k = 1
      do i = 1,n    
         ux1(i,1,1,1)=vtkxyz(k  )
         uy1(i,1,1,1)=vtkxyz(k+1)
         uz1(i,1,1,1)=vtkxyz(k+2)
         k = k+3
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine  printchecks(idump)   
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      integer     idump

      irstdumpno =idump       
      istep      =irststep

      if (nid.eq.0) write(6,*) 'restart timestep=',irststep,
     $ ', restart dumpno=',irstdumpno,', restart time=', time,'istep=',
     $    istep, 'nid',nid

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_inc                            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      !call userinc                                           

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_source                         

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real     tt, xx, yy, zz, u
      integer  ifld, e, i, j, k, n

      n = nx1*ny1*nz1*nelt

      if (.not.ifsrc) return

      call usersrc (57,reshx,reshy,reshz,resex,resey,resez)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_pmlsrf_to_face
c
c     Restrict u to faces
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'
      integer  e,ef,f
c
      call full2face(fincex, incex)
      call full2face(fincey, incey)
      call full2face(fincez, incez)
      call full2face(finchx, inchx)
      call full2face(finchy, inchy)
      call full2face(finchz, inchz)

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_charge                        

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real     tt, xx, yy, zz, u
      integer  ifld, e, i, j, k, n

      n = nx1*ny1*nz1*nelt

      call gradm1(resjx,resjy,resjz,jx)
      call gradm1(resjx,resjy,resjz,jy)
      call gradm1(resjx,resjy,resjz,jz)
      call add4(divJ,resjx,resjy,resjz,n)

      return
      end

c-----------------------------------------------------------------------
      subroutine rectangular_to_cylindrical(vr,vt,vz1,vx,vy,vz,theta)   
      real       vr,vt,vz1,vx,vy,vz,theta   

      vr= vx*cos(theta)+vy*sin(theta)
      vt=-vx*sin(theta)+vy*cos(theta)
      vz1= vz

      return
      end

c-----------------------------------------------------------------------
      subroutine cylindrical_to_rectangular(vx,vy,vz1,vr,vt,vz,theta)   
      real       vr,vt,vz1,vx,vy,vz,theta   

      vx= vr*cos(theta)-vt*sin(theta)
      vy= vr*sin(theta)+vt*cos(theta)
      vz1= vz

      return
      end
c-----------------------------------------------------------------------
      subroutine energy_field (u,v1,v2,v3,n)                 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real u(n),v1(n),v2(n),v3(n)
      real tmp

      const = 1.0

      do i=1,n 
         tmp = v1(i)**2+v2(i)**2+v3(i)**2
         tmp = const*tmp
         u(i)= sqrt (tmp)
      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine outfldrn (x,txt,kk)
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*7 txt
C
      do ie=1,nelv,4
         do iz=1,nz1,1
            if (iz.eq.1) write(6,106) txt,iz,ie,kk
            if (iz.gt.1) write(6,107) 
            i1 = ie+1
            do j=ny1,1,-1
               write(6,105) ((x(i,j,iz,ie+i0),i=1,nx1),i0=0,3)
            enddo
         enddo
      enddo
C
  107 FORMAT(' ')
  105 FORMAT(4(4f6.0,5x))
c 105 FORMAT(4f6.0,20x,4f6.0)
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',2I4,/,
     $           5X,'       X            ')
C
      return
      end

c-----------------------------------------------------------------------
      subroutine outmat (u,m,n,name6,e)

      include 'SIZE'

      real        u(m,n)
      character*6 name6
      integer     e     
c
c     Print out copies of a global matrix
c
      write(6,1) nid,m,n,name6,e
   1  format(//,3i6,'  Matrix:',2x,a6,i9,/)

      n15 = min(n,15)
      do i=1,m
         if (n15.gt.5) write(6,15) nid,name6,(u(i,j),j=1,n15)
         if (n15.le.5) write(6, 5) nid,name6,(u(i,j),j=1,n15)
      enddo
   5  format(i3,1x,a6,20(6f8.4,/,10x))
  15  format(i3,1x,a6,20(15f6.1,/,10x))
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_fc_ptr
c
c     Set up pointer to restrict u to faces ! NOTE: compact
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
c
      integer i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2

      nxyz  = nx1*ny1*nz1
      nxz   = nx1*nz1
      nface = 2*ndim
      nxzf  = nx1*nz1*nface ! red'd mod to area, unx, etc.

      k = 0

      do e=1,nelv
      do f=1,nface

         ef     = eface(f)
         js1    = skpdat(1,f)
         jf1    = skpdat(2,f)
         jskip1 = skpdat(3,f)
         js2    = skpdat(4,f)
         jf2    = skpdat(5,f)
         jskip2 = skpdat(6,f)

         i = 0
         do j2=js2,jf2,jskip2
         do j1=js1,jf1,jskip1

            i = i+1
            k = i+nxz*(ef-1)+nxzf*(e-1)           ! face   numbering
            cemface(k) = j1+nx1*(j2-1)+nxyz*(e-1) ! global numbering

            !write(6,*) 'fc',e,k,f,ef,srfout(k)!glo_num3d(cemface(k))
            ! glo_num(cemface(k))  get global numbering
            ! cemface gets the Ed's numbering at faces following pff order

         enddo
         enddo

      enddo
      enddo
      ncemface = nxzf*nelv

      return
      end
c-----------------------------------------------------------------------
      subroutine full2face(faceary, fullary)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     fullary(lx1,ly1,lz1,lelt)
      integer  i,j

      do j=1,ncemface
         i= cemface(j)
         faceary(j,1,1) = fullary(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine face2full(fullary, faceary)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     fullary(lx1,ly1,lz1,lelt)
      integer  i,j

      call rzero(fullary, lx1*ly1*lz1*lelt)

      do j= 1,ncemface
         i= cemface(j)
         fullary(i,1,1,1) = fullary(i,1,1,1)+faceary(j,1,1)
      enddo

      return
      end
crrrrrrrrrr
c-----------------------------------------------------------------------
      subroutine cem_restrict_to_face
c     Restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) call cem_set_fc_ptr
      icalld = 1

      if     (imode.eq.3) then !IF3D
        do j=1,ncemface
           i= cemface(j)
           fhx(j,1,1) = hx(i,1,1,1)
           fhy(j,1,1) = hy(i,1,1,1)
           fhz(j,1,1) = hz(i,1,1,1)
           fex(j,1,1) = ex(i,1,1,1)
           fey(j,1,1) = ey(i,1,1,1)
           fez(j,1,1) = ez(i,1,1,1)
        enddo
      elseif (imode.eq.2) then !IFTM
        do j=1,ncemface
           i= cemface(j)
           fhx(j,1,1) = hx(i,1,1,1)
           fhy(j,1,1) = hy(i,1,1,1)
           fez(j,1,1) = ez(i,1,1,1)
        enddo
      elseif (imode.eq.1) then !IFTE
        do j=1,ncemface
           i= cemface(j)
           fex(j,1,1) = ex(i,1,1,1)
           fey(j,1,1) = ey(i,1,1,1)
           fhz(j,1,1) = hz(i,1,1,1)
        enddo
      endif

c     if((ifbeam).or.(ifwake)) then
c       do j=1,ncemface
c          i= cemface(j)
c          fjz(j,1,1) = jz(i,1,1,1)
c       enddo
c     endif

c     call exitt

      return
      end

c-----------------------------------------------------------------------
      subroutine add_flux(fullfld, facefld, faceskip)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real    fullfld(lx1,ly1,lz1,lelt)
      real    facefld(*)
      integer faceskip

      k = nx1*nz1*2*ndim*nelt
      do j=1,ncemface
         i = cemface(j)
         a = area(j,1,1,1)
         fullfld(i,1,1,1) = fullfld(i,1,1,1) + a*facefld(faceskip*k+j)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_add_flux_to_res(srflx)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real srflx(1)

      k= nx1*nz1*2*ndim*nelt

      if     (imode.eq.3) then !IF3D
        do j = 1,ncemface
           i = cemface(j)
           a = area(j,1,1,1)
           reshx(i,1,1,1) = reshx(i,1,1,1) + a*srflx(0*k+j)
           reshy(i,1,1,1) = reshy(i,1,1,1) + a*srflx(1*k+j)
           reshz(i,1,1,1) = reshz(i,1,1,1) + a*srflx(2*k+j)
           resex(i,1,1,1) = resex(i,1,1,1) + a*srflx(3*k+j)
           resey(i,1,1,1) = resey(i,1,1,1) + a*srflx(4*k+j)
           resez(i,1,1,1) = resez(i,1,1,1) + a*srflx(5*k+j)
        enddo
      elseif (imode.eq.2) then !IFTM
        do j=  1,ncemface
           i = cemface(j)
           a = area(j,1,1,1)
           reshx(i,1,1,1) = reshx(i,1,1,1) + a*srflx(0*k+j)
           reshy(i,1,1,1) = reshy(i,1,1,1) + a*srflx(1*k+j)
           resez(i,1,1,1) = resez(i,1,1,1) + a*srflx(2*k+j)
        enddo
      elseif (imode.eq.1) then !IFTE
        do j = 1,ncemface
           i = cemface(j)
           a = area(j,1,1,1)
           resex(i,1,1,1) = resex(i,1,1,1) + a*srflx(0*k+j)
           resey(i,1,1,1) = resey(i,1,1,1) + a*srflx(1*k+j)
           reshz(i,1,1,1) = reshz(i,1,1,1) + a*srflx(2*k+j)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_flux(srflx)

      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      integer  icalld
      save     icalld
      data     icalld/0/
      real     srflx(1)

c...FIXME
c     if (ifsftf) then
c
c         !if (icalld.eq.0) then
c         !if (nid.eq.0) write(6,*) 'cem_outer',istep,icalld
c         !   if (ifpml) call cem_outer_pmlsrf
c         !   icalld = 1
c         !endif
c
c         if (icalld.eq.0) then
c            if (ifpml) call cem_outer_pmlsrf
c            icalld = 1
c         endif
c
c         if (if3d) then
c            call cem_flux3d_sftf (srflx)                 
c         else
c            call cem_flux2d_sftf (srflx)                 
c         endif
c
c     else 
c...FIXME

         if (if3d) then
             call cem_flux3d (srflx)
         else
             call cem_flux2d (srflx)
         endif

c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_flux2d (srflx)  

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      real*8  dclock,t0
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer nxyzf,i,j,k,e,f,nxz

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0

      if     (imode.eq.2) then !IFTM

        k = nxyzf
        do i=1,nxyzf                
           j=cemface(i)
         srflx(0*k+i)=-uny(i,1,1,1)*fEZ(i,1,1)  ! NxfHx
         srflx(1*k+i)= unx(i,1,1,1)*fEZ(i,1,1)  ! NxfHy
         srflx(2*k+i)=-unx(i,1,1,1)*fHY(i,1,1)+uny(i,1,1,1)*fHX(i,1,1) ! NxfEz
        enddo

        t0=dclock()
        call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)
        call measure_comm(t0)

        if (IFPML.or.IFPEC) call cem_flux_pec(srflx)                   

        k = nxyzf
        do i=1,nxyzf

         Y0 = Y_0(i,1,1)             
         Y1 = Y_1(i,1,1)  
         Z0 = Z_0(i,1,1)   
         Z1 = Z_1(i,1,1)    

         fu1 = uny(i,1,1,1)*srflx(2*k+i)                                
         fu2 =-unx(i,1,1,1)*srflx(2*k+i)                     
         fw3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)     

         srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)    
         srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)   
         srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)   
        enddo

      elseif (imode.eq.1) then !IFTE

        k = nxyzf
        do i=1,nxyzf                

         srflx(0*k+i)= -uny(i,1,1,1)*fHZ(i,1,1)  ! NxfHx   
         srflx(1*k+i)=  unx(i,1,1,1)*fHZ(i,1,1)  ! NxfHy   
         srflx(2*k+i)= -unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1) ! NxfEz

        enddo


 
        t0 = dclock()
        call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)
        call measure_comm(t0)

        if (IFPML.or.IFPEC) call cem_flux_pec(srflx)               

        k = nxyzf
        do i=1,nxyzf                

         Y0 = Y_0(i,1,1)
         Y1 = Y_1(i,1,1)
         Z0 = Z_0(i,1,1)
         Z1 = Z_1(i,1,1)

         fw1 = uny(i,1,1,1)*srflx(2*k+i)                            
         fw2 =-unx(i,1,1,1)*srflx(2*k+i)
         fu3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

         srflx(0*k+i) = 0.5/Z0*( Z1*srflx(0*k+i) - C0*fw1)      
         srflx(1*k+i) = 0.5/Z0*( Z1*srflx(1*k+i) - C0*fw2)       
         srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3) 

        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine measure_comm(t0)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      common  /commt/ comm_t
      real            comm_t
      integer icalld
      save    icalld
      data    icalld/0/
      real    t0
      real*8  dclock

      if (icalld.eq.0) comm_t=0.0

      dtime = dclock()-t0
      comm_t= comm_t+dtime
      icalld= icalld+1

      ifreq=5             
      if (ifexp) ifreq=larnol

      if ((mod(istep,iocomm).eq.0).and.(mod(icalld,ifreq).eq.0)) then
          if (nid.eq.0) write(6,3) dtime,comm_t  
      endif
    3 format(' gs_op_fields::',1pe11.4,' sec, ',1pe11.4,' sec')

      return
      end
cfffffff
c-----------------------------------------------------------------------
      subroutine cem_flux3d(srflx)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

c     n x dEx - c* n x n x dHx  
c
      real*8  dclock,t0
      real    srflx(1)

      nface =2*ndim
      nxyzf =nx1*nz1*nface*nelt ! face count
      nxz   =nx1*nz1
      
      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex
      k = nxyzf
      do i= 1,nxyzf
         ! - n x E
         srflx(0*k+i)=-uny(i,1,1,1)*fEZ(i,1,1)+unz(i,1,1,1)*fEY(i,1,1)
         srflx(1*k+i)=-unz(i,1,1,1)*fEX(i,1,1)+unx(i,1,1,1)*fEZ(i,1,1)
         srflx(2*k+i)=-unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1)
         ! - n x H
         srflx(3*k+i)=-uny(i,1,1,1)*fHZ(i,1,1)+unz(i,1,1,1)*fHY(i,1,1)
         srflx(4*k+i)=-unz(i,1,1,1)*fHX(i,1,1)+unx(i,1,1,1)*fHZ(i,1,1)
         srflx(5*k+i)=-unx(i,1,1,1)*fHY(i,1,1)+uny(i,1,1,1)*fHX(i,1,1)
      enddo

      t0= dclock()
      call gs_op_fields(gsh_face,srflx,nxyzf,6,1,1,0)
      call measure_comm(t0)

      !srflx now contains n x [E] and n x [H],with [F] = F^+ - F^-;(observe un^+ = - un^-)
      if (IFPEC.or.IFPML)  call cem_flux_pec(srflx)                 

      k = nxyzf
      do i=1,nxyzf

        Y0   = Y_0(i,1,1)
        Y1   = Y_1(i,1,1)
        Z0   = Z_0(i,1,1)
        Z1   = Z_1(i,1,1)
        Y02  =-0.5/Y0*Y1          
        Z02  = 0.5/Z0*Z1         
        C02Y = 0.5/Y0*C0         
        C02Z = 0.5/Z0*C0         

        ! fu = n x n x [H]
        fu1 =uny(i,1,1,1)*srflx(5*k+i)-unz(i,1,1,1)*srflx(4*k+i)
        fu2 =unz(i,1,1,1)*srflx(3*k+i)-unx(i,1,1,1)*srflx(5*k+i)
        fu3 =unx(i,1,1,1)*srflx(4*k+i)-uny(i,1,1,1)*srflx(3*k+i)

        ! fw = n x n x [E]
        fw1 =uny(i,1,1,1)*srflx(2*k+i)-unz(i,1,1,1)*srflx(1*k+i)
        fw2 =unz(i,1,1,1)*srflx(0*k+i)-unx(i,1,1,1)*srflx(2*k+i)
        fw3 =unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

        srflx(0*k+i) = Y02*srflx(0*k+i) - C02Y*fu1              
        srflx(1*k+i) = Y02*srflx(1*k+i) - C02Y*fu2
        srflx(2*k+i) = Y02*srflx(2*k+i) - C02Y*fu3          
        srflx(3*k+i) = Z02*srflx(3*k+i) - C02Z*fw1          
        srflx(4*k+i) = Z02*srflx(4*k+i) - C02Z*fw2         
        srflx(5*k+i) = Z02*srflx(5*k+i) - C02Z*fw3         

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_flux_pec(srflx)                              

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
 
      real      srflx(1)
      integer   nxyzf,nface,nxz,k,i,e,f,ef          
      character CB*3

      integer   icalld
      save      icalld
      data      icalld /0/

      common   /BDRY1/ cempec(6*lx1*lz1*2*ldim*lelt)
      common   /BDRY2/ ncempec
      integer   cempec,ncempec

      nxyzf= nx1*nz1*2*ndim*nelt
      nface= 2*ndim
      nxz  = nx1*nz1
      k    = nxyzf
      l    = 0

      if (icalld.eq.0) then

         do e= 1,nelt            
         do f= 1,nface
         do j= 1,nxz
            CB= CBC(f,e,2)
            i =(e-1)*nface*nxz+(f-1)*nxz+j
            if ((CB.eq.'PEC').or.(CB.eq.'PML')) then
              l= l+1
              cempec(l)= i
              ncempec  = l
            endif
         enddo
         enddo
         enddo
         icalld =1

      endif

      if (ncempec.eq.0) return
 
      if (if3d) then

         do j= 1,ncempec         
            i= cempec(j)

         ! here, srflx is (n x -E, n x -H)
         ! (we're on the boundary and fgs_gop left us untouched)
         ! since we want mirror conditions, i.e.
         ! n x E+ = -n x E-
         ! n x H+ =  n x H-,
         ! (cf. Hesthaven, Warburton, HO-/SM on Unst. Grids, 
         ! I. Time-Domain Solution of Maxwell's Eqns., p. 14)
         ! we have
         ! n x [E] = n x (E+ - E-) = n x -2E-
         ! n x [H] = n x (H+ - H-) = 0

         ! write(6,*) 'ii--',ncempec,l,cempec(l) 
          srflx(0*k+i)=2.0*srflx(0*k+i)
          srflx(1*k+i)=2.0*srflx(1*k+i)
          srflx(2*k+i)=2.0*srflx(2*k+i)
          srflx(3*k+i)=0
          srflx(4*k+i)=0
          srflx(5*k+i)=0
         enddo

      else
       if     (iftm) then

         do j= 1,ncempec         
            i= cempec(j)
            srflx(0*k+i)=2.0*srflx(0*k+i)
            srflx(1*k+i)=2.0*srflx(1*k+i)
            srflx(2*k+i)=0
         enddo

       elseif (ifte) then

         do j= 1,ncempec         
            i= cempec(j)
            srflx(0*k+i)=0
            srflx(1*k+i)=0
            srflx(2*k+i)=2.0*srflx(2*k+i)
         enddo

       endif
      endif
    
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_dflux_2d(srflx)
c     ffff
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'WZ'
      include 'DEALIAS'
      include 'PARALLEL'
c
      real srflx(1)

      parameter (ldd =lxd*4*lelt )
      parameter (md  =lxd        )

      common /dfaces/ fmhx(ldd),fmhy(ldd),fmhz(ldd)
     $              , fmex(ldd),fmey(ldd),fmez(ldd) 
      common /dscrn1/ srflxm (6*2*ldim*lxd*lelt)
      common /dscrn2/ srflxm0(6*2*ldim*lxd*lelt)
     $              , srflxm1(6*2*ldim*lxd*lelt)
     $              , srflxm2(6*2*ldim*lxd*lelt)
      common /weight/ zzd(lxd),wgld(lxd)
      common /invrho/ rho2i(lx1*lz1)
      real srflxm,srflxm0, srflxm1, srflxm2
      real srflx0,srflx1, srflx2                   

      integer e,f,k,i,jd,j1,nn,mm,nxyd,nxy,mx            
      integer icalld
      save    icalld
      data    icalld / 0 /

      mx    = nx1
      nxy   = nx1*nz1
      nxyd  = nxd*nzd
      nface = 2*ndim

      if (icalld.eq.0) then

         icalld = 1
         call zwgll(zzd,wgld,md)

         l = 0

         do k=1,nz1
         do i=1,nx1
            l = l+1
            rho2i(l) = 1./(wxm1(i)*wxm1(k))
         enddo
         enddo

         do e=1,nelt   ! this will need to be dynamic if we have p-refinement
         do f=1,nface

            call col3(ta,area(1,1,f,e),rho2i,nxy) ! surface Jacobian only
            call col3(wk,ta,unx(1,1,f,e),nxy)
            call intp_rsd (aud(1,1,f,e),wk,mx,md,if3d,0)  ! 1 = unx

            call col3(wk,ta,uny(1,1,f,e),nxy)
            call intp_rsd (aud(1,2,f,e),wk,mx,md,if3d,0)  ! 2 = uny

            call col3(wk,ta,unz(1,1,f,e),nxy)
            call intp_rsd (aud(1,3,f,e),wk,mx,md,if3d,0)  ! 3 = unz

            l = 0

            do k=1,nxd
            do i=1,nxd

               l   = l+1
               tmp = wgld(i)*wgld(k)
               aud(l,1,f,e) = tmp*aud(l,1,f,e)
               aud(l,2,f,e) = tmp*aud(l,2,f,e)
               aud(l,3,f,e) = tmp*aud(l,3,f,e)

            enddo
            enddo
               
         enddo
         enddo

      endif

      do e=1,nelv   
      do f=1,nface  
         call intp_rsd (fmez,fez(1,f,e),mx,md,if3d,0)                  
         call intp_rsd (fmhx,fhx(1,f,e),mx,md,if3d,0)                  
         call intp_rsd (fmhy,fhy(1,f,e),mx,md,if3d,0)                   
      enddo
      enddo

      k = 1

      myz= md*2*ndim*nelt

      do i=1,myz             

         srflxm(k  ) = aud(i,2,1,1)*fmez(i)                      ! fhx
         srflxm(k+1) =-aud(i,1,1,1)*fmez(i)
         srflxm(k+2) = aud(i,2,1,1)*fmhx(i)-aud(i,1,1,1)*fmhy(i) ! fex

         srflxm(k  ) = 0.5*srflxm(k  )                           
         srflxm(k+1) = 0.5*srflxm(k+1)                           
         srflxm(k+2) = 0.5*srflxm(k+2)                           

         k = k+3

      enddo

c     call fgs_gop_vec (gsh_face, srflxm ,'+\0',3) ! vector gs

      do i=1,myz             
         srflxm0(i)=srflxm(3*(i-1)+1)
         srflxm1(i)=srflxm(3*(i-1)+2)
         srflxm2(i)=srflxm(3*(i-1)+3)
      enddo

      do e=1,nelv   
      do f=1,nface  
         call intp_rsd (srflx0,srflxm0,mx,md,if3d,1)                    
         call intp_rsd (srflx1,srflxm1,mx,md,if3d,1)                   
         call intp_rsd (srflx2,srflxm2,mx,md,if3d,1)                   
      enddo
      enddo

      nyze = mx*2*ndim*nelt

      do i=1,nyze             
         srflx(3*(i-1)+1)=srflxm0(i)
         srflx(3*(i-1)+2)=srflxm1(i)
         srflx(3*(i-1)+3)=srflxm2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wght_dcurl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2) 

       include 'SIZE' 
       include 'TOTAL' 
       include 'EMWAVE' 
       include 'DEALIAS' 

c     ddddd
c
c     Combined weighted curl:   w12 = B*curl(u3), w3=B*curl(u1,u2)
c     with dealiasing
c
                        
      real     u1(1),u2(1),u3(1)
      real     w1(1),w2(1),w3(1)
      real     work1(1),work2(1)
      logical  ifavg
                                    
      parameter (ldd =lxd**ldim  )
      parameter (md  =lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)
 
      real           ud1r,ud1s,ud1t
     $           ,   ud2r,ud2s,ud2t
     $           ,   ud3r,ud3s,ud3t

      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)                     
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)

c     common /ddtmp5/ w3md(lxd**ldim),wgld(3*(lxd+lx1))

      real            dxmd,dxtmd,w3md,wgld,zmd

      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy,nn,mm
      save    icalld 
      data    icalld /0/
 
      mx    = nx1
      nn    = nx1-1   
      mm    = md -1

      nxyd  = md**ndim
      nxy   = mx**ndim

      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

          if (if3d) then
             l = 0
             do k=1,md
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)*wgld(k)
             enddo
             enddo
             enddo
          else
             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)              
             enddo
             enddo
          endif
                                           
          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)

      endif

      do e=1,nelt

         j  = nxyd *(e-1) 
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u1(j1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,u2(j1),mx,md,if3d, 0)
         call intp_rstd (ju3,u3(j1),mx,md,if3d, 0)

         if (if3d) then

             call local_grad3 (ud1r,ud1s,ud1t,ju1,mm,1,dxmd,dxtmd)
             call local_grad3 (ud2r,ud2s,ud2t,ju2,mm,1,dxmd,dxtmd)
             call local_grad3 (ud3r,ud3s,ud3t,ju3,mm,1,dxmd,dxtmd)

             do i=1,nxyd

                k = i+j              

                wd1(i) = ( ud3r(i)*rymd(k)
     $                   + ud3s(i)*symd(k)
     $                   + ud3t(i)*tymd(k)
     $                   - ud2r(i)*rzmd(k)
     $                   - ud2s(i)*szmd(k)
     $                   - ud2t(i)*tzmd(k) ) * w3md(i)

                wd2(i) = ( ud1r(i)*rzmd(k)
     $                   + ud1s(i)*szmd(k)
     $                   + ud1t(i)*tzmd(k) 
     $                   - ud3r(i)*rxmd(k)
     $                   - ud3s(i)*sxmd(k) 
     $                   - ud3t(i)*txmd(k) ) * w3md(i)

                wd3(i) = ( ud2r(i)*rxmd(k)
     $                   + ud2s(i)*sxmd(k) 
     $                   + ud2t(i)*txmd(k) 
     $                   - ud1r(i)*rymd(k) 
     $                   - ud1s(i)*symd(k) 
     $                   - ud1t(i)*tymd(k) ) * w3md(i)
             enddo 

         else

          call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
          call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
          call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

          do i= 1,nxyd

                k = i+j              

                wd1(i) = ( ud3r(i)*rymd(k)
     $                   + ud3s(i)*symd(k) ) * w3md(i)
 
                wd2(i) =-( ud3r(i)*rxmd(k)
     $                   + ud3s(i)*sxmd(k) ) * w3md(i)
 
                wd3(i) = ( ud2r(i)*rxmd(k)
     $                   + ud2s(i)*sxmd(k) 
     $                   - ud1r(i)*rymd(k) 
     $                   - ud1s(i)*symd(k) ) * w3md(i)

          enddo 

         endif

         call intp_rstd (w1(j1),wd1,mx,md,if3d,1) ! 1 --> backward
         call intp_rstd (w2(j1),wd2,mx,md,if3d,1) ! 1 --> backward
         call intp_rstd (w3(j1),wd3,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_grad (w1,w2,w3,u1,c1)
c
c     Weighted curl:   w = -grad(u1)
c
      include 'SIZE'
      include 'TOTAL'
c
      real w1(1),w2(1),w3(1),u1(1)                 
      integer e, k, i, j, nxyz,ntot,nn
      real c1

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)
         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) =c1*( u1r(i)*rxm1(k,1,1,1)
     $                 +u1s(i)*sxm1(k,1,1,1)
     $                 +u1t(i)*txm1(k,1,1,1) )/jacm1(i,1,1,e)
            w2(k) =c1*( u1r(i)*rym1(k,1,1,1)
     $                 +u1s(i)*sym1(k,1,1,1)
     $                 +u1t(i)*tym1(k,1,1,1) )/jacm1(i,1,1,e)
            w3(k) =c1*( u1r(i)*rzm1(k,1,1,1)
     $                 +u1s(i)*szm1(k,1,1,1)
     $                 +u1t(i)*tzm1(k,1,1,1) )/jacm1(i,1,1,e)
         enddo
         enddo
       else

         do e=1,nelt

         j = nxyz*(e-1)
         call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) =c1*(  u1r(i)*rxm1(k,1,1,1)
     $                 + u1s(i)*sxm1(k,1,1,1) )/jacm1(i,1,1,e)
            w2(k) =c1*(  u1r(i)*rym1(k,1,1,1)
     $                 + u1s(i)*sym1(k,1,1,1) )/jacm1(i,1,1,e)
c           write(6,*) 'tt',e,i,rxm1(k,1,1,1)/jacm1(i,1,1,e)
c    $                         ,rym1(k,1,1,1)/jacm1(i,1,1,e)
c           write(6,*) 'tt',e,i,u1r(i),u1s(i)                          
c           write(6,*) 'tt',e,i,xm1(k,1,1,1),ym1(k,1,1,1)             

         enddo
         enddo

       endif

       return
       end
c-----------------------------------------------------------------------
      subroutine curl_op(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)

c
c     Weighted curl:   w = B*curl(u)
c
      include 'SIZE'
      include 'TOTAL'
c
      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
      logical ifavg
      integer e, k, i, j, nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
      common /curlw/ sw1(lx1*ly1*lz1*lelt)
     $             , sw2(lx1*ly1*lz1*lelt)
     $             , sw3(lx1*ly1*lz1*lelt)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            w1(k) = (    u3r(i)*rym1(k,1,1,1)
     $                 + u3s(i)*sym1(k,1,1,1)
     $                 + u3t(i)*tym1(k,1,1,1)
     $                 - u2r(i)*rzm1(k,1,1,1)
     $                 - u2s(i)*szm1(k,1,1,1)
     $                 - u2t(i)*tzm1(k,1,1,1) ) /jacm1(i,1,1,e)
c
            w2(k) = (    u1r(i)*rzm1(k,1,1,1)
     $                 + u1s(i)*szm1(k,1,1,1)
     $                 + u1t(i)*tzm1(k,1,1,1) 
     $                 - u3r(i)*rxm1(k,1,1,1)
     $                 - u3s(i)*sxm1(k,1,1,1) 
     $                 - u3t(i)*txm1(k,1,1,1) ) /jacm1(i,1,1,e)
c
            w3(k) = (    u2r(i)*rxm1(k,1,1,1)
     $                 + u2s(i)*sxm1(k,1,1,1) 
     $                 + u2t(i)*txm1(k,1,1,1) 
     $                 - u1r(i)*rym1(k,1,1,1) 
     $                 - u1s(i)*sym1(k,1,1,1) 
     $                 - u1t(i)*tym1(k,1,1,1) ) /jacm1(i,1,1,e)

          enddo 
          enddo
    
       else

        do e= 1,nelt 

           j= nxyz*(e-1) 

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1) 
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

          k = i+j

          w1(k) = ( u3r(i)*rym1(k,1,1,1) 
     $            + u3s(i)*sym1(k,1,1,1)) /jacm1(i,1,1,e)

          w2(k) =-( u3r(i)*rxm1(k,1,1,1) 
     $            + u3s(i)*sxm1(k,1,1,1)) /jacm1(i,1,1,e)

          w3(k) = ( u2r(i)*rxm1(k,1,1,1) 
     $            + u2s(i)*sxm1(k,1,1,1) 
     $            - u1r(i)*rym1(k,1,1,1) 
     $            - u1s(i)*sym1(k,1,1,1)) /jacm1(i,1,1,e)

          enddo 

        enddo 
       endif

       return
       end
c-----------------------------------------------------------------------
c mmmmmm
      subroutine maxwell_curl(w1,w2,w3,u1,u2,u3,sgn1,sgn2)
c     Weighted curl:   dH/dt= -curl E  
c     Weighted curl:   dE/dt=  curl H  
      include 'SIZE'
      include 'TOTAL'
c
      integer mode, imode 
      real    w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      real    sgn1,sgn2  
      integer e,k,i,j,nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1

      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

c...        sgn1 is assigned "-" or "+" when calling; sgn2 is dummy; just for 3d.
c...        rym1,sym1,tym1,...jac are multiplied. w3m1 is purely weight.
c...        area(i,1,1,1) is multiplied with jacobina weight.
            u1rw = sgn1*u1r(i)*w3m1(i,1,1)
            u1sw = sgn1*u1s(i)*w3m1(i,1,1)
            u1tw = sgn1*u1t(i)*w3m1(i,1,1)
            u2rw = sgn1*u2r(i)*w3m1(i,1,1)
            u2sw = sgn1*u2s(i)*w3m1(i,1,1)
            u2tw = sgn1*u2t(i)*w3m1(i,1,1)
            u3rw = sgn1*u3r(i)*w3m1(i,1,1)
            u3sw = sgn1*u3s(i)*w3m1(i,1,1)
            u3tw = sgn1*u3t(i)*w3m1(i,1,1)

            w1(k) =      u3rw*rym1(k,1,1,1)
     $                 + u3sw*sym1(k,1,1,1)  
     $                 + u3tw*tym1(k,1,1,1)  
     $                 - u2rw*rzm1(k,1,1,1)
     $                 - u2sw*szm1(k,1,1,1)
     $                 - u2tw*tzm1(k,1,1,1)                   
 
            w2(k) =      u1rw*rzm1(k,1,1,1)
     $                 + u1sw*szm1(k,1,1,1)
     $                 + u1tw*tzm1(k,1,1,1) 
     $                 - u3rw*rxm1(k,1,1,1)
     $                 - u3sw*sxm1(k,1,1,1) 
     $                 - u3tw*txm1(k,1,1,1)                    
 
            w3(k) =      u2rw*rxm1(k,1,1,1)
     $                 + u2sw*sxm1(k,1,1,1) 
     $                 + u2tw*txm1(k,1,1,1) 
     $                 - u1rw*rym1(k,1,1,1) 
     $                 - u1sw*sym1(k,1,1,1) 
     $                 - u1tw*tym1(k,1,1,1)                 

          enddo 
          enddo

       else

        do e= 1,nelt 

           j= nxyz*(e-1) 

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1) 
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j

c...         sgn1/sgn2 is assigned "-" or "+" when calling; sgn1/sgn2=1/-1 if TM, 
c...         sgn1/sgn2 is assigned "-" or "+" when calling; sgn1/sgn2=-1/1 if TE, 

             u1rw = sgn1*u1r(i)*w3m1(i,1,1)
             u1sw = sgn1*u1s(i)*w3m1(i,1,1)
             u2rw = sgn1*u2r(i)*w3m1(i,1,1)
             u2sw = sgn1*u2s(i)*w3m1(i,1,1)

             u3rw = sgn2*u3r(i)*w3m1(i,1,1)
             u3sw = sgn2*u3s(i)*w3m1(i,1,1)

             w1(k) = ( u3rw*rym1(k,1,1,1) 
     $               + u3sw*sym1(k,1,1,1))

             w2(k) =-( u3rw*rxm1(k,1,1,1) 
     $               + u3sw*sxm1(k,1,1,1))

             w3(k) = ( u2rw*rxm1(k,1,1,1) 
     $               + u2sw*sxm1(k,1,1,1) 
     $               - u1rw*rym1(k,1,1,1) 
     $               - u1sw*sym1(k,1,1,1))

          enddo 

        enddo 
       endif

       return
       end
c-----------------------------------------------------------------------
c uuuuuu
      subroutine wght_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
c     Weighted curl:   w = B*curl(u)
      include 'SIZE'
      include 'TOTAL'
c
      logical ifavg
      real    w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
      integer e,k,i,j,nxyz,ntot,nn
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            u1rw = u1r(i)*w3m1(i,1,1)
            u1sw = u1s(i)*w3m1(i,1,1)
            u1tw = u1t(i)*w3m1(i,1,1)
            u2rw = u2r(i)*w3m1(i,1,1)
            u2sw = u2s(i)*w3m1(i,1,1)
            u2tw = u2t(i)*w3m1(i,1,1)
            u3rw = u3r(i)*w3m1(i,1,1)
            u3sw = u3s(i)*w3m1(i,1,1)
            u3tw = u3t(i)*w3m1(i,1,1)

            w1(k) =      u3rw*rym1(k,1,1,1)
     $                 + u3sw*sym1(k,1,1,1)  
     $                 + u3tw*tym1(k,1,1,1)  
     $                 - u2rw*rzm1(k,1,1,1)
     $                 - u2sw*szm1(k,1,1,1)
     $                 - u2tw*tzm1(k,1,1,1)                   
 
            w2(k) =      u1rw*rzm1(k,1,1,1)
     $                 + u1sw*szm1(k,1,1,1)
     $                 + u1tw*tzm1(k,1,1,1) 
     $                 - u3rw*rxm1(k,1,1,1)
     $                 - u3sw*sxm1(k,1,1,1) 
     $                 - u3tw*txm1(k,1,1,1)                    
 
            w3(k) =      u2rw*rxm1(k,1,1,1)
     $                 + u2sw*sxm1(k,1,1,1) 
     $                 + u2tw*txm1(k,1,1,1) 
     $                 - u1rw*rym1(k,1,1,1) 
     $                 - u1sw*sym1(k,1,1,1) 
     $                 - u1tw*tym1(k,1,1,1)                 

c           w1(k) = (    u3r(i)*rym1(k,1,1,1)
c    $                 + u3s(i)*sym1(k,1,1,1)
c    $                 + u3t(i)*tym1(k,1,1,1)
c    $                 - u2r(i)*rzm1(k,1,1,1)
c    $                 - u2s(i)*szm1(k,1,1,1)
c    $                 - u2t(i)*tzm1(k,1,1,1) ) * w3m1(i,1,1)
 
c           w2(k) = (    u1r(i)*rzm1(k,1,1,1)
c    $                 + u1s(i)*szm1(k,1,1,1)
c    $                 + u1t(i)*tzm1(k,1,1,1) 
c    $                 - u3r(i)*rxm1(k,1,1,1)
c    $                 - u3s(i)*sxm1(k,1,1,1) 
c    $                 - u3t(i)*txm1(k,1,1,1) ) * w3m1(i,1,1)
 
c           w3(k) = (    u2r(i)*rxm1(k,1,1,1)
c    $                 + u2s(i)*sxm1(k,1,1,1) 
c    $                 + u2t(i)*txm1(k,1,1,1) 
c    $                 - u1r(i)*rym1(k,1,1,1) 
c    $                 - u1s(i)*sym1(k,1,1,1) 
c    $                 - u1t(i)*tym1(k,1,1,1) ) * w3m1(i,1,1)

          enddo 
          enddo

       else

        do e= 1,nelt 

           j= nxyz*(e-1) 

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1) 
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

          k = i+j

            u1rw = u1r(i)*w3m1(i,1,1)
            u1sw = u1s(i)*w3m1(i,1,1)
            u2rw = u2r(i)*w3m1(i,1,1)
            u2sw = u2s(i)*w3m1(i,1,1)
            u3rw = u3r(i)*w3m1(i,1,1)
            u3sw = u3s(i)*w3m1(i,1,1)

          w1(k) = ( u3rw*rym1(k,1,1,1) 
     $            + u3sw*sym1(k,1,1,1))

          w2(k) =-( u3rw*rxm1(k,1,1,1) 
     $            + u3sw*sxm1(k,1,1,1))

          w3(k) = ( u2rw*rxm1(k,1,1,1) 
     $            + u2sw*sxm1(k,1,1,1) 
     $            - u1rw*rym1(k,1,1,1) 
     $            - u1sw*sym1(k,1,1,1))

          enddo 

        enddo 
       endif

       return
       end

c-----------------------------------------------------------------------
      subroutine wght_curl_2d3(w1,w2,w3,u1,u2,u3,ifavg,work1,work2) 
c-----------------------------------------------------------------------

       include 'SIZE' 
       include 'TOTAL' 
       include 'EMWAVE' 
c
c     Combined weighted curl:   w12 = B*curl(u3), w3=B*curl(u1,u2)
c
                        
       real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1) 
       logical ifavg
                                    
       integer e
       common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             ,  u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             ,  u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
       nxyz  = nx1*ny1*nz1
       ntot  = nx1*ny1*nz1*nelt
       nn    = nx1-1

       do e=1,nelt

          j = nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1)
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j 

             w1(k) = ( u3r(i)*rym1(k,1,1,1)
     $               + u3s(i)*sym1(k,1,1,1)) * w3m1(i,1,1)

             w2(k) =-( u3r(i)*rxm1(k,1,1,1) 
     $               + u3s(i)*sxm1(k,1,1,1)) * w3m1(i,1,1)

             w3(k) = ( u2r(i)*rxm1(k,1,1,1)
     $               + u2s(i)*sxm1(k,1,1,1)
     $               - u1r(i)*rym1(k,1,1,1)
     $               - u1s(i)*sym1(k,1,1,1) )* w3m1(i,1,1)

          enddo

       enddo

       return
       end

c-----------------------------------------------------------------------
      subroutine cem_div (w1,w2,w3,u1,u2,u3)
c     div:   w = div(u)/jacm1
      include 'SIZE'
      include 'TOTAL'
c
      real     w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer  e                           

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
 
      nxyz  = nx1*ny1*nz1
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(ur1,us1,ut1,u1,nn,e,dxm1,dxtm1)
         call local_grad3(ur2,us2,ut2,u2,nn,e,dxm1,dxtm1)
         call local_grad3(ur3,us3,ut3,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (  ur1(i)*rxm1(k,1,1,1)
     $               + us1(i)*sxm1(k,1,1,1)
     $               + ut1(i)*txm1(k,1,1,1))/jacm1(i,1,1,e)
            w2(k) = (  ur2(i)*rym1(k,1,1,1)
     $               + us2(i)*sym1(k,1,1,1)
     $               + ut2(i)*tym1(k,1,1,1))/jacm1(i,1,1,e)
            w3(k) = (  ur3(i)*rzm1(k,1,1,1)
     $               + us3(i)*szm1(k,1,1,1)
     $               + ut3(i)*tzm1(k,1,1,1))/jacm1(i,1,1,e)
         enddo

         enddo

      else

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
         call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (  ur1(i)*rxm1(k,1,1,1)
     $               + us1(i)*sxm1(k,1,1,1))/jacm1(i,1,1,e)
            w2(k) = (  ur2(i)*rym1(k,1,1,1)
     $               + us2(i)*sym1(k,1,1,1))/jacm1(i,1,1,e) 
            !write(6,*) '3d',e,i,u1(k)           
            !write(6,*) '3d',e,i,w1(k),w2(k)
            !write(6,*) '3d',e,i,rxm1(k,1,1,1)jacm1(i,1,1,e)
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)
C
C     -  Spectral interpolation from A to B via tensor products
C     -  scratch arrays: w(na*na*nb + nb*nb*na)
C
C     5/3/00  -- this routine replaces specmp in navier1.f, which
c                has a potential memory problem
C
C
      logical if3d
c
      real b(nb,nb,nb),a(na,na,na)
      real w(ldw)
c
      ltest = na*nb
      if (if3d) ltest = na*na*nb + nb*na*na
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine specmpnface(b,nb,a,na,ba,ab,if3d,w,ldw)
            
c     -  Spectral interpolation from A to B via tensor products
c     -  scratch arrays: w(na*na*nb + nb*nb*na)
c
c     5/3/00  -- this routine replaces specmp in navier1.f, which        
c                has a potential memory problem                           
          
      logical if3d                                                     
c          
      real b(nb,nb),a(na,na)
      real w(ldw)
c          
      ltest = na*nb
           
      if (if3d)      ltest = na*na*nb + nb*na*na
          
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c         
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine intp_rstd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      integer mx, md, idir
      logical if3d
c
      parameter       (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real            jgl,jgt,wkd
c
      parameter       (ld=2*lxd)
      common /ctmp0/  w(ld**ldim,2)

      integer   i

      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
 
      ldw = 2*(ld**ldim)
  
      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
         call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else
         call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif

      return
      end
             
c-----------------------------------------------------------------------
      subroutine intp_rsd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      integer mx, md, idir
      logical if3d
c
      parameter       (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/   d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                 jgt(ldg), wkd(lwkd)
      real             jgl,jgt,wkd
 
      parameter       (ld=2*lxd)
      common /ctmp0/   w(ld**ldim,2)

      integer   i

      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
 
      ldw = 2*(ld**ldim)
  
      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
          call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else 
          call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine gen_int(jgl,jgt,mp,np,w)
 
c     Generate interpolation from np GLL points to mp GL points
c
c        jgl  = interpolation matrix, mapping from velocity nodes to pressure
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c

      real jgl(mp,np),jgt(np*mp),w(1)
c
      iz = 1
      id = iz + np
c
      call zwgll (w(iz),jgt,np)
      call zwgll (w(id),jgt,mp)
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo
c
      call transpose(jgt,np,jgl,mp)
c
      return
      end
          
c-----------------------------------------------------------------------
      subroutine gen_dgl(dgl,dgt,mp,np,w)
c
c     Generate derivative from np GL points onto mp GL points
c
c        dgl  = interpolation matrix, mapping from velocity nodes to pressure
c        dgt  = transpose of interpolation matrix
c        w    = work array of size (3*np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
c
      real dgl(mp,np),dgt(np*mp),w(1)
c
c
      iz = 1
      id = iz + np  
c
      call zwgll  (w(iz),dgt,np)  ! GL points
      call zwgll  (w(id),dgt,mp)  ! GL points
c
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c
      n  = np -1 
      do i=1,mp   
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np  
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
         enddo
      enddo
c
      call transpose(dgt,np,dgl,mp)

c
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
      include 'SIZE'            ! need nid
      character*5  avar5,lvar5
      character*10 sub_name10
c
      if (n.gt.m) then
         write(6,*)             
         write(6,*) nid,n,m
c        write(6,*) nid,n,m,avar5,lvar5,sub_name10
c        write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i9,2(1x,a5),1x,a10)
         call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md)
c
c     Get pointer to jgl() for interpolation pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip, mx, md
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c
      if (ip.eq.0) then
c
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c
         nwrkd = mx + md

         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
c
         call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgl_ptr (ip,mx,md)
c
c     Get pointer to GL-GL interpolation dgl() for pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip, mx, md
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pdg (ij)
c
      if (ip.eq.0) then
c
         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c
         call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)

      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ifoption(ifif)        

      include 'SIZE'
      include 'EMWAVE'
      include 'PARALLEL'

      logical ifif
 
      if (np.eq.1) return

      ifip =0
      if (ifif) ifip =1

      imaxg = iglmax(ifip,1)
      if (imaxg.eq.1) ifif =.true.

      return
      end
c-----------------------------------------------------------------------
      subroutine put_face_value(face,bound0,nf)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real    bound0
      integer i,j,n,nf

      integer icalld
      save    icalld
      data    icalld /0/

      real face(nf)

      n    =nx1*ny1*nz1*nelt
      nface=2*ndim*nx1*nz1*nelt

      if (icalld.eq.0) then
          call cem_set_fc_ptr
          icalld = 1
      endif

      do j = 1,ncemface
         i = cemface(j)
         xx= xm1(i,1,1,1)
         yy= ym1(i,1,1,1)
         zz= zm1(i,1,1,1)

         rr=sqrt(xx**2+yy**2)

         if     (rr.le.bound0) then
            face(j)=-1   
         else                         
            face(j)= 1   
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zero_out(uu,n)
      integer n
      real uu(n)
    
      do i =1,n
         if (abs(uu(i)).le.1e-20) then
             uu(i)=0.0
         else
             uu(i)=uu(i)
         endif   
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine field_normalization(u,v,w,n)
      real u(n)
      real v(n)
      real w(n)

      unorm=glsc2(u,n)
      vnorm=glsc2(v,n)
      wnorm=glsc2(w,n)

      call cmult(u,unorm,n)
      call cmult(v,vnorm,n)
      call cmult(w,wnorm,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine fcsum2(xsum,asum,x,e,f)
c
c     Compute the weighted sum of X over face f of element e
c
c     x is an (NX,NY,NZ) data structure
c     f  is in the preprocessor notation 
c
c     xsum is sum (X*area)
c     asum is sum (area)

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'
      real x(lx1,ly1,lz1,1)
      integer e,f,fd

      asum = 0.
      xsum = 0.

c     Set up counters ;  fd is the dssum notation.
      call dsset(nx1,ny1,nz1)
      fd     = eface1(f)
      js1    = skpdat(1,fd)
      jf1    = skpdat(2,fd)
      jskip1 = skpdat(3,fd)
      js2    = skpdat(4,fd)
      jf2    = skpdat(5,fd)
      jskip2 = skpdat(6,fd)

      i = 0
      do j2=js2,jf2,jskip2
      do j1=js1,jf1,jskip1
         i = i+1
         xsum = xsum+area(i,1,f,e)*x(j1,j2,1,e)
         asum = asum+area(i,1,f,e)
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      real function geom_xyradius()
         
      ! return the maximum radius in the x-y plane 
      ! of the nodes in the mesh.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer n,i
      real    radius_squared,radhere_squared
      real    gop_work
      
      n   = 8*nelt
      radius_squared = 0

      do i=1,n
         radhere_squared = xc(i,1)*xc(i,1) + yc(i,1)*yc(i,1)
         if (radhere_squared.ge.radius_squared) then
             radius_squared = radhere_squared
         endif
      enddo

      call gop(radius_squared, gop_work, 'M  ', 1)
      geom_xyradius = sqrt(radius_squared)

      return
      end
c ---------------------------------------------------------
      real function final_geom_xyradius()

      ! return the maximum radius in the x-y plane 
      ! of the nodes in the mesh.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer n,i
      real radius_squared,radhere_squared
      real gop_work
      
      n = nx1*ny1*nz1*nelt
      radius_squared = 0

      do i=1,n
         radhere_squared = xm1(i,1,1,1)*xm1(i,1,1,1) 
     $                   + ym1(i,1,1,1)*ym1(i,1,1,1)
         if (radhere_squared.ge.radius_squared) then
             radius_squared = radhere_squared
         endif
      enddo

      call gop(radius_squared, gop_work, 'M  ', 1)
      final_geom_xyradius = sqrt(radius_squared)

      return
      end
c-----------------------------------------------------------------------
