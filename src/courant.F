c-----------------------------------------------------------------------
      subroutine get_dxmin ( dxmin )
c
c     Compute dxmin over domain
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      integer e


      d2m = 1.e20

      if (if3d) then
         do e=1,nelv
         do k=2,nz1-1
         do j=2,ny1-1
         do i=2,nx1-1

            km=k-1
            kp=k+1
            jm=j-1
            jp=j+1
            im=i-1
            ip=i+1

            dx = xm1(ip,j,k,e) - xm1(im,j,k,e)
            dy = ym1(ip,j,k,e) - ym1(im,j,k,e)
            dz = zm1(ip,j,k,e) - zm1(im,j,k,e)
            d2 = dx*dx + dy*dy + dz*dz
            d2m = min(d2m,d2)

            dx = xm1(i,jp,k,e) - xm1(i,jm,k,e)
            dy = ym1(i,jp,k,e) - ym1(i,jm,k,e)
            dz = zm1(i,jp,k,e) - zm1(i,jm,k,e)
            d2 = dx*dx + dy*dy + dz*dz
            d2m = min(d2m,d2)

            dx = xm1(i,j,kp,e) - xm1(i,j,km,e)
            dy = ym1(i,j,kp,e) - ym1(i,j,km,e)
            dz = zm1(i,j,kp,e) - zm1(i,j,km,e)
            d2 = dx*dx + dy*dy + dz*dz
            d2m = min(d2m,d2)

         enddo
         enddo
         enddo
         enddo

      else  ! 2D
         do e=1,nelv
         do j=2,ny1-1
         do i=2,nx1-1

            jm=j-1
            jp=j+1
            im=i-1
            ip=i+1

            dx = xm1(ip,j,1,e) - xm1(im,j,1,e)
            dy = ym1(ip,j,1,e) - ym1(im,j,1,e)
            d2 = dx*dx + dy*dy
            d2m = min(d2m,d2)

            dx = xm1(i,jp,1,e) - xm1(i,jm,1,e)
            dy = ym1(i,jp,1,e) - ym1(i,jm,1,e)
            d2 = dx*dx + dy*dy
            d2m = min(d2m,d2)

         enddo
         enddo
         enddo
      endif

      dxmin = sqrt(d2m)/2.
      dxmin = glmin(dxmin,1)  ! Take min over all processors

      !if (nid.eq.0) write(6,*) dxmin,' this is dxmin'

      return
      end

c-----------------------------------------------------------------------
      subroutine dt_adjusted_nano                
c
c     Compute dxmin over domain
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'
      INCLUDE 'EXPONE'
      INCLUDE 'NEKUSE'
      INCLUDE 'OPCTR'
      INCLUDE 'CTIMER'
      INCLUDE 'WAKEPOT'
      INCLUDE 'MOVWIN'

      integer  n_period,n_one_period

      common /timeavg/ lst_period_1,lst_period_2
      integer          lst_period_1,lst_period_2

      w_inc   = 2.0*pi*cspeed/(param(69)*nmscale)
      period  = 2.0*pi/w_inc
      n_period= param(14)  
      fintime = n_period*period   

      nsteps  = fintime/dt       

      write(6,*) 'dt 1111    ',dt,nsteps,fintime             
      dt      = fintime/(nsteps+1)
      write(6,*) 'dt 2222    ',dt,nsteps,fintime             
      
      nsteps  = fintime/dt       
      write(6,*) 'dt 3333    ',dt,nsteps,fintime             

      n_one_period= period/dt
      write(6,*) 'dt 4444    ',dt,nsteps,fintime ,n_one_period            

      lst_period_1= nsteps-n_one_period+1
      lst_period_2= nsteps          

      if (nid.eq.0) write(6,*) 'DT adjusted',dt,nsteps,n_one_period

      return
      end
