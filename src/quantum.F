c--------------------------------------------------------------------- 
      subroutine cem_quantum_parameters    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

c...  temporary parameters
      timemax = 20e-13                    ! time max
      Emin    = 1.4                       ! energy range min 
      Emax    = 1.6                       ! energy range max
      dEh     = (Emax-Emin)/(nEh+1)     
      relerr  = 1e-9                      ! relative error  
      abserr  = 1e-9                      ! absolute error 

c...  constants and units
      timeunit= 2.418884326505e-17        ! hbar/Eh SI unit 
      timemax = timemax/timeunit          !

      nqstep  = 1000 
      dt      = timemax/(nqstep-1) 
      time    = 0

      rintensity0= 0.001                  ! MW/cm**2
      rintensity= rintensity0*1e10        ! W/m**2
      energy0   = 27.450*sqrt(rintensity) ! SI(V/m) 
      energy0   = energy0/5.1421e11       ! a.u.    

      eV        = 27.21140                ! energy Eh  
      debye     = 0.3934303070      

c...  SP: surface plasmon states
      omega_s   = 1.5d0/eV
      gamma_s   = 0.2d0/eV
      s_mu      = 4095.d0*debye 

c...  QD: quantum dot
      omega_q   = 1.5d0/eV       
      gamma_q   = (1.0/1.d-9)    *timeunit
      gamma_d   = (1.d0/500.d-15)*timeunit
      q_mu      = 14.9d0*debye  
 
c...  radius in meters
      qradius   = 20e-9                  ! in meters   
      qradius   = qradius/(0.529177e-10) ! atomic units

c...  coupling constant in a.u.:
      g_couple  = 2.0/(qradius**3) 

      autime    = 2.41888e-17 

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian_setup
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk
      integer  ns,nq,nsp,nqp
      integer  nstate2
      integer  ieh,nsum,npsum
      real     tmp,tmpi


c........................................... 
c...  iqstate={0 1 0 1 0 1 0 1}
c...  isstate={0 0 1 1 2 2 3 3}
c........................................... 

      kk=0 
      do  is= 1,numsp !number of quantum states of surface plasmons
      do  iq= 1,numqd !number of quantum dots 
          kk= kk+1   
          isstate(kk)=is-1
          iqstate(kk)=iq-1
      enddo
      enddo

c........................................... 
c    (nq,nqp) 2 quantum dots ={0,1}
c........................................... 
c    (iqstate,iqstate)
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c........................................... 
c    (ns,nsp) 4 surface plasmons {0,1,2,3}
c........................................... 
c    (isstate,isstate)
c     00 , 00 , 01 , 01 , 02 , 02 , 03 , 03 
c     00 , 00 ,(01), 01 , 02 , 02 , 03 , 03 
c     10 ,(10), 11 , 11 , 12 , 12 , 13 , 13 
c     10 , 10 , 11 , 11 ,(12), 12 , 13 , 13 
c     20 , 20 , 21 ,(21), 22 , 22 , 23 , 23 
c     20 , 20 , 21 , 21 , 22 , 22 ,(23), 23 
c     30 , 30 , 31 , 31 , 32 ,(32), 33 , 33 
c     30 , 30 , 31 , 31 , 32 , 32 , 33 , 33 
c........................................... 

      nstate2 = nstate*nstate
      call rzero(hamilt_0,nstate2)
      call rzero(hamilt  ,nstate2)
      call rzero(hamilt_I,nstate2)

c...  set the diagonal terms
      do  ii = 1,nstate 
          ns = isstate(ii) 
          nq = iqstate(ii) 
          hamilt_0 (ii,ii) = ns*omega_s+nq*omega_q 
      enddo

c...  set q = 0, s couples to q' = 1, s' = s-1
      do  ii = 3,nstate,2 ! 3,5,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 

          jj = ii-1                            
          nsp= isstate(jj) 
          nqp= iqstate(jj) 
          tmpi= nsp+1
          hamilt_0(ii,jj) = -g_couple*s_mu*q_mu*sqrt(tmpi)  
      enddo

c...  q = 1, s couples to q' = 0, s' = s+1
      do  ii = 2,nstate,2 ! 2,4,6,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 

          jj = ii+1
          nsp= isstate(jj) 
          nqp= iqstate(jj) 
          tmpi= nsp
          hamilt_0(ii,jj) = -g_couple*s_mu*q_mu*dsqrt(tmpi)
      enddo

c...  if (ns==nsp) + { (nq==0)&(nqp==1)|(nq==1)&(nqp==0) }                              
      do  ii = 1,nstate,2 ! 3,5,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 
            
          nsp= isstate(ii+1) 
          nqp= iqstate(ii+1) 
          hamilt(ii  ,ii+1)=-q_mu       
          hamilt(ii+1,ii  )=-q_mu       
      enddo

c...  if (nq==nqp) + { ns==(nsp-1) }                      
      do  ii = 1,nstate     ! 1,2,3,... 
          ns = isstate(ii)
          nq = iqstate(ii)

          nsp= isstate(ii+2)
          nqp= iqstate(ii+2)
          tmpi= nsp  
          hamilt(ii,ii+2)=-sqrt(tmpi)*s_mu
      enddo

c...  if (nq==nqp) + { ns==(nsp+1) }                      
      do  ii = 3,nstate  ! 1,2,3,... 
          ns = isstate(ii)
          nq = iqstate(ii)

          nsp= isstate(ii-2)
          nqp= iqstate(ii-2)
          tmpi= nsp+1
          hamilt(ii,ii-2)=-sqrt(tmpi)*s_mu
      enddo

c... in rwa, each driving term energy0*cos(omega t) becomes an energy0/2 term:
      tmp= energy0*0.5 
      do jj= 1,nstate
      do ii= 1,nstate
         hamilt  (ii,jj)=hamilt  (ii,jj)*tmp           
         hamilt_I(ii,jj)=hamilt_0(ii,jj)+hamilt(ii,jj)
      enddo
      enddo

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
         do jj = 1,nstate
            nsp= isstate(jj)   
            nqp= iqstate(jj)   
            npsum= nsp+nqp
         do ii = 1,nstate
            ns = isstate(ii)
            nq = iqstate(ii)
            nsum= ns+nq
            frequency(ii,jj,ieh)=omega(iEh)*(nsum-npsum)
         enddo
         enddo
      enddo

c     call outmat(hamilt_0,nstate,nstate,'hamilt_0  ',nstate)
c     call outmat(hamilt  ,nstate,nstate,'hamilt    ',nstate)
c     call outmat(hamilt_I,nstate,nstate,'hamilt_I  ',nstate)
 
      return
      end

c..nnnnnnnnn
c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian(iEh)
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer    ii,jj,iEh
      real       tmprr,tmpii
      complex*16 tmpfq,tmpcc
      complex*16 ci

      ci=(0.0,1.0)

      do jj = 1,nstate
      do ii = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo

      call mxm(hamilt_I,nstate,rho_r,nstate,tmp_r1,nstate)
      call mxm(rho_r,nstate,hamilt_I,nstate,tmp_r2,nstate)
      call mxm(hamilt_I,nstate,rho_i,nstate,tmp_i1,nstate)
      call mxm(rho_i,nstate,hamilt_I,nstate,tmp_i2,nstate)

      do jj = 1,nstate
      do ii = 1,nstate
         tmpfq= frequency(ii,jj,iEh)*rho(ii,jj)        
         tmprr= tmp_r1(ii,jj)-tmp_r2(ii,jj)
         tmpii= tmp_i1(ii,jj)-tmp_i2(ii,jj)
         tmpcc= tmprr + ci*tmpii
         drho(ii,jj)= - ci*(tmpcc - tmpfq) 
         resrho_r(ii,jj)=real(drho(ii,jj) )    
         resrho_i(ii,jj)=imag(drho(ii,jj) )    
      enddo
      enddo
      
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_linblad
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  ii,jj,ns,nq,nsp,nqp, kns,knsp
      real     fac,tmp,tmp1,tmp2,tmp3

      do jj = 1,nstate
         nsp= isstate(jj)
         nqp= iqstate(jj)
      do ii = 1,nstate
         ns = isstate(ii)
         nq = iqstate(ii)
         tmp1=gamma_d*(nq+nqp-2*nq*nqp)
         tmp2=gamma_q*(nq+nqp)*0.5
         tmp3=gamma_s*(ns+nsp)*0.5
c...  QD dephasing Linblad term leads to: 
         drho(ii,jj)= drho(ii,jj)-tmp1*rho(ii,jj)
c...  QD spontaneous emission Linblad term leads to these two terms:
         drho(ii,jj)= drho(ii,jj)-tmp2*rho(ii,jj)
c...  if(iqp.eq.0.and.iq.eq.0)then  
         if ((nqp.eq.0).and.(nq.eq.0)) then
         drho(ii,jj)= drho(ii,jj)+gamma_q*rho(ii+1,jj+1)
         endif 
c...  SP spontaneous emission Linblad terms:
         drho(ii,jj)= drho(ii,jj)-tmp3*rho(ii,jj)

         kns = (ns +2)*2
         if (nq .eq.0) kns = 2*(ns +1)+1
         knsp= (nsp+2)*2
         if (nqp.eq.0) knsp= 2*(nsp+1)+1

         if ((kns.lt.nstate).and.(knsp.lt.nstate)) then
              tmp = (ns+1)*(nsp+1)
              fac = gamma_s*sqrt(tmp)
              drho(ii,jj)= drho(ii,jj)+fac*rho(kns,knsp)
         endif
         resrho_r(ii,jj)= real(drho(ii,jj))
         resrho_i(ii,jj)= imag(drho(ii,jj))
      enddo
      enddo

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_propagation_exp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh, ii  
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7


      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup_p           
 
      stime0= dclock()
      do iEh = 1,nEh
         call cem_quantum_initialize

         stime1= dclock()
         do istep = 1,nqstep
            stime2= dclock()
            call cem_quantum_exp(iEh)
            time= time+dt
            call cem_quantum_energy(iEh)     
            stime2= dclock()-stime2
         enddo
         stime1= dclock()-stime1
         call cem_quantum_post(iEh)            
      enddo
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_propagation_rk
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh, ii  
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7


      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup           
      call cem_quantum_rk_storage
 
      stime0=dclock()
      do iEh = 1,nEh
         call cem_quantum_initialize

         stime1= dclock()
         do istep= 1,nqstep
            stime2= dclock()
            do ii= 1,5
               call rk_quantum_c (ii)
               call cem_quantum_hamiltonian(iEh)         
               call cem_quantum_linblad    
               call rk_quantum_ab(ii)
            enddo
            time= time+dt
            call cem_quantum_energy(iEh)     
            stime2= dclock()-stime2
         enddo
         stime1= dclock()-stime1
         call cem_quantum_post(iEh)            
      enddo
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_propagation_rk_p
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh, ii
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7


      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup_p
      call cem_quantum_rk_storage

      stime0=dclock()
      do iEh = 1,nEh
         call cem_quantum_initialize

         stime1= dclock()
         do istep= 1,nqstep
            stime2= dclock()
            do ii= 1,5
               call rk_quantum_c (ii)
               call cem_quantum_hamiltonian_p(iEh)
               call cem_quantum_linblad_p
               call rk_quantum_ab(ii)
            enddo
            time= time+dt
            call cem_quantum_energy(iEh)
            stime2= dclock()-stime2
         enddo
         stime1= dclock()-stime1
         call cem_quantum_post(iEh)
      enddo
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_quantum_c(i)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      include 'RK5'

      integer  i

      RKtime= time+dt*rk4c(i)

      writE(6,6) dt, time, RKtime
  6   format('dt/time/RKtime:: ',1p3e20.9)

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_quantum_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      include 'RK5'
      real     ca,cb
      integer  ii,nstate2

      nstate2= nstate*nstate

      ca = rk4a(ii)
      cb = rk4b(ii)

      call rk4_upd(rho_r,k_rho_r,resrho_r,cb,ca,dt,nstate2)
      call rk4_upd(rho_i,k_rho_i,resrho_i,cb,ca,dt,nstate2)

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_initialize    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  nstate2

      nstate2= nstate*nstate

      call rzero(   rho_r,nstate2)
      call rzero(   rho_i,nstate2)
      call rzero( k_rho_r,nstate2)
      call rzero( k_rho_i,nstate2)
      call rzero(resrho_r,nstate2)
      call rzero(resrho_i,nstate2)

      rho_r(1,1)=1.0

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_energy (iEh)       
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj    
      complex*16 ci, tmpcc
      complex*16 cq_mu,cs_mu,ensp,enqd

      ci = (0.0,1.0)

      do jj=1,nstate
      do ii=1,nstate
         tmpcc     = exp(-ci*frequency(ii,jj,iEh)*time)
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))       
         rho(ii,jj)= rho(ii,jj)*tmpcc
      enddo
      enddo
      
c...  expect values of muq,mus in Debye:
      cq_mu=0.0 
      do ii=1,nstate,2 
         cq_mu=cq_mu+rho(ii,ii+1)+rho(ii+1,ii) 
      enddo 
      cq_mu=cq_mu*q_mu/0.3934303070 

      cs_mu=0.0
      do ii=1,nstate-2
         cs_mu=cs_mu+rho(ii,ii+2)+rho(ii+2,ii) 
      enddo 
      cs_mu=cs_mu*s_mu/0.3934303070 

c...  energies 
      enqd=0.0 
      do ii=1,nstate,2 
         enqd=enqd+real(rho(ii,ii))*omega_q
      enddo              

      ensp=0.0 
      do ii=1,nstate   
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! odd+even
      enddo

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_post (iEh)               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj    
      complex*16 ci,ensp
      real     rinc,sigcon,sig 
      real     qd_mu,sp_mu
      real     enqd
      real     auerg,wattocgs,xqd,xsp        

      ci = (0.0,1.0)

      do jj=1,nstate
      do ii=1,nstate
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))       
      enddo
      enddo

c...  expect values of muq,mus in Debye:
      qd_mu=0.0 
      do ii=1,nstate,2
         qd_mu=qd_mu+imag(rho(ii,ii+1))
      enddo
      qd_mu=-2.0*qd_mu*q_mu

      sp_mu=0.0
      do ii=1,nstate-2,2
         sp_mu=sp_mu+imag(rho(ii,ii+2))
      enddo
      sp_mu=-2.0*sp_mu*s_mu

      rinc  =(137.d0)*energy0/(4.0*3.14159265*omega(iEh))
      sigcon=(0.529177e-8)**2 
      sig   =(qd_mu+sp_mu)*sigcon/rinc  
      write(6,*) 

c... Also calculate avg energies times decay constants, normalized by intensity. 
c... This is possibly a crude measure of sig as well.
      enqd=0.0
      do ii=2,nstate,2
         enqd=enqd+real(rho(ii,ii))*omega_q
      enddo

      ensp=0.0
      do ii=1,nstate    
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! odd+even
      enddo

c...  convert to cross sections in cgs:
      auerg=4.35975e-11 
      wattocgs=1.7
      xqd=(gamma_q+gamma_d)*enqd*auerg/(1.d6*rintensity0*wattocgs)
      xsp=gamma_s*ensp*auerg/(1.6*rintensity0*wattocgs)        

c...  need to convert the time units from inv au to sec:
      xqd=xqd/autime     
      xsp=xsp/autime     
      write(10,10) hb_omega(iEh),sig,real(ensp),
     $            imag(ensp),enqd,real(ensp+enqd)
      write(11,11) hb_omega(iEh),sig

 10   format(1p6e20.8)
 11   format(1p2e20.8)

      return
      end        
c--------------------------------------------------------------------- 
      subroutine cem_quantum_exp(iEh)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  icalld
      save     icalld
      data     icalld/0/
      real     glsc2,beta,tmp
      integer  m,info
      integer  iEh,nstate2,n2


      qmarnol= larnol 
      m      = 5 !qmarnol

      nstate2= nstate*nstate
      n2     = 2*nstate2 

      if (icalld.eq.0) then
          call rzero(qu ,n2)
          call rzero(qv0,n2*(m+1))
          call rzero(qv ,n2*m)
          call rzero(qvt,n2*m)
          call rzero(qh0,(m+1)*m)
          call rzero(qh ,m*m)
          call rzero(qw ,m)
          call rzero(qws,2*m*(m+2))
          icalld = 1
      endif

      call copy(qu(0*nstate2+1),rho_r,nstate2)
      call copy(qu(1*nstate2+1),rho_i,nstate2)

      beta= glsc2(qu,qu,n2)
      beta= sqrt (beta)

      call quantum_arnoldi_vec (qv0,qh0,n2,m,qu,iEh)
      call quantum_arnoldi_mat (qvt,qv,qh,qv0,qh0,m,n2)

      call copy_r2z   (qzh,qh,m*m)
      call zgeev      ('N','V',m,qzh,m,qeval,qzri,m,qzr,m,qawork,
     $                  2*m,qrwork,info)
      call copy_z2z   (qzri,qzr,m*m)

      call zgetrf     (m,m,qzri,m,iqpiv,info)        !vri= LU factorization                
      call zgetri     (m,qzri,m,iqpiv,qawork,m,info) !vri= X^(-1)

      call z_exp      (qeval,dt,m)                   !eval=exp(dt*eval)
      call z_rowscale (qzri,qeval,m)                 !vri =eval*X^(-1); 
      call c_mxm      (qzr,m,qzri,m,qzh,m)           !zh =X*vri  

      call copy_z2r   (qw,qzh(1,1),m)                 !zh(i,1) (i=1:m)
      call mxm        (qv,n2,qw,m,qu,1)
      call cmult      (qu,beta,n2)

      call copy(rho_r,qu(0*nstate2+1),nstate2)
      call copy(rho_i,qu(1*nstate2+1),nstate2)


      return
      end
c--------------------------------------------------------------------- 
      subroutine quantum_arnoldi_vec(v,h,n,m,b,iEh)
      implicit none
      integer  i,j,m,n,iEh
      real     v(n,m+1),b(n)
      real     h(m+1,m)
      real     eps, s1, si, glsc2, cnt

      eps = 1.e-13

      s1 = glsc2(b,b,n)
      si = 1./sqrt(s1) ! should add a check and a comment here

      call cmult2(v(1,1),b,si,n)

      do j=1,m

         call quantum_mult (v(1,j+1),v(1,j),n,iEh)

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
      subroutine quantum_arnoldi_mat (vt,v,h,v0,h0,m,n)
      implicit none
      integer  m,n,i,j,k
      real     h0(m+1,m),h(m,m)
      real     v0(n,m+1)
      real     v (n,m)
      real     vt(m,n)

      do j=1,m
      do i=1,n
         v (i,j)= v0(i,j)
         vt(j,i)= v0(i,j)
      enddo
      enddo

      do i=1,m
      do j=1,m
         h(i,j) = h0(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine quantum_mult(w,v,n2,iEh)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM' 
      integer  iEh,n2,nstate2
      real     w(n2),v(n2)

      nstate2 = nstate*nstate

      call copy(rho_r, v(1+0*nstate2), nstate2)
      call copy(rho_i, v(1+1*nstate2), nstate2)

      call cem_quantum_hamiltonian(iEh)
      call cem_quantum_linblad

      call copy(w(1+0*nstate2), resrho_r, nstate2) 
      call copy(w(1+1*nstate2), resrho_i, nstate2) 

      return
      end

cppppppp
c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian_setup_p
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk
      integer  ns,nq,nsp,nqp
      integer  nstate3
      integer  ieh,nsum,npsum
      real     tmp,tmpi


c........................................... 
c...  iqstate={0 1 0 1 0 1 0 1}
c...  isstate={0 0 1 1 2 2 3 3}
c........................................... 

      kk=0 
      do  is= 1,numsp !number of quantum states of surface plasmons
      do  iq= 1,numqd !number of quantum dots 
          kk= kk+1   
          isstate(kk)=is-1
          iqstate(kk)=iq-1
      enddo
      enddo

c........................................... 
c    (nq,nqp) 2 quantum dots ={0,1}
c........................................... 
c    (iqstate,iqstate)
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c     00 ,(01), 00 ,(01), 00 ,(01), 00 ,(01)
c    (10), 11 ,(10), 11 ,(10), 11 ,(10), 11 
c........................................... 
c    (ns,nsp) 4 surface plasmons {0,1,2,3}
c........................................... 
c    (isstate,isstate)
c     00 , 00 , 01 , 01 , 02 , 02 , 03 , 03 
c     00 , 00 ,(01), 01 , 02 , 02 , 03 , 03 
c     10 ,(10), 11 , 11 , 12 , 12 , 13 , 13 
c     10 , 10 , 11 , 11 ,(12), 12 , 13 , 13 
c     20 , 20 , 21 ,(21), 22 , 22 , 23 , 23 
c     20 , 20 , 21 , 21 , 22 , 22 ,(23), 23 
c     30 , 30 , 31 , 31 , 32 ,(32), 33 , 33 
c     30 , 30 , 31 , 31 , 32 , 32 , 33 , 33 
c........................................... 

      nstate3 = 3*numsp*numqd
      call rzero(hamiltp_0,nstate3)
      call rzero(hamiltp  ,nstate3)
      call rzero(hamiltp_0,nstate3)

c...  set the diagonal terms
      do  ii = 1,nstate 
          ns = isstate(ii) 
          nq = iqstate(ii) 
          kk = 3*(ii-1)+1
          hamiltp_0 (kk) = ns*omega_s+nq*omega_q 
          write(6,*) 'hamp_0', kk,hamiltp_0 (kk)
      enddo

c...  q = 1, s couples to q' = 0, s' = s+1
      do  ii = 2,nstate,2 ! 2,4,6,... 
          ns  = isstate(ii  ) 
          nq  = iqstate(ii  ) 
          nsp = isstate(ii+1) 
          nqp = iqstate(ii+1) 
          tmpi= nsp
          kk  = 3*ii-1
          hamiltp_0(kk) = -g_couple*s_mu*q_mu*dsqrt(tmpi)
          write(6,*) 'hamp_0', kk,hamiltp_0 (kk)
      enddo

c...  if (ns==nsp) + { (nq==0)&(nqp==1)|(nq==1)&(nqp==0) }                              
      do  ii = 1,nstate,2 ! 1,3,5,... 
          ns = isstate(ii  ) 
          nq = iqstate(ii  ) 
          nsp= isstate(ii+1) 
          nqp= iqstate(ii+1) 
          kk = 3*(ii-1)+2
          hamiltp(kk)= -q_mu       
          write(6,*) 'hamp_0', kk,hamiltp_0 (kk)
      enddo

c...  if (nq==nqp) + { ns==(nsp-1) }                      
c...  if (nq==nqp) + { ns==(nsp+1) }                      
      do  ii = 1,nstate     ! 1,2,3,... 
          ns = isstate(ii  )
          nq = iqstate(ii  )
          nsp= isstate(ii+2)
          nqp= iqstate(ii+2)
          tmpi= nsp  
          kk  = 3*ii
          hamiltp(kk)=-sqrt(tmpi)*s_mu
      enddo


c... in rwa, each driving term energy0*cos(omega t) becomes an energy0/2 term:
      tmp= energy0*0.5 
      do ii= 1,nstate3
         hamiltp  (ii)=hamiltp(ii)*tmp           
         hamiltp_I(ii)=hamiltp_0(ii)+hamiltp(ii)
         write(6,*) 'hamp_I', ii,hamiltp_I (ii)
      enddo

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
         do jj = 1,nstate
            nsp= isstate(jj)   
            nqp= iqstate(jj)   
            npsum= nsp+nqp
         do ii = 1,nstate
            ns = isstate(ii)
            nq = iqstate(ii)
            nsum= ns+nq
            frequency(ii,jj,ieh)=omega(iEh)*(nsum-npsum)
         enddo
         enddo
      enddo

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian_p(iEh)
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer    ii,jj,iEh
      real       tmprr,tmpii
      complex*16 tmpfq,tmpcc
      complex*16 ci

      ci=(0.0,1.0)

      do jj = 1,nstate
      do ii = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo

      call cem_quantum_mxm(rho_r,tmp_r1, 1)             
      call cem_quantum_mxm(rho_r,tmp_r2,-1)             

      call cem_quantum_mxm(rho_i,tmp_i1, 1)             
      call cem_quantum_mxm(rho_i,tmp_i2,-1)             

      do jj = 1,nstate
      do ii = 1,nstate
         tmpfq= frequency(ii,jj,iEh)*rho(ii,jj)        
         tmprr= tmp_r1(ii,jj)-tmp_r2(ii,jj)
         tmpii= tmp_i1(ii,jj)-tmp_i2(ii,jj)
         tmpcc= tmprr + ci*tmpii
         drho(ii,jj)= - ci*(tmpcc - tmpfq) 
         resrho_r(ii,jj)=real(drho(ii,jj) )    
         resrho_i(ii,jj)=imag(drho(ii,jj) )    
      enddo
      enddo
      
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_mxm(rrho,tmprr,idir)             
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer    ii,jj,idir    
      real       rrho (nstate,nstate)
      real       tmprr(nstate,nstate)


c     do ii=1,nstate*3
c        write(6,*) 'h:',ii,hamiltp_I(ii)
c     enddo

c...  (idir = 1)
      if (idir.eq.1) then

      do jj = 1,nstate
      do ii = 1,1
         tmprr(ii,jj)=  hamiltp_I(3*(ii-1)+1)*rrho(ii  ,jj)
     $                + hamiltp_I(3*(ii-1)+2)*rrho(ii+1,jj)
     $                + hamiltp_I(3*(ii-1)+3)*rrho(ii+2,jj)
      enddo
      enddo

      do jj = 1,nstate
      do ii = 2,2
         tmprr(ii,jj)=  hamiltp_I(3*(ii-1)-1)*rrho(ii-1,jj)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(ii  ,jj)
     $                + hamiltp_I(3*(ii-1)+2)*rrho(ii+1,jj)
     $                + hamiltp_I(3*(ii-1)+3)*rrho(ii+2,jj)
      enddo
      enddo

      do jj = 1,nstate
      do ii = 3,nstate-2
         tmprr(ii,jj)=  hamiltp_I(3*(ii-1)-3)*rrho(ii-2,jj)
     $                + hamiltp_I(3*(ii-1)-1)*rrho(ii-1,jj)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(ii  ,jj)
     $                + hamiltp_I(3*(ii-1)+2)*rrho(ii+1,jj)
     $                + hamiltp_I(3*(ii-1)+3)*rrho(ii+2,jj)
      enddo
      enddo

      do jj = 1,nstate
      do ii = nstate-1,nstate-1
         tmprr(ii,jj)=  hamiltp_I(3*(ii-1)-3)*rrho(ii-2,jj)
     $                + hamiltp_I(3*(ii-1)-1)*rrho(ii-1,jj)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(ii  ,jj)
     $                + hamiltp_I(3*(ii-1)+2)*rrho(ii+1,jj)
      enddo
      enddo

      do jj = 1,nstate
      do ii = nstate,nstate
         tmprr(ii,jj)=  hamiltp_I(3*(ii-1)-3)*rrho(ii-2,jj)
     $                + hamiltp_I(3*(ii-1)-1)*rrho(ii-1,jj)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(ii  ,jj)
      enddo
      enddo

c...  (idir = -1)
      elseif (idir.eq.-1) then

      do jj = 1,nstate
      do ii = 1,1
         tmprr(jj,ii)=  hamiltp_I(3*(ii-1)+1)*rrho(jj,ii  )
     $                + hamiltp_I(3*(ii-1)+2)*rrho(jj,ii+1)
     $                + hamiltp_I(3*(ii-1)+3)*rrho(jj,ii+2)
      enddo
      enddo

      do jj = 1,nstate
      do ii = 2,2
         tmprr(jj,ii)=  hamiltp_I(3*(ii-1)-1)*rrho(jj,ii-1)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(jj,ii  )
     $                + hamiltp_I(3*(ii-1)+2)*rrho(jj,ii+1)
     $                + hamiltp_I(3*(ii-1)+3)*rrho(jj,ii+2)
      enddo
      enddo

      do jj = 1,nstate
      do ii = 3,nstate-2
         tmprr(jj,ii)=  hamiltp_I(3*(ii-1)-3)*rrho(jj,ii-2)
     $                + hamiltp_I(3*(ii-1)-1)*rrho(jj,ii-1)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(jj,ii  )
     $                + hamiltp_I(3*(ii-1)+2)*rrho(jj,ii+1)
     $                + hamiltp_I(3*(ii-1)+3)*rrho(jj,ii+2)
      enddo
      enddo

      do jj = 1,nstate
      do ii = nstate-1,nstate-1
         tmprr(jj,ii)=  hamiltp_I(3*(ii-1)-3)*rrho(jj,ii-2)
     $                + hamiltp_I(3*(ii-1)-1)*rrho(jj,ii-1)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(jj,ii  )
     $                + hamiltp_I(3*(ii-1)+2)*rrho(jj,ii+1)
      enddo
      enddo

      do jj = 1,nstate
      do ii = nstate,nstate
         tmprr(jj,ii)=  hamiltp_I(3*(ii-1)-3)*rrho(jj,ii-2)
     $                + hamiltp_I(3*(ii-1)-1)*rrho(jj,ii-1)
     $                + hamiltp_I(3*(ii-1)+1)*rrho(jj,ii  )
      enddo
      enddo

      else

       write(6,*) 'no idir: cem_quantum_mxm'
       call exitt

      endif

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_linblad_p
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  ii,jj,ns,nq,nsp,nqp, kns,knsp
      real     fac,tmp,tmp1,tmp2,tmp3

      do jj = 1,nstate
         nsp= isstate(jj)
         nqp= iqstate(jj)
      do ii = 1,nstate
         ns = isstate(ii)
         nq = iqstate(ii)
         tmp1=gamma_d*(nq+nqp-2*nq*nqp)
         tmp2=gamma_q*(nq+nqp)*0.5
         tmp3=gamma_s*(ns+nsp)*0.5
c...  QD dephasing Linblad term leads to: 
         drho(ii,jj)= drho(ii,jj)-tmp1*rho(ii,jj)
c...  QD spontaneous emission Linblad term leads to these two terms:
         drho(ii,jj)= drho(ii,jj)-tmp2*rho(ii,jj)
c...  if(iqp.eq.0.and.iq.eq.0)then  
         if ((nqp.eq.0).and.(nq.eq.0)) then
         drho(ii,jj)= drho(ii,jj)+gamma_q*rho(ii+1,jj+1)
         endif 
c...  SP spontaneous emission Linblad terms:
         drho(ii,jj)= drho(ii,jj)-tmp3*rho(ii,jj)

         kns = (ns +2)*2
         if (nq .eq.0) kns = 2*(ns +1)+1
         knsp= (nsp+2)*2
         if (nqp.eq.0) knsp= 2*(nsp+1)+1

         if ((kns.lt.nstate).and.(knsp.lt.nstate)) then
              tmp = (ns+1)*(nsp+1)
              fac = gamma_s*sqrt(tmp)
              drho(ii,jj)= drho(ii,jj)+fac*rho(kns,knsp)
         endif
         resrho_r(ii,jj)= real(drho(ii,jj))
         resrho_i(ii,jj)= imag(drho(ii,jj))
      enddo
      enddo

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_initialize_p    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  nstate2

      nstate2 = nstate*nstate  

      call rzero(   rho_r,nstate2)
      call rzero(   rho_i,nstate2)
      call rzero( k_rho_r,nstate2)
      call rzero( k_rho_i,nstate2)
      call rzero(resrho_r,nstate2)
      call rzero(resrho_i,nstate2)

      rho_r(1,1)=1.0

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_energy_p (iEh)       
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj    
      complex*16 ci, tmpcc
      complex*16 cq_mu,cs_mu,ensp,enqd

      ci = (0.0,1.0)

      do jj=1,nstate
      do ii=1,nstate
         tmpcc     = exp(-ci*frequency(ii,jj,iEh)*time)
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))       
         rho(ii,jj)= rho(ii,jj)*tmpcc
      enddo
      enddo
      
c...  expect values of muq,mus in Debye:
      cq_mu=0.0 
      do ii=1,nstate,2 
         cq_mu=cq_mu+rho(ii,ii+1)+rho(ii+1,ii) 
      enddo 
      cq_mu=cq_mu*q_mu/0.3934303070 

      cs_mu=0.0
      do ii=1,nstate-2
         cs_mu=cs_mu+rho(ii,ii+2)+rho(ii+2,ii) 
      enddo 
      cs_mu=cs_mu*s_mu/0.3934303070 

c...  energies 
      enqd=0.0 
      do ii=1,nstate,2 
         enqd=enqd+real(rho(ii,ii))*omega_q
      enddo              

      ensp=0.0 
      do ii=1,nstate   
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! odd+even
      enddo

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_post_p (iEh)               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj    
      complex*16 ci,ensp
      real     rinc,sigcon,sig 
      real     qd_mu,sp_mu
      real     enqd
      real     auerg,wattocgs,xqd,xsp        

      ci = (0.0,1.0)

      do jj=1,nstate
      do ii=1,nstate
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))       
      enddo
      enddo

c...  expect values of muq,mus in Debye:
      qd_mu=0.0 
      do ii=1,nstate,2
         qd_mu=qd_mu+imag(rho(ii,ii+1))
      enddo
      qd_mu=-2.0*qd_mu*q_mu

      sp_mu=0.0
      do ii=1,nstate-2,2
         sp_mu=sp_mu+imag(rho(ii,ii+2))
      enddo
      sp_mu=-2.0*sp_mu*s_mu

      rinc  =(137.d0)*energy0/(4.0*3.14159265*omega(iEh))
      sigcon=(0.529177e-8)**2 
      sig   =(qd_mu+sp_mu)*sigcon/rinc  
      write(6,*) 

c... Also calculate avg energies times decay constants, normalized by intensity. 
c... This is possibly a crude measure of sig as well.
      enqd=0.0
      do ii=2,nstate,2
         enqd=enqd+real(rho(ii,ii))*omega_q
      enddo

      ensp=0.0
      do ii=1,nstate    
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! odd+even
      enddo

c...  convert to cross sections in cgs:
      auerg=4.35975e-11 
      wattocgs=1.7
      xqd=(gamma_q+gamma_d)*enqd*auerg/(1.d6*rintensity0*wattocgs)
      xsp=gamma_s*ensp*auerg/(1.6*rintensity0*wattocgs)        

c...  need to convert the time units from inv au to sec:
      xqd=xqd/autime     
      xsp=xsp/autime     
      write(10,10) hb_omega(iEh),sig,real(ensp),
     $            imag(ensp),enqd,real(ensp+enqd)
      write(11,11) hb_omega(iEh),sig

 10   format(1p6e20.8)
 11   format(1p2e20.8)

      return
      end        
c--------------------------------------------------------------------- 
      subroutine cem_quantum_rk_storage
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RK5'


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

      return
      end
c--------------------------------------------------------------------- 
