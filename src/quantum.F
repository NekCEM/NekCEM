c--------------------------------------------------------------------- 
      subroutine cem_quantum_parameters    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

c...  temporary parameters
      timemax = 20e-13                 ! time max
      Emin    = 1.4                    ! energy range min 
      Emax    = 1.6                    ! energy range max
      dEh     = (Emax-Emin)/(nEh+1)     
      relerr  = 1e-9                   ! relative error  
      abserr  = 1e-9                   ! absolute error 

c...  constants and units
      timeunit= 2.418884326505e-17     ! hbar/Eh SI unit 
      timemax = timemax/timeunit       !

      dt      = timemax/(2000-1) 
      time    = 0

      rintensity0= 0.001                  ! MW/cm**2
      rintensity= rintensity0*1e10        ! W/m**2
      energy0   = 27.450*sqrt(rintensity)! SI(V/m) 
      energy0   = energy0/5.1421e11      ! a.u.    

      eV        = 27.21140               ! energy Eh  
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
      qradius    = 20e-9                  ! in meters   
      qradius    = qradius/(0.529177e-10) ! atomic units

c...  coupling constant in a.u.:
      g_couple  = 2.0/(qradius**3) 

      autime=2.41888e-17 

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
      call rzero(hamilt_0,nstate2)

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
      do ii= 1,nstate
      do jj= 1,nstate
         hamilt  (ii,jj)=hamilt  (ii,jj)*tmp           
         hamilt_I(ii,jj)=hamilt_0(ii,jj)+hamilt(ii,jj)
      enddo
      enddo

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh) = Emin+iEh*dEh
         hb_omega(iEh) = omega(iEh)   
         omega(iEh) = omega(iEh)/eV   
         do ii = 1,nstate
            ns = isstate(ii)
            nq = iqstate(ii)
            nsum= ns+nq
         do jj = 1,nstate
            nsp= isstate(jj)   
            nqp= iqstate(jj)   
            npsum= nsp+nqp
            frequency(ii,jj,ieh)=omega(iEh)*(nsum-npsum)
         enddo
         enddo
      enddo
 
      return
      end

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

      do ii = 1,nstate
      do jj = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo

      call mxm(hamilt_I,nstate,rho_r,nstate,tmp_r1,nstate)
      call mxm(rho_r,nstate,hamilt_I,nstate,tmp_r2,nstate)
      call mxm(hamilt_I,nstate,rho_i,nstate,tmp_i1,nstate)
      call mxm(rho_i,nstate,hamilt_I,nstate,tmp_i2,nstate)

      do ii = 1,nstate
      do jj = 1,nstate
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


      do ii = 1,nstate
         ns = isstate(ii)
         nq = iqstate(ii)
      do jj = 1,nstate
         nsp= isstate(jj)
         nqp= iqstate(jj)
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

      call outmat(resrho_r,nstate,nstate,'resrho_r  ',nstate)
      call outmat(resrho_i,nstate,nstate,'resrho_i  ',nstate)

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_propagation
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh, ii

      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup           
      call rk_storage
 
      do iEh = 1,10 
         call cem_quantum_initialize
         do istep = 1,2000
            do ii=1,5
               call rk_quantum_c (ii)
               call cem_quantum_hamiltonian(iEh)         
               call cem_quantum_linblad    
               call rk_quantum_ab(ii)
            enddo
           !call outmat(rho_r,nstate,nstate,'rho_r  ',nstate)
           !call outmat(rho_i,nstate,nstate,'rho_i  ',nstate)
            time=time+dt
            call cem_quantum_energy(iEh)     
         enddo
         call cem_quantum_post(iEh)            
      enddo

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

      writE(6,*) 'dt/time/RKtime', dt, time, RKtime

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
c     call outmat(resrho_r,nstate,nstate,'resrho_r  ',nstate)
c     call outmat(resrho_i,nstate,nstate,'resrho_i  ',nstate)
c     call outmat(k_rho_r,nstate,nstate,'k_rho_r  ',nstate)
c     call outmat(k_rho_i,nstate,nstate,'k_rho_i  ',nstate)
c     call outmat(rho_r,nstate,nstate,'rho_r  ',nstate)
c     call outmat(rho_i,nstate,nstate,'rho_i  ',nstate)
c     call exitt

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

      do ii=1,nstate
      do jj=1,nstate
         tmpcc     = exp(-ci*frequency(ii,jj,iEh)*time)
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))       
         rho(ii,jj)= rho(ii,jj)*tmpcc
      enddo
      enddo
c     call outmat(rho_r,nstate,nstate,'rho_r  ',nstate)
c     call outmat(rho_i,nstate,nstate,'rho_i  ',nstate)
c     call exitt
      
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
      do ii=ii,nstate,2 
         enqd=enqd+real(rho(ii,ii))*omega_q
      enddo              

      ensp=0.0 
      do ii=1,nstate,2 
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! odd
      enddo

      do ii=2,nstate,2
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! even
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

      do ii=1,nstate
      do jj=1,nstate
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
      do ii=ii,nstate,2
         enqd=enqd+real(rho(ii,ii))*omega_q
      enddo

      ensp=0.0
      do ii=1,nstate,2
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! odd
      enddo

      do ii=2,nstate,2
         ensp=ensp+rho(ii,ii)*omega_s*(ii-1.0)/2.0 ! even
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
