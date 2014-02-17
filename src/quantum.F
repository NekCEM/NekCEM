c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk   
      integer  ns,nq,nsp,nqp    
      integer  nstate2
      real     tmp,tmpi  

c...  temporary parameters
      timemax = 20e-13                 ! time max
      Emin    = 1.4                    ! energy range min 
      Emax    = 1.6                    ! energy range max
      dEh     = (Emax-Emin)/(nEh+1)     
      relerr  = 1e-9                   ! relative error  
      abserr  = 1e-9                   ! absolute error 

c...  constants and units
      timeunit  = 2.418884326505e-17     ! hbar/Eh SI unit 
      timemax   = timemax/timeunit       !

      rintensity= 0.001                  ! MW/cm**2
      rintensity= rintensity*1e10        ! W/m**2
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
          write(6,*) 'h0:: diag',ii,hamilt_0(ii,ii)
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
          write(6,*) 'h0::1::',ii,jj,hamilt_0(ii,jj)
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
          write(6,*) 'h0::2::',ii,jj,hamilt_0(ii,jj)
      enddo

c...  if (ns==nsp) + { (nq==0)&(nqp==1)|(nq==1)&(nqp==0) }                              
      do  ii = 1,nstate,2 ! 3,5,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 
            
          nsp= isstate(ii+1) 
          nqp= iqstate(ii+1) 
          hamilt(ii  ,ii+1)=-q_mu       
          hamilt(ii+1,ii-1)=-q_mu       
      enddo

c...  if (nq==nqp) + { ns==(nsp-1) }                      
      do  ii = 1,nstate   ! 1,2,3,... 
          ns = isstate(ii)
          nq = iqstate(ii)

          nsp= isstate(ii+2)
          nqp= iqstate(ii+2)
          tmpi= nsp+1
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

      do ii=1,nstate
      do jj=1,nstate
         write(6,*) 'h tt:: ',ii,jj,hamilt  (ii,jj)
      enddo
      enddo

c... in rwa, each driving term energy0*cos(omega t) becomes an energy0/2 term:
      tmp= energy0*0.5 
      do ii=1,nstate
      do jj=1,nstate
         hamilt  (ii,jj)=hamilt  (ii,jj)*tmp           
         hamilt_I(ii,jj)=hamilt_0(ii,jj)+hamilt(ii,jj)
      enddo
      enddo
      do ii=1,nstate
      do jj=1,nstate
         write(6,*) 'h0:: ',ii,jj,hamilt_0(ii,jj)
      enddo
      enddo
      do ii=1,nstate
      do jj=1,nstate
         write(6,*) 'h :: ',ii,jj,hamilt  (ii,jj)
      enddo
      enddo
      do ii=1,nstate
      do jj=1,nstate
         write(6,*) 'hI:: ',ii,jj,hamilt_I(ii,jj)
      enddo
      enddo


      return
      end






