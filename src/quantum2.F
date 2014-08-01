
c---------------------------------------------------------------------
      subroutine cem_quantum2_hamiltonian_setup_pls
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd)
      integer  binary_rep
      integer  iEh,nsum,npsum,nstate_min,ierr
      real*8   tmp,tmpi
      logical  q_muFlag,s_muFlag
      

      ltonq = level**numqd
c     Do the minimal binary_rep calls
      do ii=1,2*ltonq
         do kk=1,numqd
            nq_minimal(kk,ii)=binary_rep(kk,mod(ii-1,ltonq)+1)
         enddo
      enddo
c     Initially, generate upper half of s=0 part of the hamiltonian in normal way

c     the number needed to store in a single row is
c     3*2**Nqd - (2**Nqd + 1)

      nstate_min = n_in_Row*ltonq
      call rzero(ham_sys_min,nstate_min)
      call rzero(ham_qmu_min,nstate_min)
      call rzero(ham_smu_min,nstate_min)
      call rzero(ham_diag,ltonq)

c     2*level**numberOfQuantumd dots is because we want to all q possibilities
c     for the values of 0 and 1 for s
      do ii=1,ltonq
         ns=(ii-1)/ltonq
         do iq=1,numqd
            nq(iq)=nq_minimal(iq,ii)
         enddo
         do jj=1,2*ltonq
            do iq=1,numqd
               nqp(iq) = nq_minimal(iq,jj)
            enddo
            s_muFlag=.true.
            nsp = (jj-1)/ltonq
            if(jj.ge.ii) then
c to get the linearized sparse matrix position from the row/column location, 
c we take (ii+(jj-1)*n_in_row-(jj-1)
               kk = jj+(ii-1)*n_in_row-(ii-1)

                 
               do iq=1,numqd
                  q_muFlag=.true.
                  do iiq=1,numqd
c     We now have to loop over the nq to check if
c     all of the other states are equal
                     if(iq.ne.iiq.and.nq(iiq).ne.nqp(iiq)) then
                        q_muFlag=.false.
                     endif
                  enddo
                  if (q_muFlag) then

c...  q, s couples to q'=q-1, s'=s+1
                     if(nqp(iq).eq.nq(iq)-1.and.nsp.eq.ns+1)
     $                    ham_sys_min(kk) = ham_sys_min(kk)
     $                    -g_couple_m(iq)*dsqrt(dfloat((ns+1)*nq(iq)))
                  
c...  q, s couples to q'=q+1, s'=s-1
                     if(nqp(iq).eq.nq(iq)+1.and.nsp.eq.ns-1)
     $                    ham_sys_min(kk) = ham_sys_min(kk)
     $                    -g_couple_m(iq)*dsqrt(dfloat(ns*(nq(iq)+1)))
  
                  
c     now do hamilt parts
                     if(ns.eq.nsp) then
                        if(nqp(iq).eq.nq(iq)+1) then 
                           ham_qmu_min(kk) = ham_qmu_min(kk)
     $                          -q_mu_m(iq)*dsqrt(dfloat(nq(iq)+1))
                        endif
                        if(nqp(iq).eq.nq(iq)-1) then
                           ham_qmu_min(kk) = ham_qmu_min(kk)
     $                          -q_mu_m(iq)*dsqrt(dfloat(nq(iq)))
                        endif
                     endif
                  endif
                  
                  if(nq(iq).ne.nqp(iq)) then
                     s_muFlag=.false.
                  endif
               enddo
               if(s_muFlag) then
                  if(nsp.eq.(ns+1)) then
                     ham_smu_min(kk)=ham_smu_min(kk)-
     $                    dsqrt(dfloat(ns+1))*s_mu
                     

                  endif
                  
                  if(nsp.eq.(ns-1)) then
                     ham_smu_min(kk)=ham_smu_min(kk)-
     $                   dsqrt(dfloat(ns))*s_mu
                  endif
                 
               endif
            endif
         enddo
      enddo

c     set ham_diag
      do ii=1,ltonq
         do kk=1,numqd
            ham_diag(ii) = ham_diag(ii)+nq_minimal(kk,ii)*omega_q_m(kk)
         enddo
      enddo

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
      enddo
      do jj = 1,nstate/np
         npsum = (jj+c_offset-1)/ltonq
         do iq=1,numqd
            npsum = npsum + 
     $           nq_minimal(iq,mod(jj+c_offset-1,ltonq)+1)
         enddo
         
         do ii = 1,nstate
            nsum = (ii-1)/ltonq
            do iq=1,numqd
               nsum = nsum + 
     $              nq_minimal(iq,mod(ii-1,ltonq)+1)
            enddo
            state_sums(ii,jj)=(nsum-npsum)
         enddo
      enddo

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum2_initialize    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'
      integer  nstate2,IERR,ii,jj

      nstate2= nstate*nstate/np

      call czero(s_mu_frequency,nEh)
      call czero(q_mu_frequency,nEh)
      call czero(energy_t,      nEh)
      call rzero(   rho_r,nstate2)
      call rzero(   rho_i,nstate2)
      call rzero( k_rho_r,nstate2)
      call rzero( k_rho_i,nstate2)
      call rzero(resrho_r,nstate2)
      call rzero(resrho_i,nstate2)

c      do jj=1,nstate/np
c         do ii=1,nstate
c            rho_r(ii,jj) = 1
c            rho_i(ii,jj) = 1
c         enddo
c      enddo
      if(nid.eq.0) rho_r(1,1) = 1

      qctime1= 0.0
      qctime2= 0.0
      qctime3= 0.0
      qctime4= 0.0

      qrtime1= 0.0
      qrtime2= 0.0
      qrtime3= 0.0
      qrtime4= 0.0

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum_gs
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'
      integer ii,jj,kk,IERR
      

c     Setup MPI_GET things; 8 is size of real
      call MPI_WIN_CREATE(tmp_r1,nstate*nstate/np,8,MPI_INFO_NULL,
     $     MPI_COMM_WORLD,mpi_win_r,IERR)
      call MPI_WIN_CREATE(tmp_i1,nstate*nstate/np,8,MPI_INFO_NULL,
     $     MPI_COMM_WORLD,mpi_win_i,IERR)

      time=0

      call MPI_BCAST(numqd,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(numsp,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nstate,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(n_in_row,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nEh,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(nqstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(timemax,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(Emin,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(Emax,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(dEh,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(relerr,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(abserr,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(omega_s,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(gamma_s,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(s_mu,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(omega_q_m,numqd,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(gamma_d_m,numqd,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(gamma_q_m,numqd,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(g_couple_m,numqd,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(q_mu_m,numqd,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(pulse_omega_0,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(pulse_t_0,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(pulse_duration,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      call MPI_BCAST(energy0_pls,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
      nstate_p = nstate*nstate/np
      c_offset = nid*nstate/np
      dt = timemax/(nqstep-1)

      do jj=1,nstate/np
         do ii=1,nstate
            kk = jj+c_offset
            if(kk.ge.ii) then
               glo_num_q_Hv(ii+(jj-1)*nstate) = (ii-1)*nstate - 
     $              (ii)*(ii-1)/2 + kk
            else
               glo_num_q_Hv(ii+(jj-1)*nstate) = (kk-1)*nstate - 
     $              (kk-1)*(kk)/2 + ii
            endif
         enddo
      enddo

      call gs_setup(gs_handle_q_H,glo_num_q_Hv,nstate_p,
     $     MPI_COMM_WORLD,np)
 
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum2_parameters    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      common /QRDATA_I1/  i_numqd,i_numsp,i_Eh,i_nqstep
      common /QRDATA_R1/  r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      common /QRDATA_R2/  rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      common /QRDATA_R3/  r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      common /QRDATA_R4/  r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      common /QRDATA_R5/  r_pulse_duration,r_pulse_fluence
      integer             i_numqd,i_numsp,i_Eh,i_nqstep,ii
      real                r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      real                rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      real                r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      real                r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      real                r_pulse_duration,r_pulse_fluence

      real*8 fluence,tmpint
      real   one

      one = 1.0
      PI  = 4.*ATAN(one)

c      call read_quantum_input

c...  temporary parameters
c     timemax = 25e-13
c     Emin    = 1.94                       ! energy range min 
c     Emax    = 2.14                       ! energy range max
c     dEh     = (Emax-Emin)/(nEh+1)     
c     relerr  = 1e-9                      ! relative error  
c     abserr  = 1e-9                      ! absolute error 

c     numqd   = i_numqd
c     numsp   = i_numsp
      numqd   = lnumqd
      numsp   = lnumsp
      nstate  = level**numqd*numsp 
      n_in_row= 2*level**numqd - 1
      if (nid.eq.0) write(6,*) 'parameters reassigned::'
      if (nid.eq.0) write(6,*) 'level = ',level
      if (nid.eq.0) write(6,*) 'numqd = ',numqd
      if (nid.eq.0) write(6,*) 'numsp = ',numsp
      if (nid.eq.0) write(6,*) 'nstate= ',nstate

c...  temporary parameters
      nEh     = i_Eh     
      timemax = r_tmax   
      Emin    = r_Eh1                       ! energy range min 
      Emax    = r_Eh2                      ! energy range max
      dEh     = (Emax-Emin)/(nEh+1)     
      relerr  = er1                       ! relative error  
      abserr  = er2                       ! absolute error 

c...  constants and units
      timeunit= 2.418884326505e-17        ! hbar/Eh(hartree) SI unit 
      timemax = timemax/timeunit          !

c      nqstep  = 20000
      nqstep  = i_nqstep

      dt      = timemax/(nqstep-1)
      time    = 0

c     rintensity0 = 0.00001
      rintensity0 = rintense
      rintensity= rintensity0*1e10        ! W/m**2
      energy0   = 27.450*sqrt(rintensity) ! SI(V/m) 
      energy0   = energy0/5.1421e11       ! a.u.    

      eV        = 27.21140                ! energy Eh  
      debye     = 0.3934303070      

c...  SP: surface plasmon states
      omega_s   = 2.05d0/eV
      gamma_s   = 1.5d-1/eV
      s_mu      = 4.d3*debye

c...  QD: quantum dot
      omega_q   = 2.05d0/eV       
      gamma_q   = (1.9e-7)/eV
      gamma_d   = (2.0e-3)/eV
      q_mu      = 1.3d1*debye  

c... QD: Multiple Quantum Dot Parameters
c    I'm not sure that this is the proper way to do this

      if(numqd.eq.2) then
         omega_q_m(1) = r_omega_q1/eV
         omega_q_m(2) = r_omega_q2/eV
         gamma_q_m(1) = r_gamma_q1/eV
         gamma_q_m(2) = r_gamma_q2/eV
         gamma_d_m(1) = r_gamma_d1/eV
         gamma_d_m(2) = r_gamma_d2/eV
         g_couple_m(1) = r_gcouple_q1/eV
         g_couple_m(2) = r_gcouple_q2/eV
         q_mu_m(1)    = r_q_mu_1*debye
         q_mu_m(2)    = r_q_mu_1*debye
      else
         do ii=1,numqd
            omega_q_m(ii)  = omega_q
            gamma_q_m(ii)  = gamma_q
            gamma_d_m(ii)  = gamma_d
            q_mu_m(ii)     = q_mu
            g_couple_m(ii) = 1d-2/eV 
         enddo
      endif


c...  radius in meters
      qradius   = 20e-9                  ! in meters   
      qradius   = qradius/(0.529177e-10) ! atomic units

c...  coupling constant in a.u.:
      g_couple  = 2.465/qradius**3
      autime    = 2.41888e-17 

c...  pulse constants
      pulse_omega_0  = 2.05d0/eV
c    
      pulse_t_0      = r_pulse_t_0/timeunit !check units
      pulse_duration = r_pulse_duration/timeunit !check units

c     fluence information
c     Transcribed from pulse.f90 by Raman
      c_speed = 2.99792458d10 !cm/s
      d_per_au = 1.0d0 / 3.93430307d-1
      eesu_per_au = 5.14220652d17 / c_speed
      
      eps_med = 2.25
      fluence     = r_pulse_fluence*1d7              ! nJ/cm**2

      tmpint = 5e-1 * sqrt(5d-1 * pi/(2d0*log(2d0)/pulse_duration**2))
     $     * (1e0 + exp(-5e-1 * pulse_omega_0 ** 2 / 
     $     (2d0*log(2d0)/pulse_duration**2)))*timeunit
      energy0_pls =sqrt(4*pi*fluence/(c_speed*sqrt(eps_med)*tmpint))
      energy0_pls = energy0_pls / eesu_per_au

      return
      end


c--------------------------------------------------------------------- 
      subroutine cem_quantum2_propagation_exp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      integer  iEh, ii
      real*8   dclock
      real*8    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      complex  ci,sig

   
      ci = (0.0,1.0)
      
      stime0= dclock()     
      stime1= dclock()
c     Write file headers 
      write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm'
      
      call cem_quantum_parameters
      call cem_quantum_gs
      call cem_quantum2_hamiltonian_setup_pls
      call cem_quantum2_initialize
      

      do iqstep = 1,nqstep
         stime2= dclock()     
         call cem_quantum2_exp(iEh)
         time= time+dt
         call cem_quantum2_fourier_pls
     $        (s_mu_frequency,q_mu_frequency,energy_t)    
         stime2= dclock()-stime2
      enddo
      stime1= dclock()-stime1

      if(nid.eq.0) then
c     Print out results
      write(11,*)'#eV,sig_abs,sig_scatt,sig_tot'
      
      do ii=1,nEh
         rinc = 4.0*pi*(omega(ii)/timeunit)/(c_speed*eps_med**0.5)*
     $        d_per_au*1.0d-18/eesu_per_au
         sig   =rinc*(s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii)
         rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $        d_per_au*1.0d-18/eesu_per_au
         sigr = rinc * cabs((s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii))**2
         sigi  = imag(sig)
         write(11,13) hb_omega(ii),sigi,sigr,sigr+sigi 
      enddo
 13   format(1p5e20.10)
      
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      endif


      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum2_propagation_rk
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'
      integer  iEh, ii, IERR, jj
      real*8     stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      real*8   dclock
      complex  sig
      

      if(nid.eq.0) then
         call cem_quantum2_parameters
c     Write file headers 
      write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm'
      endif

      stime0= 0.0          
      stime1= 0.0          
      stime2= 0.0          
      stime3= 0.0          
      stime4= 0.0      
      stime5= 0.0      
      stime6= 0.0      
      stime7= 0.0      

      stime0= dclock()     
      stime1= dclock()
      call cem_quantum_gs
      call cem_quantum2_hamiltonian_setup_pls
c     call cem_quantum_hamiltonian_setup_pls_dense
      call cem_quantum2_initialize
      call cem_quantum_rk_storage
c     Pulsed case
      do iqstep= 1,nqstep
         stime1= dclock()
         do ii= 1,5
            call rk_quantum_c (ii)
            stime7 = dclock()
            call cem_quantum2_hamiltonian(iEh)
            stime4 = stime4 + dclock()-stime7
            stime7 = dclock()
            call cem_quantum2_linblad_m
            stime5 = stime5 + dclock()-stime7
            call rk_quantum2_ab(ii)
         enddo
         time= time+dt
         stime7= dclock()
         call cem_quantum2_fourier_pls
     $        (s_mu_frequency,q_mu_frequency,energy_t)
         stime7= dclock()-stime7
         stime1= dclock()-stime1

         qrtime1= qrtime1+stime1
         qrtime2= stime4
         qrtime3= stime5
         qrtime4= stime7
         call userchk

      enddo
      stime0 = dclock()-stime0
      
      if(nid.eq.0) then
c     Print out results
      write(6,*)'#eV,sig_abs,sig_scatt,sig_tot'
      
      do ii=1,nEh
         rinc = 4.0*pi*(omega(ii)/timeunit)/(c_speed*eps_med**0.5)*
     $        d_per_au*1.0d-18/eesu_per_au
         sig   =rinc*(s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii)
         rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $        d_per_au*1.0d-18/eesu_per_au
         sigr = rinc * cabs((s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii))**2
         sigi  = imag(sig)
         write(6,13) hb_omega(ii),sigi,sigr,sigr+sigi 
      enddo
 13   format(1p5e20.10)
      
      endif
      return
      end



c--------------------------------------------------------------------- 
      subroutine cem_quantum2_propagation_ode
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      integer  iEh, ii
      real*8  dclock,time_st
      real*8    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      complex  ci,sig

   
      ci = (0.0,1.0)
      
      stime0= dclock()     
      stime1= dclock()
c     Write file headers 
      write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm'
      
      call cem_quantum_parameters
      call cem_quantum_gs
      call cem_quantum2_hamiltonian_setup_pls
      call cem_quantum2_initialize
      
      do iqstep = 1,nqstep
         stime2= dclock() 
         time_st = time
         call cem_quantum2_ode(iEh)
         time= time_st+dt
         call cem_quantum2_fourier_pls
     $        (s_mu_frequency,q_mu_frequency,energy_t)    
         stime2= dclock()-stime2
      enddo
      stime1= dclock()-stime1

      if(nid.eq.0) then
c     Print out results
      write(11,*)'#eV,sig_abs,sig_scatt,sig_tot'
      
      do ii=1,nEh
         rinc = 4.0*pi*(omega(ii)/timeunit)/(c_speed*eps_med**0.5)*
     $        d_per_au*1.0d-18/eesu_per_au
         sig   =rinc*(s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii)
         rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $        d_per_au*1.0d-18/eesu_per_au
         sigr = rinc * cabs((s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii))**2
         sigi  = imag(sig)
         write(11,13) hb_omega(ii),sigi,sigr,sigr+sigi 
      enddo
 13   format(1p5e20.10)
      
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      endif

      return
      end




c-------------------------------------------------------------------- 
      subroutine cem_quantum2_exp(iEh)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      integer  icalld
      save     icalld
      data     icalld/0/
      real*8     glsc2,beta,tmp
      integer  m,info
      integer  iEh,nstate2,n2


      qmarnol= larnol 
      m      = 5 !qmarnol

      nstate2= nstate*nstate/np
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

      call quantum2_arnoldi_vec (qv0,qh0,n2,m,qu,iEh)
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
      subroutine quantum2_arnoldi_vec(v,h,n,m,b,iEh)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      
      integer  i,j,m,n,iEh
      real*8     v(n,m+1),b(n)
      real*8     h(m+1,m)
      real*8     eps, s1, si, glsc2, cnt

      eps = 1.e-13

      s1 = glsc2(b,b,n)
      si = 1./sqrt(s1) ! should add a check and a comment here

      call cmult2(v(1,1),b,si,n)

      do j=1,m
         call quantum2_mult (v(1,j+1),v(1,j),n,iEh)

         do i=1,j
            h(i,j) = glsc2 (v(1,j+1),v(1,i),n)
            cnt    =       -h(i,j  )
            call add2s2    (v(1,j+1),v(1,i),cnt,n)
         enddo

         s1 = glsc2(v(1,j+1),v(1,j+1),n)
         s1 = sqrt(s1)
         h(j+1,j) = s1

c         if (abs(h(j+1,j)).le.eps.and.IFPULSE.eq..false.) then
c             write(6,*) 'h', h(j+1,j)
c             call exitt
c             return   !  should add a comment here
c         else
             si = 1.0/s1
             call cmult(v(1,j+1),si,n)
c         endif

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine quantum2_mult(w,v,n2,iEh)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN' 
      integer  iEh,n2,nstate2
      real*8     w(n2),v(n2)

      nstate2 = nstate*nstate

      call copy(rho_r, v(1+0*nstate2), nstate2)
      call copy(rho_i, v(1+1*nstate2), nstate2)

      call cem_quantum2_hamiltonian(iEh)     
      call cem_quantum2_linblad_m

      call copy(w(1+0*nstate2), resrho_r, nstate2) 
      call copy(w(1+1*nstate2), resrho_i, nstate2) 

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum2_hamiltonian(iEh)
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'

      integer    ii,jj,iEh,ierr,i_one,i_zero
      real*8     tmprr,tmpii,pulseFactor,sqrtns
      real*8     tmprv(nstate*nstate/np),tmpiv(nstate*nstate/np)
      complex*16 tmpfq,tmpcc
      complex*16 ci
      real*8     dclock,ttmp

      i_one = 1
      i_zero = 0
      ci=(0.0,1.0)

c     Ensure Hermiticity of rho
      do jj = 1,nstate/np
      do ii = 1,nstate
         tmprv((jj-1)*nstate+ii) = rho_r(ii,jj)
         tmpiv((jj-1)*nstate+ii) = rho_i(ii,jj)
      enddo
      enddo

      ttmp = dclock()
      call gs_op_fields(gs_handle_q_H,tmprv,
     $     nstate_p,1,1,1,0) 
      call gs_op_fields(gs_handle_q_H,tmpiv,nstate_p,1,1,1,0)
      qctime1= qctime1+dclock()-ttmp     

      do jj=1,nstate/np
      do ii=1,nstate
         if(ii.ne.jj+c_offset) rho_r(ii,jj) = tmprv((jj-1)*nstate+ii)/2
         rho_i(ii,jj) = (2*rho_i(ii,jj) - tmpiv((jj-1)*nstate+ii))/2
         if(jj+c_offset.eq.ii) rho_i(ii,jj) = 0
      enddo
      enddo

c      if(time.gt.7500.and.np.eq.2) then
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.0) call outmat(rho_r,nstate,nstate/2,'rho_r  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.1) call outmat(rho_r,nstate,nstate/2,'rho_r  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.0) call outmat(rho_i,nstate,nstate/2,'rho_i  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.1) call outmat(rho_i,nstate,nstate/2,'rho_i  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call exitt
c      endif

c      if(time.gt.7500.and.np.eq.1) then
c         call outmat(rho_r(1,1),nstate,nstate/2,'rho_r  ',i_zero)
c         call outmat(rho_r(1,6),nstate,nstate/2,'rho_r  ',i_one)
c         call outmat(rho_i(1,1),nstate,nstate/2,'rho_i  ',i_zero)
c         call outmat(rho_i(1,6),nstate,nstate/2,'rho_i  ',i_one)
c         call exitt
c      endif

      
      do jj = 1,nstate/np
      do ii = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo

      do ii=1,(2*level**lnumqd-1)*level**lnumqd
         ham_smu_pls(ii) = ham_smu_min(ii)*pulseFactor(1)
         ham_qmu_pls(ii) = ham_qmu_min(ii)*pulseFactor(1)
      enddo
      call cem_quantum2_mxm(rho_r,tmp_r1, 1)             
      call cem_quantum2_mxm(rho_i,tmp_i1, 1)             
      
c      call mxm(hamilt_I,nstate,rho_r,nstate/np,tmp_r1,nstate)
c      call mxm(hamilt_I,nstate,rho_i,nstate/np,tmp_i1,nstate)


      do jj = 1,nstate/np
      do ii = 1,nstate
         tmprv((jj-1)*nstate+ii) = tmp_r1(ii,jj)
         tmpiv((jj-1)*nstate+ii) = tmp_i1(ii,jj)
      enddo
      enddo

      ttmp= dclock()
      call gs_op_fields(gs_handle_q_H,tmprv,
     $     nstate_p,1,1,1,0)
      call gs_op_fields(gs_handle_q_H,tmpiv,nstate_p,1,1,1,0)
      qctime2= qctime2+dclock()-ttmp     
         

      do jj = 1,nstate/np
      do ii = 1,nstate
         tmpfq= state_sums(ii,jj)*pulse_omega_0*rho(ii,jj)
         tmprr= 2*tmp_r1(ii,jj)-tmprv((jj-1)*nstate+ii)
         tmpii = tmpiv((jj-1)*nstate+ii)
         tmpcc= tmprr + ci*tmpii
         drho(ii,jj)= - ci*(tmpcc - tmpfq) 
         resrho_r(ii,jj)=real(drho(ii,jj) )    
         resrho_i(ii,jj)=imag(drho(ii,jj) )    
      enddo
      enddo

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_quantum2_mxm(rrho,tmprr,idir)             
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer    ii,jj,idir,kk,mm,ns
      real*8       rrho (lstate,lstate/lp),sqrtns
      real*8       tmprr(lstate,lstate/lp),pulseFactor

      do jj = 1,nstate/np
      do ii = 1,nstate
         tmprr(ii,jj)=0
      enddo
      enddo
      
      if(idir.eq.1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
      do jj = 1,nstate/np
      do ii = 1,1
c     Diagonal part
         ns = (ii-1)/ltonq
         tmprr(ii,jj) = tmprr(ii,jj) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(ii,jj)
         sqrtns=sqrt(real(ns) + 1.0)
         do kk=2,n_in_row
            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1
            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii+kk-1,jj)

         enddo

      enddo
      enddo
      do jj = 1,nstate/np
      do ii = 2,n_in_row-1
         ns = (ii-1)/ltonq
c     Diagonal part
         tmprr(ii,jj) = tmprr(ii,jj) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(ii,jj)

         sqrtns=sqrt(real(ns) + 1.0)
   
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=2,n_in_row
            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii+kk-1,jj)

         enddo
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,ii-1

            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii-kk,jj)

         enddo
      enddo
      enddo

c     This is the bulk of the matrix
      do jj = 1,nstate/np
      do ii = n_in_row,nstate-(n_in_row-1)
c     Diagonal part
         ns = (ii-1)/ltonq
         tmprr(ii,jj) = tmprr(ii,jj) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(ii,jj)

         sqrtns=sqrt(real(ns) + 1.0)

c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=2,n_in_row
            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii+kk-1,jj)
c            if(time.gt.5000) then
c               print*,time,(ham_smu_pls(mm)*sqrtns
c     $           +ham_qmu_pls(mm))/pulseFactor(1),mm
c            endif
c            if(time.gt.5500) then
c               print*,'new'
c               call exitt
c            endif
         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,n_in_row-1
            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii-kk,jj)

c            if(time.gt.5000) then
c               print*,time,(ham_smu_pls(mm)*sqrtns
c     $           +ham_qmu_pls(mm))/pulseFactor(1),mm
c            endif
c            if(time.gt.5500) then
c               print*,'new'
c               call exitt
c            endif


         enddo
      enddo
      enddo

c     The end parts only have the up part, not the right part.
      do jj = 1,nstate/np
      do ii = nstate-(n_in_row-1)+1,nstate-1
         ns = (ii-1)/ltonq

c     Diagonal part
         tmprr(ii,jj) = tmprr(ii,jj) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(ii,jj)

         sqrtns=sqrt(real(ns) + 1.0)

c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=2,nstate-(ii-1)
            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii+kk-1,jj)

         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii-kk,jj)

         enddo
      enddo
      enddo

c     Last element has only left of the diagonal terms
      do jj = 1,nstate/np
      do ii = nstate,nstate
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1

            tmprr(ii,jj)=tmprr(ii,jj)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(ii-kk,jj)
         enddo

c     Diagonal part
         tmprr(ii,jj) = tmprr(ii,jj) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + (ii-1)/ltonq*omega_s)*rrho(ii,jj)

      enddo
      enddo

      elseif(idir.eq.-1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
      do jj = 1,nstate/np
      do ii = 1,1
c     Diagonal part
         ns = (ii-1)/ltonq

         tmprr(jj,ii) = tmprr(jj,ii) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(jj,ii)
         sqrtns=sqrt(real(ns) + 1.0)

         do kk=2,n_in_row
            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii+kk-1)

         enddo
         
      enddo
      enddo
      do jj = 1,nstate/np
      do ii = 2,n_in_row-1
c     Diagonal part
         ns = (ii-1)/ltonq

         tmprr(jj,ii) = tmprr(jj,ii) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(jj,ii)
         sqrtns=sqrt(real(ns) + 1.0)

c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=2,n_in_row
            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii+kk-1)

         enddo
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,ii-1

            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1


            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii-kk)
  
         enddo

      enddo
      enddo
      
c     This is the bulk of the matrix
      do jj = 1,nstate/np
      do ii = n_in_row,nstate-(n_in_row-1)
c     Diagonal part
         ns = (ii-1)/ltonq

         tmprr(jj,ii) = tmprr(jj,ii) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(jj,ii)
         sqrtns=sqrt(real(ns) + 1.0)

c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=2,n_in_row

            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii+kk-1)

         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,n_in_row-1

            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1


            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii-kk)

         enddo
      enddo
      enddo

c     The end parts only have the up part, not the right part.
      do jj = 1,nstate/np
      do ii = nstate-(n_in_row-1)+1,nstate-1
c     Diagonal part
         ns = (ii-1)/ltonq

         tmprr(jj,ii) = tmprr(jj,ii) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(jj,ii)
         sqrtns=sqrt(real(ns) + 1.0)

c     Right of the diagonal; good for any row < nstate-n_in_row-1

         do kk=2,nstate-(ii-1)

            mm = mod(n_in_row*(ii-1)+kk-1,n_in_row*ltonq)+1

            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii+kk-1)

         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1

            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1


            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii-kk)

         enddo
      enddo
      enddo

c     Last element has only left of the diagonal terms
      do jj = 1,nstate/np
      do ii = nstate,nstate
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            sqrtns = sqrt(real((ii-kk-1)/ltonq)+1)
            mm = mod(n_in_row*(ii-1)-kk*(n_in_row-1),
     $           n_in_row*ltonq)+1


            tmprr(jj,ii)=tmprr(jj,ii)+(ham_sys_min(mm)*sqrtns+
     $           (ham_smu_pls(mm)*sqrtns
     $           +ham_qmu_pls(mm)))
     $           *rrho(jj,ii-kk)

         enddo
c     Diagonal part

         tmprr(jj,ii) = tmprr(jj,ii) +
     $        (ham_diag(mod((ii-1),ltonq)+1)
     $        + ns*omega_s)*rrho(jj,ii)

      enddo
      enddo
      else
         
         write(6,*) 'no idir: cem_quantum_mxm'
         call exitt
         
      endif

      
      return
      
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum2_linblad_m
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'
      integer  ii,jj,ns,iq,nsp,kns,knsp,kk,ll,mm
      integer  nqtmp1(lnumqd),nqtmp2(lnumqd)
      integer  nq(lnumqd),nqp(lnumqd),req
      integer  nstmp1,nsptmp1,target_core,ierr
      real*8     fac,tmp,tmp3,tmp1,tmp2,tmpr1,tmpi1,tmpr2,tmpi2
      complex*16  rho_passing,ci
      logical  lbFlag,qdFlag,sFlag
      integer status(MPI_STATUS_SIZE)
      integer*8 target_disp !NOT PORTABLE: MPI_ADDRESS_KIND
      real*8   dclock, ttmp
      ci = (0,1.0)
      rho_passing=0
      
      qdFlag=.false.
      sFlag=.false.
      tmpr1 = -3
      tmpi1 = -3
      tmpr2 = -3
      tmpi2 = -3

      

      do jj=1,nstate/np
      do ii=1,nstate
         rho(ii,jj) = rho_r(ii,jj) + ci*rho_i(ii,jj)
         tmp_r1(ii,jj) = real(rho(ii,jj))
         tmp_i1(ii,jj) = imag(rho(ii,jj))
      enddo
      enddo
c
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do jj = 1,nstate/np
         nsp= (jj+c_offset-1)/ltonq
         do iq=1,numqd
            nqp(iq) = nq_minimal(iq,mod(jj+c_offset-1,ltonq)+1)
         enddo
      do ii = 1,nstate
         ns = (ii-1)/ltonq
         tmp3=gamma_s*(ns+nsp)*0.5
         do iq=1,numqd
            nq(iq) = nq_minimal(iq,mod(ii-1,ltonq)+1)  
            tmp1=gamma_d_m(iq)*(nq(iq)**2+nqp(iq)**2-2*nq(iq)*nqp(iq))
            tmp2=gamma_q_m(iq)*(nq(iq)+nqp(iq))*0.5
c     Local Terms:
c...  QD dephasing Linblad term leads to:
            drho(ii,jj)= drho(ii,jj)-tmp1*rho(ii,jj)
c            if(jj+c_offset.gt.nstate/2) print*,nqp(iq),gamma_d_m(iq)
c...  QD spontaneous emission Linblad term leads to these two terms:
            drho(ii,jj)= drho(ii,jj)-tmp2*rho(ii,jj)
            !this may not be correct? I'm not certain.
            !It's certainly not the most efficient way to do this!
c     Nonlocal Terms
c     Get message, calculate term
            if(nqp(iq).eq.0.and.nq(iq).eq.0) then
               if(jj+level**(iq-1).le.nstate/np) then
c              The nonlocal term belongs to me, do the local computation
                  drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
     $              rho(ii+level**(iq-1),jj+level**(iq-1))
               else

c              The nonlocal term belongs to a different core
c              First, find the rank of the other core
                  target_core = (jj+c_offset+level**(iq-1)-1)*np/nstate         
                  target_disp = (ii+level**(iq-1))+mod(jj+level**(iq-1)
     $                 -1,nstate/np)*nstate-1

                  ttmp= dclock() 
                  call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
     $                 mpi_win_r,ierr)
                  call MPI_GET(tmpr2,1,MPI_REAL8
     $                 ,target_core,target_disp,1,MPI_REAL8,mpi_win_r
     $                 ,IERR)
                  call MPI_WIN_UNLOCK(target_core,mpi_win_r,ierr)

                  call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
     $                 mpi_win_i,ierr)
 
                  call MPI_GET(tmpi2,1,MPI_REAL8
     $                 ,target_core,target_disp,1,MPI_REAL8,mpi_win_i
     $                 ,IERR)
                  call MPI_WIN_UNLOCK(target_core,mpi_win_i,ierr)
                  qctime3= qctime3+dclock()-ttmp

                  qdFlag = .true.
               endif
            endif

            if(qdFlag) then
               rho_passing = tmpr2 + ci*tmpi2

               drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*rho_passing
            endif
            qdFlag=.false.
c     Below is for multilevel dots. It is slower than above.
c            kk = mod(ii,level**iq)
c            ll = mod(jj,level**iq)
c            if(kk.lt.level**(iq-1)+level-1.and.kk.gt.0.and.
c     $           ll.lt.level**(iq-1)+level-1.and.ll.gt.0) then
c               drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
c     $              sqrt(real((nq(iq)+1)*(nqp(iq)+1)))*
c     $              rho(ii+level**(iq-1),jj+level**(iq-1))
c            endif
         enddo

c     Receive message
         if(ii+ltonq.lt.nstate.and.
     $        jj+c_offset+ltonq.lt.nstate) then
            if(jj+ltonq.le.nstate/np) then
c              The nonlocal term belongs to me, do the local computation
               ll = ii+ltonq
               kk = jj+ltonq
               drho(ii,jj) =drho(ii,jj)+gamma_s*sqrt(real((ns+1)*(nsp+1)
     $              ))*rho(ii+ltonq,jj+ltonq)            
            else
c           The nonlocal terms belong to someone else
               target_core = (jj+c_offset+ltonq-1)*np/nstate
               target_disp = (ii+ltonq)+mod(jj+ltonq-1,
     $                 nstate/np)*nstate-1
c               print*,target_core,target_disp
               call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
     $                 mpi_win_r,ierr)

               call MPI_GET(tmpr1,1,MPI_REAL8
     $                 ,target_core,target_disp,1,MPI_REAL8,mpi_win_r
     $                 ,IERR)
               call MPI_WIN_UNLOCK(target_core,mpi_win_r,ierr)

               call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
     $                 mpi_win_i,ierr)

               call MPI_GET(tmpi1,1,MPI_REAL8
     $                 ,target_core,target_disp,1,MPI_REAL8,mpi_win_i
     $                 ,IERR)
               call MPI_WIN_UNLOCK(target_core,mpi_win_i,ierr)

               sFlag = .true.
            endif
         endif

         if(sFlag) then
            rho_passing = tmpr1 + ci*tmpi1
            drho(ii,jj)=drho(ii,jj)+gamma_s*sqrt(real((ns+1)*(nsp+1)
     $           ))*rho_passing
         endif
         sFlag=.false.
c...  SP spontaneous emission Linblad terms:
c     More local terms
         drho(ii,jj)= drho(ii,jj)-tmp3*rho(ii,jj)
         resrho_r(ii,jj)= real(drho(ii,jj))
         resrho_i(ii,jj)= imag(drho(ii,jj))
      enddo
      enddo

c      if(np.eq.2) then
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.0) call outmat(resrho_r,nstate,nstate/2,'rho_r  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.1) call outmat(resrho_r,nstate,nstate/2,'rho_r  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.0) call outmat(resrho_i,nstate,nstate/2,'rho_i  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      if(nid.eq.1) call outmat(resrho_i,nstate,nstate/2,'rho_i  ',nid)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c      call exitt
c      endif
c      if(np.eq.1) then
c         call outmat(resrho_r(1,1),nstate,nstate/2,'rho_r  ',nid)
c         call outmat(resrho_r(1,7),nstate,nstate/2,'rho_r  ',nid)
c         call outmat(resrho_i(1,1),nstate,nstate/2,'rho_i  ',nid)
c         call outmat(resrho_i(1,7),nstate,nstate/2,'rho_i  ',nid)
c         call exitt
c      endif
c     Is this barrier necesarry? 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

c--------------------------------------------------------------------
      subroutine cem_quantum2_fourier_pls(s_mu_freq,q_mu_freq,en_t)
c--------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'
      integer  iEh,ii,jj,iq,kk
      integer  nq(lnumqd),nqp,ierr
      complex*16 ci
      complex*16 tmpCMatrix(lstate,lstate),tmpCMatrix2(lstate,lstate)
      real*8     pulseFactor,fullPulse,currentPop,popSymm,popAsymm
      real*8     tmpC1,tmpC2,concurrence
      real*8     popS,popQD(lnumqd)
      real  sp_mu,qd_mu(lnumqd)
      complex  s_mu_freq(lEh),q_mu_freq(lEh),en_t(lEh)
      real     ttmp
      real*8   dclock
 
      ci = (0.0,1.0)

      do jj=1,nstate/np
      do ii=1,nstate
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))
     $        *exp(-ci*state_sums(ii,jj)*pulse_omega_0*(time-pulse_t_0))
      enddo
      enddo

c     Find the expectation values of q_mu and s_mu

      popS=0
      do iq=1,numqd
         popQD(iq) = 0
      enddo
      popSymm=0
      popAsymm=0
c     Find populatipons of quantum dots and plasmon
      do ii=1,nstate/np
         currentPop = real(rho(ii+c_offset,ii))
         popS = popS + (ii+c_offset-1)/ltonq*currentPop

         do iq=1,numqd
            popQD(iq) = popQD(iq) + 
     $           nq_minimal(iq,mod(ii+c_offset-1,ltonq)+1)
     $           *currentPop
         enddo
      enddo

c     Reduce results
      ttmp=dclock()
      if(nid.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,popS,1,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         call MPI_REDUCE(MPI_IN_PLACE,popQD,numqd,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
      else
         call MPI_REDUCE(popS,popS,1,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         call MPI_REDUCE(popQD,popQD,numqd,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
      endif
c    Find populations of symmetric an antisymmetric coupled dot states
c    Traced over plasmon degrees of freedom
c    This is NOT SCALABLE! 
c      if(numqd.eq.2) then
c         do ii=1,numsp
c            popSymm = popSymm + 0.5*real(rho(4*ii-2,4*ii-2)
c     $           + rho(4*ii-2,4*ii-1) + rho(4*ii-1,4*ii-2)
c     $           + rho(4*ii-1,4*ii-1))
c            popAsymm = popAsymm + 0.5*real(rho(4*ii-2,4*ii-2)
c     $           - rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2)
c     $           + rho(4*ii-1,4*ii-1))            
c         enddo
c     Find the concurrence; see equation 19 of He and Zhu
c         ii=1
c         concurrence = 0
c         do ii=1,numsp
c            tmpC1 = (rho(4*ii-1,4*ii-2) + rho(4*ii-2,4*ii-1))**2
c            tmpC2 = 4*imag(rho(4*ii-2,4*ii-2)-rho(4*ii-1,4*ii-1) + 
c     $           rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2))**2
c            concurrence = concurrence + sqrt(tmpC1 + tmpC2)
c         enddo
c      endif


      sp_mu = 0

      do ii=1,nstate/np
c     Calculate local terms
         if(ii+c_offset+ltonq.lt.nstate) then
         sp_mu=sp_mu+(rho(ii+c_offset+ltonq,ii))
     $        *sqrt(real((ii+c_offset+ltonq-1)/ltonq))
         endif
      enddo
c     sp_mu=s_mu*sp_mu
c     Reduce results
      if(nid.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,sp_mu,1,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         sp_mu = 2*s_mu*sp_mu
      else
         call MPI_REDUCE(sp_mu,sp_mu,1,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
      endif



c     Do local terms
c      qd_mu = 0
c      do iq=1,numqd
c         do ii=1,nstate/np,level**iq
c            do jj=0,level**(iq-1)-1
cc     Calculate local terms
c               if(ii+c_offset+level**(iq-1)+jj.lt.nstate) then
c                  qd_mu=qd_mu+(rho(ii+c_offset+level**(iq-1)+jj,ii+jj))
c     $                 *sqrt(real(nq_minimal(iq,mod(ii+jj+c_offset+
c     $                 level**(iq-1)-1,ltonq)+1)))*q_mu_m(iq)
c               endif
c            enddo
c         enddo
c      enddo

c     Will this grab garbage?
c     !! DOES HIGHER NUMQD WRONG

      do iq=1,numqd
         qd_mu(iq) = 0
         do ii=1,nstate/np,level**iq
            do jj=0,level**(iq-1)-1
               qd_mu(iq) = qd_mu(iq)+rho(ii+jj+level**(iq-1)+
     $              c_offset,ii+jj)*q_mu_m(iq)*sqrt(real(nq_minimal(iq,
     $              mod(ii+jj+c_offset+level**(iq-1)-1,ltonq)+1)))
!     $              *sqrt(real(iqstate_m(iq,ii+jj+level**(iq-1))))
            enddo
         enddo
      enddo
c      qd_mu = 2*qd_mu

c     Reduce results
      if(nid.eq.0) then
         call MPI_REDUCE(MPI_IN_PLACE,qd_mu,numqd,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
         do iq=1,numqd
            qd_mu(iq) = 2*qd_mu(iq)
         enddo
      else
         call MPI_REDUCE(qd_mu,qd_mu,numqd,MPI_REAL8,
     $        MPI_SUM,0,MPI_COMM_WORLD,IERR)
      endif

      if(nid.eq.0) then
         write(12,12)time*timeunit,sp_mu,qd_mu(1),qd_mu(2),
     $        qd_mu(1)+qd_mu(2)+sp_mu,fullPulse(1)*eesu_per_au,
     $        popS,popSymm,popAsymm,popQD(1),popQD(2),popAsymm/popSymm,
     $        concurrence
 12      format(1p20g20.8)
      
c     Now, do the fourier transform to save the value
      do ii=1,nEh
         en_t(ii) = en_t(ii)
     $        +fullPulse(1)*exp(ci*omega(ii)*(time-pulse_t_0))
         s_mu_freq(ii) = s_mu_freq(ii)
     $        +(exp(ci*omega(ii)*(time-pulse_t_0))*sp_mu)
         do iq=1,numqd
            q_mu_freq(ii) = q_mu_freq(ii)
     $           +(exp(ci*omega(ii)*(time-pulse_t_0))*qd_mu(iq))
         enddo
      enddo

      endif


      end

c-----------------------------------------------------------------------
      subroutine rk_quantum2_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'RK5'
      include 'mpif.h'
      real*8     ca,cb
      integer  ii,nstate2,ierr

      nstate2= nstate*nstate/np

      ca = rk4a(ii)
      cb = rk4b(ii)

      call rk4_upd(rho_r,k_rho_r,resrho_r,cb,ca,dt,nstate2)
      call rk4_upd(rho_i,k_rho_i,resrho_i,cb,ca,dt,nstate2)

      return
      end

c-----------------------------------------------------------------------
      subroutine f2(time_o,v_ode,w_ode)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN' 
      include 'mpif.h'
      integer  iEh,nstate2,ierr
      real*8   w_ode(2*nstate*nstate/np),v_ode(2*nstate*nstate/np)
      real*8   time_o

      nstate2 = nstate*nstate/np
      time = time_o
      call copy(rho_r, v_ode(1+0*nstate2), nstate2)
      call copy(rho_i, v_ode(1+1*nstate2), nstate2)

      call cem_quantum2_hamiltonian(iEh)     
      call cem_quantum2_linblad_m

      call copy(w_ode(1+0*nstate2), resrho_r, nstate2) 
      call copy(w_ode(1+1*nstate2), resrho_i, nstate2) 

      return
      end

c-------------------------------------------------------------------

      subroutine cem_quantum2_ode
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      external f2
      integer nstate2,iflag!,iwork2(5)
      real*8 v_ode(2*nstate*nstate/np)
c      real*8 work2(100+21*2*100)
      nstate2 = nstate*nstate/np

      iflag=1
      call copy(v_ode(1+0*nstate2), rho_r, nstate2)
      call copy(v_ode(1+1*nstate2), rho_i, nstate2)
      call ode(f2,2*nstate2,v_ode,time,time+dt,relerr,abserr,iflag,work
     $     ,iwork)
c      print*,iflag!,relerr,abserr
      call copy(rho_r,v_ode(1+0*nstate2),nstate2) 
      call copy(rho_i,v_ode(1+1*nstate2),nstate2) 


      return
      end
