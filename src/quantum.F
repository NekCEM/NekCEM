c---------------------------------------------------------------------
      subroutine read_quantum_input
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'QUANTUM'
      include 'mpif.h'

      common /QRDATA_I1/  i_numqd,i_numsp,i_Eh,i_nqstep
      common /QRDATA_R1/  r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      common /QRDATA_R2/  rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      common /QRDATA_R3/  r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      common /QRDATA_R4/  r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      common /QRDATA_R5/  r_pulse_duration,r_pulse_fluence
      integer             i_numqd,i_numsp,i_Eh,i_nqstep
      real               r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      real                rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      real                r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      real                r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      real                r_pulse_duration,r_pulse_fluence
      integer             iend


c...  this routine reads the data in *.box file for quantum casce
c...  the first integer "0" case
     

         call geti1(i_numqd,iend,7)
         call geti1(i_numsp,iend,7)
         
         call getrv(rintense,1,iend,7)
         call getr3(r_Eh1,r_Eh2,r_Eh,iend,7)
         call getr2(r_tmax,r_nqstep,iend,7)
         call getr2(er1,er2,iend,7)
         call getr2(r_omega_q1,r_omega_q2,iend,7)
         call getr2(r_gcouple_q1,r_gcouple_q2,iend,7)
         call getr2(r_gamma_q1,r_gamma_q2,iend,7)
         call getr2(r_gamma_d1,r_gamma_d2,iend,7)
         call getr2(r_q_mu_1,r_q_mu_2,iend,7)
         call getrv(r_pulse_t_0,1,iend,7)
         call getrv(r_pulse_duration,1,iend,7)
         call getrv(r_pulse_fluence,1,iend,7)
         call geti1(IFPULSE,iend,7)
         call geti1(IFDENSE,iend,7)
         call geti1(IFPARALLEL,iend,7)
         call geti1(IFRKQ,iend,7)
         
         i_Eh    =int(r_Eh)
         i_nqstep=int(r_nqstep)

c...  the box file is opened/read by nid=0

         if (nid.eq.0) then
 
         

         write(6,*) 'Read parameter from .box' 
         write(6,*) 'numqd      = ',i_numqd
         write(6,*) 'numsp      = ',i_numsp
         write(6,*) 'rintensity0= ',rintense
         write(6,*) 'Eh1,Eh2,nEh= ',r_Eh1,r_Eh2,i_Eh
         write(6,*) 'tmax,nqstep= ',r_tmax,i_nqstep
         write(6,*) 'er1,er2    = ',er1,er2
         write(6,*) 'omega_q1,omega_q2 = ',r_omega_q1,r_omega_q2
         write(6,*) 'gcouple_q1,gcouple_q2 = ',r_gcouple_q1,
     $        r_gcouple_q2
         write(6,*) 'gamma_q1,gamma_q2 = ',r_gamma_q1,r_gamma_q2
         write(6,*) 'gamma_d1,gamma_d2 = ',r_gamma_d1,r_gamma_d2
         write(6,*) 'q_mu_1,q_mu_2 = ',r_q_mu_1,r_q_mu_2
         write(6,*) 'pulse_t_0 = ',r_pulse_t_0
         write(6,*) 'pulse_duration = ',r_pulse_duration
         write(6,*) 'pulse_fluence = ',r_pulse_fluence
         write(6,*) 'IFPULSE =',IFPULSE,' IFDENSE =',IFDENSE
         write(6,*) 'IFPARALLEL =',IFPARALLEL,' IFRKQ =',IFRKQ
         
      endif
      
c      print*,'pre',ifparallel
c      if(IFPARALLEL) then

c         if(IFRKQ) then
c            call cem_quantum2_propagation_rk
c         else
c            call cem_quantum2_propagation_exp ! exponential time stepping

c         endif
c      else
c         if(IFRKQ) then
c            call cem_quantum_propagation_rk ! runge-kutta timestepping
c         else
c            call cem_quantum_propagation_exp ! exponential time stepping
c         endif
c      endif
c      call exitt
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_parameters    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

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
c      print*,'i_Eh',i_Eh
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

c--------------------------------------------------------------------
      real*8 function fullPulse(n_pls)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer n_pls
      real*8    alpha

      alpha = 2. * log(2.0) / pulse_duration**2
      fullPulse = (energy0_pls*exp(-alpha*
     $     (time-pulse_t_0)**2))
d     $     +energy0_pls*exp(-alpha*(time-1.570796*pulse_t_0)**2))
c     $     +energy0_pls*exp(-alpha*(time-2*pulse_t_0)**2))
     $     *cos(pulse_omega_0*(time-pulse_t_0))

      if(abs(fullPulse).lt.1.0e-75) fullPulse=0
c     $     *cos(pulse_omega_0*(time-pulse_t_0))
c      pulseFactor = 1.
      end

c---------------------------------------------------------------------

      real*8 function pulseFactor(n_pls)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer n_pls
      real*8    alpha
 
      alpha = 2. * log(2.0) / pulse_duration**2
      
      pulseFactor = energy0_pls/2.0*exp(-alpha*(time-pulse_t_0)**2)
c     $ +energy0_pls/2.0*exp(-alpha*(time-1.570796*pulse_t_0)**2)
c     $     +energy0_pls/2.0*exp(-alpha*(time-2*pulse_t_0)**2)
c     $     *cos(pulse_omega_0*(time-pulse_t_0))

      if(abs(pulseFactor).lt.1.0e-75) pulseFactor=0
c     $     *cos(pulse_omega_0*(time-pulse_t_0))
c      pulseFactor = 1e-10
      end

c---------------------------------------------------------------------
      
      integer function binary_rep(bit,number)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer bit,number,currentBit,maxNumber,localNumber,ii
      localNumber = number-1
      maxNumber = level**numqd
      currentBit = 0
      if(number-level**(bit-1).lt.0) then 
         binary_rep=0
         return
      endif
      do ii=1,numqd+1
         currentBit = currentBit+1
         binary_rep = mod(localNumber,level)
         if(currentBit.eq.bit) then
            return
         endif
         localNumber = int(real(localNumber/level))
      enddo
      
      end

c---------------------------------------------------------------------
      subroutine cem_quantum_hamiltonian_setup_m_p
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstate_p,binary_rep
      integer  ieh,nsum,npsum
      real     tmp,tmpi
      real     hamilt_q_mu((2*level**lnumqd-1)*lstate)
      logical  q_muFlag,s_muFlag


c     Assign isstate,iqstate
      kk=0
      do is=1,numsp !number of qunatum states of surface plasmon
         do iq=1,level**numqd !total number of quantum states for all qds
            kk=kk+1
            isstate(kk) = is-1
c           iqstate, since they are all assumed to be the same number of levels,
c           is calculated by taking the (if 2) binary representation of the number
c           and assigning its bit to the corresponding quantum dot
            do ii=1,numqd
               iqstate_m(ii,kk) = binary_rep(ii,iq)
            enddo
         enddo
      enddo

c
c     Initially, generate upper half of s=0 part of the hamiltonian in normal way

c     the number needed to store in a single row is
c     3*2**Nqd - (2**Nqd + 1))
      nstate_p = n_in_row*nstate
      call rzero(hamilt_0_m_p,nstate_p)
      call rzero(hamilt_m_p  ,nstate_p)
      call rzero(hamilt_I_m_p,nstate_p)
      call rzero(hamilt_q_mu,nstate_p)

c     2*level**numberOfQuantumd dots is because we want to all q possibilities
c     for the values of 0 and 1 for s
      do ii=1,level**numqd
         ns=isstate(ii)
         do iq=1,numqd
            nq(iq)=iqstate_m(iq,ii)
         enddo
         do jj=1,2*level**numqd
            do iq=1,numqd
               nqp(iq) = iqstate_m(iq,jj)
            enddo
            s_muFlag=.true.
            nsp = isstate(jj)
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
                  if(q_muFlag) then

c...  q=1, s couples to q'=0, s'=s+1
                     if(nqp(iq).eq.0.and.nq(iq).eq.1.and.(nsp-1).eq.ns)
     $                    hamilt_0_m_p(kk) = hamilt_0_m_p(kk) 
     $                    -g_couple*s_mu*q_mu*dsqrt(dfloat(nsp))
                  
c...  q=0, s couples to q'=1, s'=s-1
                     if(nqp(iq).eq.1.and.nq(iq).eq.0.and.(nsp+1).eq.ns)
     $                    hamilt_0_m_p(kk) = hamilt_0_m_p(kk) 
     $                    -g_couple*s_mu*q_mu*dsqrt(dfloat(nsp+1))
             
                  
c     now do hamilt parts
                     if(nq(iq).eq.0.and.nqp(iq).eq.1.or.nq(iq).eq.1.and.
     $                    nqp(iq).eq.0) then
                        if(ns.eq.nsp) then
c                           if(iq.eq.1) then
                           hamilt_q_mu(kk) =hamilt_q_mu(kk)-q_mu*energy0
     $                          *0.5
c                           endif
                        endif
                     endif
                  endif
                  
                  if(nq(iq).ne.nqp(iq)) then
                     s_muFlag=.false.
                  endif
               enddo
               if(s_muFlag) then
                  if(ns.eq.(nsp-1)) then
                     hamilt_m_p(kk)=hamilt_m_p(kk)-
     $                    dsqrt(dfloat(nsp))*s_mu 
                  endif
                  
                  if(ns.eq.(nsp+1)) then
                     hamilt_m_p(kk)=hamilt_m_p(kk)-
     $                    dsqrt(dfloat(nsp+1))*s_mu
                  endif
                  
c     if(ns.eq.nsp) then
c     hamilt_0(ii,jj)=hamilt_0(ii,jj)+omega_q*nq(iq)
c     $                     + omega_s*ns
c     endif
               endif
c     in RWA, each driving term become e0/2, so full hI = e0/2*hamt+ham0
               hamilt_m_p(kk) = hamilt_0_m_p(kk) + 
     $              hamilt_m_p(kk) * energy0*0.5
            endif
         enddo
      enddo
      
c     Now that we have the upper half of the first part of the matrix, we now propagate it through the whole matrix
c     starting after the last one we assigned, all the way to the last one to assign
      do kk=1,n_in_row*nstate
         ii = mod(kk-1,n_in_row*level**numqd)+1
          !if we are beyond the numsp+1, then we should stop adding
         !contributions from the sp transitions
         if(isstate((kk-1)/n_in_row+1)+1.ge.numsp) then
             hamilt_I_m_p(kk) = hamilt_q_mu(ii)
          else
             hamilt_I_m_p(kk) = hamilt_m_p(ii) *
     $            sqrt(real(isstate((kk-1)/n_in_row+1)+1)) + 
     $            hamilt_q_mu(ii)
          endif
      enddo

c...  set the diagonal terms
      do  ii = 1,nstate
         ns = isstate(ii)
         hamilt_I_m_p(n_in_row*(ii-1)+1) = 
     $            hamilt_I_m_p(n_in_row*(ii-1)+1) + ns*omega_s
c         Add the contribution from each quantum dot
          do kk = 1,numqd
             nq(kk) = iqstate_m(kk,ii) 
             hamilt_I_m_p(n_in_row*(ii-1)+1) = 
     $            hamilt_I_m_p(n_in_row*(ii-1)+1) + nq(kk)*omega_q_m(kk)
          enddo
      enddo


c      open(unit=123,file='matrix')
c      do ii=1,n_in_row*nstate
c         write(123,*)ii,hamilt_0_m_p(ii)
c      enddo
c      call exitt

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
         do jj = 1,nstate
            npsum = isstate(jj)
            do iq=1,numqd
               npsum = npsum + iqstate_m(iq,jj)
            enddo

            do ii = 1,nstate
               nsum = isstate(ii)
               do iq=1,numqd
                  nsum = nsum+iqstate_m(iq,ii)
               enddo
               frequency(ii,jj,ieh)=omega(iEh)*(nsum-npsum)
            enddo
         enddo
      enddo

      end
c---------------------------------------------------------------------
      subroutine cem_quantum_hamiltonian_setup_pls
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstate_p,binary_rep
      integer  iEh,nsum,npsum
      real*8   tmp,tmpi
      real*8   hamilt_q_mu(n_in_row*lstate)
      real*8   hamilt_0_pls(n_in_row*lstate)
      logical  q_muFlag,s_muFlag


c     Assign isstate,iqstate
      kk=0
      do is=1,numsp !number of qunatum states of surface plasmon
         do iq=1,level**numqd !total number of quantum states for all qds
            kk=kk+1
            isstate(kk) = is-1
c           iqstate, since they are all assumed to be the same number of levels,
c           is calculated by taking the (if 2) binary representation of the number
c           and assigning its bit to the corresponding quantum dot
            do ii=1,numqd
               iqstate_m(ii,kk) = binary_rep(ii,iq)
            enddo
         enddo
      enddo

c
c     Initially, generate upper half of s=0 part of the hamiltonian in normal way

c     the number needed to store in a single row is
c     3*2**Nqd - (2**Nqd + 1)

      nstate_p = n_in_row*nstate
      call rzero(hamilt_0_m_p,nstate_p)
      call rzero(hamilt_I_m_p,nstate_p)
      call rzero(hamilt_q_mu,nstate_p)
      call rzero(hamilt_m_p,nstate_p)


c     2*level**numberOfQuantumd dots is because we want to all q possibilities
c     for the values of 0 and 1 for s
      do ii=1,level**numqd
         ns=isstate(ii)
         do iq=1,numqd
            nq(iq)=iqstate_m(iq,ii)
         enddo
         do jj=1,2*level**numqd
            do iq=1,numqd
               nqp(iq) = iqstate_m(iq,jj)
            enddo
            s_muFlag=.true.
            nsp = isstate(jj)
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
     $                    hamilt_0_m_p(kk) = hamilt_0_m_p(kk)
     $                    -g_couple_m(iq)*dsqrt(dfloat((ns+1)*nq(iq)))
                  
c...  q, s couples to q'=q+1, s'=s-1
                     if(nqp(iq).eq.nq(iq)+1.and.nsp.eq.ns-1)
     $                    hamilt_0_m_p(kk) = hamilt_0_m_p(kk)
     $                    -g_couple_m(iq)*dsqrt(dfloat(ns*(nq(iq)+1)))
  
                  
c     now do hamilt parts
                     if(ns.eq.nsp) then
                        if(nqp(iq).eq.nq(iq)+1) then 
                           hamilt_q_mu(kk) = hamilt_q_mu(kk)
     $                          -q_mu_m(iq)*dsqrt(dfloat(nq(iq)+1))
                        endif
                        if(nqp(iq).eq.nq(iq)-1) then
                           hamilt_q_mu(kk) = hamilt_q_mu(kk)
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
                     hamilt_0_pls(kk)=hamilt_0_pls(kk)-
     $                    dsqrt(dfloat(ns+1))*s_mu
                  endif
                  
                  if(nsp.eq.(ns-1)) then
                     hamilt_0_pls(kk)=hamilt_0_pls(kk)-
     $                   dsqrt(dfloat(ns))*s_mu
                  endif
                 
c     if(ns.eq.nsp) then
c     hamilt_0(ii,jj)=hamilt_0(ii,jj)+omega_q*nq(iq)
c     $                     + omega_s*ns
c     endif
               endif
            endif
         enddo
      enddo
      
c     Now that we have the upper half of the first part of the matrix, we now propagate it through the whole matrix
c     starting after the last one we assigned, all the way to the last one to assign
      do kk=1,n_in_row*nstate
         ii = mod(kk-1,n_in_row*level**numqd)+1

         !if we are beyond the numsp+1, then we should stop adding
         !contributions from the sp transitions
         if(isstate((kk-1)/n_in_row+1)+1.ge.numsp) then
             hamilt_m_p(kk) = hamilt_q_mu(ii)

          else
             hamilt_0_m_p(kk) = hamilt_0_m_p(ii)
     $            *sqrt(real(isstate((kk-1)/n_in_row+1)+1))
             hamilt_m_p(kk) = hamilt_0_pls(ii)  
     $            *sqrt(real(isstate((kk-1)/n_in_row+1)+1)) 
     $            + hamilt_q_mu(ii)
c       We dont need extra scaling by sqrt(q), as that was taken care of in the 
c       creation of the elements. We need it for s because we only did the cases
c       of 0 and 1 for s.
          endif
      enddo


c...  set the diagonal terms
      do  ii = 1,nstate
         ns = isstate(ii)
         hamilt_0_m_p(n_in_row*(ii-1)+1) = 
     $            hamilt_0_m_p(n_in_row*(ii-1)+1) + ns*omega_s
         
c         Add the contribution from each quantum dot
          do kk = 1,numqd
             nq(kk) = iqstate_m(kk,ii) 
             hamilt_0_m_p(n_in_row*(ii-1)+1) = 
     $            hamilt_0_m_p(n_in_row*(ii-1)+1) + nq(kk)*omega_q_m(kk)
          enddo
      enddo
      
c      open(unit=123,file='matrix2')
c      do ii=1,n_in_row*nstate
c         write(123,*)ii,hamilt_0_m_p(ii)+hamilt_m_p(ii)
c      enddo
c      call exitt

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
      enddo
      do jj = 1,nstate
         npsum = isstate(jj)
         do iq=1,numqd
            npsum = npsum + iqstate_m(iq,jj)
         enddo
         
         do ii = 1,nstate
            nsum = isstate(ii)
            do iq=1,numqd
               nsum = nsum+iqstate_m(iq,ii)
            enddo
            state_sums(ii,jj)=(nsum-npsum)
         enddo
      enddo
c      call outmat(state_sums,nstate,nstate,'ss  ',nstate)
c      call exitt
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian_setup_pls_dense
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstate2,binary_rep
      integer  ieh,nsum,npsum
      real*8     tmp,tmpi
      logical  q_muFlag,s_muFlag


c........................................... 
c...  iqstate={0 1 0 1 0 1 0 1}
c...  isstate={0 0 1 1 2 2 3 3}
c........................................... 

c     Assign isstate,iqstate
      kk=0
      do is=1,numsp !number of qunatum states of surface plasmon
         do iq=1,level**numqd !total number of quantum states for all qds
            kk=kk+1
            isstate(kk) = is-1
c           iqstate, since they are all assumed to be the same number of levels,
c           is calculated by taking the n-nary representation of the number
c           and assigning its bit to the corresponding quantum dot

            do ii=1,numqd
               iqstate_m(ii,kk) = binary_rep(ii,iq)
            enddo
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
          hamilt_0 (ii,ii) = ns*omega_s
c         Add the contribution from each quantum dot
          do kk = 1,numqd
             nq(kk) = iqstate_m(kk,ii) 
             hamilt_0(ii,ii) = hamilt_0(ii,ii)+nq(kk)*omega_q_m(kk)
          enddo
      enddo

c     loop to do off diagonal terms
c     the nifty 1,3,5 structure isnt as useful in the case of multiple qdots
      do ii=1,nstate
         ns=isstate(ii)
         do iq=1,numqd
            nq(iq)=iqstate_m(iq,ii)
         enddo
         do jj=1,nstate
            !this loop may be redundant with the next
            do iq=1,numqd
               nqp(iq) = iqstate_m(iq,jj)
            enddo
c           Assign the hamiltonian
            s_muFlag=.true.
            nsp = isstate(jj)
            do iq=1,numqd
               q_muFlag=.true.
               do kk=1,numqd
c               We now have to loop over the nq to check if
c               all of the other states are equal
                  if(iq.ne.kk.and.nq(kk).ne.nqp(kk)) then
                     q_muFlag=.false.
                  endif
               enddo
c              now do hamilt parts
               if (q_muFlag) then
                  if(nqp(iq).eq.nq(iq)-1.and.nsp.eq.ns+1)
     $                 hamilt_0(ii,jj) = hamilt_0(ii,jj)
     $                 -g_couple_m(iq)*dsqrt(dfloat((ns+1)*nq(iq)))
  
c...  q, s couples to q'=q+1, s'=s-1
                  if(nqp(iq).eq.nq(iq)+1.and.nsp.eq.ns-1)
     $                 hamilt_0(ii,jj) = hamilt_0(ii,jj)
     $                 -g_couple_m(iq)*dsqrt(dfloat(ns*(nq(iq)+1)))

                  if(ns.eq.nsp) then
                     
                     if(nqp(iq).eq.nq(iq)+1) then 
                        hamilt(ii,jj) = hamilt(ii,jj)
     $                       -dsqrt(dfloat(nq(iq)+1))*q_mu_m(iq)
                      endif
                     if(nqp(iq).eq.nq(iq)-1) then
                        hamilt(ii,jj) = hamilt(ii,jj)
     $                       -dsqrt(dfloat(nq(iq)))*q_mu_m(iq)
                     endif             
                  endif
               endif
               if(nq(iq).ne.nqp(iq)) then
                  s_muFlag=.false.
               endif
            enddo
            
            if(s_muFlag) then
               if(ns.eq.(nsp-1)) then
                  hamilt(ii,jj)=hamilt(ii,jj)-
     $                 dsqrt(dfloat(ns+1))*s_mu 
               endif
               
               if(ns.eq.(nsp+1)) then
                  hamilt(ii,jj)=hamilt(ii,jj)-
     $                 dsqrt(dfloat(ns))*s_mu
               endif
            
c                  if(ns.eq.nsp) then
c                     hamilt_0(ii,jj)=hamilt_0(ii,jj)+omega_q*nq(iq)
c    $                     + omega_s*ns
c                  endif
            endif               
c      in RWA, each driving term become e0/2, so full hI = e0/2*hamt+ham0
         enddo
      enddo

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
         do jj = 1,nstate
            npsum = isstate(jj)
            do iq=1,numqd
               npsum = npsum + iqstate_m(iq,jj)
            enddo

            do ii = 1,nstate
               nsum = isstate(ii)
               do iq=1,numqd
                  nsum = nsum+iqstate_m(iq,ii)
               enddo
               state_sums(ii,jj)=(nsum-npsum)
            enddo
         enddo
      enddo

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian_setup_m
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstate2,binary_rep
      integer  ieh,nsum,npsum
      real*8     tmp,tmpi
      logical  q_muFlag,s_muFlag


c........................................... 
c...  iqstate={0 1 0 1 0 1 0 1}
c...  isstate={0 0 1 1 2 2 3 3}
c........................................... 

c     Assign isstate,iqstate
      kk=0
      do is=1,numsp !number of qunatum states of surface plasmon
         do iq=1,level**numqd !total number of quantum states for all qds
            kk=kk+1
            isstate(kk) = is-1
c           iqstate, since they are all assumed to be the same number of levels,
c           is calculated by taking the n-nary representation of the number
c           and assigning its bit to the corresponding quantum dot

            do ii=1,numqd
               iqstate_m(ii,kk) = binary_rep(ii,iq)
            enddo
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
          hamilt_0 (ii,ii) = ns*omega_s
c         Add the contribution from each quantum dot
          do kk = 1,numqd
             nq(kk) = iqstate_m(kk,ii) 
             hamilt_0 (ii,ii) = hamilt_0(ii,ii)+nq(kk)*omega_q
          enddo
      enddo

c     loop to do off diagonal terms
c     the nifty 1,3,5 structure isn't as useful in the case of multiple qdots
      do ii=1,nstate
         ns=isstate(ii)
         do iq=1,numqd
            nq(iq)=iqstate_m(iq,ii)
         enddo
         do jj=1,nstate
            !this loop may be redundant with the next
            do iq=1,numqd
               nqp(iq) = iqstate_m(iq,jj)
            enddo
c           Assign the hamiltonian
            s_muFlag=.true.
            nsp = isstate(jj)
            do iq=1,numqd
               q_muFlag=.true.
               do kk=1,numqd
c               We now have to loop over the nq to check if
c               all of the other states are equal
                  if(iq.ne.kk.and.nq(kk).ne.nqp(kk)) then
                     q_muFlag=.false.
                  endif
               enddo
c              now do hamilt parts
               if (q_muFlag) then

c...           q=1, s couples to q'=0, s'=s+1
                  if(nqp(iq).eq.0.and.nq(iq).eq.1
     $                 .and.(nsp-1).eq.ns) then
                     hamilt_0 (ii,jj) = hamilt_0(ii,jj) 
     $                    -g_couple*s_mu*q_mu*dsqrt(dfloat(nsp))
                  endif
c...  q=0, s couples to q'=1, s'=s-1
                  if(nqp(iq).eq.1.and.nq(iq).eq.0
     $                 .and.(nsp+1).eq.ns) then
                     hamilt_0 (ii,jj) = hamilt_0(ii,jj) 
     $                    -g_couple*s_mu*q_mu*dsqrt(dfloat(nsp+1))
                  endif
                  
                  if(nq(iq).eq.0.and.nqp(iq).eq.1.or.nq(iq).eq.1.and.
     $                 nqp(iq).eq.0) then
                     if(ns.eq.nsp) then
                        hamilt(ii,jj) = hamilt(ii,jj)-q_mu
                     endif
                  endif
               endif
               if(nq(iq).ne.nqp(iq)) then
                  s_muFlag=.false.
               endif
            enddo
            
            if(s_muFlag) then
               if(ns.eq.(nsp-1)) then
                  hamilt(ii,jj)=hamilt(ii,jj)-
     $                 dsqrt(dfloat(nsp))*s_mu 
               endif
               
               if(ns.eq.(nsp+1)) then
                  hamilt(ii,jj)=hamilt(ii,jj)-
     $                 dsqrt(dfloat(nsp+1))*s_mu
               endif
            
c                  if(ns.eq.nsp) then
c                     hamilt_0(ii,jj)=hamilt_0(ii,jj)+omega_q*nq(iq)
c    $                     + omega_s*ns
c                  endif
            endif               
c      in RWA, each driving term become e0/2, so full hI = e0/2*hamt+ham0
            hamilt_I(ii,jj) = hamilt_0(ii,jj)+hamilt(ii,jj)*energy0*0.5
            hamilt(ii,jj) = energy0*0.5*hamilt(ii,jj)
         enddo
      enddo


c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
         do jj = 1,nstate
            npsum = isstate(jj)
            do iq=1,numqd
               npsum = npsum + iqstate_m(iq,jj)
            enddo

            do ii = 1,nstate
               nsum = isstate(ii)
               do iq=1,numqd
                  nsum = nsum+iqstate_m(iq,ii)
               enddo
               frequency(ii,jj,ieh)=omega(iEh)*(nsum-npsum)
            enddo
         enddo
      enddo

c      call outmat(hamilt_0,nstate,nstate,'hamilt_0  ',nstate)
c      call outmat(hamilt  ,nstate,nstate,'hamilt    ',nstate)
c      call outmat(hamilt_I,nstate,nstate,'hamilt_I  ',nstate)

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_linblad_m
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  ii,jj,ns,iq,nsp,kns,knsp,kk,ll,mm
      integer  nqtmp1(lnumqd),nqtmp2(lnumqd)
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstmp1,nsptmp1
      real*8     fac,tmp,tmp3,tmp1,tmp2
      logical  lbFlag


      do jj = 1,nstate
         nsp= isstate(jj)
         do iq=1,numqd
            nqp(iq) = iqstate_m(iq,jj)
         enddo
      do ii = 1,nstate
         ns = isstate(ii)
         tmp3=gamma_s*(ns+nsp)*0.5
         do iq=1,numqd
            nq(iq) = iqstate_m(iq,ii)         
            tmp1=gamma_d_m(iq)*(nq(iq)**2+nqp(iq)**2-2*nq(iq)*nqp(iq))
c            tmp2=gamma_q_m(iq)*(nq(iq)+nqp(iq))*0.5
c...  QD dephasing Linblad term leads to:
            drho(ii,jj)= drho(ii,jj)-tmp1*rho(ii,jj)
c...  QD spontaneous emission Linblad term leads to these two terms:
            drho(ii,jj)= drho(ii,jj)-tmp2*rho(ii,jj)
            !this may not be correct? I'm not certain.
            !It's certainly not the most efficient way to do this!
           
            if(nqp(iq).eq.0.and.nq(iq).eq.0) then
               drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
     $              rho(ii+level**(iq-1),jj+level**(iq-1))
c              print*,time,ii,jj,rho(ii+level**(iq-1),jj+level**(iq-1))
c     $              ,'1'
c               if(time.gt.1000) call exitt
c               print*,'1'
c               print*,'ij',iqstate_m(iq,ii),iqstate_m(iq,jj)
c               print*,'kl',iqstate_m(iq,ll),iqstate_m(iq,kk)
            endif
            
c     Below is for multlevel qds. It is slower than above.

c            kk = mod(ii,level**iq)
c            ll = mod(jj,level**iq)
c            if(kk.lt.level**(iq-1)+level-1.and.kk.gt.0.and.
c     $           ll.lt.level**(iq-1)+level-1.and.ll.gt.0) then
c               drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
c     $              sqrt(real((nq(iq)+1)*(nqp(iq)+1)))*
c     $              rho(ii+level**(iq-1),jj+level**(iq-1))
c            endif
         enddo
         


         if(ii+level**numqd.lt.nstate.and.
     $        jj+level**numqd.lt.nstate) then
            ll = ii+level**numqd
            kk = jj+level**numqd
            drho(ii,jj) = drho(ii,jj) + gamma_s*sqrt(real((ns+1)*(nsp+1)
     $           ))*rho(ii+level**numqd,jj+
     $           level**numqd) 
c            print*,time,ii,jj,rho(ll,kk),'2'

c            print*,rho(ii+level**numqd,jj+
c     $           level**numqd),'rhopass'
         endif
         
c...  SP spontaneous emission Linblad terms:
         drho(ii,jj)= drho(ii,jj)-tmp3*rho(ii,jj)
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
      integer  iEh, iqstep,ii
      real*8  dclock
      real*8    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      complex  ci,sig

     
      ci = (0.0,1.0)

      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup
      stime0= dclock()
      
      
      stime1= dclock()
      if(IFPULSE) then

         call cem_quantum_initialize
c     Write file headers 
         write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm,concurrence,popSymm2,
     $popAsymm2'
      
         do iqstep = 1,nqstep

            stime2= dclock()            
            call cem_quantum_exp(iEh)
            time= time+dt
            call cem_quantum_fourier_pls
     $           (s_mu_frequency,q_mu_frequency,energy_t)    
            stime2= dclock()-stime2
         enddo
         stime1= dclock()-stime1
c      print the results out
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
         tmpr  = omega(ii)*real(s_mu_frequency(ii)+q_mu_frequency(ii))
         tmpi  = omega(ii)*imag(s_mu_frequency(ii)+q_mu_frequency(ii))
         write(11,13) hb_omega(ii),sigi,sigr,sigr+sigi
      enddo
        
  13   format(1p5e20.10)
      else
         do iEh = 1,nEh
            call cem_quantum_initialize
            
            stime1= dclock()
            do iqstep = 1,nqstep
               stime2= dclock()
               call cem_quantum_exp(iEh)
               time= time+dt
               call cem_quantum_energy_m(iEh)     
               stime2= dclock()-stime2
            enddo
            stime1= dclock()-stime1
            call cem_quantum_post_m(iEh)            
         enddo
      endif
         
      
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      return
      end



c--------------------------------------------------------------------- 
      subroutine cem_quantum_propagation_ode
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh, iqstep,ii
      real*8  dclock,time_st
      real*8    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      complex  ci,sig

     
      ci = (0.0,1.0)

      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup
      stime0= dclock()
      
      
      stime1= dclock()
      if(IFPULSE) then

         call cem_quantum_initialize
c     Write file headers 
c         write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
c     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
c     $popQD1,popQD2,popAsymm/popSymm,concurrence,popSymm2,
c     $popAsymm2'
      
         do iqstep = 1,nqstep

            stime2= dclock()  
            time_st = time
            call cem_quantum_ode(iEh)
c            print*,time,time_st
            time= time_st+dt
            call cem_quantum_fourier_pls
     $           (s_mu_frequency,q_mu_frequency,energy_t)    
            stime2= dclock()-stime2
         enddo
         stime1= dclock()-stime1
c      print the results out
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
         tmpr  = omega(ii)*real(s_mu_frequency(ii)+q_mu_frequency(ii))
         tmpi  = omega(ii)*imag(s_mu_frequency(ii)+q_mu_frequency(ii))
         write(11,13) hb_omega(ii),sigi,sigr,sigr+sigi 
      enddo
        
  13   format(1p5e20.10)
      else
         do iEh = 1,nEh
            call cem_quantum_initialize
            
            stime1= dclock()
            do iqstep = 1,nqstep
               stime2= dclock()
               call cem_quantum_exp(iEh)
               time= time+dt
               call cem_quantum_energy_m(iEh)     
               stime2= dclock()-stime2
            enddo
            stime1= dclock()-stime1
            call cem_quantum_post_m(iEh)            
         enddo
      endif
         
      
      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      return
      end

c-------------------------------------------------------------------

      subroutine cem_quantum_ode
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      external f
      integer nstate2,iflag,iwork(5)
      real*8 v_ode(2*nstate*nstate),work(100+21*2*nstate*nstate)
c      SAVE work, iwork
      nstate2 = nstate*nstate

      iflag=1
      call copy(v_ode(1+0*nstate2), rho_r, nstate2)
      call copy(v_ode(1+1*nstate2), rho_i, nstate2)
      call ode(f,2*nstate2,v_ode,time,time+dt,relerr,abserr,iflag,work,
     $     iwork)
      print*,iflag
      call copy(rho_r,v_ode(1+0*nstate2),nstate2) 
      call copy(rho_i,v_ode(1+1*nstate2),nstate2) 




      return
      end
c--------------------------------------------------------------------
      subroutine cem_quantum_hamiltonian_setup
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      if(IFPULSE) then
         if(IFDENSE) then
            call cem_quantum_hamiltonian_setup_pls_dense
         else
            call cem_quantum_hamiltonian_setup_pls
         endif
      else if(IFDENSE) then
         call cem_quantum_hamiltonian_setup_m
      else
         call cem_quantum_hamiltonian_setup_m_p
      endif

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_propagation_rk
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh, iqstep, ii,jj
      real*8     stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      real*8   dclock
      complex  sig
      
      call cem_quantum_parameters
      call cem_quantum_hamiltonian_setup
      call cem_quantum_rk_storage
 
      stime0=dclock()
      if(IFPULSE) then
         call cem_quantum_initialize
c     Write file headers 
         write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm'
         write(11,*)'#eV,sig_abs,sig_scatt,sig_tot'

         print*,' '
         print*,' '
         print*,' '
         print*,' '
         print*,' '
c     Pulsed case
         do iqstep= 1,nqstep
            stime2= dclock()
            do ii= 1,5
               call rk_quantum_c (ii)
               stime7 = dclock()
               call cem_quantum_hamiltonian(iEh)
               stime4 = stime4 + dclock()-stime7
               stime7 = dclock()
               call cem_quantum_linblad_m
               stime5 = stime5 + dclock()-stime7
               call rk_quantum_ab(ii)
c               if(iqstep.gt.2) call exitt
c               print*,'rk',ii
            enddo

            do jj=1,nstate
               do ii=1,nstate
                  rho_r(ii,jj) = rho_r(ii,jj)!*(1e-1)
                  rho_i(ii,jj) = rho_i(ii,jj)!*(1e-1)
               enddo
            enddo
c            print*,time,sig
c            print*,time,rho_r(1,11),rho_i(2,11)
            time= time+dt
            call cem_quantum_fourier_pls
     $           (s_mu_frequency,q_mu_frequency,energy_t)
            stime2= dclock()-stime2
         enddo
         
         stime1= dclock()-stime1
         stime3= stime2/neqn
c         open(unit=11,file='spec.out')
c     Print out results
         print*,'after write 11',IFPULSE
         do ii=1,nEh
            rinc = 4.0*pi*(omega(ii)/timeunit)/(c_speed*eps_med**0.5)*
     $           d_per_au*1.0d-18/eesu_per_au
            sig   =rinc*(s_mu_frequency(ii)+q_mu_frequency(ii))/
     $           energy_t(ii)
            rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $           d_per_au*1.0d-18/eesu_per_au
            sigr = rinc * cabs((s_mu_frequency(ii)+q_mu_frequency(ii))/
     $           energy_t(ii))**2
            sigi  = imag(sig)
            write(11,13) hb_omega(ii),sigi,sigr,sigr+sigi 
         enddo
 13      format(1p5e20.10)
         
      else 
         print*,'not pulse'
         print*,nEh
c     Continuous wave
         do iEh = 1,nEh
c            print*,iEh,nEh
            call cem_quantum_initialize
            stime1= dclock()
            time=0
            do iqstep= 1,nqstep
               stime2= dclock()
               do ii= 1,5
                  call rk_quantum_c (ii)
                  call cem_quantum_hamiltonian(iEh) 
                  call cem_quantum_linblad_m
                  call rk_quantum_ab(ii)
               enddo
c               print*,time
               time= time+dt
               call cem_quantum_fourier_pls(iEh)     
               stime2= dclock()-stime2
            enddo
            write(12,*),' '
            stime1= dclock()-stime1
            call cem_quantum_post_m(iEh)            
         enddo
      endif
      
      stime0 = dclock()-stime0
      write(6,6) stime0,stime1,stime2,stime3 
 6    format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)

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
c      writE(6,6) dt, time, RKtime
c  6   format('dt/time/RKtime:: ',1p3e20.9)

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
      real*8     ca,cb
      integer  ii,nstate2,i_zero,i_one
      i_zero = 0
      i_one = 1

      nstate2= nstate*nstate/np

      ca = rk4a(ii)
      cb = rk4b(ii)
c      if(time.gt.1500.and.time.lt.3000) then
c      call outmat(resrho_r(1,1),nstate,nstate/2,'res_r  ',i_zero)
c      call outmat(resrho_r(1,4),nstate,nstate/2,'res_r  ',i_one)
c      call outmat(resrho_i(1,1),nstate,nstate/2,'res_i  ',i_zero)
c      call outmat(resrho_i(1,4),nstate,nstate/2,'res_i  ',i_one)
c      call outmat(rho_r(1,1),nstate,nstate/2,'rho_r  ',i_zero)
c      call outmat(rho_r(1,4),nstate,nstate/2,'rho_r  ',i_one)
c      call outmat(rho_i(1,1),nstate,nstate/2,'rho_i  ',i_zero)
c      call outmat(rho_i(1,4),nstate,nstate/2,'rho_i  ',i_one)
c      endif

      call rk4_upd(rho_r,k_rho_r,resrho_r,cb,ca,dt,nstate2)
      call rk4_upd(rho_i,k_rho_i,resrho_i,cb,ca,dt,nstate2)

c      if(time.gt.1500.and.time.lt.3000) then
c      call outmat(rho_r(1,1),nstate,nstate/2,'ares_r  ',i_zero)
c      call outmat(rho_r(1,4),nstate,nstate/2,'ares_r  ',i_one)
c      call outmat(rho_i(1,1),nstate,nstate/2,'ares_i  ',i_zero)
c      call outmat(rho_i(1,4),nstate,nstate/2,'ares_i  ',i_one)
c      endif
c      if(time.gt.3000) call exitt
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_initialize    
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  nstate2,ii,jj

      nstate2= nstate*nstate

      call czero(s_mu_frequency,nEh)
      call czero(q_mu_frequency,nEh)
      call czero(energy_t,      nEh)
      call rzero(   rho_r,nstate2)
      call rzero(   rho_i,nstate2)
      call rzero( k_rho_r,nstate2)
      call rzero( k_rho_i,nstate2)
      call rzero(resrho_r,nstate2)
      call rzero(resrho_i,nstate2)

c      do jj=1,nstate
c         do ii=1,nstate
c            rho_r(ii,jj) = 1
c            if(ii.gt.jj) then
c               rho_i(ii,jj) = -1
c            else if(ii.eq.jj) then
c               rho_i(ii,jj) = 0
c            else
c               rho_i(ii,jj) = 1
c            endif
c         enddo
c      enddo
      rho_r(1,1) = 1
      return
      end

c--------------------------------------------------------------------
      subroutine cem_quantum_fourier_pls(s_mu_freq,q_mu_freq,en_t)
c--------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj,iq,kk
      integer  nq(lnumqd),nqp,info
      complex*16 ci
      complex*16 tmpCMatrix(lstate,lstate),tmpCMatrix2(lstate,lstate)
      real*8     pulseFactor,fullPulse,currentPop,popSymm,popAsymm
      real*8     tmpC1,tmpC2,concurrence,popSymm2,popAsymm2,con2
      real*8     popS,popQD(lnumqd),rho_workr(8)
      real       qd_mu(lnumqd),sp_mu,signRho
      complex    s_mu_freq(lEh),q_mu_freq(lEh),en_t(lEh),rhoOrdering
      complex    red_rho(4,4),rho_svd(4,4),rho_work(8),rhoVL(4,4)
      complex    rhoVR(4,4),rhoEval(4),spinFlip(4,4),rhoFinal(4,4)
 
      ci = (0.0,1.0)

      do jj=1,nstate
      do ii=1,nstate
         if(IFPULSE) then
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))
     $        *exp(-ci*state_sums(ii,jj)*pulse_omega_0*(time-pulse_t_0))
         else
            rho(ii,jj) = (rho_r(ii,jj)+ci*rho_i(ii,jj))
c     $        *exp(-ci*frequency(ii,jj,iEh)*time)
         endif
      enddo
      enddo
     
c     Find the expectation values of q_mu and s_mu

      popS=0
      do iq=1,numqd
         popQD(iq) = 0
      enddo
      popSymm   = 0
      popAsymm  = 0
      popSymm2  = 0
      popAsymm2 = 0
c     Find populatipons of quantum dots and plasmon
      do ii=1,nstate
         currentPop = real(rho(ii,ii))
         popS = popS + isstate(ii)*currentPop
         do iq=1,numqd
            popQD(iq) = popQD(iq)+ iqstate_m(iq,ii)*currentPop
         enddo
      enddo
c    Find populations of symmetric an antisymmetric coupled dot states
c    Traced over plasmon degrees of freedom
c    This is NOT SCALABLE! 
      if(numqd.eq.2) then
         do ii=1,numsp
            popSymm = popSymm + 0.5*real(rho(4*ii-2,4*ii-2)
     $           + rho(4*ii-2,4*ii-1) + rho(4*ii-1,4*ii-2)
     $           + rho(4*ii-1,4*ii-1))
            popAsymm = popAsymm + 0.5*real(rho(4*ii-2,4*ii-2)
     $           - rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2)
     $           + rho(4*ii-1,4*ii-1))  
            popSymm2 = popSymm2 + 0.5*real(rho(4*ii-3,4*ii-3)
     $           + rho(4*ii-3,4*ii) + rho(4*ii,4*ii-3) 
     $           + rho(4*ii,4*ii))
            popAsymm2 = popAsymm2 + 0.5*real(rho(4*ii-3,4*ii-3)
     $           - rho(4*ii-3,4*ii) - rho(4*ii,4*ii-3) 
     $           + rho(4*ii,4*ii))
         enddo
c     Find the concurrence; see equation 19 of He and Zhu
         call czero(red_rho,16)
         do ii=1,numsp
            do jj=1,4
               do kk=1,4
                  red_rho(kk,jj) = red_rho(kk,jj) + 
     $                 rho(4*ii-(4-kk),4*ii-(4-jj))
               enddo
            enddo
         enddo
         
         popSymm2 = 0.5*(red_rho(2,2) + red_rho(3,3) + 
     $        red_rho(2,3) + red_rho(3,2))
         popAsymm2 = 0.5*(red_rho(2,2) + red_rho(3,3) - 
     $        red_rho(2,3) - red_rho(3,2))
c         red_rho(2,2) = 0.5
c         red_rho(3,2) = 0.5
c         red_rho(2,3) = 0.5
c         red_rho(3,3) = 0.5
         
         call czero(spinFlip,16)
         do ii=1,4
            do jj=1,4
               rhoVL(ii,jj) = conjg(red_rho(ii,jj))
            enddo
         enddo
c         spinFlip(1,2) = -ci
c         spinFlip(2,1) = ci
c         spinFlip(3,4) = -ci
c         spinFlip(4,3) = ci
          spinFlip(1,4) = -1
          spinFlip(4,1) = -1
          spinFlip(2,3) = 1
          spinFlip(3,2) = 1
c         print*,'spinFlip'
c         call outmat(spinFlip,4,4,'spin ',nid)
c         do ii=1,4
c            write(*,"(100g15.5)") ( spinFlip(ii,jj), jj=1,4 )
c         enddo
         call zgemm('N','N',4,4,4,(1.0,0.0),
     $        spinFlip,4,rhoVL,4,(0.0,0.0),rhoVR,4)
c         call outmat(rho_svd,4,4,'svd  ',nid)
         
         call zgemm('N','N',4,4,4,(1.0,0.0),rhoVR,4,
     $        spinFlip,4,(0.0,0.0),rho_svd,4)
c         print*,'rho_svd'
c         do ii=1,4
c            write(*,"(100g15.5)") ( rho_svd(ii,jj), jj=1,4 )
c         enddo

c         call outmat(red_rho,4,4,'red  ',nid)
         call zgemm('N','N',4,4,4,(1.0,0.0),red_rho,4,rho_svd,
     $        4,(0.0,0.0),rhoFinal,4)
c         print*,'rhoFinal'
c         do ii=1,4
c            write(*,"(100g15.5)") ( rhoFinal(ii,jj), jj=1,4 )
c         enddo

c         call czero(rhoEval,4)
c         print*,'before zgeev'
         call zgeev('N','N',4,rhoFinal,4,rhoEval,rhoVL,4,rhoVR,4,
     $        rho_work,8,rho_workR,info)
c         print*,'after zgeev'
c         call exitt
         do ii=2,4
            rhoOrdering = rhoEval(ii)
            jj=ii-1
            do while(jj.ge.1.and.real(rhoEval(jj)).gt.real(rhoOrdering))
               rhoEval(jj+1) = rhoEval(jj)
               jj = jj - 1
            enddo
            rhoEval(jj+1) = rhoOrdering
         enddo

         do ii=1,4
            if(abs(real(rhoEval(ii))).lt.1e-14) rhoEval(ii)=0
         enddo
c        if(time.gt.6000) then
c           print*,'eval ',info
c           print*,rhoEval
c           print*,'red_rho'
c           do ii=1,4
c              write(*,"(100g15.5)") ( red_rho(ii,jj), jj=1,4 )
c           enddo
c           call exitt
c         endif
         concurrence = sqrt(real(rhoEval(4))) - real(sqrt(rhoEval(3))) - 
     $        sqrt(real(rhoEval(2))) - sqrt(real(rhoEval(1)))
         con2=0
         do ii=1,numsp
            tmpC1 = (rho(4*ii-1,4*ii-2) + rho(4*ii-2,4*ii-1))**2
            tmpC2 = 4*imag(rho(4*ii-2,4*ii-2)-rho(4*ii-1,4*ii-1) + 
     $           rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2))**2
            con2 = con2 + sqrt(tmpC1 + tmpC2)
         enddo
      endif
      

      sp_mu = 0
      do ii=1,nstate-level**numqd
         sp_mu=sp_mu+(rho(ii,ii+level**numqd)
     $        +rho(ii+level**numqd,ii))
     $        *sqrt(real(isstate(ii+level**numqd)))
      enddo
      sp_mu=s_mu*sp_mu

      do iq=1,numqd
         qd_mu(iq) = 0
         do ii=1,nstate,level**iq
            do jj=0,level**(iq-1)-1
               qd_mu(iq) = qd_mu(iq) + (rho(ii+jj,ii+jj+level**(iq-1))
     $              + rho(ii+jj+level**(iq-1),ii+jj))*q_mu_m(iq)
     $              *sqrt(real(iqstate_m(iq,ii+jj+level**(iq-1))))
            enddo
         enddo
      enddo

      write(12,12)time*timeunit,sp_mu,qd_mu(1),qd_mu(2),
     $     qd_mu(1)+qd_mu(2)+sp_mu,fullPulse(1)*eesu_per_au,
     $     popS,popSymm,popAsymm,popQD(1),popQD(2),popAsymm/popSymm,
     $     concurrence,con2,popSymm2,popAsymm2
 12   format(1p20g20.8)
c      print*,'after write'
c     Now, do the fourier transform to save the value
      if(IFPULSE) then
      do ii=1,nEh
         en_t(ii) = en_t(ii)
     $        +fullPulse(1)*exp(ci*omega(ii)*(time-pulse_t_0))
         s_mu_freq(ii) = s_mu_freq(ii)
     $        +(exp(ci*omega(ii)*(time-pulse_t_0))*sp_mu)
!     $        *exp(-ci*state_sums(jj,kk)*pulse_omega_0*time)
         do iq=1,numqd
            q_mu_freq(ii) = q_mu_freq(ii)
     $           +(exp(ci*omega(ii)*(time-pulse_t_0))*qd_mu(iq))
         enddo
!     $        *exp(-ci*state_sums(jj,kk)*pulse_omega_0*time)
C         write(14,12)time,pulseFactor(1),
c     $        exp(ci*omega(ii)*(time-pulse_t_0))

      enddo
      endif

      return
      end



c--------------------------------------------------------------------- 
      subroutine cem_quantum_energy_m (iEh)       
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj,iq    
      integer  nq(lnumqd),nqp
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
      do ii=1,nstate
         do iq=1,numqd
            nq(iq)=iqstate_m(iq,ii)
         enddo
         do jj=1,nstate
            do iq=1,numqd
               nqp=iqstate_m(iq,jj)
               if(nq(iq).eq.0.and.nqp.eq.1) then 
                  cq_mu=cq_mu+imag(rho(ii,jj))
               endif
            enddo
         enddo
      enddo
      cq_mu=cq_mu*q_mu/0.3934303070 

      cs_mu=0.0
      do ii=1,nstate-level**numqd
         cs_mu=cs_mu+rho(ii,ii+level**numqd)+
     $           rho(ii+level**numqd,ii) 
         
      enddo 
      cs_mu=cs_mu*s_mu/0.3934303070 

c...  energies
c     Not correct
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
      subroutine cem_quantum_post_m(iEh)               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer  iEh,ii,jj,iq,kk
      integer  nq(lnumqd),nqp,info
      complex*16 ci,ensp
      real*8     rinc,sigcon,sigi,sigr
      complex*16 tmpCMatrix(lstate,lstate),tmpCMatrix2(lstate,lstate)
      real*8     pulseFactor,fullPulse,currentPop,popSymm,popAsymm
      real*8     tmpC1,tmpC2,concurrence,popSymm2,popAsymm2,con2
      real*8     popS,popQD(lnumqd),rho_workr(8)
      real       qd_mu,sp_mu,signRho,sig
      complex    s_mu_freq(lEh),q_mu_freq(lEh),en_t(lEh),rhoOrdering
      complex    red_rho(4,4),rho_svd(4,4),rho_work(8),rhoVL(4,4)
      complex    rhoVR(4,4),rhoEval(4),spinFlip(4,4),rhoFinal(4,4)
 

      ci = (0.0,1.0)

      
      do jj=1,nstate
      do ii=1,nstate
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))
      enddo
      enddo

      sp_mu = (0,0)
      do ii=1,nstate-level**numqd
         sp_mu=sp_mu+dimag(rho(ii,ii+level**numqd))
     $        *sqrt(real(isstate(ii+level**numqd)))
      enddo
      sp_mu=-2.0*s_mu*sp_mu

      qd_mu = (0,0)
      do iq=1,numqd
         do ii=1,nstate,level**iq
            do jj=0,level**(iq-1)-1
               qd_mu = qd_mu + dimag(rho(ii+jj,ii+jj+level**(iq-1)))
            enddo
         enddo
      enddo
      qd_mu = -2.0*qd_mu*q_mu
      rinc = 4.0*pi*(omega(iEh)/timeunit)/(c_speed*eps_med**0.5)*
     $     d_per_au*1.0d-18/eesu_per_au
      sig   =rinc*(sp_mu+ qd_mu)/energy0/2.0
      rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $     d_per_au*1.0d-18/eesu_per_au
      sigr = rinc * dabs((sp_mu+qd_mu)/energy0/2)**2
      !A factor of three is needed to get it be the same as the rk_p code. It's probably
      !A units thing; they are off by almost exactly three at every points
      write(11,11) hb_omega(iEh),sig,sigr
c      print*,energy0,c_speed,eps_med,eesu_per_au
 10   format(1p6e20.8)
 11   format(1p10e20.8)

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
      real*8     glsc2,beta,tmp
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
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      
      integer  i,j,m,n,iEh
      real*8     v(n,m+1),b(n)
      real*8     h(m+1,m)
      real*8     eps, s1, si, glsc2, cnt

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

         if (abs(h(j+1,j)).le.eps.and.IFPULSE.eqv..false.) then
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
      real*8     h0(m+1,m),h(m,m)
      real*8     v0(n,m+1)
      real*8     v (n,m)
      real*8     vt(m,n)

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
      real*8     w(n2),v(n2)

      nstate2 = nstate*nstate

      call copy(rho_r, v(1+0*nstate2), nstate2)
      call copy(rho_i, v(1+1*nstate2), nstate2)

      call cem_quantum_hamiltonian(iEh)     
      call cem_quantum_linblad_m

      call copy(w(1+0*nstate2), resrho_r, nstate2) 
      call copy(w(1+1*nstate2), resrho_i, nstate2) 

      return
      end


c-----------------------------------------------------------------------
      subroutine f(time_o,v_ode,w_ode)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM' 
      integer  iEh,nstate2
      real*8   w_ode(2*nstate*nstate),v_ode(2*nstate*nstate),time_o

      nstate2 = nstate*nstate
      time = time_o
      call copy(rho_r, v_ode(1+0*nstate2), nstate2)
      call copy(rho_i, v_ode(1+1*nstate2), nstate2)

      call cem_quantum_hamiltonian(iEh)     
      call cem_quantum_linblad_m

      call copy(w_ode(1+0*nstate2), resrho_r, nstate2) 
      call copy(w_ode(1+1*nstate2), resrho_i, nstate2) 

      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian(iEh)
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer    ii,jj,iEh,i_one,i_zero
      real*8       tmprr,tmpii,pulseFactor,glamax
      complex*16 tmpfq,tmpcc
      complex*16 ci
     
      ci=(0.0,1.0)
      
      i_one = 1
      i_zero = 0
c     Symmetrize rho
      do jj=1,nstate
      do ii=1,nstate
         if(jj.gt.ii) then
            tmprr = (rho_r(ii,jj)+rho_r(jj,ii))/2
            rho_r(ii,jj) = tmprr
            rho_r(jj,ii) = tmprr
            tmpii = (rho_i(ii,jj)-rho_i(jj,ii))/2
            rho_i(ii,jj) = tmpii
            rho_i(jj,ii) = -tmpii
         endif
      enddo
      enddo

      do jj = 1,nstate
      do ii = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo
c      call outmat(resrho_r(1,1),nstate,nstate/2,'rho_r  ',i_zero)
c      call outmat(resrho_r(1,4),nstate,nstate/2,'rho_r  ',i_one)
c      call outmat(resrho_i(1,1),nstate,nstate/2,'rho_i  ',i_zero)
c      call outmat(resrho_i(1,4),nstate,nstate/2,'rho_i  ',i_one)
c      if(time.gt.3000) call exitt
      if(IFDENSE) then
         if(IFPULSE) then
            do ii=1,nstate
            do jj=1,nstate
               hamilt_I(ii,jj) = hamilt_0(ii,jj) + 
     $              pulseFactor(1)*hamilt(ii,jj)
            enddo
            enddo
         endif
         print*,'dense' 
         call mxm(hamilt_I,nstate,rho_r,nstate,tmp_r1,nstate)
         call mxm(rho_r,nstate,hamilt_I,nstate,tmp_r2,nstate)
         call mxm(hamilt_I,nstate,rho_i,nstate,tmp_i1,nstate)
         call mxm(rho_i,nstate,hamilt_I,nstate,tmp_i2,nstate)
      else
         if(IFPULSE) then
            do ii=1,nstate*n_in_row
               hamilt_I_m_p(ii) = hamilt_0_m_p(ii) + 
     $              pulseFactor(1)*hamilt_m_p(ii)
            enddo
         endif
c      open(unit=123,file='rho_r')
c      write(123,*)'rho'
c      do ii=1,nstate
c         write(123,"(100g15.5)") ( rho_r(jj,ii), jj=1,nstate )
c      enddo
c      open(unit=124,file='rho_i')
c      do ii=1,nstate
c         write(124,"(100g15.5)") ( rho_i(jj,ii), jj=1,nstate )
c      enddo

         call cem_quantum_mxm_m(rho_r,tmp_r1, 1)             
c      open(unit=123,file='rho_r')
c      write(123,*)'tmp'
c      write(124,*)''
c      do ii=1,nstate
c         write(123,"(100g15.5)") ( tmp_r1(jj,ii), jj=1,nstate )
c      enddo
c      open(unit=124,file='tmp_i1_all1')
c      do ii=1,nstate
c         write(124,"(100g15.5)") ( -tmp_i1(jj,ii), jj=1,nstate )
c      enddo
         call cem_quantum_mxm_m(rho_r,tmp_r2,-1)        
         call cem_quantum_mxm_m(rho_i,tmp_i1, 1)             
         call cem_quantum_mxm_m(rho_i,tmp_i2,-1)

      endif



c      call outmat(tmp_r2,nstate,nstate,'afterr ',nid)
c      call outmat(tmp_i2,nstate,nstate,'afteri ',nid)
c      if(time.gt.1000) call exitt

      do jj = 1,nstate
      do ii = 1,nstate
         if(IFPULSE) then
            tmpfq= state_sums(ii,jj)*pulse_omega_0*rho(ii,jj)
c            tmpfq = 0
         else
            tmpfq= frequency(ii,jj,iEh)*rho(ii,jj)
         endif
c         tmprr= tmp_r1(ii,jj)-tmp_r2(ii,jj)!+1e-17!tmp_r1(jj,ii)+1e-305
         tmprr= tmp_r1(ii,jj)-tmp_r1(jj,ii)
c         if(ii.ne.jj) tmprr=tmprr-1e-231
c         tmpii= tmp_i1(ii,jj)-tmp_i2(ii,jj)!+tmp_i1(jj,ii)
         tmpii= tmp_i1(ii,jj)+tmp_i1(jj,ii)
         tmpcc= tmprr + ci*tmpii
         drho(ii,jj)= - ci*(tmpcc - tmpfq) 
         resrho_r(ii,jj)=real(drho(ii,jj) )    
         resrho_i(ii,jj)=imag(drho(ii,jj) )    
      enddo
      enddo
c      do jj = 1,nstate
c      do ii = 1,nstate
c         tmp_r2(ii,jj) = tmp_r1(jj,ii)-tmp_r2(ii,jj)
c         tmp_i2(ii,jj) = -tmp_i1(jj,ii)-tmp_i2(ii,jj)
c      enddo
c      enddo
c      tmprr=glamax(tmp_r2,nstate*nstate)
c      tmpii=glamax(tmp_i2,nstate*nstate)
c      write(14,*) tmprr,tmpii
c      open(unit=123,file='tmp_r1_all1')
c      do ii=1,nstate
c         write(123,"(100g15.5)") ( tmp_r1(jj,ii), jj=1,nstate )
c      enddo
c      open(unit=124,file='tmp_i1_all1')
c      do ii=1,nstate
c         write(124,"(100g15.5)") ( -tmp_i1(jj,ii), jj=1,nstate )
c      enddo
c      write(123,*)''
c      write(124,*)''
c      call outmat(tmp_r2,nstate,nstate,'Hrho_r  ',i_zero)
c      call outmat(resrho_r(1,4),nstate,nstate/2,'Hrho_r  ',i_one)
c      call outmat(resrho_i(1,1),nstate,nstate,'Hrho_i  ',i_zero)
c      call outmat(resrho_i(1,4),nstate,nstate/2,'Hrho_i  ',i_one)
c      print*,time
c      call outmat(tmp_r2,nstate,nstate,'r2  ',nid)
c      call outmat(tmp_r2,nstate,nstate,'r2  ',nid)

c      call outmat(tmp_i2,nstate,nstate,'i2  ',nid)
c      call outmat(tmp_i2,nstate,nstate,'i2  ',nid)

c      call outmat(resrho_r,nstate,nstate,'res_r ',nid)
c      if(time.gt.500) call exitt
      return
      end


c--------------------------------------------------------------------- 
      subroutine cem_quantum_mxm_m(rrho,tmprr,idir)             
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer    ii,jj,idir,kk    
      real*8       rrho (lstate,lstate)
      real*8       tmprr(lstate,lstate)

      do jj = 1,nstate
      do ii = 1,nstate
         tmprr(ii,jj)=0
      enddo
      enddo
      
      if(idir.eq.1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
      do jj = 1,nstate
      do ii = 1,1
         do kk=1,n_in_row
            tmprr(ii,jj)=tmprr(ii,jj)+
     $           hamilt_I_m_p(n_in_row*(ii-1)+kk)*rrho(ii+kk-1,jj)
         enddo
         
      enddo
      enddo
      do jj = 1,nstate
      do ii = 2,n_in_row-1
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,n_in_row
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(ii+kk-1,jj)
         enddo
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,ii-1
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo

      enddo
      enddo
      
c     This is the bulk of the matrix
      do jj = 1,nstate
      do ii = n_in_row,nstate-(n_in_row-1)
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,n_in_row
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(ii+kk-1,jj)
         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,n_in_row-1
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
      enddo
      enddo

c     The end parts only have the up part, not the right part.
      do jj = 1,nstate
      do ii = nstate-(n_in_row-1)+1,nstate-1
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,nstate-(ii-1)
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(ii+kk-1,jj)
         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
      enddo
      enddo

c     Last element has only left of the diagonal terms
      do jj = 1,nstate
      do ii = nstate,nstate
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
         tmprr(ii,jj)=tmprr(ii,jj)+hamilt_I_m_p(n_in_row*(ii-1)+1)
     $        *rrho(ii,jj)
      enddo
      enddo

      elseif(idir.eq.-1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
      do jj = 1,nstate
      do ii = 1,1
         do kk=1,n_in_row
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(jj,ii+kk-1)
         enddo
         
      enddo
      enddo
      do jj = 1,nstate
      do ii = 2,n_in_row-1
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,n_in_row
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(jj,ii+kk-1)
         enddo
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,ii-1
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(jj,ii-kk)
         enddo

      enddo
      enddo
      
c     This is the bulk of the matrix
      do jj = 1,nstate
      do ii = n_in_row,nstate-(n_in_row-1)
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,n_in_row
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(jj,ii+kk-1)
         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,n_in_row-1
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(jj,ii-kk)
         enddo
      enddo
      enddo

c     The end parts only have the up part, not the right part.
      do jj = 1,nstate
      do ii = nstate-(n_in_row-1)+1,nstate-1
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,nstate-(ii-1)
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(jj,ii+kk-1)
         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(jj,ii-kk)
         enddo
      enddo
      enddo

c     Last element has only left of the diagonal terms
      do jj = 1,nstate
      do ii = nstate,nstate
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(jj,ii-kk)
         enddo
         tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+1)
     $           *rrho(jj,ii)
      enddo
      enddo
      else
         
         write(6,*) 'no idir: cem_quantum_mxm'
         call exitt
         
      endif
            
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



      subroutine ode(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
      implicit real*8(a-h,o-z)
c Shampine's differential equations integrator
c   subroutine ode integrates a system of neqn
c   first order ordinary differential equations of the form:
c             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
c             y(i) given at  t .
c   the subroutine integrates from  t  to  tout .  on return the
c   parameters in the call list are set for continuing the integration.
c   the user has only to define a new value  tout  and call  ode  again.
c
c   the differential equations are actually solved by a suite of codes
c   de ,  step , and  intrp .  ode  allocates virtual storage in the
c   arrays  work  and  iwork  and calls  de .  de  is a supervisor which
c   directs the solution.  it calls on the routines  step  and  intrp
c   to advance the integration and to interpolate at output points.
c   step  uses a modified divided difference form of the adams pece
c   formulas and local extrapolation.  it adjusts the order and step
c   size to control the local error per unit step in a generalized
c   sense.  normally each call to  step  advances the solution one step
c   in the direction of  tout .  for reasons of efficiency  de
c   integrates beyond  tout  internally, though never beyond
c   t+10*(tout-t), and calls  intrp  to interpolate the solution at
c   tout .  an option is provided to stop the integration at  tout  but
c   it should be used only if it is impossible to continue the
c   integration beyond  tout .
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c   the parameters represent:
c      f -- real subroutine f(t,y,yp) to evaluate
c                derivatives yp(i)=dy(i)/dt
c      neqn -- number of equations to be integrated (integer)
c      y(*) -- solution vector at t                 (real)
c      t -- independent variable                    (real)
c      tout -- point at which solution is desired   (real)
c      relerr,abserr -- relative and absolute error tolerances for local
c           error test (real).  at each step the code requires
c             abs(local error) .le. abs(y)*relerr + abserr
c           for each component of the local error and solution vectors
c      iflag -- indicates status of integration     (integer)
c      work(*)  (real)  -- arrays to hold information internal to
c      iwork(*) (integer)    which is necessary for subsequent calls
c
c   first call to ode --
c
c   the user must provide storage in his calling program for the arrays
c   in the call list,
c      y(neqn), work(100+21*neqn), iwork(5),
c   declare  f  in an external statement, supply the real
c   subroutine f(t,y,yp)  to evaluate
c      dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
c   and initialize the parameters:
c      neqn -- number of equations to be integrated
c      y(*) -- vector of initial conditions
c      t -- starting point of integration
c      tout -- point at which solution is desired
c      relerr,abserr -- relative and absolute local error tolerances
c      iflag -- +1,-1.  indicator to initialize the code.  normal input
c           is +1.  the user should set iflag=-1 only if it is
c           impossible to continue the integration beyond  tout .
c   all parameters except  f ,  neqn  and  tout  may be altered by the
c   code on output so must be variables in the calling program.
c
c   output from  ode  --
c
c      neqn -- unchanged
c      y(*) -- solution at  t
c      t -- last point reached in integration.  normal return has
c           t = tout .
c      tout -- unchanged
c      relerr,abserr -- normal return has tolerances unchanged.  iflag=3
c           signals tolerances increased
c      iflag = 2 -- normal return.  integration reached  tout
c            = 3 -- integration did not reach  tout  because error
c                   tolerances too small.  relerr ,  abserr  increased
c                   appropriately for continuing
c            = 4 -- integration did not reach  tout  because more than
c                   500 steps needed
c            = 5 -- integration did not reach  tout  because equations
c                   appear to be stiff
c            = 6 -- invalid input parameters (fatal error)
c           the value of  iflag  is returned negative when the input
c           value is negative and the integration does not reach  tout ,
c           i.e., -3, -4, -5.
c      work(*),iwork(*) -- information generally of no interest to the
c           user but necessary for subsequent calls.
c
c   subsequent calls to  ode --
c
c   subroutine  ode  returns with all information needed to continue
c   the integration.  if the integration reached  tout , the user need
c   only define a new  tout  and call again.  if the integration did not
c   reach  tout  and the user wants to continue, he just calls again.
c   the output value of  iflag  is the appropriate input value for
c   subsequent calls.  the only situation in which it should be altered
c   is to stop the integration internally at the new  tout , i.e.,
c   change output  iflag=2  to input  iflag=-2 .  error tolerances may
c   be changed by the user before continuing.  all other parameters must
c   remain unchanged.
c
c***********************************************************************
c*  subroutines  de  and  step  contain machine dependent constants.   *
c*  be sure they are set before using  ode .                           *
c***********************************************************************
      logical start,phase1,nornd
      dimension y(neqn),work(1),iwork(5)
      external f
      data ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart,
     1  itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/
      iyy = 100
      iwt = iyy + neqn
      ip = iwt + neqn
      iyp = ip + neqn
      iypout = iyp + neqn
      iphi = iypout + neqn
      if(abs(iflag) .eq. 1) go to 1
      start = work(istart) .gt. 0.0d0
      phase1 = work(iphase) .gt. 0.0d0
      nornd = iwork(2) .ne. -1
    1 call de(f,neqn,y,t,tout,relerr,abserr,iflag,work(iyy),
     1  work(iwt),work(ip),work(iyp),work(iypout),work(iphi),
     2  work(ialpha),work(ibeta),work(isig),work(iv),work(iw),work(ig),
     3  phase1,work(ipsi),work(ix),work(ih),work(ihold),start,
     4  work(itold),work(idelsn),iwork(1),nornd,iwork(3),iwork(4),
     5  iwork(5))
      work(istart) = -1.0d0
      if(start) work(istart) = 1.0d0
      work(iphase) = -1.0d0
      if(phase1) work(iphase) = 1.0d0
      iwork(2) = -1
      if(nornd) iwork(2) = 1
      return
      end
c
c
      subroutine de(f,neqn,y,t,tout,relerr,abserr,iflag,
     1  yy,wt,p,yp,ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold,
     2  start,told,delsgn,ns,nornd,k,kold,isnold)
      implicit real*8(a-h,o-z)
      logical stiff,crash,start,phase1,nornd
      dimension y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),
     1  ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),w(12),g(13)
      external f
c
c***********************************************************************
c*  the only machine dependent constant is based on the machine unit   *
c*  roundoff err  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
c*  in the following data statement before using  de .  the routine    *
c*  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
c*  inserted in subroutine  step1  before calling  de .
      data fouru/.888e-15/
c***********************************************************************
c
c   the constant  maxnum  is the maximum number of steps allowed in one
c   call to  de .  the user may change this limit by altering the
c   following statement
      data maxnum/10000/
c
c            ***            ***            ***
c   test for improper parameters
c
      if(neqn .lt. 1) go to 10
      if(t .eq. tout) go to 10
      if(relerr .lt. 0.0d0  .or.  abserr .lt. 0.0d0) go to 10
      eps = max(relerr,abserr)
      if(eps .le. 0.0d0) go to 10
      if(iflag .eq. 0) go to 10
      isn = sign(1,iflag)
      iflag = abs(iflag)
      if(iflag .eq. 1) go to 20
      if(t .ne. told) go to 10
      if(iflag .ge. 2  .and.  iflag .le. 5) go to 20
   10 iflag = 6
      return
c
c   on each call set interval of integration and counter for number of
c   steps.  adjust input err tolerances to define weight vector for
c   subroutine  step1
c
   20 del = tout - t
      absdel = abs(del)
      tend = t + 10.0d0*del
      if(isn .lt. 0) tend = tout
      nostep = 0
      kle4 = 0
      stiff = .false.
      releps = relerr/eps
      abseps = abserr/eps
      if(iflag .eq. 1) go to 30
      if(isnold .lt. 0) go to 30
      if(delsgn*del .gt. 0.0d0) go to 50
c
c   on start and restart also set work variables x and yy(*), store the
c   direction of integration and initialize the step size
c
   30 start = .true.
      x = t
      do 40 l = 1,neqn
   40   yy(l) = y(l)
      delsgn = 1!sign(1.0d0,del)
      h = 1!sign(max(abs(tout-x),fouru*abs(x)),tout-x)
c
c   if already past output point, interpolate and return
c
   50 if(abs(x-t) .lt. absdel) go to 60
c      if(nid.eq.0)print*,'before intrp'
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
c
c   if cannot go past output point and sufficiently close,
c   extrapolate and return
c
   60 if(isn .gt. 0  .or.  abs(tout-x) .ge. fouru*abs(x)) go to 80
      h = tout - x
      call f(x,yy,yp)
      do 70 l = 1,neqn
   70   y(l) = yy(l) + h*yp(l)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
c
c   test for too many steps
c
   80 if(nostep .lt. maxnum) go to 100
      iflag = isn*4
      if(stiff) iflag = isn*5
      do 90 l = 1,neqn
   90   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
c
c   limit step size, set weight vector and take a step
c
  100 h = sign(min(abs(h),abs(tend-x)),h)
      do 110 l = 1,neqn
  110   wt(l) = releps*abs(yy(l)) + abseps
      call step1(x,yy,f,neqn,h,eps,wt,start,
     1  hold,k,kold,crash,phi,p,yp,psi,
     2  alpha,beta,sig,v,w,g,phase1,ns,nornd)
c
c   test for tolerances too small
c
      if(.not.crash) go to 130
      iflag = isn*3
      relerr = eps*releps
      abserr = eps*abseps
      do 120 l = 1,neqn
  120   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
c
c   augment counter on number of steps and test for stiffness
c
  130 nostep = nostep + 1
      kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      go to 50
      end
c
      subroutine step1(x,y,f,neqn,h,eps,wt,start,
     1  hold,k,kold,crash,phi,p,yp,psi,
     2  alpha,beta,sig,v,w,g,phase1,ns,nornd)
      implicit real*8(a-h,o-z)
      logical start,crash,phase1,nornd
      dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
      dimension alpha(12),beta(12),sig(13),w(12),v(12),g(13),
     1  gstr(13),two(13)
      external f
c***********************************************************************
c*  the only machine dependent constants are based on the machine unit *
c*  roundoff err  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          *
c*  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling *
c*  the code.  the routine  machin  calculates  u .                    *
      data twou,fouru/.444e-15,.888e-15/
c***********************************************************************
      data two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,
     1  512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
      data gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0,
     1  0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,
     2  0.00468d0/
c
c
c       ***     begin block 0     ***
c   check if step size or err tolerance is too small for machine
c   precision.  if first step, initialize phi array and estimate a
c   starting step size.
c                   ***
c
c   if step size is too small, determine an acceptable one
c
      crash = .true.
      if(abs(h) .ge. fouru*abs(x)) go to 5
      h = sign(fouru*abs(x),h)
      return
    5 p5eps = 0.5d0*eps
c
c   if err tolerance is too small, increase it to an acceptable value
c
      round = 0.0d0
      do 10 l = 1,neqn
   10   round = round + (y(l)/wt(l))**2
      round = twou*sqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      return
   15 crash = .false.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      if(.not.start) go to 99
c
c   initialize.  compute appropriate step size for first step
c
      call f(x,y,yp)
      sum = 0.0d0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
   20   sum = sum + (yp(l)/wt(l))**2
      sum = sqrt(sum)
      absh = abs(h)
      if(eps .lt. 16.0d0*sum*h*h) absh = 0.25d0*sqrt(eps/sum)
      h = sign(max(absh,fouru*abs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0d0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
   25   phi(l,15) = 0.0d0
   99 ifail = 0
c       ***     end block 0     ***
c
c       ***     begin block 1     ***
c   compute coefficients of formulas for this step.  avoid computing
c   those quantities not changed when step size is not changed.
c                   ***
c
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
c
c   ns is the number of steps taken with size h, including the current
c   one.  when k.lt.ns, no coefficients change
c
      if(h .ne. hold) ns = 0
      if(ns.le.kold)   ns=ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
c
c   compute those components of alpha(*),beta(*),psi(*),sig(*) which
c   are changed
c
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
  105   sig(i+1) = reali*alpha(i)*sig(i)
  110 psi(k) = temp1

c
c   compute coefficients g(*)
c
c   initialize v(*) and set w(*).  g(2) is set in data statement
c
      if(ns .gt. 1) go to 120
      do 115 iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
  115   w(iq) = v(iq)
      go to 140
c
c   if order was raised, update diagonal part of v(*)
c
  120 if(k .le. kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      if(nsm2 .lt. 1) go to 130
      do 125 j = 1,nsm2
        i = k-j
  125   v(i) = v(i) - alpha(j+1)*v(i+1)
c
c   update v(*) and set w(*)
c
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
  135   w(iq) = v(iq)
      g(nsp1) = w(1)
c
c   compute the g(*) in the work vector w(*)
c
  140 nsp2 = ns + 2
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
  145     w(iq) = w(iq) - temp6*w(iq+1)
  150   g(i) = w(1)
  199   continue

c       ***     end block 1     ***
c
c       ***     begin block 2     ***
c   predict a solution p(*), evaluate derivatives using predicted
c   solution, estimate local err at order k and errs at orders k,
c   k-1, k-2 as if constant step size were used.
c                   ***
c
c   change phi to phi star
c
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
  205     phi(l,i) = temp1*phi(l,i)
  210   continue
c
c   predict solution and differences
c
  215 do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
  220   p(l) = 0.0d0
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
  225     phi(l,i) = phi(l,i) + phi(l,ip1)
  230   continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
  235   phi(l,16) = (p(l) - y(l)) - tau
      go to 250
  240 do 245 l = 1,neqn
  245   p(l) = y(l) + h*p(l)
  250 xold = x
      x = x + h
      absh = abs(h)
      call f(x,p,yp)
c
c   estimate errs at orders k,k-1,k-2
c
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        if(km2)265,260,255
  255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
  260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
  265   erk = erk + (temp4*temp3)**2
      if(km2)280,275,270
  270 erkm2 = absh*sig(km1)*gstr(km2)*sqrt(erkm2)
  275 erkm1 = absh*sig(k)*gstr(km1)*sqrt(erkm1)
  280 temp5 = absh*sqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
c
c   test if order should be lowered
c
      if(km2)299,290,285
  285 if(max(erkm1,erkm2) .le. erk) knew = km1
      go to 299
  290 if(erkm1 .le. 0.5d0*erk) knew = km1
c
c   test if step successful
c
  299 if(err .le. eps) go to 400

c       ***     end block 2     ***
c
c       ***     begin block 3     ***
c   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
c   if third consecutive failure, set order to one.  if step fails more
c   than three times, consider an optimal step size.  double err
c   tolerance and return if estimated step size is too small for machine
c   precision.
c                   ***
c
c   restore x, phi(*,*) and psi(*)
c
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
  305     phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  310   continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
  315   psi(i-1) = psi(i) - h
c
c   on third failure, set order to one.  thereafter, use optimal step
c   size
c
  320 ifail = ifail + 1
      temp2 = 0.5d0
      if(ifail - 3) 335,330,325
  325 if(p5eps .lt. 0.25d0*erk) temp2 = sqrt(p5eps/erk)
  330 knew = 1
  335 h = temp2*h
      k = knew
      if(abs(h) .ge. fouru*abs(x)) go to 340
      crash = .true.
      h = sign(fouru*abs(x),h)
      eps = eps + eps
      return
  340 go to 100

c       ***     end block 3     ***
c
c       ***     begin block 4     ***
c   the step is successful.  correct the predicted solution, evaluate
c   the derivatives using the corrected solution and update the
c   differences.  determine best order and step size for next step.
c                   ***
  400 kold = k
      hold = h
c
c   correct and evaluate
c
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do 405 l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
  405   phi(l,15) = (y(l) - p(l)) - rho
      go to 420
  410 do 415 l = 1,neqn
  415   y(l) = p(l) + temp1*(yp(l) - phi(l,1))
  420 call f(x,y,yp)
c
c   update differences for next step
c
      do 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
  425   phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      do 435 i = 1,k
        do 430 l = 1,neqn
  430     phi(l,i) = phi(l,i) + phi(l,kp1)
  435   continue
c
c   estimate err at order k+1 unless:
c     in first phase when always raise order,
c     already decided to lower order,
c     step size not constant so estimate unreliable
c
      erkp1 = 0.0d0
      if(knew .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(knew .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
  440   erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*sqrt(erkp1)
c
c   using estimated err at order k+1, determine appropriate order
c   for next step
c
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5d0*erk) go to 460
      go to 450
  445 if(erkm1 .le. min(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
c
c   here erkp1 .lt. erk .lt. max(erkm1,erkm2) else order would have
c   been lowered in block 2.  thus order is to be raised
c
c   raise order
c
  450 k = kp1
      erk = erkp1
      go to 460
c
c   lower order
c
  455 k = km1
      erk = erkm1
c
c   with new order determine appropriate step size for next step
c
  460 hnew = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hnew = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = 1!absh*max(0.5d0,min(0.9d0,r))
      hnew = sign(max(hnew,fouru*abs(x)),h)
  465 h = hnew
      return
c       ***     end block 4     ***
      end
c
      subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)
      implicit real*8(a-h,o-z)
       dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
      dimension g(13),w(13),rho(13)
      data g(1)/1.0d0/,rho(1)/1.0d0/
c
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
c
c   initialize w(*) for computing g(*)
c
      do 5 i = 1,ki
        temp1 = i
    5   w(i) = 1.0d0/temp1
      term = 0.0d0
c
c   compute g(*)
c
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
   10     w(i) = gamma*w(i) - eta*w(i+1)
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
   15   term = psijm1
c
c   interpolate
c
      do 20 l = 1,neqn
        ypout(l) = 0.0d0
   20   yout(l) = 0.0d0
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
   25     ypout(l) = ypout(l) + temp3*phi(l,i)
   30   continue
      do 35 l = 1,neqn
   35   yout(l) = y(l) + hi*yout(l)
      return
      end
