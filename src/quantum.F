c---------------------------------------------------------------------
      subroutine genbox_quantum
c---------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'ZPER'
c
      character*132 string
      character*1  string1(132)
      equivalence (string,string1)
C
c
      call blank (string,132)
 132  format(a132)
c
      if (nid.eq.0) then
c
c        Read in name of previously generated NEKTON data set.
c        (Must have same dimension and number of fields as current run)
c
         write(6,*) 'genbox_ode:: opening file'
         lfname = ltrunc(reafle,132) - 4
         call blank (string,132)
         call chcopy(string,reafle,lfname)
         call chcopy(string1(lfname+1),'.box',4)
         open (unit=7,file=string,status='old',err=999)
         call gets(string,132,iend,7)  ! This is a dummy read, for compatibility
c
c        here is where the 2d/3d determination is made....
c
         call geti1(ndim,iend,7)
         ndim = iabs(ndim)

         if (ndim.eq.0) then
             call read_quantum_input
         else
             write(6,*) 'set ndim=0 in box file'
             call exitt
         endif
         close (unit=7)

      endif

      return

 999  continue

      if (nid.eq.0) write(6,*) 'ABORT: Could not find box file ',string
      call exitt
      end

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
      common /QRDATA_R5/  r_pulse_duration,r_pulse_fluence,r_gcouple_q3
      common /QLDATA_L1/  i_IFSMALLMEMORY
      integer             i_numqd,i_numsp,i_Eh,i_nqstep
      real               r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      real                rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      real                r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      real                r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      real                r_pulse_duration,r_pulse_fluence,r_gcouple_q3
      integer             iend
      logical             i_IFSMALLMEMORY


c...  this routine reads the data in *.box file for quantum casce
c...  the first integer "0" case
     

         call geti1(i_numqd,iend,7)
         call geti1(i_numsp,iend,7)
         
         call getrv(rintense,1,iend,7)
         call getr3(r_Eh1,r_Eh2,r_Eh,iend,7)
         call getr2(r_tmax,r_nqstep,iend,7)
         call getr2(er1,er2,iend,7)
         call getr2(r_omega_q1,r_omega_q2,iend,7)
         call getr3(r_gcouple_q1,r_gcouple_q2,r_gcouple_q3,iend,7)
         call getr2(r_gamma_q1,r_gamma_q2,iend,7)
         call getr2(r_gamma_d1,r_gamma_d2,iend,7)
         call getr2(r_q_mu_1,r_q_mu_2,iend,7)
         call getrv(r_pulse_t_0,1,iend,7)
         call getrv(r_pulse_duration,1,iend,7)
         call getrv(r_pulse_fluence,1,iend,7)
         call geti1(IFPULSE,iend,7)
         call geti1(IFDENSE,iend,7)
         call geti1(i_IFSMALLMEMORY,iend,7)
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
     $        r_gcouple_q2,r_gcouple_q3
         write(6,*) 'gamma_q1,gamma_q2 = ',r_gamma_q1,r_gamma_q2
         write(6,*) 'gamma_d1,gamma_d2 = ',r_gamma_d1,r_gamma_d2
         write(6,*) 'q_mu_1,q_mu_2 = ',r_q_mu_1,r_q_mu_2
         write(6,*) 'pulse_t_0 = ',r_pulse_t_0
         write(6,*) 'pulse_duration = ',r_pulse_duration
         write(6,*) 'pulse_fluence = ',r_pulse_fluence
         write(6,*) 'IFPULSE =',IFPULSE,' IFDENSE =',IFDENSE
         write(6,*) 'IFSMALLMEMORY =',IFSMALLMEMORY,' IFRKQ =',IFRKQ
         
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
      common /QRDATA_R5/  r_pulse_duration,r_pulse_fluence,r_gcouple_q3
      integer             i_numqd,i_numsp,i_Eh,i_nqstep,ii
      real                r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      real                rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      real                r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      real                r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      real                r_pulse_duration,r_pulse_fluence,r_gcouple_q3

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
      ltonq   = level**numqd

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
      printdt = timemax/(999)             !We only need to print 1000 points
      time    = 0
      timeOld = -1

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
c    Im not sure that this is the proper way to do this

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
            g_couple_m(ii) = (1d-2)/eV 
         enddo
         g_couple_m(4) = 1d-2/eV
         g_couple_m(3) = r_gcouple_q3/ev
         g_couple_m(2) = r_gcouple_q2/ev
         g_couple_m(1) = r_gcouple_q1/ev
         if (nid.eq.0) then 
            write(6,*) 'g_couple = ',( g_couple_m(ii)*eV, ii=1,numqd )
            write(14,1234) '#g=',( g_couple_m(ii)*eV, ii=1,numqd )
            write(13,1234) '#g=',( g_couple_m(ii)*eV, ii=1,numqd )
            write(12,1234) '#g=',( g_couple_m(ii)*eV, ii=1,numqd )
            
 1234       format(' ',A3,10E9.2)
         endif
         
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
     $     (time-pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-2*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-4*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-6*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-8*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-10*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-12*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-14*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-16*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-18*pulse_t_0)**2)
c$$$     $ +energy0_pls*exp(-alpha*(time-20*pulse_t_0)**2)
     $     )*cos(pulse_omega_0*(time-pulse_t_0))

      if(abs(fullPulse).lt.1e-14) fullPulse=0

      end

c---------------------------------------------------------------------

      real*8 function pulseFactor(n_pls)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'
      integer n_pls
      real*8    alpha

      alpha = 2. * log(2.0) / (pulse_duration)**2

      pulseFactor = energy0_pls/2.0*exp(-alpha*(time-pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-2*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-4*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-6*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-8*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-10*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-12*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-14*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-16*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-18*pulse_t_0)**2)
c$$$     $ +energy0_pls/2.0*exp(-alpha*(time-20*pulse_t_0)**2)
c     $     +energy0_pls/2.0*exp(-alpha*(time-2*pulse_t_0)**2)
c     $     *cos(pulse_omega_0*(time-pulse_t_0))

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
      
c      open(unit=123,file='matrix')
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
      real     fac,tmp,tmp3,tmp1,tmp2
      logical  lbFlag

!omp parallel do
      do jj = 1,nstate
      do ii = 1,nstate
         ns = isstate(ii)
         nsp= isstate(jj)
         do iq=1,numqd
            nqp(iq) = iqstate_m(iq,jj)
         enddo
         tmp3=gamma_s*(ns+nsp)*0.5
         do iq=1,numqd
            nq(iq) = iqstate_m(iq,ii)         
            tmp1=gamma_d_m(iq)*(nq(iq)**2+nqp(iq)**2-2*nq(iq)*nqp(iq))
            tmp2=gamma_q_m(iq)*(nq(iq)+nqp(iq))*0.5 ! 10/13/14: FIXME tmp2 has to be defined
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
!omp end parallel do
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
        write(12,*)'#Time,sp_mu,qd_mu(1),qd_mu(2),
     $mu_tot,pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm,concurrence'

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
 
      stime1=0
      stime4=0
      stime5=0
      stime0=dclock()
      if(IFPULSE) then
         call cem_quantum_initialize
c     Write file headers 
        write(12,*)'#Time,sp_mu,qd_mu(1),qd_mu(2),
     $mu_tot,pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm,concurrence'
         write(11,*)'#eV,sig_abs,sig_scatt,sig_tot'
         if(numqd.eq.4) then 
            write(13,*),'#Time,4:3,4:2,3:2,4:1,2:1,3:1'
            write(14,*),'#Time,QD1,QD2,QD3,QD4,MNP'
            write(15,*),'#Time,4:3,4:2,3:2,4:1,2:1,3:1'
         endif
         if(numqd.eq.3) then
            write(13,*),'#Time,3:2,3:1,2:1'
            write(14,*),'#Time,QD1,QD2,QD3,MNP'
            write(15,*),'#Time,3:2,3:1,2:1'
         endif


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

c            print*,time,sig
c            print*,time,rho_r(1,11),rho_i(2,11)
            time= time+dt
            if(time.gt.timeOld) then
               call cem_quantum_fourier_pls
     $              (s_mu_frequency,q_mu_frequency,energy_t)
               stime2= dclock()-stime2
               timeOld = timeOld+printdt
            endif
            
         enddo
         
         stime1= dclock()-stime1
         stime3= stime2/neqn
c         open(unit=11,file='spec.out')
c     Print out results

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
               call cem_quantum_fourier_pls
     $           (s_mu_frequency,q_mu_frequency,energy_t)    
               stime2= dclock()-stime2
            enddo
            stime1= dclock()-stime1
            call cem_quantum_post_m(iEh)            
         enddo
      endif
      
      stime0 = dclock()-stime0
      write(6,6) stime0,stime1,stime2,stime3 
 6    format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      close(11)
      close(12)
      close(13)
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
c$$$      if((time.gt.10000).and.(time.lt.10500)) then
c$$$      call outmat(resrho_r,nstate,nstate,'res_r  ',nid)
c$$$      call outmat(resrho_i,nstate,nstate,'res_i  ',nid)
c$$$      call outmat(rho_r,nstate,nstate,'rho_r  ',nid)
c$$$      call outmat(rho_i,nstate,nstate,'rho_r  ',nid)
c$$$      endif
c$$$      if(time.gt.10500) call exitt

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

c$$$      do jj=1,nstate
c$$$         do ii=1,nstate
c$$$            rho_r(ii,jj) = 1
c$$$            if(ii.gt.jj) then
c$$$               rho_i(ii,jj) = -1
c$$$            else if(ii.eq.jj) then
c$$$               rho_i(ii,jj) = 0
c$$$            else
c$$$               rho_i(ii,jj) = 1
c$$$            endif
c$$$         enddo
c$$$      enddo
c     Initial conditions; change here to rho_r(2,2) for state |10>, etc
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
      integer  iEh,ii,jj,iq,kk,jjplus,kkplus
      integer  nq(lnumqd),nqp,info,twotoiq
      complex*16 ci
      complex*16 tmpCMatrix(lstate,lstate),tmpCMatrix2(lstate,lstate)
      real*8     pulseFactor,fullPulse,currentPop,popSymm,popAsymm
      real*8     tmpC1,tmpC2,concurrence,con2,conc4
      real*8     popS,popQD(lnumqd),rho_workr(2*ltonq),conc1,conc2,conc3
      real*8     popSymm1(4),popSymm2(4),popSymm3(4)
      real*8     popSymm4(4),popSymm5(4),popSymm6(4)
      real       qd_mu(lnumqd),sp_mu,signRho,rho_workr2(2*4),conc5,conc6
      complex    s_mu_freq(lEh),q_mu_freq(lEh),en_t(lEh),rhoOrdering
      complex    red_rho(ltonq,ltonq),rho_svd(ltonq,ltonq)
      complex    rho_work(2*ltonq),rhoVL(ltonq,ltonq)
      complex    rhoVR(ltonq,ltonq),rhoEval(ltonq)
      complex    spinFlip(ltonq,ltonq),rhoFinal(ltonq,ltonq)
c      complex    red_rho2(4,4),rho_svd2(4,4),rho_work2(8),rhoVL2(4,4)
c      complex    spinFlip2(4,4),rhoFinal2(4,4),rhoVR2(4,4),rhoEval2(4)
      complex    red_rho2(4,4),red_rho3(8,8),bell(4,4)
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
      if(numqd.ne.1) then
         do ii=1,numsp
            popSymm = popSymm + 0.5*real(rho(4*ii-2,4*ii-2)
     $           + rho(4*ii-2,4*ii-1) + rho(4*ii-1,4*ii-2)
     $           + rho(4*ii-1,4*ii-1))
            popAsymm = popAsymm + 0.5*real(rho(4*ii-2,4*ii-2)
     $           - rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2)
     $           + rho(4*ii-1,4*ii-1))  
         enddo
c     Find the concurrence; see Wooters
         call czero(red_rho,ltonq*ltonq)
         do ii=1,numsp
            do jj=1,ltonq
               do kk=1,ltonq
                  red_rho(kk,jj) = red_rho(kk,jj) + 
     $                 rho(ltonq*ii-(ltonq-kk),ltonq*ii-(ltonq-jj))
               enddo
            enddo
         enddo
         

         call czero(spinFlip,ltonq*ltonq)
         do ii=1,ltonq
            do jj=1,ltonq
               rhoVL(ii,jj) = conjg(red_rho(ii,jj))
            enddo
         enddo

          if(numqd.eq.2) then
          
          spinFlip(1,4) = -1
          spinFlip(4,1) = -1
          spinFlip(2,3) = 1
          spinFlip(3,2) = 1

          call czero(bell,16)
          bell(1,4) = -2*sqrt(2.0)
          bell(4,1) = -2*sqrt(2.0)
          bell(2,3) = 2*sqrt(2.0)
          bell(3,2) = 2*sqrt(2.0)

          call zgemm('N','N',4,4,4,(1.0,0.0),bell,4,
     $        red_rho,4,(0.0,0.0),rhoFinal,4)
          conc5=0
          do ii=1,4
             conc5 = conc5 + rhoFinal(ii,ii)
          enddo
          

c         Find Conditional Probabilities, according to Stephens note

          conc1 = red_rho(2,2)/(red_rho(2,2)+red_rho(4,4))
          conc2 = red_rho(4,4)/(red_rho(2,2)+red_rho(4,4))
          conc3 = red_rho(1,1)/(red_rho(1,1)+red_rho(3,3))
          conc4 = red_rho(3,3)/(red_rho(1,1)+red_rho(3,3))

          write(13,13)time*timeunit,conc1,conc2,conc3,conc4,conc5
          
          endif
          if(numqd.eq.3) then
 
c          call czero(red_rho,ltonq*ltonq)
c         W STATE
c          red_rho(2,2) = 1.0/3.0
c          red_rho(3,3) = 1.0/3.0
c          red_rho(5,5) = 1.0/3.0
c          red_rho(2,5) = 1.0/3.0
c          red_rho(5,2) = 1.0/3.0
c          red_rho(3,5) = 1.0/3.0
c          red_rho(5,3) = 1.0/3.0
c          red_rho(2,3) = 1.0/3.0
c          red_rho(3,2) = 1.0/3.0
c         MAX ENTANGLE QBIT 1          
c          red_rho(1,1) = 1.0/2.0
c          red_rho(7,7) = 1.0/2.0
c          red_rho(1,7) = 1.0/2.0
c          red_rho(7,1) = 1.0/2.0
c         MAX ENTANGLE QBIT3
c          red_rho(1,1) = 1.0/2.0
c          red_rho(4,4) = 1.0/2.0
c          red_rho(1,4) = 1.0/2.0
c          red_rho(4,1) = 1.0/2.0
c        MAX ENTANGLE QBIT2
c          red_rho(1,1) = 1.0/2.0
c          red_rho(1,6) = 1.0/2.0
c          red_rho(6,1) = 1.0/2.0
c          red_rho(6,6) = 1.0/2.0
             
c$$$             call outmat(rho_r,nstate,nstate,'rho_r',nid)
c$$$             call outmat(rho_i,nstate,nstate,'rho_i',nid)
c$$$             call outmat(red_rho,8,8,'red_rho',nid)
c$$$             call exitt()
             call cem_quantum_partial_trace_3_to_2(red_rho,red_rho2,1)
c     concurrence(3:2)
             call cem_quantum_bipartite_concurrence(red_rho2,conc1
     $            ,popSymm1)
             call cem_quantum_partial_trace_3_to_2(red_rho,red_rho2,2)
c     concurrence(3:1)
             call cem_quantum_bipartite_concurrence(red_rho2,conc2
     $            ,popSymm2)
             call cem_quantum_partial_trace_3_to_2(red_rho,red_rho2,3)
c     concurrence(2:1)
             call cem_quantum_bipartite_concurrence(red_rho2,conc3
     $            ,popSymm3)
             
             
             write(13,13)time*timeunit,conc1,conc2,conc3
             write(14,13)time*timeunit,popQD,popS
             write(15,13)time*timeunit,popSymm1,popSymm2
     $            ,popSymm3


 13          format(1p20g20.8)

             
          endif
          if(numqd.eq.4) then
c          call czero(red_rho,ltonq*ltonq)
c         W STATE
c          red_rho(2,2) = 1.0/3.0
c          red_rho(3,3) = 1.0/3.0
c          red_rho(5,5) = 1.0/3.0
c          red_rho(2,5) = 1.0/3.0
c          red_rho(5,2) = 1.0/3.0
c          red_rho(3,5) = 1.0/3.0
c          red_rho(5,3) = 1.0/3.0
c          red_rho(2,3) = 1.0/3.0
c          red_rho(3,2) = 1.0/3.0
c         MAX ENTANGLE QBIT 1          
c          red_rho(1,1) = 1.0/2.0
c          red_rho(7,7) = 1.0/2.0
c          red_rho(1,7) = 1.0/2.0
c          red_rho(7,1) = 1.0/2.0
c         MAX ENTANGLE QBIT3
c          red_rho(1,1) = 1.0/2.0
c          red_rho(4,4) = 1.0/2.0
c          red_rho(1,4) = 1.0/2.0
c          red_rho(4,1) = 1.0/2.0
c        MAX ENTANGLE QBIT2
c          red_rho(1,1) = 1.0/2.0
c          red_rho(1,6) = 1.0/2.0
c          red_rho(6,1) = 1.0/2.0
c          red_rho(6,6) = 1.0/2.0

c     First, trace out qd 1 with 4->3 quantum dots
          call cem_quantum_partial_trace_4_to_3(red_rho,red_rho3,1)
c     Now, trace 3->2 3 times, and find bipartite concurrence for each
c     First, concurrence(4:3) (trace out 2)
          call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,1)
          call cem_quantum_bipartite_concurrence(red_rho2,conc1,popSymm1
     $         )
c     concurrence(4:2) (trace out 3)
          call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,2)
          call cem_quantum_bipartite_concurrence(red_rho2,conc2,popSymm2
     $         )
c     concurrence(3:2) (trace out 4)
          call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,3)
          call cem_quantum_bipartite_concurrence(red_rho2,conc3,popSymm3
     $         )
c     trace out qd 3 from full
          call cem_quantum_partial_trace_4_to_3(red_rho,red_rho3,3)
c     concurrence(4:1)
          call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,2)
          call cem_quantum_bipartite_concurrence(red_rho2,conc4,popSymm4
     $         )
c     concurrence(2:1)
          call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,3)
          call cem_quantum_bipartite_concurrence(red_rho2,conc5,popSymm5
     $         )
c     trace out qd 4 from full
          call cem_quantum_partial_trace_4_to_3(red_rho,red_rho3,4)
c     concurrence(3:1)
          call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,3)
          call cem_quantum_bipartite_concurrence(red_rho2,conc6,popSymm6
     $         )


         write(13,13)time*timeunit,conc1,conc2,conc3,conc4,conc5,conc6
         write(14,13)time*timeunit,popQD,
     $        popS
         write(15,13)time*timeunit,popSymm1,popSymm2
     $        ,popSymm3,popSymm4,popSymm5,popSymm6


          spinFlip(1,16) = 1
          spinFlip(2,15) = -1
          spinFlip(3,14) = -1
          spinFlip(4,13) = 1
          spinFlip(5,12) = -1
          spinFlip(6,11) = 1
          spinFlip(7,10) = 1
          spinFlip(8,9)  = -1
          spinFlip(9,8)  = -1
          spinFlip(10,7) = 1
          spinFlip(11,6) = 1
          spinFlip(12,5) = -1
          spinFlip(13,4) = 1
          spinFlip(14,3) = -1
          spinFlip(15,2) = -1
          spinFlip(16,1) = 1
          endif
c         print*,'spinFlip'
c         call outmat(spinFlip,ltonq,ltonq,'spin ',nid)
c         do ii=1,ltonq
c            write(*,"(100g15.5)") ( real(spinFlip(ii,jj)), jj=1,ltonq )
c         enddo
c         print*,ltonq
c         do ii=1,8
c            write(*,"(100g15.5)") ( spinFlip(ii+8,jj+8), jj=1,8 )
c         enddo
c         call exitt
         call zgemm('N','N',ltonq,ltonq,ltonq,(1.0,0.0),
     $        spinFlip,ltonq,rhoVL,ltonq,(0.0,0.0),rhoVR,ltonq)
         call zgemm('N','N',ltonq,ltonq,ltonq,(1.0,0.0),rhoVR,ltonq,
     $        spinFlip,ltonq,(0.0,0.0),rho_svd,ltonq)
         call zgemm('N','N',ltonq,ltonq,ltonq,(1.0,0.0),red_rho,ltonq,
     $        rho_svd,ltonq,(0.0,0.0),rhoFinal,ltonq)
         call zgeev('N','N',ltonq,rhoFinal,ltonq,rhoEval,rhoVL,ltonq,
     $        rhoVR,ltonq,rho_work,2*ltonq,rho_workR,info)
         do ii=2,ltonq
            rhoOrdering = rhoEval(ii)
            jj=ii-1
            do while(jj.ge.1.and.real(rhoEval(jj)).gt.real(rhoOrdering))
               rhoEval(jj+1) = rhoEval(jj)
               jj = jj - 1
            enddo
            rhoEval(jj+1) = rhoOrdering
         enddo

         do ii=1,ltonq
            if(abs(real(rhoEval(ii))).lt.1e-14) rhoEval(ii)=0
         enddo

         concurrence=0
         do ii=1,ltonq
            concurrence = concurrence + sqrt(real(rhoEval(ii)))
         enddo
         concurrence = 2*sqrt(real(rhoEval(ltonq)))-concurrence
         con2 = sqrt(real(rhoEval(4))) - real(sqrt(rhoEval(3))) - 
     $        sqrt(real(rhoEval(2))) - sqrt(real(rhoEval(1)))
c         con2=0
c         do ii=1,numsp
c            tmpC1 = (rho(4*ii-1,4*ii-2) + rho(4*ii-2,4*ii-1))**2
c            tmpC2 = 4*imag(rho(4*ii-2,4*ii-2)-rho(4*ii-1,4*ii-1) + 
c     $           rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2))**2
c            con2 = con2 + sqrt(tmpC1 + tmpC2)
c         enddo
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
     $     concurrence,con2
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

c-------------------------------------------------------------------
      subroutine cem_quantum_bipartite_concurrence(rho_conc,bi_conc,
     $     popSymm)
c-------------------------------------------------------------------
      implicit none
      complex    rho_conc(4,4),rho_svd2(4,4),rho_work2(8),rhoVL2(4,4)
      complex    spinFlip2(4,4),rhoFinal2(4,4),rhoVR2(4,4),rhoEval2(4)
      real       rho_workr2(2*4),bi_conc,rhoOrdering,popSymm(4)
      integer    ii,jj,info

      call czero(spinFlip2,16)
      spinFlip2(1,4) = -1
      spinFlip2(4,1) = -1
      spinFlip2(2,3) = 1
      spinFlip2(3,2) = 1
      
      do ii=1,4
         do jj=1,4
            rhoVL2(ii,jj) = conjg(rho_conc(ii,jj))
         enddo
      enddo
      
      call zgemm('N','N',4,4,4,(1.0,0.0),
     $     spinFlip2,4,rhoVL2,4,(0.0,0.0),rhoVR2,4)
      call zgemm('N','N',4,4,4,(1.0,0.0),rhoVR2,4,
     $     spinFlip2,4,(0.0,0.0),rho_svd2,4)
      call zgemm('N','N',4,4,4,(1.0,0.0),rho_conc,4,
     $     rho_svd2,4,(0.0,0.0),rhoFinal2,4)
      call zgeev('N','N',4,rhoFinal2,4,rhoEval2,rhoVL2,4,
     $     rhoVR2,4,rho_work2,2*4,rho_workR2,info)

      
      do ii=2,4
         rhoOrdering = rhoEval2(ii)
         jj=ii-1
         do while(jj.ge.1.and.real(rhoEval2(jj)).gt.real(rhoOrdering))
            rhoEval2(jj+1) = rhoEval2(jj)
            jj = jj - 1
         enddo
         rhoEval2(jj+1) = rhoOrdering
      enddo
      do ii=1,4
         if(abs(real(rhoEval2(ii))).lt.1e-14) rhoEval2(ii)=0
      enddo

      popSymm(1) = 0.5*(rho_conc(2,2) + rho_conc(3,3) + 
     $     rho_conc(2,3) + rho_conc(3,2))
      popSymm(2) = 0.5*(rho_conc(2,2) + rho_conc(3,3) - 
     $     rho_conc(2,3) - rho_conc(3,2))
      popSymm(3) = 0.5*(rho_conc(1,1) + rho_conc(4,4) + 
     $     rho_conc(1,4) + rho_conc(4,1))
      popSymm(4) = 0.5*(rho_conc(1,1) + rho_conc(4,4) - 
     $     rho_conc(1,4) - rho_conc(4,1))
      


      bi_conc=0
      do ii=1,4
         bi_conc = bi_conc + sqrt(real(rhoEval2(ii)))
      enddo
      bi_conc = 2*sqrt(real(rhoEval2(4)))-bi_conc
c      print*,bi_conc


      return
      end

c--------------------------------------------------------------------
      subroutine cem_quantum_partial_trace_4_to_3(rho_16by16,rho_8by8,
     $     trace_over)
c         Reduce rho even further
c         HARD CODED IN: CONSIDER REVISING
      implicit none
      complex*16 rho_16by16(16,16),rho_8by8(8,8)
      integer trace_over,twotoiq,jjplus,kkplus,jj,kk
      
      if(trace_over.eq.1) then
         twotoiq = 1
         do jj=1,8
            do kk=1,8
c     2*jj-1 is to get the stride right.                                                          
               rho_8by8(jj,kk) =
     $               rho_16by16(2*jj-1,2*kk-1)+
     $               rho_16by16(2*jj-1+twotoiq,2*kk-1+twotoiq)
c     $               +rho_16by16(2*jj-1+twotoiq,2*kk-1)
c     $               +rho_16by16(2*jj-1,2*kk-1+twotoiq)
             enddo
          enddo
       else if(trace_over.eq.3) then
          twotoiq = 4
c        If statements to get the stride correct
          do jj=1,8
             jjplus=4
             if(jj.le.4) jjplus=0
             do kk=1,8
                kkplus=4
                if(kk.le.4) kkplus=0
                rho_8by8(jj,kk) = 
     $               rho_16by16(jj+jjplus,kk+kkplus)+
     $               rho_16by16(jj+twotoiq+jjplus,kk+twotoiq+kkplus)
c     $               +rho_16by16(jj+jjplus,kk+twotoiq+kkplus)
c     $               +rho_16by16(jj+twotoiq+jjplus,kk+kkplus)
             enddo
          enddo
       else if(trace_over.eq.4) then
c        Trace out the highest bit
          twotoiq=8
          do jj=1,8
             do kk=1,8
                rho_8by8(jj,kk) = 
     $               rho_16by16(jj,kk)+
     $               rho_16by16(jj+twotoiq,kk+twotoiq)
c     $               +rho_16by16(jj,kk+twotoiq)
c     $               +rho_16by16(jj+twotoiq,kk)
             enddo
          enddo
       endif
       
      


      return
      end

c--------------------------------------------------------------------
      subroutine cem_quantum_partial_trace_3_to_2(rho_8by8,rho_4by4,
     $     trace_over)
c--------------------------------------------------------------------
c         Reduce rho even further
c         HARD CODED IN: CONSIDER REVISING
c         NEED LOOP STRUCTURE: 11110000, 110011001100,1010101010


      implicit none
      complex*16 rho_8by8(8,8),rho_4by4(4,4)
      integer trace_over,twotoiq,jjplus,kkplus,jj,kk
c     Trace out lowest bit                                                                        
      if(trace_over.eq.1) then
         twotoiq = 1
         do jj=1,4
            do kk=1,4
c     2*jj-1 is to get the stride right.                                                          
               rho_4by4(jj,kk) =
     $               rho_8by8(2*jj-1,2*kk-1)+
     $               rho_8by8(2*jj-1+twotoiq,2*kk-1+twotoiq)
c     $               +rho_8by8(2*jj-1+twotoiq,2*kk-1)
c     $               +rho_8by8(2*jj-1,2*kk-1+twotoiq)
             enddo
          enddo
       else if(trace_over.eq.2) then
c        Trace out the second lowest bit
          twotoiq = 2
c        If statements to get the stride correct
          do jj=1,4
             jjplus=2
             if(jj.le.2) jjplus=0
             do kk=1,4
                kkplus=2
                if(kk.le.2) kkplus=0
                rho_4by4(jj,kk) = 
     $               rho_8by8(jj+jjplus,kk+kkplus)+
     $               rho_8by8(jj+twotoiq+jjplus,kk+twotoiq+kkplus)
c     $               +rho_8by8(jj+jjplus,kk+twotoiq+kkplus)
c     $               +rho_8by8(jj+twotoiq+jjplus,kk+kkplus)
 
            enddo
          enddo
       else if(trace_over.eq.3) then
c        Trace out the third lowest bit
          twotoiq=4
          do jj=1,4
             do kk=1,4
                rho_4by4(jj,kk) = 
     $               rho_8by8(jj,kk)+
     $               rho_8by8(jj+twotoiq,kk+twotoiq)
c     $               +rho_8by8(jj+twotoiq,kk)
c     $               +rho_8by8(jj,kk+twotoiq)
             enddo
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
         sp_mu=sp_mu+imag(rho(ii,ii+level**numqd))
     $        *sqrt(real(isstate(ii+level**numqd)))
      enddo
      sp_mu=-2.0*s_mu*sp_mu

      qd_mu = (0,0)
      do iq=1,numqd
         do ii=1,nstate,level**iq
            do jj=0,level**(iq-1)-1
               qd_mu = qd_mu + imag(rho(ii+jj,ii+jj+level**(iq-1)))
            enddo
         enddo
      enddo
      qd_mu = -2.0*qd_mu*q_mu
      rinc = 4.0*pi*(omega(iEh)/timeunit)/(c_speed*eps_med**0.5)*
     $     d_per_au*1.0d-18/eesu_per_au
      sig   =rinc*(sp_mu+ qd_mu)/energy0/2.0
      rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $     d_per_au*1.0d-18/eesu_per_au
      sigr = rinc * abs((sp_mu+qd_mu)/energy0/2)**2
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
c$$$      open(unit=123,file='matrix')
c$$$      do ii=1,n_in_row*nstate
c$$$         write(123,*)ii,hamilt_0_m_p(ii)+hamilt_m_p(ii)
c$$$      enddo
c      call exitt

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
!$omp parallel do collapse(1)
      do jj = 1,nstate
      do ii = 1,nstate
         tmprr(ii,jj)=0
      enddo
      enddo

      if(idir.eq.1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
!$omp parallel do collapse(2)
      do jj = 1,nstate
      do ii = 1,1
         do kk=1,n_in_row
            tmprr(ii,jj)=tmprr(ii,jj)+
     $           hamilt_I_m_p(n_in_row*(ii-1)+kk)*rrho(ii+kk-1,jj)
         enddo
      enddo
      enddo
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do
!$omp parallel do collapse(1)      
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
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do

      elseif(idir.eq.-1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
!$omp parallel do collapse(2)
      do jj = 1,nstate
      do ii = 1,1
         do kk=1,n_in_row
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(jj,ii+kk-1)
         enddo
      enddo
      enddo
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do
!$omp parallel do collapse(1)
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
!$omp end parallel do

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

