
c---------------------------------------------------------------------
      subroutine cem_quantum2_hamiltonian_setup_pls
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstate_p,binary_rep
      integer  iEh,nsum,npsum,nstate_min
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

      nstate_p   = n_in_row*nstate
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
                 
c     if(ns.eq.nsp) then
c     hamilt_0(ii,jj)=hamilt_0(ii,jj)+omega_q*nq(iq)
c     $                     + omega_s*ns
c     endif
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
         do jj = 1,nstate
            npsum = (jj-1)/ltonq
            do iq=1,numqd
               npsum = npsum + nq_minimal(iq,mod(jj-1,ltonq)+1)
            enddo

            do ii = 1,nstate
               nsum = (ii-1)/ltonq
               do iq=1,numqd
                  nsum = nsum + nq_minimal(iq,mod(ii-1,ltonq)+1)
               enddo
               state_sums(ii,jj)=(nsum-npsum)
            enddo
         enddo
      enddo
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_quantum2_propagation_exp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      integer  iEh, iqstep,ii
      real*8  dclock
      real*8    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      complex  ci,sig

     
      ci = (0.0,1.0)

      call cem_quantum_parameters
      call cem_quantum2_hamiltonian_setup_pls
      stime0= dclock()
      
      
      stime1= dclock()
      call cem_quantum_initialize
      
c     Write file headers 
      write(12,*)'#Time,Re(sp_mu),Im(sp_mu),re(qd_mu),im(qd_mu),
     $re(mu_tot),im(mu_tot),pulse,popSp,popSymm,popAsymm,
     $popQD1,popQD2,popAsymm/popSymm'
      
      do iqstep = 1,nqstep
         stime2= dclock()     
         call cem_quantum2_exp(iEh)
         time= time+dt
         call cem_quantum2_fourier_pls
     $        (s_mu_frequency,q_mu_frequency,energy_t)    
         stime2= dclock()-stime2
      enddo
      stime1= dclock()-stime1
c     print the results out
      write(11,*)'#eV,sig_abs,sig_scatt,sig_tot'
      
      do ii=1,nEh
         rinc = 4.0*pi*(omega(ii)/timeunit)/(c_speed*eps_med**0.5)*
     $        d_per_au*1.0d-18/eesu_per_au
         sig   =rinc*(s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii)
         rinc = 8.0d0*pi*(omega(ii)/timeunit)**4/(3.0d0*c_speed**4)*
     $        d_per_au*1.0d-18/eesu_per_au
         sigr = rinc * cdabs((s_mu_frequency(ii)+q_mu_frequency(ii))/
     $        energy_t(ii))**2
         sigi  = imag(sig)
         tmpr  = omega(ii)*real(s_mu_frequency(ii)+q_mu_frequency(ii))
         tmpi  = omega(ii)*imag(s_mu_frequency(ii)+q_mu_frequency(ii))
         write(11,13) hb_omega(ii),sigi,sigr,sigr+sigi 
      enddo

 13   format(1p5e20.10)

      stime0= dclock()-stime0
      stime3= stime2/neqn

      write(6,6) stime0,stime1,stime2,stime3 
  6   format('stime0/stime1/stime2/stime2_pts:: ',1p4e20.9)
      return
      end


c--------------------------------------------------------------------- 
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

         if (abs(h(j+1,j)).le.eps.and.IFPULSE.eq..false.) then
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

      integer    ii,jj,iEh
      real*8     tmprr,tmpii,pulseFactor,sqrtns
      complex*16 tmpfq,tmpcc
      complex*16 ci
     
      ci=(0.0,1.0)

      do jj = 1,nstate
      do ii = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo

      do ii=1,(2*level**lnumqd-1)*level**lnumqd
         ham_smu_pls(ii) = ham_smu_min(ii)*pulseFactor(1)
         ham_qmu_pls(ii) = ham_qmu_min(ii)*pulseFactor(1)
      enddo         

      call cem_quantum2_mxm(rho_r,tmp_r1, 1)             
      call cem_quantum2_mxm(rho_r,tmp_r2,-1)        
      call cem_quantum2_mxm(rho_i,tmp_i1, 1)             
      call cem_quantum2_mxm(rho_i,tmp_i2,-1)
      

      do jj = 1,nstate
      do ii = 1,nstate
         tmpfq= state_sums(ii,jj)*pulse_omega_0*rho(ii,jj)
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
      subroutine cem_quantum2_mxm(rrho,tmprr,idir)             
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer    ii,jj,idir,kk,mm,ns
      real*8       rrho (lstate,lstate),sqrtns
      real*8       tmprr(lstate,lstate),pulseFactor

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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      do jj = 1,nstate
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
      integer  ii,jj,ns,iq,nsp,kns,knsp,kk,ll,mm
      integer  nqtmp1(lnumqd),nqtmp2(lnumqd)
      integer  nq(lnumqd),nqp(lnumqd)
      integer  nstmp1,nsptmp1
      real*8     fac,tmp,tmp3,tmp1,tmp2
      logical  lbFlag


      do jj = 1,nstate
         nsp= (jj-1)/ltonq
         do iq=1,numqd
            nqp(iq) = nq_minimal(iq,mod(jj-1,ltonq)+1)
         enddo
      do ii = 1,nstate
         ns = (ii-1)/ltonq
         tmp3=gamma_s*(ns+nsp)*0.5
         do iq=1,numqd
            nq(iq) = nq_minimal(iq,mod(ii-1,ltonq)+1)  
            tmp1=gamma_d_m(iq)*(nq(iq)**2+nqp(iq)**2-2*nq(iq)*nqp(iq))
            tmp2=gamma_q_m(iq)*(nq(iq)+nqp(iq))*0.5
c...  QD dephasing Linblad term leads to:
            drho(ii,jj)= drho(ii,jj)-tmp1*rho(ii,jj)
c...  QD spontaneous emission Linblad term leads to these two terms:
            drho(ii,jj)= drho(ii,jj)-tmp2*rho(ii,jj)
            !this may not be correct? I'm not certain.
            !It's certainly not the most efficient way to do this!
           
            if(nqp(iq).eq.0.and.nq(iq).eq.0) then
               drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
     $              rho(ii+level**(iq-1),jj+level**(iq-1))
c               print*,'1'
c               print*,'ij',iqstate_m(iq,ii),iqstate_m(iq,jj)
c               print*,'kl',iqstate_m(iq,ll),iqstate_m(iq,kk)
            endif

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
         


         if(ii+ltonq.lt.nstate.and.
     $        jj+ltonq.lt.nstate) then
            ll = ii+ltonq
            kk = jj+ltonq
            drho(ii,jj) = drho(ii,jj) + gamma_s*sqrt(real((ns+1)*(nsp+1)
     $           ))*rho(ii+ltonq,jj+
     $           ltonq)            
         endif
         
c...  SP spontaneous emission Linblad terms:
         drho(ii,jj)= drho(ii,jj)-tmp3*rho(ii,jj)
         resrho_r(ii,jj)= real(drho(ii,jj))
         resrho_i(ii,jj)= imag(drho(ii,jj))

      enddo
      enddo
      return
      end

c--------------------------------------------------------------------
      subroutine cem_quantum2_fourier_pls(s_mu_freq,q_mu_freq,en_t)
c--------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      integer  iEh,ii,jj,iq,kk
      integer  nq(lnumqd),nqp
      complex*16 ci
      complex*16 tmpCMatrix(lstate,lstate),tmpCMatrix2(lstate,lstate)
      real*8     pulseFactor,fullPulse,currentPop,popSymm,popAsymm
      real*8     tmpC1,tmpC2,concurrence
      real*8     popS,popQD(lnumqd)
      complex  qd_mu,sp_mu
      complex  s_mu_freq(lEh),q_mu_freq(lEh),en_t(lEh)
 
      ci = (0.0,1.0)

      do jj=1,nstate
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
      do ii=1,nstate
         currentPop = real(rho(ii,ii))
         popS = popS + (ii-1)/ltonq*currentPop
         do iq=1,numqd
            popQD(iq) = popQD(iq) + 
     $           nq_minimal(iq,mod(ii-1,ltonq)+1)*currentPop
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
         enddo
c     Find the concurrence; see equation 19 of He and Zhu
         ii=1
         concurrence = 0
         do ii=1,numsp
            tmpC1 = (rho(4*ii-1,4*ii-2) + rho(4*ii-2,4*ii-1))**2
            tmpC2 = 4*imag(rho(4*ii-2,4*ii-2)-rho(4*ii-1,4*ii-1) + 
     $           rho(4*ii-2,4*ii-1) - rho(4*ii-1,4*ii-2))**2
            concurrence = concurrence + sqrt(tmpC1 + tmpC2)
         enddo
      endif


      sp_mu = 0
      do ii=1,nstate-ltonq
         sp_mu=sp_mu+(rho(ii,ii+ltonq)
     $        +rho(ii+ltonq,ii))
     $        *sqrt(real((ii+ltonq-1)/ltonq))
      enddo
      sp_mu=s_mu*sp_mu

      qd_mu = (0,0)
      do iq=1,numqd
         do ii=1,nstate,level**iq
            do jj=0,level**(iq-1)-1
               qd_mu = qd_mu + (rho(ii+jj,ii+jj+level**(iq-1))
     $              + rho(ii+jj+level**(iq-1),ii+jj))*q_mu_m(iq)
     $              *sqrt(real(nq_minimal(iq,mod(ii+jj+level**(iq-1)-1,
     $              ltonq)+1)))
            enddo
         enddo
      enddo

      write(12,12)time*timeunit,sp_mu,qd_mu,qd_mu+sp_mu,
     $     fullPulse(1)*eesu_per_au,
     $     popS,popSymm,popAsymm,popQD(1),popQD(2),popAsymm/popSymm,
     $     concurrence
 12   format(1p20g20.8)
c     Now, do the fourier transform to save the value
      do ii=1,nEh
         en_t(ii) = en_t(ii)
     $        +fullPulse(1)*exp(ci*omega(ii)*(time-pulse_t_0))
         s_mu_freq(ii) = s_mu_freq(ii)
     $        +(exp(ci*omega(ii)*(time-pulse_t_0))*sp_mu)
         q_mu_freq(ii) = q_mu_freq(ii)
     $        +(exp(ci*omega(ii)*(time-pulse_t_0))*qd_mu)

      enddo


      end
