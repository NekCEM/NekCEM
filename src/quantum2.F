c--------------------------------------------------------------------- 
      subroutine cem_quantum2_mxm_m(rrho,tmprr,idir)             
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer    ii,jj,idir,kk    
      real*8       rrho (lstate,lstate/lp)
      real*8       tmprr(lstate,lstate/lp),tmp_acc

c!$acc data present(tmprr)

c!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
c      do jj = 1,nstate/np
c      do ii = 1,nstate
c         tmprr(ii,jj)=0
c      enddo
c      enddo
c!$omp end parallel do
c!$acc end parallel loop
c!$acc end data      
      if(idir.eq.1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist

ccc FIXME Reductions here?
!$acc data present(tmprr,rrho,hamilt_I_m_p)                                                          
!$acc parallel loop collapse(2) 
c!$omp parallel do collapse(1)
      do jj = 1,nstate/np
      do ii = 1,1
         tmp_acc=0
         do kk=1,n_in_row
            tmp_acc=tmp_acc+
     $           hamilt_I_m_p(n_in_row*(ii-1)+kk)*rrho(ii+kk-1,jj)
         enddo
         tmprr(ii,jj)=tmp_acc
         
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
      do jj = 1,nstate/np
      do ii = 2,n_in_row-1
         tmp_acc=0
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,n_in_row
            tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(ii+kk-1,jj)
         enddo
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

         do kk=1,ii-1
            tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
         tmprr(ii,jj)=tmp_acc
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
c     This is the bulk of the matrix
c!$acc parallel loop gang private(tmp_acc)
      do jj = 1,nstate/np
c!$acc loop worker private(tmp_acc)
      do ii = n_in_row,nstate-(n_in_row-1)
         tmp_acc=0
c     Right of the diagonal; good for any row < nstate-n_in_row-1
c!$acc loop vector reduction(+:tmp_acc)
         do kk=1,n_in_row-1
            tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(ii+kk-1,jj)+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
         
         tmprr(ii,jj) = tmp_acc+hamilt_I_m_p(n_in_row*ii)
     $           *rrho(ii+n_in_row-1,jj)

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1

c         do kk=1,n_in_row-1
c           tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)-
c     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
c         enddo
c         tmprr(ii,jj)=tmp_acc

      enddo
      enddo
c!$omp end parallel do
c!$acc end parallel loop

!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
c     The end parts only have the up part, not the right part.
      do jj = 1,nstate/np
      do ii = nstate-(n_in_row-1)+1,nstate-1
         tmp_acc=0
c     Right of the diagonal; good for any row < nstate-n_in_row-1
         do kk=1,nstate-(ii-1)
            tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(ii+kk-1,jj)
         enddo

c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
         tmprr(ii,jj)=tmp_acc
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop
!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
c     Last element has only left of the diagonal terms
      do jj = 1,nstate/np
      do ii = nstate,nstate
         tmp_acc=0
c     Left of the diagonal; instead of incrementing to the right, we
c     increment 'up', because  hamiltonian is symmetric
c     good for any row >n_in_row-1
         do kk=1,n_in_row-1
            tmp_acc=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)-
     $           kk*(n_in_row-1)+1)*rrho(ii-kk,jj)
         enddo
         tmprr(ii,jj)=tmp_acc+hamilt_I_m_p(n_in_row*(ii-1)+1)
     $        *rrho(ii,jj)
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

!$acc end data

      elseif(idir.eq.-1) then
c     The beginning n_in_row-1 rows need special treatment; 
c     only the diagonal and the terms to the right exist
      do jj = 1,nstate/np
      do ii = 1,1
         do kk=1,n_in_row
            tmprr(jj,ii)=tmprr(jj,ii)+hamilt_I_m_p(n_in_row*(ii-1)+kk)
     $           *rrho(jj,ii+kk-1)
         enddo
         
      enddo
      enddo
      do jj = 1,nstate/np
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
      do jj = 1,nstate/np
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
      do jj = 1,nstate/np
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
      do jj = 1,nstate/np
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
      
      subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)

c*********************************************************************72
c
cc DNSCSR converts Dense to Compressed Row Sparse format.
c
c converts a densely stored matrix into a row orientied
c compactly sparse matrix. ( reverse of csrdns )
c Note: this routine does not check whether an element
c is small. It considers that a(i,j) is zero if it is exactly
c equal to zero: see test below.
c
c on entry:
c
c
c nrow      = row-dimension of a
c ncol      = column dimension of a
c nzmax = maximum number of nonzero elements allowed. This
c         should be set to be the lengths of the arrays a and ja.
c dns   = input nrow x ncol (dense) matrix.
c ndns      = first dimension of dns.
c
c on return:
c
c
c a, ja, ia = value, column, pointer  arrays for output matrix
c
c ierr      = integer error indicator:
c         ierr .eq. 0 means normal retur
c         ierr .eq. i means that the the code stopped while
c         processing row number i, because there was no space left in
c         a, and ja (as defined by parameter nzmax).
c
      double precision dns(ndns,*),a(*)
      integer ia(*),ja(*)

      ierr = 0
      next = 1
      ia(1) = 1
      do 4 i=1,nrow
         do 3 j=1,ncol
            if (dns(i,j) .eq. 0.0) goto 3
            if (next .gt. nzmax) then
               ierr = i
               return
            end if
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
 3                 continue
         ia(i+1) = next
 4           continue
      return
      end

      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr)

c*********************************************************************72
c
cc CSRDNS converts Compressed Sparse Row to Dense format.
c
c converts a row-stored sparse matrix into a densely stored one
c
c On entry:
c
c
c nrow      = row-dimension of a
c ncol      = column dimension of a
c a,
c ja,
c ia    = input matrix in compressed sparse row format.
c         (a=value array, ja=column array, ia=pointer array)
c dns   = array where to store dense matrix
c ndns      = first dimension of array dns
c
c on return:
c
c dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
c
c ierr  = integer error indicator.
c         ierr .eq. 0  means normal return
c         ierr .eq. i  means that the code has stopped when processing
c         row number i, because it found a column number .gt. ncol.
c
      double precision dns(ndns,*),a(*)
      integer ja(*),ia(*)
      
      ierr = 0
      do 1 i=1, nrow
         do 2 j=1,ncol
          dns(i,j) = 0.0
 2     continue
 1    continue
      
      do 4 i=1,nrow
         do 3 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (j .gt. ncol) then
               ierr = i
               return
            end if
            dns(i,j) = a(k)
 3       continue
 4    continue
      return
      end
      
      
c-------------------------------------------
c Turn the banded hamiltonian to csr format
c Based off of SparseKit, http://people.sc.fsu.edu/~jburkardt/f77_src/sparsekit/sparsekit.f
c-------------------------------------------
c--------------------------------------------------------------------
      subroutine cem_quantum2_hamiltonian_to_csr
c--------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      
      integer ii,jj,kk

c     FIXME: 12/11/14: If numqd=1, this breaks! It overwrite nEh
      
c     First row
      ham_nnz=1
      ham_ia(1)=1
      do ii=1,1
         do kk=1,n_in_row
             if(hamilt_m_p(n_in_row*(ii-1)+kk)+
     $           hamilt_0_m_p(n_in_row*(ii-1)+kk).ne.0.0) then
               ham_ja(ham_nnz) = ii+kk-1
               ham_m_p_a(ham_nnz) = hamilt_m_p(n_in_row*(ii-1)+kk)
               ham_0_m_p_a(ham_nnz)  = hamilt_0_m_p(n_in_row*(ii-1)+kk)
               ham_nnz = ham_nnz+1
            endif
         enddo
         ham_ia(ii+1) = ham_nnz
      enddo

c     Next few rows
      do ii=2,n_in_row-1
         do kk=1,n_in_row
             if(hamilt_m_p(n_in_row*(ii-1)+kk)+
     $           hamilt_0_m_p(n_in_row*(ii-1)+kk).ne.0.0) then
               ham_ja(ham_nnz) = ii+kk-1
               ham_m_p_a(ham_nnz) = hamilt_m_p(n_in_row*(ii-1)+kk)
               ham_0_m_p_a(ham_nnz)  = hamilt_0_m_p(n_in_row*(ii-1)+kk)
               ham_nnz = ham_nnz+1
            endif
         enddo
         do kk=1,ii-1
            if(hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)+
     $           hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1).ne.0.0) 
     $           then
               ham_ja(ham_nnz) = ii-kk
               ham_m_p_a(ham_nnz) = 
     $              hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_0_m_p_a(ham_nnz)  = 
     $              hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_nnz = ham_nnz+1
            endif
         enddo
         ham_ia(ii+1)=ham_nnz
      enddo

c     Bulk
      do ii = n_in_row,nstate-(n_in_row-1)
         do kk = 1,n_in_row
             if(hamilt_m_p(n_in_row*(ii-1)+kk)+
     $           hamilt_0_m_p(n_in_row*(ii-1)+kk).ne.0.0) then
               ham_ja(ham_nnz) = ii+kk-1
               ham_m_p_a(ham_nnz) = hamilt_m_p(n_in_row*(ii-1)+kk)
               ham_0_m_p_a(ham_nnz)  = hamilt_0_m_p(n_in_row*(ii-1)+kk)
               ham_nnz = ham_nnz+1
            endif
         enddo
         do kk=1,n_in_row-1
            if(hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)+
     $           hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1).ne.0.0) 
     $           then
               ham_ja(ham_nnz) = ii-kk
               ham_m_p_a(ham_nnz) = 
     $              hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_0_m_p_a(ham_nnz)  = 
     $              hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_nnz = ham_nnz+1
            endif
         enddo
         ham_ia(ii+1)=ham_nnz
      enddo

c last few rows
      do ii=nstate-(n_in_row-1)+1,nstate-1
         do kk=1,nstate-(ii-1)
             if(hamilt_m_p(n_in_row*(ii-1)+kk)+
     $           hamilt_0_m_p(n_in_row*(ii-1)+kk).ne.0.0) then
               ham_ja(ham_nnz) = ii+kk-1
               ham_m_p_a(ham_nnz) = hamilt_m_p(n_in_row*(ii-1)+kk)
               ham_0_m_p_a(ham_nnz)  = hamilt_0_m_p(n_in_row*(ii-1)+kk)
               ham_nnz = ham_nnz+1
            endif
         enddo
         do kk=1,n_in_row-1
            if(hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)+
     $           hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1).ne.0.0) 
     $           then
               ham_ja(ham_nnz) = ii-kk
               ham_m_p_a(ham_nnz) = 
     $              hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_0_m_p_a(ham_nnz)  = 
     $              hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_nnz = ham_nnz+1
            endif
         enddo
         ham_ia(ii+1)=ham_nnz
      enddo

c     last row
      do ii=nstate,nstate
         do kk=1,n_in_row-1
            if(hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)+
     $           hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1).ne.0.0) 
     $           then
               ham_ja(ham_nnz) = ii-kk
               ham_m_p_a(ham_nnz) = 
     $              hamilt_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_0_m_p_a(ham_nnz)  = 
     $              hamilt_0_m_p(n_in_row*(ii-1)-kk*(n_in_row-1)+1)
               ham_nnz = ham_nnz+1
            endif
         enddo
         ham_ja(ham_nnz)=ii
         ham_m_p_a(ham_nnz) = 
     $        hamilt_m_p(n_in_row*(ii-1)+1)
         ham_0_m_p_a(ham_nnz)  = 
     $        hamilt_0_m_p(n_in_row*(ii-1)+1)
         ham_nnz = ham_nnz+1
         
         ham_ia(ii+1)=ham_nnz
      enddo

      return
      end


c-------------------------------------------
c csr mxm
c Based off of SparseKit, http://people.sc.fsu.edu/~jburkardt/f77_src/sparsekit/sparsekit.f
c and Mills Vectorized sparse mxm for csr
c-------------------------------------------
c---------------------------------------------------------------------
      subroutine cem_quantum2_csr_mxm(rrho,tmprr)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer ii,jj,kk,istart,iend
      real    sum
      real*8       rrho (lstate,lstate/lp)
      real*8       tmprr(lstate,lstate/lp),tmp_acc

!$acc data present(ham_ia,ham_ja,ham_a,rrho,tmprr)
!$acc parallel loop private(istart,iend)
      do jj=1,nstate/np
!$acc loop vector private(sum)
      do ii=1,nstate
         istart = ham_ia(ii)
         iend   = ham_ia(ii+1)-1
         sum = 0
         do kk=istart,iend
            sum = sum + ham_a(kk)*rrho(ham_ja(kk),jj)
         enddo
         tmprr(ii,jj) = sum
      enddo
      enddo
!$acc end data
      return
      end



c---------------------------------------------------------------------
      subroutine cem_quantum2_hamiltonian_setup_pls_f
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer  iq,is,ii,jj,kk,iiq
      integer  ns,nsp,nstate_ham
      integer  nq(lnumqd),nqp(lnumqd)
      integer  binary_rep
      integer  iEh,nsum,npsum
      real*8   tmp,tmpi
      real*8   hamilt_q_mu(n_in_row*lstate)
      real*8   hamilt_0_pls(n_in_row*lstate)
      logical  q_muFlag,s_muFlag

 

c!$acc data present(hamilt_0_m_p,hamilt_m_p,nq_minimal)
c!$acc&     pcreate(hamilt_q_mu,hamilt_0_pls)
c!$acc&     present(omega_q_m,gamma_q_m,gamma_d_m,q_mu_m,g_couple_m)
c!$acc&     present(isstate,iqstate_m)
c!$acc&     present(state_sums)

ccc FIXME can't accelerate because call to binary rep?
      do ii=1,2*ltonq
         do kk=1,numqd
            nq_minimal(kk,ii)=binary_rep(kk,mod(ii-1,ltonq)+1)
         enddo
      enddo


c!$acc update device(nq_minimal)

      kk=0
c     Assign isstate,iqstate
ccc FIXME Have to update device? doesnt support function calls (binary_rep) 
      do is=1,numsp !number of qunatum states of surface plasmon
         do iq=1,level**numqd !total number of quantum states for all qds
            kk=kk+1
c            kk=(is-1)*level**numqd+iq
            isstate(kk) = is-1
c           iqstate, since they are all assumed to be the same number of levels,
c           is calculated by taking the (if 2) binary representation of the number
c           and assigning its bit to the corresponding quantum dot
            do ii=1,numqd
               iqstate_m(ii,kk) = binary_rep(ii,iq)
            enddo
         enddo
      enddo
 
c!$acc update device(isstate,iqstate_m)

c
c     Initially, generate upper half of s=0 part of the hamiltonian in normal way

c     the number needed to store in a single row is
c     3*2**Nqd - (2**Nqd + 1)

      nstate_ham = n_in_row*nstate
ccc   Declare here?
c#ifdef _OPENACC
c      call rzero_acc(hamilt_0_m_p,nstate_ham)
c      call rzero_acc(hamilt_I_m_p,nstate_ham)
c      call rzero_acc(hamilt_m_p,nstate_ham)
c      call rzero_acc(hamilt_0_pls,nstate_ham)
c      call rzero_acc(hamilt_q_mu,nstate_ham)
c#else
      call rzero(hamilt_0_m_p,nstate_ham)
      call rzero(hamilt_I_m_p,nstate_ham)
      call rzero(hamilt_q_mu,nstate_ham)
      call rzero(hamilt_m_p,nstate_ham)
      call rzero(hamilt_0_pls,nstate_ham)
c#endif


c     2*level**numberOfQuantumd dots is because we want to all q possibilities
c     for the values of 0 and 1 for s

c!$acc parallel loop collapse(2)
      do ii=1,level**numqd
         do jj=1,2*level**numqd
            ns=isstate(ii)
            do iq=1,numqd
               nq(iq)=iqstate_m(iq,ii)
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
c!$acc end parallel loop
      
c     Now that we have the upper half of the first part of the matrix, we now propagate it through the whole matrix
c     starting after the last one we assigned, all the way to the last one to assign
c!$acc parallel loop
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
c!$acc end parallel loop

c...  set the diagonal terms
c!$acc parallel loop
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
c!$acc end parallel loop

      
c      open(unit=123,file='matrix2')
c      do ii=1,n_in_row*nstate
c         write(123,*)ii,hamilt_0_m_p(ii)+hamilt_m_p(ii)
c      enddo
c      call exitt
c...  define the frequencies for different quantum states at each energy level iEh.
ccc This data is only for the host?
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
      enddo


c!$acc parallel loop
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
c!$acc end parallel loop
c!$acc end data

c      open(unit=123,file='matrix2')
c      do ii=1,n_in_row*nstate
c         write(123,*)ii,hamilt_0_m_p(ii)+hamilt_m_p(ii)
c      enddo
c      call exitt

      
      return
      end


c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian_ndot_dense
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'

      integer  iq,is,ii,jj,kk,iiq,ll
      integer  ns,nsp
      integer  nq(lnumqd),nqp(lnumqd),nq2,nqp2
      integer  nstate2,binary_rep
      integer  ieh,nsum,npsum
      real     tmp,tmpi,hamilt_0(nstate,nstate),hamilt(nstate,nstate)
      real     hamilt_I(nstate,nstate),ham(nstate,nstate)
      real     ham_tmp((2*level**lnumqd-1)*lstate)
      logical  q_muFlag,s_muFlag


c     Assign isstate,iqstate
      kk=0
      do iq=1,level**numqd      !total number of quantum states for all qds
         kk=kk+1
c     iqstate, since they are all assumed to be the same number of levels,
c     is calculated by taking the n-nary representation of the number
c     and assigning its bit to the corresponding quantum dot
         do ii=1,numqd
            iqstate_m(ii,kk) = binary_rep(ii,iq)
         enddo
      enddo
      
      nstate  = level**numqd
      nstate2 = nstate*nstate
      call rzero(hamilt_0,nstate2)
      call rzero(hamilt  ,nstate2)
      call rzero(hamilt_I,nstate2)

c...  set the diagonal terms
      do  ii = 1,nstate 
c         Add the contribution from each quantum dot
          do kk = 1,numqd
             nq(kk) = iqstate_m(kk,ii) 
             hamilt_0(ii,ii) = hamilt_0(ii,ii)+nq(kk)*omega_q_m(kk)
          enddo
      enddo

c     loop to do off diagonal terms
c     the nifty 1,3,5 structure isnt as useful in the case of multiple qdots
      do ii=1,nstate
         do iq=1,numqd
            nq(iq)  = iqstate_m(iq,ii)
         enddo
         do jj=1,nstate
            !this loop may be redundant with the next
            do iq=1,numqd
               nqp(iq) = iqstate_m(iq,jj)
            enddo
c           Assign the hamiltonian
            do iq=1,numqd
               q_muFlag=.true.
               do kk=1,numqd
c               We now have to loop over the nq to check if
c               all of the other states are equal
                  if(iq.ne.kk.and.nq(kk).ne.nqp(kk)) then
                     q_muFlag=.false.
                  endif
               enddo
               if(q_muFlag.eqv..true.) then
                     if(nqp(iq).eq.nq(iq)+1) then 
                        hamilt(ii,jj) = hamilt(ii,jj)
     $                       -dsqrt(dfloat(nq(iq)+1))*q_mu_m(iq)
                      endif
                     if(nqp(iq).eq.nq(iq)-1) then
                        hamilt(ii,jj) = hamilt(ii,jj)
     $                       -dsqrt(dfloat(nq(iq)))*q_mu_m(iq)
                     endif             
                  endif
               do kk=iq+1,numqd
c              now do hamilt parts
                  q_muFlag=.true.
c                 Loop through to check to see if the other quantum dots are
c                 in different states; if they are, we don't include the interaction
c                 term in this element
                  do ll=1,numqd
                     if(ll.ne.kk.and.ll.ne.iq.and.nq(ll).ne.nqp(ll))then
                        q_muFlag=.false.
                     endif
                  enddo
                  if(q_muFlag.eqv..true.) then

                     nq2  = iqstate_m(kk,ii)
                     nqp2 = iqstate_m(kk,jj)

c... Coupling between dot ii and jj is described g_ndot(ii,jj)
                     if(nqp(iq).eq.nq(iq)-1.and.nqp2.eq.nq2+1) then
                        hamilt_0(ii,jj) = hamilt_0(ii,jj)
     $                       -g_ndot(iq,kk)*
     $                       dsqrt(dfloat((nq2+1)*nq(iq)))
                     endif
                     
                     if(nqp(iq).eq.nq(iq)+1.and.nqp2.eq.nq2-1) then
                        hamilt_0(ii,jj) = hamilt_0(ii,jj)
     $                       -g_ndot(iq,kk)*
     $                       dsqrt(dfloat(nq2*(nq(iq)+1)))
                     endif
                  endif
               enddo
            enddo
         enddo
      enddo
     
      !We have to temporarily construct the full matrix to get the numbering right
      do jj=1,nstate
         do ii=1,nstate
            ham(ii,jj) = hamilt_0(ii,jj) + hamilt(ii,jj)
         enddo
      enddo

      !Get csr format of the full matrix
      call dnscsr(nstate,nstate,n_in_row*nstate,ham,nstate,
     $     ham_tmp,ham_ja,ham_ia,ii)

      ham_nnz = 1
      !Subtract off parts to get the time independent and time dependent versions
      do ii=1,nstate
         do kk=ham_ia(ii),ham_ia(ii+1)-1
            jj              = ham_ja(kk)
            ham_0_m_p_a(kk) = ham_tmp(kk) - hamilt(ii,jj)
            ham_m_p_a(kk)   = ham_tmp(kk) - hamilt_0(ii,jj)
            ham_nnz         = ham_nnz + 1
         enddo
      enddo

c$$$      call csrdns(nstate,nstate,ham_0_m_p_a,ham_ja,ham_ia,ham,nstate,ii)
c$$$      call outmat(ham,nstate,nstate,'ham_0     ',nid)
c$$$      call csrdns(nstate,nstate,ham_m_p_a,ham_ja,ham_ia,ham,nstate,ii)
c$$$      call outmat(ham,nstate,nstate,'ham       ',nid)
c$$$
c$$$      call exitt()

c...  define the frequencies for different quantum states at each energy level iEh.
      do iEh= 1,nEh   
         omega(iEh)   = Emin+iEh*dEh
         hb_omega(iEh)= omega(iEh)   
         omega(iEh)   = omega(iEh)/eV   
         do jj = 1,nstate
            npsum = 0
            do iq=1,numqd
               npsum = npsum + iqstate_m(iq,jj)
            enddo

            do ii = 1,nstate
               nsum = 0
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

      nstate2         = nstate*nstate/np
      maxPartialTrace = level**numqd

      call czero(s_mu_frequency,nEh)
      call czero(q_mu_frequency,nEh)
      call czero(energy_t,      nEh)
#ifdef GPU
#ifdef _OPENACC
      call rzero_acc(   rho_r,nstate2)
      call rzero_acc(   rho_i,nstate2)
      call rzero_acc( k_rho_r,nstate2)
      call rzero_acc( k_rho_i,nstate2)
      call rzero_acc(resrho_r,nstate2)
      call rzero_acc(resrho_i,nstate2)
#endif
#else
      call rzero(   rho_r,nstate2)
      call rzero(   rho_i,nstate2)
      call rzero( k_rho_r,nstate2)
      call rzero( k_rho_i,nstate2)
      call rzero(resrho_r,nstate2)
      call rzero(resrho_i,nstate2)

#endif
      call rzero(rho_r,nstate2)


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


c!$acc parallel
c      if(nid.eq.0) rho_r(1,1) = 1
c!$acc end parallel

      if(nid.eq.0) rho_r(1,1) = 1
!$acc update device(rho_r)
c!$acc update device(rho_i)
c!$acc update device(k_rho_r)
c!$acc update device(k_rho_i)
c!$acc update device(resrho_r)
c!$acc update device(resrho_i)

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
      integer ii,jj,kk,IERR,gloJ,gloI,locJ
      integer*8 size_type
      

c     Setup MPI_GET things; 8 is size of real
c      call MPI_WIN_CREATE(tmp_r1,nstate*nstate/np,8,MPI_INFO_NULL,
c     $     MPI_COMM_WORLD,mpi_win_r,IERR)
c      call MPI_WIN_CREATE(tmp_i1,nstate*nstate/np,8,MPI_INFO_NULL,
c     $     MPI_COMM_WORLD,mpi_win_i,IERR)

c!$acc enter data create(tmp_c1)
c!$acc data present(tmp_c1)
c!$acc host_data use_device(tmp_c1)
      size_type = nstate*nstate/np
      call MPI_WIN_CREATE(tmp_c1,size_type,16,MPI_INFO_NULL,
     $     MPI_COMM_WORLD,mpi_win_c,IERR)
c!$acc end host_data
c!$acc end data
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
      call MPI_BCAST(eV,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
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
      call MPI_BCAST(IFSMALLMEMORY,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

      nstate_p = nstate*nstate/np
      c_offset = nid*nstate/np
      dt      = timemax/(nqstep-1)
      printdt = timemax/(999)             !We only need to print 1000 points
      time    = 0
      timeOld = -1

      ltonq = level**numqd
      maxPartialTrace=level**numqd

ccc This doesn't need to be done on the GPU, because it is for host to communicate
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

      do ii=1,numsp
         do jj=1,ltonq
            do kk=1,ltonq
               gloJ = ltonq*ii-(ltonq-jj)
               gloI = ltonq*ii-(ltonq-kk)

               if((gloJ.gt.nid*nstate/np).and.
     $              (gloJ.le.(nid+1)*nstate/np)) then
                  locJ = gloJ - c_offset
                  glo_num_q_concV(gloI+(locJ-1)*nstate) =
     $                 kk+(jj-1)*maxPartialTrace
               endif

c$$$               red_rho(kk,jj) = red_rho(kk,jj) + 
c$$$     $              rho(ltonq*ii-(ltonq-kk),ltonq*ii-(ltonq-jj))

            enddo
         enddo
      enddo
      
c$$$      open(unit=123,file='matrix')
c$$$      do ii=1,nstate*nstate
c$$$         gloJ = ii/nstate
c$$$         gloI = mod(ii,nstate)
c$$$         write(123,*)ii,gloI,gloJ,glo_num_q_concV(ii)
c$$$      enddo
c$$$      call outmat(k_rho_r,nstate,nstate,'glo_num ',nid)

      call gs_setup(gs_handle_q_H,glo_num_q_Hv,nstate_p,
     $     MPI_COMM_WORLD,np)
      call gs_setup(gs_handle_q_conc,glo_num_q_concV,nstate_p,
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
      common /QRDATA_R5/  r_pulse_duration,r_pulse_fluence,r_gcouple_q3
      common /QLDATA_L1/  i_IFSMALLMEMORY
      integer             i_numqd,i_numsp,i_Eh,i_nqstep,ii,jj
      real                r_Eh1,r_Eh2,r_Eh,r_tmax,r_nqstep,er1,er2
      real                rintense,r_omega_q1,r_omega_q2,r_gcouple_q1
      real                r_gcouple_q2,r_gamma_q1,r_gamma_q2,r_gamma_d1
      real                r_gamma_d2,r_pulse_t_0,r_q_mu_1,r_q_mu_2
      real                r_pulse_duration,r_pulse_fluence,r_gcouple_q3
      logical             i_IFSMALLMEMORY
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
      printdt = timemax/(999)
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
      if(numsp.ne.1) then
         omega_s   = 2.05d0/eV
         gamma_s   = 1.5d-1/eV
         s_mu      = 4.d3*debye
      else
         omega_s   = 0
         gamma_s   = 0
         s_mu      = 0
      endif

c...  QD: quantum dot
      omega_q   = 2.05d0/eV       
      gamma_q   = (1.9e-7)/eV
      gamma_d   = (2.0e-3)/eV
      q_mu      = 1.3d1*debye  

      IFSMALLMEMORY=i_IFSMALLMEMORY

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
      if(numqd.gt.2) then
         g_couple_m(4) = r_gcouple_q3/eV
         g_couple_m(3) = r_gcouple_q3/ev
         g_couple_m(2) = r_gcouple_q2/ev
         g_couple_m(1) = r_gcouple_q1/ev
      endif
      if (nid.eq.0) then 
         write(6,*) 'g_couple = ',( g_couple_m(ii)*eV, ii=1,numqd )
         write(14,1234) '#g=',( g_couple_m(ii)*eV, ii=1,numqd )
         write(13,1234) '#g=',( g_couple_m(ii)*eV, ii=1,numqd )
         write(12,1234) '#g=',( g_couple_m(ii)*eV, ii=1,numqd )

 1234 format(' ',A3,10E9.2)
      endif

c...  radius in meters
      qradius   = 20e-9                  ! in meters   
      qradius   = qradius/(0.529177e-10) ! atomic units


c... Coupling between dot ii and jj is described g_ndot(ii,jj)
      do ii = 1,numqd
         do jj = ii+1,numqd
            !Assuming the dots are all in a line, their distance is
            !equal to the difference of their labels
           g_ndot(ii,jj) =  2.465/((jj-ii)*qradius)**3
         enddo
      enddo

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
      include 'RK5'
      include 'mpif.h'
      integer  iEh, ii, IERR, jj
      real*8     stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7
      real*8     sigr,sigi,tmpr,tmpi,rinc,sigcon,pulseFactor
      real*8   dclock
      complex  sig

      if(nid.eq.0) then
         call cem_quantum2_parameters

c     Write file headers 
         write(12,*)'#Time,sp_mu,qd_mu(1),qd_mu(2),
     $mu_tot,pulse,popSp'
         if(numqd.eq.4) write(13,*)'#Time,4:3,4:2,3:2,4:1,2:1,3:1'
         if(numqd.eq.3) write(13,*)'#Time,3:2,3:1,2:1'
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
      call cem_quantum_rk_storage
c      call cem_quantum_hamiltonian_ndot_dense
c    *_setup_pls does the minimum version (folded up)
c    *_setup_pls_f does the full, banded version
      call cem_quantum2_hamiltonian_setup_pls
      call cem_quantum2_hamiltonian_setup_pls_f
      call cem_quantum2_hamiltonian_to_csr

ccc Setup acc stuff

!$acc data copyin(hamilt_I_m_p,hamilt_0_m_p,hamilt_m_p)
!$acc&     copyin(ham_a,ham_ia,ham_ja,ham_0_m_p_a,ham_m_p_a)
!$acc&     pcreate(rho_i,rho_r,rho,drho,resrho_r,resrho_i)
!$acc&     copyin(isstate,iqstate_m,state_sums,nq_minimal)
!$acc&     pcreate(tmp_r1,tmp_i1,tmp_c1,k_rho_r,k_rho_i)
!$acc&     copyin(rk4a,rk4b,rk4c)
ccc parameters for the GPU to use - scalars automatically imported
c!$acc&     copyin(numqd,numsp,nstate,n_in_row,nEh,nqstep,timemax)
c!$acc&     copyin(Emin,Emax,dEh,relerr,abserr)
c!$acc&     copyin(omega_s,gamma_s,s_mu)
!$acc&     copyin(omega_q_m,gamma_q_m,gamma_d_m,q_mu_m,g_couple_m)
c!$acc&     copyin(pulse_omega_0,pulse_t_0,pulse_duration,energy0_pls)
c!$acc&     copyin(nstate_p,c_offset,dt,nid,ltonq,IFSMALLMEMORY)

      call cem_quantum2_initialize

c      Pulsed case
c      do iqstep= 1,nqstep
       iqstep=1
       do while(time.lt.timemax)
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
c         if(abs(time-pulse_t_0).lt.3*pulse_duration) then
c             dt=dt_small
c         else
c             dt=dt_large
c         endif
         time= time+dt
         stime7= dclock()
         if(time.gt.timeOld) then
            call cem_quantum2_fourier_pls
     $           (s_mu_frequency,q_mu_frequency,energy_t)
            timeOld = timeOld+printdt
            stime7= dclock()-stime7
         endif
         stime1= dclock()-stime1
         qrtime1= qrtime1+stime1
         qrtime2= stime4
         qrtime3= stime5
         qrtime4= stime7
         call userchk
         iqstep=iqstep+1
      enddo
      stime0 = dclock()-stime0
      if(nid.eq.0) write(6,*)'stime0:',stime0 
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
      close(12)
      close(11)
      endif

!$acc end data

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

      integer    ii,jj,iEh,ierr,i_one,i_zero,i0,j0
      real*8     tmprr,tmpii,pulseFactor_loc,sqrtns
      real*8     tmprv(lstate*lstate/lp),alpha
      real*8     tmpiv(lstate*lstate/lp),pulseFactor
      real*8     tmppp(2*lstate*lstate/lp)
      complex*16 tmpfq,tmpcc
      complex*16 ci,tmp_gs_op(nstate*nstate/np)
      real*8     dclock,ttmp

      i_one = 1
      i_zero = 0
      ci=(0.0,1.0)


!$acc data pcreate(tmprv,tmpiv)
!$acc&     present(rho_r,rho_i,tmp_r1,tmp_i1,resrho_i,resrho_r,drho)
!$acc&     present (hamilt_0_m_p,hamilt_I_m_p,hamilt_m_p,rho,state_sums)
!$acc&     present(ham_a,ham_0_m_p_a,ham_m_p_a)

c     Ensure Hermiticity of rho
!$acc parallel loop collapse(2) private(j0)
c!$omp parallel do collapse(2)
      do jj = 1,nstate/np
      do ii = 1,nstate
         j0= (jj-1)*nstate+ii
         tmprv(j0) = rho_r(ii,jj)
         tmpiv(j0) = rho_i(ii,jj)
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

      ttmp = dclock()

      call nekgsync()!#FIXME WHY IS THIS NECESSARY?
      call gs_op_fields(gs_handle_q_H,tmprv,
     $     nstate_p,1,1,1,0) 
      call gs_op_fields(gs_handle_q_H,tmpiv,nstate_p,1,1,1,0)
      qctime1= qctime1+dclock()-ttmp     

!$acc parallel loop collapse(2) private(j0)
c!$omp parallel do collapse(1)
      do jj=1,nstate/np
      do ii=1,nstate
         j0= (jj-1)*nstate+ii
         if(ii.ne.jj+c_offset) rho_r(ii,jj) = tmprv(j0)/2
         rho_i(ii,jj) = (2*rho_i(ii,jj) - 
     $        tmpiv(j0))/2
         if(jj+c_offset.eq.ii) rho_i(ii,jj) = 0
c         rho_r(ii,jj) = tmprv(j0)
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

c!$acc update host(rho_r)
c      call outmat(rho_r,nstate,nstate/4,'rho_r  ',nid)
c      call exitt
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

!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1) 
      do jj = 1,nstate/np
      do ii = 1,nstate
         rho(ii,jj)= rho_r(ii,jj)+ci*rho_i(ii,jj)      
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

ccc TEMP!!! Disabled if(IFSMALLMEMORY) so that compiler would compile acc

c      if(IFSMALLMEMORY) then
c!$acc parallel loop
c!$omp parallel do 
c         do ii=1,(2*level**lnumqd-1)*level**lnumqd
c            ham_smu_pls(ii) = ham_smu_min(ii)*pulseFactor(1)
c            ham_qmu_pls(ii) = ham_qmu_min(ii)*pulseFactor(1)
c         enddo
c!$omp end parallel do
c!$acc end parallel loop
c         call cem_quantum2_mxm(rho_r,tmp_r1, 1)             
c         call cem_quantum2_mxm(rho_i,tmp_i1, 1)             
c      else

      alpha = 2. * log(2.0) / pulse_duration**2

      pulseFactor_loc = energy0_pls/2.0*exp(-alpha*(time-pulse_t_0)**2)
      pulseFactor_loc = pulseFactor(1)
!$acc parallel loop
c!$omp parallel do
         do ii=1,nstate*n_in_row
            hamilt_I_m_p(ii) = hamilt_0_m_p(ii) + 
     $           pulseFactor_loc*hamilt_m_p(ii)
         enddo

c      call exitt

c!$omp end parallel do
!$acc end parallel loop

!$acc parallel loop
      do ii=1,ham_nnz
         ham_a(ii) = ham_0_m_p_a(ii) + pulseFactor_loc*ham_m_p_a(ii)
      enddo
!$acc end parallel loop


c!$acc update host(hamilt_I_m_p)
c         call cem_quantum2_mxm_m(rho_r,tmp_r1,1)
c         call cem_quantum2_mxm_m(rho_i,tmp_i1,1)
         call cem_quantum2_csr_mxm(rho_r,tmp_r1)
         call cem_quantum2_csr_mxm(rho_i,tmp_i1)
c     endif

!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
      do jj = 1,nstate/np
      do ii = 1,nstate
         tmprv((jj-1)*nstate+ii) = tmp_r1(ii,jj)
         tmpiv((jj-1)*nstate+ii) = tmp_i1(ii,jj)
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

c      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      ttmp= dclock()
      call gs_op_fields(gs_handle_q_H,tmprv,
     $     nstate_p,1,1,1,0)
      call gs_op_fields(gs_handle_q_H,tmpiv,nstate_p,1,1,1,0)
c      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
      qctime1= qctime1+dclock()-ttmp     

!$acc parallel loop collapse(2) private(tmpfq,tmprr,tmpii,tmpcc)         
c!$omp parallel do collapse(1)
      do jj = 1,nstate/np
      do ii = 1,nstate
         tmpfq= state_sums(ii,jj)*pulse_omega_0*rho(ii,jj)
         tmprr= 2*tmp_r1(ii,jj)-tmprv((jj-1)*nstate+ii)
         if(ii.eq.jj+c_offset) tmprr = 0
         tmpii = tmpiv((jj-1)*nstate+ii)
         if(ii.eq.jj+c_offset) tmpii = 2*tmp_i1(ii,jj)
         tmpcc= tmprr + ci*tmpii
         drho(ii,jj)= - ci*(tmpcc - tmpfq) 
         resrho_r(ii,jj)=real(drho(ii,jj) )    
         resrho_i(ii,jj)=imag(drho(ii,jj) )    
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

!$acc end data

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
      integer  ii,jj,iq,ns,nsp,kns,knsp,kk,ll,mm,rrequest(nstate/np)
      integer  nqtmp1(lnumqd),nqtmp2(lnumqd),req,srequest(nstate/np)
      real*8   nq,nqp
      integer  nstmp1,nsptmp1,target_core,ierr,ltoiq
      real*8   fac,tmp,tmp3,tmp1,tmp2,tmpr1,tmpi1
      real*8   tmpr2(nstate),tmpi2(nstate)
      complex*16  rho_passing,ci,tmpc2(nstate),tmpsnd(nstate)
      logical  lbFlag,qdFlag,sFlag
      integer status(MPI_STATUS_SIZE),get_counter
      integer*8 target_disp !NOT PORTABLE: MPI_ADDRESS_KIND
      real*8   dclock, ttmp
      ci = (0,1.0)
      rho_passing=0
      get_counter=0
      qdFlag=.false.
      sFlag=.false.

!$acc data present(rho_r,rho_i,rho,drho,resrho_i,resrho_r,tmp_c1)
!$acc&     present(gamma_d_m,gamma_q_m,nq_minimal)
!$acc&     pcreate(tmpr2,tmpi2,tmpc2,tmpsnd)


c      do jj=1,nstate/np
c         do iq=1,numqd
c            ltoiq = level**iq
c            nqp = nq_minimal(iq,mod(jj+c_offset-1-ltoiq,ltonq)+1)
c            target_core = (jj+c_offset-ltoiq-1)*np/nstate
c            if(target_core.ge.0.and.target_core.ne.nid.and.nqp.eq.0)
c     $           then
c               !pack snd buffer
c               do ii=1,nstate
c                  tmpsnd(ii) = rho(ii,jj)
c               enddo
cc               print*,'send ',nid,jj+c_offset,target_core,
cc     $              jj+c_offset-ltoiq
c               call MPI_ISEND(tmpsnd,nstate,MPI_DOUBLE_COMPLEX,
c     $              target_core,0,MPI_COMM_WORLD,
c     $              srequest(jj),ierr)
c            endif
c         enddo
c      enddo



!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)
      do jj=1,nstate/np
      do ii=1,nstate
c         rho(ii,jj) = 1 + ci*1
         rho(ii,jj) = rho_r(ii,jj) + ci*rho_i(ii,jj)
         tmp_c1(ii,jj) = rho(ii,jj)
c         tmp_r1(ii,jj) = real(rho(ii,jj))
c         tmp_i1(ii,jj) = imag(rho(ii,jj))
 
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop
c
c     Is this barrier needed?
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!$acc parallel loop collapse(2) private(nsp,ns,tmp3,nq,nqp,tmp1)
c!$omp parallel do collapse(1)
      do jj = 1,nstate/np
      do ii = 1,nstate
            nsp  = (jj+c_offset-1)/ltonq
            ns   = (ii-1)/ltonq
         do iq=1,numqd
            nq   = nq_minimal(iq,mod(ii-1,ltonq)+1)
            nqp  = nq_minimal(iq,mod(jj+c_offset-1,ltonq)+1)
            tmp1 = gamma_d_m(iq)*(nq-nqp)*(nq-nqp) + 
     $           gamma_q_m(iq)*(nq+nqp)*0.5
c     Local Terms:
c...  QD dephasing Linblad term leads to:
c...  QD spontaneous emission Linblad term leads to these two terms:
c...  SP spontaneous emission Linblad terms:
            drho(ii,jj)= drho(ii,jj)-tmp1*rho(ii,jj)
         enddo
         tmp3 = gamma_s*(ns+nsp)*0.5
         drho(ii,jj) = drho(ii,jj) - tmp3*rho(ii,jj)
      enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop

c     Nonlocal terms, gamma_q

ccc Synchronization here? Everyone needs to be done before we can start
ccc MPI_GETing from the other cores

!$acc update host(rho,drho,tmp_c1)
!      if(nid.eq.0) print*,'before host data
      do jj=1,nstate/np
         ltoiq = 1
         target_core = (jj+c_offset+ltoiq-1)*np/nstate  
         do iq=1,numqd
            nqp = nq_minimal(iq,mod(jj+c_offset-1,ltonq)+1)
            if(target_core.ne.nid.and.target_core.lt.np.and.nqp.eq.0)
     $           then
               target_disp = (jj+c_offset+ltoiq-target_core
     $              *nstate/np-1)*nstate
               ttmp= dclock() 
c               call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
c     $              mpi_win_r,ierr)

c               call MPI_GET(tmpr2,nstate,MPI_REAL8
c     $              ,target_core,target_disp,nstate,MPI_REAL8,mpi_win_r
c     $              ,IERR)

c               call MPI_WIN_UNLOCK(target_core,mpi_win_r,ierr)

c               call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
c     $              mpi_win_i,ierr)
               
c               call MPI_GET(tmpi2,nstate,MPI_REAL8
c     $              ,target_core,target_disp,nstate,MPI_REAL8,mpi_win_i
c     $              ,IERR)
c               call MPI_WIN_UNLOCK(target_core,mpi_win_i,ierr)
c!$acc host_data use_device(tmpc2)

               call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,
     $              0,mpi_win_c,ierr)
               call MPI_GET(tmpc2,nstate,MPI_DOUBLE_COMPLEX
     $              ,target_core,target_disp,nstate,MPI_DOUBLE_COMPLEX
     $              ,mpi_win_c,IERR)
               call MPI_WIN_UNLOCK(target_core,mpi_win_c,ierr)

c               call MPI_IRECV(tmpc2,nstate,MPI_DOUBLE_COMPLEX,
c     $              target_core,0,MPI_COMM_WORLD,rrequest(jj),ierr)
c               print*,'recv ',target_core,jj+c_offset+ltoiq,
c     $              nid,jj+c_offset

c!$acc end host_data
               qctime3= qctime3+dclock()-ttmp
            endif
            
c!$acc kernels 
            do ii=1,nstate
               nq = nq_minimal(iq,mod(ii-1,ltonq)+1)
               if(nqp.eq.0.and.nq.eq.0) then
                  if(target_core.eq.nid) then
c     The nonlocal term belongs to me, do the local computation
                     drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
     $                    rho(ii+ltoiq,jj+ltoiq)
c                     print*,'local',ii+level**(iq-1),jj+level**(iq-1)
c     $                    ,nid
                  else
c     The nonlocal term belongs to a different core; it is the buffer
c                     drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
c     $                    (tmpr2(ii+level**(iq-1))+
c     $                    ci*tmpi2(ii+level**(iq-1)))
c                     call MPI_WAIT(rrequest(jj),status,ierr)
c!$acc update device(tmpc2)
                     drho(ii,jj) = drho(ii,jj)+gamma_q_m(iq)*
     $                    (tmpc2(ii+ltoiq))
c                     print*,'nonlocal',ii+level**(iq-1),jj+level**(iq-1)
c                     print*,'rho_r',tmpr2(ii+level**(iq-1))
c                     print*,'rho_i',tmpi2(ii+level**(iq-1))
c                     do kk=1,nstate
c                        print*,'r',kk,tmpr2(kk)
c                        print*,'i',kk,tmpi2(kk)
c                     enddo
                  endif
               endif
            enddo
c!$acc end kernels
            ltoiq = ltoiq*level
         enddo
      enddo



c     Nonlocal terms, gamma_s
      do jj=1,nstate/np
         nsp= (jj+c_offset-1)/ltonq
         
         target_core = (jj+c_offset+ltonq-1)*np/nstate
         
         if(target_core.ne.nid.and.target_core.lt.np) then
            target_disp = (jj+c_offset+ltonq-target_core
     $           *nstate/np-1)*nstate
            ttmp = dclock()
c            call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
c     $           mpi_win_r,ierr)
            
c            call MPI_GET(tmpr2,nstate,MPI_REAL8
c     $           ,target_core,target_disp,nstate,MPI_REAL8,mpi_win_r
c     $           ,IERR)
            
c            call MPI_WIN_UNLOCK(target_core,mpi_win_r,ierr)
            
c            call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
c     $           mpi_win_i,ierr)
            
c            call MPI_GET(tmpi2,nstate,MPI_REAL8
c     $           ,target_core,target_disp,nstate,MPI_REAL8,mpi_win_i
c     $           ,IERR)
c            call MPI_WIN_UNLOCK(target_core,mpi_win_i,ierr)
c!$acc host_data use_device(tmpc2,tmp_c1)

            call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
     $           mpi_win_c,ierr)

            call MPI_GET(tmpc2,nstate,MPI_DOUBLE_COMPLEX
     $           ,target_core,target_disp,nstate,MPI_DOUBLE_COMPLEX
     $           ,mpi_win_c,IERR)

            call MPI_WIN_UNLOCK(target_core,mpi_win_c,ierr)
c!$acc end host_data

            qctime3= qctime3+dclock()-ttmp
         endif
         
c!$acc kernels 
         do ii=1,nstate
            if(ii+ltonq.lt.nstate.and.jj+c_offset+ltonq.lt.nstate) then
               ns = (ii-1)/ltonq
               
               if(target_core.eq.nid) then
c     Nonlocal terms belong to me
                  drho(ii,jj) =drho(ii,jj)+gamma_s*sqrt(real((ns+1)*
     $                 (nsp+1)))*rho(ii+ltonq,jj+ltonq)
               else
c     Nonlocal terms belong to someone else
c                  drho(ii,jj) = drho(ii,jj)+gamma_s*sqrt(real((ns+1)*
c     $                 (nsp+1)))*(tmpr2(ii+ltonq)+ci*tmpi2(ii+ltonq))
                  drho(ii,jj) = drho(ii,jj)+gamma_s*sqrt(real((ns+1)*
     $                 (nsp+1)))*tmpc2(ii+ltonq)
               endif
               
            endif
         enddo
c!$acc end kernels
      enddo

!$acc update device(drho)

!$acc parallel loop collapse(2)
c!$omp parallel do collapse(1)         
      do jj=1,nstate/np
         do ii=1,nstate
            resrho_r(ii,jj)= real(drho(ii,jj))
            resrho_i(ii,jj)= imag(drho(ii,jj))
         enddo
      enddo
c!$omp end parallel do
!$acc end parallel loop


!$acc end data


c      do iq=1,numqd
c         ltoiq = level**(iq-1)
c         target_core = nid+level**(iq-1)
c         target_disp = 0
c         if(level**(iq-1).ge.nstate/np) then
c           The terms are on another core. Go get them.
c            ttmp= dclock() 
c            call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
c     $           mpi_win_r,ierr)
c            call MPI_GET(tmpr2,nstate*nstate/np,MPI_REAL8
c     $           ,target_core,target_disp,1,MPI_REAL8,mpi_win_r
c     $           ,IERR)
c            call MPI_WIN_UNLOCK(target_core,mpi_win_r,ierr)
            
c            call MPI_WIN_LOCK(MPI_LOCK_SHARED,target_core,0,
c     $           mpi_win_i,ierr)
            
c            call MPI_GET(tmpi2,nstate*nstate/np,MPI_REAL8
c     $           ,target_core,target_disp,1,MPI_REAL8,mpi_win_i
c     $           ,IERR)
c            call MPI_WIN_UNLOCK(target_core,mpi_win_i,ierr)
c            qctime3= qctime3+dclock()-ttmp
            
c            do jj=1,nstate/np




c            enddo
c         else
c        Local terms


c      if(time.gt.1000) then
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
c      if(time.gt.20) call exitt
c      endif
cc      endif
c      if(np.eq.1) then
c         call outmat(resrho_r(1,1),nstate,nstate/2,'rho_r  ',nid) 
c         call outmat(resrho_r(1,7),nstate,nstate/2,'rho_r  ',nid)
c         call outmat(resrho_i(1,1),nstate,nstate/2,'rho_i  ',nid)
c         call outmat(resrho_i(1,7),nstate,nstate/2,'rho_i  ',nid)
c         call exitt
c      endif
c     Is this barrier necesarry? 
c      call MPI_BARRIER(MPI_COMM_WORLD,IERR)
c      call exitt
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

!$acc data present(state_sums,rho_r,rho_i,rho)
!$acc parallel loop collapse(2)
      do jj=1,nstate/np
      do ii=1,nstate
         rho(ii,jj)= (rho_r(ii,jj)+ci*rho_i(ii,jj))
     $        *exp(-ci*state_sums(ii,jj)*pulse_omega_0*(time-pulse_t_0))
      enddo
      enddo
!$acc end parallel loop
!$acc update host(rho)
!$acc end data
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

      sp_mu = 0

c     FIXME: 12/11/14 Matt- calculates sp_mu wrong (but rho and pop are right)

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
     $        popS,popQD(1),popQD(2)
c     $        concurrence

         write(14,12)time*timeunit,popQD,popS
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

      call quantum2_concurrence()

      end

c**************************
c
c subroutine quantum2_concurrence finds and prints the concurrence
c using wooters formula. 
c
c
c**************************

c-----------------------------------------------------------------------
      subroutine quantum2_concurrence()
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUMMIN'
      include 'mpif.h'
      complex red_rho2(4,4),red_rho3(8,8),ci
      real*8  conc1,conc2,conc3,conc4,conc5,conc6
      real*8  popSymm1(4),popSymm2(4),popSymm3(4),popSymm4(4)
      real*8  popSymm5(4),popSymm6(4)
      real*8  tmprv(lstate*lstate/lp)
      real*8  tmpiv(lstate*lstate/lp)
      integer ii,jj,j0

      ci=(0.0,1.0)

c     Partial trace to maxPartialTrace size

c     Prepare communcation vector
c!$acc parallel loop collapse(2) private(j0)
c!$omp parallel do collapse(2)
      do jj = 1,nstate/np
      do ii = 1,nstate
         j0= (jj-1)*nstate+ii
         tmprv(j0) = rho_r(ii,jj)
         tmpiv(j0) = rho_i(ii,jj)
      enddo
      enddo
c!$omp end parallel do
c!$acc end parallel loop

      call gs_op_fields(gs_handle_q_conc,tmprv,
     $     nstate_p,1,1,1,0) 
      call gs_op_fields(gs_handle_q_conc,tmpiv,
     $     nstate_p,1,1,1,0) 

c     store result in red_rho
      do jj=1,maxPartialTrace
         do ii=1,maxPartialTrace
            j0= (jj-1)*nstate+ii
            red_rho(ii,jj) = (tmprv(j0) + ci*tmpiv(j0))
     $        *exp(-ci*state_sums(ii,jj)*pulse_omega_0*(time-pulse_t_0))
         enddo
      enddo

c$$$      call outmat(rho_r,nstate,nstate,'rho_r',nid)
c$$$      call outmat(rho_i,nstate,nstate,'rho_i',nid)
c$$$      if(nid.eq.0) call outmat(red_rho,8,8,'red_rho',nid)
c$$$      call exitt()
      if((maxPartialTrace.eq.4).and.(nid.eq.0)) then
         call cem_quantum_bipartite_concurrence(red_rho,conc1,popSymm1)
         write(13,13)time*timeunit,conc1
 13      format(1p20g20.8)
         write(15,13)time*timeunit,popSymm1         
      endif



      if((maxPartialTrace.eq.8).and.(nid.eq.0)) then
         call cem_quantum_partial_trace_3_to_2(red_rho,red_rho2,1)
c     concurrence(3:2)
         call cem_quantum_bipartite_concurrence(red_rho2,conc1
     $        ,popSymm1)
         call cem_quantum_partial_trace_3_to_2(red_rho,red_rho2,2)
c     concurrence(3:1)
         call cem_quantum_bipartite_concurrence(red_rho2,conc2
     $        ,popSymm2)
         call cem_quantum_partial_trace_3_to_2(red_rho,red_rho2,3)
c     concurrence(2:1)
         call cem_quantum_bipartite_concurrence(red_rho2,conc3
     $        ,popSymm3)
         
         write(13,13)time*timeunit,conc1,conc2,conc3
         write(15,13)time*timeunit,popSymm1,popSymm2
     $        ,popSymm3

      endif

      if((maxPartialTrace.eq.16).and.(nid.eq.0)) then
c     First, trace out qd 1 with 4->3 quantum dots
         call cem_quantum_partial_trace_4_to_3(red_rho,red_rho3,1)
c     Now, trace 3->2 3 times, and find bipartite concurrence for each
c     First, concurrence(4:3) (trace out 2)
         call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,1)
         call cem_quantum_bipartite_concurrence(red_rho2,conc1
     $        ,popSymm1)
c     concurrence(4:2) (trace out 3)
         call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,2)
         call cem_quantum_bipartite_concurrence(red_rho2,conc2
     $        ,popSymm2)
c     concurrence(3:2) (trace out 4)
         call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,3)
         call cem_quantum_bipartite_concurrence(red_rho2,conc3
     $        ,popSymm3)
c     trace out qd 3 from full
         call cem_quantum_partial_trace_4_to_3(red_rho,red_rho3,3)
c     concurrence(4:1)
         call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,2)
         call cem_quantum_bipartite_concurrence(red_rho2,conc4
     $        ,popSymm4)
c     concurrence(2:1)
         call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,3)
         call cem_quantum_bipartite_concurrence(red_rho2,conc5
     $        ,popSymm5)
c     trace out qd 4 from full
         call cem_quantum_partial_trace_4_to_3(red_rho,red_rho3,4)
c     concurrence(3:1)
         call cem_quantum_partial_trace_3_to_2(red_rho3,red_rho2,3)
         call cem_quantum_bipartite_concurrence(red_rho2,conc6
     $        ,popSymm6)
         write(13,13)time*timeunit,conc1,conc2,conc3,
     $        conc4,conc5,conc6
         write(15,13)time*timeunit,popSymm1,popSymm2
     $        ,popSymm3,popSymm4,popSymm5,popSymm6
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
c$$$      if((time.gt.10000).and.(time.lt.10500)) then
c$$$      call outmat(resrho_r,nstate,nstate,'res_r  ',nid)
c$$$      call outmat(resrho_i,nstate,nstate,'res_i  ',nid)
c$$$      call outmat(rho_r,nstate,nstate,'rho_r  ',nid)
c$$$      call outmat(rho_i,nstate,nstate,'rho_r  ',nid)
c$$$      endif
c$$$      if(time.gt.10500) call exitt
#ifdef GPU
#ifdef _OPENACC
      call rk4_upd_acc(rho_r,k_rho_r,resrho_r,cb,ca,dt,nstate2)
      call rk4_upd_acc(rho_i,k_rho_i,resrho_i,cb,ca,dt,nstate2)
#endif
#else
      call rk4_upd(rho_r,k_rho_r,resrho_r,cb,ca,dt,nstate2)
      call rk4_upd(rho_i,k_rho_i,resrho_i,cb,ca,dt,nstate2)
#endif

      return
      end
