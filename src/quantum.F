c--------------------------------------------------------------------- 
      subroutine cem_quantum_hamiltonian
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'QUANTUM'

      integer  is,iq,ii,jj   

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

c...  set the diagonal terms
      do  ii = 1,nstate 
          ns = isstate(ii) 
          nq = iqstate(ii) 
          hamilt_0 (ii,ii) = ns*oms+nq*omq 
      enddo

c...  set q = 0, s couples to q' = 1, s' = s-1
      do  ii = 3,nstate,2 ! 3,5,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 

          jj = ii-1                            
          nsp= isstate(jj) 
          nqp= iqstate(jj) 
          hamilt_0(ii,jj) = -g_coupling*mus*muq*sqrt(nsp+1) 
      enddo

c...  q = 1, s couples to q' = 0, s' = s+1
      do  ii = 2,nstate,2 ! 2,4,6,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 

          jj = ii+1
          nsp= isstate(jj) 
          nqp= iqstate(jj) 
          hamilt_0(ii,jj) = -g_coupling*mus*muq*sqrt(nsp) 
      enddo

c...  if (ns==nsp) + { (nq==0)&(nqp==1)|(nq==1)&(nqp==0) }                              
      do  ii = 3,nstate,2 ! 3,5,... 
          ns = isstate(ii) 
          nq = iqstate(ii) 
            
          nsp= isstate(ii+1) 
          nqp= iqstate(ii+1) 
          hamilt(ii  ,ii+1)=-muq       
          hamilt(ii+1,ii-1)=-muq       
      enddo

c...  if (nq==nqp) + { ns==(nsp-1) }                      
      do  ii = 1,nstate   ! 1,2,3,... 
          ns = isstate(ii)
          nq = iqstate(ii)

          nsp= isstate(ii+2)
          nqp= iqstate(ii+2)
          hamilt(ii,ii+2)=-sqrt(nsp)*mus
      enddo

c...  if (nq==nqp) + { ns==(nsp+1) }                      
      do  ii = 3,nstate  ! 1,2,3,... 
          ns = isstate(ii)
          nq = iqstate(ii)

          nsp= isstate(ii-2)
          nqp= iqstate(ii-2)
          hamilt(ii,ii-2)=-sqrt(nsp+1)*mus
      enddo

c... in rwa, each driving term energy0*cos(omega t) becomes an energy0/2 term:
      tmp= energy0*0.5 
      do ii=1,nstate2
         hamilt  (ii,1)=hamilt(ii,1)*tmp           
         hamilt_I(ii,1)=hamilt_0(ii,1)+hamilt(ii,1)
      enddo




