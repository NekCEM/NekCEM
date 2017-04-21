c--------------------------------------------------------------------- 
c...  This file includes the CUDA routines for solving Maxwell eqs.
c...  Globally defined variable: npts, npts3, nxyz, nxzf,nxzfl, nfaces
c--------------------------------------------------------------------- 

c...  preproessor compile options: CPU/OMP/CUDA/OpenACC
#define GPU_VALIDATE 0  ! testing purpose for cuda 
#ifdef GPU

c-----------------------------------------------------------------------
      subroutine ck_gpu_val(expected,actual,k,msg)
c-----------------------------------------------------------------------
      implicit none
      real tol,expected,actual
      integer k
      character msg*(*)
      
      tol = 1.0E-8
      if (abs(expected-actual).ge.tol) then
        write(6,'(A,A, I10, E30.16, E30.16)')
     $    'validation error: ',msg,k,expected,actual
      endif
      end

c-----------------------------------------------------------------------
      subroutine maxwell_wght_curl_cuda(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      include 'SIZE'
      include 'TOTAL'
c
c     real    w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      real    w1(lpts),w2(lpts),w3(lpts),u1(lpts),u2(lpts),u3(lpts)
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      real    rxmn_k, sxmn_k, txmn_k
      real    rymn_k, symn_k, tymn_k
      real    rzmn_k, szmn_k, tzmn_k
      integer e,k,i,j,ntot,nn,szrl,sznpts,szlpts1
      integer icalld
      save    icalld
      data    icalld/0/
      real*8  dclock    
      real    varw
      real    wght_curl_time               
      real    gpu_time1,gpu_time2
      real    cpu_time1,cpu_time2
      real    cpu_time1_total,cpu_time2_total
c.............
      common /ctmp1/ u1r(lpts1),u1s(lpts1),u1t(lpts1) 
     $             , u2r(lpts1),u2s(lpts1),u2t(lpts1)
     $             , u3r(lpts1),u3s(lpts1),u3t(lpts1)
      real u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t    
c      real u1r(lpts1),u1s(lpts1),u1t(lpts1) 
c     $   , u2r(lpts1),u2s(lpts1),u2t(lpts1)
c     $   , u3r(lpts1),u3s(lpts1),u3t(lpts1)
      structure /memptr/
         integer sync
         integer sz
         integer host
         integer dev
         character vname*10
      end structure
      record /memptr/ mp_d, mp_dt ! D and Dt memptr's 
      record /memptr/ mp_w3mn
      record /memptr/ mp_rxmn,mp_sxmn,mp_txmn,
     $                mp_rymn,mp_symn,mp_tymn,
     $                mp_rzmn,mp_szmn,mp_tzmn,
     $                mp_u1r,mp_u1s,mp_u1t,
     $                mp_u2r,mp_u2s,mp_u2t,
     $                mp_u3r,mp_u3s,mp_u3t,
     $                mp_u1,mp_u2,mp_u3,
     $                mp_w1,mp_w2,mp_w3
#ifdef GPU_VALIDATE
      real cu1r(lxyzm),cu1s(lxyzm),cu1t(lxyzm)
     $    ,cu2r(lxyzm),cu2s(lxyzm),cu2t(lxyzm)
     $    ,cu3r(lxyzm),cu3s(lxyzm),cu3t(lxyzm)
      real cw1(1),cw2(1),cw3(1)
#endif
      
c.............
      nn    = nx1-1

      if (icalld.eq.0) then
         write(6,*) 'tt== lpts1',lpts1,npts,nn,nelt,nxyz
         write(6,*) 'tt== lx1',lx1,ly1,lz1,lelt

         szrl   = 8
         sznpts = szrl*npts
         szlpts1= szrl*lpts1
         ! initialize D memptr
         mp_d.sync    = 3
         mp_d.sz      = szrl*nx1*nx1
         mp_d.host    = LOC(dxm1)
         ! initialize Dt memptr
         mp_dt.sync   = 3
         mp_dt.sz     = szrl*nx1*nx1
         mp_dt.host   = LOC(dxtm1)
         mp_w3mn.sync = 3
         mp_w3mn.sz   = szrl*nxyz
         mp_w3mn.host = LOC(w3mn)
         mp_rxmn.sync = 3
         mp_rxmn.sz   = sznpts
         mp_rxmn.host = LOC(rxmn)
         mp_rymn.sync = 3
         mp_rymn.sz   = sznpts
         mp_rymn.host = LOC(rymn)
         mp_rzmn.sync = 3
         mp_rzmn.sz   = sznpts
         mp_rzmn.host = LOC(rzmn)
         mp_sxmn.sync = 3
         mp_sxmn.sz   = sznpts
         mp_sxmn.host = LOC(sxmn)
         mp_symn.sync = 3
         mp_symn.sz   = sznpts
         mp_symn.host = LOC(symn)
         mp_szmn.sync = 3
         mp_szmn.sz   = sznpts
         mp_szmn.host = LOC(szmn)
         mp_txmn.sync = 3
         mp_txmn.sz   = sznpts
         mp_txmn.host = LOC(txmn)
         mp_tymn.sync = 3
         mp_tymn.sz   = sznpts
         mp_tymn.host = LOC(tymn)
         mp_tzmn.sync = 3
         mp_tzmn.sz   = sznpts
         mp_tzmn.host = LOC(tzmn)

         mp_u1r.sync  = 1
         mp_u1r.sz    = sznpts
         mp_u1r.host  = LOC(u1r)
         mp_u2r.sync  = 1
         mp_u2r.sz    = sznpts
         mp_u2r.host  = LOC(u2r)
         mp_u3r.sync  = 1
         mp_u3r.sz    = sznpts
         mp_u3r.host  = LOC(u3r)
         mp_u1s.sync  = 1
         mp_u1s.sz    = sznpts
         mp_u1s.host  = LOC(u1s)
         mp_u2s.sync  = 1
         mp_u2s.sz    = sznpts
         mp_u2s.host  = LOC(u2s)
         mp_u3s.sync  = 1
         mp_u3s.sz    = sznpts
         mp_u3s.host  = LOC(u3s)
         mp_u1t.sync  = 1
         mp_u1t.sz    = sznpts
         mp_u1t.host  = LOC(u1t)
         mp_u2t.sync  = 1
         mp_u2t.sz    = sznpts
         mp_u2t.host  = LOC(u2t)
         mp_u3t.sync  = 1
         mp_u3t.sz    = sznpts
         mp_u3t.host  = LOC(u3t)

         mp_u1.sync   = 1
         mp_u2.sync   = 1
         mp_u3.sync   = 1
         mp_u1.sz     = szlpts1
         mp_u2.sz     = szlpts1
         mp_u3.sz     = szlpts1

         mp_w1.sync   = 1
         mp_w2.sync   = 1
         mp_w3.sync   = 1
         mp_w1.sz     = szlpts1
         mp_w2.sz     = szlpts1
         mp_w3.sz     = szlpts1

         icalld=1
      endif

      ! CPU->GPU
      mp_u1.sync=mp_u1.sync.OR.2
      mp_u2.sync=mp_u2.sync.OR.2
      mp_u3.sync=mp_u3.sync.OR.2
      mp_u1.host=LOC(u1)
      mp_u2.host=LOC(u2)
      mp_u3.host=LOC(u3)
      ! GPU->CPU
      mp_w1.sync=mp_w1.sync.OR.4
      mp_w2.sync=mp_w2.sync.OR.4
      mp_w3.sync=mp_w3.sync.OR.4
      mp_w1.host=LOC(w1)
      mp_w2.host=LOC(w2)
      mp_w3.host=LOC(w3)

      wght_curl_time=dclock()
      gpu_time1     =dclock()

      call local_grad3_gpu (mp_u1r,mp_u1s,mp_u1t,
     $                      mp_u2r,mp_u2s,mp_u2t,
     $                      mp_u3r,mp_u3s,mp_u3t,
     $                      mp_u1 ,mp_u2 ,mp_u3 ,
     $                      mp_d,  mp_dt, nx1,nelt,lpts1,nid)

      gpu_time1=dclock()-gpu_time1

      gpu_time2=dclock()

      call curl_gpu(mp_u1r, mp_u1s, mp_u1t,
     $              mp_u2r, mp_u2s, mp_u2t,
     $              mp_u3r, mp_u3s, mp_u3t,
     $              mp_rxmn,mp_sxmn,mp_txmn,
     $              mp_rymn,mp_symn,mp_tymn,
     $              mp_rzmn,mp_szmn,mp_tzmn,
     $              mp_w1,  mp_w2,  mp_w3,mp_w3mn,nxyz,nelt,lpts1)

      gpu_time2=dclock()-gpu_time2
      wght_curl_time=dclock()-wght_curl_time

#ifdef GPU_VALIDATE
      do e=1,nelt
         j = nxyz*(e-1)
         call local_grad3(cu1r,cu1s,cu1t,u1,nn,e,dxm1,dxtm1)        
         call local_grad3(cu2r,cu2s,cu2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(cu3r,cu3s,cu3t,u3,nn,e,dxm1,dxtm1)
         do i=1,nxyz
            k = i+j
c            call ck_gpu_val(cu1r(i),u1r(k),k,'u1r')
c            call ck_gpu_val(cu2r(i),u2r(k),k,'u2r')
c            call ck_gpu_val(cu3r(i),u3r(k),k,'u3r')
c            call ck_gpu_val(cu1s(i),u1s(k),k,'u1s')
c            call ck_gpu_val(cu2s(i),u2s(k),k,'u2s')
c            call ck_gpu_val(cu3s(i),u3s(k),k,'u3s')
c            call ck_gpu_val(cu1t(i),u1t(k),k,'u1t')
c            call ck_gpu_val(cu2t(i),u2t(k),k,'u2t')
c            call ck_gpu_val(cu3t(i),u3t(k),k,'u3t')

            varw = w3mn(i)
            cw1(k)= cu3r(i)*varw*rymn(k)
     $            + cu3s(i)*varw*symn(k)
     $            + cu3t(i)*varw*tymn(k)
     $            - cu2r(i)*varw*rzmn(k)
     $            - cu2s(i)*varw*szmn(k)
     $            - cu2t(i)*varw*tzmn(k)
            call ck_gpu_val(cw1(k),w1(k),k,'w1')

            cw2(k)= cu1r(i)*varw*rzmn(k)
     $            + cu1s(i)*varw*szmn(k)
     $            + cu1t(i)*varw*tzmn(k)
     $            - cu3r(i)*varw*rxmn(k)
     $            - cu3s(i)*varw*sxmn(k)
     $            - cu3t(i)*varw*txmn(k)
            call ck_gpu_val(cw2(k),w2(k),k,'w2')

            cw3(k)= cu2r(i)*varw*rxmn(k)
     $            + cu2s(i)*varw*sxmn(k)
     $            + cu2t(i)*varw*txmn(k)
     $            - cu1r(i)*varw*rymn(k)
     $            - cu1s(i)*varw*symn(k)
     $            - cu1t(i)*varw*tymn(k)
            call ck_gpu_val(cw3(k),w3(k),k,'w3')
         enddo
      enddo
#endif

#ifdef TIMER 
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
      write(6,*) '----------------------------------::'
      write(6,*) 'maxwell_wght_curl_cuda       (sec)::',wght_curl_time
      write(6,*) '        local_grad3_gpu::gpu_time1::',gpu_time1     
      write(6,*) '                         gpu_time2::',gpu_time2
      write(6,*) '----------------------------------::'
      endif
      endif
#endif

       return
       end
#endif
