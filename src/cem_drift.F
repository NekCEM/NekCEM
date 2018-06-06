c---------------------------------------------------------------------
      subroutine cem_drift_init
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'ZPER'
      include 'BCS'
      include 'NEWTON' ! flag

      integer  i, npts3, nxzfl3
      integer  npts21,lpts13

      real dummy(lx1*ly1*lz1*lelt)

      npts3 = 3*npts
      lpts13 = 3*lpts1
      nxzfl3= 3*nxzfl
      npts21 = 21*npts

c...  these to be deleted all when IFHEX/IFTET are done

      if (nid.eq.0) write(6,*) 'set variables: cem_drift_temporary'
      call rzero (xmn ,npts)
      call rzero (ymn ,npts)
      call rzero (zmn ,npts)
      call rzero (bmn ,npts)
      call rzero (rxmn,npts)
      call rzero (rymn,npts)
      call rzero (rzmn,npts)
      call rzero (sxmn,npts)
      call rzero (symn,npts)
      call rzero (szmn,npts)
      call rzero (txmn,npts)
      call rzero (tymn,npts)
      call rzero (tzmn,npts)
      call rzero (jacm,npts)

      call rzero (cN   ,lpts1*lcdim)
      call rzero (cN1  ,npts)
      call rzero (cP1  ,npts)
      call rzero (cE1  ,npts)
      call rzero (scN  ,lpts1*lcdim)
      call rzero (gradp,lpts13)
      call rzero (cNJ  ,lpts13*lcdim)
      call rzero (rescN ,lpts1*lcdim)
      call rzero (rescNJ,lpts13*lcdim)

      call rzero (charge,npts)
      call rzero (potent,npts)
      call rzero (potent1,npts)
      call rzero (rhs   ,npts)
      call rzero (spotent,npts)
      call rzero (epotent,npts)

      call rzero (errcN,lpts1*lcdim)
      call rzero (kcN  ,lpts1*lcdim)
      call rzero (qbm1 ,npts)

      call rzero (fcN  ,nxzfl)
      call rzero (fcNJ ,nxzfl3)

      call rzero (bc_flux,lxzfl*ldimt)
      call rzero (bc_conv,lxzfl*ldimt)

      call rzero (w3mn ,nxyz )
      call rzero (unxm ,nxzfl)
      call rzero (unym ,nxzfl)
      call rzero (unzm ,nxzfl)
      call rzero (aream,nxzfl)
      call rzero (approx_o,npts21)

      do i=1,npts
         xmn (i)= xm1 (i,1,1,1)
         ymn (i)= ym1 (i,1,1,1)
         zmn (i)= zm1 (i,1,1,1)
         bmn (i)= bm1 (i,1,1,1)
         rxmn(i)= rxm1(i,1,1,1)
         rymn(i)= rym1(i,1,1,1)
         rzmn(i)= rzm1(i,1,1,1)
         sxmn(i)= sxm1(i,1,1,1)
         symn(i)= sym1(i,1,1,1)
         szmn(i)= szm1(i,1,1,1)
         txmn(i)= txm1(i,1,1,1)
         tymn(i)= tym1(i,1,1,1)
         tzmn(i)= tzm1(i,1,1,1)
         jacm(i)= jacm1(i,1,1,1)
      enddo

      do i=1,nxyz
         w3mn(i)= w3m1 (i,1,1)
      enddo

      do i=1,nxzfl
         unxm (i)= unx (i,1,1,1)
         unym (i)= uny (i,1,1,1)
         unzm (i)= unz (i,1,1,1)
         aream(i)= area(i,1,1,1)
      enddo
      if (nid.eq.0) write(6,*) 'done: temporary variables for geom'

c...  get global numbering on face and set material properties
      call cem_set_fc_ptr       !global numbering index on face
      call cem_drift_uvp        !set parameters
      call cem_drift_bcs_fc_ptr !bounday settings

c...  inverse mass matrix including material constants
      call invers2(qbm1,bmn,npts)

      if (nid.eq.0) write(6,*) 'done: inverse mass matrix'

      if (ifrestart) then
         call restart_swap
      else
         call userini(time,cN)
      endif

c.... setup masks for multi-fields for drift
      do i=1,nfield
        call rone (dmask(1,i),npts)
        call setdmask (dmask(1,i),npts,i)
      enddo

c.... setup multiplicity 
      call rone    (mult,npts)
      call dssum   (mult,lx1,ly1,lz1)
      call invcol1 (mult,npts)

      ifparamNT = .false.
      if (abs(param(1)).ne.0) then   ! todo restr Newton later
        call cem_drift_newton
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_uvp
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'

      call uservp(1,1,1,1)

      return
      end
c---------------------------------------------------------------------
c.... 1st/2nd order semi-implicit BDF scheme
      subroutine cem_drift_op_bdf
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'      

      integer ic
      real tmp_cN(lpts1,lcdim),tmp_po(lpts1)  ! Rishardson
      real r2_cN(lpts1,lcdim),r2_po(lpts1)!

      if     (ifbdf1) then    
        
        call cem_drift_source  
        call cem_drift_poisson
        if (ifsteric) then
          call cem_drift_sem_bdf1_steric
        else 
          call cem_drift_sem_bdf1
        endif

      elseif (ifbdf1r) then ! Richardson on BDF1

        do ic = 1,lcdim ! copy c_{n-1}
          call copy(tmp_cN(1,ic),cN(1,ic),npts)
        enddo

C       half time step and do twice: c_{n-1} -> c_{n-1/2} -> c_{n}
        dt = 0.5*dt

        call cem_drift_source
        call cem_drift_poisson 
        call copy(tmp_po,potent,npts) ! re-use
        if (ifsteric) then
          call cem_drift_sem_bdf1_steric
        else 
          call cem_drift_sem_bdf1
        endif

        time = time+dt! time = t_{n-1/2}

        call cem_drift_source
        call cem_drift_poisson 
        if (ifsteric) then
          call cem_drift_sem_bdf1_steric
        else 
          call cem_drift_sem_bdf1
        endif

        time = time-dt ! time = t_{n-1}

        do ic = 1,lcdim
          call copy(r2_cN(1,ic),cN(1,ic),npts)
        enddo
c        call copy(r2_po,potent,npts)
        do ic = 1,lcdim
          call copy(cN(1,ic),tmp_cN(1,ic),npts)
        enddo
        call copy(potent,tmp_po,npts)

C       origin BDF1: c_{n-1} -> c_{n}
        dt = 2.0*dt

        call cem_drift_source
c        call cem_drift_poisson ! unnecessary
        if (ifsteric) then
          call cem_drift_sem_bdf1_steric
        else 
          call cem_drift_sem_bdf1
        endif

C       Richarson extrapolation: cN = 2* r2half_cN - bdf(cN,dt)
        do ic = 1,lcdim
          call cmult(r2_cN(1,ic),2.0,npts)
          call add2s1(cN(1,ic),r2_cN(1,ic),-1.0,npts)
        enddo
c        call cmult(r2_po,2.0,npts)
c        call add2s1(potent,r2_po,-1.0,npts)

      endif

      if (iffilter) call q_filter(0.01)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'    

      real dummy(lx1*ly1*lz1*lelt)

      call usersrc (time,rhs,rescN,dummy,dummy,dummy,dummy)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'

      real wa(lx1*ly1*lz1*lelt)
      real wb(lx1*ly1*lz1*lelt)

      integer icalled
      save icalled
      data icalled/0/

      integer isd,imsh,maxit
      real    tol

c...  setup helmholtz: [h1 A + h2 B] u = rhs where u = u_0+u_b
      if (icalled.eq.0) then
        call copy (mask,dmask(1,1),npts)
        call copy (h1,d_permit,npts)
        icalled = 1
      endif
      call rzero (h2,npts)

      isd = 1
      imsh = 1
      maxit = 5000
      tol = param(22)

c...  setup helmholtz: [h1 A + h2 B] u_0 = rhs - h1 A u_b
      call bcneusc (wa,npts,1,-1)    
      call col2 (wa,d_permit,npts)
      call add2 (h2,wa,npts)
      call bcdirsc (potent,npts,1)    
      call axhelm3 (wa,potent,h1,h2,1,1) 
      call sub3 (wb,rhs,wa,npts) ! wb = rhs-wa       
      call bcneusc (wa,npts,1,1)    
      call col2 (wa,d_permit,npts)
      call add2 (wb,wa,npts) 

c...  solve u_0 after assemble and mask for boundary on right hand side  
      call dssum (wb,lx1,ly1,lz1)
      call col2 (wb,mask,npts)

      if (ifsemg) then
        call copy      (wa,wb,npts)
        call hmh_gmres_mg (wa,h1,h2,mult,mask,maxit)
      else
        call cggo2 (wa,wb,h1,h2,mult,mask,imsh,isd,npts,tol,maxit)
      endif


c...  add u = u_0 + u_b (u_0= potent, u_b=potent_b if applicable)  
      call add2(potent,wa,npts)

c...  Compute electric fields cEE=-\grad (potent)
      call cem_drift_grad(gradp(1,1),gradp(1,2),gradp(1,3),potent)

      call chsign(gradp(1,1),npts)
      call chsign(gradp(1,2),npts)
      call chsign(gradp(1,3),npts)

      return
      end
c----------------------------------------------------------------------
      subroutine cem_drift_grad(w1,w2,w3,u1)
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real w1(1),w2(1),w3(1),u1(1)
      real rxmnk,sxmnk,txmnk,rymnk,symnk,tymnk,rzmnk,szmnk,tzmnk
      integer e,k,i,j,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real u1r,u1s,u1t

      nn = lx1-1

      if (if3d) then

        do e = 1,nelt

          j = nxyz*(e-1)
          call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

          do i = 1,nxyz

            k = i+j
            rxmnk= rxmn(k)/jacm(k)
            sxmnk= sxmn(k)/jacm(k)
            txmnk= txmn(k)/jacm(k)
            rymnk= rymn(k)/jacm(k)
            symnk= symn(k)/jacm(k)
            tymnk= tymn(k)/jacm(k)
            rzmnk= rzmn(k)/jacm(k)
            szmnk= szmn(k)/jacm(k)
            tzmnk= tzmn(k)/jacm(k)
            w1(k)= u1r(i)*rxmnk+u1s(i)*sxmnk+u1t(i)*txmnk
            w2(k)= u1r(i)*rymnk+u1s(i)*symnk+u1t(i)*tymnk
            w3(k)= u1r(i)*rzmnk+u1s(i)*szmnk+u1t(i)*tzmnk

          enddo

        enddo

      else
         do e = 1,nelt

           j= nxyz*(e-1)
           call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

           do i = 1,nxyz

             k = i+j
             rxmnk= rxmn(k)/jacm(k)
             sxmnk= sxmn(k)/jacm(k)
             rymnk= rymn(k)/jacm(k)
             symnk= symn(k)/jacm(k)
             w1(k)= u1r(i)*rxmnk+u1s(i)*sxmnk
             w2(k)= u1r(i)*rymnk+u1s(i)*symnk
             w3(k) =0.0

          enddo
          enddo

      endif

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'

      integer ic
      npts= lx1*ly1*lz1*nelt 

      call cem_drift_sem_bdf1_init

      do ic = 1,lcdim
        call cem_drift_rhs_bdf1(ic)
        call cem_drift_lhs_bdf1(ic)
      enddo

      return
      end

c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_init
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'

      integer ic 

      npts= lx1*ly1*lz1*nelt 
                             
      do ic = 1,lcdim
        call cmult2 (h1c(1,ic),diff_n(1,ic),dt,npts) ! h1n= h1n * (diff_n*dt)
        call rone (h2c(1,ic),npts)             ! h2n= 1.0
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_lhs_bdf1(ic)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'

      integer  isd,imsh,maxit
      real     tol

      integer ic
      real wa (lx1*ly1*lz1*lelt)
      real wb (lx1*ly1*lz1*lelt)

      call bcneusc0 (rescn(1,ic),mu_n(1,ic),cn(1,ic),npts,ic+1,1)

c...  solve cn
      
      call bcneusc (wa,npts,ic+1,-1)
      call col2 (wa,h1c(1,ic),npts)
      call add2 (h2c(1,ic),wa,npts)
      call bcdirsc (cn(1,ic),npts,ic+1)
      call axhelm3 (wa,cn(1,ic),h1c(1,ic),h2c(1,ic),1,1)
      call sub3 (wb,rescN(1,ic),wa,npts) 
      call bcneusc (wa,npts,ic+1,1)
      call col2 (wa,h1c(1,ic),npts) 
      call add2 (wb,wa,npts) 

c...  set cggo parameters

      isd = 1
      imsh = 1
      maxit = 5000
      tol = param(22)

      call dssum (wb,lx1,ly1,lz1)
      call col2 (wb,dmask(1,ic+1),npts)
      call cggo2(wa,wb,h1c(1,ic),h2c(1,ic)
     $     ,mult,dmask(1,ic+1),imsh,isd,npts,tol,maxit)
      call add2 (cn(1,ic),wa,npts)

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_rhs_bdf1(ic)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'
      real wn(lpts), wp(lpts)
      real w_h1(lpts),w_h2(lpts)
      integer ic

      call rzero(w_h2,npts)

c      do ic = 1,lcdim
        call col3 (wn,cn(1,ic),bm1,npts)
        call add2s1 (rescN(1,ic),wn,dt,npts) ! resN=cN+dt*(resN) with weights

c...  resN=resN -zN wN: (wN= -cn\grad\Phi,\grad v), h1= mu_n*cN, h2=0
c...  resP=resP -zP wP: (wP=  cp\grad\Phi,\grad v), h1= mu_p*cP, h2=0
c...  compute the nonlinear terms: cN*(grad(phi),cP*(grad(phi)
c...  and add it to the right-hand-side
        call col3(w_h1,cn(1,ic),mu_n(1,ic),npts) ! wk_h1n= mu_n * cN
        call axhelm3(wn,potent,w_h1,w_h2,1,1)
        call add2s2 (rescN(1,ic),wn,-1.0*zvalence(ic)*dt,npts)
c      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_div(wk0,wkx,wky,wkz)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      real    wkx(1),wky(1),wkz(1),wk0(1)

      if (if3d) then
      call drift_grad(wk1,wk2,wk3,wkx)
      call copy(wk0,wk1,npts)
      call drift_grad(wk1,wk2,wk3,wky)
      call add2s2(wk0,wk2,1.0,npts)
      call drift_grad(wk1,wk2,wk3,wky)
      call add2s2(wk0,wk3,1.0,npts)
      else
      call drift_grad(wk1,wk2,wk3,wkx)
      call copy(wk0,wk1,npts)
      call drift_grad(wk1,wk2,wk3,wky)
      call add2s2(wk0,wk2,1.0,npts)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine drift_wght_div_d(w0,w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DEALIAS'
c     ddddd
c     Combined weighted div:   w0 = B*div(u1,u2,u3) on a finer mesh
c     with dealiasing
c
      real     w0(1)
      real     u1(1),u2(1),u3(1)
      real     w1(1),w2(1),w3(1)

      integer    ldd,md
      parameter (ldd = lxd**ldim  )
      parameter ( md = lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)

      real         ud1r,ud1s,ud1t
     $           , ud2r,ud2s,ud2t
     $           , ud3r,ud3s,ud3t

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw

      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd0(ldd),wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd0,wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)
      real            dxmd,dxtmd,w3md,wgld,zmd

      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy, nn, mm, mx
      save    icalld
      data    icalld /0/

      mx    = lx1
      nn    = lx1-1
      mm    = md -1
      nxyd  = md**ndim
      nxy   = mx**ndim

      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

          if (if3d) then
             l = 0
             do k=1,md
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)*wgld(k)
             enddo
             enddo
             enddo
          else
             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)
             enddo
             enddo
          endif

          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)

      endif

      do e=1,nelt

         j  = nxyd *(e-1)
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u1(j1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,u2(j1),mx,md,if3d, 0)
         call intp_rstd (ju3,u3(j1),mx,md,if3d, 0)

         if (if3d) then

            call local_grad3 (ud1r,ud1s,ud1t,ju1,mm,1,dxmd,dxtmd)
            call local_grad3 (ud2r,ud2s,ud2t,ju2,mm,1,dxmd,dxtmd)
            call local_grad3 (ud3r,ud3s,ud3t,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j

            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u1tw = ud1t(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)
            u2tw = ud2t(i)*w3md(i)
            u3rw = ud3r(i)*w3md(i)
            u3sw = ud3s(i)*w3md(i)
            u3tw = ud3t(i)*w3md(i)

            wd1(i) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k)
     $               + u1tw*txmd(k))

            wd2(i) = ( u2rw*rymd(k)
     $               + u2sw*symd(k)
     $               + u2tw*tymd(k))

            wd3(i) = ( u3rw*rzmd(k)
     $               + u3sw*szmd(k)
     $               + u3tw*tzmd(k))

            wd0(i) = wd1(i)+wd2(i)+wd3(i)

            enddo

        else

            call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
            call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
            call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j

            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)

            wd1(i) = (  u1rw*rxmd(k)
     $               + u1sw*sxmd(k))

            wd2(i) = (  u2rw*rymd(k)
     $               + u2sw*symd(k))

            wd0(i) = wd1(i)+wd2(i)

           enddo

         endif

         call intp_rstd (w0(j1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_div_grad_d(w0,w1,w2,w3,u1,u2,u3,u4)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DEALIAS'
c     ddddd
c     Combined weighted div(u4*(u1,u2,u3)^t):   w0 = B*div(u4*(u1,u2,u3)^t) on a finer mesh
c     with dealiasing
c
      real     w0(1)
      real     u1(1),u2(1),u3(1),u4(1)
      real     w1(1),w2(1),w3(1)

      integer    ldd,md
      parameter (ldd = lxd**ldim  )
      parameter ( md = lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)

      real         ud1r,ud1s,ud1t
     $           , ud2r,ud2s,ud2t
     $           , ud3r,ud3s,ud3t

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      real    ju4(ldd)
      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd0(ldd),wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd0,wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)
      real            dxmd,dxtmd,w3md,wgld,zmd

      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy, nn, mm, mx
      save    icalld
      data    icalld /0/

      mx    = lx1
      nn    = lx1-1
      mm    = md -1
      nxyd  = md**ndim
      nxy   = mx**ndim

      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

          if (if3d) then
             l = 0
             do k=1,md
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)*wgld(k)
             enddo
             enddo
             enddo
          else
             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)
             enddo
             enddo
          endif

          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)

      endif

      do e=1,nelt

         j  = nxyd *(e-1)
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u1(j1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,u2(j1),mx,md,if3d, 0)
         call intp_rstd (ju3,u3(j1),mx,md,if3d, 0)
         call intp_rstd (ju4,u4(j1),mx,md,if3d, 0)
         call col2(ju1,ju4,nxyd)
         call col2(ju2,ju4,nxyd)
         call col2(ju3,ju4,nxyd)

         if (if3d) then

            call local_grad3 (ud1r,ud1s,ud1t,ju1,mm,1,dxmd,dxtmd)
            call local_grad3 (ud2r,ud2s,ud2t,ju2,mm,1,dxmd,dxtmd)
            call local_grad3 (ud3r,ud3s,ud3t,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j

            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u1tw = ud1t(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)
            u2tw = ud2t(i)*w3md(i)
            u3rw = ud3r(i)*w3md(i)
            u3sw = ud3s(i)*w3md(i)
            u3tw = ud3t(i)*w3md(i)

            wd1(i) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k)
     $               + u1tw*txmd(k))

            wd2(i) = ( u2rw*rymd(k)
     $               + u2sw*symd(k)
     $               + u2tw*tymd(k))

            wd3(i) = ( u3rw*rzmd(k)
     $               + u3sw*szmd(k)
     $               + u3tw*tzmd(k))

            wd0(i) = wd1(i)+wd2(i)+wd3(i)

            enddo

        else

            call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
            call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
            call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

            do i=1,nxyd

            k = i+j

            u1rw = ud1r(i)*w3md(i)
            u1sw = ud1s(i)*w3md(i)
            u2rw = ud2r(i)*w3md(i)
            u2sw = ud2s(i)*w3md(i)

            wd1(i) = ( u1rw*rxmd(k)
     $               + u1sw*sxmd(k))

            wd2(i) = ( u2rw*rymd(k)
     $               + u2sw*symd(k))

            wd0(i) = wd1(i)+wd2(i)

           enddo

         endif

         call intp_rstd (w0(j1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_intp_u_grad_v(w0,u,v)
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DEALIAS'
c     ddddd
c     Combined weighted div(u4*(u1,u2,u3)^t):   w0 = B*div(u4*(u1,u2,u3)^t) on a finer mesh
c     with dealiasing
c
      real     w0(lx1,ly1,lz1,1)
      real     u(lx1,ly1,lz1,1)
      real     v(lx1,ly1,lz1,1)

      integer    ldd,md
      parameter (ldd = lxd**ldim  )
      parameter ( md = lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)

      real         ud1r,ud1s,ud1t
     $           , ud2r,ud2s,ud2t
     $           , ud3r,ud3s,ud3t

      common /t1t2/dudrd(lxd,lyd,lzd)
     $            ,dudsd(lxd,lyd,lzd)
     $            ,tmp1d(lxd,lyd,lzd)
     $            ,tmp2d(lxd,lyd,lzd)

      real         dudrd,dudsd,tmp1d,tmp2d
      real         tm1(lxd,lyd,lzd)
      real         tm2(lxd,lyd,lzd)

      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      real    ju4(ldd)
      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd0(ldd),wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd0,wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)
      real            dxmd,dxtmd,w3md,wgld,zmd
      common /ddtmp6/ g1m1d(lxd**ldim*lelt)
     $               ,g2m1d(lxd**ldim*lelt)
     $               ,g4m1d(lxd**ldim*lelt)
     $               ,jacmd(lxd**ldim*lelt)
     $               ,wjd(lxd**ldim*lelt)
      real            g1m1d,g2m1d,g4m1d,jacmd,wjd
      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy, nn, mm, mx,nptsd
      save    icalld
      data    icalld /0/
      EQUIVALENCE    (DUDRd,TM1),(DUDSd,TM2)
      mx    = lx1
      nn    = lx1-1
      mm    = md -1
      nxyd  = md**ndim
      nxy   = mx**ndim
      nptsd = nxyd*nelt
      if (icalld.eq.0) then

          icalld = 1

          call zwgll(zmd,wgld,md)

             l = 0
             do j=1,md
             do i=1,md
                l=l+1
                w3md(l) = wgld(i)*wgld(j)
             enddo
             enddo

          do e=1,nelt   ! this will need to be dynamic if we have p-refinement

             j1 = nxy *(e-1) + 1
             jd = nxyd*(e-1) + 1

             call intp_rstd (rxmd(jd),rxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (sxmd(jd),sxm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (txmd(jd),txm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rymd(jd),rym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (symd(jd),sym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tymd(jd),tym1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (rzmd(jd),rzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (szmd(jd),szm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (tzmd(jd),tzm1(j1,1,1,1),mx,md,if3d, 0)
             call intp_rstd (jacmd(jd),jacm1(j1,1,1,1),mx,md,if3d, 0)
          CALL INVERS2 (WJd(jd),JACMd(jd),nxyd)
         CALL VDOT2(G1M1d(jd),RXMd(jd),RYMd(jd),RXMd(jd),RYMd(jd),nxyd)
         CALL VDOT2(G2M1d(jd),sXMd(jd),sYMd(jd),sXMd(jd),sYMd(jd),nxyd)
         CALL VDOT2(G4M1d(jd),RXMd(jd),RYMd(jd),SXMd(jd),SYMd(jd),nxyd)
          CALL COL2  (G1M1d(jd),WJd(jd),nxyd)
          CALL COL2  (G2M1d(jd),WJd(jd),nxyd)
          CALL COL2  (G4M1d(jd),WJd(jd),nxyd)
          CALL COL2 (G1M1d(jd),W3Md,nxyd)
          CALL COL2 (G2M1d(jd),W3Md,nxyd)
          CALL COL2 (G4M1d(jd),W3Md,nxyd)

          enddo

          call gen_dgl(dxmd,dxtmd,md,md,wgld)
      endif


      call rzero(w0,npts)
      do e=1,nelt
         call rzero(wd0,nxyd)
         j  = nxyd *(e-1) + 1
         j1 = nxy  *(e-1) + 1

         call intp_rstd (ju1,u(j1,1,1,1),mx,md,if3d, 0)  ! 0 ==> forward
         call intp_rstd (ju2,v(j1,1,1,1),mx,md,if3d, 0)

         CALL MXM  (DXMd,md,ju1,md,DUDRd,md)
         CALL MXM  (ju1,md,DxTMd,md,DUDSd,md)
         CALL COL3 (TMP1d,DUDRd,G1M1d(j),nxyd)
         CALL COL3 (TMP2d,DUDSd,G2M1d(j),nxyd)

         CALL ADDCOL3 (TMP1d,DUDSd,G4M1d(j),NXYd)
         CALL ADDCOL3 (TMP2d,DUDRd,G4M1d(j),NXYd)
         CALL COL2 (TMP1d,ju2,NXYd)
         CALL COL2 (TMP2d,ju2,NXYd)
         CALL MXM  (DXTMd,md,TMP1d,md,TM1,md)
         CALL MXM  (TMP2d,md,DxMd,md,TM2,md)
         CALL ADD2 (wd0,TM1,NXYd)
         CALL ADD2 (wd0,TM2,NXYd)


         call intp_rstd (w0(j1,1,1,1),wd0,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres2(phi,res,h1,h2,wt,mask,isd,imsh,n,tol)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,outer,isd,imsh
      real     phi(lpts),res(lpts),h1(lpts),h2(lpts),wt(lpts),mask(lpts)
      real     tol,alpha,l,temp
      real*8   etime1,dnekclock

      if (nid.eq.0) write(6,*) 'start: hmh_gmres2'

      n = lx1*ly1*lz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      call rzero(x,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.5000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r,res,n)                  ! r = res
         else
            !update residual
            call copy   (r,res,n)                  ! r = res
            call axhelm3(w,x,h1,h2,imsh,isd)       ! w = A x
            call dssum  (w,lx1,ly1,lz1)
            call col2   (w,mask,n)
            call add2s2 (r,w,-1.,n)                ! r = r - w
         endif

         gamma(1) = glsc3(r,r,wt,n)                ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call axhelm3(w,v(1,j),h1,h2,imsh,isd)  ! w = A v
            call dssum  (w,lx1,ly1,lz1)
            call col2   (w,mask,n)

c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i

               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))!*norm_fac
c           ratio = rnorm/div0
            if ((nid.eq.0).and.(istep.le.5))
     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' Divergence')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,v(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

       call copy(phi,x,n)

c      call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nid.eq.0) then
          if ((mod(istep,iocomm).eq.0).or.(istep.le.10))
     $    write(6,9999) istep,iter,tolpss,etime_p,etime1
      endif
 9999 format(' ',' ',i9,' gmres   : iteration#',i5,1p3e12.4)

      return
      end
c=======================================================================
      subroutine fast_poisson_test(mode)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'ZPER'  ! lelx,lely,lelz in nek5000

      common /rwork/ ff(lpts),uu(lpts),wwk(lpts),wwk1(lpts),wwk2(lpts)
      real    ff,uu,wwk,wwk1,wwk2,wmax,glamax
      real    kwave2
      real    pi2,xx,yy
      integer mode,i

      nelx = param(116)+0.001
      nely = param(117)+0.001
      nelz = param(118)+0.001

      if (nelx.gt.lelx) call exitti('lelx too small:$',lelx,1)
      if (nely.gt.lely) call exitti('lely too small:$',lely,1)
      if (nelz.gt.lelz) call exitti('lelz too small:$',lelz,1)


      pi2 = pi/2.
      kwave2 =0.0 !100 !0.0 ! pi*pi

      call gfdm_init(lx1,ly1,lz1,.false.,0.0)

      if (mode.eq.1) then
        do i=1,npts
         xx = xm1(i,1,1,1) !*i
         yy = ym1(i,1,1,1) !*i
         ff(i) = pi*pi2*cos(pi2*xx)*cos(pi2*yy)
     $         - kwave2*cos(pi2*xx)*cos(pi2*yy)
         uu(i) = cos(pi2*xx)*cos(pi2*yy)
         enddo
      else
        do i=1,npts
         xx = xm1(i,1,1,1) !*i
         yy = ym1(i,1,1,1) !*i
         ff(i) = pi*pi2*sin(pi2*xx)*sin(pi2*yy)
     $         - kwave2*sin(pi2*xx)*sin(pi2*yy)
         uu(i) = sin(pi2*xx)*sin(pi2*yy)
         enddo
      endif

      call col2(ff,bm1,npts)

      call gfdm_pres_solv(wwk,ff,wwk1,wwk2,.false.,0.0) ! (A - kwave2*B)z = r

      do i=1,npts
         write(10,10) i,uu(i),wwk(i)
      enddo
  10  format(i8,2e25.15)

      call sub3(wwk1,uu,wwk,npts)

      wmax=glamax(wwk1,npts)

      if (nid.eq.0) write(6,*) 'gdfm:: maxerr',wmax

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bcs_fc_ptr
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'BCS'
      integer  ifld,i,i0,i1,j,l1,l2,l3,e,f
      character CB*3
      integer iglsum
      integer glo_ncemface_pec
      integer glo_ncemface_nmn
      integer glo_ncemface_rob

      i0= 1
      i1= nfield

      do ifld= i0,i1
         l1=0
         l2=0
         l3=0
      do e= 1,nelt
      do f= 1,nfaces
      do j= 1,nxzf
         CB= CBC(f,e,ifld)
         i =(e-1)*nfaces*nxzf+(f-1)*nxzf+j !global face numbering
         if  (CB.eq.'N  ') then
            l1= l1+1
            cemface_nmn(l1,ifld)= i
         elseif  (CB.eq.'PEC') then
            l2= l2+1
            cemface_pec(l2,ifld)= i
         elseif  (CB.eq.'R  ') then
            l3= l3+1
            cemface_rob(l3,ifld)= i
         endif
      enddo
      enddo
      enddo

      ncemface_nmn(ifld)  = l1
      ncemface_pec(ifld)  = l2
      ncemface_rob(ifld)  = l3
      enddo

      do ifld=i0,i1    

        glo_ncemface_pec= iglsum(ncemface_pec(ifld),1)
        glo_ncemface_nmn= iglsum(ncemface_nmn(ifld),1)
        glo_ncemface_rob= iglsum(ncemface_rob(ifld),1)

        if (nid.eq.0) then
         write(6,*) 'tot#pts: BCs ncemface_pec',ifld,glo_ncemface_pec
         write(6,*) 'tot#pts: BCs ncemface_nmn',ifld,glo_ncemface_nmn
         write(6,*) 'tot#pts: BCs ncemface_rob',ifld,glo_ncemface_rob
        endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bd_nmn(bdin,bdout,ifld)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'BCS'

      real    bdin(1),bdout(1)
      integer ifld,i,i0,j

      do i=1,ncemface_nmn(ifld)
         j=cemface_nmn(i,ifld)
         i0=cemface(j)
         bdout(i0)=bdout(i0)+bdin(i)*aream(j)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bd_pec(bdin,bdout,ifld)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'BCS'
      real    bdin(1),bdout(1)
      integer ifld,i,i0,j

      do i=1,ncemface_pec(ifld)
         j =cemface_pec(i,ifld)
         i0=cemface(j)
         bdout(i0)=bdout(i0)+bdin(i)*aream(j)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_bd_rob  (bdin,bdout,ifld)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'BCS'
      real    bdin(1),bdout(1)
      integer ifld,i,i0,j

      do i=1,ncemface_rob (ifld)
         j =cemface_rob (i,ifld)
         i0=cemface(j)
         bdout(i0)=bdout(i0)+bdin(i)*aream(j)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_scaling(flag)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      integer flag
      integer ic
      if (flag.eq.1) then
         do ic=1,lcdim
           call cmult(cN(1,ic),scale_concent,npts)
         enddo
         call cmult(potent,scale_potent,npts)
         call cmult(gradp,scale_potent,lpts1*3)
      else

      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine setdmask (dmask,n,ifld)
c-----------------------------------------------------------------------
C     set dmask for poisson problem
      include 'SIZE'
      include 'TOTAL'
      integer   e,f,i,j,n,ifld
      real      dmask(1)
      character CB*3

      do e=1,nelt
      do f=1,nfaces

         CB =CBC(f,e,ifld)
c        Homogeneous Dirichlet boundary conditions
         if (CB.eq.'PEC') call facev (dmask,e,f,0.0,lx1,ly1,lz1)
         if (CB.eq.'PML') call facev (dmask,e,f,0.0,lx1,ly1,lz1) ! fischer/mmin: 1/23/07

      enddo
      enddo

      call dsop(dmask,'MUL',lx1,ly1,lz1)

c     if (nid.eq.0) write(6,*) ' done: setdmaskp'

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_op_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'

      if (nid.eq.0) write(6,*) 'compute eigenvalues'

c     call cem_drift_eig        !call drift_eig_arpack
c     call cem_end

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_grad(w1,w2,w3,u1)
c---------------------------------------------------------------------
      implicit none
c     Non-weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real    u1r,u1s,u1t

      nn    = lx1-1

      if (if3d) then

         do e=1,nelt
         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (    u1r(i)*rxmn(k)
     $                 + u1s(i)*sxmn(k)
     $                 + u1t(i)*txmn(k))/jacm(k)

            w2(k) = (    u1r(i)*rymn(k)
     $                 + u1s(i)*symn(k)
     $                 + u1t(i)*tymn(k))/jacm(k)

            w3(k) = (    u1r(i)*rzmn(k)
     $                 + u1s(i)*szmn(k)
     $                 + u1t(i)*tzmn(k))/jacm(k)

          enddo

          enddo

       else

c        ...2d grad
          do e= 1,nelt

           j= nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j

             w1(k) = ( u1r(i)*rxmn(k)
     $               + u1s(i)*sxmn(k))/jacm(k)

             w2(k) =(  u1r(i)*rymn(k)
     $               + u1s(i)*symn(k))/jacm(k)
             w3(k)= 0.0
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_wght_grad(w1,w2,w3,u1)
c-----------------------------------------------------------------------
      implicit none
c     Weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
 
      integer mode        
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real    u1r,u1s,u1t
c
      nn    = lx1-1

      if (ldim.eq.3) then 

         do e=1,nelt
         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
c...        rym1,sym1,tym1,...jac are multiplied. w3m1 is purely weight.
c...        area(i,1,1,1) is multiplied with jacobina weight.
            u1rw = u1r(i)*w3mn(i)
            u1sw = u1s(i)*w3mn(i)
            u1tw = u1t(i)*w3mn(i)
            w1(k) =      u1rw*rxmn(k)
     $                 + u1sw*sxmn(k)
     $                 + u1tw*txmn(k)

            w2(k) =      u1rw*rymn(k)
     $                 + u1sw*symn(k)
     $                 + u1tw*tymn(k)

            w3(k) =      u1rw*rzmn(k)
     $                 + u1sw*szmn(k)
     $                 + u1tw*tzmn(k)

          enddo

          enddo

       else

c        ...2d grad
          do e= 1,nelt

           j= nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j


             u1rw = u1r(i)*w3mn(i)
             u1sw = u1s(i)*w3mn(i)

c             u1rw = u1r(i) !*w3mn(i)
c             u1sw = u1s(i) !*w3mn(i)

             w1(k) = ( u1rw*rxmn(k)
     $               + u1sw*sxmn(k))!/jacm(k)

             w2(k) =(  u1rw*rymn(k)
     $               + u1sw*symn(k))!/jacm(k)
             w3(k)=w2(k)
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine drift_wght_div (w0,w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     div:   w = div(u)/jacm
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer  k,i,j,e,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real    ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw

      nn    = lx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(ur1,us1,ut1,u1,nn,e,dxm1,dxtm1)
         call local_grad3(ur2,us2,ut2,u2,nn,e,dxm1,dxtm1)
         call local_grad3(ur3,us3,ut3,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            u1rw = ur1(i)*w3mn(i)
            u1sw = us1(i)*w3mn(i)
            u1tw = ut1(i)*w3mn(i)
            u2rw = ur2(i)*w3mn(i)
            u2sw = us2(i)*w3mn(i)
            u2tw = ut2(i)*w3mn(i)
            u3rw = ur3(i)*w3mn(i)
            u3sw = us3(i)*w3mn(i)
            u3tw = ut3(i)*w3mn(i)

            w1(k) = (  u1rw*rxmn(k)
     $               + u1sw*sxmn(k)
     $               + u1tw*txmn(k))

            w2(k) = (  u2rw*rymn(k)
     $               + u2sw*symn(k)
     $               + u2tw*tymn(k))

            w3(k) = (  u3rw*rzmn(k)
     $               + u3sw*szmn(k)
     $               + u3tw*tzmn(k))

            w0(k) = w1(k)+w2(k)+w3(k)
         enddo

         enddo

      else

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
         call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            u1rw = ur1(i)*w3mn(i)
            u1sw = us1(i)*w3mn(i)
            u2rw = ur2(i)*w3mn(i)
            u2sw = us2(i)*w3mn(i)

            w1(k) = (  u1rw*rxmn(k)
     $               + u1sw*sxmn(k))

            w2(k) = (  u2rw*rymn(k)
     $               + u2sw*symn(k))

            w0(k) = w1(k)+w2(k)
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bcdirsc (wa,n,ifld)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'   

      integer n,ifld,itype
      real wa(lpts)
      real wtmp(lx1*ly1*lz1*lelt)

c.... fixme: 7/16/17 misun, this is still inactive

      call rzero(wtmp,n)

      if (ncemface_pec(ifld).gt.0) then
c       do i = 1,ncemface_pec(ifld)
c         j = cemface_pec(i,ifld)
c         k = cemface(j)
c         wtmp(k) = wtmp(k) + wa(k)*bm1(k,1,1,1)
c       enddo
c       call col2(wa,dmask(1,ifld),n)
c       call add2(wa,wtmp,n)                
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bcneusc (wa,n,ifld,itype) 
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'    

      integer n,ifld,itype
      real wa(lpts)

      call rzero(wa,n)

      if (itype.eq.-1) then

        if (ncemface_rob(ifld).gt.0) then
          do i = 1,ncemface_rob(ifld)
            j = cemface_rob(i,ifld)
            k = cemface(j)
            wa(k) = wa(k) + bc_conv(j,ifld)*area(j,1,1,1)/bm1(k,1,1,1)
          enddo
        endif

      else

        if (ncemface_rob(ifld).gt.0) then
          do i = 1,ncemface_rob(ifld)
            j = cemface_rob(i,ifld)
            k = cemface(j)
            wa(k) = wa(k) + bc_flux(j,ifld)*area(j,1,1,1)      
          enddo
        endif

        if (ncemface_nmn(ifld).gt.0) then
          do i = 1,ncemface_nmn(ifld)
            j = cemface_nmn(i,ifld)
            k = cemface(j)
            wa(k) = wa(k) + bc_flux(j,ifld)*area(j,1,1,1)
          enddo
        endif

      endif
      
      return
      end
c-----------------------------------------------------------------------
c     this should be aready included in bc_flux as current flux 
c     this routine is tempory
      subroutine bcneusc0 (wa,wb,wc,n,ifld,itype) 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'
      integer n,ifld,ifld_c,itype
      real wa(lpts),wb(lpts),wc(lpts)
      real tmp
      integer i,j,k

      ifld_c = ifld-1

      if (itype.eq.1) then

        if (ncemface_rob(ifld).gt.0) then
          do i = 1,ncemface_rob(ifld)
            j = cemface_rob(i,ifld)
            k = cemface(j)
            tmp = gradp(k,1)*unxm(j)+gradp(k,2)*unym(j)
     $          + gradp(k,3)*unzm(j)
            wa(k) = wa(k) 
     $            - tmp*zvalence(ifld_c)*dt*wc(k)*wb(k)*aream(j)
          enddo
        endif

        if (ncemface_nmn(ifld).gt.0) then
          do i = 1,ncemface_nmn(ifld)
            j = cemface_nmn(i,ifld)
            k = cemface(j)
            tmp = gradp(k,1)*unxm(j)+gradp(k,2)*unym(j)
     $          + gradp(k,3)*unzm(j)
            wa(k) = wa(k) 
     $            - tmp*zvalence(ifld_c)*dt*wc(k)*wb(k)*aream(j)
          enddo
        endif

      endif
      
      return
      end
c-----------------------------------------------------------------------
