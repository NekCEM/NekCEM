c-----------------------------------------------------------------------
      subroutine gfdm_pres_solv(z,r,ug,wg,ifemat,kwave2) ! (A - kwave2*B)z = r

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      include 'ZPER'

      real z (lx2*ly2*lz2*lelv),r (lx2*ly2*lz2*lelv)
      real ug(lx2*ly2*lz2*lelv),wg(lx2*ly2*lz2*lelv)
      real kwave2
      logical ifemat

      mfld = 1
      if (ifmhd.and.ifield.gt.1) mfld = 2

c     Set up diagonal for solver

      ii= 1          ! we'll need to change this for multiple flds
      l = ngfdm_p(1)
      m = ngfdm_p(2)
      n = ngfdm_p(3)
      i = mlp(1,mfld)
      j = mlp(2,mfld)
      k = mlp(3,mfld)

      if (nid.eq.0) write(6,*) 'gfdm_pres_solv: mfld:',mfld
      if (nid.eq.0) write(6,*) 'l,m,n,i,j,k:',l,m,n,i,j,k
       
      call gfdm_set_diagp(wavep(ii),tpn1,mcex
     $         ,eigp(i),l,eigp(j),m,eigp(k),n,ifemat,kwave2)

      mp   = ngfdm_p(pst2lex(1))
      ms   = ngfdm_p(pst2lex(2))
      mt   = ngfdm_p(pst2lex(3))

      msfp = msp(pst2lex(1),mfld)
      msfs = msp(pst2lex(2),mfld)
      msft = msp(pst2lex(3),mfld)

      ntot  = nx2*ny2*nz2*nelv

      m     = ntot/mp
      nwave = mcex
      mpt   = mcex/(ms*mt)

      dtbdi = 1   
      if (ifemat.and..not.ifmhd) dtbdi=bd(1)/dt ! scale by BDF coefficient/DT

c     write(6,*) 'msfp:',msfp,msfs,msft,mp,ms,mt,sp(msfp)
      call outmat(sp(msfp),mp,mp,'S  x  ',istep)
      call outmat(sp(msfs),ms,ms,'S  y  ',ifield)
      call outmat(sp(msft),mt,mt,'S  z  ',mt)
      call outmat(wavep   ,mp,ms,'WAVE  ',mt)
c     stop

      if (if3d) then
         call copy    (ug,r,ntot)
c        write(6,*) nid,'SWAP_IP 1:',ntot
         call swap_ip (ug,tpn2,ntot)
         call mxm     (ug,m,sp (msfp),mp,wg,mp)
         call cexr    (ug,wg,m,mp,part_out,part_in,msg_id,wdsize,nid,np)
c        write(6,*) nid,'SWAP_IP 2:',ntot
         call swap_ip (ug,ind23,mcex)
         call mxm     (spt(msfs),ms,ug,ms,wg,mpt*mt)
         call mxm     (wg,ms*mpt,sp (msft),mt,ug,mt)
         call col2s2  (ug,wavep,dtbdi,mcex)
         call mxm     (ug,ms*mpt,spt(msft) ,mt,wg,mt)
         call mxm     (sp (msfs),ms,wg,ms,ug,mpt*mt)
         call swapt_ip(ug,ind23,mcex)
         call cextr   (wg,m,mp,ug,part_out,part_in,msg_id,wdsize,nid,np)
         call mxm     (wg,m,spt(msfp),mp,z,mp)
         call swapt_ip(z,tpn2,ntot)
      else
         call copy    (ug,r,ntot)
         call swap_ip (ug,tpn2,ntot)
         call mxm     (ug,m,sp (msfp),mp,wg,mp)
         call cexr    (ug,wg,m,mp,part_out,part_in,msg_id,wdsize,nid,np)
         call swap_ip (ug,ind23,mcex)
         call mxm     (spt(msfs),ms,ug,ms,wg,mpt)  ! wg= St*r*S
         call col2s2  (wg,wavep,dtbdi,mcex)        ! wg= wg*wavep*dtbdi 
         call mxm     (sp (msfs),ms,wg,ms,ug,mpt)  ! ug= S*wg 
         call swapt_ip(ug,ind23,mcex)
         call cextr   (wg,m,mp,ug,part_out,part_in,msg_id,wdsize,nid,np)
         call mxm     (wg,m,spt(msfp),mp,z,mp)     ! z = wg* St
         call swapt_ip(z,tpn2,ntot)
      endif

      return
      end
c-----------------------------------------------------------------------
