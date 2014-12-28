c-----------------------------------------------------------------------
c    BEGIN  Projection Method
c-----------------------------------------------------------------------
      subroutine hmhzpf(name,u,r,h1,h2,mask,mult,imesh,tli,maxit,isd,bi)
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'FDMH1'
      include 'CTIMER'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           MASK (LX1,LY1,LZ1,1)
      REAL           MULT (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      COMMON /CTMP0/ W1   (LX1,LY1,LZ1,LELT)
     $ ,             W2   (LX1,LY1,LZ1,LELT)
c
      etime1=dnekclock()
 
      npts = nx1*ny1*nz1*nelt

c     IF (IMESH.EQ.1) NTOT = NX1*NY1*NZ1*NELV
c     IF (IMESH.EQ.2) NTOT = NX1*NY1*NZ1*NELT
 
      npts = nx1*ny1*nz1*nelt
      tol  = param(22) 

c     CALL CHKTCG1 (TOL,R,H1,H2,MASK,MULT,IMESH,ISD)
 
 
c     Set flags for overlapping Schwarz preconditioner (pff 11/12/98)
c
                          kfldfdm = -1
c     if (name.eq.'TEMP') kfldfdm =  0
c     if (name.eq.'VELX') kfldfdm =  1
c     if (name.eq.'VELY') kfldfdm =  2
c     if (name.eq.'VELZ') kfldfdm =  3
      if (name.eq.'PRES') kfldfdm =  ndim+1
c
      call cggo
     $      (u,r,h1,h2,mask,mult,imesh,tol,maxit,isd,bi,name)
      thmhz=thmhz+(dnekclock()-etime1)
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dd_hsolve(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd
     $                     ,approx,napprox,bi)
c
c     Either std. Helmholtz solve, or a projection + Helmholtz solve
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           vmk  (LX1,LY1,LZ1,1)
      REAL           vml  (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      REAL           approx (1)
      integer        napprox(1)
c     common /ctmp2/ pp1   (lx1,ly1,lz1,lelt)
c     common /ctmp3/ pp2   (2+2*mxprev)

      logical ifstdh
      character*4  cname
      character*6  name6

      logical ifwt,ifvec

      call chcopy(cname,name,4)
      call capit (cname,4)
      call blank (name6,6)
      call chcopy(name6,name,4)

      ifwt  = .true.
      ifvec = .false.

      call col2   (r,vmk,npts)  ! mask
      call dssum  (r,nx1,ny1,nz1)

      if (nid.eq.0) write(6,*) 'call project1',ifield,ifwt,ifvec

      call project1(r,n,approx,napprox,h1,h2,vmk,vml,ifwt,ifvec,name6)

      call cggo (u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)

      call project2(u,r,n,approx,napprox,h1,h2,vmk,vml,ifwt,ifvec,name6)

      return
      end
c-----------------------------------------------------------------------
      subroutine hsolve(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd
     $                 ,approx,napprox,bi)
c
c     Either std. Helmholtz solve, or a projection + Helmholtz solve
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           vmk  (LX1,LY1,LZ1,1)
      REAL           vml  (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      REAL           approx (1)
      integer        napprox(1)
c     common /ctmp2/ pp1   (lx1,ly1,lz1,lelt)
c     common /ctmp3/ pp2   (2+2*mxprev)

      logical ifstdh
      character*4  cname
      character*6  name6

      logical ifwt,ifvec

      call chcopy(cname,name,4)
      call capit (cname,4)


      p945 = param(94)
      if (cname.eq.'PRES') p945 = param(95)

                          ifstdh = .false.
      if (param(93).eq.0) ifstdh = .true.
      if (p945.eq.0)      ifstdh = .true.
      if (istep.lt.p945)  ifstdh = .true.

      if (ifstdh) then
         call hmholtz(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd)
      else

         n = nx1*ny1*nz1*nelfld(ifield)

         call col2   (r,vmk,n)
         call dssum  (r,nx1,ny1,nz1)

c        call projh  (r,h1,h2,bi,vml,vmk,approx,napprox,pp1,pp2,name)
c        call hmhzpf (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)
c        call gensh  (u,h1,h2,vml,vmk,approx,napprox,pp1,pp2,name)

         call blank (name6,6)
         call chcopy(name6,name,4)
         ifwt  = .true.
         ifvec = .false.

         write(6,*) 'call project1',ifield,ifwt,ifvec
         call project1
     $       (r,n,approx,napprox,h1,h2,vmk,vml,ifwt,ifvec,name6)

         call hmhzpf (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)

         call project2
     $       (u,r,n,approx,napprox,h1,h2,vmk,vml,ifwt,ifvec,name6)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine project1(b,n,rvar,ivar,h1,h2,msk,w,ifwt,ifvec,name6)

c     1. Compute the projection of x onto X

c     2. Re-orthogonalize the X basis set and corresponding B=A*X
c        vectors if A has changed.

c     Output:  b = b - projection of b onto B

c     Input:   n     = length of field (or multifields, when ifvec=true)
c              rvar  = real array of field values, including old h1,h2, etc.
c              ivar  = integer array of pointers, etc.
c              h1    = current h1, for Axhelm(.,.,h1,h2,...)
c              h2    = current h2
c              msk   = mask for Dirichlet BCs
c              w     = weight for inner products (typ. w=vmult, tmult, etc.)
c              ifwt  = use weighted inner products when ifwt=.true.
c              ifvec = are x and b vectors, or scalar fields?
c              name6 = discriminator for action of A*x

c     The idea here is to have one pair of projection routines for 
c     constructing the new rhs (project1) and reconstructing the new
c     solution (x = xbar + dx) plus updating the approximation space.
c     The latter functions are done in project2.
c
c     The approximation space X and corresponding right-hand sides,
c     B := A*X are stored in rvar, as well as h1old and h2old and a
c     couple of other auxiliary arrays.

c     In this new code, we retain both X and B=A*X and we re-orthogonalize
c     at each timestep (with no extra matrix-vector products, but O(n m^2)
c     work.   The idea is to retain fresh vectors by injecting the most 
c     recent solution and pushing the oldest off the stack, hopefully 
c     keeping the number of vectors, m, small.


      include 'SIZE'   ! For nid/nio
      include 'TSTEP'  ! For istep

      real b(n),rvar(n,1),h1(n),h2(n),w(n),msk(n)
      integer ivar(1)
      character*6 name6
      logical ifwt,ifvec

      nn = n
      if (ifvec) nn = n*ndim

      call proj_get_ivar
     $   (m,mmx,ixb,ibb,ix,ib,ih1,ih2,ivar,n,ifvec,name6)
c     if (nid.eq.0) write(6,*) m,mmx,' IVAR2'
 
      if (m.eq.0) return

      ireset=iproj_chk_h1h2(rvar(ih1,1),rvar(ih2,1),h1,h2,n) ! Updated matrix?

      bb4 = glsc3(b,w,b,n)
      bb4 = sqrt(bb4)


c     Re-orthogonalize basis set w.r.t. new vectors if space has changed.

      if (ireset.eq.1) then

         do j=0,m-1         ! First, set B := A*X
            jb = ib+j*nn
            jx = ix+j*nn
            call proj_matvec (rvar(jb,1),rvar(jx,1),n,h1,h2,msk,name6)
         enddo

         if ((mod(istep,iocomm).eq.0).or.(istep.le.10)) then
              if (nid.eq.0) write(6,*) 'Reorthogonalize Basis:'
         endif

         call proj_ortho    ! Orthogonalize X & B basis sets
     $      (rvar(ix,1),rvar(ib,1),n,m,w,ifwt,ifvec,name6)

      endif

c     ixb is pointer to xbar,  ibb is pointer to bbar := A*xbar

      call project1_a(rvar(ixb,1),rvar(ibb,1),b,rvar(ix,1),rvar(ib,1)
     $               ,n,m,w,ifwt,ifvec)

      baf = glsc3(b,w,b,n)
      baf = sqrt(baf)
      ratio = bb4/baf

      if ((mod(istep,iocomm).eq.0).or.(istep.le.10)) then
           if (nid.eq.0) write(6,1) istep,bb4,baf,ratio,m,name6
      endif
    1 format(i8,1p3e14.5,i4,1x,a6,' PROJECT')

      return
      end
c-----------------------------------------------------------------------
      subroutine project1_a(xbar,bbar,b,xx,bb,n,m,w,ifwt,ifvec)

c     xbar is best fit in xx, bbar = A*xbar
c     b <-- b - bbar

      include 'SIZE'
      real xbar(n),bbar(n),b(n),xx(n,m),bb(n,m),w(n)
      logical ifwt,ifvec

      parameter (lxprev=20)
      real alpha(lxprev),work(lxprev)


      if (m.eq.0) return

      if (ifwt) then
         do j=1,m
            alpha(j)=vlsc3(xx(1,j),w,b,n)
         enddo
      else
         do j=1,m
            alpha(j)=vlsc2(xx(1,j),b,n)
         enddo
      endif
      call gop(alpha,work,'+  ',m)

      call cmult2(xbar,xx(1,1),alpha(1),n)
      call cmult2(bbar,bb(1,1),alpha(1),n)

      do j=1,m
         call add2s2(xbar,xx(1,j),alpha(j),n)
         call add2s2(bbar,bb(1,j),alpha(j),n)
      enddo

      call sub2(b,bbar,n)

      return
      end
c-----------------------------------------------------------------------
      function iproj_chk_h1h2(h1old,h2old,h1,h2,n)

c     Matrix has changed if h1/h2 differ from old values

      real h1(n),h2(n),h1old(n),h2old(n)

      dh1 = 0.
      dh2 = 0.
      do i=1,n
         dh1 = max(dh1,abs(h1(i)-h1old(i)))
         dh2 = max(dh2,abs(h2(i)-h2old(i)))
      enddo
      dh = max(dh1,dh2)
      dh = glmax(dh,1)  ! Max across all processors

      iproj_chk_h1h2 = 0

      if (dh.gt.0) then

         call copy(h1old,h1,n)   ! Save old h1 / h2 values
         call copy(h2old,h2,n)

         iproj_chk_h1h2 = 1      ! Force re-orthogonalization of basis

      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine proj_matvec(b,x,n,h1,h2,msk,name6)
      include 'SIZE'
      include 'TOTAL'
      real b(n),x(n),h1(n),h2(n),msk(n)
      character*6 name6

c     This is the default matvec for nekcem.

c     The code can later be updated to support different matvec
c     implementations, which would be discriminated by the character
c     string "name6"

      imsh = 1
      isd  = 1
      call axhelm3  (b,x,h1,h2,imsh,isd)       ! b = A x
      call dssum   (b,nx1,ny1,nz1)
      call col2    (b,msk,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine proj_ortho(xx,bb,n,m,w,ifwt,ifvec,name6)
      include 'SIZE'
      real xx(n,1),bb(n,1),w(n)
      character*6 name6
      logical ifwt,ifvec

      if (ifwt) then
         alpha = glsc3(xx(1,m),w,bb(1,m),n)
c        if (nid.eq.0) write(6,*) alpha,'alpha 3'
      else
         alpha = glsc2(xx(1,m),bb(1,m),n)
c        if (nid.eq.0) write(6,*) alpha,'alpha 2'
      endif
      scale = 1./sqrt(alpha)
      call cmult(xx(1,m),scale,n)
      call cmult(bb(1,m),scale,n)

      do k=m-1,1,-1  ! Reorthogonalize, starting with latest solution

         do j=m,k+1,-1
            alpha = 0.
            if (ifwt) then
               alpha = alpha + .5*(vlsc3(xx(1,j),w,bb(1,k),n)
     $                       +     vlsc3(bb(1,j),w,xx(1,k),n))
            else
               alpha = alpha + .5*(vlsc2(xx(1,j),bb(1,k),n)
     $                       +     vlsc2(bb(1,j),xx(1,k),n))
            endif
            scale = -glsum(alpha,1)
            call add2s2(xx(1,k),xx(1,j),scale,n)
            call add2s2(bb(1,k),bb(1,j),scale,n)
         enddo

         if (ifwt) then
            alpha = glsc3(xx(1,k),w,bb(1,k),n)
         else
            alpha = glsc2(xx(1,k),bb(1,k),n)
         endif
         scale = 1./sqrt(alpha)
         call cmult(xx(1,k),scale,n)
         call cmult(bb(1,k),scale,n)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine project2(x,b,n,rvar,ivar,h1,h2,msk,w,ifwt,ifvec,name6)
      real x(n),b(n),rvar(n,1),h1(n),h2(n),w(n),msk(n)
      integer ivar(1)
      character*6 name6
      logical ifwt,ifvec

      call proj_get_ivar(m,mmx,ixb,ibb,ix,ib,ih1,ih2,ivar,n,ifvec,name6)

c     ix  is pointer to X,     ib  is pointer to B
c     ixb is pointer to xbar,  ibb is pointer to bbar := A*xbar

      call project2_a(x,rvar(ixb,1),rvar(ix,1),rvar(ib,1),
     $                n,m,mmx,h1,h2,msk,w,ifwt,ifvec,name6)

      ivar(2) = m ! Update number of saved vectors

      return
      end
c-----------------------------------------------------------------------
      subroutine project2_a
     $      (x,xbar,xx,bb,n,m,mmx,h1,h2,msk,w,ifwt,ifvec,name6)

      include 'SIZE'
      real x(n),xbar(n),xx(n,1),bb(n,1),h1(n),h2(n),w(n),msk(n)
      character*6 name6
      logical ifwt,ifvec

      nn = n
      if (ifvec) nn=ndim*n

      call add2        (x,xbar,n)      ! Restore desired solution

      if (m.eq.mmx) then ! Push old vector off the stack
         do k=2,mmx
            call copy     (xx(1,k-1),xx(1,k),nn)
            call copy     (bb(1,k-1),bb(1,k),nn)
         enddo
      endif

      m = min(m+1,mmx)
c     if (nid.eq.0) write(6,*) m,mmx,' this is mmx'

      call copy        (xx(1,m),x,nn)   ! Update (X,B)
      call proj_matvec (bb(1,m),xx(1,m),n,h1,h2,msk,name6)
    
      call proj_ortho  (xx,bb,n,m,w,ifwt,ifvec,name6) ! w=mult array

      return
      end
c-----------------------------------------------------------------------
      subroutine proj_get_ivar
     $    (m,mmx,ixb,ibb,ix,ib,ih1,ih2,ivar,n,ifvec,name6)

      include 'SIZE'
      include 'TSTEP'
      include 'DRIFT'

      logical ifvec
      character*6 name6

c     integer ivar(10)
      integer ivar(1) ! FIXME misun 10/16/2014 

      integer icalld
      save    icalld
      data    icalld/0/

      if (icalld.eq.0.or.icalld.eq.istep) then
         ivar(2)=0
         icalld=istep
      endif
 
      m    = ivar(2)
      mmx  = ivar(1)

      mmx = (mxprev-4)/2

      nn = n
      if (ifvec) nn = n*ndim  ! Number of entries in a vector


      ih1  = 1
      ih2  = ih1 + n
      ixb  = ih2 + n      ! pointer to xbar
      ibb  = ixb + nn     !    "    to bbar
      ix   = ibb + nn     !    "    to B
      ib   = ix  + nn*mmx !    "    to X

c     if (nid.eq.0) write(6,*) istep,m,mmx,n,' IVAR'

      return
      end
c-----------------------------------------------------------------------
c     BEGIN DG part 
      subroutine cem_drift_dg_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      common /scrndd/ srflxd(12*2*ldim*lx1*lz1*lelt)        
      real            srflxd
      call cem_drift_dg
      call cem_drift_restrict_to_face
      call cem_drift_flux (srflxd)
      call cem_drift_add_flux_to_res (srflxd)
      call cem_drift_source_dg   
      call cem_drift_invqmass_dg
      
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      integer i
      call cem_drift_sem1
      call cem_drift_sem_source1   
      call col2(cN,dmask(1,2),npts)
      call col2(cP,dmask(1,3),npts)
      call col2(cE,dmask(1,4),npts)
      call cem_drift_invqmass_sem
      call col2(cN,dmask(1,2),npts)
      call col2(cP,dmask(1,3),npts)
      call col2(cE,dmask(1,4),npts)
c      call cem_drift_poisson1 ! return potent and electric feild cEE by constant coeff poisson equation
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_sem1               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      real    wk0(lpts4),wk1(lpts4),wk2(lpts4),wk3(lpts4)
      real    wkxyz(lpts4*3),wkx(lpts4),wky(lpts4),wkz(lpts4)
      integer  i,isd,ifld,imsh,maxit
      real    xx,yy,aa,bb,cc,tt

c...  2d case for now
      isd   = 1
      imsh  = 1

c...  pre-computation for DtN stuffs
      do ifld=1,nfield
      call rone     (dmask(1,ifld),npts)      ! for grids                         
      call dd_setdmaskp(dmask(1,ifld),npts,ifld)
      enddo

c...  added by Ying Oct 6, 2014
c     check value for nfield, if nfield=1, set dmask(:,1)=dmask(:,2)=dmask(:,3)=dmask(:,4)
      if (nfield.eq.1) then
         do i=2,LDIMT
            call copy(dmask(1,i),dmask(1,1),npts)
         enddo
      endif
      call rone    (mult,npts)
      call dssum   (mult,nx1,ny1,nz1)
      call invcol1 (mult,npts)

      do i=1,npts
         h1n(i)=-diff_n(i)/charge(i)
         h1p(i)=-diff_p(i)/charge(i)
      enddo

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call axhelm3(resN,cN,h1n,h2n,imsh,isd)!compute AU=-h1\delta u
      call axhelm3(resP,cP,h1p,h2p,imsh,isd)!compute AU=-h1\delta u

      if (ifdealias) then

      call drift_intp_u_grad_v(wk0,potent,cN)
      call add2s2 (resN,wk0,+1.0,npts)
      call drift_intp_u_grad_v(wk0,potent,cP)
      call add2s2 (resP,wk0,-1.0,npts)

      else

      call axhelm3(wk0,potent,cN,h2n,imsh,isd)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0,+1.0,npts)
      call axhelm3(wk0,potent,cP,h2p,imsh,isd)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0,npts)

      endif

      call dssum  (resN,nx1,ny1,nz1)
      call dssum  (resP,nx1,ny1,nz1)
      if (ifexct) then
      do i=1,npts
         h1e(i)=-diff_e(i)/charge(i)
      enddo
      call rzero(h2e,npts)
      call axhelm3(resE,cE,h1e,h2e,imsh,isd)!compute AU=-h1\delta u
      call dssum  (resE,nx1,ny1,nz1)
      endif

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_dg                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i
      real xx,yy,tt,cc,glamax!fixme delte later 
      
      call drift_wght_div  (resN,wNJ(1,1),wNJ(1,2),wNJ(1,3),
     $                      cNJ(1,1),cNJ(1,2),cNJ(1,3))
      call drift_wght_div  (resP,wPJ(1,1),wPJ(1,2),wPJ(1,3),
     $                      cPJ(1,1),cPJ(1,2),cPJ(1,3))
      call drift_wght_div  (resE,wEJ(1,1),wEJ(1,2),wEJ(1,3),
     $                      cEJ(1,1),cEJ(1,2),cEJ(1,3))
      
      call drift_wght_grad (resNQ(1,1),resNQ(1,2),resNQ(1,3),cN)
      call drift_wght_grad (resPQ(1,1),resPQ(1,2),resPQ(1,3),cP)
      call drift_wght_grad (resEQ(1,1),resEQ(1,2),resEQ(1,3),cE)
      call drift_wght_grad (resEN(1,1),resEN(1,2),resEN(1,3),potent)
      if (if3d) then
      do i=1,npts
         resNQ(i,1)=diff_n(i)/charge(i)*resNQ(i,1)
     $              -mu_n(i)*cN(i)*resEN(i,1)
         resNQ(i,2)=diff_n(i)/charge(i)*resNQ(i,2)
     $              -mu_n(i)*cN(i)*resEN(i,2)
         resNQ(i,3)=diff_n(i)/charge(i)*resNQ(i,3)
     $              -mu_n(i)*cN(i)*resEN(i,3)
         resPQ(i,1)=diff_p(i)/charge(i)*resPQ(i,1)
     $              +mu_p(i)*cP(i)*resEN(i,1)
         resPQ(i,2)=diff_p(i)/charge(i)*resPQ(i,2)
     $              +mu_p(i)*cP(i)*resEN(i,2)
         resPQ(i,3)=diff_p(i)/charge(i)*resPQ(i,3)
     $              +mu_p(i)*cP(i)*resEN(i,3)
         resEQ(i,1)=diff_e(i)/charge(i)*resEQ(i,1)
         resEQ(i,2)=diff_e(i)/charge(i)*resEQ(i,2)
         resEQ(i,3)=diff_e(i)/charge(i)*resEQ(i,3)
      enddo
      else
      do i=1,npts
         resNQ(i,1)=diff_n(i)/charge(i)*resNQ(i,1)
     $              -mu_n(i)*cN(i)*resEN(i,1)
         resNQ(i,2)=diff_n(i)/charge(i)*resNQ(i,2)
     $              -mu_n(i)*cN(i)*resEN(i,2)
         resPQ(i,1)=diff_p(i)/charge(i)*resPQ(i,1)
     $              +mu_p(i)*cP(i)*resEN(i,1)
         resPQ(i,2)=diff_p(i)/charge(i)*resPQ(i,2)
     $              +mu_p(i)*cP(i)*resEN(i,2)
         resEQ(i,1)=diff_e(i)/charge(i)*resEQ(i,1)
         resEQ(i,2)=diff_e(i)/charge(i)*resEQ(i,2)
      enddo
      endif

c...  FIXME
c      call cem_drift_force 

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_restrict_to_face
c-----------------------------------------------------------------------
      implicit none
c     Restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT'
      include 'POISSON'
      integer  i,j

      do j=1,ncemface
          i= cemface(j)
           fcNJ(j,1)= cNJ(i,1)
           fcNJ(j,2)= cNJ(i,2)
           fcNJ(j,3)= cNJ(i,3)
           fcPJ(j,1)= cPJ(i,1)
           fcPJ(j,2)= cPJ(i,2)
           fcPJ(j,3)= cPJ(i,3)
           fcEJ(j,1)= cEJ(i,1)
           fcEJ(j,2)= cEJ(i,2)
           fcEJ(j,3)= cEJ(i,3)
           fcN(j)   = cN(i)*diff_n(i)/charge(i)
           fcP(j)   = cP(i)*diff_p(i)/charge(i)
           fcE(j)   = cE(i)*diff_e(i)/charge(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_add_flux_to_res(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'

      real     srflx(1), a
      integer  i,j,k

      k   = nxzfl            

      if (if3d) then
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resN (i)   = resN (i)   + a*srflx( 0*k+j)
         resP (i)   = resP (i)   + a*srflx( 1*k+j)
         resE (i)   = resE (i)   + a*srflx( 2*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx( 3*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx( 4*k+j)
         resNQ(i,3) = resNQ(i,3) + a*srflx( 5*k+j)
         resPQ(i,1) = resPQ(i,1) + a*srflx( 6*k+j)
         resPQ(i,2) = resPQ(i,2) + a*srflx( 7*k+j)
         resPQ(i,3) = resPQ(i,3) + a*srflx( 8*k+j)
         resEQ(i,1) = resEQ(i,1) + a*srflx( 9*k+j)
         resEQ(i,2) = resEQ(i,2) + a*srflx(10*k+j)
         resEQ(i,3) = resEQ(i,3) + a*srflx(11*k+j)
         enddo
      else
         do j = 1,ncemface
         i = cemface(j)
         a = aream(j)
         resN (i)   = resN (i)   + a*srflx( 0*k+j)
         resP (i)   = resP (i)   + a*srflx( 1*k+j)
         resE (i)   = resE (i)   + a*srflx( 2*k+j)
         resNQ(i,1) = resNQ(i,1) + a*srflx( 3*k+j)
         resNQ(i,2) = resNQ(i,2) + a*srflx( 4*k+j)
         resPQ(i,1) = resPQ(i,1) + a*srflx( 5*k+j)
         resPQ(i,2) = resPQ(i,2) + a*srflx( 6*k+j)
         resEQ(i,1) = resEQ(i,1) + a*srflx( 7*k+j)
         resEQ(i,2) = resEQ(i,2) + a*srflx( 8*k+j)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source_dg    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'RK5'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      real    dumm1,dummy1(lpts),dummy2(lpts),dummy3(lpts)
      integer i
      real tt,xx,yy,aa,bb,cc !delete later fixme
      if (.not.ifsrc) return
     
      call usersrc (dumm1,wk1,wk2,wk3,dummy1,dummy2,dummy3)
c      call col2    (wk1,bmn,npts) !fixme for dg should including the boundary info or not.. 
      call add2s2  (resN,wk1,1.,npts) 
c      call col2    (wk2,bmn,npts) 
      call add2s2  (resP,wk2,1.,npts) 
c      call col2    (wk3,bmn,npts) 
      call add2s2  (resE,wk3,1.,npts) 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_source1  
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'
      include 'RK5'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      real    dummy1(lpts)
      real    dummy2(lpts)
      real    dummy3(lpts)
      integer i,dumm1
      if (.not.ifsrc) return

      call usersrc (dumm1,wk1,wk2,wk3,dummy1,dummy2,dummy3)
      call dssum   (wk1,nx1,ny1,nz1)
      call dssum   (wk2,nx1,ny1,nz1)
      call add2s2  (resN,wk1,1.,npts) 
      call add2s2  (resP,wk2,1.,npts)
      if(ifexct) then
      call dssum   (wk3,nx1,ny1,nz1)
      call add2s2  (resE,wk3,1.,npts)
      endif 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_invqmass_sem
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'
      call col2(resN,qbm1,npts) !resN=resN*qbm1
      call col2(resP,qbm1,npts) !resP=resP*qbm1
      if (ifexct) call col2(resE,qbm1,npts) !resE=resE*qbm1
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_invqmass_dg
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      call col2(resN      ,qbm1,npts)       !resN =resN*qbm1
      call col2(resNQ(1,1),qbm1,npts) !resNQ=resNQ*qbm1
      call col2(resNQ(1,2),qbm1,npts) !array for copy current
      call col2(resNQ(1,3),qbm1,npts) ! not used for RK45
      call col2(resP      ,qbm1,npts) !resP =resP*qbm1
      call col2(resPQ(1,1),qbm1,npts) !resPQ=resPQ*qbm1
      call col2(resPQ(1,2),qbm1,npts) 
      call col2(resPQ(1,3),qbm1,npts) 
      call col2(resE      ,qbm1,npts) !resE =resE*qbm1
      call col2(resEQ(1,1),qbm1,npts) !resEQ=resEQ*qbm1
      call col2(resEQ(1,2),qbm1,npts) 
      call col2(resEQ(1,3),qbm1,npts) 
      call copy(cNJ,resNQ,npts*3)
      call copy(cPJ,resPQ,npts*3)
      call copy(cEJ,resEQ,npts*3)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_flux (srflx)  
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'DRIFT'
      include 'RK5'   
     
      real*8  dclock,t0
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer i,j,k,e,f,i0
      real    unx0,uny0,unz0
      real    tmp,tmpi1,tmpi2,tmpi3,tmpr1,tmpr2,tmpr3
      real    tmpur,tmpui,grad_r

      integer icalld
      save    icalld
      data    icalld/0/
      character CB*3

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0
      !call usersol(RKtime,sQr,sQi,sQi,sUr,sUi,sUi) !FIXME

      k = nxzfl

      if (if3d) then
        do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
           unz0=unzm(i)
           srflx( 0*k+i)=-(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)+unz0*fcNJ(i,3))
           srflx( 1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)+unz0*fcPJ(i,3))
           srflx( 2*k+i)=-(unx0*fcEJ(i,1)+uny0*fcEJ(i,2)+unz0*fcEJ(i,3))
           srflx( 3*k+i)=-unx0*fcN(i) 
           srflx( 4*k+i)=-uny0*fcN(i) 
           srflx( 5*k+i)=-unz0*fcN(i) 
           srflx( 6*k+i)=-unx0*fcP(i) 
           srflx( 7*k+i)=-uny0*fcP(i) 
           srflx( 8*k+i)=-unz0*fcP(i) 
           srflx( 9*k+i)=-unx0*fcE(i) 
           srflx(10*k+i)=-uny0*fcE(i) 
           srflx(11*k+i)=-unz0*fcE(i) 

         enddo

        else
         do i=1,nxzfl                
           unx0=unxm(i)
           uny0=unym(i)
C.........FIXME figure out the right sign +/-
           srflx(0*k+i)=-(unx0*fcNJ(i,1)+uny0*fcNJ(i,2)) 
c           srflx(0*k+i)=srflx(0*k+i)+fcN(i)*(unx0+uny0) 
           srflx(1*k+i)=-(unx0*fcPJ(i,1)+uny0*fcPJ(i,2)) 
           srflx(2*k+i)=-(unx0*fcEJ(i,1)+uny0*fcEJ(i,2)) 
           srflx(3*k+i)=-unx0*fcN(i)  ! NxfHx
           srflx(4*k+i)=-uny0*fcN(i)  ! NxfHx
           srflx(5*k+i)=-unx0*fcP(i)  ! NxfHx
           srflx(6*k+i)=-uny0*fcP(i)  ! NxfHx
           srflx(7*k+i)=-unx0*fcE(i)  ! NxfHx
           srflx(8*k+i)=-uny0*fcE(i)  ! NxfHx

         enddo
        endif

        t0=dclock()

        if (if3d) then
           call gs_op_fields(gsh_face,srflx,nxzfl,12,1,1,0)
        else
           call gs_op_fields(gsh_face,srflx,nxzfl,9,1,1,0)
        endif
        call measure_comm(t0)
        
        if (IFPEC) call cem_drift_flux_pec(srflx)
         
        k = nxzfl

        if (if3d) then
         do i=1,nxzfl
         srflx( 0*k+i)= 0.5*srflx( 0*k+i)
         srflx( 1*k+i)= 0.5*srflx( 1*k+i) 
         srflx( 2*k+i)= 0.5*srflx( 2*k+i) 
         srflx( 3*k+i)= 0.5*srflx( 3*k+i) 
         srflx( 4*k+i)= 0.5*srflx( 4*k+i) 
         srflx( 5*k+i)= 0.5*srflx( 5*k+i) 
         srflx( 6*k+i)= 0.5*srflx( 6*k+i) 
         srflx( 7*k+i)= 0.5*srflx( 7*k+i) 
         srflx( 8*k+i)= 0.5*srflx( 8*k+i) 
         srflx( 9*k+i)= 0.5*srflx( 9*k+i) 
         srflx(10*k+i)= 0.5*srflx(10*k+i) 
         srflx(11*k+i)= 0.5*srflx(11*k+i) 
         enddo
        else

         do i=1,nxzfl
         srflx( 0*k+i)= 0.5*srflx( 0*k+i)
         srflx( 1*k+i)= 0.5*srflx( 1*k+i) 
         srflx( 2*k+i)= 0.5*srflx( 2*k+i) 
         srflx( 3*k+i)= 0.5*srflx( 3*k+i) 
         srflx( 4*k+i)= 0.5*srflx( 4*k+i) 
         srflx( 5*k+i)= 0.5*srflx( 5*k+i) 
         srflx( 6*k+i)= 0.5*srflx( 6*k+i) 
         srflx( 7*k+i)= 0.5*srflx( 7*k+i) 
         srflx( 8*k+i)= 0.5*srflx( 8*k+i) 
         enddo
        endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_flux_pec(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'DRIFT'  
      include 'PARALLEL'
      include 'NEKUSE'
      include 'RK5'    
 
      real     srflx(1)
      integer  nxyzf,nface,nxz,k,i0,i,j,e,f,ef          
      real     unx0,uny0,unz0
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      k    = nxzfl
       call usersol(RKtime,scn,scp,sce,scn,scp,sce)

      if (if3d) then
         do e= 1,nelt
         do f= 1,nfaces
         do j= 1,nxzf

         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,1)
         if (CB.eq.'PEC') then

             unx0=unxm(i)
             uny0=unym(i)
             unz0=unzm(i)
             srflx( 0*k+i)= 0.d0
             srflx( 1*k+i)= 0.d0
             srflx( 2*k+i)= 0.d0
             srflx( 3*k+i)= 2.d0*(srflx( 3*k+i)+unx0*scN(i0))
             srflx( 4*k+i)= 2.d0*(srflx( 4*k+i)+uny0*scN(i0))
             srflx( 5*k+i)= 2.d0*(srflx( 5*k+i)+unz0*scN(i0))
             srflx( 6*k+i)= 2.d0*(srflx( 6*k+i)+unx0*scP(i0))
             srflx( 7*k+i)= 2.d0*(srflx( 7*k+i)+uny0*scP(i0))
             srflx( 8*k+i)= 2.d0*(srflx( 8*k+i)+unz0*scP(i0))
             srflx( 9*k+i)= 2.d0*(srflx( 9*k+i)+unx0*scE(i0))
             srflx(10*k+i)= 2.d0*(srflx(10*k+i)+uny0*scE(i0))
             srflx(11*k+i)= 2.d0*(srflx(11*k+i)+unz0*scE(i0))

         endif
         enddo
         enddo
         enddo
      else
         do e= 1,nelt
         do f= 1,nfaces
         do j= 1,nxzf

         i = (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering
         i0= cemface(i)

         CB = CBC(f,e,1)
         if (CB.eq.'PEC') then

             unx0=unxm(i)
             uny0=unym(i)
             srflx(0*k+i)= 0.d0
c             srflx(0*k+i)= 2.d0*(srflx( 0*k+i)-(unx0+uny0)*scN(i0))
             srflx(1*k+i)= 0.d0
             srflx(2*k+i)= 0.d0
             srflx(3*k+i)= 2.d0*(srflx(3*k+i)+unx0*scN(i0))
             srflx(4*k+i)= 2.d0*(srflx(4*k+i)+uny0*scN(i0))
             srflx(5*k+i)= 2.d0*(srflx(5*k+i)+unx0*scP(i0))
             srflx(6*k+i)= 2.d0*(srflx(6*k+i)+uny0*scP(i0))
             srflx(7*k+i)= 2.d0*(srflx(7*k+i)+unx0*scE(i0))
             srflx(8*k+i)= 2.d0*(srflx(8*k+i)+uny0*scE(i0))

         endif
           
         enddo
         enddo
         enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine rk_drift_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'
      include 'DRIFT'

      integer  ii   
      real     ca,cb
c     kHX=rk4a(i)*kHX

      ca = rk4a(ii)
      cb = rk4b(ii)
      call rk4_upd(cN,kcN,resN,cb,ca,dt,npts)
      call rk4_upd(cP,kcP,resP,cb,ca,dt,npts)
      if (ifexct) call rk4_upd(cE,kcE,resE,cb,ca,dt,npts)

      return
      end
c----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs2
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
      real wkn(lpts),wkp(lpts),wke(lpts)
      real wk1(lpts),wk2(lpts),wk3(lpts)
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      tt=time+1.0
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call rone (h1n,npts)
      call rzero(h2n,npts)
      if (ifdealias) then
      call drift_intp_u_grad_v(wkn,cN,h1n)
      call drift_intp_u_grad_v(wkp,cP,h1n)
      else
      call axhelm2(wkn,cN,h1n,h2n,1,1)       ! w = A x
      call axhelm2(wkp,cP,h1n,h2n,1,1)       ! w = A x
      endif
      
      do i=1, npts
c below define rhs for potential, which should be given in source FIXME
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(pi*xx)**2*sin(tt)**2   ! u
         cc=-2.0*pi**2*sin(tt)**2*cos(2*pi*xx)
         ee=cc-charge(i)*aa
         rhs(i)=cP(i)-cN(i)
     $          +dt*(wk2(i)-wk1(i))
         rhs(i)=charge(i)*rhs(i)
         rhs(i)=rhs(i)+ee
      enddo       
      call col2(rhs,bm1,npts)
      call add2s2(rhs,wkn, 1.0*dt,npts)
      call add2s2(rhs,wkp,-1.0*dt,npts)
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_source_semi    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk1(lpts),wk2(lpts),wk3(lpts)
      integer i
      real tt,xx,yy,aa,bb,cc,a1,a2,dd !delete later fixme
c      if (.not.ifsrc) return
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call copy    (resN,cN,npts)
      call add2s2  (resN,wk1,1.*dt,npts) 
      return
      end

c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres(phi,res,dpv,h1,h2,wt,mask,isd,imsh
     $           ,n,tol,iflag)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      integer  n,outer,isd,imsh,iflag
      real     phi(n),res(n),dpv(n,3),h1(n),h2(n),wt(n),mask(n)
      real     tol,alpha,l,temp
      real*8   etime1,dnekclock

      if (nid.eq.0) write(6,*) 'start: hmh_gmres2'

      n = nx1*ny1*nz1*nelv

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
            call axhelm4(w,x,dpv,h1,h2,imsh,isd,iflag)       ! w = A x
            call dssum  (w,nx1,ny1,nz1)
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

            call axhelm4(w,v(1,j),dpv,h1,h2,imsh,isd,iflag)  ! w = A v
            call dssum  (w,nx1,ny1,nz1)
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
c            ratio = rnorm/div0
            if (nid.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,istep
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
 9999 format('\n',' ',i9,' gmres   : iteration#',i5,1p3e12.4)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_sem_bdf2_124               
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'DRIFT'
      integer  i,isd,ifld,imsh,maxit
      real tol,tt
      integer  icalled
      save     icalled
      data     icalled/0/

    
      if (icalled.eq.0) then
         do ifld=1,nfield
         call rone    (dmask(1,ifld),npts )      ! for grids                         
         call dd_setdmaskp(dmask(1,ifld),npts,ifld)
         enddo
         call rone    (mult,npts)
         call dssum   (mult,nx1,ny1,nz1)
         call invcol1 (mult,npts)
         icalled=1
      endif
c.... I just keep this temporarily
      if (istep.eq.1) then
        tt=time+dt
        call usersol(tt, scn, scp, sce, scn, scp, sce)
        call copy(cN1,cN,npts)
        call copy(cP1,cP,npts)
        call copy(potent1,potent,npts)
        call copy(cN,scN,npts)
        call copy(cP,scP,npts)
        return 
      endif     
 
c---- rhs=rhs*dt+ u^n
      call cem_drift_source_bdf2_124   
      call dssum   (resN,nx1,ny1,nz1)
      call col2    (resN,dmask(1,2),npts)  
      call dssum   (resP,nx1,ny1,nz1)
      call col2    (resP,dmask(1,3),npts)  

c...  elliptic: 2nd + constant+ 1st
      do i=1,npts
         h1n(i)=diff_n(i)/charge(i)*dt
         h2n(i)=1.5 ! mu_n*\delta\phi*dt+1.0
      enddo
      do i=1,npts
         h1p(i)=diff_p(i)/charge(i)*dt
         h2p(i)=1.5 !- mu_p*\delta\phi
      enddo

      call copy(cN1,cN,npts)
      call copy(cP1,cP,npts)
      call copy(potent1,potent,npts)

      isd   = 1
      imsh  = 1
      tol   = param(22)
      maxit=2000

      if (IFCG) then
         call cggo2 (cN,resN,h1n,h2n,mult,dmask(1,2)
     $              ,imsh,isd,npts,tol,maxit)
         call cggo2 (cP,resP,h1p,h2p,mult,dmask(1,3)
     $              ,imsh,isd,npts,tol,maxit)
      endif
      if (IFGMRES) then
         call drift_hmh_gmres2(cN,resN,h1n,h2n
     $                   ,mult,dmask(1,2),isd,imsh,npts,tol)
         call drift_hmh_gmres2(cP,resP,h1p,h2p
     $                   ,mult,dmask(1,3),isd,imsh,npts,tol)
         if (nid.eq.0) write(6,*) 'done: elliptic_gmres3'
      endif
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_drift_source_bdf2_124    
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'DRIFT'
      include 'POISSON'

      COMMON /SCRLB/ WNX   (LX1*LY1*LZ1*LELT)
     $,              WNY   (LX1*LY1*LZ1*LELT)
     $,              WNZ   (LX1*LY1*LZ1*LELT)
      real    wnx,wny,wnz
      real    wk0(lpts*3),wk1(lpts),wk2(lpts),wk3(lpts)
      integer i,npts3
      real tt,xx,yy,aa,bb,cc,a1,a2,dd !delete later fixme

      if (param(17).eq.-3) then
      call rzero(resN,npts)
      call add2s2(resN,cN,2.0,npts)
      call add2s2(resN,cN1,-0.5,npts)
      call rzero(resP,npts)
      call add2s2(resP,cP,2.0,npts)
      call add2s2(resP,cP1,-0.5,npts)
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call add2s2 (resN,wk1,1.0*dt,npts)
      call add2s2 (resP,wk2,1.0*dt,npts)
      call col2(resN,bm1,npts)
      call col2(resP,bm1,npts)

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call rzero  (wk1 ,npts)
      call add2s2 (wk1 ,cN , 2.0,npts)
      call add2s2 (wk1 ,cN1,-1.0,npts)
      call rzero  (wk2 ,npts)
      call add2s2 (wk2 ,cP , 2.0,npts)
      call add2s2 (wk2 ,cP1,-1.0,npts)
      call rzero  (wk3 ,npts)
      call add2s2 (wk3 ,potent  , 1.0,npts)
      call axhelm3(wk0 ,wk3,wk1,h2n,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0, 1.0*dt,npts)
      call axhelm3(wk0 ,wk3,wk2,h2p,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0*dt,npts)
      elseif (param(17).eq.-4) then
      call rzero(resN,npts)
      call add2s2(resN,cN,2.0,npts)
      call add2s2(resN,cN1,-0.5,npts)
      call rzero(resP,npts)
      call add2s2(resP,cP,2.0,npts)
      call add2s2(resP,cP1,-0.5,npts)
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call add2s2 (resN,wk1,1.0*dt,npts)
      call add2s2 (resP,wk2,1.0*dt,npts)
      call col2(resN,bm1,npts)
      call col2(resP,bm1,npts)

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call rzero  (wk1 ,npts)
      call add2s2 (wk1 ,cN , 2.0,npts)
      call add2s2 (wk1 ,cN1,-1.0,npts)
      call rzero  (wk2 ,npts)
      call add2s2 (wk2 ,cP , 2.0,npts)
      call add2s2 (wk2 ,cP1,-1.0,npts)
      call rzero  (wk3 ,npts)
      call add2s2 (wk3 ,potent  , 1.0,npts)
      call axhelm3(wk0 ,wk3,wk1,h2n,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0, 1.0*dt,npts)
      call axhelm3(wk0 ,wk3,wk2,h2p,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0*dt,npts)
      elseif (param(17).eq.-5) then
      call rzero(resN,npts)
      call add2s2(resN,cN,2.0,npts)
      call add2s2(resN,cN1,-0.5,npts)
      call rzero(resP,npts)
      call add2s2(resP,cP,2.0,npts)
      call add2s2(resP,cP1,-0.5,npts)
      call usersrc (57,wk1,wk2,wk3,wk1,wk2,wk3)
      call add2s2 (resN,wk1,1.0*dt,npts)
      call add2s2 (resP,wk2,1.0*dt,npts)
      call col2(resN,bm1,npts)
      call col2(resP,bm1,npts)

      call rzero(h2n,npts)
      call rzero(h2p,npts)
      call rzero  (wk1 ,npts)
      call add2s2 (wk1 ,cN , 2.0,npts)
      call add2s2 (wk1 ,cN1,-1.0,npts)
      call rzero  (wk2 ,npts)
      call add2s2 (wk2 ,cP , 2.0,npts)
      call add2s2 (wk2 ,cP1,-1.0,npts)
      call rzero  (wk3 ,npts)
      call add2s2 (wk3 ,potent  , 2.0,npts)
      call add2s2 (wk3 ,potent1 ,-1.0,npts)
      call axhelm3(wk0 ,wk3,wk1,h2n,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resN,wk0, 1.0*dt,npts)
      call axhelm3(wk0 ,wk3,wk2,h2p,1,1)!wk0=(cn\grad\Phi,\grad v)
      call add2s2 (resP,wk0,-1.0*dt,npts)
      else
      call exitt
      endif
      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_drift_poisson_bdf2_24
c----------------------------------------------------------------------
c...  implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT' 
      include 'POISSON'
      integer i
      real    tt !fixme delete finally
      integer icalled
      save    icalled
      data    icalled/0/

c     if (nid.eq.0) write(6,*) '\n call cem_drift_poissonq'
      npts =  nx1*ny1*nz1*nelt

      if (icalled.eq.0) then
          do ifld=1,nfield
          call rone    (dmask(1,ifld),npts )      ! for grids                         
          call dd_setdmaskp(dmask(1,ifld),npts,ifld)
          enddo
          call rone    (mult,npts)
          call dssum   (mult,nx1,ny1,nz1)
          call invcol1 (mult,npts)

          call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
          call rzero   (h2,npts)
          icalled = 1
      endif

      call cem_drift_poisson_rhs24  
      call col2    (rhs,bm1,npts)            ! B*rhs
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,dmask(1,1),npts)

      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = param(22) 

      if (IFCG) then
        call cggo2 (potent,rhs,h1,h2,mult,dmask(1,1)
     $             ,imsh,isd,npts,tol,maxit)
      endif
      if (IFgmres) then
         call drift_hmh_gmres2(potent,rhs,h1,h2
     $                   ,mult,dmask(1,1),isd,imsh,npts,tol)
      endif

      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent) !cEE=\grad potent 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs24
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      tt=time+dt+1.0
      if (istep.eq. 0) tt=time+1.0
      do i=1,npts
         rhs(i)=charge(i)*(cP(i)-cN(i))*2.0
         rhs(i)=rhs(i)-charge(i)*(cP1(i)-cN1(i))
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(pi*xx)**2*sin(tt)**2   ! u
         cc=-2.0*pi**2*sin(tt)**2*cos(2*pi*xx)
         ee=cc-charge(i)*aa
         rhs(i)=rhs(i)+ee
      enddo
     
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson1
c----------------------------------------------------------------------
c...  implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DRIFT' 
      include 'POISSON'
      include 'ZPER'   
      integer i
      real    tt !fixme delete finally
      integer icalled
      save    icalled
      data    icalled/0/
      real    wk0(lpts),wk1(lpts),wk2(lpts)

c     if (nid.eq.0) write(6,*) '\n call cem_drift_poissonq'
      npts =  nx1*ny1*nz1*nelt

      if (icalled.eq.0) then
         call rone    (dmask(1,1),npts*(nfield) )      ! for grids                         
         call dd_setdmaskp(dmask(1,1),npts,1)
         call dd_setdmaskp(dmask(1,2),npts,2)
         call dd_setdmaskp(dmask(1,3),npts,3)
         call dd_setdmaskp(dmask(1,4),npts,4)

          call rone    (mult,npts)
          call dssum   (mult,nx1,ny1,nz1)
          call invcol1 (mult,npts)

          call rone    (h1,npts)        ! Solve  [ h1*A + h2*B ]  phi = B*rhs
          call rzero   (h2,npts)

          if (IFGFDMDD) call gfdm_init(nx1,ny1,nz1,.false.,0.0)      

          icalled = 1
      endif

c      dt1=rktime1-rktime
      dt1=dt
      call cem_drift_poisson_rhs1_rk  
      call usersrc (57,wk0,wk1,wk2,wk0,wk1,wk2)
      call add2s2(rhs,wk0,-1.0*dt1,npts)
      call add2s2(rhs,wk1, 1.0*dt1,npts)
      call col2    (rhs,bm1,npts)            ! B*rhs
      call axhelm3(wk0,cN,h1,h2,1,1)
      call add2s2(rhs,wk0, 1.0*dt1,npts)
      call axhelm3(wk0,cP,h1,h2,1,1)
      call add2s2(rhs,wk0,-1.0*dt1,npts)
      
      call dssum   (rhs,nx1,ny1,nz1)
      call col2    (rhs,dmask(1,1),npts)
      do i=1,npts
         h1(i)=1.0+charge(i)*dt1*(mu_n(i)*cN(i)+mu_p(i)*cP(i))
      enddo
      isd   = 1
      imsh  = 1
      maxit = 2000
      tol   = param(22) 

      if (IFPRECON) then
         call gfdm_pres_solv(potent,rhs,wk1,wk2,.false.,0.0) 
      endif
      
      if (IFCG) then
        call cggo2 (potent,rhs,h1,h2,mult,dmask(1,1)
     $             ,imsh,isd,npts,tol,maxit)
      endif
      if (IFgmres) then
         call drift_hmh_gmres2(potent,rhs,h1,h2
     $                   ,mult,dmask(1,1),isd,imsh,npts,tol)
      endif
   
      call cem_drift_grad(cEE(1,1),cEE(1,2),cEE(1,3),potent) !cEE=\grad potent 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drift_poisson_rhs1_rk
c----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' 
      include 'POISSON'
      include 'EMWAVE'
      include 'RK5'
      integer i
      real xx,yy,tt,aa,bb,cc,ee
C.... define the right hand side 
C.... rhs=charge*(cP-cN)
c.... for solve -\div \permit \grad u = f 
      if (ifrk) then
         tt=rktime1+1.0
      endif
      if (ifbdf) then
         tt=time+dt+1.0
         if (istep.eq. 0) tt=time+1.0
      endif
      do i=1,npts
         rhs(i)=charge(i)*(cP(i)-cN(i))
         xx=xm1(i,1,1,1)
         yy=ym1(i,1,1,1)
         aa=sin(pi*xx)**2*sin(tt)**2   ! u
         cc=-2.0*pi**2*sin(tt)**2*cos(2*pi*xx)
         ee=cc-charge(i)*aa
         rhs(i)=rhs(i)+ee
      enddo
      return
      end
c-----------------------------------------------------------------------
! end of file Complex.F  
