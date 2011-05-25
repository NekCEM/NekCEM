!> \file cem_dg.F 
!!
!! \brief CEM main subroutines          
!!
!! A more detailed description of the file


!> \brief Brief description of the subroutine
!!
!! A more detailed descripton goes here. This is an example
!! of a detailed description of a method
!!
!! @param param1 a description of the first parameter
!! @param param2 a description of the second parameter
!! @todo things to do can be called out here
!! @return the return value goes here
!! @see cem_op()
!! @see READATW()
!! @see RDMESH()
!! @see RDMESHW()

c--------------------------------------------------------------------- 
      subroutine cem (igeom)           
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'EXPONE'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)

      if (ifeig) then
          if (param(97).ne.0) then
             call cem_eig_arpack
          else
             call cem_eig
          endif
          call exitt
      endif

      if     (IFRK4) then
          
          do ii=1,5
             call rk_c (ii)
             call cem_op
             call rk_ab(ii)
          enddo
          
C        !if (ifwake  ) call cem_divergence
          if (iffilter) call q_filter(0.05) 
C         if (ifscat  ) call total_field    
C         if (param(14).gt.0) call cem_energy        
C         if (ifwake  ) call cem_wakepotential
                                                  
      elseif (ifexp) then
          call zexp_op 
         ! call zexp_op1 
      else
          write(6,*) ' time stepping: not defined'
          call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_eig
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

#if  !defined(CEM_EIG_MEMORY_SMALL) &&    ! defined(CEM_EIG_MEMORY_BIG)
      parameter  (lelm = 1, lq1= 1  )     ! For IFRK4 runs
      parameter  (mm =lq1*lq1*6*lelm)     ! For IFEIG runs: 2D
      parameter  (mm2=mm*mm         )     ! For IFEIG runs: 2D
#endif

#ifdef CEM_EIG_MEMORY_SMALL
     !parameter  (lelm = 18,lq1=lx1 )     ! For IFEIG runs (make sure lelm > nelv for ev test)
      parameter  (lelm =  1, lq1=  1)     ! For IFEIG runs (make sure lelm > nelv for ev test)
      parameter  (mm =  1)!q1*lq1*6*lelm)     ! For IFEIG runs: 2D
      parameter  (mm2=  1)!mm*mm         )     ! For IFEIG runs: 2D
#endif

#ifdef CEM_EIG_MEMORY_BIG
     !parameter  (lelm = 3,lq1=lx1)       ! For IFEIG runs  3D (make sure lelm > nelv for ev test)
     !parameter  (lelm = 1,lq1=lx1)       ! For IFEIG runs  3D (make sure lelm > nelv for ev test)
      parameter  (lelm =  1, lq1=  1)     ! For IFEIG runs (make sure lelm > nelv for ev test)
      parameter  (mm =lq1*lq1*lq1*6*lelm) ! For IFEIG runs: 3D
      parameter  (mm2=mm*mm*mm)           ! For IFEIG runs: 3D
#endif
                                    
      common /ecem/   res( mm2)
     $            ,   ev ( mm2)
     $            ,   er ( mm )
     $            ,   ei ( mm )
     $            ,   rw1(4*mm)
      common /scrns/  srflx(6*2*ldim*lx1*lz1*lelt)

      real            res, evecres, evalresr, evalresi,rw1        

      integer         info,ii,i,jj,j
      integer         k,k1,n,n3,n6,nn,nn2   
      integer         icalld
      save            icalld
      data            icalld /0/


      n   = nx1*ny1*nz1*nelv
      n3  = nx1*ny1*nz1*nelv*3
      n6  = nx1*ny1*nz1*nelv*6

      if (n.gt.mm) then
         write(6,*) 'insufficient memory in cem_eig, abort',mm,n
         call exitt
      endif
      
      if (if3D)  then
          nn  = n6
      else
          nn  = n3
      endif

      nn2 = nn*nn 

      call rzero(res,nn2) 

      if     (IF3D) then

         k1 = 0

         do jj=1,6                  
         do j =1,n

         k1 =  k1+1
         k  = (k1-1)*nn+1  

         call rzero(hx,n) 
         call rzero(hy,n) 
         call rzero(hz,n) 
         call rzero(ex,n) 
         call rzero(ey,n) 
         call rzero(ez,n) 

         call rzero(reshx,n) 
         call rzero(reshy,n) 
         call rzero(reshz,n) 
         call rzero(resex,n) 
         call rzero(resey,n) 
         call rzero(resez,n) 

         if (jj.eq.1) hx(j,1,1,1) = 1.  ! jth col of identity matrix
         if (jj.eq.2) hy(j,1,1,1) = 1.
         if (jj.eq.3) hz(j,1,1,1) = 1.
         if (jj.eq.4) ex(j,1,1,1) = 1.                         
         if (jj.eq.5) ey(j,1,1,1) = 1.
         if (jj.eq.6) ez(j,1,1,1) = 1.

         call cem_curl
         call cem_restrict_to_face
         call cem_flux (srflx)
         call cem_add_flux_to_res (srflx)

         do i=1,n

            if (abs(reshx(i,1,1,1)).lt.1.e-13) reshx(i,1,1,1)=0.
            if (abs(reshy(i,1,1,1)).lt.1.e-13) reshy(i,1,1,1)=0.
            if (abs(reshz(i,1,1,1)).lt.1.e-13) reshz(i,1,1,1)=0.
            if (abs(resex(i,1,1,1)).lt.1.e-13) resex(i,1,1,1)=0.
            if (abs(resey(i,1,1,1)).lt.1.e-13) resey(i,1,1,1)=0.
            if (abs(resez(i,1,1,1)).lt.1.e-13) resez(i,1,1,1)=0.

         enddo

         call copy(res(k+0*n),reshx,n)
         call copy(res(k+1*n),reshy,n)
         call copy(res(k+2*n),reshz,n)
         call copy(res(k+3*n),resex,n)
         call copy(res(k+4*n),resey,n)
         call copy(res(k+5*n),resez,n)

         enddo
         enddo

      else

         if     (IFTM) then

         k1 = 0

         do jj=1,3                  

         do j =1,n

         k1 =  k1+1
         k  = (k1-1)*nn+1  

         call rzero(hx,n) 
         call rzero(hy,n) 
         call rzero(ez,n) 

         call rzero(reshx,n) 
         call rzero(reshy,n) 
         call rzero(resez,n) 

         if (jj.eq.1) hx(j,1,1,1) = 1.  ! jth col of identity matrix
         if (jj.eq.2) hy(j,1,1,1) = 1.
         if (jj.eq.3) ez(j,1,1,1) = 1.

         call cem_curl
         call cem_restrict_to_face
         call cem_flux (srflx)
         call cem_add_flux_to_res (srflx)

         do i=1,n
            if (abs(reshx(i,1,1,1)).lt.1.e-13) reshx(i,1,1,1)=0.
            if (abs(reshy(i,1,1,1)).lt.1.e-13) reshy(i,1,1,1)=0.
            if (abs(resez(i,1,1,1)).lt.1.e-13) resez(i,1,1,1)=0.
         enddo

         call copy(res(k+0*n),reshx,n)
         call copy(res(k+1*n),reshy,n)
         call copy(res(k+2*n),resez,n)

         enddo
         enddo
 
         elseif (IFTE) then

         k1 = 0

         do jj=1,3                  
         do j =1,n

         k1 =  k1+1
         k  = (k1-1)*nn+1  

         call rzero(ex,n) 
         call rzero(ey,n) 
         call rzero(hz,n) 

         call rzero(resex,n) 
         call rzero(resey,n) 
         call rzero(reshz,n) 

         if (jj.eq.1) ex(j,1,1,1) = 1.
         if (jj.eq.2) ey(j,1,1,1) = 1.  ! jth col of identity matrix
         if (jj.eq.3) hz(j,1,1,1) = 1.

         call cem_curl
         call cem_restrict_to_face
         call cem_flux (srflx)
         call cem_add_flux_to_res (srflx)

         do i=1,n
            if (abs(resex(i,1,1,1)).lt.1.e-13) resex(i,1,1,1)=0.
            if (abs(resey(i,1,1,1)).lt.1.e-13) resey(i,1,1,1)=0.
            if (abs(reshz(i,1,1,1)).lt.1.e-13) reshz(i,1,1,1)=0.
         enddo

         call copy (res(k+0*n),resex,n)
         call copy (res(k+1*n),resey,n)
         call copy (res(k+2*n),reshz,n)


         enddo
         enddo

         else
            write(6,*) 'IFTM IFTE not defined'
            call exitt
         endif

      endif 

      if (nelv.eq.1) call outmat(res,nx1*3,nx1*3,'res   ',1)

      do i=1,nn2
         write(60,*) res(i)                  
      enddo

      write(6,*) '                                      '
      write(6,*) 'IFEIG T :: computing eigenvalues      '
      write(6,*) '                                      '

      call dgeev('N','V',nn,res,nn,er,ei,ev,nn,ev,nn,rw1,4*nn,info)

      do i=1,nn
         write(59,*)  er(i),ei(i)                  
      enddo

    1 format(1p1e19.12)

      if (icalld.eq.0) then
        write(6,*) ' IF3D =',if3d,', TM=',iftm, ', CENTRAL=',
     $       ifcentral ,', PEC=',ifpec,', RK4=',IFRK4, 
     $       ', DEALIAS=',ifdealias,', SCAT=',ifscat,', IFNM =',ifnm
        write(6,*) ' IFSRC=', ifsrc,', TE=', ifte,', UPWIND =', 
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude,
     $       ', LORENTZ=',ifloren

        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=', 
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter 
        icalld= 1     
      endif
      return
      end

c-----------------------------------------------------------------------
      subroutine gather_fields(dest, myhx, myhy, myhz, myex, myey, myez)

      ! copy separate fields to contiguous storage at dest

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real myhx(*), myhy(*), myhz(*), myex(*), myey(*), myez(*)
      real dest(*)

      integer dpf

      dpf = nx1*ny1*nz1*nelv ! dof per field

      if (hxactive) call copy(dest(1+(fldidx(1)-1)*dpf), myhx,dpf)
      if (hyactive) call copy(dest(1+(fldidx(2)-1)*dpf), myhy,dpf)
      if (hzactive) call copy(dest(1+(fldidx(3)-1)*dpf), myhz,dpf)
      if (exactive) call copy(dest(1+(fldidx(4)-1)*dpf), myex,dpf)
      if (eyactive) call copy(dest(1+(fldidx(5)-1)*dpf), myey,dpf)
      if (ezactive) call copy(dest(1+(fldidx(6)-1)*dpf), myez,dpf)

      return
      end

C----------------------------------------------------------------------
      subroutine scatter_fields(myhx, myhy, myhz, myex, myey, myez, src)

      ! copy contiguous storage at src to separate fields

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real myhx(*), myhy(*), myhz(*), myex(*), myey(*), myez(*)
      real src(*)

      integer dpf

      dpf   = nx1*ny1*nz1*nelv ! dof per field

      if (hxactive) call copy(myhx, src(1+(fldidx(1)-1)*dpf), dpf)
      if (hyactive) call copy(myhy, src(1+(fldidx(2)-1)*dpf), dpf)
      if (hzactive) call copy(myhz, src(1+(fldidx(3)-1)*dpf), dpf)
      if (exactive) call copy(myex, src(1+(fldidx(4)-1)*dpf), dpf)
      if (eyactive) call copy(myey, src(1+(fldidx(5)-1)*dpf), dpf)
      if (ezactive) call copy(myez, src(1+(fldidx(6)-1)*dpf), dpf)

      return
      end
c-----------------------------------------------------------------------
      subroutine vec_energy2(a,b,c,d,n)
      real     a(1),b(1),c(1)
      real     d
      integer  n

      do i=1,n
         a(i)= a(i)+d*(b(i)**2+c(i)**2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine vec_energy3(a,b,c,d,e,n)
      real     a(1),b(1),c(1),d(1)
      real     e
      integer  n
      
      do i = 1,n
         a(i)= a(i)+e*(b(i)**2+c(i)**2+d(i)**2)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine l2_of_field(fullfld, l2norm)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'

      integer dpf
      real    l2norm

      dpf    = nx1*ny1*nz1*nelv ! dof per field
      tmp    = glsc3(fullfld, bm1, fullfld, dpf)/volvm1
      if (tmp.gt.0) tmp=sqrt(tmp)
      l2norm = tmp                                                
      return
      end
c-----------------------------------------------------------------------
      subroutine integral_of_field(fullfld, integral)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'

      integer  dpf
      real     integral
      real     fullfld(lx1,ly1,lz1,lelt)

      n        = nx1*ny1*nz1*nelv ! dof per field
      integral = glsc2(fullfld,bm1,n)/volvm1

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_eig_arpack
C----------------------------------------------------------------------
#if USE_ARPACK
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      ! calculation temporaries
      real       srflx(6*2*ldim*lx1*lz1*lelt)

      ! problem size parameters
      integer    ldofperfld, ldof
      parameter (ldofperfld   = lx1*ly1*lz1*lelv)
      parameter (ldof = ldofperfld*6)

      integer    ndofperfld, ndof

      integer    eigvcount,ritzvcount, lworkl
      parameter (eigvcount = 10)
      parameter (ritzvcount = eigvcount*2+20)
      parameter (lworkl = 3*ritzvcount*ritzvcount + 6*ritzvcount)

      ! spectral shift parameters
      real       shift
      parameter (shift = 1.)

      ! ARPACK dnaupd arguments
      character  abmat*1, awhich*2
      integer    aido, ainfo, aldv, alworkl, an, ancv, anev
      real       atol

      integer    aiparam(11), aipntr(14)
      real       aresid(ldof), av(ldof,ritzvcount)
      real       aworkd(3*ldof), aworkl(lworkl)

      ! ARPACK dneupd extra arguments
      logical   arvec
      character ahowmny*1
      logical   aselect(ritzvcount)
      real      adr(eigvcount+1)
      real      adi(eigvcount+1)
      real      az(ldof, eigvcount+1)
      integer   aldz
      real      asigmar, asigmai
      real      aworkev(3*ritzvcount)

      ! loop variable
      integer   i, loopcnt

      ! start of code

      ndofperfld = nx1*ny1*nz1*nelv

      if (if3d) then
          ndof = 6*ndofperfld
      else
          ndof = 3*ndofperfld
      endif

      write(*,*) 'entering arpack eigenvalue computation'
      write(*,*) ndofperfld, ' dofs per field'
      write(*,*) ndof, ' dofs total'

      ! set parameters for arpack
      aido    =  0    ! initial call to RCI
      abmat   = 'I'   ! standard eigenvalue problem
      an      =  ndof
      awhich  = 'LR'  ! want evalues with Largest Real part
      anev    =  eigvcount ! how many eigenvalues
      atol    =  1e-5 ! permitted relative error
      ainfo   =  0    ! use random initial residual
      ancv    =  ritzvcount ! how many ritz vectors (.ge. arpnev+2)
      aldv    =  ldof ! leading dimension of arpv

      aiparam(1) = 1  ! use exact shifts
      aiparam(3) = 100 * ritzvcount ! maximum # of iterations
      aiparam(4) = 1
      aiparam(7) = 1  ! mode 1: Ax = lambda x

      alworkl = lworkl

      loopcnt = 0
      do while (aido.ne.99)
        loopcnt = loopcnt + 1
        if (mod(loopcnt,100).eq.0) then
          write (*,*) 'arpack iteration', loopcnt
        endif
        call dnaupd
     &  (aido, abmat, an, awhich, anev, atol, aresid, 
     &  ancv, av, aldv, aiparam, aipntr,
     &  aworkd, aworkl, alworkl, ainfo)

        if (ainfo.ne.0) then
          write (*,*) 'arpack dnaupd error number:', ainfo
          write (*,*) 'converged eigenvalues     :', aiparam(5)
          call exitt
        endif

        if ((aido.ne.99).and.((aido.eq.1).or.(aido.eq.-1))) then
          call rzero(reshx, ndofperfld) 
          call rzero(reshy, ndofperfld) 
          call rzero(reshz, ndofperfld) 
          call rzero(resex, ndofperfld) 
          call rzero(resey, ndofperfld) 
          call rzero(resez, ndofperfld) 

          call scatter_fields(hx, hy, hz, ex, ey, ez, aworkd(aipntr(1)))

          call cem_curl
          call cem_restrict_to_face
          call cem_flux (srflx)
          call cem_add_flux_to_res (srflx)

          call gather_fields(aworkd(aipntr(2)), 
     &         reshx, reshy, reshz, resex, resey, resez)
          if (shift.ne.0) then
              do i=0,ndof-1
                 aworkd(aipntr(2)+i) = aworkd(aipntr(2)+i) +
     &                                 shift*aworkd(aipntr(1)+i)
              enddo
          endif
        endif
      enddo

      arvec = .false. ! no eigenvectors
      ahowmny = 'R' ! Ritz vectors, but doesn't matter (see arvec)
      aldz = ldof

      call dneupd
     &(arvec, ahowmny, aselect, adr, adi, az, aldz, asigmar, asigmai, 
     & aworkev, abmat,
     & an, awhich, anev, atol, aresid, ancv, av, aldv, aiparam, aipntr, 
     & aworkd, aworkl, alworkl, ainfo)
      if (ainfo.ne.0) then
        write (*,*) 'arpack dneupd error number:', ainfo
        call exitt
      endif

      do i=1,eigvcount
         write(59,*)  adr(i)-shift, adi(i)                  
      enddo
      write (*,*) 'eigenvalue computation complete.'
      write (*,*) 'eigenvalue with max real part:'
      write (*,*) adr(1)-shift, adi(1)
#endif // USE_ARPACK

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_op
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'

      integer n
      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)

      n= nx1*ny1*nz1*nelt
    
c... FIXME
c     if (ifsrc) then
c        call cem_inc
c        call cem_pmlsrf_to_face ! not used                 
c     endif
c... FIXME

      call cem_curl
      call cem_restrict_to_face
      call cem_flux (srflx)
      call cem_add_flux_to_res (srflx)

c... FIXME
c     if (ifpml) call pml_step
c     if (ifsrc) then 
c        if (ifsftf) then
c            call cem_source ! FIXME check out: nano pml, aps pml     
c        else
c            call cem_source ! FIXME check out: nano pml, aps pml     
c        endif  
c     endif
c
c     if (ifdrude) call cem_drude
c... FIXME
 
      call cem_invqmass

      return
      end

c-----------------------------------------------------------------------
      subroutine zexp_op1
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      integer icalled
      save    icalled
      data    icalled/0/

      common /arn7/ chkv(larnol,larnol),chkw(larnol)
      real          chkv,chkw
      common /arn8/ crV  (lx1*ly1*lz1*lelt*(3+3*(ldim-2)),larnol)
     $            , crVt (larnol,lx1*ly1*lz1*lelt*(3+3*(ldim-2)))
     $            , ceye (larnol,larnol)
     $            , czh  (larnol,larnol)
      common /arn9/ wru1 (lx1*ly1*lz1*lelt*(3+3*(ldim-2)))
     $            , wru2 (larnol)
      real          wru1, wru2         
      common /ar10/ zru  (larnol)
      complex*16    zru 
      common /ar11/ wzh1 (larnol)
     $            , wzh2 (larnol)
      complex*16    wzh1, wzh2 

      n6 = nx1*ny1*nz1*nelt*(3+3*(ndim-2))
      n  = nx1*ny1*nz1*nelt
      m  = larnol

      if (icalld.eq.0) then
          call rzero(ru ,n6)
          call rzero(rv0,n6*(m+1))
          call rzero(rv ,n6*m)
          call rzero(rvt,n6*m)
          call rzero(rh0,(m+1)*m)
          call rzero(rh ,m*m)
          call rzero(rw ,m)
          call rzero(wsp,2*m*(m+2))
          icalld = 1
      endif

      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,1)
      beta= glsc2(ru,ru,n6)
      beta= sqrt (beta)
      beta2=1.0/beta
      call copy (wru1,ru,n6)
      call cmult(wru1,beta2,n6)

      call arnoldi_vec (rv0,rh0,n6,m,ru)
      call arnoldi_mat (rvt,rv,rh,rv0,rh0,m,n6)          

      !test-start: orthonormality check
      call copy (crvt,rvt,n6*m)
      call copy (crv ,rv ,n6*m)
      call mxm  (crvt,m,crv,n6,czh,m)
      call mxm  (rvt,m,wru1,n6,wru2,1)
       
      gg=anorm(czh,m,m)
      write(6,*) 'arnorm',nx1,m,gg-1.0 
      call exitt
      call copy_r2z   (zru,wru2,m)
      !test-end  : orthonormality check

      call copy_r2z   (zh,rh,m*m)
      call zgeev  ('N','V',m,zh,m,eval,vri,m,vr,m,awork,2*m,rwork,info)
      call copy_z2z   (vri,vr,m*m)

      call zgetrf     (m,m,vri,m,ipiv,info)      !vri= LU factorization                
      call zgetri     (m,vri,m,ipiv,awork,m,info)!vri= X^(-1)

      call c_mxm      (vri,m,zru,m,wzh1,1)   !zh  = X*vri  
c     do i=1,m
c        write(6,*) 'w-1', i,real(wzh1(i)),real(vri(i,1))
c     enddo
c     do i=1,m
c        write(6,*) 'w-2', i,imag(wzh1(i)),imag(vri(i,1))
c     enddo
      do i=1,m
         write(6,*) 'wru', i,wru2(i)
      enddo
c     call exitt
      call z_exp      (eval,dt,m)           !eval= exp(dt*eval)
      call z_col2     (wzh1,eval,m)
      call c_mxm      (vr,m,wzh1,m,wzh2,1)    !zh  =X*vri  
      call copy_z2r   (rw,wzh2,m)       !zh(i,1) (i=1:m)

      call mxm        (rv,n6,rw,m,ru,1)      
      call cmult      (ru,beta,n6)      
      call flds_to_vec(ru,hx,hy,hz,ex,ey,ez,n6,n,-1)  

      return
      end
c-----------------------------------------------------------------------
      subroutine zexp_op
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      integer icalled
      save    icalled
      data    icalled/0/


      n6 = nx1*ny1*nz1*nelt*(3+3*(ndim-2))
      n  = nx1*ny1*nz1*nelt
      m  = larnol

      if (icalld.eq.0) then
          call rzero(ru ,n6)
          call rzero(rv0,n6*(m+1))
          call rzero(rv ,n6*m)
          call rzero(rvt,n6*m)
          call rzero(rh0,(m+1)*m)
          call rzero(rh ,m*m)
          call rzero(rw ,m)
          call rzero(wsp,2*m*(m+2))
          icalld = 1
      endif

      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,1)
      beta= glsc2(ru,ru,n6)
      beta= sqrt (beta)

      call arnoldi_vec (rv0,rh0,n6,m,ru)
      call arnoldi_mat (rvt,rv,rh,rv0,rh0,m,n6)          

      call copy_r2z    (zh,rh,m*m)
      call zgeev  ('N','V',m,zh,m,eval,vri,m,vr,m,awork,2*m,rwork,info)
      call copy_z2z    (vri,vr,m*m)

      call zgetrf      (m,m,vri,m,ipiv,info)      !vri= LU factorization                
      call zgetri      (m,vri,m,ipiv,awork,m,info)!vri= X^(-1)

      call z_rowscale2 (zh2,vri(1,1),eval,dt,m)     !vri =eval*X^(-1); vri=X^(-1)
      call c_mxm       (vr,m,zh2,m,zh3,1)!zh =X*vri  
      call copy_z2r    (rw,zh3,m)        !zh(i,1) (i=1:m)

      call mxm         (rv,n6,rw,m,ru,1)      
      call cmult       (ru,beta,n6)      
      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,-1)  

      return
      end
c-----------------------------------------------------------------------
      subroutine exp_op 
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'

      integer icalled
      save    icalled
      data    icalled/0/

      n6 = nx1*ny1*nz1*nelt*(3+3*(ndim-2))
      n  = nx1*ny1*nz1*nelt
      m  = larnol
      
      if (icalld.eq.0) then
          call rzero(ru, n6)
          call rzero(rv0,n6*(m+1))
          call rzero(rh0,(m+1)*m )
          call rzero(rv ,n6*m)
          call rzero(rvt,n6*m)
          call rzero(rh ,m*m)
          call rzero(rw ,m)
          call rzero(wsp,2*m*(m+2))
          icalld = 1
      endif

      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,1)  
      call arnoldi_vec (rv0,rh0,n6,m,ru)   
      call arnoldi_mat (rvt,rv,rh,rv0,rh0,m,n6)  
      call mxm         (rvt,m,ru,n6,rw,1) 
      call dnchbv      (m,dt,rh,m,rw,wsp)   
      call mxm         (rv,n6,rw,m,ru,1)  
      call flds_to_vec (ru,hx,hy,hz,ex,ey,ez,n6,n,-1)        

      return
      end
c-----------------------------------------------------------------------
      subroutine dtcontrol(stepopt,step,u,u1,u2,n,m,tol)

      integer i,n,m
      real    error,const,tol
      real    stepopt,newstep
      real    step,errchk   
      real    u(n),u1(n),u2(n)
   
      call sub3 (u,u1,u2,n) ! difference u=u1-u2
      error = glamax(u,n)   ! get maximun error

      const0= (step**m*(1.0-1.0/2.0**m))
      const =  error/const0                             

      newstep= step       
      errchk = const*newstep**m 

      if    (error.lt.tol) then
          do i = 1,10000000
             errchk = const*newstep**i 
             if (errchk .lt. tol ) then
                 newstep = newstep*2.0         
                 errchk = const*newstep**i 
                 if (errchk .gt. tol ) then
                 newstep=newstep                 
                 write(6,*) '(dt/2)^',i-1,', DT=',newstep,', tol=',tol 
                 write(6,*) 'STEP SIZE INCREASED !!!',errchk           
                 stepopt = newstep
                 return
                 endif
             endif
          enddo
      elseif (error.gt.tol ) then
          do i = 1,10000000
             errchk = const*newstep**i 
             if (errchk .gt. tol ) then
                 newstep = newstep/2.0         
                 errchk = const*newstep**i 
                if (errchk .lt. tol ) then
                 newstep=newstep                     
                 write(6,*) '(dt*2)^',i+1,', DT=',newstep
                 write(6,*) 'STEP SIZE REDUCED !!!',errchk,'DT=',newstep
                 stepopt = newstep
                 return
                endif
             endif
          enddo
      else
          write(6,*) 'ERROR:: OUT OF RANGE'  ! FIXME                          
          call exitt
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine flds_to_vec(u,ux1,uy1,uz1,ux2,uy2,uz2,n6,n,dir)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  n6,n,dir
      real     u(n6)
      real     ux1(n),uy1(n),uz1(n)
      real     ux2(n),uy2(n),uz2(n)

      if   (dir.eq.1) then

        if        (IF3D) then
           call copy (u(1 + 0*n),ux1,n)
           call copy (u(1 + 1*n),uy1,n)
           call copy (u(1 + 2*n),uz1,n)
           call copy (u(1 + 3*n),ux2,n)
           call copy (u(1 + 4*n),uy2,n)
           call copy (u(1 + 5*n),uz2,n)
        else
         if     (IFTM) then
           call copy (u(1 + 0*n),ux1,n)
           call copy (u(1 + 1*n),uy1,n)
           call copy (u(1 + 2*n),uz2,n)
         elseif (IFTE) then
           call copy (u(1 + 0*n),ux2,n)
           call copy (u(1 + 1*n),uy2,n)
           call copy (u(1 + 2*n),uz1,n)
         endif
        endif

      elseif (dir.eq.-1) then 

        if       (IF3D) then
           call copy (ux1,u(1 + 0*n),n)
           call copy (uy1,u(1 + 1*n),n)
           call copy (uz1,u(1 + 2*n),n)
           call copy (ux2,u(1 + 3*n),n)
           call copy (uy2,u(1 + 4*n),n)
           call copy (uz2,u(1 + 5*n),n)
       else
        if     (IFTM) then
           call copy (ux1,u(1 + 0*n),n)
           call copy (uy1,u(1 + 1*n),n)
           call copy (uz2,u(1 + 2*n),n)
        elseif (IFTE) then
           call copy (ux2,u(1 + 0*n),n)
           call copy (uy2,u(1 + 1*n),n)
           call copy (uz1,u(1 + 2*n),n)
        endif
       endif

      else
         write(6,*) 'ERROR:: define forward=1, backward=-1 '
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine arnoldi_mat (vt,v,h,v0,h0,m,n)
      integer  m,n, i,j,k
      real     h0(m+1,m),h(m,m)
      real     v0(n,m+1)
      real     v (n,m)
      real     vt(m,n)     

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
      subroutine arnoldi_vec(v,h,n,m,b)
      integer  e,ii,jj, i,j,k
      real     h(m+1,m)         
      real     v(n,m+1),b(n)
      real     cnt, eps, s1, si

      eps = 1.e-13

      s1 = glsc2(b,b,n)  
      si = 1./sqrt(s1) ! should add a check and a comment here

      call cmult2(v(1,1),b,si,n)

      do j=1,m 

         call amult (v(1,j+1),v(1,j),n)                              

         do i=1,j
            h(i,j) = glsc2 (v(1,j+1),v(1,i),n)
            cnt    =       -h(i,j  )
            call add2s2    (v(1,j+1),v(1,i),cnt,n)
         enddo

         s1 = glsc2(v(1,j+1),v(1,j+1),n)
         s1 = sqrt(s1)                                
         h(j+1,j) = s1                                      

         if (abs(h(j+1,j)).le. eps) then
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
      subroutine arnoldi_vec1(v,h,n,m,b)
      integer  e,ii,jj, i,j,k
      real     h(m+1,m)
      real     v(n,m+1),b(n)
      real     cnt, eps, s1, si

      eps = 1.e-13

      do j=1,m

         call amult (v(1,j+1),v(1,j),n)

         do i=1,j
            h(i,j) = glsc2 (v(1,j+1),v(1,i),n)
            cnt    =       -h(i,j  )
            call add2s2    (v(1,j+1),v(1,i),cnt,n)
         enddo

         s1 = glsc2(v(1,j+1),v(1,j+1),n)
         s1 = sqrt(s1)
         h(j+1,j) = s1

         if (abs(h(j+1,j)).le. eps) then
             write(6,*) 'h', h(j+1,j)
             return   !  should add a comment here
         else
             si = 1.0/s1
             call cmult(v(1,j+1),si,n)
         endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine amult(w,v,n6)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      common  /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real     w(n6),v(n6)
      integer  n                 

      n = nx1*ny1*nz1*nelt 

      if        (IF3D) then
            call copy(hx,v(1 + 0*n),n)
            call copy(hy,v(1 + 1*n),n)
            call copy(hz,v(1 + 2*n),n)
            call copy(ex,v(1 + 3*n),n)
            call copy(ey,v(1 + 4*n),n)
            call copy(ez,v(1 + 5*n),n)
      else
         if     (IFTM) then
            call copy(hx,v(1 + 0*n),n)
            call copy(hy,v(1 + 1*n),n)
            call copy(ez,v(1 + 2*n),n)
         elseif (IFTE) then
            call copy(ex,v(1 + 0*n),n)
            call copy(ey,v(1 + 1*n),n)
            call copy(hz,v(1 + 2*n),n)
         else
            write(6,*) 'ERROR::IFTM or IFTE not defined'
            call exitt
         endif
      endif

      call cem_curl
      call cem_restrict_to_face
      call cem_flux (srflx)
      call cem_add_flux_to_res (srflx)
      call cem_invqmass

      if        (IF3D) then
            call copy(w(1 + 0*n),reshx,n)
            call copy(w(1 + 1*n),reshy,n)
            call copy(w(1 + 2*n),reshz,n)
            call copy(w(1 + 3*n),resex,n)
            call copy(w(1 + 4*n),resey,n)
            call copy(w(1 + 5*n),resez,n)
      else
         if     (IFTM) then
            call copy(w(1 + 0*n),reshx,n)
            call copy(w(1 + 1*n),reshy,n)
            call copy(w(1 + 2*n),resez,n)
         elseif (IFTE) then
            call copy(w(1 + 0*n),resex,n)
            call copy(w(1 + 1*n),resey,n)
            call copy(w(1 + 2*n),reshz,n)
         else
            write(6,*) 'ERROR::IFTM or IFTE not defined'
            call exitt
         endif
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine copy_r2z(z,r,n)
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         z(i)=r(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2r(r,z,n)
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         r(i)=real(z(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2z(z2,z1,n)
      complex*16 z1(n),z2(n)

      do i=1,n
         z2(i)=z1(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_col2(a,b,n)
      complex*16 a(n),b(n)

      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c_mxm(a,n1,b,n2,c,n3) ! complex mxm
      include 'SIZE'
      include 'PARALLEL'   
c
      real a(1),b(1),c(1)   
      real c_one(2),c_zero(2)

      save c_one   ,c_zero
      data c_one  / 1. , 0. /
      data c_zero / 0. , 0. /

      !c=n1xn3  b=n2xn3  a=n1xn2

       call zgemm( 'N','N',n1,n3,n2,c_one,a,n1,b,n2,c_zero,c,n1) !(ifblas)
      !call cgemm( 'N','N',n1,n3,n2,c_one,a,n1,b,n2,c_zero,c,n1)

      return
      end
c-----------------------------------------------------------------------
      real function anorm(a,n,m)
      real a(n,m)

      anorm =0.0
      do j=1,m
         tmp=0.0
      do i=1,n
         tmp=tmp+abs(a(i,j))
      enddo      
      if (tmp.ge.anorm) anorm=tmp
      enddo      

      return
      end
c-----------------------------------------------------------------------
      subroutine z_rowscale2(a,b,c,dt,n)
      complex*16 a(n),b(n),c(n)
      real       dt

      do i=1,n
         a(i)=b(i)*exp(c(i)*dt)
      enddo
      return
      end 
c-----------------------------------------------------------------------
      subroutine z_rowscale(a,s,n)
      complex*16 a(n,n),s(n)

      do i=1,n
      do j=1,n
         a(i,j) = s(i)*a(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_exp(f,dt,n)
      real*8     dt
      complex*16 f(n),a

      do i=1,n
         a    = dt*f(i)
         f(i) = exp(a)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine z_exp_arnoldi(eh,h,dt,m,ch,eval,evec,eveci,ipiv,w1,w2)

      include 'SIZE'
      include 'PARALLEL'

      real     eh(m,m),h(m,m)       
      real     ch(2,m,m),eval(2*m),evec(2*m*m),eveci(2*m*m)
      real     w1(4*m*m),w2(4*m*m)                  
      real     dt
      integer  ipiv(m) 
      integer  info,i,j
      
      if (wdsize.eq.4) then
          write(6,*) 'real4 not supported in z_exp_arnoldi'
          call exitt
      endif

      do j=1,m
      do i=1,m
         ch(1,i,j) = h(i,j)
         ch(2,i,j) = 0.
      enddo
      enddo

      call zgeev('N','V',m,ch,m,eval,evec,m,evec,m,w1,2*m,w2,info)

      if (info.ne.0) then
         write(6,*) 'exp_arnoldi fail:',info,m
      endif

      call copy       (eveci,evec,2*m*m)            
      call zgetri     (m,eveci,m,ipiv,w2,2*m,info) ! inverse of evec = X^(-1)
      call z_exp      (eval,dt,m )                 ! dt*lambda
      call z_rowscale (eveci,eval,m )              ! eveci=exp(dt*eval)*X^(-1)
      call c_mxm      (evec,m,eveci,m,ch,m)        ! ch=X*exp(dt*lambda)*X^(-1)  

      do j=1,m
      do i=1,m
         eh(i,j) = ch(1,i,j)
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine sizetest(ch,n) ! complex mxm
      complex*16 ch(n,n)
      integer i,j

      write(6,*) 'ttt'

      do j=1,n
      do i=1,n
         write(6,*) 'complex',i,j,ch(i,j)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_uvp2
C----------------------------------------------------------------------

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_uvp           
C----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'EMWAVE'

      integer e

      do e=1,nelt
         call nekuvp(e)
      enddo

      call cem_mater1 
      call cem_mater2
      call cem_mater3

      return
      end
C----------------------------------------------------------------------
      subroutine cem_energy         
C----------------------------------------------------------------------
c...  this works OK with param(103)=3: needs to be updated
c...  by completely replacing old subroutine cem_out_fields
c...  09/11/2006 misun

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'POISSON'
      
      logical    outex, outey, outez, outhx, outhy, outhz
      integer    n

      integer dumpno
      integer icalld
      save    icalld,dumpno
      data    icalld /0/, dumpno /0/

      common /timeavg/ lst_period_1,lst_period_2
      integer          lst_period_1,lst_period_2

      n = nx1*ny1*nz1*nelt

c...  will be meaning ful only for pure-scattered

      tmp= dt*cspeed/(param(69)*nmscale)
      period_length= (param(69)*nmscale)/cspeed

      if (iostep.lt.0)  iostep= nopstep  ! auto assign if negative
      
      if (istep.le.nopstep) then
          nperiod= 1
          if (mod(istep,iostep).eq.1) then
          if (nid.eq.0) write(6,*) '.'
          if (nid.eq.0) write(6,*) 'First period = ',nperiod,nopstep
          endif                                 
      else
          if (mod(istep,nopstep).eq.1) then
              nperiod= nperiod+1
          if (nperiod.le.n_period) then
          if (nid.eq.0) write(6,*) 'Start new period = ',nperiod    
          endif
          endif
      endif

      if (mod(istep,nopstep).eq.1) then
          call rzero(avgtSx,n)
          call rzero(avgtSy,n)
          call rzero(avgtSz,n)
          call rzero(avgSx,n)
          call rzero(avgSy,n)
          call rzero(avgSz,n)
          call rzero(avgISx,n)
          call rzero(avgISy,n)
          call rzero(avgISz,n)
          call rzero(engyE,n)
          call rzero(engyH,n)
          call rzero(engyIE,n)
          call rzero(engyIH,n)
      endif

      if      (ifscat) then
          call vcross (poytSx,poytSy,poytSz,
     $                 totex,totey,totez,tothx,tothy,tothz,n)
          call vcross (poySx,poySy,poySz,
     $                 ex,ey,ez,hx,hy,hz,n)
          call vcross (poyISx,poyISy,poyISz,
     $                 incex,incey,incez,inchx,inchy,inchz,n)

      elseif (ifsftf) then
          call vcross (poySx,poySy,poySz,ex,ey,ez,hx,hy,hz,n)
          call vcross (poyISx,poyISy,poyISz,
     $                 incex,incey,incez,inchx,inchy,inchz,n)

      endif
 
      call add2s2(avgtSx,poytSx,tmp,n)
      call add2s2(avgtSy,poytSy,tmp,n)
      call add2s2(avgtSx,poytSz,tmp,n)
      call add2s2(avgSx ,poySx ,tmp,n)
      call add2s2(avgSy ,poySy ,tmp,n)
      call add2s2(avgSz ,poySz ,tmp,n)
      call add2s2(avgISx,poyISx,tmp,n)
      call add2s2(avgISy,poyISy,tmp,n)
      call add2s2(avgISz,poyISz,tmp,n)
                    
      call vec_energy3(engyE,ex,ey,ez,tmp,n)
      call vec_energy3(engyH,hx,hy,hz,tmp,n)
      call vec_energy3(engyIE,incex,incey,incez,tmp,n)
      call vec_energy3(engyIH,inchx,inchy,inchz,tmp,n)

      return
      end

c ...mpi file routines starts here...
C----------------------------------------------------------------------
      subroutine cem_out_fields6
C----------------------------------------------------------------------
c...  this works OK with param(103)=6,-6,5,8:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-07-08 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'POISSON'

      logical outex, outey, outez, outhx, outhy, outhz
      integer n
      integer numgroups,numfields,maxnumnodes
      integer maxnelt      
      save    maxnelt 

      integer icalld,dumpno
      save    icalld,dumpno
      data    icalld /0/,dumpno /0/

      n = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          maxnelt = iglmax(nelt,1)
          icalld  = 1
      endif

c ... VTK FORMAT BINARY

      numgroups = param(82)             ! param(82) = # groups 
      numfields = param(83)             ! param(83) = max # fields gonna write
      maxnumnodes = nx1*ny1*nz1*maxnelt ! max # nodes per processor

      call initrbio (numgroups, numfields, maxnumnodes )        !rbio_nekcem.c
      call openfile6 (dumpno,nid)     ! rbio_nekcem.c
      call vtk_dump_header6()        ! write nodes and cells

      outex = exactive
      outey = eyactive
      outez = ezactive
      outhx = hxactive
      outhy = hyactive
      outhz = hzactive

      if (ifsol) then

        if (outhx) call sub3(errhx, hx, shx, n)
        if (outhy) call sub3(errhy, hy, shy, n)
        if (outhz) call sub3(errhz, hz, shz, n)
        if (outex) call sub3(errex, ex, sex, n)
        if (outey) call sub3(errey, ey, sey, n)
        if (outez) call sub3(errez, ez, sez, n)

        call vtk_dump_field6(1,hx,hy,hz,istep,time) !11 errH
        call vtk_dump_field6(2,ex,ey,ez,istep,time) !12 errE
        call vtk_dump_field6(11,errhx,errhy,errhz,istep,time) !11 errH
        call vtk_dump_field6(12,errex,errey,errez,istep,time) !12 errE
        call vtk_dump_field6(13,shx,shy,shz,istep,time)       !13 solhx
        call vtk_dump_field6(14,sex,sey,sez,istep,time)       !14 solex
        call vtk_dump_field6(101,phi,phi,phi,istep,time)      !14 potE 
        call vtk_dump_field6(102,sphi,sphi,sphi,istep,time)   !14 spotE  
      else

       write(6,*) 'fd 6----'
       call vtk_dump_field6(21,tothx,tothy,tothz,istep,time)    !21 totH
       call vtk_dump_field6(22,totex,totey,totez,istep,time)    !22 totE
       call vtk_dump_field6(23,hx,hy,hz,istep,time)             !23 scatH
       call vtk_dump_field6(24,ex,ey,ez,istep,time)             !24 scatE
       call vtk_dump_field6(25,inchx,inchy,inchz,istep,time)    !25 incH
       call vtk_dump_field6(26,incex,incey,incez,istep,time)    !26 incE
       write(6,*) 'fd 6----'

       if (ifdrude) then
       call vtk_dump_field6(27,engyE,engyE,engyE,istep,time)    !27 engyH
       call vtk_dump_field6(28,engyIE,engyIE,engyIE,istep,time) !27 engyE
       call vtk_dump_field6(29,poySX,poySY,poySZ,istep,time)    !27 engyH
       call vtk_dump_field6(30,poytSX,poytSY,poytSZ,istep,time) !27 engyH
      !call vtk_dump_field6(31,poyISX,poyISY,poyISZ,istep,time) !27 !engyH
       call vtk_dump_field6(32,avgSX,avgSY,avgSZ,istep,time)    !27 engyH
       call vtk_dump_field6(33,avgtSX,avgtSY,avgtSZ,istep,time) !27 engyH
      !call vtk_dump_field6(34,avgISX,avgISY,avgISZ,istep,time) !27 !engyH
       endif

       call vtk_dump_field6(201,faceif,faceif,faceif,istep,time)!51 region
       call vtk_dump_field6(202,faceid,faceid,faceid,istep,time)!52 face     

       call vtk_dump_field6(203,faceb0,faceb0,faceb0,istep,time)!53 no bdry       
       call vtk_dump_field6(204,faceb1,faceb1,faceb1,istep,time)!54 pec         
       call vtk_dump_field6(205,pmlsigma(1,1,1,1,1),pmlsigma(1,1,1,1,2),
     $                          pmlsigma(1,1,1,1,3),istep,time) !55 pml       
       call vtk_dump_field6(206,faceb3,faceb3,faceb3,istep,time)!56 p              

      endif


      call closefile6()  ! vtkbin.c

c      call pvtk_xml(dumpno)
c     if(param(103).eq.8) then
c       if(nid.eq.0) call pvtk_nmm(dumpno)
c     elseif(param(103).eq.5) then
c       if(nid.eq.0) call pvtk_nm(dumpno)
c     endif

      dumpno = dumpno +1

      if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
        write(6,*)
        write(6,*) ' IFSRC=',ifsrc, ', CENTRAL=',ifcentral ,
     $    ', PEC=',ifpec,', RK4=',IFRK4,', DEALIAS=',ifdealias,
     $    ', SCAT=',ifscat,', SFTF =',ifsftf,', IFNM =',ifnm
        write(6,*) ' TM=', iftm,', TE=', ifte,', UPWIND =',
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude
        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=',
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter

        if (IFDEALIAS) write(6,33) 'Dealias_M  =',lxd
        if (IFEXP    ) write(6,33) 'Krylov_Dim =',marnol
        icalld =1
        write(6,*) ' Direct dumping into VTK p(103) = 6 or -6:::::::'

      endif

  33  format(a14,i2)

      return
      end

C----------------------------------------------------------------------
      subroutine cem_out_fields4
C----------------------------------------------------------------------
c...  this works OK with param(103)=4:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-06-28 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'POISSON'

      logical    outex, outey, outez, outhx, outhy, outhz
      integer    n

      integer dumpno
      integer icalld
      save    icalld,dumpno
      data    icalld /0/, dumpno /0/

      common /timeavg/ lst_period_1,lst_period_2
      integer          lst_period_1,lst_period_2

      n = nx1*ny1*nz1*nelt

c ... VTK FORMAT BINARY

      dumpno= dumpno + irstdumpno   ! added restart dumpno misun 5/25/11

      call openfile4 (dumpno,nid)   ! vtkbin.c

c...  want to pass (istep,dt,time,dumpno) for header
      call vtk_dump_header4()       ! write nodes and cells 

      outex = exactive
      outey = eyactive
      outez = ezactive
      outhx = hxactive
      outhy = hyactive
      outhz = hzactive
      if (ifsol) then

        if (outhx) call sub3(errhx, hx, shx, n)
        if (outhy) call sub3(errhy, hy, shy, n)
        if (outhz) call sub3(errhz, hz, shz, n)
        if (outex) call sub3(errex, ex, sex, n)
        if (outey) call sub3(errey, ey, sey, n)
        if (outez) call sub3(errez, ez, sez, n)

        call vtk_dump_field4(1,hx,hy,hz,istep,time) !11 errH
        call vtk_dump_field4(2,ex,ey,ez,istep,time) !12 errE
        call vtk_dump_field4(11,errhx,errhy,errhz,istep,time) !11 errH
        call vtk_dump_field4(12,errex,errey,errez,istep,time) !12 errE
        call vtk_dump_field4(13,shx,shy,shz,istep,time)       !13 solhx
        call vtk_dump_field4(14,sex,sey,sez,istep,time)       !14 solex
        call vtk_dump_field4(101,phi,phi,phi,istep,time)      !14 potE 
        call vtk_dump_field4(102,sphi,sphi,sphi,istep,time)   !14 spotE  
      else

       call vtk_dump_field4(21,tothx,tothy,tothz,istep,time)    !21 totH
       call vtk_dump_field4(22,totex,totey,totez,istep,time)    !22 totE
       call vtk_dump_field4(23,hx,hy,hz,istep,time)             !23 scatH
       call vtk_dump_field4(24,ex,ey,ez,istep,time)             !24 scatE
       call vtk_dump_field4(25,inchx,inchy,inchz,istep,time)    !25 incH
       call vtk_dump_field4(26,incex,incey,incez,istep,time)    !26 incE

       if (ifdrude) then
       call vtk_dump_field4(27,engyE,engyE,engyE,istep,time)    !27 engyH
       call vtk_dump_field4(28,engyIE,engyIE,engyIE,istep,time) !27 engyE
       call vtk_dump_field4(29,poySX,poySY,poySZ,istep,time)    !27 engyH
       call vtk_dump_field4(30,poytSX,poytSY,poytSZ,istep,time) !27 engyH
      !call vtk_dump_field4(31,poyISX,poyISY,poyISZ,istep,time) !27 !engyH
       call vtk_dump_field4(32,avgSX,avgSY,avgSZ,istep,time)    !27 engyH
       call vtk_dump_field4(33,avgtSX,avgtSY,avgtSZ,istep,time) !27 engyH
      !call vtk_dump_field4(34,avgISX,avgISY,avgISZ,istep,time) !27 !engyH
       endif

       call vtk_dump_field4(201,faceif,faceif,faceif,istep,time)!51 region
       call vtk_dump_field4(202,faceid,faceid,faceid,istep,time)!52 face     

       call vtk_dump_field4(203,faceb0,faceb0,faceb0,istep,time)!53 no bdry       
       call vtk_dump_field4(204,faceb1,faceb1,faceb1,istep,time)!54 pec         
       call vtk_dump_field4(205,pmlsigma(1,1,1,1,1),pmlsigma(1,1,1,1,2),
     $                          pmlsigma(1,1,1,1,3),istep,time) !55 pml       
       call vtk_dump_field4(206,faceb3,faceb3,faceb3,istep,time)!56 p              

      endif


      call closefile4()  ! vtkbin.c

c      call pvtk_xml(dumpno)

      dumpno = dumpno +1

      if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
        write(6,*)
        write(6,*) ' IFSRC=',ifsrc, ', CENTRAL=',ifcentral ,
     $    ', PEC=',ifpec,', RK4=',IFRK4,', DEALIAS=',ifdealias,
     $    ', SCAT=',ifscat,', SFTF =',ifsftf,', IFNM =',ifnm
        write(6,*) ' TM=', iftm,', TE=', ifte,', UPWIND =',
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude
        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=',
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter

        if (IFDEALIAS) write(6,33) 'Dealias_M  =',lxd
        if (IFEXP    ) write(6,33) 'Krylov_Dim =',marnol
        icalld =1
        write(6,*) ' Direct dumping into VTK Binary p(103) = 4 :::::::'

      endif

  33  format(a14,i2)

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields3
C----------------------------------------------------------------------
c...  this works OK with param(103)=3: needs to be updated
c...  by completely replacing old subroutine cem_out_fields
c...  09/11/2006 misun

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'POISSON'
      
      logical    outex, outey, outez, outhx, outhy, outhz
      integer    n

      integer dumpno
      integer icalld
      save    icalld,dumpno
      data    icalld /0/, dumpno /0/

      common /timeavg/ lst_period_1,lst_period_2
      integer          lst_period_1,lst_period_2

      n = nx1*ny1*nz1*nelt

c ... VTK FORMAT BINARY
      
      call openfile (dumpno,nid)     ! vtkbin.c
      call vtk_dump_header3()        ! write nodes and cells

      outex = exactive
      outey = eyactive
      outez = ezactive
      outhx = hxactive
      outhy = hyactive
      outhz = hzactive

      if (ifsol) then

        if (outhx) call sub3(errhx, hx, shx, n)
        if (outhy) call sub3(errhy, hy, shy, n)
        if (outhz) call sub3(errhz, hz, shz, n)
        if (outex) call sub3(errex, ex, sex, n)
        if (outey) call sub3(errey, ey, sey, n)
        if (outez) call sub3(errez, ez, sez, n)

        call vtk_dump_field3(1,hx,hy,hz,istep,time) !11 errH
        call vtk_dump_field3(2,ex,ey,ez,istep,time) !12 errE
        call vtk_dump_field3(11,errhx,errhy,errhz,istep,time) !11 errH
        call vtk_dump_field3(12,errex,errey,errez,istep,time) !12 errE
        call vtk_dump_field3(13,shx,shy,shz,istep,time)       !13 solhx
        call vtk_dump_field3(14,sex,sey,sez,istep,time)       !14 solex
        call vtk_dump_field3(101,phi,phi,phi,istep,time)      !14 potE 
        call vtk_dump_field3(102,sphi,sphi,sphi,istep,time)   !14 spotE  
      else

       call vtk_dump_field3(21,tothx,tothy,tothz,istep,time)    !21 totH
       call vtk_dump_field3(22,totex,totey,totez,istep,time)    !22 totE
       call vtk_dump_field3(23,hx,hy,hz,istep,time)             !23 scatH
       call vtk_dump_field3(24,ex,ey,ez,istep,time)             !24 scatE
       call vtk_dump_field3(25,inchx,inchy,inchz,istep,time)    !25 incH
       call vtk_dump_field3(26,incex,incey,incez,istep,time)    !26 incE

       if (ifdrude) then
       call vtk_dump_field3(27,engyE,engyE,engyE,istep,time)    !27 engyH
       call vtk_dump_field3(28,engyIE,engyIE,engyIE,istep,time) !27 engyE
       call vtk_dump_field3(29,poySX,poySY,poySZ,istep,time)    !27 engyH
       call vtk_dump_field3(30,poytSX,poytSY,poytSZ,istep,time) !27 engyH
      !call vtk_dump_field3(31,poyISX,poyISY,poyISZ,istep,time) !27 engyH
       call vtk_dump_field3(32,avgSX,avgSY,avgSZ,istep,time)    !27 engyH
       call vtk_dump_field3(33,avgtSX,avgtSY,avgtSZ,istep,time) !27 engyH
      !call vtk_dump_field3(34,avgISX,avgISY,avgISZ,istep,time) !27 engyH
       endif

       call vtk_dump_field3(201,faceif,faceif,faceif,istep,time)!51 region
       call vtk_dump_field3(202,faceid,faceid,faceid,istep,time)!52 face     

       call vtk_dump_field3(203,faceb0,faceb0,faceb0,istep,time)!53 no bdry       
       call vtk_dump_field3(204,faceb1,faceb1,faceb1,istep,time)!54 pec         
       call vtk_dump_field3(205,pmlsigma(1,1,1,1,1),pmlsigma(1,1,1,1,2),
     $                          pmlsigma(1,1,1,1,3),istep,time) !55 pml       
       call vtk_dump_field3(206,faceb3,faceb3,faceb3,istep,time)!56 p              

      endif


      call closefile()  ! vtkbin.c

      if(nid.eq.0) call pvtk_xml3(dumpno)

      dumpno = dumpno +1

      if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
        write(6,*)
        write(6,*) ' IFSRC=',ifsrc, ', CENTRAL=',ifcentral ,
     $    ', PEC=',ifpec,', RK4=',IFRK4,', DEALIAS=',ifdealias,
     $    ', SCAT=',ifscat,', SFTF =',ifsftf,', IFNM =',ifnm 
        write(6,*) ' TM=', iftm,', TE=', ifte,', UPWIND =', 
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude,
     $       ', LORENTZ=',ifloren

        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=', 
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter 

        if (IFDEALIAS) write(6,33) 'Dealias_M  =',lxd
        if (IFEXP    ) write(6,33) 'Krylov_Dim =',marnol
        icalld =1 
        write(6,*) ' Direct dumping into VTK Binary format :::::::'

      endif

  33  format(a14,i2)

      return 
      end

C----------------------------------------------------------------------
      subroutine cem_out_fields2
C----------------------------------------------------------------------
c...  this works OK with param(103)=2: needs to be updated
c...  by completely replacing old subroutine cem_out_fields
c...  09/11/2006 misun

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      
      integer    fnum
      parameter (fnum=111)
      logical    outex, outey, outez, outhx, outhy, outhz
      integer    n

      character*33 filename
      data         filename /'./vtk/ascii-NN-p000000-t00000.vtk'/
      CHARACTER*1  fnarray(33)
      equivalence  (filename,fnarray)

      integer   dumpno
      integer   icalld
      save      icalld,dumpno
      data      icalld /0/, dumpno /0/

      character*1  numrl(0:9)
      data numrl /'0','1','2','3','4','5','6','7','8','9'/

      fnarray(29)= numrl(mod(dumpno,10))
      fnarray(28)= numrl(mod(dumpno,100)/10)
      fnarray(27)= numrl(mod(dumpno,1000)/100)
      fnarray(26)= numrl(mod(dumpno,10000)/1000)
      fnarray(25)= numrl(mod(dumpno,100000)/10000)

      fnarray(22)= numrl(mod(nid,10))
      fnarray(21)= numrl(mod(nid,100)/10)
      fnarray(20)= numrl(mod(nid,1000)/100)
      fnarray(19)= numrl(mod(nid,10000)/1000)
      fnarray(18)= numrl(mod(nid,100000)/10000)
      fnarray(17)= numrl(mod(nid,1000000)/100000)

      n = nx1*ny1*nz1*nelt

c ... VTK FORMAT: ASCII

      open(fnum, file=filename)      
      call vtk_dump_header2(fnum)

      outex = exactive
      outey = eyactive
      outez = ezactive
      outhx = hxactive
      outhy = hyactive
      outhz = hzactive

      if (ifscat) then

        call vtk_dump_field2(fnum,'totH  ',tothx,tothy,tothz,istep,time)
        call vtk_dump_field2(fnum,'totE  ',totex,totey,totez,istep,time)
        call vtk_dump_field2(fnum,'scatH ',hx,hy,hz,istep,time)
        call vtk_dump_field2(fnum,'scatE ',ex,ey,ez,istep,time)
        call vtk_dump_field2(fnum,'incH  ',inchx,inchy,inchz,istep,time)
        call vtk_dump_field2(fnum,'incE  ',incex,incey,incez,istep,time)

      else

        call vtk_dump_field2(fnum,'H     ',hx,hy,hz,istep,time)
        call vtk_dump_field2(fnum,'E     ',ex,ey,ez,istep,time)

      endif

      call vtk_dump_field2
     $              (fnum,'region',faceif,faceif,faceif,istep,time)!51 region
      call vtk_dump_field2
     $              (fnum,'face  ',faceid,faceid,faceid,istep,time)!52 face line
   
      if (ifsol) then

        if (outhx) call sub3(errhx, hx, shx, n)
        if (outhy) call sub3(errhy, hy, shy, n)
        if (outhz) call sub3(errhz, hz, shz, n)
        if (outex) call sub3(errex, ex, sex, n)
        if (outey) call sub3(errey, ey, sey, n)
        if (outez) call sub3(errez, ez, sez, n)

       call vtk_dump_field2(fnum,'errH  ',errhx,errhy,errhz,istep,time)
       call vtk_dump_field2(fnum,'errE  ',errex,errey,errez,istep,time)
       call vtk_dump_field2(fnum,'solH  ',shx,shy,shz,istep,time)
       call vtk_dump_field2(fnum,'solE  ',sex,sey,sez,istep,time)
          
      endif

      close(fnum)

      if(nid.eq.0) call pvtk_xml2(dumpno)   

      dumpno = dumpno +1

      if ((nid.eq.0).or.(icalld.eq.0)) then
        write(6,*)
        write(6,*) ' IF3D =',if3d,', TM=',iftm, ', CENTRAL=',
     $       ifcentral ,', PEC=',ifpec,', RK4=',IFRK4, 
     $       ', DEALIAS=',ifdealias,', SCAT=',ifscat,', SFTF =',ifsftf,
     $       ', IFNM =',ifnm 
        write(6,*) ' IFSRC=', ifsrc,', TE=', ifte,', UPWIND =', 
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude,
     $       ', LORENTZ=',ifloren

        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=', 
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter 

        if (IFDEALIAS) write(6,33) 'Dealias_M  =',lxd
        if (IFEXP    ) write(6,33) 'Krylov_Dim =',marnol
        icalld =1 
        write(6,*) ' Direct dumping into VTK ASCII format :::::::'
      endif

  33  format(a14,i2)


      return 
      end

C----------------------------------------------------------------------
      subroutine cem_out_fields
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      
      logical outex, outey, outez, outhx, outhy, outhz
      integer n
      integer icalld
      save    icalld
      data    icalld /0/


      ifxyo=.true.
      if (istep.gt.iostep) ifxyo = .false.  ! dump geom for 1st file only
      ifvo  = .true.
      ifpo  = .false.
      ifto  = .false.

      n = nx1*ny1*nz1*nelt

        if (NSTEPS.eq.0) then
            call outpost(hx,hy,hz,pr,hx,'hxy')
            call outpost(ex,ey,ez,pr,ey,'exy')
            return
        endif

        if (ifsol) then

          call outpost(hx,hy,hz,pr,hx,'hxy')
          call outpost(ex,ey,ez,pr,ey,'exy')

        else 

          if (ifscat) then

              call outpost(tothx,tothy,tothz,pr,hx,'h-t')
              call outpost(totex,totey,totez,pr,ey,'e-t')
              call outpost(hx,hy,hz,pr,hx,'hxy')
              call outpost(ex,ey,ez,pr,ey,'exy')
              call outpost(engyh,engyh,engyh,pr,ey,'enh')
              call outpost(engye,engye,engye,pr,ex,'ene')

          else

             call outpost(hx,hy,hz,pr,hx,'hxy')
             call outpost(ex,ey,ez,pr,ey,'exy')

          endif ! ifscat

        endif ! ifsol


      if ((nid.eq.0).or.(icalld.eq.0)) then
        write(6,*)
        write(6,*) ' IF3D =',if3d,', TM=',iftm, ', CENTRAL=',
     $       ifcentral ,', PEC=',ifpec,', RK4=',IFRK4, 
     $       ', DEALIAS=',ifdealias,', SCAT=',ifscat,', SFTF =',ifsftf,
     $       ', IFNM =',ifnm 
        write(6,*) ' IFSRC=', ifsrc,', TE=', ifte,', UPWIND =', 
     $       ifupwind, ', PML=',ifpml,', EXP=', IFEXP,', DTCTL=',
     $       ifdtctl,', DIELEC=',ifdielec,', DRUDE=',ifdrude,
     $       ', LORENTZ=',ifloren
        write(6,*) ' IFSOL=',ifsol,', POISSON=',ifpoisson,', VARDT=', 
     $              ifvardt,', WAKE=',ifwake,', IFNAPOLY=',ifnapy,
     $              ', BEAM=', ifbeam,', FILTER=',iffilter 

        if (IFDEALIAS) write(6,33) 'Dealias_M  =',lxd
        if (IFEXP    ) write(6,33) 'Krylov_Dim =',marnol
        icalld =0 
      endif

      if (ifvtk) then
          write(6,*) ' Dumping into .dat format .......'
      else
          write(6,*) ' Dumping into postx format >>>>>>>'
      endif 


  33  format(a14,i2)

      return 
      end

C----------------------------------------------------------------------
      subroutine cem_out
C----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'DEALIAS'
      include 'EXPONE'

      n=nx1*ny1*nz1*nelt

      if (param(96).ne.0) then
        return
      endif

      if (iostep.gt.0) then

         if (mod(istep,iostep).eq.0) then
           call starttiming
           if     (param(103).eq.8) then  ! added by Jing Fu 2010-07-08
              call set_ascii_nmm
              call cem_out_fields6
	   elseif (param(103).eq.5) then  ! added by Jing Fu 2010-08-12
              call set_ascii_nm
              call cem_out_fields6
           elseif (param(103).eq.6) then  ! added by Jing Fu 2010-07-08
              call cem_out_fields6
           elseif (param(103).eq.-6) then ! added by Jing Fu 2010-07-16
              call set_ascii_true
              call cem_out_fields6
           elseif (param(103).eq.4) then  ! added by Jing Fu 2010-06-28
              call cem_out_fields4
           elseif (param(103).eq.3) then  ! temporary for testing
              call cem_out_fields3
           elseif (param(103).eq.2) then  ! temporary for testing
              call cem_out_fields2
           elseif (param(103).eq.1) then  ! temporary for testing
              call cem_out_fields
           elseif (param(103).eq.0) then  ! temporary for testing
              call cem_out_fields4
              call cem_out_fields3
              call cem_out_fields2
           endif
           call endtiming
           iparam=param(103)
           call writeiotrace(iparam, istep)
         endif


      elseif (istep.eq.nsteps) then

	call starttiming
           if     (param(103).eq.8) then  ! added by Jing Fu 2010-07-08
              call set_ascii_nmm
              call cem_out_fields6
	      call free_rbio_buffer
           elseif (param(103).eq.5) then  ! added by Jing Fu 2010-08-12
              call set_ascii_nm
              call cem_out_fields6
	      call free_rbio_buffer
           elseif (param(103).eq.6) then  ! added by Jing Fu 2010-07-08
              call cem_out_fields6
	      call free_rbio_buffer
           elseif (param(103).eq.-6) then ! added by Jing Fu 2010-07-16
              call set_ascii_true
              call cem_out_fields6
	      call free_rbio_buffer
           elseif (param(103).eq.4) then  ! added by Jing Fu 2010-06-28
              call cem_out_fields4
           elseif (param(103).eq.3) then  ! temporary for testing
              call cem_out_fields3
           elseif (param(103).eq.2) then  ! temporary for testing
              call cem_out_fields2
           elseif (param(103).eq.1) then  ! temporary for testing
              call cem_out_fields
           elseif (param(103).eq.0) then  ! temporary for testing
              call cem_out_fields4
              call cem_out_fields3
              call cem_out_fields2
           endif
           call endtiming
           iparam=param(103)
           call writeiotrace(iparam, istep)

      endif
c
      return
      end

c---------------------------------------------------------------
      subroutine cem_geombound
c---------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'GEOMBOUND'

      n = nx1*ny1*nz1*nelv
      xmin = glmin(xm1,n)
      xmax = glmax(xm1,n)
      ymin = glmin(ym1,n)
      ymax = glmax(ym1,n)
      zmin = glmin(zm1,n)
      zmax = glmax(zm1,n)

      xsize = xmax-xmin
      ysize = ymax-ymin
      zsize = zmax-zmin

      xcenter = 0.5*(xmin+xmax)
      ycenter = 0.5*(ymin+ymax)
      zcenter = 0.5*(zmin+zmax)

      return
      end
c---------------------------------------------------------------
      subroutine cem_init
c---------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'

      n    = nx1*ny1*nz1*nelv
      n6   = nx1*ny1*nz1*nelv*(3+3*(ndim-2))
      nsrf = nx1*nz1*nelv*2*ndim
      m    = marnol

      call rzero(shx,n)
      call rzero(shy,n)
      call rzero(shz,n)
      call rzero(sex,n)
      call rzero(sey,n)
      call rzero(sez,n)

      call rzero(hx,n)
      call rzero(hy,n)
      call rzero(hz,n)
      call rzero(ex,n)
      call rzero(ey,n)
      call rzero(ez,n)
      call rzero(jx,n)
      call rzero(jy,n)
      call rzero(jz,n)

      call rzero(tothx,n)
      call rzero(tothy,n)
      call rzero(tothz,n)
      call rzero(totex,n)
      call rzero(totey,n)
      call rzero(totez,n)

      call rzero(inchx,n)
      call rzero(inchy,n)
      call rzero(inchz,n)
      call rzero(incex,n)
      call rzero(incey,n)
      call rzero(incez,n)

      call rzero(poytSX,n)
      call rzero(poytSY,n)
      call rzero(poytSZ,n)
      call rzero(avgtSX,n)
      call rzero(avgtSY,n)
      call rzero(avgtSZ,n)

      call rzero(poySX,n)
      call rzero(poySY,n)
      call rzero(poySZ,n)
      call rzero(avgSX,n)
      call rzero(avgSY,n)
      call rzero(avgSZ,n)

      call rzero(poyISX,n)
      call rzero(poyISY,n)
      call rzero(poyISZ,n)
      call rzero(avgISX,n)
      call rzero(avgISY,n)
      call rzero(avgISZ,n)

      call rzero(engyE,n)
      call rzero(engyH,n)
      call rzero(engyIE,n)
      call rzero(engyIH,n)

      call rzero(ww1,NSRF)
      call rzero(ww2,NSRF)
      call rzero(ww3,NSRF)
      call rzero(uu1,NSRF)
      call rzero(uu2,NSRF)
      call rzero(uu3,NSRF)

      call rzero(Fhx,NSRF)
      call rzero(Fhy,NSRF)
      call rzero(Fhz,NSRF)
      call rzero(Fex,NSRF)
      call rzero(Fey,NSRF)
      call rzero(Fez,NSRF)
      call rzero(Fjx,NSRF)
      call rzero(Fjy,NSRF)
      call rzero(Fjz,NSRF)
      call rzero(fincex,NSRF)
      call rzero(fincey,NSRF)
      call rzero(fincez,NSRF)
      call rzero(finchx,NSRF)
      call rzero(finchy,NSRF)
      call rzero(finchz,NSRF)
      call rzero(srftmp,6*NSRF)
      call rzero(srfinc,6*NSRF)
      call rzero(srfin ,6*NSRF)
      call rzero(srfout,6*NSRF)
      call rzero(srflxi,6*NSRF)
      call rzero(srflxo,6*NSRF)

      call rzero(Ghx,NSRF)
      call rzero(Ghy,NSRF)
      call rzero(Ghz,NSRF)
      call rzero(Gex,NSRF)
      call rzero(Gey,NSRF)
      call rzero(Gez,NSRF)
      call rzero(Gjx,NSRF)
      call rzero(Gjy,NSRF)
      call rzero(Gjz,NSRF)

      call rzero(Khx,n)
      call rzero(Khy,n)
      call rzero(Khz,n)
      call rzero(Kex,n)
      call rzero(Key,n)
      call rzero(Kez,n)
      call rzero(Kjx,n)
      call rzero(Kjy,n)
      call rzero(Kjz,n)

      call rzero(Kbx,n)
      call rzero(Kby,n)
      call rzero(Kbz,n)
      call rzero(Kdx,n)
      call rzero(Kdy,n)
      call rzero(Kdz,n)

      call rzero(REShx,n)
      call rzero(REShy,n)
      call rzero(REShz,n)
      call rzero(RESex,n)
      call rzero(RESey,n)
      call rzero(RESez,n)
      call rzero(RESjx,n)
      call rzero(RESjy,n)
      call rzero(RESjz,n)

      call rzero(px,3*n)
      call rzero(py,3*n)
      call rzero(pz,3*n)
      call rzero(qx,3*n)
      call rzero(qy,3*n)
      call rzero(qz,3*n)

      call rzero(Kpx,3*n)
      call rzero(Kpy,3*n)
      call rzero(Kpz,3*n)
      call rzero(Kqx,3*n)
      call rzero(Kqy,3*n)
      call rzero(Kqz,3*n)

      call rzero(RESpx,3*n)
      call rzero(RESpy,3*n)
      call rzero(RESpz,3*n)
      call rzero(RESqx,3*n)
      call rzero(RESqy,3*n)
      call rzero(RESqz,3*n)

      call rzero(lorena,3*n)
      call rzero(lorenb,3*n)
      call rzero(lorenc,3*n)

      !call izero(if_in_region,nelt)
      call rzero(epsilo,n)
      call rzero(mu,    n)
      call rzero(impede,n)
      call rzero(conduc,n)
      call rzero(drudea,n)
      call rzero(drudeb,n)

      call rzero(scat1,n)
      call rzero(scat2,n)
      call rzero(baseinc1,n)
      call rzero(baseinc2,n)
      call rzero(base11,n)
      call rzero(base12,n)
      call rzero(base21,n)
      call rzero(base22,n)

      call rzero(errex,n)
      call rzero(errey,n)
      call rzero(errez,n)
      call rzero(errhx,n)
      call rzero(errhy,n)
      call rzero(errhz,n)

      nxyzl=nx1*ny1*nz1*nelt*ndim 
      call rzero(pmlsigma,nxyzl)
      call rzero(pmldx,n)
      call rzero(pmldy,n)
      call rzero(pmldz,n)
      call rzero(pmlbx,n)
      call rzero(pmlby,n)
      call rzero(pmlbz,n)

      call rzero(respmldx,n)                      
      call rzero(respmldy,n)                      
      call rzero(respmldz,n)                      
      call rzero(respmlbx,n)                      
      call rzero(respmlby,n)                      
      call rzero(respmlbz,n)                      

      call rzero(respmlex,n)                      
      call rzero(respmley,n)                      
      call rzero(respmlez,n)                      
      call rzero(respmlhx,n)                      
      call rzero(respmlhy,n)                      
      call rzero(respmlhz,n)                      

      call rzero(kpmlbx,n)
      call rzero(kpmlby,n)
      call rzero(kpmlbz,n)
      call rzero(kpmldx,n)
      call rzero(kpmldy,n)
      call rzero(kpmldz,n)

      call izero(pmlptr,nelt)

c     call rzero(arnolh0,(m+1)*m)
c     call rzero(arnolh ,    m*m)
c      call rzero(arnolrh,    m*m)
c      call rzero(arnolch,  2*m*m)
c     call rzero(arnole,m)
c     call rzero(arnolf,m)

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_mater2
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer i, n
      real    tmpmu, tmpep, tmpt1, tmpt2

      n   = nx1*ny1*nz1*nelt

      do i=1,n 

         tmpmu = MU      (i,1,1,1)                      
         tmpep = EPSILO  (i,1,1,1)
         tmpt1 = sqrt(tmpmu/tmpep)
         tmpt2 = sqrt(tmpep/tmpmu)

         IMPEDE(i,1,1,1) = tmpt1                
         CONDUC(i,1,1,1) = tmpt2                                       

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_mater3 
c
c     Restrict u to faces
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
 
      integer nxyzf
      integer icalld
      save    icalld
      data    icalld  /0/

      real const
      character CB*3
      integer   e,f,i,j

      nxyzf = nx1*nz1*2*ndim*nelt
      nface = 2*ndim
      nxz   = nx1*nz1

      if (icalld.eq.0) then
         call cem_set_fc_ptr
         icalld = 1
      endif

      nxyzf = nx1*nz1*2*ndim*nelt

      do j=1,ncemface
         i              = cemface(j)
         ZIMPEDE(j,1,1) = IMPEDE(i,1,1,1)
         YCONDUC(j,1,1) = CONDUC(i,1,1,1)
         FACEF  (j,1,1) = FACEIF(i,1,1,1) 
      enddo

      call copy (Z_0,ZIMPEDE,nxyzf)
      call copy (Y_0,YCONDUC,nxyzf)

      call gs_op_fields(gsh_face,Z_0,nxyzf,1,1,1,0)
      call gs_op_fields(gsh_face,Y_0,nxyzf,1,1,1,0)

      do e= 1,nelt            
      do f= 1,nface
      do j= 1,nxz

         i   = (e-1)*nxz*nface + nxz*(f-1) + j  ! face numbering

         CB = CBC(f,e,2)

         if (CB.eq.'PEC') then
             Y_0(i,1,1) = 2.0*Y_0(i,1,1)
             Y_1(i,1,1) = 2.0*Y_1(i,1,1)
             Z_0(i,1,1) = 2.0*Z_0(i,1,1)
             Z_1(i,1,1) = 2.0*Z_1(i,1,1)
         endif

      enddo
      enddo
      enddo

      call sub3 (Z_1,Z_0,ZIMPEDE,NXYZF)
      call sub3 (Y_1,Y_0,YCONDUC,NXYZF)

      const = 0.5
 
      call cmult (Z_0,const,nxyzf)
      call cmult (Y_0,const,nxyzf)

      return
      end
c-----------------------------------------------------------------------
c parallel checked!
      subroutine beam_setup 
      include 'SIZE'
      include 'TOTAL'   
      include 'EMWAVE'

      ! set arrays for curarry, rarry0, and zarry0              
      call arrays_setup  

      if (nid.eq.0) write(6,*) '.'
      if (nid.eq.0) write(6,7)    

      ! this is preset for wake zone including napoly curve zone
      ! use parameters
      !rsigma = param(88)
      !zsigma = param(88)
      rsigma = param(70)
      zsigma = param(71)
      zwidth = param(72)*zsigma
      zfrst  = zrr0 (2 )       !=glmin(zarry0(2),1): global min        
      ztmp1  = zfrst+zwidth

      ! this is preset for zposition, zcut                          
      icalld1=0
      icalld2=0

      do i = 1,nzrr0   
         if (zrr0(i).ge.ztmp1) then
             if (icalld1.EQ.0) then
                 zposition = zrr0(i)
                 izposition= i
                 icalld1=1
             endif
         endif
      enddo

      ztmp2= zposition+zwidth
      do i = 1,nzrr0   
         if (zrr0(i).ge.ztmp2) then
             if (icalld2.EQ.0) then
                 zcut = zrr0(i)
                 izcut= i
                 icalld2=1
             endif
         endif
      enddo

      if (nid.eq.0) write(6,10) zposition
      if (nid.eq.0) write(6,20) zcut                         

      ! this is to divide into three zones                                                    
      ! wnapy1 ---- wzone1 ---- wzone2 ---- wnapy2                                               
      ! local info
      iwake1 =izcut+1
      iwake2 =nzrr0-izcut
      wnapy1 =zrr0(iwake1-1) ! see wake.F
      wnapy2 =zrr0(iwake2+1)

      if (.not.(ifmovw)) then 
      if (iwake1-1.ge.iwake2+1) then
        if (nid.eq.0) then
            write(6,*) 'ERR: wake zone is zero: increase z levels' 
            write(6,*) 'ERR: or reduce the size of the bunch' 
            write(6,*) 'ERR: wnapy1, wnapy2',wnapy1,wnapy2 
        endif
        call exitt
      endif
      endif

      wzone1 =zrr0(iwake1)   ! see wake.F
      wzone2 =zrr0(iwake2)

      !self-checking for num of elts with z <= zcut
      !call radial_elts_num_chk  FIXME: this isn't correct chk          
      !but do I need this nel_zarry ??? not really so far.                                
       
      !err chks to be done: size of tubes, radial num of elts 08/17/07
      if (nid.eq.0) write(6,30) wzone1,wzone2                        
      if (nid.eq.0) write(6,40) wnapy1,wzone1,wzone2,wnapy2            

      ! find elets in napoly curve zone
      call find_napoly_curve_elts
      call get_napoly_arrays       

      ! assign radius for straight path
      if (IFWAKET) then
          wakerad1=glmin(rarry0(2),1)
          wakerad2=glmin(rarry0(3),1)
      else
         if (IFNAPY) then
             wakerad1=glmin(rarry0(2),1)
         else
             wakerad1=glmax(wnprad1,2*nelt)
             !call get_rad_no_napoly(wakerad) : temporarily blocked for 
             !moving window test -- 04/27/09 
             wakerad1=max(wakerad,wakerad1)       
             wakerad1=glmin(rarry0(2),1)  !TEMP
         endif
      endif

      wakerad1=param(74)
      if (nid.eq.0) write(6,9) wakerad1,wakerad2

   7  format ( ' Beam and Wakepotential setups::')
   8  format ( ' wakepothen path at radius=',f15.7)
   9  format ( ' wakepoten path at radius=',f15.7,' and radius=',f15.7)
  10  format ( ' assigned zposition =',f15.7)
  20  format ( ' assigned zcut      =',f15.7)
  30  format ( ' assigned wake mid zone =[',f15.7,',',f15.7,']')
  40  format ( ' assigned napoly curve zone =[' 
     $           ,f15.7,',',f15.7,'] :: [',f15.7,',',f15.7,']')

      return
      end
c-----------------------------------------------------------------------
c this is not used: not correct checking
      subroutine radial_elts_num_chk               
      include 'SIZE'
      include 'TOTAL'   
      include 'EMWAVE'
      include 'WAKEPOT'

      nxyz=nx1*ny1*nz1

      lnum=nel_zarry(1)
      do i=1,nzarry0
         if (zarry0(i).le.wnapy1) then
             l=nel_zarry(i)
             if (l.eq.lnum) then
             else
                 if (nid.eq.0) write(6,*) 'rad elt mismatch 1',l,lnum
             endif
         endif
      enddo

      lnum=nel_zarry(nzarry0-1)
      do i=1,nzarry0-1
         if (zarry0(i).ge.wzone2) then
             l=nel_zarry(i)
             if (l.eq.lnum) then
             else
                 if (nid.eq.0) write(6,*) 'rad elt mismatch 2',l,lnum
             endif
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
c ... parallel checked!
c     BEAM, WAKE and NAPOLY-CONTOUR SETUP

      subroutine arrays_setup 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  i,e,f

      call find_curve_elts    ! curvelt=1 (1-side), -1 (2-side), 0(none) 
      call get_rarrays        ! radius including one-side curve cylinder         
      call get_curarrays      ! radius only two-sided curve for cylinder
      call get_zarrays        ! z values at interfaces for each proc.    
      call get_glo_ptr_zarrays! zarry zz0 in 1D globally assigned.          
      call get_glo_nel_zarrays! get num of transverse elts zarry0[k,k+1]
      call chk_direction           

      return
      end
c-----------------------------------------------------------------------
c ... this determines direction of each element with two curved sides.
c     This works only for cylindrical meshes (from n2to3). Thus face 5
c     and 6 is not considered. If iedg with curves = (2,6) (4,8), node 
c     is assigned in a radial direction 'R'. If edge with curves = (1,5)
c     (3,7), nodes are assigned in theta direction 'T'
c...  this is for interpolation in x-y direction for wakepotentials.

      subroutine chk_direction
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'

      integer e

      do e=1,nelt
         eltdir(e)='NOCUR'             
         if (curvelt(e).eq.-1) then ! check elt with two-curve sides  
            do iedg=1,8
               if (mod(iedg,2).EQ.0) then
                   if (ccurve(iedg,e).EQ.'C') then
                   eltdir(e)='RADIS'             
                   !write(6,*) e,eltdir(e)
                   endif
               else
                   if (ccurve(iedg,e).EQ.'C') then
                   eltdir(e)='THETA'             
                   !write(6,*) e,eltdir(e)
                   endif
               endif
            enddo

         endif

      enddo

      return
      end
c-----------------------------------------------------------------------
c     parallel
c     nel_zarry (k): 
c     The number of transverse elts within [zarry0(k),zarry0(k+1)] 
c     izelev(e): level (or location) of the elt 'e' within zarry0 frame.
 
      subroutine get_glo_nel_zarrays                                    
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'
    
      integer e

      call izero(nel_zarry_glo,nlevg)

      nxyz = nx1*ny1*nz1
      do e=1,nelt
         zavg= vlsum (zm1(1,1,1,e),nxyz)
         zavg= zavg/nxyz
         do k= 1,nzarry0-1
            if (zarry0(k).le.zavg.and.zavg.lt.zarry0(k+1)) then
                kg= izptr(k)
                nel_zarry_glo(kg)= nel_zarry_glo(kg)+1
                izelev(e)= kg    ! assign the level of e
                !write(6,*) 'nel',k,nel_zarry_glo(k),e,izelev(e)      
            endif
         enddo
      enddo

      do k=1,nlevg 
         nel_zarry_glo(k) = iglsum(nel_zarry_glo(k),1)
        !if (nid.eq.3) write(6,*) 'nel_zarry_glo',k,nel_zarry_glo(k)
        !write(6,*) 'nel_zarry_glo',node,k,nel_zarry_glo(k)
      enddo

      return
      end
c-----------------------------------------------------------------------
c     parallel
c     nel_zarry (k): 
c     The number of transverse elts within [zarry0(k),zarry0(k+1)] 
c     izelev(e): level (or location) of the elt 'e' within zarry0 frame.
c
c     FIXME : this is not correct implementation

      subroutine get_nel_zarray                                
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'
    
      integer e

      nxyz = nx1*ny1*nz1

      !write(6,*) 'nzarry0=',nzarry0
      call izero(nel_zarry,nzarry0)

      do e=1,nelt
         zavg = vlsum (zm1(1,1,1,e),nxyz)
         zavg = zavg/nxyz
         !write(6,*) 'zavg',e,zavg,zavg/nxyz

         do k=1,nzarry0-1
            if (zarry0(k).le.zavg.and.zavg.lt.zarry0(k+1)) then
                nel_zarry(k)= nel_zarry(k)+1
                izelev (e)  = k     ! assign the level of e
                !write(6,*) 'nel',k,nel_zarry(k),e,izelev(e)          
            endif
         enddo

      enddo

      !do k=1,nzarry0-1
      !   nel_zarry(k) = iglsum(nel_zarry(k),1)
      !   !write(6,*) 'nel_zarry',k,nel_zarry(k)
      !enddo

      return
      end
c-----------------------------------------------------------------------
c ... parallel
      subroutine get_napoly_arrays      
      include 'SIZE'
      include 'TOTAL'   
      include 'EMWAVE'
      integer  e

      !if (.NOT.IFNAPY) return
  
      rmax=-1.e+20
      rmin= 1.e+20
      tmax=-1.e+20
      tmin= 1.e+20

      if (nid.eq.0) write(6,*) '.'

      nxyz = nx1*ny1*nz1
      do k= 1,nwnpel(1)
         e= wnpel(1,k)        
         do i=1,nxyz
            xx=xm1(i,1,1,e)
            yy=ym1(i,1,1,e)
            rr=xx**2+yy**2  
            theta=atan2(yy,xx)
            if (rr.gt.0) rr=sqrt(rr)
            rmax= max(   rr,rmax)
            rmin= min(   rr,rmin)
            tmax= max(theta,tmax)
            tmin= min(theta,tmin)
         enddo
         wnprad1(1,e)=rmin
         wnprad1(2,e)=rmax
         wnptht1(1,e)=tmin
         wnptht1(2,e)=tmax
         if (nid.eq.0) write(6,10) e,rmin,rmax,tmax,tmin
      enddo

      rmax=-1.e+20
      rmin= 1.e+20
      tmax=-1.e+20
      tmin= 1.e+20
      nxyz = nx1*ny1*nz1
      do k= 1,nwnpel(2)
         e= wnpel(2,k)        
         do i=1,nxyz
            xx=xm1(i,1,1,e)
            yy=ym1(i,1,1,e)
            rr=xx**2+yy**2  
            theta=atan2(yy,xx)
            if (rr.gt.0) rr=sqrt(rr)
            rmax= max(   rr,rmax)
            rmin= min(   rr,rmin)
            tmax= max(theta,tmax)
            tmin= min(theta,tmin)
         enddo
         wnprad2(1,e)=rmin
         wnprad2(2,e)=rmax
         wnptht2(1,e)=tmin
         wnptht2(2,e)=tmax
         if (nid.eq.0) write(6,20) e,rmin,rmax,tmin,tmax
      enddo

      rmax=-1.e+20
      rmin= 1.e+20
      tmax=-1.e+20
      tmin= 1.e+20
      nxyz = nx1*ny1*nz1
      do k= 1,nwnpel(3)
         e= wnpel(3,k)        
         do i=1,nxyz
            xx=xm1(i,1,1,e)
            yy=ym1(i,1,1,e)
            rr=xx**2+yy**2  
            theta=atan2(yy,xx)
            if (rr.gt.0) rr=sqrt(rr)
            rmax= max(   rr,rmax)
            rmin= min(   rr,rmin)
            tmax= max(theta,tmax)
            tmin= min(theta,tmin)
         enddo
         wnprad2(2,e)=rmax
         wnptht2(1,e)=tmin
         wnptht2(2,e)=tmax
         if (nid.eq.0) write(6,20) e,rmin,rmax,tmin,tmax
      enddo

      if (nid.eq.0) write(6,*) '.'

  10  format ('  rmin:rmax:tmin:tmax in napoly elts 1=',i8,4f15.7)
  20  format ('  rmin:rmax:tmin:tmax in napoly elts 2=',i8,4f15.7)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_rad_no_napoly(rmax)             
      include 'SIZE'
      include 'TOTAL'   
      include 'EMWAVE'
      real     rmax
      integer  e

      rmax=-1.e+20
      rmin= 1.e+20
      tmax=-1.e+20
      tmin= 1.e+20
      nxyz = nx1*ny1*nz1
      do k= 1,nwnpel(3)
         e= wnpel(3,k)        
         do i=1,nxyz
            xx=xm1(i,1,1,e)
            yy=ym1(i,1,1,e)
            rr=xx**2+yy**2  
            theta=atan2(yy,xx)
            if (rr.gt.0) rr=sqrt(rr)
            rmax= max(   rr,rmax)
         enddo
      enddo
      rmax =glmax(rmax,1)
      rmax0=rmax

      rmax=-1.e+20
      rmin= 1.e+20
      tmax=-1.e+20
      tmin= 1.e+20
      nxyz = nx1*ny1*nz1
      do k= 1,nwnpel(4)
         e= wnpel(4,k)        
         do i=1,nxyz
            xx=xm1(i,1,1,e)
            yy=ym1(i,1,1,e)
            rr=xx**2+yy**2  
            theta=atan2(yy,xx)
            if (rr.gt.0) rr=sqrt(rr)
            rmax= max(   rr,rmax)
         enddo
      enddo
      rmax=glmax(rmax,1)


      if (abs(rmax-rmax0).GT.1.e-5) then
          write(6,*) 'Err: one-sided curve in napoly',rmax0,rmax
          ! should be only one-sided curve and no napoly element
          call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
c ... parallel
      subroutine find_napoly_curve_elts
      include 'SIZE'
      include 'TOTAL'   
      include 'EMWAVE'

      integer  radmin


      ! if (.NOT.IFNAPY) return

      ! self-check if there are minimum two curve radius.
      if (IFNAPY) then
          nelt2=2*nelt
          radmin=iglmin(curvelt,nelt2) 
          if (radmin.gt. -1) then
              if (nid.eq.0) then
              write(6,*) 'ERR:: mesh should have two inner-cylinders 
     $                          to use Napoly contours'
              endif
              call exitt
          endif
      endif

      ! find two-sided curve element within napoly curve zone  
      if (nid.eq.0) write(6,*) '.'
      if (nid.eq.0) write(6,*) 'find two-sided curve elements: '

      k1=0
      k2=0
      k3=0
      k4=0
      nxyz=nx1*ny1*nz1
      do iel=1,nelt
         if (curvelt(iel).EQ.-1) then       ! two-sided curve
             zlmax = vlmax(zm1(1,1,1,iel),nxyz)
             zlmin = vlmin(zm1(1,1,1,iel),nxyz)
             zlavg = 0.5*(zlmax+zlmin)
             if     ((zlavg.GT.wnapy1).AND.(zlavg.LT.wzone1)) then
                  k1=k1+1
                  wnpel(1,k1)=iel
                  nwnpel(1)=k1
                  if (nid.eq.0) write(6,10) k1,iel
             elseif ((zlavg.LT.wnapy2).AND.(zlavg.GT.wzone2)) then
                  k2=k2+1
                  wnpel(2,k2)=iel
                  nwnpel(2)=k2
                  if (nid.eq.0) write(6,20) k2,iel
             endif
         elseif (curvelt(iel).EQ. 1) then   ! one-sided curve
             zlmax = vlmax(zm1(1,1,1,iel),nxyz)
             zlmin = vlmin(zm1(1,1,1,iel),nxyz)
             zlavg = 0.5*(zlmax+zlmin)
             if     ((zlavg.GT.wnapy1).AND.(zlavg.LT.wzone1)) then
                  k3=k3+1
                  wnpel(3,k3)=iel
                  nwnpel(3)=k3
                  if (nid.eq.0) write(6,30) k3,iel
             elseif ((zlavg.LT.wnapy2).AND.(zlavg.GT.wzone2)) then
                  k4=k4+1
                  wnpel(4,k4)=iel
                  nwnpel(4)=k4
                  if (nid.eq.0) write(6,40) k4,iel
             endif
         endif
      enddo
      if (nid.eq.0) write(6,*) '.'

      ! nwnpel can be different for each node 
      !write(6,*) 'nwnpel--1',nwnpel(1)
      !write(6,*) 'nwnpel--2',nwnpel(2)
      !write(6,*) 'nwnpel--3',nwnpel(3)
      !write(6,*) 'nwnpel--4',nwnpel(4)
      !call exitt

  10  format ('  napoly curve zone 1: elt=',2i8)
  20  format ('  napoly curve zone 2: elt=',2i8)
  30  format ('  no napoly curve zone 1: elt=',2i8)
  40  format ('  no napoly curve zone 2: elt=',2i8)

      if (nid.eq.0) write(6,*) '.'

      return
      end
c-----------------------------------------------------------------------
c     NAPOLY-CONTOUR SETUP
c
c ... Find elements with 2 curve-sided:            
c     To support napoly path for wake potential calculation
c     Works only for n2to3 mesh (not for spherical curve)

      subroutine find_curve_elts
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'EMWAVE'
      real     rtmp0, rtmp1
    
      call izero(curvelt,nelt)

      do iel=1,nelt
         do iedg=1,4    ! caution the index is from 1 to 2 
            if (ccurve(  iedg,iel).EQ.'C') then 
                if     (iedg.EQ.1) then 
                    rtmp11=curve(1,1,iel) 
                    rtmp21=curve(1,3,iel) 
                elseif (iedg.EQ.2) then 
                    rtmp11=curve(1,2,iel) 
                    rtmp21=curve(1,4,iel) 
                elseif (iedg.EQ.3) then 
                    rtmp11=curve(1,3,iel) 
                    rtmp21=curve(1,1,iel) 
                elseif (iedg.EQ.4) then 
                    rtmp11=curve(1,4,iel) 
                    rtmp21=curve(1,2,iel) 
                endif
                if     ((rtmp11*rtmp21.LT.0)) then ! cylindrically 2-side curve          
                         curvelt(iel)=-1   
                elseif ((rtmp11*rtmp21.GT.0)) then !                                     
                      if (nid.eq.0) write(6,*) 'ERR: no such curve elts'
                else
                    if ((rtmp11.GT.0).OR.(rtmp21.GT.0)) then ! cylindrically 1-side                
                         curvelt(iel)= 1   
                    endif
                endif
            endif
         enddo

      enddo

      return
      end

c-----------------------------------------------------------------------
c     NAPOLY-CONTOUR SETUP
c
      subroutine get_rarrays
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      real     rmax, rmin  
      integer  e, nxyz

      nelt2=2*nelt
      call rzero (zin0  ,nelt2)
      call izero (iwork1,nelt2)

      nxyz=nx1*ny1*nz1

      do e=1,nelt           
         rmax=-1.0e+20
         rmin= 1.0e+20
         if     (curvelt(e).EQ. 1) then
            do i=1,nxyz
               xx=xm1(i,1,1,e)
               yy=ym1(i,1,1,e)
               rr=xx**2+yy**2
               if (rr.gt.0) rr=sqrt(rr)
               rmax=max(rr,rmax)
               rmin=min(rr,rmin)
            enddo
            rin0(1,e)=0   
            rin0(2,e)=rmax
         elseif (curvelt(e).EQ.-1) then
            do i=1,nxyz
               xx=xm1(i,1,1,e)
               yy=ym1(i,1,1,e)
               rr=xx**2+yy**2
               if (rr.gt.0) rr=sqrt(rr)
               rmax=max(rr,rmax)
               rmin=min(rr,rmin)
            enddo
            rin0(1,e)=rmin
            rin0(2,e)=rmax
         else
            rin0(1,e)=0      
            rin0(2,e)=0    
         endif
      enddo

      nsort = 2*nelt
      tol   = 1.e-5       ! FIXME 
      call reduced_sorts(rarry0,nrarry0,rin0,nsort,iwork1,tol)

      if (nid.eq.0) write(6,*) '.'              
      do i=1,nrarry0
         ! rarry(1)=0 : dummy
         if (nid.eq.0) write(6,10) i,rarry0(i)
      enddo
      if (nid.eq.0) write(6,*) '.'              

  10  format(' all curve radius array(',i5,')=',f15.7)

      return
      end
c-----------------------------------------------------------------------
c     NAPOLY-CONTOUR SETUP
c
      subroutine get_curarrays
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      real     rmax, rmin  
      integer  e, nxyz

      nxyz=nx1*ny1*nz1

      do e=1,nelt           
         rmax=-1.0e+20
         rmin= 1.0e+20
         if     (curvelt(e).EQ.-1) then
            do i=1,nxyz
               xx=xm1(i,1,1,e)
               yy=ym1(i,1,1,e)
               rr=xx**2+yy**2
               if (rr.gt.0) rr=sqrt(rr)
               rmax=max(rr,rmax)
               rmin=min(rr,rmin)
            enddo
            curveltr(1,e)=rmin
            curveltr(2,e)=rmax
            icalld = 1 
           !write(6,*) 'radius',iel,curveltr(1,e),curveltr(2,e)
         else
            curveltr(1,e)=0      
            curveltr(2,e)=0    
         endif
      enddo

      nsort = 2*nelt
      tol   = 1.e-5 
      call reduced_sorts(curarry,ncurarry,curveltr,nsort,icurwork,tol)

      ! curarry(1)=0 : dummy
      if (ncurarry.eq.1) then
         if (nid.eq.0)  write(6,9)                  
      else
         do i=1,ncurarry
           if (nid.eq.0)  write(6,10) i,curarry(i)
         enddo
      endif

   9  format(' no double-side curve radius array')
  10  format(' double-side curve radius array(',i5,')=',f15.7)

      return
      end
c-----------------------------------------------------------------------
      subroutine reduced_sorts(uout,n1,uin,n,iwork,tol)           
      integer    n1,n,i,kmax  
      integer    iwork(1)
      real       uout (1),uin(1)                    
      real       umin,tol

      call sorts(uin,uin,iwork,n)

      kmax   = 1
      umin   = uin(1)     
      uout(1)= umin
      n1     = kmax

      do i=2,n
         if ((uin(i)-uout(kmax)).gt.tol) then
            kmax=kmax+1
            uout(kmax)=uin(i)
            n1=kmax
         endif       
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_swap_across_faces(fld)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real fld    (LX1*LZ1,2*LDIM,LELT)
      real fldcpy (LX1*LZ1,2*LDIM,LELT)
      real origfld(LX1*LZ1,2*LDIM,LELT)

      call copy(origfld, fld, ncemface)
      call copy( fldcpy, fld, ncemface)
      call gs_op(gsh_face,fldcpy,1,1,0)

      ! the additional copy was necessary to avoid aliasing fld 
      ! fld = fldcpy - origfld

      call sub3 (fld, fldcpy, origfld, ncemface)

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_verify_comm

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real rndfld  (lx1*lz1,2*ldim,lelt)
      real rndcopy (lx1*lz1,2*ldim,lelt)
      real facexfld(lx1*lz1,2*ldim,lelt)
      real faceyfld(lx1*lz1,2*ldim,lelt)
      real facezfld(lx1*lz1,2*ldim,lelt)
      real fullxfld(lx1,ly1,lz1,lelt)
      real fullyfld(lx1,ly1,lz1,lelt)
      real fullzfld(lx1,ly1,lz1,lelt)
      real samefld (lx1,ly1,lz1,lelt)

      integer n,j,i

      n = nx1*ny1*nz1*nelv

      ! fill field with random garbage
      do j=1,ncemface
         rndfld (j,1,1) = zbqlu01(0.)
         rndcopy(j,1,1) = rndfld(j,1,1)
      enddo

      ! perform communication
      call gs_op(gsh_face,rndfld,1,1,0)
      !call cem_swap_across_faces(rndfld)

      ! compare to original, flag identical values
      call rzero(samefld, n)
      do j=1,ncemface
         i              = cemface(j)
         if (abs(rndfld(j,1,1)-rndcopy(j,1,1)).le.1e-12) then
           samefld(i,1,1,1) = 1
         endif
      enddo

      ! fill fields of coordinates of faces
      do j=1,ncemface
         i              = cemface(j)
         facexfld(j,1,1) = xm1(i,1,1,1)
         faceyfld(j,1,1) = ym1(i,1,1,1)
         facezfld(j,1,1) = zm1(i,1,1,1)
      enddo

      ! swap them
      call cem_swap_across_faces(facexfld)
      call cem_swap_across_faces(faceyfld)
      call cem_swap_across_faces(facezfld)

      call rzero(fullxfld, n)
      call rzero(fullyfld, n)
      call rzero(fullzfld, n)

      ! compute the difference, put in full-space field
      do j=1,ncemface
         i              = cemface(j)
         fullxfld(i,1,1,1) = facexfld(j,1,1) - xm1(i,1,1,1)
         fullyfld(i,1,1,1) = faceyfld(j,1,1) - ym1(i,1,1,1)
         fullzfld(i,1,1,1) = facezfld(j,1,1) - zm1(i,1,1,1)
      enddo

      ! and dump for visualization
      call vtk_dump_geometry
      call vtk_dump_field(fullxfld, 'CMVX  ', 0, 0.)
      call vtk_dump_field(fullyfld, 'CMVY  ', 0, 0.)
      call vtk_dump_field(fullzfld, 'CMVZ  ', 0, 0.)
      call vtk_dump_field(samefld, 'SMFX  ', 0, 0.)
      call vtk_dump_field(samefld, 'SMFY  ', 0, 0.)
      call vtk_dump_field(samefld, 'SMFZ  ', 0, 0.)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_mater1
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real     sum1, sum2, sum3, sum4, sum5, sum6
      real     old1, old2, old3, old4, old5, old6
      real     new1, new2, new3, new4, new5, new6
      integer  n,i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2
      integer icalld                                               
      save    icalld 
      data    icalld /0/

      n   = nx1*ny1*nz1

      do e=1,nelv

         do i=1,n
            FACEIF(i,1,1,e)= if_in_region(e)                     
         enddo
                
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE CEM_CURL2        
c---------------------------------------------------------------

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'EMWAVE'

      COMMON /SCRMG/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)

      call curl_3d(resEX,resEY,resEZ,HX,HY,HZ,.false.,w1,w2) ! check done
      call curl_3d(resHX,resHY,resHZ,EX,EY,EZ,.false.,w1,w2) ! check done
 
      n = nx1*ny1*nz1*nelt

      call chsign(resHX,n)
      call chsign(resHY,n)
      call chsign(resHZ,n)

      RETURN
      END

c-----------------------------------------------------------------------
      subroutine curl_3d(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
c
      include 'SIZE'
      include 'TOTAL'
c
      logical ifavg
c
      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
c
      ntot  = nx1*ny1*nz1*nelt

      call dudxyz(work1,u3,rym1,sym1,tym1,jacm1,1,2) ! work1=dw/dy ; work2=dv/dz
      call dudxyz(work2,u2,rzm1,szm1,tzm1,jacm1,1,3)
      call sub3  (w1,work1,work2,ntot)

      call dudxyz(work1,u1,rzm1,szm1,tzm1,jacm1,1,3) ! work1=du/dz ; work2=dw/dx
      call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
      call sub3  (w2,work1,work2,ntot)

      call dudxyz(work1,u2,rxm1,sxm1,txm1,jacm1,1,1) ! work1=dv/dx ; work2=du/dy
      call dudxyz(work2,u1,rym1,sym1,tym1,jacm1,1,2)
      call sub3  (w3,work1,work2,ntot)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_divergence  ! pointwise                     

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      n = nx1*ny1*nz1*nelt

      call cem_div(divex,divey,divez,ex,ey,ez)
      call cem_div(divhx,divhy,divhz,hx,hy,hz)

      do i=1,n
         divE(i,1,1,1)= divex(i,1,1,1)
     $                 +divey(i,1,1,1)
     $                 +divez(i,1,1,1)
         divH(i,1,1,1)= divhx(i,1,1,1)
     $                 +divhy(i,1,1,1)
     $                 +divhz(i,1,1,1)
         
         errDivE(i,1,1,1)=abs(divE(i,1,1,1)+chgq(i,1,1,1))
         errDivH(i,1,1,1)=abs(divH(i,1,1,1))
      enddo

      if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
        divEmax =glmax(errDivE,n)
        divHmax =glmax(errDivH,n)
        write(6,*) '.'
        write(6,*) ' divEmax',divEmax,', divHmax',divHmax 
        write(6,*) '.'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_curl ! weighted form of curl

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      COMMON /SCRMG/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)

      common /curlw/ swx(lx1,ly1,lz1,lelt)
     $             , swy(lx1,ly1,lz1,lelt)
     $             , swz(lx1,ly1,lz1,lelt)


      n = nx1*ny1*nz1*nelt
                        
      if (ifdealias) then

        if (if3d) then

          call wght_dcurl (resEX,resEY,resEZ,HX,HY,HZ,.false.,w1,w2)
          call wght_dcurl (resHX,resHY,resHZ,EX,EY,EZ,.false.,w1,w2)
          call chsign (resHX,n)
          call chsign (resHY,n)
          call chsign (resHZ,n)

        else

          if     (iftm) then
            call wght_dcurl(resHX,resHY,resEZ,HX,HY,EZ,.false.,w1,w2)
            call chsign (resHX,n)
            call chsign (resHY,n)
          elseif (ifte) then
            call wght_dcurl(resEX,resEY,resHZ,EX,EY,HZ,.false.,w1,w2)
            call chsign (resHZ,n)
          else
            write(6,*) 'IFTE or IFTM:: not defined'
            call exitt
          endif

        endif

      else 

        if (if3d) then

          call wght_curl  (resEX,resEY,resEZ,HX,HY,HZ,.false.,w1,w2)
          call wght_curl  (resHX,resHY,resHZ,EX,EY,EZ,.false.,w1,w2)
          call chsign (resHX,n)
          call chsign (resHY,n)
          call chsign (resHZ,n)

        else

          if     (iftm) then
            call wght_curl (resHX,resHY,resEZ,HX,HY,EZ,.false.,w1,w2)
            call chsign (resHX,n)
            call chsign (resHY,n)
          elseif (ifte) then
            call wght_curl (resEX,resEY,resHZ,EX,EY,HZ,.false.,w1,w2)
            call chsign (resHZ,n)
          else
            write(6,*) 'IFTE or IFTM:: not defined'
            call exitt
          endif

        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drude 
C---------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

c     integer  n
c     real     tmpex(lx1*ly1*lz1*lelt)
c     real     tmpey(lx1*ly1*lz1*lelt)
c     real     tmpez(lx1*ly1*lz1*lelt)
c
      N = NX1*NY1*NZ1*NELT

      call SUBCOL3(resEX,JX,BM1,N)   ! (curl H)x - Jx
      call SUBCOL3(resEY,JY,BM1,N)   ! (curl H)y - Jy
      call SUBCOL3(resEZ,JZ,BM1,N)   ! (curl H)z - Jz

      if (IFSCAT) then
          call add3 (tmpex,EX,incEX,N)     
          call add3 (tmpey,EY,incEY,N)
          call add3 (tmpez,EZ,incEZ,N)
      else
          call copy (tmpex,EX,N)     
          call copy (tmpey,EY,N)
          call copy (tmpez,EZ,N)
      endif

      ! no multiplication with bm1, thus no invmass for resJ
      call VDOT2(resJX,JX,tmpex,DRUDEA,DRUDEB,N) ! Jx_t= a*Jx+b*Ex
      call VDOT2(resJY,JY,tmpey,DRUDEA,DRUDEB,N) ! Jy_t= a*Jy+b*Ey 
      call VDOT2(resJZ,JZ,tmpez,DRUDEA,DRUDEB,N) ! Jz_t= a*Jz+b*Ez

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_lorentz
C---------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

c     !if (istep.eq.1) then
c       aa1= glmax(reshx,n)
c       aa2= glmax(reshy,n)
c       aa3= glmax(reshz,n)
c       bb1= glmax(resex,n)
c       bb2= glmax(resey,n)
c       bb3= glmax(resez,n)
c       cc1= glmax(resjx,n)
c       cc2= glmax(resjy,n)
c       cc3= glmax(resjz,n)
c       if (nid.eq.0) write(6,*) 'aa,bb',istep,aa1,aa2,aa3,
c    $                             bb1,bb2,bb3,cc1,cc2,cc3
c     !endif

      integer  n

      N = NX1*NY1*NZ1*NELT

      call SUBCOL3(resEX,PX(1,1,1,1,1),BM1,N)   ! (curl H)x - Jx-P1-P2
      call SUBCOL3(resEY,PY(1,1,1,1,1),BM1,N)   ! (curl H)y - Jy
      call SUBCOL3(resEZ,PZ(1,1,1,1,1),BM1,N)   ! (curl H)z - Jz

      call SUBCOL3(resEX,PX(1,1,1,1,2),BM1,N)   ! (curl H)x - Jx-P1-P2
      call SUBCOL3(resEY,PY(1,1,1,1,2),BM1,N)   ! (curl H)y - Jy
      call SUBCOL3(resEZ,PZ(1,1,1,1,2),BM1,N)   ! (curl H)z - Jz

      call SUBCOL3(resEX,PX(1,1,1,1,3),BM1,N)   ! (curl H)x - Jx-P1-P2
      call SUBCOL3(resEY,PY(1,1,1,1,3),BM1,N)   ! (curl H)y - Jy
      call SUBCOL3(resEZ,PZ(1,1,1,1,3),BM1,N)   ! (curl H)z - Jz

      if (IFSCAT) then  ! check this out
          call add3 (tmpex,EX,incEX,N)
          call add3 (tmpey,EY,incEY,N)
          call add3 (tmpez,EZ,incEZ,N)
      else
          call copy (tmpex,EX,N)
          call copy (tmpey,EY,N)
          call copy (tmpez,EZ,N)
      endif

      ! no multiplication with bm1, thus no invmass for resQ, resP
      call vdot3(resPX(1,1,1,1,1),PX(1,1,1,1,1),QX(1,1,1,1,1),tmpEX,
     $     LORENA(1,1,1,1,1),LORENB(1,1,1,1,1),LORENC(1,1,1,1,1),n)
      call vdot3(resPY(1,1,1,1,1),PY(1,1,1,1,1),QY(1,1,1,1,1),tmpEY,
     $     LORENA(1,1,1,1,1),LORENB(1,1,1,1,1),LORENC(1,1,1,1,1),n)
      call vdot3(resPZ(1,1,1,1,1),PZ(1,1,1,1,1),QZ(1,1,1,1,1),tmpEZ,
     $     LORENA(1,1,1,1,1),LORENB(1,1,1,1,1),LORENC(1,1,1,1,1),n)
      call copy(resQX(1,1,1,1,1),PX(1,1,1,1,1),n)
      call copy(resQY(1,1,1,1,1),PY(1,1,1,1,1),n)
      call copy(resQZ(1,1,1,1,1),PZ(1,1,1,1,1),n)

      call vdot3(resPX(1,1,1,1,2),PX(1,1,1,1,2),QX(1,1,1,1,2),tmpEX,
     $     LORENA(1,1,1,1,2),LORENB(1,1,1,1,2),LORENC(1,1,1,1,2),n)
      call vdot3(resPY(1,1,1,1,2),PY(1,1,1,1,2),QY(1,1,1,1,2),tmpEY,
     $     LORENA(1,1,1,1,2),LORENB(1,1,1,1,2),LORENC(1,1,1,1,2),n)
      call vdot3(resPZ(1,1,1,1,2),PZ(1,1,1,1,2),QZ(1,1,1,1,2),tmpEZ,
     $     LORENA(1,1,1,1,2),LORENB(1,1,1,1,2),LORENC(1,1,1,1,2),n)
      call copy(resQX(1,1,1,1,2),PX(1,1,1,1,2),n)
      call copy(resQY(1,1,1,1,2),PY(1,1,1,1,2),n)
      call copy(resQZ(1,1,1,1,2),PZ(1,1,1,1,2),n)

      call vdot3(resPX(1,1,1,1,3),PX(1,1,1,1,3),QX(1,1,1,1,3),tmpEX,
     $     LORENA(1,1,1,1,3),LORENB(1,1,1,1,3),LORENC(1,1,1,1,3),n)
      call vdot3(resPY(1,1,1,1,3),PY(1,1,1,1,3),QY(1,1,1,1,3),tmpEY,
     $     LORENA(1,1,1,1,3),LORENB(1,1,1,1,3),LORENC(1,1,1,1,3),n)
      call vdot3(resPZ(1,1,1,1,3),PZ(1,1,1,1,3),QZ(1,1,1,1,2),tmpEZ,
     $     LORENA(1,1,1,1,3),LORENB(1,1,1,1,3),LORENC(1,1,1,1,3),n)
      call copy(resQX(1,1,1,1,3),PX(1,1,1,1,3),n)
      call copy(resQY(1,1,1,1,3),PY(1,1,1,1,3),n)
      call copy(resQZ(1,1,1,1,3),PZ(1,1,1,1,3),n)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_invqmass
C---------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer icalld
      save    icalld
      data    icalld  /0/

      integer n 

      n = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call COL3(QBM1,EPSILO,BM1,n)
          call COL3(QBM2,MU    ,BM1,n)
          call invcol1(QBM1,n)
          call invcol1(QBM2,n)
          call invers2(QBM3,BM1,n)
          icalld = 1
      endif

      if (if3d) then
          call col2 (reshx,qbm2,n)
          call col2 (reshy,qbm2,n)
          call col2 (reshz,qbm2,n)
          call col2 (resex,qbm1,n)
          call col2 (resey,qbm1,n)
          call col2 (resez,qbm1,n)
      else
         if     (IFTE) then
           call col2 (reshz,qbm2,n)
           call col2 (resex,qbm1,n)
           call col2 (resey,qbm1,n)
         elseif (IFTM) then
           call col2 (reshx,qbm2,n)
           call col2 (reshy,qbm2,n)
           call col2 (resez,qbm1,n)
         endif
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_l2err
C---------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer  n
      real     l2err

      N = NX1*NY1*NZ1*NELT

      l2err =0.0
#if 0
      write(6,*) 'test 11'
      l2err = l2err+ GLSC3(resHX,QBM2,HX,N)
      l2err = l2err+ GLSC3(resHY,QBM2,HY,N)
      l2err = l2err+ GLSC3(resHZ,QBM2,HZ,N)
      l2err = l2err+ GLSC3(resEX,QBM1,EX,N)
      write(6,*) 'test 22'
      l2err = l2err+ GLSC3(resEY,QBM1,EY,N)
      l2err = l2err+ GLSC3(resEZ,QBM1,EZ,N)
      write(6,*) 'test 33'

      l2err = l2err+ GLSC3(HX,BM2,HX,N)
      l2err = l2err+ GLSC3(HY,BM2,HY,N)
      l2err = l2err+ GLSC3(HZ,BM2,HZ,N)
      l2err = l2err+ GLSC3(EX,BM1,EX,N)
      write(6,*) 'test 44'
      l2err = l2err+ GLSC3(EY,BM1,EY,N)
      l2err = l2err+ GLSC3(EZ,BM1,EZ,N)

      l2err = sqrt(l2err)
      energy=l2err
#endif
      return
      end


c-----------------------------------------------------------------------
      subroutine rk_c(i)
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      integer  i

      RKtime=time+dt*rk4c(i)
      if (i.ne.5) then
        nextsubdt = rk4c(i+1)
      else
        nextsubdt = rk4c(1)
      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine rk_ab(i)
c-----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'RK5'
      include 'PML'

      integer i

      n=nx1*ny1*nz1*nelt

c     kHX=rk4a(i)*kHX

      ca = rk4a(i)
      cb = rk4b(i)

      if (if3d) then

         call rk4_upd(hx,khx,reshx,cb,ca,dt,n)
         call rk4_upd(hy,khy,reshy,cb,ca,dt,n)
         call rk4_upd(hz,khz,reshz,cb,ca,dt,n)
         call rk4_upd(ex,kex,resex,cb,ca,dt,n)
         call rk4_upd(ey,key,resey,cb,ca,dt,n)
         call rk4_upd(ez,kez,resez,cb,ca,dt,n)

         if (ifdrude) then
         call rk4_upd(jx,kjx,resjx,cb,ca,dt,n)
         call rk4_upd(jy,kjy,resjy,cb,ca,dt,n)
         call rk4_upd(jz,kjz,resjz,cb,ca,dt,n)
         endif

         if (ifloren) then
         call rk4_upd(px,kpx,respx,cb,ca,dt,3*n)
         call rk4_upd(py,kpy,respy,cb,ca,dt,3*n)
         call rk4_upd(pz,kpz,respz,cb,ca,dt,3*n)
         call rk4_upd(qx,kqx,resqx,cb,ca,dt,3*n)
         call rk4_upd(qy,kqy,resqy,cb,ca,dt,3*n)
         call rk4_upd(qz,kqz,resqz,cb,ca,dt,3*n)
         endif

         if (ifpml) then
         call rk4_upd(pmlbx,kpmlbx,respmlbx,cb,ca,dt,n)
         call rk4_upd(pmlby,kpmlby,respmlby,cb,ca,dt,n)
         call rk4_upd(pmlbz,kpmlbz,respmlbz,cb,ca,dt,n)
         call rk4_upd(pmldx,kpmldx,respmldx,cb,ca,dt,n)
         call rk4_upd(pmldy,kpmldy,respmldy,cb,ca,dt,n)
         call rk4_upd(pmldz,kpmldz,respmldz,cb,ca,dt,n)
         endif

      else

        if     (IFTE) then
           call rk4_upd(hz,khz,reshz,cb,ca,dt,n)
           call rk4_upd(ex,kex,resex,cb,ca,dt,n)
           call rk4_upd(ey,key,resey,cb,ca,dt,n)
           if (ifdrude) then
               call rk4_upd(jx,kjx,resjx,cb,ca,dt,n)
               call rk4_upd(jy,kjy,resjy,cb,ca,dt,n)
           endif

           if (ifloren) then
              call rk4_upd(px,kpx,respx,cb,ca,dt,3*n)
              call rk4_upd(py,kpy,respy,cb,ca,dt,3*n)
              call rk4_upd(qx,kqx,resqx,cb,ca,dt,3*n)
              call rk4_upd(qy,kqy,resqy,cb,ca,dt,3*n)
           endif

           if (ifpml) then
               call rk4_upd(pmlbz,kpmlbz,respmlbz,cb,ca,dt,n)
               call rk4_upd(pmldx,kpmldx,respmldx,cb,ca,dt,n)
               call rk4_upd(pmldy,kpmldy,respmldy,cb,ca,dt,n)
           endif
        elseif (IFTM) then
           call rk4_upd(hx,khx,reshx,cb,ca,dt,n)
           call rk4_upd(hy,khy,reshy,cb,ca,dt,n)
           call rk4_upd(ez,kez,resez,cb,ca,dt,n)
           if (ifdrude) then
              call rk4_upd(jx,kjx,resjx,cb,ca,dt,n)
              call rk4_upd(jy,kjy,resjy,cb,ca,dt,n)
           endif

           if (ifloren) then
              call rk4_upd(px,kpx,respx,cb,ca,dt,3*n)
              call rk4_upd(py,kpy,respy,cb,ca,dt,3*n)
              call rk4_upd(qx,kqx,resqx,cb,ca,dt,3*n)
              call rk4_upd(qy,kqy,resqy,cb,ca,dt,3*n)
           endif

           if (ifpml) then
               call rk4_upd(pmlbx,kpmlbx,respmlbx,cb,ca,dt,n)
               call rk4_upd(pmlby,kpmlby,respmlby,cb,ca,dt,n)
               call rk4_upd(pmldz,kpmldz,respmldz,cb,ca,dt,n)
           endif

        endif

      endif

      return
      end


c-----------------------------------------------------------------------
      subroutine rk4_upd(h,kh,resh,cb,ca,dt,n)
c-----------------------------------------------------------------------
      real h(1),kh(1),resh(1)
c
      nresid = iand(n,3)

      do i=1,n-nresid,4

         kh(i  ) = ca*kh(i  ) + dt*resh(i  )
         kh(i+1) = ca*kh(i+1) + dt*resh(i+1)
         kh(i+2) = ca*kh(i+2) + dt*resh(i+2)
         kh(i+3) = ca*kh(i+3) + dt*resh(i+3)

         h (i  ) =     h(i  ) + cb*kh  (i  )
         h (i+1) =     h(i+1) + cb*kh  (i+1)
         h (i+2) =     h(i+2) + cb*kh  (i+2)
         h (i+3) =     h(i+3) + cb*kh  (i+3)

      enddo

      if (nresid.eq.0) then

         return

      elseif (nresid.eq.1) then
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n  ) =     h(n  ) + cb*kh  (n  )
      elseif (nresid.eq.2) then
         kh(n-1) = ca*kh(n-1) + dt*resh(n-1)
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n-1) =     h(n-1) + cb*kh  (n-1)
         h (n  ) =     h(n  ) + cb*kh  (n  )
      else
         kh(n-2) = ca*kh(n-2) + dt*resh(n-2)
         kh(n-1) = ca*kh(n-1) + dt*resh(n-1)
         kh(n  ) = ca*kh(n  ) + dt*resh(n  )
         h (n-2) =     h(n-2) + cb*kh  (n-2)
         h (n-1) =     h(n-1) + cb*kh  (n-1)
         h (n  ) =     h(n  ) + cb*kh  (n  )
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_storage

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
c-----------------------------------------------------------------------
      subroutine restart2
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      common  vtkuxyz1(3*lx1*ly1*lz1*lelt)
      common  vtkuxyz2(3*lx1*ly1*lz1*lelt)
      common  vtkuxyz3(3*lx1*ly1*lz1*lelt)
      real    vtkuxyz1,vtkuxyz2,vtkuxyz3
      integer dumpno,n
      real*8  dclock,t0

      n = nx1*ny1*nz1*nelt

      irstdumpno=param(94) ! declared in EMWAVE
      !irststep  = xxx   (timestep from vtk file header)
      !time      = xxx   (time from vtk file header)
      !dt        = xxx   (dt from vtk file header)
      !idumpno    = xxx   (dumpno from vtk file header)

      t0=dclock()
      if (nid.eq.0) write(6,*) 'start restart: param(94)=',param(94)
      
      if     (param(24).eq.4) then   ! added by Jing Fu 2010-06-28
         if  (nid.eq.0) write(6,*) 'start reading xxx.vtk format' 

         call openfile_read4
         call readheader4(irststep,time,dt,idump)
         call checkdumpno(idump,irstdumpno)
         call readnodes4(vtkuxyz1,n) 
         call readfield4(vtkuxyz2,n) 
         call readfield4(vtkuxyz3,n) 

         call closefile_read4

      elseif (param(24).eq. 5) then  
         if  (nid.eq.0) write(6,*) 'start reading xxx.vtk format' 

      elseif (param(24).eq. 6) then  
         if  (nid.eq.0) write(6,*) 'start reading xxx.vtk format' 

      elseif (param(24).eq.-6) then  
         if  (nid.eq.0) write(6,*) 'start reading xxx.vtk format' 

      elseif (param(24).eq. 8) then 
         if  (nid.eq.0) write(6,*) 'start reading xxx.vtk format' 

      endif

      if (nid.eq.0) write(6,*) 'end restart',dclock()-t0

      call save2vectors(xm1,ym1,zm1,vtkuxyz1,n)
      call save2vectors(hx,hy,hz,vtkuxyz2,n)
      call save2vectors(ex,ey,ez,vtkuxyz3,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine save2vectors (ux1,uy1,uz1,vtkxyz,n)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer i,k,n
      real    ux1(lx1,ly1,lz1,lelt)
      real    uy1(lx1,ly1,lz1,lelt)
      real    uz1(lx1,ly1,lz1,lelt)
      real    vtkxyz(3*lx1*ly1*lz1*lelt)

      !n=nx1*ny1*nz1*nelt

      k = 1
      do i = 1,n    
         ux1(i,1,1,1)=vtkxyz(k  )
         uy1(i,1,1,1)=vtkxyz(k+1)
         uz1(i,1,1,1)=vtkxyz(k+2)
         k = k+3
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine  checkdumpno(idump,irstdumpno)
      include 'SIZE'
      include 'TOTAL'
      
      integer     idump, irstdumpno

      if (idump.ne.irstdumpno) then
          if (nid.eq.0) write(6,*)
     $       'non matching dumpno for restart: set param(94)=',idump
          call exitt
      else
          if (nid.eq.0) write(6,*) 'irststep, rst_time, rst_dt, 
     $                  irstdumpno',irststep, time, dt, irstdumpno
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_inc                            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      call userinc                                           

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_source                         

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real     tt, xx, yy, zz, u
      integer  ifld, e, i, j, k, n

      n = nx1*ny1*nz1*nelt

      call usersrc (57,reshx,reshy,reshz,resex,resey,resez)

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_pmlsrf_to_face
c
c     Restrict u to faces
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'
      integer  e,ef,f
c
      call full2face(fincex, incex)
      call full2face(fincey, incey)
      call full2face(fincez, incez)
      call full2face(finchx, inchx)
      call full2face(finchy, inchy)
      call full2face(finchz, inchz)

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_charge                        

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real     tt, xx, yy, zz, u
      integer  ifld, e, i, j, k, n

      n = nx1*ny1*nz1*nelt

      call gradm1(resjx,resjy,resjz,jx)
      call gradm1(resjx,resjy,resjz,jy)
      call gradm1(resjx,resjy,resjz,jz)
      call add4(divJ,resjx,resjy,resjz,n)

      return
      end

c-----------------------------------------------------------------------
      subroutine rectangular_to_cylindrical(vr,vt,vz1,vx,vy,vz,theta)   
      real       vr,vt,vz1,vx,vy,vz,theta   

      vr= vx*cos(theta)+vy*sin(theta)
      vt=-vx*sin(theta)+vy*cos(theta)
      vz1= vz

      return
      end

c-----------------------------------------------------------------------
      subroutine cylindrical_to_rectangular(vx,vy,vz1,vr,vt,vz,theta)   
      real       vr,vt,vz1,vx,vy,vz,theta   

      vx= vr*cos(theta)-vt*sin(theta)
      vy= vr*sin(theta)+vt*cos(theta)
      vz1= vz

      return
      end
c-----------------------------------------------------------------------
      subroutine poynting_field                                

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      const = 1.0

      n=nx1*ny1*nz1*nelt

      do i=1,n 

       
         exx=ex(i,1,1,1)
         eyy=ey(i,1,1,1)
         ezz=ez(i,1,1,1)
         hxx=hx(i,1,1,1)
         hyy=hy(i,1,1,1)
         hzz=hz(i,1,1,1)

         poytSX(i,1,1,1)=const*(eyy*bzz-ezz*byy)
         poytSY(i,1,1,1)=const*(ezz*bxx-exx*bzz)
         poytSZ(i,1,1,1)=const*(exx*byy-eyy*bxx)

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine energy_field (u,v1,v2,v3,n)                 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real u(n),v1(n),v2(n),v3(n)
      real tmp

      const = 1.0

      do i=1,n 
         tmp = v1(i)**2+v2(i)**2+v3(i)**2
         tmp = const*tmp
         u(i)= sqrt (tmp)
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine total_field                           
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      n= nx1*ny1*nz1*nelt
      call sub3(totex,ex,incex,n)
      call sub3(totey,ey,incey,n)
      call sub3(totez,ez,incez,n)
      call sub3(tothx,hx,inchx,n)
      call sub3(tothy,hy,inchy,n)
      call sub3(tothz,hz,inchz,n)

      return
      end

c-----------------------------------------------------------------------
      subroutine scat_field (u,v,n)                 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real u(n)
      real v(n)

      integer nxyz,e,i,j

      nxyz = nx1*ny1*nz1

      do e = 1, nelt

         ! if (pmlel(e).eq.1)  then
             ! tmp = 0.0
             ! do i = 1,nxyz
                ! i2=nxyz*(e-1)+i
                ! u(i2) = u(i2)+tmp*v(i2)
             ! enddo
         ! else
             tmp = 1.0
             do i = 1,nxyz
                i2=nxyz*(e-1)+i
                u(i2) = u(i2)+tmp*v(i2)
             enddo
         ! endif

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine outfldrn (x,txt,kk)
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*7 txt
C
      do ie=1,nelv,4
         do iz=1,nz1,1
            if (iz.eq.1) write(6,106) txt,iz,ie,kk
            if (iz.gt.1) write(6,107) 
            i1 = ie+1
            do j=ny1,1,-1
               write(6,105) ((x(i,j,iz,ie+i0),i=1,nx1),i0=0,3)
            enddo
         enddo
      enddo
C
  107 FORMAT(' ')
  105 FORMAT(4(4f6.0,5x))
c 105 FORMAT(4f6.0,20x,4f6.0)
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',2I4,/,
     $           5X,'       X            ')
C
      return
      end

c-----------------------------------------------------------------------
      subroutine outmat (u,m,n,name6,e)

      include 'SIZE'

      real        u(m,n)
      character*6 name6
      integer     e     
c
c     Print out copies of a global matrix
c
      write(6,1) nid,m,n,name6,e
   1  format(//,3i6,'  Matrix:',2x,a6,i9,/)

      n15 = min(n,15)
      do i=1,m
         if (n15.gt.5) write(6,15) nid,name6,(u(i,j),j=1,n15)
         if (n15.le.5) write(6, 5) nid,name6,(u(i,j),j=1,n15)
      enddo
   5  format(i3,1x,a6,20(6f8.4,/,10x))
  15  format(i3,1x,a6,20(15f6.1,/,10x))
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_set_fc_ptr
c
c     Set up pointer to restrict u to faces ! NOTE: compact
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON2D'
c
      integer i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2

      nxyz  = nx1*ny1*nz1
      nxz   = nx1*nz1
      nface = 2*ndim
      nxzf  = nx1*nz1*nface ! red'd mod to area, unx, etc.

      k = 0

      do e=1,nelv
      do f=1,nface

         ef     = eface(f)
         js1    = skpdat(1,f)
         jf1    = skpdat(2,f)
         jskip1 = skpdat(3,f)
         js2    = skpdat(4,f)
         jf2    = skpdat(5,f)
         jskip2 = skpdat(6,f)

         i = 0
         do j2=js2,jf2,jskip2
         do j1=js1,jf1,jskip1

            i = i+1
            k = i+nxz*(ef-1)+nxzf*(e-1)           ! face   numbering
            cemface(k) = j1+nx1*(j2-1)+nxyz*(e-1) ! global numbering

            !write(6,*) 'fc',e,k,f,ef,srfout(k)!glo_num3d(cemface(k))
            ! glo_num(cemface(k))  get global numbering
            ! cemface gets the Ed's numbering at faces following pff order

         enddo
         enddo

      enddo
      enddo
      ncemface = nxzf*nelv

      return
      end
c-----------------------------------------------------------------------
      subroutine full2face(faceary, fullary)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     fullary(lx1,ly1,lz1,lelt)
      integer  i,j

      do j=1,ncemface
         i= cemface(j)
         faceary(j,1,1) = fullary(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine face2full(fullary, faceary)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     fullary(lx1,ly1,lz1,lelt)
      integer  i,j

      call rzero(fullary, lx1*ly1*lz1*lelt)

      do j= 1,ncemface
         i= cemface(j)
         fullary(i,1,1,1) = fullary(i,1,1,1)+faceary(j,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_restrict_to_face
c
c     Restrict u to faces
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) then
          call cem_set_fc_ptr
          icalld = 1
      endif
      if     (if3d) then
        do j=1,ncemface
          i= cemface(j)
          fhx(j,1,1) = hx(i,1,1,1)
          fhy(j,1,1) = hy(i,1,1,1)
          fhz(j,1,1) = hz(i,1,1,1)
          fex(j,1,1) = ex(i,1,1,1)
          fey(j,1,1) = ey(i,1,1,1)
          fez(j,1,1) = ez(i,1,1,1)
        enddo
      else
        if     (iftm) then
        do j=1,ncemface
           i= cemface(j)
           fhx(j,1,1) = hx(i,1,1,1)
           fhy(j,1,1) = hy(i,1,1,1)
           fez(j,1,1) = ez(i,1,1,1)
        enddo
        elseif (ifte) then
        do j=1,ncemface
           i= cemface(j)
           fex(j,1,1) = ex(i,1,1,1)
           fey(j,1,1) = ey(i,1,1,1)
           fhz(j,1,1) = hz(i,1,1,1)
        enddo
        endif
      endif

      if((ifbeam).or.(ifwake)) then
        do j=1,ncemface
           i= cemface(j)
           fjz(j,1,1) = jz(i,1,1,1)
        enddo
      endif

c     call exitt

      return
      end

c-----------------------------------------------------------------------
      subroutine add_flux(fullfld, facefld, faceskip)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real    fullfld(lx1,ly1,lz1,lelt)
      real    facefld(*)
      integer faceskip

      k = nx1*nz1*2*ndim*nelt
      do j=1,ncemface
         i = cemface(j)
         a = area(j,1,1,1)
         fullfld(i,1,1,1) = fullfld(i,1,1,1) + a*facefld(faceskip*k+j)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_add_flux_to_res(srflx)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real srflx(1)

      k= nx1*nz1*2*ndim*nelt

      if (if3d) then
        do j=1,ncemface
           i = cemface(j)
           a = area(j,1,1,1)
           reshx(i,1,1,1) = reshx(i,1,1,1) + a*srflx(0*k+j)
           reshy(i,1,1,1) = reshy(i,1,1,1) + a*srflx(1*k+j)
           reshz(i,1,1,1) = reshz(i,1,1,1) + a*srflx(2*k+j)
           resex(i,1,1,1) = resex(i,1,1,1) + a*srflx(3*k+j)
           resey(i,1,1,1) = resey(i,1,1,1) + a*srflx(4*k+j)
           resez(i,1,1,1) = resez(i,1,1,1) + a*srflx(5*k+j)
        enddo
      else
        if     (iftm) then
        do j=1,ncemface
           i = cemface(j)
           a = area(j,1,1,1)
           reshx(i,1,1,1) = reshx(i,1,1,1) + a*srflx(0*k+j)
           reshy(i,1,1,1) = reshy(i,1,1,1) + a*srflx(1*k+j)
           resez(i,1,1,1) = resez(i,1,1,1) + a*srflx(2*k+j)
        enddo
        elseif (ifte) then
        do j=1,ncemface
           i = cemface(j)
           a = area(j,1,1,1)
           resex(i,1,1,1) = resex(i,1,1,1) + a*srflx(0*k+j)
           resey(i,1,1,1) = resey(i,1,1,1) + a*srflx(1*k+j)
           reshz(i,1,1,1) = reshz(i,1,1,1) + a*srflx(2*k+j)
        enddo
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_flux(srflx)

      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      integer  icalld
      save     icalld
      data     icalld/0/
      real     srflx(1)

c...FIXME
c     if (ifsftf) then
c
c         !if (icalld.eq.0) then
c         !if (nid.eq.0) write(6,*) 'cem_outer',istep,icalld
c         !   if (ifpml) call cem_outer_pmlsrf
c         !   icalld = 1
c         !endif
c
c         if (icalld.eq.0) then
c            if (ifpml) call cem_outer_pmlsrf
c            icalld = 1
c         endif
c
c         if (if3d) then
c            call cem_flux3d_sftf (srflx)                 
c         else
c            call cem_flux2d_sftf (srflx)                 
c         endif
c
c     else 
c...FIXME

         if (if3d) then
             call cem_flux3d (srflx)
         else
             call cem_flux2d (srflx)
         endif

c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_outer_pmlsrf               
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'PML'
c
      integer nxyzf,i,j,k,e,f,ef,nxz

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1

c.....basic setting for inner pml surface
      k = nxyzf
      do e=1,nelt
      do f=1,2*ndim
         ef=eface(f)
         if (if3d) then
            do i=1,nx1*nz1
               j=(e-1)*nx1*nz1*2*ndim+(f-1)*nxz+i 
               if (pmlsrf(e,f).eq.1) then
                  srftmp(0*k+j)=1.0     
                  srftmp(1*k+j)=1.0     
                  srftmp(2*k+j)=1.0     
                  srftmp(3*k+j)=1.0     
                  srftmp(4*k+j)=1.0     
                  srftmp(5*k+j)=1.0     
                 if (incflag.eq.1) then
                  if (f.eq.5) then ! FIXME: mesh restriction *n2to3 type 
                  srftmp(0*k+j)=1.0     
                  srftmp(1*k+j)=1.0     
                  srftmp(2*k+j)=1.0     
                  srftmp(3*k+j)=1.0     
                  srftmp(4*k+j)=1.0     
                  srftmp(5*k+j)=1.0     
                  endif
                 elseif (incflag.eq.-1) then
                  if (f.eq.6) then ! FIXME: mesh restriction *n2to3 type 
                  srftmp(0*k+j)=1.0     
                  srftmp(1*k+j)=1.0     
                  srftmp(2*k+j)=1.0     
                  srftmp(3*k+j)=1.0     
                  srftmp(4*k+j)=1.0     
                  srftmp(5*k+j)=1.0     
                  endif
                 endif
               else
                  srftmp(0*k+j)=1.0     
                  srftmp(1*k+j)=1.0     
                  srftmp(2*k+j)=1.0     
                  srftmp(3*k+j)=1.0     
                  srftmp(4*k+j)=1.0     
                  srftmp(5*k+j)=1.0     
               endif
            enddo
         else
            do i=1,nx1*nz1
               j=(e-1)*nx1*nz1*2*ndim+(f-1)*nxz+i 
               if (pmlsrf(e,f).eq.1) then
                  srftmp(0*k+j)=1.0     
                  srftmp(1*k+j)=1.0     
                  srftmp(2*k+j)=1.0     
                 if (incflag.eq.1) then
                  if (f.eq.1) then ! FIXME: mesh restriction *n1to2 type 
                  srfinc(0*k+j)=1.0     
                  srfinc(1*k+j)=1.0     
                  srfinc(2*k+j)=1.0     
                  endif
                 elseif (incflag.eq.-1) then
                  if (f.eq.3) then ! FIXME: mesh restriction *n1to2 type 
                  srfinc(0*k+j)=1.0     
                  srfinc(1*k+j)=1.0     
                  srfinc(2*k+j)=1.0     
                  endif
                 endif
               else
                  srftmp(0*k+j)=1.0     
                  srftmp(1*k+j)=1.0     
                  srftmp(2*k+j)=1.0     
               endif
            enddo
         endif
      enddo
      enddo


      if (if3d) then
      call gs_op_fields(gsh_face,srftmp,nxyzf,6,1,1,0)
      else
      call gs_op_fields(gsh_face,srftmp,nxyzf,3,1,1,0)
      endif

c.....adding 32 to the inner pml surface 
      k = nxyzf
      do e=1,nelt
      do f=1,2*ndim
         ef= eface(f)
         if (if3d) then
            do i=1,nx1*nz1
               j=(e-1)*nx1*nz1*2*ndim+(f-1)*nxz+i 
               if (pmlsrf(e,f).eq.1) then
               srftmp(0*k+j)=srftmp(0*k+j)+32.0
               srftmp(1*k+j)=srftmp(1*k+j)+32.0
               srftmp(2*k+j)=srftmp(2*k+j)+32.0
               srftmp(3*k+j)=srftmp(3*k+j)+32.0
               srftmp(4*k+j)=srftmp(4*k+j)+32.0
               srftmp(5*k+j)=srftmp(5*k+j)+32.0
               else
               srftmp(0*k+j)=srftmp(0*k+j)
               srftmp(1*k+j)=srftmp(1*k+j)
               srftmp(2*k+j)=srftmp(2*k+j)
               srftmp(3*k+j)=srftmp(3*k+j)
               srftmp(4*k+j)=srftmp(4*k+j)
               srftmp(5*k+j)=srftmp(5*k+j)
               endif
            enddo
         else
            do i=1,nx1*nz1
               j=(e-1)*nx1*nz1*2*ndim+(f-1)*nxz+i 
               if (pmlsrf(e,f).eq.1) then    ! inner pml    
               srftmp(0*k+j)=srftmp(0*k+j)+32.0
               srftmp(1*k+j)=srftmp(1*k+j)+32.0
               srftmp(2*k+j)=srftmp(2*k+j)+32.0
               else
               srftmp(0*k+j)=srftmp(0*k+j)
               srftmp(1*k+j)=srftmp(1*k+j)
               srftmp(2*k+j)=srftmp(2*k+j)
               endif
            enddo
         endif
      enddo
      enddo

      ichk1=0
      ichk2=0

      nsrf = nx1*nz1*2*ndim*nelt

      call rzero(srfout,6*nsrf)
      call rzero(srfin ,6*nsrf)

      k= nxyzf
      if (if3d) then

         do i=1,nx1*nz1*2*ndim*nelt

         if     (srftmp(0*k+i).eq. 1) then
                 srfout(0*k+i)=1 !srftmp(0*k+i)
                 srfout(1*k+i)=1 !srftmp(1*k+i)
                 srfout(2*k+i)=1 !srftmp(2*k+i)
                 srfout(3*k+i)=1 !srftmp(3*k+i)
                 srfout(4*k+i)=1 !srftmp(4*k+i)
                 srfout(5*k+i)=1 !srftmp(4*k+i)

                 ichk1= ichk1+1

         elseif (srftmp(0*k+i).eq.33) then

                 srfin (0*k+i)=1 !srftmp(0*k+i)
                 srfin (1*k+i)=1 !srftmp(1*k+i)
                 srfin (2*k+i)=1 !srftmp(2*k+i)
                 srfin (3*k+i)=1 !srftmp(3*k+i)
                 srfin (4*k+i)=1 !srftmp(4*k+i)
                 srfin (5*k+i)=1 !srftmp(4*k+i)
                 ichk2= ichk2+1

         endif  

         enddo

      else

         do i=1,nx1*nz1*2*ndim*nelt

         if     (srftmp(0*k+i).eq. 1) then

                 srfout(0*k+i)=1 !srftmp(0*k+i)
                 srfout(1*k+i)=1 !srftmp(1*k+i)
                 srfout(2*k+i)=1 !srftmp(2*k+i)
                 ichk1= ichk1+1

         elseif (srftmp(0*k+i).eq.33) then

                 srfin (0*k+i)=1 !srftmp(0*k+i)
                 srfin (1*k+i)=1 !srftmp(1*k+i)
                 srfin (2*k+i)=1 !srftmp(2*k+i)
                 ichk2= ichk2+1

         endif  

         enddo
      endif

      ichk1= iglsum(ichk1,1)
      ichk2= iglsum(ichk2,1)

      if (ichk1.eq.ichk2) then
c       if (nid.eq.0) write(6,*) 'pml face equal    ',ichk1,ichk2    
      else
        if (nid.eq.0) write(6,*) 'pml face not equal',ichk1,ichk2     
         call exitt
      endif

      k= nxyzf
      if (if3d) then
         l=0
         do i=1,nx1*nz1*2*ndim*nelt
            if  (srfinc(0*k+1).eq.1) then
                 j=cemface(i)
                 xminval=xm1(j,1,1,1)
                 yminval=ym1(j,1,1,1)
                 zminval=zm1(j,1,1,1)
                 xmaxval=xm1(j,1,1,1)
                 ymaxval=ym1(j,1,1,1)
                 zmaxval=zm1(j,1,1,1)
                 l=l+1
            endif
         enddo
      else
         l=0
         do i=1,nx1*nz1*2*ndim*nelt
            if  (srfinc(0*k+i).eq.1) then
                 j=cemface(i)
                 xminval=xm1(j,1,1,1)
                 yminval=ym1(j,1,1,1)
                 zminval=zm1(j,1,1,1)
                 xmaxval=xm1(j,1,1,1)
                 ymaxval=ym1(j,1,1,1)
                 zmaxval=zm1(j,1,1,1)
                 l=l+1
            endif
         enddo

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_flux2d_sftf (srflx)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'PML'
 
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer nxyzf,i,j,k,e,f,nxz

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0

      if     (IFTM) then

      k= nxyzf
      do i= 1,nxyzf 
         j= cemface(i)

         if (srfin(0*k+i).eq.1) then
         srflxi(0*k+i)=-uny(i,1,1,1)*(incEZ(j,1,1,1))           ! NxfHx
         srflxi(1*k+i)= unx(i,1,1,1)*(incEZ(j,1,1,1))           ! NxfHy
         srflxi(2*k+i)=-unx(i,1,1,1)*(incHY(j,1,1,1))                    
     $                 +uny(i,1,1,1)*(incHX(j,1,1,1))           ! NxfEz
         endif
      enddo
 
      do i= 1,nxyzf 
         j= cemface(i)
         if (srfout (0*k+i).eq.1) then
         srflxo(0*k+i)=-uny(i,1,1,1)*(incEZ(j,1,1,1))                ! NxfHx
         srflxo(1*k+i)= unx(i,1,1,1)*(incEZ(j,1,1,1))                ! NxfHy
         srflxo(2*k+i)=-unx(i,1,1,1)*(incHY(j,1,1,1))                    
     $                 +uny(i,1,1,1)*(incHX(j,1,1,1))                ! NxfEz
         endif
      enddo

      ll=0
      do i= 1,nxyzf
         j= cemface(i)

         if      (srfinc(0*k+i).eq.1) then
          srflx(0*k+i)=-uny(i,1,1,1)*(fEZ(i,1,1)-incez(j,1,1,1))! NxfHx
          srflx(1*k+i)= unx(i,1,1,1)*(fEZ(i,1,1)-incez(j,1,1,1))! NxfHy
          srflx(2*k+i)=-unx(i,1,1,1)*(fHY(i,1,1)-inchy(j,1,1,1))         
     $                 +uny(i,1,1,1)*(fHX(i,1,1)-inchx(j,1,1,1))! NxfEz
          ll=ll+1
         else
          srflx(0*k+i)=-uny(i,1,1,1)*(fEZ(i,1,1))               ! NxfHx
          srflx(1*k+i)= unx(i,1,1,1)*(fEZ(i,1,1))               ! NxfHy
          srflx(2*k+i)=-unx(i,1,1,1)*(fHY(i,1,1))                    
     $               +uny(i,1,1,1)*(fHX(i,1,1))               ! NxfEz
         endif
      enddo
 
      call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)

      if (IFPML.or.IFPEC) call cem_flux_pec(srflx)

      k= nxyzf
      do i=1,nxyzf

         Y0 = Y_0(i,1,1)
         Y1 = Y_1(i,1,1)
         Z0 = Z_0(i,1,1)
         Z1 = Z_1(i,1,1)

         fu1 = uny(i,1,1,1)*srflx(2*k+i)
         fu2 =-unx(i,1,1,1)*srflx(2*k+i)
         fw3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

         srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)
         srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)
         srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)

      enddo

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  1. check out the blocked srflx
c  2. check out the location of finchz
c     currently not at the surface location
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      elseif (IFTE) then

      k = nxyzf
      do i=1,nxyzf
         j=cemface(i)

         if      (srfinc(0*k+i).eq.1) then
       
         srflx(0*k+i)=-uny(i,1,1,1)*(fHZ(i,1,1)-inchx(j,1,1,1))  ! NxfHx   
         srflx(1*k+i)= unx(i,1,1,1)*(fHZ(i,1,1)-inchz(j,1,1,1))  ! NxfHy   
         srflx(2*k+i)=-unx(i,1,1,1)*(fEY(i,1,1)-incey(j,1,1,1))
     $                +uny(i,1,1,1)*(fEX(i,1,1)-incex(j,1,1,1)) !  NxfEz
         else 
         srflx(0*k+i)=-uny(i,1,1,1)*fHZ(i,1,1)  ! NxfHx   
         srflx(1*k+i)= unx(i,1,1,1)*fHZ(i,1,1)  ! NxfHy   
         srflx(2*k+i)=-unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1) !  NxfEz
         endif

      enddo

      call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)

      if (IFPML.or.IFPEC) call cem_flux_pec(srflx)                

      do i=1,nxyzf

         Y0 = Y_0(i,1,1)
         Y1 = Y_1(i,1,1)
         Z0 = Z_0(i,1,1)
         Z1 = Z_1(i,1,1)

         fw1 = uny(i,1,1,1)*srflx(2*k+i)
         fw2 =-unx(i,1,1,1)*srflx(2*k+i)
         fu3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

         srflx(0*k+i) = 0.5/Z0*( Z1*srflx(0*k+i) - C0*fw1)
         srflx(1*k+i) = 0.5/Z0*( Z1*srflx(1*k+i) - C0*fw2)
         srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3)

      enddo

      else
         write(6,*) 'ERROR:: IFTM or IFTE'
         call exitt
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_flux2d_sftf1 (srflx)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'PML'
c
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer nxyzf,i,j,k,e,f,nxz

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0
      
      if     (IFTM) then

c...  to be substracted inside: srfin
      k = nxyzf
      do i=1,nxyzf 
         j=cemface(i)
         if (srfin(0*k+i).eq.1) then
         srflxi(0*k+i)=-uny(i,1,1,1)*(incEZ(j,1,1,1))                ! NxfHx
         srflxi(1*k+i)= unx(i,1,1,1)*(incEZ(j,1,1,1))                ! NxfHy
         srflxi(2*k+i)=-unx(i,1,1,1)*(incHY(j,1,1,1))                    
     $                 +uny(i,1,1,1)*(incHX(j,1,1,1))                ! NxfEz
         endif
      enddo

c...  to be added outside: srfout  
      do i=1,nxyzf 
         j=cemface(i)
         if (srfout (0*k+i).eq.1) then
         srflxo(0*k+i)=-uny(i,1,1,1)*(incEZ(j,1,1,1))                ! NxfHx
         srflxo(1*k+i)= unx(i,1,1,1)*(incEZ(j,1,1,1))                ! NxfHy
         srflxo(2*k+i)=-unx(i,1,1,1)*(incHY(j,1,1,1))                    
     $                 +uny(i,1,1,1)*(incHX(j,1,1,1))                ! NxfEz
c        write(6,*) 'ttt',i,unx(i,1,1,1),srflxo(k+1)
         endif
      enddo
      
      kk= 0
      k = nxyzf
      do i=1,nxyzf
         j=cemface(i)
c ...    substract incident field inside
         if     ((srfin(0*k+i).eq.1).and.(srfout(0*k+i).eq.0)) then
          srflx(0*k+i)=-uny(i,1,1,1)*(fEZ(i,1,1)-incez(j,1,1,1)) ! NxfHx
          srflx(1*k+i)= unx(i,1,1,1)*(fEZ(i,1,1)-incez(j,1,1,1)) ! NxfHy
          srflx(2*k+i)=-unx(i,1,1,1)*(fHY(i,1,1)-inchy(j,1,1,1))          
     $                 +uny(i,1,1,1)*(fHX(i,1,1)-inchx(j,1,1,1)) ! NxfEz
c ...    add incident field outside
         elseif ((srfin(0*k+i).eq.0).and.(srfout(0*k+i).eq.1)) then
          srflx(0*k+i)=-uny(i,1,1,1)*(fEZ(i,1,1)+incez(j,1,1,1)) ! NxfHx
          srflx(1*k+i)= unx(i,1,1,1)*(fEZ(i,1,1)+incez(j,1,1,1)) ! NxfHy
          srflx(2*k+i)=-unx(i,1,1,1)*(fHY(i,1,1)+inchy(j,1,1,1))
     $                 +uny(i,1,1,1)*(fHX(i,1,1)+inchx(j,1,1,1)) ! NxfEz
c ...    regular                      
         elseif ((srfin(0*k+i).eq.0).and.(srfout(0*k+i).eq.0)) then
          srflx(0*k+i)=-uny(i,1,1,1)*fEZ(i,1,1)!-incez(j,1,1,1)) ! NxfHx
          srflx(1*k+i)= unx(i,1,1,1)*fEZ(i,1,1)!-incez(j,1,1,1)) ! NxfHy
          srflx(2*k+i)=-unx(i,1,1,1)*fHY(i,1,1)
     $                 +uny(i,1,1,1)*fHX(i,1,1)!-inchx(j,1,1,1)) ! NxfEz
         else
            call exitt
         endif
      enddo
 
      call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)

      k = nxyzf
      do i=1,nxyzf
c        srflx(0*k+i)=srflx(0*k+i)+srflxi(k  )-srflxo(k  )             
c        srflx(k+1)=srflx(k+1)+srflxi(k+1)-srflxo(k+1)                
c        srflx(k+2)=srflx(k+2)+srflxi(k+2)-srflxo(k+2)                
c        write(6,22) '  ',i,srflxi(k+1),srflxo(k+1),srflx(k+1)
      enddo
  22  format(a2,i3,3f10.2)

      if (IFPML.or.IFPEC) call cem_flux_pec(srflx)

      k = nxyzf
      do i=1,nxyzf

         Y0 = Y_0(i,1,1)
         Y1 = Y_1(i,1,1)
         Z0 = Z_0(i,1,1)
         Z1 = Z_1(i,1,1)

         fu1= uny(i,1,1,1)*srflx(2*k+i)
         fu2=-unx(i,1,1,1)*srflx(2*k+i)
         fw3= unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

         srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)
         srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)
         srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)

      enddo

c     k =  1
c     do i=1,nxyzf
c        write(6,22) 'f-',i,srflxi(k+1),srflxo(k+1),srflx(k+1)
c        k = k+3
c     enddo

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c  1. check out the blocked srflx
c  2. check out the location of finchz
c     currently not at the surface location
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      elseif (IFTE) then

      k = nxyzf
      do i=1,nxyzf

         srflx(0*k+i)=-uny(i,1,1,1)*fHZ(i,1,1)  ! NxfHx   
         srflx(1*k+i)= unx(i,1,1,1)*fHZ(i,1,1)  ! NxfHy   
         srflx(2*k+i)=-unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1) !  NxfEz

      enddo

      call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)

      if (IFPML.or.IFPEC) call cem_flux_pec(srflx)                

      k = nxyzf
      do i=1,nxyzf

         Y0 = Y_0(i,1,1)
         Y1 = Y_1(i,1,1)
         Z0 = Z_0(i,1,1)
         Z1 = Z_1(i,1,1)

         fw1 = uny(i,1,1,1)*srflx(2*k+i)
         fw2 =-unx(i,1,1,1)*srflx(2*k+i)
         fu3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

         srflx(0*k+i) = 0.5/Z0*( Z1*srflx(0*k+i) - C0*fw1)
         srflx(1*k+i) = 0.5/Z0*( Z1*srflx(1*k+i) - C0*fw2)
         srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3)

      enddo

      else
         write(6,*) 'ERROR:: IFTM or IFTE'
         call exitt
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_flux2d (srflx)  

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
c
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer nxyzf,i,j,k,e,f,nxz
      real*8  dclock,t0

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0

      if     (IFTM) then

      k = nxyzf
      do i=1,nxyzf                
         srflx(0*k+i)=-uny(i,1,1,1)*fEZ(i,1,1)  ! NxfHx
         srflx(1*k+i)= unx(i,1,1,1)*fEZ(i,1,1)  ! NxfHy
         srflx(2*k+i)=-unx(i,1,1,1)*fHY(i,1,1)+uny(i,1,1,1)*fHX(i,1,1) ! NxfEz
      enddo

      t0=dclock()
      call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)
      if (mod(istep,iocomm).eq.0) then
      if (nid.eq.0)  write(6,2)  dclock()-t0
      endif

      if (IFPML.or.IFPEC) call cem_flux_pec(srflx)                   

      k = nxyzf
      do i=1,nxyzf

         Y0 = Y_0(i,1,1)             
         Y1 = Y_1(i,1,1)  
         Z0 = Z_0(i,1,1)   
         Z1 = Z_1(i,1,1)    

         fu1 = uny(i,1,1,1)*srflx(2*k+i)                                
         fu2 =-unx(i,1,1,1)*srflx(2*k+i)                     
         fw3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)     

         srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)    
         srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)   
         srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)   
      enddo

      elseif (IFTE) then

      k = nxyzf
      do i=1,nxyzf                

         srflx(0*k+i)= -uny(i,1,1,1)*fHZ(i,1,1)  ! NxfHx   
         srflx(1*k+i)=  unx(i,1,1,1)*fHZ(i,1,1)  ! NxfHy   
         srflx(2*k+i)= -unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1) ! NxfEz

      enddo

      t0=dclock()
      call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)
      if (mod(istep,iocomm).eq.0) then
      if (nid.eq.0)  write(6,2)  dclock()-t0
      endif
    2 format(' gs_op_fields time::',1pe11.4,' seconds ')

      if (IFPML.or.IFPEC) call cem_flux_pec(srflx)               

      k = nxyzf
      do i=1,nxyzf                

         Y0 = Y_0(i,1,1)
         Y1 = Y_1(i,1,1)
         Z0 = Z_0(i,1,1)
         Z1 = Z_1(i,1,1)

         fw1 = uny(i,1,1,1)*srflx(2*k+i)                            
         fw2 =-unx(i,1,1,1)*srflx(2*k+i)
         fu3 = unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

         srflx(0*k+i) = 0.5/Z0*( Z1*srflx(0*k+i) - C0*fw1)      
         srflx(1*k+i) = 0.5/Z0*( Z1*srflx(1*k+i) - C0*fw2)       
         srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3) 

      enddo

      else
         write(6,*) 'ERROR:: IFTM or IFTE'
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_flux3d_sftf(srflx)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'

c      n x dEx - c* n x n x dHx  
c
      real srflx(1)
      real const
      integer e, f 
      character CB*3 

      nxyzf =nx1*nz1*2*ndim*nelt ! face count
      nface =2*ndim
      nxz   =nx1*nz1
      
      if (IFCENTRAL) C0 = 0.0
      if (IFUPWIND ) C0 = 1.0

      k = nxyzf
      do i=1,nxyzf
         j=cemface(i)
       if (srfin (0*k+i).eq.1) then
c...     to be substracted inside
         srflxi(0*k+i)=-uny(i,1,1,1)*incEZ(j,1,1,1)
     $                 +unz(i,1,1,1)*incEY(j,1,1,1)
         srflxi(1*k+i)=-unz(i,1,1,1)*incEX(j,1,1,1)
     $                 +unx(i,1,1,1)*incEZ(j,1,1,1)
         srflxi(2*k+i)=-unx(i,1,1,1)*incEY(j,1,1,1)
     $                 +uny(i,1,1,1)*incEX(j,1,1,1)
         srflxi(3*k+i)=-uny(i,1,1,1)*incHZ(j,1,1,1)
     $                 +unz(i,1,1,1)*incHY(j,1,1,1)
         srflxi(4*k+i)=-unz(i,1,1,1)*incHX(j,1,1,1)
     $                 +unx(i,1,1,1)*incHZ(j,1,1,1)
         srflxi(5*k+i)=-unx(i,1,1,1)*incHY(j,1,1,1)
     $                 +uny(i,1,1,1)*incHX(j,1,1,1)
         else
         srflxi(0*k+i)= 0
         srflxi(1*k+i)= 0
         srflxi(2*k+i)= 0
         srflxi(3*k+i)= 0
         srflxi(4*k+i)= 0
         srflxi(5*k+i)= 0
         endif
      enddo

      do i=1,nxyzf
         j=cemface(i)
       if (srfout (0*k+i).eq.1) then
c...     to be substracted inside
         srflxo(0*k+i)=-uny(i,1,1,1)*incEZ(j,1,1,1)
     $                 +unz(i,1,1,1)*incEY(j,1,1,1)
         srflxo(1*k+i)=-unz(i,1,1,1)*incEX(j,1,1,1)
     $                 +unx(i,1,1,1)*incEZ(j,1,1,1)
         srflxo(2*k+i)=-unx(i,1,1,1)*incEY(j,1,1,1)
     $                 +uny(i,1,1,1)*incEX(j,1,1,1)
         srflxo(3*k+i)=-uny(i,1,1,1)*incHZ(j,1,1,1)
     $                 +unz(i,1,1,1)*incHY(j,1,1,1)
         srflxo(4*k+i)=-unz(i,1,1,1)*incHX(j,1,1,1)
     $                 +unx(i,1,1,1)*incHZ(j,1,1,1)
         srflxo(5*k+i)=-unx(i,1,1,1)*incHY(j,1,1,1)
     $                 +uny(i,1,1,1)*incHX(j,1,1,1)
         else
         srflxo(0*k+i)= 0
         srflxo(1*k+i)= 0
         srflxo(2*k+i)= 0
         srflxo(3*k+i)= 0
         srflxo(4*k+i)= 0
         srflxo(5*k+i)= 0
         endif
      enddo

      k = nxyzf
      do i=1,nxyzf
         j=cemface(i)
         if     (srfinc(0*k+i).eq.1) then

         srflx(0*k+i)=-uny(i,1,1,1)*(fEZ(i,1,1)-incez(j,1,1,1))
     $                +unz(i,1,1,1)*(fEY(i,1,1)-incey(j,1,1,1))
         srflx(1*k+i)=-unz(i,1,1,1)*(fEX(i,1,1)-incex(j,1,1,1))
     $                +unx(i,1,1,1)*(fEZ(i,1,1)-incez(j,1,1,1))
         srflx(2*k+i)=-unx(i,1,1,1)*(fEY(i,1,1)-incey(j,1,1,1))
     $                +uny(i,1,1,1)*(fEX(i,1,1)-incex(j,1,1,1))
         srflx(3*k+i)=-uny(i,1,1,1)*(fHZ(i,1,1)-inchz(j,1,1,1))
     $                +unz(i,1,1,1)*(fHY(i,1,1)-inchy(j,1,1,1))
         srflx(4*k+i)=-unz(i,1,1,1)*(fHX(i,1,1)-inchx(j,1,1,1))
     $                +unx(i,1,1,1)*(fHZ(i,1,1)-inchz(j,1,1,1))
         srflx(5*k+i)=-unx(i,1,1,1)*(fHY(i,1,1)-inchy(j,1,1,1))
     $                +uny(i,1,1,1)*(fHX(i,1,1)-inchx(j,1,1,1))

         else   

         srflx(0*k+i)=-uny(i,1,1,1)*fEZ(i,1,1)+unz(i,1,1,1)*fEY(i,1,1)
         srflx(1*k+i)=-unz(i,1,1,1)*fEX(i,1,1)+unx(i,1,1,1)*fEZ(i,1,1)
         srflx(2*k+i)=-unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1)
         srflx(3*k+i)=-uny(i,1,1,1)*fHZ(i,1,1)+unz(i,1,1,1)*fHY(i,1,1)
         srflx(4*k+i)=-unz(i,1,1,1)*fHX(i,1,1)+unx(i,1,1,1)*fHZ(i,1,1)
         srflx(5*k+i)=-unx(i,1,1,1)*fHY(i,1,1)+uny(i,1,1,1)*fHX(i,1,1)

         endif

      enddo

      call gs_op_fields(gsh_face,srflx,nxyzf,6,1,1,0)

      ! srflx now contains n x [E] and n x [H],
      ! with [F] = F^+ - F^-
      ! (observe un^+ = - un^-)

      if (IFPEC.or.IFPML)  call cem_flux_pec(srflx)                    

      k = nxyzf
      do i=1,nxyzf

        Y0 = Y_0(i,1,1)
        Y1 = Y_1(i,1,1)
        Z0 = Z_0(i,1,1)
        Z1 = Z_1(i,1,1)

        ! fu = n x n x [H]
        fu1 =uny(i,1,1,1)*srflx(5*k+i)-unz(i,1,1,1)*srflx(4*k+i)
        fu2 =unz(i,1,1,1)*srflx(3*k+i)-unx(i,1,1,1)*srflx(5*k+i)
        fu3 =unx(i,1,1,1)*srflx(4*k+i)-uny(i,1,1,1)*srflx(3*k+i)

        ! fw = n x n x [E]
        fw1 =uny(i,1,1,1)*srflx(2*k+i)-unz(i,1,1,1)*srflx(1*k+i)
        fw2 =unz(i,1,1,1)*srflx(0*k+i)-unx(i,1,1,1)*srflx(2*k+i)
        fw3 =unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

        srflx(0*k+i) = 0.5/Y0*(-Y1*srflx(0*k+i) - C0*fu1)                  
        srflx(1*k+i) = 0.5/Y0*(-Y1*srflx(1*k+i) - C0*fu2)
        srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3)         
        srflx(3*k+i) = 0.5/Z0*( Z1*srflx(3*k+i) - C0*fw1)         
        srflx(4*k+i) = 0.5/Z0*( Z1*srflx(4*k+i) - C0*fw2)        
        srflx(5*k+i) = 0.5/Z0*( Z1*srflx(5*k+i) - C0*fw3)        

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_flux3d(srflx)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

c      n x dEx - c* n x n x dHx  
c
      real srflx(1)
      real const
      integer e, f 
      character CB*3 
      real*8 dclock,t0

      nxyzf =nx1*nz1*2*ndim*nelt ! face count
      nface =2*ndim
      nxz   =nx1*nz1
      
      if (IFCENTRAL) C0 = 0.0
      if (IFUPWIND ) C0 = 1.0

      !  srflx(1)-> reshx
      !  srflx(4)-> resex
      k = nxyzf
      do i=1,nxyzf
         ! - n x E
         srflx(0*k+i)=-uny(i,1,1,1)*fEZ(i,1,1)+unz(i,1,1,1)*fEY(i,1,1)
         srflx(1*k+i)=-unz(i,1,1,1)*fEX(i,1,1)+unx(i,1,1,1)*fEZ(i,1,1)
         srflx(2*k+i)=-unx(i,1,1,1)*fEY(i,1,1)+uny(i,1,1,1)*fEX(i,1,1)
         ! - n x H
         srflx(3*k+i)=-uny(i,1,1,1)*fHZ(i,1,1)+unz(i,1,1,1)*fHY(i,1,1)
         srflx(4*k+i)=-unz(i,1,1,1)*fHX(i,1,1)+unx(i,1,1,1)*fHZ(i,1,1)
         srflx(5*k+i)=-unx(i,1,1,1)*fHY(i,1,1)+uny(i,1,1,1)*fHX(i,1,1)
      enddo

      t0=dclock()
      call gs_op_fields(gsh_face,srflx,nxyzf,6,1,1,0)
      if (mod(istep,iocomm).eq.0) then
      if (nid.eq.0)  write(6,3)  dclock()-t0
      endif
    3 format(' gs_op_fields::',1pe11.4,' seconds ')

      ! srflx now contains n x [E] and n x [H],
      ! with [F] = F^+ - F^-
      ! (observe un^+ = - un^-)

      if (IFPEC.or.IFPML)  call cem_flux_pec(srflx)                 

      k = nxyzf
      do i=1,nxyzf

        Y0 = Y_0(i,1,1)
        Y1 = Y_1(i,1,1)
        Z0 = Z_0(i,1,1)
        Z1 = Z_1(i,1,1)

        ! fu = n x n x [H]
        fu1 =uny(i,1,1,1)*srflx(5*k+i)-unz(i,1,1,1)*srflx(4*k+i)
        fu2 =unz(i,1,1,1)*srflx(3*k+i)-unx(i,1,1,1)*srflx(5*k+i)
        fu3 =unx(i,1,1,1)*srflx(4*k+i)-uny(i,1,1,1)*srflx(3*k+i)

        ! fw = n x n x [E]
        fw1 =uny(i,1,1,1)*srflx(2*k+i)-unz(i,1,1,1)*srflx(1*k+i)
        fw2 =unz(i,1,1,1)*srflx(0*k+i)-unx(i,1,1,1)*srflx(2*k+i)
        fw3 =unx(i,1,1,1)*srflx(1*k+i)-uny(i,1,1,1)*srflx(0*k+i)

        srflx(0*k+i) = 0.5/Y0*(-Y1*srflx(0*k+i) - C0*fu1)              
        srflx(1*k+i) = 0.5/Y0*(-Y1*srflx(1*k+i) - C0*fu2)
        srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3)         
        srflx(3*k+i) = 0.5/Z0*( Z1*srflx(3*k+i) - C0*fw1)         
        srflx(4*k+i) = 0.5/Z0*( Z1*srflx(4*k+i) - C0*fw2)        
        srflx(5*k+i) = 0.5/Z0*( Z1*srflx(5*k+i) - C0*fw3)        

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine assign_bdry_output                            

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
 
      common /scrnsbd/ srflx_bd0(6*2*ldim*lx1*lz1*lelt)
     $               , srflx_bd1(6*2*ldim*lx1*lz1*lelt)
     $               , srflx_bd2(6*2*ldim*lx1*lz1*lelt)
     $               , srflx_bd3(6*2*ldim*lx1*lz1*lelt)

      integer  nxyzf,nface,nxz,k,i,e,f,ef          
      character CB*3

      nxyzf= nx1*nz1*2*ndim*nelt
      nface= 2*ndim
      nxz  = nx1*nz1

      k0 = 0
      k1 = 0
      k2 = 0
      k3 = 0

      do e= 1,nelt            
      do f= 1,nface
      do j= 1,nxz

       CB = CBC(f,e,2)

       if     ((CB.eq.'PEC').or.(CB.eq.'pec')) then

         if(hxactive) srflx_bd1(k1+fldidx(1))=10                      
         if(hyactive) srflx_bd1(k1+fldidx(2))=10
         if(hzactive) srflx_bd1(k1+fldidx(3))=10                      
         if(exactive) srflx_bd1(k1+fldidx(4))=10
         if(eyactive) srflx_bd1(k1+fldidx(5))=10
         if(ezactive) srflx_bd1(k1+fldidx(6))=10
         !write(6,*) 'pec:: ',CB

       elseif ((CB.eq.'PML').or.(CB.eq.'pml')) then

         if(hxactive) srflx_bd2(k2+fldidx(1))=20                      
         if(hyactive) srflx_bd2(k2+fldidx(2))=20
         if(hzactive) srflx_bd2(k2+fldidx(3))=20                      
         if(exactive) srflx_bd2(k2+fldidx(4))=20
         if(eyactive) srflx_bd2(k2+fldidx(5))=20
         if(ezactive) srflx_bd2(k2+fldidx(6))=20
         !write(6,*) 'pml:: ',CB

       elseif ((CB.eq.'P  ').or.(CB.eq.'p  ')) then

         if(hxactive) srflx_bd3(k3+fldidx(1))=30                      
         if(hyactive) srflx_bd3(k3+fldidx(2))=30
         if(hzactive) srflx_bd3(k3+fldidx(3))=30                      
         if(exactive) srflx_bd3(k3+fldidx(4))=30
         if(eyactive) srflx_bd3(k3+fldidx(5))=30
         if(ezactive) srflx_bd3(k3+fldidx(6))=30
         !write(6,*) 'p  :: ',CB

       elseif (CB.eq.'E  ') then

         if(hxactive) srflx_bd0(k0+fldidx(1))=1                     
         if(hyactive) srflx_bd0(k0+fldidx(2))=1 
         if(hzactive) srflx_bd0(k0+fldidx(3))=1                     
         if(exactive) srflx_bd0(k0+fldidx(4))=1 
         if(eyactive) srflx_bd0(k0+fldidx(5))=1 
         if(ezactive) srflx_bd0(k0+fldidx(6))=1 
         !write(6,*) 'e  :: ',CB

       endif

       k0=k0+fldcnt
       k1=k1+fldcnt
       k2=k2+fldcnt
       k3=k3+fldcnt

      enddo
      enddo
      enddo

      if (icalld.eq.0) then
        call cem_set_fc_ptr
        icalld = 1
      endif

      n = nx1*ny1*nz1*nelt

      call rzero(FACEB0,n)
      call rzero(FACEB1,n)
      call rzero(FACEB2,n)
      call rzero(FACEB3,n)

      k0 = 1
      k1 = 1
      k2 = 1
      k3 = 1

      do j=1,ncemface

         i = cemface(j)

         FACEB0(i,1,1,1) = srflx_bd0(k0)                        
         FACEB1(i,1,1,1) = srflx_bd1(k1)                        
         FACEB2(i,1,1,1) = srflx_bd2(k2)                        
         FACEB3(i,1,1,1) = srflx_bd3(k3)                        

         k0 = k0 + fldcnt   
         k1 = k1 + fldcnt   
         k2 = k2 + fldcnt   
         k3 = k3 + fldcnt   

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_flux_pec(srflx)                              

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
 
      real     srflx(1)
      integer  nxyzf,nface,nxz,k,i,e,f,ef          
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      common  /BDRY1/ cempec(6*lx1*lz1*2*ldim*lelt)
      common  /BDRY2/ ncempec
      integer  cempec,ncempec

      nxyzf= nx1*nz1*2*ndim*nelt
      nface= 2*ndim
      nxz  = nx1*nz1
      k    = nxyzf
      l    = 0

      if (icalld.eq.0) then

         do e= 1,nelt            
         do f= 1,nface
         do j= 1,nxz
            CB= CBC(f,e,2)
            i =(e-1)*nface*nxz+(f-1)*nxz+j
            if ((CB.eq.'PEC').or.(CB.eq.'PML')) then
              l= l+1
              cempec(l)= i
              ncempec  = l
            endif
         enddo
         enddo
         enddo
         icalld =1

      endif

      if (ncempec.eq.0) return
 
      if (if3d) then

         do j= 1,ncempec         
            i= cempec(j)

         ! here, srflx is (n x -E, n x -H)
         ! (we're on the boundary and fgs_gop left us untouched)
         ! since we want mirror conditions, i.e.
         ! n x E+ = -n x E-
         ! n x H+ =  n x H-,
         ! (cf. Hesthaven, Warburton, HO-/SM on Unst. Grids, 
         ! I. Time-Domain Solution of Maxwell's Eqns., p. 14)
         ! we have
         ! n x [E] = n x (E+ - E-) = n x -2E-
         ! n x [H] = n x (H+ - H-) = 0

         ! write(6,*) 'ii--',ncempec,l,cempec(l) 
          srflx(0*k+i)=2.0*srflx(0*k+i)
          srflx(1*k+i)=2.0*srflx(1*k+i)
          srflx(2*k+i)=2.0*srflx(2*k+i)
          srflx(3*k+i)=0
          srflx(4*k+i)=0
          srflx(5*k+i)=0
         enddo

      else
       if     (iftm) then

         do j= 1,ncempec         
            i= cempec(j)
            srflx(0*k+i)=2.0*srflx(0*k+i)
            srflx(1*k+i)=2.0*srflx(1*k+i)
            srflx(2*k+i)=0
         enddo

       elseif (ifte) then

         do j= 1,ncempec         
            i= cempec(j)
            srflx(0*k+i)=0
            srflx(1*k+i)=0
            srflx(2*k+i)=2.0*srflx(2*k+i)
         enddo

       endif
      endif
    
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_dflux_2d(srflx)
c     ffff
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'WZ'
      include 'DEALIAS'
      include 'PARALLEL'
c
      real srflx(1)

      parameter (ldd =lxd*4*lelt )
      parameter (md  =lxd        )

      common /dfaces/ fmhx(ldd),fmhy(ldd),fmhz(ldd)
     $              , fmex(ldd),fmey(ldd),fmez(ldd) 
      common /dscrn1/ srflxm (6*2*ldim*lxd*lelt)
      common /dscrn2/ srflxm0(6*2*ldim*lxd*lelt)
     $              , srflxm1(6*2*ldim*lxd*lelt)
     $              , srflxm2(6*2*ldim*lxd*lelt)
      common /weight/ zzd(lxd),wgld(lxd)
      common /invrho/ rho2i(lx1*lz1)
      real srflxm,srflxm0, srflxm1, srflxm2
      real srflx0,srflx1, srflx2                   

      integer e,f,k,i,jd,j1,nn,mm,nxyd,nxy,mx            
      integer icalld
      save    icalld
      data    icalld / 0 /

      mx    = nx1
      nxy   = nx1*nz1
      nxyd  = nxd*nzd
      nface = 2*ndim

      if (icalld.eq.0) then

         icalld = 1
         call zwgll(zzd,wgld,md)

         l = 0

         do k=1,nz1
         do i=1,nx1
            l = l+1
            rho2i(l) = 1./(wxm1(i)*wxm1(k))
         enddo
         enddo

         do e=1,nelt   ! this will need to be dynamic if we have p-refinement
         do f=1,nface

            call col3(ta,area(1,1,f,e),rho2i,nxy) ! surface Jacobian only
            call col3(wk,ta,unx(1,1,f,e),nxy)
            call intp_rsd (aud(1,1,f,e),wk,mx,md,if3d,0)  ! 1 = unx

            call col3(wk,ta,uny(1,1,f,e),nxy)
            call intp_rsd (aud(1,2,f,e),wk,mx,md,if3d,0)  ! 2 = uny

            call col3(wk,ta,unz(1,1,f,e),nxy)
            call intp_rsd (aud(1,3,f,e),wk,mx,md,if3d,0)  ! 3 = unz

            l = 0

            do k=1,nxd
            do i=1,nxd

               l   = l+1
               tmp = wgld(i)*wgld(k)
               aud(l,1,f,e) = tmp*aud(l,1,f,e)
               aud(l,2,f,e) = tmp*aud(l,2,f,e)
               aud(l,3,f,e) = tmp*aud(l,3,f,e)

            enddo
            enddo
               
         enddo
         enddo

      endif

      do e=1,nelv   
      do f=1,nface  
         call intp_rsd (fmez,fez(1,f,e),mx,md,if3d,0)                  
         call intp_rsd (fmhx,fhx(1,f,e),mx,md,if3d,0)                  
         call intp_rsd (fmhy,fhy(1,f,e),mx,md,if3d,0)                   
      enddo
      enddo

      k = 1

      myz= md*2*ndim*nelt

      do i=1,myz             

         srflxm(k  ) = aud(i,2,1,1)*fmez(i)                      ! fhx
         srflxm(k+1) =-aud(i,1,1,1)*fmez(i)
         srflxm(k+2) = aud(i,2,1,1)*fmhx(i)-aud(i,1,1,1)*fmhy(i) ! fex

         srflxm(k  ) = 0.5*srflxm(k  )                           
         srflxm(k+1) = 0.5*srflxm(k+1)                           
         srflxm(k+2) = 0.5*srflxm(k+2)                           

         k = k+3

      enddo

c     call fgs_gop_vec (gsh_face, srflxm ,'+\0',3) ! vector gs

      do i=1,myz             
         srflxm0(i)=srflxm(3*(i-1)+1)
         srflxm1(i)=srflxm(3*(i-1)+2)
         srflxm2(i)=srflxm(3*(i-1)+3)
      enddo

      do e=1,nelv   
      do f=1,nface  
         call intp_rsd (srflx0,srflxm0,mx,md,if3d,1)                    
         call intp_rsd (srflx1,srflxm1,mx,md,if3d,1)                   
         call intp_rsd (srflx2,srflxm2,mx,md,if3d,1)                   
      enddo
      enddo

      nyze = mx*2*ndim*nelt

      do i=1,nyze             
         srflx(3*(i-1)+1)=srflxm0(i)
         srflx(3*(i-1)+2)=srflxm1(i)
         srflx(3*(i-1)+3)=srflxm2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine wght_dcurl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2) 

       include 'SIZE' 
       include 'TOTAL' 
       include 'EMWAVE' 
       include 'DEALIAS' 

c     ddddd
c
c     Combined weighted curl:   w12 = B*curl(u3), w3=B*curl(u1,u2)
c     with dealiasing
c
                        
      real     u1(1),u2(1),u3(1)
      real     w1(1),w2(1),w3(1)
      real     work1(1),work2(1)
      logical  ifavg
                                    
      parameter (ldd =lxd**ldim  )
      parameter (md  =lxd        )

      common /t1t/ ud1r(ldd),ud1s(ldd),ud1t(ldd)
     $           , ud2r(ldd),ud2s(ldd),ud2t(ldd)
     $           , ud3r(ldd),ud3s(ldd),ud3t(ldd)
 
      real           ud1r,ud1s,ud1t
     $           ,   ud2r,ud2s,ud2t
     $           ,   ud3r,ud3s,ud3t

      common /ddtmp2/ ju1(ldd),ju2(ldd),ju3(ldd)
      real            ju1,ju2,ju3

      common /ddtmp3/ wd1(ldd),wd2(ldd),wd3(ldd)
      real            wd1,wd2,wd3
      common /ddtmp4/ dxmd(lxd,lxd),dxtmd(lxd,lxd)                     
      common /ddtmp5/ w3md(lxd**ldim),wgld(lxd),zmd(lxd)

c     common /ddtmp5/ w3md(lxd**ldim),wgld(3*(lxd+lx1))

      real            dxmd,dxtmd,w3md,wgld,zmd

      integer icalld, e, l, i, j, k, j1, jd, nxyd, nxy,nn,mm
      save    icalld 
      data    icalld /0/
 
      mx    = nx1
      nn    = nx1-1   
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

                wd1(i) = ( ud3r(i)*rymd(k)
     $                   + ud3s(i)*symd(k)
     $                   + ud3t(i)*tymd(k)
     $                   - ud2r(i)*rzmd(k)
     $                   - ud2s(i)*szmd(k)
     $                   - ud2t(i)*tzmd(k) ) * w3md(i)

                wd2(i) = ( ud1r(i)*rzmd(k)
     $                   + ud1s(i)*szmd(k)
     $                   + ud1t(i)*tzmd(k) 
     $                   - ud3r(i)*rxmd(k)
     $                   - ud3s(i)*sxmd(k) 
     $                   - ud3t(i)*txmd(k) ) * w3md(i)

                wd3(i) = ( ud2r(i)*rxmd(k)
     $                   + ud2s(i)*sxmd(k) 
     $                   + ud2t(i)*txmd(k) 
     $                   - ud1r(i)*rymd(k) 
     $                   - ud1s(i)*symd(k) 
     $                   - ud1t(i)*tymd(k) ) * w3md(i)
             enddo 

         else

          call local_grad2 (ud1r,ud1s,ju1,mm,1,dxmd,dxtmd)
          call local_grad2 (ud2r,ud2s,ju2,mm,1,dxmd,dxtmd)
          call local_grad2 (ud3r,ud3s,ju3,mm,1,dxmd,dxtmd)

          do i= 1,nxyd

                k = i+j              

                wd1(i) = ( ud3r(i)*rymd(k)
     $                   + ud3s(i)*symd(k) ) * w3md(i)
 
                wd2(i) =-( ud3r(i)*rxmd(k)
     $                   + ud3s(i)*sxmd(k) ) * w3md(i)
 
                wd3(i) = ( ud2r(i)*rxmd(k)
     $                   + ud2s(i)*sxmd(k) 
     $                   - ud1r(i)*rymd(k) 
     $                   - ud1s(i)*symd(k) ) * w3md(i)

          enddo 

         endif

         call intp_rstd (w1(j1),wd1,mx,md,if3d,1) ! 1 --> backward
         call intp_rstd (w2(j1),wd2,mx,md,if3d,1) ! 1 --> backward
         call intp_rstd (w3(j1),wd3,mx,md,if3d,1) ! 1 --> backward

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_grad (w1,w2,w3,u1,c1)
c
c     Weighted curl:   w = -grad(u1)
c
      include 'SIZE'
      include 'TOTAL'
c
      real w1(1),w2(1),w3(1),u1(1)                 
      integer e, k, i, j, nxyz,ntot,nn
      real c1

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)
         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) =c1*( u1r(i)*rxm1(k,1,1,1)
     $                 +u1s(i)*sxm1(k,1,1,1)
     $                 +u1t(i)*txm1(k,1,1,1) )/jacm1(i,1,1,e)
            w2(k) =c1*( u1r(i)*rym1(k,1,1,1)
     $                 +u1s(i)*sym1(k,1,1,1)
     $                 +u1t(i)*tym1(k,1,1,1) )/jacm1(i,1,1,e)
            w3(k) =c1*( u1r(i)*rzm1(k,1,1,1)
     $                 +u1s(i)*szm1(k,1,1,1)
     $                 +u1t(i)*tzm1(k,1,1,1) )/jacm1(i,1,1,e)
         enddo
         enddo
       else

         do e=1,nelt

         j = nxyz*(e-1)
         call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) =c1*(  u1r(i)*rxm1(k,1,1,1)
     $                 + u1s(i)*sxm1(k,1,1,1) )/jacm1(i,1,1,e)
            w2(k) =c1*(  u1r(i)*rym1(k,1,1,1)
     $                 + u1s(i)*sym1(k,1,1,1) )/jacm1(i,1,1,e)
c           write(6,*) 'tt',e,i,rxm1(k,1,1,1)/jacm1(i,1,1,e)
c    $                         ,rym1(k,1,1,1)/jacm1(i,1,1,e)
c           write(6,*) 'tt',e,i,u1r(i),u1s(i)                          
c           write(6,*) 'tt',e,i,xm1(k,1,1,1),ym1(k,1,1,1)             

         enddo
         enddo

       endif

       return
       end
c-----------------------------------------------------------------------
      subroutine curl_op(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)

c
c     Weighted curl:   w = B*curl(u)
c
      include 'SIZE'
      include 'TOTAL'
c
      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
      logical ifavg
      integer e, k, i, j, nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
      common /curlw/ sw1(lx1*ly1*lz1*lelt)
     $             , sw2(lx1*ly1*lz1*lelt)
     $             , sw3(lx1*ly1*lz1*lelt)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            w1(k) = (    u3r(i)*rym1(k,1,1,1)
     $                 + u3s(i)*sym1(k,1,1,1)
     $                 + u3t(i)*tym1(k,1,1,1)
     $                 - u2r(i)*rzm1(k,1,1,1)
     $                 - u2s(i)*szm1(k,1,1,1)
     $                 - u2t(i)*tzm1(k,1,1,1) ) /jacm1(i,1,1,e)
c
            w2(k) = (    u1r(i)*rzm1(k,1,1,1)
     $                 + u1s(i)*szm1(k,1,1,1)
     $                 + u1t(i)*tzm1(k,1,1,1) 
     $                 - u3r(i)*rxm1(k,1,1,1)
     $                 - u3s(i)*sxm1(k,1,1,1) 
     $                 - u3t(i)*txm1(k,1,1,1) ) /jacm1(i,1,1,e)
c
            w3(k) = (    u2r(i)*rxm1(k,1,1,1)
     $                 + u2s(i)*sxm1(k,1,1,1) 
     $                 + u2t(i)*txm1(k,1,1,1) 
     $                 - u1r(i)*rym1(k,1,1,1) 
     $                 - u1s(i)*sym1(k,1,1,1) 
     $                 - u1t(i)*tym1(k,1,1,1) ) /jacm1(i,1,1,e)

          enddo 
          enddo
    
       else

        do e= 1,nelt 

           j= nxyz*(e-1) 

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1) 
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

          k = i+j

          w1(k) = ( u3r(i)*rym1(k,1,1,1) 
     $            + u3s(i)*sym1(k,1,1,1)) /jacm1(i,1,1,e)

          w2(k) =-( u3r(i)*rxm1(k,1,1,1) 
     $            + u3s(i)*sxm1(k,1,1,1)) /jacm1(i,1,1,e)

          w3(k) = ( u2r(i)*rxm1(k,1,1,1) 
     $            + u2s(i)*sxm1(k,1,1,1) 
     $            - u1r(i)*rym1(k,1,1,1) 
     $            - u1s(i)*sym1(k,1,1,1)) /jacm1(i,1,1,e)

          enddo 

        enddo 
       endif

       return
       end
c-----------------------------------------------------------------------
      subroutine wght_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)

c
c     Weighted curl:   w = B*curl(u)
c
      include 'SIZE'
      include 'TOTAL'
c
      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
      logical ifavg
      integer e, k, i, j, nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
c
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            w1(k) = (    u3r(i)*rym1(k,1,1,1)
     $                 + u3s(i)*sym1(k,1,1,1)
     $                 + u3t(i)*tym1(k,1,1,1)
     $                 - u2r(i)*rzm1(k,1,1,1)
     $                 - u2s(i)*szm1(k,1,1,1)
     $                 - u2t(i)*tzm1(k,1,1,1) ) * w3m1(i,1,1)
c
            w2(k) = (    u1r(i)*rzm1(k,1,1,1)
     $                 + u1s(i)*szm1(k,1,1,1)
     $                 + u1t(i)*tzm1(k,1,1,1) 
     $                 - u3r(i)*rxm1(k,1,1,1)
     $                 - u3s(i)*sxm1(k,1,1,1) 
     $                 - u3t(i)*txm1(k,1,1,1) ) * w3m1(i,1,1)
c
            w3(k) = (    u2r(i)*rxm1(k,1,1,1)
     $                 + u2s(i)*sxm1(k,1,1,1) 
     $                 + u2t(i)*txm1(k,1,1,1) 
     $                 - u1r(i)*rym1(k,1,1,1) 
     $                 - u1s(i)*sym1(k,1,1,1) 
     $                 - u1t(i)*tym1(k,1,1,1) ) * w3m1(i,1,1)

          enddo 
          enddo

       else

        do e= 1,nelt 

           j= nxyz*(e-1) 

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1) 
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

          k = i+j

          w1(k) = ( u3r(i)*rym1(k,1,1,1) 
     $            + u3s(i)*sym1(k,1,1,1)) * w3m1(i,1,1)

          w2(k) =-( u3r(i)*rxm1(k,1,1,1) 
     $            + u3s(i)*sxm1(k,1,1,1)) * w3m1(i,1,1) 

          w3(k) = ( u2r(i)*rxm1(k,1,1,1) 
     $            + u2s(i)*sxm1(k,1,1,1) 
     $            - u1r(i)*rym1(k,1,1,1) 
     $            - u1s(i)*sym1(k,1,1,1)) * w3m1(i,1,1) 

          enddo 

        enddo 
       endif

       return
       end

c-----------------------------------------------------------------------
      subroutine wght_curl_2d3(w1,w2,w3,u1,u2,u3,ifavg,work1,work2) 
c-----------------------------------------------------------------------

       include 'SIZE' 
       include 'TOTAL' 
       include 'EMWAVE' 
c
c     Combined weighted curl:   w12 = B*curl(u3), w3=B*curl(u1,u2)
c
                        
       real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1) 
       logical ifavg
                                    
       integer e
       common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             ,  u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             ,  u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
       nxyz  = nx1*ny1*nz1
       ntot  = nx1*ny1*nz1*nelt
       nn    = nx1-1

       do e=1,nelt

          j = nxyz*(e-1)

          call local_grad2(u1r,u1s,u1,nn,e,dxm1,dxtm1)
          call local_grad2(u2r,u2s,u2,nn,e,dxm1,dxtm1)
          call local_grad2(u3r,u3s,u3,nn,e,dxm1,dxtm1)

          do i=1,nxyz

             k = i+j 

             w1(k) = ( u3r(i)*rym1(k,1,1,1)
     $               + u3s(i)*sym1(k,1,1,1)) * w3m1(i,1,1)

             w2(k) =-( u3r(i)*rxm1(k,1,1,1) 
     $               + u3s(i)*sxm1(k,1,1,1)) * w3m1(i,1,1)

             w3(k) = ( u2r(i)*rxm1(k,1,1,1)
     $               + u2s(i)*sxm1(k,1,1,1)
     $               - u1r(i)*rym1(k,1,1,1)
     $               - u1s(i)*sym1(k,1,1,1) )* w3m1(i,1,1)

          enddo

       enddo

       return
       end

c-----------------------------------------------------------------------
      subroutine cem_div (w1,w2,w3,u1,u2,u3)
c     div:   w = div(u)/jacm1
      include 'SIZE'
      include 'TOTAL'
c
      real     w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer  e                           

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
 
      nxyz  = nx1*ny1*nz1
      nn    = nx1-1
c
      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(ur1,us1,ut1,u1,nn,e,dxm1,dxtm1)
         call local_grad3(ur2,us2,ut2,u2,nn,e,dxm1,dxtm1)
         call local_grad3(ur3,us3,ut3,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (  ur1(i)*rxm1(k,1,1,1)
     $               + us1(i)*sxm1(k,1,1,1)
     $               + ut1(i)*txm1(k,1,1,1))/jacm1(i,1,1,e)
            w2(k) = (  ur2(i)*rym1(k,1,1,1)
     $               + us2(i)*sym1(k,1,1,1)
     $               + ut2(i)*tym1(k,1,1,1))/jacm1(i,1,1,e)
            w3(k) = (  ur3(i)*rzm1(k,1,1,1)
     $               + us3(i)*szm1(k,1,1,1)
     $               + ut3(i)*tzm1(k,1,1,1))/jacm1(i,1,1,e)
         enddo

         enddo

      else

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad2(ur1,us1,u1,nn,e,dxm1,dxtm1)
         call local_grad2(ur2,us2,u2,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j
            w1(k) = (  ur1(i)*rxm1(k,1,1,1)
     $               + us1(i)*sxm1(k,1,1,1))/jacm1(i,1,1,e)
            w2(k) = (  ur2(i)*rym1(k,1,1,1)
     $               + us2(i)*sym1(k,1,1,1))/jacm1(i,1,1,e) 
            !write(6,*) '3d',e,i,u1(k)           
            !write(6,*) '3d',e,i,w1(k),w2(k)
            !write(6,*) '3d',e,i,rxm1(k,1,1,1)jacm1(i,1,1,e)
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)
C
C     -  Spectral interpolation from A to B via tensor products
C     -  scratch arrays: w(na*na*nb + nb*nb*na)
C
C     5/3/00  -- this routine replaces specmp in navier1.f, which
c                has a potential memory problem
C
C
      logical if3d
c
      real b(nb,nb,nb),a(na,na,na)
      real w(ldw)
c
      ltest = na*nb
      if (if3d) ltest = na*na*nb + nb*na*na
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine specmpnface(b,nb,a,na,ba,ab,if3d,w,ldw)
            
c     -  Spectral interpolation from A to B via tensor products
c     -  scratch arrays: w(na*na*nb + nb*nb*na)
c
c     5/3/00  -- this routine replaces specmp in navier1.f, which        
c                has a potential memory problem                           
          
      logical if3d                                                     
c          
      real b(nb,nb),a(na,na)
      real w(ldw)
c          
      ltest = na*nb
           
      if (if3d)      ltest = na*na*nb + nb*na*na
          
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif
c         
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine intp_rstd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      integer mx, md, idir
      logical if3d
c
      parameter       (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real            jgl,jgt,wkd
c
      parameter       (ld=2*lxd)
      common /ctmp0/  w(ld**ldim,2)

      integer   i

      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
 
      ldw = 2*(ld**ldim)
  
      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
         call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else
         call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif

      return
      end
             
c-----------------------------------------------------------------------
      subroutine intp_rsd(ju,u,mx,md,if3d,idir)
c
c     GLL interpolation from mx to md.
c
c     If idir ^= 0, then apply transpose operator  (md to mx)
c
      include 'SIZE'
c
      real    ju(1),u(1)
      integer mx, md, idir
      logical if3d
c
      parameter       (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/   d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                 jgt(ldg), wkd(lwkd)
      real             jgl,jgt,wkd
 
      parameter       (ld=2*lxd)
      common /ctmp0/   w(ld**ldim,2)

      integer   i

      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
 
      ldw = 2*(ld**ldim)
  
      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
          call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else 
          call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine gen_int(jgl,jgt,mp,np,w)
 
c     Generate interpolation from np GLL points to mp GL points
c
c        jgl  = interpolation matrix, mapping from velocity nodes to pressure
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c

      real jgl(mp,np),jgt(np*mp),w(1)
c
      iz = 1
      id = iz + np
c
      call zwgll (w(iz),jgt,np)
      call zwgll (w(id),jgt,mp)
c
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo
c
      call transpose(jgt,np,jgl,mp)
c
      return
      end
          
c-----------------------------------------------------------------------
      subroutine gen_dgl(dgl,dgt,mp,np,w)
c
c     Generate derivative from np GL points onto mp GL points
c
c        dgl  = interpolation matrix, mapping from velocity nodes to pressure
c        dgt  = transpose of interpolation matrix
c        w    = work array of size (3*np+mp)
c
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c
c
c
      real dgl(mp,np),dgt(np*mp),w(1)
c
c
      iz = 1
      id = iz + np  
c
      call zwgll  (w(iz),dgt,np)  ! GL points
      call zwgll  (w(id),dgt,mp)  ! GL points
c
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c
      n  = np -1 
      do i=1,mp   
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np  
            dgl(i,j) = dgt(np+j)                       ! Derivative matrix
         enddo
      enddo
c
      call transpose(dgt,np,dgl,mp)

c
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_chk(n,m,avar5,lvar5,sub_name10)
      include 'SIZE'            ! need nid
      character*5  avar5,lvar5
      character*10 sub_name10
c
      if (n.gt.m) then
         write(6,*)             
         write(6,*) nid,n,m
c        write(6,*) nid,n,m,avar5,lvar5,sub_name10
c        write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i9,2(1x,a5),1x,a10)
         call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md)
c
c     Get pointer to jgl() for interpolation pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip, mx, md
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c
      if (ip.eq.0) then
c
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c
         nwrkd = mx + md

         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
c
         call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgl_ptr (ip,mx,md)
c
c     Get pointer to GL-GL interpolation dgl() for pair (mx,md)
c
      include 'SIZE'
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip, mx, md
c
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c
      ij = md + ld*(mx-1)
      ip = pdg (ij)
c
      if (ip.eq.0) then
c
         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)
c
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c
         call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)

      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rstd(ur,us,ut,u,mx,md,if3d,ju)
c
      include 'SIZE'
      include 'DXYZ'
c
      real    ur(1),us(1),ut(1),u(1),ju(1)
      logical if3d
c
      parameter (ldg=lxd**3,lwkd=4*lx1)
      common /dgrad/  d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),
     $                jgt(ldg), wkd(lwkd)
      real jgl,jgt
      integer ip,mx,md,m0
c
      m0 = md-1
      if (mx.ne.md) then
c
         call intp_rstd(ju,u,mx,md,if3d,0) ! 0 = forward
c
         call get_dgl_ptr (ip,md,md)

         if (if3d) then
            call local_grad3(ur,us,ut,ju,m0,1,dg(ip),dgt(ip))
         else
            call local_grad2(ur,us   ,ju,m0,1,dg(ip),dgt(ip))
         endif
c
      else
         m0 = nx1-1
         if (if3d) then
            call local_grad3(ur,us,ut,u,m0,1,dxm1,dxtm1)
         else
            call local_grad2(ur,us   ,u,m0,1,dxm1,dxtm1)
         endif
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_options2     

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer k, nxyzf
      nxyzf = nx1*nz1*2*ndim*nelt

      IFVARDT  =.false.
      IFDIELEC =.false.
      IFIC     =.false.
      IFPOISSON=.false.
      IFDTCTL  =.false.
      IFSCAT   =.false.
      IFDRUDE  =.false.
      IFBEAM   =.false.
      IFWAKE   =.false.
      IFNAPY   =.false.
      IFSOL    =.false.
      IFVTK    =.false.

      if (param( 29).ne.0)  IFVARDT  =.true.
      if (param( 93).ne.0)  IFDIELEC =.true.
      if (param( 94).ne.0)  IFIC     =.true.
      if (param( 95).ne.0)  IFPOISSON=.true.
      if (param( 99).ne.0)  IFDTCTL  =.true.
      if (param(100).ne.0)  IFSCAT   =.true.
      if (param(101).eq.1)  IFDRUDE  =.true.
      if (param(101).eq.2)  IFBEAM   =.true.
      if (param(101).eq.3)  IFWAKE   =.true.

      if (param(101).eq.4)  then
          IFWAKE =.true.
          IFNAPY =.true.
      endif

      if (param(102).ne.0)  IFSOL    =.true.
      if (param(103).ne.0)  IFVTK    =.true.

      write(6,*) 'IFDRUDE     ',IFDRUDE
      write(6,*) 'IFBEAM      ',IFBEAM 
      write(6,*) 'IFWAKE      ',IFWAKE 
      write(6,*) 'IFNAPY      ',IFNAPY 
      write(6,*) 'IFVARDT     ',IFVARDT
      write(6,*) 'IFPBG       ',IFPBG                     
      write(6,*) 'IFDIELEC    ',IFDIELEC 
      write(6,*) 'IFPOISSON   ',IFPOISSON 
      write(6,*) 'IFIC        ',IFIC            
      write(6,*) 'IFSOL       ',IFSOL                      
      write(6,*) 'IFSRC       ',IFSRC
      write(6,*) 'IFCENTRAL   ',IFCENTRAL
      write(6,*) 'IFUPWIND    ',IFUPWIND
      write(6,*) 'IFDEALIAS   ',IFDEALIAS 
      write(6,*) 'IFTM        ',IFTM
      write(6,*) 'IFTE        ',IFTE
      write(6,*) 'IFEIG       ',IFEIG
      write(6,*) 'IFRK4       ',IFRK4 
      write(6,*) 'IFEXP       ',IFEXP
      write(6,*) 'IFDTCTL     ',IFDTCTL
      write(6,*) 'IFNM        ',IFNM
      write(6,*) 'IFVTK       ',IFVTK

      if (IFNM) then
          cspeed  = param(1)
          permit_0= param(2)
          permea_0= param(3)
          permit_r= param(4)
          permea_r= param(5) ! not used
          gamma_p = param(6)
          freq_p  = param(7)
          xyzscale= 0.0
          if (IFDRUDE) xyzscale=1.e+9 ! for visual tool
          if (IFBEAM ) xyzscale=1.0  
          if (IFWAKE ) xyzscale=1.0 
      else
          cspeed  = 1.0          
          permit_0= 1.0         
          permea_0= 1.0             
          permit_r= param(4)
          permea_r= param(5) ! not used
          gamma_p = param(8)
          freq_p  = param(9)
          xyzscale= 1.0 
      endif

      write(6,*) ' '             
      write(6,*) 'cspeed   =',cspeed
      write(6,*) 'permit_0 =',permit_0
      write(6,*) 'permea_0 =',permea_0
      write(6,*) 'permit_r =',permit_r
      write(6,*) 'permea_r =',permea_r
      write(6,*) ' '             
      write(6,*) 'fscale   =',param(59) 
      write(6,*) ' '             

      if     ((IFCENTRAL).and.(IFUPWIND))  then
             write(6,*) 'ERROR:: not allowed IFCENTRAL=T & IFUPWIND= T'
             call exitt
      elseif ((.not.(IFCENTRAL)).and.(.not.(IFUPWIND))) then
             write(6,*) 'ERROR:: not allowed IFCENTRAL=F & IFUPWIND= F'
             call exitt
      endif

      if     (.not.IF3D) then
         if     ((IFTM).and.(IFTE))  then
             write(6,*) 'ERROR:: not allowed IFTM= T & IFTE= T'
             call exitt
         elseif  ((.not.(IFTM)).and.(.not.(IFTE))) then
             write(6,*) 'ERROR:: not allowed IFTM= F & IFTE= F'
             call exitt
         endif
      endif 

      if ((param(29).ne.0).and.(param(84).eq.0)) call exitt

      write(6,*) ' '             
      write(6,*) 'cspeed   =',cspeed
      write(6,*) 'permit_0 =',permit_0
      write(6,*) 'permea_0 =',permea_0
      write(6,*) 'permit_r =',permit_r
      write(6,*) 'permea_r =',permea_r
      write(6,*) ' '             
      write(6,*) 'fscale   =',param(59) 
      write(6,*) ' '             

      call ifoption (IFPML) ! globalize logical value to all processors
      call ifoption (IFPEC)
      call assign_bdry_output

      if (if3d) then
        hxactive = .true.
        hyactive = .true.
        hzactive = .true.
        exactive = .true.
        eyactive = .true.
        ezactive = .true.
        
        fldidx(1) = 1
        fldidx(2) = 2
        fldidx(3) = 3
        fldidx(4) = 4
        fldidx(5) = 5
        fldidx(6) = 6
        fldcnt = 6
      elseif (iftm) then
        hxactive = .true.
        hyactive = .true.
        hzactive = .false.
        exactive = .false.
        eyactive = .false.
        ezactive = .true.

        fldidx(1) = 1 !hx
        fldidx(2) = 2
        fldidx(3) = 0 !hz
        fldidx(4) = 0 !ex
        fldidx(5) = 0
        fldidx(6) = 3 !ez
        fldcnt = 3
      elseif (ifte) then
        hxactive = .false.
        hyactive = .false.
        hzactive = .true.
        exactive = .true.
        eyactive = .true.
        ezactive = .false.

        fldidx(1) = 0 !hx
        fldidx(2) = 0
        fldidx(3) = 3 !hz
        fldidx(4) = 1 !ex
        fldidx(5) = 2
        fldidx(6) = 0 !ez
        fldcnt = 3
      endif

      if     (IFEIG) then
         write(6,*) '                                      '
         write(6,*) 'IFEIG T :: start eigenvalue computation'
         write(6,*) '                                      '
      else
         if     ((IFRK4).and.(IFEXP)) then
              write(6,*) 'ERROR :: not allowed  IFRK4= T & IFEXP= T'
              call exitt
         elseif ((.not.(IFRK4)).and.(.not.(IFEXP))) then
              write(6,*) 'ERROR :: not allowed  IFRK4 F &  IFEXP F'
              call exitt
         endif
      endif

      return
      end  
c-----------------------------------------------------------------------
      subroutine pout(u,n)        
      integer n
      real u(n)

      do i=1,n
         write(6,*) i,u(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine pout2(u1,u2,n)        
      integer n
      real u1(n)
      real u2(n)

      do i=1,n
         write(6,*) i,u1(i),u2(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_option_chk2
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer k, nxyzf
      nxyzf = nx1*nz1*2*ndim*nelt

      call assign_bdry_output                            
      !call face_value_chk   

      return
      end
c-----------------------------------------------------------------------
      subroutine ifoption(ifif)        

      include 'SIZE'
      include 'EMWAVE'
      include 'PARALLEL'

      logical ifif
 
      if (np.eq.1) return

      ifip =0
      if (ifif) ifip =1

      imaxg = iglmax(ifip,1)
      if (imaxg.eq.1) ifif =.true.

      return
      end
c-----------------------------------------------------------------------
c global z array pointer for parallel 
c izptr: computes 1D global z interface pointer index 
c zrr0: computes 1D global z interface arrays on GLL
c wrr0: computes 1D global z interface redundancy arrays on GLL
c zz0 : computes 1D global z arrays on GLL
c zw0 : computes 1D global z arrays on equidistant points
      subroutine get_glo_ptr_zarrays            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'WAKEPOT'
      real     w(1)
      integer  e

      nxyz = nx1*ny1*nz1
      zlmax= 1.e+20
      zlmin= 1.e+20
      do e=1,nelt
         zlmax0= vlmax(zm1(1,1,1,e),nxyz) ! local max on each node
         zlmin0= vlmin(zm1(1,1,1,e),nxyz) ! zlmin= zarry0(1)
         zlmax = min(zlmax,zlmax0)
         zlmin = min(zlmin,zlmin0)
      enddo

      n= nx1*ny1*nz1*nelt
      zgmax= glmax(zm1,n)        ! global max on global nodes
      zgmin= glmin(zm1,n)

c ... zgap = zarry0(2)-zarry0(1) ! constant: restriction to be fixed.
      zgap = zlmax-zlmin         ! constant: restriction to be fixed.
      zdist= (zlmin-zgmin)/zgap  ! zlmin>zgmin
      ndist= int(zdist+0.1)      ! distance from global minimum z: 0 1 2 3 ...
      zlevg= (zgmax-zgmin)/zgap  ! total levels in global
      nlevg= int(zlevg+0.1)      ! global
      nzrr0= nlevg+1             

      call rzero(zrr0 ,nzrr0)
      call rzero(wrr0 ,nzrr0)
      call rzero(work0,nzrr0)

c ... get global level array: zrr0
      do i=1,nzarry0
         izptr(i)= ndist+i        ! allows overlap izptr for nodes.
         zrr0(izptr(i))= zarry0(i)! z values on each processor
         wrr0(izptr(i))= 1.0      ! will use as redundant number
      enddo
      call gop(wrr0,work0,'+  ',nzrr0) !
      call gop(zrr0,work0,'+  ',nzrr0) ! redundancy
      call invcol2   (zrr0,wrr0,nzrr0) ! remove reduandancy

      nzlev= nz1*nlevg
      call zwgll(zi,wi,nz1)       ! GLL 
      do i = 1,nzarry0-1          ! redundancy chk array
         ind1= izptr(i)           ! global pointer
      do k = 1,nz1
         ind2= k+(ind1-1)*nz1
         zind(ind2)= 1.0
      enddo
      enddo
      call gop (zind,zwork,'+  ',nzlev) ! this is c

      nzz0= (nz1-1)*nlevg+1       ! no redundancy
      do i= 1,nlevg
      do k= 1,nz1-1
         ind= k+(i-1)*(nz1-1)
         zz0(ind)= (zrr0(i+1)-zrr0(i))/2.*(zi(k)+1.)+zrr0(i) !GLL
      enddo
      enddo
      zz0(nzz0)= zgmax  ! not used

      nuniform = luniform 
      nzw0= nuniform*nlevg        ! total nuniform1*nlevels in global
      do i= 1,nuniform
         zwh  = 2.0/nuniform      ! for equi-distant grids 
         zw(i)=-1.+(i-1)*zwh      ! not using the end point for rhs
      enddo

      do i= 1,nlevg
      do k= 1,nuniform
         ind= k+(i-1)*nuniform
         zw0(ind)=(zrr0(i+1)-zrr0(i))/2.*(zw(k)+1.)+zrr0(i)
         !write(6,*) 'tt-',i,k,ind,zw0(ind)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_zarrays           
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      ! FIXME  3/7/07 misun  parallel implementation needed
      integer  e

      ! base zarry with equidistant
      nelt2= 2*nelt
      call rzero (zin0  ,nelt2)
      call izero (iwork1,nelt2)

      nxyz=nx1*ny1*nz1
      do e=1,nelt   
         zlmin  =vlmin(zm1(1,1,1,e),nxyz)
         zlmax  =vlmax(zm1(1,1,1,e),nxyz)
         zin0(1,e)=zlmin 
         zin0(2,e)=zlmax 
      enddo

      nsort = 2*nelt
      tol   = 1.e-5 
      call reduced_sorts(zarry0,nzarry0,zin0,nsort,iwork1,tol)

      ! GLL points zarry with varibale dz
      n   =nx1*ny1*nz1*nelt
      call rzero (zin    ,n)
      call izero (workvec,n)
      call copy  (zin,zm1,n)

      nsort = n          
      tol   = 1.e-8 
      call reduced_sorts(zarry,nzarry,zin,nsort,workvec,tol)

      return
      end
c ---------------------------------------------------------
      subroutine put_face_value(face,bound0,nf)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real    bound0
      integer i,j,n,nf

      integer icalld
      save    icalld
      data    icalld /0/

      real face(nf)

      n    =nx1*ny1*nz1*nelt
      nface=2*ndim*nx1*nz1*nelt

      if (icalld.eq.0) then
          call cem_set_fc_ptr
          icalld = 1
      endif

      do j = 1,ncemface
         i = cemface(j)
         xx= xm1(i,1,1,1)
         yy= ym1(i,1,1,1)
         zz= zm1(i,1,1,1)

         rr=sqrt(xx**2+yy**2)

         if     (rr.le.bound0) then
            face(j)=-1   
         else                         
            face(j)= 1   
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine zero_out(uu,n)
      integer n
      real uu(n)
    
      do i =1,n
         if (abs(uu(i)).le.1e-20) then
             uu(i)=0.0
         else
             uu(i)=uu(i)
         endif   
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine field_normalization(u,v,w,n)
      real u(n)
      real v(n)
      real w(n)

      unorm=glsc2(u,n)
      vnorm=glsc2(v,n)
      wnorm=glsc2(w,n)

      call cmult(u,unorm,n)
      call cmult(v,vnorm,n)
      call cmult(w,wnorm,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine vertex_scale (a,b)
      include 'SIZE'
      include 'TOTAL'

c     a   : separate radius (if a=0, no separate radius)
c     b   : scaling factor

      real     a,b   
      integer  e


      rinner = a 
      ro_new = b                   
      router = 0.

        do e=1,nelv
         do i=1,8
            r = xc(i,e)**2 + yc(i,e)**2
            router = max(router,r)
         enddo
        enddo
        router = sqrt(router)

        do e=1,nelv
         ravg = 0
         do i=1,8
            r = xc(i,e)**2 + yc(i,e)**2
            if (r.gt.0) r = sqrt(r)
            ravg = ravg + r
         enddo
         ravg = ravg /8.

         if (ravg.gt.rinner) then ! we're outside the inner cyl
           do i=1,8
              r = xc(i,e)**2 + yc(i,e)**2
              if (r.gt.0) r = sqrt(r)
              s = (r-rinner)/(router-rinner)
              rnew    = rinner + s*(ro_new-rinner)
              theta   = atan2(yc(i,e),xc(i,e))
              xc(i,e) = rnew*cos(theta)
              yc(i,e) = rnew*sin(theta)
           enddo
           do i=1,8
           if ((ccurve(i,e).eq.'C').and.(curve(1,i,e).gt.rinner)) then
              ! curve (1,..): the first argument indicates radius
                r = curve(1,i,e)
                s = (r-rinner)/(router-rinner)
                rnew = rinner + s*(ro_new-rinner)
                curve(1,i,e) = rnew
           endif
           enddo
         endif
        enddo

        if (nid.eq.0) write(6,*) 'vertex scaled by factor=',b

      return
      end
c-----------------------------------------------------------------------
      subroutine vertex_theta ()
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer n,i
      real    radius_squared,radhere_squared
      real    gop_work
      
      n   = 4*(ndim-1)*nelt

      do i=1,n
         theta   = atan2(yc(i,1),xc(i,1))
         tin0(i) = theta
      enddo

      call sorts(tout0,tin0,iwork3,n)
      k = 1
      tarry0(1)=tout0(1)
      do i=2,n     
         if ((tout0(i)-tarry0(k)).gt.1e-14) then
             k=k+1
             tarry0(k)=tout0(i)
         endif
      enddo
      ntarry0=k

      write(6,* ) '.'           
      do i=1,ntarry0
         write(6,10) i,tarry0(i)
      enddo
   10 format(' theta  array:',i5,f15.7)
      write(6,* ) '.'           

      return
      end
c-----------------------------------------------------------------------
      subroutine vertex_adjust()
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      real xsize, ysize, zsize
      real xmin,xmax,ymin,ymax,zmin,zmax
      real radius

      real e1,e2

      do ir=1,nrarry0

      radius = rarry0(ir)
      radius = int(10*radius+0.1)/10.
      rarry0(ir) = radius

      e1 = radius*radius - 1e-4
      e2 = radius*radius + 1e-4

      do ie=1,nelt
      do i=1,8
         rr = xc(i,ie)*xc(i,ie) + yc(i,ie)*yc(i,ie)
         if ((rr.gt.e1).and.(rr.le.e2)) then
            if (rr.gt.0) then
                rr = sqrt(rr)
                rn = radius/rr
                xc(i,ie) = rn*xc(i,ie)
                yc(i,ie) = rn*yc(i,ie)
            endif
         endif
      enddo
      enddo

      enddo

      do i=1,nrarry0
         write(6,10) i,rarry0(i)
      enddo
      write(6,* ) '.'           
   10 format(' final rarray:',i5,f15.7)

      return
      end

c-----------------------------------------------------------------------
      subroutine fcsum2(xsum,asum,x,e,f)
c
c     Compute the weighted sum of X over face f of element e
c
c     x is an (NX,NY,NZ) data structure
c     f  is in the preprocessor notation 
c
c     xsum is sum (X*area)
c     asum is sum (area)

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'
      real x(lx1,ly1,lz1,1)
      integer e,f,fd

      asum = 0.
      xsum = 0.

c     Set up counters ;  fd is the dssum notation.
      call dsset(nx1,ny1,nz1)
      fd     = eface1(f)
      js1    = skpdat(1,fd)
      jf1    = skpdat(2,fd)
      jskip1 = skpdat(3,fd)
      js2    = skpdat(4,fd)
      jf2    = skpdat(5,fd)
      jskip2 = skpdat(6,fd)

      i = 0
      do j2=js2,jf2,jskip2
      do j1=js1,jf1,jskip1
         i = i+1
         xsum = xsum+area(i,1,f,e)*x(j1,j2,1,e)
         asum = asum+area(i,1,f,e)
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      real function geom_xyradius()
         
      ! return the maximum radius in the x-y plane 
      ! of the nodes in the mesh.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer n,i
      real    radius_squared,radhere_squared
      real    gop_work
      
      n   = 8*nelt
      radius_squared = 0

      do i=1,n
         radhere_squared = xc(i,1)*xc(i,1) + yc(i,1)*yc(i,1)
         if (radhere_squared.ge.radius_squared) then
             radius_squared = radhere_squared
         endif
      enddo

      call gop(radius_squared, gop_work, 'M  ', 1)
      geom_xyradius = sqrt(radius_squared)

      return
      end
c ---------------------------------------------------------
      real function final_geom_xyradius()

      ! return the maximum radius in the x-y plane 
      ! of the nodes in the mesh.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer n,i
      real radius_squared,radhere_squared
      real gop_work
      
      n = nx1*ny1*nz1*nelt
      radius_squared = 0

      do i=1,n
         radhere_squared = xm1(i,1,1,1)*xm1(i,1,1,1) 
     $                   + ym1(i,1,1,1)*ym1(i,1,1,1)
         if (radhere_squared.ge.radius_squared) then
             radius_squared = radhere_squared
         endif
      enddo

      call gop(radius_squared, gop_work, 'M  ', 1)
      final_geom_xyradius = sqrt(radius_squared)

      return
      end


