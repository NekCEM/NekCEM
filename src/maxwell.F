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
      subroutine cem_maxwell_temprary                
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'TMPINPUT'
       
c...  these to be deleted all when IFHEX/IFTET are done

          nxyz=nx1*ny1*nz1      
          npts=nx1*ny1*nz1*nelt
          npts3=nx1*ny1*nz1*nelt*3
          nxzf=nx1*nz1*2*ndim*nelt

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
          call rzero (shn,npts3)
          call rzero (sen,npts3)
          call rzero (hn,npts3)
          call rzero (en,npts3)
          call rzero (permittivity,npts)
          call rzero (permeability,npts)

          call rzero (w3mn,nxyz)
          call rzero (unxm,nxzf)
          call rzero (unym,nxzf)
          call rzero (unzm,nxzf)
          call rzero (aream,nxzf)
          call rzero (yy_0,nxzf)
          call rzero (yy_1,nxzf)
          call rzero (zz_0,nxzf)
          call rzero (zz_1,nxzf)

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
             aream(i)= area(i,1,1,1)
          enddo

          do i=1,nxyz
             w3mn(i)= w3m1 (i,1,1)
          enddo

          do i=1,nxzf
             unxm(i)=unx(i,1,1,1)
             unym(i)=uny(i,1,1,1)
             unzm(i)=unz(i,1,1,1)
             yy_0(i)=y_0(i,1,1)
             yy_1(i)=y_1(i,1,1)
             zz_0(i)=z_0(i,1,1)
             zz_1(i)=z_1(i,1,1)
          enddo

c...      be careful before deleting  
          call cem_set_fc_ptr

          call copy (permittivity,epsilo,npts)
          call copy (permeability,  mu  ,npts)
          call col3 (ebm1,permittivity,bmn,npts)
          call col3 (hbm1,permeability,bmn,npts)
          call invcol1(ebm1,npts)
          call invcol1(hbm1,npts)

          if (nid.eq.0) write(6,*) 'vvvvvvvvv---'
          call usersol(time,sHN(1,1),sHN(1,2),sHN(1,3),
     $                         sEN(1,1),sEN(1,2),sEN(1,3))
             call usersol(time, HN(1,1), HN(1,2), HN(1,3),
     $                          EN(1,1), EN(1,2), EN(1,3))

          gg=glmax(hn(1,3),npts)
          hh=glmax(en(1,1),npts)
          gg1=glmax(shn(1,3),npts)
          hh1=glmax(sen(1,1),npts)
          aa=glmax(jacm,npts)
          bb=glmax(bmn,npts)
          xxm=glmax(xmn,npts)

          if (nid.eq.0) write(6,*) 'tt---',gg,hh,gg1,hh1, time
          if (nid.eq.0) write(6,*) 'oo---',aa,bb,xxm, time
      return
      end

c--------------------------------------------------------------------- 
      subroutine cem_maxwell_op_rk                   
c--------------------------------------------------------------------- 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'RK5'
      integer  icalld
      save     icalld
      data     icalld/0/

c...  this is temporary: to be deleted all when IFHEX/IFTET are done
      if (icalld.eq.0)  call cem_maxwell_temprary
      icalld= 1

c ... 5-stage 4th-order RK
      do ii=1,5             
         call rk_c (ii)
         call cem_maxwell_op   
         call rk_maxwell_ab(ii)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)

      call cem_maxwell      
      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux (srflx)
      call cem_maxwell_add_flux_to_res (srflx)
c     call cem_source
c     call cem_apps    
      call cem_maxwell_invqmass

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell 
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'

      npts = nx1*ny1*nz1*nelt
      npts3= nx1*ny1*nz1*nelt

c...  sacrify 2D optimization for simplicity::
      if     (imode.eq.3) then !IF3D

      call maxwell_wght_curl (resEN(1,1),resEN(1,2),resEN(1,3),HN(1,1),
     $                                                 HN(1,2),HN(1,3))
      call maxwell_wght_curl (resHN(1,1),resHN(1,2),resHN(1,3),EN(1,1),
     $                                                 EN(1,2),EN(1,3))
      call chsign (resHN,npts3)

      elseif (imode.eq.2) then !IFTE

      call maxwell_wght_curl (resHN(1,1),resHN(1,2),resEN(1,3),HN(1,1),
     $                                                 HN(1,2),EN(1,3))
      call chsign (resHN(1,1),npts)
      call chsign (resHN(1,2),npts)

      elseif (imode.eq.1) then !IFTE
      call maxwell_wght_curl (resEN(1,1),resEN(1,2),resHN(1,3),EN(1,1),
     $                                                 EN(1,2),HN(1,3))
      call chsign (resHN(1,3),npts)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_restrict_to_face
c-----------------------------------------------------------------------
c     restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'

c...  sacrify efficiency for 2D for simplicity
      do j=1, ncemface
           i= cemface(j)
           fHN(j,1) = HN(i,1)
           fHN(j,2) = HN(i,2)
           fHN(j,3) = HN(i,3)
           fEN(j,1) = EN(i,1)
           fEN(j,2) = EN(i,2)
           fEN(j,3) = EN(i,3)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_add_flux_to_res(srflx)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'TMPINPUT'

      real     srflx(1)

      nxzf= nx1*nz1*2*ndim*nelt
      k   = nxzf               

      if     (imode.eq.3) then !IF3D
        do j = 1,ncemface
           i = cemface(j)
           a = aream(j)
           resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
           resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
           resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
           resEN(i,1) = resEN(i,1) + a*srflx(3*k+j)
           resEN(i,2) = resEN(i,2) + a*srflx(4*k+j)
           resEN(i,3) = resEN(i,3) + a*srflx(5*k+j)
        enddo
      elseif (imode.eq.2) then !IFTM
        do j=  1,ncemface
           i = cemface(j)
           a = aream(j)
           resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
           resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
           resEN(i,3) = resEN(i,3) + a*srflx(2*k+j)
        enddo
      elseif (imode.eq.1) then !IFTE
        do j = 1,ncemface
           i = cemface(j)
           a = aream(j)
           resEN(i,1) = resEN(i,1) + a*srflx(0*k+j)
           resEN(i,2) = resEN(i,2) + a*srflx(1*k+j)
           resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux(srflx)
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'
      include 'MAXWELL'

      real     srflx(1)

      if (if3d) then
         call cem_maxwell_flux3d (srflx)
      else
         call cem_maxwell_flux2d (srflx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux2d (srflx)  
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'PARALLEL'
      include 'TMPINPUT'

      real*8  dclock,t0
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1
      integer nxyzf,i,j,k,e,f,nxz

      nxyzf = nx1*nz1*2*ndim*nelt
      nxz   = nx1*nz1

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0

      if     (imode.eq.2) then !IFTM

        k = nxyzf
        do i=1,nxyzf                
           j=cemface(i)
         srflx(0*k+i)=-unym(i)*fEN(i,3)  ! NxfHx
         srflx(1*k+i)= unxm(i)*fEN(i,3)  ! NxfHy
         srflx(2*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1) ! NxfEz
        enddo

        t0=dclock()
        call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)
        call measure_comm(t0)

        if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)       

        k = nxyzf
        do i=1,nxyzf

         Y0 = YY_0(i) 
         Y1 = YY_1(i)  
         Z0 = ZZ_0(i)   
         Z1 = ZZ_1(i)    

         fu1 = unym(i)*srflx(2*k+i)                                
         fu2 =-unxm(i)*srflx(2*k+i)                     
         fw3 = unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)     

         srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)    
         srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)   
         srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)   
        enddo

      elseif (imode.eq.1) then !IFTE

        k = nxyzf
        do i=1,nxyzf                

         srflx(0*k+i)= -unym(i)*fHN(i,3)  ! NxfHx   
         srflx(1*k+i)=  unxm(i)*fHN(i,3)  ! NxfHy   
         srflx(2*k+i)= -unxm(i)*fEN(i,2)+unym(i)*fEN(i,1) ! NxfEz

        enddo

           
        t0 = dclock()
        call gs_op_fields(gsh_face,srflx,nxyzf,3,1,1,0)
        call measure_comm(t0)

        if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)        

        k = nxyzf
        do i=1,nxyzf                

         Y0 = YY_0(i)
         Y1 = YY_1(i)
         Z0 = ZZ_0(i)
         Z1 = ZZ_1(i)

         fw1 = unym(i)*srflx(2*k+i)                            
         fw2 =-unxm(i)*srflx(2*k+i)
         fu3 = unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

         srflx(0*k+i) = 0.5/Z0*( Z1*srflx(0*k+i) - C0*fw1)      
         srflx(1*k+i) = 0.5/Z0*( Z1*srflx(1*k+i) - C0*fw2)       
         srflx(2*k+i) = 0.5/Y0*(-Y1*srflx(2*k+i) - C0*fu3) 

        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux3d(srflx)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'PARALLEL'
      include 'TMPINPUT'

c      n x dEx - c* n x n x dHx  
c
      real*8  dclock,t0
      real    srflx(1)

      nxyzf =nx1*nz1*2*ndim*nelt ! face count
      nxz   =nx1*nz1
      
      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex
      k = nxyzf
      do i= 1,nxyzf
         ! - n x E
         srflx(0*k+i)=-unym(i)*fEN(i,3)+unzm(i)*fEN(i,2)
         srflx(1*k+i)=-unzm(i)*fEN(i,1)+unxm(i)*fEN(i,3)
         srflx(2*k+i)=-unxm(i)*fEN(i,2)+unym(i)*fEN(i,1)
         ! - n x H
         srflx(3*k+i)=-unym(i)*fHN(i,3)+unzm(i)*fHN(i,2)
         srflx(4*k+i)=-unzm(i)*fHN(i,1)+unxm(i)*fHN(i,3)
         srflx(5*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
      enddo

      t0= dclock()
      call gs_op_fields(gsh_face,srflx,nxyzf,6,1,1,0)
      call measure_comm(t0)

      !srflx now contains n x [E] and n x [H],with [F] = F^+ - F^-;(observe un^+ = - un^-)
      if (IFPEC.or.IFPML)  call cem_maxwell_flux_pec(srflx)            

      k = nxyzf
      do i=1,nxyzf

        Y0   = YY_0(i)
        Y1   = YY_1(i)
        Z0   = ZZ_0(i)
        Z1   = ZZ_1(i)
        Y02  =-0.5/Y0*Y1          
        Z02  = 0.5/Z0*Z1         
        C02Y = 0.5/Y0*C0         
        C02Z = 0.5/Z0*C0         

        ! fu = n x n x [H]
        fu1 =unym(i)*srflx(5*k+i)-unzm(i)*srflx(4*k+i)
        fu2 =unzm(i)*srflx(3*k+i)-unxm(i)*srflx(5*k+i)
        fu3 =unxm(i)*srflx(4*k+i)-unym(i)*srflx(3*k+i)

        ! fw = n x n x [E]
        fw1 =unym(i)*srflx(2*k+i)-unzm(i)*srflx(1*k+i)
        fw2 =unzm(i)*srflx(0*k+i)-unxm(i)*srflx(2*k+i)
        fw3 =unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)

        srflx(0*k+i) = Y02*srflx(0*k+i) - C02Y*fu1              
        srflx(1*k+i) = Y02*srflx(1*k+i) - C02Y*fu2
        srflx(2*k+i) = Y02*srflx(2*k+i) - C02Y*fu3          
        srflx(3*k+i) = Z02*srflx(3*k+i) - C02Z*fw1          
        srflx(4*k+i) = Z02*srflx(4*k+i) - C02Z*fw2         
        srflx(5*k+i) = Z02*srflx(5*k+i) - C02Z*fw3         

      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux_pec(srflx)                          
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'PARALLEL'
      include 'NEKUSE'
      include 'TMPINPUT'
 
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
      subroutine cem_maxwell_apps 
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'

c     if (ifdrude)  call cem_maxwell_drude
c     if (ifloren)  call cem_maxwell_lorentz

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_uvp           
C----------------------------------------------------------------------

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'

      integer e

      do e=1,nelt
         call nekuvp(e)
      enddo

c     call cem_mater1 
c     call cem_mater2
c     call cem_mater3

      return
      end
c-----------------------------------------------------------------------
      subroutine maxwell_wght_curl(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      include 'SIZE'
      include 'TOTAL'
      include 'TMPINPUT' !temporary
c
      real    w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      integer e,k,i,j,nxyz,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
      nxyz  = nx1*ny1*nz1
      ntot  = nx1*ny1*nz1*nelt
      nn    = nx1-1


      if (if3d) then

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            u1rw = u1r(i)*w3mn(i)
            u1sw = u1s(i)*w3mn(i)
            u1tw = u1t(i)*w3mn(i)
            u2rw = u2r(i)*w3mn(i)
            u2sw = u2s(i)*w3mn(i)
            u2tw = u2t(i)*w3mn(i)
            u3rw = u3r(i)*w3mn(i)
            u3sw = u3s(i)*w3mn(i)
            u3tw = u3t(i)*w3mn(i)

            w1(k) =      u3rw*rymn(k)     
     $                 + u3sw*symn(k)  
     $                 + u3tw*tymn(k)  
     $                 - u2rw*rzmn(k)
     $                 - u2sw*szmn(k)
     $                 - u2tw*tzmn(k)                   
 
            w2(k) =      u1rw*rzmn(k)    
     $                 + u1sw*szmn(k)     
     $                 + u1tw*tzmn(k)     
     $                 - u3rw*rxmn(k)
     $                 - u3sw*sxmn(k) 
     $                 - u3tw*txmn(k)                    
 
            w3(k) =      u2rw*rxmn(k)
     $                 + u2sw*sxmn(k) 
     $                 + u2tw*txmn(k) 
     $                 - u1rw*rymn(k) 
     $                 - u1sw*symn(k) 
     $                 - u1tw*tymn(k)                 

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

c...         sgn1/sgn2 is assigned "-" or "+" when calling; sgn1/sgn2=1/-1 if TM, 
c...         sgn1/sgn2 is assigned "-" or "+" when calling; sgn1/sgn2=-1/1 if TE, 

             u1rw = u1r(i)*w3m1(i,1,1)
             u1sw = u1s(i)*w3m1(i,1,1)
             u2rw = u2r(i)*w3m1(i,1,1)
             u2sw = u2s(i)*w3m1(i,1,1)

             u3rw = u3r(i)*w3m1(i,1,1)
             u3sw = u3s(i)*w3m1(i,1,1)

             w1(k) = ( u3rw*rym1(k,1,1,1)
     $               + u3sw*sym1(k,1,1,1))

             w2(k) =-( u3rw*rxm1(k,1,1,1)
     $               + u3sw*sxm1(k,1,1,1))

             w3(k) = ( u2rw*rxm1(k,1,1,1)
     $               + u2sw*sxm1(k,1,1,1)
     $               - u1rw*rym1(k,1,1,1)
     $               - u1sw*sym1(k,1,1,1))

          enddo
          enddo

       endif

       return
       end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_invqmass
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'

      npts = nx1*ny1*nz1*nelt

      if     (imode.eq.3) then !IF3D

          do i=1,npts
             resHN(i,1)=resHN(i,1)*hbm1(i)
             resHN(i,2)=resHN(i,2)*hbm1(i)
             resHN(i,3)=resHN(i,3)*hbm1(i)
             resEN(i,1)=resEN(i,1)*ebm1(i)
             resEN(i,2)=resEN(i,2)*ebm1(i)
             resEN(i,3)=resEN(i,3)*ebm1(i)
          enddo

      elseif (imode.eq.2) then !IFTM

          do i=1,npts
             resHN(i,1)=reshN(i,1)*hbm1(i) 
             resHN(i,2)=reshN(i,2)*hbm1(i)
             resEN(i,3)=reseN(i,3)*ebm1(i)
          enddo
 
      elseif (imode.eq.1) then !IFTE

           do i=1,npts
             resHN(i,3)=resHN(i,3)*hbm1(i)
             resEN(i,1)=resEN(i,1)*ebm1(i)
             resEN(i,2)=resEN(i,2)*ebm1(i)
           enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_maxwell_ab(ii)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'MAXWELL'
      include 'RK5'

      integer ii

      npts = nx1*ny1*nz1*nelt
      npts3= nx1*ny1*nz1*nelt*3

c     kHX=rk4a(i)*kHX
      ca = rk4a(ii)
      cb = rk4b(ii)

      if     (imode.eq.3) then

      call rk4_upd(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
      call rk4_upd(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
      call rk4_upd(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)
      call rk4_upd(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
      call rk4_upd(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)
      call rk4_upd(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)

      elseif (imod.eq.2) then

      call rk4_upd(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
      call rk4_upd(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
      call rk4_upd(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)

      elseif (imod.eq.1) then

      call rk4_upd(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)
      call rk4_upd(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
      call rk4_upd(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)

      endif

      return
      end
c-----------------------------------------------------------------------
