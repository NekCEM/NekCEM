!> \file maxwell.F
!!
!! \brief Maxwell solver main subroutines          
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
c...  Globally defined variable: npts, npts3, nxyz, nxzf,nxzfl, nfaces
#define GPU_VALIDATE 0
c--------------------------------------------------------------------- 
      subroutine cem_maxwell_temporary                
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'
      integer  i,npts3

      npts3= 3*npts

c...  these to be deleted all when IFHEX/IFTET are done

          if (nid.eq.0) write(6,*) 'start: initialize geom and fields'
          !coordinates and geom
          call rzero (xmn ,npts)
          call rzero (ymn ,npts)
          call rzero (zmn ,npts)
          call rzero (bmn ,npts)
          call rzero (bminv,npts)
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

          call rzero (w3mn ,nxyz)

          call rzero (unxm ,nxzfl)
          call rzero (unym ,nxzfl)
          call rzero (unzm ,nxzfl)
          call rzero (aream,nxzfl)

          !lorentz-drude model
          call rzero(lorena,npts3)
          call rzero(lorenb,npts3)
          call rzero(lorenc,npts3)

          !material properties 
          call rzero(impede,npts)
          call rzero(conduc,npts)
          call rzero(drudea,npts)
          call rzero(drudeb,npts)
          call rzero(faceif,npts)
          call izero(if_in_region,nelt)

          !pml related arrays
          call rzero(pmlsigma,npts3)
          call rzero(pmldn,npts3)
          call rzero(pmlbn,npts3)
          call rzero(respmldn,npts3)
          call rzero(respmlbn,npts3)
          call rzero(respmlen,npts3)
          call rzero(respmlhn,npts3)
          call rzero(kpmlbn,npts3)
          call rzero(kpmldn,npts3)
          call izero(pmlptr,nelt)

          do i=1,npts
             xmn (i) = xm1  (i,1,1,1)
             ymn (i) = ym1  (i,1,1,1)
             zmn (i) = zm1  (i,1,1,1)
             bmn (i) = bm1  (i,1,1,1)
             rxmn(i) = rxm1 (i,1,1,1)
             rymn(i) = rym1 (i,1,1,1)
             rzmn(i) = rzm1 (i,1,1,1)
             sxmn(i) = sxm1 (i,1,1,1)
             symn(i) = sym1 (i,1,1,1)
             szmn(i) = szm1 (i,1,1,1)
             txmn(i) = txm1 (i,1,1,1)
             tymn(i) = tym1 (i,1,1,1)
             tzmn(i) = tzm1 (i,1,1,1)
             jacm(i) = jacm1(i,1,1,1)
             bminv(i)= 1.0/bmn(i)
          enddo

          do i=1,nxyz
             w3mn(i) = w3m1 (i,1,1)
          enddo

          do i=1,nxzfl
             unxm (i)= unx (i,1,1,1)
             unym (i)= uny (i,1,1,1)
             unzm (i)= unz (i,1,1,1)
             aream(i)= area(i,1,1,1)
          enddo
          if (nid.eq.0) write(6,*) 'done: initialization geom/fields'

c...      get global numbering on face and set material properties           
          call cem_set_fc_ptr        !global numbering index on face
          call cem_maxwell_uvp       !returns material properties
          call cem_maxwell_materials !other material related setting

c...      inverse mass matrix including material constants 
          call col3 (ebm1,permittivity,bmn,npts)
          call col3 (hbm1,permeability,bmn,npts)
          call invcol1 (ebm1,npts)
          call invcol1 (hbm1,npts)
          if (nid.eq.0) write(6,*) 'done: inverse mass matrix'

      return
      end
c--------------------------------------------------------------------- 
      subroutine cem_maxwell_uvp
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  ie, npts3
      real     glmax,glmin
      real     permit_max,permit_min
      real     permea_max,permea_min
      real     drudea_max,drudea_min
      real     drudeb_max,drudeb_min
      real     drudec_max,drudec_min
      real     tlorena_max,tlorena_min
      real     tlorenb_max,tlorenb_min
      real     tlorenc_max,tlorenc_min

      npts3 = 3*npts

c...  nekuvp returns material properites defined in uservp
      do ie=1,nelt
         call nekuvp(ie)
      enddo

c...  check permittivity/permeability setting 
      permit_max=glmax(permittivity,npts)
      permit_min=glmin(permittivity,npts)
      permea_max=glmax(permeability,npts)
      permea_min=glmin(permeability,npts)

      if (nid.eq.0) write(6,*) 'npts',npts,nx1,ny1,nz1,nelt
      if (nid.eq.0) write(6,*) 'permit_max=',permit_max
      if (nid.eq.0) write(6,*) 'permit_min=',permit_min
      if (nid.eq.0) write(6,*) 'permea_max=',permea_max
      if (nid.eq.0) write(6,*) 'permea_min=',permea_min

      if ((permit_max.le.0).or.(permit_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define permit  in uservp (.usr)'
          call exitt
      endif
      if ((permea_max.le.0).or.(permea_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define permea  in uservp (.usr)'
          call exitt
      endif

c...  drude parameters       
      if (IFDRUDE) then

      drudea_max=glmax(drudea,npts)
      drudea_min=glmin(drudea,npts)
      drudeb_max=glmax(drudeb,npts)
      drudeb_min=glmin(drudeb,npts)

      if (nid.eq.0) write(6,*) 'drudea_max=',drudea_max
      if (nid.eq.0) write(6,*) 'drudea_min=',drudea_min
      if (nid.eq.0) write(6,*) 'drudeb_max=',drudeb_max
      if (nid.eq.0) write(6,*) 'drudeb_min=',drudeb_min

      if ((drudea_max.le.0).or.(drudea_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define drudea  in uservp (.usr)'
          call exitt
      endif
      if ((drudeb_max.le.0).or.(drudeb_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define drudeb  in uservp (.usr)'
          call exitt
      endif

      endif

c...  drude-lorentz parameters       
      if (IFLOREN) then

      tlorena_max=glmax(lorena,npts3)
      tlorena_min=glmin(lorena,npts3)
      tlorenb_max=glmax(lorenb,npts3)
      tlorenb_min=glmin(lorenb,npts3)
      tlorenc_max=glmax(lorenc,npts3)
      tlorenc_min=glmin(lorenc,npts3)

      if (nid.eq.0) write(6,*) 'tlorena_max=',tlorena_max
      if (nid.eq.0) write(6,*) 'tlorena_min=',tlorena_min
      if (nid.eq.0) write(6,*) 'tlorenb_max=',tlorenb_max
      if (nid.eq.0) write(6,*) 'tlorenb_min=',tlorenb_min
      if (nid.eq.0) write(6,*) 'tlorenc_max=',tlorenc_max
      if (nid.eq.0) write(6,*) 'tlorenc_min=',tlorenc_min

      if ((tlorena_max.le.0).or.(tlorena_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define lorena  in uservp (.usr)'
          call exitt
      endif
      if ((tlorenb_max.le.0).or.(tlorenb_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define lorenb  in uservp (.usr)'
          call exitt
      endif
      if ((tlorenc_max.le.0).or.(tlorenc_min.le.0)) then
          if (nid.eq.0) write(6,*) 'define lorenc  in uservp (.usr)'
          call exitt
      endif

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_materials
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer   e,f,i,j
      character CB*3
      real      tmpep,tmpmu,tmpt1,tmpt2,const

      if (nid.eq.0) write(6,*) 'start: cem_maxwell_materials'
c...  temporary use for checking some values on volume:
      do e= 1,nelt
         do i=1,nxyz
            faceif(i+nxyz*(e-1)) = if_in_region(e)             
         enddo
      enddo

c...  impedence/conductance
      do i=1,npts
         tmpmu = permittivity (i)  
         tmpep = permeability (i)   
         tmpt1 = sqrt(tmpmu/tmpep)
         tmpt2 = sqrt(tmpep/tmpmu)

         impede(i) = tmpt1
         conduc(i) = tmpt2
      enddo

c...  restrict impedence/conductance to face
      do j = 1,ncemface
         i = cemface(j)
         zimpede(j) = impede(i)
         yconduc(j) = conduc(i)
         facef  (j) = faceif(i) 
      enddo

      call copy (Z_0,zimpede,nxzfl)
      call copy (Y_0,yconduc,nxzfl)

c...  addition of the face values 
c...  Z_0=Z_0^+ + Z_0^-, Y_0=Y_0^+ + Y_0^-            
      call gs_op_fields(gsh_face,Z_0,nxzfl,1,1,1,0)
      call gs_op_fields(gsh_face,Y_0,nxzfl,1,1,1,0)
c     do j = 1,ncemface
c        write(6,*) 'zz--',j,Z_0(j),Y_0(j)      
c     enddo

c...  special setting for PEC boundary 
c...  Z_0= 2*Z_0^-, Y_0= 2*Y_0^-            
      do e= 1,nelt            
      do f= 1,nfaces
      do j= 1,nxzf

         i= (e-1)*nxzf*nfaces + nxzf*(f-1) + j  ! face numbering
         CB = CBC(f,e,2)

         if (CB.eq.'PEC') then
             Y_0(i) = 2.0*Y_0(i)
             Y_1(i) = 2.0*Y_1(i)
             Z_0(i) = 2.0*Z_0(i)
             Z_1(i) = 2.0*Z_1(i)
         endif

      enddo
      enddo
      enddo

c...  Z_1 = Z_0 + zimpede
      call sub3 (Z_1,Z_0,zimpede,nxzfl)
      call sub3 (Y_1,Y_0,yconduc,nxzfl)

      const = 0.5
 
c...  Z_0 = 0.5 * Z_0
      call cmult (Z_0,const,nxzfl)
      call cmult (Y_0,const,nxzfl)

      if (nid.eq.0) write(6,*) 'done: cem_maxwell_materials'

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_rk                   
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  ii

      real*8  dclock
      real    stime0,stime5
c ... 5-stage 4th-order RK
      stime5=dclock()
      do ii=1,5             

         stime0=dclock()
         call rk_c (ii)
         call cem_maxwell_op   
         call rk_maxwell_ab(ii)
         stime0=dclock()-stime0

      enddo
      stime5=dclock()-stime5

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
          write(6,*) '----------------------------------::'
          write(6,*) 'cem_maxwell_op_rk            (sec)::',stime0
          write(6,*) 'cem_maxwell_op_rk per istep  (sec)::',stime5
          write(6,*) '----------------------------------::'
      endif
      endif
#endif


      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx
      real*8  dclock
      real    stime0,stime1,stime2,stime3,stime4,stime5,stime6,stime7

      stime1=dclock()
      call cem_maxwell      
      stime1=dclock()-stime1       

      stime2=dclock()
      call cem_maxwell_restrict_to_face
      stime2=dclock()-stime2       

      stime3=dclock()
      call cem_maxwell_flux (srflx)
      stime3=dclock()-stime3       

      stime4=dclock()
      call cem_maxwell_add_flux_to_res (srflx)
      stime4=dclock()-stime4       

      stime5=dclock()
      call cem_maxwell_apps    
      stime5=dclock()-stime5       

      stime6=dclock()
      call cem_maxwell_source
      stime6=dclock()-stime6       

      stime7=dclock()
      call cem_maxwell_invqmass
      stime7=dclock()-stime7       
      stime0=stime1+stime2+stime3+stime4+stime5+stime6+stime7

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then 
          write(6,*) 'cem_maxwell                  (sec)::',stime1 
          write(6,*) 'cem_maxwell_restrict_to_face (sec)::',stime2 
          write(6,*) 'cem_maxwell_flux             (sec)::',stime3 
          write(6,*) 'cem_maxwell_add_flux_to_res  (sec)::',stime4 
          write(6,*) 'cem_maxwell_apps             (sec)::',stime5 
          write(6,*) 'cem_maxwell_source           (sec)::',stime6 
          write(6,*) 'cem_maxwell_invqmass         (sec)::',stime7 
          write(6,*) 'cem_maxwell_op               (sec)::',stime0 
      endif
      endif
#endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

c...  sacrify 2D optimization for simplicity::
      if     (imode.eq.3) then !IF3D

      call maxwell_wght_curl (resEN(1,1),resEN(1,2),resEN(1,3),HN(1,1),
     $                                                 HN(1,2),HN(1,3))
      call maxwell_wght_curl (resHN(1,1),resHN(1,2),resHN(1,3),EN(1,1),
     $                                                 EN(1,2),EN(1,3))
      call chsign (resHN(1,1),npts)
      call chsign (resHN(1,2),npts)
      call chsign (resHN(1,3),npts)

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
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  i,j,k
      real     srflx(1),a

      k = nxzfl               

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
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

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
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

      real*8  dclock,t0
      real    srflx(1)
      real    C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k,e,f

      if  (IFCENTRAL)  C0 = 0.0
      if  (IFUPWIND )  C0 = 1.0

      if     (imode.eq.2) then !IFTM

        k = nxzfl
        do i=1,nxzfl                
           j=cemface(i)
           srflx(0*k+i)=-unym(i)*fEN(i,3)  ! NxfHx
           srflx(1*k+i)= unxm(i)*fEN(i,3)  ! NxfHy
           srflx(2*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1) ! NxfEz
        enddo

        t0=dclock()
        call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)
        call measure_comm(t0)

        if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)       

        k = nxzfl
        do i=1,nxzfl

         Y0 = Y_0(i) 
         Y1 = Y_1(i)  
         Z0 = Z_0(i)   
         Z1 = Z_1(i)    

         fu1 = unym(i)*srflx(2*k+i)                                
         fu2 =-unxm(i)*srflx(2*k+i)                     
         fw3 = unxm(i)*srflx(1*k+i)-unym(i)*srflx(0*k+i)     

         srflx(0*k+i)= 0.5/Y0*(-Y1*srflx(0*k+i)-C0*fu1)    
         srflx(1*k+i)= 0.5/Y0*(-Y1*srflx(1*k+i)-C0*fu2)   
         srflx(2*k+i)= 0.5/Z0*( Z1*srflx(2*k+i)-C0*fw3)   
        enddo

      elseif (imode.eq.1) then !IFTE

        k = nxzfl
        do i=1,nxzfl                

         srflx(0*k+i)= -unym(i)*fHN(i,3)  ! NxfHx   
         srflx(1*k+i)=  unxm(i)*fHN(i,3)  ! NxfHy   
         srflx(2*k+i)= -unxm(i)*fEN(i,2)+unym(i)*fEN(i,1) ! NxfEz

        enddo

           
        t0 = dclock()
        call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)
        call measure_comm(t0)

        if (IFPML.or.IFPEC) call cem_maxwell_flux_pec(srflx)        

        k = nxzfl
        do i=1,nxzfl                

         Y0 = Y_0(i)
         Y1 = Y_1(i)
         Z0 = Z_0(i)
         Z1 = Z_1(i)

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
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'

c      n x dEx - c* n x n x dHx  
c
      real*8  dclock,t0
      real    srflx(1)
      real    C0,y0,y1,z0,z1,y02,z02,c02y,c02z,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k

      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex
      k = nxzfl
      do i= 1,nxzfl
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
      call gs_op_fields(gsh_face,srflx,nxzfl,6,1,1,0)
      call measure_comm(t0)

      !srflx now contains n x [E] and n x [H],with [F] = F^+ - F^-;(observe un^+ = - un^-)
      if (IFPEC.or.IFPML)  call cem_maxwell_flux_pec(srflx)            

      k = nxzfl
      do i=1,nxzfl

        Y0   = Y_0(i)
        Y1   = Y_1(i)
        Z0   = Z_0(i)
        Z1   = Z_1(i)
c       if ((nid.eq.0).and.(i.eq.1)) write(6,*) 'tt====y0',Y0,Z1
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
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'EMWAVE'
      include 'PARALLEL'
      include 'NEKUSE'
 
      real      srflx(1)
      integer   k,i,j,l,e,f,ef          
      character CB*3

      integer  icalld
      save     icalld
      data     icalld /0/

      common  /BDRY1/ cempec(6*lx1*lz1*2*ldim*lelt)
      common  /BDRY2/ ncempec
      integer  cempec,ncempec

      k = nxzfl
      l = 0

      if (icalld.eq.0) then

         do e= 1,nelt            
         do f= 1,nfaces
         do j= 1,nxzf
            CB= CBC(f,e,2)
            i =(e-1)*nfaces*nxzf+(f-1)*nxzf+j
            if ((CB.eq.'PEC').or.(CB.eq.'PML')) then
              l= l+1
              cempec(l)= i
              ncempec  = l
             !write(6,*) 'PML: flux',ncempec,CB
            endif
         enddo
         enddo
         enddo
         icalld =1

      endif

      if (ncempec.eq.0) return
 
      k = nxzfl
      if     (imode.eq.3) then !IF3D

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

          srflx(0*k+i)=2.0*srflx(0*k+i)
          srflx(1*k+i)=2.0*srflx(1*k+i)
          srflx(2*k+i)=2.0*srflx(2*k+i)
          srflx(3*k+i)=0
          srflx(4*k+i)=0
          srflx(5*k+i)=0
         enddo

      elseif (imode.eq.2) then !IFTM

         do j= 1,ncempec         
            i= cempec(j)
            srflx(0*k+i)=2.0*srflx(0*k+i)
            srflx(1*k+i)=2.0*srflx(1*k+i)
            srflx(2*k+i)=0
         enddo

      elseif (imode.eq.1) then !IFTE

         do j= 1,ncempec         
            i= cempec(j)
            srflx(0*k+i)=0
            srflx(1*k+i)=0
            srflx(2*k+i)=2.0*srflx(2*k+i)
         enddo

      endif
    
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_apps 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

      if (ifpml )  call pml_step
c     if (ifdrude)  call cem_maxwell_drude
c     if (ifloren)  call cem_maxwell_lorentz

      return
      end

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
      subroutine maxwell_wght_curl(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      include 'SIZE'
      include 'TOTAL'
c
      real    w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
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
#ifdef GPU
c      common /ctmp1/ u1r(lpts1),u1s(lpts1),u1t(lpts1) 
c     $             , u2r(lpts1),u2s(lpts1),u2t(lpts1)
c     $             , u3r(lpts1),u3s(lpts1),u3t(lpts1)
c      real u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t    
      real u1r(lpts1),u1s(lpts1),u1t(lpts1) 
     $   , u2r(lpts1),u2s(lpts1),u2t(lpts1)
     $   , u3r(lpts1),u3s(lpts1),u3t(lpts1)
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
#if GPU_VALIDATE
      real cu1r(nxyz),cu1s(nxyz),cu1t(nxyz)
     $    ,cu2r(nxyz),cu2s(nxyz),cu2t(nxyz)
     $    ,cu3r(nxyz),cu3s(nxyz),cu3t(nxyz)
      real cw1(1),cw2(1),cw3(1)
#endif
#else 
      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)
      real u1r,u2r,u3r,u1s,u2s,u3s,u1t,u2t,u3t    
#endif

c.............
#ifdef GPU
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

#if GPU_VALIDATE
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
      write(6,*) 'maxwell_wght_curl            (sec)::',wght_curl_time
      write(6,*) '        local_grad3_gpu::gpu_time1::',gpu_time1     
      write(6,*) '                         gpu_time2::',gpu_time2
      write(6,*) '----------------------------------::'
      endif
      endif
#endif

      return
#endif
c.............

      nn    = nx1-1

      wght_curl_time=dclock()

      if (if3d) then

         cpu_time1_total=0.0
         cpu_time2_total=0.0

         do e=1,nelt

         j = nxyz*(e-1)

         cpu_time1=dclock()           
         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)        
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)
         cpu_time1=dclock()-cpu_time1
         cpu_time1_total=cpu_time1_total+cpu_time1

         cpu_time2=dclock()           
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
            rxmn_k= rxmn(k)
            sxmn_k= sxmn(k)
            txmn_k= txmn(k)
            rymn_k= rymn(k)
            symn_k= symn(k)
            tymn_k= tymn(k)
            rzmn_k= rzmn(k)
            szmn_k= szmn(k)
            tzmn_k= tzmn(k)

         w1(k)= u3rw*rymn_k
     $        + u3sw*symn_k
     $        + u3tw*tymn_k
     $        - u2rw*rzmn_k
     $        - u2sw*szmn_k
     $        - u2tw*tzmn_k

         w2(k)= u1rw*rzmn_k
     $        + u1sw*szmn_k
     $        + u1tw*tzmn_k
     $        - u3rw*rxmn_k
     $        - u3sw*sxmn_k
     $        - u3tw*txmn_k

         w3(k)= u2rw*rxmn_k
     $        + u2sw*sxmn_k
     $        + u2tw*txmn_k
     $        - u1rw*rymn_k
     $        - u1sw*symn_k
     $        - u1tw*tymn_k

          enddo 
          cpu_time2      = dclock()-cpu_time2
          cpu_time2_total= cpu_time2_total+cpu_time2
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

             u1rw = u1r(i)*w3mn(i)
             u1sw = u1s(i)*w3mn(i)
             u2rw = u2r(i)*w3mn(i)
             u2sw = u2s(i)*w3mn(i)

             u3rw = u3r(i)*w3mn(i)
             u3sw = u3s(i)*w3mn(i)

             rxmn_k= rxmn(k)
             sxmn_k= sxmn(k)
             rymn_k= rymn(k)
             symn_k= symn(k)

             w1(k) = ( u3rw*rymn_k
     $               + u3sw*symn_k )

             w2(k) =-( u3rw*rxmn_k
     $               + u3sw*sxmn_k )

             w3(k) = ( u2rw*rxmn_k 
     $               + u2sw*sxmn_k 
     $               - u1rw*rymn_k 
     $               - u1sw*symn_k )

          enddo
          enddo

       endif

      wght_curl_time=dclock()-wght_curl_time

#ifdef TIMER 
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
      write(6,*) '----------------------------------::'
      write(6,*) 'maxwell_wght_curl            (sec)::',wght_curl_time
      write(6,*) '     3*local_grad3*nelt::cpu_time1::',cpu_time1_total
      write(6,*) '                         cpu_time2::',cpu_time2_total
      write(6,*) '----------------------------------::'
      endif
      endif
#endif

       return
       end
c-----------------------------------------------------------------------
      subroutine maxwell_curl(w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     weighted curl:   (w1,w2,w3) = curl x (u1,u2,u3)
      include 'SIZE'
      include 'TOTAL'
c
      real    w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      real    u1rw,u1sw,u1tw
      real    u2rw,u2sw,u2tw
      real    u3rw,u3sw,u3tw
      integer e,k,i,j,ntot,nn
      integer icalld
      save    icalld
      data    icalld/0/
      real*8  dclock    
      real    varw
      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
     $             , u2r(lx1*ly1*lz1),u2s(lx1*ly1*lz1),u2t(lx1*ly1*lz1)
     $             , u3r(lx1*ly1*lz1),u3s(lx1*ly1*lz1),u3t(lx1*ly1*lz1)

      real  u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t

      nn    = nx1-1

         do e=1,nelt

         j = nxyz*(e-1)

         call local_grad3(u1r,u1s,u1t,u1,nn,e,dxm1,dxtm1)        
         call local_grad3(u2r,u2s,u2t,u2,nn,e,dxm1,dxtm1)
         call local_grad3(u3r,u3s,u3t,u3,nn,e,dxm1,dxtm1)

         do i=1,nxyz

            k = i+j

            u1rw = u1r(i)/jacm(k)         
            u1sw = u1s(i)/jacm(k)
            u1tw = u1t(i)/jacm(k)
            u2rw = u2r(i)/jacm(k)
            u2sw = u2s(i)/jacm(k)
            u2tw = u2t(i)/jacm(k)
            u3rw = u3r(i)/jacm(k)
            u3sw = u3s(i)/jacm(k)
            u3tw = u3t(i)/jacm(k)

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

      return
      end
c-----------------------------------------------------------------------
      subroutine maxwell_grad(w1,w2,w3,u1)
c-----------------------------------------------------------------------
      implicit none
c     Non-weighted grad:  w1=du1/dx, w2=du1/dy, w3=du1/dz
      include 'SIZE'
      include 'TOTAL'
c
      integer mode, imode
      real    w1(1),w2(1),w3(1),u1(1)
      real    u1rw,u1sw,u1tw
      integer e,k,i,j,ntot,nn

      common /ctmp1/ u1r(lx1*ly1*lz1),u1s(lx1*ly1*lz1),u1t(lx1*ly1*lz1)
      real    u1r,u1s,u1t
c
      nn    = nx1-1

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

             w2(k) = ( u1r(i)*rymn(k)
     $               + u1s(i)*symn(k))/jacm(k)
             w3(k) = 0    
          enddo
          enddo


      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine maxwell_div (w0,w1,w2,w3,u1,u2,u3)
c-----------------------------------------------------------------------
      implicit none
c     div:     w = div(u)/jacm
      include 'SIZE'
      include 'TOTAL'
c
      real     w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1)
      integer  e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real  ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3

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

            w1(k) = (  ur1(i)*rxmn(k)
     $               + us1(i)*sxmn(k)
     $               + ut1(i)*txmn(k) )/jacm(k)

            w2(k) = (  ur2(i)*rymn(k)
     $               + us2(i)*symn(k)
     $               + ut2(i)*tymn(k) )/jacm(k)

            w3(k) = (  ur3(i)*rzmn(k)
     $               + us3(i)*szmn(k)
     $               + ut3(i)*tzmn(k) )/jacm(k)

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

            w1(k) = (  ur1(i)*rxmn(k)
     $               + us1(i)*sxmn(k) )/jacm(k)
            w2(k) = (  ur2(i)*rymn(k)
     $               + us2(i)*symn(k) )/jacm(k)
         
            w0(k) = w1(k)+w2(k)       
         enddo

         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_source
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      if (.not.ifsrc) return

      call usersrc (57,reshn(1,1),reshn(1,2),reshn(1,3),
     $                 resen(1,1),resen(1,2),resen(1,3))

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_invqmass
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'
      integer  i

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
             resEN(i,1)=resEN(i,1)*ebm1(i)
             resEN(i,2)=resEN(i,2)*ebm1(i)
             resHN(i,3)=resHN(i,3)*hbm1(i)
           enddo

      endif

      if  (ifpml) then

           do i=1,npts
             respmlBN(i,1)=respmlBN(i,1)*bminv(i)   
             respmlBN(i,2)=respmlBN(i,2)*bminv(i)  
             respmlBN(i,3)=respmlBN(i,3)*bminv(i)   
             respmlDN(i,1)=respmlDN(i,1)*bminv(i)   
             respmlDN(i,2)=respmlDN(i,2)*bminv(i)  
             respmlDN(i,3)=respmlDN(i,3)*bminv(i)   
           enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_maxwell_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'
      include 'PML'
      real     ca,cb
      integer  ii
      real*8 dclock
      real   stime0

      stime0=dclock()

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

      elseif (imode.eq.2) then

      call rk4_upd(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
      call rk4_upd(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
      call rk4_upd(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)

      elseif (imode.eq.1) then

      call rk4_upd(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
      call rk4_upd(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)
      call rk4_upd(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)

      endif

      if (ifpml) then

      call rk4_upd(pmlbn(1,1),kpmlbn(1,1),respmlbn(1,1),cb,ca,dt,npts)
      call rk4_upd(pmlbn(1,2),kpmlbn(1,2),respmlbn(1,2),cb,ca,dt,npts)
      call rk4_upd(pmlbn(1,3),kpmlbn(1,3),respmlbn(1,3),cb,ca,dt,npts)
      call rk4_upd(pmldn(1,1),kpmldn(1,1),respmldn(1,1),cb,ca,dt,npts)
      call rk4_upd(pmldn(1,2),kpmldn(1,2),respmldn(1,2),cb,ca,dt,npts)
      call rk4_upd(pmldn(1,3),kpmldn(1,3),respmldn(1,3),cb,ca,dt,npts)

      endif

      stime0=dclock()-stime0

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
          write(6,*) 'rk_maxwell_ab                (sec)::',stime0
      endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      if (nid.eq.0) write(6,*) 'compute eigenvalues'

      call cem_maxwell_eig        !call schrod_eig_arpack
      call cem_end

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_eig
C----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer       lpts_eig2
      parameter    (lpts_eig2=lpts_eig*lpts_eig)

      common /EIGC/ res(  lpts_eig2)
     $            , ev (  lpts_eig2)
     $            , er (  lpts_eig)
     $            , ei (  lpts_eig)
     $            , rew(4*lpts_eig)
      real          res,ev,er,ei,rew
      integer       info,i,j,jj,k,k1,nn,nn2,npts2,npts3,npts4,npts6


      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real           srflx

      nfield=3*(ndim-1)

      npts2= npts*2               
      npts3= npts*3               
      npts4= npts*4               
      npts6= npts*6               

      if (nfield*npts.gt.lpts_eig) then
         write(6,*) 'lpts_eig',lpts_eig, nfield*npts 
         write(6,*) 'insufficient memory in cem_maxwell_eig, abort'    
         call exitt
      endif

      if (nid.eq.0) write(6,*) 'start cem_maxwell_eig',imode
      nn  = npts*nfield
      nn2 = nn*nn
      call rzero(res,nn2)

c... if3d
      if (imode.eq.3) then

      k1 = 0
      do jj= 1,nfield
      do j = 1,npts

         k1 =  k1+1
         k  = (k1-1)*nn+1

         call rzero(HN,npts3)
         call rzero(EN,npts3)

         if (jj.eq.1) HN(j,1)= 1.  ! jth col of identity matrix
         if (jj.eq.2) HN(j,2)= 1.  ! jth col of identity matrix
         if (jj.eq.3) HN(j,3)= 1.  ! jth col of identity matrix
         if (jj.eq.4) EN(j,1)= 1.  ! jth col of identity matrix
         if (jj.eq.5) EN(j,2)= 1.  ! jth col of identity matrix
         if (jj.eq.6) EN(j,3)= 1.  ! jth col of identity matrix

         call cem_maxwell
         call cem_maxwell_restrict_to_face
         call cem_maxwell_flux (srflx)
         call cem_maxwell_add_flux_to_res (srflx)
        !call cem_maxwell_invqmass
         
         call copy(res(k+0*npts),resHN(1,1),npts)
         call copy(res(k+1*npts),resHN(1,2),npts)
         call copy(res(k+2*npts),resHN(1,3),npts)
         call copy(res(k+3*npts),resEN(1,1),npts)
         call copy(res(k+4*npts),resEN(1,2),npts)
         call copy(res(k+5*npts),resEN(1,3),npts)

      enddo
      enddo

      elseif (imode.eq.2) then

      k1 = 0
      do jj = 1,nfield
      do j  = 1,npts

         k1 =  k1+1
         k  = (k1-1)*nn+1

         call rzero(HN,npts3)
         call rzero(EN,npts3)

         if (jj.eq.1) HN(j,1)= 1.  ! jth col of identity matrix
         if (jj.eq.2) HN(j,2)= 1.  ! jth col of identity matrix
         if (jj.eq.3) EN(j,3)= 1.  ! jth col of identity matrix

         call cem_maxwell
         call cem_maxwell_restrict_to_face
         call cem_maxwell_flux (srflx)
         call cem_maxwell_add_flux_to_res (srflx)
        !call cem_maxwell_invqmass

         call copy(res(k+0*npts),resHN(1,1),npts)
         call copy(res(k+1*npts),resHN(1,2),npts)
         call copy(res(k+2*npts),resEN(1,3),npts)

      enddo
      enddo

      elseif (imode.eq.1) then

      k1 = 0
      do jj = 1,nfield
      do j  = 1,npts

         k1 =  k1+1
         k  = (k1-1)*nn+1

         call rzero(EN(1,1),npts)
         call rzero(EN(1,2),npts)
         call rzero(HN(1,3),npts)

         if (jj.eq.1) EN(j,1)= 1.  ! jth col of identity matrix
         if (jj.eq.2) EN(j,2)= 1.  ! jth col of identity matrix
         if (jj.eq.3) HN(j,3)= 1.  ! jth col of identity matrix

         call cem_maxwell
         call cem_maxwell_restrict_to_face
         call cem_maxwell_flux (srflx)
         call cem_maxwell_add_flux_to_res (srflx)
        !call cem_maxwell_invqmass

         call copy(res(k+0*npts),resEN(1,1),npts)
         call copy(res(k+1*npts),resEN(1,2),npts)
         call copy(res(k+2*npts),resHN(1,3),npts)
        
      enddo
      enddo

      endif

      if (nid.eq.0) write(6,*) 'obtain spatial operator'

      do i=1,nn2  
         if (abs(res(i)).lt.1.e-13) res(i)=0.0
      enddo

      open(50,file='spatial_matrix.dat')

      do i=1,nn2
         write(50,*) res(i)
      enddo
      close(50)


      if (nid.eq.0) then
          write(6,*) '                                      '
          write(6,*) 'IFEIG T :: computing eigenvalues      '
          write(6,*) '                                      '
      endif

      call dgeev('N','V',nn,res,nn,er,ei,ev,nn,ev,nn,rew,4*nn,info)

      open(60,file='eigenvalue.dat')
      do i=1,nn
        write(60,*)  er(i),ei(i)
      enddo
      close(60)

    1 format(1p1e19.12)

      if (nid.eq.0) write(6,*) 'eigenvalue computation completed!!!'
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_exp
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      integer  icalld
      save     icalld
      data     icalld/0/
      real     glsc2,beta,tmp
      integer  i,j,k,n6,n,m,info

      n6 = npts*(3+3*(ndim-2))
      n  = npts
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

      if     (imode.eq.3) then
      call copy(ru(0*npts+1),hn(1,1),npts)
      call copy(ru(1*npts+1),hn(1,2),npts)
      call copy(ru(2*npts+1),hn(1,3),npts)
      call copy(ru(3*npts+1),en(1,1),npts)
      call copy(ru(4*npts+1),en(1,2),npts)
      call copy(ru(5*npts+1),en(1,3),npts)
      elseif (imode.eq.2) then
      call copy(ru(0*npts+1),hn(1,1),npts)
      call copy(ru(1*npts+1),hn(1,2),npts)
      call copy(ru(2*npts+1),en(1,3),npts)
      elseif (imode.eq.1) then
      call copy(ru(0*npts+1),en(1,1),npts)
      call copy(ru(1*npts+1),en(1,2),npts)
      call copy(ru(2*npts+1),hn(1,3),npts)
      endif

      beta= glsc2(ru,ru,n6)
      beta= sqrt (beta)

      call arnoldi_vec (rv0,rh0,n6,m,ru)
      call arnoldi_mat (rvt,rv,rh,rv0,rh0,m,n6)          

      call copy_r2z    (zh,rh,m*m)
      call zgeev  ('N','V',m,zh,m,eval,vri,m,vr,m,awork,2*m,rwork,info)
      call copy_z2z    (vri,vr,m*m)

      call zgetrf      (m,m,vri,m,ipiv,info)      !vri= LU factorization                
      call zgetri      (m,vri,m,ipiv,awork,m,info)!vri= X^(-1)

      call z_exp       (eval,dt,m)      !eval=exp(dt*eval)
      call z_rowscale  (vri,eval,m)     !vri =eval*X^(-1); 
      call c_mxm       (vr,m,vri,m,zh,m)         !zh =X*vri  

      call copy_z2r    (rw,zh(1,1),m)                 !zh(i,1) (i=1:m)
      call mxm         (rv,n6,rw,m,ru,1)      
      call cmult       (ru,beta,n6)      

      if     (imode.eq.3) then
      call copy(hn(1,1),ru(0*npts+1),npts)
      call copy(hn(1,2),ru(1*npts+1),npts)
      call copy(hn(1,3),ru(2*npts+1),npts)
      call copy(en(1,1),ru(3*npts+1),npts)
      call copy(en(1,2),ru(4*npts+1),npts)
      call copy(en(1,3),ru(5*npts+1),npts)
      elseif (imode.eq.2) then
      call copy(hn(1,1),ru(0*npts+1),npts)
      call copy(hn(1,2),ru(1*npts+1),npts)
      call copy(en(1,3),ru(2*npts+1),npts)
      elseif (imode.eq.1) then
      call copy(en(1,1),ru(0*npts+1),npts)
      call copy(en(1,2),ru(1*npts+1),npts)
      call copy(hn(1,3),ru(2*npts+1),npts)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine arnoldi_mat (vt,v,h,v0,h0,m,n)
      implicit none
      integer  m,n,i,j,k
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
      implicit none
      integer  i,j,m,n
      real     v(n,m+1),b(n)
      real     h(m+1,m)         
      real     eps, s1, si, glsc2, cnt

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
      subroutine amult(w,v,n6)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer  n6
      real     w(n6),v(n6)
      common  /scrns/ srflx(6*2*ldim*lx1*lz1*lelt)
      real            srflx


      if        (imode.eq.3) then
           call copy(hn(1,1), v(1+0*npts), npts)
           call copy(hn(1,2), v(1+1*npts), npts)
           call copy(hn(1,3), v(1+2*npts), npts)
           call copy(en(1,1), v(1+3*npts), npts)
           call copy(en(1,2), v(1+4*npts), npts)
           call copy(en(1,3), v(1+5*npts), npts)
      elseif    (imode.eq.2) then
           call copy(hn(1,1), v(1+0*npts), npts)
           call copy(hn(1,2), v(1+1*npts), npts)
           call copy(en(1,3), v(1+2*npts), npts)
      elseif    (imode.eq.1) then
           call copy(en(1,1), v(1+0*npts), npts)
           call copy(en(1,2), v(1+1*npts), npts)
           call copy(hn(1,3), v(1+2*npts), npts)
      endif

      call cem_maxwell
      call cem_maxwell_restrict_to_face
      call cem_maxwell_flux (srflx)
      call cem_maxwell_add_flux_to_res (srflx)
      call cem_maxwell_apps
      call cem_maxwell_source
      call cem_maxwell_invqmass

      if        (imode.eq.3) then
           call copy(w(1+0*npts), reshn(1,1), npts)
           call copy(w(1+1*npts), reshn(1,2), npts)
           call copy(w(1+2*npts), reshn(1,3), npts)
           call copy(w(1+3*npts), resen(1,1), npts)
           call copy(w(1+4*npts), resen(1,2), npts)
           call copy(w(1+5*npts), resen(1,3), npts)
      elseif    (imode.eq.2) then
           call copy(w(1+0*npts), reshn(1,1), npts)
           call copy(w(1+1*npts), reshn(1,2), npts)
           call copy(w(1+2*npts), resen(1,3), npts)
      elseif    (imode.eq.1) then
           call copy(w(1+0*npts), resen(1,1), npts)
           call copy(w(1+1*npts), resen(1,2), npts)
           call copy(w(1+2*npts), reshn(1,3), npts)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_r2z(z,r,n)
      implicit none 
      integer    n,i
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         z(i)=r(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2r(r,z,n)
      implicit none 
      integer    n,i
      complex*16 z(n)
      real*8     r(n)

      do i=1,n
         r(i)=real(z(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_z2z(z2,z1,n)
      implicit none 
      integer    n,i
      complex*16 z1(n),z2(n)

      do i=1,n
         z2(i)=z1(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c_mxm(a,n1,b,n2,c,n3) ! complex mxm
      implicit none 
      integer    n1,n2,n3,i
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
      subroutine z_rowscale(a,s,n)
      implicit none
      integer    n,i,j
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
      implicit none
      integer    n,i
      real*8     dt
      complex*16 f(n),a

      do i=1,n
         a    = dt*f(i)
         f(i) = exp(a)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine z_rowscale2(a,b,c,dt,n)
      implicit none
      integer    n,i
      complex*16 a(n),b(n),c(n)
      real       dt

      do i=1,n
         a(i)=b(i)*exp(c(i)*dt)
      enddo

      return
      end 
c-----------------------------------------------------------------------
