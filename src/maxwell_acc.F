C========  OpenACC

#ifdef GPU

#ifdef _OPENACC
c-----------------------------------------------------------------------
      subroutine rk_c_acc(i)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      integer  i

      RKtime=time+dt*rk4c(i)


      return
      end
c-----------------------------------------------------------------------
      subroutine rk_storage_acc
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RK5'

      if     (IFRK45) then

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

      elseif (IFRK22) then

      rk4a(1) =   0.0 
      rk4a(2) =  -1.0                         

      rk4b(1) =   1.0 
      rk4b(2) =   1.0/2.0                       

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_rk_acc                   
c--------------------------------------------------------------------- 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'RTIMER'
      integer  ii

      real*8  dclock
      real    stime0,stime5
c ... 5-stage 4th-order RK
      stime5=dclock()
      do ii=1,5             

         stime0=dclock()
         call rk_c_acc (ii)
         call cem_maxwell_op_acc   
         call rk_maxwell_ab_acc(ii)
         stime0=dclock()-stime0

      enddo
      stime5= dclock()-stime5
      acc_t = acc_t   +stime5     

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
          write(6,*) '-------------------------------------::'
          write(6,*) 'cem_maxwell_op_rk_acc           (sec)::',stime0
          write(6,*) 'cem_maxwell_op_rk_acc per istep (sec)::',stime5
          write(6,*) '-------------------------------------::'
      endif
      endif
#endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_op_acc
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

      call cem_maxwell_acc
      stime1=dclock()-stime1       

      stime2=dclock()
      call cem_maxwell_restrict_to_face_acc
      stime2=dclock()-stime2       

      stime3=dclock()
      call cem_maxwell_flux_acc (srflx)
      stime3=dclock()-stime3       

      stime4=dclock()

      call cem_maxwell_add_flux_to_res_acc (srflx)
      stime4=dclock()-stime4       

      stime5=dclock()
      call cem_maxwell_apps_acc
      stime5=dclock()-stime5       

      stime6=dclock()
      call cem_maxwell_source_acc
      stime6=dclock()-stime6       

      stime7=dclock()
      call cem_maxwell_invqmass_acc

      stime7=dclock()-stime7       
      stime0=stime1+stime2+stime3+stime4+stime5+stime6+stime7

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then 
          write(6,*) '-------------------------------------::'
          write(6,*) 'cem_maxwell_acc                 (sec)::',stime1 
          write(6,*) 'cem_maxwell_restrict_to_face_acc(sec)::',stime2 
          write(6,*) 'cem_maxwell_flux_acc            (sec)::',stime3 
          write(6,*) 'cem_maxwell_add_flux_to_res_acc (sec)::',stime4 
          write(6,*) 'cem_maxwell_apps_acc            (sec)::',stime5 
          write(6,*) 'cem_maxwell_source_acc          (sec)::',stime6 
          write(6,*) 'cem_maxwell_invqmass_acc        (sec)::',stime7 
          write(6,*) 'cem_maxwell_op_acc              (sec)::',stime0 
          write(6,*) '-------------------------------------::'
      endif
      endif
#endif
      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_acc 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

c...  sacrify 2D optimization for simplicity::
      if (ifdealias) then

       write(6,*) "OpenACC is not implemented for ifdealias=true"
       call exitt

       if     (imode.eq.3) then !IF3D

       call maxwell_wght_dcurl(resEN(1,1),resEN(1,2),resEN(1,3),HN(1,1),
     $                                                 HN(1,2),HN(1,3))
       call maxwell_wght_dcurl(resHN(1,1),resHN(1,2),resHN(1,3),EN(1,1),
     $                                                 EN(1,2),EN(1,3))
       call chsign (resHN(1,1),npts)
       call chsign (resHN(1,2),npts)
       call chsign (resHN(1,3),npts)

       elseif (imode.eq.2) then !IFTE

       call maxwell_wght_dcurl(resHN(1,1),resHN(1,2),resEN(1,3),HN(1,1),
     $                                                  HN(1,2),EN(1,3))
       call chsign (resHN(1,1),npts)
       call chsign (resHN(1,2),npts)

       elseif (imode.eq.1) then !IFTE
       call maxwell_wght_dcurl(resEN(1,1),resEN(1,2),resHN(1,3),EN(1,1),
     $                                                  EN(1,2),HN(1,3))
       call chsign (resHN(1,3),npts)

       endif

      else

!$ACC DATA PRESENT(resEN,resHN,HN,EN)
       if     (imode.eq.3) then !IF3D

       call maxwell_wght_curl_acc (resEN(1,1),resEN(1,2),resEN(1,3)
     $         ,                      HN(1,1),   HN(1,2),   HN(1,3))

       call maxwell_wght_curl_acc (resHN(1,1),resHN(1,2),resHN(1,3)
     $      ,                         EN(1,1),   EN(1,2),   EN(1,3))
       call chsign_acc (resHN(1,1),npts)
       call chsign_acc (resHN(1,2),npts)
       call chsign_acc (resHN(1,3),npts)

       elseif (imode.eq.2) then !IFTE

       call maxwell_wght_curl_acc (resHN(1,1),resHN(1,2),resEN(1,3)
     $         ,                      HN(1,1),   HN(1,2),   EN(1,3))
       call chsign_acc (resHN(1,1),npts)
       call chsign_acc (resHN(1,2),npts)

       elseif (imode.eq.1) then !IFTE
       call maxwell_wght_curl_acc (resEN(1,1),resEN(1,2),resHN(1,3)
     $         ,                      EN(1,1),   EN(1,2),   HN(1,3))
       call chsign_acc (resHN(1,3),npts)
       endif

!$ACC END DATA

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_restrict_to_face_acc
c-----------------------------------------------------------------------
c     restrict u to faces
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

c...  sacrify efficiency for 2D for simplicity
!$ACC DATA PRESENT (fHN,fEN,HN,EN,cemface)
!$ACC PARALLEL LOOP
      do j=1, ncemface
           i= cemface(j)
           fHN(j,1) = HN(i,1)
           fHN(j,2) = HN(i,2)
           fHN(j,3) = HN(i,3)
           fEN(j,1) = EN(i,1)
           fEN(j,2) = EN(i,2)
           fEN(j,3) = EN(i,3)
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux_acc(srflx)
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'EMWAVE'

      real     srflx(1)

      if (if3d) then
         call cem_maxwell_flux3d_acc (srflx)
      else
         call cem_maxwell_flux2d_acc (srflx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux2d_acc (srflx)  
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

!$ACC  DATA PRESENT(srflx,cemface,unxm,unym,fEN,fHN)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)
      if     (imode.eq.2) then !IFTM

        k = nxzfl
!$ACC PARALLEL LOOP
        do i=1,nxzfl                
           j=cemface(i)
           srflx(0*k+i)=-unym(i)*fEN(i,3)  ! NxfHx
           srflx(1*k+i)= unxm(i)*fEN(i,3)  ! NxfHy
           srflx(2*k+i)=-unxm(i)*fHN(i,2)+unym(i)*fHN(i,1) ! NxfEz
        enddo
!$ACC END PARALLEL LOOP

        t0=dclock()
        call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)

!$ACC WAIT
        call measure_comm_acc(t0)

        if (IFPML.or.IFPEC) then
           write(*,*) "OpenACC is not implemented for (IFPML.or.IFPEC)"
           stop
           call cem_maxwell_flux_pec(srflx)       
        endif

        k = nxzfl
!$ACC PARALLEL LOOP
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
!ACC END PARALLEL LOOP

      elseif (imode.eq.1) then !IFTE

        k = nxzfl
!$ACC PARALLEL LOOP
        do i=1,nxzfl                
         srflx(0*k+i)= -unym(i)*fHN(i,3)  ! NxfHx   
         srflx(1*k+i)=  unxm(i)*fHN(i,3)  ! NxfHy   
         srflx(2*k+i)= -unxm(i)*fEN(i,2)+unym(i)*fEN(i,1) ! NxfEz
        enddo
!$ACC END PARALLEL LOOP

        t0 = dclock()
        call gs_op_fields_acc(gsh_face,srflx,nxzfl,3,1,1,0)

!$ACC WAIT
        call measure_comm_acc(t0)

        if (IFPML.or.IFPEC) then
           write(*,*) "OpenACC is not implemented for (IFPML.or.IFPEC)"
           stop
           call cem_maxwell_flux_pec(srflx)        
        endif

        k = nxzfl
!$ACC PARALLEL LOOP
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
!$ACC END PARALLEL LOOP 
      endif
!$ACC END DATA

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_maxwell_flux3d_acc(srflx)
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
      real    srflx(6*2*ldim*lx1*lz1*lelt)
      real    C0,y0,y1,z0,z1,y02,z02,c02y,c02z,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k

      if (ifcentral) C0 = 0.0
      if (ifupwind ) C0 = 1.0

      !srflx(1)-> reshx
      !srflx(4)-> resex

!$ACC  DATA PRESENT(srflx,fEN,fHN,unxm,unym,unzm)
!$ACC& PRESENT(Y_0,Y_1,Z_0,Z_1)

      k = nxzfl

!$ACC PARALLEL LOOP
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
!$ACC END PARALLEL LOOP

      t0= dclock()

      call gs_op_fields_acc(gsh_face,srflx,nxzfl,6,1,1,0)

!$ACC WAIT
      call measure_comm_acc(t0)

      !srflx now contains n x [E] and n x [H],with [F] = F^+ - F^-;(observe un^+ = - un^-)
      if (IFPEC.or.IFPML) then
         write(*,*) "OpenACC is not implemented for (IFPEC.or.IFPML)"
         stop
         call cem_maxwell_flux_pec(srflx)            
      endif

      k = nxzfl

!$ACC PARALLEL LOOP
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
!$ACC END PARALLEL LOOP

!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_add_flux_to_res_acc(srflx)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      integer  i,j,k,p,l,nptr
      real    srflx(6*2*ldim*lx1*lz1*lelt),a

      COMMON /cemfce_acc/ cemface2(2*ldim*lx1*lz1*lelt)
     $     ,              cemface_ptr(0:2*ldim*lx1*lz1*lelt)
     $     ,              cemface_ind(2*ldim*lx1*lz1*lelt)
      integer cemface2, cemface_ptr, cemface_ind

      k = nxzfl               
      nptr = cemface_ptr(0)

!$ACC DATA PRESENT(resHN,resEN,aream,srflx,cemface2,cemface_ptr)
      if     (imode.eq.3) then !IF3D
!$ACC PARALLEL LOOP 
        do p = 1,nptr
!$ACC LOOP SEQ
           do l = cemface_ptr(p),cemface_ptr(p+1)-1
              i = cemface_ind(l)
              j = cemface2(l)
              a = aream(j)
              resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
              resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
              resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
              resEN(i,1) = resEN(i,1) + a*srflx(3*k+j)
              resEN(i,2) = resEN(i,2) + a*srflx(4*k+j)
              resEN(i,3) = resEN(i,3) + a*srflx(5*k+j)
           enddo
        enddo

      elseif (imode.eq.2) then !IFTM
!!$ACC PARALLEL LOOP
        do j=  1,ncemface
           i = cemface(j)
           a = aream(j)
           resHN(i,1) = resHN(i,1) + a*srflx(0*k+j)
           resHN(i,2) = resHN(i,2) + a*srflx(1*k+j)
           resEN(i,3) = resEN(i,3) + a*srflx(2*k+j)
        enddo
!!$ACC END PARALLEL LOOP
      elseif (imode.eq.1) then !IFTE
!!$ACC PARALLEL LOOP
        do j = 1,ncemface
           i = cemface(j)
           a = aream(j)
           resEN(i,1) = resEN(i,1) + a*srflx(0*k+j)
           resEN(i,2) = resEN(i,2) + a*srflx(1*k+j)
           resHN(i,3) = resHN(i,3) + a*srflx(2*k+j)
        enddo
!!$ACC END PARALLEL LOOP
      endif
!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_apps_acc 
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'

c     if (ifdrude) call cem_maxwell_drude_acc

      if (ifloren .or.ifpml) then
         write(*,*) "OpenACC is not implemented for (ifloren .or.ifpml)"
         stop
      endif

c     if (ifloren) call cem_maxwell_lorentz
c     if (ifhydro) call cem_maxwell_hydraulic    
c     if (ifpml  ) call pml_step

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_drude_acc 
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      integer  i
      real     bmv
      real     ca1, cb1, cc1
   
!$ACC DATA PRESENT(resEN,JN1,resJN1,EN,resPN1,drudea,drudeb,bmn)

!$ACC PARALLEL LOOP
      do i=1,npts

         bmv = bmn(i)
         ca1 = drudea(i)
         cb1 = drudeb(i)
      
         resEN(i,1)= resEN(i,1)-JN1(i,1)*bmv
         resEN(i,2)= resEN(i,2)-JN1(i,2)*bmv
         resEN(i,3)= resEN(i,3)-JN1(i,3)*bmv

         resJN1(i,1)= ca1*JN1(i,1) + cb1*EN(i,1)                
         resJN1(i,2)= ca1*JN1(i,2) + cb1*EN(i,2)              
         resJN1(i,3)= ca1*JN1(i,3) + cb1*EN(i,3)                  
 
         resPN1(i,1)= JN1(i,1)
         resPN1(i,2)= JN1(i,2)                      
         resPN1(i,3)= JN1(i,3)                    
      enddo
!$ACC END PARALLEL LOOP

!$ACC END DATA
      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_invqmass_acc
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'PML'
      integer  i

!$ACC DATA PRESENT(resHN,resEN,hbm1,ebm1)
      if     (imode.eq.3) then !IF3D
!$ACC PARALLEL LOOP
          do i=1,npts
             resHN(i,1)=resHN(i,1)*hbm1(i)
             resHN(i,2)=resHN(i,2)*hbm1(i)
             resHN(i,3)=resHN(i,3)*hbm1(i)
             resEN(i,1)=resEN(i,1)*ebm1(i)
             resEN(i,2)=resEN(i,2)*ebm1(i)
             resEN(i,3)=resEN(i,3)*ebm1(i)
          enddo
!$ACC END PARALLEL LOOP
      elseif (imode.eq.2) then !IFTM
!$ACC PARALLEL LOOP  
          do i=1,npts
             resHN(i,1)=reshN(i,1)*hbm1(i) 
             resHN(i,2)=reshN(i,2)*hbm1(i)
             resEN(i,3)=reseN(i,3)*ebm1(i)
          enddo
!$ACC END PARALLEL LOOP  
      elseif (imode.eq.1) then !IFTE
!$ACC PARALLEL LOOP
           do i=1,npts
             resEN(i,1)=resEN(i,1)*ebm1(i)
             resEN(i,2)=resEN(i,2)*ebm1(i)
             resHN(i,3)=resHN(i,3)*hbm1(i)
           enddo
!$ACC END PARALLEL LOOP 
      endif

      if  (ifpml) then
         write(*,*) "OpenACC is not implemented for ifpml"
         stop
           do i=1,npts
             respmlBN(i,1)=respmlBN(i,1)*bminv(i)   
             respmlBN(i,2)=respmlBN(i,2)*bminv(i)  
             respmlBN(i,3)=respmlBN(i,3)*bminv(i)   
             respmlDN(i,1)=respmlDN(i,1)*bminv(i)   
             respmlDN(i,2)=respmlDN(i,2)*bminv(i)  
             respmlDN(i,3)=respmlDN(i,3)*bminv(i)   
           enddo

      endif
!$ACC END DATA

      return
      end

c-----------------------------------------------------------------------
      subroutine rk_maxwell_ab_acc(ii)
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

      call rk4_upd_acc(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)
      call rk4_upd_acc(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)

      elseif (imode.eq.2) then

      call rk4_upd_acc(HN(1,1),kHN(1,1),resHN(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(HN(1,2),kHN(1,2),resHN(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(EN(1,3),kEN(1,3),resEN(1,3),cb,ca,dt,npts)

      elseif (imode.eq.1) then

      call rk4_upd_acc(EN(1,1),kEN(1,1),resEN(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(EN(1,2),kEN(1,2),resEN(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(HN(1,3),kHN(1,3),resHN(1,3),cb,ca,dt,npts)

      endif

      if (ifdrude) then
      call rk4_upd_acc(JN1(1,1),kJN1(1,1),resJN1(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(JN1(1,2),kJN1(1,2),resJN1(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(JN1(1,3),kJN1(1,3),resJN1(1,3),cb,ca,dt,npts)
      call rk4_upd_acc(PN1(1,1),kPN1(1,1),resPN1(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(PN1(1,2),kPN1(1,2),resPN1(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(PN1(1,3),kPN1(1,3),resPN1(1,3),cb,ca,dt,npts)
      endif

      if (ifloren) then
      call rk4_upd_acc(JN1(1,1),kJN1(1,1),resJN1(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(JN1(1,2),kJN1(1,2),resJN1(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(JN1(1,3),kJN1(1,3),resJN1(1,3),cb,ca,dt,npts)
      call rk4_upd_acc(JN2(1,1),kJN2(1,1),resJN2(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(JN2(1,2),kJN2(1,2),resJN2(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(JN2(1,3),kJN2(1,3),resJN2(1,3),cb,ca,dt,npts)
      call rk4_upd_acc(JN3(1,1),kJN3(1,1),resJN3(1,1),cb,ca,dt,npts)
      call rk4_upd_acc(JN3(1,2),kJN3(1,2),resJN3(1,2),cb,ca,dt,npts)
      call rk4_upd_acc(JN3(1,3),kJN3(1,3),resJN3(1,3),cb,ca,dt,npts)
      endif

      if (ifpml) then

         call rk4_upd_acc(pmlbn(1,1),kpmlbn(1,1),respmlbn(1,1)
     $        ,           cb,ca,dt,npts)
         call rk4_upd_acc(pmlbn(1,2),kpmlbn(1,2),respmlbn(1,2)
     $        ,           cb,ca,dt,npts)
         call rk4_upd_acc(pmlbn(1,3),kpmlbn(1,3),respmlbn(1,3)
     $        ,           cb,ca,dt,npts)
         call rk4_upd_acc(pmldn(1,1),kpmldn(1,1),respmldn(1,1)
     $        ,           cb,ca,dt,npts)
         call rk4_upd_acc(pmldn(1,2),kpmldn(1,2),respmldn(1,2)
     $        ,           cb,ca,dt,npts)
         call rk4_upd_acc(pmldn(1,3),kpmldn(1,3),respmldn(1,3)
     $        ,           cb,ca,dt,npts)

      endif

      stime0=dclock()-stime0

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
      if (nid.eq.0) then
          write(6,*) 'rk_maxwell_ab_acc            (sec)::',stime0
      endif
      endif
#endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_maxwell_source_acc
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      if (.not.ifsrc) return

      if (ifsrc) then
         write(*,*) "OpenACC is not implemented for ifsrc=.true."
         stop
      endif

      call usersrc (57,reshn(1,1),reshn(1,2),reshn(1,3),
     $                 resen(1,1),resen(1,2),resen(1,3))

      return
      end

c-----------------------------------------------------------------------
      subroutine local_grad3_acc ( u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,
     $                          u1,u2,u3,d,dtranse,nn,nelti,lpts1i,nidi)
c-----------------------------------------------------------------------
      include 'SIZE'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt) 
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt) 
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)

      real u1 (lx1,ly1,lz1,lelt) 
      real u2 (lx1,ly1,lz1,lelt)
      real u3 (lx1,ly1,lz1,lelt)

      real d(lx1,ly1), dtrans(lx,ly1)

      real tmpr1,tmps1,tmpt1
      real tmpr2,tmps2,tmpt2
      real tmpr3,tmps3,tmpt3
      integer i,j,k,l,e 

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)
!$ACC&     PRESENT(u1,u2,u3,d,dtrans)

!$ACC KERNELS ASYNC
!$ACC LOOP GANG
      do e = 1,nelt
!$ACC LOOP VECTOR(LZ1) 
         do k = 1,nz1
!$ACC LOOP VECTOR(LY1) 
         do j = 1,ny1
!$ACC LOOP VECTOR(LX1) 
         do i = 1,nx1
            tmpr1 = 0.0
            tmpr2 = 0.0
            tmpr3 = 0.0
            do l=1,nx1
               tmpr1=tmpr1+d(i,l)*u1(l,j,k,e)
               tmpr2=tmpr2+d(i,l)*u2(l,j,k,e)
               tmpr3=tmpr3+d(i,l)*u3(l,j,k,e)
            enddo
            u1r(i,j,k,e) = tmpr1
            u2r(i,j,k,e) = tmpr2
            u3r(i,j,k,e) = tmpr3
         enddo
         enddo
         enddo
      enddo
!$ACC END KERNELS

!$ACC KERNELS ASYNC
!$ACC LOOP GANG   
      do e = 1,nelt
!$ACC LOOP VECTOR(LZ1)
         do k = 1,nz1
!$ACC LOOP VECTOR(LY1)
         do j = 1,ny1
!$ACC LOOP VECTOR(LX1)
         do i = 1,nx1
            tmps1 = 0.0
            tmps2 = 0.0
            tmps3 = 0.0
            do l=1,nx1
               tmps1=tmps1+d(j,l)*u1(i,l,k,e)
               tmps2=tmps2+d(j,l)*u2(i,l,k,e)
               tmps3=tmps3+d(j,l)*u3(i,l,k,e)
            enddo
            u1s(i,j,k,e) = tmps1
            u2s(i,j,k,e) = tmps2
            u3s(i,j,k,e) = tmps3
         enddo
         enddo
         enddo
      enddo
!$ACC END KERNELS

!$ACC KERNELS ASYNC
!$ACC LOOP GANG   
      do e = 1,nelt
!$ACC LOOP VECTOR(LZ1) 
         do k = 1,nz1
!$ACC LOOP VECTOR(LY1) 
         do j = 1,ny1
!$ACC LOOP VECTOR(LX1) 
         do i = 1,nx1
            tmpt1 = 0.0
            tmpt2 = 0.0
            tmpt3 = 0.0
            do l=1,nx1
               tmpt1=tmpt1+d(k,l)*u1(i,j,l,e)
               tmpt2=tmpt2+d(k,l)*u2(i,j,l,e)
               tmpt3=tmpt3+d(k,l)*u3(i,j,l,e)
            enddo
            u1t(i,j,k,e) = tmpt1
            u2t(i,j,k,e) = tmpt2
            u3t(i,j,k,e) = tmpt3
         enddo
         enddo
         enddo
      enddo
!$ACC END KERNELS

!$ACC WAIT
!$ACC END DATA

      return
      end

c----------------------------------------------------------------------
      subroutine curl_acc(u1r, u1s, u1t, u2r, u2s, u2t, u3r, u3s, u3t,
     $                    rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn,
     $                    w1,  w2,  w3,  w3mn,nxyzi,nelti,lpts1i)
c-----------------------------------------------------------------------
      implicit none

      include 'SIZE'

      real u1r(lx1,ly1,lz1,lelt)
      real u1s(lx1,ly1,lz1,lelt)
      real u1t(lx1,ly1,lz1,lelt)
      real u2r(lx1,ly1,lz1,lelt)
      real u2s(lx1,ly1,lz1,lelt)
      real u2t(lx1,ly1,lz1,lelt)
      real u3r(lx1,ly1,lz1,lelt)
      real u3s(lx1,ly1,lz1,lelt)
      real u3t(lx1,ly1,lz1,lelt)
      real w1 (lx1,ly1,lz1,lelt)
      real w2 (lx1,ly1,lz1,lelt)
      real w3 (lx1,ly1,lz1,lelt)

      real w3mn(lx1,ly1,lz1)
      real rxmn(lx1,ly1,lz1,lelt)
      real sxmn(lx1,ly1,lz1,lelt)
      real txmn(lx1,ly1,lz1,lelt)
      real rymn(lx1,ly1,lz1,lelt)
      real symn(lx1,ly1,lz1,lelt)
      real tymn(lx1,ly1,lz1,lelt)
      real rzmn(lx1,ly1,lz1,lelt)
      real szmn(lx1,ly1,lz1,lelt)
      real tzmn(lx1,ly1,lz1,lelt)

      real tmp
      integer nxyzi,nelti,lpts1i
      integer i,j,k,l,e

!$ACC DATA PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,w1,w2,w3)
!$ACC&     PRESENT(w3mn,rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn)

!$ACC KERNELS     
!$ACC LOOP GANG   
      do e = 1,nelt
!$ACC LOOP VECTOR(LZ1)
         do k = 1,nz1
!$ACC LOOP VECTOR(LY1)
         do j = 1,ny1
!$ACC LOOP VECTOR(LX1)
         do i = 1,nx1
            u1r(i,j,k,e) = u1r(i,j,k,e) * w3mn(i,j,k)
            u1s(i,j,k,e) = u1s(i,j,k,e) * w3mn(i,j,k)
            u1t(i,j,k,e) = u1t(i,j,k,e) * w3mn(i,j,k)
            u2r(i,j,k,e) = u2r(i,j,k,e) * w3mn(i,j,k)
            u2s(i,j,k,e) = u2s(i,j,k,e) * w3mn(i,j,k)
            u2t(i,j,k,e) = u2t(i,j,k,e) * w3mn(i,j,k)
            u3r(i,j,k,e) = u3r(i,j,k,e) * w3mn(i,j,k)
            u3s(i,j,k,e) = u3s(i,j,k,e) * w3mn(i,j,k)
            u3t(i,j,k,e) = u3t(i,j,k,e) * w3mn(i,j,k)
         enddo
         enddo
         enddo
      enddo

!$ACC LOOP GANG   
      do e = 1,nelt
!$ACC LOOP VECTOR(LZ1)    
         do k = 1,nz1
!$ACC LOOP VECTOR(LY1)    
         do j = 1,ny1
!$ACC LOOP VECTOR(LX1)
         do i = 1,nx1
            w1(i,j,k,e)= u3r(i,j,k,e)*rymn(i,j,k,e)
     $                 + u3s(i,j,k,e)*symn(i,j,k,e)
     $                 + u3t(i,j,k,e)*tymn(i,j,k,e)
     $                 - u2r(i,j,k,e)*rzmn(i,j,k,e)
     $                 - u2s(i,j,k,e)*szmn(i,j,k,e)
     $                 - u2t(i,j,k,e)*tzmn(i,j,k,e)

            w2(i,j,k,e)= u1r(i,j,k,e)*rzmn(i,j,k,e)
     $                 + u1s(i,j,k,e)*szmn(i,j,k,e)
     $                 + u1t(i,j,k,e)*tzmn(i,j,k,e)
     $                 - u3r(i,j,k,e)*rxmn(i,j,k,e)
     $                 - u3s(i,j,k,e)*sxmn(i,j,k,e)
     $                 - u3t(i,j,k,e)*txmn(i,j,k,e)

            w3(i,j,k,e)= u2r(i,j,k,e)*rxmn(i,j,k,e)
     $                 + u2s(i,j,k,e)*sxmn(i,j,k,e)
     $                 + u2t(i,j,k,e)*txmn(i,j,k,e)
     $                 - u1r(i,j,k,e)*rymn(i,j,k,e)
     $                 - u1s(i,j,k,e)*symn(i,j,k,e)
     $                 - u1t(i,j,k,e)*tymn(i,j,k,e)
         enddo
         enddo
         enddo
      enddo
!$ACC END KERNELS  
!$ACC END DATA

      return 
      end
c------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine maxwell_wght_curl_acc(w1,w2,w3,u1,u2,u3)
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
c.............

      wght_curl_time=dclock()
      gpu_time1     =dclock()

!$ACC DATA PRESENT(w1,w2,w3,u1,u2,u3)
!$ACC&     PRESENT(w3mn,rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn)
!$ACC&     PRESENT(dxm1,dxtm1)
!$ACC&     PRESENT(u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t)

      call local_grad3_acc (u1r,u1s,u1t,u2r,u2s,u2t,u3r,u3s,u3t,
     $                      u1 ,u2 ,u3 ,dxm1,dxtm1, nx1,nelt,lpts1,nid)

      gpu_time1=dclock()-gpu_time1

      gpu_time2=dclock()

      call curl_acc(u1r, u1s, u1t, u2r, u2s, u2t, u3r, u3s, u3t,
     $              rxmn,sxmn,txmn,rymn,symn,tymn,rzmn,szmn,tzmn,
     $              w1,  w2,  w3,  w3mn,nxyz,nelt,lpts1)
!$ACC END DATA

      gpu_time2=dclock()-gpu_time2
      wght_curl_time=dclock()-wght_curl_time


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
       end

#endif 

#endif

