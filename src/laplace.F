c-----------------------------------------------------------------------
      subroutine optimize_mesh
c
c     3D Laplace solver
c
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'

      if (nid.eq.0) then
         write(6,*) '.'
         write(6,*) 'Start 2D Poisson in 3D >>>>>> '
      endif
     
      call presetup_poisson2d
      call poisson2d_init                
      call get_elt_ptr                                
      call set_grs_2d
      call set_charge_2d
      call cg_2d          

      n=nx1*ny1*nz1*nelt
      hxmah=glmax(hy,n)
      if (nid.eq.0) write(6,*) 'hxmax',hxmax
      call put_soln_2d            
      hxmax=glmax(hy,n)
      if (nid.eq.0) write(6,*) 'hxmax',hxmax

      return
      end
c-----------------------------------------------------------------------
c parallelized
      subroutine cem_set_face_ptr_poisson                      
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON2D'
      include 'SCRCT'
      include 'BEAM'
c
      integer i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2

      nxyz  = nx1*ny1*nz1
      nxz   = nx1*nz1
      nface = 2*ndim
      nxzf  = nx1*nz1*nface ! red'd mod to area, unx, etc.

c ..  display arrays for z

      if (nid.eq.0) write(6,7) nzrr0  
   7  format(' total levels +1 =',i4)
 
      if (nid.eq.0) then
      do ii=1,nzrr0
         write(6,8) ii,zrr0(ii),zposition
      enddo
   8  format(' zarry0/zposition::',i4,2f15.8)
 
      write(6,*) '.'
      write(6,9) face,zposition,zcut
   9  format (' face=',i4,': zposition=',f10.5,': zcut=',f10.5)
      write(6,*) '.'
      endif

      tol = 1.e-12
c ..  get global numbering for 2d slice
      k1=0
      do ie= 1,nelt         
         ieg= lglel(ie,node)
         zmax= vlmax(zm1(1,1,1,ie),nxyz)
      do f = face,face    ! face =5
      do j = 1,nxz 
         ef= eface(f)
         k = j+nxz*(ef-1)+nxzf*(ie-1)
         j1= cemface(k)                        
         xx= xm1(j1,1,1,1)
         yy= ym1(j1,1,1,1)
         zz= zm1(j1,1,1,1)

         if (zmax.le.zcut+tol) then      ! miniumn z elts=2 < zcut   ?
         if ((abs(zz-zposition)).le.1.e-12) then
            k1=k1+1
            glo_num2d(k1)=glo_num3d_face(k)  ! global number
            !if (nid.eq.0) write(111,*) 'glo-',node,glo_num2d(k1)
            !if (nid.eq.1) write(112,*) 'glo-',node,glo_num2d(k1)
            !if (nid.eq.2) write(113,*) 'glo-',node,glo_num2d(k1)
            !if (nid.eq.3) write(114,*) 'glo-',node,glo_num2d(k1)
            !poissonid(k1)=j1
            !write(6,*) 'poissid',k1,j1,cemface(k)
            !write(6,*) k1,f,glo_num2d(k1),glo_num3d_face(k)
            !glo_num3d_face(k)=glo_num3d(cemface(k))
            !Extract only the global numbers on face=5 
            !of 3d mesh, where it satisfies zz=zposition
            !glo_num3d has volume information
            !glo_num3d_face has face information
         endif
         endif
      enddo
      enddo
      enddo

c ..  assign size of the numbering for 2d slice
      num2d = k1   

      return
      end
c-----------------------------------------------------------------------
      subroutine presetup_poisson2d
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      include 'SCRCT'
      include 'EMWAVE'
      include 'POISSON2D'
          
      common /scrmg/  glo_num(4*lx1*ly1*lz1*lelv)
      common /scrns/  vertex (7*lx1*ly1*lz1*lelv)
      common /work/   u(lx1*ly1*lz1*lelt)
      integer glo_num,vertex
      integer fgs_init  !check if this is ok for bgl
           
      ncrnr= 2**ndim
      mz0  = 1
      mz1  = 1

      if (if3d) mz0=1
      if (if3d) mz1=nz1

c ..  3d volume
      ntot=nx1*ny1*nz1*nelt

c ..  save glo_num3d
      call icopy (glo_num3d,glo_num,ntot) 
                           
c ..  keep one array glo_num3d_face to extract face info from glo_num3d
      call icopy (glo_num3d_face,glo_num,ntot) 
      call icopy (vertex,glo_num3d_face,ntot ) 
                               
c ..  change structure of glo_num3d_face as 'face array' of glo_num3d                           
      call iface_vert2d (glo_num3d_face,vertex,mz0,mz1,nelt)   

c ..  get 2d volume info 'glo_num2d' for 2dslice extracted 
c ..  from glo_num3d_face; get 'num2d'= size of glo_num2d           
      call cem_set_face_ptr_poisson  

c ..  2d slice volume
      ntot = num2d ! in parallel: num2d=nx1*nz1*2*ndim*nel2d         

c ..  get gop handle info through fgs_init

      gsh_hndl_2d = fgs_init(glo_num2d,ntot,np)

      if (nid.eq.0) 
     $ write(6,*) 'fgs_init: gsh_hndl_2d of glo_num2d=',gsh_hndl_2d

c ..  setup mult on 2d volume 
      call rone     (mult2d,ntot)
      call fgs_gop  (gsh_hndl_2d,mult2d,'+\0') ! dssum
      call invcol1  (mult2d,ntot)

      if (nid.eq.0)
     $    write(6,*) 'fgs_gop: dssum mult2d',ntot,np,gsh_hndl_2d

c ..  setup mask on 2d volume 
      call rone     (mask2d,ntot)
      call setmask2d(mask2d,ntot)                   
      if (nid.eq.0)
     $    write(6,*) 'setup: mask2d',ntot,np                     

      return
      end
c-----------------------------------------------------------------------
      subroutine setmask2d(smask,n) 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON2D'
      include 'SCRCT'
      include 'BEAM'

      integer   n
      real      smask(n)
      
      common   /nekcb/ cb
      character cb*3
      integer   i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2

      nxyz  = nx1*ny1*nz1
      nxz   = nx1*nz1
      nface = 2*ndim
      nxzf  = nx1*nz1*nface ! red'd mod to area, unx, etc.

c     write(6,9) face,zposition,zcut
c  9  format (' face=',i4,': zposition=',f10.5,': zcut=',f10.5)

c     call cem_set_fc_ptr  ! get cemface necessary ?

      do ie=1,nelt         
      do f =1,nface               
         ef=eface(f)
         cb =cbc(ef,ie,2)  ! ifield=2 for CEM
         !write(6,*) 'cb=',cb,', f=',f,', face=',face,', e=',ie
         if ((f.ne.5).and.(f.ne.6)) then ! face =5
             do j =1,nx1*nz1 
                ef=eface(f)
                k = j+nxz*(ef-1)+nxzf*(ie-1)
                j1=cemface(k)                        
                itmp1=glo_num3d_face(k)
                !write(6,*) 'itmp1',itmp1
                do ii=1,n           
                   itmp2=glo_num2d(ii)
                   if (itmp2.eq.itmp1) then
                      if ((cb.eq.'PEC').or.(cb.eq.'PML')) then
                        ! write(6,*) 'mask',itmp2,f,cb 
                          smask(ii)=0
                        ! search global numbers on face=1,2,3,4 (no 5,6).
                        ! thus face 5 and 6 can be pec or pml even with
                        ! 1 elt in z direction.        
                        ! glo_num matches with glo_num2d for 2d slice volume.
                        ! then check the boundary conditions.
                        ! set zero mask for the pec and pml boundaries.
                      endif
                   endif
                enddo
              enddo
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine poisson2d_init                    
      include 'SIZE'
      include 'EMWAVE'
      include 'POISSON2D'

      integer  e,i,j,l

      nel2d = num2d/(nx1*nz1)
      nelg2d= iglsum(nel2d,1)
      if (lelg2d.lt.nelg2d) then
          if (nid.eq.0) write(6,*) 'lelg2d < nelg2d',lelg2d,nelg2d    
      endif

      if (nid.eq.0) then
      write(6,*) '.'                              
      write(6,7) nelg2d,num2d            
   7  format(' nelg2d=',i5,', num2d=',i5)

      write(6,8) nelg2d*nx1*ny1
   8  format(' check (num2d.eq.nel2d*nx1*ny1)?:  nelg2d*nx1*ny1=',i8)
      endif

      if (num2d.ne.nel2d*nx1*nz1) then
          if (nid.eq.0) then
          write(6,*) 'ERR: num2d, nel2d*nx1*nz1',num2d,nel2d*nx1*nz1
          endif
          call exitt              
      endif

      n = nx1*ny1*nel2d        

      call rzero ( EX2d,n)
      call rzero ( EY2d,n)
      call rzero ( HX2d,n)
      call rzero ( HY2d,n)
      call rzero (pot2d,n)
      call rzero (chg2d,n)
      call rzero (wgt2d,n)

      call zwgll (z1d,wgt1d,nx1)

      do e=1,nel2d
         l=0
      do j=1,ny1
      do i=1,nx1
         l=l+1
         wgt2d(l,e)=wgt1d(i)*wgt1d(j)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_elt_ptr            
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON2D'
      include 'EMWAVE'
      include 'BEAM'

      integer  e,ei,e2d,ef,f,e3d,f3d,eg,el2d,elg2d
      integer  icalld
      save     icalld
      data     icalld /0/

c      call cem_set_fc_ptr ! cemface

      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7

      nxz = nx1*nz1
      nxyz= nx1*ny1*nz1
      nxzf= 2*ndim*nx1*nz1
 
      tol = 1.e-12
      k   = 0
      do ie=1,nelt         
         eg  = lglel(ie,node)
         zmax= vlmax(zm1(1,1,1,ie),nxyz)

         do i =face,face ! restrict to the case with f or ef =5
         do j =1,nx1*nz1 

         ef=eface(i)
         jj=cemface(j+nxz*(ef-1)+nxzf*(ie-1))
         xx=xm1(jj,1,1,1)
         yy=ym1(jj,1,1,1)
         zz=zm1(jj,1,1,1)

         if (zmax.le.zcut+tol) then ! miniumn z elts=2 < zcut ?
             el2d = mod1(eg,nelg2d) ! write(6,*) el2d, eg,e
             if ((abs(zz-zposition)).le.eps) then
                  e3dto2d (ie  )= el2d  ! assign 2d elt
                  e2dto3d (el2d)= ie    ! assign 3d elt
                  e3d2d   (eg  )= el2d  ! assign 2d elt
                  e2d3d   (el2d)= eg    ! assign 3d elt
                  k=k+1
                  cemslice(k)= jj      ! assign 3d index
                  ncemslice  = nx1*nz1*nel2d
                  !write(6,*) 'k',k,cemslice(k),ncemslice
                  if (k.gt.ncemslice) then
                      write(6,*) 'ERR:: k is gt ncemslice'
                      call exitt
                  endif
             endif 
         endif
         enddo
         enddo
      enddo

c ... global 
      call igop (e2d3d,iwrk2d,'+  ',nelg2d) 
      call igop (e3d2d,iwrk3d,'+  ',nelgt )

c ... local to global
      l = 0
      do ie= 1,nelt
         if (e3dto2d(ie).gt.0) then
             l=l+1
             el2d= e3dto2d(ie)
             lgel2d(l)= el2d  ! el2d=e3dto2d(eg)
         endif 
      enddo

      if (.not.(l.eq.nel2d)) then
           write(6,*) 'nel2d not matching max glle2d',l,nel2d
           call exitt
      endif


      if (nid.eq.0) write(6,*) '.'
      do i=1,nelg2d
         if (nid.eq.0) write(6,11) i, e2d3d(i)
      enddo
      if (nid.eq.0) write(6,*) '.'

      do i=1,nelgt 
         if (e3d2d(i).ge.1) then
             if (nid.eq.0) write(6,12) i, e3d2d(i)  
         else
            !write(6,13) FIXME                    
            !call exitt
         endif
      enddo
      if (nid.eq.0) write(6,*) '.'

   9  format(' nel2d =',i4)
  10  format(' zarry0/zposition::',i4,2f15.8)
  11  format(' 3D element  used::',2i4)
  12  format(' 2D element found::',2i4)
  13  format(' ERR:: No element having z=zposition')

      return
      end
c-----------------------------------------------------------------------
      subroutine set_charge_2d
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'POISSON2D'
      include 'NEKUSE'
      include 'EMWAVE'

      integer e,e2d,e3d,f3d

      do e=1,nel2d
         e2d = lgel2d (e  )
         e3d = e2dto3d(e2d)  ! get global element number
         f3d = 5           ! Assume, for now, only dealin w/ face 5
         call facexs(  x2d (1,1,e),   xm1(1,1,1,e3d),f3d,0)
         call facexs(  y2d (1,1,e),   ym1(1,1,1,e3d),f3d,0)
         call facexs(eps2d (1,1,e),epsilo(1,1,1,e3d),f3d,0)
      enddo

      n   = nx1*ny1*nel2d
      pi  = 4.0*atan(1.0)
      mode= param(35)     !FIXME temporary:: check out if any confilicts
      !what was this? x2dmax=glmax(x2d,nx1*ny1*nel2d)

      do i=1,n    
         xx = x2d(i,1,1)  
         yy = y2d(i,1,1)  
         rr = sqrt(xx**2+yy**2)
         theta= atan2(yy,xx)

         if     (mode.eq.0) then
            chg2d (i,1,1)= bnch_rdensity(rr)/permit           
            spot2d(i,1,1)= 0.0                         
            sex2d (i,1,1)= 0.0                
            sey2d (i,1,1)= 0.0                 
         elseif (mode.eq.1) then
            chg2d (i,1,1)= 4                             ! cylpec  
            spot2d(i,1,1)= (rr**2-1)                  
            ee1=-2*rr
            ee2= 0
            ee3= 0
            sex2d(i,1,1) = cos(theta)*ee1-sin(theta)*ee2
            sey2d(i,1,1) = sin(theta)*ee1+cos(theta)*ee2 
         elseif (mode.eq.2) then
            chg2d(i,1,1) =-2*pi*pi*cos(pi*xx)*cos(pi*yy) ! recper  
            spot2d(i,1,1)= cos(pi*xx)*cos(pi*yy)
            sex2d(i,1,1) = pi*sin(pi*xx)*cos(pi*yy)
            sey2d(i,1,1) = pi*cos(pi*xx)*sin(pi*yy)
         elseif (mode.eq.3) then
            chg2d(i,1,1) = -2*(2-xx**2-yy**2)             ! recpec  
            spot2d(i,1,1)= (1-xx**2)*(1-yy**2)
            sex2d(i,1,1) = 2*xx*(1-yy**2)
            sey2d(i,1,1) = 2*yy*(1-xx**2)
         elseif (mode.eq.4) then
            chg2d (i,1,1)= bnch_rdensity(rr) !exp(-0.5*(rr/0.1)**2)! no minus gaussian
            spot2d(i,1,1)= 0.0                         
            sex2d (i,1,1)= 0.0                
            sey2d (i,1,1)= 0.0                 
         elseif (mode.eq.5) then
            if (abs(rr).le.1e-14) then
            chg2d(i,1,1) = 1.0 ! bnch_rdensity(rr) !exp(-0.5*(rr/0.1)**2)! gaussian
            else
            chg2d(i,1,1) = 0.0 
            endif
            spot2d(i,1,1)= 0.0                         
            sex2d(i,1,1) = 0.0                
            sey2d(i,1,1) = 0.0                 
         elseif (mode.eq.6) then
            offaxis=0.3
            rr = sqrt((xx-offaxis)**2+yy**2)
            chg2d (i,1,1)= bnch_rdensity(rr)!exp(-0.5*(rr/0.1)**2) ! why minus ? FIXME               
            spot2d(i,1,1)= 0.0                         
            sex2d (i,1,1)= 0.0                
            sey2d (i,1,1)= 0.0                 
         else
            if (nid.eq.0) then
            write(6,*) 'set mode=0 (gaussian),1
     6      (cylpec),2(recper),3(recpec),6(off-axis)'
            write(6,*) 'otherwise define charge'
            call exitt
            endif
         endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_grs_2d
c
c    2D Poisson solver
c
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'POISSON2D'

      common /srcsf/ grs(3,lx1,ly1,lelt)
      integer        e,e2d,e3d,f3d
      common /jacinv/ jacmi(lx1,ly1,lz1,lelt)
      real    jacmi

      call invers2(jacmi,jacm1,nx1*ny1*nz1*nelt) 
      ! do inverse of jacm1 first
      ! otherwise it gets stuck if do "1/j3d"    

      do e=1,nel2d
         e2d = lgel2d (e  )
         e3d = e2dto3d(e2d)  ! get global element number
         f3d = 5           ! Assume, for now, only dealin w/ face 5

         call facexs(j3d (1,e),jacmi(1,1,1,e3d),f3d,0)  
         call facexs(rx2d(1,e),rxm1 (1,1,1,e3d),f3d,0)  ! rx*jac
         call facexs(ry2d(1,e),rym1 (1,1,1,e3d),f3d,0)
         call facexs(sx2d(1,e),sxm1 (1,1,1,e3d),f3d,0)
         call facexs(sy2d(1,e),sym1 (1,1,1,e3d),f3d,0)

         do i=1,nx1*ny1

          rx2d(i,e)=rx2d(i,e)*j3d(i,e)   ! This is necessary because nekton
          ry2d(i,e)=ry2d(i,e)*j3d(i,e)   ! stores Jac*rx in place of rx.
          sx2d(i,e)=sx2d(i,e)*j3d(i,e)
          sy2d(i,e)=sy2d(i,e)*j3d(i,e)
         
          grs(1,i,1,e)=
     $    (rx2d(i,e)*rx2d(i,e)+ry2d(i,e)*ry2d(i,e))*area(i,1,f3d,e3d)! G1M1
          grs(2,i,1,e)=
     $    (rx2d(i,e)*sx2d(i,e)+ry2d(i,e)*sy2d(i,e))*area(i,1,f3d,e3d)! G4M1
          grs(3,i,1,e)=
     $    (sx2d(i,e)*sx2d(i,e)+sy2d(i,e)*sy2d(i,e))*area(i,1,f3d,e3d)! G2M1

          bm2d(i,1,e)=area(i,1,f3d,e3d)             
         
         enddo
       
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine put_soln_2d          
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'POISSON'
      include 'POISSON2D'
      include 'BEAM'
      include 'NEKUSE'

      integer e,eg,el2d,elg2d

      nxyz  = nx1*ny1*nz1
      n     = nx1*ny1*nz1*nelt
      nxy2d = nx1*ny1*nel2d
      nxyg2d= nx1*ny1*nelg2d

      tol = 1.e-12
      rimpede  =sqrt(permea/permit)

      do el2d  = 1,nel2d
         elg2d = lgel2d(el2d)                                
         do i=1,nx1*ny1 
         gex2d (i,1,elg2d)= ex2d(i,1,el2d)    
         gey2d (i,1,elg2d)= ey2d(i,1,el2d)    
         ghx2d (i,1,elg2d)= hx2d(i,1,el2d)/rimpede    
         ghy2d (i,1,elg2d)= hy2d(i,1,el2d)/rimpede    
         enddo
      enddo

      call gop (gex2d,wrk2d,'+  ',nxyg2d)
      call gop (gey2d,wrk2d,'+  ',nxyg2d)
      call gop (ghx2d,wrk2d,'+  ',nxyg2d)
      call gop (ghy2d,wrk2d,'+  ',nxyg2d)

      do e= 1,nelt
         eg  = lglel(e,node)
         zmax= vlmax(zm1(1,1,1,e),nxyz)
        
         if (zmax.le.zcut+tol) then       !miniumn z elts=2 < zcut ?
             elg2d = mod1(eg,nelg2d)  

             do iz=1,nz1
             do i =1,nx1*ny1

                 xx = xm1(i,1,iz,e)
                 yy = ym1(i,1,iz,e)
                 zz = zm1(i,1,iz,e)
                 zzz= zm1(i,1,iz,e)
                 rr = xx**2+yy**2
                 if (rr.gt.0)  rr= sqrt(rr)            

                 if (param(30).eq.0) zz = bnch_z(rktime,zz)
                 if (param(30).eq.1) zz = bnch_z(  time,zz)

                 rdensity = bnch_rdensity(rr)                   
                 zdensity = bnch_zdensity(zz)              !*geps2d(i,1,elg2d)
                 
                 EX (i,1,iz,e) = zdensity* gEX2D(i,1,elg2d)!zdensity         
                 EY (i,1,iz,e) = zdensity* gEY2D(i,1,elg2d)!zdensity 
                 EZ (i,1,iz,e) = 0.0                                 
                 HX (i,1,iz,e) = zdensity* gHX2D(i,1,elg2d)!zdensity 
                 HY (i,1,iz,e) = zdensity* gHY2D(i,1,elg2d)!zdensity 
                 HZ (i,1,iz,e) = 0.0                               

             enddo
             enddo
        endif
      enddo 
      exmax=glmax(ex,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,u,n)   ! this is ax for 2D poisson

      include 'SIZE'
      include 'DXYZ'
      include 'POISSON2D'

       
      integer  n
      real     w(n),u(n)

      common /srcsf/ grs(3,lx1,ly1,lelt)
      common /srccg/ ur(lx1*ly1),us(lx1*ly1)
     $             , wr(lx1*ly1),ws(lx1*ly1)
      integer e

      nxy= nx1*ny1

      do e=1,nel2d

         call mxm(dxm1  ,nx1,u(1+nxy*(e-1)),nx1,ur,nx1)
         call mxm(u(1+nxy*(e-1)),nx1,dytm1 ,nx1,us,nx1)
         do i=1,nx1*ny1
            wr(i) = grs(1,i,1,e)*ur(i) + grs(2,i,1,e)*us(i)
            ws(i) = grs(2,i,1,e)*ur(i) + grs(3,i,1,e)*us(i)
         enddo

         call mxm(dxtm1,nx1,  wr,nx1,w(1+nxy*(e-1)),nx1)
         call mxm(ws   ,nx1,dym1,nx1,wr            ,nx1)

         do i=1,nx1*ny1
            w(i+nxy*(e-1)) = w(i+nxy*(e-1)) + wr(i)
         enddo

      enddo
      
      call fgs_gop(gsh_hndl_2d,w,'+\0') ! dssum
      call col2   (w,mask2d,n)          ! bdry

      return
      end
c-----------------------------------------------------------------------
c ggggg
      subroutine cg2d(x,f,mult,mask,r,w,p,z,n,niter)
      include 'SIZE'

c     Solve Ax=f where A is SPD and is invoked by the routine ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided routine ax(w,z,n) returns  w := Az,  
c
c     User-provided routine solveM(z,r,n) ) returns  z := M^-1 r,  
c
c     User-provided array wght() is used to scale inner-products
c     of the form (p,z) = p'*wght*z
c
      !common /srcsf/ grs(3,lx1,ly1,lelt)
      !common /srccg/ ur(lx1*ly1),us(lx1*ly1)
      !$             , wr(lx1*ly1),ws(lx1*ly1)

      integer n, niter
      real    eps, rhsnorm
      real    x(n),f(n),mult(n),mask(n),r(n),w(n),p(n),z(n)

      logical ifmcor


c     set machine tolerances
 
      one = 1.
      eps = 1.e-20
      rhsnorm = sqrt(glsc2(f,f,n))
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7
      eps= eps*rhsnorm

      !tube-9: just remain above: no additional setting for eps

      if (nid.eq.0) write(6,5) eps,n,niter
   5  format(' eps=',e15.7,', n=',i6, ', iter=',i6)

      pap = 0.0

      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)

      ifmcor= .false.
      skmin =  glmin(mask,n) 

      if (skmin.gt.0) ifmcor = .true.

      if (ifmcor) then
          smean = -1./glsum(mult,n)
          rmean = smean*glsc2(r,mult,n)
          call cadd(r,rmean,n)
      endif

      rnorm= sqrt(glsc3(r,mult,r,n)) ! inner product

      iter = 0
      if (nid.eq.0) write(6,*) 'iter;rnorm',iter,rnorm,ifmcor
      if (nid.eq.0) write(6,*) '.'                            

      miter = niter

      rtz1=1.0
      do 1000 iter=1,miter

         call solveM(z,r,n)         ! Invert preconditioner here

         if (ifmcor) then
             rmean = smean*glsc2(z,mult,n)
             call cadd(z,rmean,n)
         endif

         rtz2 = rtz1
         rtz1 = glsc3(r,mult,z,n)   ! parallel inner product
         beta = rtz1/rtz2

         if (iter.eq.1) beta=0.0
 
         call add2s1 (p,z,beta,n)    ! p=p*beta+z
         call ax     (w,p,n)

         pap = glsc3 (w,mult,p,n)

         alpha=  rtz1/pap
         alphm= -alpha

         call add2s2(x,p,alpha,n)
         call add2s2(r,w,alphm,n)

         rtr = glsc2(r,r,n)

         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)

         if (nid.eq.0)  write(6,6) iter,rnorm,alpha,beta,pap 
c        if (iter.le.10.or.mod(iter,10).eq.0) write(6,6) iter,rnorm
    6    format(' cg:',i8,1p4e12.4)

         if (rtr.le.rlim2) goto 1001

 1000 continue
 1001 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine solveM(z,r,n)
      real z(n),r(n)

      call copy(z,r,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine cg_2d
c
c     2D Poisson solver
c
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON2D'

      common /srcsf/ grs(3,lx1,ly1,lelv)

      common /srcuz/ r(lx1*ly1*lelv)  , w(lx1*ly1*lelv)
     $             , p(lx1*ly1*lelv)  , z(lx1*ly1*lelv)
     $             , u(lx1*ly1*lelv)  , f(lx1*ly1*lelv)
c    $             , wgt2d(lx1*ly1*lelv) 

      integer e,e3d,f3d

      real xx,yy,tmp

      n     = nx1*ny1*nel2d
      niter = 100000 ! temp
      pi    = 4.0*atan(1.0)

      call col3    (rhs2d,chg2d,bm2d,n)
      call chsign  (rhs2d,n)
      call fgs_gop (gsh_hndl_2d,rhs2d,'+\0') ! dssum
      call col2    (rhs2d,mask2d,n)

      call cg2d    (pot2d,rhs2d,mult2d,mask2d,r,w,p,z,n,niter)

      if (nid.eq.0) write(6,*) '.'             
      if (nid.eq.0) write(6,*) 'POISSON Solver Done!'

      call copy( phi2d, pot2d,n) ! copy solution   
      call copy(sphi2d,spot2d,n) ! copy solution   

      do i=1,n
         ephi2d(i)=abs(phi2d(i)-sphi2d(i))
      enddo

      phimax = glmax(ephi2d,n)
      if (nid.eq.0) write(6,*) 'poisson 2d in 3d: phimax=',phimax
      if (nid.eq.0) write(6,*) '.'             

      ! gradm1: caution the sign (const) !!! 

      if (param(35).eq.0) then
          const =  1.0   ! beam and wake problem
      else
          const = -1.0   ! always minus for examples: mode= 1,2,3,4
      endif
      call cem_grad2d (ex2d,ey2d,pot2d,const)

      if (nid.eq.0) write(6,*) 'Electrostatic Solution!'

      do i=1,n
         xx=x2d(i,1,1)
         yy=y2d(i,1,1)
         theta = atan2(yy,xx)
         e1=ex2d(i,1,1)
         e2=ey2d(i,1,1)
         e3=0.0
         call rectangular_to_cylindrical(ee1,ee2,ee3,e1,e2,e3,theta)
         hh1= 0.0
         hh2= ee1                                             
         hh3= 0.0
         call cylindrical_to_rectangular(h1,h2,h3,hh1,hh2,hh3,theta)
         hx2d(i,1,1)=h1
         hy2d(i,1,1)=h2
      enddo

      ! divergence

      call cem_div2d(divEX2d,divEY2d,ex2d,ey2d)  ! checked 
      call cem_div2d(divHX2d,divHY2d,hx2d,hy2d)

      do i=1,n
          divE2d(i,1,1)=divEX2d(i,1,1)+divEY2d(i,1,1)
          divH2d(i,1,1)=divHX2d(i,1,1)+divHY2d(i,1,1)
      enddo

      if (nid.eq.0) write(6,*) 'Divergence Check!'

      do i=1,n
         errEX2d  (i,1,1)= abs(sex2d (i,1,1)-ex2d (i,1,1))
         errEY2d  (i,1,1)= abs(sey2d (i,1,1)-ey2d (i,1,1))
         errDivE2d(i,1,1)= abs(divE2d(i,1,1)+chg2d(i,1,1)) ! chg2d=rho/epsilo
         errDivH2d(i,1,1)= abs(divH2d(i,1,1))
      enddo

      exmax   = glmax(errEX2d  ,n)
      exmax   = glmax(errEX2d  ,n)
      divEmax = glmax(errDivE2d,n)
      divHmax = glmax(errDivH2d,n)

      if (nid.eq.0) then
      write(6,*) '.'
      write(6,*) 'poisson 2d in 3d: exmax  =',exmax                    
      write(6,*) 'poisson 2d in 3d: eymax  =',eymax
      write(6,*) 'poisson 2d in 3d: divEmax=',divEmax
      write(6,*) 'poisson 2d in 3d: divHmax=',divHmax
      write(6,*) '.'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine cem_grad2d(w1,w2,u1,const)
      include 'SIZE'
      include 'DXYZ'
      include 'POISSON2D'
      real     w1(1),w2(1),u1(1)
      real     const

      integer          e                  

      nn = nx1-1

      do e=1,nel2d
         j = nx1*ny1*(e-1)
         call local_grad2(ur2d,us2d,u1,nn,e,dxm1,dxtm1)
         do i=1,nx1*ny1
            k=i+j
            w1(k)=const*(ur2d(i)*rx2d(k,1)+us2d(i)*sx2d(k,1))
            w2(k)=const*(ur2d(i)*ry2d(k,1)+us2d(i)*sy2d(k,1))
            !write(6,*) 'ex',e,i,rx2d(k,1),ry2d(k,1)
            !write(6,*) 'ex',i,ur2d(i),us2d(i)
            !write(6,*) 'ex',e,i,w1(k),w2(k)                 
            !write(6,*) 'ex',e,i,x2d(i,1,e),y2d(i,1,e)            
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_div2d(w1,w2,u1,u2)
      include 'SIZE'
      include 'DXYZ'
      include 'POISSON2D'
      real     w1(1),w2(1),u1(1),u2(1)
      common  /div2d/ ur2d1(lx1*ly1*lelt),us2d1(lx1*ly1*lelt)
     $              , ur2d2(lx1*ly1*lelt),us2d2(lx1*ly1*lelt)
      integer  e                  

      nn = nx1-1

      do e=1,nel2d
         j = nx1*ny1*(e-1)
         call local_grad2(ur2d1,us2d1,u1,nn,e,dxm1,dxtm1)
         call local_grad2(ur2d2,us2d2,u2,nn,e,dxm1,dxtm1)
         do i=1,nx1*ny1
            k=i+j
            w1(k)=(ur2d1(i)*rx2d(k,1)+us2d1(i)*sx2d(k,1))
            w2(k)=(ur2d2(i)*ry2d(k,1)+us2d2(i)*sy2d(k,1))
           ! write(6,*) 'ex',e,i,u1(k)                  
           !write(6,*) 'ex',e,i,rx2d(k,1),ry2d(k,1)
           !write(6,*) 'ex',i,ur2d(i),us2d(i)
           ! write(6,*) 'ex',e,i,w1(k),w2(k)                 
           !write(6,*) 'ex',e,i,x2d(i,1,e),y2d(i,1,e)            
         enddo
      enddo

      return
      end
