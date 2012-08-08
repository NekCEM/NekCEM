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

      if   (IFSCHROD ) then

          if (ifrk ) call cem_schrod_op_rk
          if (ifeig) call cem_schrod_op_eig

      else

          if (ifrk ) call cem_maxwell_op_rk
          if (ifeig) call cem_maxwell_op_eig

      endif

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields6(dumpno)
C----------------------------------------------------------------------
c...  this works OK with param(81)=6,-6,5,8:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-07-08 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'

      integer n,dumpno
      integer maxnelt,numgroups,numfields,maxnumnodes
      save    maxnelt,numgroups,numfields,maxnumnodes

      integer icalld        
      save    icalld          
      data    icalld /0/               

      common  ncell,ncellsize,ncelltype     
      integer ncell,ncellsize,ncelltype

      common /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real            vtkxyz
      common /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer         fpcnt
      common /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      if (icalld.eq.0) then
          maxnelt = iglmax(nelt,1)
          call bcast(maxnelt,isize)
          maxnumnodes = nx1*ny1*nz1*maxnelt ! max # nodes per processor

          if (nid.eq.0) then
              numfields   = param(80)           ! this number has to be exact.                            
              numgroups   = param(82)           ! param(82) = # groups 
          endif
          call bcast(numfields,isize)
          call bcast(numgroups,isize)

          call vtk_dump_header_nonswap()!save nodes and cells 
          icalld=1
      endif

      npts = nx1*ny1*nz1*nelt

c ... VTK FORMAT BINARY

      io_option= 6
      call print_options(io_option)

      call initrbio    (numgroups,numfields,maxnumnodes) !rbio_nekcem.c
      call openfile6   (dumpno,nid)                     
      call writeheader6()
      call writenodes6 (vtkxyz,npts)

      if  (if3d) then
          call write3dcells6(fpcnt,nelt,ncell,npts)
      else
          call write2dcells6(fpcnt,nelt,ncell,npts)
      endif

c...  number of fields are four: 
      if (ifschrod) then
      call writefield6 (11,vtkout1,npts) 
      call writefield6 (12,vtkout2,npts) 
      call writefield6 (13,vtkout3,npts)
      call writefield6 (14,vtkout4,npts)
      else
      call writefield6 (1,vtkout1,npts) 
      call writefield6 (2,vtkout2,npts) 
      call writefield6 (3,vtkout3,npts)
      call writefield6 (4,vtkout4,npts)
      endif

      inumfields=4

      if (inumfields.ne.numfields) then
          if (nid.eq.0) write(6,*) 'ERR: param(80) should be set as 4'
          call exitt
      endif 

      call closefile6  ()

c     call pvtk_xml(dumpno)
c     if(param(81).eq.8) then
c       if(nid.eq.0) call pvtk_nmm(dumpno)
c     elseif(param(81).eq.18) then
c       if(nid.eq.0) call pvtk_nmm(dumpno)
c     elseif(param(81).eq.5) then
c       if(nid.eq.0) call pvtk_nm(dumpno)
c     endif

      return
      end

C----------------------------------------------------------------------
      subroutine cem_restart_out(dumpno)    
C----------------------------------------------------------------------
c...  this works OK with param(81)=4:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-06-28 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      include 'SCHROD'

      integer icalld
      save    icalld            
      data    icalld /0/
      integer n,dumpno,idump

      common  ncell,ncellsize,ncelltype     
      integer ncell,ncellsize,ncelltype
      common  /s_rvtk/ vtkxyz_s(3*lx1*ly1*lz1*lelt)
      real             vtkxyz_s
      common  /s_ivtk/ fpcnt_s (8*lx1*ly1*lz1*lelt)
      integer          fpcnt_s
      common  /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               , vtkout2(3*lx1*ly1*lz1*lelt)
     $               , vtkout3(3*lx1*ly1*lz1*lelt)
     $               , vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      npts = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
         call vtk_dump_header_swap()     !save nodes and cells 
         icalld=1
      endif

      call openfile_restart (dumpno,nid) !vtkbin.c

c...  calling c-routines
      call writeheader4(istep,dumpno,nx1,np,time,dt)
      call writenodes4 (vtkxyz_s,npts)
      if (if3d) then
          call write3dcells4_swap(fpcnt_s,nelt,ncell,npts)
      else
          call write2dcells4_swap(fpcnt_s,nelt,ncell,npts)  
      endif

      if (ifschrod) then

        if (ifdouble) then
          call writefield4_double (11,vtkout1,npts)           
          call writefield4_double (12,vtkout2,npts) 
          call writefield4_double (13,vtkout3,npts)           
          call writefield4_double (14,vtkout4,npts) 
        else
          call writefield4 (11,vtkout1,npts)           
          call writefield4 (12,vtkout2,npts) 
          call writefield4 (13,vtkout3,npts)           
          call writefield4 (14,vtkout4,npts) 
        endif

      else

        if (ifdouble) then
          call writefield4_double (1,vtkout1,npts)           
          call writefield4_double (2,vtkout2,npts) 
          call writefield4_double (3,vtkout3,npts)           
          call writefield4_double (4,vtkout4,npts) 
        else
          call writefield4 (1,vtkout1,npts)           
          call writefield4 (2,vtkout2,npts) 
          call writefield4 (3,vtkout3,npts)           
          call writefield4 (4,vtkout4,npts) 
        endif

      endif

      call closefile_restart () 

c     if ((nid.eq.0).and.(mod(istep,iostep).eq.0)) then
c       write(6,*) ' restart-file-write completed  ::::'
c     endif

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields4(dumpno)
C----------------------------------------------------------------------
c...  this works OK with param(81)=4:  
c...  this produces a single file using MPI collective write function
c...  Jing Fu 2010-06-28 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'SCHROD'

      integer icalld
      save    icalld         
      data    icalld /0/                 
      integer n,dumpno,idump

      common  ncell,ncellsize,ncelltype     
      integer ncell,ncellsize,ncelltype

      common /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real            vtkxyz
      common /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer         fpcnt
      common /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      io_option= 4
      call print_options(io_option)

      npts = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call vtk_dump_header_nonswap() !save nodes and cells 
          icalld=1
      endif

c...  calling c-routines
      call openfile4   (dumpno,nid)   ! vtkbin.c
      call writeheader4(istep,dumpno,nx1,np,time,dt)
      call writenodes4 (vtkxyz,npts)
      if (if3d) then
          call write3dcells4(fpcnt,nelt,ncell,npts)
      else
          call write2dcells4(fpcnt,nelt,ncell,npts)  
      endif

      if (ifschrod) then

      call writefield4 (11,vtkout1,npts) 
      call writefield4 (12,vtkout2,npts) 
      call writefield4 (13,vtkout3,npts)
      call writefield4 (14,vtkout4,npts)

      else
      call writefield4 (1,vtkout1,npts) 
      call writefield4 (2,vtkout2,npts) 
      call writefield4 (3,vtkout3,npts)
      call writefield4 (4,vtkout4,npts)
      endif

      call closefile4  () 

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields3(dumpno)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'SCHROD'

      integer icalld
      save    icalld         
      data    icalld /0/
      integer n,dumpno,idump

      common  ncell,ncellsize,ncelltype
      integer ncell,ncellsize,ncelltype

      common /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real            vtkxyz
      common /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer         fpcnt
      common /outvtk/ vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)

      io_option= 3
      call print_options(io_option)

      npts = nx1*ny1*nz1*nelt

      if (icalld.eq.0) then
          call vtk_dump_header_nonswap() !save nodes and cells
          icalld=1
      endif

c...  calling c-routines
      call openfile   (dumpno,nid)     ! vtkbin.c
      call writeheader()
      call writenodes (vtkxyz,npts)
      if (if3d) then
          call write3dcells(fpcnt,nelt,ncell,npts)
      else
          call write2dcells(fpcnt,nelt,ncell,npts)
      endif

      if (ifschrod) then
      call writefield (11,vtkout1,npts) 
      call writefield (12,vtkout2,npts) 
      call writefield (13,vtkout3,npts)
      call writefield (14,vtkout4,npts)
      else
      call writefield (1,vtkout1,npts) 
      call writefield (2,vtkout2,npts) 
      call writefield (3,vtkout3,npts)
      call writefield (4,vtkout4,npts)
      endif

      call closefile()  ! vtkbin.c

      if(nid.eq.0) call pvtk_xml3(dumpno)

      return 
      end
C----------------------------------------------------------------------
      subroutine cem_out_fields2(dumpno)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'SCHROD'
      
      integer    fnum, dumpno
      parameter (fnum=111)
      logical    outex, outey, outez, outhx, outhy, outhz
      integer    n

      character*33 filename
      data         filename /'./vtk/ascii-NN-p000000-t00000.vtk'/
      CHARACTER*1  fnarray(33)
      equivalence  (filename,fnarray)

      integer   icalld
      save      icalld          
      data      icalld /0/             

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

      io_option= 2
      call print_options(io_option)

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

      return 
      end

C----------------------------------------------------------------------
      subroutine cem_out_fields(dumpno)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'EXPONE'
      include 'NEKUSE'
      include 'PML'
      include 'SCHROD'
      
      logical outex, outey, outez, outhx, outhy, outhz
      integer n,dumpno
      integer icalld
      save    icalld
      data    icalld /0/


      ifxyo=.true.
      if (istep.gt.iostep) ifxyo = .false.  ! dump geom for 1st file only
      ifvo  = .true.
      ifpo  = .false.
      ifto  = .false.

      io_option= 0
      call print_options(io_option)

      n = nx1*ny1*nz1*nelt

        if (NSTEPS.eq.0) then
C           call outpost(hx,hy,hz,pr,hx,'hxy')
C           call outpost(ex,ey,ez,pr,ey,'exy')
            return
        endif

        if (ifsol) then

c         call outpost(hx,hy,hz,pr,hx,'hxy')
c         call outpost(ex,ey,ez,pr,ey,'exy')

        else 

          if (ifscat) then

c             call outpost(tothx,tothy,tothz,pr,hx,'h-t')
c             call outpost(totex,totey,totez,pr,ey,'e-t')
c             call outpost(hx,hy,hz,pr,hx,'hxy')
c             call outpost(ex,ey,ez,pr,ey,'exy')
c             call outpost(engyh,engyh,engyh,pr,ey,'enh')
c             call outpost(engye,engye,engye,pr,ex,'ene')

          else

c            call outpost(hx,hy,hz,pr,hx,'hxy')
c            call outpost(ex,ey,ez,pr,ey,'exy')

          endif ! ifscat

        endif ! ifsol

      return 
      end

C----------------------------------------------------------------------
      subroutine print_options(io_option)
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      
      integer  io_option

      if ((mod(istep,iostep).eq.0).or.(istep.eq.nsteps)) then
      if (nid.eq.0) then
        write(6,*)
        write(6,34) ifsrc,iftm,ifte,ifcentral,ifupwind,ifpml,ifpec,
     $  ifeig,ifexp,ifrk45,ifrk44,ifrk33,ifrk22,ifsim5,ifsim4,ifsim3,
     $  ifsim2,ifsim1                     

        write(6,35) ifdealias,ifscat,ifsftf,ifnm,ifdielec,iffilter,
     $  ifsol,ifpoisson,ifdrude,ifloren,ifwake,ifnapy,ifbeam,lxd,
     $  marnol

c...   not printint due to check if (io_option.eq.0) for now, it ok
        if (io_option.lt.0) return
        if (io_option.eq.0) write(6,3) istep,io_option
        if (io_option.eq.0) write(6,4) istep,io_option
        if (io_option.eq.0) write(6,5) istep,io_option
        if (io_option.eq.0) write(6,6) istep,io_option
        if (io_option.eq.0) write(6,7) istep,io_option
      endif

      endif

   33 format(a12,i10)
   34 format(' SRC=',L1,', TM=',L1,', TE=',L1,', CENTRAL=',L1,
     $ ', UPWIND=',L1,', PML=',L1,', PEC=',L1,' EIG=',L1,', EXP=',L1,
     $ ', RK45=',L1,', RK4=',L1,', RK3=',L1,', RK2=',L1,'SIM5=',L1,
     $ ', SIM4=',L1,', SIM3=',L1,', SIM2=',L1,', SIM1=',L1)
   35 format(' DEALIAS=',L1,', SCAT=',L1,', SFTF=',L1,', NM=',L1,
     $ ', DIELEC=',L1,', FILTER=',L1,', SOL=',L1,', POISSON=',L1,
     $ ', DRUDE=',L1,', LORENTZ=',L1,', WAKE=',L1,', NAPOLY=',L1,
     $ ', BEAM=',L1,', Dealias_M=',i3,', Krylov_Dim=',i3,'\n')     

   3  format(' ---> dumping postx format, istep=',i8,', iooption=',i3)
   4  format(' ---> dumping ascii  vtk, istep=',i8,', iooption=',i3)
   5  format(' ---> dumping binary vtk, istep=',i8,', iooption=',i3)
   6  format(' ---> dumping binary vtk, istep=',i8,', iooption=',i3)
   7  format(' ---> dumping binary vtk, istep=',i8,', iooption=',i3)

      return
      end
C----------------------------------------------------------------------
      subroutine cem_out
C----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD'
      include 'PML'    

      common  /outvtk/vtkout1(3*lx1*ly1*lz1*lelt)
     $               ,vtkout2(3*lx1*ly1*lz1*lelt)
     $               ,vtkout3(3*lx1*ly1*lz1*lelt)
     $               ,vtkout4(3*lx1*ly1*lz1*lelt)
      real    vtkout1,vtkout2,vtkout3,vtkout4

      common /ccpu  / cpu_t,cpu_dtime,cpu_chk
      common /commt / comm_t
      common /iotime/ cpu_t_io, cpu_io_step, swaptime
      real            cpu_t,cpu_dtime,cpu_chk,comm_t 
      real            cpu_t_io, cpu_io_step, swaptime

      common /iocnt/  dumpno
      integer         dumpno
      integer iooption,numfiles,irestart,idouble
      save    iooption,numfiles,irestart,idouble
      integer icomptrace,iotrace
      save    icomptrace,iotrace        
      integer iparam1,iparam2,iparam3,iparam4,iparam5,iparam6

      integer icalld
      save    icalld
      data    icalld/0/
      real*8  dclock

c...  dumpno for restarting case comes from previous step
      if (icalld.eq.0) then
         cpu_t_io= 0.0

         if (nid.eq.0) then
             iparam1=param(81)
             iparam2=param(82)
             iparam3=param(83)      !io-frequency for restart files     
             iparam4=0              !restart dumpno             
             iparam5=param(85)      !write compute trace                  
             iparam6=param(86)      !write io trace                       
             iparam7=param(87)      !write io trace                       
         else
             iparam1=0             
             iparam2=0            
             iparam3=0              !io-frequency for restart files     
             iparam4=0                 
             iparam5=0              !write compute trace                  
             iparam6=0              !write io trace                       
             iparam7=0              !write io trace                       
         endif
         if (ifrestart) then
            if (nid.eq.0) then
             iparam4= param(84)+1  !restart dumpno update
            else
             iparam4= 0                     
            endif
         endif
c...     iglsum is not doing smart bcast'ing, revert back to bcast -- jingfu
         call bcast(iparam1,isize) 
         call bcast(iparam2,isize) 
         call bcast(iparam3,isize) 
         call bcast(iparam4,isize) 
         call bcast(iparam5,isize) 
         call bcast(iparam6,isize) 
         call bcast(iparam7,isize) 

         iooption    =iparam1      
         numfiles    =iparam2        
         irestart    =iparam3                                       
         dumpno      =iparam4
         icomptrace  =iparam5
         iotrace     =iparam6
         idouble     =iparam7

         IFDOUBLE=.true.  
         if (idouble.ne.0) IFDOUBLE=.false. 
         if (nid.eq.0) write(6,*) 'IFDOUBLE=',IFDOUBLE
         if (nid.eq.0) write(6,*) 'param(81)=',iparam1
         if (nid.eq.0) write(6,*) 'param(82)=',iparam2
         if (nid.eq.0) write(6,*) 'param(83)=',irestart            
         if (nid.eq.0) write(6,*) 'dumpno   =',dumpno

         icalld  = 1
         call pass_io_params(iooption,numfiles)
      endif

c...  param(85) decide if compute trace should be generated
      if (icomptrace.gt.0) call writecomputetrace(istep,cpu_dtime,cpu_t)

c...  no output/restart files will be generated
      if (iooption.eq.0) return 

c...  this produces output for (io-option=iparam) with timing info
      if (iostep.gt.0) then

      if ((mod(istep,iostep).eq.0).or.(istep.eq.nsteps)) then
        
         if (ifschrod) then

         call vtk_nonswap_field( Ur, Ui, Ui,vtkout1)
         call vtk_nonswap_field( Qr(1,1),Qr(1,2),Qr(1,3),vtkout2)
         call vtk_nonswap_field( Qi(1,1),Qi(1,2),Qi(1,3),vtkout3)
         call vtk_nonswap_field( sUr,sQr(1,1),sQi(1,1),vtkout4)

         else

         call vtk_nonswap_field( HN(1,1), HN(1,2), HN(1,3),vtkout1)
         call vtk_nonswap_field( EN(1,1), EN(1,2), EN(1,3),vtkout2)
         call vtk_nonswap_field(sHN(1,1),sHN(1,2),sHN(1,3),vtkout3)
c        call vtk_nonswap_field(sEN(1,1),sEN(1,2),sEN(1,3),vtkout4)

         call vtk_nonswap_field( pmlsigma(1,1), pmlsigma(1,2), 
     $                           pmlsigma(1,3), vtkout4)

         endif

         pretime = dclock()
         call starttiming
         call checkpoint_write(iooption,dumpno)
c...     maybe need a thread_join function to ensure the last io step finished
         call endtiming
         endtime = dclock()-pretime
         cpu_t_io= cpu_t_io+endtime

         call printio(iooption,istep)

c...  produce restart files    
         if (irestart.gt.0) then !io-frequency for restart files >= nonzero

          if (mod(istep,irestart*iostep).eq.0) then

#ifdef NO_MPI
           write(6,*) 'non-mpi run: no restart file produced'
#else
           swaptime=dclock() 

           if (ifschrod) then
               call vtk_swap_field( Ur, Ui, Ui,vtkout1)
               call vtk_swap_field( Qr(1,1),Qr(1,2),Qr(1,3),vtkout2)
               call vtk_swap_field( Qi(1,1),Qi(1,2),Qi(1,3),vtkout3)
               call vtk_swap_field( sUr,sQr(1,1),sQi(1,1),vtkout4)
           else
               call vtk_swap_field(HN(1,1),HN(1,2),HN(1,3),vtkout1)
               call vtk_swap_field(EN(1,1),EN(1,2),EN(1,3),vtkout2)
               call vtk_swap_field(sHN(1,1),sHN(1,2),sHN(1,3),vtkout3)
               call vtk_swap_field(sEN(1,1),sEN(1,2),sEN(1,3),vtkout4)
           endif

           swaptime=dclock()-swaptime 
           swaptime=swaptime/4.0
           i99=99

          !For threading, restart frequency should be zero on BG/P: fixme -- misun/jing 3/23/2012 
          !To be tested the same case on XK6: not likely to cause the same situation.                      
           if (nid.eq.0) write(6,3) dumpno,irestart*iostep,swaptime
           call checkpoint_write(i99,dumpno)
#endif

          endif     
         endif     
    3    format(' ---> write restart files at dumpno=',i10,
     $          ', irestart*iostep=',i10,', swaptime=',e11.4,'\n')
         dumpno = dumpno + 1  

      endif
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine cem_drude 
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      n = nx1*ny1*nz1*nelt    

c...  simpler with explicit writing than calling mat1 routines many times
      do i=1,n
        !cdrude1 = DRUDEA(i)!*bm1(i,1,1,1)
        !cdrude2 = DRUDEB(i)!*bm1(i,1,1,1)
 
        !resJN(i,1)= cdrude1*JN(i,1)+cdrude2*EN(i,1)
        !resJN(i,2)= cdrude1*JN(i,2)+cdrude2*EN(i,2)
        !resJN(i,3)= cdrude1*JN(i,3)+cdrude2*EN(i,3)
 
        !resJX(i,1,1,1,1)= resJx(i,1,1,1,1)+incEx(i,1,1,1)
        !resJY(i,1,1,1,1)= resJy(i,1,1,1,1)+incEy(i,1,1,1)
        !resJZ(i,1,1,1,1)= resJz(i,1,1,1,1)+incEz(i,1,1,1)

        !resEN(i,1) = resEN(i,1)-JN(i,1)*bmn(i)
        !resEN(i,2) = resEN(i,2)-JN(i,2)*bmn(i)
        !resEN(i,3) = resEN(i,3)-JN(i,3)*bmn(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_lorentz
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      n = nx1*ny1*nz1*nelt    

      do i=1,n

c        clorenA1 = LORENA(i,1)!*bm1(i,1,1,1)
c        clorenB1 = LORENB(i,1)!*bm1(i,1,1,1)
c        clorenC1 = LORENC(i,1)!*bm1(i,1,1,1)
c
c        clorenA1_Jx= clorenA1*JN(i,1)                
c        clorenA1_Jy= clorenA1*JN(i,2)                
c        clorenA1_Jz= clorenA1*JN(i,3)                
c
c        clorenB1_Px= clorenB1*PN(i,1)                
c        clorenB1_Py= clorenB1*PN(i,2)                
c        clorenB1_Pz= clorenB1*PN(i,3)                
c
c        clorenC1_Ex= clorenC1*EN(i)                
c        clorenC1_Ey= clorenC1*EN(i)                
c        clorenC1_Ez= clorenC1*EN(i)                
c
c        resJN(i,1)= clorenA1_Jx+clorenB1_Px+clorenC1_Ex 
c        resJN(i,2)= c1orenA1_Jy+clorenB1_Py+clorenC1_Ey
c        resJN(i,3)= clorenA1_Jz+clorenB1_Pz+clorenC1_Ez
c
c        clorenA2 = LORENA(i,2)!*bm1(i,1,1,1)
c        clorenB2 = LORENB(i,2)!*bm1(i,1,1,1)
c        clorenC2 = LORENC(i,2)!*bm1(i,1,1,1)
c
c        clorenA2_Jx= clorenA2*JN(i,2)                
c        clorenA2_Jy= clorenA2*JN(i,2)                
c        clorenA2_Jz= clorenA2*JN(i,2)                
c
c        clorenB2_Px= clorenB2*PN(i,2)                
c        clorenB2_Py= clorenB2*PN(i,2)                
c        clorenB2_Pz= clorenB2*PN(i,2)                
c
c        clorenC2_Ex= clorenC2*EN(i)                
c        clorenC2_Ey= clorenC2*EN(i)                
c        clorenC2_Ez= clorenC2*EN(i)                
c
c        resJx(i,2)= clorenA2_Jx+clorenB2_Px+clorenC2_Ex 
c        resJy(i,2)= c1orenA2_Jy+clorenB2_Py+clorenC2_Ey
c        resJz(i,2)= clorenA2_Jz+clorenB2_Pz+clorenC2_Ez
c
c        clorenA3 = LORENA(i,3)!*bm1(i,1,1,1)
c        clorenB3 = LORENB(i,3)!*bm1(i,1,1,1)
c        clorenC3 = LORENC(i,3)!*bm1(i,1,1,1)

c        clorenA3_Jx= clorenA3*JN(i,3)                
c        clorenA3_Jy= clorenA3*JN(i,3)                
c        clorenA3_Jz= clorenA3*JN(i,3)                
c
c        clorenB3_Px= clorenB3*PN(i,3)                
c        clorenB3_Py= clorenB3*PN(i,3)                
c        clorenB3_Pz= clorenB3*PN(i,3)                
c
c        clorenC3_Ex= clorenC3*EN(i)                
c        clorenC3_Ey= clorenC3*EN(i)                
c        clorenC3_Ez= clorenC3*EN(i)                
c
c        resJN(i,3)= clorenA3_Jx+clorenB3_Px+clorenC3_Ex 
c        resJN(i,3)= c1orenA3_Jy+clorenB3_Py+clorenC3_Ey
c        resJN(i,3)= clorenA3_Jz+clorenB3_Pz+clorenC3_Ez
c
c        cJx=(JN(i,1)+JN(i,2)+JN(i,3))*bmn(i)
c        cJy=(JN(i,1)+JN(i,2)+JN(i,3))*bmn(i)
c        cJz=(JN(i,1)+JN(i,2)+JN(i,3))*bmn(i)
c
c        resEN(i)= resEN(i)-cJx                 
c        resEN(i)= resEN(i)-cJy                 
c        resEN(i)= resEN(i)-cJz            
c
c        resPN(i,1)= JN(i,1)                     
c        resPN(i,1)= JN(i,1)                     
c        resPN(i,1)= JN(i,1)                
c        resPN(i,2)= JN(i,2)                     
c        resPN(i,2)= JN(i,2)                     
c        resPN(i,2)= JN(i,2)                
c        resPN(i,3)= JN(i,3)                     
c        resPN(i,3)= JN(i,3)                     
c        resPN(i,3)= JN(i,3)                

      enddo

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


      return
      end
c-----------------------------------------------------------------------
      subroutine rk4_upd(h,kh,resh,cb,ca,dt,n)
c-----------------------------------------------------------------------
      real h(1),kh(1),resh(1)

c     do i=1,n          
c        kh(i  ) = ca*kh(i  ) + dt*resh(i  )
c        h (i  ) =     h(i  ) + cb*kh  (i  )
c     enddo
c     return


c     below routine calculates 4 consecutive ones at one time
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

c     rk4a(1) =   0.0 
c     rk4a(2) =  -11/15 !756391.0/934407.0
c     rk4a(3) =  -5/3   !36441873.0/15625003.0
c     rk4a(4) =  -1     !1953125.0/1085297.0
c     !rk4a(5) =   0.0
c
c     rk4b(1) =   1/3  !8.0/141.0 
c     rk4b(2) =   5/6  !6627.0/2000.0
c     rk4b(3) =   3/5  !609375.0/1085297.0
c     rk4b(4) =   1/4  !198961.0/526383.0
c     !rk4b(5) =   0.0
c
c     rk4c(1) =   0.0
c     rk4c(2) =   1.0/3.0
c     rk4c(3) =   5.0/9.0
c     rk4c(4) =   8.0/9.0
c     !rk4c(5) =   1.0
c
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine restart_swap
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'SCHROD' 

      common  /RSTTMP/  xmr1(lx1,ly1,lz1,lelt)
     $                , ymr1(lx1,ly1,lz1,lelt)
     $                , zmr1(lx1,ly1,lz1,lelt)
     $                , exr (lx1,ly1,lz1,lelt)
     $                , eyr (lx1,ly1,lz1,lelt)
     $                , ezr (lx1,ly1,lz1,lelt)
     $                , hxr (lx1,ly1,lz1,lelt)
     $                , hyr (lx1,ly1,lz1,lelt)
     $                , hzr (lx1,ly1,lz1,lelt)
      common  /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real             vtkxyz
      common  /c_ivtk/ fpcnt (8*lx1*lx1*lx1*lelt)
      integer          fpcnt

      !write in float --> read in double: doesn't work
      !write in float --> read in float : recover up to 7-8 digits
      !write in double--> read in double: recover up to 16 digits

      integer n,dumpno,idump,idummy
      integer iooption,numfiles     
      real*8  dclock,t0

      integer e,eg

      integer icalld
      save    icalld
      data    icalld/0/

      npts = nx1*ny1*nz1*nelt !                            
      nxyz = nx1*ny1*nz1      !                               
      nxyz3= nx1*ny1*nz1*3    !

      t0= dclock()

      if (nid.eq.0) then
          iooption    = param(81)
          numfiles    = param(82)
          irstdumpno  = param(84) 
          idouble     = param(87) 
      else
          iooption    = 0           
          numfiles    = 0        
          irstdumpno  = 0         
          idouble     = 0         
      endif
 
      call bcast(iooption  ,isize)
      call bcast(numfiles  ,isize)
      call bcast(irstdumpno,isize)
      call bcast(idouble   ,isize)

      IFDOUBLE=.true.                    ! zero: double
      if (idouble.ne.0) IFDOUBLE=.false. ! nonzero: float
      if (nid.eq.0) write(6,*) 'IFDOUBLE=',IFDOUBLE
      if (nid.eq.0) write(6,*) 'restart format param(87)=',idouble  
      if (nid.eq.0) write(6,*) 'restart param(84)=',irstdumpno

      t0= dclock()
      call pass_io_params  (iooption,numfiles)
      call openfile_restart(irstdumpno, nid)
      call readheader4     (irststep,idump,nn1,nn2,restart_time,dtvtk)

      if (nid.eq.0) then
          time = time + restart_time
      else
          time = 0.0
      endif
      time = glsum(time,1)
      if (nid.eq.0) write(6,*) 'restart at time=time+rst_time',time

      !call printchecks  (idump)
      call readnodes4   (vtkxyz,npts) ! dummy reading but necessary     
      call swap_real_backward(vtkxyz,nxyz3) 
      call save2vectors (xmr1,ymr1,zmr1,vtkxyz,npts)
      
      if (if3d) then
           call read3dcells4(idummy,idummy,idummy,npts)  !no return
      else
           call read2dcells4(idummy,idummy,idummy,npts)  !no return
      endif

      if (ifschrod) then

          if (ifdouble) then
          call readfield4_double  (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (ur,ui,ui,vtkxyz,npts)
          call readfield4_double  (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qr(1,1),qr(1,2),qr(1,3),vtkxyz,npts)
          call readfield4_double  (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qi(1,1),qi(1,2),qi(1,3),vtkxyz,npts)
          else
          call readfield4         (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (ur,ui,ui,vtkxyz,npts)
          call readfield4         (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qr(1,1),qr(1,2),qr(1,3),vtkxyz,npts)
          call readfield4         (idummy,vtkxyz,npts)    
          call swap_real_backward (vtkxyz,nxyz3)
          call save2vectors       (qi(1,1),qi(1,2),qi(1,3),vtkxyz,npts)
          endif

      else

          if (ifdouble) then
          call readfield4_double   (idummy,vtkxyz,npts)   
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (hn(1,1),hn(1,2),hn(1,3),vtkxyz,npts)
          call readfield4_double   (idummy,vtkxyz,npts)   
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (en(1,1),en(1,2),en(1,3),vtkxyz,npts)
          else
          call readfield4          (idummy,vtkxyz,npts)   
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (hn(1,1),hn(1,2),hn(1,3),vtkxyz,npts)
          call readfield4          (idummy,vtkxyz,npts)    
          call swap_real_backward  (vtkxyz,nxyz3)
          call save2vectors        (en(1,1),en(1,2),en(1,3),vtkxyz,npts)
          endif

      endif
      call closefile_restart

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
      subroutine  printchecks(idump)   
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      
      integer     idump

      irstdumpno =idump       
      istep      =irststep

      if (nid.eq.0) write(6,*) 'restart timestep=',irststep,
     $ ', restart dumpno=',irstdumpno,', restart time=', time,'istep=',
     $    istep, 'nid',nid

      return
      end
c-----------------------------------------------------------------------
      subroutine cem_inc                            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      !call userinc                                           

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

      if (.not.ifsrc) return

      call usersrc (57,reshn(1,1),reshn(1,2),reshn(1,3),
     $                 resen(1,1),resen(1,2),resen(1,3))

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
c
      integer i,e,f,ef,js1,jf1,jskip1,js2,jf2,jskip2,j1,j2

      nxyz  = nx1*ny1*nz1
      nxz   = nx1*nz1
      nface = 2*ndim
      nxzf  = nx1*nz1*nface ! red'd mod to area, unx, etc.

      if (nid.eq.0) then
         write(6,*) 'start: cem_set_fc_ptr, volume index j=cemface(i)'
      endif

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

      if (nid.eq.0) then
         write(6,*) 'done: cem_set_fc_ptr, ncemface= ',ncemface
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine measure_comm(t0)
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'EMWAVE'
      include 'PARALLEL'
      common  /commt/ comm_t
      real            comm_t
      integer icalld
      save    icalld
      data    icalld/0/
      real    t0
      real*8  dclock

      if (icalld.eq.0) comm_t=0.0

      dtime = dclock()-t0
      comm_t= comm_t+dtime
      icalld= icalld+1

      ifreq=5             
      if (ifexp) ifreq=larnol

      if ((mod(istep,iocomm).eq.0).and.(mod(icalld,ifreq).eq.0)) then
          if (nid.eq.0) write(6,3) dtime,comm_t  
      endif
    3 format(' gs_op_fields::',1pe11.4,' sec, ',1pe11.4,' sec')

      return
      end
c-----------------------------------------------------------------------
      subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)
c-----------------------------------------------------------------------
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
      subroutine ifoption(ifif)        
      include 'SIZE'
      include 'EMWAVE'
      include 'PARALLEL'

      integer ifip
      logical ifif
 
      ifip =0
      if (ifif) ifip =1

      imaxg = iglmax(ifip,1)
      if (imaxg.eq.1) ifif =.true.
   
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
c-----------------------------------------------------------------------
      subroutine gh_face_extend(x,zg,n,gh_type,e,v)
c-----------------------------------------------------------------------
      include 'SIZE'

      real x(1),zg(1),e(1),v(1)
      integer gh_type

      if (ndim.eq.2) then
         call gh_face_extend_2d(x,zg,n,gh_type,e,v)
      else
         call gh_face_extend_3d(x,zg,n,gh_type,e,v)
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine gh_face_extend_2d(x,zg,n,gh_type,e,v)
c-----------------------------------------------------------------------
c
c     Extend 2D faces into interior via gordon hall
c
c     gh_type:  1 - vertex only
c               2 - vertex and faces
c
c
      real x(n,n)
      real zg(n)
      real e(n,n)
      real v(n,n)
      integer gh_type
c
c     Build vertex interpolant
c
      ntot=n*n
      call rzero(v,ntot)
      do jj=1,n,n-1
      do ii=1,n,n-1
         do j=1,n
         do i=1,n
            si     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            sj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            v(i,j) = v(i,j) + si*sj*x(ii,jj)
         enddo
         enddo
      enddo
      enddo
      if (gh_type.eq.1) then
         call copy(x,v,ntot)
         return
      endif


c     Extend 4 edges
      call rzero(e,ntot)
c
c     x-edges
c
      do jj=1,n,n-1
         do j=1,n
         do i=1,n
            hj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            e(i,j) = e(i,j) + hj*(x(i,jj)-v(i,jj))
         enddo
         enddo
      enddo
c
c     y-edges
c
      do ii=1,n,n-1
         do j=1,n
         do i=1,n
            hi     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            e(i,j) = e(i,j) + hi*(x(ii,j)-v(ii,j))
         enddo
         enddo
      enddo

      call add3(x,e,v,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine gh_face_extend_3d(x,zg,n,gh_type,e,v)
c-----------------------------------------------------------------------
c
c     Extend faces into interior via gordon hall
c
c     gh_type:  1 - vertex only
c               2 - vertex and edges
c               3 - vertex, edges, and faces
c
c
      real x(n,n,n)
      real zg(n)
      real e(n,n,n)
      real v(n,n,n)
      integer gh_type
c
c     Build vertex interpolant
c
      ntot=n*n*n
      call rzero(v,ntot)
      do kk=1,n,n-1
      do jj=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            si       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            sj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            sk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            v(i,j,k) = v(i,j,k) + si*sj*sk*x(ii,jj,kk)
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
      if (gh_type.eq.1) then
         call copy(x,v,ntot)
         return
      endif
c
c
c     Extend 12 edges
      call rzero(e,ntot)
c
c     x-edges
c
      do kk=1,n,n-1
      do jj=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            e(i,j,k) = e(i,j,k) + hj*hk*(x(i,jj,kk)-v(i,jj,kk))
         enddo
         enddo
         enddo
      enddo
      enddo
c
c     y-edges
c
      do kk=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            e(i,j,k) = e(i,j,k) + hi*hk*(x(ii,j,kk)-v(ii,j,kk))
         enddo
         enddo
         enddo
      enddo
      enddo
c
c     z-edges
c
      do jj=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            e(i,j,k) = e(i,j,k) + hi*hj*(x(ii,jj,k)-v(ii,jj,k))
         enddo
         enddo
         enddo
      enddo
      enddo
c
      call add2(e,v,ntot)
c
      if (gh_type.eq.2) then
         call copy(x,e,ntot)
         return
      endif
c
c     Extend faces
c
      call rzero(v,ntot)
c
c     x-edges
c
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            v(i,j,k) = v(i,j,k) + hi*(x(ii,j,k)-e(ii,j,k))
         enddo
         enddo
         enddo
      enddo
c
c     y-edges
c
      do jj=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            v(i,j,k) = v(i,j,k) + hj*(x(i,jj,k)-e(i,jj,k))
         enddo
         enddo
         enddo
      enddo
c
c     z-edges
c
      do kk=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            v(i,j,k) = v(i,j,k) + hk*(x(i,j,kk)-e(i,j,kk))
         enddo
         enddo
         enddo
      enddo
c
      call add2(v,e,ntot)
      call copy(x,v,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine map_to_crs(a,na,b,nb,if3d,w,ldw)
c-----------------------------------------------------------------------
c
c     Input:   b
c     Output:  a
c
      real a(1),b(1),w(1)
      logical if3d
c
      parameter(lx=40)
      real za(lx),zb(lx)
c
      real iba(lx*lx),ibat(lx*lx)
      save iba,ibat
c
      integer nao,nbo
      save    nao,nbo
      data    nao,nbo  / -9, -9/
c
      if (na.gt.lx.or.nb.gt.lx) then
         write(6,*)'ERROR: increase lx in map_to_crs to max:',na,nb
         call exitt
      endif
c
      if (na.ne.nao  .or.   nb.ne.nbo) then
         nao = na
         nbo = nb
         call zwgll(za,w,na)
         call zwgll(zb,w,nb)
         call igllm(iba,ibat,zb,za,nb,na,nb,na)
      endif
c
      call specmpn(a,na,b,nb,iba,ibat,if3d,w,ldw)

      return
      end
c-----------------------------------------------------------------------
