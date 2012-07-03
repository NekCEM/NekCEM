      
      subroutine pvtk_xml2(dumpno)

      include 'SIZE'
      include 'EMWAVE'
      include 'PARALLEL'

      integer      unitno1
      parameter   (unitno1=211)
      integer      i, j, k, e, dumpno

      character*25 filename
      data         filename /'./vtk/ascii-NN-00000.pvtk'/
      CHARACTER*1  fnarray(25)
      equivalence (filename,fnarray)

      CHARACTER*29 vtkname1
      data         vtkname1 /'"ascii-NN-p000000-t00000.vtk"'/
      CHARACTER*1  vtkarray1(29)
      equivalence (vtkname1,vtkarray1)

      CHARACTER*29 vtkname2
      data         vtkname2 /'"ascii-NN-p000000-t00000.vtk"'/
      CHARACTER*1  vtkarray2(29)
      equivalence (vtkname2,vtkarray2)

      character*1  numrl(0:9)
      data         numrl /'0','1','2','3','4','5','6','7','8','9'/

      fnarray(20) = numrl(mod(dumpno,10))
      fnarray(19) = numrl(mod(dumpno,100)/10)
      fnarray(18) = numrl(mod(dumpno,1000)/100)
      fnarray(17) = numrl(mod(dumpno,10000)/1000)
      fnarray(16) = numrl(mod(dumpno,100000)/10000)

c     vtkarray1(18) = numrl(mod(dumpno,10))
c     vtkarray1(17) = numrl(mod(dumpno,100)/10)
c     vtkarray1(16) = numrl(mod(dumpno,1000)/100)
c     vtkarray1(15) = numrl(mod(dumpno,10000)/1000)
c     vtkarray1(14) = numrl(mod(dumpno,100000)/10000)
c     vtkarray1(11) = numrl(mod(nid,10))
c     vtkarray1(10) = numrl(mod(nid,100)/10)
c     vtkarray1( 9) = numrl(mod(nid,1000)/100)
c     vtkarray1( 8) = numrl(mod(nid,10000)/1000)
c     vtkarray1( 7) = numrl(mod(nid,100000)/10000)
c     vtkarray1( 6) = numrl(mod(nid,1000000)/100000)

      open (unitno1, file=filename)
      write(unitno1,*) '<File version="pvtk-1.0"'
      write(unitno1,*) ' dataType="vtkUnstructuredGrid"'
      write(unitno1,1)  np            

      do i=0,np-1    
         vtkarray2(24) = numrl(mod(dumpno,10))
         vtkarray2(23) = numrl(mod(dumpno,100)/10)
         vtkarray2(22) = numrl(mod(dumpno,1000)/100)
         vtkarray2(21) = numrl(mod(dumpno,10000)/1000)
         vtkarray2(20) = numrl(mod(dumpno,100000)/10000)
         vtkarray2(17) = numrl(mod(i,10))
         vtkarray2(16) = numrl(mod(i,100)/10)
         vtkarray2( 15) = numrl(mod(i,1000)/100)
         vtkarray2( 14) = numrl(mod(i,10000)/1000)
         vtkarray2( 13) = numrl(mod(i,100000)/10000)
         vtkarray2( 12) = numrl(mod(i,1000000)/100000)
         write(unitno1,2) vtkname2 
      enddo

      write(unitno1,3) 

  1   format ('   numberOfPieces=" ',i7,' " >')
  2   format ('   <Piece fileName=',A29,'/>')
  3   format (' </File>')
      
      close(unitno1)

      return
      end
c     -------------------------------------------------------------

      
      subroutine pvtk_xml3(dumpno)

      include 'SIZE'
      include 'EMWAVE'
      include 'PARALLEL'

      integer      unitno1
      parameter   (unitno1=211)
      integer      i, j, k, e, dumpno

      character*26 filename
      data         filename /'./vtk/binary-NN-00000.pvtk'/
      CHARACTER*1  fnarray(26)
      equivalence (filename,fnarray)

      CHARACTER*30 vtkname1
      data         vtkname1 /'"binary-NN-p000000-t00000.vtk"'/
      CHARACTER*1  vtkarray1(30)
      equivalence (vtkname1,vtkarray1)

      CHARACTER*30 vtkname2
      data         vtkname2 /'"binary-NN-p000000-t00000.vtk"'/
      CHARACTER*1  vtkarray2(30)
      equivalence (vtkname2,vtkarray2)

      character*1  numrl(0:9)
      data         numrl /'0','1','2','3','4','5','6','7','8','9'/

      fnarray(21) = numrl(mod(dumpno,10))
      fnarray(20) = numrl(mod(dumpno,100)/10)
      fnarray(19) = numrl(mod(dumpno,1000)/100)
      fnarray(18) = numrl(mod(dumpno,10000)/1000)
      fnarray(17) = numrl(mod(dumpno,100000)/10000)

c     vtkarray1(18) = numrl(mod(dumpno,10))
c     vtkarray1(17) = numrl(mod(dumpno,100)/10)
c     vtkarray1(16) = numrl(mod(dumpno,1000)/100)
c     vtkarray1(15) = numrl(mod(dumpno,10000)/1000)
c     vtkarray1(14) = numrl(mod(dumpno,100000)/10000)
c     vtkarray1(11) = numrl(mod(nid,10))
c     vtkarray1(10) = numrl(mod(nid,100)/10)
c     vtkarray1( 9) = numrl(mod(nid,1000)/100)
c     vtkarray1( 8) = numrl(mod(nid,10000)/1000)
c     vtkarray1( 7) = numrl(mod(nid,100000)/10000)
c     vtkarray1( 6) = numrl(mod(nid,1000000)/100000)

      open (unitno1, file=filename)
      write(unitno1,*) '<File version="pvtk-1.0"'
      write(unitno1,*) ' dataType="vtkUnstructuredGrid"'
      write(unitno1,1)  np            

      do i=0,np-1    
         vtkarray2(25) = numrl(mod(dumpno,10))
         vtkarray2(24) = numrl(mod(dumpno,100)/10)
         vtkarray2(23) = numrl(mod(dumpno,1000)/100)
         vtkarray2(22) = numrl(mod(dumpno,10000)/1000)
         vtkarray2(21) = numrl(mod(dumpno,100000)/10000)
         vtkarray2(18) = numrl(mod(i,10))
         vtkarray2(17) = numrl(mod(i,100)/10)
         vtkarray2( 16) = numrl(mod(i,1000)/100)
         vtkarray2( 15) = numrl(mod(i,10000)/1000)
         vtkarray2( 14) = numrl(mod(i,100000)/10000)
         vtkarray2( 13) = numrl(mod(i,1000000)/100000)
         write(unitno1,2) vtkname2 
      enddo

      write(unitno1,3) 

  1   format ('   numberOfPieces=" ',i7,' " >')
  2   format ('   <Piece fileName=',A30,'/>')
  3   format (' </File>')
      
      close(unitno1)

      return
      end
c     -------------------------------------------------------------
      subroutine vtk_dump_geometry()
      include 'SIZE'
      include 'TOTAL'

      integer unitno
      parameter (unitno=17)
      integer i, j, k, e

      character*22 filename
      data filename /'out/vtkgeom-n00000.dat'/
      CHARACTER*1  fnarray(22)
      equivalence (filename,fnarray)

      character*1  numrl(0:9)
      data numrl /'0','1','2','3','4','5','6','7','8','9'/

      fnarray(18) = numrl(mod(nid,10))
      fnarray(17) = numrl(mod(nid,100)/10)
      fnarray(16) = numrl(mod(nid,1000)/100)
      fnarray(15) = numrl(mod(nid,10000)/1000)
      fnarray(14) = numrl(mod(nid,100000)/10000)

      open  (unitno, file=filename)
      write (unitno,*) nelt, nx1, ny1, nz1
      
      do e=1,nelt
      do k=1, nz1
      do j=1, ny1
      do i=1, nx1
        
        write (unitno,10) xm1(i,j,k,e), ym1(i,j,k,e), zm1(i,j,k,e)

      enddo
      enddo
      enddo
      enddo

10    format (3d20.10)

      close(unitno)

      return

      end
c     -------------------------------------------------------------
      subroutine vtk_dump_field(field, id, step, tt)
      include 'SIZE'
      include 'TOTAL'

      real        field(lx1,ly1,lz1,lelt)
      character*6 id
      integer     step
      real        tt

      integer     unitno 
      parameter  (unitno=39)

      logical     icalld
      save        icalld            
      data        icalld /.false./

      integer     dumpno
      save        dumpno
      data        dumpno /0/

      character*26 filename
      data filename /'out/vtk-p000000-t00000.dat'/
      CHARACTER*1  fnarray(26)
      equivalence (filename,fnarray)
      
      integer i,j,k,e

      character*1  numrl(0:9)
      data numrl /'0','1','2','3','4','5','6','7','8','9'/

      if (.not.icalld) then
        call vtk_dump_geometry
        icalld = .true.
      endif
      
      fnarray(22) = numrl(mod(dumpno,10))
      fnarray(21) = numrl(mod(dumpno,100)/10)
      fnarray(20) = numrl(mod(dumpno,1000)/100)
      fnarray(19) = numrl(mod(dumpno,10000)/1000)
      fnarray(18) = numrl(mod(dumpno,100000)/10000)

      fnarray(15) = numrl(mod(nid,10))
      fnarray(14) = numrl(mod(nid,100)/10)
      fnarray(13) = numrl(mod(nid,1000)/100)
      fnarray(12) = numrl(mod(nid,10000)/1000)
      fnarray(11) = numrl(mod(nid,100000)/10000)
      fnarray(10) = numrl(mod(nid,1000000)/100000)

      open(unitno, file=filename)
      write (unitno,20) id, step, tt
20    format (A,i10,d20.10)

      do e=1, nelt
      do k=1, nz1
      do j=1, ny1
      do i=1, nx1
        
        write (unitno,30) field(i,j,k,e)
30      format (d20.10)

      enddo
      enddo
      enddo
      enddo

      dumpno = dumpno + 1

      close(unitno)

      return

      end

c     -------------------------------------------------------------
      subroutine vtk_dump_geometry2(unitno)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer  unitno
      integer  i, j, k, e, f

      character*18 filename
      data         filename /'./vtk/em-cells.vtk'/
      CHARACTER*1  fnarray(18)
      equivalence (filename,fnarray)

      integer nfpcnt (6,0:lxi*2*ndim)
      integer nfpcnt0(6,0:lxi*2*ndim)

      if (IF3D) then
        npts  = nx1*ny1*nz1*nelt
        nfpts = nx1*2+(nx1-2)*2
        nfgap = nx1-1
        nfline= nelt*6
      else
        npts  = nx1*ny1*nz1*nelt
        nfpts = nx1*2+(nx1-2)*2
        nfgap = nx1-1
        nfline= nelt*1
      endif                               

c ... writing in VTK format

      open (unitno, file=filename)

      write(unitno,1)
      write(unitno,2)
      write(unitno,3)
      write(unitno,4)
      write(unitno,5) npts                     

  1   format ('# vtk DataFile Version 2.0')
  2   format ('Cell geometry')
  3   format ('ASCII')
  4   format ('DATASET UNSTRUCTURED_GRID')
c This gives errors w/ visit and paraview, when having the values
c overranging the velues E-45 or something like that 
c 5   format ('POINTS ',i7,' float') 
  5   format ('POINTS ',i7,'double')

      const = 1.0 !xyzscale
      do i=1,nx1*ny1*nz1*nelt
        xx = xm1(i,1,1,1)*const 
        yy = ym1(i,1,1,1)*const 
        zz = zm1(i,1,1,1)*const 
        write (unitno,6) xx, yy, zz 
      enddo

  6   format (3e15.10)

      write(unitno,7) nfline,(nfpts+1)*nfline                    

      
      if (IF3D) then

          do i=1,nfgap 

             ! face 1
             nfpcnt(1,0*nfgap+i)=i                        
             nfpcnt(1,1*nfgap+i)=nfgap+(nfgap+1)*i        
             nfpcnt(1,2*nfgap+i)=(nfgap+1)*(nfgap+1)-1-i        
             nfpcnt(1,3*nfgap+i)=(nfgap+1)*(nfgap)-(nfgap+1)*i      

             ! face 6 
             nfgaptmp = (nfgap+1)*(nfgap+1)*nfgap        
             nfpcnt(2,0*nfgap+i)=nfpcnt(1,0*nfgap+i)+nfgaptmp          
             nfpcnt(2,1*nfgap+i)=nfpcnt(1,1*nfgap+i)+nfgaptmp        
             nfpcnt(2,2*nfgap+i)=nfpcnt(1,2*nfgap+i)+nfgaptmp        
             nfpcnt(2,3*nfgap+i)=nfpcnt(1,3*nfgap+i)+nfgaptmp

             ! face 2 
             nfpcnt(3,0*nfgap+i)=(nfgap+1)*i                  
             nfpcnt(3,1*nfgap+i)=(nfgap+1)*(nfgap+1)*i+(nfgap+1)*nfgap 
             nfpcnt(3,2*nfgap+i)=(nfgap+1)*(nfgap+1)*nfgap
     $                          +(nfgap-i)*(nfgap+1)     
             nfpcnt(3,3*nfgap+i)=(nfgap+1)*(nfgap+1)*(nfgap-i)

             ! face 3 
             nfpcnt(4,0*nfgap+i)=nfpcnt(3,0*nfgap+i)+nfgap             
             nfpcnt(4,1*nfgap+i)=nfpcnt(3,1*nfgap+i)+nfgap      
             nfpcnt(4,2*nfgap+i)=nfpcnt(3,2*nfgap+i)+nfgap        
             nfpcnt(4,3*nfgap+i)=nfpcnt(3,3*nfgap+i)+nfgap

             ! face 4 
             nfpcnt(5,0*nfgap+i)=i                  
             nfpcnt(5,1*nfgap+i)=(nfgap+1)*(nfgap+1)*i+nfgap 
             nfpcnt(5,2*nfgap+i)=(nfgap+1)*(nfgap+1)*nfgap+(nfgap-i) 
             nfpcnt(5,3*nfgap+i)=(nfgap+1)*(nfgap+1)*(nfgap-i)

             ! face 5 
             nfgaptmp = (nfgap+1)*nfgap        
             nfpcnt(6,0*nfgap+i)=nfpcnt(5,0*nfgap+i)+nfgaptmp         
             nfpcnt(6,1*nfgap+i)=nfpcnt(5,1*nfgap+i)+nfgaptmp   
             nfpcnt(6,2*nfgap+i)=nfpcnt(5,2*nfgap+i)+nfgaptmp     
             nfpcnt(6,3*nfgap+i)=nfpcnt(5,3*nfgap+i)+nfgaptmp

          enddo

          nfpcnt(1,0)=nfpcnt(1,nfgap*4)
          nfpcnt(2,0)=nfpcnt(2,nfgap*4)
          nfpcnt(3,0)=nfpcnt(3,nfgap*4)
          nfpcnt(4,0)=nfpcnt(4,nfgap*4)
          nfpcnt(5,0)=nfpcnt(5,nfgap*4)
          nfpcnt(6,0)=nfpcnt(6,nfgap*4)

          do f=1,6
          do i=0,4*nfgap-1 
             nfpcnt0(f,i)=nfpcnt(f,i)                    
          enddo
          enddo

          do e = 1,nelt

             nfgaptmp =(nfgap+1)*(nfgap+1)*(nfgap+1)*(e-1)
             do f=1,6
             do i =0,4*nfgap-1 
                nfpcnt(f,i)=nfpcnt0(f,i)+nfgaptmp                   
             enddo
             enddo

             write(unitno,8) nfpts,(nfpcnt(1,i),i=0,4*nfgap-1)      
             write(unitno,8) nfpts,(nfpcnt(2,i),i=0,4*nfgap-1)        
             write(unitno,8) nfpts,(nfpcnt(3,i),i=0,4*nfgap-1)        
             write(unitno,8) nfpts,(nfpcnt(4,i),i=0,4*nfgap-1)        
             write(unitno,8) nfpts,(nfpcnt(5,i),i=0,4*nfgap-1)        
             write(unitno,8) nfpts,(nfpcnt(6,i),i=0,4*nfgap-1)        

          enddo

      else

        do e=1,nelt
        do i=1,nfgap 
           i2=(e-1)*nx1*nx1          
           nfpcnt(1,0*nfgap+i)=i2+i                        
           nfpcnt(1,1*nfgap+i)=i2+nfgap+(nfgap+1)*i        
           nfpcnt(1,2*nfgap+i)=i2+(nfgap+1)*(nfgap+1)-1-i        
           nfpcnt(1,3*nfgap+i)=i2+(nfgap+1)*(nfgap)-(nfgap+1)*i        
        enddo
           nfpcnt(1,0)=nfpcnt(1,nfgap*4)
           write(unitno,8) nfpts, (nfpcnt(1,i),i=0,4*nfgap-1)           
        enddo

      endif

      write(unitno,9) nfline

  7   format ('CELLS',2i7)
  8   format (13i7)
  9   format ('CELL_TYPES',i7)

      ncelldata = 7

      write (unitno,10) (ncelldata,i=1,nfline)
  10  format (70000i3)        

      close(unitno)

c...  end writing in VTK format

      return
      end

c     -------------------------------------------------------------
      subroutine vtk_dump_header2(unitno)                            
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer e, i, j, k, f
      integer unitno
      integer npts,ncell,ncellsize,ncelltype
      integer nfpcnt(1:2+2*ndim)

      logical icalld
      save    icalld                
      data    icalld /.false./

      if (IF3D) then
          npts =  nx1*ny1*nz1*nelt
          ncell= (nx1-1)*(ny1-1)*(nz1-1)*nelt
          ncellsize= ncell*(4*(ndim-1)+1)
          ncelltype= 12           
      else
          npts =  nx1*ny1*nelt
          ncell= (nx1-1)*(ny1-1)*nelt
          ncellsize= ncell*(4*(ndim-1)+1)
          ncelltype= 9 
      endif

c...  writing in VTK format
     
      write(unitno,1)
      write(unitno,2)
      write(unitno,3)
      write(unitno,4)
      write(unitno,5) npts                       

  1   format ('# vtk DataFile Version 2.0')
  2   format ('EM data')
  3   format ('ASCII')
  4   format ('DATASET UNSTRUCTURED_GRID')
  5   format ('POINTS ',i12,' float')

      const= 1.0 !xyzscale
      do i=1,nx1*ny1*nz1*nelt
         xx=const*xm1(i,1,1,1)
         yy=const*ym1(i,1,1,1)
         zz=const*zm1(i,1,1,1)
         write (unitno,30) xx, yy, zz
      enddo
 30   format (3e18.8)
 
      write(unitno, 6)  ncell, ncellsize                  
  6   format ('CELLS      ',2i12)

      nfnum = 4*(ndim-1)
      nfgap = nx1  
      nfcnt = 0

      if (IF3D) then

        do e=1,nelt
        do j=1,(nfgap-1)                                  
           nfcnt = (j-1)*nfgap*nfgap+(e-1)*nfgap*nfgap*nfgap            
           do i=1,(nfgap-1)*(nfgap-1)                  
             nfpcnt(1)=nfcnt
             nfpcnt(2)=nfcnt+1
             nfpcnt(3)=nfcnt+nfgap+1
             nfpcnt(4)=nfcnt+nfgap     
             nfpcnt(5)=nfpcnt(1)+nfgap*nfgap     
             nfpcnt(6)=nfpcnt(2)+nfgap*nfgap    
             nfpcnt(7)=nfpcnt(3)+nfgap*nfgap     
             nfpcnt(8)=nfpcnt(4)+nfgap*nfgap     
             write (unitno,7) nfnum, (nfpcnt(k),k=1,nfnum)  
             nfcnt = nfcnt+1 
             if (mod(nfcnt+1,nfgap).eq.0) nfcnt = nfcnt+1 
           enddo
        enddo
        enddo

      else

        do e=1,nelt                
           nfcnt = (e-1)*nfgap*nfgap
           do i=1,(nfgap-1)*(nfgap-1)
             nfpcnt(1)=nfcnt
             nfpcnt(2)=nfcnt+1
             nfpcnt(3)=nfcnt+nfgap+1
             nfpcnt(4)=nfcnt+nfgap     
             write (unitno,7) nfnum, (nfpcnt(k),k=1,nfnum)  
             nfcnt = nfcnt+1 
             if (mod(nfcnt+1,nfgap).eq.0) nfcnt = nfcnt+1 
           enddo
        enddo

      endif

      write(unitno, 8)  ncell                                
      write(unitno, 9) (ncelltype,i=1,ncell)
      write(unitno,10)  npts                          
  7   format ( 9i7 )
  8   format ('CELL_TYPES ',i7  )
  9   format (70000i3 )
 10   format ('POINT_DATA ',i7 )

      return
      end

c     -------------------------------------------------------------
      subroutine vtk_dump_field2(nunit,field1,ux1,uy1,uz1,step,tt)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer nunit
      real    ux1(lx1,ly1,lz1,lelt)
      real    uy1(lx1,ly1,lz1,lelt)
      real    uz1(lx1,ly1,lz1,lelt)

      real    tt

      character*6 field1
      integer     step

      write (nunit,10) field1              
10    format ('VECTORS ',A,' float')
      
      nxyz = nx1*ny1*nz1*nelt

      do i=1,nxyz                
         write (nunit,30) ux1(i,1,1,1),uy1(i,1,1,1),uz1(i,1,1,1) 
      enddo
30    format (3e18.8)

      return

      end
c     -------------------------------------------------------------
      subroutine vtk_dump_header_swap()
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer idumpno
      integer e, i, j, k, f, eg
      integer unitno
      common  ncell,ncellsize,ncelltype
      integer ncell,ncellsize,ncelltype    

      common  /s_rvtk/ vtkxyz_s(3*lx1*ly1*lz1*lelt)
      real             vtkxyz_s
      common  /s_ivtk/ fpcnt_s (8*lx1*ly1*lz1*lelt)
      integer          fpcnt_s 

      npts  = nx1*ny1*nz1*nelt
      nxyz  = nx1*ny1*nz1

      !call writeheader4(istep,idumpno,nx1,np,time,dt)

      if (IF3D) then 
          ncell= (nx1-1)*(ny1-1)*(nz1-1)*nelt
          ncellsize= ncell*(4*(ndim-1)+1)
          ncelltype= 12   
      else
          ncell= (nx1-1)*(ny1-1)*nelt
          ncellsize= ncell*(4*(ndim-1)+1)
          ncelltype= 9
      endif

c...  writing in VTK format

      const = 1.0 !xyzscale

      k = 1
      do i = 1, npts
         vtkxyz_s (k  )=const*xm1(i,1,1,1)
         vtkxyz_s (k+1)=const*ym1(i,1,1,1)
         vtkxyz_s (k+2)=const*zm1(i,1,1,1)
         k = k+3
      enddo

      m = 3*nxyz
      call swap_real ( vtkxyz_s , m )

      !call writenodes4(vtkxyz,npts)     

      nfnum = 4*(ndim-1)
      nfgap = nx1
      nfcnt = 0

      if (IF3D) then

        k = 0
        do e=1,nelt
         eg = lglel(e)
         do j=1,(nfgap-1)
           nfcnt = (j-1)*nfgap*nfgap+(eg-1)*nfgap*nfgap*nfgap
           do i=1,(nfgap-1)*(nfgap-1)
             fpcnt_s(k+1)=nfcnt
             fpcnt_s(k+2)=nfcnt+1
             fpcnt_s(k+3)=nfcnt+nfgap+1
             fpcnt_s(k+4)=nfcnt+nfgap
             fpcnt_s(k+5)=fpcnt_s(k+1)+nfgap*nfgap
             fpcnt_s(k+6)=fpcnt_s(k+2)+nfgap*nfgap
             fpcnt_s(k+7)=fpcnt_s(k+3)+nfgap*nfgap
             fpcnt_s(k+8)=fpcnt_s(k+4)+nfgap*nfgap
             nfcnt = nfcnt+1
             if (mod(nfcnt+1,nfgap).eq.0) nfcnt = nfcnt+1
             k = k+8
           enddo
         enddo
        enddo

        if (k.gt.lx1*lx1*lx1*8*lelt) then
          write(6,*) 'ERR:: fpcnt size is small in vtk_dump_header'
          call exitt
        endif

        !if (nid.eq.0) write(6,*) 'dumpvtk.F: b4 swap_even_ints'
        m = 8
        call swap_even_ints(fpcnt_s, m) ! m must be even 
        !call write3dcells4(fpcnt,nelt,ncell,npts)                  

      else ! 2D

        k = 0
        do e=1,nelt
           eg = lglel(e)
           nfcnt = (eg-1)*nfgap*nfgap
           do i=1,(nfgap-1)*(nfgap-1)
             fpcnt_s(k+1)=nfcnt
             fpcnt_s(k+2)=nfcnt+1
             fpcnt_s(k+3)=nfcnt+nfgap+1
             fpcnt_s(k+4)=nfcnt+nfgap
             nfcnt = nfcnt+1
             if (mod(nfcnt+1,nfgap).eq.0) nfcnt = nfcnt+1
             k = k+4
           enddo
        enddo

        if (k.gt.lx1*lx1*4*lelt) then
          write(6,*) 'ERR:: fpcnt size is small in vtk_dump_header'
          call exitt
        endif

        m = 4
        call swap_even_ints(fpcnt_s, m) ! m must be even 
        k=1

        !call write2dcells4(fpcnt,nelt,ncell,npts)                   

      endif

      return
      end
c     -------------------------------------------------------------
      subroutine vtk_dump_header_nonswap()
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer e, i, j, k, f
      integer unitno

      common  ncell,ncellsize,ncelltype
      integer ncell,ncellsize,ncelltype

      common  /c_rvtk/ vtkxyz(3*lx1*ly1*lz1*lelt)
      real             vtkxyz
      common  /c_ivtk/ fpcnt (8*lx1*ly1*lz1*lelt)
      integer          fpcnt

      nxyz  = nx1*ny1*nz1*nelt

c...  calculate sizes of arrays for total points and cells
      if (IF3D) then
          npts =  nx1*ny1*nz1*nelt
          ncell= (nx1-1)*(ny1-1)*(nz1-1)*nelt
          ncellsize= ncell*(4*(ndim-1)+1)
          ncelltype= 12
      else
          npts =  nx1*ny1*nelt
          ncell= (nx1-1)*(ny1-1)*nelt
          ncellsize= ncell*(4*(ndim-1)+1)
          ncelltype= 9
      endif

c...  save coordinate
      const = xyzscale
      k = 1
      do i = 1, nxyz
         vtkxyz (k  )=const*xm1(i,1,1,1)
         vtkxyz (k+1)=const*ym1(i,1,1,1)
         vtkxyz (k+2)=const*zm1(i,1,1,1)
         k = k+3
      enddo

      nfnum = 4*(ndim-1)
      nfgap = nx1
      nfcnt = 0

c...  compute cell numbering
      if (IF3D) then

        k = 0
        do e=1,nelt
        do j=1,(nfgap-1)
           nfcnt = (j-1)*nfgap*nfgap+(e-1)*nfgap*nfgap*nfgap
           do i=1,(nfgap-1)*(nfgap-1)
             fpcnt(k+1)=nfcnt
             fpcnt(k+2)=nfcnt+1
             fpcnt(k+3)=nfcnt+nfgap+1
             fpcnt(k+4)=nfcnt+nfgap
             fpcnt(k+5)=fpcnt(k+1)+nfgap*nfgap
             fpcnt(k+6)=fpcnt(k+2)+nfgap*nfgap
             fpcnt(k+7)=fpcnt(k+3)+nfgap*nfgap
             fpcnt(k+8)=fpcnt(k+4)+nfgap*nfgap
             nfcnt = nfcnt+1
             if (mod(nfcnt+1,nfgap).eq.0) nfcnt = nfcnt+1
             k = k+8
           enddo
        enddo
        enddo

        if (k.gt.lx1*lx1*lx1*8*lelt) then
          write(6,*) 'ERR:: fpcnt size is small in vtk_dump_header'
          call exitt
        endif

       !call write3dcells6(fpcnt,nelt,ncell,npts)                      

      else
        k = 0
        do e=1,nelt
           nfcnt = (e-1)*nfgap*nfgap
           do i=1,(nfgap-1)*(nfgap-1)
             fpcnt(k+1)=nfcnt
             fpcnt(k+2)=nfcnt+1
             fpcnt(k+3)=nfcnt+nfgap+1
             fpcnt(k+4)=nfcnt+nfgap
             nfcnt = nfcnt+1
             if (mod(nfcnt+1,nfgap).eq.0) nfcnt = nfcnt+1
             k = k+4
           enddo
        enddo

        if (k.gt.lx1*lx1*4*lelt) then
          write(6,*) 'ERR:: fpcnt size is small in vtk_dump_header'
          call exitt
        endif

       !call write2dcells6(fpcnt,nelt,ncell,npts)                      

      endif

      return
      end
c     -------------------------------------------------------------
      subroutine vtk_swap_field(ux1,uy1,uz1,vtkout)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real ux1(lx1,ly1,lz1,lelt)
      real uy1(lx1,ly1,lz1,lelt)
      real uz1(lx1,ly1,lz1,lelt)
      real vtkout(3*lx1*ly1*lz1*lelt)


      npts = nx1*ny1*nz1*nelt
      nxyz = nx1*ny1*nz1
      nxyz3= nx1*ny1*nz1*3

      k = 1
      do i = 1, npts
         vtkout(k  )=ux1(i,1,1,1)
         vtkout(k+1)=uy1(i,1,1,1)
         vtkout(k+2)=uz1(i,1,1,1)
         k = k+3
      enddo

      call swap_real   (vtkout, nxyz3)

      return
      end
c     -------------------------------------------------------------
      subroutine vtk_nonswap_field(ux1,uy1,uz1,vtkout)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      real ux1(lx1,ly1,lz1,lelt)
      real uy1(lx1,ly1,lz1,lelt)
      real uz1(lx1,ly1,lz1,lelt)
      real vtkout(3*lx1*ly1*lz1*lelt)

      integer i, npts
      npts = nx1*ny1*nz1*nelt

      k = 1
      do i = 1, npts
         vtkout(k  )=ux1(i,1,1,1)
         vtkout(k+1)=uy1(i,1,1,1)
         vtkout(k+2)=uz1(i,1,1,1)
         k = k+3
      enddo

      return
      end
c     -------------------------------------------------------------
      subroutine tswap_ip(x,m,p,n,w) ! in-place permutation: x' = x(p)
      real    x(m,n),w(m)
      integer p(n)
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            call copy(w,x(1,k),m)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! tswap_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  call copy(x(1,last),w,m)
                  p(last) = -p(last)
                  goto 10
               else
                  call copy(x(1,last),x(1,next),m)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine get_target_proc_list(list3)

      include 'SIZE'
      include 'TOTAL'
      common /nekmpi/ nidx,npx,nekcomm,nekgroup,nekreal

      integer list3(3,lelt)
      integer e,e0,e1,eg,eg_desired,el,cr_t

      common /rstswap/ plist(2,lelt),qlist(2,lelt),cr_s
      integer plist,qlist,cr_s

      integer icalld
      save    icalld
      data    icalld /0/

      if (icalld.eq.1) return
      icalld = 1

      e1 = igl_running_sum(nelt)
      e0 = e1 - nelt

      call crystal_new(cr_t,nekcomm,np)

c     Algorithm:  three cr calls:

c     1&2: Swap current element list and desired element list

c       3: Match lists and return target processor



c     FIRST CALL:
      do e=1,nelt
         eg = lglel(e)
         plist(1,e) = mod(eg,np) ! target proc for swap 1
         plist(2,e) = eg
      enddo

      kp = 1      ! Row # that holds the target processor #
      n1 = nelt
      call crystal_ituple_transfer(cr_h,plist,2,n1,lelt,kp)

c     SECOND CALL:
      do e=1,nelt
         eg_desired = e0+e  ! Consecutive range associated with proc nid.
         qlist(1,e) = mod(eg_desired,np) ! target proc for swap 1
         qlist(2,e) = eg_desired
      enddo

      kp = 1      ! Row # that holds the target processor #
      n2 = nelt
      call crystal_ituple_transfer(cr_h,qlist,2,n2,lelt,kp)

c     locally sort columns list(:,1...n) in ascending order
      key  = 2
      nkey = 1
      call crystal_ituple_sort(cr_h,plist,2,n1,key,nkey)
      call crystal_ituple_sort(cr_h,qlist,2,n2,key,nkey)

      do e=1,n2
         list3(1,e) = plist(1,e)  ! Target processor id
         list3(2,e) = qlist(1,e)  ! Target processor id
         list3(3,e) = qlist(2,e)  ! Target processor el
      enddo


      kp = 1      ! Row # that holds the target processor #
      n2 = nelt
      call crystal_ituple_transfer(cr_h,list3,3,n1,lelt,kp)

      do e=1,nelt
         mid = list3(2,e)
         eg  = list3(3,e)
         el  = gllel(eg)
         plist(1,el) = mid ! Target processor 
         plist(2,el) = eg  ! Target processor 
      enddo

      !call crystal_done(cr_t)

      return
      end
c     -------------------------------------------------------------
      subroutine swap_even_ints ( idata , m )

c     ! Only swaps tuples of even integers m

      include 'SIZE'
      include 'TOTAL'
      common /nekmpi/ nidx,npx,nekcomm,nekgroup,nekreal

      integer idata(m,lelt)

      common /rstswap/ plist(2,lelt),qlist(2,lelt),cr_s
      common /rstswap2/ qlist2(lelt)
      integer plist,qlist,cr_s,qlist2

      parameter(ldw=lx1*ly1*lz1*lelt)
      common /ctmp0/ wk1(ldw),wk2(lelt)

      integer icalld
      save    icalld
      data    icalld /0/

      if (m.gt.ldw) call exitti('Increase ldw in swap_real$',m)
      if (mod(m,2).gt.0) 
     $     call exitti('ERROR: swap_even_ints, m not even.$',m)

      if (icalld.eq.0) then
         icalld = 1
         call get_target_proc_list( wk1 )
         call crystal_new(cr_s,nekcomm,np)
      endif

      kp = 1  ! Key for target processor list
      n  = nelt
      call icopy(qlist,plist,2*n)

      m2 = m
      if (wdsize.eq.8) m2=m/2
      call crystal_transfer
     $             (cr_s,n,lelt,qlist,2,dummy,0,idata,m2,kp)


c     NEXT:  Sort data into element-ordered list.

      do i=1,n
         qlist(i,1) = qlist(2,i)         ! compact qlist
      enddo

      call isort   (qlist,qlist2,n)
      call tswap_ip(idata,m2,qlist2,n,wk1) ! in-place permutation: x' = x(p)

      return
      end
c     -------------------------------------------------------------
      subroutine swap_real ( rdata , m )

      include 'SIZE'
      include 'TOTAL'
      common /nekmpi/ nidx,npx,nekcomm,nekgroup,nekreal

      real rdata(m,lelt)

      common /rstswap/ plist(2,lelt),qlist(2,lelt),cr_s
      common /rstswap2/ qlist2(lelt)
      integer plist,qlist,cr_s,qlist2
      integer e,eg

      parameter(ldw=lx1*ly1*lz1*lelt)
      common /ctmp0/ wk1(ldw),wk2(lelt)

      integer icalld
      save    icalld
      data    icalld /0/

      if (m.gt.ldw) call exitti('Increase ldw in swap_real$',m)

      if (icalld.eq.0) then
         icalld = 1
         call get_target_proc_list( wk1 )
         call crystal_new(cr_s,nekcomm,np)
      endif

      kp = 1  ! Key for target processor list
      n  = nelt
      call icopy(qlist,plist,2*n)

      call crystal_transfer
     $             (cr_s,n,lelt,qlist,2,dummy,0,rdata,m,kp)

c     NEXT:  Sort data into element-ordered list.

      do i=1,n
         qlist(i,1) = qlist(2,i)         ! compact qlist
      enddo

      call isort   (qlist,qlist2,n)
      call tswap_ip(rdata,m,qlist2,n,wk1) ! in-place permutation: x' = x(p)

      return
      end
c     -------------------------------------------------------------
      subroutine swap_real_backward ( rdata , m )
      include 'SIZE'
      include 'TOTAL'
      common /nekmpi/ nidx,npx,nekcomm,nekgroup,nekreal

      real    rdata(m,lelt)

      common /rstswap/ plist(2,lelt),qlist(2,lelt),cr_s
      common /rstswap2/ qlist2(lelt)
      integer plist,qlist,cr_s,qlist2

      parameter(ldw=lx1*ly1*lz1*lelt)
      common /ctmp0/ wk1(ldw),wk2(lelt)

      integer e0,e1,eg,e

      if (m.gt.ldw) call exitti('Increase ldw in swap_real$',m)

      call crystal_new(cr_s,nekcomm,np)

      e1 = igl_running_sum(nelt)
      e0 = e1 - nelt

      do e=1,nelt
         eg=e0+e
         qlist(1,e) = gllnid(eg)
         qlist(2,e) = gllel (eg)
      enddo

      kp= 1      ! Row # that holds the target processor #
      n = nelt
      call crystal_transfer
     $             (cr_s,n,lelt,qlist,2,dummy,0,rdata,m,kp)

      do i=1,n
         qlist(i,1) = qlist(2,i)         ! compact qlist
      enddo

      call isort    (qlist,qlist2,n)
      call tswap_ip (rdata,m,qlist2,n,wk1) ! in-place permutation: x' = x(p)

      !if (nid.eq.0) write(6,*) 'swap_real_backward: done'

      call crystal_done(cr_s)              

      return
      end
c     -------------------------------------------------------------
