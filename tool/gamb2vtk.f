    
      PROGRAM SIMPLEX
    
      integer      lxi,ld,lbc
      parameter   (lxi= 10)
      parameter   (ld=1000)
      parameter   (lbc= 20)
      character*80 meshname
      character*80 meshname_new
      character*80 meshname_rea
      character*35 CBC          
      integer   numnp,nelem,ngrps,nbsets,ndfcd,ndfvl
      integer   nface1,nface2
      real      point(ld,3)
      integer   mapv (ld,4)
      integer   ibc1,ibc2,ibc3,ibc4
      integer   BCE(ld,lbc)
      integer   BCF(ld,lbc)

      call blank(meshname,80)
      call blank(meshname_new,80)

      write(*,*) 'type: name of gambit (neu) file'
      read(*,80)  meshname
      write(*,*) 'type: new filename for vtk'
      read(*,80)  meshname_new
c     write(*,*) 'type: new filename for rea'
c     read(*,80)  meshname_rea
  80  format(a80)

c... read neu file
      open(unit=80,file=meshname,status='old')

      read(80,*)              
      read(80,*)              
      read(80,*)             
      read(80,*)           
      read(80,*)               
      read(80,*)            
      read(80,*) numnp,nelem,ngrps,nbsets,ndfcd,ndfvl

      read(80,*)               
      read(80,*)            
      do i=1,numnp
         read(80,*) inum,point(i,1),point(i,2),point(i,3)
      enddo

      read(80,*)               
      read(80,*)            

      do i=1,nelem 
         read(80,*) inum,nface1,nface2,
     $              mapv(i,1),mapv(i,2),mapv(i,3),mapv(i,4)
      enddo
      read(80,*)               
      read(80,*)               
      read(80,*)            
      read(80,*)            
      read(80,*)            

      nelmode=10
      m=nelem/nelmode
      do i=1,m+1   
         read(80,*) 
      enddo

      do j=1,nbsets
      read(80,*)            
      read(80,*)            
      read(80,*)  CBC,ibc1,ibc2,ibc3,ibc4             

      do i=1,ibc2
         read(80,*) BCE(i,j),iedge,BCF(i,j)
      enddo
      read(80,*)            

      enddo

      close(80)
c... end of reading neu file

      open(unit=81,file=meshname_new,status='unknown')
      write(81,90) 
      write(81,91)
      write(81,92) 
      write(81,93)
      write(81,94) numnp

   90 format('# vtk DataFile Version 2.0')
   91 format('Unstructure Grid for TRI/TET')
   92 format('ASCII')                  
   93 format('DATASET UNSTRUCTURED_GRID')
   94 format('POINTS ',i8,' float')

      do i=1,numnp
         write(81,*) point(i,1),point(i,2),point(i,3)
      enddo

      write(81,95) nelem,nelem*5
   95 format('CELLS ',i10,i10)

      do i=1,nelem 
         write(81,*) 4,mapv(i,1)-1,mapv(i,2)-1,mapv(i,3)-1,mapv(i,4)-1
      enddo

      write(81,96) nelem                
   96 format('CELL_TYPES ',i10)

      ncelltype=10
      write(81,15) (ncelltype,i=1,nelem)
  15  format(70000i3)

      write(81,97)  numnp          
   97 format('POINT_DATA ',i10)

      write(81,98)                              
   98 format('VECTORS E float')            

      do i=1,numnp
         write(81,*) 1.0, 2.0, 3.0 
      enddo

      close(81)
      end

      subroutine blank(s,n)
      character*1 s(1)
      do i=1,n
        s(i)=' '
      enddo
      return
      end


      
