    
      PROGRAM SIMPLEX
    
      integer      lxi,ld,lbc
      parameter   (lxi= 10)
      parameter   (ld=1000)
      parameter   (lbc= 20)
      character*80 meshname
      character*80 meshname_map
      integer   numnp,nelem,ngrps,nbsets,ndfcd,ndfvl
      integer   nface1,nface2
      real      point(ld,3)
      integer   mapv (ld,4)
      integer   ibc1,ibc2,ibc3,ibc4
      integer   BCE(ld,lbc)
      integer   BCF(ld,lbc)
      character*35 RDB          
      character*3  BCC(ld,lbc)
      character*3  CBC(lbc)          
      character*3  BC          
      integer   ie,ifc,ibc

      call blank(meshname,80)
      call blank(meshname_map,80)

      write(*,*) 'type: name of gambit (neu) file'
      read(*,70)  meshname
      write(*,*) 'type: new filename for map'
      read(*,70)  meshname_map
  70  format(a80)

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

      close(80)
c... end of reading neu file

      open(unit=81,file=meshname_map,status='unknown')
      write(81,94) nelem, nelem, nelem,nelem,nelem,nelem,0  
   94 format(7i8)

      do i=1,nelem
         i1=mapv(i,1)
         i2=mapv(i,2)
         i3=mapv(i,3)
         i4=mapv(i,4)
         write(81,95) i, i1,i2,i3,i4 
      enddo
   95 format(5i10)

      close(81)
      end

      subroutine blank(s,n)
      character*1 s(1)
      do i=1,n
        s(i)=' '
      enddo
      return
      end


      
