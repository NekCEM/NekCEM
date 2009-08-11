    
      PROGRAM SIMPLEX
    
      integer      lxi,ld,lbc
      parameter   (lxi= 10)
      parameter   (ld=1000)
      parameter   (lbc= 20)
      character*80 meshname
      character*80 meshname_new
      character*80 meshname_rea
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
      call blank(meshname_new,80)
      call blank(CBC(1),3)
      call blank(CBC(2),3)
      call blank(CBC(3),3)

      write(*,*) 'type: name of gambit (neu) file'
      read(*,70)  meshname
      write(*,*) 'type: new filename for rea'
      read(*,70)  meshname_rea
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

      nelmode=10
      m= nelem/nelmode
      do i=1,m+1   
         read(80,*) 
      enddo

c ... set boundary conditions
      write(*,*) 'type BC1'
      read (*,71)  CBC(1)             
      write(*,*) 'type BC2'
      read (*,71)  CBC(2)             
      write(*,*) 'type BC3'
      read (*,71)  CBC(3)             
  71  format(a3)
      !write(6,*) 'CBC',CBC(1),CBC(2),CBC(3)

c ... set=1st boundary 
      ib1=1
      ib2=2
      ib3=3
      read(80,*)            
      read(80,*)            
      read(80,*)  RDB,ibc1,ibc2,ibc3,ibc4             
      write(*,*)  RDB,ibc1,ibc2,ibc3,ibc4             
      do i=1,ibc2
         read(80,*) BCE(i,ib1),iedge,BCF(i,ib1)
         BCC(i,ib1)= CBC(ib1)
      enddo
c ... set=2nd boundary 
      read(80,*)            
      read(80,*)            
      read(80,*)  RDB,ibc1,ibc2,ibc3,ibc4             
      write(*,*)  RDB,ibc1,ibc2,ibc3,ibc4             
      do i=1,ibc2
         read(80,*) BCE(i,ib2),iedge,BCF(i,ib2)
         BCC(i,ib2)= CBC(ib2)
      enddo
c ... set=3rd boundary 
      read(80,*)            
      read(80,*)            
      read(80,*)  RDB,ibc1,ibc2,ibc3,ibc4             
      write(*,*)  RDB,ibc1,ibc2,ibc3,ibc4             
      do i=1,ibc2
         read(80,*) BCE(i,ib3),iedge,BCF(i,ib3)
         BCC(i,ib3)= CBC(ib3)
      enddo

      close(80)
c... end of reading neu file

      open(unit=81,file=meshname_rea,status='unknown')
      write(81,94) nelem, 3, nelem 
   94 format(3i8,'   NEL,NDIM,NELV')

      do i=1,nelem
         i1=mapv(i,1)
         i2=mapv(i,2)
         i3=mapv(i,3)
         i4=mapv(i,4)
         write(81,95) i, i 
         write(81,101) point(i1,1),point(i2,1),point(i3,1),point(i4,1)
         write(81,101) point(i1,2),point(i2,2),point(i3,2),point(i4,2)
         write(81,101) point(i1,3),point(i2,3),point(i3,3),point(i4,3)
      enddo
   95 format('     ELEMENT   ',i10,' [ ',i10,' ] ')

      
      write(81,96) 
      write(81,97) 
      write(81,98) 
      write(81,99) 
      write(81,100) 

   96 format('***** CURVED SIDE DATA *****')
   97 format(' 0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')
   98 format('***** BOUNDARY CONDITIONS *****')
   99 format('***** NO FLUID   BOUNDARY CONDITIONS *****')
  100 format('***** THERMAL BOUNDARY CONDITIONS *****')
  101 format(4e18.6)

      nface = 4
      do ie = 1, nelem
         do ifc= 1, nface 
           BC = 'E  '
           do ibc= 1,nbsets
              if (BCF(ie,ibc).eq.ifc) then
                BC = BCC(ie,ibc)
              endif
            enddo
         write(81,102)  BC, ie, ifc, 0, 0, 0.0, 0.0, 0.0
      enddo
      enddo
  102 format(1x,a3,4i8,3e17.6)


      close(81)
      end

      subroutine blank(s,n)
      character*1 s(1)
      do i=1,n
        s(i)=' '
      enddo
      return
      end


      
