c-----------------------------------------------------------------------
      program trans
c
c     read and translate (replicate) 2d data set
c     Paul F. Fischer  (pff@cfm.brown.edu)  9/6/95.
c
      character*80 file
      character*1  file1(80)
      equivalence (file1,file)
      character*80 fout
      character*1  fout1(80)
      equivalence (fout1,fout)
c
      parameter(nelm=80000)
      common /array/ x(4,nelm),y(4,nelm),bc(5,4,nelm)
c
c     Get file name
      write(6,*) 'Input old (source) file name:'      
      call blank(file,80)
      read(5,80) file
      len = ltrunc(file,80)
   80 format(a80)
c
c     Get file name
      write(6,*) 'Input new (output) file name:'      
      call blank(fout,80)
      read(5,80) fout
      lou = ltrunc(fout,80)
c
      nrep = 1
      write(6,*) 'input number of reps: (0, 1, 2, 3,... etc.?):'
      read (5,*) nrep
c
c     Translate .rea data (for 2d only at present)
c
      call ccopy(file1(len+1),'.rea',4)
      call ccopy(fout1(lou+1),'.rea',4)
c
      open(unit=10, file=file)
      open(unit=11, file=fout)
      call transrea2(nrep,neln)
c
      write(6,*)
      write(6,6) neln,(fout1(k),k=1,len+4)
    6 format(i4,' elements written to ',40a1)
c
      close (unit=10)
      close (unit=11)
c
c     Translate .fld data, if available
c
      call ccopy(file1(len+1),'.fld',4)
      call ccopy(fout1(lou+1),'.fld',4)
      open(unit=10, file=file,status='old',err=999)
      open(unit=11, file=fout)
c
      do idmp=1,100
         call transfld(nrep,neln)
      enddo
c
      write(6,*)
      write(6,6) neln,(fout1(k),k=1,len+4)
c
      close (unit=10)
      close (unit=11)
  999 continue
      stop
      end
c-----------------------------------------------------------------------
      subroutine transfld(nrep,neln)
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
      character*4 s4
c
c     Nekton stuff
      parameter(nelm=800)
      parameter(nxm=9)
      common /array2/ data(6*nelm*nxm*nxm)
c
      integer nel,nx1,ny1
      save    nel,nx1,ny1
c
      if (nel.eq.0) then
         read(10,*,end=999,err=999) nel,nx1,ny1
         rewind(10)
      endif
c
      ndat = 0
      ndim = 2
      read(10,80,end=999,err=999) string
      if (indx1(string,' X ',3).ne.0) ndat = ndat + ndim
      if (indx1(string,' U ',3).ne.0) ndat = ndat + ndim
      if (indx1(string,' P ',3).ne.0) ndat = ndat + 1
      if (indx1(string,' T ',3).ne.0) ndat = ndat + 1
c
      write(s4,'(i4)') neln
      call ccopy(string,s4,4)
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
      read (10,'(6g11.4)') (data(ie),ie=1,nel)
      write(11,'(6g11.4)') (data(1),ie=1,neln)
c
      j1 = 1
      j2 = ndat
      do ie=1,nel
         do k=1,nx1*ny1
            read (10,*) (data(j),j=j1,j2)
            write(11,'(1p5e14.6)') (data(j),j=j1,j2)
            j1 = j2+1
            j2 = j2+ndat
         enddo
      enddo
c
      do irep = 1,nrep
      j1 = 1
      j2 = ndat
      do ie=1,nel
         do k=1,nx1*ny1
            write(11,'(1p5e14.6)') (data(j),j=j1,j2)
            j1 = j2+1
            j2 = j2+ndat
         enddo
      enddo
      enddo
c
   80 format(a80)
   81 format(80a1)
  999 continue
      return
      end
c
c-----------------------------------------------------------------------
      subroutine transrea(nrep,neln)
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
      character*4  string20
      equivalence (string20,string1(21))
c
c     Nekton stuff
      parameter(nelm=80000)
      common /array/ x(4,nelm),y(4,nelm),bc(5,4,nelm)
      character*3 cbc(4,nelm)
      character*1 ans
c
      logical ifflow,ifheat
c
c
c     Read parameters
c
      call readwrite2(string,'IFCOMP',6,'IFCHAR',6)
c
c     Read screen coordinates
c
      scale = nrep+1
c
      read (10,'(4e14.6)') (x(k,1),k=1,4)
      call cmult(x,scale,4)
      write(11,'(4e14.6)') (x(k,1),k=1,4)
      call readwrite(string,'MESH DATA',9)
c
      read (10,*) nel,ndim
      neln = nel+nrep*nel
      write(11,11) neln,ndim,neln
   11 format(3i6,11x,'NEL,NDIM,NELV')
c
c
c     Read & write xy data 
      do ie=1,nel      
         read (10,80) string
c        write(6 ,80) string
         read (10,*)   (x(k,ie),k=1,4)
         read (10,*)   (y(k,ie),k=1,4)
         len = ltrunc(string,80)
c        write(11,81) (string1(k),k=1,len)
c        write(11,90)  (x(k,ie),k=1,4)
c        write(11,90)  (y(k,ie),k=1,4)
      enddo
   90 format(4e14.6)
      len = ltrunc(string,80)
c
c     Translate Write xy data 
      nel4 = 4*nel
      xmin = glmin(x,nel4)
      xmax = glmax(x,nel4)
      write(6,*) 'rescale x/y? (scale = 1 ==> No.)',xmin,xmax
      read (5,*) rescale
      call cmult    (x,rescale,nel4)
      call cmult    (y,rescale,nel4)
c
      xmin = glmin(x,nel4)
      xmax = glmax(x,nel4)
      dx   = xmax-xmin
      write(6,*) 'xmin/max:',xmin,xmax
c
c     Hardwire translation distance to be 10.
c     dx   = 10.0
c
      do irep=0,nrep
         do ie=1,nel      
            ien = ie + nel*irep
            write(string20,'(i4)') ien
            write(11,81) (string1(k),k=1,len)
            write(11,90)  (x(k,ie),k=1,4)
            write(11,90)  (y(k,ie),k=1,4)
         enddo
         call cadd(x,dx,nel4)
      enddo
c
c     New. Curved sides added 4/25/97.  pff.
c
      read (10,80) string
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
      read (10, *) ncurve
      ncnew = (1+nrep)*ncurve
      write(11,96) ncnew
   96 format(i6,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')
      write(6,*) 'this is ncurve:',ncurve,ncnew,nel,neln,nrep
c
      do ic = 1,ncurve
         if (nel.lt.1000) then
            read(10,60) iedg,ieg,r1,r2,r3,r4,r5,ans
         else
            read(10,61) iedg,ieg,r1,r2,r3,r4,r5,ans
         endif
         r1 = rescale*r1
         r2 = rescale*r2
         r3 = rescale*r3
         r4 = rescale*r4
         r5 = rescale*r5
c
         do irep=0,nrep
           ien = ieg + irep*nel
           if (neln.lt.1000) then
             write(11,60) iedg,ien,r1,r2,r3,r4,r5,ans
c            write( 6,60) iedg,ien,r1,r2,r3,r4,r5,ans
           else
             write(11,61) iedg,ien,r1,r2,r3,r4,r5,ans
c            write( 6,61) iedg,ien,r1,r2,r3,r4,r5,ans
           endif
         enddo
c
      enddo
   60 format(i3,i3,5g14.6,1x,a1)
   61 format(i2,i6,5g14.6,1x,a1)
c
c     Boundary info....
c
      read (10,80) string
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
c
c     Fluid bc's
      read (10,80) string
c
      ifflow = .true.
      if (indx1(string,'NO',2).gt.0) ifflow=.false.
c
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
      if (ifflow) then
         do ie = 1,nel
         do  k = 1,4
            read (10,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
c           write( 6,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            if (cbc(k,ie).eq.'P  '.or.cbc(k,ie).eq.'E  ') 
     $          cbc(k,ie)='   '
            if (ie.lt.1000) then
               write(11,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            else
               write(11,21) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            endif
   20       FORMAT(1x,A3,2I3,5G14.6)
   21       FORMAT(1x,a3,i5,i1,5g14.6)
c  21       FORMAT(1x,a3,i5,i1,1p5e14.6)
         enddo
         enddo
c
         do irep = 1,nrep
         do ie = 1,nel
            id = ie+nel*irep
            do  k = 1,4
               if (id.lt.1000) then
                  write(11,20) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               else
                  write(11,21) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               endif
            enddo
         enddo
         enddo
      endif
c
c     Thermal bc's
      read (10,80) string
c
      ifheat = .true.
      if (indx1(string,'NO',2).gt.0) ifheat=.false.
c
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c     write( 6,81) (string1(k),k=1,len)
c     write( 6,*) 'ifheat:',ifheat
c
      if (ifheat) then
         do ie = 1,nel
         do  k = 1,4
            read (10,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            if (cbc(k,ie).eq.'P  '.or.cbc(k,ie).eq.'E  ') 
     $          cbc(k,ie)='   '
            if (ie.lt.1000) then
               write(11,20) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
            else
               write(11,21) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
            endif
         enddo
         enddo
c
         do irep = 1,nrep
         do ie = 1,nel
            id = ie+nel*irep
            do  k = 1,4
               if (id.lt.1000) then
                  write(11,20) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               else
                  write(11,21) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               endif
            enddo
         enddo
         enddo
      endif
c
      call readwrite(string,'endendend',9)
c
   80 format(a80)
   81 format(80a1)
      return
      end
c-----------------------------------------------------------------------
      subroutine blank(s,n)
      character*1 s(1)
      do i=1,n
        s(i)=' '
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      function ltrunc(s,n)
      character*1 s(1)
      ltrunc = 0
      do j=n,1,-1
         if (s(j).ne.' ') then
            ltrunc = j 
            return
         endif
      enddo
      return
      end
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*80 S1,S2
C
      N1=80-L2+1
      INDX1=0
      IF (N1.LT.1) RETURN
C
      DO 300 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            RETURN
         ENDIF
300   CONTINUE
C
      RETURN
      END
c
c-----------------------------------------------------------------------
      subroutine readwrite(sout,key,nk)
      character*80 sout,key
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      do i=1,90000      
         call blank(string,80)
         read (10,80,end=100,err=100) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)
         if (indx1(string,key,nk).ne.0) return
      enddo
  100 continue
c
   80 format(a80)
   81 format(80a1)
      return
      end
c
c-----------------------------------------------------------------------
      subroutine readwrite2(sout,key1,nk1,key2,nk2)
      character*80 sout,key1,key2
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      do i=1,90000      
         call blank(string,80)
         read (10,80,end=100,err=100) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)
         write(6 ,81) (string1(k),k=1,len)
         if (indx1(string,key1,nk1).ne.0) return
         if (indx1(string,key2,nk2).ne.0) return
      enddo
  100 continue
   80 format(a80)
   81 format(80a1)
      return
      end
c
c-----------------------------------------------------------------------
      function glmax(a,n)
      real a(1)
      tmax=-99.0E20
      do 100 i=1,n
         tmax=max(tmax,a(i))
  100 continue
      glmax=tmax
      return
      end
c
c-----------------------------------------------------------------------
      function glmin(a,n)
      real a(1)
      tmin=99.0E20
      do 100 i=1,n
         tmin=min(tmin,a(i))
  100 continue
      glmin = tmin
      return
      end
c-----------------------------------------------------------------------
      subroutine copy(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cadd(x,c,n)
      real x(1)
      do i=1,n
         x(i) = x(i)+c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cmult(x,c,n)
      real x(1)
      do i=1,n
         x(i) = x(i)*c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ccopy(x,y,n)
      character*1 x(1),y(1)
      do i=1,n
         x(i) = y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine transrea2(nrep,neln)
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
      character*4  string20
      equivalence (string20,string1(21))
      character*1  ans,ans2
c
c     Nekton stuff
      parameter(nelm=80000)
      common /array/ x(4,nelm),y(4,nelm),bc(5,4,nelm)
      character*3 cbc(4,nelm)
c
      logical ifflow,ifheat
c
c
c     Read parameters
c
      call readwrite2(string,'IFCOMP',6,'IFNM  ',6) !changed for CEM  
      !call readwrite2(string,'IFCOMP',6,'IFCHAR',6) !original          
c
c     Read screen coordinates
c
      scale = nrep+1
c
      read (10,'(4e14.6)') (x(k,1),k=1,4)
      call cmult(x,scale,4)
      write(11,'(4e14.6)') (x(k,1),k=1,4)
      call readwrite(string,'MESH DATA',9)
c
      read (10,*) nel,ndim
      neln = nel+nrep*nel
      write(11,11) neln,ndim,neln
   11 format(3i6,11x,'NEL,NDIM,NELV')
c
c
c     Read & write xy data 
      do ie=1,nel      
         read (10,80) string
c        write(6 ,80) string
         read (10,*)   (x(k,ie),k=1,4)
         read (10,*)   (y(k,ie),k=1,4)
         len = ltrunc(string,80)
      enddo
   90 format(4e14.6)
      len = ltrunc(string,80)
c
c     Translate Write xy data 
      nel4 = 4*nel
      xmin = glmin(x,nel4)
      xmax = glmax(x,nel4)
c
      write(6,*) 'rescale x/y? (scale = 1 ==> No.)',xmin,xmax
      read (5,*) rescale
c
      call cmult  (x,rescale,nel4)
      call cmult  (y,rescale,nel4)
      write(6,*) 'translate in x, y, or rotate? (x,y,r):'
      read (5,*) ans
      if (ans.eq.'x'.or.ans.eq.'X') then
         dy = 0
         dr = 0
         xmin = glmin(x,nel4)
         xmax = glmax(x,nel4)
         dx   = xmax-xmin
         write(6,*)'dx/xmin/max:',dx,xmin,xmax,' Keep default dx? (y/n)'
         read (5,*) ans
         if (ans.eq.'n'.or.ans.eq.'N') then
            write(6,*) 'Input new dx:'
            read (5,*) dx
         endif
      elseif (ans.eq.'y'.or.ans.eq.'Y') then
         dx = 0
         dr = 0
         ymin = glmin(y,nel4)
         ymax = glmax(y,nel4)
         dy   = ymax-ymin
         write(6,*)'dy/ymin/max:',dy,ymin,ymax,' Keep default dy? (y/n)'
         read (5,*) ans
         if (ans.eq.'n'.or.ans.eq.'N') then
            write(6,*) 'Input new dy:'
            read (5,*) dy
         endif
      else
         dx = 0
         dy = 0
         write(6,*)'input angle to rotate (degrees):'
         read (5,*) dr
         one = 1.
         pi  = 4.*atan(one)
         dr  = pi*dr/180
      endif
c
      dxk = dx
      dyk = dy
      drk = dr
      ctk = cos(drk)
      stk = sin(drk)
      do irep=1,nrep
         do ie=1,nel      
            ien = ie + nel*irep
            call copy(x(1,ien),x(1,ie),4)
            call copy(y(1,ien),y(1,ie),4)
            do j=1,4
               xx = ctk*x(j,ien) - stk*y(j,ien)
               yy = stk*x(j,ien) + ctk*y(j,ien)
               x(j,ien) = xx
               y(j,ien) = yy
            enddo
            call cadd(x(1,ien),dxk,4)
            call cadd(y(1,ien),dyk,4)
         enddo
         dxk = dxk + dx
         dyk = dyk + dy
         drk = drk + dr
         ctk = cos(drk)
         stk = sin(drk)
      enddo
      neln4 = neln*4
c
c     Manipulation loop:
c
   10 continue
      xmin = glmin(x,neln4)
      xmax = glmax(x,neln4)
      ymin = glmin(y,neln4)
      ymax = glmax(y,neln4)
      write(6,*) 'Translate/Rotate/Shift/Exit? (t,r,e,s)? '
     $           ,xmin,xmax,ymin,ymax
      read (5,*) ans
      if (ans.eq.'t'.or.ans.eq.'T') then
         write(6,*) 'Input x_shift, y_shift:'
         read (5,*) xtrans,ytrans
         call cadd(x,xtrans,neln4)
         call cadd(y,ytrans,neln4)
         goto 10
      elseif (ans.eq.'s'.or.ans.eq.'S') then
         write(6,*) 'Shift x or y ?'
         read (5,*) ans
         write(6,*) 'Input location separating shift section (x or y):'
         read (5,*) spoint
         write(6,*) 'Shift region? (<, >, or = to abort)'
         read (5,*) ans2
         write(6,*) 'Input amount to shift:'
         read (5,*) shift
         nshift = 0
         do i=1,neln4
            if (ans.eq.'x'.or.ans.eq.'X') then
               if (ans2.eq.'>'.and.x(i,1).gt.spoint) then
                  x(i,1) = x(i,1)+shift
                  nshift = nshift + 1
               elseif (ans2.eq.'<'.and.x(i,1).lt.spoint) then
                  x(i,1) = x(i,1)+shift
                  nshift = nshift + 1
               endif
            elseif (ans.eq.'y'.or.ans.eq.'Y') then
               if (ans2.eq.'>'.and.y(i,1).gt.spoint) then
                  y(i,1) = y(i,1)+shift
                  nshift = nshift + 1
               elseif (ans2.eq.'<'.and.y(i,1).lt.spoint) then
                  y(i,1) = y(i,1)+shift
                  nshift = nshift + 1
               endif
            endif
         enddo
         goto 10
      elseif (ans.eq.'r'.or.ans.eq.'R') then
         write(6,*) 'Input rotation angle (degrees):'
         read (5,*) angle
         one = 1.
         pi  = 4.*atan(one)
         angle = pi*angle/180.
         ca = cos(angle)
         sa = sin(angle)
         do i=1,neln4
            xx = x(i,1)
            yy = y(i,1)
            x(i,1) = ca*xx - sa*yy
            y(i,1) = sa*xx + ca*yy
         enddo
         goto 10
      endif
c
c     Output
c
      do ie=1,neln
         write(string20,'(i4)') ie
         write(11,81) (string1(k),k=1,len)
         write(11,90)  (x(k,ie),k=1,4)
         write(11,90)  (y(k,ie),k=1,4)
      enddo
c
c     New. Curved sides added 4/25/97.  pff.
c
      read (10,80) string
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
      read (10, *) ncurve
      ncnew = (1+nrep)*ncurve
      write(11,96) ncnew
   96 format(i6,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')
      write(6,*) 'this is ncurve:',ncurve,ncnew,nel,neln,nrep
c
      do ic = 1,ncurve
         if (nel.lt.1000) then
            read(10,60) iedg,ieg,r1,r2,r3,r4,r5,ans
         else
            read(10,61) iedg,ieg,r1,r2,r3,r4,r5,ans
         endif
         r1 = rescale*r1
         r2 = rescale*r2
         r3 = rescale*r3
         r4 = rescale*r4
         r5 = rescale*r5
c
         do irep=0,nrep
           ien = ieg + irep*nel
           if (neln.lt.1000) then
             write(11,60) iedg,ien,r1,r2,r3,r4,r5,ans
c            write( 6,60) iedg,ien,r1,r2,r3,r4,r5,ans
           else
             write(11,61) iedg,ien,r1,r2,r3,r4,r5,ans
c            write( 6,61) iedg,ien,r1,r2,r3,r4,r5,ans
           endif
         enddo
c
      enddo
   60 format(i3,i3,5g14.6,1x,a1)
   61 format(i2,i6,5g14.6,1x,a1)
c
c     Boundary info....
c
      read (10,80) string
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
c
c     Fluid bc's
      read (10,80) string
c
      ifflow = .true.
      if (indx1(string,'NO',2).gt.0) ifflow=.false.
c
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
      if (ifflow) then
         do ie = 1,nel
         do  k = 1,4
            read (10,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
c           write( 6,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            if (cbc(k,ie).eq.'P  '.or.cbc(k,ie).eq.'E  ') 
     $          cbc(k,ie)='   '
            if (ie.lt.1000) then
               write(11,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            else
               write(11,21) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            endif
   20       FORMAT(1x,A3,2I3,5G14.6)
   21       FORMAT(1x,a3,i5,i1,5g14.6)
c  21       FORMAT(1x,a3,i5,i1,1p5e14.6)
         enddo
         enddo
c
         do irep = 1,nrep
         do ie = 1,nel
            id = ie+nel*irep
            do  k = 1,4
               if (id.lt.1000) then
                  write(11,20) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               else
                  write(11,21) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               endif
            enddo
         enddo
         enddo
      endif
c
c     Thermal bc's
      read (10,80) string
c
      ifheat = .true.
      if (indx1(string,'NO',2).gt.0) ifheat=.false.
c
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c     write( 6,81) (string1(k),k=1,len)
c     write( 6,*) 'ifheat:',ifheat
c
      if (ifheat) then
         do ie = 1,nel
         do  k = 1,4
            read (10,20) cbc(k,ie),id,jd,(bc(j,k,ie),j=1,5)
            if (cbc(k,ie).eq.'P  '.or.cbc(k,ie).eq.'E  ') 
     $          cbc(k,ie)='   '
            if (ie.lt.1000) then
               write(11,20) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
            else
               write(11,21) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
            endif
         enddo
         enddo
c
         do irep = 1,nrep
         do ie = 1,nel
            id = ie+nel*irep
            do  k = 1,4
               if (id.lt.1000) then
                  write(11,20) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               else
                  write(11,21) cbc(k,ie),id,k,(bc(j,k,ie),j=1,5)
               endif
            enddo
         enddo
         enddo
      endif
c
      call readwrite(string,'endendend',9)
c
   80 format(a80)
   81 format(80a1)
      return
      end
c-----------------------------------------------------------------------
