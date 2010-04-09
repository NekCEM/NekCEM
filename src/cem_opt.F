c-----------------------------------------------------------------------
      subroutine set_options 

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'

      integer k, nxyzf
      nxyzf = nx1*nz1*2*ndim*nelt

      IFCEM    =.TRUE.
      IFVARDT  =.false.
      IFFILTER =.false.
      IFDIELEC =.false.
      IFIC     =.false.
      IFPOISSON=.false.
      IFDTCTL  =.false.
      IFSCAT   =.false.
      IFSFTF   =.false.
      IFDRUDE  =.false.
      IFLOREN  =.false.
      IFBEAM   =.false.
      IFWAKE   =.false.
      IFNAPY   =.false.
      IFWAVEGUIDE=.false.
      IFPBG    =.false.
      IFWAKET  =.false.
      IFSOL    =.false.
      IFVTK    =.false.

      if (param( 17).eq.1)  IFVARDT  =.true.
      if (param( 85).ne.0)  IFFILTER =.true.
      if (param( 93).ne.0)  IFDIELEC =.true.
      if (param( 94).ne.0)  IFIC     =.true.
      if (param( 95).ne.0)  IFPOISSON=.true.
      if (param( 99).ne.0)  IFDTCTL  =.true.
      if (param(100).eq.1)  IFSCAT   =.true.
      if (param(100).eq.2)  IFSFTF   =.true.
      if (param(101).eq.1)  IFDRUDE  =.true.
      if (param(101).eq.12) IFLOREN  =.true.
      if (param(101).eq.2)  IFBEAM   =.true.

      if ((param(66).lt.0).or.(param(67).lt.0).or.(param(68).lt.0)) then
          incflag = -1
      else
          incflag=1
      endif

      if (param(101).eq.3)  then
             IFBEAM   =.true.
             IFWAKE   =.true.
      endif
      if (param(101).eq.4)  then
             IFBEAM   =.true.
             IFWAKE   =.true.
             IFNAPY   =.true.
      endif
      if (param(101).eq.5)  then
             IFBEAM   =.true.
             IFWAKE   =.true.
             IFNAPY   =.true.
             IFWAKET  =.true.
      endif
      if (param(102).ne.0)  IFSOL    =.true.
      if (param(103).ne.0)  IFVTK    =.true.

      if (nid.eq.0) then

      write(6,*) 'IFDRUDE     ',IFDRUDE
      write(6,*) 'IFLOREN     ',IFLOREN
      write(6,*) 'IFBEAM      ',IFBEAM 
      write(6,*) 'IFWAKE      ',IFWAKE 
      write(6,*) 'IFMOVW      ',IFMOVW 
      write(6,*) 'IFFILTER    ',IFFILTER
      write(6,*) 'IFNAPY      ',IFNAPY 
      write(6,*) 'IFWAKET     ',IFWAKET
      write(6,*) 'IFVARDT     ',IFVARDT
      write(6,*) 'IFPBG       ',IFPBG                     
      write(6,*) 'IFDIELEC    ',IFDIELEC 
      write(6,*) 'IFPOISSON   ',IFPOISSON 
      write(6,*) 'IFIC        ',IFIC            
      write(6,*) 'IFSOL       ',IFSOL                      
      write(6,*) 'IFSRC       ',IFSRC
      write(6,*) 'IFSCAT      ',IFSCAT
      write(6,*) 'IFSFTF      ',IFSFTF
      write(6,*) 'IFCENTRAL   ',IFCENTRAL
      write(6,*) 'IFUPWIND    ',IFUPWIND
      write(6,*) 'IFDEALIAS   ',IFDEALIAS 
      write(6,*) 'IFTM        ',IFTM
      write(6,*) 'IFTE        ',IFTE
      write(6,*) 'IFEIG       ',IFEIG
      write(6,*) 'IFRK4       ',IFRK4 
      write(6,*) 'IFEXP       ',IFEXP
      write(6,*) 'IFDTCTL     ',IFDTCTL
      write(6,*) 'IFNM        ',IFNM
      write(6,*) 'IFVTK       ',IFVTK

      endif

      if (IFEXP) write(6,*) 'Krylov_Dim =',larnol

      if (IFNM) then

          cspeed  = param(1)
          permit_0= param(2)
          permea_0= param(3)
          permit_r= param(4)
          permea_r= param(5) ! not used
          gamma_p = param(6)
          freq_p  = param(7)

          xyzscale= 1.0      
          if (IFDRUDE) xyzscale=1.e+9 ! for visual tool
          if (IFBEAM ) xyzscale=1.0  
          if (IFWAKE ) xyzscale=1.0 

          nmscale = 1.0   
          if (IFDRUDE) nmscale =1.e-9 ! for visual tool
          if (IFBEAM ) nmscale =1.0 
          if (IFWAKE ) nmscale =1.0 
           
      else
          cspeed  = 1.0          
          permit_0= 1.0         
          permea_0= 1.0             
          permit_r= param(4)
          permea_r= param(5) ! not used
          gamma_p = param(8)
          freq_p  = param(9)

          xyzscale= 1.0 
          nmscale = 1.0
      endif

      if (nid.eq.0) then

      write(6,*)                  
      write(6,*) 'cspeed   =',cspeed
      write(6,*) 'permit_0 =',permit_0
      write(6,*) 'permea_0 =',permea_0
      write(6,*) 'permit_r =',permit_r
      write(6,*) 'permea_r =',permea_r
      write(6,*)                  
      write(6,*) 'nmscale  =',nmscale   
      write(6,*) 'fscale   =',param(59) 
      write(6,*)                  

      endif

      if     ((IFCENTRAL).and.(IFUPWIND))  then
             write(6,*) 'ERROR:: not allowed IFCENTRAL=T & IFUPWIND= T'
             call exitt
      elseif ((.not.(IFCENTRAL)).and.(.not.(IFUPWIND))) then
             write(6,*) 'ERROR:: not allowed IFCENTRAL=F & IFUPWIND= F'
             call exitt
      endif

      if     (.not.IF3D) then
         if     ((IFTM).and.(IFTE))  then
             write(6,*) 'ERROR:: not allowed IFTM= T & IFTE= T'
             call exitt
         elseif  ((.not.(IFTM)).and.(.not.(IFTE))) then
             write(6,*) 'ERROR:: not allowed IFTM= F & IFTE= F'
             call exitt
         endif
      endif 

      call ifoption (IFPML) ! globalize logical value to all processors
      call ifoption (IFPEC)
      call ifoption (IFPP )
      call assign_bdry_output

      if (if3d) then
        hxactive = .true.
        hyactive = .true.
        hzactive = .true.
        exactive = .true.
        eyactive = .true.
        ezactive = .true.
        
        fldidx(1) = 1
        fldidx(2) = 2
        fldidx(3) = 3
        fldidx(4) = 4
        fldidx(5) = 5
        fldidx(6) = 6
        fldcnt = 6
      elseif (iftm) then
        hxactive = .true.
        hyactive = .true.
        hzactive = .false.
        exactive = .false.
        eyactive = .false.
        ezactive = .true.

        fldidx(1) = 1 !hx
        fldidx(2) = 2
        fldidx(3) = 0 !hz
        fldidx(4) = 0 !ex
        fldidx(5) = 0
        fldidx(6) = 3 !ez
        fldcnt = 3
      elseif (ifte) then
        hxactive = .false.
        hyactive = .false.
        hzactive = .true.
        exactive = .true.
        eyactive = .true.
        ezactive = .false.

        fldidx(1) = 0 !hx
        fldidx(2) = 0
        fldidx(3) = 3 !hz
        fldidx(4) = 1 !ex
        fldidx(5) = 2
        fldidx(6) = 0 !ez
        fldcnt = 3
      endif

      if     (IFEIG) then
         if (nid.eq.0) write(6,*) 
         if (nid.eq.0) 
     $       write(6,*) 'IFEIG T :: start eigenvalue computation'
         if (nid.eq.0) write(6,*) 
      else
         if     ((IFRK4).and.(IFEXP)) then
            if (nid.eq.0)  
     $       write(6,*) 'ERROR :: not allowed  IFRK4= T & IFEXP= T'
             call exitt
         elseif ((.not.(IFRK4)).and.(.not.(IFEXP))) then
            if (nid.eq.0) 
     $       write(6,*) 'ERROR :: not allowed  IFRK4 F &  IFEXP F'
             call exitt
         endif
      endif

      return
      end  
