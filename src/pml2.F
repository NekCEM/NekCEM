      
       subroutine pml_faces(faceary)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'

      ! Given the values on all the faces (in faceary),
      ! this routine will increment the two opposing faces
      ! if at least one of them has a non-zero value.
      !
      ! If you think about it, this will look a little bit
      ! like marching through opposing faces.
      ! 
      ! (This procedure applies to all values on a face,
      ! so you need to be careful with edges (which are shared
      ! between faces). Edge values may "march off" parallel to
      ! the face.)

      real faceary(lx1*lz1,2*ldim,lelt)

      integer nfaces
      integer axis, elt, idx
      integer posface, negface

      nfaces= 2*ndim   
      npmlsrf= nelt*nfaces
      call izero(pmlsrf,npmlsrf)

      do  ie= 1,nelv
        do axis=0,ldim-1
           posface= eface(2+2*axis)  ! Ed's numbering
           negface= eface(1+2*axis)
          do idx = 1,lx1*lz1  ! check out for 3d
            itmp1= int(faceary(idx,posface ,ie))
            itmp2= int(faceary(idx,negface ,ie))

            if (   ((itmp1.eq.1).and.(itmp2.eq.0))
     $        .or. ((itmp1.eq.0).and.(itmp2.eq.1)) ) then

               if (pmltag(ie).eq.0) then
                   pmlsrf(ie,posface)=itmp1 ! Ed's numbering
                   pmlsrf(ie,negface)=itmp2
c                  write(6,14) ie,idx,posface,negface,      
c    $             faceary(idx,posface,ie),faceary(idx,negface,ie)
               endif

            endif

          enddo
        enddo
      enddo

      return
      end

      subroutine march_faces(faceary)

      ! Given the values on all the faces (in faceary),
      ! this routine will increment the two opposing faces
      ! if at least one of them has a non-zero value.
      !
      ! If you think about it, this will look a little bit
      ! like marching through opposing faces.
      ! 
      ! (This procedure applies to all values on a face,
      ! so you need to be careful with edges (which are shared
      ! between faces). Edge values may "march off" parallel to
      ! the face.)

      include 'SIZE'
      include 'TOTAL'

      real faceary(lx1*lz1,2*ldim,lelt)

      integer nfaces
      integer axis, elt, idx
      integer posface, negface
      integer posface1, negface1

      real sum_here

      nxyzf =lx1*lz1*2*ldim*lelt

      ! implement marching inside each element
      do ie= 1,nelv
        do axis=0,ldim-1
           posface= eface(2+2*axis)  ! Ed's numbering
           negface= eface(1+2*axis)
        do idx = 1,lx1*lz1
           sum_here= faceary(idx,posface,ie)+faceary(idx,negface,ie)
           if (sum_here.ne.0) then
               faceary(idx,posface,ie)=faceary(idx,posface,ie)+1
               faceary(idx,negface,ie)=faceary(idx,negface,ie)+1
           endif
c          write(6,13) ie,idx,posface,negface,
c    $     faceary(idx,posface,ie),faceary(idx,negface,ie) 
        enddo
        enddo
      enddo

 13   format(4i5,1p2e11.2) 

      !maximum across neighboring elements

      call gs_op_fields(gsh_face,faceary,nxyzf,1,1,1,0)

      ! changes from P to PML in rea contains the connectivity
      ! info from the original built in prex. so that it still
      ! recognizes neighboring element with periodic connectivity.
      ! so that maximum with worng neighboring faces gives
      ! incorrect info. 
      ! BTW at this moment it's not well-understood why fgs_gop is 
      ! it doesn't seem necessary to have this. FIXME misun 10/23/08

      return
      end
c-----------------------------------------------------------------------
      subroutine intelt2full(fullary, eltary)

      ! Convert an integer per-element array to a
      ! visualizable per-point float array.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'

      integer eltary(lelt)
      real fullary(lx1,ly1,lz1,lelt)

      integer elt,i

      call rzero(fullary, lx1*ly1*lz1*lelt)

      do elt=1,nelt
        do i=1,nx1*ny1*nz1
          fullary(i,1,1,elt) = eltary(elt)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pml_fill_faceary(faceary, thick)

      ! Fill a face array with numbers such that the faces
      ! at the innermost layer of the PML (and only those) have 
      ! a value of 1. Values increase (by 1 per layer) towards
      ! the outside.

      include 'SIZE'
      include 'TOTAL'
      include 'PML'

      real    faceary(lx1*lz1,2*ldim,lelt)
      integer thick
      integer elt, face, i, axis, ix, iz
      character*3 cb           

      call rzero(faceary, lx1*lz1*2*ldim*lelt)

      ! produce a face array that has ones at the *interior*
      ! of all PML faces.
      ! (We need to avoid the edges since those will march off 
      ! parallel to the face, which we do not want.)

      do ie=1,nelv
        do ifc=1,2*ndim
           cb = cbc(ifc,ie,2)

          if  (cb.eq.'PML' .or. cb.eq.'pml') then
            if (if3d) then
              do ix=1,nx1-2
              do iz=1,nz1-2
                faceary(ix+iz*lx1+1, ifc, ie) = 1
                ! interior points of the face are 1
              enddo
              enddo
            else
              do ix=2,nx1-1
                faceary(ix, ifc, ie) = 1
                ! interior points of the face are 1
              enddo
            endif
          endif
        enddo
      enddo

      ! feed the face array to march_faces to
      ! produce the desired thickness of the pml.

      do i=1,thick
        call march_faces(faceary)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dir_local_to_global(globdir, elt, dir)

      ! Translate the direction `dir' from local to global,
      ! w.r.t. element number `el'.

      ! (A direction is essentially is essentially a face
      ! index in pff/sym convention, indicating -x,+x,-y,+y...
      ! by 1,2,3,4,...)

      ! This assumes that the Jacobian of element number `el'
      ! is a permutation matrix, or at least very nearly so.
      ! (The PML will fail if that's not the case anyway...)

      ! NOTE: It'd be wicked nice if we didn't need this routine.

      include 'SIZE'
      include 'TOTAL'

      integer globdir, elt, dir, axis
      real locvec(3), globvec(3)

      integer i, argmax
      real biggest

      call rzero(locvec, 3)

      ! create a local unit vector that points in `dir'.
      axis = (dir-1)/2+1
      if (mod(dir-1,2).eq.0) then
        locvec(axis) = -1
      else
        locvec(axis) = 1
      endif

      ! I'm using only one point to decide this. Am I crazy?

      ! A learning experience is one of those things that says, 
      ! 'You know that thing you just did? Don't do that.'
      !                 -- Douglas Adams

      globvec(1) =   rxm1(1,1,1,elt)*locvec(1)
     $             + sxm1(1,1,1,elt)*locvec(2)
     $             + txm1(1,1,1,elt)*locvec(3)
      globvec(2) =   rym1(1,1,1,elt)*locvec(1)
     $             + sym1(1,1,1,elt)*locvec(2)
     $             + tym1(1,1,1,elt)*locvec(3)
      globvec(3) =   rzm1(1,1,1,elt)*locvec(1)
     $             + szm1(1,1,1,elt)*locvec(2)
     $             + tzm1(1,1,1,elt)*locvec(3)

      ! determine argmax_i |globvec(i)|
      biggest = 0
      argmax = 0
      do i=1,3
        if (abs(globvec(i)).ge.biggest) then
          biggest = abs(globvec(i))
          argmax = i
        endif
      enddo

      ! from the argmax, determine the global direction
      globdir = (argmax-1)*2+1 ! default to neg. direction
      if (globvec(argmax).ge.0) then ! if pos. direction
        globdir = globdir + 1 ! say so.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pml_extent_and_tags(inner, outer, eltag, faceary)

      include 'SIZE'
      include 'TOTAL'

      ! Fill inner, outer and eltag as described below, based
      ! on faceary.

      ! inner and outer boundary coordinates of the pml, i.e.
      ! 
      ! |PMLPMLPML|........computational domain........|PMLPMLPML|
      ! 1         2                                    3         4
      !
      ! +-----> x
      ! 
      ! outer(1) = 1
      ! inner(1) = 2
      ! inner(2) = 3
      ! outer(2) = 4

      real inner(2*ldim),outer(2*ldim)

      ! see pmltag (src/PML, common block PML) for a description
      ! of what eltag is

      integer eltag(lelt)

      ! faceary is as described in pml_fill_faceary
      real faceary(lx1*lz1,2*ldim,lelt)

      ! END PARAMETERS

      ! The "indicative point index"
      ! Since we are using interior points of the faces (see march_faces),
      ! this gives one (arbitrary) point index in the 
      ! interior of a face that we can use to check if that face
      ! is part of a PML.

      integer  ind_ptidx 
      real     far_here, far_opp
      real     pmlinf
      parameter (pmlinf = 1e20)

      ! table of opposing face indices, symm/pff convention
      integer oppface(6) 
      data    oppface /2,1,4,3,6,5/

      integer elt, face, axis, globface, globaxis
      integer globsign
      real    mincoord, maxcoord
      real    gop_work
      real    minx, maxx, miny, maxy, minz, maxz

      ! BEGIN CODE

      ! determine the indicative point index (see above)

      if (if3d) then
        ind_ptidx = 1+1*lx1+1
      else
        ind_ptidx = 2
      endif

      ! initialize outer, inner to "infinities"

      do axis=0,ndim-1
         outer((axis*2)+1) =  pmlinf
         inner((axis*2)+1) = -pmlinf
         inner((axis*2)+2) =  pmlinf
         outer((axis*2)+2) = -pmlinf
      enddo

      ! fill outer, inner, eltag for local processor

      do elt=1,nelv
        eltag(elt) = 0
      enddo

      do elt=1,nelv
        do axis=1,ndim
          face = (axis-1)*2 + 1 ! pff/sym convention

          ! face vs. oppface(face): pff ordering
          ! eface: ed's ordering
          far_here= faceary(ind_ptidx, eface(face), elt)
          far_opp = faceary(ind_ptidx, eface(oppface(face)), elt)

          if ((far_here.ne.0).and.(far_opp.ne.0)) then
            ! element number `elt' is inside the PML

            if (far_here.eq.far_opp) then
              if (nid.eq.0) then
              write (*,*) 'No "gradient" in PML indicators.'
              write (*,*) 'This should not happen, unless your PMLs'
              write (*,*) 'are colliding in the middle.'
              write (*,*) 'Aborting.'
              endif
              call exitt
            endif

            call dir_local_to_global(globface, elt, face)
            globaxis = (globface-1)/2+1
            globsign = mod(globface-1,2)*2 - 1

            ! force globface to be the negative-side face
            globface = (globaxis-1)*2 + 1

            ! obtain mincoord, maxcoord of element
            if (globaxis.eq.1) then
              mincoord = vlmin(xm1(1,1,1,elt),lx1*ly1*lz1)
              maxcoord = vlmax(xm1(1,1,1,elt),lx1*ly1*lz1)
            elseif (globaxis.eq.2) then
              mincoord = vlmin(ym1(1,1,1,elt),lx1*ly1*lz1)
              maxcoord = vlmax(ym1(1,1,1,elt),lx1*ly1*lz1)
            elseif (globaxis.eq.3) then
              mincoord = vlmin(zm1(1,1,1,elt),lx1*ly1*lz1)
              maxcoord = vlmax(zm1(1,1,1,elt),lx1*ly1*lz1)
            endif

            ! update inner, outer

            if (globsign*(far_here-far_opp).gt.0) then
              ! Everything taken into account, we're on the
              ! positive side of the PML.
              ! outer grows towards +
              ! inner grows towards -

              inner(oppface(globface)) = 
     $            min(inner(oppface(globface)), mincoord)
              outer(oppface(globface)) = 
     $            max(outer(oppface(globface)), maxcoord)

              eltag(elt) = ior(eltag(elt), 
     $          ishft(1, oppface(globface)-1))
            endif

            if (globsign*(far_here-far_opp).lt.0) then
              ! Everything taken into account, we're on the
              ! negative side of the PML.
              ! outer grows towards -
              ! inner grows towards +

              inner(globface) = max(inner(globface), maxcoord)
              outer(globface) = min(outer(globface), mincoord)

              eltag(elt) = ior(eltag(elt), ishft(1, globface-1))
            endif
          endif
        enddo
      enddo

      ! fill outer, inner globally

      call gop(inner(1), gop_work, 'M  ', 1)
      call gop(inner(2), gop_work, 'm  ', 1)
      call gop(inner(3), gop_work, 'M  ', 1)
      call gop(inner(4), gop_work, 'm  ', 1)

      call gop(outer(1), gop_work, 'm  ', 1)
      call gop(outer(2), gop_work, 'M  ', 1)
      call gop(outer(3), gop_work, 'm  ', 1)
      call gop(outer(4), gop_work, 'M  ', 1)

      if (if3d) then
        call gop(inner(5), gop_work, 'M  ', 1)
        call gop(inner(6), gop_work, 'm  ', 1)

        call gop(outer(5), gop_work, 'm  ', 1)
        call gop(outer(6), gop_work, 'M  ', 1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine pml_calc_sigma(inner, outer, eltag, order, referr)

      ! Calculate the profile of the PML parameter sigma, and store it
      ! in the `PML' common block variable pmlsigma.

      ! `order' gives the order of the polynomial grading of the
      ! profile, `referr' is the desired reflection error for waves
      ! orthogonally incident on the PML.

      ! See Taflove book, (7.55), (7.57).

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'

      real     inner(2*ldim),outer(2*ldim)
      integer  eltag(lelt)
      real     order, referr

      integer  elt,e, i, face, axis
      real     point(3), zero2one
      real     eta, sigmamax
      real     width

      integer icalld
      save    icalld
      data    icalld/0/

      if (icalld.eq.0) then
          l=0
          do ie=1,nelt
             if (pmltag(ie).ne.0) then
                 l=l+1
                 pmlptr(l)=ie       !pml pointer                             
                !write(6,*) 'pml',nid,ie,pmltag(ie)! check out meaning
             endif
          enddo
          maxpml = l                !local total pml elts           
          maxpmlg= iglsum(maxpml,1) !global total pml elts

          if (nid.eq.0) write(6,*) 'local  total: PML elts=',maxpml
          if (nid.eq.0) write(6,*) 'global total: PML elts=',maxpmlg
          icalld = 1
      endif

      n = nx1*ny1*nz1*nelt*ndim
      call rzero(pmlsigma,n)

      do e=1,maxpml
         elt = pmlptr(e) 

          do face=1,2*ndim ! pff/sym convention
             axis = (face-1)/2 + 1
             width = abs(outer(face)-inner(face))

             if (iand(eltag(elt),ishft(1,face-1)).ne.0) then 
             ! there is a PML in the direction of `face' here

              do i=1,nx1*ny1*nz1
                j=i+nx1*ny1*nz1*(elt-1)
                eta  = sqrt(permeability(j)/permittivity(j))
                sigmamax = -(order+1)*log(referr)/(2*eta*width)

                point(1) = xm1(i,1,1,elt)
                point(2) = ym1(i,1,1,elt)
                point(3) = zm1(i,1,1,elt)

                zero2one = (point(axis)-inner(face))/
     $                     (outer(face)-inner(face))

                constant = sigmamax * zero2one**order
                pmlsigma(j,axis) = constant 

               enddo
              endif
          enddo
      enddo
  11  format('zero--',2i3,2e17.7)
   
      if (if3d) then
      else
        do e=1,nelt       
          do i=1,nx1*ny1*nz1
             j=i+nx1*ny1*nz1*(e-1)
             pmlsigma(j,3) = 0                     
          enddo
        enddo
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine pml_setup

      include 'SIZE'
      include 'TOTAL'
      include 'PML'
      include 'SCRATCH'

      integer  pmlthick
      real     pmlorder
      real     pmlreferr 
      ! how much of the incident signal is permitted to be reflected
      ! by the PML. See the term R(0) in formula (7.57) of the Taflove
      ! book.

      pmlthick = int(param(77))
      pmlorder = param(78)
      pmlreferr= param(79)

      if ((pmlthick.lt.1).or.(pmlthick.gt.10)) then
        write (*,*) 'Something is wrong with your pmlthick setting'
        write (*,*) 'in your .rea file. It needs to be between one'
        write (*,*) 'and ten.'
        call exitt
      endif
      if ((pmlorder.lt.1).or.(pmlorder.gt.10)) then
        write (*,*) 'Something is wrong with your pmlorder setting'
        write (*,*) 'in your .rea file. It needs to be between one'
        write (*,*) 'and ten.'
        call exitt
      endif
      if ((pmlreferr.lt.0).or.(pmlreferr.gt.1)) then
        write (*,*) 'Something is wrong with your pmlreferr setting'
        write (*,*) 'in your .rea file. It needs to be between zero'
        write (*,*) 'and one.'
        call exitt
      endif

      call pml_fill_faceary   (faceary, pmlthick)
      call pml_extent_and_tags(pmlinner, pmlouter, pmltag, faceary)
      call pml_calc_sigma     (pmlinner, pmlouter, pmltag, 
     $                         pmlorder, pmlreferr)
      call pml_faces          (faceary)
      return
      end
c-----------------------------------------------------------------------
      subroutine find_peak(field)

      ! This routine helps diagnose localized blowups by printing
      ! the element and index where the current maximum of the full
      ! field `field' is located.

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'

      real field(lx1,ly1,lz1,lelt)
      integer elt, i
      real absmax, abshere
      integer maxelt, maxi

      absmax = 0
      maxelt = 0
      maxi = 0

      do elt=1,nelt
        do i=1,lx1*ly1*lz1
          abshere = abs(field(i,1,1,elt))
          if (abshere.gt.absmax) then
            absmax = abshere
            maxelt = elt
            maxi = i
          endif
        enddo
      enddo

      write (*,*) 'max',absmax,' reached at elt', maxelt,' index', maxi

      return
      end

c-----------------------------------------------------------------------
      subroutine pml_step

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'

      real     sigx, sigy, sigz
      real     inv_mass
      integer  e,ie


      n   = nx1*ny1*nz1 
      nxyz= nx1*ny1*nz1 

      if (if3d) then
       do ie=1,maxpml
          e=pmlptr(ie) 
         do i=1,n                   

         j=i+nxyz*(e-1)

         bminv= 1./bm1(i,1,1,e)
         sigx = pmlsigma(j,1)
         sigy = pmlsigma(j,2)
         sigz = pmlsigma(j,3)

         respmlbn(j,1)=reshn(j,1)*bminv-sigy/permittivity(j)*pmlbn(j,1)
         respmlbn(j,2)=reshn(j,2)*bminv-sigz/permittivity(j)*pmlbn(j,2)
         respmlbn(j,3)=reshn(j,3)*bminv-sigx/permittivity(j)*pmlbn(j,3)
         respmldn(j,1)=resen(j,1)*bminv-sigy/permittivity(j)*pmldn(j,1)
         respmldn(j,2)=resen(j,2)*bminv-sigz/permittivity(j)*pmldn(j,2)
         respmldn(j,3)=resen(j,3)*bminv-sigx/permittivity(j)*pmldn(j,3)

         respmlhn(j,1)=-sigy/permittivity(j)*pmlbn(j,1) 
     $                 +sigx/permittivity(j)*pmlbn(j,1)
     $                 -sigz*permeability(j)/permittivity(j)*hn(j,1)
         respmlhn(j,2)=-sigz/permittivity(j)*pmlbn(j,2) 
     $                 +sigy/permittivity(j)*pmlbn(j,2)
     $                 -sigx*permeability(j)/permittivity(j)*hn(j,2)
         respmlhn(j,3)=-sigx/permittivity(j)*pmlbn(j,3) 
     $                 +sigz*permeability(j)*pmlbn(j,3)
     $                 -sigy*permeability(j)/permittivity(j)*hn(j,3)
         respmlen(j,1)=-sigy/permittivity(j)*pmldn(j,1) 
     $                 +sigx/permittivity(j)*pmldn(j,1)
     $                 -sigz*en(j,1)
         respmlen(j,2)=-sigz/permittivity(j)*pmldn(j,2) 
     $                 +sigy/permittivity(j)*pmldn(j,2)
     $                 -sigx*en(j,2)
         respmlen(j,3)=-sigx/permittivity(j)*pmldn(j,3) 
     $                 +sigz/permittivity(j)*pmldn(j,3)
     $                 -sigy*en(j,3)

         reshn(j,1)=reshn(j,1)+respmlhn(j,1)*bm1(i,1,1,e) 
         reshn(j,2)=reshn(j,2)+respmlhn(j,2)*bm1(i,1,1,e) 
         reshn(j,3)=reshn(j,3)+respmlhn(j,3)*bm1(i,1,1,e) 
         resen(j,1)=resen(j,1)+respmlen(j,1)*bm1(i,1,1,e) 
         resen(j,2)=resen(j,2)+respmlen(j,2)*bm1(i,1,1,e) 
         resen(j,3)=resen(j,3)+respmlen(j,3)*bm1(i,1,1,e) 
 
         enddo
       enddo
        
      else

       if     (IFTE) then

        do ie=1,maxpml
         e =pmlptr(ie) 
         do i = 1,n                

         j=i+nxyz*(e-1)

         bminv= 1./bm1(i,1,1,e)
         sigx = pmlsigma(j,1)
         sigy = pmlsigma(j,2)
         sigz = 0

         respmlbn(j,3)=reshn(j,3)*bminv-sigx/permittivity(j)*pmlbn(j,3)
         respmldn(j,1)=resen(j,1)*bminv-sigy/permittivity(j)*pmldn(j,1)
         respmldn(j,2)=resen(j,2)*bminv-sigz/permittivity(j)*pmldn(j,2)

         respmlhn(j,3)=-sigx/permittivity(j)*pmlbn(j,3) 
     $                 +sigz*permeability(j)*pmlbn(j,3)
     $                 -sigy*permeability(j)/permittivity(j)*hn(j,3)
         respmlen(j,1)=-sigy/permittivity(j)*pmldn(j,1) 
     $                 +sigx/permittivity(j)*pmldn(j,1)-sigz*en(j,1)
         respmlen(j,2)=-sigz/permittivity(j)*pmldn(j,2) 
     $                 +sigy/permittivity(j)*pmldn(j,2)-sigx*en(j,2)

         reshn(j,3)=reshn(j,3)+respmlhn(j,3)*bm1(i,1,1,e) 
         resen(j,1)=resen(j,1)+respmlen(j,1)*bm1(i,1,1,e) 
         resen(j,2)=resen(j,2)+respmlen(j,2)*bm1(i,1,1,e) 

         enddo
        enddo

       elseif (IFTM) then

        do ie=1,maxpml
          e =pmlptr(ie) 

         do i = 1,n                

         j=i+nxyz*(e-1)
         bminv= 1./bm1(i,1,1,e)
         sigx = pmlsigma(j,1)
         sigy = pmlsigma(j,2)
         sigz = 0

         respmlbn(j,1)=reshn(j,1)*bminv-sigy/permittivity(j)*pmlbn(j,1)
         respmlbn(j,2)=reshn(j,2)*bminv-sigz/permittivity(j)*pmlbn(j,2)
         respmldn(j,3)=resen(j,3)*bminv-sigx/permittivity(j)*pmldn(j,3)


         respmlhn(j,1)=-sigy/permittivity(j)*pmlbn(j,1)
     $                 +sigx/permittivity(j)*pmlbn(j,1)
     $                 -sigz*permeability(j)/permittivity(j)*hn(j,1)
         respmlhn(j,2)=-sigz/permittivity(j)*pmlbn(j,2)
     $                 +sigy/permittivity(j)*pmlbn(j,2)
     $                 -sigx*permeability(j)/permittivity(j)*hn(j,2)


         respmlen(j,3)=-sigx/permittivity(j)*pmldn(j,3)
     $                 +sigz/permittivity(j)*pmldn(j,3)-sigy*en(j,3)

         reshn(j,1)=reshn(j,1)+respmlhn(j,1)*bm1(i,1,1,e)
         reshn(j,2)=reshn(j,2)+respmlhn(j,2)*bm1(i,1,1,e)
         resen(j,3)=resen(j,3)+respmlen(j,3)*bm1(i,1,1,e)

         enddo
        enddo
       else
         if (nid.eq.0) write(6,*) 'pml_step: IFTE or IFTM:not defined'
         call exitt
       endif

      endif

      return
      end
