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

      real sum_here

      ! implement marching inside each element
      do elt=1,nelv
        do axis=0,ldim-1
          posface = eface(2+2*axis)
          negface = eface(1+2*axis)
          do idx = 1,lx1*lz1
            sum_here = faceary(idx,posface,elt)+faceary(idx,negface,elt)
            if (sum_here.ne.0) then
              faceary(idx,posface,elt) = faceary(idx,posface,elt) + 1
              faceary(idx,negface,elt) = faceary(idx,negface,elt) + 1
            endif
          enddo
        enddo
      enddo

      ! maximum across neighboring elements
      call fgs_gop(gsh_face, faceary, 'M\0')

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

      real faceary(lx1*lz1,2*ldim,lelt)
      integer thick

      integer elt, face, i, axis, ix, iz

      character*3 bc_here

      call rzero(faceary, lx1*lz1*2*ldim*lelt)

      ! produce a face array that has ones at the *interior*
      ! of all PML faces.
      ! (We need to avoid the edges since those will march off 
      ! parallel to the face, which we do not want.)

      do elt=1,nelv
        do face=1,2*ndim
          bc_here = cbc(face, elt, 2)

          if  (bc_here.eq.'PML' .or. bc_here.eq.'pml') then
            if (if3d) then
              do ix=1,lx1-2
              do iz=1,lz1-2
                faceary(ix+iz*lx1+1, face, elt) = 1
              enddo
              enddo
            else
              do ix=2,lx1-1
                faceary(ix, face, elt) = 1
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
      integer ind_ptidx 

      real far_here, far_opp

      real pmlinf
      parameter (pmlinf = 1e20)

      ! table of opposing face indices, symm/pff convention
      integer oppface(6) 
      data oppface /2,1,4,3,6,5/

      integer elt, face, axis, globface, globaxis
      integer globsign
      real mincoord, maxcoord
      real gop_work

      real minx, maxx, miny, maxy, minz, maxz

      ! BEGIN CODE

      ! determine the indicative point index (see above)

      if (if3d) then
        ind_ptidx = 1+1*lx1+1
      else
        ind_ptidx = 2
      endif

      ! initialize outer, inner to "infinities"

      do axis=0,ndim-1
        outer((axis*2)+1) = pmlinf
        inner((axis*2)+1) = -pmlinf
        inner((axis*2)+2) = pmlinf
        outer((axis*2)+2) = -pmlinf
      enddo

      ! fill outer, inner, eltag for local processor

      do elt=1,nelv
        eltag(elt) = 0
      enddo

      do elt=1,nelv
        do axis=1,ndim
          face = (axis-1)*2 + 1 ! pff/sym convention

          far_here = faceary(ind_ptidx, eface(face), elt)
          far_opp = faceary(ind_ptidx, eface(oppface(face)), elt)

          if ((far_here.ne.0).and.(far_opp.ne.0)) then
            ! element number `elt' is inside the PML

            if (far_here.eq.far_opp) then
              write (*,*) 'No "gradient" in PML indicators.'
              write (*,*) 'This should not happen, unless your PMLs'
              write (*,*) 'are colliding in the middle.'
              write (*,*) 'Aborting.'
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

      real inner(2*ldim),outer(2*ldim)
      integer eltag(lelt)
      real order, referr

      integer elt, i, face, axis
      real point(3), zero2one
      real eta, sigmamax
      real width

      call rzero(pmlsigma, lx1*ly1*lz1*lelt*ldim)

      do elt=1,nelt
        if (eltag(elt).ne.0) then
          ! there is some PML here

          do face=1,2*ndim ! pff/sym convention
            axis = (face-1)/2 + 1

            width = abs(outer(face)-inner(face))

            if (iand(eltag(elt),ishft(1,face-1)).ne.0) then
              ! there is a PML in the direction of `face' here
              do i=1,nx1*ny1*nz1

                ss  = mu(i,1,1,elt)/epsilo(i,1,1,elt)
                ss  = sqrt(ss)
                eta = ss           !1 or  FIXME: should be sqrt(mu/epsilon)
               !eta = 1.0 !or  FIXME: should be sqrt(mu/epsilon)
                sigmamax = -(order+1)*log(referr)/(2*eta*width)

                point(1) = xm1(i,1,1,elt)
                point(2) = ym1(i,1,1,elt)
                point(3) = zm1(i,1,1,elt)

                zero2one = (point(axis)-inner(face))/
     $                     (outer(face)-inner(face))

                pmlsigma(i,1,1,elt,axis) = sigmamax * zero2one**order
              enddo
            endif
          enddo
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine pml_zero_fields
      include 'SIZE'
      include 'TOTAL'
      include 'PML'

      integer n

      n = lx1*ly1*lz1*lelt

      call rzero(pmlbx, n)
      call rzero(pmlby, n)
      call rzero(pmlbz, n)
      call rzero(pmldx, n)
      call rzero(pmldy, n)
      call rzero(pmldz, n)

      return
      end
c-----------------------------------------------------------------------
      subroutine pml_setup

      include 'SIZE'
      include 'TOTAL'
      include 'PML'
      include 'SCRATCH'

      integer pmlthick
      real pmlorder
      real pmlreferr 
      ! how much of the incident signal is permitted to be reflected
      ! by the PML. See the term R(0) in formula (7.57) of the Taflove
      ! book.

      pmlthick = int(param(77))
      pmlorder = param(78)
      pmlreferr = param(79)

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
      call pml_zero_fields

      !write(6,*) 'pml',pmlinner(1)
      !write(6,*) 'pml',pmlinner(2)
      !write(6,*) 'pml',pmlinner(3)
      !write(6,*) 'pml',pmlinner(4)
      !write(6,*) 'pml',pmlinner(5)
      !write(6,*) 'pml',pmlinner(6)
      !write(6,*) 'pml',pmlinner(7)
      !call exitt
          
      ! for debugging:

      !call face2full(fullary, faceary)
      !call vtk_dump_field(fullary, 'PML   ', 0, 0.)

      ! call intelt2full(fullary, pmltag)
      ! call vtk_dump_field(fullary, 'PMLTAG', 0, 0.)
      ! call exitt

      ! call vtk_dump_field(pmlsigma(1,1,1,1,1), 'PMLSI1', 0, 0.)
      ! call vtk_dump_field(pmlsigma(1,1,1,1,2), 'PMLSI2', 0, 0.)
      ! if (if3d) then
        ! call vtk_dump_field(pmlsigma(1,1,1,1,3), 'PMLSI3', 0, 0.)
      ! endif

      ! call exitt
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

      write (*,*) 'max',absmax,' reached at element', maxelt,
     $ ' index', maxi

      return
      end

c-----------------------------------------------------------------------
      subroutine pml_step

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'PML'

      integer elt, i

      real sigx, sigy, sigz
      real inv_mass
      real x,y,z

      ! FIXME This PML does not even attempt to worry about
      ! spatially varying material coefficients.
      ! permit_0 and permea_0 are taken into account, though.

      ! FIXME this loop might cause a good bit of cache thrashing.
      ! If this turns out to be a perf bottleneck, just break it
      ! up into a sequence of smaller loops or *yuck* even use mat1.f
      ! functions.

      ! call find_peak(hx)

      do elt=1,nelt
        if (pmltag(elt).ne.0) then ! is there a PML here?
          do i=1,lx1*ly1*lz1
            sigx = pmlsigma(i,1,1,elt,1)
            sigy = pmlsigma(i,1,1,elt,2)
            if (if3d) then
              sigz = pmlsigma(i,1,1,elt,3)
            else
              sigz = 0
            endif

            x = xm1(i,1,1,elt)
            y = ym1(i,1,1,elt)
            z = zm1(i,1,1,elt)

            ! see doc/pml.tm to make sense of this.

            ! Note that we assume that at the point where we get called, 
            ! rese[xyz] = curl H
            ! resh[xyz] = - curl E

            inv_mass = 1./bm1(i,1,1,elt)

            if(hxactive) respmlbx(i,1,1,elt) = 
     $        inv_mass*reshx(i,1,1,elt)
     $        -sigy/permit_0*pmlbx(i,1,1,elt) 
            ! cem_invqmass for only reshx term
            if(hyactive) respmlby(i,1,1,elt) = 
     $        inv_mass*reshy(i,1,1,elt)
     $        -sigz/permit_0*pmlby(i,1,1,elt)
            if(hzactive) respmlbz(i,1,1,elt) = 
     $        inv_mass*reshz(i,1,1,elt)
     $        -sigx/permit_0*pmlbz(i,1,1,elt)

            if(exactive) respmldx(i,1,1,elt) = 
     $        inv_mass*resex(i,1,1,elt)
     $        -sigy/permit_0*pmldx(i,1,1,elt)
            if(eyactive) respmldy(i,1,1,elt) = 
     $        inv_mass*resey(i,1,1,elt)
     $        -sigz/permit_0*pmldy(i,1,1,elt)
            if(ezactive) respmldz(i,1,1,elt) = 
     $        inv_mass*resez(i,1,1,elt)
     $        -sigx/permit_0*pmldz(i,1,1,elt)

            if(hxactive) respmlhx(i,1,1,elt) = 
     $        -sigy*permit_0/permea_0*pmlbx(i,1,1,elt) 
     $        +sigx/permit_0/permea_0*pmlbx(i,1,1,elt)
     $        -sigz/permit_0*hx(i,1,1,elt)
            if(hyactive) respmlhy(i,1,1,elt) = 
     $        -sigz*permit_0/permea_0*pmlby(i,1,1,elt) 
     $        +sigy/permit_0/permea_0*pmlby(i,1,1,elt)
     $        -sigx/permit_0*hy(i,1,1,elt)
            if(hzactive) respmlhz(i,1,1,elt) = 
     $        -sigx*permit_0/permea_0*pmlbz(i,1,1,elt) 
     $        +sigz/permit_0/permea_0*pmlbz(i,1,1,elt)
     $        -sigy/permit_0*hz(i,1,1,elt)

            if(exactive) respmlex(i,1,1,elt) = 
     $         -sigy/permit_0**2 * pmldx(i,1,1,elt)
     $         +sigx/permit_0**2 * pmldx(i,1,1,elt)
     $         -sigz/permit_0*ex(i,1,1,elt)
            if (eyactive) respmley(i,1,1,elt) = 
     $         -sigz/permit_0**2 * pmldy(i,1,1,elt)
     $         +sigy/permit_0**2 * pmldy(i,1,1,elt)
     $         -sigx/permit_0*ey(i,1,1,elt)
            if(ezactive) respmlez(i,1,1,elt) = 
     $         -sigx/permit_0**2 * pmldz(i,1,1,elt)
     $         +sigz/permit_0**2 * pmldz(i,1,1,elt)
     $         -sigy/permit_0*ez(i,1,1,elt)

          enddo
        endif
      enddo

      !call vtk_dump_field(pmlbx(1,1,1,1), 'PMLBX ', istep, time)
      !call vtk_dump_field(pmlby(1,1,1,1), 'PMLBY ', istep, time)
      !call vtk_dump_field(pmlbz(1,1,1,1), 'PMLBZ ', istep, time)
      !call vtk_dump_field(pmldx(1,1,1,1), 'PMLDX ', istep, time)
      !call vtk_dump_field(pmldy(1,1,1,1), 'PMLDY ', istep, time)
      !call vtk_dump_field(pmldz(1,1,1,1), 'PMLDZ ', istep, time)
      
      return
      end
