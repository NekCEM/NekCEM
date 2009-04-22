c-----------------------------------------------------------------------
      real function bnch_z (tt,zz)
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'BEAM'

      bnch_z = zz-cspeed*tt
      bnch_z = bnch_z-zposition

      return
      end
c-----------------------------------------------------------------------
      real function bnch_rdensity(rr)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'BEAM'

      bnch_normal  = 1.0/(sqrt(2.0*pi)*rsigma)
      bnch_rdensity= exp(-0.5*(rr/rsigma)**2)*bnch_normal

      return
      end
c-----------------------------------------------------------------------
      real function bnch_zdensity(zz)

      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'
      include 'NEKUSE'
      include 'BEAM'
      integer  icalld
      save     icalld
      data     icalld/0/

      real     bnch_normal
      save     bnch_normal

      bnch_normal  = 1.0/(sqrt(2.0*pi)*zsigma)
      zz  = zz-cspeed*rktime-zposition                    
      bnch_zdensity= exp(-0.5*(zz/zsigma)**2)*bnch_normal
      !write(6,*) 'bnch---',zz,bnch_zdensity,bnch_normal ,(zz/zsigma)**2

      return
      end
c-----------------------------------------------------------------------
