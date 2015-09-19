      subroutine pbath(Psys)
c $Id: pbath.f,v 1.3 2005/06/09 14:02:00 hal Exp $
c
c     bath: weak coupling method for pressure control, see
c           Berendsen et al., J Chem Phys 81(1984), 3684

      implicit none

      double precision Psys

      integer i
      double precision mu

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'units.inc'
      include 'press.inc'



c     apply the coupling bath
      mu = 1.0D0 - (dt * Pcomp * (Popt - Psys) ) / (3 * taup)

c     scale box edges and recalculate volume
      boxx = mu * boxx
      boxy = mu * boxy
      boxz = mu * boxz
      Volume = boxx * boxy * boxz / 1.0D30

c     scale positions
      do i = 1, nwc
        x(i) = mu * x(i)
        y(i) = mu * y(i)
        z(i) = mu * z(i)
      enddo

      end
