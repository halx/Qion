      subroutine bath(Tsys)
c $Id: bath.f,v 1.5 2005/06/21 11:48:13 hal Exp $
c
c     bath: weak coupling method for temperature control, see
c           Berendsen et al., J Chem Phys 81(1984), 3684

      implicit none

      double precision Tsys

      integer i
      double precision lambda

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'units.inc'



c     apply the coupling bath
      lambda = sqrt(1.0D0 + dt / taut * (Topt / Tsys - 1.0D0))

      do i = 1, nwc
        vx(i) = lambda * vx(i)
        vy(i) = lambda * vy(i)
        vz(i) = lambda * vz(i)
      enddo

      end
