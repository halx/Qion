      subroutine corr
c $Id: corr.f,v 1.6 2005/06/10 10:46:41 hal Exp $
c
c     corr: corrector method of Adams-Bashforth (sure?)
c
c     (this routine was written by Lutz Schaefer and copied by
c      Ecki, 1.4.1990)
c

      implicit none

      integer i

      include 'params.inc'
      include 'sizes.inc'
      include 'eqmot.inc'
      include 'rvf.inc'



      do i = 1, nwc
        vx(i) = vxo(i) + corp3(i) * (fxo(i) + fx(i))
        vy(i) = vyo(i) + corp3(i) * (fyo(i) + fy(i))
        vz(i) = vzo(i) + corp3(i) * (fzo(i) + fz(i))

        x(i) = xo(i) + 0.5D0 * corp1 * (vx(i) + vxo(i) - corp2(i) *
     $      (fx(i) - fxo(i)))
        y(i) = yo(i) + 0.5D0 * corp1 * (vy(i) + vyo(i) - corp2(i) *
     $      (fy(i) - fyo(i)))
        z(i) = zo(i) + 0.5D0 * corp1 * (vz(i) + vzo(i) - corp2(i) *
     $      (fz(i) - fzo(i)))

        if (iwrap) call pbcmic (x(i), y(i), z(i))
      enddo

      end
