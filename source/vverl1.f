      subroutine vverl1
c $Id: vverl1.f,v 1.3 2005/05/24 10:49:24 hal Exp $
c
c     velocity Verlet part 1, see [AT87; pp. 87]
c

      implicit none

      integer i

      include 'params.inc'
      include 'sizes.inc'
      include 'eqmot.inc'
      include 'rvf.inc'



      do i = 1, nwc
        x(i) = x(i) + dtv * vx(i) + dtsq2 * fx(i) / mass(i)
        y(i) = y(i) + dtv * vy(i) + dtsq2 * fy(i) / mass(i)
        z(i) = z(i) + dtv * vz(i) + dtsq2 * fz(i) / mass(i)

        vx(i) = vx(i) + dt2i * fx(i) / mass(i)
        vy(i) = vy(i) + dt2i * fy(i) / mass(i)
        vz(i) = vz(i) + dt2i * fz(i) / mass(i)
      enddo

      end
