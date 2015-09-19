      subroutine vverl2
c $Id: vverl2.f,v 1.6 2005/06/10 10:46:41 hal Exp $
c
c     velocity Verlet part 2, see [AT87; pp. 87]
c

      implicit none

      integer i

      include 'params.inc'
      include 'sizes.inc'
      include 'eqmot.inc'
      include 'rvf.inc'



      do i = 1, nwc
        vx(i) = vx(i) + dt2i * fx(i) / mass(i)
        vy(i) = vy(i) + dt2i * fy(i) / mass(i)
        vz(i) = vz(i) + dt2i * fz(i) / mass(i)

        if (iwrap) call pbcmic(x(i), y(i), z(i))
      enddo

      end
