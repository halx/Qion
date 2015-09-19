      subroutine pbcmic (x, y, z)
c $Id: pbcmic.f,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
c
c     pbcmic: apply periodic boundary conditions/minimum image convention
c

      implicit none

      double precision x, y, z

      include 'params.inc'



      x = x - boxx * anint (x / boxx)
      y = y - boxy * anint (y / boxy)
      z = z - boxz * anint (z / boxz)

      end
