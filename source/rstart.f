      subroutine rstart
c $Id: rstart.f,v 1.3 2005/06/07 13:43:26 hal Exp $
c
c     rstart: write the restart coordinates, velocities, and forces
c             write coordinates in xyz format
c

      implicit none

      integer i

      include 'params.inc'
      include 'stpcnt.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'units.inc'
      include 'files.inc'



      open (ustart, err = 500, file = frst, status = 'unknown')
      rewind (ustart)

      write (ustart, 8000) nfi, boxx, boxy, boxz, no, nh, nc
      write (ustart, 8010) (xo(i), yo(i), zo(i), i = 1, nwc)
      write (ustart, 8020) (vxo(i), vyo(i), vzo(i), i = 1, nwc)
      write (ustart, 8020) (fxo(i), fyo(i), fzo(i), i = 1, nwc)
      write (ustart, *)
      write (ustart, 8010) (x(i), y(i), z(i), i = 1, nwc)
      write (ustart, 8020) (vx(i), vy(i), vz(i), i = 1, nwc)
      write (ustart, 8020) (fx(i), fy(i), fz(i), i = 1, nwc)

      close (ustart)

      return


  500 write (stderr, *) 'Could not open: ', frst
      call exit (1)

 8000 format (1x, i9, 3(1x, f10.7), 3(1x, i4) )
 8010 format (3(1x, f20.16) )
 8020 format (3(1x, f24.13) )

      end
