      subroutine wrtraj
c $Id: wrtraj.f,v 1.7 2005/06/21 11:48:13 hal Exp $
c
c     wrtraj: write trajectories and velocities if requested
c

      implicit none

      integer i

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'stpcnt.inc'
      include 'units.inc'
      include 'files.inc'



c     write coordinates
      if (ftraj .ne. ' ') then
        write (utraj, 8000) (x(i), y(i), z(i), i = 1, nwc)
        write (utraj, 8000) boxx, boxy, boxz

        call flush (utraj)
      endif

c     write velocities
      if (fveloc .ne. ' ') then
        write (uveloc, 8010) (vx(i), vy(i), vz(i), i = 1, nwc)
        call flush (uveloc)
      endif

c     write forces along trajectory
      if (fftraj .ne. ' ') then
        write (uftraj, 8010) (fx(i), fy(i), fz(i), i = 1, nwc)
        call flush (uftraj)
      endif


 8000 format (10(1x, f7.3))
 8010 format (6(1x, f12.3))

      end
