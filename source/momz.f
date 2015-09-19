      subroutine momz
c $Id: momz.f,v 1.2 2005/05/26 11:54:42 hal Exp $
c
c     momz:  set total linear momentum of system to zero
c

      implicit none

      external ssum
      double precision ssum

      integer i
      double precision pxsum, pysum, pzsum, vxcorr, vycorr, vzcorr, mtot

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'



c     calculate total linear momentum
      pxsum = amo * ssum(vx, 1, no)
      pysum = amo * ssum(vy, 1, no)
      pzsum = amo * ssum(vz, 1, no)

      pxsum = pxsum + amh * ssum(vx, no + 1, nw)
      pysum = pysum + amh * ssum(vy, no + 1, nw)
      pzsum = pzsum + amh * ssum(vz, no + 1, nw)

      pxsum = pxsum + amion * ssum(vx, nw + 1, nwc)
      pysum = pysum + amion * ssum(vy, nw + 1, nwc)
      pzsum = pzsum + amion * ssum(vz, nw + 1, nwc)

c     calculate total mass
      mtot = dble(no) * amo + dble(nh) * amh + dble(nc) * amion

c     calculate velocity vector for total system
      vxcorr = pxsum / mtot
      vycorr = pysum / mtot
      vzcorr = pzsum / mtot

c     set total momentum to zero
      do i = 1, nwc
        vx(i) = vx(i) - vxcorr
        vy(i) = vy(i) - vycorr
        vz(i) = vz(i) - vzcorr

        vxo(i) = vxo(i) - vxcorr
        vyo(i) = vyo(i) - vycorr
        vzo(i) = vzo(i) - vzcorr
      enddo

      end
