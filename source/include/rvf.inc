c $Id: rvf.inc,v 1.2 2005/05/26 11:54:42 hal Exp $
c
c     coordinates, velocities and forces
c

      double precision x, y, z, vx, vy, vz, fx, fy, fz
      double precision xo, yo, zo, vxo, vyo, vzo, fxo, fyo, fzo

      common /rvf/ x(mxatom), y(mxatom), z(mxatom),
     $    vx(mxatom), vy(mxatom), vz(mxatom),
     $    fx(mxatom), fy(mxatom), fz(mxatom),
     $    xo(mxatom), yo(mxatom), zo(mxatom),
     $    vxo(mxatom), vyo(mxatom), vzo(mxatom),
     $    fxo(mxatom), fyo(mxatom), fzo(mxatom)

      save /rvf/


c     local variables:
c     mode: fortran
c     end:
