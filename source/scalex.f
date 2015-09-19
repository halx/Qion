      subroutine scalex (Toxy, Thydr)
c $Id: scalex.f,v 1.2 2005/05/26 11:54:42 hal Exp $
c
c     scalex: "hard" temperature scaling for H2O intended for the
c             beginning (first nscale steps) of a simulation of a
c             highly strained system
c

      implicit none

      double precision Toxy, Thydr

      integer i
      logical scaled
      double precision scafac, dTmax
      parameter (dTmax = 10.0D0)

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'units.inc'



      scaled = .false.

c     scale the oxygen temperature
      if (Toxy .gt. (Topt + dTmax) ) then
        scafac = sqrt (Topt / Toxy)

        do i = 1, no
          vx(i) = scafac * vx(i)
          vy(i) = scafac * vy(i)
          vz(i) = scafac * vz(i)

          vxo(i) = scafac * vxo(i)
          vyo(i) = scafac * vyo(i)
          vzo(i) = scafac * vzo(i)
        enddo

        scaled = .true.
      endif

c     scale the hydrogen temperature
      if (Thydr .gt. (Topt + dTmax) ) then
        scafac = sqrt (Topt / Thydr)

        do i = no + 1, nw
          vx(i) = scafac * vx(i)
          vy(i) = scafac * vy(i)
          vz(i) = scafac * vz(i)

          vxo(i) = scafac * vxo(i)
          vyo(i) = scafac * vyo(i)
          vzo(i) = scafac * vzo(i)
        enddo

        scaled = .true.
      endif

c     set total momentum of system to zero, if scaling was performed
      if (scaled) call momz

      end
