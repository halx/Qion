      subroutine pf3bdc (IdxLst, NoNbrs, Epot)
c $Id: pf3bdc.f,v 1.11 2005/06/13 10:38:34 hal Exp $
c
c     pf3bdc: calculate three-body correction forces between ion and water
c

      implicit none

      integer IdxLst(*), NoNbrs
      double precision Epot

      integer i, j, l, k
      double precision xion1, yion1, zion1, xO2, yO2, zO2, xO3, yO3, zO3
      double precision r12, r13, r23, r12kp1, r13kp1, r23kp1
      double precision x12, y12, z12, xe12, ye12, ze12
      double precision x13, y13, z13, xe13, ye13, ze13
      double precision x23, y23, z23, xe23, ye23, ze23
      double precision V, Vrep, Vattr, V2fsh1, V2fsh2, fshift
      double precision fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'qmmm.inc'
      include 'press.inc'



      xion1 = x(nwc)
      yion1 = y(nwc)
      zion1 = z(nwc)

c     iterate over the ion's oxygen neighbours
      do k = 1, NoNbrs - 1
        i = IdxLst(k)

        xO2 = x(i)
        yO2 = y(i)
        zO2 = z(i)

c       ion1-O2 distance
        x12 = xion1 - xO2
        y12 = yion1 - yO2
        z12 = zion1 - zO2

        call pbcmic(x12, y12, z12)

        r12 = sqrt(x12**2 + y12**2 + z12**2)

        if (r12 .le. r3bd) then
          xe12 = x12 / r12
          ye12 = y12 / r12
          ze12 = z12 / r12

          do l = k + 1, NoNbrs
            j = IdxLst(l)

            if (iqmmm .and. inqm(i) .and. inqm(j)) goto 100

            xO3 = x(j)
            yO3 = y(j)
            zO3 = z(j)

c           ion1-O3 distance
            x13 = xion1 - xO3
            y13 = yion1 - yO3
            z13 = zion1 - zO3

            call pbcmic(x13, y13, z13)

            r13 = sqrt(x13**2 + y13**2 + z13**2)

            if (r13 .le. r3bd) then

c             O2-O3 distance
              x23 = xO2 - xO3
              y23 = yO2 - yO3
              z23 = zO2 - zO3

              call pbcmic(x23, y23, z23)

              r23 = sqrt(x23**2 + y23**2 + z23**2)

              xe13 = x13 / r13
              ye13 = y13 / r13
              ze13 = z13 / r13

              xe23 = x23 / r23
              ye23 = y23 / r23
              ze23 = z23 / r23

c             potential energy
              Vrep = a3bd *
     $            exp (-b3bd * r12) *
     $            exp (-b3bd * r13) *
     $            exp (-c3bd * r23)

              Vattr = -d3bd / r12**k3bd - d3bd / r13**k3bd -
     $            e3bd / r23**k3bd

              V = Vrep + Vattr


c             helper variables
              fshift = (r3bd - r12)**2 * (r3bd - r13)**2

              V2fsh1 = 2 * V * (r3bd - r12)**2 * (r3bd - r13)
              V2fsh2 = 2 * V * (r3bd - r12)    * (r3bd - r13)**2

              r12kp1 = r12**(k3bd+1)
              r13kp1 = r13**(k3bd+1)
              r23kp1 = r23**(k3bd+1)


c             shifted potential
              Epot = Epot + V * fshift

c             force exerted on ion1
              fx1 = (b3bd * Vrep * (xe12 + xe13) -
     $            kD3bd * (xe12 / r12kp1 + xe13 / r13kp1)) * fshift +
     $            V2fsh2 * xe12 + V2fsh1 * xe13
              fy1 = (b3bd * Vrep * (ye12 + ye13) -
     $            kD3bd * (ye12 / r12kp1 + ye13 / r13kp1)) * fshift +
     $            V2fsh2 * ye12 + V2fsh1 * ye13
              fz1 = (b3bd * Vrep * (ze12 + ze13) -
     $            kD3bd * (ze12 / r12kp1 + ze13 / r13kp1)) * fshift +
     $            V2fsh2 * ze12 + V2fsh1 * ze13

              fx(nwc) = fx(nwc) + fx1
              fy(nwc) = fy(nwc) + fy1
              fz(nwc) = fz(nwc) + fz1


c             force exerted on O2
              fx2 = (-Vrep * (b3bd * xe12 - c3bd * xe23) +
     $            kD3bd * xe12 / r12kp1 - kE3bd * xe23 / r23kp1) *
     $            fshift - V2fsh2 * xe12
              fy2 = (-Vrep * (b3bd * ye12 - c3bd * ye23) +
     $            kD3bd * ye12 / r12kp1 - kE3bd * ye23 / r23kp1) *
     $            fshift - V2fsh2 * ye12
              fz2 = (-Vrep * (b3bd * ze12 - c3bd * ze23) +
     $            kD3bd * ze12 / r12kp1 - kE3bd * ze23 / r23kp1) *
     $            fshift - V2fsh2 * ze12

              fx(i) = fx(i) + fx2
              fy(i) = fy(i) + fy2
              fz(i) = fz(i) + fz2


c             force exerted on O3
c             because f1 + f2 + f3 = 0:
              fx3 = -fx1 - fx2
              fy3 = -fy1 - fy2
              fz3 = -fz1 - fz2

              fx(j) = fx(j) + fx3
              fy(j) = fy(j) + fy3
              fz(j) = fz(j) + fz3

c             need two forces and two distances NOT joining these force
c             centers for virial
              virial = virial +
     $            fx1 * x13 + fy1 * y13 + fz1 * z13 +
     $            fx2 * x23 + fy2 * y23 + fz2 * z23
            endif

  100       continue
          enddo
        endif
      enddo

      end
