      subroutine pfintr
c $Id: pfintr.f,v 1.10 2005/06/13 10:38:34 hal Exp $
c
c     pfintr: calculate the intra-molecular potential energy and forces
c             of H2O, see
c             P. Bopp et. al, Chem. Phys. Letters 98(1983), 129--133
c
c             coded(?) by J. Seitz-Beywl 07.06.89,
c                         Josef Boecker  09.03.94
c

      implicit none

      integer i, j, k1, k2
      double precision xo1, yo1, zo1, xh1, yh1, zh1, xh2, yh2, zh2
      double precision xoh1, xoh2, yoh1, yoh2, zoh1, zoh2
      double precision sqdst1, sqdst2, dist1, dist2
      double precision xhh, yhh, zhh, sqdthh, cosang, sinang, angle
      double precision rho1, rho2, dalpha
      double precision xeoh1, yeoh1, zeoh1, xeoh2, yeoh2, zeoh2
      double precision st2x, st2y, st2z, st3x, st3y
      double precision st3z, ffan, ffh1, ffh2
      double precision fxh1, fyh1, fzh1, fxh2, fyh2, fzh2
      double precision fxo1, fyo1, fzo1

      include 'params.inc'
      include 'consts.inc'
      include 'sizes.inc'
      include 'intra.inc'
      include 'rvf.inc'
      include 'qmmm.inc'
      include 'energy.inc'
      include 'press.inc'



      do i = 1, no
        if (iqmmm .and. inqm(i)) goto 100

        k1 = no + 2 * i - 1
        k2 = k1 + 1

        ffh1 = 0.0D0
        ffh2 = 0.0D0
        ffan = 0.0D0

        xo1 = x(i)
        yo1 = y(i)
        zo1 = z(i)

        xh1 = x(k1)
        yh1 = y(k1)
        zh1 = z(k1)

        xh2 = x(k2)
        yh2 = y(k2)
        zh2 = z(k2)

c       calculate the distance components between O and H1
        xoh1 = xo1 - xh1
        yoh1 = yo1 - yh1
        zoh1 = zo1 - zh1

        call pbcmic(xoh1, yoh1, zoh1)

c       calculate the distance between O and H1
        sqdst1 = xoh1**2 + yoh1**2 + zoh1**2
        dist1 = sqrt(sqdst1)

c       calculate the distance components between O and H2
        xoh2 = xo1 - xh2
        yoh2 = yo1 - yh2
        zoh2 = zo1 - zh2

        call pbcmic(xoh2, yoh2, zoh2)

c       calculate the distance between O and H2
        sqdst2 = xoh2**2 + yoh2**2 + zoh2**2
        dist2 = sqrt(sqdst2)

c       calculate the distance components between H1 and H2
        xhh = xh1 - xh2
        yhh = yh1 - yh2
        zhh = zh1 - zh2

        call pbcmic(xhh, yhh, zhh)

c       calculate the square distance between H1 and H2
        sqdthh = xhh**2 + yhh**2 + zhh**2

c       cosine law used to calculate the current H-O-H angle,
c       we assume that sqdthh is always the hypothenuse
        cosang = (sqdthh - sqdst1 - sqdst2) / (-2.0D0 * dist1 * dist2)
        angle = acos(cosang)
        sinang = sin(angle)

        rho1 = (dist1 - doh0) / dist1
        rho2 = (dist2 - doh0) / dist2
        dalpha = angle - alpha0

c       calculate the energy and the forces
        do j = 1, npotp
          Eintr = Eintr + ppintr(j) * rho1**irho1(j) *
     $        rho2**irho2(j) * dalpha**ialpha(j)

          ffh1 = ffh1 - ppintr(j) * rho2**irho2(j) * dalpha**ialpha(j) *
     $        irho1(j) * rho1**(irho1(j) - 1) * (doh0 / sqdst1)
          ffh2 = ffh2 - ppintr(j) * rho1**irho1(j) * dalpha**ialpha(j) *
     $        irho2(j) * rho2**(irho2(j) - 1) * (doh0 / sqdst2)
          ffan = ffan - ppintr(j) * rho1**irho1(j) * rho2**irho2(j) *
     $        ialpha(j) * dalpha**(ialpha(j) - 1)
        enddo


c
c       calculate the B-matrix elements (Wilson, Decius & Cross p. 57)
c

c       unity vector r1(O-H)
        xeoh1 = xoh1 / dist1
        yeoh1 = yoh1 / dist1
        zeoh1 = zoh1 / dist1

c       unity vector r2(O-H)
        xeoh2 = xoh2 / dist2
        yeoh2 = yoh2 / dist2
        zeoh2 = zoh2 / dist2

        st2x = (xeoh2 - cosang * xeoh1) / (sinang * dist1)
        st2y = (yeoh2 - cosang * yeoh1) / (sinang * dist1)
        st2z = (zeoh2 - cosang * zeoh1) / (sinang * dist1)

        st3x = (xeoh1 - cosang * xeoh2) / (sinang * dist2)
        st3y = (yeoh1 - cosang * yeoh2) / (sinang * dist2)
        st3z = (zeoh1 - cosang * zeoh2) / (sinang * dist2)

        fxh1 = -ffh1 * xeoh1 + ffan * st2x
        fyh1 = -ffh1 * yeoh1 + ffan * st2y
        fzh1 = -ffh1 * zeoh1 + ffan * st2z

        fxh2 = -ffh2 * xeoh2 + ffan * st3x
        fyh2 = -ffh2 * yeoh2 + ffan * st3y
        fzh2 = -ffh2 * zeoh2 + ffan * st3z

c       f1 + f2 + f3 = 0 and hence:
        fxo1 = -fxh1 - fxh2
        fyo1 = -fyh1 - fyh2
        fzo1 = -fzh1 - fzh2

c       store three-body forces
        fx(i) = fx(i) + fxo1
        fy(i) = fy(i) + fyo1
        fz(i) = fz(i) + fzo1

        fx(k1) = fx(k1) + fxh1
        fy(k1) = fy(k1) + fyh1
        fz(k1) = fz(k1) + fzh1

        fx(k2) = fx(k2) + fxh2
        fy(k2) = fy(k2) + fyh2
        fz(k2) = fz(k2) + fzh2

c       need two forces and two distances NOT joining these force
c       centers for virial
        virial = virial +
     $      fxo1 * xoh2 + fyo1 * yoh2 + fzo1 * zoh2 +
     $      fxh1 * xhh  + fyh1 * yhh  + fzh1 * zhh

  100   continue
      enddo

      end
