      subroutine pfwthh
c $Id: pfwthh.f,v 1.12 2005/06/21 11:48:13 hal Exp $
c
c     pfwthh: calculate the potential energy and forces between
c             hydrogen atoms of CF2 water, see
c             F. H. Stillinger and A. Rahman, JPC 68(1978), 666--670
c

      implicit none

      integer i, j, k
      double precision xi, yi, zi, rxij, ryij, rzij, rij, rsqij
      double precision fxij, fyij, fzij, fijr, Force
      double precision VnonC, VCoul, VCsum, VnCsum
      double precision hlpHH4, hlpHH7, hexHH1, hexHH2

      include 'params.inc'
      include 'potfor.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'
      include 'energy.inc'
      include 'press.inc'
      include 'qmmm.inc'



      VCsum = 0.0D0
      VnCsum = 0.0D0

      do k = 1, nonHH
        if (nbrHH(k) .lt. 0) then
          i = -nbrHH(k)

          xi = x(i)
          yi = y(i)
          zi = z(i)
        else
          j = nbrHH(k)

          if (iqmmm .and. inqm(i) .and. inqm(j)) goto 100

c         calculate the distance components between i and j particle
          rxij = xi - x(j)
          ryij = yi - y(j)
          rzij = zi - z(j)

          call pbcmic(rxij, ryij, rzij)

c         calculate the distance between i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         calculate Coulomb potential and forces for distances < rcut
          if (rsqij .le. rsqcut) then
            rij = sqrt(rsqij)

c           except in the case of point charges
            if (ipntch .and. (inqm(i) .neqv. inqm(j)) ) goto 110

            VCoul = parHH(1) / rij
            Force = parHH(1) / rsqij

c           calculate the shifted-force potential [AT87; pp. 145] and
c           the shifted force
            VCsum = VCsum + (VCoul - vccbhh - dccbhh * (rij - rcut) )
            fijr = (Force + dccbhh) / rij
            virial = virial + fijr * rsqij

c           force components
            fxij = fijr * rxij
            fyij = fijr * ryij
            fzij = fijr * rzij

c           sum force components for current atom i
            fx(i) = fx(i) + fxij
            fy(i) = fy(i) + fyij
            fz(i) = fz(i) + fzij

c           subtract force components from neighbour atom j
            fx(j) = fx(j) - fxij
            fy(j) = fy(j) - fyij
            fz(j) = fz(j) - fzij

  110       continue

c           calculate non-Coulomb potential and forces for dists < rHHnC
            if (rij .le. rHHnC) then

c             help parameters
              hlpHH4 = rij - parHH(4)
              hlpHH7 = rij - parHH(7)
              hexHH1 = 1.0D0 + exp (parHH(3) * hlpHH4)
              hexHH2 = exp (-parHH(6) * hlpHH7**2)

              VnonC = parHH(2) / hexHH1 - parHH(5) * hexHH2
              Force = parHH(2) * parHH(3) * exp (parHH(3) * hlpHH4) /
     $            hexHH1**2 - parHH(5) * parHH(6) * 2.0D0 * hlpHH7 *
     $            hexHH2

              VnCsum = VnCsum +
     $            (VnonC - vcshhh - dcshhh * (rij - rHHnC) )
              fijr = (Force + dcshhh) / rij
              virial = virial + fijr * rsqij

              fxij = fijr * rxij
              fyij = fijr * ryij
              fzij = fijr * rzij

              fx(i) = fx(i) + fxij
              fy(i) = fy(i) + fyij
              fz(i) = fz(i) + fzij

              fx(j) = fx(j) - fxij
              fy(j) = fy(j) - fyij
              fz(j) = fz(j) - fzij
            endif
          endif
        endif

  100   continue
      enddo

      ECou = ECou + VCsum
      EnonC = EnonC + VnCsum
      EpotHH = VCsum + VnCsum

      end
