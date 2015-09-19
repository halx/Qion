      subroutine pfwtOH
c $Id: pfwtoh.f,v 1.12 2005/06/21 11:48:13 hal Exp $
c
c     pfwtOH: calculate the potential energy and forces between
c             oxygen and hydrogen atoms of CF2 water, see
c             F. H. Stillinger and A. Rahman, JPC 68(1978), 666--670
c

      implicit none

      integer i, j, k
      double precision xi, yi, zi, rxij, ryij, rzij, rij, rsqij
      double precision Force, fxij, fyij, fzij, fijr
      double precision VnonC, VCoul, VCsum, VnCsum
      double precision hlpOH6, hlpOH9, hexOH1, hexOH2

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

      do k = 1, nonOH
        if (nbrOH(k) .lt. 0) then
          i = -nbrOH(k)

          xi = x(i)
          yi = y(i)
          zi = z(i)
        else
          j = nbrOH(k)

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

            VCoul = parOH(1) / rij
            Force = parOH(1) / rsqij

c           calculate the shifted-force potential [AT87; pp. 145] and
c           the shifted force
            VCsum = VCsum + (VCoul - vccboh - dccboh * (rij - rcut) )
            fijr = (Force + dccboh) / rij
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

c           calculate non-Coulomb potential and forces for dists < rOHnC
            if (rij .le. rOHnC) then
              hlpOH6 = rij - parOH(6)
              hlpOH9 = rij - parOH(9)
              hexOH1 = 1.0D0 + exp (parOH(5) * hlpOH6)
              hexOH2 = 1.0D0 + exp (parOH(8) * hlpOH9)

              VnonC = parOH(2) / rij**parOH(3) - parOH(4) / hexOH1 -
     $            parOH(7) / hexOH2

              Force = parOH(3) * parOH(2) / rij**(parOH(3) + 1.0D0) -
     $            parOH(4) * parOH(5) * exp (parOH(5) * hlpOH6) /
     $            hexOH1**2 -
     $            parOH(7) * parOH(8) * exp (parOH(8) * hlpOH9) /
     $            hexOH2**2

              VnCsum = VnCsum +
     $            (VnonC - vcshoh - dcshoh * (rij - rOHnC) )
              fijr = (Force + dcshoh) / rij
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
      EpotOH = VCsum + VnCsum

      end
