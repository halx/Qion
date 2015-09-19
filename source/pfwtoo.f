      subroutine pfwtOO
c $Id: pfwtoo.f,v 1.13 2005/06/21 11:48:13 hal Exp $
c
c     pfwtOO: calculate the potential energy and forces between
c             oxygen atoms of CF2 water, see
c             F. H. Stillinger and A. Rahman, JPC 68(1978), 666--670
c

      implicit none

      integer i, j, k
      double precision xi, yi, zi, rxij, ryij, rzij, rsqij, rij
      double precision FCoul, FnonC, fijr, fijrC, fxij, fyij, fzij
      double precision VnonC, VCoul, VnCsum, VCsum, hlpOO6, hlpOO9

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
      VCoul = 0.0D0

      do k = 1, nonOO
        if (nbrOO(k) .lt. 0) then
          i = -nbrOO(k)

          xi = x(i)
          yi = y(i)
          zi = z(i)
        else
          j = nbrOO(k)

c         cycle if both atoms are in the QM zone
          if (iqmmm .and. inqm(i) .and. inqm(j)) goto 100

c         calculate the distance components between i and j particle
          rxij = xi - x(j)
          ryij = yi - y(j)
          rzij = zi - z(j)

          call pbcmic(rxij, ryij, rzij)

c         calculate the distance between i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         calculate potential and forces for distances < rcut
          if (rsqij .le. rsqcut) then
            rij = sqrt (rsqij)

c           helper parameters
            hlpOO6 = rij - parOO(6)
            hlpOO9 = rij - parOO(9)

c           Coulomb potential and forces
            if (ipntch .and. (inqm(i) .neqv. inqm(j)) ) then
              fijrC = 0.0D0
            else
              VCoul = parOO(1) / rij
              FCoul = parOO(1) / rsqij
c             calculate the shifted-force potential [AT87; pp. 145] and
c             the shifted force
              VCsum = VCsum + (VCoul - vccboo - dccboo * (rij - rcut) )
              fijrC = (FCoul + dccboo) / rij
            endif

c           non-Coulom potential and forces
            VnonC =
     $          parOO(2) / rij**parOO(3) -
     $          parOO(4) * exp (-parOO(5) * hlpOO6**2) -
     $          parOO(7) * exp (-parOO(8) * hlpOO9**2)

            FnonC =
     $          parOO(3)*parOO(2) / rij**(parOO(3) + 1.0D0) -
     $          parOO(4) * parOO(5) * 2.0D0 * hlpOO6 *
     $          exp(-parOO(5) * hlpOO6**2) -
     $          parOO(7) * parOO(8) * 2.0D0 * hlpOO9 *
     $          exp(-parOO(8) * hlpOO9**2)

c           calculate the shifted-force potential
            VnCsum = VnCsum + (VnonC - vcshoo - dcshoo * (rij - rcut) )
            fijr = fijrC + (FnonC + dcshoo) / rij
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
          endif
        endif

  100   continue
      enddo

      ECou = ECou + VCsum
      EnonC = EnonC + VnCsum
      EpotOO = VCsum + VnCsum

      end
