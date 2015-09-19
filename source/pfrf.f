      subroutine pfrf (IdxLst, NoNbrs, RFcon, VRFcut, dRFcut)
c $Id: pfrf.f,v 1.6 2005/06/10 10:47:12 hal Exp $
c
c     pfrf: reaction-field long-range correction
c           (originally written by Teerakiat Kerdcharoen)

      implicit none

      integer IdxLst(*), NoNbrs
      double precision RFcon, VRFcut, dRFcut

      integer i, j, k
      double precision xi, yi, zi, rxij, ryij, rzij, rij, rsqij
      double precision fxij, fyij, fzij, fijr, Force, VRF

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'energy.inc'
      include 'press.inc'



      do k = 1, NoNbrs
        if (IdxLst(k) .lt. 0) then
          i = -IdxLst(k)

          xi = x(i)
          yi = y(i)
          zi = z(i)
        else
          j = IdxLst(k)

c         calculate the distance components between i and j particle
          rxij = xi - x(j)
          ryij = yi - y(j)
          rzij = zi - z(j)

          call pbcmic(rxij, ryij, rzij)

c         calculate the distance between i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         calculate potential and forces for distances less than rcut
          if (rsqij .le. rsqcut) then
            rij = sqrt(rsqij)

            VRF = RFcon * rsqij / rcut3
            Force = -2.0D0 * RFcon * rij / rcut3

c           calculate the shifted-force potential [AT87; pp. 145] and
c           the shifted force
            EnrRF = EnrRF + (VRF - VRFcut - dRFcut * (rij - rcut) )
            fijr = (Force + dRFcut) / rij
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
      enddo

      end
