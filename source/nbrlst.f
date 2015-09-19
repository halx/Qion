      subroutine nbrlst (IdxLst, NoNbrs, LBnd, UBnd, cutoff)
c $Id: nbrlst.f,v 1.5 2005/06/02 10:20:03 hal Exp $
c
c     nbrlst: generate neighbour list (general routine for ion)
c

      implicit none

      integer IdxLst(*)
      integer NoNbrs, LBnd, UBnd
      double precision cutoff

      integer j
      double precision xi, yi, zi, rxij, ryij, rzij, rsqn, rsqij

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'



c     reset the neighbour list counter
      NoNbrs = 0

c     suare distance of neighbour list radius
      rsqn = (cutoff + rplus)**2

c     coordinates of the ion
      xi = x(nwc)
      yi = y(nwc)
      zi = z(nwc)

      do j = LBnd, UBnd
c       calculate distance vector components for i and j particle
        rxij = xi - x(j)
        ryij = yi - y(j)
        rzij = zi - z(j)

        call pbcmic (rxij, ryij, rzij)

c       calculate square distance vector for i and j particle
        rsqij = rxij**2 + ryij**2 + rzij**2

c       find all neighbours within rsqn
        if (rsqij .le. rsqn) then
          NoNbrs = NoNbrs + 1
          IdxLst(NoNbrs) = j
        endif
      enddo

      end
