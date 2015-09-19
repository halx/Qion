      subroutine nbrion
c $Id: nbrion.f,v 1.5 2005/06/02 10:20:03 hal Exp $
c
c     nbrion: create ion-O and ion-H neighbour lists
c

      implicit none

      integer i
      double precision xi, yi, zi, rxij, ryij, rzij
      double precision rsqn, rsqij, rsq3bd

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'



c     set neighbour list pointers to zero
      noniO = 0
      noniH = 0
      NoN3bd = 0

c     suare distance of neighbour list radius
      rsqn = (rcut + rplus)**2
      rsq3bd = r3bd**2

c     coordinates of the ion
      xi = x(nwc)
      yi = y(nwc)
      zi = z(nwc)

      do i = 1, no
c       calculate distance vector components for i and j particle
        rxij = xi - x(i)
        ryij = yi - y(i)
        rzij = zi - z(i)

        call pbcmic(rxij, ryij, rzij)

c       calculate square distance vector for i and j particle
        rsqij = rxij**2 + ryij**2 + rzij**2

c       find all neighbours within rsqn
        if (rsqij .le. rsqn) then
          noniO = noniO + 1
          nbriO(noniO) = i

          noniH = noniH + 2
          nbriH(noniH-1) = no + 2 * i - 1
          nbriH(noniH) = no + 2 * i

c         create neighbour list for 3-body corrections
          if (ithree .and. rsqij .le. rsq3bd) then
            NoN3bd = NoN3bd + 1
            Idx3bd(NoN3bd) = i
          endif
        endif
      enddo

      end
