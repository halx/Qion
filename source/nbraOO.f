      subroutine nbraOO
c $Id: nbraOO.f,v 1.3 2005/06/02 10:20:03 hal Exp $
c
c     nebroo: create atom based O-O neighbour list
c

      implicit none

      integer i, j
      double precision xi, yi, zi, rxij, ryij, rzij, rsqij, rsqn

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'



c     reset the neighbour list counter
      nonOO = 0

c     square distance of neighbour list radius
      rsqn = (rcut + rplus)**2

c     find the O neighbours of all oxygen atoms
      do i = 1, no - 1
        xi = x(i)
        yi = y(i)
        zi = z(i)

        nonOO = nonOO + 1
        nbrOO(nonOO) = -i

        do j = i + 1, no
c         calculate distance vector components for i and j particle
          rxij = xi - x(j)
          ryij = yi - y(j)
          rzij = zi - z(j)

          call pbcmic (rxij, ryij, rzij)

c         calculate square distance vector for i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         find all neighbours within the cutoff radius
          if (rsqij .le. rsqn) then
            nonOO = nonOO + 1
            nbrOO(nonOO) = j
          endif
        enddo
      enddo

      end
