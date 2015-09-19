      subroutine nbrwt
c $Id: nbrwt.f,v 1.6 2005/06/02 10:20:03 hal Exp $
c
c     nbrwt: create O-O, O-H, and H-H neighbour lists
c

      implicit none

      integer i, j
      double precision xi, yi, zi, rxij, ryij, rzij, rsqij, rsqn

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'



c     reset the neighbour list pointers
      nonOO = 0
      nonOH = 0
      nonHH = 0

c     square distance of neighbour list radius
      rsqn = (rcut + rplus)**2

c     iterate over all O-O atom pairs
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

          call pbcmic(rxij, ryij, rzij)

c         calculate square distance vector for i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         find all neighbours within the cutoff radius
          if (rsqij .le. rsqn) then
            nonOO = nonOO + 1
            nbrOO(nonOO) = j

c           store j hydrogen neighbours of oxygen i
            nonOH = nonOH + 1
            nbrOH(nonOH) = -i
            nonOH = nonOH + 2
            nbrOH(nonOH-1) = no + 2 * j - 1
            nbrOH(nonOH)   = no + 2 * j

c           store i hydrogen neighbours of oxygen j
            nonOH = nonOH + 1
            nbrOH(nonOH) = -j
            nonOH = nonOH + 2
            nbrOH(nonOH-1) = no + 2 * i - 1
            nbrOH(nonOH)   = no + 2 * i

c           store j hydrogen neighbours of first hydrogen i
            nonHH = nonHH + 1
            nbrHH(nonHH) = -(no + 2 * i - 1)
            nonHH = nonHH + 2
            nbrHH(nonHH-1) = no + 2 * j - 1
            nbrHH(nonHH)   = no + 2 * j

c           store i hydrogen neighbours of second hydrogen i
            nonHH = nonHH + 1
            nbrHH(nonHH) = -(no + 2 * i)
            nonHH = nonHH + 2
            nbrHH(nonHH-1) = no + 2 * j - 1
            nbrHH(nonHH)   = no + 2 * j
          endif
        enddo
      enddo

      end
