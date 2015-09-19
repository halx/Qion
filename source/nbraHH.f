      subroutine nbraHH
c $Id: nbraHH.f,v 1.3 2005/06/02 10:20:03 hal Exp $
c
c     nebrhh: create atom based H-H neighbour list
c

      implicit none

      integer i, j, idx, lower, upper
      double precision xi, yi, zi, rxij, ryij, rzij, rsqij, rsqn

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'



c     set the index of the O-O neighbour list to zero
      nonHH = 0

c     square distance of neighbour list radius
      rsqn = (rcut + rplus)**2

      if (intra) then
        upper = nw - 2
      else
        upper = nw - 1
      endif

c     find the H neighbours of all hydrogen atoms
      do i = no + 1, upper
        idx = i - no

        xi = x(i)
        yi = y(i)
        zi = z(i)

        nonHH = nonHH + 1
        nbrHH(nonHH) = -i

        if (intra) then
c         make sure that hydrogens on the same water molecule are *not*
c         counted as neighbours
          lower = i + 1 + mod(idx, 2)
        else
          lower = i + 1
        endif

        do j = lower, nw
c         calculate distance vector components for i and j particle
          rxij = xi - x(j)
          ryij = yi - y(j)
          rzij = zi - z(j)

          call pbcmic (rxij, ryij, rzij)

c         calculate square distance vector for i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         find all neighbours within the cutoff radius
          if (rsqij .le. rsqn) then
            nonHH = nonHH + 1
            nbrHH(nonHH) = j
          endif
        enddo
      enddo

      end
