      subroutine nbraOH
c $Id: nbraOH.f,v 1.4 2005/06/02 10:20:03 hal Exp $
c
c     nebroh: create atom based O-H neighbour list
c

      implicit none

      integer i, j
      double precision xi, yi, zi, rxij, ryij, rzij, rsqij, rsqn


      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'



c     set the index of the O-H neighbour lists to zero
      nonOH = 0

c     square distances of neighbour list radius
      rsqn = (rcut + rplus)**2

c     find the H neighbours of all oxygen atoms
      do i = 1, no
        xi = x(i)
        yi = y(i)
        zi = z(i)

        nonOH = nonOH + 1
        nbrOH(nonOH) = -i

        do j = no + 1, nw
c         make sure hydrogens bound to current oxygen i are not counted
          if (intra .and.
     $        (j .eq. (no + 2*i - 1) .or. j .eq. (no + 2*i))) goto 100

c         calculate distance vector components for i and j particle
          rxij = xi - x(j)
          ryij = yi - y(j)
          rzij = zi - z(j)

          call pbcmic (rxij, ryij, rzij)

c         calculate square distance vector for i and j particle
          rsqij = rxij**2 + ryij**2 + rzij**2

c         find all neighbours within the cutoff radius
          if (rsqij .le. rsqn) then
            nonOH = nonOH + 1
            nbrOH(nonOH) = j
          endif

  100     continue
        enddo
      enddo

      end
