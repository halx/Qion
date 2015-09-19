      subroutine qmlist
c $Id: qmlist.f,v 1.2 2005/06/02 10:20:03 hal Exp $
c
c     qmlist: create a list of oxygen atoms within the QM zone
c

      implicit none

      integer i, start
      double precision xc, yc, zc, rxij, ryij, rzij, rsqn, rsqij

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'qmmm.inc'



      if (ion) then
        start = 1

        xc = x(nwc)
        yc = y(nwc)
        zc = z(nwc)
      else
        start = 2

        xc = x(1)
        yc = y(1)
        zc = z(1)
      endif

      rsqn = roff**2

      do i = start, no
        rxij = xc - x(i)
        ryij = yc - y(i)
        rzij = zc - z(i)

        call pbcmic (rxij, ryij, rzij)

        rsqij = rxij**2 + ryij**2 + rzij**2
        rhlist(i) = rsqij

        if (rsqij .le. rsqn) then
          inqm(i) = .true.
          inqm(no+2*i-1) = .true.
          inqm(no+2*i) = .true.
        else
          inqm(i) = .false.
          inqm(no+2*i-1) = .false.
          inqm(no+2*i) = .false.
        endif
      enddo

      end
