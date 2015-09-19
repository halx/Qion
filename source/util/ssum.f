      double precision function ssum (a, l, u)
c $Id: ssum.f,v 1.1.1.1 2003/06/17 10:15:23 hal Exp $
c
c     ssum: calculate the sum of elements of a one dimensional array,
c           starting from a lower bound `l' up to an upper bound `u'
c
c           in: a ... one dimensional array of assumed size
c               l ... lower index
c               u ... upper index
c

      implicit none

      integer l, u
      double precision a(*)

      integer i



      ssum = 0.0D0

      do i = l, u
        ssum = ssum + a(i)
      enddo

      end

