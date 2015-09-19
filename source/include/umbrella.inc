c $Id: umbrella.inc,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
c
c     data structures for umbrella sampling
c

      integer mxumb, maxust
      parameter (mxumb = 100, maxust = 1000)

      integer winlen, wincnt

      double precision uhisto, kumb, d0umb

      common /umbprm/ uhisto(maxust), kumb(mxumb), d0umb(mxumb),
     $    winlen, wincnt

      save /umbprm/


c     local variables:
c     mode: fortran
c     end:
