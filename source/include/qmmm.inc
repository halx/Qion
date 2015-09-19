c $Id: qmmm.inc,v 1.2 2005/05/26 11:28:42 hal Exp $
c
c     QM/MM lists
c

      logical inqm
      double precision rhlist

      common /qmmmpa/ rhlist(mxno), inqm(mxatom)

      save /qmmmpa/



c     local variables:
c     mode: fortran
c     end:
