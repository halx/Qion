c $Id: intra.inc,v 1.3 2005/06/03 11:26:29 hal Exp $
c
c     intra-molecular water potential
c     (P. Bopp et. al., Chem. Phys. Letters 98(1983), 129--133)
c

c     maximum number of intra-molecular potential parameters
      integer npotp
      parameter (npotp = 18)

      integer irho1, irho2, ialpha
      double precision ppintr

      common /pint1/ ppintr(npotp), irho1(npotp), irho2(npotp),
     $    ialpha(npotp)


      double precision doh0, alpha0
      common /pint2/ doh0, alpha0


      save /pint1/, /pint2/



c     local variables:
c     mode: fortran
c     end:
