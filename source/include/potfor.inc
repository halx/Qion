c $Id: potfor.inc,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
c
c     potential and force parameters and constants
c

c     maximum number of potential parameters for the CF2 water model
c     (F. H. Stillinger and A. Rahman, JPC 68(1978), 666--670)
      integer nparOO, nparOH, nparHH
      parameter (nparOO = 9, nparOH = 9, nparHH = 7)

c     maximum number of potential parameters for ion-O and ion-H
      integer mxpiO, mxpiH, mxppiO, mxppiH
      parameter (mxpiO = 6, mxpiH = 6, mxppiO = 4, mxppiH = 4)


      double precision vcshoo, vccboo, dcshoo, dccboo
      double precision vcshoh, vccboh, dcshoh, dccboh
      double precision vcshhh, vccbhh, dcshhh, dccbhh
      double precision vcshco, vccbco, dcshco, dccbco
      double precision vcshch, vccbch, dcshch, dccbch

      common /pfcon/
     $    vcshoo, vccboo, dcshoo, dccboo,
     $    vcshoh, vccboh, dcshoh, dccboh,
     $    vcshhh, vccbhh, dcshhh, dccbhh,
     $    vcshco, vccbco, dcshco, dccbco,
     $    vcshch, vccbch, dcshch, dccbch


      double precision qO, qH, qionO, qionH
      common /charge/ qO, qH, qionO, qionH


      integer ppionO, ppionH
      double precision parOO, parOH, parHH, pionO, pionH, rOHnC, rHHnC
      common /pfprm/ parOO(nparOO), parOH(nparOH), parHH(nparHH),
     $    pionO(mxpiO), pionH(mxpiH), ppionO(mxppiO), ppionH(mxppiH)
      common /cutoff/ rOHnC, rHHnC


      save /pfcon/, /charge/, /pfprm/, /cutoff/



c     local variables:
c     mode: fortran
c     end:
