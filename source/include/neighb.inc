c $Id: neighb.inc,v 1.5 2005/05/25 13:07:40 hal Exp $
c
c     neighbour list tables
c

      integer nonOO, nonOH, nonHH, noniO, noniH
      integer nbrOO, nbrOH, nbrHH, nbriO, nbriH
      integer Idx3bd, NoN3bd

      common /neigwt/ nbrOO(mxOO), nbrOH(mxOH), nbrHH(mxHH), nonOH,
     $    nonOO, nonHH
      common /neigiO/ nbriO(mxiono), nbriH(mxionh), Idx3bd(mxiono),
     $    noniO, noniH, NoN3bd

      save /neigwt/, /neigiO/


c     local variables:
c     mode: fortran
c     end:
