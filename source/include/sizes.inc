c $Id: sizes.inc,v 1.3 2005/06/03 11:26:29 hal Exp $
c
c     number of atoms determining static sizes of arrays
c

c     maximum atom numbers
      integer mxatom, mxno, mxnh
      parameter (mxatom = 3001, mxno = 1000, mxnh = 2000)

c     neighbour list sizes for ~3000 atoms and a cutoff radius of ~15.0A
      integer mxoo, mxoh, mxhh, mxiono, mxionh
      parameter (mxoo = 300000, mxoh = 1200000, mxhh = 1200000)
      parameter (mxiono = 1000, mxionh = 2000)

c     maximum numbers of oxygens and hydrogens in the QM zone
      integer mxOqm, mxHqm
      parameter (mxOqm = 50, mxHqm = 100)

c     current number of atoms
      integer no, nh, nc, nw, nwc
      common /sizes/ no, nh, nc, nw, nwc


      save /sizes/



c     local variables:
c     mode: fortran
c     end:
