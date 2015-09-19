      subroutine nbratb
c $Id: nbratb.f,v 1.3 2005/05/25 13:07:40 hal Exp $
c
c     nebr: create atom based neighbour lists
c

      implicit none

      include 'params.inc'
      include 'sizes.inc'
      include 'neighb.inc'



      call nbraOO
      call nbraOH
      call nbraHH

      if (ion) then
        call nbrlst(nbriO, noniO,      1, no, rcut)
        call nbrlst(nbriH, noniH, no + 1, nw, rcut)
      endif

      if (ithree) call nbrlst(Idx3bd, NoN3bd, 1, no, r3bd)

      end
