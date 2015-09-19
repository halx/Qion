c $Id: energy.inc,v 1.3 2005/06/21 11:48:15 hal Exp $
c
c     energy variables
c

      integer natqm

      double precision EpotOO, EpotOH, EpotHH, EpotiO, EpotiH, EnonC
      double precision ECou, Enrrf, EQM, Eintr, Ebias, V3bd, V0

      common /energy/ EpotOO, EpotOH, EpotHH, EpotiO, EpotiH, EnonC,
     $    ECou, Enrrf, EQM, Eintr, Ebias, V3bd, V0, natqm

      save /energy/


c     local variables:
c     mode: fortran
c     end:
