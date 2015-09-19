c $Id: rfpot.inc,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
c
c     reaction field potentials and forces
c

      double precision VRFcOO, VRFcOH, VRFcHH, VRFciO, VRFciH
      double precision dRFcOO, dRFcOH, dRFcHH, dRFciO, dRFciH
      double precision cRFOO, cRFOH, cRFHH, cRFiO, cRFiH

      common /rfpfOO/ VRFcOO, dRFcOO, cRFOO
      common /rfpfOH/ VRFcOH, dRFcOH, cRFOH
      common /rfpfHH/ VRFcHH, dRFcHH, cRFHH
      common /rfpfCO/ VRFciO, dRFciO, cRFiO
      common /rfpfCH/ VRFciH, dRFciH, cRFiH


      save /rfpfOO/, /rfpfOH/, /rfpfHH/, /rfpfCO/, /rfpfCH/



c     local variables:
c     mode: fortran
c     end:
