      subroutine rfinit
c $Id: rfinit.f,v 1.2 2005/06/02 10:20:03 hal Exp $
c
c     rfinit: calculate potential and first derivative at the cutoff
c             distance for later computation of shifted-force potential
c             (reaction field)
c

      implicit none

      double precision eH2O, term

      include 'params.inc'
      include 'potfor.inc'
      include 'rfpot.inc'
      include 'consts.inc'



c     dielectric constant (permittivity) of water at chosen temperature
c     e_H2O(T) = a + bT + cT**2, fitted within 273.3--372.2 K
      eH2O = 0.24921D3 - 0.79069D0 * Topt + 0.72997D-03 * Topt**2

      term = (eH2O - 1.0D0) * coucnv / (2.0D0 * eH2O + 1.0D0)

c     constants
      cRFOO = qO * qO * term
      cRFOH = qO * qH * term
      cRFHH = qH * qH * term
      cRFiO = qion * qO * term
      cRFiH = qion * qH * term

c     RF potential at O-O cutoff
      VRFcOO = cRFOO / rcut

c     derivative of RF potential at O-O cutoff
      dRFcOO = 2.0D0 * cRFOO / rsqcut

      VRFcOH = cRFOH / rcut
      dRFcOH = 2.0D0 * cRFOH / rsqcut

      VRFcHH = cRFHH / rcut
      dRFcHH = 2.0D0 * cRFHH / rsqcut

      VRFciO = cRFiO / rcut
      dRFciO = 2.0D0 * cRFiO / rsqcut

      VRFciH = cRFiH / rcut
      dRFciH = 2.0D0 * cRFiH / rsqcut

      end
