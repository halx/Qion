      subroutine potfor
c $Id: potfor.f,v 1.8 2005/06/21 11:48:13 hal Exp $
c
c     potfor: sets energies and forces to zero and calls all the power
c             subroutines
c             pfwtxx  - calculate the shortrange forces and energy of water
c             pfintr  - takes care of the intra-molecular energy of H2O
c             pf2bd   - calculate ion-O and ion-H pair potential
c             pf3bdc  - calculate O-ion-O three body correction
c             pfrf    - calculate the reaction field approximation
c             qmmm    - calculates QM/MM forces and energies
c             umbpot  - umbrella biasing potential
c

      implicit none

      integer i

      include 'params.inc'
      include 'potfor.inc'
      include 'rfpot.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'neighb.inc'
      include 'energy.inc'
      include 'press.inc'



c     set forces to zero
      do i = 1, nwc
        fx(i) = 0.0D0
        fy(i) = 0.0D0
        fz(i) = 0.0D0
      enddo

c     reaction-field energy
      EnrRF = 0.0D0

c     Coulomb and non-Coulomb energy
      ECou = 0.0D0
      EnonC = 0.0D0

c     potential energies
      EpotOO = 0.0D0
      EpotOH = 0.0D0
      Eintr = 0.0D0
      EpotHH = 0.0D0
      EpotiO = 0.0D0
      EpotiH = 0.0D0

c     energy from external force program
      EQM = 0.0D0

c     three-body energy
      V3bd = 0.0D0

c     virial
      virial = 0.0D0

c     call various potential subroutines
      call pfwtoo
      call pfwtoh
      call pfwthh

      if (ion) then
        call pf2bd(nbriO, noniO, pionO, ppionO, qionO,
     $      vccbco, vcshco, dccbco, dcshco, ECou, EnonC, EpotiO)
        call pf2bd(nbriH, noniH, pionH, ppionH, qionH,
     $      vccbch, vcshch, dccbch, dcshch, ECou, EnonC, EpotiH)
      endif

      if (intra) call pfintr

      if (irf) then
        call pfrf(nbrOO, nonOO, cRFOO, VRFcOO, dRFcOO)
        call pfrf(nbrOH, nonOH, cRFOH, VRFcOH, dRFcOH)
        call pfrf(nbrHH, nonHH, cRFHH, VRFcHH, dRFcHH)
        call pfrf(nbriO, noniO, cRFiO, VRFciO, dRFciO)
        call pfrf(nbriH, noniH, cRFiH, VRFciH, dRFciH)
      endif

      if (ithree) call pf3bdc(Idx3bd, NoN3bd, V3bd)

      if (iqmmm) call qmmm

      if (iumb) call umbpot(Ebias)

      end
