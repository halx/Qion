      subroutine pinit
c $Id: pinit.f,v 1.3 2005/06/02 10:29:04 hal Exp $
c
c     pinit: calculate potential and first derivative at the cutoff
c            distance for later computation of shifted-force potential
c            initialize some useful constants
c

      implicit none

      double precision hlpOO6, hlpOO9, hlpOH6, hlpOH9, hexOH1, hexOH2
      double precision hlpHH4, hlpHH7, hexHH1, hexHH2

      include 'params.inc'
      include 'potfor.inc'
      include 'consts.inc'



c     helper parameters
      hlpOO6 = rcut - parOO(6)
      hlpOO9 = rcut - parOO(9)
      hlpOH6 = rOHnC - parOH(6)
      hlpOH9 = rOHnC - parOH(9)
      hexOH1 = 1.0D0 + exp (parOH(5) * hlpOH6)
      hexOH2 = 1.0D0 + exp (parOH(8) * hlpOH9)
      hlpHH4 = rHHnC - parHH(4)
      hlpHH7 = rHHnC - parHH(7)
      hexHH1 = 1.0D0 + exp (parHH(3) * hlpHH4)
      hexHH2 = exp (-parHH(6) * hlpHH7**2)

c     Coulomb term: q_1 * q_2 * conversion factor
      qionO = qO * qion * coucnv
      qionH = qH * qion * coucnv

c     potential at cutoff O-O
      vccboo = parOO(1) / rcut
      vcshoo = parOO(2) / rcut**parOO(3) -
     $    parOO(4) * exp (-parOO(5) * hlpOO6**2) -
     $    parOO(7) * exp (-parOO(8) * hlpOO9**2)

c     derivative at cutoff O-O
      dccboo = -parOO(1) / rsqcut
      dcshoo = -parOO(3) * parOO(2) / rcut**(parOO(3) + 1.0D0) +
     $    parOO(4) * parOO(5) * 2.0D0 * hlpOO6 *
     $    exp (-parOO(5) * hlpOO6**2)+
     $    parOO(7) * parOO(8) * 2.0D0 * hlpOO9 *
     $    exp (-parOO(8) * hlpOO9**2)

c     potential at cutoff O-H
      vccboh = parOH(1) / rcut
      vcshoh = parOH(2) / rOHnC**parOH(3) - parOH(4) / hexOH1 -
     $    parOH(7) / hexOH2

c     derivative at cutoff O-H
      dccboh = -parOH(1) / rsqcut
      dcshoh = -parOH(3) * parOH(2) / rOHnC**(parOH(3) + 1.0D0) +
     $    parOH(4) * parOH(5) * exp (parOH(5) * hlpOH6) / hexOH1**2 +
     $    parOH(7) * parOH(8) * exp (parOH(8) * hlpOH9) / hexOH2**2

c     potential at cutoff H-H
      vccbhh = parHH(1) / rcut
      vcshhh = parHH(2) / hexHH1 - parHH(5) * hexHH2

c     derivative at cutoff H-H
      dccbhh = -parHH(1) / rsqcut
      dcshhh = -parHH(2) * parHH(3) *
     $    exp (parHH(3) * hlpHH4) / hexHH1**2 +
     $    parHH(5) * parHH(6) * 2.0D0 * hlpHH7 * hexHH2

c     potential at cutoff ion-O
      vccbco = qionO / rcut
      vcshco = pionO(1) / rcut**ppionO(1) +
     $    pionO(2) / rcut**ppionO(2) +
     $    pionO(3) / rcut**ppionO(3) +
     $    pionO(4) / rcut**ppionO(4) +
     $    pionO(5) * exp (-pionO(6) * rcut)

c     derivative at cutoff ion-O
      dccbco = -qionO / rsqcut
      dcshco = -ppionO(1) * pionO(1) / rcut**(ppionO(1) + 1) -
     $    ppionO(2) * pionO(2) / rcut**(ppionO(2) + 1) -
     $    ppionO(3) * pionO(3) / rcut**(ppionO(3) + 1) -
     $    ppionO(4) * pionO(4) / rcut**(ppionO(4) + 1) -
     $    pionO(6) * pionO(5) * exp (-pionO(6) * rcut)

c     potential at cutoff ion-H
      vccbch = qionH / rcut
      vcshch = pionH(1) / rcut**ppionH(1) +
     $    pionH(2) / rcut**ppionH(2) +
     $    pionH(3) / rcut**ppionH(3) +
     $    pionH(4) / rcut**ppionH(4) +
     $    pionH(5) * exp (-pionH(6) * rcut)

c     derivative at cutoff ion-H
      dccbch = -qionH / rsqcut
      dcshch = -ppionH(1) * pionH(1) / rcut**(ppionH(1) + 1) -
     $    ppionH(2) * pionH(2) / rcut**(ppionH(2) + 1) -
     $    ppionH(3) * pionH(3) / rcut**(ppionH(3) + 1) -
     $    ppionH(4) * pionH(4) / rcut**(ppionH(4) + 1) -
     $    pionH(6) * pionH(5) * exp (-pionH(6) * rcut)

      end
