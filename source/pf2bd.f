      subroutine pf2bd (IdxLst, NoNbrs, param, pparam,
     $    q1q2, VCcut, VnCcut, dCcut, dnCcut, Eqq, EvdW, Epot)
c $Id: pf2bd.f,v 1.10 2005/06/21 11:48:13 hal Exp $
c
c     pf2bd: calculate pair potential energies and forces
c
c     The general formula for the pair potential is:
c
c     q_1*q_2    A      B      C      D
c     ------- + ---- + ---- + ---- + ---- + E*exp(-F*r)
c        r      r**a   r**b   r**c   r**d
c
c     A--E are energy parameters, F is the exponential parameter,
c     a--d are integer powers
c
c     The general formula of the force is:
c
c     q_1*q_2      A          B          C          D
c     ------- + -------- + -------- + -------- + -------- + F*E*exp(-F*r)
c      r**2     r**(a+1)   r**(b+1)   r**(c+1)   r**(a+1)

      implicit none

      integer IdxLst(*), NoNbrs, pparam(*)
      double precision q1q2, param(*)
      double precision VCcut, VnCcut, dCcut, dnCcut
      double precision Eqq, EvdW, Epot

      integer j, k
      double precision xi, yi, zi, rxij, ryij, rzij, rij, rsqij
      double precision FCoul, FnonC, fijr, fijrC, fxij, fyij, fzij
      double precision VCoul, VnonC, VCsum, VnCsum

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'press.inc'
      include 'qmmm.inc'



      VCsum = 0.0D0
      VnCsum = 0.0D0
      VCoul = 0.0D0

c     coordinates of the ion
      xi = x(nwc)
      yi = y(nwc)
      zi = z(nwc)

      do k = 1, NoNbrs
        j = IdxLst(k)

        if (iqmmm .and. inqm(j)) goto 100

c       calculate the distance components between i and j particle
        rxij = xi - x(j)
        ryij = yi - y(j)
        rzij = zi - z(j)

        call pbcmic(rxij, ryij, rzij)

c       calculate the distance between i and j particle
        rsqij = rxij**2 + ryij**2 + rzij**2

c       calculate potential and forces for distances < rcut
        if (rsqij .le. rsqcut) then
          rij = sqrt(rsqij)

c         Coulomb potential and forces
          if (ipntch .and. .not. inqm(j)) then
            fijrC = 0.0D0
          else
            VCoul = q1q2 / rij
            FCoul = q1q2 / rsqij
c           calculate the shifted-force potential [AT87; pp. 145] and
c           the shifted force
            VCsum = VCsum + (VCoul - VCcut - dCcut * (rij - rcut) )
            fijrC = (FCoul + dCcut) / rij
          endif

c         non-Coulomb potential and forces
          VnonC =
     $        param(1) / rij**pparam(1) +
     $        param(2) / rij**pparam(2) +
     $        param(3) / rij**pparam(3) +
     $        param(4) / rij**pparam(4) +
     $        param(5) * exp (-param(6) * rij)

          FnonC =
     $        pparam(1) * param(1) / rij**(pparam(1) + 1) +
     $        pparam(2) * param(2) / rij**(pparam(2) + 1) +
     $        pparam(3) * param(3) / rij**(pparam(3) + 1) +
     $        pparam(4) * param(4) / rij**(pparam(4) + 1) +
     $        param(6) * param(5) * exp (-param(6) * rij)

c         shifted-force potential
          VnCsum = VnCsum + (VnonC - VnCcut - dnCcut * (rij - rcut) )
          fijr = fijrC + (FnonC + dnCcut) / rij
          virial = virial + fijr * rsqij

c         force components
          fxij = fijr * rxij
          fyij = fijr * ryij
          fzij = fijr * rzij

c         sum force components for current atom i
          fx(nwc) = fx(nwc) + fxij
          fy(nwc) = fy(nwc) + fyij
          fz(nwc) = fz(nwc) + fzij

c         subtract force components from neighbour atom j
          fx(j) = fx(j) - fxij
          fy(j) = fy(j) - fyij
          fz(j) = fz(j) - fzij
        endif

  100   continue
      enddo

      Eqq = Eqq + VCsum
      EvdW = EvdW + VnCsum
      Epot = VCsum + VnCsum

      end
