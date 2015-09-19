      subroutine umbpot (Epot)
c $Id: umbpot.f,v 1.6 2005/06/03 11:26:26 hal Exp $
c
c     umbpot: apply biasing energy and forces along reaction coordinate
c

      implicit none

      double precision Epot

      integer st

      double precision rxij, ryij, rzij, rij
      double precision fijr, fxij, fyij, fzij, Force

      include 'params.inc'
      include 'umbrella.inc'
      include 'stpcnt.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'units.inc'
      include 'consts.inc'



      st = mod (cstep, winlen)

      if (st .eq. 0) wincnt = wincnt + 1

c     special case: counting starts at 1
c     special case: last step is divisible by winlen
      if ( (st .eq. 0 .or. cstep .eq. 1) .and. cstep .lt. nstep) then
        write (uout, *)
        write (uout, '(1x, a, i2)') 'starting window no. ', wincnt
        write (uout, '(1x, a, f6.2, a, f5.2)') 'k = ',
     $      kumb(wincnt) * encnv, ', d0 = ', d0umb(wincnt)
        write (uout, '(1x, a, i7)')
     $      'switching to equilibration at step ', cstep

        call flush(uout)
      endif


c     calculate distance along reaction coordinate
      rxij = x(nuidx1) - x(nuidx2)
      ryij = y(nuidx1) - y(nuidx2)
      rzij = z(nuidx1) - z(nuidx2)

      call pbcmic(rxij, ryij, rzij)

      rij = sqrt(rxij**2 + ryij**2 + rzij**2)

c     calculate biasing energy and forces
      Epot  = 0.5D0 * kumb(wincnt) * (rij - d0umb(wincnt) )**2
      Force =        -kumb(wincnt) * (rij - d0umb(wincnt) )

      fijr = Force / rij

      fxij = fijr * rxij
      fyij = fijr * ryij
      fzij = fijr * rzij

      fx(nuidx1) = fx(nuidx1) + fxij
      fy(nuidx1) = fy(nuidx1) + fyij
      fz(nuidx1) = fz(nuidx1) + fzij

      fx(nuidx2) = fx(nuidx2) - fxij
      fy(nuidx2) = fy(nuidx2) - fyij
      fz(nuidx2) = fz(nuidx2) - fzij

      end
