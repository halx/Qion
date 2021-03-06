      subroutine umbst
c $Id: umbst.f,v 1.13 2005/06/21 11:48:13 hal Exp $
c
c     umbst: calculate statistics generated by umbrella potential
c

      implicit none

      integer st, bin, i
      double precision rxij, ryij, rzij, rij

      include 'params.inc'
      include 'umbrella.inc'
      include 'stpcnt.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'energy.inc'
      include 'files.inc'
      include 'units.inc'



      st = numbeq + (wincnt - 1) * winlen

      if (cstep .ge. st) then
        if (cstep .eq. st) then
          write (uout, *)
          write (uout, '(1x, a, i7)')
     $        'switching to data collection at step ', cstep

          call flush(uout)
        endif

        rxij = x(nuidx1) - x(nuidx2)
        ryij = y(nuidx1) - y(nuidx2)
        rzij = z(nuidx1) - z(nuidx2)

        call pbcmic(rxij, ryij, rzij)

        rij = sqrt(rxij**2 + ryij**2 + rzij**2)

c       collect distances in histogram
        bin = nint(rij / hdel)
        uhisto(bin) = uhisto(bin) + 1

        if (fumbtr .ne. ' ') then
          write (uumbtr, '(1x, i8, 1x, f8.5, 1x, f12.5)') cstep, rij, V0
        endif
      endif

c     write histogram
      if (mod(cstep, winlen) .eq. 0) then
        do i = 1, maxust
          if (uhisto(i) .gt. 0) then
            write (uumbst, '(1x, f7.4, 1x, i6)')
     $          dble(i) * hdel, uhisto(i)
          endif

          uhisto(i) = 0
        enddo

        write (uumbst, *)
        call flush(uumbst)

        if (fumbtr .ne. ' ') then
          write (uumbtr, *)
          call flush(uumbtr)
        endif
      endif

      end
