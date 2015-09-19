      subroutine finish (error)
c $Id: finish.f,v 1.3 2003/06/18 09:25:21 hal Exp $
c
c     finish:  write final information
c

      implicit none

      integer error

      intrinsic hostnm, etime
      external lentrm
      integer lentrm, hostnm
      real etime

      real runtim, dummy(2)
      integer mins, hours, days, hnerr, length
      character date*24, hname*60

      include 'params.inc'
      include 'units.inc'
      include 'files.inc'



c     runtime in seconds
      runtim = etime (dummy)

      days = int (runtim / 86400.0)
      runtim = runtim - (days * 86400)

      hours = int (runtim / 3600.0)
      runtim = runtim - (hours * 3600)

      mins = int (runtim / 60.0)
      runtim = runtim - (mins * 60)


      write (uout, *)

      write (uout, '(1x, 2a)')  '-------------------------------------',
     $    '-----------------------------------------'

      write (uout, '(/)')

      write (uout, '(1x, a, i3, 2(a, i2), a, f4.1, a)')
     $    'CPU cycles wasted: ', days, ' days ', hours, ' hours ', mins,
     $    ' minutes ', runtim, ' seconds'
      write (uout, *)

c     current date and hostname
      call fdate (date)
      hnerr = hostnm (hname)

      if (hnerr .eq. 0) then
        length = lentrm (hname)

        if (error .ne. 0) then
          write (uout, '(1x, 4a)') 'ABEND ', date, ' on ',
     $        hname(1:length)
        else
          write (uout, '(1x, 4a)') 'program terminated normally ', date,
     $        ' on ', hname(1:length)
        endif
      else
        if (error .ne. 0) then
          write (uout, *) 'ABEND'
        else
          write (uout, *) 'program terminated normally'
        endif
      endif


      close (uout)
      close (uenout)

      if (ftraj .ne. ' ') close (utraj)
      if (fveloc .ne. ' ') close (uveloc)
      if (iumb) then
        close (uumbst)
        if (fumbtr .ne. ' ') close (uumbtr)
      endif

      call exit (error)

      end
