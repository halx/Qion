      subroutine prslin (unit, buff, status)
c $Id: prslin.f,v 1.1.1.1 2003/06/17 10:15:23 hal Exp $
c
c     prslin: parse `buff' and discard comment lines
c
c     in:  unit - I/O unit number
c     out: buff - buffer to which contents of file is written
c          status - error flag
c

      implicit none

      external lentrm
      integer lentrm

      integer unit, status
      character*(*) buff

      integer i, j, start, length
      character cmtchr*5, line*80

      data cmtchr /'#!*%'/



      status = 0
      start = 1

c     preparse file and read complete file into temporary buffer
  100 continue
      read (unit, '(A80)', end = 120, err = 130) line

      length = lentrm (line)

c     ignore empty lines and comment lines
      if (line .eq. ' ') goto 100

      do i = 1, length
        if (line(i:i) .ne. ' ') then
          do j = 1, len (cmtchr)
            if (line(i:i) .eq. cmtchr(j:j)) goto 100
          enddo

          goto 110
        endif
      enddo

  110 continue

c     fill buffer with contents of file
      buff(start:start+length) = line(1:length)
      start = start + length + 1

      goto 100

  120 continue

      return


  130 continue
      status = 1

      end
