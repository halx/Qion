      integer function lentrm (string)
c $Id: lentrm.f,v 1.1.1.1 2003/06/17 10:15:23 hal Exp $
c
c     lentrm: calculate the length of a string without trailing spaces
c
c     in:  string of which length should be determined
c     out: length of trimmed string, 0 if string contains spaces only
c

      implicit none

      integer i
      character*(*) string



      lentrm = 0

      do i = len (string), 1, -1
        if (string(i:i) .ne. ' ') then
          lentrm = i
          return
        endif
      enddo

      end
