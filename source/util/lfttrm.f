      integer function lfttrm (string)
c $Id: lfttrm.f,v 1.1.1.1 2003/06/17 10:15:23 hal Exp $
c
c     lfttrm: calculate the position of the first non-blank in a string
c
c     in:  string of which first non-blank should be determined
c     out: position of the first non-blank in the string

      implicit none

      integer i
      character*(*) string



      lfttrm = 0

      do i = 1, len (string)
        if (string(i:i) .ne. ' ') then
          lfttrm = i
          return
        endif
      enddo

      end
