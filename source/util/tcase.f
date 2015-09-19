      subroutine tcase (string, tchar)
c $Id: tcase.f,v 1.1.1.1 2003/06/17 10:15:23 hal Exp $
c
c     tcase: transform case of `string' according to `tchar'
c
c     inout: string to be transform,
c     in:    tchar 'u' for uppercase, tchar 'l' for lowercase
c

      implicit none

      character*(*) string
      character tchar

      integer i, j

      character*26 lower, upper
      data lower /'abcdefghijklmnopqrstuvwxyz'/
      data upper /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/


      
      if (tchar .eq. 'l') then
        do i = 1, len (string)
          j = index (upper, string(i:i))
          if (j .gt. 0) string(i:i) = lower(j:j)
        enddo
      elseif (tchar .eq. 'u') then
        do i = 1, len (string)
          j = index (lower, string(i:i))
          if (j .gt. 0) string(i:i) = upper(j:j)
        enddo
      endif

      end
