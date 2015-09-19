      double precision function getmas (elem)
c $Id: getmas.f,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
c
c     getmas: return the atom mass pertaining to element `elem',
c             (data taken from http://www.webelements.com/)
c
c     in: element name `elem'
c     out: corresponding mass
c

      implicit none

      character*2 elem

      integer i
      character*2 cmpstr

      integer nentr
      parameter (nentr = 69)

      double precision ramass(nentr)
      character*2 ptab(nentr)

c     (partial) periodic table
      data ptab /'h ', 'd ', 'he', 'li', 'be', 'b', 'c ', 'n ', 'o ',
     $    'f ', 'ne', 'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar',
     $    'k ', 'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni',
     $    'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr',
     $    'y ', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd',
     $    'in', 'sn', 'sb', 'te', 'i ', 'xe', 'cs', 'ba', 'hf', 'ta',
     $    'w ', 're', 'os', 'ir', 'pt', 'au', 'hg', 'tl', 'pb', 'bi'/
      data ramass /1.00794D0, 2.014D0, 4.002602D0, 6.941D0, 9.012182D0,
     $    10.811D0, 12.0107D0, 14.00674D0, 15.9994D0, 18.9984032D0,
     $    20.1797D0, 22.998977D0, 24.305D0, 26.981538D0, 28.0855D0,
     $    30.973761D0, 32.066D0, 35.4527D0, 39.948D0, 39.0983D0,
     $    40.078D0, 44.95591D0, 47.867D0, 50.9415D0, 51.9961D0,
     $    54.93804D0, 55.845D0, 58.9332D0, 58.6934D0, 63.546D0, 65.39D0,
     $    69.723D0, 72.61D0, 74.9216D0, 78.96D0, 79.904D0, 83.8D0,
     $    85.4678D0, 87.62D0, 88.90585D0, 91.224D0, 92.90638D0, 95.94D0,
     $    98.0D0, 101.07D0, 102.9055D0, 106.42D0, 107.8682D0, 112.411D0,
     $    114.818D0, 118.71D0, 121.76D0, 127.6D0, 126.90447D0, 131.29D0,
     $    132.90545D0, 137.327D0, 178.49D0, 180.9479D0, 183.83D0,
     $    186.207D0, 190.23D0, 192.217D0, 195.078D0, 196.96655D0,
     $    200.59D0, 204.3833D0, 207.2D0, 208.98038D0/



      getmas = 0.0D0

c     copy element name to temporary string since tcase modifies its
c     first argument
      cmpstr = elem
      call tcase (cmpstr, 'l')

      do i = 1, nentr
        if (cmpstr .eq. ptab(i)) then
          getmas = ramass(i)
          return
        endif
      enddo

      end
