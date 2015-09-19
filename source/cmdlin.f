      subroutine cmdlin
c $Id: cmdlin.f,v 1.1.1.1 2003/06/17 10:15:22 hal Exp $
c
c     cmdlin: parse command line
c

      implicit none

      intrinsic IArgC
      external lentrm
      integer IArgC, lentrm

      integer i, arg, nargs, start, length
      character*128 argv, prgnam

      include 'units.inc'
      include 'files.inc'



      start = 1

c     get number of arguments
      nargs = IArgC ()

      if (nargs .gt. 0) then
        arg = 0

  100   continue

        arg = arg + 1

        call GetArg (arg, argv)

        if (argv(1:2) .eq. '-O') then
          mode = 'unknown'
        else if (argv(1:2) .eq. '-w') then
          arg = arg + 1
          call GetArg (arg, fwater)
        else if (argv(1:2) .eq. '-h') then
          call GetArg (0, prgnam)

c         find last occurence of '/'
          do i = len (prgnam), 1, -1
            if (prgnam(i:i) .eq. '/') then
              start = i + 1
              goto 110
            endif
          enddo

  110     continue

          length = lentrm (prgnam)

          write (*,*) prgnam(start:length), 
     $        ' [-w file] [-O] [-h] [<] input_file'
          write (*,*) 'perform QM/MM-MD calculations'
          write (*,*)
          write (*,*) '-w file   user supplied CF2 water parameter file'
          write (*,*) '-O        overwrite mode'
          write (*,*) '-h        this help'
          write (*,*)

          call exit (0)
        else
          fin = argv
        endif

        if (arg .lt. nargs) goto 100


        if (fin .ne. '*stdin*') then
          length = lentrm (fin)
          open (uin, err = 500, file = fin, status = 'old')
        endif
      endif

      return


  500 write (stderr, *) 'Could not open/read from: ', fin(1:length)
      call exit (1)

      end
