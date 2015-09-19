      subroutine qmmm
c $Id: qmmm.f,v 1.3 2005/05/26 11:28:42 hal Exp $
c
c     qmmm: perform QM calculations by passing coordinates and
c           point charges to external program and read computed
c           forces
c

      implicit none

      external lentrm
      integer lentrm

      character hnchng

      integer i, j, k, start
      integer cntO, ocntO, cntpnt
      integer length, status, errno
      double precision smfac
      double precision xc, yc, zc, rxijO, ryijO, rzijO, rijO
      double precision rxijH, ryijH, rzijH
      double precision fxabin, fyabin, fzabin

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'potfor.inc'
      include 'intra.inc'
      include 'qmmm.inc'
      include 'energy.inc'
      include 'atoms.inc'
      include 'neighb.inc'
      include 'consts.inc'
      include 'units.inc'
      include 'files.inc'

c     store indices
      integer pntlst(mxno), listO(mxOqm)

      save ocntO
      data ocntO /0/



c     initialization
      cntO = 0
      cntpnt = 0
      hnchng = '0'

c
c     create coordinates for external force program from QM zone
c

      open (ucoord, err = 500, file = fcoord, status = 'unknown')
      rewind (ucoord)

      if (ipntch) then
        open (upntch, err = 510, file = fpntch, status = 'unknown')
        rewind (upntch)
      endif

      if (ion) then
        start = 1

c       ion is central atom
        xc = x(nwc)
        yc = y(nwc)
        zc = z(nwc)

        write (ucoord, 8000) ionnam, 0.0, 0.0, 0.0
      else
        start = 2

c       first oxygen is central atom
        xc = x(1)
        yc = y(1)
        zc = z(1)

        write (ucoord, 8000) 'O ', 0.0, 0.0, 0.0

        cntO = cntO + 1
        listO(cntO) = 1

        do j = no + 1, no + 2
          rxijH = x(j) - xc
          ryijH = y(j) - yc
          rzijH = z(j) - zc

          call pbcmic(rxijH, ryijH, rzijH)

          write (ucoord, 8000) 'H ', rxijH, ryijH, rzijH
        enddo
      endif

      do i = start, no
c       distance O-central atom
        rxijO = x(i) - xc
        ryijO = y(i) - yc
        rzijO = z(i) - zc

        call pbcmic (rxijO, ryijO, rzijO)

        if (inqm(i)) then
          write (ucoord, 8000) 'O ', rxijO, ryijO, rzijO

          cntO = cntO + 1
          listO(cntO) = i

          do j = no + 2 * i - 1, no + 2 * i
            rxijH = x(j) - xc
            ryijH = y(j) - yc
            rzijH = z(j) - zc

            call pbcmic(rxijH, ryijH, rzijH)

            write (ucoord, 8000) 'H ', rxijH, ryijH, rzijH
          enddo
        else if (ipntch) then
          write (upntch, 8010) rxijO, ryijO, rzijO, qO

          cntpnt = cntpnt + 1
          pntlst(cntpnt) = i

          do j = no + 2 * i - 1, no + 2 * i
            rxijH = x(j) - xc
            ryijH = y(j) - yc
            rzijH = z(j) - zc

            call pbcmic(rxijH, ryijH, rzijH)

            write (upntch, 8010) rxijH, ryijH, rzijH, qH
          enddo
        endif
      enddo

      close (ucoord)
      if (ipntch) close (upntch)

      natqm = cntO


c
c     call external force program
c

c     check if number of oxygens has changed in QM zone
      if (cntO .ne. ocntO) hnchng = '1'
      ocntO = cntO

c     check if user asks external force program to *not* perform the
c     initial guess
      if (iguess) then
        hnchng = '0'
        iguess = .false.
      endif

c     arguments to external force program are:
c
c     arg #1: file containing atomnames+coordinates
c     arg #2: file into which external force program writes forces
c     arg #3: file containing point charges+coordinates
c     arg #4: flag indicating changes in QM zone

      call fsys(frcpth // ' ' // fcoord // ' ' // fforce // ' ' //
     $    fpntch // ' ' // hnchng, status, errno)

      if (status .ne. 0 .or. errno .ne. 0) then
        length = lentrm (frcpth)

        write (stderr, *) 'An error occured while executing/running ',
     $      frcpth(1:length), '.  Aborting!'

        write (uout, *)

        write (uout, *) 'An error occured while executing/running ',
     $      frcpth(1:length), '.  Aborting!'
        write (uout, *)

        call finish(3)
      endif


c
c     read energy and forces calculated by external force program and apply
c     smoothing function to forces between ron and roff, otherwise use
c     QM forces (r <= ron) or MD forces (r > roff)
c

      open (uforce, err = 520, file = fforce, status = 'old')
      rewind (uforce)

c     energy in hartree
      read (uforce, *) EQM
      EQM = EQM * hart2k

c     read force components of ion (always use full QM forces)
      if (ion) then
        read (uforce, *) fxabin, fyabin, fzabin

        fx(nwc) = fx(nwc) + fxabin * frccnv
        fy(nwc) = fy(nwc) + fyabin * frccnv
        fz(nwc) = fz(nwc) + fzabin * frccnv
      endif

      do i = 1, cntO
        j = listO(i)
        rijO = sqrt(rhlist(j))

c       smoothing factor is 1.0 (= full QM forces) for particles
c       inside ron, between ron/roff the switching function is used
        if (rijO .le. ron) then
          smfac = 1.0D0
        else
          smfac = (roff**2 - rijO**2)**2 *
     $        (roff**2 + 2.0D0 * rijO**2 - 3.0D0 * ron**2) /
     $        (roff**2 - ron**2)**3
        endif

c       read O force components and smooth forces if necessary
        read (uforce, *) fxabin, fyabin, fzabin

        fx(j) = fx(j) + smfac * (fxabin * frccnv)
        fy(j) = fy(j) + smfac * (fyabin * frccnv)
        fz(j) = fz(j) + smfac * (fzabin * frccnv)

c       read H force components and smooth forces if necessary
        do k = no + 2 * j - 1, no + 2 * j
          read (uforce, *) fxabin, fyabin, fzabin

          fx(k) = fx(k) + smfac * (fxabin * frccnv)
          fy(k) = fy(k) + smfac * (fyabin * frccnv)
          fz(k) = fz(k) + smfac * (fzabin * frccnv)
        enddo
      enddo

      if (ipntch) then
        do i = 1, cntpnt
          j = pntlst(i)

c         read O force components
          read (uforce, *) fxabin, fyabin, fzabin

          fx(j) = fx(j) + fxabin * frccnv
          fy(j) = fy(j) + fyabin * frccnv
          fz(j) = fz(j) + fzabin * frccnv

c         read H force components
          do k = no + 2 * j - 1, no + 2 * j
            read (uforce, *) fxabin, fyabin, fzabin

            fx(k) = fx(k) + fxabin * frccnv
            fy(k) = fy(k) + fyabin * frccnv
            fz(k) = fz(k) + fzabin * frccnv
          enddo
        enddo
      endif

      close (uforce)

      return


  500 write (stderr, *) 'Could not open: ', fcoord
      call exit (1)
  510 write (stderr, *) 'Could not open: ', fpntch
      call exit (1)
  520 write (stderr, *) 'Could not open: ', fforce
      call exit (1)


 8000 format (1x, a2, 3(2x, f14.8) )
 8010 format (1x, 3(2x, f14.8), 2x, f8.5)

      end
