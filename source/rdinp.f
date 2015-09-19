      subroutine rdinp
c $Id: rdinp.f,v 1.24 2005/06/13 10:38:34 hal Exp $
c
c     rdinp: read control parameters and start coordinates, velocities,
c            and forces for the water molecules and the ion from an
c            input file or stdin
c

      implicit none

      external lentrm
      integer lentrm

      integer i, length, ierror

      character qmfpr*5, qmdir*128, buffer*1000

      include 'sizes.inc'
      include 'params.inc'
      include 'potfor.inc'
      include 'qmmm.inc'
      include 'umbrella.inc'
      include 'intra.inc'
      include 'press.inc'
      include 'rvf.inc'
      include 'eqmot.inc'
      include 'atoms.inc'
      include 'stpcnt.inc'
      include 'consts.inc'
      include 'units.inc'
      include 'files.inc'

      namelist /cntrl/ ion, init, intra, irf, ithree, iqmmm, iguess,
     $    iumb, ipntch, iatmbc, nscale, ncom, nstep, nstat, nout, nrst,
     $    noutcv, naver, nwin, numbeq, nuidx1, nuidx2, hdel, Topt, taut,
     $    dt, dens, rcut, ron, roff, r3bd, rplus, amo, amh, amion, qion,
     $    a3bd, b3bd, c3bd, d3bd, e3bd, k3bd, ionnam, fstart,frst, fout,
     $    finfo, fenout, ftraj, fveloc, fftraj, fqmen, fumbst, fumbtr,
     $    iNPT, Popt, Pcomp, iwrap, taup, eqmot, debug

      data buffer /' '/



c     read namelist
      read (uin, cntrl)

c     sanity checks
      if (nstat .lt. 2 .or. nrst .lt. 2 .or. noutcv .lt. 2) then
        write (stderr, *) 'nstat, nrst, and noutcv must be >= 2'
        call exit (2)
      endif

      if (nout .lt. 1) then
        write (stderr, *) 'nout must be >= 1'
        call exit (2)
      endif

      if (iqmmm) then
        if (ron .eq. 0.0D0 .or. roff .eq. 0.0D0) then
          write (stderr, 8000)
     $        'since a QM/MM calculation should be done, the ',
     $        'following parameters must be ',
     $        'supplied explicitly: ron, roff'
          call exit (2)
        endif

c       arbitrary distance but should prevent stupid errors
        if (roff .lt. 2.5D0) then
          write (stderr, *) 'QM/MM radius is too small (< 2.5)'
          call exit (2)
        endif

        if (ron .gt. roff) then
          write (stderr, *)
     $        'setting ron greater than roff doesn''t make sense'
          call exit (2)
        endif
      endif

      if (ipntch .and. .not. iqmmm) then
        write (stderr, *) 'point charges requested but makes only ',
     $      'sense in connection with QM/MM calculation'
        call exit (2)
      endif

      if (.not. iatmbc .and. .not. intra) then
        write (stderr, *) 'molecule based cutoff only supported ',
     $      'for intra molecular potential'
        call exit (2)
      endif

      if (ithree .and.  (a3bd .eq. 0.0D0 .or. b3bd .eq. 0.0D0 .or.
     $                   c3bd .eq. 0.0D0)) then
        write (stderr, 8000)
     $      'since a three-body potential should be included, these ',
     $      'parameters must be ',
     $      'supplied explicitly: a3bd, b3bd, c3bd'
        call exit (2)
      endif

      if (topt .gt. 5000.0D0) then
        write (stderr, *)
     $      'due to numerical reasons the maximum temperature is ',
     $      'restricted to 5000 K'
        call exit (2)
      endif

      if (eqmot .lt. 1 .or. eqmot .gt. 2) then
        write (stderr, *)
     $      'equation of motion must be either 1 or 2'
        call exit (2)
      endif

c     ion-O and ion-H potential parameters
      if (ion) then
        if (ionnam .eq. '  ' .or. qion .eq. 0.0D0 .or. rcut .le. 0.0D0)
     $      then
          write (stderr, 8000)
     $        'at least one of the following parameters is missing, ',
     $        'but must be supplied ',
     $        'explicitly: ionnam, qion, rcut'
          call exit (2)
        endif

        read (uin, *) (pionO(i), i = 1, mxpiO)
        read (uin, *) (ppionO(i), i = 1, mxppiO)
        read (uin, *) (pionH(i), i = 1, mxpiH)
        read (uin, *) (ppionH(i), i = 1, mxppiH)
      endif

c     k, d0 values for umbrella sampling: Ebias = kumb * (d - d0umb)**2
      if (iumb) then
        if (nuidx1 .le. 0 .or. nuidx2 .le. 0) then
          write (stderr, 8000)
     $        'since an umbrella  potential should be included, these ',
     $        'parameters must be ',
     $        'supplied explicitly: nuidx1, nuidx2'
          call exit (3)
        endif

        if (mod (nstep, nwin) .ne. 0) then
          write (stderr, *) 'nstep must be a multiple of nwin'
          call exit (3)
        endif

        winlen = nstep / nwin

        if (numbeq .ge. winlen) then
          write (stderr, *) 'number of data collection steps is zero'
          call exit (3)
        endif

        do i = 1, nwin
          read (uin, *) kumb(i), d0umb(i)
          kumb(i) = kumb(i) / encnv
        enddo

c       zero out histogram
        do i = 1, maxust
          uhisto(i) = 0
        enddo

        wincnt = 1
      endif

      close (uin)


c     read environement variables
      call GetEnv ('QION_FPROG', qmfpr)

      if (qmfpr .eq. ' ') qmfpr = 'tm'

      call GetEnv ('QION_ROOT', qmdir)

      if (qmdir .eq. ' ') then
        write (stderr, *) 'Enviroment variable QION_ROOT set?'
        call exit (2)
      endif

      length = lentrm (qmdir)

c     create complete path for calc_forces
      frcpth = qmdir
      frcpth(length+1:) = '/' // frcexe // '.' // qmfpr

c     create complete path for file containing water model description
      if (fwater .eq. ' ') then
        fwater = qmdir
        fwater(length+1:) = '/' // fwtCF2
      endif

c     read parameters for water model
      length = lentrm (fwater)
      open (uwater, err = 500, file = fwater, status = 'old')
      rewind (uwater)

      read (uwater, '(a)') wattit

      call prslin (uwater, buffer, ierror)

      if (ierror .ne. 0) then
        write (stderr, *) 'An error occured while reading ', uwater
        call exit (1)
      endif

c     get the water parameters from the temporary buffer
      read (buffer, *)
     $    qO, qH,
     $    rOHnC, rHHnC,
     $    (parOO(i), i = 1, nparOO),
     $    (parOH(i), i = 1, nparOH),
     $    (parHH(i), i = 1, nparHH),
     $    doh0, alpha0,
     $    (irho1(i), irho2(i), ialpha(i), ppintr(i), i = 1, npotp)

      close (uwater)


      length = lentrm (fstart)
      open (ustart, err = 510, file = fstart, status = 'old')
      rewind (ustart)

c     read step number, box sizes and number of atoms
      read (ustart, *) nfi, boxx, boxy, boxz, no, nh, nc

c     some sanity check
      if (no .gt. mxno .or. nh .gt. mxnh) then
        write (stderr, '(2a, i4, a)') 'number of atoms exceeds ',
     $      'internal maximum of ', mxatom, ' atoms'
        call exit (2)
      endif

      if (mod (nh, no) .ne. 0) then
        write (stderr, *)
     $      'there is something wrong with the number of atoms:',
     $      'number of hydrogens not twice as much as oxygens'
        call exit (2)
      endif

      if (.not. ion .and. nc .gt. 0) then
        write (stderr, '(2a)')
     $      'the input coordinates obviously contain an ion, ',
     $      'but a simulation without an ion was requested'
        call exit (2)
      endif

c     number of waters and number of total atoms
      nw = no + nh
      nwc = nw + nc

      read (ustart, *) (xo(i), yo(i), zo(i), i = 1, nwc)
      read (ustart, *) (vxo(i), vyo(i), vzo(i), i = 1, nwc)
      read (ustart, *) (fxo(i), fyo(i), fzo(i), i = 1, nwc)
      read (ustart, *) (x(i), y(i), z(i), i = 1, nwc)
      read (ustart, *) (vx(i), vy(i), vz(i), i = 1, nwc)
      read (ustart, *) (fx(i), fy(i), fz(i), i = 1, nwc)

      close (ustart)

      return


  500 write (stderr, *) 'Could not open/read from: ', fwater(1:length)
      call exit (1)
  510 write (stderr, *) 'Could not open/read from: ', fstart(1:length)
      call exit (1)

 8000 format (1x, 2a, /, 1x, a)

      end
