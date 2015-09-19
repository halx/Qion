      subroutine wrout
c $Id: _wrout.f,v 1.34 2005/06/21 11:48:13 hal Exp $
c
c     wrout: write information banner
c

      implicit none

      intrinsic HostNm
      external lentrm, lfttrm
      integer HostNm, lentrm, lfttrm

      integer i, err, start, length, len2, lenion, icharg
      character chstr*2, date*24, hname*60, user*20, buffer*80

      include 'sizes.inc'
      include 'params.inc'
      include 'potfor.inc'
      include 'umbrella.inc'
      include 'intra.inc'
      include 'press.inc'
      include 'atoms.inc'
      include 'stpcnt.inc'
      include 'consts.inc'
      include 'units.inc'
      include 'files.inc'



c     write current date, hostname and user name to output file
      call FDate(date)
      call GetLog(user)
      length = lentrm(user)
      err = HostNm (hname)

      if (err .eq. 0) then
        len2 = lentrm(hname)

        write (uout, '(1x, 6a)') 'program started ', date, ' on ',
     $      hname(1:len2), ' by ', user(1:length)
      else
        write (uout, '(1x, 4a)') 'program started ', date, ' by ',
     $      user(1:length)
      endif

      write (uout, *)

c     write banner
      write (uout, '(1x, 2a)') '+-------------------------------------',
     $    '---------------------------------------+'
      write (uout, 8000) '|', '|'
      write (uout, '(1x, a, t31, a, t79, a)')
     $    '|', 'Q I O N   (QM/MM-MD)', '|'
      write (uout, 8000) '|', '|'
      write (uout, 8000)
     $    '|   program version: @VERSION@',
     $    '|'
      write (uout, 8000)
     $    '|   build system:    @KERNEL@/@MACHINE@',
     $    '|'
      write (uout, 8000)
     $    '|   build host:      @HOST@',
     $    '|'
      write (uout, 8000)
     $    '|   build date:      @DATE@',
     $    '|'
      write (uout, 8000) '|', '|'
      write (uout, '(1x, a, t49, a, t79, a)')
     $    '|', '(c) 2000-2005 Hannes Loeffler', '|'
      write (uout, '(1x, 2a)') '+-------------------------------------',
     $    '---------------------------------------+'
      write (uout, '(2(/))')


c     name of ion and its charge
      lenion = lentrm(ionnam)
      icharg = nint(qion)

      if (icharg .gt. 0) then
        write (chstr, '(i1, a)') icharg, '+'
      endif

      if (icharg .lt. 0) then
        write (chstr, '(i1, a)') -icharg, '-'
      endif

      if (ion) then
        write (uout, *) 'ion to be investigated: ', ionnam(1:lenion),
     $      '(',chstr, ')'
      else
        write (uout, *) 'simulation of neat liquid'
      endif

      write (uout, *)

c     filenames
      write (uout, *) 'input files: '

      length = lentrm(fin)
      write (uout, 8010) 'control file:   ', fin(1:length)

      length = lentrm(fstart)
      write (uout, 8010) 'start rvf:      ', fstart(1:length)

      length = lentrm(fwater)
      write (uout, 8010) 'CF2 parameters: ', fwater(1:length)
      write (uout, *)

      if (mode .eq. 'unknown') then
        write (uout, *) 'output files (overwrite mode!):'
      else
        write (uout, *) 'output files:'
      endif

      length = lentrm(fout)
      write (uout, 8010) 'log file:       ', fout(1:length)

      length = lentrm(fenout)
      write (uout, 8010) 'energy file:    ', fenout(1:length)

      if (fqmen .ne. ' ') then
        length = lentrm(fqmen)
        write (uout, 8010) 'QM/MM energies: ', fqmen(1:length)
      else
        write (uout, 8010) 'QM/MM energies: *not written*'
      endif

      length = lentrm(finfo)
      write (uout, 8010) 'info file:      ', finfo(1:length)

      length = lentrm(frst)
      write (uout, 8010) 'restart rvf:    ', frst(1:length)

      if (ftraj .ne. ' ') then
        length = lentrm(ftraj)
        write (uout, 8010) 'trajectory:     ', ftraj(1:length)
      else
        write (uout, 8010) 'trajectory:     *not written*'
      endif

      if (fveloc .ne. ' ') then
        length = lentrm(fveloc)
        write (uout, 8010) 'velocities:     ', fveloc(1:length)
      else
        write (uout, 8010) 'velocities:     *not written*'
      endif

      if (fftraj .ne. ' ') then
        length = lentrm(fftraj)
        write (uout, 8010) 'forces:         ', fftraj(1:length)
      else
        write (uout, 8010) 'forces:         *not written*'
      endif

      if (iumb) then
        length = lentrm(fumbst)
        write (uout, 8010) 'umbrella histo: ', fumbst(1:length)

        if (fumbtr .ne. ' ') then
          length = lentrm(fumbtr)
          write (uout, 8010) 'umbrella trace: ', fumbtr(1:length)
        endif
      endif

      write (uout, '(/)')


c     echo the input file
      write (uout, *) 'this is the input file:'
      write (uout, '(2a)') '>-----------------------------------------',
     $    '-------------------------------------'

      length = lentrm(fin)
      open (uin, err = 500, file = fin, status = 'old')

  100 continue
      read (uin, 8015, end = 110) buffer
      length = lentrm(buffer)
      write (uout, 8015) buffer(1:length)
      goto 100
  110 continue

      close (uin)

      write (uout, '(2a)') '>-----------------------------------------',
     $    '-------------------------------------'
      write (uout, '(/)')


c     information about parameters for QM/MM calculation
      if (iqmmm) then
        write (uout, *) 'QM/MM calculation will be performed'

c       find last occurence of '/'
        do i = len (frcpth), 1, -1
          if (frcpth(i:i) .eq. '/') then
            start = i + 1
            goto 120
          endif
        enddo

  120   continue

        length = lentrm(frcpth)

        if (ipntch) then
          write (uout, *) '  point charges will be included'
        endif

        write (uout, *) '  external force program "',
     $      frcpth(start:length), '" will be used'

        if (iguess) then
          write (uout, *) '  user supplied initial guess will be used'
        else
          write (uout, *)
     $        '  initial guess will be calculated from scratch'
        endif

        write (uout, '(1x, a, f5.2, a)') '  QM/MM zone is ',
     $      2 * roff, ' A in diameter'
        write (uout,'(1x, 2(a, f5.2), a)')
     $      '  smoothing function will be applied between ', ron,
     $      ' A and ', roff, ' A'
        write (uout, *)
      endif

c     should three-body potentials be included?
      if (ithree) then
        write (uout, *) 'three-body potentials will be included'
      endif

c     should the intra-molecular potential be included?
      if (intra) then
        write (uout, *) 'intra-molecular potential will be applied'
      endif

c     should umbrella sampling be performed?
      if (iumb) then
        write (uout, *) 'umbrella biasing potential will be applied'
      endif

c     should reaction field be applied?
      if (irf) then
        write (uout, *) 'reaction field will be applied'
      endif

      write (uout, *)

c     title of parameter file (CF2 water model)
      length = lentrm(wattit)
      start = lfttrm(wattit)

      write (uout, '(1x, 2a)') 'the CF2 parameter file has the ',
     $    'title:'
      write (uout, *) '> ', wattit(start:length), ' <'
      write (uout, *)
      write (uout, *) 'CF2 charges and geometry:'
      write (uout, 8020) 'q_O = ', qO, ' a.u.'
      write (uout, 8020) 'q_H = ', qH, ' a.u.'
      write (uout, 8020) 'r_0 = ', doh0, ' A'
      write (uout, 8020) 'a_0 = ', alpha0 * todeg, ' degrees'
      write (uout, '(/)')

c     user supplied parameters
      if (ion) then
        write (uout, *) 'two-body functions:'
        write (uout, *)
        write (uout, '(7x, a)') 'q_1*q_2    A      B      C      D'
        write (uout, '(3x, 2a)') 'V = ------- + ---- + ---- + ---- +',
     $      ' ---- + E * exp(-F*r)'
        write (uout, '(10x, a)') 'r      r**a   r**b   r**c   r**d'
        write (uout, *)
        write (uout, '(10x, 2a, 6x, 2a)') ionnam(1:lenion),
     $      '-O two-body parameters:', ionnam(1:lenion),
     $      '-H two-body parameters:'

        do i = 1, mxpiO - 2
          write (uout, 8030) char (64 + i), '/', char (96 + i), ': ',
     $        pionO(i) * encnv, ' A**', ppionO(i), '*kcal/mole',
     $        pionH(i) * encnv, ' A**', ppionH(i), '*kcal/mole'
        enddo

        write (uout, 8040) 'E: ',
     $        pionO(mxpiO-1) * encnv, ' kcal/mole      ',
     $        pionH(mxpiH-1) * encnv, ' kcal/mole      '
        write (uout, 8040) 'F: ',
     $        pionO(mxpiO), ' A**(-1)        ', pionH(mxpiH), ' A**(-1)'
        write (uout, '(/)')
      endif

      if (ithree) then
        write (uout, *) 'three-body correction function:'
        write (uout, '(t56, a)') 'D        D        E'
        write (uout, '(3x, 2a)') 'V = A * exp(-B*r12) * exp(-B*r13)',
     $      ' * exp(-C*r23) - ------ - ------ - ------'
        write (uout, '(t54, a)') 'r12**k   r13**k   r23**k'
        write (uout, *)
        write (uout, '(10x, 2a)') ionnam(1:lenion),
     $      '-water three-body parameters:'
        write (uout, 8050) 'A = ', a3bd * encnv, ' kcal/mole'
        write (uout, 8050) 'B = ', b3bd, ' A**(-1)'
        write (uout, 8050) 'C = ', c3bd, ' A**(-1)'
        write (uout, 8060) 'D = ', d3bd * encnv, ' A*kcal/mole'
        write (uout, 8060) 'E = ', e3bd * encnv, ' A*kcal/mole'
        write (uout, '(3x, a, i2)') 'k = ', k3bd
        write (uout, '(/)')
      endif

      if (iumb) then
        write (uout, '(1x, a, 2(a, i5))') 'umbrella biasing potential ',
     $      'applied to distance ', nuidx1, ' -- ', nuidx2
        write (uout, '(1x, a, i7, a)')
     $     'first ', numbeq, ' steps of each window are equilibration'

        if (nwin .gt. 1) then
          write (uout, '(1x, a, i2, a, i8, a)') 'using ', nwin,
     $        ' windows of length ', winlen, ' each:'
        else
          write (uout, '(1x, a, i8)')
     $        'using a single window of length ', winlen
        endif

        write (uout, *)
        write (uout, '(1x, a)') '       k    d0'

        do i = 1, nwin
          write (uout, '(3x, f6.2, 1x, f5.2)') kumb(i) * encnv, d0umb(i)
        enddo

        write (uout, '(/)')
      endif


c     equation of motion
      if (eqmot .eq. 1) then
        write (uout, *) 'velocity Verlet algorithm is used'
      else if (eqmot .eq. 2) then
        write (uout, *) 'predictor/corrector algorithm is used'
      endif

      write (uout, *)

c     timesteps
      if (init) then
        write (uout, 8065)
     $      'timestep counting will start at step:', 1
      else
        write (uout, 8065)
     $      'timestep counting will start at step:', nfi
      endif

      write (uout, '(1x, a, 18x, i8)') 'number of timesteps:', nstep
      write (uout, '(1x, a, 19x, f6.4, 1x, a)') 'length of a timestep:',
     $    dt / 1.0D-15, 'fs'
      write (uout, '(1x, a, 2x, f10.4, 1x, a)')
     $    'total simulation time of this run:', nstep * dt / tops,
     $    'ps'
      write (uout, *)

c     temperature regulation
      write (uout, 8070) 'target temperature: ', Topt, ' K'
      write (uout, 8070) 'relaxation time:    ', taut / tops, ' ps'

c     pressure regulation
      if (iNPT) then
        write (uout, 8070) 'target pressure:    ', Popt, ' bar'
        write (uout, 8070) 'relaxation time:    ', taup / tops, ' ps'
        write (uout, '(1x, a, f9.3, a)')
     $      'compressibility:    ', Pcomp * 1.0D5, '*D-5 bar**-1'
      endif

      if (nscale .gt. 0) then
        write (uout, '(1x, 2a, i5, a)') 'solvent molecules will be ',
     $      'hard-scaled within the first ', nscale, ' steps'
      endif

      write (uout, *)

c     frequency of output
      write (uout, 8080) 'output updated every                   ',
     $    nout * nstat, ' steps'
      write (uout, 8080) 'statistics evaluated every             ',
     $     nstat, ' steps'

      if (naver .gt. 0) then
        write (uout, 8080) 'averages evaluated every               ',
     $      naver, ' steps'
      endif

      write (uout, 8080) 'restart coordinates updated every      ',
     $     nrst, ' steps'
      if (iwrap) write (uout, '(3x, a)') '(coordinates wrapped)'

      write (uout, 8080) 'translational CoM motion removed every ',
     $    ncom, ' steps'

      if (ftraj .ne. ' ') then
        write (uout, 8080) 'trajectory is updated every            ',
     $      noutcv, ' steps'
        if (iwrap) write (uout, '(3x, a)') '(coordinates wrapped)'
      endif

      if (fveloc .ne. ' ') then
        write (uout, 8080) 'velocities are updated every           ',
     $      noutcv, ' steps'
      endif

      if (fftraj .ne. ' ') then
        write (uout, 8080) 'force trajectory updated every         ',
     $      noutcv, ' steps'
      endif

      write (uout, *)

c     number of particles
      write (uout, '(1x, a, i5)') 'number of molecules: ', no + nc

      if (ion) then
        write (uout, '(3x, a, 8x, i5, 3a, f9.5)')
     $      'ions:      ', nc, ', m(', ionnam(1:lenion), ') = ', amion
      endif

      if (lenion .gt. 1) then
        write (uout, 8090) 'oxygens:   ', no, ', m(O)  = ', amo
        write (uout, 8090) 'hydrogens: ', nh, ', m(H)  = ', amh
      else
        write (uout, 8090) 'oxygens:   ', no, ', m(O) = ', amo
        write (uout, 8090) 'hydrogens: ', nh, ', m(H) = ', amh
      endif

      write (uout, *)

c     density and box dimensions
      if (iNPT) then
        write (uout, *) 'initial box dimensions:'
      else
        write (uout, '(1x, a, f7.5, a)') 'density of the solution: ',
     $      dens,' g/cm**3'
        write (uout, *)
        write (uout, *) 'box dimensions:'
      endif

      write (uout, 8070) '  x = ', boxx, ' A'
      write (uout, 8070) '  y = ', boxy, ' A'
      write (uout, 8070) '  z = ', boxz, ' A'
      write (uout, 8070) '  V = ', Volume * 1.0D30, ' A**3'
      write (uout, *)

c     cutoffs
      if (iatmbc) then
        write (uout, *) 'using atom based cutoffs'
      else
        write (uout, *) 'using molecule based cutoffs'
      endif

      write (uout, *)
      write (uout, '(4(1x, a, f6.3, a, /))')
     $    'general cutoffs:     r_cut      = ', rcut, ' A',
     $    'non-Coulomb cutoffs: r_cut(O-H) = ', rOHnC, ' A',
     $    '                     r_cut(H-H) = ', rHHnC, ' A',
     $    'overhang radius:     rplus      = ', rplus, ' A'

      if (ithree) then
        write (uout, '(1x, a, f6.3, a, /)')
     $      'three-body limit:    r_cut(3bd) = ', r3bd, ' A'
      endif

      write (uout, *)
      write (uout, *)
      write (uout, *) 'starting simulation...'
      write (uout, *)
      write (uout, '(1x, 2a)')  '-----------------------------------',
     $    '-----------------------------------------'

      call flush (uout)

      return


  500 write (stderr, *) 'Could not open/read from: ', fin(1:length)

 8000 format (1x, a, t79, a)
 8010 format (3x, 2a)
 8015 format (a)
 8020 format (3x, a, f9.5, a)
 8030 format (3x, 4a, f11.2, a, bz, i2, a, 4x, f11.2, a, i2, a)
 8040 format (3x, a, 2x, f11.2, a, 4x, f11.2, a)
 8050 format (3x, a, f10.7, a)
 8060 format (3x, a, f18.15, a)
 8065 format (1x, a, 1x, i8)
 8070 format (1x, a, f9.3, a)
 8080 format (1x, a, i7, a)
 8090 format (3x, a, 8x, i5, a, f9.5)

      end
