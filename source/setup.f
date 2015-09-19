      subroutine setup
c $Id: setup.f,v 1.3 2005/06/21 11:48:13 hal Exp $
c
c     setup: open files and pre-set constants
c

      implicit none

      external lentrm, getmas
      integer lentrm
      double precision getmas

      integer i, length
      double precision t1, t2

      include 'sizes.inc'
      include 'params.inc'
      include 'potfor.inc'
      include 'qmmm.inc'
      include 'intra.inc'
      include 'press.inc'
      include 'eqmot.inc'
      include 'atoms.inc'
      include 'consts.inc'
      include 'units.inc'
      include 'files.inc'



c     general output file
      length = lentrm (fout)
      open (uout, err = 500, file = fout, status = mode)
      rewind (uout)

c     energy
      length = lentrm(fenout)
      open (uenout, err = 510, file = fenout, status = mode)
      rewind uenout


c     retrieve atom masses from the internal periodic table
      if (amh .eq. 0.0D0) amh = getmas ('H ')
      if (amo .eq. 0.0D0) amo = getmas ('O ')

      if (ion .and. amion .eq. 0.0D0) then
        amion = getmas (ionnam)

        if (amion .eq. 0.0D0) then
          write (stderr, *) 'the atom ', ionnam, ' is unknown'
          call exit (2)
        endif
      endif

c     calculate volume and box lengths of the cubic box
      if (.not. iNPT) then
        Volume = (dble(no) * amo + dble(nh) * amh + dble(nc) * amion) /
     $      (NL * dens * 1.0D6)

        boxx = (Volume * 1.0D30)**(1.0D0/3.0D0)

c       we assume a cubic box
        boxy = boxx
        boxz = boxx
      else
        Volume = boxx * boxy * boxz / 1.0D30
      endif


c     trajectory, and velocity files
      if (ftraj .ne. ' ') then
        length = lentrm(ftraj)

        open (utraj, err = 520, file = ftraj, status = mode)
        rewind utraj

c       write number of atoms and box dimensions to trajectory file
        if (ion) then
          write (utraj, 8000)
     $        nwc, '3', 'O_CF2', no, 'H_CF2', nh, ionnam, nc
        else
          write (utraj, 8010)
     $        nwc, '2', 'O_CF2', no, 'H_CF2', nh
        endif
      endif

      if (fveloc .ne. ' ') then
        length = lentrm(fveloc)

        open (uveloc, err = 530, file = fveloc, status = mode)
        rewind uveloc

        write (uveloc, 8020) nwc
      endif

      if (fftraj .ne. ' ') then
        length = lentrm(fftraj)

        open (uftraj, err = 540, file = fftraj, status = mode)
        rewind uftraj

        write (uftraj, 8020) nwc
      endif

c     QM/MM energy file
      if (iqmmm .and. fqmen .ne. ' ') then
        length = lentrm(fqmen)

        open (uqmen, err = 550, file = fqmen, status = mode)
        rewind uqmen
      endif

c     umbrella statistics file
      if (iumb) then
        length = lentrm(fumbst)

        open (uumbst, err = 560, file = fumbst, status = mode)
        rewind uumbst

        length = lentrm(fumbtr)

        open (uumbtr, err = 570, file = fumbtr, status = mode)
        rewind uumbtr
      endif


c     initialize QM list
      do i = 1, nwc
        inqm(i) = .false.
      enddo


      if (ion) then
c       energy conversion from [kcal/mol] into 1.0D23 * [J/particle]
c       last (exponential) parameter is not an energy parameter!
        do i = 1, mxpiO - 1
          pionO(i) = pionO(i) / encnv
        enddo

c       last (exponential) parameter is not an energy parameter!
        do i = 1, mxpiH - 1
          pionH(i) = pionH(i) / encnv
        enddo
      endif

c     first three-body parameter
      if (ithree) then
        a3bd = a3bd / encnv
        d3bd = d3bd / encnv
        e3bd = e3bd / encnv

        kD3bd = d3bd * k3bd
        kE3bd = e3bd * k3bd
      endif

c     squared and cubed cutoff radius
      rsqcut = rcut**2
      rcut3  = rcut**3

c     relaxation times and timestep input in ps
      taut = taut * tops
      taup = taup * tops
      dt = dt * tops

c     equilibrium H-O-H angle used for intra-molecular potential
      alpha0 = alpha0 * torad

c     constants for velocity verlet
      dt2 = dt / 2.0D0
      dtsq2 = dt * dt2
      dtv = dt * vcnv
      dt2i = dt2 * vcnvi

c     cutoffs for CF2 water model (Coulomb all, non-Coulomb O-O) and
c     cutoffs for ion-O and ion-H (Coulomb and non-Coulomb)
      if (rcut .gt. boxx / 2.0D0) rcut = boxx / 2.0D0

      if (roff .gt. rcut) then
        write (stderr, *)
     $      'it doesn''t make sense to set roff greater then rcut'
        call exit(2)
      endif

c     calculate coefficients for the predictor and corrector subroutines
      predp1 = 2.0D0 * dt * vcnv
      corp1 = dt * vcnv

c     oxygens
      t1 = dt * vcnvi * NL / (2.0D0 * amo)
      t2 = t1 / 3.0D0

      do i = 1, no
        predp2(i) = t2
        corp2(i) = t2
        predp3(i) = t1
        corp3(i) = t1
        mass(i) = amo / NL
      enddo

c     hydrogens
      t1 = dt * vcnvi * NL / (2.0D0 * amh)
      t2 = t1 / 3.0D0

      do i = no + 1, nw
        predp2(i) = t2
        corp2(i) = t2
        predp3(i) = t1
        corp3(i) = t1
        mass(i) = amh / NL
      enddo

c     ions
      if (ion) then
        t1 = dt * vcnvi * NL / (2.0D0 * amion)
        t2 = t1 / 3.0D0

        do i = nw + 1, nwc
          predp2(i) = t2
          corp2(i) = t2
          predp3(i) = t1
          corp3(i) = t1
          mass(i) = amion / NL
        enddo
      endif

      return


  500 write (stderr, *) 'Could not open/write to: ', fout(1:length)
      call exit (1)
  510 write (stderr, *) 'Could not open/write to: ', fenout(1:length)
      call exit (1)
  520 write (stderr, *) 'Could not open/write to: ', ftraj(1:length)
      call exit (1)
  530 write (stderr, *) 'Could not open/write to: ', fveloc(1:length)
      call exit (1)
  540 write (stderr, *) 'Could not open/write to: ', fftraj(1:length)
      call exit (1)
  550 write (stderr, *) 'Could not open/write to: ', fqmen(1:length)
      call exit (1)
  560 write (stderr, *) 'Could not open/write to: ', fumbst(1:length)
      call exit (1)
  570 write (stderr, *) 'Could not open/write to: ', fumbtr(1:length)
      call exit (1)

 8000 format (1x, i4, 1x, a, 3(1x, a, 1x, i4))
 8010 format (1x, i4, 1x, a, 2(1x, a, 1x, i4))
 8020 format (1x, i4)

      end
