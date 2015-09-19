      program main
c $Id: main.f,v 1.16 2005/06/18 08:25:29 hal Exp $
c
c
c     Qion: classical MD or QM/MM-MD with a single ion in CF2 water or pure CF2 water
c
c
c     authors (all known to date):
c       (probably) started by Heinzinger's group in Mainz
c       Klaus R. Liedl & Teerakiat Kerdcharoen (original hotspot implementation)
c       Anan Tongraar
c       Hannes H. Loeffler:
c         o major rewrite into a "real" program (mostly new program now)
c         o simulation of pure water
c         o umbrella sampling
c         o point charges for QM/MM
c         o complete rewrite of QM/MM code
c         o complete rewrite of neighbour list code (molecule based possible)
c         o many extensions
c         o bug fixes!
c
c
c     The following references (Molecular Dynamics) have been used:
c
c     o [AT87]  M. P. Allen and D. J. Tildesley, "Computer Simulation
c               of Liquids", Oxford University Press, Oxford, 1987
c     o [RAP95] D. C. Rapaport, "The Art of Molecular Dynamics
c               Simulation", Cambridge University Press, 1995
c



      implicit none

      logical nbrnow

      include 'params.inc'
      include 'stpcnt.inc'



c
c     initialization
c

c     initialize neighbour list flag
      nbrnow = .true.

c     parse command line
      call cmdlin

c     read input files and perform some sanity checks
      call rdinp

c     open files and pre-set constants
      call setup

c     initialize shifted potential and force calculation
      call pinit

c     initialize RF potential and RF force calculation
      if (irf) call rfinit

c     time step counter
      if (init) nfi = 0

c     write parameters and values
      call wrout


c
c     start the simulation loop
c
      do cstep = 1, nstep
        nfi = nfi + 1

c       predict new configuration
        if (eqmot .eq. 1) then
          call vverl1
        else if (eqmot .eq. 2) then
          call pred
        endif

c       recreate the neighbour list if necessary
        if (nbrnow) then
          nbrnow = .false.

          if (iatmbc) then
            call nbratb
          else
            call nbrwt
            call nbrion
          endif
        endif

c       create O list of QM zone
        if (iqmmm) call qmlist

c       calculate the potentials and forces
        call potfor

c       correct configuration
        if (eqmot .eq. 1) then
          call vverl2
        else if (eqmot .eq. 2) then
          call corr
        endif

c       remove translational center-of-mass motion
        if (mod(cstep, ncom) .eq. 0) call momz

c       write statistics about current and recent configurations
        call wrstat(nbrnow)

c       calculate umbrella sampling statistics
        if (iumb) call umbst

c       write restart configuration
        if (mod(cstep, nrst) .eq. 0) call rstart

c       write trajectory, velocities and forces if requested
        if (mod(nfi, noutcv) .eq. 0) call wrtraj

      enddo

c
c     happy end...
c
c     output calculation time and terminate program
c
      call finish(0)


      end
