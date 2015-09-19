      block data prmini
c $Id: prmini.f,v 1.15 2005/06/10 10:47:12 hal Exp $
c
c     initial values of input parameters, potential parameters and
c     filenames
c

      implicit none

      include 'params.inc'
      include 'potfor.inc'
      include 'umbrella.inc'
      include 'sizes.inc'
      include 'atoms.inc'
      include 'files.inc'


c     optional input parameters with reasonable defaults
      data ion, init, intra, irf /4 * .true./
      data ithree, iqmmm, iguess, iumb, ipntch /5 * .false./
      data iatmbc, iwrap /2 * .true./
      data nstep, nout, wincnt /3 * 1/
      data nstat, nrst, noutcv /3 * 10/
      data nscale, nwin, numbeq, nuidx1, nuidx2 /5 * 0/
      data eqmot /1/
      data naver /100/
      data ncom /10000/
      data Topt /298.16D0/, taut /0.5D0/, dt /0.0002D0/, dens /0.997D0/
      data Popt /1.0D0/, taup /1.0D0/, Pcomp /4.6D-5/, iNPT /.false./
      data boxx, boxy, boxz /3 * 0.0D0/
      data r3bd /6.0D0/, hdel /0.1D0/

c     atom masses retrieved from a look-up table by default
      data amh, amo, amion /3 * 0.0D0/

c     these parameters *must* be set in the input file
      data qion, rcut, roff, ron /4 * 0.0D0/
      data a3bd, b3bd, c3bd, d3bd, e3bd /5 * 0.0D0/
      data k3bd /0/
      data ionnam /'  '/

c     ion-O and ion-H potentials read from the input file
      data pionO /mxpiO * 0.0D0/, ppionO /mxppiO * 0/
      data pionH /mxpiH * 0.0D0/, ppionH /mxppiH * 0/

c     additional radius for neighbour lists
      data rplus /0.4D0/

c     default access mode and default filenames
      data mode /'new'/, fin /'*stdin*'/, fstart /'qion.rvf'/,
     $    frst /'qion.rst'/, fout /'qion.out'/,
     $    finfo /'qion.info'/, fenout /'qion.en'/,
     $    ftraj /' '/, fveloc/' '/, fqmen /' '/, fwater /' '/,
     $    fftraj /' '/, fumbst /'qion_umb.hst'/, fumbtr /' '/

c     debugging flag
      data debug /0/

      end
