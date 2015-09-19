c $Id: params.inc,v 1.13 2005/06/10 10:45:24 hal Exp $
c
c     input parameters
c

      logical ion, init, intra, irf, ithree, iqmmm, iguess
      logical iumb, ipntch, iatmbc, iNPT, iwrap
      integer nscale, ncom, nstep, nstat, nout, nrst, noutcv
      integer naver, nwin, numbeq, nuidx1, nuidx2
      integer k3bd, eqmot, debug
      double precision Topt, taut, dt, dens, rcut, ron, roff
      double precision r3bd, amo, amh, amion, qion, a3bd, b3bd, c3bd
      double precision d3bd, e3bd, hdel, Popt, Pcomp, taup

      common /inpprm/ dt, Topt, taut, Pcomp, taup, Popt, dens, rcut,
     $    ron, roff, r3bd, amo,amh, amion, qion, a3bd, b3bd, c3bd, d3bd,
     $    e3bd, hdel, iNPT, ion, init, intra, irf, ithree, iqmmm,
     $    iguess, iumb, ipntch, iatmbc, iwrap, nscale, ncom, nstep,
     $    nstat, nout, nrst, noutcv, naver, nwin, numbeq,
     $    nuidx1, nuidx2, k3bd, eqmot, debug


      double precision rsqcut, rcut3, rplus, boxx, boxy, boxz

      common /prm2/ rsqcut, rcut3, rplus, boxx, boxy, boxz


c     constants for three-body potential
      double precision kD3bd, kE3bd

      common /pot3bd/ kD3bd, kE3bd


      save /inpprm/, /prm2/, /pot3bd/



c     local variables:
c     mode: fortran
c     end:
