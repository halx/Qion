c $Id: files.inc,v 1.7 2005/06/07 10:12:30 hal Exp $
c
c     static and dynamic file names
c

      character*(*) fwtCF2, fcoord, fforce, frcexe, fpntch

c     fwater and frcexe are searched below QION_ROOT
      parameter (fwtCF2 = 'data/CF2_water.dat',
     $    frcexe = 'scripts/calc_forces',
     $    fcoord = 'qmmm_coords',
     $    fforce = 'qmmm_forces',
     $    fpntch = 'qmmm_pntch')


      character*7 mode
      character*80 wattit
      character*128 frcpth, fin, fstart, frst, fout, finfo, fenout,
     $    ftraj, fveloc, fftraj, fqmen, fwater, fumbst, fumbtr

      common /files/ frcpth, fin, fstart, frst, fout, finfo, fenout,
     $    ftraj, fveloc, fftraj, fqmen, fwater, fumbst, fumbtr,
     $    wattit, mode 

      save /files/


c     Local Variables:
c     mode: fortran
c     End:
