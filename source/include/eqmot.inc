c $Id: eqmot.inc,v 1.3 2005/05/24 10:49:24 hal Exp $
c
c     predictor-corrector coefficients
c
c     predp? and corp? are needed for the predictor/corrector method
c     dt* and mass are needed for the velocity Verlet algorithm
c

      double precision predp1, predp2, predp3, corp1, corp2, corp3, mass
      double precision dt2, dtsq2, dtv, dt2i

      common /eqmot/ predp2(mxatom), predp3(mxatom), corp2(mxatom),
     $    corp3(mxatom), predp1, corp1, mass(mxatom), dt2, dtsq2, dtv,
     $    dt2i

      save /eqmot/



c     local variables:
c     mode: fortran
c     end:
