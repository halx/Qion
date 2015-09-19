      subroutine pred
c $Id: pred.f,v 1.4 2005/05/26 11:54:42 hal Exp $
c
c     pred: predictor method of Adams-Bashforth (sure?)
c
c     (this routine was written by Lutz Schaefer and copied by
c      Ecki, 1.4.1990)
c

      implicit none

      integer i

      include 'sizes.inc'
      include 'eqmot.inc'
      include 'rvf.inc'

      double precision xi, yi, zi
      double precision xold, yold, zold
      double precision vxold, vyold, vzold



      do i = 1, nwc
c       predict new positions
        xold = xo(i)
        yold = yo(i)
        zold = zo(i)

        call pbcmic(xold, yold, zold)

        xi = x(i)
        yi = y(i)
        zi = z(i)

        call pbcmic(xi, yi, zi)

        xo(i) = xi
        yo(i) = yi
        zo(i) = zi

        x(i) = xold + predp1 * (vxo(i) + 2.0D0 * predp2(i) * (2.0D0 *
     $      fx(i) + fxo(i)))
        y(i) = yold + predp1 * (vyo(i) + 2.0D0 * predp2(i) * (2.0D0 *
     $      fy(i) + fyo(i)))
        z(i) = zold + predp1 * (vzo(i) + 2.0D0 * predp2(i) * (2.0D0 *
     $      fz(i) + fzo(i)))

c       predict new velocities
        vxold = vxo(i)
        vyold = vyo(i)
        vzold = vzo(i)

        vxo(i) = vx(i)
        vyo(i) = vy(i)
        vzo(i) = vz(i)

        vx(i) = vxold + 4.0D0 * predp3(i) * fxo(i)
        vy(i) = vyold + 4.0D0 * predp3(i) * fyo(i)
        vz(i) = vzold + 4.0D0 * predp3(i) * fzo(i)

        fxo(i) = fx(i)
        fyo(i) = fy(i)
        fzo(i) = fz(i)
      enddo

      end
