      subroutine wrstat (nbrnow)
c $Id: wrstat.f,v 1.19 2005/06/21 11:48:13 hal Exp $
c
c     wrstat: write energy averages (in kcal/mol) and temperature,
c             perform temperature coupling, check if neighbour list
c             should be updated
c
c     in: nbrnow - flag indicating if neighbour list should be recreated
c

      implicit none

      logical nbrnow

      external ssum
      double precision ssum

      integer i

      double precision vsqmax, maxdis
      double precision Tsys, Press, Vol, Etot, EMM, Epot, Ekin, EkO, EkH
      double precision Ekion
      double precision sEpOO, sEpOH, sEpHH, sEintr, sEpiO, sPress, sVol
      double precision sEpiH, sE3bd, sEbias, sEnrRF, sEnonC, sECou
      double precision sEkO, sEkH, sEkion, sEtot, sEQM
      double precision zEpOO, zEpOH, zEpHH, zEintr, zEpiO, zPress, zVol
      double precision zEpiH, zE3bd, zEbias, zEnrRF, zEnonC, zECou
      double precision zEkO, zEkH, zEkion, zEtot, zEQM
      double precision qEpOO, qEpOH, qEpHH, qEintr, qEpiO, qPress, qVol
      double precision qEpiH, qE3bd, qEbias, qEnrRF, qEnonC, qECou
      double precision qEkO, qEkH, qEkion, qEMM, qEtot, qEQM, qEpot
      double precision qEkin, qT
      double precision aEpot, aEpOO, aEpOH, aEpHH, aEintr, aEpiO, aEpiH
      double precision aEbias, aEnonC, aECou, aEnrRF, aEMM, aE3bd
      double precision aPress, aVol, aTH, aTO, aTion, aTtot
      double precision aEkin, aEkO, aEkH, aEkion, aEtot, aEQM

      save sEpOO, sEpOH, sEpHH, sEintr, sEpiO, sEpiH, sEnrRF, sEnonC
      save sECou, sEkO, sEkH, sEkion, sE3bd, sEbias, sPress, sVol
      save sEtot, sEQM
      save zEpOO, zEpOH, zEpHH, zEintr, zEpiO, zEpiH, zEnrRF, zEnonC
      save zECou, zEkO, zEkH, zEkion, zE3bd, zEbias, zPress, zVol
      save zEtot, zEQM
      save qEpOO, qEpOH, qEpHH, qEintr, qEpiO, qEpiH, qEnrRF, qEnonC
      save qECou, qEkO, qEkH, qEkion, qE3bd, qEbias, qPress, qVol
      save qEMM, qEtot, qEQM, qEpot, qEkin, qT
      save maxdis

      include 'params.inc'
      include 'sizes.inc'
      include 'rvf.inc'
      include 'energy.inc'
      include 'press.inc'
      include 'stpcnt.inc'
      include 'consts.inc'
      include 'units.inc'
      include 'files.inc'

      double precision vsq(mxatom)

      data sEnrRF, sEintr, sECou, sEnonC, sEpOO, sEpOH /6 * 0.0D0/
      data sEpHH, sEpiO, sEpiH, sEkO, sEkH, sEkion /6 * 0.0D0/
      data sE3bd, sEbias, sPress, sVol, sEtot, sEQM /6 * 0.0D0/
      data zEnrRF, zEintr, zECou, zEnonC, zEpOO, zEpOH /6 * 0.0D0/
      data zEpHH, zEpiO, zEpiH, zEkO, zEkH, zEkion /6 * 0.0D0/
      data zE3bd, zEbias, zPress, zVol, zEtot, zEQM /6 * 0.0D0/
      data qEnrRF, qEintr, qECou, qEnonC, qEpOO, qEpOH /6 * 0.0D0/
      data qEpHH, qEpiO, qEpiH, qEkO, qEkH, qEkion /6 * 0.0D0/
      data qE3bd, qEbias, qPress, qVol, qEMM, qEtot /6 * 0.0D0/
      data qEQM, qEpot, qEkin, qT /4 * 0.0D0/
      data maxdis /0.0D0/



c     sum the potential energies, in 1.0D23 * [J/particle]
      sEnrRF = sEnrRF + EnrRF
      sEintr = sEintr + Eintr
      sEnonC = sEnonC + EnonC
      sECou = sECou + ECou
      sEpOO = sEpOO + EpotOO
      sEpOH = sEpOH + EpotOH
      sEpHH = sEpHH + EpotHH
      sEpiO = sEpiO + EpotiO
      sEpiH = sEpiH + EpotiH
      sE3bd = sE3bd + V3bd
      sEbias = sEbias + Ebias

      zEnrRF = zEnrRF + EnrRF
      zEintr = zEintr + Eintr
      zEnonC = zEnonC + EnonC
      zECou = zECou + ECou
      zEpOO = zEpOO + EpotOO
      zEpOH = zEpOH + EpotOH
      zEpHH = zEpHH + EpotHH
      zEpiO = zEpiO + EpotiO
      zEpiH = zEpiH + EpotiH
      zE3bd = zE3bd + V3bd
      zEbias = zEbias + Ebias

      qEnrRF = qEnrRF + (EnrRF * encnv)**2
      qEintr = qEintr + (Eintr * encnv)**2
      qEnonC = qEnonC + (EnonC * encnv)**2
      qECou = qECou + (ECou * encnv)**2
      qEpOO = qEpOO + (EpotOO * encnv)**2
      qEpOH = qEpOH + (EpotOH * encnv)**2
      qEpHH = qEpHH + (EpotHH * encnv)**2
      qEpiO = qEpiO + (EpotiO * encnv)**2
      qEpiH = qEpiH + (EpotiH * encnv)**2
      qE3bd = qE3bd + (V3bd * encnv)**2
      qEbias = qEbias + (Ebias * encnv)**2

c     dynamic neighbour list update, see [RAP; pp. 52]
      vsqmax = 0.0D0

      do i = 1, nwc
        vsq(i) = vx(i)**2 + vy(i)**2 + vz(i)**2
        if (vsq(i) .gt. vsqmax) vsqmax = vsq(i)
      enddo

c     calculate maximum displacement
      maxdis = maxdis + sqrt (vsqmax) * dt * vcnv

      if (maxdis .gt. 0.5D0 * rplus) then
        nbrnow = .true.
        maxdis = 0.0D0
      endif

c     calculate the kinetic energies, in 1.0D23 * [J/particle]
      EkO = 0.5D0 * amo * ssum (vsq, 1, no) * vcnv2 / NL
      EkH = 0.5D0 * amh * ssum (vsq, no + 1, nw) * vcnv2 / NL
      Ekion = 0.5D0 * amion * ssum (vsq, nw + 1, nwc) * vcnv2 / NL
      Ekin = EkO + EkH + Ekion

c     accumulate kinetic energies
      sEkO = sEkO + EkO
      sEkH = sEkH + EkH
      sEkion = sEkion + Ekion

      zEkO = zEkO + EkO
      zEkH = zEkH + EkH
      zEkion = zEkion + Ekion

      qEkO = qEkO + (EkO * encnv)**2
      qEkH = qEkH + (EkH * encnv)**2
      qEkion = qEkion + (Ekion * encnv)**2
      qEkin = qEkin + (Ekin * encnv)**2

c     calculate instantaneous temperature and perform coupling
      Tsys = 2.0D0 * Ekin / (3.0D23 * dble(Nwc) * kB)
      qT = qT + Tsys**2
      call bath(Tsys)

c     calculate instantaneous pressure, energies in 1.0D23 * [J]
c     [N/m**2] = 1.0D5 * [Pa] -> [bar]
      Press = (2.0D0 * Ekin + virial) / (3.0D28 * Volume)
      if (iNPT) call pbath(Press)
      Vol = Volume * 1.0D30

c     accumulate pressure
      sPress = sPress + Press
      sVol = sVol + Vol
      zPress = zPress + Press
      zVol = zVol + Vol
      qPress = qPress + Press**2
      qVol = qVol + Vol**2

c     instantenous MM energies
      Epot = EpotOO + EpotOH + EpotHH + Eintr + EpotiO + EpotiH + EnrRF
     $    + V3bd + Ebias
      V0 = (EpotOO + EpotOH + EpotHH + Eintr + EpotiO + EpotiH + EnrRF
     $    + V3bd) * encnv
      EMM = (Epot + Ekin) * encnv

c     instantenous total energy
      Etot = EQM + EMM

      sEtot = sEtot + Etot
      sEQM = sEQM + EQM
      zEtot = zEtot + Etot
      zEQM = zEQM + EQM

      qEMM = qEMM + EMM**2
      qEQM = qEQM + EQM**2
      qEtot = qEtot + Etot**2
      qEpot = qEpot + (Epot * encnv)**2

c     write QM/MM energies _EVERY_ step
      if (iqmmm .and. fqmen .ne. ' ') then
        write (uqmen, 8000) nfi, natqm, Etot, EQM, EMM
      endif

c     compute average values and write them out
      if (mod(cstep, nstat) .eq. 0) then
        aEnrRF = encnv * sEnrRF / dble(nstat)
        aEintr = encnv * sEintr / dble(nstat)
        aEnonC = encnv * sEnonC / dble(nstat)
        aECou = encnv * sECou / dble(nstat)
        aEpOO = encnv * sEpOO / dble(nstat)
        aEpOH = encnv * sEpOH / dble(nstat)
        aEpHH = encnv * sEpHH / dble(nstat)
        aEpiO = encnv * sEpiO / dble(nstat)
        aEpiH = encnv * sEpiH / dble(nstat)
        aE3bd = encnv * sE3bd / dble(nstat)
        aEbias = encnv * sEbias / dble(nstat)

        aEtot = sEtot / dble(nstat)
        aEQM = sEQM / dble(nstat)

        aEkO = sEkO / dble(nstat)
c       "temperature" of Os
        aTO = aEkO * 2.0D-23 / (3.0D0 * dble(no) * kB)
        aEkO = aEkO * encnv

        aEkH = sEkH / dble(nstat)
c       "temperature" of Hs
        aTH = aEkH * 2.0D-23 / (3.0D0 * dble(nh) * kB)
        aEkH = aEkH * encnv

        if (ion) then
          aEkion = sEkion / dble(nstat)
c         "temperature" of ion
          aTion = aEkion * 2.0D-23 / (3.0D0 * dble(nc) * kB)
          aEkion = aEkion * encnv
        else
          aTion = 0.0D0
        endif

        aTtot = (dble(no) * aTO + dble(nh) * aTH +
     $      dble(nc) * aTion) / dble(nwc)

c       average pressure and volume
        aPress = sPress / dble(nstat)
        aVol = sVol / dble(nstat)

c       average potential energy
        aEpot = aEpOO + aEpOH + aEpHH + aEintr + aEpiO + aEpiH + aEnrRF
     $      + aE3bd + aEbias

c       average kinetic energy
        aEkin = aEkO + aEkH + aEkion

c       average MM energy
        aEMM = aEpot + aEkin


c       write info file
        open (uinfo, err = 500, file = finfo, status = 'unknown')
        rewind uinfo

        write (uinfo, 8010) 'step# = ', nfi, '|*NQM  = ', natqm,
     $      '|T     = ', aTtot, '|p     = ', aPress
        write (uinfo, 8020) 'E_MM  = ', aEMM, '|Etot  = ', aEtot,
     $      '|E_QM  = ', aEQM, '|V     = ', aVol
        write (uinfo, 8020) 'E_pot = ', aEpot, '|E_O-O = ', aEpOO,
     $      '|E_O-H = ', aEpOH, '|E_H-H = ', aEpHH
        write (uinfo, 8020) 'E_int = ', aEintr, '|E_i-O = ', aEpiO,
     $      '|E_i-H = ', aEpiH, '|E_bia = ', aEbias
        write (uinfo, 8020) 'E_kin = ', aEkin, '|E_k_O = ', aEkO,
     $      '|E_k_H = ', aEkH, '|E_k_i = ', aEkion
        write (uinfo, 8020) 'E_RF  = ', aEnrRF, '|E_Cou = ',
     $      aECou, '|E_nC  = ', aEnonC, '|E_3bd = ', aE3bd

        close (uinfo)

c       write to output file
        if (mod(cstep, nout * nstat) .eq. 0) then
          write (uout, *)
          write (uout, 8010) 'step# = ', nfi, '|*NQM  = ', natqm,
     $        '|T     = ', aTtot, '|p     = ', aPress
          write (uout, 8020) 'E_MM  = ', aEMM, '|Etot  = ', aEtot,
     $        '|E_QM  = ', aEQM, '|V     = ', aVol
          write (uout, 8020) 'E_pot = ', aEpot, '|E_O-O = ', aEpOO,
     $        '|E_O-H = ', aEpOH, '|E_H-H = ', aEpHH
          write (uout, 8020) 'E_int = ', aEintr, '|E_i-O = ', aEpiO,
     $        '|E_i-H = ', aEpiH, '|E_bia = ', aEbias
          write (uout, 8020) 'E_kin = ', aEkin, '|E_k_O = ', aEkO,
     $        '|E_k_H = ', aEkH, '|E_k_i = ', aEkion
          write (uout, 8020) 'E_RF  = ', aEnrRF, '|E_Cou = ',
     $        aECou, '|E_nC  = ', aEnonC, '|E_3bd = ', aE3bd
        endif

c       write to energy file
        write (uenout, 8030) nfi, natqm, aTtot, aPress, aVol, aEMM,
     $      aEpot, aEpOO, aEpOH, aEpHH, aEintr, aEpiO, aEpiH, aEbias,
     $      aEkin, aEkO, aEkH, aEkion, aEnrRF, aECou, aEnonC, aE3bd

        call flush(uenout)


        if (mod(cstep, naver) .eq. 0) then
          aEnrRF = encnv * zEnrRF / dble(naver)
          aEintr = encnv * zEintr / dble(naver)
          aEnonC = encnv * zEnonC / dble(naver)
          aECou = encnv * zECou / dble(naver)
          aEpOO = encnv * zEpOO / dble(naver)
          aEpOH = encnv * zEpOH / dble(naver)
          aEpHH = encnv * zEpHH / dble(naver)
          aEpiO = encnv * zEpiO / dble(naver)
          aEpiH = encnv * zEpiH / dble(naver)
          aE3bd = encnv * zE3bd / dble(naver)
          aEbias = encnv * zEbias / dble(naver)
          aEtot = zEtot / dble(naver)
          aEQM = zEQM / dble(naver)

          qEnrRF = sqrt((qEnrRF/dble(naver)) - aEnrRF**2)
          qEintr = sqrt((qEintr/dble(naver)) - aEintr**2)
          qEnonC = sqrt((qEnonC/dble(naver)) - aEnonC**2)
          qECou = sqrt((qECou/dble(naver)) - aECou**2)
          qEpOO = sqrt((qEpOO/dble(naver)) - aEpOO**2)
          qEpOH = sqrt((qEpOH/dble(naver)) - aEpOH**2)
          qEpHH = sqrt((qEpHH/dble(naver)) - aEpHH**2)
          qEpiO = sqrt((qEpiO/dble(naver)) - aEpiO**2)
          qEpiH = sqrt((qEpiH/dble(naver)) - aEpiH**2)
          qE3bd = sqrt((qE3bd/dble(naver)) - aE3bd**2)
          qEbias = sqrt((qEbias/dble(naver)) - aEbias**2)
          qEQM = sqrt((qEQM/dble(naver)) - aEQM**2)

          aEkO = zEkO / dble(naver)
c         "temperature" of Os
          aTO = aEkO * 2.0D-23 / (3.0D0 * dble(no) * kB)
          aEkO = aEkO * encnv
          qEkO = sqrt((qEkO/dble(naver)) - aEkO**2)

          aEkH = zEkH / dble(naver)
c         "temperature" of Hs
          aTH = aEkH * 2.0D-23 / (3.0D0 * dble(nh) * kB)
          aEkH = aEkH * encnv
          qEkH = sqrt((qEkH/dble(naver)) - aEkH**2)

          if (ion) then
            aEkion = zEkion / dble(naver)
c           "temperature" of ion
            aTion = aEkion * 2.0D-23 / (3.0D0 * dble(nc) * kB)
            aEkion = aEkion * encnv
            qEkion = sqrt((qEkion/dble(naver)) - aEkion**2)
          else
            aTion = 0.0D0
          endif

          aTtot = (dble(no) * aTO + dble(nh) * aTH +
     $        dble(nc) * aTion) / dble(nwc)
          qT = sqrt((qT/dble(naver)) - aTtot**2)

c         average pressure and volume
          aPress = zPress / dble(naver)
          aVol = zVol / dble(naver)

          qPress = sqrt((qPress/dble(naver)) - aPress**2)
          
          if (iNPT) then
            qVol = sqrt((qVol/dble(naver)) - aVol**2)
          else
            qVol = 0.0d0
          endif

c         average potential energy
          aEpot = aEpOO + aEpOH + aEpHH + aEintr + aEpiO + aEpiH
     $        + aEnrRF + aE3bd + aEbias
          qEpot = sqrt((qEpot/dble(naver)) - aEpot**2)

c         average kinetic energy
          aEkin = aEkO + aEkH + aEkion
          qEkin = sqrt((qEkin/dble(naver)) - aEkin**2)

c         average MM energy
          aEMM = aEpot + aEkin
          qEMM = sqrt((qEMM/dble(naver)) - aEMM**2)

          qEQM = sqrt((qEQM/dble(naver)) - aEQM**2)
          qEtot = sqrt((qEtot/dble(naver)) - aEtot**2)

          write (uout, *)
          write (uout, '(21x, a, i6, a)') 'averages over the last ',
     $        naver, ' steps:'
          write (uout, 8010) 'AVERAGES', nfi, '|*NQM  = ', natqm,
     $        '|T     = ', aTtot, '|p     = ', aPress
          write (uout, 8020) 'E_MM  = ', aEMM, '|Etot  = ', aEtot,
     $        '|E_QM  = ', zEQM, '|V     = ', aVol
          write (uout, 8020) 'E_pot = ', aEpot, '|E_O-O = ', aEpOO,
     $        '|E_O-H = ', aEpOH, '|E_H-H = ', aEpHH
          write (uout, 8020) 'E_int = ', aEintr, '|E_i-O = ', aEpiO,
     $        '|E_i-H = ', aEpiH, '|E_bia = ', aEbias
          write (uout, 8020) 'E_kin = ', aEkin, '|E_k_O = ', aEkO,
     $        '|E_k_H = ', aEkH, '|E_k_i = ', aEkion
          write (uout, 8020) 'E_RF  = ', aEnrRF, '|E_Cou = ',
     $        aECou, '|E_nC  = ', aEnonC, '|E_3bd = ', aE3bd
          write (uout, *)

          write (uout, '(17x, 2a, i6, a)') 'RMS fluctuations over ',
     $        'the last ', naver, ' steps:'
          write (uout, 8010) 'FLUCTS  ', nfi, '|*NQM  = ', natqm,
     $        '|T     = ', qT, '|p     = ', qPress
          write (uout, 8020) 'E_MM  = ', qEMM, '|Etot  = ', qEtot,
     $        '|E_QM  = ', qEQM, '|V     = ', qVol
          write (uout, 8020) 'E_pot = ', qEpot, '|E_O-O = ', qEpOO,
     $        '|E_O-H = ', qEpOH, '|E_H-H = ', qEpHH
          write (uout, 8020) 'E_int = ', qEintr, '|E_i-O = ', qEpiO,
     $        '|E_i-H = ', qEpiH, '|E_bia = ', qEbias
          write (uout, 8020) 'E_kin = ', qEkin, '|E_k_O = ', qEkO,
     $        '|E_k_H = ', qEkH, '|E_k_i = ', qEkion
          write (uout, 8020) 'E_RF  = ', qEnrRF, '|E_Cou = ',
     $        qECou, '|E_nC  = ', qEnonC, '|E_3bd = ', qE3bd
          write (uout, *)

c         reset sums for the next cycle
          zEnrRF = 0.0D0
          zEintr = 0.0D0
          zEnonC = 0.0D0
          zECou = 0.0D0
          zEpOO = 0.0D0
          zEpOH = 0.0D0
          zEpHH = 0.0D0
          zEpiO = 0.0D0
          zEpiH = 0.0D0
          zEkO = 0.0D0
          zEkH = 0.0D0
          zEkion = 0.0D0
          zE3bd = 0.0D0
          zEbias = 0.0D0
          zPress = 0.0D0
          zVol = 0.0D0
          zEtot = 0.0D0
          zEQM = 0.0D0

          qEnrRF = 0.0D0
          qEintr = 0.0D0
          qEnonC = 0.0D0
          qECou = 0.0D0
          qEpOO = 0.0D0
          qEpOH = 0.0D0
          qEpHH = 0.0D0
          qEpiO = 0.0D0
          qEpiH = 0.0D0
          qEkO = 0.0D0
          qEkH = 0.0D0
          qEkion = 0.0D0
          qE3bd = 0.0D0
          qEbias = 0.0D0
          qPress = 0.0D0
          qVol = 0.0D0
          qEtot = 0.0D0
          qEQM = 0.0D0
          qEMM = 0.0D0
          qEpot = 0.0D0
          qT = 0.0D0
          qEkin = 0.0D0
      endif

        call flush(uout)

c       check whether the H2O velocities should be hard-scaled
        if (cstep .le. nscale) then
          call scalex(aTO, aTH)
        endif

c       we perform the flush every nstat step and not in each step
        if (iqmmm .and. fqmen .ne. ' ') call flush(uqmen)

c       reset sums for the next cycle
        sEnrRF = 0.0D0
        sEintr = 0.0D0
        sEnonC = 0.0D0
        sECou = 0.0D0
        sEpOO = 0.0D0
        sEpOH = 0.0D0
        sEpHH = 0.0D0
        sEpiO = 0.0D0
        sEpiH = 0.0D0
        sEkO = 0.0D0
        sEkH = 0.0D0
        sEkion = 0.0D0
        sE3bd = 0.0D0
        sEbias = 0.0D0
        sPress = 0.0D0
        sVol = 0.0D0
        sEtot = 0.0D0
        sEQM = 0.0D0
      endif

      return


  500 write (stderr, *) 'Could not open: ', finfo
      call exit(1)

 8000 format (2(1x, i8), 3(1x, E20.10))
 8010 format (2(a, i11), 2(a, f11.2))
 8020 format (4(a, f11.2))
 8030 format (2(1x, i8), 21(1x, E20.10))

      end
