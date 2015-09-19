c $Id: consts.inc,v 1.2 2005/05/22 08:10:09 hal Exp $
c
c     physical constants taken from CRC 77th edition
c

      double precision pi, NL, Rgas, kB, qe, e0
      double precision coucnv, encnv, frccnv, J2cal, bohr2A
      double precision hart2J, hart2k, A2m, tops, torad, todeg
      double precision vcnv, vcnv2, vcnvi

c     pi
      parameter (pi = 3.141592653589793D0)

c     Loschmidt's constant (particles/mole)
      parameter (NL = 6.0221367D+23)

c     picoseconds
      parameter (tops = 1.0D-12)

c     angle unit conversion
      parameter (torad = pi / 180.0D0, todeg = 180.0D0 / pi)

c     gas constant [J/(K*mol)] (used only here internally)
      parameter (Rgas = 8.31451D0)

c     Boltzmann constant [J/K]
      parameter (kB = Rgas / NL)

c     unit charge [C] (used only here internally)
      parameter (qe = 1.60217733D-19)

c     permittivity of vacuum [F/m] (used only here internally)
      parameter (e0 = 8.854187817D-12)

c     Ångström into meter (used only here internally)
      parameter (A2m = 1.0D-10)

c     bohr into Ångström (used only here internally)
      parameter (bohr2A = 0.529177249D0)

c     Joule into calories (used only here internally)
      parameter (J2cal = 0.239006D0)

c     hartree into Joule/mol (used only here internally)
      parameter (hart2J = 2625500.02D0)

c     hartree into kcal/mol
      parameter (hart2k = 627.51D0)

c     Coulomb energy conversion factor: q_1 * q_2 * coucnv
c     [J/particle] -> 1.0D23 * [J/(Å*particle)] (input energy parameters
c     are in 1.0D23 * [J])
c     q_1 and q_2 are partial atom charges
      parameter (coucnv = qe**2 / (4.0D0 * pi * e0) * (1.0D23 / A2m))

c     energy conversion factor: 1.0D23 * [J/particle] -> [kcal/mol]
      parameter (encnv = NL * J2cal / (1.0D3 * 1.0D23))

c     force conversion factor: [hartree/bohr] -> 1.0D23 * [J/Å]
      parameter (frccnv = 1.0D23 * hart2J / (bohr2A * NL))

c     velocity conversion factor
      parameter (vcnv = 1.0D8, vcnv2 = vcnv**2, vcnvi = 1 / vcnv)


c     Local Variables:
c     mode: fortran
c     End:
