BEGIN_PROVIDER [double precision, au2D]
&BEGIN_PROVIDER [double precision, au2eV]
&BEGIN_PROVIDER [double precision, Ha2nm]

  implicit none
  double precision :: Ha2J

  !BEGIN_DOC
  ! Atomic units to Debye
  ! Atomic units to eV
  ! Hartree to wavelength in nanometer
  !END_DOC

  ! 1 Ha 27.21138602d0 eV
  ! 1 au = 2.5415802529d0 D
  !planck_cte = 6.62606957d-34
  !light_speed = 2.99792458d10
  Ha2J = 4.35974434d-18

  au2D = 2.5415802529d0
  au2eV = 27.21138602d0
  Ha2nm = 1d9 * (planck_cte * light_speed) / Ha2J

END_PROVIDER

