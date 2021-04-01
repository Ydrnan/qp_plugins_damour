subroutine dn_rho_model(rho,nb_iter)

  include 'constants.h'

  implicit none

  !=============================================================================
  ! Compute the ratio : rho = (prev_energy - energy) / (prev_energy - e_modele)
  !=============================================================================

  !===========
  ! Intent in
  !===========
  double precision, intent(out) :: rho
  integer, intent(in)           :: nb_iter
  ! rho : double precision, ratio :
  ! rho = (prev_energy - energy) / (prev_energy - e_modele)
  ! nb_iter : integer, number of iteration

  !==========
  ! Internal
  !==========
  double precision :: energy, prev_energy, e_model
  ! energy      : double precision, energy of the actual step
  ! prev_energy : double precision, energy of the previous step
  ! e_model     : double precision, predicted energy for the actual step

  !open(unit=10,file='energy.dat')
  !read(10,*) energy
  !close(10)

  !=============
  ! Calculation
  !=============

  if (debug) then
    print*,'Enter in dn_rho_model'
  endif

  ! Energy of the actual step
  energy = sum(ci_energy(1:N_states) / dble(N_states))

  print*, 'ci_energy', energy

  if ( nb_iter >=1) then

    open(unit=11,file='prev_energy.dat')
      read(11,*) prev_energy
    close(11)

    open(unit=12,file='e_model.dat')
      read(12,*) e_model
    close(12)

    rho = (prev_energy - energy) / (prev_energy - e_model)

    print*, 'prev_energy :', prev_energy
    print*, 'e_model :', e_model
    print*, 'energy :', energy
    print*, 'prev_energy - energy :', prev_energy - energy
    print*, 'e_model - energy :', prev_energy - e_model
    print*, 'Rho :', rho

  elseif (nb_iter == -1) then

    ! nb_iter = -1 corresponds to a cancellation of the previous step
    rho = 0.5d0

    ! the previous energy doesn't change
    open(unit=11,file='prev_energy.dat')
      read(11,*) prev_energy
    close(11)

    print*, 'prev_energy :', prev_energy

  else

    ! nb_iter = 0 corresponds to the first step after initialization
    rho = 0.5d0

  endif

  if (rho >= 0.1d0 .and. nb_iter >= 0) then

    ! Replacement of the previous energy
    open(unit=11,file='prev_energy.dat')
      write(11,*) energy
    close(11)

    print*, 'energy -> prev_energy :', energy

  endif

  if (debug) then
    print*,'Leaves dn_rho_model'
  endif

end subroutine
