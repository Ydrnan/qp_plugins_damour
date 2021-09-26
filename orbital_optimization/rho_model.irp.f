subroutine rho_model(prev_energy,e_model,rho)

  include 'constants.h'

  implicit none

  !=============================================================================
  ! Compute the ratio : rho = (prev_energy - energy) / (prev_energy - e_modele)
  !=============================================================================

  !====
  ! in
  !====
  double precision, intent(inout)  :: prev_energy
  double precision, intent(in)  :: e_model
  ! nb_iter : integer, number of iteration
  ! prev_energy : double precision, energy of the previous step
  ! e_model     : double precision, predicted energy for the actual step

  !=====
  ! out
  !=====
  double precision, intent(out) :: rho
  ! rho : double precision, representing the quality of the model

  !==========
  ! Internal
  !==========
  double precision :: energy
  integer :: i
  ! energy      : double precision, energy of the actual step

  !=============
  ! Calculation
  !=============

  !if (debug) then
    print*,''
    print*,'---Enter in rho_model---'
  !endif

  ! Energy of the actual step
  !energy = sum(ci_energy(1:N_states) / dble(N_states))

  energy = 0d0
  do i = 1, N_states
    energy = energy + ci_energy(i) * state_average_weight(i)
  enddo
  energy = energy / DBLE(N_states)

  print*, 'ci_energy :', energy

  rho = (prev_energy - energy) / (prev_energy - e_model)

  print*, 'previous energy, prev_energy :', prev_energy
  print*, 'predicted energy, e_model :', e_model
  print*, 'real energy, energy :', energy
  print*, 'prev_energy - energy :', prev_energy - energy
  print*, 'e_model - energy :', prev_energy - e_model
  print*, 'Rho :', rho

  ! Modification of prev_energy in function of rho
  if (rho < 0.1) then 
    print*, 'Rho < 0.1, the previous energy does not change'
    print*, 'prev_energy :', prev_energy  
  else
    prev_energy = energy
    print*, 'Rho >= 0.1, energy -> prev_energy :', energy
  endif

  !if (debug) then
    print*,'---Leave rho_model---'
    print*,''
  !endif

end subroutine
