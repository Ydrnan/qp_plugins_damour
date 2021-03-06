subroutine rho_model(rho,nb_iter,prev_energy,e_model,cancel_step)

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
  double precision, intent(inout)  :: prev_energy
  double precision, intent(in)  :: e_model
  logical, intent(in) :: cancel_step
  ! rho : double precision, ratio :
  ! rho = (prev_energy - energy) / (prev_energy - e_modele)
  ! nb_iter : integer, number of iteration

  !==========
  ! Internal
  !==========
  double precision :: energy!, e_model
  integer :: i
  ! energy      : double precision, energy of the actual step
  ! prev_energy : double precision, energy of the previous step
  ! e_model     : double precision, predicted energy for the actual step

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
    energy = energy + ci_energy(i)*state_average_weight(i)
  enddo
  energy = energy / DBLE(N_states)

  print*, 'ci_energy', energy

  if ( nb_iter >= 1  .and. .not.(cancel_step)) then

    rho = (prev_energy - energy) / (prev_energy - e_model)

    print*, 'previous energy, prev_energy :', prev_energy
    print*, 'predicted energy, e_model :', e_model
    print*, 'real energy, energy :', energy
    print*, 'prev_energy - energy :', prev_energy - energy
    print*, 'e_model - energy :', prev_energy - e_model
    print*, 'Rho :', rho

  elseif (cancel_step) then

    print*,'cancel_step = ', cancel_step
    rho = 0.5d0
    print*,'rho initialized :', rho

    ! the previous energy doesn't change since the step is cancelled

    print*, 'prev_energy :', prev_energy

  else
 
    print*,'Iteration 0, initialization'
    ! nb_iter = 0 corresponds to the first step after initialization
    print*, 'energy -> prev_energy :', energy
    prev_energy = energy
    rho = 0.5d0 
    print*,'Rho initialized :', rho

  endif

  if (rho >= 0.1d0 .and. nb_iter >= 0 .and. .not. cancel_step) then

!    ! Replacement of the previous energy
     prev_energy = energy

    print*, 'energy -> prev_energy :', energy

  endif

  !if (debug) then
    print*,'---Leave rho_model---'
    print*,''
  !endif

end subroutine
