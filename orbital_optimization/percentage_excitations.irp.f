program percentage_excitations

  implicit none

  double precision :: percentage
  integer :: state, excitation_degree

 
  do state = 1, n_states
    print*,'State:', state
    do excitation_degree = 1, 2
      call percentage_excitation(excitation_degree, state, percentage)
      write (*, "(A12,I1,F10.4)") 'Percentage T',excitation_degree,percentage
    enddo
    print*,''
  enddo

end

subroutine percentage_excitation(excitation_degree, state, percentage)

  implicit none

  integer :: i,j
  integer :: degree
  integer, intent(in) :: excitation_degree
  integer, intent(in) :: state
  double precision, intent(out) :: percentage

  ! state 1 -> ground state

  percentage = 0d0

  do i = 2, n_det
    ! Debug
    !call print_det(psi_det_sorted(n_int,:,1),n_int)
    !call print_det(psi_det_sorted(n_int,:,i),n_int)

    ! Excitation degree between the ref and the det |i>
    call get_excitation_degree(psi_det_sorted(n_int,1,1),psi_det_sorted(n_int,1,i),degree,n_int)

    ! Debug
    !print*,i,degree

    ! Sum the corresponding psi_coef**2
    if (degree == excitation_degree) then
      percentage = percentage + psi_coef_sorted(i,state)**2
    endif
  enddo

end
