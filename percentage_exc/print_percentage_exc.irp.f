!program percentage_excitations
!
!  implicit none
!
!  double precision, allocatable :: percentage(:), T(:,:)
!  double precision :: accu
!  integer :: state, excitation_degree, max_exc_degree
!
!  max_exc_degree = 10
!
!  allocate(percentage(n_states), T(n_states, max_exc_degree))
!
!  do excitation_degree = 1, max_exc_degree
!    call percentage_excitation(excitation_degree, percentage)
!    do state = 1, n_states
!      T(state, excitation_degree) = percentage(state)
!    enddo
!  enddo
!
!  do state = 1, n_states
!    print*,''
!    write(*,'(A7,I2)') 'State:', state  
!    print*,'Exc      %Exc        sum%     100%-sum%'
!    accu = 100d0*psi_coef_sorted(1,state)**2
!    write (*, '(A2,I2,1pE12.4,1pE12.4,1pE12.4)') '%T', 0, 100d0*psi_coef_sorted(1,state)**2, accu, 100d0 - accu
!    do excitation_degree = 1, max_exc_degree 
!      accu = accu + T(state, excitation_degree)
!      write (*, '(A2,I2,1pE12.4,1pE12.4,1pE12.4)') '%T', excitation_degree, T(state, excitation_degree), accu, 100d0-accu
!    enddo
!  enddo
!  deallocate(percentage,T)
!
!end
!
!subroutine percentage_excitation(excitation_degree, percentage)
!
!  implicit none
!
!  integer :: i,j
!  integer :: degree
!  integer, intent(in) :: excitation_degree
!  integer :: state
!  double precision, intent(out) :: percentage(n_states)
!  double precision :: t1, t2
!
!  ! state 1 -> ground state
!
!  call wall_time(t1)
!  percentage = 0d0
!
!  do i = 2, n_det
!    ! Debug
!    !call print_det(psi_det_sorted(n_int,:,1),n_int)
!    !call print_det(psi_det_sorted(n_int,:,i),n_int)
!
!    ! Excitation degree between the ref and the det |i>
!    call get_excitation_degree(psi_det_sorted(n_int,1,1),psi_det_sorted(n_int,1,i),degree,n_int)
!
!    ! Debug
!    !print*,i,degree
!
!    ! Sum the corresponding psi_coef**2
!    do state = 1, n_states
!      if (degree == excitation_degree) then
!        percentage(state) = percentage(state) + psi_coef_sorted(i,state)**2
!      endif
!    enddo
!  enddo
!  percentage = percentage * 100d0
!  call wall_time(t2)
!  write(*,'(A12,I2,1pE12.4)') 'Total time T', excitation_degree, t2-t1
!
!end

program print_percentage_exc

  implicit none

  read_wf = .True.
  TOUCH read_wf

  call run_print_percentage_exc

end
