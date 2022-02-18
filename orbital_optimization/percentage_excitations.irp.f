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

program print_percentage_T

  implicit none

  integer :: i,s, max_exc_degree
  integer, allocatable :: list_states(:)
  double precision, allocatable :: percentage(:,:), accu(:)

  max_exc_degree = min(2*mo_num - elec_alpha_num - elec_beta_num, elec_alpha_num + elec_beta_num)
  allocate(percentage(max_exc_degree+1, n_states), accu(n_states), list_states(n_states))

  call percentage_excitations(max_exc_degree, percentage)

  do s = 1, n_states
   list_states(s) = s
  enddo   

  print*,''
  print*,'Percentage of the excitations per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  do i = 1, 5
    !write (exc, '(i4)') i-1
    write (*, '(A2,I2,10(1pE12.4))') '%T', i-1, percentage(i,:)    
  enddo

  print*,''
  print*,'Sum of the contributions per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  accu = 0d0
  do i = 1, 5
    do s = 1, n_states
      accu(s) = accu(s) + percentage(i,s)
    enddo
    write (*, '(A2,I2,10(F12.4))') '%T', i-1, accu(:)
  enddo

  print*,''
  print*,'Missing contributions per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  accu = 0d0
  do i = 1, 5
    do s = 1, n_states
      accu(s) = accu(s) + percentage(i,s)
    enddo
    write (*, '(A2,I2,10(1pE12.4))') '%T', i-1, 100d0-accu(:)        
  enddo

  deallocate(percentage, accu, list_states)

end

subroutine percentage_excitations(max_exc_degree, percentage)

  implicit none

  ! in 
  integer, intent(in) :: max_exc_degree
  ! max_exc_degree = min(2*mo_num - elec_alpha_num - elec_beta_num, elec_alpha_num + elec_beta_num)
  ! + 1 to have the %T0
  ! out
  double precision, intent(out) :: percentage(max_exc_degree + 1, n_states) 

  ! internal
  integer :: i, s, exc_degree

  percentage = 0d0

  ! %T(n,s_state) = \sum_i psi_coef_sorted(i,s)**2 s.t. excitation_degree(|HF>,|i>) = n

  ! Ref
  do s = 1, n_states
    percentage(1,s) = psi_coef_sorted(1,s)**2
  enddo

  ! Others
  do i = 2, n_det
    call get_excitation_degree(psi_det_sorted(n_int,1,1), psi_det_sorted(n_int,1,i), exc_degree, n_int)
    do s = 1, n_states
      percentage(exc_degree+1, s) = percentage(exc_degree+1, s) + psi_coef_sorted(i,s)**2
    enddo
  enddo

  percentage = percentage *100d0

end
