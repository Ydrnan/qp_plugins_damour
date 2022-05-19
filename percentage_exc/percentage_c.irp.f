subroutine run_print_percentage_c

  implicit none

  integer :: i,s
  integer, allocatable :: list_states(:)
  double precision, allocatable :: percentage(:,:), accu(:)
  character(len=2) :: exc

  allocate(percentage(max_exc_degree+1, n_states), accu(n_states), list_states(n_states))

  call percentage_c(percentage)
  
  do s = 1, n_states
   list_states(s) = s
  enddo   

  print*,''
  print*,'Percentage of the excitations per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  if (percentage_in_exp) then
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%C', adjustl(exc), percentage(i,:)
    enddo
  else
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%C', adjustl(exc), percentage(i,:)
    enddo
  endif

  print*,''
  print*,'Percentage of the excitations'
  print*,'in intermediate normalization, %C0=1:'
  write(*,'(A4,10(I12))') '', 1
  if (percentage_in_exp) then
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%C', adjustl(exc), (percentage(i,:))*1d0/(percentage(1,:))
    enddo
  else
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%C', adjustl(exc), (percentage(i,:))*1d0/(percentage(1,:))
    enddo
  endif

  print*,''
  print*,'Sum of the contributions per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  accu = 0d0
  if (percentage_in_exp) then
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%C', adjustl(exc), accu(:)
    enddo
  else
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%C', adjustl(exc), accu(:)
    enddo
  endif

  print*,''
  print*,'Missing contributions per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  if (percentage_in_exp) then
    accu = 0d0
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%C', adjustl(exc), 100d0-accu(:)        
    enddo
  else
    accu = 0d0
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%C', adjustl(exc), 100d0-accu(:)        
    enddo
  endif

  deallocate(percentage, accu, list_states)

end

subroutine percentage_c(percentage)

  implicit none

  ! max_exc_degree = min(2*mo_num - elec_alpha_num - elec_beta_num, elec_alpha_num + elec_beta_num)
  ! + 1 to have the %C_0
  ! out
  double precision, intent(out) :: percentage(max_exc_degree + 1, n_states) 

  ! internal
  integer :: i, s, exc_degree

  percentage = 0d0

  ! %C(n,s_state) = \sum_i psi_coef(i,s)**2 s.t. excitation_degree(|HF>,|i>) = n

  ! Contribution of the ref
  do s = 1, n_states
    percentage(1,s) = psi_coef(1,s)**2
  enddo
  
  ! Others determinantss
  do i = 2, n_det
    call get_excitation_degree(psi_det(1,1,1), psi_det(1,1,i), exc_degree, n_int)
    do s = 1, n_states
      percentage(exc_degree+1, s) = percentage(exc_degree+1, s) + psi_coef(i,s)**2
    enddo
  enddo

  !percentage = dsqrt(percentage)

  percentage = percentage *100d0

end
