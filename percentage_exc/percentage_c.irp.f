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
      write (*, '(A2,A2,10(1pE12.4))') '%C', adjustl(exc), percentage(i,:)/percentage(1,:)
    enddo
  else
    do i = 1, min(max_exc_degree+1,nb_max_percentage)
      write (exc,'(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%C', adjustl(exc), percentage(i,:)/percentage(1,:)
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
  integer :: i, s, exc_degree, idx_hf

  percentage = 0d0

  ! %C(n,s_state) = \sum_i psi_coef(i,s)**2 s.t. excitation_degree(|HF>,|i>) = n

  ! Contribution of HF det
  call find_hf_in_wf(psi_det,N_det,N_int,idx_hf)
  do s = 1, n_states
    percentage(1,s) = psi_coef(idx_hf,s)**2
  enddo
  
  ! Others determinants
  do i = 1, n_det
    call get_excitation_degree(HF_bitmask, psi_det(1,1,i), exc_degree, n_int)
    if (exc_degree == 0) then
      cycle
    endif
    do s = 1, n_states
      percentage(exc_degree+1, s) = percentage(exc_degree+1, s) + psi_coef(i,s)**2
    enddo
  enddo

  percentage = percentage *100d0

end

subroutine find_hf_in_wf(psidet,Ndet,Nint,idx_hf)

  implicit none

  ! in
  integer(bit_kind), intent(in) :: psidet(Nint,2,Ndet) 
  integer, intent(in) :: Ndet                          
  integer, intent(in) :: Nint                          

  ! out
  integer :: idx_hf

  ! internal
  integer :: i_int,j
  logical :: is_hf

  idx_hf = 0
  ! %T_0 
  do j = 1, Ndet
    is_hf = .True.
    do i_int = 1, Nint
      if (psidet(i_int,1,j) /= HF_bitmask(i_int,1) .or. &
          psidet(i_int,2,j) /= HF_bitmask(i_int,2)) then
        is_hf = .False.
        exit
      endif
      if (.not. is_hf) then
        exit
      endif
    enddo
    if (.not. is_hf) then
      cycle
    endif
    idx_hf = j
    exit
  enddo

  if (idx_hf == 0) then
    print*,'Hf determinant not found, abort'
    call abort
  endif

end
