program true_fci

  implicit none

  integer, allocatable :: occ_a(:,:), vir_a(:,:), det_a(:,:)
  integer, allocatable :: beg_degree_occ_a(:), n_degree_occ_a(:), end_degree_occ_a(:)
  integer, allocatable :: beg_degree_vir_a(:), n_degree_vir_a(:), end_degree_vir_a(:)
  integer, allocatable :: beg_degree_det_a(:), n_degree_det_a(:), end_degree_det_a(:)
  integer :: N_int_occ_a, N_int_vir_a, n_det_a
  integer :: n_unique_occ_a, n_unique_vir_a, n_unique_det_a
  integer :: i,j,k,l,m,n,p,q
  integer :: binom_coef
  integer :: degree, idx_degree, n_e
  integer :: max_exc

  max_exc = 5

  allocate(beg_degree_occ_a(max_exc+1), n_degree_occ_a(max_exc+1), end_degree_occ_a(max_exc+1))
  allocate(beg_degree_vir_a(max_exc+1), n_degree_vir_a(max_exc+1), end_degree_vir_a(max_exc+1))

  ! Number of configuration per excitation degree for occ/vir, a/b parts

  ! The i-th element indicates the total number of configurations
  ! with a degree inferior or equal to i-1

  ! occ_a
  !degree_occ_a(1) = 1
  !do n_e = elec_alpha_num-1, elec_alpha_num - max_exc, -1
  !  i = elec_alpha_num - n_e + 1
  !  degree_occ_a(i) = degree_occ_a(i-1) + binom_coef(n_e,elec_alpha_num)
  !enddo

  beg_degree_occ_a(1) = 1
  n_degree_occ_a(1) = 1
  end_degree_occ_a(1) = 1
  do n_e = elec_alpha_num-1, elec_alpha_num - max_exc, -1
    i = elec_alpha_num - n_e + 1
    beg_degree_occ_a(i) = end_degree_occ_a(i-1) + 1 
    n_degree_occ_a(i) = binom_coef(n_e,elec_alpha_num)
    end_degree_occ_a(i) = beg_degree_occ_a(i) + n_degree_occ_a(i) - 1
  enddo

  n_unique_occ_a = end_degree_occ_a(max_exc + 1)
  print*,'buffer size occ_a', n_unique_occ_a

  ! vir a
  !degree_vir_a(1) = 1
  !do n_e = 1, max_exc
  !  i = n_e + 1
  !  degree_vir_a(i) = degree_vir_a(i-1) + binom_coef(n_e,mo_num-elec_alpha_num)
  !enddo

  beg_degree_vir_a(1) = 1
  n_degree_vir_a(1) = 1
  end_degree_vir_a(1) = 1
  do n_e = 1, max_exc
    i = n_e + 1
    beg_degree_vir_a(i) = end_degree_vir_a(i-1) + 1 
    n_degree_vir_a(i) = binom_coef(n_e,mo_num-elec_alpha_num)
    end_degree_vir_a(i) = beg_degree_vir_a(i) + n_degree_vir_a(i) - 1
  enddo

  n_unique_vir_a = end_degree_vir_a(max_exc + 1)
  print*,'buffer size vir_a', n_unique_vir_a

  allocate(occ_a(elec_alpha_num,n_unique_occ_a),vir_a(mo_num-elec_alpha_num,n_unique_vir_a))

  do degree = 1, max_exc
  !  call gen_new_det_occ(degree)
  enddo

  ! Occupied alpha spin-orbitals
  print*,'### Occ alpha ###'
  occ_a = 0
  ! HF
  occ_a(:,1) = 1
  ! Excitations
  do n_e = elec_alpha_num-1, elec_alpha_num - max_exc, -1
    i = elec_alpha_num - n_e + 1
    !idx_degree = degree_occ_a(i-1) + 1
    idx_degree = beg_degree_occ_a(i)
    call gen_k_in_n(n_e,elec_alpha_num,idx_degree,n_unique_occ_a,occ_a)
  enddo
  print*,''  

  ! Virtual alpha spin-orbitals
  print*,'### Vir alpha ###'
  vir_a = 1
  ! HF
  vir_a(:,1) = 0
  ! Excitations
  do n_e = 1, max_exc
    i = n_e + 1
    !idx_degree = degree_vir_a(i-1) + 1
    idx_degree = beg_degree_vir_a(i) 
    call gen_k_in_n(n_e,mo_num-elec_alpha_num,idx_degree,n_unique_vir_a,vir_a)
  enddo

  print*,'Total number of configuration for occupied alpha:'
  print*,n_unique_occ_a
  print*,'Total number of configuration for virtual alpha:'
  print*,n_unique_vir_a
  print*,'Total number of alpha occ+vir possible for each excitation degree:'

  allocate(n_degree_det_a(max_exc + 1),beg_degree_det_a(max_exc+1),end_degree_det_a(max_exc+1))

  beg_degree_det_a(1) = 1
  end_degree_det_a(1) = 1
  n_degree_det_a(1) = 1
  n_det_a = 1
  do i = 2, max_exc + 1
    beg_degree_det_a(i) = end_degree_det_a(i-1) + 1
    n_degree_det_a(i) =  n_degree_occ_a(i) * n_degree_vir_a(i)
    end_degree_det_a(i) = beg_degree_det_a(i) + n_degree_det_a(i) - 1
    n_det_a = n_det_a + n_degree_det_a(i)
    print*, i-1, n_degree_det_a(i)
  enddo
  print*,'Total:',n_det_a

  ! Alpha determinants
  ! Generation of all the unique alpha determinants by combining
  ! the different occ and vir parts for each degree of excitation
  allocate(det_a(mo_num,n_det_a))

  det_a = 0
  j = 1
  do i = 1, max_exc + 1
    do while (j <= end_degree_det_a(i))
      do k = 1, n_degree_occ_a(i)
        do l = 1, n_degree_vir_a(i)
          do p = 1, elec_alpha_num
            det_a(p,j) = occ_a(p,beg_degree_occ_a(i) + k - 1)
          enddo
          do p = elec_alpha_num + 1, mo_num
            det_a(p,j) = vir_a(p-elec_alpha_num,beg_degree_vir_a(i) + l - 1)
          enddo
          write(*,'(100(I1))') det_a(:,j)
          j = j + 1
        enddo
      enddo
    enddo
  enddo

  !deallocate(occ_a,vir_a,degree_occ_a,degree_vir_a)

end

subroutine gen_k_in_n(k,n,idx_base,size_buffer,buffer_det)

  implicit none

  integer, intent(in) :: k, n, idx_base, size_buffer
  integer, intent(inout) :: buffer_det(n,size_buffer)
  integer, allocatable :: det(:)
  logical :: is_complete
  integer :: i,m
  integer :: ndet_kn

  if (k > n) then
    print*, 'k>n, canot generate any configuration, abort'
    call abort
  endif

  allocate(det(n))

  ! Generate the first configuration
  det = 0
  do i = 1, k
     det(i) = 1
  enddo
  
  print*,'================================================'
  write(*,'(I4,A7,I4,A14)') k, ' e- in ', n, ' spin-orbitals'
  write(*,'(100(I1))') det(:)

  ! First configuration
  m = 1
  buffer_det(:,idx_base + m - 1) = det(:)

  ! Loop to generate all the configuration of
  ! k e- in n spin-orbitals
  is_complete = .False.
  do while (.not. is_complete)
    ! Generate the next configuration
    call  next_det(det,n,is_complete)
    ! To avoid printing 2 times the last one
    if (is_complete) then
      exit
    endif
    write(*,'(100(I1))') det(:)

    ! Include the configuration in the buff
    m = m + 1
    buffer_det(:,idx_base + m - 1) = det(:) 
  enddo

  write(*,'(A34,I4,A7,I4,A15)') 'Total number of configurations for', k, ' e- in ',n,' spin-orbitals:'
  write(*,'(I10)') m
  print*,''

  deallocate(det)

end

subroutine gen_new_det_occ(degree)

  implicit none

  integer, intent(in) :: degree
  integer, allocatable :: det(:)
  integer :: i,j,k,l,m,n
  integer :: size_det, ndet_tot
  logical :: is_complete

  size_det = elec_alpha_num
  allocate(det(size_det))

  det = 0
  do i = 1, elec_alpha_num - degree
    det(i) = 1
  enddo
  print*,'Degree:',degree
  write(*,'(100(I1))') det(:)

  ndet_tot = 1
  is_complete = .False.
  do while (.not. is_complete)
    call  next_det(det,size_det,is_complete)
    ! To avoid printing 2 times the last one
    if (is_complete) then
      exit
    endif
    write(*,'(100(I1))') det(:)
    ndet_tot = ndet_tot + 1
  enddo
  print*,'Total number:', ndet_tot, degree
  print*,''
 
  deallocate(det)
  
end
 
subroutine next_det(det,size_det,is_complete)

  implicit none

  integer, intent(in) :: size_det
  integer, intent(inout) :: det(size_det)
  logical, intent(inout) :: is_complete
  integer :: z,o
  integer :: i,j

  call find_first_one(det,size_det,o,is_complete)
  call find_first_zero(det,size_det,o,z,is_complete)

  if (z == 0 .or. o == 0) then
    print*,'End'
    return
  endif

  det(z) = 1   ! 0 becomes 1 
  det(z-1) = 0 ! 1 becomes 0 

  ! Put all the 1 before z at the beginning of det
  j = 1
  do i= 1, z-2
    if (det(i) == 1) then
      det(i) = 0
      det(j) = 1
      j = j+1
    endif
  enddo
  
end

subroutine find_first_one(det,size_det,first_one,is_complete)

  implicit none

  integer, intent(in) :: size_det, det(size_det)
  integer, intent(out) :: first_one
  logical, intent(inout) :: is_complete
  integer :: i

  first_one = 0
  do i = 1, size_det
    if (det(i) == 1) then
      first_one = i
      return
    endif
  enddo

  if (first_one == 0) then
    is_complete = .True.
  endif

end

subroutine find_first_zero(det,size_det,one_idx,first_zero,is_complete)

  implicit none

  integer, intent(in) :: size_det, det(size_det), one_idx
  integer, intent(out) :: first_zero
  logical, intent(inout) :: is_complete
  integer :: i
  
  first_zero = 0
  do i = one_idx + 1, size_det
    if (det(i) == 0) then
      first_zero = i
      return
    endif
  enddo
   
  if (first_zero == 0) then
    is_complete = .True.
  endif

end

function binom_coef(k,n)

  implicit none

  integer,intent(in) :: k,n
  integer :: factorial
  integer :: binom_coef
  
  binom_coef = factorial(n)/(factorial(k)*factorial(n-k))

end

function factorial(n)

  implicit none

  integer, intent(in) :: n
  integer :: factorial

  integer :: i

  factorial = 1

  if (n <= 1) then 
    return
  endif

  do i = 1, n
    factorial = factorial * i
  enddo

end
