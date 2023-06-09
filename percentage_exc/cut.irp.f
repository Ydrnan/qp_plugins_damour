program cut

  use bitmasks

  implicit none

  integer :: i,j,k,l,m,n,na,nb,nab,n_max,beg,sze,idx
  integer :: Na_save,Nb_save,N_det_save
  integer, allocatable :: iorder(:)
  double precision :: norm,  energy, e
  double precision, allocatable :: norm_sort(:)
  integer(bit_kind), allocatable :: psi_det_save(:,:,:), psi_a_save(:,:), psi_b_save(:,:), cp_det(:,:,:)
  double precision, allocatable :: psi_coef_save(:,:), psi_bilinear_matrix_save(:,:,:)
  double precision, allocatable :: e_line(:,:), e_column(:,:), cp_coef(:,:)
  double precision :: e_target
  double precision, allocatable :: e_tmp(:)
  integer, allocatable :: tmp_order(:)
  integer :: n_col, n_line, cp_ndet
  logical :: must_exit, allowed_to_exit
  double precision :: cost0, cost

  PROVIDE mo_two_e_integrals_in_map mo_one_e_integrals

  if (.not. read_wf) then
    stop 'Please set read_wf to true'
  endif

  e_target = target_energy_qmc

  if (e_target == 0d0) then
    stop 'Please set target_energy_qmc to the targeted energy'
  endif

  ! Initial cut to reduce the number of alpha/beta parts
  print*,'Initial cut...'
  n_max = nb_line_max_qmc

  call reduce_psi_bilinear(2*n_max)
  print*,'Done'

  print*,'N_det', N_det
  N_det_save = N_det
  Na_save = N_det_alpha_unique
  Nb_save = N_det_beta_unique

  ! Copy after initial cut
  allocate(psi_det_save(N_int,2,N_det),psi_coef_save(N_det,N_states))
  allocate(psi_bilinear_matrix_save(N_det_alpha_unique,N_det_beta_unique,N_states))
  psi_det_save = psi_det
  psi_coef_save = psi_coef
  if (.True.) then
    psi_bilinear_matrix_save = psi_bilinear_matrix 
  endif

  ! Add columns/lines with largest PT2 ?
  allocate(psi_a_save(N_int,Na_save),psi_b_save(N_int,Nb_save))
  allocate(e_line(Na_save,N_states),e_column(Nb_save,N_states))
  psi_a_save = psi_det_alpha_unique
  psi_b_save = psi_det_beta_unique

  N_det = 1
  TOUCH psi_det psi_coef N_det
  !print*,'norm1',dsqrt(sum(psi_coef(1:N_det,1)**2))
  
  allocate(e_tmp(Na_save+Nb_save))
  allocate(tmp_order(Na_save+Nb_save))
 
  n_max = 2
  must_exit = .False.
 
  do while (.True.) 

    ! pt2 contrib of each line/col
    call pt2_by_line(Na_save,psi_a_save,Nb_save,psi_b_save,psi_bilinear_matrix_save,e_line,e_column)

    ! Counts the non-zero contrib and sorts them
    n_line = 0
    n_col = 0

    do i = 1, Na_save
      tmp_order(i) = i
      e_tmp(i) = e_line(i,1)

      if (e_tmp(i) /= 0d0) then
        n_line = n_line + 1
      endif
    enddo
    do i = 1, Nb_save
      tmp_order(Na_save+i) = -i
      e_tmp(Na_save+i) = e_column(i,1)
      if (e_tmp(Na_save+i) /= 0d0) then
        n_col = n_col + 1
      endif
    enddo
    call dsort(e_tmp,tmp_order,Na_save+Nb_save)

    ! Copy the wf
    allocate(cp_coef(N_det,1),cp_det(N_int,2,N_det))   
    cp_coef = psi_coef
    cp_det = psi_det
    cp_ndet = N_det
 
    allowed_to_exit = .True.
    do while (.True.)
      
      call add_line_col_from_pt2(n_max,n_col,n_line,Na_save,psi_a_save, &
         Nb_save,psi_b_save,psi_bilinear_matrix_save,tmp_order)

      na = N_det_alpha_unique
      nb = N_det_beta_unique
      cost0 = elec_alpha_num**3 + elec_beta_num**3
      cost = (na-1) * elec_alpha_num**2 + &
             (nb-1) * elec_beta_num**2 + &
              elec_alpha_num**3 + elec_beta_num**3
      cost = cost/cost0

     if (.True.) then
      !print*,'E:',n_max,N_det,psi_energy+nuclear_repulsion
      write(*,'(I10,F10.2,F20.10,F14.8)') N_det,cost,psi_energy+nuclear_repulsion,dsqrt(sum(psi_coef(1:N_det,:)**2))
     endif
      
     !print*,''
     !print*,'new'
     !do i = 1, N_det
     !  print*,psi_det(:,:,i)
     !enddo
     if (n_max < 2) then
       must_exit = .True.
       exit        
     endif

     ! Remove the added determinant
     if (.True.) then
        energy = psi_energy(1)+nuclear_repulsion
     endif
     if (energy < e_target) then
       n_max = n_max / 2
       psi_coef = cp_coef    
       psi_det = cp_det
       N_det = cp_ndet
       TOUCH psi_coef psi_det N_det
       allowed_to_exit = .False.
       cycle
     !else
     !  n_max = n_max * 2
     endif

     if (allowed_to_exit) then
       n_max = n_max * 2
       exit
     else
       deallocate(cp_coef,cp_det)
       allocate(cp_coef(N_det,1),cp_det(N_int,2,N_det))   
       cp_coef = psi_coef
       cp_det = psi_det
       cp_ndet = N_det
     endif
    enddo
    

    deallocate(cp_coef,cp_det)
    
    if (N_det == N_det_save) exit
    if (must_exit) exit
    !if (psi_energy(1)+nuclear_repulsion > e_target) cycle

  enddo

  ! Remove columns/lines by norm
  
end

subroutine add_line_col_from_pt2(n_max,n_col,n_line,Na_save,psi_a_save, &
        Nb_save,psi_b_save,psi_bilinear_matrix_save,tmp_order)

  implicit none

  integer, intent(in) :: Na_save,Nb_save,n_line,n_col,n_max
  integer(bit_kind), intent(in) :: psi_a_save(N_int,Na_save), psi_b_save(N_int,Nb_save)
  integer, intent(in) :: tmp_order(n_line+n_col)
  double precision, intent(in) :: psi_bilinear_matrix_save(Na_save,Nb_save,N_states)

  integer, allocatable :: idx_a(:), idx_b(:)
  integer(bit_kind), allocatable :: buff_a(:,:), buff_b(:,:), buff_det(:,:,:)
  double precision, allocatable :: buff_coef(:)

  integer :: i,j,k,l,a,b,idx_i,idx_j,ia,ib,na,nb
  integer, external :: get_index_in_psi_det_alpha_unique, get_index_in_psi_det_beta_unique
  logical :: is_eq
  logical, external :: is_in_wavefunction, is_in_psi_det_alpha_unique, is_in_psi_det_beta_unique

  allocate(idx_a(n_line),idx_b(n_line))
  allocate(buff_a(N_int,n_line),buff_b(N_int,n_col))

  na = N_det_alpha_unique
  nb = n_det_beta_unique

  ia = 1
  ib = 1
  buff_a = 0
  buff_b = 0
  idx_a = 0
  idx_b = 0

  ! New a/b det
  do i = 1, n_col+n_line
    j = tmp_order(i)
    if (j > 0) then
      if (is_in_psi_det_alpha_unique(psi_a_save(1,j),N_int)) cycle
      buff_a(:,ia) = psi_a_save(:,j)
      idx_a(ia) = j
      ia = ia + 1
    else
      if (is_in_psi_det_beta_unique(psi_b_save(1,-j),N_int)) cycle
      buff_b(:,ib) = psi_b_save(:,-j)
      idx_b(ib) = -j
      ib = ib + 1
    endif
    if (ia+ib-2 == n_max) exit
  enddo
  ia = ia - 1
  ib = ib - 1
  
  allocate(buff_det(N_int,2,na*ib+ia*nb+ia*ib),buff_coef(na*ib+ia*nb+ia*ib))

  ! New_a - old_b
  k = 1
  do j = 1, N_det_beta_unique
    ! Index beta in psi_b_save
    do b = 1, Nb_save
      is_eq = .True.
      do l = 1, N_int
        if (psi_det_beta_unique(l,j) /= psi_b_save(l,b)) then
          is_eq = .False.
          exit
        endif
      enddo
      if (is_eq) exit
    enddo
    do i = 1, ia
      if (psi_bilinear_matrix_save(idx_a(i),b,1) == 0d0) cycle
      buff_det(:,1,k) = buff_a(:,i)
      buff_det(:,2,k) = psi_b_save(:,b)
      !if (is_in_wavefunction(buff_det(1,1,k),N_int)) cycle
      buff_coef(k) = psi_bilinear_matrix_save(idx_a(i),b,1)
      k = k + 1
    enddo
  enddo
  
  ! Old_a- New_b
  do j = 1, ib
    do i = 1, N_det_alpha_unique
      ! Index alpha in psi_a_save
      do a = 1, Na_save
        is_eq = .True.
        do l = 1, N_int
          if (psi_det_alpha_unique(l,i) /= psi_a_save(l,a)) then
            is_eq = .False.
            exit
          endif
        enddo
        if (is_eq) exit
      enddo
      if (psi_bilinear_matrix_save(a,idx_b(j),1) == 0d0) cycle
      buff_det(:,1,k) = psi_a_save(:,a)
      buff_det(:,2,k) = buff_b(:,j)
      !if (is_in_wavefunction(buff_det(1,1,k),N_int)) cycle
      buff_coef(k) = psi_bilinear_matrix_save(a,idx_b(j),1)
      k = k + 1
    enddo
  enddo

  ! new_a - new_b
  do j = 1, ib
    do i = 1, ia
      if (psi_bilinear_matrix_save(idx_a(i),idx_b(j),1) == 0d0) cycle
      buff_det(:,1,k) = buff_a(:,i)
      buff_det(:,2,k) = buff_b(:,j)
      if (is_in_wavefunction(buff_det(1,1,k),N_int)) cycle
      buff_coef(k) = psi_bilinear_matrix_save(idx_a(i),idx_b(j),1)
      k = k + 1
    enddo
  enddo
  k = k - 1
  call add_det_wf(buff_det,buff_coef,k)
  deallocate(buff_det,buff_coef,buff_a,buff_b,idx_a,idx_b)

end

subroutine reduce_psi_bilinear(N)

  implicit none
 
  integer, intent(in) :: N

  integer :: i,j,k,l,m,na,nb,nab
  integer, allocatable :: iorder(:)
  double precision, allocatable :: norm_sort(:)

  nab = n_det_alpha_unique+n_det_beta_unique
  allocate ( norm_sort(0:nab), iorder(0:nab) )

  norm_sort(0) = 0.d0
  iorder(0) = 0
  do i=1,n_det_alpha_unique
   norm_sort(i) = det_alpha_norm(i)
   iorder(i) = i
  enddo

  do i=1,n_det_beta_unique
   norm_sort(i+n_det_alpha_unique) = det_beta_norm(i)
   iorder(i+n_det_alpha_unique) = -i
  enddo

  call dsort(norm_sort(1),iorder(1),nab)

  na = n_det_alpha_unique
  nb = n_det_beta_unique
  do j=1,nab
    if (na+nb <= N) then
      exit
    endif
    i = iorder(j)
    if ((i<0).and.(nb>1)) then
      nb -= 1
      do k=1,n_det
        if (psi_bilinear_matrix_columns(k) == -i) then
          psi_bilinear_matrix_values(k,1) = 0.d0
        endif
      enddo
    else if ((i>0).and.(na>1)) then
      na -= 1
      do k=1,n_det
        if (psi_bilinear_matrix_rows(k) ==  i) then
          psi_bilinear_matrix_values(k,1) = 0.d0
        endif
      enddo
    endif
  enddo

  do k=1,N_states
    psi_coef(1:N_det,k) = psi_bilinear_matrix_values(1:N_det,k)
    call dset_order(psi_coef(1,k),psi_bilinear_matrix_order_reverse,N_det)
  enddo
  TOUCH psi_det psi_coef

  psi_det = psi_det_sorted
  psi_coef = psi_coef_sorted

  do m=1,n_det
   if (psi_coef_sorted(m,1) == 0.d0) exit
  enddo

  N_det = m-1
  TOUCH psi_det psi_coef N_det

end

subroutine add_det_wf(l_det,l_coef,N)

  implicit none

  integer, intent(in) :: N
  integer(bit_kind), intent(in) ::l_det(N_int,2,N)
  double precision, intent(in) :: l_coef(N,N_states)
  integer :: i,j,k,s,new
  logical :: is_in,is_eq
  double precision :: norm

  !print*,'new',new,N
  !if (new /= N) call abort
  !do i = 1, new
  !  !call print_det(cp_det(1,1,i),N_int)
  !  print*,cp_det(:,:,i)
  !enddo

  !print*,'N_det',N_det
  !do i = 1, N_det
  !!  call print_det(psi_det(1,1,i),N_int)
  !  print*,psi_det(:,:,i),psi_coef(i,1)
  !enddo
  !print*,'N',N
  !do i = 1, N
  !  call print_det(l_det(1,1,i),N_int)
  !  !print*,l_det(:,:,i),l_coef(i,1)
  !enddo

  ! Add det in the wf 
  call fill_H_apply_buffer_no_selection(N,l_det,N_int,0)
  norm = dsqrt(sum(psi_coef(1:N_det,:)**2))
  call copy_H_apply_buffer_to_wf
  psi_coef = psi_coef * norm
  TOUCH psi_det psi_coef
  !print*,N_det,new
  !print*,'N_det',N_det
  !do i = 1, N_det
  !  !call print_det(psi_det(1,1,i),N_int)
  !  print*,psi_det(:,:,i)
  !enddo

  ! Update psi_coef from l_coef
  do s = 1, N_states
    do i = 1, N
      psi_coef(N_det-N+i,s) = l_coef(i,s)
    enddo
  enddo
  do i = 1, N
    psi_det(:,:,N_det-N+i) = l_det(:,:,i)
  enddo
  TOUCH psi_det psi_coef
  !print*,'N_det',N_det
  !do i = 1, N_det
  !  call print_det(psi_det(1,1,i),N_int)
  !!  print*,psi_det(:,:,i),psi_coef(i,1)
  !enddo

  ! Sort by coef
  psi_det = psi_det_sorted
  psi_coef = psi_coef_sorted
  TOUCH psi_det psi_coef

end

subroutine pt2_by_line(Na_save,psi_a_save,Nb_save,psi_b_save,psi_bilinear_matrix_save,e_line,e_column)

  implicit none
 
  integer, intent(in) :: Na_save,Nb_save
  integer(bit_kind), intent(in) :: psi_a_save(N_int,Na_save), psi_b_save(N_int,Nb_save)
  double precision, intent(in) :: psi_bilinear_matrix_save(Na_save,Nb_save,N_states)
  double precision, intent(out) :: e_line(Na_save,N_states), e_column(Nb_save,N_states)
  double precision :: ihpsi(N_states), e(N_states)
  double precision :: ihi,c
  integer :: i,j,k,l,s,a,b
  integer(bit_kind) :: det(N_int,2)
  logical :: is_eq, is_in
  logical, external :: is_in_wavefunction
  logical, external :: is_in_psi_det_alpha_unique, is_in_psi_det_beta_unique

  e_column = 0d0
  e_line = 0d0

  !$omp parallel &
  !$omp shared(Na_save,Nb_save,psi_a_save,psi_b_save,psi_det,psi_coef,&
  !$omp psi_bilinear_matrix_save,N_det,N_states,psi_energy,&
  !$omp e_line,e_column,N_int) &
  !$omp private(det,e,ihpsi,ihi,i,j,c,s,is_in) &
  !$omp default(none)

  !$omp do
  do j = 1, Nb_save
    det(:,2) = psi_b_save(:,j)
    do i = 1, Na_save
      det(:,1) = psi_a_save(:,i)
      is_in = is_in_psi_det_alpha_unique(det(1,1),N_int)
      if (is_in) cycle
      c = 0d0
      do s = 1, N_states
        c = c + abs(psi_bilinear_matrix_save(i,j,s))
      enddo
      if (c==0d0) cycle
      if (is_in_wavefunction(det,N_int)) cycle
      call i_H_psi(det,psi_det,psi_coef,N_int,N_det,N_det,N_states,ihpsi)   
      call i_h_j(det,det,N_int,ihi)  
      do s = 1, N_states
        e = -abs(ihpsi(s)**2/(psi_energy(s) - ihi))
        !print*,ihpsi(s)**2
        !print*,(psi_energy(s) - ihi)
        e_line(i,s) = e_line(i,s) + e(s)
        !print*,i,e(s)!,psi_bilinear_matrix_save(i,j,s)
      enddo
    enddo
  enddo
  !$omp end do nowait

  !$omp do
  do j = 1, Nb_save
    det(:,2) = psi_b_save(:,j)
    is_in = is_in_psi_det_beta_unique(det(1,2),N_int)
    if (is_in) cycle
    do i = 1, Na_save
      det(:,1) = psi_a_save(:,i)
      c = 0d0
      do s = 1, N_states
        c = c + abs(psi_bilinear_matrix_save(i,j,s))
      enddo
      if (c==0d0) cycle
      if (is_in_wavefunction(det,N_int)) cycle
      call i_H_psi(det,psi_det,psi_coef,N_int,N_det,N_det,N_states,ihpsi)   
      call i_h_j(det,det,N_int,ihi)  
      do s = 1, N_states
        e = -abs(ihpsi(s)**2/(psi_energy(s) - ihi))
        !print*,ihpsi(s)**2
        !print*,(psi_energy(s) - ihi)
        e_column(j,s) = e_column(j,s) + e(s)
        !print*,i,e(s)!,psi_bilinear_matrix_save(i,j,s)
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  
end
