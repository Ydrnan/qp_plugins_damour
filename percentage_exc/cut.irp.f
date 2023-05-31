program cut

  use bitmasks

  implicit none

  integer :: i,j,k,l,m,n,na,nb,nab,n_max
  integer :: Na_save,Nb_save,N_det_save
  integer, allocatable :: iorder(:)
  double precision :: norm,  energy, e
  double precision, allocatable :: norm_sort(:)
  integer(bit_kind), allocatable :: psi_det_save(:,:,:), psi_a_save(:,:), psi_b_save(:,:)
  double precision, allocatable :: psi_coef_save(:,:), psi_bilinear_matrix_save(:,:,:)
  double precision, allocatable :: e_line(:,:), e_column(:,:)

  PROVIDE mo_two_e_integrals_in_map mo_one_e_integrals

  if (.not. read_wf) then
    stop 'Please set read_wf to true'
  endif

  ! Initial cut to reduce the number of alpha/beta parts
  print*,'Initial cut...'
  n_max = 10000

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
  psi_bilinear_matrix_save = psi_bilinear_matrix 

  ! Add columns/lines with largest PT2 ?
  allocate(psi_a_save(N_int,Na_save),psi_b_save(N_int,Nb_save))
  allocate(e_line(Na_save,N_states),e_column(Nb_save,N_states))
  psi_a_save = psi_det_alpha_unique
  psi_b_save = psi_det_beta_unique

  N_det = 1
  TOUCH psi_det psi_coef N_det
  print*,'norm1',dsqrt(sum(psi_coef(1:N_det,1)**2))
  
  double precision, allocatable :: e_tmp(:), buff_coef(:)
  integer, allocatable :: tmp_order(:),idx_a(:),idx_b(:)
  integer(bit_kind), allocatable :: buff_det(:,:,:), buff_a(:,:), buff_b(:,:)
  integer :: ia,ib,a,b,n1,n2,n3
  logical :: is_eq
  logical, external :: is_in_wavefunction
  
  allocate(e_tmp(Na_save+Nb_save))
  allocate(tmp_order(Na_save+Nb_save))
  allocate(buff_a(N_int,Na_save),buff_b(N_int,Nb_save))
  allocate(idx_a(Na_save),idx_b(Nb_save))

  
  do while (.True.)
    na = N_det_alpha_unique
    nb = N_det_beta_unique
    !print*,'ici'
    call pt2_by_line(Na_save,psi_a_save,Nb_save,psi_b_save,psi_bilinear_matrix_save,e_line,e_column)
     if (N_det == 9) exit
     do i = 1, Na_save
       tmp_order(i) = i
       e_tmp(i) = e_line(i,1)
     enddo
     do i = 1, Nb_save
       tmp_order(Na_save+i) = -i
       e_tmp(Na_save+i) = e_column(i,1)
     enddo
     call dsort(e_tmp,tmp_order,Na_save+Nb_save)
     !print*,tmp_order

     e = psi_energy(1) + nuclear_repulsion
     ia = 1
     ib = 1
     !print*,e_tmp(1:20)
     !print*,tmp_order(1:20)
     do i = 1, min(Na_save+Nb_save,2*N_det)
       e = e + e_tmp(i)
       j = tmp_order(i)
       if (e_tmp(i) == 0d0) exit
       if (j > 0) then
         buff_a(:,ia) = psi_a_save(:,j)
         idx_a(ia) = j
         ia = ia + 1
       else
         buff_b(:,ib) = psi_b_save(:,-j)
         idx_b(ib) = -j
         ib = ib + 1
       endif
       !print*,j
     enddo
     ia = ia - 1
     ib = ib - 1
     !print*,'ia/ib',ia,ib
     !print*,'e',e_tmp(1:ia+ib)
     !print*,'idx_a',idx_a(1:ia)
     !print*,'idx_b',idx_b(1:ib)
     
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
         !print*,'na',i,buff_a(:,i)
         buff_det(:,2,k) = psi_b_save(:,b)
         !if (is_in_wavefunction(buff_det(1,1,k),N_int)) then
         !   cycle
         !   print*,'p1'
         !   print*,'bilinear'
         !   do l = 1, Na_save
         !     write(*,'(100(E12.3))') psi_bilinear_matrix_save(l,:,1)
         !   enddo
         !   call abort
         !endif
         buff_coef(k) = psi_bilinear_matrix_save(idx_a(i),b,1)
         !print*,'1',idx_a(i),b,buff_coef(k),buff_det(:,:,k)
         k = k + 1
       enddo
     enddo
     !print*,'p1'
     !n1 = k-1
     !do i = 1, n1
     !   print*,buff_det(:,:,i)
     !!   !call print_det(buff_det(1,1,i),N_int)
     !enddo
   
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
         !print*,'nb',buff_b(:,j)
         !if (is_in_wavefunction(buff_det(1,1,k),N_int)) then
         !   cycle
         !   print*,'p2'
         !   call abort
         !endif
         buff_coef(k) = psi_bilinear_matrix_save(a,idx_b(j),1)
         !print*,'2',a,idx_b(j),buff_coef(k),buff_det(:,:,k)
         k = k + 1
       enddo
     enddo
     !print*,'p2'
     !n2 = k-1
     !do i = n1+1, n2
     !   print*,buff_det(:,:,i)
     !   !call print_det(buff_det(1,1,i),N_int)
     !enddo

     ! new_a - new_b
     do j = 1, ib
       do i = 1, ia
         !print*,idx_a(i),idx_b(j),psi_bilinear_matrix_save(idx_a(i),idx_b(j),1)
         if (psi_bilinear_matrix_save(idx_a(i),idx_b(j),1) == 0d0) cycle
         buff_det(:,1,k) = buff_a(:,i)
         buff_det(:,2,k) = buff_b(:,j)
         !if (is_in_wavefunction(buff_det(1,1,k),N_int)) then
         !   cycle
         !   print*,'p3'
         !   call abort
         !endif
         buff_coef(k) = psi_bilinear_matrix_save(idx_a(i),idx_b(j),1)
         !print*,'3',idx_a(i),idx_b(j),buff_coef(k),buff_det(:,:,k)
         k = k + 1
       enddo
     enddo
     !print*,'p3'
     !n3 = k-1
     !do i = n2+1, n3
     !   print*,buff_det(:,:,i)
     !   !call print_det(buff_det(1,1,i),N_int)
     !enddo
     !print*,N_det,k-1
     !prit*,k-1
     !print*,'k',k
     !do i = 1, k-1
     !  print*,'coef',buff_coef(i),buff_det(:,:,i)
     !!  call print_det(buff_det(1,1,i),N_int)
     !enddo
     call add_det_wf(buff_det,buff_coef,k-1)
     deallocate(buff_det,buff_coef)
     !call print_det(psi_det(1,1,1),N_int)
     print*,'E',N_det,psi_energy+nuclear_repulsion
     !print*,'bilinear'
     !do l = 1, N_det_alpha_unique
     !  write(*,'(100(E12.3))') psi_bilinear_matrix(l,:,1)
     !enddo
     if (N_det == N_det_save) exit
     !do i = 1, N_det
     !  print*,psi_coef(i,1),psi_det(:,:,i)
     !!  !print*,psi_det(:,:,i)
     !!  !call print_det(psi_det(1,1,i),N_int)
     !enddo
     !print*,'bilinear'
     !print*,idx_a(1:ia)
     !print*,idx_b(1:ib)
     !do i = 1, N_det_alpha_unique
     !  write(*,'(100(E12.3))') psi_bilinear_matrix(i,:,1)
     !enddo
     !if (N_det > 15)  call abort
    
  enddo

     !print*,'bilinear'
     !print*,'norm',dsqrt(sum(psi_bilinear_matrix_save(:,:,1)**2))
     !do l = 1, Na_save
     !  write(*,'(100(E12.3))') psi_bilinear_matrix_save(l,:,1)
     !enddo
  ! Remove columns/lines by norm
  
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

  !print*,psi_det_alpha_unique(:,1:N_det_alpha_unique)
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
      !if (e(1) /= 0d0) then
      !  print*,'col',det,is_in
      !endif
    enddo
  enddo

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
      !if (e(1) /= 0d0) then
      !  print*,'lin',det
      !endif
    enddo
  enddo
  
  !print*,'c',e_column
  !print*,'l',e_line
  
  !do i = 1, Na_save
  !  call print_det((/psi_a_save(1,i),psi_a_save(1,i)/),N_int)
  !  print*,'l',i,e_line(i,1)
  !enddo
  !do j = 1, Nb_save
  !  call print_det((/psi_b_save(1,j),psi_b_save(1,j)/),N_int)
  !  print*,'c',j,e_column(j,1)
  !enddo

end
