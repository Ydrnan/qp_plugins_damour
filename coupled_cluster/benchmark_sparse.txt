program benchmarck_sparse

  use mkl_spblas
  use iso_c_binding, only :c_int, c_double

  implicit none

  double precision, allocatable :: A(:,:), B(:,:), C(:,:), list_A(:), tmp_list_A(:)
  integer, allocatable :: list_2d_key(:,:), list_key(:)
  integer, allocatable ::  tmp_list_row(:), list_row(:), nb_val_row(:), list_index_row(:), list_index_last(:)
  double precision :: t1,t2, epsilon
  integer :: i,j,k,l,n,p, nb_A, info, nb_row, row

  type(sparse_matrix_t)         :: sp_A
  type(matrix_descr)            :: descr

  n = 1000
  epsilon = 0.99d0

  allocate(A(n,n),B(n,n),C(n,n), tmp_list_A(n*n), list_key(n*n))

  call random_number(A)
  call random_number(B)

  ! dense
  call wall_time(t1)
  call dgemm('N','N',n,n,n, 1d0, A, size(A,1), B, size(B,1), 0d0, C, size(C,1))
  call wall_time(t2)
  print*,'dense',t2-t1
  !do i = 1, n
  !  write(*,'(100(1pE14.5))') C(i,:)
  !enddo

  ! Sparse   
  k = 1
  do j = 1, n
    do i = 1, n
      tmp_list_A(k) = A(i,j) 
      list_key(k) = k
      k = k + 1
    enddo
  enddo
  
  tmp_list_A = - tmp_list_A

  ! Sort by ascending order
  call dsort(tmp_list_A, list_key, n*n)
  tmp_list_A = -tmp_list_A

  k = 1
  do while ((k <= n*n) .and. tmp_list_A(min(k, n*n)) >= epsilon)
    k = k + 1
  enddo
  nb_A = k - 1  
  print*,'nb_A', nb_A, DBLE(nb_A)/DBLE(n*n)*100d0

  allocate(list_A(nb_A),list_2d_key(nb_A,2), tmp_list_row(nb_A))

  ! Row indexes
  do k = 1, nb_A
    l = list_key(k)
    call index_1d_to_2d(n, l, i, j)
    list_2d_key(k,1) = i
  enddo

  ! Column indexes
  do k = 1, nb_A
    l = list_key(k)
    call index_1d_to_2d(n, l, i, j)
    list_2d_key(k,2) = j
  enddo

  do k = 1, nb_A
     i = list_2d_key(k,1)
     j = list_2d_key(k,2)
     list_A(k) =  A(i,j) 
  enddo

  ! Sort the t2 by ascending order of indexes
  call sort_2d_key(n, nb_A, list_2d_key, list_A)

  !! ### TEST ###
  ! list of rows with non-zero elements
  list_row = 0
  do p = 1, nb_A
    i = list_2d_key(p,1)
    tmp_list_row(i) = i
  enddo

  nb_row = 0
  do i = 1, n
    if (tmp_list_row(i) /= 0) then
       nb_row = nb_row + 1
    endif
  enddo

  allocate(list_row(nb_row), nb_val_row(nb_row), list_index_row(nb_row), list_index_last(nb_row))

  j = 1
  do i = 1, n
    if (tmp_list_row(i) /= 0) then
      list_row(j) = tmp_list_row(i)
      j = j + 1
    endif
  enddo

  ! Number of non-zero elements per row
  i = 1
  nb_val_row = 0
  do p = 1, nb_A
    if (list_2d_key(p,1) == list_row(i)) then
      nb_val_row(i) = nb_val_row(i) + 1
    else
      i = i + 1
      nb_val_row(i) = nb_val_row(i) + 1
    endif
  enddo

  ! Index of each row in the vector
  list_index_row(1) = 1
  do i = 2, nb_row
    list_index_row(i) = list_index_row(i-1) + nb_val_row(i-1)
  enddo

  ! Index of the last non-zero element per row in the vector
  do i = 1, nb_row
    list_index_last(i) = list_index_row(i) + nb_val_row(i) - 1
  enddo

  !do i = 1, n 
  ! print*, 'A', A(i,:)
  !enddo

  !print*,'list_A', list_A(:)
  !print*,'row', list_row(:)
  !print*,'start', list_index_row(:)
  !print*,'end', list_index_last(:)
  

  ! vec -> sparse CSR
  ! info = mkl_sparse_d_create_csr(A, indexing, rows, cols, rows_start, rows_end, col_indx, values)
  !info = mkl_sparse_d_create_csr(sp_A, SPARSE_INDEX_BASE_ONE, n, n, list_index_row, list_index_last, ***, list_A)

  ! vec -> sparse COO
  info = mkl_sparse_d_create_coo(sp_A, SPARSE_INDEX_BASE_ONE, n, n, nb_A, list_2d_key(:,1), list_2d_key(:,2), list_A)
  print*,'info:',info
  descr%type = SPARSE_MATRIX_TYPE_GENERAL

  call wall_time(t1)
  ! Sparse-dense matrix multiplication
  info = mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE, 1d0, sp_A, descr, SPARSE_LAYOUT_ROW_MAJOR, B, n, size(B,1), 0d0, C, size(C,1))  
  call wall_time(t2)
  print*, 'sparse', t2 -t1
  !do i = 1, n
  !  write(*,'(100(1pE14.5))') C(i,:)
  !enddo

  C = 0d0
  call wall_time(t1)
  !do j = 1, n
  !  do i = 1, n
  !    do k = 1, n
  !      C(i,j) = C(i,j) + A(k,i) * C(k,j)
  !    enddo
  !  enddo
  !enddo

  do j = 1, n
    do i = 1, nb_row
      row = list_row(i)
      do k = 1, nb_val_row(i)
        l = list_index_row(i) + k - 1 ! index in list_A
        p = list_2d_key(l,2)          ! row index in B
        C(i,j) = C(i,j) + list_A(l) * B(p,j)
        !print*, 'C',i,j
        !print*, 'A',list_2d_key(l,1),list_2d_key(l,2)
        !print*, 'B', p,j 
      enddo
    enddo
  enddo
       
      
  call wall_time(t2)
  print*, 'test', t2 -t1
  !do i = 1, n
  !  write(*,'(100(1pE14.5))') C(i,:)
  !enddo

  deallocate(tmp_list_A,list_A,list_key,list_2d_key)

end

subroutine sort_2d_key(n, nb_t2, list_2d_key, list_t2)

  implicit none

  integer, intent(in)    :: nb_t2, n
  integer, intent(inout) :: list_2d_key(nb_t2,2)
  double precision, intent(inout) :: list_t2(nb_t2)
  integer, allocatable   :: tmp_list(:,:), index(:), key(:)
  double precision, allocatable :: tmp_t2(:)
  integer :: i,j,k,l

  allocate(tmp_list(nb_t2,2), index(nb_t2), key(nb_t2), tmp_t2(nb_t2))

  ! index to sort by (i,j) by ascending order (1,1), (1,2),...,(2,1), (2,2), ..., (n,n)
  do k = 1, nb_t2
    index(k) = (list_2d_key(k,1)-1) * n + list_2d_key(k,2)
  enddo

  ! sort
  call isort(index, key, nb_t2)

  ! tmp array contening the sorted key 
  do k = 1, nb_t2
    l = key(k)
    tmp_list(k,1) = list_2d_key(l,1)
    tmp_list(k,2) = list_2d_key(l,2)
    tmp_t2(k) = list_t2(l)
  enddo

  ! and put them in the array
  do k = 1, nb_t2
    list_2d_key(k,1) = tmp_list(k,1)  
    list_2d_key(k,2) = tmp_list(k,2)
    list_t2(k) = tmp_t2(k)
  enddo
  
  deallocate(tmp_list,index,key,tmp_t2)
  
end

subroutine index_1d_to_2d(n,k,i,j)

  implicit none

  integer, intent(in) :: n,k
  integer, intent(out) :: i,j

  ! k index in the list, list ordered column
  ! 1  p   ...
  ! 2  p+1 ...
  ! 3  p+2 ...
  ! :  :   ...  

  j = ((k-1)/n) + 1
  i = modulo((k-1),n) + 1
  
end
