program index
  implicit none
  integer :: i,j,k
  integer, allocatable :: mat(:,:)

  allocate(mat(4,4))
 
  do i = 1, 4
    do j = 1, 4
       call mat_to_vec_index(i,j,k)
       mat(i,j) = k
    enddo
  enddo
  
  print*,'(i,j) -> k'
  do i = 1, 4
    print*, mat(i,:)
  enddo

  print*,''
  print*,'k -> i j'
  do k = 1, 4*3/2
    call vec_to_mat_index(k,i,j)
    print*,k,i,j
  enddo 

end
