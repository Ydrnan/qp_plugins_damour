subroutine dm_newton_test(R)
  implicit none
  double precision, intent(in) :: R(mo_num,mo_num)

  double precision, allocatable :: new_mos(:,:)
  integer :: i
 
  !============
  ! Allocation
  !============
  allocate(new_mos(mo_num,mo_num))
  
  !=============
  ! Calculation
  !=============
  
  call dgemm('N','T',ao_num,mo_num,mo_num,1d0,mo_coef,size(mo_coef,1),R,size(R,1),0d0,new_mos,size(new_mos,1)) 

  !=========
  ! Storage
  !=========

  print*,'Save MOs...'
  
  mo_coef = new_mos
  print*, mo_coef(:,:)
  call save_mos
  
  print*,'Done, MOs saved'
 
  !==============
  ! Deallocation
  !==============

  deallocate(new_mos)

end subroutine
