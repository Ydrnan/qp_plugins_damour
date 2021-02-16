subroutine dm_newton_test(R)
  implicit none
  double precision, intent(in) :: R(mo_num,mo_num)

  double precision, allocatable :: new_mos(:,:)

  allocate(new_mos(mo_num,mo_num))

  call dgemm('N','N',mo_num,mo_num,mo_num,1d0,R,mo_num,mo_coef,mo_num,0d0,new_mos,mo_num)
  print*,'saaaaave mos'
  mo_coef = new_mos
  call save_mos
end subroutine
