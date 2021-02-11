program test_plugin
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  double precision :: get_two_e_integral
  double precision :: e, nuc
  integer :: i,j

  print *, 'Hello world'
  print*, mo_num

  e = 0d0
  do i=1,mo_num
    do j=1,mo_num
       e = e + mo_one_e_integrals(j,i) &
       * (one_e_dm_mo_beta(j,i,1) + one_e_dm_mo_alpha(j,i,1))
    enddo
  enddo 
  
  print*,e 
 
  nuc = nuclear_repulsion
 
  ! subroutine de verification
  ! call print_energy_components

end
