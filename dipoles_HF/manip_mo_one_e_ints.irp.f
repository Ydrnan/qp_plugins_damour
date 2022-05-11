subroutine run

  implicit none

  integer :: i,j,p,q
  double precision :: lambda

  lambda = 0d0

  print*,mo_one_e_integrals(1,1)

  call pert_mo_one_e_ints(10d0,1)
 ! print*,mo_one_e_integrals(1,1)

!  call init_mo_one_e_ints

!  call ao_to_mo(ao_one_e_integrals,size(ao_one_e_integrals,1),mo_one_e_integrals,size(mo_one_e_integrals,1))
!
  !print*,mo_one_e_integrals(1,1)
  !do j = 1, mo_num
  !  do i = 1, mo_num
  !    mo_one_e_integrals(i,j) = 1.2d0!mo_one_e_integrals(i,j) + lambda * mo_dipole_x(i,j)
  !  enddo
  !enddo

  !call ezfio_set_mo_one_e_ints_mo_one_e_integrals(mo_one_e_integrals)
  !call ezfio_set_mo_one_e_ints_io_mo_one_e_integrals('Read')
  !call ezfio_set_mo_one_e_ints_io_mo_integrals_kinetic('None')
  !call ezfio_set_mo_one_e_ints_io_mo_integrals_pseudo('None')
  !call ezfio_set_mo_one_e_ints_io_mo_integrals_n_e('None')
!
!  open(unit=10,file='H_mo.qp')
!  do j = 1, mo_num
!    do i = 1, mo_num
!      write(10,*)i,j,mo_one_e_integrals(i,j)
!    enddo
!  enddo
!  close(10)
  
end

subroutine pert_mo_one_e_ints(lambda,dir)

  implicit none

  double precision, intent(in) :: lambda
  integer, intent(in) :: dir

  integer :: i,j

  PROVIDE mo_one_e_integrals mo_dipole_x mo_dipole_y mo_dipole_z 

  call init_mo_one_e_ints

  if (dir == 1) then
  do j = 1, mo_num
    do i = 1, mo_num
      mo_one_e_integrals(i,j) = mo_one_e_integrals(i,j) + lambda * mo_dipole_x(i,j)
    enddo
  enddo
  elseif (dir == 2) then
  do j = 1, mo_num
    do i = 1, mo_num
      mo_one_e_integrals(i,j) = mo_one_e_integrals(i,j) + lambda * mo_dipole_y(i,j)
    enddo
  enddo
  elseif (dir == 3) then
  do j = 1, mo_num
    do i = 1, mo_num
      mo_one_e_integrals(i,j) = mo_one_e_integrals(i,j) + lambda * mo_dipole_z(i,j)
    enddo
  enddo
  else
   print*,'bad dir value in pert_mo_one_e_ints:',dir
   call abort
  endif

  call save_mo_one_e_ints(mo_one_e_integrals)

end

subroutine save_mo_one_e_ints(A)

  implicit none

  double precision, intent(in) :: A(mo_num,mo_num)

  call ezfio_set_mo_one_e_ints_mo_one_e_integrals(A)
  !call ezfio_set_mo_one_e_ints_io_mo_one_e_integrals('Read')
  call ezfio_set_mo_one_e_ints_io_mo_one_e_integrals('None')
  call ezfio_set_mo_one_e_ints_io_mo_integrals_kinetic('None')
  call ezfio_set_mo_one_e_ints_io_mo_integrals_pseudo('None')
  call ezfio_set_mo_one_e_ints_io_mo_integrals_n_e('None')

end

subroutine init_mo_one_e_ints

  implicit none

  PROVIDE ao_one_e_integrals

  ! AOs to MOs
  call ao_to_mo(ao_one_e_integrals,size(ao_one_e_integrals,1),mo_one_e_integrals,size(mo_one_e_integrals,1))

  call save_mo_one_e_ints(mo_one_e_integrals)

end
