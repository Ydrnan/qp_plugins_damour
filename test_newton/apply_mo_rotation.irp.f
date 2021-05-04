subroutine apply_mo_rotation(R,prev_mos,new_mos)
  
  include 'constants.h'

  implicit none
  
  !=================================================================
  ! Compute the new MOs with the previous MOs and a rotation matrix
  !=================================================================

  !===========
  ! Variables
  !===========

  !===========
  ! intent in
  !===========
  double precision, intent(in)  :: R(mo_num,mo_num)
  ! R       : double precision mo_num by mo_num double precision matrix, rotation matrix

  ! out 
  double precision, intent(out) :: prev_mos(ao_num,mo_num)
  

  !==========
  ! internal
  !==========
  double precision,intent(in)  :: new_mos(ao_num,mo_num)
  integer                       :: i,j
  ! new_mos : ao_num by mo_num double precision matrix, new coefficients of the MOs 
  ! i       : integer, index 
 
  ! Provided 
  ! mo_num  : number of MOs
  ! ao_num  : number of AOs
  ! mo_coef : ao_num by mo_num double precision matrix, contains the coefficients of the MOs
 
  !============
  ! Allocation
  !============

  !=============
  ! Calculation
  !=============
 
  if (debug) then
    print*,'Enter in apply_mo_rotation'
  endif

  ! Product of old MOs (mo_coef) by Rotation matrix (R) 

  call dgemm('N','T',ao_num,mo_num,mo_num,1d0,mo_coef,size(mo_coef,1),R,size(R,1),0d0,new_mos,size(new_mos,1)) 
  
  !=========
  ! Storage
  !=========

  print*,'Save MOs...'
  
  prev_mos = mo_coef
  mo_coef = new_mos

  if (debug) then  
    print*,'New mo_coef : '
    do i=1,mo_num
      write(*,'(100(F10.5))') mo_coef(i,:)
    enddo
  endif

  ! Save the new MOs and change the label
  mo_label = 'MCSCF'
  call save_mos
  call ezfio_set_determinants_mo_label(mo_label)
  
  print*,'Done, MOs saved'
 
  !==============
  ! Deallocation
  !==============

  if (debug) then
    print*,'Leave apply_mo_rotation'
  endif

end subroutine
