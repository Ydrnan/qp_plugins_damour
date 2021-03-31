subroutine in_mat_vec_index(i,p,q)

  include 'constants.h'

  implicit none

  !==================================================
  ! Compute the indexes p,q of a matrix element with
  ! the vector index i
  !
  ! Vector -> lower/upper diagonal matrix
  !==================================================

  !===========
  ! Variables
  !===========

  !====
  ! in
  !====
  integer,intent(in) :: i
  
  !=====
  ! out
  !=====
  integer, intent(out) :: p,q
  
  !==========
  ! internal 
  !==========
  integer :: a,b
  double precision :: da

  !==============
  ! Calculations
  !==============

  if (debug) then
    print*,'Enter in in_mat_vec_index'
  endif
  
  da = 0.5d0*(1+ sqrt(1d0+8d0*DBLE(i)))
  a = INT(da) 
  if ((a*(a-1))/2==i) then
    p = a-1
  else
    p = a
  endif
  b = p*(p-1)/2
 
  ! Matrix element indexes
  p = p + 1
  q = i - b 

  if (debug) then
    print*,'Leaves in_mat_vec_index'
  endif
end subroutine
