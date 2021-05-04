subroutine mat_to_vec_index(p,q,i)

  include 'constants.h'

  implicit none
  
  !==============================================
  ! Compute the index i of a vector element from
  ! the indexes p,q of a matrix element
  !
  ! Lower/upper diagonal matrix -> vector
  !==============================================

  !===========
  ! Variables
  !===========
  
  !====
  ! in
  !====
  integer, intent(in) :: p,q
  
  !===== 
  ! out
  !=====
  integer, intent(out) :: i 

  !==========
  ! internal
  !==========
  integer :: a,b
  double precision :: da

  !=============
  ! Calculation
  !=============
 
  if (debug) then
    print*, 'Enter in mat_to_vec_index'
  endif

  a = p-1
  b = a*(a-1)/2
  
  i = q+b

  if (debug) then
    print*, 'Leave mat_to_vec_index'
  endif
end subroutine
