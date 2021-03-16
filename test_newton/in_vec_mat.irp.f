subroutine in_vec_to_mat(p,q,i)
  implicit none
  
  integer, intent(in) :: p,q
  integer, intent(out) :: i 

  integer :: a,b
  double precision :: da

  a = p-1
  b = a*(a-1)/2
  
  i = q+b
end subroutine
