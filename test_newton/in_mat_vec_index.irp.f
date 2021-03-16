subroutine in_mat_vec_index(i,p,q)
  implicit none
  integer,intent(in) :: i
  integer, intent(out) :: p,q
  integer :: a,b
  double precision :: da
  
  da = 0.5d0*(1+ sqrt(1d0+8d0*DBLE(i)))
  a = INT(da) 
  if ((a*(a-1))/2==i) then
    p = a-1
  else
    p = a
  endif
  b = p*(p-1)/2
 
  p = p + 1
  q = i - b 
end subroutine
