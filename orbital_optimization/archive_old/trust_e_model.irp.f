subroutine trust_e_model(n,v_grad,H,x,prev_energy,e_model)
   
  include 'constants.h' 

  implicit none
 
  !===================================================
  ! Compute the energy predicted by the Taylor series
  !===================================================

  integer, intent(in)          :: n
  double precision, intent(in) :: v_grad(n),H(n,n),x(n)
  double precision, intent(inout) :: prev_energy
  double precision, intent(out) :: e_model
  ! n      : integer, n = mo_num*(mo_num-1)/2
  ! v_grad : double precision vector of size n containing the gradient
  ! H      : n by n double precision matrix containing the hessian
  ! x      : double precision vector of size n containg the elements
  !          to compute the next step

  ! E(x) = E(0) + v_grad^T.x + 1/2 . x^T.H.x

  !double precision              :: e_model
  double precision              :: part_1, part_2
  double precision, allocatable :: part_2a(:)
  ! e_model     : double precision, predicted energy after the actual step
  ! prev_energy : double precision, energy after the previous step
  ! part_1      : double precision -> v_grad^T.x
  ! part_2      : double precision -> 1/2 . x^T.H.x 
  ! part_2a     : double precision vector of size n, temporary vector
  !               containing the result of H.x

  integer :: i,j
  ! i,j : integer, indexes

  !Function
  double precision              :: ddot
  ! ddot : double precision Blas function, dot product

  !============
  ! Allocation
  !============

  allocate(part_2a(n))

  !=============
  ! Calculation
  !=============

!  if (debug) then
    print*,''
    print*,'---Enter in trust_e_model---'
!  endif

  ! v_grad.x
  part_1 = ddot(n,v_grad,1,x,1)
 
  if (debug) then
    print*,'g^T.x : ', part_1
  endif  

  ! H.x
  call dgemv('N',n,n,1d0,H,size(H,1),x,1,0d0,part_2a,1)
  
  ! 1/2 . x^T.H.x
  part_2 = 0.5d0 * ddot(n,x,1,part_2a,1)
  
  if (debug) then
    print*,'1/2 x^T.H.x : ', part_2 
  endif

  ! Verif, pourquoi part_1 et 2 sont positifs ??
  ! Positif car le produit des MOs.R se fait dans dm_newton avec R^T
  ! au lieu de R ... donc ça change des signes 

  e_model = prev_energy - part_1 - part_2

  ! Writing the predicted energy
  print*, 'e_model for the next iteration : ', e_model

  !==============
  ! Deallocation
  !==============

  deallocate(part_2a)

  !if (debug) then
    print*,'---Leave trust_e_model---'
    print*,''
  !endif 
 
end subroutine 
