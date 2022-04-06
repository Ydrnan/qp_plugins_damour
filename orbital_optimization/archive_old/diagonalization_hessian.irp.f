subroutine diagonalization_hessian(n,H,e_val,w)

  include 'constants.h'

  implicit none

  !===========
  ! Variables
  !===========
 

  !====
  ! in
  !====

  integer, intent(in) :: n
  double precision, intent(in) :: H(n,n)

  !=====
  ! out
  !=====

  double precision, intent(out) :: e_val(n), w(n,n)

  !==========
  ! internal
  !==========

  double precision, allocatable :: work(:,:)
  integer                       :: info,lwork!, nb_iter
  integer                       :: i
  integer                       :: nb_negative_vp


  !============
  ! Allocation
  !============

  lwork=3*n-1

  allocate(work(lwork,n))

  !=============
  ! Calculation
  !=============

  ! Copy the hessian matrix, the eigenvectors will be store in W
  W=H

  ! Diagonalization of the hessian
  call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info)

  if (info /= 0) then
      print*, 'Error diagonalization : trust_region'
      call ABORT
  endif

  if (debug) then
    print *, 'vp Hess:'
    write(*,'(100(F10.5))')  real(e_val(:))
  endif

  ! Number of negative eigenvalues
  nb_negative_vp = 0
  do i = 1, n
    if (e_val(i) < -1d-12) then
      nb_negative_vp = nb_negative_vp + 1
      print*,'e_val < 0 :', e_val(i)
    endif
  enddo
  print*,'Number of negative eigenvalues :', nb_negative_vp

  !==============
  ! Deallocation
  !==============

  deallocate(work)

end subroutine
