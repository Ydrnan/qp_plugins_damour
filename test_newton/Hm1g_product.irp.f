subroutine Hm1g_product(n,Hm1,v_grad,m_Hm1g,Hm1g)
  
  include 'constants.h'  
  
  !==============================
  ! Compute the product H^-1 . g
  ! 
  ! The result is a vector but we
  ! do the transformation into a
  ! matrix for the next step
  !==============================

  implicit none
 
  !===========
  ! Variables
  !=========== 

  !====
  ! in
  !====
  integer, intent(in)          :: n 
  double precision, intent(in) :: Hm1(n,n),v_grad(n)
  ! n      : integer, n = mo_num*(mo_num-1)/2
  ! Hm1    : n by n double precision matrix, the inverse of the hessian
  ! v_grad : double precision vector of length n containing the gradient
 
  !=====
  ! out
  !=====
  double precision, intent(out) :: m_Hm1g(mo_num,mo_num),Hm1g(n) 
  ! Hm1g   : double precision vector of length n, result of the product Hm1.g
  ! m_Hm1g : mo_num by mo_num double precision matrix builds from the Hm1g vector  

  integer                       :: p,q,k
  ! p,q,k  : integer, indexes

  ! degmv : Blas routine, matrix vector product

  !=============
  ! Calculation
  !=============
 
  if (debug) then
    print*,'Enter in Hm1g_product'
  endif
  
  ! Product Hm1^T.g
  call dgemv('T',n,n,1d0,Hm1,size(Hm1,1),v_grad,1,0d0,Hm1g,1)

  ! Debug 
  if (debug) then
    print*,'Vector of size n, Hm1.g :'
    write(*,'(100(F10.5))') Hm1g(:)
  endif

  ! Vector with n element -> mo_num by mo_num matrix (lower diagonal)
  do q=1,mo_num
    do p=1,mo_num
      if (p>q) then
        call mat_to_vec_index(p,q,k)
        m_Hm1g(p,q) = Hm1g(k)
      else
        m_Hm1g(p,q)=0d0
      endif
    enddo
  enddo

  ! Antisymmetrization of the previous matrix
  do q=1,mo_num
    do p=1,mo_num
      if (p<q) then
        m_Hm1g(p,q) = - m_Hm1g(q,p)
      endif
    enddo
  enddo

  ! Debug
  if (debug) then  
    print*,'mo_num by mo_num matrix from the vector Hm1.g :'
    do p=1,mo_num
      write(*,'(100(F10.5))') m_Hm1g(p,:)
    enddo
  endif

  if (debug) then
    print*,'Leave Hm1g_product'
  endif

end subroutine

