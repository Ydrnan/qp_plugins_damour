subroutine dm_inversion(method,n,H,Hm1)
  
  include 'constants.h'
 
  implicit none
  
  !=====================================
  ! Matrix inversion by diagonalization 
  !=====================================
 
  !===========
  ! Variables
  !===========

  ! in
  integer, intent(in)           :: method
  integer, intent(in)           :: n
  double precision, intent(in)  :: H(n,n)
  ! method   : integer, full hessian -> 1, or diagonal hessain -> 2
  ! n        :  integer, n = mo_num*(mo_num-1)/2
  ! H        : n by n double precision matrix containing the hessian

  ! out
  double precision, intent(out) :: Hm1(n,n)
  ! Hm1      :  n by n double precision matrix, the inverse matrix of the hessian 

  ! internal
  double precision, allocatable :: e_val(:),W(:,:),Hm1_tmpr(:,:)
  double precision, allocatable :: work(:,:)
  integer                       :: info,lwork
  integer                       :: i
  double precision              :: t1,t2,t3
  ! e_val    : double precision vector containing the n eigenvalues of the H matrix
  ! W        : n by n double precision matrix containg the eigenvectors from the diagonalization of H
  ! Hm1_tmpr : n by n double precision temporary matrix to invert H 
  ! work     : lwork by n double precision working matrix for the diagonalization
  ! info     : integer  :
                ! if info = 0, the execution is successful
                ! if info = k, the k-th parameter has an illegal value
                ! if info = -k, the algorithm failed
  ! lwork    : integer for the working matrix work
  ! i        : integer, index for the matrices
  ! t1,t2,t3 : double precision, t3 = t2 - t1, time to compute the inverse matrix 
 
  CALL CPU_TIME(t1)

  if (method == 1) then

    ! For full hessian    

    !============
    ! Allocation
    !============
  
    lwork=3*n-1
  
    allocate(W(n,n),Hm1_tmpr(n,n))
    allocate(work(lwork,n),e_val(n))
  
    !=============
    ! Calculation
    !=============
  
    ! Diagonalization
    
    W = H
  
    call dsyev('V','U',n,W,size(W,1),e_val,work,lwork,info)
    
    if (info /= 0) then
        print*, 'Error diagonalization : dm_inversion'
        call ABORT
    endif
  
    Hm1=0d0
  
    if (debug) then
      print *, 'Eigenvalues of the n by n Hessian matrix (dm_inversion) :'
      write(*,'(100(F10.5))')  real(e_val(:))
    endif
  
    ! Inversion of the eigenvalues
    do i=1,n
  
        if ( (e_val(i)>0.d0)) then
              Hm1(i,i)=1d0/max(1.d-4,e_val(i))
        else
              Hm1(i,i)=0d0  !1d0/min(-1.d-4,e_val(i))
        endif
  
    enddo
  
    ! Back transformation to the initial basis : Hm1 = W.(1/e_val).W^T
    call dgemm('N','T',n,n,n,1d0,Hm1,size(Hm1,1),W,size(W,1),0d0,Hm1_tmpr,size(Hm1_tmpr,1))
    call dgemm('N','N',n,n,n,1d0,W,size(W,1),Hm1_tmpr,size(Hm1_tmpr,1),0d0,Hm1,size(Hm1,1))
   
    !==============
    ! Deallocation
    !==============
    
    deallocate(Hm1_tmpr)
    deallocate(work,e_val)

  else

    ! For diagonal Hessian
    Hm1=0d0
    do i=1,n
      if ( (H(i,i)>0.d0)) then    
       Hm1(i,i) = 1d0/max(1.d-4,H(i,i))
      else 
        Hm1(i,i) = 0d0
      endif
    enddo
  
  endif
   
  CALL CPU_TIME(t2)
  
  t2 = t2-t1
  print*, 'Time to invert the Hessian (dm_inversion) : ', t2  

end subroutine
