subroutine diis_cc(err,all_err,t,all_t,sze,m,iter)

  implicit none

  ! https://hal.archives-ouvertes.fr/hal-02492983/document
  ! Maxime Chupin, Mi-Song Dupuy, Guillaume Legendre, Eric Séré. Convergence analysis of adaptive
  ! DIIS algorithms witerh application to electronic ground state calculations. 
  ! ESAIM: Mathematical Modelling and Numerical Analysis, EDP Sciences, 2021, 55 (6), pp.2785 - 2825. 10.1051/m2an/2021069ff.ffhal-02492983v5
  ! 
  ! t_{k+1} = g(t_k)
  ! err_k = f(t_k)
  !
  ! m_k = min(m,k)
  ! m maximal depth
  ! t_{k+1} = \sum_{i=0}^{m_k} c_i^k g(t_{k-m_k+i})
  ! \sum_{i=0}^{m_k} c_i^k = 1
  ! 
  ! b_{ij}^k = < err^{k-m_k+j}, err^{k-m_k+i} >
  ! 
  ! (b   -1) ( c^k    ) = (  0 )
  ! (-1   0) ( \lambda)   ( -1 ) 
  ! 
  ! In: t_0, m
  ! err_0 = g(t_0)
  ! k = 0
  ! m_k = 0
  ! while ||r_k|| > CC
  !   A.x=b
  !   t_{k+1} = \sum_{i=0}^{m_k} c_i^k g(t_{k-m_k+i})
  !   err_{k+1} = f(t_{k+1})
  !   m_{k+1} = min(m_k+1,m)
  !   k = k +1
  ! end

  ! {err_i}_{i=1}^{m_it} -> B -> c
  ! {t_i}_{i=1}^{m_it}, c, {err_i}_{i=1}^{m_it} -> t_{m_it+1}

  integer, intent(in)             :: m,iter,sze
  double precision, intent(in)    :: err(sze), all_err(sze,m)
  double precision, intent(in)    :: all_t(sze,m)
  double precision, intent(out)   :: t(sze)
  double precision, allocatable   :: B(:,:), c(:), zero(:)
  integer                         :: m_iter
  integer                         :: i,j,k
  integer                         :: info
  integer, allocatable            :: ipiv(:)
  
  m_iter = min(m,iter)
  print*,'m_iter',m_iter
  allocate(B(m_iter+1,m_iter+1), c(m_iter), zero(m_iter+1))
  allocate(ipiv(m+1))

  ! B(i,j) = < err(iter-m_iter+j),err(iter-m_iter+i) > ! iter-m_iter will be zero for us
  B = 0d0
  do j = 1, m_iter
    do i = 1, m_iter
      do k = 1, sze
        B(i,j) = B(i,j) + all_err(k,m+1-i) * all_err(k,m+1-j)
      enddo
    enddo
  enddo
  do i = 1, m_iter
    B(i,m_iter+1) = -1
  enddo
  do j = 1, m_iter
    B(m_iter+1,j) = -1
  enddo
  !print*,'B'
  !do i = 1, m_iter+1
  !  print*,B(i,:)
  !enddo

  ! (0 0 .... 0 -1)
  zero = 0d0
  zero(m_iter+1) = -1d0

  ! Solve B.c = zero
  call dgesv(m_iter+1, 1, B, size(B,1), ipiv, zero, size(zero,1), info)
  if (info /= 0) then
    print*,'Error dgesv:', info
    call abort
  endif
  c = zero(1:m_iter)
  !print*,'c',c
  !print*,'all_t' 
  !do i = 1, m
  !  print*,all_t(:,i)
  !enddo
  !print*,'all_r' 
  !do i = 1, m
  !  print*,all_err(:,i)
  !enddo

  ! update T 
  t(:) = 0d0
  do i = 1, m_iter
    t(:) = t(:) + c(i) * (all_t(:,m+1-i) + all_err(:,m+1-i))
  enddo

  !print*,'new t',t

  deallocate(ipiv,B,c,zero)

end

subroutine update_all_err(err,all_err,sze,m,iter)

  implicit none

  BEGIN_DOC
  ! Shift all the err vectors of the previous iterations to add the new one
  END_DOC

  integer, intent(in)             :: m, iter, sze
  double precision, intent(in)    :: err(sze)
  double precision, intent(inout) :: all_err(sze,m)
  integer                         :: i,j
  integer                         :: m_iter

  m_iter = min(m,iter)

  ! Shift
  do i = 1, m-1
    all_err(:,i) = all_err(:,i+1)
  enddo
  !print*,'shift err'
  !do i = 1, m
  !  print*,i, all_err(:,i)
  !enddo

  ! New
  all_err(:,m) = err(:)

  !print*,'Updated err'
  !do i = 1, m
  !  print*,i, all_err(:,i)
  !enddo

end

subroutine update_all_t(t,all_t,sze,m,iter)

  implicit none

  BEGIN_DOC
  ! Shift all the t vectors of the previous iterations to add the new one
  END_DOC

  integer, intent(in)             :: m, iter, sze
  double precision, intent(in)    :: t(sze)
  double precision, intent(inout) :: all_t(sze,m)
  integer                         :: i,j
  integer                         :: m_iter

  m_iter = min(m,iter)

  ! Shift
  do i = 1, m-1
    all_t(:,i) = all_t(:,i+1)
  enddo

  ! New
  all_t(:,m) = t(:)

  !print*,'Updated t'
  !do i = 1, m
  !  print*,i, all_t(:,i)
  !enddo

end

