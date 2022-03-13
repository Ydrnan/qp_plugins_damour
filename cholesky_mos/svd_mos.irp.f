program svd_mos
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  END_DOC
  print *, 'Hello world'
  call ao_decomposition()
end

subroutine ao_decomposition()

  implicit none

  double precision :: get_ao_two_e_integral
  integer :: mu,nu,la,si
  integer :: i,j,k,l,p,q,r,s,alpha
  integer :: info
  double precision, allocatable :: A(:,:,:,:), tmp_A(:,:,:,:), B(:,:,:,:), U(:,:,:,:), Vt(:,:,:,:), e(:), Ld(:,:,:,:)
  double precision, allocatable :: L_svd(:,:,:), tmp_U(:,:,:), diag_e(:,:)
  integer :: nb_error
  double precision :: max_error
  integer :: lwork
  double precision, allocatable :: work(:)
  double precision :: thresh
  integer :: n_e

  lwork = 1
  
  allocate(A(ao_num, ao_num, ao_num, ao_num), tmp_A(ao_num,ao_num,ao_num,ao_num), B(ao_num, ao_num, ao_num, ao_num))
  allocate(Ld(ao_num, ao_num, ao_num, ao_num))
  allocate(U(ao_num, ao_num, ao_num, ao_num), Vt(ao_num, ao_num, ao_num, ao_num), e(ao_num*ao_num))
  allocate(work(lwork))  

  PROVIDE ao_integrals_map

  do si = 1, ao_num
    do la = 1, ao_num
      do nu = 1, ao_num
        do mu = 1, ao_num
          A(mu,nu,la,si) = get_ao_two_e_integral(mu,nu,la,si,ao_integrals_map)
        enddo
      enddo
    enddo
  enddo

  print*,'AOs two e ints'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') A(mu,nu,:,:)
    enddo
  enddo 
  print*,''

  tmp_A = A

  ! optimal work size
  lwork = -1
  work(1) = 0d0
  call dgesvd('A', 'A', ao_num*ao_num, ao_num*ao_num, &
              tmp_A, size(tmp_A,1)*size(tmp_A,2), e, U, size(U,1)*size(U,2), &
              Vt, size(Vt,1)*size(Vt,2), work, lwork, info)

  !print*,work(1)
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))

  ! svd
  call dgesvd('A', 'A', ao_num*ao_num, ao_num*ao_num, & 
              tmp_A, size(tmp_A,1)*size(tmp_A,2), e, U, size(U,1)*size(U,2), & 
              Vt, size(Vt,1)*size(Vt,2), work, lwork, info)


  if (info /= 0) then
    print*,info
    print*,'Error in dgesvd'
    call abort
  endif

  print*,'U before removing singular vectors'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') U(mu,nu,:,:)
    enddo
  enddo

  print*,'e'
  print*,e(:)
  
  thresh = 0d0
  n_e = 0d0
  do i = 1, ao_num*ao_num
    if (e(i) > thresh) then
      e(i) = dsqrt(e(i))
      n_e = n_e + 1
    else
      e(i) = 0d0
    endif
  enddo

  allocate(tmp_U(ao_num,ao_num,n_e), diag_e(n_e,n_e), L_svd(ao_num,ao_num,n_e))

  ! diag_e
  diag_e = 0d0
  do alpha = 1, n_e
    diag_e(alpha,alpha) = e(alpha)
  enddo

  ! tmp_U
  do alpha = 1, n_e
    si = (alpha-1)/ao_num + 1 ! ao_num = n_la
    la = alpha - (si-1) * ao_num
    do nu = 1, ao_num
      do mu = 1, ao_num
        tmp_U(mu,nu,alpha) =  U(mu,nu,la,si)
      enddo
    enddo
  enddo

  ! L_svd = tmp_U . diag_e
  call dgemm('N','N', ao_num*ao_num, n_e, n_e, &
             1d0, tmp_U, size(tmp_U,1) * size(tmp_U,2), diag_e, size(diag_e,1), &
             0d0, L_svd, size(L_svd,1) * size(L_svd,2))

  print*,'L_svd'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') L_svd(mu,nu,:)
    enddo
  enddo

  ! L.L^T should be equal to the initial matrix 
  call dgemm('N','T',ao_num*ao_num, ao_num*ao_num, n_e, &
             1d0, L_svd, size(L_svd,1)*size(L_svd,2), L_svd, size(L_svd,3), &
             0d0, B, size(B,1)*size(B,2))

  print*,'L_svd.L_svd^T ?= A'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') B(mu,nu,:,:)
    enddo
  enddo

  ! Check Ld.Ld^T = A
  print*, 'Check L_svd.L_svd^T = A'
  nb_error = 0
  max_error = 0d0
  do si = 1, ao_num
    do la = 1, ao_num
      do nu = 1, ao_num
        do mu = 1, ao_num
          if (dabs(B(mu,nu,la,si) - A(mu,nu,la,si)) > 1d-12) then
            nb_error = nb_error + 1
            if (dabs(B(mu,nu,la,si) - A(mu,nu,la,si)) > max_error) then
              max_error = dabs(B(mu,nu,la,si) - A(mu,nu,la,si))
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  print*, 'Nb error:', nb_error
  print*, 'Max_error:', max_error

  double precision, allocatable :: L_pnu(:,:,:), L_pq(:,:,:), L_las(:,:,:), L_rs(:,:,:)
  double precision, allocatable :: A_mo(:,:,:,:)
  allocate(L_pnu(mo_num,ao_num,ao_num*ao_num), L_pq(mo_num,mo_num,ao_num*ao_num))
  allocate(L_las(ao_num,mo_num,ao_num*ao_num), L_rs(mo_num,mo_num,ao_num*ao_num))
  allocate(A_mo(mo_num,mo_num,mo_num,mo_num))

  ! mu,nu,la,si -> mu,q,la,si = mu,q,alpha
  L_pnu = 0d0
  do alpha = 1, n_e
    si = (alpha-1)/ao_num + 1 ! ao_num = n_la 
    la = alpha - (si-1) * ao_num
    do nu = 1, ao_num
      do p = 1, mo_num
        do mu = 1, ao_num
          L_pnu(p,nu,alpha) = L_pnu(p,nu,alpha) + mo_coef(mu,p) * Ld(mu,nu,la,si)
        enddo
      enddo
    enddo
  enddo
  
  call dgemm('T','N', mo_num, ao_num*n_e, ao_num, &
             1d0, mo_coef, size(mo_coef,2), L_svd, size(L_svd,1), &
             0d0, L_pnu, size(L_pnu,1))
  
  ! mu,q,alpha -> p,q,alpha
  L_pq = 0d0
  do alpha = 1, ao_num*ao_num
    do q = 1, mo_num
      do p = 1, mo_num
        do nu = 1, ao_num
          L_pq(p,q,alpha) = L_pq(p,q,alpha) + mo_coef(nu,q) * L_pnu(p,nu,alpha)
        enddo
      enddo
    enddo
  enddo

  ! Not working
  !do p = 1, mo_num
  !  call dgemm('T','N',mo_num, ao_num*ao_num, ao_num, &
  !             1d0, mo_coef, size(mo_coef,2), L_pnu(p,1,1), size(L_pnu,2), &
  !             0d0, L_pq(p,1,1), size(L_pq,2))
  !enddo

  ! ### Just to check ###
  ! (mu,nu,la,si)^T -> la,si,mu,nu -> la,s,mu,nu = la,s,alpha
  L_las = 0d0
  do alpha = 1, n_e
    nu = (alpha-1)/ao_num + 1 ! ao_num = n_la 
    mu = alpha - (nu-1) * ao_num
    do s = 1, mo_num
      do la = 1, ao_num
        do si = 1, ao_num
          L_las(la,s,alpha) = L_las(la,s,alpha) + mo_coef(si,s) * L_svd(la,si,alpha)
        enddo
      enddo
    enddo
  enddo

  ! la,s,alpha -> r,s,alpha
  L_rs = 0d0
  do alpha = 1, n_e
    do s = 1, mo_num
      do r = 1, mo_num
        do la = 1, ao_num
          L_rs(r,s,alpha) = L_rs(r,s,alpha) + mo_coef(la,r) * L_las(la,s,alpha)
        enddo
      enddo
    enddo
  enddo

  ! L-pq and L_rs should be equals
  print*,'L_pq'
  do j= 1, mo_num
    do i = 1, mo_num
      write(*,'(100F14.6)') L_pq(i,j,:)
    enddo
  enddo
 
  print*,'L_rs'
  do j= 1, mo_num
    do i = 1, mo_num
      write(*,'(100F14.6)') L_rs(i,j,:)
    enddo
  enddo
 
  ! Check
  print*,'Check L_pq = L_rs'
  nb_error = 0
  max_error = 0d0
  do alpha = 1, n_e
    do j = 1, mo_num
      do i = 1, mo_num
        if (dabs(L_pq(i,j,alpha) - L_rs(i,j,alpha)) > 1d-12) then
          nb_error = nb_error + 1
          if (dabs(L_pq(i,j,alpha) - L_rs(i,j,alpha)) > max_error) then
            max_error = dabs(L_pq(i,j,alpha) - L_rs(i,j,alpha))
          endif
        endif
      enddo
    enddo
  enddo 
  print*,'nb_error:',nb_error
  print*,'max error:', max_error

  ! Transfo to MO ints
  !A_mo =0d0
  !do s = 1, mo_num
  !  do r = 1, mo_num
  !    do q = 1, mo_num 
  !      do p = 1, mo_num
  !        do alpha = 1, ao_num*ao_num
  !          A_mo(p,q,r,s) = A_mo(p,q,r,s) + L_pq(p,q,alpha) *  L_rs(r,s,alpha)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','T',mo_num*mo_num, mo_num*mo_num, n_e, &
             1d0, L_pq, size(L_pq,1) * size(L_pq,2), &
             L_pq, size(L_pq,3), 0d0, A_mo, size(A_mo,1)*size(A_mo,2))

  print*,'Int in MO basis'
  do j= 1, mo_num
    do i = 1, mo_num
      write(*,'(100F14.6)') A_mo(i,j,:,:)
    enddo
  enddo

  PROVIDE mo_integrals_map

  double precision :: get_two_e_integral
  double precision, allocatable :: mo_ints(:,:,:,:)
  allocate(mo_ints(mo_num,mo_num,mo_num,mo_num))

  do s = 1, mo_num
    do r = 1, mo_num
      do q = 1, mo_num
        do p = 1, mo_num
          mo_ints(p,q,r,s) = get_two_e_integral(p,q,r,s,mo_integrals_map)
        enddo
      enddo
    enddo
  enddo

  print*,'Expected result'
  do j= 1, mo_num
    do i = 1, mo_num
      write(*,'(100F14.6)') mo_ints(i,j,:,:)
    enddo
  enddo

  ! Check
  print*,'Result check'
  nb_error = 0
  max_error = 0d0
  do s = 1, mo_num
    do r = 1, mo_num
      do p = 1, mo_num
        do q = 1, mo_num
          if (dabs(mo_ints(p,q,r,s) - A_mo(p,q,r,s)) > 1d-12) then
            nb_error = nb_error + 1
            if (dabs(mo_ints(p,q,r,s) - A_mo(p,q,r,s)) > max_error) then
              max_error = dabs(mo_ints(p,q,r,s) - A_mo(p,q,r,s))
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  print*, 'Nb error:', nb_error
  print*, 'Max_error:', max_error

  deallocate(A, tmp_A, B, Ld, A_mo, L_pq, L_rs, L_pnu, L_las, mo_ints)  

end
