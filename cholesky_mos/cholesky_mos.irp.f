program cholesky_mos
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
  integer :: mu,nu,la,si, n_alpha
  integer :: i,j,k,l,p,q,r,s,alpha
  integer :: info
  double precision, allocatable :: A(:,:,:,:), tmp_A(:,:,:,:), B(:,:,:,:), Ld(:,:,:,:)
  integer :: nb_error
  double precision :: max_error

  allocate(A(ao_num, ao_num, ao_num, ao_num), tmp_A(ao_num,ao_num,ao_num,ao_num), B(ao_num, ao_num, ao_num, ao_num))
  allocate(Ld(ao_num, ao_num, ao_num, ao_num))

  PROVIDE ao_integrals_map

  do si = 1, ao_num
    do la = 1, ao_num
      do nu = 1, ao_num
        do mu = 1, ao_num
          !A(mu,nu,la,si) = get_ao_two_e_integral(mu,nu,la,si,ao_integrals_map)
          A(mu,nu,la,si) = get_ao_two_e_integral(mu,la,nu,si,ao_integrals_map)
        enddo
      enddo
    enddo
  enddo
  do si = 1, ao_num
    do la = 1, ao_num
      do nu = 1, ao_num
        do mu = 1, ao_num
          if (dabs(A(mu,nu,la,si) - A(la,si,mu,nu)) > 1d-16) then
            call abort
          endif
          if (dabs(A(mu,nu,la,si) - A(nu,mu,la,si)) > 1d-16) then
            call abort
          endif
          if (dabs(A(mu,nu,la,si) - A(mu,nu,si,la)) > 1d-16) then
            call abort
          endif
          if (dabs(A(mu,nu,la,si) - A(nu,mu,si,la)) > 1d-16) then
            call abort
          endif
        enddo
      enddo
    enddo
  enddo

  if (debug) then
    print*,'AOs two e ints'
    do nu = 1, ao_num
      do mu = 1, ao_num
        write(*,'(100F10.6)') A(mu,nu,:,:)
      enddo
    enddo 
    print*,''
  endif

  ! copy for cholesky decomposition
  Ld = A

  call dpotrf('L', ao_num*ao_num, Ld, size(Ld,1)*size(Ld,2), info)

  if (info /= 0) then
    print*,info
    print*,'Error in dpotrf'
    call abort
  endif

  if (debug) then
    print*,'Ld before removing upper diag elements'
    do nu = 1, ao_num
      do mu = 1, ao_num
        write(*,'(100F10.6)') Ld(mu,nu,:,:)
      enddo
    enddo
  endif

  ! Put to zero the upper diagonal elements of L
  do q = 2, ao_num*ao_num
    do p = 1, q-1
      nu = (p-1)/ao_num + 1 ! ao_num = n_mu
      mu = p - (nu-1) * ao_num
   
      si = (q-1)/ao_num + 1 ! ao_num = n_la
      la = q - (si-1) * ao_num
      Ld(mu,nu,la,si) = 0d0
    enddo
  enddo

  if (debug) then
    print*,'Ld'
    do nu = 1, ao_num
      do mu = 1, ao_num
        write(*,'(100F10.6)') Ld(mu,nu,:,:)
      enddo
    enddo
  endif

  ! L.L^T should be equal to the initial matrix 
  call dgemm('N','T',ao_num*ao_num, ao_num*ao_num,ao_num*ao_num, &
             1d0, Ld, size(Ld,1)*size(Ld,2), Ld, size(Ld,3)*size(Ld,4), &
             0d0, B, size(B,1)*size(B,2))

  if (debug) then
    print*,'Ld.Ld^T ?= A'
    do nu = 1, ao_num
      do mu = 1, ao_num
        write(*,'(100F10.6)') B(mu,nu,:,:)
      enddo
    enddo
  endif
 
  ! Check Ld.Ld^T = A
  print*, 'Check Ld.Ld^T = A'
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


  double precision, allocatable :: L_pq(:,:,:)
  double precision, allocatable :: L_pnu(:,:,:)
  double precision, allocatable :: A_mo(:,:,:,:)

  n_alpha = ao_num**2

  allocate(L_pq(mo_num,mo_num,n_alpha))
  allocate(L_pnu(mo_num,ao_num,n_alpha))
  allocate(A_mo(mo_num,mo_num,mo_num,mo_num))

  ! mu,nu,la,si -> mu,q,la,si = mu,q,alpha
  !L_pnu = 0d0
  !do alpha = 1, n_alpha
  !  si = (alpha-1)/ao_num + 1 ! ao_num = n_la 
  !  la = alpha - (si-1) * ao_num
  !  do nu = 1, ao_num
  !    do p = 1, mo_num
  !      do mu = 1, ao_num
  !        L_pnu(p,nu,alpha) = L_pnu(p,nu,alpha) + mo_coef(mu,p) * Ld(mu,nu,la,si)
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call dgemm('T','N', mo_num, ao_num*n_alpha, ao_num, &
             1d0, mo_coef, size(mo_coef,2), Ld, size(Ld,1), &
             0d0, L_pnu, size(L_pnu,1))
  
  ! mu,q,alpha -> p,q,alpha
  !L_pq = 0d0
  !do alpha = 1, n_alpha
  !  do q = 1, mo_num
  !    do p = 1, mo_num
  !      do nu = 1, ao_num
  !        L_pq(p,q,alpha) = L_pq(p,q,alpha) + mo_coef(nu,q) * L_pnu(p,nu,alpha)
  !      enddo
  !    enddo
  !  enddo
  !enddo

  do alpha = 1, n_alpha
    call dgemm('N','N',mo_num, ao_num, ao_num, &
               1d0, L_pnu(1,1,alpha), size(L_pnu,1), mo_coef, size(mo_coef,1), &
               0d0, L_pq(1,1,alpha), size(L_pq,1))
  enddo

  ! Transfo to MO ints
  !A_mo =0d0
  !do s = 1, mo_num
  !  do r = 1, mo_num
  !    do q = 1, mo_num 
  !      do p = 1, mo_num
  !        do alpha = 1, ao_num*ao_num
  !          A_mo(p,q,r,s) = A_mo(p,q,r,s) + L_pq(p,q,alpha) *  L_pq(r,s,alpha)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo

  call dgemm('N','T',mo_num*mo_num, mo_num*mo_num, n_alpha, &
             1d0, L_pq, size(L_pq,1) * size(L_pq,2), &
             L_pq, size(L_pq,1)*size(L_pq,2), 0d0, A_mo, size(A_mo,1)*size(A_mo,2))

  if (debug) then
    print*,'Int in MO basis'
    do j= 1, mo_num
      do i = 1, mo_num
        write(*,'(100F10.6)') A_mo(i,j,:,:)
      enddo
    enddo
  endif

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

  if (debug) then
    print*,'Expected result'
    do j= 1, mo_num
      do i = 1, mo_num
        write(*,'(100F10.6)') mo_ints(i,j,:,:)
      enddo
    enddo
  endif

  ! Check
  print*,'Result check'
  nb_error = 0
  max_error = 0d0
  do s = 1, mo_num
    do r = 1, mo_num
      do p = 1, mo_num
        do q = 1, mo_num
          if (dabs(mo_ints(p,q,r,s) - A_mo(p,q,r,s)) > 1d-14) then
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

  deallocate(A, tmp_A, B, Ld, A_mo, L_pq, mo_ints)  

end
