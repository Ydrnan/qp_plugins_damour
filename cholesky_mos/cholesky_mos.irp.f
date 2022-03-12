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
  integer :: mu,nu,la,si
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

  ! copy for cholesky decomposition
  Ld = A

  call dpotrf('L', ao_num*ao_num, Ld, size(Ld,1)*size(Ld,2), info)

  if (info /= 0) then
    print*,info
    print*,'Error in dpotrf'
    call abort
  endif

  print*,'Ld before removing upper diag elements'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') Ld(mu,nu,:,:)
    enddo
  enddo

  ! Put the lower diag elements in a tmp matrix
  !tmp_A = 0d0
  !do q = 1, ao_num*ao_num
  !  do p = q, ao_num*ao_num
  !    j = (p-1)/ao_num + 1 ! ao_num = n_i
  !    i = p - (j-1) * ao_num
  !  
  !    l = (q-1)/ao_num + 1 ! ao_num = n_k
  !    k = q - (l-1) * ao_num
  !    tmp_A(i,j,k,l) = A(i,j,k,l)
  !    print*, i,j,k,l, A(i,j,k,l)
  !  enddo
  !enddo
  !A = tmp_A

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

  print*,'Ld'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') Ld(mu,nu,:,:)
    enddo
  enddo

  ! L.L^T should be equal to the initial matrix 
  call dgemm('N','T',ao_num*ao_num, ao_num*ao_num,ao_num*ao_num, 1d0, Ld, size(Ld,3)*size(Ld,4), Ld, size(Ld,1)*size(Ld,2), 0d0, B, size(B,1)*size(B,2))

  print*,'Ld.Ld^T ?= A'
  do nu = 1, ao_num
    do mu = 1, ao_num
      write(*,'(100F10.6)') B(mu,nu,:,:)
    enddo
  enddo
 
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


  double precision, allocatable :: L_muq(:,:,:), L_pq(:,:,:), L_las(:,:,:), L_rs(:,:,:)
  double precision, allocatable :: A_mo(:,:,:,:)
  allocate(L_muq(ao_num,mo_num,ao_num*ao_num), L_pq(mo_num,mo_num,ao_num*ao_num))
  allocate(L_las(ao_num,mo_num,ao_num*ao_num), L_rs(mo_num,mo_num,ao_num*ao_num))
  allocate(A_mo(mo_num,mo_num,mo_num,mo_num))

  ! mu,nu,la,si -> mu,q,la,si = mu,q,alpha
  L_muq = 0d0
  do alpha = 1, ao_num*ao_num
    si = (alpha-1)/ao_num + 1 ! ao_num = n_la 
    la = alpha - (si-1) * ao_num
    do q = 1, mo_num
      do mu = 1, ao_num
        do nu = 1, ao_num
          L_muq(mu,q,alpha) = L_muq(mu,q,alpha) + mo_coef(nu,q) * Ld(mu,nu,la,si)
        enddo
      enddo
    enddo
  enddo

  ! mu,q,alpha -> p,q,alpha
  L_pq = 0d0
  do alpha = 1, ao_num*ao_num
    do q = 1, mo_num
      do p = 1, mo_num
        do mu = 1, ao_num
          L_pq(p,q,alpha) = L_pq(p,q,alpha) + mo_coef(mu,p) * L_muq(mu,q,alpha)
        enddo
      enddo
    enddo
  enddo

  ! (mu,nu,la,si)^T -> la,si,mu,nu -> la,s,mu,nu = la,s,alpha
  L_las = 0d0
  do alpha = 1, ao_num*ao_num
    nu = (alpha-1)/ao_num + 1 ! ao_num = n_la 
    mu = alpha - (nu-1) * ao_num
    do s = 1, mo_num
      do la = 1, ao_num
        do si = 1, ao_num
          L_las(la,s,alpha) = L_las(la,s,alpha) + mo_coef(si,s) * Ld(la,si,mu,nu)
        enddo
      enddo
    enddo
  enddo

  ! la,s,alpha -> r,s,alpha
  L_rs = 0d0
  do alpha = 1, ao_num*ao_num
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
  do alpha = 1, ao_num*ao_num
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
  A_mo =0d0
  do s = 1, mo_num
    do r = 1, mo_num
      do q = 1, mo_num 
        do p = 1, mo_num
          do alpha = 1, ao_num*ao_num
            A_mo(p,q,r,s) = A_mo(p,q,r,s) + L_pq(p,q,alpha) *  L_rs(r,s,alpha)
          enddo
        enddo
      enddo
    enddo
  enddo

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

  deallocate(A, tmp_A, B, Ld, A_mo, L_pq, L_rs, L_muq, L_las, mo_ints)  

end
