subroutine run_print_percentage_t

  implicit none

  integer :: i,s, nb_T
  integer, allocatable :: list_states(:)
  double precision, allocatable :: percentage(:,:), accu(:)
  character(len=2) :: exc

  nb_T = 3

  if (nb_T < 1) then
    print*,'nb_T should be >= 1, abort'
    call abort
  endif


  allocate(percentage(nb_T, n_states), accu(n_states), list_states(n_states))

  call percentage_t(nb_T,percentage)
  
  do s = 1, n_states
   list_states(s) = s
  enddo   

  print*,''
  print*,'Percentage of the excitations per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  if (percentage_in_exp) then
    do i = 1, nb_T
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%T', adjustl(exc), percentage(i,:)    
    enddo
  else
    do i = 1, nb_T
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%T', adjustl(exc), percentage(i,:)
    enddo
  endif

  print*,''
  print*,'Percentage of the excitations per state in'
  print*,'intermediate normalization, %T0=1:'
  write(*,'(A4,10(I12))') '', list_states(:)
  if (percentage_in_exp) then
    do i = 1, nb_T
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%T', adjustl(exc), percentage(i,1)/percentage(1,:)
    enddo
  else
    do i = 1, nb_T
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%T', adjustl(exc), percentage(i,1)/percentage(1,:)
    enddo
  endif

  print*,''
  print*,'Sum of the contributions per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  accu = 0d0
  if (percentage_in_exp) then
    do i = 1, nb_T
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%T', adjustl(exc), accu(:)
    enddo
  else
    do i = 1, nb_T
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%T', adjustl(exc), accu(:)
    enddo
  endif

  print*,''
  print*,'Missing contributions per state:'
  write(*,'(A4,10(I12))') '', list_states(:)
  accu = 0d0
  if (percentage_in_exp) then
    do i = 1, nb_T
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(1pE12.4))') '%T', adjustl(exc), 100d0-accu(:)        
    enddo
  else
    do i = 1, nb_T
      do s = 1, n_states
        accu(s) = accu(s) + percentage(i,s)
      enddo
      write (exc, '(I2)') i-1
      write (*, '(A2,A2,10(F12.4))') '%T', adjustl(exc), 100d0-accu(:)        
    enddo
  endif

  deallocate(percentage, accu, list_states)

end

subroutine percentage_t(nb_T,percentage)

  implicit none

  ! in
  integer :: nb_T

  ! out
  double precision, intent(out) :: percentage(nb_T, n_states) 

  ! internal
  integer :: i,j, s, n_occ, n_vir, nb_t1, i_int, idx_hf
  double precision, allocatable :: t1_amplitude(:,:,:), t2_amplitude(:,:,:,:,:)
  double precision, allocatable :: c1_coef(:,:)
  integer(bit_kind), allocatable :: c1_det(:,:,:)

  ! function
  double precision :: dnrm2

  n_occ = elec_alpha_num + elec_beta_num
  n_vir = 2*mo_num - n_occ

  allocate(t1_amplitude(n_occ,n_vir,N_states))
  allocate(t2_amplitude(n_occ,n_occ,n_vir,n_vir,N_states))

  percentage = 0d0

  ! %T_0 
  call find_det(HF_bitmask,idx_hf)
  do s = 1, n_states
    percentage(1,s) = psi_coef(idx_hf,s)**2
  enddo
  
  ! %T1
  ! t_i^a = c_i^a, %T1 = \sum t1^2
  if (nb_T >= 2) then
  call ci_coef_to_t1(n_occ, n_vir, nb_t1, t1_amplitude)
    do s = 1, N_states
      percentage(2,s) = dnrm2(n_occ * n_vir, t1_amplitude(1,1,s), 1)**2
    enddo
  endif
  
  allocate(c1_coef(nb_t1,N_states), c1_det(N_int,2,nb_t1))
  !call t1_to_c1(n_occ,n_vir,t1_amplitude,nb_t1,c1_coef,c1_det)

  ! %T2
  if (nb_T >= 3) then
    call ci_coef_to_t2(n_occ, n_vir, t1_amplitude, t2_amplitude)
    do s = 1, N_states
      percentage(3,s) = dnrm2(n_occ*n_occ * n_vir*n_vir, t2_amplitude(1,1,1,1,s), 1)**2
    enddo
  endif

  ! %T3
  !if (nb_T >= 4 .and. elec_alpha_num + elec_beta_num >= 3) then
  !  call ci_coef_to_t3(n_occ, n_vir, t1_amplitude, t2_amplitude, t3_amplitude)
  !  do s = 1, N_states
  !    percentage(4,s) = dnrm2(n_occ*n_occ*n_occ*n_vir*n_vir*n_vir, t3_amplitude(1,1,1,1,1,1,s), 1)**2
  !  enddo
  !endif

  ! Frac to %
  percentage = percentage * 100d0

  deallocate(t1_amplitude,t2_amplitude)
  deallocate(c1_coef,c1_det)

end

subroutine ci_coef_to_t1(n_occ, n_vir, nb_t1, t1_amplitude)

  implicit none

  ! in 
  integer, intent(in) :: n_occ, n_vir

  ! out 
  double precision, intent(out) :: t1_amplitude(n_occ,n_vir,N_states)
  integer, intent(out) :: nb_t1

  ! internal
  integer :: i,a,u,h1,p1,s1,h2,p2,s2,s
  integer :: degree, exc(0:2,2,2)
  double precision :: phase

  t1_amplitude = 0d0

  nb_t1 = 0
  do u = 1, N_det
    call get_excitation(HF_bitmask,psi_det(1,1,u),exc,degree,phase,N_int)
    if (degree == 1) then
      call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
      do s = 1, N_states
        ! alpha
        if (s1==1) then
          i = h1
          a = p1-elec_alpha_num
        ! beta
        else
          i = h1+elec_alpha_num
          a = p1+mo_num-2*elec_alpha_num
        endif
        if (i > 2*elec_alpha_num .or. a > 2*mo_num - 2*elec_alpha_num) then
          print*,'bug'
          print*,u
          print*,i,a
          print*,h1,h2,p1,p2,s1,s2
          call print_det(HF_bitmask,N_int)
          call print_det(psi_det(1,1,u),N_int)
          call abort
        endif
        t1_amplitude(i,a,s) = psi_coef(u,s)*phase 
      enddo

      ! Count
      do s = 1, N_states
        if (dabs(psi_coef(u,s)) > 1d-15) then
          nb_t1 = nb_t1 + 1
          exit
        endif
      enddo

    endif    
  enddo

end

!subroutine t1_to_c1(n_occ,n_vir,t1_amplitude,nb_c1,c1_coef,c1_det)
!
!  implicit none
!
!  ! in 
!  integer, intent(in) :: n_occ, n_vir, nb_c1
!  double precision, intent(in) :: t1_amplitude(n_occ,n_vir,N_states)
!
!  ! out 
!  double precision, intent(out) :: c1_coef(nb_c1,N_states)
!  integer(bit_kind), intent(out) :: c1_det(N_int,2,nb_c1)
!
!  ! internal 
!  integer :: i,j,a,b,u,sa,si,s
!  integer :: sigma,h1,p1,s1,h2,p2,s2
!  integer(bit_kind), allocatable :: res(:,:)
!  logical :: ok, non_zero
!
!  allocate(res(N_int,2))
!
!  print*,'nb_c1', nb_c1
!
!  do u = 1, nb_c1
!    c1_det(:,:,u) = psi_det(:,:,1)
!  enddo
!
!  u = 1
!  do a = 1, n_vir
!    if (a > mo_num - elec_alpha_num) then 
!      p1 = a - mo_num+2*elec_alpha_num
!      sa = 2
!    else 
!      p1 = a + elec_alpha_num
!      sa = 1
!    endif
!    do i = 1, n_occ
!      if (i > elec_alpha_num) then 
!        h1 = i - elec_alpha_num
!        si = 2
!      else
!        h1 = i
!        si = 1
!      endif
!
!      ! Non-existing amplitudes
!      if (sa /= si) then
!        cycle
!      endif
!
!      s1 = sa
!
!      do s = 1, N_states
!        non_zero = .False.          
!        if (dabs(t1_amplitude(i,a,s)) > 1d-15) then
!          non_zero = .True.  
!        endif
!        if (non_zero) then
!          exit
!        endif
!      enddo
!
!      if (.not. non_zero) then 
!        cycle
!      endif
!
!      print*,u
!      !print*,i,a
!      !print*,h1,p1,s1
!      call apply_hole(c1_det(1,1,u), s1, h1, res, ok, N_int)
!      !call print_det(res,N_int) 
!      c1_det(:,:,u) = res
!      call apply_particle(c1_det(1,1,u), s1, p1, res, ok, N_int)
!      c1_det(:,:,u) = res
!      call print_det(c1_det(1,1,u),N_int) 
!
!      ! c1 to t1
!      do s = 1, N_states 
!        c1_coef(u,s) = t1_amplitude(i,a,s)
!      enddo
!      print*,t1_amplitude(i,a,:)
!
!      u = u + 1
!    enddo
!  enddo  
!
!  deallocate(res)
!
!end

subroutine ci_coef_to_t2(n_occ, n_vir, t1_amplitude, t2_amplitude)

  implicit none

  ! in 
  integer, intent(in) :: n_occ, n_vir
  double precision, intent(in) :: t1_amplitude(n_occ,n_vir,N_states)

  ! out 
  double precision,intent(out) :: t2_amplitude(n_occ,n_occ,n_vir,n_vir,N_states)

  ! internal
  integer :: i,j,a,b,s,u
  integer :: h1,p1,s1,h2,p2,s2
  integer :: degree, exc(0:2,2,2)
  double precision :: phase, accu
  double precision, allocatable :: c2_coef(:,:,:,:,:)

  allocate(c2_coef(n_occ,n_occ,n_vir,n_vir,N_states))

  c2_coef = 0d0
  t2_amplitude = 0d0

  do u = 1, N_det
    call get_excitation(HF_bitmask,psi_det(1,1,u),exc,degree,phase,N_int)
    if (degree == 2) then
      call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
      do s = 1, N_states
        ! alpha alpha
        if (s1==1 .and. s2==1) then
          i = h1
          j = h2
          a = p1-elec_alpha_num
          b = p2-elec_alpha_num
        ! alpha beta
        elseif (s1==1 .and. s2==2) then
          i = h1
          j = h2+elec_alpha_num
          a = p1-elec_alpha_num
          b = p2+mo_num-2*elec_alpha_num
        ! beta alpha
        elseif (s1==2 .and. s2==1) then
          i = h1+elec_alpha_num
          j = h2
          a = p1+mo_num-2*elec_alpha_num
          b = p2-elec_alpha_num
        ! beta beta
        else
          i = h1+elec_alpha_num
          j = h2+elec_alpha_num
          a = p1+mo_num-2*elec_alpha_num
          b =  p2+mo_num-2*elec_alpha_num
        endif
        c2_coef(i,j,a,b,s) = psi_coef(u,s)*phase
        !print '(6(I3))', h1,h2,p1,p2,s1,s2
        !print '(4(I3))', i,j,a,b
      enddo 
    endif
  enddo

  !do s = 1, N_states
  !  print*,'%C2',dnrm2(n_occ*n_occ*n_vir*n_vir, c2_coef(1,1,1,1,s),1)**2 * 100d0
  !enddo

  ! c2 to t2
  accu = 0d0
  do s = 1, N_states
    do b = 1, n_vir
      do a = 1, n_vir
        do j = 1, n_occ
          do i = 1, n_occ
            t2_amplitude(i,j,a,b,s) = c2_coef(i,j,a,b,s) - t1_amplitude(i,a,s)*t1_amplitude(j,b,s) + t1_amplitude(i,b,s)*t1_amplitude(j,a,s)
          enddo
        enddo
      enddo
    enddo
  enddo
  !do s = 1, N_states
  !  print*,'%T2',dnrm2(n_occ*n_occ*n_vir*n_vir, t2_amplitude(1,1,1,1,s),1)**2 * 100d0
  !enddo

  deallocate(c2_coef)

end

!subroutine ci_coef_to_t3(n_occ, n_vir, t1_amplitude, t2_amplitude, t3_amplitude)
!
!  implicit none
!
!  ! in 
!  integer, intent(in) :: n_occ, n_vir
!  double precision, intent(in) :: t1_amplitude(n_occ,n_vir,N_states)
!  double precision,intent(in) :: t2_amplitude(n_occ,n_occ,n_vir,n_vir,N_states)
!
!  ! out 
!  double precision,intent(out) :: t3_amplitude(n_occ,n_occ,n_occ,n_vir,n_vir,n_vir,N_states)
!
!  ! internal
!  integer :: i,j,a,b,s,u
!  integer :: h1,p1,s1,h2,p2,s2
!  integer :: degree, exc(0:2,2,2)
!  double precision :: phase
!  double precision, allocatable :: c3_coef(:,:,:,:,:,:,:)
!
!  allocate(c3_coef(n_occ,n_occ,n_occ,n_vir,n_vir,n_vir,N_states))
!
!  c3_coef = 0d0
!  t3_amplitude = 0d0
!
!  do u = 1, N_det
!    call get_triple_excitation(HF_bitmask,psi_det(1,1,u),exc,degree,phase,N_int)
!    if (degree == 3) then
!      call decode_triple_exc(exc,degree,h1,p1,h2,p2,h3,p3,s1,s2,s3)
!      do s = 1, N_states
!        ! alpha alpha alpha
!        if (s1==1 .and. s2==1 .and. s3==1) then
!          i = h1
!          j = h2
!          k = h3
!          a = p1-elec_alpha_num
!          b = p2-elec_alpha_num
!          c = p3-elec_alpha_num
!        ! alpha alpha beta
!        elseif (s1==1 .and. s2==1 .and. s3==2) then
!          i = h1
!          j = h2
!          k = h3+elec_alpha_num
!          a = p1-elec_alpha_num
!          b = p2-elec_alpha_num
!          c = p3+mo_num-2*elec_alpha_num
!        ! alpha beta alpha
!        elseif (s1==2 .and. s2==2 .and. s3==1) then
!          i = h1
!          j = h2+elec_alpha_num
!          k = h3
!          a = p1-elec_alpha_num
!          b = p2+mo_num-2*elec_alpha_num
!          c = p3-elec_alpha_num
!        ! beta alpha alpha 
!        elseif (s1==2 .and. s2==1 .and. s3==1) then
!          i = h1+elec_alpha_num
!          j = h2
!          k = h3
!          a = p1+mo_num-2*elec_alpha_num
!          b = p2-elec_alpha_num
!          c = p3-elec_alpha_num
!        ! beta beta alpha 
!        elseif (s1==2 .and. s2==2 .and. s3==1) then
!          i = h1+elec_alpha_num
!          j = h2+elec_alpha_num
!          k = h3
!          a = p1+mo_num-2*elec_alpha_num
!          b = p2+mo_num-2*elec_alpha_num
!          c = p3-elec_alpha_num
!        ! beta alpha beta 
!        elseif (s1==2 .and. s2==1 .and. s3==2) then
!          i = h1+elec_alpha_num
!          j = h2
!          k = h3+elec_alpha_num
!          a = p1+mo_num-2*elec_alpha_num
!          b = p2-elec_alpha_num
!          c = p3+mo_num-2*elec_alpha_num
!        ! alpha beta beta
!        elseif (s1==1 .and. s2==2 .and. s3==2) then
!          i = h1
!          j = h2+elec_alpha_num
!          k = h3+elec_alpha_num
!          a = p1-elec_alpha_num
!          b = p2+mo_num-2*elec_alpha_num
!          c = p3+mo_num-2*elec_alpha_num
!        ! beta beta beta
!        else
!          i = h1+elec_alpha_num
!          j = h2+elec_alpha_num
!          k = h3+elec_alpha_num
!          a = p1+mo_num-2*elec_alpha_num
!          b = p2+mo_num-2*elec_alpha_num
!          c = p3+mo_num-2*elec_alpha_num
!        endif
!        c3_coef(i,j,k,a,b,c,s) = psi_coef(u,s)*phase
!      enddo
!    endif
!  enddo
!
!  ! c2 to t2
!  accu = 0d0
!  do s = 1, N_states
!    do c = 1, n_vir
!      do b = 1, n_vir
!        do a = 1, n_vir
!          do k = 1, n_occ
!            do j = 1, n_occ
!              do i = 1, n_occ
!                t3_amplitude(i,j,k,a,b,c,s) = c3(i,j,k,a,b,c,s) + &
!   t1(i,c,s) * t1(j,b,s) * t1(k,a,s) - t1(i,b,s) * t1(j,c,s) * t1(k,a,s) - & 
!   t1(i,c,s) * t1(j,a,s) * t1(k,b,s) + t1(i,a,s) * t1(j,c,s) * t1(k,b,s) + &
!   t1(i,b,s) * t1(j,a,s) * t1(k,c,s) - t1(i,a,s) * t1(j,b,s) * t1(k,c,s) - &
!   t1(k,c,s) * t2(i,j,a,b,s) + t1(k,b,s) t2(i,j,a,c,s) - &
!   t1(k,a,s) * t2(i,j,b,c,s) + t1(j,c,s) t2(i,k,a,b,s) - &
!   t1(j,b,s) * t2(i,k,a,c,s) + t1(j,a,s) t2(i,k,b,c,s) - &
!   t1(i,c,s) * t2(j,k,a,b,s) + t1(i,b,s) t2(j,k,a,c,s) - &
!   t1(i,a,s) * t2(j,k,b,c,s)
!              enddo
!            enddo
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!
!  deallocate(c3_coeff)
!
!end

subroutine oned_to_2d(k,dim_i,i,j)

  implicit none

  ! in
  integer, intent(in) :: k,dim_i

  ! out 
  integer, intent(out) :: i,j

  i = (k-1)/dim_i + i
  j = mod(k,dim_i+1) + 1

end

subroutine twod_to_1d(dim_i,i,j,k)

  implicit none

  ! in
  integer, intent(in) :: i,j,dim_i

  ! out
  integer, intent(out) :: k

  k = (j-1) * dim_i + i

end

