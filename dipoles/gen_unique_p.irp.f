!subroutine pert_dipole_moments
!
!  use bitmasks
!
!  implicit none
!
!  integer :: n_P
!  integer(bit_kind), allocatable :: psi_P(:,:,:)
!  double precision, allocatable  :: c_P(:)
!  double precision :: dip1,dip1_x,dip1_y,dip1_z,dip,dip_x,dip_y,dip_z
!
!  n_P = 100
!
!  allocate(c_P(n_P),psi_P(N_int,2,n_P))
!
!  ! Generates the first order wf (only the mono excitations)
!  call gen_unique_p(n_P,c_P,psi_P)
!
!  ! Dipole moments
!  call comp_dip(dip,dip_x,dip_y,dip_z)
!
!  ! First order dipole moments
!  !call first_order_dip_w_dm(n_P,psi_P,c_P,dip1_x,dip1_y,dip1_z,dip1)
!  call first_order_dip(n_P,psi_P,c_P,dip1_x,dip1_y,dip1_z,dip1)
!
!  print*,'0-order:', dip
!  print*,'1-order correction:',dip1
!  print*,'Total:',dsqrt((dip_x + dip1_x)**2 + (dip_y + dip1_y)**2+ (dip_z + dip1_z)**2)
!
!  deallocate(c_P,psi_P)
!
!end

subroutine gen_unique_p

  use bitmasks

  implicit none

  integer                        :: n_P
  integer(bit_kind), allocatable :: psi_P(:,:,:)
  double precision, allocatable  :: c_P(:,:)

  integer(bit_kind), allocatable :: cp_psi_P(:,:,:)
  integer                        :: nb_generators
  integer(bit_kind)              :: det(N_int,2), res(N_int,2)
  integer                        :: degree, exc(0:2,2,2)
  integer                        :: i,j,l,k,r,s,u,istate,n_P_max
  integer                        :: sigma,h1,p1,s1,h2,p2,s2,p,q
  logical                        :: ok, in_wf
  double precision               :: phase, norm, percentage
  double precision               :: t1,t2
  double precision               :: dip1(N_states),dip1_x(N_states),dip1_y(N_states),dip1_z(N_states)
  double precision               :: dip(N_states),dip_x(N_states),dip_y(N_states),dip_z(N_states)
  double precision, allocatable :: gamma_1(:,:,:), gamma_1_test(:,:,:)

  allocate(gamma_1(mo_num,mo_num,N_states),gamma_1_test(mo_num,mo_num,N_states))

  ! Generation of all the perturbed determinants
  n_P_max = min(100,size(psi_det,3)*2*(elec_alpha_num*(mo_num-elec_alpha_num)))

  allocate(psi_P(1,2,n_P_max))

  n_P = 1
  percentage = 0d0

  ! Nb of generators for the mono excitations 
  i = 1
  norm = 0d0
  do while (norm <= threshold_generators .and. i < size(psi_det,3))
    norm = norm + psi_coef(i,1)**2
    i = i + 1
  enddo
  nb_generators = i-1
  print*,'nb_generators',nb_generators

  call wall_time(t1)
  do i = 1, nb_generators!N_det

    if (dble(i)/dble(N_det)*100d0 >= percentage) then
      call wall_time(t2)
      print*,int(percentage), t2-t1
      percentage = percentage + 10d0
      call wall_time(t1)
    endif

    ! For alpha or beta excitation
    do sigma = 1, 2
      ! r should be occupied
      do r = n_core_orb + 1, mo_num
        ! s should be virtual 
        do s = n_core_orb + 1, mo_num
          det(:,:) = psi_det(:,:,i)
          
          ! Create P
          s1 = sigma
          h1 = r
          p1 = s
          call apply_hole(det, s1, h1, res, ok, N_int)
          if (.not. ok) then
            cycle
          endif
          det = res
          call apply_particle(det, s1, p1, res, ok, N_int)
          if (.not. ok) then
            cycle
          endif
          det = res
  
          ! Check if P is already in the wf
          in_wf = .False.
          do j = 1, N_det
            call get_excitation_degree(det,psi_det(1,1,j),degree,N_int)
            if (degree == 0) then
              in_wf = .True.
              exit
            endif
          enddo

          ! Already in the wf
          if (in_wf) then
            cycle
          endif

          ! Check if P is already in the perturbed wf
          in_wf = .False.
          if (n_P > 1) then
            do l = 1, n_P
              call get_excitation_degree(det,psi_P(1,1,l),degree,N_int) 
              if (degree == 0) then
                in_wf = .True.
                exit
              endif
            enddo   
          endif

          ! Already in the perturbed wf
          if (in_wf) then
            cycle
          endif

          psi_P(:,:,n_P) = det(:,:)
          n_P = n_P + 1 
          
          ! Build a bigger array if necessary
          if (n_P > n_P_max) then
            n_P_max = n_P * 2
            allocate(cp_psi_P(1,2,n_P_max))

            cp_psi_P = 0
            do l = 1, n_P-1
              cp_psi_P(:,:,l) = psi_P(:,:,l)
            enddo

            deallocate(psi_P)
            allocate(psi_P(1,2,n_P_max))

            psi_P = cp_psi_P
            deallocate(cp_psi_P)
          endif

        enddo
      enddo
    enddo
  enddo

  n_P = n_P - 1
  allocate(cp_psi_P(1,2,n_P))
  cp_psi_P = 0
  do l = 1, n_P
    cp_psi_P(:,:,l) = psi_P(:,:,l)
  enddo
  deallocate(psi_P)
  allocate(psi_P(1,2,n_P))
  psi_P = cp_psi_P
  deallocate(cp_psi_P)

  print*,'Number of single excitations generated:',n_P

  allocate(c_P(n_P,N_states))

  call first_order_coef(n_P,psi_P,c_P)

  ! Dipole moments
  do istate = 1, N_states
    dip(istate) = multi_s_dipole_moment(istate,istate)
    dip_x(istate) = multi_s_x_dipole_moment(istate,istate)
    dip_y(istate) = multi_s_y_dipole_moment(istate,istate)
    dip_z(istate) = multi_s_z_dipole_moment(istate,istate)
  enddo
  print*,'0-order:', dip(:)

  ! First order dipole moments
  ! w dm
  call first_order_dm(n_P,psi_P,c_P,gamma_1)
  call first_order_dip_w_dm(gamma_1,dip1_x,dip1_y,dip1_z,dip1)
  print*,'w dm'
  print*,'1-order correction:',dip1(:)
  ! without dm
  call first_order_dip(n_P,psi_P,c_P,dip1_x,dip1_y,dip1_z,dip1)
  print*,'without dm'
  print*,'1-order correction:',dip1(:)

  ! check
  call test_first_order_dm(n_P,psi_P,c_P,gamma_1_test)
  call first_order_dip_w_dm(gamma_1_test,dip1_x,dip1_y,dip1_z,dip1)

  deallocate(c_P,psi_P,gamma_1)

end

subroutine first_order_coef(n_P,psi_P,c_P)

  implicit none
 
  integer :: n_P 
  double precision, intent(out) :: c_P(n_P,N_states)
  integer(bit_kind), intent(in) :: psi_P(1,2,n_P)
  integer :: l, istate
  double precision :: norm
  double precision :: hpp,ihpsi(N_states)
  
  c_P = 0d0
  do l = 1, n_P
    call i_H_j(psi_P(1,1,l),psi_P(1,1,l),N_int,hpp)
    call i_H_psi(psi_P(1,1,l),psi_det,psi_coef,N_int,N_det,N_det,N_states,ihpsi)
    do istate = 1, N_states
      c_P(l,istate) = ihpsi(istate) / (ci_energy(istate) - nuclear_repulsion - hpp)
    enddo
  enddo

  do istate = 1, N_states
    norm = 0d0
    do l = 1, n_P
      norm = norm + c_P(l,istate)**2
    enddo
    print*,'norm',istate,norm
  enddo

end 

subroutine first_order_dm(n_P,psi_P,c_P,gamma_1)

  implicit none

  integer, intent(in) :: n_P
  integer(bit_kind), intent(in) :: psi_P(N_int,2,n_P)
  double precision, intent(in) :: c_P(n_P,N_states)
  double precision, intent(out) :: gamma_1(mo_num,mo_num,N_states)
  integer :: i,l,p,q,h1,p1,h2,p2,s1,s2,degree, exc(0:2,2,2)
  double precision :: phase, t1, t2, percentage
  integer :: istate
  double precision, allocatable :: gamma_tilde_1(:,:,:)

  PROVIDE mo_dipole_x mo_dipole_y mo_dipole_z

  ! first order 1-dm
  allocate(gamma_tilde_1(mo_num,mo_num,N_states))

  gamma_tilde_1 = 0d0

  percentage = 0d0
  call wall_time(t1)
  do l = 1, n_P
    if (dble(l)/dble(n_P)*100d0 >= percentage) then
      call wall_time(t2)
      print*,int(percentage), t2-t1
      percentage = percentage + 10d0
      call wall_time(t1)
    endif
    do i = 1, N_det
      call get_excitation(psi_P(1,1,l),psi_det(1,1,i),exc,degree,phase,N_int)
      if (degree == 1) then
        call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
        p = h1
        q = p1
        do istate = 1, N_states
          gamma_tilde_1(p,q,istate) = gamma_tilde_1(p,q,istate) + c_P(l,istate) * psi_coef(i,istate) * phase
        enddo
      endif
    enddo
  enddo

  do istate = 1, N_states
    do q = 1, mo_num
      do p = 1, mo_num
        gamma_1(p,q,istate) = gamma_tilde_1(p,q,istate) + gamma_tilde_1(q,p,istate)
      enddo
    enddo
    !print*,gamma_1(:,:,istate)
  enddo

  deallocate(gamma_tilde_1)

end

subroutine test_first_order_dm(n_P,psi_P,c_P,gamma_1)

  implicit none

  integer, intent(in) :: n_P
  integer(bit_kind), intent(in) :: psi_P(N_int,2,n_P)
  double precision, intent(in) :: c_P(n_P,N_states)
  double precision, intent(out) :: gamma_1(mo_num,mo_num,N_states)
  integer :: i,l,p,q,h1,p1,h2,p2,s1,s2,degree, exc(0:2,2,2)
  double precision :: phase, t1, t2, percentage
  integer :: istate
  double precision, allocatable :: gamma_tilde_1(:,:,:), i_aa_psi_array(:,:,:)

  PROVIDE mo_dipole_x mo_dipole_y mo_dipole_z

  ! first order 1-dm
  allocate(gamma_tilde_1(mo_num,mo_num,N_states))
  allocate(i_aa_psi_array(mo_num,mo_num,N_states))

  gamma_tilde_1 = 0d0

  percentage = 0d0
  do l = 1, n_P
    call i_op_aa_psi(psi_P(1,1,l),psi_det,psi_coef,N_int,N_det,N_det,N_states,i_aa_psi_array)
    do istate = 1, N_states
      gamma_tilde_1(:,:,istate) = gamma_tilde_1(:,:,istate) + i_aa_psi_array(:,:,istate) * c_P(l,istate)
    enddo
  enddo

  do istate = 1, N_states
    do q = 1, mo_num
      do p = 1, mo_num
        gamma_1(p,q,istate) = gamma_tilde_1(p,q,istate) + gamma_tilde_1(q,p,istate)
      enddo
    enddo
  enddo

  deallocate(gamma_tilde_1)

end

subroutine first_order_dip_w_dm(gamma_1,dip1_x,dip1_y,dip1_z,dip1)

  implicit none

  double precision, intent(in) :: gamma_1(mo_num,mo_num,N_states)
  double precision, intent(out) :: dip1_x(N_states),dip1_y(N_states),dip1_z(N_states),dip1(N_states)
  integer :: i,l,p,q,h1,p1,h2,p2,s1,s2,degree, exc(0:2,2,2)
  double precision :: phase, t1, t2, percentage
  integer :: ii,istate

  PROVIDE mo_dipole_x mo_dipole_y mo_dipole_z

  dip1_x = 0d0
  dip1_y = 0d0
  dip1_z = 0d0
  do istate = 1, N_states
    do q = 1, mo_num
      do p = 1, mo_num
        dip1_x(istate) =  dip1_x(istate) - gamma_1(p,q,istate) * mo_dipole_x(p,q)
        dip1_y(istate) =  dip1_y(istate) - gamma_1(p,q,istate) * mo_dipole_y(p,q)
        dip1_z(istate) =  dip1_z(istate) - gamma_1(p,q,istate) * mo_dipole_z(p,q)
      enddo
    enddo
  enddo

  print*,'Dip1_x', dip1_x(:) * au2D
  print*,'Dip1_y', dip1_y(:) * au2D
  print*,'Dip1_z', dip1_z(:) * au2D
 
  do istate = 1, N_states
    dip1(istate) = dsqrt(dip1_x(istate)**2 + dip1_y(istate)**2 + dip1_z(istate)**2) * au2D
  enddo 

end

subroutine first_order_dip(n_P,psi_P,c_P,dip1_x,dip1_y,dip1_z,dip1)

  implicit none

  integer, intent(in) :: n_P
  integer(bit_kind), intent(in) :: psi_P(1,2,n_P)
  double precision, intent(in) :: c_P(n_P,N_states)
  double precision, intent(out) :: dip1_x(N_states),dip1_y(N_states),dip1_z(N_states),dip1(N_states)

  double precision :: phase,t1,t2,percentage
  integer :: k,l,i,p,q,h1,p1,h2,p2,s1,s2,degree,ii, exc(0:2,2,2), istate

  PROVIDE mo_dipole_x mo_dipole_y mo_dipole_z

  ! Without density matrix
  dip1_x = 0d0
  dip1_y = 0d0
  dip1_z = 0d0
  percentage = 0d0

  call wall_time(t1)
  do l = 1, n_P
    ii = l
    if (dble(ii)/dble(n_P)*100d0 >= percentage) then
      call wall_time(t2)
      print*,int(percentage), t2-t1
      percentage = percentage + 10d0
      call wall_time(t1)
    endif
    do i = 1, N_det
      call get_excitation(psi_P(1,1,l),psi_det(1,1,i),exc,degree,phase,N_int)
      if (degree == 1) then
        call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
        p = h1
        q = p1
        do istate = 1, N_states
          dip1_x(istate) = dip1_x(istate) - c_P(l,istate) * psi_coef(i,istate) * phase * (mo_dipole_x(p,q) + mo_dipole_x(q,p))
          dip1_y(istate) = dip1_y(istate) - c_P(l,istate) * psi_coef(i,istate) * phase * (mo_dipole_y(p,q) + mo_dipole_y(q,p))
          dip1_z(istate) = dip1_z(istate) - c_P(l,istate) * psi_coef(i,istate) * phase * (mo_dipole_z(p,q) + mo_dipole_z(q,p))
        enddo
      endif
    enddo
  enddo
 
  ! Error, the terms < P_k | r | P_k > are missing

  !call wall_time(t1)
  !double precision :: dip2_z
  !integer :: list_mo_a(:), list_mo_b(:)
  !allocate(list_mo_a(elec_alpha_num),list_mo_b(elec_beta_num)) 
  
  !dip2_z = 0d0
  !percentage = 0d0
  !do k = 1, n_P
  !   call bitstring_to_list(psi_P(1,1,k), list_mo_a, elec_alpha_num, N_int)
  !   call bitstring_to_list(psi_P(1,2,k), list_mo_b, elec_beta_num, N_int)
  !   do u = 1, elec_alpha_num
  !      p = list_mo_a(u) 
  !     dip2_z = dip2_z - c_P(k) * c_P(k) * mo_dipole_z(p,p)
  !   enddo
  !   do u = 1, elec_beta_num
  !      p = list_mo_b(u) 
  !     dip2_z = dip2_z - c_P(k) * c_P(k) * mo_dipole_z(p,p)
  !   enddo
  !  do l = k+1, n_P
  !    ii = k
  !    if (dble(ii)/dble(n_P)*100d0 >= percentage) then
  !      call wall_time(t2)
  !      print*,int(percentage), t2-t1
  !      percentage = percentage + 10d0
  !      call wall_time(t1)
  !    endif
  !    call get_excitation(psi_P(1,1,k),psi_P(1,1,l),exc,degree,phase,N_int)
  !    if (degree == 1) then
  !      call decode_exc(exc,degree,h1,p1,h2,p2,s1,s2)
  !      p = h1
  !      q = p1
  !      dip2_z = dip2_z - c_P(k) * c_P(l) * phase * (mo_dipole_z(p,q) + mo_dipole_z(q,p))
  !    endif
  !  enddo
  !enddo 

  print*,'Using gen_unique_p:'
  print*,'Dip1_x', dip1_x(:) * au2D
  print*,'Dip1_y', dip1_y(:) * au2D
  print*,'Dip1_z', dip1_z(:) * au2D
  !print*,'Dip2_z', dip2_z * au2D

  do istate = 1, N_states
    dip1(istate) = dsqrt(dip1_x(istate)**2 + dip1_y(istate)**2 + dip1_z(istate)**2) * au2D
  enddo

end
