subroutine i_r_j(key_i,key_j,Nint,i_x_j,i_y_j,i_z_j)

  use bitmasks

  implicit none

  BEGIN_DOC
  ! Returns $\langle I| x |J \rangle$, $\langle I| y |J \rangle$ and $\langle I| z |J \rangle$,
  ! where $I$ and $J$ are determinants.
  END_DOC

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: i_x_j,i_y_j,i_z_j

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  integer                        :: m,n,p,q,a,b
  integer                        :: i,j,k
  !integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: occ_a(Nint*bit_kind_size)
  integer                        :: occ_b(Nint*bit_kind_size)
  !double precision               :: diag_H_mat_elem
  double precision               :: phase
  !integer                        :: n_occ_ab(2)

  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  i_x_j = 0d0
  i_y_j = 0d0
  i_z_j = 0d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    ! < I | r | J >
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      if (exc(0,1,1) == 1) then
        ! Single alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      else
        ! Single beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      i_x_j = (mo_dipole_x(m,p) + mo_dipole_x(p,m)) * phase
      i_y_j = (mo_dipole_y(m,p) + mo_dipole_y(p,m)) * phase
      i_z_j = (mo_dipole_z(m,p) + mo_dipole_z(p,m)) * phase

    ! < I | r | I >
    case (0)
      ! alpha part
      call bitstring_to_list(key_i, occ_a, elec_alpha_num, Nint)
      do a = 1, elec_alpha_num
        p = occ_a(a)
        i_x_j = i_x_j + mo_dipole_x(p,p)
        i_y_j = i_y_j + mo_dipole_y(p,p)
        i_z_j = i_z_j + mo_dipole_z(p,p)
      enddo
      ! beta part
      call bitstring_to_list(key_i, occ_b, elec_beta_num, Nint)
      do b = 1, elec_beta_num
        p = occ_b(b)
        i_x_j = i_x_j + mo_dipole_x(p,p)
        i_y_j = i_y_j + mo_dipole_y(p,p)
        i_z_j = i_z_j + mo_dipole_z(p,p)
      enddo
  end select
end

subroutine i_r_psi(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_x_psi_array, i_y_psi_array, i_z_psi_array)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle I|x|Psi \rangle  = \sum_J c_J \langle I | x | J \rangle$.
  ! Computes $\langle I|y|Psi \rangle  = \sum_J c_J \langle I | y | J \rangle$.
  ! Computes $\langle I|z|Psi \rangle  = \sum_J c_J \langle I | z | J \rangle$.
  !
  ! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|I \rangle$
  ! is connected.
  ! The i_H_psi_minilist is much faster but requires to build the
  ! minilists.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_x_psi_array(Nstate), i_y_psi_array(Nstate), i_z_psi_array(Nstate)

  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: i_x_j,i_y_j,i_z_j
  integer, allocatable           :: idx(:)

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))

  i_x_psi_array = 0.d0
  i_y_psi_array = 0.d0
  i_z_psi_array = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      !call i_H_j(keys(1,1,i),key,Nint,hij)
      !i_H_psi_array(1) = i_H_psi_array(1) + coef(i,1)*hij
      call i_r_j(keys(1,1,i),key,Nint,i_x_j,i_y_j,i_z_j)
      i_x_psi_array(1) = i_x_psi_array(1) + i_x_j * coef(i,1)
      i_y_psi_array(1) = i_y_psi_array(1) + i_y_j * coef(i,1)
      i_z_psi_array(1) = i_z_psi_array(1) + i_z_j * coef(i,1)
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      !call i_H_j(keys(1,1,i),key,Nint,hij)
      call i_r_j(keys(1,1,i),key,Nint,i_x_j,i_y_j,i_z_j)
      do j = 1, Nstate
        !i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
        i_x_psi_array(j) = i_x_psi_array(j) + i_x_j * coef(i,j)
        i_y_psi_array(j) = i_y_psi_array(j) + i_y_j * coef(i,j)
        i_z_psi_array(j) = i_z_psi_array(j) + i_z_j * coef(i,j)
      enddo
    enddo

  endif

end

subroutine i_monoe_op_j(key_i,key_j,Nint,monoe_ints,i_op_j)

  use bitmasks

  implicit none

  BEGIN_DOC
  ! Returns $\langle I| a mono_electronic operator |J \rangle$ where $I$ and $J$ are determinants.
  END_DOC

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(in)   :: monoe_ints(mo_num,mo_num)
  double precision, intent(out)  :: i_op_j

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  integer                        :: m,n,p,q,a,b
  integer                        :: i,j,k
  !integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: occ_a(Nint*bit_kind_size)
  integer                        :: occ_b(Nint*bit_kind_size)
  !double precision               :: diag_H_mat_elem
  double precision               :: phase
  !integer                        :: n_occ_ab(2)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)


  i_op_j = 0d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    ! < I | op | J >
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      if (exc(0,1,1) == 1) then
        ! Single alpha
        m = exc(1,1,1)
        p = exc(1,2,1)
        spin = 1
      else
        ! Single beta
        m = exc(1,1,2)
        p = exc(1,2,2)
        spin = 2
      endif
      i_op_j = (monoe_ints(m,p) + monoe_ints(p,m)) * phase

    case (0)
      ! < I | op | I >
      ! alpha part
      call bitstring_to_list(key_i, occ_a, elec_alpha_num, Nint)
      do a = 1, elec_alpha_num
        p = occ_a(a)
        i_op_j = i_op_j + monoe_ints(p,p)
      enddo
      ! beta part
      call bitstring_to_list(key_i, occ_b, elec_beta_num, Nint)
      do b = 1, elec_beta_num
        p = occ_b(b)
        i_op_j = i_op_j + monoe_ints(p,p)
      enddo
  end select
end

subroutine i_monoe_op_psi(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,monoe_ints,i_op_psi_array)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle I|a mono electronic operator |Psi \rangle
  !  = \sum_J c_J \langle I | a mono electronic operator | J \rangle$.
  !
  ! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|I \rangle$
  ! is connected.
  ! The i_H_psi_minilist is much faster but requires to build the
  ! minilists.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_op_psi_array(Nstate)

  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: i_op_j
  integer, allocatable           :: idx(:)

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))

  i_op_psi_array = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      call i_monoe_op_j(keys(1,1,i),key,Nint,monoe_ints,i_op_j)
      i_op_psi_array(1) = i_op_psi_array(1) + i_op_j * coef(i,1)
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      call i_monoe_op_j(keys(1,1,i),key,Nint,monoe_ints,i_op_j)
      do j = 1, Nstate
        i_op_psi_array(j) = i_op_psi_array(j) + i_op_j * coef(i,j)
      enddo
    enddo

  endif

end

subroutine i_aa_j(key_i,key_j,Nint,i_aa_j)

  use bitmasks

  implicit none

  BEGIN_DOC
  ! Returns $\langle I| a_a^\dagger a_i | J \rangle$ with the phase, where $I$ and $J$ are determinants.
  ! a_a^\dagger creates an electron in th MO $\phi_a$
  ! a_i annihilates an electron in th MO $\phi_i$
  END_DOC

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: i_aa_j

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  integer                        :: m,n,p,q,a,b
  integer                        :: i,j,k
  !integer                        :: occ(Nint*bit_kind_size,2)
  integer                        :: occ_a(Nint*bit_kind_size)
  integer                        :: occ_b(Nint*bit_kind_size)
  !double precision               :: diag_H_mat_elem
  double precision               :: phase
  !integer                        :: n_occ_ab(2)
  PROVIDE mo_two_e_integrals_in_map mo_integrals_map big_array_exchange_integrals

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  i_aa_j = 0d0
  !DIR$ FORCEINLINE
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    ! < I | a_a^\dagger a_i | J >
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
      !DIR$ FORCEINLINE
      i_aa_j = phase

    ! < I | a_i^\dagger a_i | I >
    case (0)
      i_aa_j = 1

  end select
end

subroutine i_aa_psi(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_aa_psi)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle I| a_a^\dagger a_i |Psi \rangle
  !  = \sum_J c_J \langle I | a_a^\dagger a_i  | J \rangle$.
  !
  ! a_a^\dagger creates an electron in th MO $\phi_a$
  ! a_i annihilates an electron in th MO $\phi_i$
  !
  ! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|I \rangle$
  ! is connected.
  ! The i_H_psi_minilist is much faster but requires to build the
  ! minilists.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_aa_psi

  integer                        :: i, ii,j
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: i_op_j
  integer, allocatable           :: idx(:)

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))

  i_aa_psi = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)

  do ii=1,idx(0)
    i = idx(ii)
    !DIR$ FORCEINLINE
    call i_aa_j(keys(1,1,i),key,Nint,i_a_j)
    i_aa_psi = i_op_psi + i_aa_j * coef(i,1)
  enddo

end
