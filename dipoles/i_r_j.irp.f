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
  integer                        :: occ_a(Nint*bit_kind_size)
  integer                        :: occ_b(Nint*bit_kind_size)
  double precision               :: phase

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  i_x_j = 0d0
  i_y_j = 0d0
  i_z_j = 0d0

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
      i_x_j = mo_dipole_x(m,p) * phase
      i_y_j = mo_dipole_y(m,p) * phase
      i_z_j = mo_dipole_z(m,p) * phase

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
      call i_r_j(keys(1,1,i),key,Nint,i_x_j,i_y_j,i_z_j)
      i_x_psi_array(1) = i_x_psi_array(1) + i_x_j * coef(i,1)
      i_y_psi_array(1) = i_y_psi_array(1) + i_y_j * coef(i,1)
      i_z_psi_array(1) = i_z_psi_array(1) + i_z_j * coef(i,1)
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      call i_r_j(keys(1,1,i),key,Nint,i_x_j,i_y_j,i_z_j)
      do j = 1, Nstate
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
  integer                        :: occ_a(Nint*bit_kind_size)
  integer                        :: occ_b(Nint*bit_kind_size)
  double precision               :: phase

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)


  i_op_j = 0d0
  call get_excitation_degree(key_i,key_j,degree,Nint)
  integer :: spin
  select case (degree)
    ! < I | op | J >
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
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

    ! < I | op | I >
    case (0)
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
  ! Computes $\langle I|a mono electronic operator |Psi \rangle$
  ! $ = \sum_J c_J \langle I | a mono electronic operator | J \rangle$.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate), monoe_ints(mo_num,mo_num)
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
      call i_monoe_op_j(keys(1,1,i),key,Nint,monoe_ints,i_op_j)
      i_op_psi_array(1) = i_op_psi_array(1) + i_op_j * coef(i,1)
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      call i_monoe_op_j(keys(1,1,i),key,Nint,monoe_ints,i_op_j)
      do j = 1, Nstate
        i_op_psi_array(j) = i_op_psi_array(j) + i_op_j * coef(i,j)
      enddo
    enddo

  endif

end

subroutine i_op_aa_j(key_i,key_j,Nint,i_aa_j,idx_h,idx_p)

  use bitmasks

  implicit none

  BEGIN_DOC
  ! Returns $\langle I| a_a^\dagger a_i |J \rangle$ where $I$ and $J$ are determinants.
  ! $a_a^\dagger$ creates an electron in the orbital $\psi_a$.
  ! $a_i^$ annihilates an electron in the orbital $\psi_i$.
  ! idx_h = i or 0 if I = J
  ! idx_p = a or 0 if I = J
  END_DOC

  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: key_i(Nint,2), key_j(Nint,2)
  double precision, intent(out)  :: i_aa_j
  integer, intent(out)           :: idx_h, idx_p

  integer                        :: exc(0:2,2,2)
  integer                        :: degree
  integer                        :: m,n,p,q,a,b
  double precision               :: phase
  integer :: spin

  ASSERT (Nint > 0)
  ASSERT (Nint == N_int)
  ASSERT (sum(popcnt(key_i(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_i(:,2))) == elec_beta_num)
  ASSERT (sum(popcnt(key_j(:,1))) == elec_alpha_num)
  ASSERT (sum(popcnt(key_j(:,2))) == elec_beta_num)

  i_aa_j = 0d0

  call get_excitation_degree(key_i,key_j,degree,Nint)

  idx_h = -1
  idx_p = -1
  select case (degree)
    ! < I | a_a^\dagger a_i | J >
    case (1)
      call get_single_excitation(key_i,key_j,exc,phase,Nint)
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
      i_aa_j = phase
      idx_h = m
      idx_p = p

    ! < I | a_a^\dagger a_i | I >
    case (0)
      i_aa_j = 1d0
      idx_h = 0
      idx_p = 0

  end select
end

subroutine i_op_aa_psi(key,keys,coef,Nint,Ndet,Ndet_max,Nstate,i_aa_psi_array)
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Computes $\langle i|a_a^\dagger a_i |Psi \rangle$.
  ! $a_a^\dagger$ creates an electron in the orbital $\psi_a$.
  ! $a_i^$ annihilates an electron in the orbital $\psi_i$.
  END_DOC
  integer, intent(in)            :: Nint, Ndet,Ndet_max,Nstate
  integer(bit_kind), intent(in)  :: keys(Nint,2,Ndet)
  integer(bit_kind), intent(in)  :: key(Nint,2)
  double precision, intent(in)   :: coef(Ndet_max,Nstate)
  double precision, intent(out)  :: i_aa_psi_array(mo_num,mo_num,Nstate)

  integer                        :: i, ii,jstate,p
  double precision               :: phase
  integer                        :: exc(0:2,2,2)
  double precision               :: i_aa_j
  integer                        :: idx_h,  idx_p
  integer, allocatable           :: idx(:)

  ASSERT (Nint > 0)
  ASSERT (N_int == Nint)
  ASSERT (Nstate > 0)
  ASSERT (Ndet > 0)
  ASSERT (Ndet_max >= Ndet)
  allocate(idx(0:Ndet))

  i_aa_psi_array = 0.d0

  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i = idx(ii)
      call i_op_aa_j(keys(1,1,i),key,Nint,i_aa_j,idx_h,idx_p)
      ! if I = J
      if (idx_h == 0) then
        do p = 1, mo_num
          i_aa_psi_array(p,p,1) = i_aa_psi_array(p,p,1) + coef(i,1) * i_aa_j
        enddo
      elseif (idx_h == -1) then
        ! nothing
        cycle
      else
        i_aa_psi_array(idx_h,idx_p,1) = i_aa_psi_array(idx_h,idx_p,1) + coef(i,1) * i_aa_j
      endif
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      call i_op_aa_j(keys(1,1,i),key,Nint,i_aa_j,idx_h,idx_p)
      do jstate = 1, Nstate
        ! if I = J
        if (idx_h == 0) then
          do p = 1, mo_num
            i_aa_psi_array(p,p,jstate) = i_aa_psi_array(p,p,jstate) + coef(i,jstate) * i_aa_j
          enddo
        elseif (idx_h == -1) then
          ! nothing
          cycle
        else
          i_aa_psi_array(idx_h,idx_p,jstate) = i_aa_psi_array(idx_h,idx_p,jstate) + coef(i,jstate) * i_aa_j
        endif
      enddo
    enddo

  endif

end 
