subroutine i_r_j(key_i,key_j,Nint,i_x_j,i_y_j,i_z_j)

  use bitmasks

  implicit none

  BEGIN_DOC
  ! Returns $\langle i|H|j \rangle$ where $i$ and $j$ are determinants.
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

    case (0)
      ! < I | r | I >
      call bitstring_to_list(key_i, occ_a, elec_alpha_num, Nint)
      do a = 1, elec_alpha_num
        p = occ_a(a)
        i_x_j = i_x_j + mo_dipole_x(p,p)
        i_y_j = i_y_j + mo_dipole_y(p,p)
        i_z_j = i_z_j + mo_dipole_z(p,p)
      enddo
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
! Computes $\langle i|H|Psi \rangle  = \sum_J c_J \langle i | H | J \rangle$.
!
! Uses filter_connected_i_H_psi0 to get all the $|J \rangle$ to which $|i \rangle$
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

  ! call filter_connected_i_r_psi0()
  call filter_connected_i_H_psi0(keys,key,Nint,Ndet,idx)
  if (Nstate == 1) then

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      !call i_H_j(keys(1,1,i),key,Nint,hij)
      !i_H_psi_array(1) = i_H_psi_array(1) + coef(i,1)*hij
      call i_r_j(keys(1,1,i),key,Nint,i_x_j,i_y_j,i_z_j)
      i_x_psi_array = i_x_psi_array + i_x_j * coef(i,1)
      i_y_psi_array = i_y_psi_array + i_y_j * coef(i,1)
      i_z_psi_array = i_z_psi_array + i_z_j * coef(i,1)
    enddo

  else

    do ii=1,idx(0)
      i = idx(ii)
      !DIR$ FORCEINLINE
      !call i_H_j(keys(1,1,i),key,Nint,hij)
      call i_r_j(keys(1,1,i),key,Nint,i_x_j,i_y_j,i_z_j)
      do j = 1, Nstate
        !i_H_psi_array(j) = i_H_psi_array(j) + coef(i,j)*hij
        i_x_psi_array = i_x_psi_array + i_x_j * coef(i,j)
        i_y_psi_array = i_y_psi_array + i_y_j * coef(i,j)
        i_z_psi_array = i_z_psi_array + i_z_j * coef(i,j)
      enddo
    enddo

  endif

end
