BEGIN_PROVIDER [ double precision, spin_fock_matrix_mo, (spin_mo_num_nc, spin_mo_num_nc) ]
 implicit none
 BEGIN_DOC
 ! Fock matrix in the spin-orbital basis
 END_DOC
 integer :: p,q
 double precision :: F(mo_num,mo_num)

 spin_fock_matrix_mo = 0.d0
 call get_fock_matrix_alpha(hf_bitmask,F)
 do q=1,spin_mo_num_nc,2
  do p=1,spin_mo_num_nc,2
    spin_fock_matrix_mo(p,q) = F((p+1)/2+n_core_orb, (q+1)/2+n_core_orb)
  enddo
 enddo

 call get_fock_matrix_beta (hf_bitmask,F)
 do q=2,spin_mo_num_nc,2
  do p=2,spin_mo_num_nc,2
    spin_fock_matrix_mo(p,q) = F((p+1)/2+n_core_orb, (q+1)/2+n_core_orb)
  enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, spin_fock_matrix_diag_mo, (spin_mo_num_nc) ]
 implicit none
 BEGIN_DOC
 ! Diagonal of the Fock matrix in the spin-orbital basis
 END_DOC
 integer :: p

 do p=1,spin_mo_num_nc
    spin_fock_matrix_diag_mo(p) = spin_fock_matrix_mo(p,p)
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, spin_fock_matrix_mo_oo, (spin_occ_nc_num, spin_occ_nc_num) ]
 implicit none
 BEGIN_DOC
 ! Occupied-Occupied block of the Fock matrix in the spin-orbital basis 
 END_DOC
 integer :: p,q

 do q=1,spin_occ_nc_num
  do p=1,spin_occ_nc_num
    spin_fock_matrix_mo_oo(p,q) = spin_fock_matrix_mo(p,q) 
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, spin_fock_matrix_mo_vv, (spin_vir_num, spin_vir_num) ]
 implicit none
 BEGIN_DOC
 ! Virtual-Virtual block of the Fock matrix in the spin-orbital basis 
 END_DOC
 integer :: p,q

 do q=1,spin_vir_num
  do p=1,spin_vir_num
    spin_fock_matrix_mo_vv(p,q) = spin_fock_matrix_mo(p+spin_occ_nc_num, q+spin_occ_nc_num) 
  enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, spin_fock_matrix_mo_ov, (spin_occ_nc_num, spin_vir_num) ]
 implicit none
 BEGIN_DOC
 ! Occupied-Virtual block of the Fock matrix in the spin-orbital basis 
 END_DOC
 integer :: p,q

 do q=1,spin_vir_num
  do p=1,spin_occ_nc_num
    spin_fock_matrix_mo_ov(p,q) = spin_fock_matrix_mo(p, q+spin_occ_nc_num) 
  enddo
 enddo

END_PROVIDER



subroutine get_fock_matrix_alpha(det,F)
  implicit none
  BEGIN_DOC
! Returns the alpha Fock matrix in MO basis associated with the determinant given as input
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  double precision, intent(out) :: F(mo_num,mo_num)

  integer :: i,j,k

  F(:,:) = fock_operator_closed_shell_ref_bitmask(:,:)

end    



subroutine get_fock_matrix_beta(det,F)
  implicit none
  BEGIN_DOC
! Returns the beta Fock matrix in MO basis associated with the determinant given as input
  END_DOC
  integer(bit_kind), intent(in) :: det(N_int,2)
  double precision, intent(out) :: F(mo_num,mo_num)

  integer :: i,j,k

  F(:,:) = fock_operator_closed_shell_ref_bitmask(:,:)

end    



