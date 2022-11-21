BEGIN_PROVIDER [ integer, spin_mo_num ]
&BEGIN_PROVIDER [ integer, spin_mo_num_nc ]
 implicit none
 BEGIN_DOC
 ! Number of spin-orbitals
 END_DOC
 spin_mo_num = mo_num * 2
 spin_mo_num_nc = mo_num * 2 - 2 * n_core_orb
END_PROVIDER


 BEGIN_PROVIDER [ integer, spin_occ_num ]
&BEGIN_PROVIDER [ integer, spin_occ_nc_num ]
&BEGIN_PROVIDER [ integer, spin_vir_num ]
 implicit none
 BEGIN_DOC
 ! Number of occupied/virtual spin-orbitals
 END_DOC
 !spin_occ_num = elec_num
 !spin_vir_num = spin_mo_num - spin_occ_num
 spin_occ_num = elec_num
 spin_occ_nc_num = elec_num - 2*n_core_orb
 spin_vir_num = spin_mo_num - elec_num
END_PROVIDER


