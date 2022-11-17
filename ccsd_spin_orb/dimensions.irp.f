BEGIN_PROVIDER [ integer, spin_mo_num ]
 implicit none
 BEGIN_DOC
 ! Number of spin-orbitals
 END_DOC
 spin_mo_num = mo_num * 2
END_PROVIDER


 BEGIN_PROVIDER [ integer, spin_occ_num ]
&BEGIN_PROVIDER [ integer, spin_vir_num ]
 implicit none
 BEGIN_DOC
 ! Number of occupied/virtual spin-orbitals
 END_DOC
 spin_occ_num = elec_num
 spin_vir_num = spin_mo_num - spin_occ_num
END_PROVIDER


