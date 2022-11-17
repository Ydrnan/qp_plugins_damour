 BEGIN_PROVIDER [ double precision, eO, (1:spin_occ_num) ]
&BEGIN_PROVIDER [ double precision, eV, (1:spin_vir_num) ]
 implicit none
 BEGIN_DOC
 ! Eigenvalues of the spin-orbitals
 END_DOC
 eO(:) = spin_fock_matrix_diag_mo(1:spin_occ_num)
 eV(:) = spin_fock_matrix_diag_mo(spin_occ_num+1:spin_mo_num)
END_PROVIDER



BEGIN_PROVIDER [ double precision, delta_OV, (spin_occ_num, spin_vir_num) ]
 implicit none
 BEGIN_DOC
 ! Energy denominator for CC
 END_DOC
 integer :: i, a

 do a=1,spin_vir_num
   do i=1,spin_occ_num
     delta_OV(i,a) = eO(i) - eV(a)
!     delta_OV(i,a) = min(eO(i) - eV(a), -1.d-3)
   enddo
 enddo

END_PROVIDER


BEGIN_PROVIDER [ double precision, delta_OOVV, (spin_occ_num,spin_occ_num,spin_vir_num,spin_vir_num) ]
 implicit none
 BEGIN_DOC
 ! Energy denominator for CC
 END_DOC
 integer :: i,j,a,b

 do b=1,spin_vir_num
  do a=1,spin_vir_num
   do j=1,spin_occ_num
    do i=1,spin_occ_num
     delta_OOVV(i,j,a,b) = delta_OV(i,a) + delta_OV(j,b)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER


BEGIN_PROVIDER [ double precision, delta_OOOVVV, (spin_occ_num,spin_occ_num,spin_occ_num,spin_vir_num,spin_vir_num,spin_vir_num) ]
 implicit none
 BEGIN_DOC
 ! Energy denominator for CC
 END_DOC
 integer :: i,j,k,a,b,c

 do c=1,spin_vir_num
  do b=1,spin_vir_num
   do a=1,spin_vir_num
    do k=1,spin_occ_num
     do j=1,spin_occ_num
      do i=1,spin_occ_num
        delta_OOOVVV(i,j,k,a,b,c) = delta_OV(i,a) + delta_OV(j,b) + delta_OV(k,c)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
END_PROVIDER


