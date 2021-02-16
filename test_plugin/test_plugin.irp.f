program test_plugin
  implicit none
  BEGIN_DOC
! TODO : Put the documentation of the program here
  
  ! Check if read_wf = true, else :
  ! qp set determinant read_wf true 
 
  END_DOC
  double precision :: get_two_e_integral
  double precision :: mono_e, bi_e, nuc
  integer :: i,j,k,l

  print *, 'Hello world'
  print*, mo_num

  ! Mono e- part
  mono_e = 0d0
  
  do i=1,mo_num
    do j=1,mo_num
       mono_e = mono_e + mo_one_e_integrals(j,i) &
        *(one_e_dm_mo_beta(j,i,1) + one_e_dm_mo_alpha(j,i,1))
    enddo
  enddo 
  
  print*,'mono_e', mono_e 
 
  ! Nuclear repulsion  
  nuc = nuclear_repulsion
  
  print*,'nuc', nuc

  ! Bi e- part
  ! 1/2 factor is include in full_occ_2_rdm_spin_trace_mo
  bi_e = 0d0 
    
  do i=1,mo_num
    do j=1,mo_num
      do k=1,mo_num
        do l=1,mo_num
          bi_e = bi_e + get_two_e_integral(l,k,j,i,1) &
            * full_occ_2_rdm_spin_trace_mo(l,k,j,i,1)
        enddo
      enddo
    enddo
  enddo
  
  print*,'bi_e', bi_e

  print*, 'Total energy : ', mono_e + bi_e + nuc
  
  ! subroutine de verification
  ! call print_energy_components

end
