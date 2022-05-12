program debug_density_matrix

  implicit none

  integer :: i,j,istate,jstate
  integer :: p,q,tmp_p
  integer :: degree, exc(0:2,2,2)
  integer :: h1,p1,h2,p2,s1,s2
  double precision :: phase, max_error
  integer :: nb_error
  !double precision, allocatable :: dm(:,:,:,:)
  integer, allocatable :: list_occ_a(:), list_occ_b(:)

  PROVIDE dm multi_s_dm

  allocate(list_occ_a(elec_alpha_num), list_occ_b(elec_beta_num))
  !allocate(dm(mo_num,mo_num,N_states,N_states))

  !dm = 0d0

  do istate = 1, N_states
    do jstate = 1, N_states  
    
      do i = 1, N_det
        do j = 1, N_det 
          if (i == j) then

            ! diag part
            call bitstring_to_list(psi_det_alpha(N_int,i), list_occ_a, elec_alpha_num, N_int)
            call bitstring_to_list(psi_det_beta(N_int,i), list_occ_b, elec_beta_num, N_int)
            
            ! alpha
            do p = 1, elec_alpha_num
              tmp_p = list_occ_a(p)
              dm(tmp_p,tmp_p,istate,jstate) = dm(tmp_p,tmp_p,istate,jstate) + psi_coef(i,istate) * psi_coef(i,jstate)
            enddo
       
            ! beta
            do p = 1, elec_beta_num
              tmp_p = list_occ_b(p)
              dm(tmp_p,tmp_p,istate,jstate) = dm(tmp_p,tmp_p,istate,jstate) + psi_coef(i,istate) * psi_coef(i,jstate)
            enddo
  
          else
        
            ! extra diag part
            call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
            if (degree == 1) then
              exc = 0
              call get_single_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,phase,N_int)
              call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
              dm(h1,p1,istate,jstate) = dm(h1,p1,istate,jstate) + psi_coef(i,istate) * psi_coef(j,jstate) * phase
            endif
          endif
        enddo
      enddo
      
      nb_error = 0
      max_error = 0d0
      print*,'States:',istate,jstate
      do q = 1, mo_num
        do p = 1, mo_num
         print*,dm(p,q,istate,jstate), multi_s_dm(p,q,istate,jstate),dm(p,q,istate,jstate) - multi_s_dm(p,q,istate,jstate)
          !if (dabs(dm(p,q,istate,jstate) - multi_s_dm(p,q,istate,jstate)) > 1d-12) then
            !nb_error = nb_error + 1
            !if (dabs(dm(p,q,istate,jstate) - multi_s_dm(p,q,istate,jstate)) > max_error) then
            !  max_error = dabs(dm(p,q,istate,jstate))
            !endif
          !endif
        enddo
      enddo
      print*,'nb_error:',nb_error
      print*,'max_error:',max_error
      print*,''

    enddo
  enddo

  deallocate(list_occ_a,list_occ_b)!,dm)

end
