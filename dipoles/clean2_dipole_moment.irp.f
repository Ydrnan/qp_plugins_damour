subroutine clean2_dipole_moment

  implicit none

  integer :: i,j,istate,jstate
  integer :: p,q,tmp_p
  integer :: degree, exc(0:2,2,2)
  integer :: h1,p1,h2,p2,s1,s2
  double precision :: phase, d, nuclei_part_x, nuclei_part_y, nuclei_part_z
  double precision :: d_x, d_y, d_z
  integer, allocatable :: list_occ_a(:), list_occ_b(:)

  allocate(list_occ_a(elec_alpha_num), list_occ_b(elec_beta_num))

  print*,''
  print*,'************ Clean dipole *********************'
  
  ! nucl part
  nuclei_part_x = 0.d0
  nuclei_part_y = 0.d0
  nuclei_part_z = 0.d0
  do i = 1,nucl_num
    nuclei_part_x += nucl_charge(i) * nucl_coord(i,1)
    nuclei_part_y += nucl_charge(i) * nucl_coord(i,2)
    nuclei_part_z += nucl_charge(i) * nucl_coord(i,3)
  enddo

  do istate = 1, 1! N_states
    do jstate = istate, N_states  
      d_x = 0d0
      d_y = 0d0
      d_z = 0d0
    
      do i = 1, N_det
        do j = 1, N_det 
          if (i == j) then

            ! diag part
            call bitstring_to_list(psi_det_alpha(N_int,i), list_occ_a, elec_alpha_num, N_int)
            call bitstring_to_list(psi_det_beta(N_int,i), list_occ_b, elec_beta_num, N_int)
            
            ! alpha
            do p = 1, elec_alpha_num
              tmp_p = list_occ_a(p)
              d_x = d_x + psi_coef(i,istate) * psi_coef(i,jstate) * mo_dipole_x(tmp_p,tmp_p)
              d_y = d_y + psi_coef(i,istate) * psi_coef(i,jstate) * mo_dipole_y(tmp_p,tmp_p)
              d_z = d_z + psi_coef(i,istate) * psi_coef(i,jstate) * mo_dipole_z(tmp_p,tmp_p)
            enddo
       
            ! beta
            do p = 1, elec_beta_num
              tmp_p = list_occ_b(p)
              d_x = d_x + psi_coef(i,istate) * psi_coef(i,jstate) * mo_dipole_x(tmp_p,tmp_p)
              d_y = d_y + psi_coef(i,istate) * psi_coef(i,jstate) * mo_dipole_y(tmp_p,tmp_p)
              d_z = d_z + psi_coef(i,istate) * psi_coef(i,jstate) * mo_dipole_z(tmp_p,tmp_p)
            enddo
  
          else
        
            ! extra diag part
            call get_excitation_degree(psi_det(1,1,i),psi_det(1,1,j),degree,N_int)
            if (degree == 1) then
              exc = 0
              call get_single_excitation(psi_det(1,1,i),psi_det(1,1,j),exc,phase,N_int)
              call decode_exc(exc,1,h1,p1,h2,p2,s1,s2)
              !print*,'indexes',h1,p1,h2,p2,s1,s2
              d_x = d_x + psi_coef(i,istate) * psi_coef(j,jstate) * phase * mo_dipole_x(h1,p1)
              d_y = d_y + psi_coef(i,istate) * psi_coef(j,jstate) * phase * mo_dipole_y(h1,p1)
              d_z = d_z + psi_coef(i,istate) * psi_coef(j,jstate) * phase * mo_dipole_z(h1,p1)
            endif
          endif
        enddo
      enddo
  
      d_x = - d_x
      d_y = - d_y
      d_z = - d_z

      if (istate == jstate) then
        d_x = d_x + nuclei_part_x
        d_y = d_y + nuclei_part_y
        d_z = d_z + nuclei_part_z
      endif

      d = dsqrt(d_x**2 + d_y**2 + d_z**2)
     
      write(*,'(A6,I2,A4,I2,A10,F12.6,A3)') 'State:', istate, ' -> ', jstate, ' Delta_E = ', (ci_energy(jstate) - ci_energy(istate))/0.0367502d0, ' eV'
      write(*,'(A6,I1,A3,A4,I1,A4,1pE16.8,A5)') ' <Psi_',istate,'|r|','Psi_',jstate,'> = ', d, ' (au)'
      write(*,'(A6,I1,A3,A4,I1,A4,1pE16.8,A5)') ' <Psi_',istate,'|r|','Psi_',jstate,'> = ', d*au2D, ' (D)'
      write(*,'(A3,I1,A1,I2,A3,1pE16.8)') ' f_',istate,'^',jstate,' = ', d*d*(2d0/3d0) * (ci_energy(jstate) - ci_energy(istate))
      print*,''

    enddo
  enddo
  print*,'***********************************************'
  print*,''

  deallocate(list_occ_a,list_occ_b)

end
