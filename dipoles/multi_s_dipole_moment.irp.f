 BEGIN_PROVIDER [double precision, multi_s_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_x_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_dipole_moment, (N_states, N_states)]

  implicit none

  BEGIN_DOC
  ! Providers for :
  ! <Psi_m|µ_x|Psi_n>
  ! <Psi_m|µ_y|Psi_n>
  ! <Psi_m|µ_z|Psi_n>
  ! ||µ|| = sqrt(µ_x^2 + µ_y^2 + µ_z^2)
  !
  ! \bra{\Psi^n} x \bra{\Psi^m} = \sum_p \gamma_{pp}^{nm} \bra{\phi_p} x \ket{\phi_p} 
  !   + \sum_{pq} \gamma_{pq}^{nm} \bra{\phi_p} x \ket{\phi_q}
  ! \Psi: wf
  ! n,m indexes for the states
  ! p,q: general spatial MOs 
  ! gamma^{nm}: density matrix \bra{\Psi^n} a^{\dagger}_a a_i \ket{\Psi^m}
  END_DOC

  ! TODO: Better loops and OMP

  integer          :: istate,jstate ! States
  integer          :: i,j           ! general spatial MOs
  double precision :: nuclei_part_x, nuclei_part_y, nuclei_part_z
 
  multi_s_x_dipole_moment = 0.d0
  multi_s_y_dipole_moment = 0.d0
  multi_s_z_dipole_moment = 0.d0
 
  do jstate = 1, N_states
    do istate = 1, N_states
 
      do i = 1, mo_num  
        ! Diag part
        multi_s_x_dipole_moment(istate,jstate) -= multi_s_dm(i,i,istate,jstate) * mo_dipole_x(i,i)
        multi_s_y_dipole_moment(istate,jstate) -= multi_s_dm(i,i,istate,jstate) * mo_dipole_y(i,i)
        multi_s_z_dipole_moment(istate,jstate) -= multi_s_dm(i,i,istate,jstate) * mo_dipole_z(i,i)
 
        do j = 1, mo_num  
          if (i == j) then
           cycle
          endif
          ! Extra diag part
          multi_s_x_dipole_moment(istate,jstate) -= multi_s_dm(j,i,istate,jstate) * mo_dipole_x(j,i)  
          multi_s_y_dipole_moment(istate,jstate) -= multi_s_dm(j,i,istate,jstate) * mo_dipole_y(j,i) 
          multi_s_z_dipole_moment(istate,jstate) -= multi_s_dm(j,i,istate,jstate) * mo_dipole_z(j,i) 
        enddo
      enddo
 
    enddo
  enddo
 
  ! Nuclei part
  nuclei_part_x = 0.d0
  nuclei_part_y = 0.d0
  nuclei_part_z = 0.d0
 
  do i = 1,nucl_num 
    nuclei_part_x += nucl_charge(i) * nucl_coord(i,1) 
    nuclei_part_y += nucl_charge(i) * nucl_coord(i,2) 
    nuclei_part_z += nucl_charge(i) * nucl_coord(i,3) 
  enddo
 
  ! Only if istate = jstate, otherwise 0 by the orthogonality of the states
  do istate = 1, N_states
    multi_s_x_dipole_moment(istate,istate) += nuclei_part_x
    multi_s_y_dipole_moment(istate,istate) += nuclei_part_y
    multi_s_z_dipole_moment(istate,istate) += nuclei_part_z
  enddo
 
  ! d = <Psi|r|Psi>
  do jstate = 1, N_states
    do istate = 1, N_states
      multi_s_dipole_moment(istate,jstate) = &
        dsqrt(multi_s_x_dipole_moment(istate,jstate)**2 & 
            + multi_s_y_dipole_moment(istate,jstate)**2 &
            + multi_s_z_dipole_moment(istate,jstate)**2) 
    enddo
  enddo

END_PROVIDER

