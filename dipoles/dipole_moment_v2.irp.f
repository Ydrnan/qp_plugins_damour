 BEGIN_PROVIDER [double precision, multi_s_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_x_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_y_dipole_moment, (N_states, N_states)]
&BEGIN_PROVIDER [double precision, multi_s_z_dipole_moment, (N_states, N_states)]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,jstate,i,j 
 double precision :: weight, r(3)  
 double precision :: nuclei_part_z,nuclei_part_y,nuclei_part_x

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

! print*,'mat'
! do i = 1, mo_num
!   write(*,'(40E10.3)') dm(i,:,1,2) - multi_s_dm(i,:,1,2)
! enddo
! print*,'mat2'
! do i = 1, mo_num
!   write(*,'(40E10.3)') multi_s_dm(i,:,1,2)
! enddo
! call abort
!
! do jstate = 1, N_states
!   do istate = 1, N_states
!
!     do i = 1, mo_num
!       do j = 1, mo_num
!         if (dabs(dm(j,i,istate,jstate) - multi_s_dm(j,i,jstate,istate)) > 1d-12) then
!           print*,dm(j,i,istate,jstate),  multi_s_dm(j,i,jstate,istate)
!           !call abort
!         endif
!       enddo
!     enddo
!
!   enddo
! enddo

 ! Nuclei part
 nuclei_part_x = 0.d0
 nuclei_part_y = 0.d0
 nuclei_part_z = 0.d0

 do i = 1,nucl_num 
   nuclei_part_x += nucl_charge(i) * nucl_coord(i,1) 
   nuclei_part_y += nucl_charge(i) * nucl_coord(i,2) 
   nuclei_part_z += nucl_charge(i) * nucl_coord(i,3) 
 enddo

 ! Only if istate = jstate
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

!BEGIN_PROVIDER [double precision, au2D]
!  implicit none
!  !BEGIN_DOC
!  ! Atomic units to Debye
!  !END_DOC
!  au2D = 2.541765d0
!END_PROVIDER

subroutine print_dipole_moment_xyz_v2

  implicit none

  integer :: istate
  double precision, allocatable :: d(:), d_x(:), d_y(:), d_z(:)

  allocate(d(N_states),d_x(N_states),d_y(N_states),d_z(N_states))
 
  do istate = 1, N_states 
    d_x(istate) = multi_s_x_dipole_moment(istate,istate)
    d_y(istate) = multi_s_y_dipole_moment(istate,istate)
    d_z(istate) = multi_s_z_dipole_moment(istate,istate)
    d(istate) = multi_s_dipole_moment(istate,istate)
  enddo

  print*, ''
  print*,  '****************************************'
  print*, '# Dipole moment :'
  write(*,'(A17)',advance='no') ' State :         '
  do istate = 1,N_states
    write(*,'(i16)',advance='no') istate
  end do
  write(*,*) ''
  write(*,'(A20,100(1pE16.8))') ' <Psi|x|Psi> (au) = ', d_x(:)
  write(*,'(A20,100(1pE16.8))') ' <Psi|x|Psi> (D)  = ', d_x(:) * au2D
  write(*,'(A20,100(1pE16.8))') ' <Psi|y|Psi> (au) = ', d_y(:)
  write(*,'(A20,100(1pE16.8))') ' <Psi|y|Psi> (D)  = ', d_y(:) * au2D
  write(*,'(A20,100(1pE16.8))') ' <Psi|z|Psi> (au) = ', d_z(:)
  write(*,'(A20,100(1pE16.8))') ' <Psi|z|Psi> (D)  = ', d_z(:) * au2D
  write(*,'(A20,100(1pE16.8))') ' <Psi|r|Psi> (au) = ', d(:)
  write(*,'(A20,100(1pE16.8))') ' <Psi|r|Psi> (D)  = ', d(:) * au2D
  print*,  '****************************************'
  print*,''

  deallocate(d,d_x,d_y,d_z)

 end

subroutine print_oscillator_strength

  implicit none
  
  integer :: istate,jstate
  double precision :: f,d

  print*,''
  print*,'****************************************'
  print*, '# Oscillator strength :'

  do jstate = 1, 1 !N_states
    do istate = jstate + 1, N_states
      d = multi_s_dipole_moment(istate,jstate)
      print*,'d',d
      f = 2d0/3d0 * d * d * (ci_energy(istate) - ci_energy(jstate))
      write(*,'(A6,I2,A4,I2,A10,F12.6,A9,1pE16.8)') 'State:', jstate, ' -> ', istate, ' Delta_E = ', (ci_energy(istate) - ci_energy(jstate))/0.0367502d0, ' eV, f = ', f
    enddo
  enddo

  print*,'****************************************'
  print*,''

end
