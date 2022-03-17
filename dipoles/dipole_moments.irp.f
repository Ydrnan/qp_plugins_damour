 BEGIN_PROVIDER [double precision, my_z_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, my_y_dipole_moment, (N_states)]
&BEGIN_PROVIDER [double precision, my_x_dipole_moment, (N_states)]
 implicit none
 BEGIN_DOC
 ! blablabla 
 END_DOC
 integer :: ipoint,istate,i,j 
 double precision :: weight, r(3)  
 double precision :: cpu0,cpu1,nuclei_part_z,nuclei_part_y,nuclei_part_x

! call cpu_time(cpu0)
 my_z_dipole_moment = 0.d0
 my_y_dipole_moment = 0.d0
 my_x_dipole_moment = 0.d0
 do istate = 1, N_states
  do i = 1, mo_num  
   my_z_dipole_moment(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_z(i,i)
   my_y_dipole_moment(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_y(i,i)
   my_x_dipole_moment(istate) += -(one_e_dm_mo_alpha(i,i,istate)+one_e_dm_mo_beta(i,i,istate)) *  mo_dipole_x(i,i)
   do j = i+1, mo_num  
    my_z_dipole_moment(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_z(j,i) 
    my_y_dipole_moment(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_y(j,i) 
    my_x_dipole_moment(istate) += - 2.d0 * (one_e_dm_mo_alpha(j,i,istate)+one_e_dm_mo_beta(j,i,istate)) *  mo_dipole_x(j,i) 
   enddo
  enddo
 enddo
 
! print*,'electron part for z_dipole = ',z_dipole_moment
! print*,'electron part for y_dipole = ',y_dipole_moment
! print*,'electron part for x_dipole = ',x_dipole_moment
! 
 nuclei_part_z = 0.d0
 nuclei_part_y = 0.d0
 nuclei_part_x = 0.d0
 do i = 1,nucl_num 
  nuclei_part_z += nucl_charge(i) * nucl_coord(i,3) 
  nuclei_part_y += nucl_charge(i) * nucl_coord(i,2) 
  nuclei_part_x += nucl_charge(i) * nucl_coord(i,1) 
 enddo
! print*,'nuclei   part for z_dipole = ',nuclei_part_z
! print*,'nuclei   part for y_dipole = ',nuclei_part_y
! print*,'nuclei   part for x_dipole = ',nuclei_part_x
!
 do istate = 1, N_states
  my_z_dipole_moment(istate) += nuclei_part_z
  my_y_dipole_moment(istate) += nuclei_part_y
  my_x_dipole_moment(istate) += nuclei_part_x
 enddo

! call cpu_time(cpu1)
! print*,'Time to provide the dipole moment :',cpu1-cpu0
END_PROVIDER

BEGIN_PROVIDER [double precision, au2D]
  implicit none
  !BEGIN_DOC
  ! Atomic units to Debye
  !END_DOC
  au2D = 2.541765d0
END_PROVIDER

 subroutine print_dipole_moments_xyz
 implicit none
  integer :: i
  print*, ''
  print*, ''
  print*,  '****************************************'
  print*, '# Dipole moments :'
  write(*,'(A17)',advance='no') ' State :         '
  do i = 1,N_states
    write(*,'(i16)',advance='no') i
  end do
  write(*,*) ''
  write(*,'(A24,100(1pE16.8))') ' x_dipole_moment (au) = ',my_x_dipole_moment
  write(*,'(A24,100(1pE16.8))') ' x_dipole_moment (D)  = ',my_x_dipole_moment * au2D
  write(*,'(A24,100(1pE16.8))') ' y_dipole_moment (au) = ',my_y_dipole_moment
  write(*,'(A24,100(1pE16.8))') ' y_dipole_moment (D)  = ',my_y_dipole_moment * au2D
  write(*,'(A24,100(1pE16.8))') ' z_dipole_moment (au) = ',my_z_dipole_moment
  write(*,'(A24,100(1pE16.8))') ' z_dipole_moment (D)  = ',my_z_dipole_moment * au2D
  !print*, 'x_dipole_moment                  = ',x_dipole_moment
  !print*, 'y_dipole_moment                  = ',y_dipole_moment
  !print*, 'z_dipole_moment                  = ',z_dipole_moment
  print*,  '****************************************'
 end
