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

  print*,''
  print*,'# Dipoles:'
  print*,'=============================================='
  print*,' Dipole moments (au)'
  print*,' State      X           Y           Z         ||µ||' 

  do istate = 1, N_states 
    write(*,'(I5,4(F12.6))') istate, d_x(istate), d_y(istate), d_z(istate), d(istate)
  enddo

  print*,''
  print*,' Dipole moments (D)'
  print*,' State      X           Y           Z         ||µ||' 

  do istate = 1, N_states 
    write(*,'(I5,4(F12.6))') istate, d_x(istate)*au2D, d_y(istate)*au2D, d_z(istate)*au2D, d(istate)*au2D
  enddo

  print*,'=============================================='
  print*,''

  deallocate(d,d_x,d_y,d_z)

 end

subroutine print_transition_dipole_moment

  implicit none

  integer :: istate,jstate
  double precision :: f, d, d_x, d_y, d_z

  print*,''
  print*,'# Transition dipoles:'
  print*,'=============================================='
  print*,' Transition dipole moments (au)'
  write(*,'(A90)') '   #  Transition        X           Y           Z         ||µ||     Dip. str.   Osc. str.'
 
  do jstate = 1, 1 !N_states
    do istate = jstate + 1, N_states
      d_x = multi_s_x_dipole_moment(istate,jstate)
      d_y = multi_s_y_dipole_moment(istate,jstate)
      d_z = multi_s_z_dipole_moment(istate,jstate)
      d = multi_s_dipole_moment(istate,jstate)
      f = 2d0/3d0 * d * d * (ci_energy(istate) - ci_energy(jstate))
      write(*,'(I4,I4,A5,I3,6(F12.6))') (istate-1), jstate, '  -> ', istate, d_x, d_y, d_z, d, d_x**2 + d_y**2 + d_z**2, f
    enddo
  enddo

  print*,''
  print*,' Transition dipole moments (D)'
  write(*,'(A90)') '   #  Transition        X           Y           Z         ||µ||     Dip. str.   Osc. str.'
  
  do jstate = 1, 1 !N_states
    do istate = jstate + 1, N_states
      d_x = multi_s_x_dipole_moment(istate,jstate) * au2D
      d_y = multi_s_y_dipole_moment(istate,jstate) * au2D
      d_z = multi_s_z_dipole_moment(istate,jstate) * au2D
      d = multi_s_dipole_moment(istate,jstate) 
      f = 2d0/3d0 * d * d * (ci_energy(istate) - ci_energy(jstate))
      d = multi_s_dipole_moment(istate,jstate) * au2D
      write(*,'(I4,I4,A5,I3,6(F12.6))') (istate-1), jstate, '  -> ', istate, d_x, d_y, d_z, d, d_x**2 + d_y**2 + d_z**2, f
    enddo
  enddo
  print*,'=============================================='
  print*,''

end

subroutine print_oscillator_strength

  implicit none
  
  integer :: istate,jstate,det1,det2,det3,k
  double precision ::coef1,coef2,coef3
  double precision :: f,d,d_x,d_y,d_z

  if (N_states == 1 .or. N_det == 1) then
    return
  endif

  print*,''
  print*,'# Oscillator strength:'
  print*,'=============================================='
  !print*, 'Transition n.    Exc (eV)  f  <S^2>  %T1  %T2'
  do jstate = 1, 1 !N_states
    do istate = jstate + 1, N_states
      d = multi_s_dipole_moment(istate,jstate)
      f = 2d0/3d0 * d * d * (ci_energy(istate) - ci_energy(jstate))

      write(*,'(A13,I3,A1,F12.6,A7,F12.6)') 'Transition n.', (istate-1), ':', (ci_energy(istate) - ci_energy(jstate))/0.0367502d0, ' eV  f=',f
      write(*,'(A7,I4,A4,I4,A6,F8.4,A6,F8.4)') '       ', istate, '  ->', jstate, ', %T1=', percent_exc(2,istate), ', %T2=',percent_exc(3,istate)
      write(*,'(A6,F6.2,A2,F6.2)') '<S^2> ', s2_values(istate),'  ', s2_values(jstate) 
   
      if (print_det_state) then
        call print_state(jstate)
        call print_state(istate)
      endif

    enddo
  enddo

  print*,'=============================================='
  print*,''

end

subroutine print_state(istate)

  implicit none
  
  integer, intent(in) :: istate
  integer :: degree, h1,p1,h2,p2, exc(0:2,2)
  double precision :: phase
  integer, allocatable :: det(:)
  double precision, allocatable :: coef(:)
  integer :: i
  
  n_det_print = min(N_det,n_det_print)

  if (N_det == 1) then
    return
  endif

  allocate(det(N_det), coef(N_det))

  ! Key for dsort
  do i = 1, N_det
    det(i) = i
  enddo

  ! -abs to sort in the right order
  do i = 1, N_det
    coef(i) = -dabs(psi_coef(i,istate))
  enddo
 
  call dsort(coef, det, N_det)

  ! I don't like minus sign
  do i = 1, N_det
    coef(i)  = dabs(coef(i))
  enddo

  print*,''
  write(*,'(A8,I3)') '# State ', istate
  do i = 1, n_det_print
    print*,''
    write(*,'(A6,I10,A8,I4,A11,1pE14.6)') 'Det n. ', i, ', state n.', istate, ', |coef| = ', coef(i)

    ! Exc
    call get_excitation_degree(psi_det(N_int,1,1),psi_det(N_int,1,det(i)),degree,N_int)
    write(*,'(A20,I3)') 'Excitation degree = ', degree

    ! Print exc
    if (degree <= 4 .and. degree > 0) then
      call get_excitation_spin(psi_det_alpha(N_int,1),psi_det_alpha(N_int,det(i)),exc,degree,phase,N_int) 
      call decode_exc_spin(exc,h1,p1,h2,p2)
      if (degree == 1 ) then
        write(*,'(A4,I4,A3,I4)') 'Exc:', h1, ' ->', p1
      elseif (degree == 2) then
        write(*,'(A4,I4,A3,I4,A1,I4,A3,I4)') 'Exc:', h1, ' ->', p1, ',', h2, ' ->', p2
      endif

      call get_excitation_spin(psi_det_beta(N_int,1),psi_det_beta(N_int,det(i)),exc,degree,phase,N_int)
      call decode_exc_spin(exc,h1,p1,h2,p2)
      if (degree == 1) then
        write(*,'(A4,I4,A3,I4)') 'Exc:', h1, ' ->', p1
      elseif (degree == 2) then
        write(*,'(A4,I4,A3,I4,A1,I4,A3,I4)') 'Exc:', h1, ' ->', p1, ',', h2, ' ->', p2
      endif
    endif

    ! Print det
    call print_det(psi_det(N_int,1,det(i)),N_int)
  enddo
  print*,''

  deallocate(det, coef)

end
