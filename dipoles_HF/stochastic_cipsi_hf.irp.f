subroutine run_stochastic_cipsi_hf
  use selection_types
  implicit none
  BEGIN_DOC
! Selected Full Configuration Interaction with Stochastic selection and PT2.
  END_DOC
  integer                        :: i,j,k
  double precision, allocatable  :: zeros(:)
  integer                        :: to_select
  type(pt2_type)                 :: pt2_data, pt2_data_err
  logical, external              :: qp_stop

  double precision :: lambda, nuclei_part_x, nuclei_part_y, nuclei_part_z, nuclear_part
  double precision, allocatable :: e_lambda(:,:,:)
  integer :: s, dip_axis
  dip_axis = 3
  lambda = 1d-6
  allocate(e_lambda(N_states,2,5))
  e_lambda = 0d0

  double precision :: rss
  double precision, external :: memory_of_double
  PROVIDE H_apply_buffer_allocated distributed_davidson mo_two_e_integrals_in_map

  N_iter = 1
  threshold_generators = 1.d0
  SOFT_TOUCH threshold_generators

  rss = memory_of_double(N_states)*4.d0
  call check_mem(rss,irp_here)

  allocate (zeros(N_states))
  call pt2_alloc(pt2_data, N_states)
  call pt2_alloc(pt2_data_err, N_states)

  double precision               :: hf_energy_ref
  logical                        :: has
  double precision               :: relative_error

  relative_error=PT2_relative_error

  zeros = 0.d0
  pt2_data % pt2   = -huge(1.e0)
  pt2_data % rpt2  = -huge(1.e0)
  pt2_data % overlap= 0.d0
  pt2_data % variance = huge(1.e0)

  if (s2_eig) then
    call make_s2_eigenfunction
  endif
  call diagonalize_CI
  call save_wavefunction

  call ezfio_has_hartree_fock_energy(has)
  if (has) then
    call ezfio_get_hartree_fock_energy(hf_energy_ref)
  else
    hf_energy_ref = ref_bitmask_energy
  endif

  if (N_det > N_det_max) then
    psi_det = psi_det_sorted
    psi_coef = psi_coef_sorted
    N_det = N_det_max
    soft_touch N_det psi_det psi_coef
    if (s2_eig) then
      call make_s2_eigenfunction
    endif
    call diagonalize_CI
    call save_wavefunction
  endif

  double precision :: correlation_energy_ratio

  correlation_energy_ratio = 0.d0

  do while (                                                         &
        (N_det < N_det_max) .and.                                    &
        (sum(abs(pt2_data % pt2(1:N_states)) * state_average_weight(1:N_states)) > pt2_max) .and.               &
        (sum(abs(pt2_data % variance(1:N_states)) * state_average_weight(1:N_states)) > variance_max) .and.     &
        (correlation_energy_ratio <= correlation_energy_ratio_max)   &
        )
      write(*,'(A)')  '--------------------------------------------------------------------------------'


    to_select = int(sqrt(dble(N_states))*dble(N_det)*selection_factor)
    to_select = max(N_states_diag, to_select)


    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,to_select) ! Stochastic PT2 and selection

    correlation_energy_ratio = (psi_energy_with_nucl_rep(1) - hf_energy_ref)  /     &
                    (psi_energy_with_nucl_rep(1) + pt2_data % rpt2(1) - hf_energy_ref)
    correlation_energy_ratio = min(1.d0,correlation_energy_ratio)

    call write_double(6,correlation_energy_ratio, 'Correlation ratio')
    call print_summary(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)

    call save_iterations(psi_energy_with_nucl_rep(1:N_states),pt2_data % rpt2,N_det)
    call print_extrapolated_energy()

    call print_dipole_moment_xyz_v2

    ! Nuclei part
    nuclei_part_x = 0.d0
    nuclei_part_y = 0.d0
    nuclei_part_z = 0.d0
  
    do i = 1,nucl_num
      nuclei_part_x += nucl_charge(i) * nucl_coord(i,1)
      nuclei_part_y += nucl_charge(i) * nucl_coord(i,2)
      nuclei_part_z += nucl_charge(i) * nucl_coord(i,3)
    enddo

    if (dip_axis == 1) then
      nuclear_part = nuclei_part_x
    elseif (dip_axis == 2) then 
      nuclear_part = nuclei_part_y
    else
      nuclear_part = nuclei_part_z
    endif

    do s = 1, N_states
      e_lambda(s,1,3) = ci_energy(s)
      e_lambda(s,2,3) = pt2_data % pt2(s)
    enddo
    print*,''
    print*,'# E(lambda) at',N_det
    do s = 1, N_states
      print*,'State:',s
      print*,'E                  ',e_lambda(s,1,:)
      print*,'E(la) - E          ',e_lambda(s,1,:) - e_lambda(s,1,3)
      print*,'dE(la)             ',e_lambda(s,1,2) - e_lambda(s,1,4)
      print*,'PT2                ',e_lambda(s,2,:)
      print*,'PT2(la) -PT2       ',e_lambda(s,2,:) - e_lambda(s,2,3)
      print*,'dPT2(la)           ',e_lambda(s,2,2) - e_lambda(s,2,4)
      print*,'E(la)+PT2(la)-E+PT2',e_lambda(s,1,:) - e_lambda(s,1,3) + e_lambda(s,2,:) - e_lambda(s,2,3)
      print*,'d(E(la)+PT2(la))   ',e_lambda(s,1,2) - e_lambda(s,1,4) + e_lambda(s,2,2) - e_lambda(s,2,4)
      print*,'Dip2', (-((e_lambda(s,1,4) - e_lambda(s,1,2)) / (2d0 * lambda)) + nuclear_part)*au2D
      print*,'Dip2pt2', (-(((e_lambda(s,1,4) + e_lambda(s,2,4)) - (e_lambda(s,1,2) + e_lambda(s,2,2))) / (2d0 * lambda)) +nuclear_part)*au2D
      !print*,'Dip4', (-((e_lambda(s,1,1) &
      !              - 8d0 * e_lambda(s,1,2) &
      !              + 8d0 * e_lambda(s,1,4) &
      !              - e_lambda(s,1,5)) / (12d0 * lambda)) +nuclear_part)*au2D
      !print*,'Dip4pt2', (-(((e_lambda(s,1,1) + e_lambda(s,2,1)) &
      !              - 8d0 * (e_lambda(s,1,2) + e_lambda(s,2,2)) &
      !              + 8d0 * (e_lambda(s,1,4) + e_lambda(s,2,4)) &
      !              - (e_lambda(s,1,5) + e_lambda(s,2,5))) / (12d0 * lambda)) +nuclear_part)*au2D
    enddo
    print*,''
    
    open(unit=11,file='fci_hf.dat',access='append')
    write(11,'(A10,I10)') adjustl('N_det'), N_det
    write(11,'(A10,I2)')  adjustl('N_states'), N_states
    write(11,'(A16,10(F16.8))') adjustl('E(-2lambda)'),   e_lambda(:,1,1)
    write(11,'(A16,10(F16.8))') adjustl('E(-1lambda)'),   e_lambda(:,1,2)
    write(11,'(A16,10(F16.8))') adjustl('E(2lambda)'),    e_lambda(:,1,4)
    write(11,'(A16,10(F16.8))') adjustl('E(2lambda)'),    e_lambda(:,1,5)
    write(11,'(A16,10(F16.8))') adjustl('PT2(-2lambda)'), e_lambda(:,2,1)
    write(11,'(A16,10(F16.8))') adjustl('PT2(-1lambda)'), e_lambda(:,2,2)
    write(11,'(A16,10(F16.8))') adjustl('PT2(1lambda)'),  e_lambda(:,2,4)
    write(11,'(A16,10(F16.8))') adjustl('PT2(2lambda)'),  e_lambda(:,2,5)
    close(11)

    N_iter += 1
    e_lambda = 0d0

    if (qp_stop()) exit

    ! Add selected determinants
    call copy_H_apply_buffer_to_wf()
    if (save_wf_after_selection) then
      call save_wavefunction
    endif

    PROVIDE  psi_coef
    PROVIDE  psi_det
    PROVIDE  psi_det_sorted

    ! mo_one_e_ints - 2 lambda
    !call pert_mo_one_e_ints(-2*lambda,dip_axis)
    !TOUCH mo_one_e_integrals
    !print*,''
    !print*,'======================================'
    !print*,'# -2lambda'
    !print*,'======================================'
    !call diagonalize_ci
    !call pt2_dealloc(pt2_data)
    !call pt2_dealloc(pt2_data_err)
    !call pt2_alloc(pt2_data, N_states)
    !call pt2_alloc(pt2_data_err, N_states)
    !call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,0) ! Stochastic PT2, no selection
    !call print_summary(psi_energy_with_nucl_rep, &
    !   pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

    !do s = 1, N_states
    !  e_lambda(s,1,1) = ci_energy(s)
    !  e_lambda(s,2,1) = pt2_data % pt2(s)
    !enddo

    ! mo_one_e_ints - lambda
    call pert_mo_one_e_ints(-lambda,dip_axis)
    TOUCH mo_one_e_integrals
    print*,''
    print*,'======================================'
    print*,'# -lambda'
    print*,'======================================'
    call diagonalize_ci
    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,0) ! Stochastic PT2, no selection
    call print_summary(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

    do s = 1, N_states
      e_lambda(s,1,2) = ci_energy(s)
      e_lambda(s,2,2) = pt2_data % pt2(s)
    enddo


    ! mo_one_e_ints + lambda
    call pert_mo_one_e_ints(lambda,dip_axis)
    TOUCH mo_one_e_integrals
    print*,''
    print*,'======================================'
    print*,'# +lambda'
    print*,'======================================'
    call diagonalize_ci
    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,0) ! Stochastic PT2, no selection
    call print_summary(psi_energy_with_nucl_rep, &
       pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

    do s = 1, N_states
      e_lambda(s,1,4) = ci_energy(s)
      e_lambda(s,2,4) = pt2_data % pt2(s)
    enddo

    ! mo_one_e_ints + 2 lambda
    !call pert_mo_one_e_ints(2*lambda,dip_axis)
    !TOUCH mo_one_e_integrals
    !print*,''
    !print*,'======================================'
    !print*,'# +2lambda'
    !print*,'======================================'
    !call diagonalize_ci
    !call pt2_dealloc(pt2_data)
    !call pt2_dealloc(pt2_data_err)
    !call pt2_alloc(pt2_data, N_states)
    !call pt2_alloc(pt2_data_err, N_states)
    !call ZMQ_pt2(psi_energy_with_nucl_rep,pt2_data,pt2_data_err,relative_error,0) ! Stochastic PT2, no selection
    !call print_summary(psi_energy_with_nucl_rep, &
    !   pt2_data, pt2_data_err, N_det,N_configuration,N_states,psi_s2)

    !do s = 1, N_states
    !  e_lambda(s,1,5) = ci_energy(s)
    !  e_lambda(s,2,5) = pt2_data % pt2(s)
    !enddo

    ! Back to the initial ones
    call init_mo_one_e_ints
    TOUCH mo_one_e_integrals

    ! Normal diagonalization
    print*,''
    print*,'======================================'
    print*,'# Normal'
    print*,'======================================'
    call diagonalize_CI
    call save_wavefunction
    call save_energy(psi_energy_with_nucl_rep, zeros)
    if (qp_stop()) exit
  enddo

  call print_dipole_moment_xyz_v2

  ! Nuclei part
  nuclei_part_x = 0.d0
  nuclei_part_y = 0.d0
  nuclei_part_z = 0.d0

  do i = 1,nucl_num
    nuclei_part_x += nucl_charge(i) * nucl_coord(i,1)
    nuclei_part_y += nucl_charge(i) * nucl_coord(i,2)
    nuclei_part_z += nucl_charge(i) * nucl_coord(i,3)
  enddo

  if (dip_axis == 1) then
    nuclear_part = nuclei_part_x
  elseif (dip_axis == 2) then
    nuclear_part = nuclei_part_y
  else
    nuclear_part = nuclei_part_z
  endif

  do s = 1, N_states
    e_lambda(s,1,3) = ci_energy(s)
    e_lambda(s,2,3) = pt2_data % pt2(s)
  enddo
  print*,''
  print*,'# E(lambda) at',N_det
  do s = 1, N_states
    print*,'State:',s
    print*,'E                  ',e_lambda(s,1,:)
    print*,'E(la) - E          ',e_lambda(s,1,:) - e_lambda(s,1,3)
    print*,'dE(la)             ',e_lambda(s,1,2) - e_lambda(s,1,4)
    print*,'PT2                ',e_lambda(s,2,:)
    print*,'PT2(la) -PT2       ',e_lambda(s,2,:) - e_lambda(s,2,3)
    print*,'dPT2(la)           ',e_lambda(s,2,2) - e_lambda(s,2,4)
    print*,'E(la)+PT2(la)-E+PT2',e_lambda(s,1,:) - e_lambda(s,1,3) + e_lambda(s,2,:) - e_lambda(s,2,3)
    print*,'d(E(la)+PT2(la))   ',e_lambda(s,1,2) - e_lambda(s,1,4) + e_lambda(s,2,2) - e_lambda(s,2,4)
    print*,'Dip2', (-((e_lambda(s,1,4) - e_lambda(s,1,2)) / (2d0 * lambda)) + nuclear_part)*au2D
    print*,'Dip2pt2', (-(((e_lambda(s,1,4) + e_lambda(s,2,4)) - (e_lambda(s,1,2) + e_lambda(s,2,2))) / (2d0 * lambda)) +nuclear_part)*au2D
  enddo

  if (.not.qp_stop()) then
    if (N_det < N_det_max) then
        call diagonalize_CI
        call save_wavefunction
        call save_energy(psi_energy_with_nucl_rep, zeros)
    endif

    call pt2_dealloc(pt2_data)
    call pt2_dealloc(pt2_data_err)
    call pt2_alloc(pt2_data, N_states)
    call pt2_alloc(pt2_data_err, N_states)
    call ZMQ_pt2(psi_energy_with_nucl_rep, pt2_data, pt2_data_err, relative_error, 0) ! Stochastic PT2

    call save_energy(psi_energy_with_nucl_rep, pt2_data % pt2)
    call print_summary(psi_energy_with_nucl_rep, &
       pt2_data , pt2_data_err, N_det, N_configuration, N_states, psi_s2)
    call save_iterations(psi_energy_with_nucl_rep(1:N_states),pt2_data % rpt2,N_det)
    call print_extrapolated_energy()
  endif
  call pt2_dealloc(pt2_data)
  call pt2_dealloc(pt2_data_err)

end
