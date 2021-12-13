program opti
  implicit none

  read_wf = .true. ! must be True for the orbital optimization !!!
  TOUCH read_wf
  call algo_trust2(n_mo_dim_act, dim_list_act_orb, list_act)

end

subroutine algo_trust(tmp_n, tmp_list_size,tmp_list)

  implicit none

  ! Variables

  ! In
  integer, intent(in) :: tmp_n, tmp_list_size,tmp_list(tmp_list_size)
  ! condition if tmp_list_size==0 -> cartesian ?

  ! Out
  ! Rien ou un truc pour savoir si ça c'est bien passé
  
  ! Internal
  double precision, allocatable :: e_val(:), W(:,:), tmp_R(:,:), R(:,:), tmp_x(:),tmp_m_x(:,:)
  double precision, allocatable :: prev_mos(:,:)
  double precision              :: criterion, prev_criterion, criterion_model
  double precision              :: delta, rho
  logical                       :: not_converged, cancel_step, must_exit,enforce_step_cancellation
  integer                       :: nb_iter, info
  integer                       :: i,j,tmp_i,tmp_j

  allocate(W(tmp_n,tmp_n),e_val(tmp_n),tmp_x(tmp_n),tmp_m_x(tmp_list_size,tmp_list_size))
  allocate(tmp_R(tmp_list_size,tmp_list_size), R(mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num))

  ! Necessary ?
  PROVIDE my_hessian_list_opt my_gradient_list_opt my_st_av_energy 

  ! Initialization
  delta = 0d0 
  nb_iter = 0 ! Must starts at 0 !!!
  rho = 0.5d0 ! Must starts at 0.5 
  not_converged = .True. ! Must be true

  ! Compute the criterion before the loop
  prev_criterion = my_st_av_energy

  do while (not_converged)

      if (nb_iter > 0) then
        PROVIDE my_hessian_list_opt my_gradient_list_opt
      endif

      ! Diagonalization of the hessian 
      call diagonalization_hessian(tmp_n,my_hessian_list_opt,e_val,W)

      cancel_step = .True. ! To enter in the loop just after 

      ! Loop to Reduce the trust radius until the criterion decreases and rho >= thresh_rho
      do while (cancel_step)

          ! Hessian,gradient,Criterion -> x 
          call trust_region_step_w_expected_e(tmp_n, my_hessian_list_opt ,W, e_val, my_gradient_list_opt, prev_criterion, rho, nb_iter, delta, criterion_model, tmp_x, must_exit)
          
          if (must_exit) then
              ! if step_in_trust_region sets must_exit on true for numerical reasons
              print*,'algo_trust1 sends the message : Exit'
              exit
          endif

          ! 1D tmp -> 2D tmp 
          call vec_to_mat_v2(tmp_n,tmp_list_size,tmp_x,tmp_m_x)

          ! Rotation submatrix (square matrix tmp_list_size by tmp_list_size)
          call rotation_matrix(tmp_m_x,tmp_list_size,tmp_R,tmp_list_size,tmp_list_size,info, enforce_step_cancellation)

          ! In the case where the error in the rotation matrix is greater thazn 1d-12
          if (enforce_step_cancellation) then
           print*, 'Forces the step cancellation, too big error in the rotation matrix'
           rho = 0d0
           cycle
         endif


          ! TODO subroutine ?
          ! tmp_R to R, subspace to full space
          R = 0d0
          do i = 1, mo_num
              R(i,i) = 1d0 ! 1 on the diagonal because it is a rotation matrix, 1 = nothing change for the corresponding orbital
          enddo
          do tmp_j = 1, tmp_list_size
              j = tmp_list(tmp_j)
              do tmp_i = 1, tmp_list_size
                  i = tmp_list(tmp_i)
                  R(i,j) = tmp_R(tmp_i,tmp_j)
              enddo
          enddo

          ! Rotation of the MOs
          call apply_mo_rotation(R, prev_mos)

          ! touch mo_coef 
          print*,'On TOUCH mo_coef my_hessian_list_opt my_gradient_list_opt my_CC2_opt'
          TOUCH mo_coef my_hessian_list_opt my_gradient_list_opt my_CC2_opt my_st_av_energy

          ! To update the other parameters if needed
          !call update_parameters()

          ! New criteriion
          print*,'On update les parametres'
          call clear_mo_map
          TOUCH mo_coef psi_det psi_coef
          call diagonalize_ci
          call save_wavefunction_unsorted

          print*,'on reprovide E'
          PROVIDE my_st_av_energy
          criterion = my_st_av_energy

          ! Criterion -> step accepted or rejected 
          call trust_region_is_step_cancelled(nb_iter,prev_criterion, criterion, criterion_model,rho,cancel_step)

          print*,'STEP CANCELLED ???:', cancel_step

          ! Cancel the previous step (mo_coef = prev_mos if you keep them...)
          if (cancel_step) then
              ! Replacement by the previous MOs
              mo_coef = prev_mos

              ! Avoid the recomputation of the hessian and the gradient
              print*,'On TOUCH my_hessian_list_opt my_gradient_list_opt my_CC2_opt'
              TOUCH my_hessian_list_opt my_gradient_list_opt my_CC2_opt
          endif      

      enddo

      call save_mos() !### depend of the time for 1 iteration

      ! To exit the external loop if must_exit = .True.
      if (must_exit) then
          exit
      endif 

      ! Step accepted, nb iteration + 1
      nb_iter = nb_iter + 1

      ! To invalid the gradient and the hessian
      FREE my_hessian_list_opt my_gradient_list_opt

      ! Unnecessary
      PROVIDE my_CC2_opt

      ! To exit
      if (dabs(my_CC2_opt) < 1d-6) then
        not_converged = .False.
      endif

      if (nb_iter > 10) then
        not_converged = .False.
      endif 

  enddo
  
 deallocate(e_val, W, tmp_R, R, tmp_x, prev_mos)

end


subroutine algo_trust2(tmp_n, tmp_list_size, tmp_list)

  implicit none

  ! Variables

  ! In
  integer, intent(in) :: tmp_n, tmp_list_size, tmp_list(tmp_list_size)

  ! Out
  ! Rien ou un truc pour savoir si ça c'est bien passé
  
  ! Internal
  double precision, allocatable :: e_val(:), W(:,:), tmp_R(:,:), R(:,:), tmp_x(:), tmp_m_x(:,:)
  double precision, allocatable :: prev_mos(:,:)
  double precision              :: criterion, prev_criterion, criterion_model
  double precision              :: delta, rho
  logical                       :: not_converged, cancel_step, must_exit
  integer                       :: nb_iter, info
  integer                       :: i,j,tmp_i,tmp_j

  allocate(W(tmp_n,tmp_n),e_val(tmp_n),tmp_x(tmp_n),tmp_m_x(tmp_list_size,tmp_list_size))
  allocate(tmp_R(tmp_list_size,tmp_list_size), R(mo_num,mo_num))
  allocate(prev_mos(ao_num,mo_num))

  ! Necessary ?
  PROVIDE my_hessian_list_opt my_gradient_list_opt my_st_av_energy 

  ! Initialization
  delta = 0d0 
  nb_iter = 0 ! Must starts at 0 !!!
  rho = 0.5d0 ! Must starts at 0.5 
  not_converged = .True. ! Must be true

  ! Compute the criterion before the loop
  prev_criterion = my_st_av_energy

  do while (not_converged)

      if (nb_iter > 0) then
        PROVIDE my_hessian_list_opt my_gradient_list_opt
      endif

      ! Diagonalization of the hessian 
      call diagonalization_hessian(tmp_n, my_hessian_list_opt, e_val, W)

      cancel_step = .True. ! To enter in the loop just after 

      ! Loop to Reduce the trust radius until the criterion decreases and rho >= thresh_rho
      do while (cancel_step)

          ! Hessian,gradient,Criterion -> x 
          call trust_region_step_w_expected_e(tmp_n, my_hessian_list_opt, W, e_val, my_gradient_list_opt, &
               prev_criterion, rho, nb_iter, delta, criterion_model, tmp_x, must_exit)

          if (must_exit) then
              ! if step_in_trust_region sets must_exit on true for numerical reasons
              print*,'trust_region_step_w_expected_e sent the message : Exit'
              exit
          endif

          ! 1D tmp -> 2D tmp 
          call vec_to_mat_v2(tmp_n, tmp_list_size, tmp_x, tmp_m_x)

          ! Rotation submatrix (square matrix tmp_list_size by tmp_list_size)
          call rotation_matrix(tmp_m_x, tmp_list_size, tmp_R, tmp_list_size, tmp_list_size, info)

          ! TODO subroutine ?
          ! tmp_R to R, subspace to full space
          R = 0d0
          do i = 1, mo_num
              R(i,i) = 1d0 ! 1 on the diagonal because it is a rotation matrix, 1 = nothing change for the corresponding orbital
          enddo
          do tmp_j = 1, tmp_list_size
              j = tmp_list(tmp_j)
              do tmp_i = 1, tmp_list_size
                  i = tmp_list(tmp_i)
                  R(i,j) = tmp_R(tmp_i, tmp_j)
              enddo
          enddo

          ! Rotation of the MOs
          call apply_mo_rotation(R, prev_mos)

          ! touch mo_coef 
          TOUCH mo_coef my_hessian_list_opt my_gradient_list_opt my_st_av_energy my_CC2_opt

          ! To update the other parameters if needed
          call update_parameters()

          ! New criterion
          PROVIDE my_st_av_energy
          criterion = my_st_av_energy

          ! Criterion -> step accepted or rejected 
          call trust_region_is_step_cancelled(nb_iter, prev_criterion, criterion, criterion_model, rho, cancel_step)

          ! Cancel the previous step (mo_coef = prev_mos if you keep them...)
          if (cancel_step) then
              ! Replacement by the previous MOs
              mo_coef = prev_mos

              ! Avoid the recomputation of the hessian and the gradient
              TOUCH mo_coef my_hessian_list_opt my_gradient_list_opt my_CC2_opt 
          endif      

      enddo

      call save_mos() !### depend of the time for 1 iteration

      ! To exit the external loop if must_exit = .True.
      if (must_exit) then
          exit
      endif 

      ! Step accepted, nb iteration + 1
      nb_iter = nb_iter + 1

      ! To invalid the gradient and the hessian
      FREE my_hessian_list_opt my_gradient_list_opt

      PROVIDE my_CC2_opt

      ! To exit
      if (dabs(my_CC2_opt) < thresh_opt_max_elem_grad) then
        not_converged = .False.
      endif

      if (nb_iter > optimization_max_nb_iter) then
        not_converged = .False.
      endif

      if (delta < thresh_delta) then
        not_converged = .False.
      endif

  enddo
  
 deallocate(e_val, W, tmp_R, R, tmp_x, prev_mos)

end
