program test_gp

  implicit none

  call algo_trust_cartesian_template(2)

end

subroutine algo_trust_cartesian_template(tmp_n)

  implicit none

  ! Variables

  ! In
  integer, intent(in) :: tmp_n
  
  ! Out
  ! Rien ou un truc pour savoir si ça c'est bien passé
  
  ! Internal
  double precision, allocatable :: e_val(:), W(:,:), tmp_x(:)
  double precision              :: criterion, prev_criterion, criterion_model
  double precision              :: delta, rho
  logical                       :: not_converged, cancel_step, must_exit
  integer                       :: nb_iter
  integer                       :: i,j

  allocate(W(tmp_n,tmp_n),e_val(tmp_n),tmp_x(tmp_n))

  PROVIDE df2_him df_him f_him x_position cc_him
  print*,'H',df2_him
  print*,'g',df_him
  print*,'pos',x_position

  ! Initialization
  delta = 0d0 
  nb_iter = 0 ! Must starts at 0 !!!
  rho = 0.5d0 ! Must starts at 0.5 
  not_converged = .True. ! Must be true

  ! Compute the criterion before the loop
  prev_criterion = f_him

  do while (not_converged)

      if (nb_iter > 0) then
        PROVIDE df2_him df_him
      endif

      ! Diagonalization of the hessian 
      call diagonalization_hessian(tmp_n,df2_him,e_val,W)

      cancel_step = .True. ! To enter in the loop just after 

      ! Loop to Reduce the trust radius until the criterion decreases and rho >= thresh_rho
      do while (cancel_step)

          ! Hessian,gradient,Criterion -> x 
          call trust_region_step_w_expected_e(tmp_n,df2_him, W, e_val,df_him, &
               prev_criterion, rho, nb_iter, delta, criterion_model, tmp_x, must_exit)

          if (must_exit) then
              ! if step_in_trust_region sets must_exit on true for numerical reasons
              print*,'trust_region_step_w_expected_e sent the message : Exit'
              exit
          endif

          ! New coordinates, check the sign 
          print*,'prev_pos', x_position
          x_position = x_position + tmp_x
          print*,'dx', tmp_x
          print*,'new pos:',x_position

          ! touch x_position
          TOUCH x_position

          ! New criterion
          PROVIDE f_him
          criterion = f_him

          ! Criterion -> step accepted or rejected 
          call trust_region_is_step_cancelled(nb_iter,prev_criterion, criterion, criterion_model,rho,cancel_step)

          ! Cancel the previous step
          if (cancel_step) then
              ! Replacement by the previous coordinates, check the sign 
              x_position = x_position - tmp_x

              ! Avoid the recomputation of the hessian and the gradient
              TOUCH x_position df2_him df_him f_him cc_him
          endif      

      enddo

      ! To exit the external loop if must_exit = .True.
      if (must_exit) then
          exit
      endif 

      ! Step accepted, nb iteration + 1
      nb_iter = nb_iter + 1

      ! To invalid the gradient and the hessian
      FREE df2_him df_him

      ! Unnecessary
      PROVIDE cc_him

      ! To exit
      if (dabs(cc_him) < thresh_opt_max_elem_grad) then
        not_converged = .False.
      endif

      if (nb_iter > optimization_max_nb_iter) then
        not_converged = .False.
      endif

      if (delta < thresh_delta) then
        not_converged = .False.
      endif
      print*,cc_him
      print*,(dabs(cc_him) < thresh_opt_max_elem_grad)
      print*,(nb_iter > optimization_max_nb_iter)
      print*,(delta < thresh_delta)
      print*,'not_converged',not_converged
  enddo
  
 deallocate(e_val, W, tmp_x)

end

