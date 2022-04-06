subroutine run_wrapper_cipsi

  implicit none

  if (.not.is_zmq_slave) then
    PROVIDE psi_det psi_coef mo_two_e_integrals_in_map

    if (do_pt2) then
      call run_stochastic_cipsi
    else
      call run_cipsi
    endif

  else
    PROVIDE mo_two_e_integrals_in_map pt2_min_parallel_tasks

    call run_slave_cipsi

  endif

end
