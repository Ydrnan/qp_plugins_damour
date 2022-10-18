program print_percentage_t

  implicit none

  BEGIN_DOC
  ! %T_k = % of contribution of the T_k amplitudes
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  call run_print_percentage_t

end
