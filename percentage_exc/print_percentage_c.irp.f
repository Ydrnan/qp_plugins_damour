program print_percentage_c

  implicit none

  BEGIN_DOC
  ! %C_k = % of contribution of the k excited determinants to the norm
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  call run_print_percentage_c

end
