program s_ci0_cipsi

  implicit none

  BEGIN_DOC
  ! Selected DOCI
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  excitation_max = -1
  TOUCH excitation_max

  seniority_max = 0
  TOUCH seniority_max

  call run_wrapper_cipsi
 
end
