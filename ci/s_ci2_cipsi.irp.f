program s_ci2_cipsi

  implicit none

  BEGIN_DOC
  ! Selected seniority 2 ci
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  excitation_max = -1
  TOUCH excitation_max

  seniority_max = 2
  TOUCH seniority_max

  call run_wrapper_cipsi
 
end
