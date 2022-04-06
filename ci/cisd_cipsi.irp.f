program cisd_cipsi

  implicit none

  BEGIN_DOC
  ! Selected cisd
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  excitation_max = 2
  TOUCH excitation_max

  seniority_max = -1
  TOUCH seniority_max

  call run_wrapper_cipsi
 
end
