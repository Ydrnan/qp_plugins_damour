program hci_1_cipsi

  implicit none

  BEGIN_DOC
  ! Hierarchy CI 1
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  twice_hierarchy_max = 2
  TOUCH twice_hierarchy_max

  call run_wrapper_cipsi
 
end
