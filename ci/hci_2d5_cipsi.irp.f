program hci_2d5_cipsi

  implicit none

  BEGIN_DOC
  ! Hierarchy CI 2.5
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  twice_hierarchy_max = 5
  TOUCH twice_hierarchy_max

  call run_wrapper_cipsi
 
end
