program hci_3_cipsi

  implicit none

  BEGIN_DOC
  ! Hierarchy CI 3
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  twice_hierarchy_max = 6
  TOUCH twice_hierarchy_max

  call run_wrapper_cipsi
 
end
