program hci_2_cipsi

  implicit none

  BEGIN_DOC
  ! Hierarchy CI 2
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  twice_hierarchy_max = 4
  TOUCH twice_hierarchy_max

  call run_wrapper_cipsi
 
end
