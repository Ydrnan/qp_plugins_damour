program cisdt

  implicit none

  BEGIN_DOC
  ! True CISDT without cipsi selection
  END_DOC

  excitation_max = 3
  TOUCH excitation_max

  call gen_fci_wf
  call diagonalize_ci
  call save_wavefunction

end
