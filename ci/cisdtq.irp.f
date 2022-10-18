program cisdtq

  implicit none

  BEGIN_DOC
  ! True CISDTQ without cipsi selection
  END_DOC

  excitation_max = 4
  TOUCH excitation_max

  call gen_fci_wf
  call diagonalize_ci
  call save_wavefunction

end
