program true_fci

  implicit none
  
  ! Generate the FCI wave function if 
  ! Cipsi:
  ! excitation_max = -1
  ! else generate the wave function with all the excitations 
  ! <= excitation_max

  read_wf = .False.
  TOUCH read_wf

  call gen_fci_wf

end
