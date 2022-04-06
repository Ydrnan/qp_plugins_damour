program h_ci1_cipsi

  implicit none

  BEGIN_DOC
  ! Selected hierarchy ci 1, hCI1
  ! Hierarchy Configuration Interaction: Combining Seniority Number and Excitation Degree
  ! Fábris Kossoski, Yann Damour, Pierre-François Loos
  ! arXiv:2203.06154
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  excitation_max = 1
  TOUCH excitation_max

  seniority_max = 2
  TOUCH seniority_max

  call run_wrapper_cipsi
 
end
