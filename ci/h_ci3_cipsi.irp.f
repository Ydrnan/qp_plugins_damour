program h_ci3_cipsi

  implicit none

  BEGIN_DOC
  ! Selected hierarchy ci 3, hCI3
  ! Hierarchy Configuration Interaction: Combining Seniority Number and Excitation Degree
  ! Fábris Kossoski, Yann Damour, Pierre-François Loos
  ! arXiv:2203.06154
  END_DOC

  read_wf = .True.
  TOUCH read_wf

  excitation_max = 4
  TOUCH excitation_max

  seniority_max = 6
  TOUCH seniority_max

  call run_wrapper_cipsi
 
end
