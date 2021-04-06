program run_optimization
  implicit none
 
  N_det_max = 10
 
  call run_cipsi
  
  N_det_max = 30
 TOUCH N_det_max 
 print*,'@@@@@', N_det_max
 ! call run_cipsi

end program
