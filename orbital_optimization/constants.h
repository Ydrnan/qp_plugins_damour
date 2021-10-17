  logical, parameter :: debug=.False.
  
  double precision, parameter :: pi = 3.1415926535897932d0
  
  ! Global variable for orbital optimization
  integer, parameter :: method = 1
  ! method = 1 -> full hessian
  ! method = 2 -> diagonal hessian
  
  ! Global variables for trust region
  double precision, parameter :: thresh_rho = 0.1d0
  double precision :: thresh_eig = 1d-12
  double precision, parameter :: thresh_rho_2 = 0.1d0
  double precision, parameter :: thresh_cc = 1d-6
  integer, parameter :: nb_it_max_lambda = 100
  double precision, parameter :: thresh_wtg = 1d-6
  double precision, parameter :: thresh_wtg2 = 1d-6
  logical, parameter :: absolute_eig = .False.
  logical, parameter :: avoid_saddle = .True.
  integer, parameter :: version_lambda_search = 2
  
  ! thresh_rho: threshold for the step cancellation
  !             if rho < thresh_rho, the step is cancelled
  ! thresh_eig: threshold for DABS(e_val(i)+lambda) in
  !             in order to avoid division by zero
  ! thresh_rho_2: threshold for the step cancellation
  !               if rho_2 < thresh_rho_2 in the research
  !               if the optimal lambda
  ! thresh_cc: threshold for the research of lambda,
  !            to leave the loop when
  !            dabs(1d0-||x||^2/delta^2) > thresh_cc
  ! nb_it_max_lambda: maximal number of iteration to find
  !                   the optimal lambda
  ! thresh_wtg: threshold to considere when the dot product
  !             of the eigenvector w by the gradient v_grad
  !             is 0
  ! absolute_eig: if True the trust region use the absolute
  !               value of the eigenvalues
  ! avoid_saddle: test in order to avoid saddle pioints
  ! version_lambda_search : Research of the optimal lambda 
  !                         by solving:
  !                         * ||x||^2 - delta^2 = 0
  !                         * 1/||x||^2 - 1/delta^2 = 0
