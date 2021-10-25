  logical, parameter :: debug=.False.
  
  double precision, parameter :: pi = 3.1415926535897932d0
  
  ! Global variables for orbital optimization
  integer, parameter :: method = 2
  integer, parameter :: nb_iter_orb_opt = 20
  integer, parameter :: nb_cancel_max = 100
  integer, parameter :: nb_cancel_tot_max = 1000
  ! method = 1 -> full hessian
  ! method = 2 -> diagonal hessian
  ! nb_iter_orb_opt: Maximal number of steps in the orbital optimization
  ! nb_cancel_max: maximal number of cancel steps
  ! nb_cancel_tot_max: maximal number of cancel step (total)

  ! Global variables for orbital localization
  integer, parameter :: nb_iter_loc = 10000
  logical, parameter :: sort_mos_by_e = .False.
  logical, parameter :: kick_in_mos = .True.
  double precision, parameter :: angle_pre_rot = 0.1d0
  integer, parameter :: method_loc = 2
  integer, parameter :: nb_iter_max_loc = 10000

  
  ! Global variables for trust region
  double precision, parameter :: thresh_rho = 0.1d0
  double precision :: thresh_eig = 1d-12
  double precision, parameter :: thresh_rho_2 = 0.1d0
  double precision, parameter :: thresh_cc = 1d-6
  integer, parameter :: nb_it_max_lambda = 100
  integer, parameter :: nb_it_max_pre_search = 20
  double precision, parameter :: thresh_wtg = 1d-6
  double precision, parameter :: thresh_wtg2 = 1d-6
  logical, parameter :: absolute_eig = .False.
  logical, parameter :: avoid_saddle = .False.
  integer, parameter :: version_avoid_saddle = 2
  integer, parameter :: version_lambda_search = 2
  double precision, parameter :: thresh_model =  1d-12
  double precision, parameter :: thresh_model_2 =  1d-12
  double precision, parameter :: thresh_delta = 1d-10
  
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
  ! avoid_saddle: test in order to avoid saddle points
  ! version_avoid_saddle: cf trust_region
  ! version_lambda_search: Research of the optimal lambda 
  !                        by solving:
  !                         * ||x||^2 - delta^2 = 0
  !                         * 1/||x||^2 - 1/delta^2 = 0
  ! thresh_model: if ABS(criterion - criterion_model) < thresh
  !               the program exit to avoid numerical problem
  ! thresh_model_2: if ABS(criterion - criterion_model) < thresh
  !                 during the research of the optimal lambda
  !                 it just prints a warning
