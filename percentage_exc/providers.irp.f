BEGIN_PROVIDER [ integer, max_exc_degree]

  !BEGIN_DOC
  ! Maximal degree of excitation possible
  !END_DOC

  implicit none

  max_exc_degree = min(2*mo_num - elec_alpha_num - elec_beta_num, elec_alpha_num + elec_beta_num)

END_PROVIDER

BEGIN_PROVIDER [ double precision, percent_c, (max_exc_degree + 1, N_states)]

  !BEGIN_DOC
  ! Percentage of each kind of excitation
  !END_DOC

  implicit none

  call percentage_c(percent_c)

END_PROVIDER

