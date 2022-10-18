program diag

  implicit none

  call diagonalize_ci
  TOUCH mo_coef mo_integrals_map
  call diagonalize_ci
  TOUCH mo_coef mo_integrals_map
  call diagonalize_ci
  TOUCH mo_coef mo_integrals_map
  call diagonalize_ci
  TOUCH mo_coef mo_integrals_map
  call diagonalize_ci

end
