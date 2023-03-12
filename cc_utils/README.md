# cc_utils

Utils for the CC modules.

## Keywords
- cc_thresh_conv: Threshold for the convergence of the residual equations. Default: 1e-6.
- cc_max_iter: Maximum number of iterations. Default: 100.
- cc_diis_depth: Diis depth. Default: 8.
- cc_level_shift: Level shift for the CC. Default: 0.0.
- cc_level_shift_guess: Level shift for the MP guess of the amplitudes. Default: 0.0.
- cc_update_method: Method used to update the CC amplitudes. none -> normal, diis -> with diis. Default: diis.
- cc_guess_t1: Guess used to initialize the T1 amplitudes. none -> 0, MP -> perturbation theory, read -> read from disk. Default: MP.
- cc_guess_t2: Guess used to initialize the T2 amplitudes. none -> 0, MP -> perturbation theory, read -> read from disk. Default: MP.
- cc_write_t1: If true, it will write on disk the T1 amplitudes at the end of the calculation. Default: False.
- cc_write_t2: If true, it will write on disk the T2 amplitudes at the end of the calculation. Default: False.
- cc_par_t: If true, the CCSD(T) will be computed.
- cc_ref: Index of the reference determinant in psi_det for CC calculation. Default: 1.


