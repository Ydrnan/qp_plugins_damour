#!/bin/bash

source $QP_ROOT/quantum_package.rc

FILE=
EZFIO=${FILE}.ezfio

qp set_file $EZFIO
qp set orbital_optimization optimization_method diag
qp set orbital_optimization normalized_st_av_weight true
qp set orbital_optimization start_from_wf true
#qp set orbital_optimization n_det_start 5
qp set orbital_optimization n_det_max_opt 2e5
qp set orbital_optimization targeted_accuracy_cipsi 1e-5
qp set orbital_optimization optimization_max_nb_iter 20
qp set orbital_optimization thresh_opt_max_elem_grad 1e-4

qp set determinants read_wf true
qp set determinants n_states 2
qp set davidson n_states_diag 4
qp set determinants n_det_max 2e7
qp run cis > $FILE.cis.out
qp run optimization > $FILE.opt.out

qp reset -d
tar zcf $EZFIO ${EZFIO}_save_opt_mos.tar.gz

qp run cis > $FILE.opt_cis.out
qp run fci > $FILE.oipt_fci.out

tar zcf ${EZFIO}_save_opt_res.tar.gz $EZFIO
