#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run() {
  thresh1=1e-4
  thresh2=1e-5
  test_exe scf || skip
  qp set_file $1
  qp edit --check
  qp set_frozen_core
  qp set mol_properties calc_dipole_moment true
  qp set mol_properties calc_tr_dipole_moment true
  qp set mol_properties calc_osc_str true
  qp set determinants n_states 2
  qp set cipsi excitation_max 2
  qp reset -d
  file="$(echo $1 | sed 's/.ezfio//g')"
  qp run fci | tee $file.fci_dip.out
  qp set determinants n_states 1
  qp set cipsi excitation_max -1
  qp reset -d
  dip1="$(grep -A 2 'Dipole moments (D)' ${file}.fci_dip.out | tail -n 1 | awk '{print $5}')"
  dip2="$(grep -A 3 'Dipole moments (D)' ${file}.fci_dip.out | tail -n 1 | awk '{print $5}')"
  fl="$(grep 'Transition n.' ${file}.fci_dip.out | tail -n 1 | awk '{print $12}')"
  fl="${fl::-1}"
  fv="$(grep 'Transition n.' ${file}.fci_dip.out | tail -n 1 | awk '{print $14}')"
  fv="${fv::-1}"
  fm="$(grep 'Transition n.' ${file}.fci_dip.out | tail -n 1 | awk '{print $16}')"
  fm="${fm::-1}"
  eq $dip1 $2 $thresh1
  eq $dip2 $3 $thresh1
  eq $fl $4 $thresh2
  eq $fv $5 $thresh2
  eq $fm $6 $thresh2
}

@test "clf" {
run  clf.ezfio 1.170762 0.556674 0.000591 0.003908 0.001519
}

@test "clo" {
run  clo.ezfio 1.106317 3.260396 0.060811 0.008104 0.022200 
}

@test "h2o" {
run  h2o.ezfio 1.961689 0.716766 0.114905 0.131365 0.122860
}

@test "h2s" {
run  h2s.ezfio 1.264523 1.352781 0.195914 0.135379 0.162858 
}

@test "lif" {
run  lif.ezfio 6.315678 3.644930 0.008834 0.004596 0.006372
}

@test "oh" {
run  oh.ezfio 1.636729 1.710037 0.000000 0.000000 0.000000
}
