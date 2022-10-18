#!/bin/bash

source $QP_ROOT/quantum_package.rc

rm -r reference_results
rm reference_results.tar.gz
rm -r *.ezfio

### H2 ###
qp create_ezfio h2.xyz -b 6-31g -o h2_doci_save.ezfio
qp set_file h2_doci_save.ezfio
qp set determinants read_wf true
qp set cipsi seniority_max 0
qp run scf > h2.scf.out
qp run fci > h2.fci.out

cp -r h2_doci_save.ezfio  h2_doci_diag_save.ezfio
qp set_file h2_doci_diag_save.ezfio
qp set orbital_optimization optimization_method diag

### N2 ###
qp create_ezfio n2.xyz -b cc-pvdz -o n2_save.ezfio
qp set_file n2_save.ezfio
qp set determinants read_wf true
qp set determinants n_det_max 200
qp run scf > n2.scf.out
qp run fci > n2.fci.out

cp -r n2_save.ezfio n2_diag_save.ezfio
qp set_file n2_diag_save.ezfio
qp set orbital_optimization optimization_method diag

#FC
qp create_ezfio n2.xyz -b cc-pvdz -o n2_fc_save.ezfio
qp set_file n2_fc_save.ezfio
qp set determinants read_wf true
qp set determinants n_det_max 200
qp run scf > n2_fc.scf.out
qp set_frozen_core > n2_fc.fc.out
qp run fci > n2_fc.fci.out

cp -r n2_fc_save.ezfio n2_fc_diag_save.ezfio
qp set_file n2_fc_diag_save.ezfio
qp set orbital_optimization optimization_method diag

### H2CO ###
qp create_ezfio h2co.xyz -b cc-pvdz -o h2co_save.ezfio
qp set_file h2co_save.ezfio
qp set determinants read_wf true
qp set determinants n_det_max 200
qp run scf > h2co.scf.out
qp run fci > h2co.fci.out

cp -r h2co_save.ezfio h2co_diag_save.ezfio
qp set_file h2co_diag_save.ezfio
qp set orbital_optimization optimization_method diag

#FC
qp create_ezfio h2co.xyz -b cc-pvdz -o h2co_fc_save.ezfio
qp set_file h2co_fc_save.ezfio
qp set determinants read_wf true
qp set determinants n_det_max 200
qp run scf > h2co_fc.scf.out
qp set_frozen_core > h2co_fc.fc.out
qp run fci > h2co_fc.fci.out

cp -r h2co_fc_save.ezfio h2co_fc_diag_save.ezfio
qp set_file h2co_fc_diag_save.ezfio
qp set orbital_optimization optimization_method diag

mkdir reference_results
cp -r *_save.ezfio reference_results/.
tar zcf reference_results.tar.gz reference_results

list_ezfio=$(ls -d *.ezfio)

echo "" > data_opt.txt

for ezfio in $list_ezfio
do
    new=$(echo $ezfio | sed "s/_save//")
    file=$(echo $new | sed "s/.ezfio//")
    cp -r $ezfio $new
    echo $new
    qp set_file $new
    qp run orb_opt_trust_v2 > $file.opt.out
    E=$(grep "Energy of state" $file.opt.out | tail -n 1 | awk '{printf $6}')
    echo $new $E >> data_opt.txt
done


