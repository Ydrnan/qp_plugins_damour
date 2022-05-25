#!/bin/bash

source $QP_ROOT/quantum_package.rc

rm -r reference_results
rm reference_results.tar.gz
rm -r *.ezfio

### BH ###
qp create_ezfio -a bh.xyz -b aug-cc-pvdz -o bh_cis_save.ezfio
qp set_file bh_cis_save.ezfio
qp set determinants read_wf true
qp set determinants n_states 2
qp set determinants n_det_max 2000
qp run scf > bh.scf.out
qp set_frozen_core > bh.fc.out
qp run cis_w_dipoles > bh.cis.out

cp -r bh_cis_save.ezfio bh_cisd_save.ezfio
qp set_file bh_cisd.ezfio
qp run cisd_w_dipoles > bh.cisd.out

cp -r bh_cis_save.ezfio bh_fci_save.ezfio
qp set_file bh_fci.ezfio
qp run cis > /dev/null
qp run fci_w_dipoles > bh.fci.out

### HCl ###
qp create_ezfio -a hcl.xyz -b aug-cc-pvdz -o hcl_cis_save.ezfio
qp set_file hcl_cis_save.ezfio
qp set determinants read_wf true
qp set determinants n_states 2
qp set determinants n_det_max 2000
qp run scf > hcl.scf.out
qp set_frozen_core > hcl.fc.out
qp run cis_w_dipoles > hcl.cis.out

cp -r hcl_cis_save.ezfio hcl_cisd_save.ezfio
qp set_file hcl_cisd.ezfio
qp run cisd_w_dipoles > hcl.cisd.out

cp -r hcl_cis_save.ezfio hcl_fci_save.ezfio
qp set_file hcl_fci.ezfio
qp run cis > /dev/null
qp run fci_w_dipoles > hcl.fci.out

### H2O ###
qp create_ezfio -a h2o.xyz -b aug-cc-pvdz -o h2o_cis_save.ezfio
qp set_file h2o_cis_save.ezfio
qp set determinants read_wf true
qp set determinants n_states 2
qp set determinants n_det_max 2000
qp run scf > h2o.scf.out
qp set_frozen_core > h2o.fc.out
qp run cis_w_dipoles > h2o.cis.out

cp -r h2o_cis_save.ezfio h2o_cisd_save.ezfio
qp set_file h2o_cisd.ezfio
qp run cisd_w_dipoles > h2o.cisd.out

cp -r h2o_cis_save.ezfio h2o_fci_save.ezfio
qp set_file h2o_fci.ezfio
qp run cis > /dev/null
qp run fci_w_dipoles > h2o.fci.out


mkdir reference_results
cp -r *.ezfio reference_results/.
tar zcf reference_results.tar.gz reference_results

list_mol=$(ls *.xyz)

echo "" > data_dip.txt
echo "" > data_osc.txt

for mol in $list_mol
do
    mol=$(echo $mol | sed "s/.xyz//")
    echo $mol
    
    dip0_cis=$(grep -A 4  "Dipole moments (D)" $mol.cis.out  | grep " 0 " | awk '{printf $5}')
    dip1_cis=$(grep -A 4  "Dipole moments (D)" $mol.cis.out  | grep " 1 " | awk '{printf $5}')
    dip0_cisd=$(grep -A 4 "Dipole moments (D)" $mol.cisd.out | grep " 0 " | awk '{printf $5}')
    dip1_cisd=$(grep -A 4 "Dipole moments (D)" $mol.cisd.out | grep " 1 " | awk '{printf $5}')
    dip0_fci=$(grep -A 4  "Dipole moments (D)" $mol.fci.out  | grep " 0 " | tail -n 1 | awk '{printf $5}')
    dip1_fci=$(grep -A 4  "Dipole moments (D)" $mol.fci.out  | grep " 1 " | tail -n 1 | awk '{printf $5}')
    echo $mol " cis  " $dip0_cis  "   " $dip1_cis  >> data_dip.txt
    echo $mol " cisd " $dip0_cisd "   " $dip1_cisd >> data_dip.txt
    echo $mol " fci  " $dip0_fci  "   " $dip1_fci  >> data_dip.txt

    fl=$(grep "#  Transition n." $mol.fci.out | grep " 1:" | tail -n 1 | awk '{printf $12}' | sed "s/,//" )
    fv=$(grep "#  Transition n." $mol.fci.out | grep " 1:" | tail -n 1 | awk '{printf $14}' | sed "s/,//" )
    fm=$(grep "#  Transition n." $mol.fci.out | grep " 1:" | tail -n 1 | awk '{printf $16}' | sed "s/,//" )
    echo $mol " " $fl " " $fv " " $fm >> data_osc.txt

done


