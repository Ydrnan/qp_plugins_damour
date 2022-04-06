#!/bin/bash

source $QP_ROOT/quantum_package.rc

FILES="h2o_sto_3g h2o_6_31g"

for FILE in $FILES
do
    echo $FILE
    EZFIO=$FILE.ezfio
    qp set_file $EZFIO
    qp run scf > ${FILE}.scf.out
    
    qp set coupled_cluster sccd_method bi_int
    qp run spin_orb_ccd > ${FILE}_bi_int.out
    
    qp set coupled_cluster sccd_method guess_mp2
    qp run spin_orb_ccd > ${FILE}_guess_mp2.out
    
    qp set coupled_cluster sccd_method estimated_e
    qp run spin_orb_ccd > ${FILE}_estimated_e.out
    
    grep "Result" ${FILE}_bi_int.out | cut -d ":" -f2 > ${FILE}_bi_int.out.dat
    grep "Result" ${FILE}_guess_mp2.out | cut -d ":" -f2 > ${FILE}_guess_mp2.out.dat
    grep "Result" ${FILE}_estimated_e.out | cut -d ":" -f2 > ${FILE}_estimated_e.out.dat

done

