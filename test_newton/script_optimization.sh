#!/bin/bash
#SBATCH -p xeonv3 -N 1 -n 1 -c 24 --exclusive
source /home/ydamour/qp2/quantum_package.rc
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=benzene
BASIS=HF_opt_S2
DIR=${XYZ}_${BASIS}.ezfio
FILE=${XYZ}_${BASIS}
INITIAL_DIR

cp -r ${INITIAL_DIR} ${DIR}

qp set_file $DIR
qp reset -a
qp reset -d
qp run scf > ${DIR}/${FILE}.scf

# frozen core ?
qp set_frozen_core > ${DIR}/${FILE}.frozen

qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF

# S^2 true or false
qp set determinants s2_eig true
#qp set determinants s2_eig false

# Number of determinants
Ndet=5
Nb_max_det = 200000

while [ ${Ndet} -lt ${Nb_max_det} ]
do
    echo ${Ndet}
    qp set determinants n_det_max ${Ndet}
    qp run fci > ${DIR}/${FILE}_${Ndet}.fci
 
    grep "Summary at N_det = " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_n_det.dat
    grep "# E   " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_energy.dat
    grep "# PT2   " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_pt2.dat
    grep "# rPT2   " ${DIR}/${FILE}_${Ndet}.fci | tail -1 >> ${DIR}/${FILE}_rpt2.dat

    qp run orb_opt_trust > ${DIR}/${FILE}_opt_trash_${Ndet}.dat
    grep "Max element in gardient :" ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_grad_${Ndet}.dat
    grep "Number of negative eigenvalues :"  ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_eval_${Ndet}.dat        
    grep "e_val < 0 :"  ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_eval_value_${Ndet}.dat
    grep "Energy of state    1" ${DIR}/${FILE}_opt_trash_${Ndet}.dat > ${DIR}/${FILE}_opt_energy_${Ndet}.dat

    qp run pt2 > ${DIR}/${FILE}_${Ndet}.pt2 

    grep "Summary at N_det = " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_n_det.dat
    grep "# E   " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_energy.dat
    grep "# PT2   " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_pt2.dat
    grep "# rPT2   " ${DIR}/${FILE}_${Ndet}.pt2 | tail -1 >> ${DIR}/${FILE}_rpt2.dat

Ndet=$[${Ndet}*2]    
done

paste ${DIR}/${FILE}_n_det_opt.dat ${DIR}/${FILE}_energy_opt.dat > ${DIR}/${FILE}_result_opt.dat
paste ${DIR}/${FILE}_n_det.dat ${DIR}/${FILE}_energy.dat > ${DIR}/${FILE}_tmp1.dat
paste ${DIR}/${FILE}_tmp1.dat ${DIR}/${FILE}_pt2.dat > ${DIR}/${FILE}_tmp2.dat
paste ${DIR}/${FILE}_tmp2.dat ${DIR}/${FILE}_rpt2.dat > ${DIR}/${FILE}_result.dat

# CIPSI calculation
DIR2=${FILE}_cipsi_S2.ezfio
FILE2=${FILE}_cipsi_S2
cp -r ${DIR} ${DIR2}

qp set_file ${DIR2}
qp reset -d
qp set determinants read_wf true
qp set determinants mo_label MCSCF
qp set mo_basis mo_label MCSCF

# S^2 true or false
qp set determinants s2_eig true
#qp set determinants s2_eig false

# Number of determinants
qp set determinants n_det_max 3e6

qp run fci > ${DIR2}/${FILE2}.fci

grep "Summary at N_det = " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_n_det.dat
grep "# E   " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_energy.dat
grep "# PT2   " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_pt2.dat
grep "# rPT2   " ${DIR2}/${FILE2}.fci >> ${DIR2}/${FILE2}_rpt2.dat

paste ${DIR2}/${FILE2}_n_det_opt.dat ${DIR2}/${FILE2}_energy_opt.dat > ${DIR2}/${FILE2}_result_opt.dat

paste ${DIR2}/${FILE2}_n_det.dat ${DIR2}/${FILE2}_energy.dat > ${DIR2}/${FILE2}_tmp1.dat
paste ${DIR2}/${FILE2}_tmp1.dat ${DIR2}/${FILE2}_pt2.dat > ${DIR2}/${FILE2}_tmp2.dat
paste ${DIR2}/${FILE2}_tmp2.dat ${DIR2}/${FILE2}_rpt2.dat > ${DIR2}/${FILE2}_result.dat

#       -p xeonv1 -N 1 -n 1 -c 16 --exclusive
#       -p xeonv2 -N 1 -n 1 -c 20 --exclusive
#       -p xeonv3 -N 1 -n 1 -c 24 --exclusive
#       -p xeonv4 -N 1 -n 1 -c 28 --exclusive
#       -p xeonv5 -N 1 -n 1 -c 32 --exclusive
#       -p xeonv6 -N 1 -n 1 -c 32 --exclusive

