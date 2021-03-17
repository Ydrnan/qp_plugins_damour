#!/bin/bash

# for the cluster :
# Change the source !!!!!!!!!
#SBATCH -p xeonv3 -N 1 -n 1 -c 24 --exclusive
source /home/ydamour/qp2/quantum_package.rc
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=CN
BASIS=cc_pvdz
EXTRA=1e4_h_diag

IT=_it_

DIR=${XYZ}_${BASIS}_${EXTRA}.ezfio
FILE=${XYZ}_${BASIS}_${EXTRA}
PATH_CIPSI=../../../y_calculs
PATH_OPT=../plugins/qp_plugins_damour/test_newton

# path between 
# /qp2/plugin/qp_plugin_damour/test_newton the directory for optimizations
# and 
# /qp2/y_calculs/ the directory for the CIPSI calculations

# Go to the directory for a first CIPSI calculation
cd $PATH_CIPSI
qp_run fci ${DIR} > ${DIR}/${FILE}${IT}.fci
echo $(echo 0) "   "  $(grep "E               =" ${DIR}/${FILE}${IT}.fci | tail -1) >> ${DIR}/optimization.dat

# Optimization
for ((i=1 ; 100 - $i ; i++))
do
		cd ${PATH_OPT}
		qp_run orb_opt ${PATH_CIPSI}/${DIR} > ${PATH_CIPSI}/${DIR}/orb_trash${IT}${i}.dat
         
		echo $i

		cd ${PATH_CIPSI}
		qp_run diagonalize_h ${DIR} > ${DIR}/${FILE} > ${DIR}/${FILE}${IT}${i}.fci

		echo $(echo $i) "   "  $(grep "E               =" ${DIR}/${FILE}${IT}${i}.fci | tail -1) "   " $(grep "Gradient norm :" ${DIR}/orb_trash${IT}${i}.dat)  >> ${DIR}/optimization.dat
done


#       -p xeonv1 -N 1 -n 1 -c 16 --exclusive
#       -p xeonv2 -N 1 -n 1 -c 20 --exclusive
#       -p xeonv3 -N 1 -n 1 -c 24 --exclusive
#       -p xeonv4 -N 1 -n 1 -c 28 --exclusive
#       -p xeonv5 -N 1 -n 1 -c 32 --exclusive
#       -p xeonv6 -N 1 -n 1 -c 32 --exclusive

