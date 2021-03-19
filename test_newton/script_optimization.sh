#!/bin/bash

# for the cluster :
# Change the source !!!!!!!!!
#SBATCH -p xeonv3 -N 1 -n 1 -c 24 --exclusive
source /home/ydamour/qp2/quantum_package.rc
module load intel/2019.0
module load python/3.7.6-gcc-9.2.0
module load gcc/8.2.0

XYZ=copie_trust2  #CN
BASIS=cc_pvtz
EXTRA=DOCI

IT=_it_

#DIR=${XYZ}_${BASIS}_${EXTRA}.ezfio
DIR=${XYZ}.ezfio
#FILE=${XYZ}_${BASIS}_${EXTRA}
FILE=${XYZ}
PATH_CIPSI=../../../y_calculs
PATH_OPT=../plugins/qp_plugins_damour/test_newton

# path between 
# /qp2/plugin/qp_plugin_damour/test_newton the directory for optimizations
# and 
# /qp2/y_calculs/ the directory for the CIPSI calculations

echo ${DIR}

#Initialisation du trust region
qp_run init_nb_iteration ${PATH_CIPSI}/${DIR}

# Go to the directory for a first CIPSI calculation
cd $PATH_CIPSI
#qp_run fci ${DIR} > ${DIR}/${FILE}${IT}.fci  # normalement deja fait dans le repertoire
#echo $(echo 0) "   "  $(grep "E               =" ${DIR}/${FILE}${IT}.fci | tail -1) >> ${DIR}/optimization.dat

# Optimization
for ((i=1 ; 20 - $i ; i++))
do
		cd ${PATH_OPT}
		qp_run orb_opt ${PATH_CIPSI}/${DIR} > ${PATH_CIPSI}/${DIR}/orb_trash${IT}${i}.dat
         
		echo $i

		cd ${PATH_CIPSI}
		qp_run diagonalize_h ${DIR} > ${DIR}/${FILE} > ${DIR}/${FILE}${IT}${i}.diagonalize

		echo $i >> ${DIR}/iteration.dat
		grep "N_det =" ${DIR}/${FILE}${IT}${i}.diagonalize >> ${DIR}/nb_det.dat
	    grep "* Energy of state    1" ${DIR}/${FILE}${IT}${i}.diagonalize >> ${DIR}/energy.dat
	    grep "Gradient norm :" ${DIR}/orb_trash${IT}${i}.dat  >> ${DIR}/norm_grad.dat
done

paste ${DIR}/iteration.dat ${DIR}/energy.dat > ${DIR}/tmp_opt.dat
paste ${DIR}/tmp_opt.dat ${DIR}/norm_grad.dat > ${DIR}/tmp2_opt.dat
paste ${DIR}/tmp2_opt.dat ${DIR}/nb_det.dat >> ${DIR}/optimization.dat

cd ${PATH_OPT}

#       -p xeonv1 -N 1 -n 1 -c 16 --exclusive
#       -p xeonv2 -N 1 -n 1 -c 20 --exclusive
#       -p xeonv3 -N 1 -n 1 -c 24 --exclusive
#       -p xeonv4 -N 1 -n 1 -c 28 --exclusive
#       -p xeonv5 -N 1 -n 1 -c 32 --exclusive
#       -p xeonv6 -N 1 -n 1 -c 32 --exclusive

