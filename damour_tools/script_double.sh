#!/bin/bash
#SBATCH -p xeonv6 -N 1 -n 1 -c 32 --exclusive

set -e

export OMP_PROC_BIND=false
# Put the right path to source qp2
source $QP_ROOT/quantum_package.rc
#module list
echo "Hostename" $HOSTNAME
echo "OMP_PROC_BIND" $OMP_PROC_BIND
#module load intel/2019.3
#module load gcc/9.2.0
#module load python/3.7.6-gcc-9.2.0

MOL=h2 #tetrazine
BASIS=cc-pvqz  #6-31+g_star
CHARGE=0
MULTIPLICITY=1

STATE_0=0
STATE_1=1
METHOD=none #break_spatial_sym # break_spatial_sym, fb_loc, pm_loc
OPT_METHOD=diag          # diag, full
SELECTED_STATES=false    # if the states are close
N_MAX_STATES=4           # In this case, how many states ?
N_DET_MAX=2e7             

FILE1=${MOL}_hf
FILE2=${MOL}_no
FILE3=${MOL}_oo
DIR=${MOL}
EZFIO1=${FILE1}.ezfio
EZFIO2=${FILE2}.ezfio
EZFIO3=${FILE3}.ezfio

mkdir $DIR
cd $DIR
cp ../${MOL}.xyz .

qp_create_ezfio ${MOL}.xyz -b ${BASIS} -c ${CHARGE} -m ${MULTIPLICITY} -o ${EZFIO1}

echo ${EZFIO1}
qp set_file ${EZFIO1}
qp set determinants read_wf true

# SCF and frozen_core
qp run scf > ${FILE1}.scf.out
qp set_frozen_core > ${FILE1}.fc.out

# Save HFOs
tar zcf ${FILE1}_save_mos.ezfio.tar.gz ${EZFIO1}

# Selection of the states
if [[ ${SELECT_STATES} == true ]]
then
	qp set determinants n_states ${N_MAX_STATES}
	qp set determinants n_det_max 1e5
	qp run fci > ${FILE1}.pre_fci.out
	
	qp edit -s [${STATE_0},${STATE_1}]
fi

# cipsi w HFOs
qp set determinants n_states 2
qp set determinants n_det_max 1e7
qp run fci > ${FILE1}.fci.out

tar zcf ${FILE1}_save_cispi_res.tar.gz ${EZFIO1}

# NOs
cp -r ${EZFIO1} ${EZFIO2}
qp set_file ${EZFIO2}
qp run save_natorb > ${FILE2}.save_natorb.out
tar zcf ${FILE2}_save_mos.tar.gz ${EZFIO2}
cp -r ${EZFIO2} ${EZFIO3}

# cipsi w NOs
qp reset -d
qp set determinants n_states 2
qp set determinants n_det_max 1e7
qp run fci > ${FILE2}.nofci.out
tar zcf ${FILE2}_save_cipsi_res.tar.gz ${EZFIO2}

qp set_file ${EZFIO3}

if [[ ${METHOD} == break_spatial_sym ]] 
then
	qp set orbital_optimization security_mo_class false
	qp set orbital_optimization angle_pre_rot 1e-3
	qp set_mo_class -d [] -a [] -v []
	qp run break_spatial_sym > ${FILE3}.break_sym.out

elif [[ ${METHOD} == fb_loc ]]
then
        qp set orbital_optimization localization_method boys
        qp set orbital_optimization localization_max_nb_iter 1e4
	qp set orbital_optimization angle_pre_rot 1e-3
        qp set_mo_class -d [] -a [] -v []
        qp run localization > ${FILE3}.loc.out
        qp set_mo_class -c [] -a []

	tar zcf ${FILE3}_save_fb_lo_mos.tar.gz ${EZFIO3}


elif [[ ${METHOD} == pm_loc ]]
then
	qp set orbital_optimization localization_method pipek
	qp set orbital_optimization localization_max_nb_iter 1e4
	qp set orbital_optimization angle_pre_rot 1e-3
	qp set_mo_class -d [] -a [] -v []
	qp run localization > ${FILE3}.loc.out
	qp set_mo_class -c [] -a []

	tar zcf ${FILE3}_save_pm_lo_mos.tar.gz ${EZFIO3}
else
	echo ""
fi

# OOs
## Selecition of the states
if [[ ${SELECT_STATES} == true ]]
then
	qp set determinants n_states ${N_MAX_STATES}
	qp set determinants n_det_max 1e5
	qp run fci > ${FILE3}.pre_fci.out
	
	qp edit -s [${STATE_0},${STATE_1}]
fi

qp set determinants n_states 2

qp set orbital_optimization optimization_method ${OPT_METHOD}
qp set orbital_optimization normalized_st_av_weight true
qp set orbital_optimization start_from_wf true
qp set orbital_optimization n_det_start 5
qp set orbital_optimization n_det_max_opt 1e5
qp set orbital_optimization targeted_accuracy_cipsi 1e-5
qp set orbital_optimization optimization_max_nb_iter 10
qp set orbital_optimization thresh_opt_max_elem_grad 1e-4

qp run optimization > ${FILE3}.opt.out

## Save OOs
qp reset -d
tar zcf ${FILE3}_save_mos.tar.gz ${EZFIO3}

## Selection of the states
if [[ ${SELECT_STATES} == true ]]
then
	qp set determinants n_states ${N_MAX_STATES}
	qp set determinants n_det_max 1e4
	qp run fci > ${FILE3}.pre_opt_fci.out
	
	qp edit -s [${STATE_0},${STATE_1}]
fi

# cispi w OOs
qp set determinants n_states 2
qp set determinants n_det_max ${N_DET_MAX}
qp run fci > ${FILE3}.opt_fci.out

tar zcf ${FILE3}_save_cipsi_res.tar.gz ${EZFIO3}

# Data extraction and extrapolation of the fci/excitation energies
## with OOs
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extract_E_cipsi.py ${FILE2}.nofci.out
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extrapolation_fci.py ${FILE2}.nofci.out.dat > ${FILE2}.extrapolation_fci.dat
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/cipsi_error.py ${FILE2}.nofci.out.dat > ${FILE2}.cipsi_error.dat

## with OOs
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extract_E_cipsi.py ${FILE3}.opt_fci.out
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/extrapolation_fci.py ${FILE3}.opt_fci.out.dat > ${FILE3}.extrapolation_fci.dat
python3 $QP_ROOT/plugins/qp_plugins_damour/damour_tools/cipsi_error.py ${FILE3}.opt_fci.out.dat > ${FILE3}.cipsi_error.dat

